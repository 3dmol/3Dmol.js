global.URL.createObjectURL = function () {};
const $3Dmol = require('../../build/3Dmol.js');

const roundtrip = state => JSON.parse(JSON.stringify(state));
const xyz = '2\nfirst\nC -1 0 0\nO 1 0 0\n2\nsecond\nC -1 1 0\nO 1 1 0\n';
const dx = 'object 1 class gridpositions counts 2 2 2\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 2 class gridconnections counts 2 2 2\nobject 3 class array type double rank 0 items 8 data follows\n0 1 1 2 1 2 2 3\n';
let viewers;

beforeEach(() => {
    viewers = [];
    jest.spyOn($3Dmol.GLViewer.prototype, 'setupRenderer').mockImplementation(function () {
        const canvas = document.createElement('canvas');
        this.renderer = {
            domElement: canvas, getCanvas: () => canvas, getAspect: () => 1,
            setSize() {}, setClearColorHex: jest.fn(), setViewport() {}, supportsVolumetric: () => true,
            supportedExtensions: () => ({}), render: jest.fn(), dispose: jest.fn()
        };
    });
    $3Dmol.setSyncSurface(true);
});
afterEach(() => {
    viewers.forEach(viewer => { viewer.destroy(); viewer.container?.remove(); });
    document.body.replaceChildren();
    jest.restoreAllMocks();
});
function viewer() {
    const element = document.createElement('div'); document.body.append(element);
    const result = new $3Dmol.GLViewer(element, { nomouse: true }); viewers.push(result); return result;
}

function parts(state) { return roundtrip({ shapes: state.shapes, labels: state.labels, surfaces: state.surfaces, view: state.view }); }

test('JSON roundtrip preserves compound/custom shapes, labels, transformed masked volumes and surface geometry', async () => {
    const source = viewer(), target = viewer();
    source.addModel('2\natoms\nC -1 0 0\nO 1 0 0\n', 'xyz');
    source.setStyle({}, { sphere: { radius: 0.4 } });
    const compound = source.addShape({ color: 'red' });
    compound.addSphere({ center: { x: 0, y: 2, z: 0 }, radius: 0.5 });
    compound.addLine({ start: { x: -1, y: 2, z: 0 }, end: { x: 1, y: 2, z: 0 }, color: 'red' });
    compound.updateStyle({ color: 'blue', opacity: 0.4, frame: 1 });
    source.addCustom({ vertexArr: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }],
        faceArr: [0, 1, 2], color: ['red', 'green', 'blue'] });
    const volume = new $3Dmol.VolumeData(dx, 'dx');
    volume.matrix = new $3Dmol.Matrix4(); volume.matrix.elements[12] = 2;
    source.addIsosurface(volume, { isoval: 1.5, color: 'purple' });
    const density = source.addVolumetricRender(volume, { coords: [], subsamples: 3,
        transferfn: [{ value: 0, color: 'red', opacity: 0 }, { value: 3, color: 'blue', opacity: 0.8 }] });
    const annotation = source.addLine({ start: { x: -1, y: -2, z: 0 }, end: { x: 1, y: -2, z: 0 }, opacity: 0.25 });
    const label = source.addLabel('Classroom', { position: { x: 3, y: 2, z: 0 }, fontColor: new $3Dmol.Color('red') });
    label.hide();
    const surfaceId = await source.addSurface($3Dmol.SurfaceType.VDW, { color: 'green', opacity: 0.4 }, { elem: 'C' }, {});
    source.setBackgroundColor('navy', 0.3); source.setProjection('orthographic'); source.setSlab(-5, 10);
    source.render();
    const snapshot = roundtrip(source.getInternalState());
    expect(snapshot.shapes).toHaveLength(5);
    expect(snapshot.surfaces[0].parts[0].geometry.groups[0].vertices).toBeGreaterThan(0);
    target.setInternalState(snapshot);
    expect(parts(target.getInternalState())).toEqual(parts(snapshot));
    const restoredVolume = target.getShape(density.shapePosition);
    expect(restoredVolume.data.data).toBeInstanceOf(Float32Array);
    expect(Array.from(restoredVolume.data.data).every(value => value === Infinity)).toBe(true);
    expect(restoredVolume.texmatrix).toBeInstanceOf($3Dmol.Matrix4);
    expect(target.getShape(compound.shapePosition).opacity).toBe(0.4);
    expect(target.selectedAtoms({})).toHaveLength(2);
    const line = target.getShape(annotation.shapePosition).renderedShapeObj.children[0];
    expect(line.material.opacity).toBe(0.25);
    expect(line.material.transparent).toBe(true);
    expect(line.material.wireframe).toBe(true);
    restoredVolume.hidden = true; target.render();
    expect(restoredVolume.renderedShapeObj).toBeNull();
    restoredVolume.hidden = false; target.render();
    expect(restoredVolume.renderedShapeObj).not.toBeNull();
    source.setSurfaceMaterialStyle(surfaceId, { color: 'yellow', opacity: 0.8, wireframe: true });
    target.setInternalState({ surfaces: roundtrip(source.getInternalState()).surfaces }, { partial: true });
    expect(target.getSurface(surfaceId)[0].mat.opacity).toBe(0.8);
    expect(target.getSurface(surfaceId)[0].mat.wireframe).toBe(true);
    expect(parts(target.getInternalState())).toEqual(parts(source.getInternalState()));
}, 20000);

test('compact model styles restore gradient objects and materialize custom color functions', () => {
    const source = viewer(), target = viewer();
    source.addModel('2\natoms\nC -1 0 0\nO 1 0 0\n', 'xyz');
    source.setStyle({}, { sphere: { colorscheme: { prop: 'x', gradient: new $3Dmol.Gradient.CustomLinear(-1, 1, ['red', 'blue']) } },
        stick: { colorfunc: atom => atom.elem === 'C' ? 0x00ff00 : 0xffff00 } });
    target.setInternalState(roundtrip(source.getInternalState({ compact: true })));
    const atoms = target.selectedAtoms({});
    expect(atoms[0].style.sphere.colorscheme.gradient.valueToHex(-1)).toBe(0xff0000);
    expect(atoms[0].style.sphere.colorscheme.gradient.valueToHex(1)).toBe(0x0000ff);
    expect(atoms.map(atom => atom.style.stick.color)).toEqual([0x00ff00, 0xffff00]);
    expect(atoms[0].style.stick.colorfunc).toBeUndefined();
});

test('compact appearance updates leave coordinates, frames, model identity and unrelated geometry intact', async () => {
    const source = viewer(), target = viewer();
    source.addModelsAsFrames(xyz, 'xyz');
    source.setStyle({}, { stick: { color: 'red' } });
    source.addSphere({ center: { x: 4, y: 0, z: 0 }, radius: 1 });
    const first = roundtrip(source.getInternalState({ compact: true }));
    target.setInternalState(first);
    const model = target.getModel(), shape = target.getShape(0), frames = model.getInternalState().frames;
    source.setStyle({}, { sphere: { color: 'blue' } });
    source.getModel().hide();
    const next = roundtrip(source.getInternalState({ compact: true }));
    expect(next.models).toEqual(first.models);
    expect(next.modelStyles).not.toEqual(first.modelStyles);
    target.setInternalState({ modelStyles: next.modelStyles }, { partial: true });
    expect(target.getModel()).toBe(model);
    expect(model.getInternalState().frames).toBe(frames);
    expect(target.getShape(0)).toBe(shape);
    expect(model.isHidden()).toBe(true);
    expect(model.selectedAtoms({})[0].style).toEqual({ sphere: { color: 'blue' } });
    await source.setFrame(1);
    expect(roundtrip(source.getInternalState({ compact: true })).models).toEqual(next.models);
});

test('partial overlays remove replaced objects without clearing models or evaluating callback strings', () => {
    const source = viewer(), target = viewer();
    source.addSphere({ center: { x: 0, y: 0, z: 0 }, radius: 1 });
    target.addModel('1\natom\nC 0 0 0\n', 'xyz');
    const model = target.getModel();
    const snapshot = roundtrip(source.getInternalState());
    snapshot.shapes[0].style.callback = 'function () { throw new Error("must not execute"); }';
    target.setInternalState({ shapes: snapshot.shapes }, { partial: true });
    expect(target.getModel()).toBe(model);
    expect(target.getShape(0).callback).toBeUndefined();
    expect(target.getShape(0).clickable).toBe(false);
    target.setInternalState({ shapes: [], labels: [], surfaces: [] }, { partial: true });
    expect(target.getShape(0)).toBeUndefined();
    expect(target.getModel()).toBe(model);
    target.setInternalState({ models: [] });
    expect(target.selectedAtoms({})).toHaveLength(0);
});

test('remote snapshots retain the timeline without requesting frames or wrapping overlays early', async () => {
    const source = viewer(), target = viewer();
    source.addModel('1\natom\nC 3 0 0\n', 'xyz');
    const internal = source.getModel().getInternalState();
    internal.frames = [internal.atoms];
    Object.assign(internal.frames, { url: 'https://trajectory.example', path: 'classroom', numFrames: 20 });
    source.getModel().setInternalState(internal);
    source.addSphere({ center: { x: 0, y: 2, z: 0 }, radius: 0.3, frame: 5 });
    const snapshot = roundtrip(source.getInternalState({ compact: true }));
    const request = jest.spyOn(global, 'fetch');
    target.setInternalState(snapshot);
    await target.setFrame(11); target.render();
    expect(target.getNumFrames()).toBe(20);
    expect(target.getFrame()).toBe(11);
    expect(target.getShape(0).renderedShapeObj).toBeNull();
    expect(target.selectedAtoms({})[0].x).toBe(3);
    expect(request).not.toHaveBeenCalled();
    expect(roundtrip(target.getInternalState({ compact: true })).models[0].snapshotFrameCount).toBe(20);
});

test('lightweight render subscriptions defer snapshots and unsubscribe cleanly', () => {
    const source = viewer(), listener = jest.fn();
    const capture = jest.spyOn(source, 'getInternalState');
    const unsubscribe = source.subscribeStateChanges(listener);
    source.render();
    expect(listener).toHaveBeenCalledTimes(1);
    expect(capture).not.toHaveBeenCalled();
    unsubscribe(); source.render();
    expect(listener).toHaveBeenCalledTimes(1);
});
