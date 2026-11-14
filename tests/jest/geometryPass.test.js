global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function camera() {
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"), getAspect: () => 1,
            setSize() {}, setClearColorHex() {}, setViewport() {}, supportedExtensions: () => ({})
        };
    });
    try {
        return new $3Dmol.GLViewer(document.createElement("div"), { nomouse: true, autoDestroy: false }).camera;
    } finally {
        setup.mockRestore();
    }
}

function objectWithGroups(count, line = false) {
    let root;
    const parent = { add(value) { root = value; }, remove() {} };
    if (line) {
        const model = new $3Dmol.GLModel(0);
        model.addAtoms([
            { x: 0, y: 0, z: 0, elem: "C", bonds: [1], bondOrder: [1] },
            { x: 1, y: 0, z: 0, elem: "C", bonds: [0], bondOrder: [1] }
        ]);
        model.setStyle({}, { line: {} });
        model.globj(parent, {});
    } else {
        const shape = new $3Dmol.GLShape({});
        shape.addCustom({
            vertexArr: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }],
            faceArr: [0, 1, 2]
        });
        shape.globj(parent);
    }
    const object = root.children[0];
    const geometry = object.geometry, template = geometry.geometryGroups[0];
    // Small, independent real GeometryGroups avoid maximum-capacity fixture allocations.
    geometry.geometryGroups = Array.from({ length: count }, (_, i) => {
        const group = new template.constructor(i);
        for (const key of ["vertexArray", "colorArray", "normalArray", "faceArray", "lineArray"])
            group[key] = template[key] ? new template[key].constructor(template[key]) : null;
        for (const key of ["vertices", "faceidx", "lineidx", "useOffset"])
            group[key] = template[key];
        return group;
    });
    geometry.groups = count;
    return object;
}

function fixture(objects) {
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: {}, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer, bindings = new Map(), uploads = [];
    let constant = 100;
    const gl = new Proxy({
        isContextLost: () => false,
        createBuffer: jest.fn(() => ({})), deleteBuffer: jest.fn(),
        bindBuffer: (target, buffer) => bindings.set(target, buffer),
        bufferData: jest.fn((target, data) => uploads.push({ full: true, buffer: bindings.get(target), data })),
        bufferSubData: jest.fn((target, offset, data) => uploads.push({ full: false, buffer: bindings.get(target), data }))
    }, { get(target, key) {
        if (key in target) return target[key];
        return target[key] = /^[A-Z0-9_]+$/.test(key) ? constant++ : jest.fn();
    } });
    renderer._gl = gl;
    renderer._webglversion = 1;
    renderer.autoUpdateScene = false;
    renderer.sprites.render = jest.fn();
    renderer.renderBuffer = jest.fn(); // Keep actual record classification and draw-pass iteration.
    const scene = {
        __objectsAdded: [...objects], __objectsRemoved: [], __webglObjects: [], __webglSprites: [],
        __lights: [], fog: {}
    };
    renderer.initWebGLObjects(scene);
    uploads.length = 0;
    gl.bufferData.mockClear(); gl.bufferSubData.mockClear();
    return { renderer, gl, scene, uploads };
}

function watchGroupReads(geometry) {
    let reads = 0;
    geometry.geometryGroups = new Proxy(geometry.geometryGroups, { get(target, key, receiver) {
        if (/^\d+$/.test(String(key))) reads++;
        return Reflect.get(target, key, receiver);
    } });
    return () => reads;
}

function expectMatrices(object, viewCamera) {
    const expected = new $3Dmol.Matrix4().multiplyMatrices(viewCamera.matrixWorldInverse, object.matrixWorld);
    expect(Array.from(object._modelViewMatrix.elements)).toEqual(Array.from(expected.elements));
    const normal = new $3Dmol.Matrix3().getInverse(expected).transpose();
    expect(Array.from(object._normalMatrix.elements)).toEqual(Array.from(normal.elements));
}

describe("geometry-group preparation", () => {
    test.each([false, true])("clean groups are never rescanned or uploaded (line: %s)", line => {
        const object = objectWithGroups(100, line), env = fixture([object]);
        const reads = watchGroupReads(object.geometry);
        object.geometry.buffersNeedUpdate = true;
        env.renderer.initWebGLObjects(env.scene);
        expect(reads()).toBe(0);
        expect(env.uploads).toEqual([]);
        expect(object.geometry.buffersNeedUpdate).toBe(false);
        expect(env.scene.__webglObjects).toHaveLength(100);
    });

    test("coordinate-only updates visit each group once and preserve other buffers", () => {
        const object = objectWithGroups(100), { renderer, gl, scene, uploads } = fixture([object]);
        const reads = watchGroupReads(object.geometry);
        object.geometry.verticesNeedUpdate = object.geometry.normalsNeedUpdate = true;
        renderer.initWebGLObjects(scene);
        expect(reads()).toBe(100);
        expect(gl.bufferData).not.toHaveBeenCalled();
        expect(uploads).toHaveLength(200);
        expect(new Set(uploads.map(upload => upload.buffer)).size).toBe(200);
        for (const group of object.geometry.geometryGroups) {
            expect(uploads).toContainEqual({ full: false, buffer: group.__webglVertexBuffer, data: group.vertexArray });
            expect(uploads).toContainEqual({ full: false, buffer: group.__webglNormalBuffer, data: group.normalArray });
        }
        const before = reads();
        renderer.initWebGLObjects(scene);
        expect(reads()).toBe(before);
        expect(uploads).toHaveLength(200);
    });

    test.each(["colorsNeedUpdate", "elementsNeedUpdate"])("%s keeps the full-upload path", flag => {
        const object = objectWithGroups(4), { renderer, gl, scene, uploads } = fixture([object]);
        const reads = watchGroupReads(object.geometry);
        object.geometry.verticesNeedUpdate = object.geometry[flag] = true;
        renderer.initWebGLObjects(scene);
        expect(reads()).toBe(4);
        expect(gl.bufferSubData).not.toHaveBeenCalled();
        expect(uploads).toHaveLength(20); // Vertex, color, normal, face and wireframe buffers per group.
        expect(new Set(uploads.map(upload => upload.buffer)).size).toBe(20);
        expect(object.geometry[flag]).toBe(false);
    });

    test("disposed shared geometry is recreated and uploaded before its records are reused", () => {
        const first = objectWithGroups(3), second = first.clone();
        const { renderer, gl, scene, uploads } = fixture([first, second]);
        const oldBuffers = first.geometry.geometryGroups.map(group => group.__webglVertexBuffer);
        first.geometry.dispose();
        gl.createBuffer.mockClear();
        renderer.initWebGLObjects(scene);
        expect(gl.createBuffer).toHaveBeenCalledTimes(15);
        expect(uploads).toHaveLength(15);
        expect(scene.__webglObjects).toHaveLength(6);
        for (const group of first.geometry.geometryGroups) expect(oldBuffers).not.toContain(group.__webglVertexBuffer);
        renderer.initWebGLObjects(scene);
        expect(uploads).toHaveLength(15);
    });

    test("shared geometry uploads once but each instance gets its own matrices and draw records", () => {
        const first = objectWithGroups(4), second = first.clone();
        first.matrixWorld.setPosition(new $3Dmol.Vector3(2, 0, 0));
        second.matrixWorld.setPosition(new $3Dmol.Vector3(-3, 1, 0));
        const { renderer, scene, uploads } = fixture([first, second]);
        const matrices = jest.spyOn(renderer, "setupMatrices"), unroll = jest.spyOn(renderer, "unrollBufferMaterial");
        first.geometry.verticesNeedUpdate = true;
        const viewCamera = camera();
        renderer.render(scene, viewCamera);
        expect(uploads).toHaveLength(4);
        expect(matrices.mock.calls.map(call => call[0])).toEqual([first, second]);
        expect(unroll).toHaveBeenCalledTimes(8);
        expect(renderer.renderBuffer).toHaveBeenCalledTimes(8);
        expectMatrices(first, viewCamera); expectMatrices(second, viewCamera);
    });

    test("matrix reuse is limited to one render and recomputes after object or camera movement", () => {
        const object = objectWithGroups(100), { renderer, scene } = fixture([object]);
        const matrices = jest.spyOn(renderer, "setupMatrices"), viewCamera = camera();
        renderer.render(scene, viewCamera);
        expect(matrices).toHaveBeenCalledTimes(1);
        expect(renderer.renderBuffer).toHaveBeenCalledTimes(100);
        object.matrixWorld.setPosition(new $3Dmol.Vector3(5, 2, 1));
        viewCamera.position.z += 20;
        renderer.render(scene, viewCamera);
        expect(matrices).toHaveBeenCalledTimes(2);
        expectMatrices(object, viewCamera);
        const otherCamera = camera();
        otherCamera.position.x = 4;
        renderer.render(scene, otherCamera);
        expect(matrices).toHaveBeenCalledTimes(3);
        expectMatrices(object, otherCamera);
    });

    test("interleaved records and invisible objects retain correct preparation and drawing", () => {
        const a = objectWithGroups(2), b = objectWithGroups(1), hidden = objectWithGroups(1);
        hidden.visible = false;
        const { renderer, scene } = fixture([a, b, hidden]);
        const [a0, a1, b0, h0] = scene.__webglObjects;
        scene.__webglObjects = [a0, b0, h0, a1];
        const matrices = jest.spyOn(renderer, "setupMatrices"), unroll = jest.spyOn(renderer, "unrollBufferMaterial");
        const viewCamera = camera();
        renderer.render(scene, viewCamera);
        expect(matrices.mock.calls.map(call => call[0])).toEqual([a, b, a]);
        expect(unroll).toHaveBeenCalledTimes(3);
        expect(renderer.renderBuffer).toHaveBeenCalledTimes(3);
        expect(h0.render).toBe(false);
        expectMatrices(a, viewCamera); expectMatrices(b, viewCamera);
    });
});
