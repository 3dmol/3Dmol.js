/** @jest-environment node */
// Exercise surface geometry without a WebGL context or native canvas renderer.
jest.mock('canvas', () => ({}));
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!doctype html><html><body></body></html>');
global.window = dom.window;
global.document = dom.window.document;
global.URL.createObjectURL = jest.fn(() => 'blob:surface-test');
dom.window.URL.createObjectURL = jest.fn(() => 'blob:surface-test');
Object.defineProperty(global, '$3Dmol', {
    configurable: true,
    get: () => dom.window.$3Dmol,
    set: value => { dom.window.$3Dmol = value; }
});
const $3Dmol = require('../../build/3Dmol.js');
const originalSyncSurface = $3Dmol.syncSurface;
beforeAll(() => $3Dmol.setSyncSurface(true));
afterAll(() => {
    $3Dmol.setSyncSurface(originalSyncSurface);
    dom.window.close();
    delete global.window;
    delete global.document;
    delete global.$3Dmol;
});

function model(id, atoms) {
    const result = new $3Dmol.GLModel(id);
    result.addAtoms(atoms.map(atom => ({ ...atom, bonds: [] })));
    return result;
}

function viewer(models) {
    return Object.assign(Object.create($3Dmol.GLViewer.prototype), {
        models, surfaces: {}, render: jest.fn()
    });
}

async function geometry(instance, shown, calculated) {
    const id = await instance.addSurface($3Dmol.SurfaceType.VDW,
        { color: 'blue' }, shown, calculated);
    return Array.from(instance.getSurface(id)).flatMap(part => part.geo.geometryGroups.map(group => ({
        vertices: Array.from(group.vertexArray.slice(0, group.vertices * 3)),
        faces: Array.from(group.faceArray.slice(0, group.faceidx)),
        normals: Array.from(group.normalArray.slice(0, group.vertices * 3)),
        colors: Array.from(group.colorArray.slice(0, group.vertices * 3))
    })));
}

describe('surface calculation selections', () => {
    const carbon = { elem: 'C', x: 0, y: 0, z: 0 };
    const oxygen = { elem: 'O', x: 2.4, y: 0, z: 0 };

    test('surface chunks distinguish matching atom indices in different models', () => {
        const instance = viewer([model(0, [carbon]), model(1, [oxygen])]);
        // addSurface copies the display and calculation atoms independently.
        const calculated = instance.selectedAtoms({}).map(atom => ({ ...atom }));
        const shown = instance.selectedAtoms({ model: 0 }).map(atom => ({ ...atom }));
        expect(calculated.map(atom => atom.index)).toEqual([0, 0]);
        const chunks = instance.carveUpExtent($3Dmol.getExtent(shown, true), calculated, shown);
        expect(chunks).toHaveLength(1);
        expect(chunks[0].atoms).toEqual([0, 1]);
        expect(chunks[0].toshow).toEqual([0]);
    });

    test('single-model geometry is identical with implicit or explicit calculation atoms', async () => {
        const instance = viewer([model(0, [carbon, oxygen])]);
        const implicit = await geometry(instance, {});
        const explicit = await geometry(instance, {}, {});
        expect(implicit.some(group => group.faces.length > 0)).toBe(true);
        expect(explicit).toEqual(implicit);
    });

    test('displaying one model against all atoms matches the equivalent single-model surface', async () => {
        const combined = viewer([model(0, [carbon, oxygen])]);
        const separate = viewer([model(0, [carbon]), model(1, [oxygen])]);
        const expected = await geometry(combined, { index: 0 }, {});
        const actual = await geometry(separate, { model: 0 }, {});
        expect(expected.some(group => group.faces.length > 0)).toBe(true);
        expect(actual).toEqual(expected);
    });
});
