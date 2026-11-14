global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function fixture() {
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: {}, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer, uniforms = {};
    let constant = 100;
    const gl = new Proxy({
        uniformMatrix4fv: jest.fn((location, transpose, elements) => { uniforms[location] = Array.from(elements); })
    }, { get(target, key) {
        if (!(key in target)) target[key] = /^[A-Z0-9_]+$/.test(key) ? constant++ : jest.fn();
        return target[key];
    } });
    renderer._gl = gl;
    renderer._fb = renderer._screenshader = renderer._blurshader = {};
    renderer._AOshader = { uniforms: { projectionMatrix: "projection", projinv: "inverse" } };
    renderer.renderObjects = jest.fn();
    return { renderer, gl, uniforms };
}
function record(x = 0, properties = {}) {
    const material = { clone: jest.fn(() => ({ opacity: 1, dispose: jest.fn() })) };
    const matrix = new $3Dmol.Matrix4();
    matrix.elements[12] = x;
    return { render: true, opaque: material, object: { _modelViewMatrix: matrix }, ...properties };
}
const scene = records => ({ __lights: [], fog: {}, __webglObjects: records });
const camera = () => ({ projectionMatrix: new $3Dmol.Matrix4() });
const flush = async () => { await Promise.resolve(); await Promise.resolve(); };
afterEach(async () => { await flush(); jest.restoreAllMocks(); });

test("empty material cleanup clears stale AO flags without creating filtered lists", () => {
    const { renderer } = fixture();
    const base = { dispose: jest.fn(), map: { dispose: jest.fn() } };
    const entry = { object: { material: base }, hasAO: true };
    const filter = jest.spyOn(Array.prototype, "filter");
    let filters;
    try {
        for (let i = 0; i < 100; i++) renderer.releaseBufferMaterials(entry);
        filters = filter.mock.calls.length;
    } finally {
        filter.mockRestore();
    }
    expect(filters).toBe(0);
    expect(entry.hasAO).toBe(false);
    expect([entry.opaqueDepth, entry.transparentDepth, entry.opaqueShaded]).toEqual([null, null, null]);
    expect(base.dispose).not.toHaveBeenCalled();
    expect(base.map.dispose).not.toHaveBeenCalled();
});

test("owned pass materials still retire once at the next microtask while base resources survive", async () => {
    const { renderer } = fixture();
    const base = { dispose: jest.fn(), map: { dispose: jest.fn() } };
    const shared = { dispose: jest.fn() }, transparent = { dispose: jest.fn() };
    const entry = { object: { material: base }, opaqueDepth: shared, opaqueShaded: shared, transparentDepth: transparent, hasAO: true };
    renderer.releaseBufferMaterials(entry);
    renderer.releaseBufferMaterials(entry);
    expect(shared.dispose).not.toHaveBeenCalled();
    expect(transparent.dispose).not.toHaveBeenCalled();
    expect(entry.hasAO).toBe(false);
    expect([entry.opaqueDepth, entry.transparentDepth, entry.opaqueShaded]).toEqual([null, null, null]);
    await flush();
    expect(shared.dispose).toHaveBeenCalledTimes(1);
    expect(transparent.dispose).toHaveBeenCalledTimes(1);
    expect(base.dispose).not.toHaveBeenCalled();
    expect(base.map.dispose).not.toHaveBeenCalled();
});

test("AO uses the first matching object and prepares depth materials for every later match", () => {
    const { renderer, uniforms } = fixture();
    const hidden = record(100, { render: false }), nonmatching = record(200, { opaque: null });
    const first = record(3), second = record(7), third = record(11);
    first.opaque.wireframe = true;
    const records = [hidden, nonmatching, first, second, third];
    const view = camera(), world = scene(records);
    renderer.setShading(world, view, "opaque");
    expect(uniforms.projection).toEqual(Array.from(first.object._modelViewMatrix.elements));
    expect(hidden.opaque.clone).not.toHaveBeenCalled();
    expect(first.opaque.clone).not.toHaveBeenCalled();
    for (const entry of [second, third]) {
        expect(entry.opaque.clone).toHaveBeenCalledTimes(1);
        expect(entry.opaqueDepth.opacity).toBe(0);
    }
    expect(renderer.renderObjects).toHaveBeenCalledWith(records, true, "opaqueDepth", view, world.__lights, world.fog, false);
    renderer.setShading(world, view, "opaque");
    for (const entry of [second, third]) expect(entry.opaque.clone).toHaveBeenCalledTimes(1);
});

test("AO with no matching objects performs no pass or matrix upload", () => {
    const { renderer, gl } = fixture();
    renderer.setShading(scene([]), camera(), "opaque");
    renderer.setShading(scene([record(1, { render: false }), record(2, { opaque: null })]), camera(), "opaque");
    expect(renderer.renderObjects).not.toHaveBeenCalled();
    expect(gl.framebufferTexture2D).not.toHaveBeenCalled();
    expect(gl.uniformMatrix4fv).not.toHaveBeenCalled();
    expect(gl.drawArrays).not.toHaveBeenCalled();
});

test("AO reuses matrix storage while updating values for different cameras and first objects", () => {
    const { renderer, uniforms } = fixture();
    const forward = renderer._fullProjModelMatrix, inverse = renderer._fullProjModelMatrixInv;
    const forwardData = forward.elements, inverseData = inverse.elements;
    const first = record(3), next = record(-9), view = camera(), world = scene([first, next]);
    for (let pass = 0; pass < 3; pass++) {
        first.render = pass === 0;
        view.projectionMatrix.elements[0] = pass + 2;
        view.projectionMatrix.elements[5] = pass + 3;
        view.projectionMatrix.elements[9] = pass / 4;
        const object = pass === 0 ? first : next;
        const expected = new $3Dmol.Matrix4().multiplyMatrices(view.projectionMatrix, object.object._modelViewMatrix);
        const expectedInverse = new $3Dmol.Matrix4().getInverse(expected);
        renderer.setShading(world, view, "opaque");
        expect(renderer._fullProjModelMatrix).toBe(forward);
        expect(renderer._fullProjModelMatrixInv).toBe(inverse);
        expect(forward.elements).toBe(forwardData);
        expect(inverse.elements).toBe(inverseData);
        expect(uniforms.projection).toEqual(Array.from(expected.elements));
        expect(uniforms.inverse).toEqual(Array.from(expectedInverse.elements));
    }
});

test("a singular AO projection resets the reused inverse rather than retaining its previous values", () => {
    const { renderer, uniforms } = fixture();
    const view = camera(), world = scene([record(7)]);
    renderer.setShading(world, view, "opaque");
    expect(uniforms.inverse[12]).toBe(-7);
    view.projectionMatrix.elements.fill(0);
    jest.spyOn(console, "warn").mockImplementation(() => {});
    renderer.setShading(world, view, "opaque");
    expect(uniforms.inverse).toEqual(Array.from(new $3Dmol.Matrix4().elements));
    expect(uniforms.projection).toEqual(Array(16).fill(0));
});
