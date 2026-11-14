global.URL.createObjectURL = function () {};
const M = require("../../build/3Dmol.js");
const viewers = [];
function viewer(orthographic = false) {
    const setup = jest.spyOn(M.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        const canvas = document.createElement("canvas");
        this.renderer = { domElement: canvas, getCanvas: () => canvas, getAspect: () => 1,
            getXRatio: () => 1, getYRatio: () => 1, setSize() {}, setClearColorHex() {}, dispose() {} };
    });
    let v;
    try { v = new M.GLViewer(document.createElement("div"), { nomouse: true, autoDestroy: false, orthographic }); }
    finally { setup.mockRestore(); }
    viewers.push(v);
    v.WIDTH = v.HEIGHT = 400;
    v.show = jest.fn(() => { v.setSlabAndFog(); v.scene.updateMatrixWorld(); v.camera.updateMatrixWorld(); });
    v.show();
    return v;
}
afterEach(() => { viewers.splice(0).forEach(v => v.destroy()); jest.restoreAllMocks(); jest.useRealTimers(); });
function oldProjection(v, point) {
    const t = new M.Vector3(point.x, point.y, point.z).applyMatrix4(v.modelGroup.matrixWorld);
    v.projector.projectVector(t, v.camera);
    return { x: 200 * (t.x + 1), y: -200 * (t.y - 1) };
}
function oldDistance(v, screen, model) {
    const p = new M.Vector3(model.x, model.y, model.z).applyMatrix4(v.modelGroup.matrixWorld), world = p.clone();
    v.projector.projectVector(p, v.camera);
    const t = new M.Vector3(screen.x * 2 / 400 - 1, screen.y * 2 / -400 + 1, p.z);
    v.projector.unprojectVector(t, v.camera);
    return t.distanceTo(world);
}
test.each([false, true])("batched projections reuse one vector and preserve independent results (ortho: %s)", ortho => {
    const v = viewer(ortho);
    v.modelGroup.position.set(2, -3, 4); v.rotate(37, "y");
    const coords = Array.from({ length: 100 }, (_, i) => ({ x: i / 3, y: Math.sin(i), z: i % 7 }));
    const expected = coords.map(point => oldProjection(v, point));
    const seen = new Set(), original = M.Vector3.prototype.applyMatrix4;
    const spy = jest.spyOn(M.Vector3.prototype, "applyMatrix4").mockImplementation(function (matrix) { seen.add(this); return original.call(this, matrix); });
    const result = v.modelToScreen(coords);
    spy.mockRestore();
    expect(result).toEqual(expected); expect(seen.size).toBe(1);
    expect(new Set(result).size).toBe(coords.length);
    result[0].x = 12345;
    expect(v.modelToScreen(coords[1])).toEqual(expected[1]);
    expect(v.modelToScreen([])).toEqual([]);
    const unusual = [{ x: NaN, y: -0 }, { x: Infinity, y: 2, z: 3 }];
    expect(v.modelToScreen(unusual)).toEqual(unusual.map(point => oldProjection(v, point)));
});
test("nested projection calls retain their own working vector", () => {
    const v = viewer();
    const point = { get x() { expect(v.modelToScreen({ x: 8, y: 9, z: 1 })).toEqual(oldProjection(v, { x: 8, y: 9, z: 1 })); return 2; }, y: 3, z: 4 };
    expect(v.modelToScreen([{ x: 1, y: 1, z: 1 }, point])).toEqual([oldProjection(v, { x: 1, y: 1, z: 1 }), oldProjection(v, { x: 2, y: 3, z: 4 })]);
});
test.each([false, true])("screen distance preserves original arithmetic and fresh public offsets (ortho: %s)", ortho => {
    const v = viewer(ortho); v.rotate(17, "x");
    const points = [{ x: 1, y: 2, z: 3 }, { x: NaN, y: -0, z: undefined }, { x: Infinity, y: 0, z: 2 }];
    for (const model of points) for (const screen of [{ x: 210, y: 187 }, { x: NaN, y: -0 }]) {
        const expected = oldDistance(v, screen, model), clone = jest.spyOn(M.Vector3.prototype, "clone");
        expect(Object.is(v.screenToModelDistance(screen, model), expected)).toBe(true);
        expect(clone).not.toHaveBeenCalled(); clone.mockRestore();
    }
    const first = v.screenOffsetToModel(10, 20), saved = first.clone();
    const second = v.screenOffsetToModel(30, 40);
    expect(first).not.toBe(second); expect(first).toEqual(saved);
});
test("rotation reuses its quaternion without retaining it in animation targets or callbacks", () => {
    const v = viewer(), scratch = v.dq;
    for (const axis of [{ x: 1, y: 2, z: 3 }, { x: NaN, y: -0, z: 1 }]) {
        const current = v.rotationGroup.quaternion.clone(), s = Math.sin(Math.PI * 17 / 180 / 2), c = Math.cos(Math.PI * 17 / 180 / 2);
        const expected = current.multiply(new M.Quaternion(axis.x * s, axis.y * s, axis.z * s, c).normalize());
        v.rotate(17, axis);
        expect(v.rotationGroup.quaternion).toEqual(expected); expect(v.dq).toBe(scratch);
    }
    const animate = jest.spyOn(v, "animateMotion").mockImplementation(() => {});
    v.rotate(25, "z", 60, true);
    const target = animate.mock.calls[0][4], saved = target.clone();
    expect(target).not.toBe(scratch);
    v.rotate(19, "x"); expect(target).toEqual(saved);
    let nested = false;
    v.show.mockImplementation(() => { if (!nested) { nested = true; v.rotate(11, "y"); } });
    const expected = v.rotationGroup.quaternion.clone().multiply(new M.Quaternion(0, 0, Math.sin(Math.PI * 7 / 180 / 2), Math.cos(Math.PI * 7 / 180 / 2)).normalize()).multiply(new M.Quaternion(0, Math.sin(Math.PI * 11 / 180 / 2), 0, Math.cos(Math.PI * 11 / 180 / 2)).normalize());
    v.rotate(7, "z");
    expect(v.rotationGroup.quaternion).toEqual(expected);
});
test("surface normals match per-face vector calculations without retaining face temporaries", () => {
    const vertices = [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }, { x: 0, y: 0, z: 2 }];
    const faces = [0, 1, 2, 0, 3, 1, 1, 3, 2, 0, 0, 1];
    const seen = new Set(), original = M.Vector3.prototype.cross;
    const cross = jest.spyOn(M.Vector3.prototype, "cross").mockImplementation(function (v) { seen.add(this); return original.call(this, v); });
    const mesh = M.GLViewer.generateSurfaceMesh([], { vertices, faces }, { vertexColors: true });
    cross.mockRestore();
    const geo = mesh.geometry.geometryGroups[0], expected = new Float32Array(geo.normalArray.length);
    for (let i = 0; i < faces.length; i += 3) {
        const [a, b, c] = faces.slice(i, i + 3).map(index => new M.Vector3(...geo.vertexArray.slice(index * 3, index * 3 + 3)));
        c.subVectors(c, b); a.subVectors(a, b); c.cross(a).normalize();
        for (let j = 0; j < 3; j++) {
            const offset = faces[i + j] * 3;
            expected[offset] += c.x; expected[offset + 1] += c.y; expected[offset + 2] += c.z;
        }
    }
    expect(geo.normalArray).toEqual(expected); expect(seen.size).toBe(1);
    mesh.geometry.dispose();
});
