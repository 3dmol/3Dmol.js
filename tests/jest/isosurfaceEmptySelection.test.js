global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function volume() {
    return {
        size: { x: 3, y: 3, z: 3 },
        origin: new $3Dmol.Vector3(0, 0, 0),
        unit: new $3Dmol.Vector3(1, 1, 1),
        data: Float32Array.from({ length: 27 }, (_, i) => Math.floor(i / 9) - 1)
    };
}

function viewer() {
    const result = Object.create($3Dmol.GLViewer.prototype);
    result.models = [];
    result.shapes = [];
    const model = new $3Dmol.GLModel(0);
    model.addAtoms([{ x: 1, y: 1, z: 1, elem: "C", bonds: [], chain: "A" }]);
    result.models.push(model);
    return result;
}

const vertexCount = shape => shape.geo.geometryGroups.reduce((sum, group) => sum + group.vertices, 0);

describe("isosurfaces restricted to empty regions", () => {
    test.each([
        ["explicit coordinates", () => ({ coords: [] })],
        ["legacy selectedRegion", () => ({ selectedRegion: [] })],
        ["unmatched atom selection", () => ({ selection: { chain: "missing" } })]
    ])("%s completes with no geometry and finite bounds", (_name, makeSpec) => {
        const callback = jest.fn();
        const v = viewer();
        const shape = v.addIsosurface(volume(), { ...makeSpec(), opacity: 0.4 }, callback);
        expect(v.shapes).toEqual([shape]);
        expect(vertexCount(shape)).toBe(0);
        expect(shape.boundingSphere.center).toEqual(new $3Dmol.Vector3(0, 0, 0));
        expect(shape.boundingSphere.radius).toBe(0);
        expect(shape.opacity).toBe(0.4);
        expect(callback).toHaveBeenCalledTimes(1);
        expect(callback).toHaveBeenCalledWith();
    });

    test("an empty addition preserves existing geometry and bounds while updating style", () => {
        const shape = new $3Dmol.GLShape({});
        shape.addSphere({ center: { x: 4, y: 5, z: 6 }, radius: 2 });
        const vertices = vertexCount(shape);
        const center = shape.boundingSphere.center.clone();
        const radius = shape.boundingSphere.radius;
        const callback = jest.fn(() => expect(shape.opacity).toBe(0.4));
        shape.addIsosurface(volume(), { coords: [], opacity: 0.4 }, callback);
        expect(vertexCount(shape)).toBe(vertices);
        expect(shape.boundingSphere.center).toEqual(center);
        expect(shape.boundingSphere.radius).toBe(radius);
        expect(callback).toHaveBeenCalledTimes(1);
    });

    test("explicit empty coordinates take precedence over a matching atom selection", () => {
        const shape = viewer().addIsosurface(volume(), { coords: [], selection: { chain: "A" } });
        expect(vertexCount(shape)).toBe(0);
    });

    test.each([
        ["unrestricted", {}],
        ["matching atom selection", { selection: { chain: "A" }, seldist: 5 }]
    ])("%s still creates an isosurface", (_name, spec) => {
        const callback = jest.fn();
        const shape = viewer().addIsosurface(volume(), { ...spec, isoval: 0.5, smoothness: 0 }, callback);
        expect(vertexCount(shape)).toBeGreaterThan(0);
        expect(Number.isFinite(shape.boundingSphere.radius)).toBe(true);
        expect(callback).toHaveBeenCalledTimes(1);
    });
});
