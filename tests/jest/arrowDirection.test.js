global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function arrow(spec) {
    const shape = new $3Dmol.GLShape({});
    shape.addArrow(spec);
    return shape;
}

const xyz = point => ({ x: point.x, y: point.y, z: point.z });
const vertices = shape => shape.geo.geometryGroups.flatMap(group =>
    Array.from(group.vertexArray.slice(0, group.vertices * 3)));

describe("arrow direction coordinates", () => {
    test.each([3, -2, 0])("plain directions match Vector3 geometry for length %s", length => {
        const dir = Object.freeze({ x: 2, y: -1, z: 0 });
        const start = Object.freeze({ x: 1, y: 2, z: 3 });
        const spec = { start, dir, length, end: { x: 40, y: 50, z: 60 } };
        const expected = { x: 1 + 2 * length, y: 2 - length, z: 3 };
        const shape = arrow(spec);
        const vectorShape = arrow({ start, dir: new $3Dmol.Vector3(2, -1, 0), length });
        expect(xyz(spec.end)).toEqual(expected);
        expect(spec.dir).toBe(dir);
        expect(dir).toEqual({ x: 2, y: -1, z: 0 });
        const actual = vertices(shape);
        expect(actual.length).toBeGreaterThan(0);
        expect(actual.every(Number.isFinite)).toBe(true);
        expect(actual).toEqual(vertices(vectorShape));
        expect(Number.isFinite(shape.boundingSphere.radius)).toBe(true);
    });

    test("existing Vector3 inputs retain their components and magnitude", () => {
        const dir = Object.freeze(new $3Dmol.Vector3(0, 2, -1));
        const spec = { dir, length: 3 };
        arrow(spec);
        expect(xyz(spec.start)).toEqual({ x: 0, y: 0, z: 0 });
        expect(xyz(spec.end)).toEqual({ x: 0, y: 6, z: -3 });
        expect(spec.dir).toBe(dir);
        expect(xyz(dir)).toEqual({ x: 0, y: 2, z: -1 });
    });

    test("plain directions can omit both start and end", () => {
        const spec = { dir: { x: 0, y: 1, z: 0 }, length: 4 };
        arrow(spec);
        expect(xyz(spec.start)).toEqual({ x: 0, y: 0, z: 0 });
        expect(xyz(spec.end)).toEqual({ x: 0, y: 4, z: 0 });
    });

    test.each([
        ["missing length", { dir: { x: 0, y: 1, z: 0 } }],
        ["missing direction", { length: 4 }]
    ])("%s preserves explicit endpoints and the default endpoint", (_name, options) => {
        const explicit = { ...options, end: { x: 7, y: 8, z: 9 } };
        arrow(explicit);
        expect(xyz(explicit.end)).toEqual({ x: 7, y: 8, z: 9 });
        const defaults = { ...options };
        arrow(defaults);
        expect(xyz(defaults.end)).toEqual({ x: 3, y: 0, z: 0 });
    });

    test("viewer resolves an atom selection before applying the plain direction", () => {
        const viewer = Object.create($3Dmol.GLViewer.prototype);
        viewer.models = [];
        viewer.shapes = [];
        const model = new $3Dmol.GLModel(0);
        model.addAtoms([
            { x: 2, y: 4, z: 6, elem: "C", bonds: [], chain: "A" },
            { x: 4, y: 8, z: 10, elem: "C", bonds: [], chain: "A" }
        ]);
        viewer.models.push(model);
        const spec = { start: { chain: "A" }, dir: { x: 0, y: 2, z: -1 }, length: 3 };
        const shape = viewer.addArrow(spec);
        expect(xyz(spec.start)).toEqual({ x: 3, y: 6, z: 8 });
        expect(xyz(spec.end)).toEqual({ x: 3, y: 12, z: 5 });
        expect(viewer.shapes).toEqual([shape]);
        expect(vertices(shape).length).toBeGreaterThan(0);
        expect(vertices(shape).every(Number.isFinite)).toBe(true);
    });
});
