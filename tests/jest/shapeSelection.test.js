global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function viewer() {
    // Exercise actual model selection and shape geometry without a WebGL context.
    const result = Object.create($3Dmol.GLViewer.prototype);
    result.models = [];
    result.shapes = [];
    const model = new $3Dmol.GLModel(0);
    model.addAtoms([
        { x: 2, y: 4, z: 6, chain: "A" },
        { x: 8, y: 10, z: 12, chain: "A" },
        { x: 3, y: 5, z: 7, chain: "A" },
        { x: 50, y: 60, z: 70, chain: "B" }
    ].map(atom => ({ ...atom, elem: "C", bonds: [] })));
    result.models.push(model);
    return result;
}

const selection = () => ({ x: [2, 8, 3], y: [4, 10, 5], z: [6, 12, 7] });
const center = { x: 5, y: 7, z: 9 };
const endpoint = { x: 14, y: 16, z: 18 };
const origin = { x: 0, y: 0, z: 0 };
const xyz = point => ({ x: point.x, y: point.y, z: point.z });

function expectGeometry(shape) {
    let vertices = 0;
    for (const group of shape.geo.geometryGroups) {
        expect(group.vertexArray.slice(0, group.vertices * 3).every(Number.isFinite)).toBe(true);
        vertices += group.vertices;
    }
    expect(vertices).toBeGreaterThan(0);
    expect(Number.isFinite(shape.boundingSphere.radius)).toBe(true);
}

describe("shape positions from atom selections", () => {
    test("coordinate arrays select the bounding-box center, not the arithmetic mean", () => {
        const input = selection();
        expect(viewer().getSelectionCenter(input)).toEqual(center);
        expect(input).toEqual(selection());
    });

    test("mixed coordinate filters and additional selection fields are respected", () => {
        const v = viewer();
        expect(v.getSelectionCenter({ x: [2, 8], y: 4, z: 6 })).toEqual({ x: 2, y: 4, z: 6 });
        expect(v.getSelectionCenter({ ...selection(), predicate: atom => atom.x < 4 }))
            .toEqual({ x: 2.5, y: 4.5, z: 6.5 });
        expect(v.getSelectionCenter({ chain: "A" })).toEqual(center);
    });

    test("empty selections return origin while an explicit empty filter selects all", () => {
        const v = viewer();
        expect(v.getSelectionCenter({ x: [], y: [4], z: [6] })).toEqual(origin);
        expect(v.getSelectionCenter({})).toEqual({ x: 26, y: 32, z: 38 });
        v.models = [];
        expect(v.getSelectionCenter(selection())).toEqual(origin);
    });

    test.each([
        ["plain coordinates", () => ({ x: 0, y: -2, z: 4 })],
        ["Vector3", () => new $3Dmol.Vector3(2, 4, 6)],
        ["atom coordinates", () => ({ x: 2, y: 4, z: 6, elem: "C", chain: "B" })],
        ["null-prototype coordinates", () => Object.assign(Object.create(null), { x: 2, y: 4, z: 6 })]
    ])("%s remain explicit positions with their identity intact", (_name, makePoint) => {
        const point = makePoint();
        expect(viewer().getSelectionCenter(point)).toBe(point);
    });

    test("sphere selection centers generate finite geometry", () => {
        const shape = viewer().addSphere({ center: selection(), radius: 0.5 });
        expect(xyz(shape.intersectionShape.sphere[0].center)).toEqual(center);
        expectGeometry(shape);
    });

    test.each(["center", "corner"])("box %s accepts coordinate-array selections", field => {
        const shape = viewer().addBox({ [field]: selection(), dimensions: { w: 2, h: 4, d: 6 } });
        expect(xyz(shape.boundingSphere.center)).toEqual(field === "center" ? center : { x: 6, y: 9, z: 12 });
        expectGeometry(shape);
    });

    test("cylinder selection endpoints generate finite geometry", () => {
        const shape = viewer().addCylinder({ start: selection(), end: endpoint, radius: 0.2 });
        const cylinder = shape.intersectionShape.cylinder[0];
        expect(xyz(cylinder.c1)).toEqual(center);
        expect(xyz(cylinder.c2)).toEqual(endpoint);
        expectGeometry(shape);
    });

    test("line selection endpoints are written as numeric vertices", () => {
        const shape = viewer().addLine({ start: selection(), end: endpoint });
        expect(Array.from(shape.geo.geometryGroups[0].vertexArray.slice(0, 6))).toEqual([5, 7, 9, 14, 16, 18]);
        expectGeometry(shape);
    });

    test("arrow selections resolve both endpoints before geometry generation", () => {
        const spec = { start: selection(), end: endpoint, radius: 0.2 };
        const shape = viewer().addArrow(spec);
        expect(xyz(spec.start)).toEqual(center);
        expect(xyz(spec.end)).toEqual(endpoint);
        expectGeometry(shape);
    });

    test.each([undefined, null])("missing sphere and line positions (%s) retain shape defaults", position => {
        const v = viewer();
        const sphere = v.addSphere({ center: position, radius: 0.5 });
        const line = v.addLine({ start: position, end: position });
        expect(xyz(sphere.intersectionShape.sphere[0].center)).toEqual(origin);
        expect(Array.from(line.geo.geometryGroups[0].vertexArray.slice(0, 6))).toEqual([0, 0, 0, 3, 0, 0]);
        expectGeometry(sphere);
        expectGeometry(line);
    });

    test("an arrow with direction and length can omit its endpoint", () => {
        const spec = { start: selection(), dir: new $3Dmol.Vector3(0, 1, 0), length: 3, radius: 0.2 };
        const shape = viewer().addArrow(spec);
        expect(xyz(spec.start)).toEqual(center);
        expect(xyz(spec.end)).toEqual({ x: 5, y: 10, z: 9 });
        expectGeometry(shape);
    });

    test("omitted arrow and cylinder coordinates use their existing defaults", () => {
        const v = viewer(), spec = { radius: 0.2 };
        const arrow = v.addArrow(spec);
        expect(xyz(spec.start)).toEqual(origin);
        expect(xyz(spec.end)).toEqual({ x: 3, y: 0, z: 0 });
        const cylinder = v.addCylinder({ end: endpoint, radius: 0.2 });
        expect(xyz(cylinder.intersectionShape.cylinder[0].c1)).toEqual(origin);
        expect(xyz(cylinder.intersectionShape.cylinder[0].c2)).toEqual(endpoint);
        expectGeometry(arrow);
        expectGeometry(cylinder);
    });
});
