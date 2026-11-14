global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const vec = (x, y, z) => new $3Dmol.Vector3(x, y, z);
let raycaster;

beforeAll(() => {
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = { domElement: document.createElement("canvas"), getAspect: () => 1,
            setSize() {}, setClearColorHex() {} };
    });
    try {
        raycaster = new $3Dmol.GLViewer(document.createElement("div"), { nomouse: true, autoDestroy: false }).raycaster;
    } finally { setup.mockRestore(); }
});
afterEach(() => jest.restoreAllMocks());

function triangle(z = 5, x = 0) {
    return new $3Dmol.Triangle(vec(x - 1, -1, z), vec(x + 1, -1, z), vec(x, 1, z));
}
function object(triangles = [], other = {}) {
    return { intersectionShape: { triangle: triangles, sphere: [], cylinder: [], line: [], ...other } };
}
function pick(objects, matrix = new $3Dmol.Matrix4(), origin = vec(0, 0, 0)) {
    raycaster.set(origin, vec(0, 0, 1));
    return raycaster.intersectObjects({ matrixWorld: matrix }, objects);
}

// These checks exercise scratch reuse both within a shape and between separate queries.
test("repeated hits and misses preserve distances, ordering, and source vertices", () => {
    const far = object([triangle(9)]), near = object([triangle(3)]), miss = object([triangle(4, 10)]);
    const before = JSON.stringify([far, near, miss]);
    for (let i = 0; i < 4; i++) {
        const hits = pick([far, miss, near]);
        expect(hits.map(hit => hit.clickable)).toEqual([near, far]);
        expect(hits.map(hit => hit.distance)).toEqual([3, 9]);
        expect(pick([near, far], undefined, vec(20, 0, 0))).toEqual([]);
    }
    expect(JSON.stringify([far, near, miss])).toBe(before);
});

test("winding, behind-camera triangles, and degenerate triangles retain their culling", () => {
    const reversed = triangle();
    [reversed.b, reversed.c] = [reversed.c, reversed.b];
    const degenerate = new $3Dmol.Triangle(vec(-1, 0, 5), vec(0, 0, 5), vec(1, 0, 5));
    const valid = object([triangle()]);
    const hits = pick([object([reversed, triangle(-5), degenerate]), valid]);
    expect(hits).toEqual([{ clickable: valid, distance: 5 }]);
    // A second culling pass must not read a normal left by the valid triangle.
    expect(pick([object([reversed, degenerate])])).toEqual([]);
});

test.each([[0, -1], [-1, -1], [0, 0]])("accepts boundary/interior ray at (%s, %s)", (x, y) => {
    const target = object([triangle()]);
    expect(pick([target], undefined, vec(x, y, 0))).toEqual([{ clickable: target, distance: 5 }]);
});

test("transformed triangles coexist with cylinders, lines, and spheres using the same scratch vectors", () => {
    const matrix = new $3Dmol.Matrix4();
    Object.assign(matrix.elements, { 0: 2, 5: 3, 10: 4, 12: 5, 13: -2, 14: 1 });
    const mixed = object([triangle(3)], {
        sphere: [new $3Dmol.Sphere(vec(0, 0, 1), 0.25)],
        cylinder: [new $3Dmol.Cylinder(vec(-1, 0, 4), vec(1, 0, 4), 0.25)],
        line: [vec(-1, 0, 2), vec(1, 0, 2)]
    });
    const last = object([triangle(5)]), before = JSON.stringify([mixed, last]);
    const hits = pick([mixed, last], matrix, vec(5, -2, 0));
    expect(hits.map(hit => hit.distance)).toEqual([4, 9, 13, 16, 21]);
    expect(hits.map(hit => hit.clickable)).toEqual([mixed, mixed, mixed, mixed, last]);
    expect(JSON.stringify([mixed, last])).toBe(before);
    matrix.elements[0] = -2;
    expect(pick([last], matrix, vec(5, -2, 0))).toEqual([]);
});

test("testing ten thousand triangles allocates no cloned normal vectors, including successful hits", () => {
    const triangles = Array.from({ length: 10000 }, () => triangle(5, 10));
    const target = object(triangles), group = { matrixWorld: new $3Dmol.Matrix4() };
    raycaster.set(vec(0, 0, 0), vec(0, 0, 1));
    const clone = jest.spyOn($3Dmol.Vector3.prototype, "clone");
    let misses, hits, clones;
    try {
        misses = raycaster.intersectObjects(group, [target]);
        triangles.push(triangle(7), triangle(2));
        hits = raycaster.intersectObjects(group, [target]);
        clones = clone.mock.calls.length;
    } finally { clone.mockRestore(); }
    expect(misses).toEqual([]);
    expect(hits.map(hit => hit.distance)).toEqual([2, 7]);
    expect(hits.every(hit => hit.clickable === target)).toBe(true);
    expect(clones).toBe(0);
});

test("the public triangle normal remains a fresh unit vector independent of picking scratch space", () => {
    const face = triangle(), before = JSON.stringify(face);
    const first = face.getNormal(), second = face.getNormal();
    expect(first).not.toBe(second);
    expect(first.length()).toBe(1);
    expect([first.x, first.y, first.z]).toEqual([0, 0, -1]);
    first.set(8, 9, 10);
    pick([object([face])]);
    expect([second.x, second.y, second.z]).toEqual([0, 0, -1]);
    expect([first.x, first.y, first.z]).toEqual([8, 9, 10]);
    expect(JSON.stringify(face)).toBe(before);
});
