global.$ = require("jquery");
global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const transferfn = [
    { value: 0, color: "black", opacity: 0 },
    { value: 2, color: "white", opacity: 1 }
];
const matrix = (...values) => new $3Dmol.Matrix4(...values);

function volume(transform = null, unit = { x: 1, y: 1, z: 1 }) {
    const data = new $3Dmol.VolumeData("", "cube");
    Object.assign(data, {
        matrix: transform, unit, size: { x: 5, y: 5, z: 5 },
        data: new Float32Array(125).fill(1)
    });
    return data;
}

function boundedLookup(data) {
    const getIndex = data.getIndex.bind(data);
    let calls = 0;
    data.getIndex = (...args) => {
        // The original zero/negative increment must fail this test instead of hanging Jest.
        if (++calls > 1000) throw new Error("Volumetric selection did not finish sampling.");
        return getIndex(...args);
    };
    return () => calls;
}

function checkSpacing(data, expected) {
    boundedLookup(data);
    const center = new $3Dmol.Vector3(2, 2, 2);
    if (data.matrix) center.applyMatrix4(data.matrix);
    else {
        center.x *= data.unit.x;
        center.y *= data.unit.y;
        center.z *= data.unit.z;
    }
    const render = new $3Dmol.GLVolumetricRender(data, {
        transferfn, coords: [center], seldist: expected / 4
    });
    expect(render.hidden).toBe(false);
    expect(render.minunit).toBeCloseTo(expected, 6);
    expect(render.data.getVal(center.x, center.y, center.z)).toBe(1);

    const group = { add: jest.fn(), remove: jest.fn() };
    render.globj(group);
    expect(group.add).toHaveBeenCalledTimes(1);
    // The same positive spacing must reach the shader's sampling material.
    expect(group.add.mock.calls[0][0].children[0].material.unit).toBeCloseTo(expected, 6);
}

describe("volumetric sampling spacing", () => {
    const a = Math.SQRT1_2;
    test.each([
        ["axis-aligned grid", () => volume(), 1],
        ["90-degree rotation", () => volume(matrix(
            0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
        )), 1],
        ["45-degree rotation", () => volume(matrix(
            a, -a, 0, 0, a, a, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
        )), 1],
        ["translated anisotropic rotation", () => volume(matrix(
            0, -3, 0, 12, 2, 0, 0, -8, 0, 0, 4, 5, 0, 0, 0, 1
        )), 2],
        ["skew grid", () => volume(matrix(
            1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
        )), 1],
        ["negative axis spacing", () => volume(null, { x: -1, y: 2, z: 3 }), 1],
        ["matrix reflection", () => volume(matrix(
            -1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1
        )), 1]
    ])("%s samples with a positive distance", (_name, makeVolume, expected) => {
        checkSpacing(makeVolume(), expected);
    });

    test("rotated CUBE data samples with the parsed voxel-axis length", () => {
        const cube = [
            "rotated fixture", "angstrom grid", "0 0 0 0",
            "-5 0 1 0", "5 -1 0 0", "5 0 0 1",
            Array(125).fill("1").join(" ")
        ].join("\n");
        checkSpacing(new $3Dmol.VolumeData(cube, "cube"), 1);
    });

    test.each([
        ["zero axis", () => volume(null, { x: 0, y: 1, z: 1 })],
        ["NaN axis", () => volume(null, { x: NaN, y: 1, z: 1 })],
        ["infinite axis", () => volume(null, { x: Infinity, y: 1, z: 1 })],
        ["zero matrix axis", () => volume(matrix(
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
        ))],
        ["NaN matrix axis", () => {
            const transform = new $3Dmol.Matrix4();
            transform.elements[0] = NaN;
            return volume(transform);
        }],
        ["infinite matrix axis", () => {
            const transform = new $3Dmol.Matrix4();
            transform.elements[0] = Infinity;
            return volume(transform);
        }]
    ])("%s warns and skips rendering without modifying voxels", (_name, makeVolume) => {
        const data = makeVolume();
        const before = Array.from(data.data);
        const calls = boundedLookup(data);
        const warn = jest.spyOn(console, "warn").mockImplementation(() => {});
        try {
            const render = new $3Dmol.GLVolumetricRender(data, {
                transferfn, coords: [{ x: 0, y: 0, z: 0 }], seldist: 1
            });
            expect(warn).toHaveBeenCalledTimes(1);
            expect(render.hidden).toBe(true);
            expect(calls()).toBe(0);
            expect(Array.from(data.data)).toEqual(before);
            const group = { add: jest.fn(), remove: jest.fn() };
            render.globj(group);
            render.removegl(group);
            expect(group.add).not.toHaveBeenCalled();
        } finally {
            warn.mockRestore();
        }
    });
});
