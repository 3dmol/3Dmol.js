global.$ = require("jquery");
global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function volume() {
    const data = new $3Dmol.VolumeData("", "cube");
    data.size = { x: 2, y: 3, z: 4 };
    data.data = Float32Array.from({ length: 24 }, (_, i) => i + 0.25);
    return data;
}

function checkRoundTrip(data) {
    for (let index = 0; index < data.data.length; index++) {
        const point = data.getCoordinates(index);
        expect(data.getIndex(point.x, point.y, point.z)).toBe(index);
        expect(data.getVal(point.x, point.y, point.z)).toBe(data.data[index]);
    }
}

describe("volume grid coordinates", () => {
    test("unequal grid dimensions use z-fastest row and plane strides", () => {
        const data = volume();
        const samples = [
            [0, { x: 0, y: 0, z: 0 }],
            [1, { x: 0, y: 0, z: 1 }],
            [3, { x: 0, y: 0, z: 3 }],
            [4, { x: 0, y: 1, z: 0 }],
            [11, { x: 0, y: 2, z: 3 }],
            [12, { x: 1, y: 0, z: 0 }],
            [23, { x: 1, y: 2, z: 3 }]
        ];
        for (const [index, expected] of samples) {
            expect(data.getCoordinates(index)).toEqual(expected);
        }
    });

    test("every coordinate round trips to its original index and value", () => {
        checkRoundTrip(volume());
    });

    test("origin and signed anisotropic units locate grid samples", () => {
        const data = volume();
        data.origin = { x: -3, y: 7, z: 2 };
        data.unit = { x: 2, y: -0.5, z: 1.25 };
        expect(data.getCoordinates(0)).toEqual({ x: -3, y: 7, z: 2 });
        expect(data.getCoordinates(4)).toEqual({ x: -3, y: 6.5, z: 2 });
        expect(data.getCoordinates(23)).toEqual({ x: -1, y: 6, z: 5.75 });
        checkRoundTrip(data);
    });

    test("matrix rotation, reflection and translation take precedence over origin and unit", () => {
        const data = volume();
        data.matrix = new $3Dmol.Matrix4(
            0, -2, 0, 10,
            3, 0, 0, -4,
            0, 0, -0.5, 7,
            0, 0, 0, 1
        );
        data.origin = { x: 100, y: 200, z: 300 };
        data.unit = { x: 8, y: -9, z: 12 };
        expect(data.getCoordinates(0)).toEqual({ x: 10, y: -4, z: 7 });
        expect(data.getCoordinates(12)).toEqual({ x: 10, y: -1, z: 7 });
        expect(data.getCoordinates(23)).toEqual({ x: 6, y: -1, z: 5.5 });
        checkRoundTrip(data);
    });

    test("parsed skewed CUBE coordinates include the basis vectors and origin", () => {
        const cube = [
            "coordinate fixture", "skewed angstrom grid", "0 1 2 -3",
            "-2 2 0.25 -0.5", "3 0.5 3 0.25", "4 -0.25 0.5 4",
            Array.from({ length: 24 }, (_, i) => i + 0.25).join(" ")
        ].join("\n");
        const data = new $3Dmol.VolumeData(cube, "cube");
        expect(data.data.length).toBe(24);
        expect(data.matrix).toBeTruthy();
        expect(data.getCoordinates(0)).toEqual({ x: 1, y: 2, z: -3 });
        expect(data.getCoordinates(1)).toEqual({ x: 0.75, y: 2.5, z: 1 });
        expect(data.getCoordinates(4)).toEqual({ x: 1.5, y: 5, z: -2.75 });
        expect(data.getCoordinates(12)).toEqual({ x: 3, y: 2.25, z: -3.5 });
        expect(data.getCoordinates(23)).toEqual({ x: 3.25, y: 9.75, z: 9 });
        checkRoundTrip(data);
    });

    test("each result is an independent plain coordinate object", () => {
        const data = volume();
        const first = data.getCoordinates(0);
        const second = data.getCoordinates(0);
        expect(Object.getPrototypeOf(first)).toBe(Object.prototype);
        expect(Object.keys(first).sort()).toEqual(["x", "y", "z"]);
        expect(first).not.toBe(second);
        first.x = 100;
        expect(second).toEqual({ x: 0, y: 0, z: 0 });
        expect(data.getCoordinates(0)).toEqual(second);
    });

    test("singleton axes preserve the remaining grid strides", () => {
        const data = volume();
        data.size = { x: 1, y: 3, z: 4 };
        expect(data.getCoordinates(11)).toEqual({ x: 0, y: 2, z: 3 });
        data.size = { x: 3, y: 1, z: 4 };
        expect(data.getCoordinates(9)).toEqual({ x: 2, y: 0, z: 1 });
        data.size = { x: 3, y: 4, z: 1 };
        expect(data.getCoordinates(9)).toEqual({ x: 2, y: 1, z: 0 });
    });

    test("large valid indices are not truncated to 32-bit integers", () => {
        const data = volume();
        // Coordinate lookup needs only metadata, not a multi-gigabyte voxel array.
        data.size = { x: 50000, y: 50000, z: 2 };
        expect(data.getCoordinates(4999999999)).toEqual({ x: 49999, y: 49999, z: 1 });
    });
});
