global.$ = require("jquery");
global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const center = { x: 0, y: 0, z: 0 };
const values = [1, 2, 3, 4, 5];
const selected = [1, 2, Infinity, Infinity, Infinity];
const transferfn = [
    { value: 0, color: "black", opacity: 0 },
    { value: 5, color: "white", opacity: 1 }
];

function volume(matrix = null) {
    const data = new $3Dmol.VolumeData("", "cube");
    Object.assign(data, {
        matrix, size: { x: 5, y: 1, z: 1 }, unit: { x: 1, y: 1, z: 1 },
        origin: center, data: new Float32Array(values)
    });
    return data;
}

function render(spec, data = volume(), viewer) {
    const original = data.data;
    const result = new $3Dmol.GLVolumetricRender(data, { transferfn, ...spec }, viewer);
    expect(data.data).toBe(original);
    expect(Array.from(data.data)).toEqual(values);
    return result;
}

describe("volumetric selection distance", () => {
    test("omitted distance matches the documented default of 2", () => {
        const spec = { coords: [center] };
        const implicit = render(spec);
        const explicit = render({ ...spec, seldist: 2 });
        expect(implicit.hidden).toBe(false);
        expect(Array.from(implicit.data.data)).toEqual(selected);
        expect(implicit.data.data).toEqual(explicit.data.data);
        expect(spec.seldist).toBeUndefined();
    });

    test("an explicitly undefined distance uses the default", () => {
        expect(Array.from(render({ coords: [center], seldist: undefined }).data.data)).toEqual(selected);
    });

    test("an atom selection also uses the default", () => {
        const selection = { chain: "A" };
        const viewer = { selectedAtoms: jest.fn(() => [center]) };
        const result = render({ selection }, volume(), viewer);
        expect(viewer.selectedAtoms).toHaveBeenCalledWith(selection);
        expect(Array.from(result.data.data)).toEqual(selected);
    });

    test.each([
        ["empty coordinates", { coords: [] }, undefined],
        ["empty atom selection", { selection: {} }, { selectedAtoms: () => [] }]
    ])("%s masks the whole volume", (_name, spec, viewer) => {
        expect(Array.from(render(spec, volume(), viewer).data.data)).toEqual(Array(5).fill(Infinity));
    });

    test.each([
        [0, [Infinity, Infinity, Infinity, Infinity, Infinity]],
        [3, [1, 2, 3, Infinity, Infinity]]
    ])("explicit distance %s is preserved", (seldist, expected) => {
        expect(Array.from(render({ coords: [center], seldist }).data.data)).toEqual(expected);
    });

    test.each([null, -1, NaN, Infinity])("invalid distance %s still warns and skips rendering", seldist => {
        const warn = jest.spyOn(console, "warn").mockImplementation(() => {});
        try {
            const result = render({ coords: [center], seldist });
            expect(warn).toHaveBeenCalledTimes(1);
            expect(result.hidden).toBe(true);
            const group = { add: jest.fn(), remove: jest.fn() };
            result.globj(group);
            result.removegl(group);
            expect(group.add).not.toHaveBeenCalled();
        } finally {
            warn.mockRestore();
        }
    });

    test("without coordinates or a selection the volume remains unmasked", () => {
        const data = volume();
        expect(render({}, data).data).toBe(data);
        expect(render({ seldist: 2 }, data).data).toBe(data);
    });

    test("coordinates take precedence over an atom selection", () => {
        const viewer = { selectedAtoms: jest.fn(() => []) };
        const result = render({ coords: [center], selection: {} }, volume(), viewer);
        expect(viewer.selectedAtoms).not.toHaveBeenCalled();
        expect(Array.from(result.data.data)).toEqual(selected);
    });

    test("the default distance is measured in world coordinates on a transformed grid", () => {
        const matrix = new $3Dmol.Matrix4(
            0, -1, 0, 10, 1, 0, 0, -5, 0, 0, 1, 2, 0, 0, 0, 1
        );
        const result = render({ coords: [{ x: 10, y: -5, z: 2 }] }, volume(matrix));
        expect(Array.from(result.data.data)).toEqual(selected);
    });
});
