global.$ = require("jquery");
global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const point = (value, color = "black", opacity = 0) => ({ value, color, opacity });

function volume() {
    const data = new $3Dmol.VolumeData("", "cube");
    Object.assign(data, {
        size: { x: 2, y: 2, z: 2 }, unit: { x: 1, y: 1, z: 1 },
        origin: { x: 0, y: 0, z: 0 }, data: new Float32Array(8).fill(1)
    });
    return data;
}

describe("volumetric transfer functions", () => {
    beforeEach(() => {
        const interpolate = $3Dmol.GLVolumetricRender.interpolateArray;
        jest.spyOn($3Dmol.GLVolumetricRender, "interpolateArray").mockImplementation((data, count) => {
            // Overflow must fail the old implementation promptly instead of hanging Jest.
            if (!Number.isFinite(count) || count < 0 || count > 256) {
                throw new Error("Invalid transfer-function interpolation size.");
            }
            return interpolate(data, count);
        });
    });

    afterEach(() => jest.restoreAllMocks());

    test.each([
        ["missing specification", () => undefined],
        ["missing transfer function", () => ({})],
        ["null transfer function", () => ({ transferfn: null })],
        ["array-like object", () => ({ transferfn: { 0: point(0), 1: point(1), length: 2 } })],
        ["empty array", () => ({ transferfn: [] })],
        ["one point", () => ({ transferfn: [point(0)] })],
        ["sparse array", () => ({ transferfn: [point(0), , point(1)] })],
        ["null point", () => ({ transferfn: [point(0), null] })],
        ["missing value", () => ({ transferfn: [point(0), { color: "white", opacity: 1 }] })],
        ["empty value", () => ({ transferfn: [point(0), point("")] })],
        ["nonnumeric value", () => ({ transferfn: [point(0), point("invalid")] })],
        ["object value", () => ({ transferfn: [point(0), point({ toString: () => "1" })] })],
        ["symbol value", () => ({ transferfn: [point(0), point(Symbol())] })],
        ["NaN value", () => ({ transferfn: [point(0), point(NaN)] })],
        ["infinite value", () => ({ transferfn: [point(0), point(Infinity)] })],
        ["overflowing numeric string", () => ({ transferfn: [point(0), point("1e309")] })],
        ["constant values", () => ({ transferfn: [point(2), point("2.0")] })],
        ["overflowing range", () => ({ transferfn: [point(-1e308), point(1e308)] })],
        ["overflowing texture position", () => ({ transferfn: [point(0), point(1e308)] })]
    ])("%s warns and safely skips rendering", (_name, makeSpec) => {
        const data = volume();
        const before = Array.from(data.data);
        const warn = jest.spyOn(console, "warn").mockImplementation(() => {});
        const render = new $3Dmol.GLVolumetricRender(data, makeSpec());
        expect(warn).toHaveBeenCalledTimes(1);
        expect(render.hidden).toBe(true);
        expect(render.geo).toBeUndefined();
        expect(Array.from(data.data)).toEqual(before);
        const group = { add: jest.fn(), remove: jest.fn() };
        render.globj(group);
        render.removegl(group);
        expect(group.add).not.toHaveBeenCalled();
        expect(group.remove).not.toHaveBeenCalled();
    });

    test.each([
        ["sorted points", () => [point(0), point(1, "white", 1)]],
        ["unsorted points", () => [point(1, "white", 1), point(0)]],
        ["numeric strings", () => [point("0", "black", "0"), point("1", "white", "1")]],
        ["duplicate endpoints", () => [point(0), point(0), point(1, "white", 1), point(1, "white", 1)]]
    ])("%s produces the expected complete RGBA ramp", (_name, makePoints) => {
        const warn = jest.spyOn(console, "warn").mockImplementation(() => {});
        const render = new $3Dmol.GLVolumetricRender(volume(), { transferfn: makePoints() });
        expect(warn).not.toHaveBeenCalled();
        expect(render.hidden).toBe(false);
        expect(render.min).toBe(0);
        expect(render.max).toBe(1);
        // A black/transparent to white/opaque ramp contains every byte value once per channel.
        const expected = Uint8ClampedArray.from(Array.from({ length: 1024 }, (_, i) => Math.floor(i / 4)));
        expect(render.transferfunctionbuffer).toEqual(expected);
        const group = { add: jest.fn(), remove: jest.fn() };
        render.globj(group);
        expect(group.add).toHaveBeenCalledTimes(1);
        render.removegl(group);
        expect(group.remove).toHaveBeenCalledWith(group.add.mock.calls[0][0]);
    });
});
