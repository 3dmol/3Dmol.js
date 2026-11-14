global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const point = (value, color, opacity) => ({ value, color, opacity });
const stops = () => [point("1", "white", "1"), point("0", "black", "0")];

function render(transferfn) {
    const data = new $3Dmol.VolumeData("", "cube");
    Object.assign(data, {
        size: { x: 2, y: 2, z: 2 }, unit: { x: 1, y: 1, z: 1 },
        origin: { x: 0, y: 0, z: 0 }, data: new Float32Array(8).fill(1)
    });
    return new $3Dmol.GLVolumetricRender(data, { transferfn });
}

function expectRamp(result) {
    expect(result.hidden).toBe(false);
    expect(result.min).toBe(0);
    expect(result.max).toBe(1);
    // Every channel must contain the full black/transparent to white/opaque ramp.
    expect(result.transferfunctionbuffer).toEqual(
        Uint8ClampedArray.from({ length: 1024 }, (_, i) => Math.floor(i / 4))
    );
}

describe("transfer-function input ownership", () => {
    let warn;
    beforeEach(() => { warn = jest.spyOn(console, "warn").mockImplementation(() => {}); });
    afterEach(() => warn.mockRestore());

    test("normalization retains caller values, order, and point identities", () => {
        const transferfn = stops(), [high, low] = transferfn;
        expectRamp(render(transferfn));
        expect(transferfn[0]).toBe(high);
        expect(transferfn[1]).toBe(low);
        expect(transferfn).toEqual(stops());
        expect(warn).not.toHaveBeenCalled();
    });

    test.each([
        ["frozen string points", () => stops().map(Object.freeze)],
        ["frozen numeric points", () => [Object.freeze(point(0, "black", 0)), Object.freeze(point(1, "white", 1))]],
        ["frozen array and points", () => Object.freeze(stops().map(Object.freeze))],
        ["getter-only values", () => [0, 1].map(value => ({
            get value() { return String(value); }, color: value ? "white" : "black", opacity: value
        }))],
        ["inherited read-only fields", () => [0, 1].map(value => Object.create(Object.freeze(
            point(String(value), value ? "white" : "black", value)
        )))]]
    )("%s render without changing caller-owned properties", (_name, makePoints) => {
        const transferfn = makePoints();
        const before = transferfn.map(value => Object.getOwnPropertyDescriptors(value));
        expectRamp(render(transferfn));
        expect(transferfn.map(value => Object.getOwnPropertyDescriptors(value))).toEqual(before);
        expect(warn).not.toHaveBeenCalled();
    });

    test("duplicate shared points can be reused by multiple renderers", () => {
        const [high, low] = stops(), transferfn = [high, low, low, high];
        const first = render(transferfn), second = render(transferfn);
        expectRamp(first);
        expectRamp(second);
        expect(transferfn[0]).toBe(transferfn[3]);
        expect(transferfn[1]).toBe(transferfn[2]);
        expect([high, low]).toEqual(stops());
        expect(warn).not.toHaveBeenCalled();
    });

    test("normalization does not read unrelated point metadata", () => {
        const transferfn = stops();
        for (const value of transferfn) Object.defineProperty(value, "metadata", {
            enumerable: true, get() { throw new Error("Unrelated metadata was read."); }
        });
        expectRamp(render(transferfn));
        expect(transferfn.map(value => value.value)).toEqual(["1", "0"]);
        expect(warn).not.toHaveBeenCalled();
    });

    test.each([["2", "2.0"], ["-1e308", "1e308"]])(
        "invalid frozen range %s to %s warns without mutation or throwing", (low, high) => {
            const transferfn = Object.freeze([
                Object.freeze(point(low, "black", 0)), Object.freeze(point(high, "white", 1))
            ]);
            const result = render(transferfn);
            expect(result.hidden).toBe(true);
            expect(result.geo).toBeUndefined();
            expect(transferfn.map(value => value.value)).toEqual([low, high]);
            expect(warn).toHaveBeenCalledTimes(1);
        }
    );
});
