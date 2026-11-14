global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const G = $3Dmol.Gradient;

function model() {
    const result = new $3Dmol.GLModel(0);
    result.addAtoms([10, 20, 30].map((score, x) => ({
        x, y: 0, z: 0, elem: "C", bonds: [], properties: { score }
    })));
    return result;
}

const colors = m => m.selectedAtoms({}).map(atom => atom.color);

describe("gradient validation", () => {
    let warn;
    beforeEach(() => { warn = jest.spyOn(console, "warn").mockImplementation(() => {}); });
    afterEach(() => warn.mockRestore());

    test.each([
        ["null", null], ["undefined", undefined], ["empty description", {}],
        ["unknown name", { gradient: "missing" }], ["numeric name", { gradient: 7 }],
        ["constructor", { gradient: "constructor" }], ["toString", { gradient: "toString" }],
        ["__proto__", { gradient: "__proto__" }],
        ["noncallable methods", { range: 5, valueToHex: "bad" }]
    ])("%s warns and produces a usable fallback", (_name, input) => {
        const fallback = G.getGradient(input);
        expect(fallback).toBeInstanceOf(G.RWB);
        expect(fallback.range()).toEqual([-1, 1]);
        expect([-1, 0, 1].map(value => fallback.valueToHex(value))).toEqual([0xff0000, 0xffffff, 0x0000ff]);
        expect(warn).toHaveBeenCalledTimes(1);
    });

    test("valid builtin aliases and custom linear descriptions retain their colors", () => {
        for (const name of ["rwb", "RWB", "roygb", "ROYGB", "sinebow", "linear"]) {
            const spec = { gradient: name, min: -3, max: 5 };
            if (name === "linear") spec.colors = ["red", "yellow", "blue"];
            const actual = G.getGradient(spec);
            const expected = name === "linear" ? new G.CustomLinear(-3, 5, spec.colors) : new G.builtinGradients[name](-3, 5);
            expect(actual.range()).toEqual([-3, 5]);
            const samples = [-4, -3, -1, 1, 3, 5, 6];
            expect(samples.map(value => actual.valueToHex(value))).toEqual(samples.map(value => expected.valueToHex(value)));
        }
        const shorthand = G.getGradient({ gradient: "linear_red_yellow_blue", min: 0, max: 10 });
        expect([0, 5, 10].map(value => shorthand.valueToHex(value))).toEqual([0xff0000, 0xffff00, 0x0000ff]);
        expect(warn).not.toHaveBeenCalled();
    });

    test("structural custom gradients retain identity and control property coloring", () => {
        const custom = {
            gradient: "custom", range: () => [10, 20],
            valueToHex: jest.fn(() => 0x123456)
        };
        expect(G.getGradient(custom)).toBe(custom);
        const actual = $3Dmol.getColorFromStyle({ elem: "C", properties: { score: 15 } }, {
            colorscheme: { prop: "score", gradient: custom }
        });
        expect(actual.getHex()).toBe(0x123456);
        expect(custom.valueToHex).toHaveBeenCalledWith(15, [10, 20]);
        expect(warn).not.toHaveBeenCalled();
    });

    test.each([{}, { gradient: "constructor" }])("invalid label background %j still renders", input => {
        const label = new $3Dmol.Label("Gradient", { backgroundGradient: input });
        try {
            label.setContext();
            expect(label.canvas.width).toBeGreaterThan(0);
            expect(label.canvas.height).toBeGreaterThan(0);
            expect(label.sprite.material.map.image).toBe(label.canvas);
            expect(warn).toHaveBeenCalledTimes(1);
        } finally {
            label.dispose();
        }
    });

    test("a structural label gradient supplies its own range and colors", () => {
        const samples = [];
        const custom = { range: () => [0, 100], valueToHex: value => { samples.push(value); return 0x123456; } };
        const label = new $3Dmol.Label("Custom", { backgroundGradient: custom });
        try {
            label.setContext();
            expect(samples).toHaveLength(11);
            expect(samples[0]).toBe(0);
            expect(samples.at(-1)).toBeCloseTo(100);
            expect(warn).not.toHaveBeenCalled();
        } finally {
            label.dispose();
        }
    });

    test("volume styles normalize invalid frozen descriptions without mutating them", () => {
        const original = Object.freeze({ gradient: "missing", min: 0, max: 50 });
        const style = { volscheme: original };
        $3Dmol.adjustVolumeStyle(style);
        expect(style.volscheme).toBeInstanceOf(G.RWB);
        expect(style.volscheme.valueToHex(0)).toBe(0xffffff);
        expect(original).toEqual({ gradient: "missing", min: 0, max: 50 });
        expect(warn).toHaveBeenCalledTimes(1);
    });

    test.each(["missing", null, 2])("invalid property gradient %s uses the fallback", gradient => {
        const actual = $3Dmol.getColorFromStyle({ elem: "C", color: 0, properties: { score: 0 } }, {
            colorscheme: { prop: "score", gradient }
        });
        expect(actual.getHex()).toBe(0xffffff);
        expect(warn).toHaveBeenCalledTimes(1);
    });

    test.each(["missing", "constructor", "toString", "__proto__"])("unknown model gradient %s retains automatic data ranges", name => {
        const m = model();
        m.setColorByProperty({}, "score", name);
        expect(colors(m)).toEqual([0xff0000, 0xffffff, 0x0000ff]);
        expect(warn).toHaveBeenCalledTimes(1);
    });

    test("model fallback respects an explicit range", () => {
        const m = model();
        m.setColorByProperty({}, "score", "missing", [0, 100]);
        const expected = new G.RWB(0, 100);
        expect(colors(m)).toEqual([10, 20, 30].map(value => expected.valueToHex(value)));
        expect(warn).toHaveBeenCalledTimes(1);
    });

    test("valid model names still infer the range and custom objects accept an explicit range", () => {
        const m = model();
        m.setColorByProperty({}, "score", "rwb");
        expect(colors(m)).toEqual([0xff0000, 0xffffff, 0x0000ff]);
        const custom = { range: () => [10, 30], valueToHex: jest.fn(() => 0xabcdef) };
        m.setColorByProperty({}, "score", custom, [0, 100]);
        expect(colors(m)).toEqual([0xabcdef, 0xabcdef, 0xabcdef]);
        expect(custom.valueToHex.mock.calls).toEqual([[10, [0, 100]], [20, [0, 100]], [30, [0, 100]]]);
        expect(warn).not.toHaveBeenCalled();
    });
});
