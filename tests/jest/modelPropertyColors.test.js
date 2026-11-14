global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function model(atoms) {
    const result = new $3Dmol.GLModel(0);
    result.addAtoms(atoms.map((atom, x) => ({ x, y: 0, z: 0, elem: "C", bonds: [], ...atom })));
    return result;
}

const colors = m => m.selectedAtoms({}).map(atom => atom.color);

describe("model property colors", () => {
    test.each([undefined, {}])("built-in atom fields work with properties %j", properties => {
        const m = model([10, 20, 30].map(b => ({ b, properties })));
        m.setColorByProperty({}, "b", "rwb");
        expect(colors(m)).toEqual([0xff0000, 0xffffff, 0x0000ff]);
    });

    test("nested properties override direct fields and missing nested values fall back", () => {
        const m = model([
            { b: 100, properties: { b: 10 } },
            { b: 20, properties: {} },
            { b: 30 }
        ]);
        m.setColorByProperty({}, "b", "rwb");
        expect(colors(m)).toEqual([0xff0000, 0xffffff, 0x0000ff]);
    });

    test("resolved numeric strings use the requested range and missing values keep their colors", () => {
        const m = model([{ b: "10.5" }, { properties: { b: "20.5" } }, { color: 0x123456 }]);
        const gradient = { range: () => [0, 1], valueToHex: jest.fn(() => 0xabcdef) };
        m.setColorByProperty({}, "b", gradient, [0, 30]);
        expect(gradient.valueToHex.mock.calls).toEqual([[10.5, [0, 30]], [20.5, [0, 30]]]);
        expect(colors(m)).toEqual([0xabcdef, 0xabcdef, 0x123456]);
    });

    test("nested properties retain their existing coloring", () => {
        const m = model([10, 20, 30].map(score => ({ properties: { score } })));
        m.setColorByProperty({}, "score", "rwb");
        expect(colors(m)).toEqual([0xff0000, 0xffffff, 0x0000ff]);
    });
});
