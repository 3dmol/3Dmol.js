global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const fs = require("fs");
const path = require("path");

const data = fs.readFileSync(path.resolve(__dirname, "../test_structs/1ubq.cif"), "utf8");
const secondaryStructure = model => model.selectedAtoms({}).map(atom => atom.ss);
const reference = new $3Dmol.GLModel(0);
reference.addMolData(data, "cif");
const computed = secondaryStructure(reference);
const originalFetch = global.fetch;

function viewer(parse = true) {
    const model = new $3Dmol.GLModel(0);
    if (!parse) model.addMolData = jest.fn();
    return { addModel: jest.fn(() => model), zoomTo: jest.fn(), render: jest.fn(), model };
}

beforeEach(() => {
    global.fetch = jest.fn(async () => ({
        ok: true, text: async () => data, arrayBuffer: async () => new ArrayBuffer(0)
    }));
});
afterEach(() => { global.fetch = originalFetch; });

describe("downloaded CIF secondary structure", () => {
    test.each(["cif", "mcif", "CIF"])("PDB download in %s format keeps computing by default", async format => {
        const v = viewer();
        const model = await $3Dmol.download("pdb:1ubq", v, { format });
        expect(model).toBe(v.model);
        expect(secondaryStructure(model)).toEqual(computed);
        expect(computed).toContain("h");
        expect(computed).toContain("s");
        expect(global.fetch).toHaveBeenCalledWith(`https://files.rcsb.org/view/1UBQ.${format}`);
        expect(v.zoomTo).toHaveBeenCalledTimes(1);
        expect(v.render).toHaveBeenCalledTimes(1);
    });

    test.each([false, true])("explicit noComputeSecondaryStructure=%s is honored", async skip => {
        const v = viewer();
        const options = { format: "cif", noComputeSecondaryStructure: skip };
        const model = await $3Dmol.download("pdb:1ubq", v, options);
        expect(secondaryStructure(model)).toEqual(skip ? computed.map(() => "c") : computed);
        expect(options.noComputeSecondaryStructure).toBe(skip);
    });

    test("callback downloads keep their return value and computed default", async () => {
        const v = viewer();
        const callback = jest.fn();
        let immediate;
        const loaded = new Promise(resolve => {
            immediate = $3Dmol.download("pdb:1ubq", v, { format: "cif" }, model => {
                callback(model);
                resolve(model);
            });
        });
        expect(immediate).toBe(v.model);
        const model = await loaded;
        expect(model).toBe(immediate);
        expect(secondaryStructure(model)).toEqual(computed);
        expect(callback).toHaveBeenCalledTimes(1);
        expect(callback).toHaveBeenCalledWith(model);
    });

    test.each([
        ["pdb:1ubq", { format: "pdb" }, true],
        ["pdb:1ubq", {}, true],
        ["bcif:1ubq", {}, true],
        ["pdb:1ubq", { format: "cif.gz" }, false]
    ])("%s with %j passes the appropriate parser default", async (query, options, expected) => {
        const v = viewer(false);
        await $3Dmol.download(query, v, options);
        expect(v.model.addMolData).toHaveBeenCalledTimes(1);
        expect(v.model.addMolData.mock.calls[0][2].noComputeSecondaryStructure).toBe(expected);
    });
});
