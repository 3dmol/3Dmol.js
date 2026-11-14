global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const fs = require("fs");
const path = require("path");

const data = fs.readFileSync(path.resolve(__dirname, "../test_structs/1ubq.cif"), "utf8");

function withoutSecondaryStructure(models) {
    return models.map(atoms => atoms.map(atom => {
        const { ss, ssbegin, ssend, ...rest } = atom;
        return rest;
    }));
}

const terBlock = `data_first
loop_
_atom_site.id
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 ATOM N N ALA A 1 0 0 0
2 TER C CA ALA A 1 1 0 0
3 ATOM C CA ALA A 1 2 0 0
#
`;

describe("CIF secondary structure options", () => {
    test("default and explicitly enabled computation retain helices and sheets", () => {
        const defaultModels = $3Dmol.Parsers.cif(data, {});
        const enabled = $3Dmol.Parsers.cif(data, { noComputeSecondaryStructure: false });
        expect(enabled).toEqual(defaultModels);
        expect(enabled[0].some(atom => atom.ss === "h")).toBe(true);
        expect(enabled[0].some(atom => atom.ss === "s")).toBe(true);
    });

    test.each(["cif", "mcif"])("%s skips secondary structure computation when requested", format => {
        const models = $3Dmol.Parsers[format](data, { noComputeSecondaryStructure: true });
        expect(models[0].length).toBeGreaterThan(0);
        expect(models[0].every(atom => atom.ss === "c")).toBe(true);
        expect(models[0].every(atom => !atom.ssbegin && !atom.ssend)).toBe(true);
    });

    test("disabling computation preserves coordinates, bonds, atom metadata and model metadata", () => {
        const baseline = $3Dmol.Parsers.cif(data, {});
        const skipped = $3Dmol.Parsers.cif(data, { noComputeSecondaryStructure: true });
        expect(withoutSecondaryStructure(skipped)).toEqual(withoutSecondaryStructure(baseline));
        expect(skipped.modelData).toEqual(baseline.modelData);
        expect(skipped[0].some(atom => atom.bonds.length > 0)).toBe(true);
    });

    test("disabling both bond assignment and computation keeps contiguous atom indices", () => {
        const models = $3Dmol.Parsers.cif(data, {
            assignBonds: false, noComputeSecondaryStructure: true
        });
        expect(models[0].every((atom, index) => atom.index === index)).toBe(true);
        expect(models[0].every(atom => atom.bonds.length === 0 && atom.bondOrder.length === 0)).toBe(true);
        expect(models[0].every(atom => atom.ss === "c")).toBe(true);
    });

    test("model selections by atom index work without bond or secondary structure computation", () => {
        const model = new $3Dmol.GLModel(0);
        model.addMolData(data, "cif", { assignBonds: false, noComputeSecondaryStructure: true });
        const atoms = model.selectedAtoms({});
        const indices = [0, 10, atoms.length - 1];
        expect(model.selectedAtoms({ index: indices })).toEqual(indices.map(index => atoms[index]));
    });

    test.each([1, 2])("indices remain contiguous after skipped TER rows in %i data blocks", count => {
        const input = count === 1 ? terBlock : terBlock + terBlock.replace("data_first", "data_second");
        const models = $3Dmol.Parsers.cif(input, {
            assignBonds: false, noComputeSecondaryStructure: true
        });
        expect(models).toHaveLength(count);
        for (const atoms of models) {
            expect(atoms.map(atom => atom.index)).toEqual([0, 1]);
            expect(atoms.map(atom => atom.serial)).toEqual([0, 2]);
        }
    });
});
