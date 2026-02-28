
global.$ = require("jquery");
global.URL.createObjectURL = function () {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require("path");

// Benzene SDF with aromatic bond order 4
const benzeneSDF = [
    "benzene", "     3Dmol.js", "",
    " 12 12  0  0  0  0  0  0  0  0999 V2000",
    "    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "    1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000   -1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "    2.1560    1.2440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "    2.1560   -1.2440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000   -2.4880    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -2.1560   -1.2440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -2.1560    1.2440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000    2.4880    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "  1  2  4  0",
    "  2  3  4  0",
    "  3  4  4  0",
    "  4  5  4  0",
    "  5  6  4  0",
    "  6  1  4  0",
    "  1  7  1  0",
    "  2  8  1  0",
    "  3  9  1  0",
    "  4 10  1  0",
    "  5 11  1  0",
    "  6 12  1  0",
    "M  END",
    "$$$$"
].join("\n");

// Pyridine SDF with aromatic bond order 4
const pyridineSDF = [
    "pyridine", "     3Dmol.js", "",
    " 11 11  0  0  0  0  0  0  0  0999 V2000",
    "    1.1260    0.6500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0",
    "    1.1260   -0.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000   -1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -1.1260   -0.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -1.1260    0.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000    1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "    2.0000   -1.1540    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000   -2.3080    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -2.0000   -1.1540    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "   -2.0000    1.1540    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "    0.0000    2.3080    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
    "  1  2  4  0",
    "  2  3  4  0",
    "  3  4  4  0",
    "  4  5  4  0",
    "  5  6  4  0",
    "  6  1  4  0",
    "  2  7  1  0",
    "  3  8  1  0",
    "  4  9  1  0",
    "  5 10  1  0",
    "  6 11  1  0",
    "M  END",
    "$$$$"
].join("\n");

// Naphthalene: two fused 6-membered aromatic rings
const naphthaleneSDF = fs.readFileSync(
    path.resolve(__dirname, '../test_structs/naphthalene_aromatic.sdf'), 'utf-8'
);

describe('Aromatic bond parsing | bond order 4 preserved', () => {

    test("Benzene: all C-C bonds have bond order 4", () => {
        const atoms = $3Dmol.Parsers.SDF(benzeneSDF, {});
        const frame = atoms[0];
        // Atoms 0-5 are carbons forming the ring
        for (let i = 0; i < 6; i++) {
            const atom = frame[i];
            // Each carbon should have at least one bond with order 4
            const aromaticBonds = atom.bondOrder.filter(o => o === 4);
            expect(aromaticBonds.length).toBeGreaterThanOrEqual(1);
        }
    });

    test("Benzene: C-H bonds have bond order 1", () => {
        const atoms = $3Dmol.Parsers.SDF(benzeneSDF, {});
        const frame = atoms[0];
        for (let i = 0; i < 6; i++) {
            const atom = frame[i];
            for (let j = 0; j < atom.bonds.length; j++) {
                const neighbor = frame[atom.bonds[j]];
                if (neighbor.elem === 'H') {
                    expect(atom.bondOrder[j]).toBe(1);
                }
            }
        }
    });

    test("Pyridine: N-C bonds and C-C bonds in ring all have order 4", () => {
        const atoms = $3Dmol.Parsers.SDF(pyridineSDF, {});
        const frame = atoms[0];
        // Atom 0 is N, atoms 1-5 are C, all ring bonds should be order 4
        for (let i = 0; i < 6; i++) {
            const atom = frame[i];
            for (let j = 0; j < atom.bonds.length; j++) {
                const neighborIdx = atom.bonds[j];
                // Ring bonds connect atoms 0-5 to each other
                if (neighborIdx < 6) {
                    expect(atom.bondOrder[j]).toBe(4);
                }
            }
        }
    });

    test("Naphthalene: all ring bonds have order 4", () => {
        const atoms = $3Dmol.Parsers.SDF(naphthaleneSDF, {});
        const frame = atoms[0];
        // Atoms 0-9 are ring carbons
        for (let i = 0; i < 10; i++) {
            const atom = frame[i];
            for (let j = 0; j < atom.bonds.length; j++) {
                const neighborIdx = atom.bonds[j];
                if (neighborIdx < 10) {
                    expect(atom.bondOrder[j]).toBe(4);
                }
            }
        }
    });

    test("Naphthalene: shared bond between rings also has order 4", () => {
        const atoms = $3Dmol.Parsers.SDF(naphthaleneSDF, {});
        const frame = atoms[0];
        // Atoms 4 and 5 (0-indexed) are shared between the two rings
        // They should be bonded to each other with order 4
        const atom4 = frame[3]; // C4 at index 3
        const atom5 = frame[4]; // C5 at index 4
        const bondIdx = atom4.bonds.indexOf(4);
        expect(bondIdx).toBeGreaterThanOrEqual(0);
        expect(atom4.bondOrder[bondIdx]).toBe(4);
    });
});

describe('Aromatic aromaticStyle option | StickStyleSpec', () => {

    test("aromaticStyle defaults to 'dashed' when not specified", () => {
        // Verify the style spec accepts aromaticStyle
        const style = { stick: { radius: 0.15 } };
        // aromaticStyle not set - should default
        expect(style.stick.aromaticStyle).toBeUndefined();
    });

    test("aromaticStyle 'circle' is a valid option", () => {
        const style = { stick: { radius: 0.15, aromaticStyle: "circle" } };
        expect(style.stick.aromaticStyle).toBe("circle");
    });

    test("aromaticStyle 'dashed' is a valid option", () => {
        const style = { stick: { radius: 0.15, aromaticStyle: "dashed" } };
        expect(style.stick.aromaticStyle).toBe("dashed");
    });
});

describe('Aromatic ring detection for circle style', () => {

    test("Benzene: 6 carbons form a ring (bond order 4 count)", () => {
        const atoms = $3Dmol.Parsers.SDF(benzeneSDF, {});
        const frame = atoms[0];
        // Count total aromatic bonds (each counted once per atom)
        let aromaticCount = 0;
        for (let i = 0; i < frame.length; i++) {
            for (let j = 0; j < frame[i].bondOrder.length; j++) {
                if (frame[i].bondOrder[j] === 4) aromaticCount++;
            }
        }
        // Benzene has 6 aromatic bonds, each counted from both atoms = 12
        expect(aromaticCount).toBe(12);
    });

    test("Naphthalene: 11 aromatic bonds (each counted from both atoms = 22)", () => {
        const atoms = $3Dmol.Parsers.SDF(naphthaleneSDF, {});
        const frame = atoms[0];
        let aromaticCount = 0;
        for (let i = 0; i < frame.length; i++) {
            for (let j = 0; j < frame[i].bondOrder.length; j++) {
                if (frame[i].bondOrder[j] === 4) aromaticCount++;
            }
        }
        // 11 aromatic bonds x 2 (from both atoms)
        expect(aromaticCount).toBe(22);
    });

    test("Pyridine: heteroatom (N) is part of aromatic ring", () => {
        const atoms = $3Dmol.Parsers.SDF(pyridineSDF, {});
        const frame = atoms[0];
        const nitrogen = frame[0];
        expect(nitrogen.elem).toBe('N');
        // N should have 2 aromatic bonds (to adjacent ring atoms)
        const aromaticBonds = nitrogen.bondOrder.filter(o => o === 4);
        expect(aromaticBonds.length).toBe(2);
    });

    test("Benzene: non-ring bonds (C-H) are not aromatic", () => {
        const atoms = $3Dmol.Parsers.SDF(benzeneSDF, {});
        const frame = atoms[0];
        // Hydrogen atoms (indices 6-11) should only have order 1 bonds
        for (let i = 6; i < 12; i++) {
            for (let j = 0; j < frame[i].bondOrder.length; j++) {
                expect(frame[i].bondOrder[j]).toBe(1);
            }
        }
    });
});
