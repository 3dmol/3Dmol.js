
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require("path");

describe('Function SDF | parseV2000 |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/aromaticsdf.sdf'), 'utf-8')
    // undefined keepH
    let atoms = $3Dmol.Parsers.SDF(data, {});
    let atomCount = atoms[0].length;
    
    test("Atoms is empty when input data is not greater than 3 lines",()=>{
        let test_data1 = "line1\nline2\nline3";
        let test_atoms1 = $3Dmol.Parsers.SDF(test_data1, {});
        expect(test_atoms1.length).toBe(1);
        expect(test_atoms1[0]).toEqual([]);
        expect(test_atoms1.modelData).toEqual([{}]);
    });

    test("Atoms is empty when length of line 4 of input data is less than 38",()=>{
        let test_data2 = "line1\nline2\nline3\nlessthan38";
        let test_atoms2 = $3Dmol.Parsers.SDF(test_data2, {});
        expect(test_atoms2.length).toBe(1);
        expect(test_atoms2[0]).toEqual([]);
        expect(test_atoms2.modelData).toEqual([{}]);
    });

    test("Atoms is empty when atomCount is not a number", ()=>{
        let test_data3 = "line1\nline2\nline3\natomCount";
        let test_atoms3 = $3Dmol.Parsers.SDF(test_data3, {});
        expect(test_atoms3.length).toBe(1);
        expect(test_atoms3[0]).toEqual([]);
        expect(test_atoms3.modelData).toEqual([{}]);
    });

    test("Atoms is empty when atomCount is not greater than 0", ()=>{
        let test_data4 = "line1\nline2\nline3\n -1";
        let test_atoms4 = $3Dmol.Parsers.SDF(test_data4, {});
        expect(test_atoms4.length).toBe(1);
        expect(test_atoms4[0]).toEqual([]);
        expect(test_atoms4.modelData).toEqual([{}]);
    });

    test("Atoms is empty when the number of lines of input is less than 4 + atomCount + bondCount", ()=>{
        let test_data5 = "line1\nline2\nline3\n  1 0 ";// atomCount=0, bondCount=1
        let test_atoms5 = $3Dmol.Parsers.SDF(test_data5, {});
        expect(test_atoms5.length).toBe(1);
        expect(test_atoms5[0]).toEqual([]);
        expect(test_atoms5.modelData).toEqual([{}]);
    });

    test("Length of atoms is 1", ()=>{
        expect(atoms.length).toBe(1); 
    });
    
    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

// keepH = false; noH = true;
describe('Function SDF | parseV2000 | options: keepH |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/aromaticsdf.sdf'), 'utf-8')
    let atoms = $3Dmol.Parsers.SDF(data, {keepH:false});

    test("Atoms should match the snapshot when keepH is false", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    
});

// 3DMOL_STYLE property parsing
describe('Function SDF | parseV2000 | 3DMOL_STYLE |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/styled_benzene.sdf'), 'utf-8')
    let atoms = $3Dmol.Parsers.SDF(data, {});

    test("Bond style from 3DMOL_STYLE is applied to correct atom bondStyles", ()=>{
        // Bond 1-2 (serials 1,2 → indices 0,1). Lower index is 0.
        // Atom 0 should have bondStyles on the bond to atom 1.
        let atom0 = atoms[0][0];
        let bondToAtom1 = atom0.bonds.indexOf(1);
        expect(bondToAtom1).toBeGreaterThanOrEqual(0);
        expect(atom0.bondStyles).toBeDefined();
        expect(atom0.bondStyles[bondToAtom1]).toBeDefined();
        expect(atom0.bondStyles[bondToAtom1].color1).toBe("red");
        expect(atom0.bondStyles[bondToAtom1].color2).toBe("blue");
    });

    test("solidColor and dashedColor from 3DMOL_STYLE are applied per-bond", ()=>{
        // Bond 3-4 (serials 3,4 → indices 2,3). Lower index is 2.
        let atom2 = atoms[0][2];
        let bondToAtom3 = atom2.bonds.indexOf(3);
        expect(bondToAtom3).toBeGreaterThanOrEqual(0);
        expect(atom2.bondStyles).toBeDefined();
        expect(atom2.bondStyles[bondToAtom3]).toBeDefined();
        expect(atom2.bondStyles[bondToAtom3].dashedBondConfig).toBeDefined();
        expect(atom2.bondStyles[bondToAtom3].dashedBondConfig.solidColor).toBe("green");
        expect(atom2.bondStyles[bondToAtom3].dashedBondConfig.dashedColor).toBe("orange");
    });

    test("Model-level style from 3DMOL_STYLE is stored in modelData", ()=>{
        expect(atoms.modelData).toBeDefined();
        expect(atoms.modelData[0]).toBeDefined();
        expect(atoms.modelData[0].style).toBeDefined();
        expect(atoms.modelData[0].style.stick).toBeDefined();
        expect(atoms.modelData[0].style.stick.radius).toBe(0.12);
    });

    test("Bonds without 3DMOL_STYLE overrides have no bondStyles", ()=>{
        // Atom 4 (index 4) has no bond style overrides
        let atom4 = atoms[0][4];
        // bondStyles may be undefined or have no entries for this atom's bonds
        if (atom4.bondStyles) {
            atom4.bonds.forEach((_, idx) => {
                expect(atom4.bondStyles[idx]).toBeUndefined();
            });
        }
    });
});

// 3DMOL_STYLE with reversed bond serials
describe('Function SDF | parseV2000 | 3DMOL_STYLE reversed serials |', ()=>{
    // Build a minimal SDF with reversed bond serial order in 3DMOL_STYLE
    const sdf = [
        "test", "", "",
        "  3  2  0  0  0  0  0  0  0  0999 V2000",
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
        "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
        "    3.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0",
        "  1  2  1  0",
        "  2  3  1  0",
        "M  END",
        '> <3DMOL_STYLE>',
        '{"bonds":{"2-1":{"color1":"green","color2":"purple"}}}',
        "",
        "$$$$"
    ].join("\n");

    let atoms = $3Dmol.Parsers.SDF(sdf, {});

    test("Reversed serial order swaps color1/color2", ()=>{
        // Bond 2-1 means serial2=2(idx1), serial1=1(idx0). idx1 > idx0, so swap.
        // Lower index atom is 0, which holds the bond to atom 1.
        let atom0 = atoms[0][0];
        let bondToAtom1 = atom0.bonds.indexOf(1);
        expect(bondToAtom1).toBeGreaterThanOrEqual(0);
        expect(atom0.bondStyles).toBeDefined();
        expect(atom0.bondStyles[bondToAtom1]).toBeDefined();
        // Swapped: original color1=green becomes color2, original color2=purple becomes color1
        expect(atom0.bondStyles[bondToAtom1].color1).toBe("purple");
        expect(atom0.bondStyles[bondToAtom1].color2).toBe("green");
    });
});

// V3000 with 3DMOL_STYLE
describe('Function SDF | parseV3000 | 3DMOL_STYLE |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/styled_pyridine_v3000.sdf'), 'utf-8')
    let atoms = $3Dmol.Parsers.SDF(data, {});

    test("V3000 parses correct number of atoms", ()=>{
        expect(atoms[0].length).toBe(11);
    });

    test("V3000 parses aromatic bond orders", ()=>{
        // Atom 0 (N) bonds to atoms 1 and 5 with order 1.5
        let atom0 = atoms[0][0];
        let bondToAtom1 = atom0.bonds.indexOf(1);
        let bondToAtom5 = atom0.bonds.indexOf(5);
        expect(bondToAtom1).toBeGreaterThanOrEqual(0);
        expect(bondToAtom5).toBeGreaterThanOrEqual(0);
        expect(atom0.bondOrder[bondToAtom1]).toBe(1.5);
        expect(atom0.bondOrder[bondToAtom5]).toBe(1.5);
    });

    test("V3000 3DMOL_STYLE bond colors are applied", ()=>{
        // Bond 2-3 (serials 2,3 → indices 1,2). Lower index is 1.
        let atom1 = atoms[0][1];
        let bondToAtom2 = atom1.bonds.indexOf(2);
        expect(bondToAtom2).toBeGreaterThanOrEqual(0);
        expect(atom1.bondStyles).toBeDefined();
        expect(atom1.bondStyles[bondToAtom2]).toBeDefined();
        expect(atom1.bondStyles[bondToAtom2].color1).toBe("red");
        expect(atom1.bondStyles[bondToAtom2].color2).toBe("blue");
    });

    test("V3000 3DMOL_STYLE solidColor/dashedColor are applied", ()=>{
        // Bond 5-6 (serials 5,6 → indices 4,5). Lower index is 4.
        let atom4 = atoms[0][4];
        let bondToAtom5 = atom4.bonds.indexOf(5);
        expect(bondToAtom5).toBeGreaterThanOrEqual(0);
        expect(atom4.bondStyles).toBeDefined();
        expect(atom4.bondStyles[bondToAtom5]).toBeDefined();
        expect(atom4.bondStyles[bondToAtom5].dashedBondConfig).toBeDefined();
        expect(atom4.bondStyles[bondToAtom5].dashedBondConfig.solidColor).toBe("#00FF00");
        expect(atom4.bondStyles[bondToAtom5].dashedBondConfig.dashedColor).toBe("#FF8800");
    });

    test("V3000 3DMOL_STYLE model-level style is stored", ()=>{
        expect(atoms.modelData).toBeDefined();
        expect(atoms.modelData[0]).toBeDefined();
        expect(atoms.modelData[0].style).toBeDefined();
        expect(atoms.modelData[0].style.stick).toBeDefined();
        expect(atoms.modelData[0].style.stick.radius).toBe(0.12);
    });
});

// multimodel
describe('Function SDF | parseV2000 | options: multimodel, onemol |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/aromaticsdf.sdf'), 'utf-8')

    test("Atoms should match the snapshot when multimodel is true and onemol is false/undefinedd", ()=>{
        let atoms = $3Dmol.Parsers.SDF(data, {multimodel:true, onemol:false}); 
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot when multimodel is true and onemol is true", ()=>{
        let atoms = $3Dmol.Parsers.SDF(data, {multimodel:true, onemol:true}); 
        expect(atoms).toMatchSnapshot();
    });

});
