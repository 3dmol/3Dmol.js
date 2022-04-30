/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../build/3Dmol.js");
const fs = require('fs');


describe('Function VASP | Input: CONTCAR |', () => {
    const data = fs.readFileSync('tests/test_structs/CONTCAR', 'utf-8');
    let atoms = $3Dmol.Parsers.VASP(data);
    
    test("Atom is empty if input has less than 3 lines", () => {
        let str = "data";
        let result = $3Dmol.Parsers.VASP(str);
        expect(result).toEqual([[]]);
    });

    test("Atom is empty if second line of input is not a number", () => {
        let str = "line\nline\nline";
        let result = $3Dmol.Parsers.VASP(str);
        expect(result).toEqual([[]]);
    });

    test("Atoms is empty if second line of input is negative", () => {
        let str = "line1\n-1\nline3";
        let result = $3Dmol.Parsers.VASP(str);
        expect(result).toEqual([[]]);
    });
    
    test("Atoms length is 1", () => {
        expect(atoms.length).toBe(1);
    });
    
    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});


describe('Function CUBE', ()=>{ 
    // CUBE returns 'atoms'
    const data = fs.readFileSync('tests/auto/data/1fas.cube', 'utf-8');
    let atoms = $3Dmol.Parsers.CUBE(data, {});

    test("Atom is empty if input has less than 6 lines", () => {
        let str = "line1\nline2\nline3";
        let result = $3Dmol.Parsers.CUBE(str);
        expect(result).toEqual([[]]); // if lines.length < 6, atoms is empty
    });

    test("Length of atoms is equal to 1", ()=>{
        expect(atoms.length).toBe(1);
    });

    test("Length of atoms[0] is 913", ()=>{
        let len = atoms[0].length;
        expect(atoms[0].length).toBe(len);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function CUBE | assignBonds: false |', ()=>{ 
    const data = fs.readFileSync('tests/auto/data/1fas.cube', 'utf-8');
    let atoms = $3Dmol.Parsers.CUBE(data, {assignBonds:false});

    test("All bonds are empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].bonds).toEqual([]);
        }
    });
    
    test("All bondOrder are empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].bondOrder).toEqual([]);
        }
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function XYZ | Input: C111tiny.xyz | options:{} |', ()=>{
    const data = fs.readFileSync('tests/auto/data/C111tiny.xyz', 'utf-8')
    let atomCount = 17411;
    let atoms = $3Dmol.Parsers.XYZ(data, {}); //assignbonds

    test("Atom is empty when input data has less than or equal to 3 lines.", ()=>{
        let test_data = '1\nline2';
        let empty_atom = $3Dmol.Parsers.XYZ(test_data, {});
        expect(empty_atom).toEqual([[]]); //if lines.length < 3, atoms is empty
    });
    
    test("Atom is empty when the atom count is less than 0", ()=>{ //first line (atom count)
        let test_data = '0\n20\n30\n';
        let empty_atom = $3Dmol.Parsers.XYZ(test_data, {});
        expect(empty_atom).toEqual([[]]);
    });

    test("Atom is empty when the atom count is not a number.", ()=>{
        let test_data = 'line1\nline2\nline3\n';
        let empty_atom = $3Dmol.Parsers.XYZ(test_data, {});
        expect(empty_atom).toEqual([[]]);
    });

    test("Atom is empty when input data has less than 'atom count + 2' lines", ()=>{
        let test_data = '2\nline2\n'; // atom count = 2
        let empty_atom = $3Dmol.Parsers.XYZ(test_data, {});
        expect(empty_atom).toEqual([[]]);
    });

    test("Length of atoms is 1", ()=>{
        expect(atoms.length).toBe(1); 
    });
    
    test("Length of atoms[0] is 50 (atom count)", ()=>{
        let atomCount = 50;
        expect(atoms[0].length).toBe(atomCount); // atom count 
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });   

});
 
/* onemol gives no change */

describe('Function XYZ | Input: C111tiny.xyz | options: assignBonds', ()=>{
    const data = fs.readFileSync('tests/auto/data/C111tiny.xyz', 'utf-8');
    let atoms = $3Dmol.Parsers.XYZ(data, {assignBonds:false}); //assignbonds

    // assignBonds affects bonds and bondOrder
    test("Bonds is empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].bonds).toEqual([]);
        }
    });

    test("BondOrder is empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].bondOrder).toEqual([]);
        }
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });   

});

describe('Function XYZ | Input: C111tiny.xyz | options: multimodel, onemol |', ()=>{
    const data = fs.readFileSync('tests/auto/data/C111tiny.xyz', 'utf-8');

    test("Atoms should match the snapshot when multimodel is true", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {multimodel:true}); 
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot when onemol is true", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {onemol:true}); 
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot when both multimodel and onemol are true", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {multimodel:true, onemol:true}); 
        expect(atoms).toMatchSnapshot();
    });

});


describe('Function SDF | parseV2000 |', ()=>{
    const data = fs.readFileSync('tests/test_structs/aromaticsdf.sdf', 'utf-8')
    // undefined keepH
    let atoms = $3Dmol.Parsers.SDF(data, {});
    let atomCount = atoms[0].length;
    
    test("Atoms is empty when input data is not greater than 3 lines",()=>{
        let test_data1 = "line1\nline2\nline3";
        let test_atoms1 = $3Dmol.Parsers.SDF(test_data1, {});
        expect(test_atoms1).toEqual([[]]);
    });

    test("Atoms is empty when length line 4 of input data is less than 38",()=>{
        let test_data2 = "line1\nline2\nline3\nlessthan38";
        let test_atoms2 = $3Dmol.Parsers.SDF(test_data2, {});
        expect(test_atoms2).toEqual([[]]);
    });

    test("Atoms is empty when atomCount is not a number", ()=>{
        let test_data3 = "line1\nline2\nline3\natomCount";
        let test_atoms3 = $3Dmol.Parsers.SDF(test_data3, {});
        expect(test_atoms3).toEqual([[]]);
    });

    test("Atoms is empty when atomCount is not greater than 0", ()=>{
        let test_data4 = "line1\nline2\nline3\n -1";
        let test_atoms4 = $3Dmol.Parsers.SDF(test_data4, {});
        expect(test_atoms4).toEqual([[]]);
    });

    test("Atoms is empty when the number of data lines < 4 + atomCount + bondCount", ()=>{
        let test_data5 = "line1\nline2\nline3\n  1 0 ";// atomCount=0, bondCount=1
        let test_atoms5 = $3Dmol.Parsers.SDF(test_data5, {});
        expect(test_atoms5).toEqual([[]]);
    });

    test("Length of atoms is 1", ()=>{
        expect(atoms.length).toBe(1); 
    });

    test("Length of atoms[0] is equal to atomCount (19)", ()=>{
        expect(atoms[0].length).toBe(atomCount); // atom count 
    });
    
    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

// keepH = false; noH = true;
describe('Function SDF | parseV2000 | options: keepH', ()=>{
    const data = fs.readFileSync('tests/test_structs/aromaticsdf.sdf', 'utf-8')
    let atoms = $3Dmol.Parsers.SDF(data, {keepH:false});

    test("Atoms should match the snapshot when keepH is false", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    
});

// multimodel
describe('Function SDF | parseV2000 | options: multimodel, onemol |', ()=>{
    const data = fs.readFileSync('tests/test_structs/aromaticsdf.sdf', 'utf-8')

    test("Atoms should match the snapshot when multimodel is true and onemol is false/undefinedd", ()=>{
        let atoms = $3Dmol.Parsers.SDF(data, {multimodel:true, onemol:false}); 
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot when multimodel is true and onemol is true", ()=>{
        let atoms = $3Dmol.Parsers.SDF(data, {multimodel:true, onemol:true}); 
        expect(atoms).toMatchSnapshot();
    });

});


describe('Function json', ()=>{
    const data = '{"m":[{"a":[{"x":85,"y":144},{"x":102.32050807568878,"y":134},{"x":67.67949192431124,"y":134},{"x":85,"y":164}],"b":[{"b":0,"e":1},{"b":0,"e":2},{"b":0,"e":3}]},{"a":[{"x":213,"y":131},{"x":230.32050807568876,"y":121},{"x":247.64101615137753,"y":131},{"x":264.9615242270663,"y":121}],"b":[{"b":0,"e":1},{"b":1,"e":2},{"b":2,"e":3}]}]}'
    let atoms = $3Dmol.Parsers.json(data, {});

     test("Atom is not empty ", ()=>{
        expect(atoms).not.toBe([]);
    });

});


describe('Function cif | Input: multiple.cif | option: assignBonds', ()=>{
    const data = fs.readFileSync('tests/test_structs/multiple.cif', 'utf-8')
    let atoms = $3Dmol.Parsers.cif(data, {});
    let atoms2 = $3Dmol.Parsers.cif(data, {assignBonds:false});
    
    test("Atoms is not empty", ()=>{
        expect(atoms).not.toBe([]);
    });
    
    test("Atoms length is 98", ()=>{
        expect(atoms.length).toBe(98);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms bonds are empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms2.length; i++){
            for(let j = 0; j < atoms2[0].length; j++){
                expect(atoms2[i][j].bonds).toEqual([]);
            }
        }
    });

    test("Atoms bondOrder are empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms2.length; i++){
            for(let j = 0; j < atoms2[0].length; j++){
                expect(atoms2[i][j].bondOrder).toEqual([]);
            }
        }
    });

    test("Atoms should match the snapshot (assignBonds: false)", ()=>{
        expect(atoms2).toMatchSnapshot();
    });

});

describe('Function cif | Input: multiple.cif | option: duplicateAssemblyAtoms, dontConnectDuplicatedAtoms', ()=>{
    const data = fs.readFileSync('tests/test_structs/multiple.cif', 'utf-8')
    let atoms = $3Dmol.Parsers.cif(data, {duplicateAssemblyAtoms:true, dontConnectDuplicatedAtoms:false});

    test("Atoms should match the snapshot (duplicateAssemblyAtoms:true, dontConnectDuplicatedAtoms:false)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});


describe('Function MOL2 | Input: multiple.mol2 |', ()=>{
    const data = fs.readFileSync('tests/test_structs/multiple.mol2', 'utf-8')
    let atoms = $3Dmol.Parsers.MOL2(data, {});// cannot use multimodel
    
    test("Atoms is not empty", ()=>{
        expect(atoms).not.toBe([[]]);
    });

    test("Atoms length is 1", ()=>{
        expect(atoms.length).toBe(1);
    });
    
    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    
});

//noH
describe('Function MOL2 | Input: multiple.mol2 | options: keepH |', ()=>{
    const data = fs.readFileSync('tests/test_structs/multiple.mol2', 'utf-8')
    let atoms = $3Dmol.Parsers.MOL2(data, {keepH:false}); // noH = true
    
    test("Atoms should match the snapshot (keepH: false)", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    
});


describe('Function PDBQT', ()=>{
    const data = fs.readFileSync('tests/test_structs/1B5S.pdb', 'utf-8')
    let atoms = $3Dmol.Parsers.PDBQT(data, {});

    test("Length of atoms should be 1", ()=>{
        expect(atoms.length).toBe(1);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function PDBQT | option: multimodel, onemol', ()=>{
    const data = fs.readFileSync('tests/test_structs/1B5S.pdb', 'utf-8')
    
    test("Atoms should match the snapshot (multimodel:true, onemol:false/undefined)", ()=>{
        let atoms = $3Dmol.Parsers.PDBQT(data, {multimodel:true, onemol:false});
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot (multimodel:true, onemol:true)", ()=>{
        let atoms = $3Dmol.Parsers.PDBQT(data, {multimodel:true, onemol:true});
        expect(atoms).toMatchSnapshot();
    });

});


describe('Function PQR', ()=>{
    const data = fs.readFileSync('tests/auto/data/1fas.pqr', 'utf-8')
    let atoms = $3Dmol.Parsers.PQR(data, {});

    test("Length of atoms should be 1", ()=>{
        expect(atoms.length).toBe(1);
    });

    /*
    test("snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    */ 

});

/*
describe('Function MMTF', ()=>{
    const data = fs.readFileSync('tests/test_structs/1B5S.pdb', 'utf-8')
    //let atoms = $3Dmol.Parsers.mmtf(data, {});
    global.$3Dmol.download("pdb:2V0E", global.$3Dmol.viewers,{multimodel:true, frames:true})
        .then(function(m1){
            global.$3Dmol.download("mmtf:4HHB", global.$3Dmol.viewers, {multimodel:true, frames:true})
                .then(function(m2){

                });
        });
    
    test("mmtf input", ()=>{
        //expect(atoms).toBeDefined();
    });
});
*/
/*
describe('function PRMTOP', ()=>{
    const data = fs.readFileSync('tests/auto/data/model1.prmtop', 'utf-8')
    let atoms = $3Dmol.Parsers.PRMTOP(data);
    
    test("PRMTOP input", ()=>{
        expect(atoms).toBeDefined();
    });
});
*/