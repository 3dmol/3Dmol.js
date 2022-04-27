/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../build/3Dmol.js");
const fs = require('fs');

describe('Function VASP\nInput: CONTCAR', () => {
    const data = fs.readFileSync('tests/test_structs/CONTCAR', 'utf-8');
    let atoms = $3Dmol.Parsers.VASP(data);
    /*
    test("test VASP input", ()=>{
        expect(atoms).toBeDefined();
    });

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
    */
    test("Atoms is not empty", () => {
        expect(atoms).not.toEqual([[]]);
    });

    test("Atoms length is 1", () => {
        expect(atoms.length).toBe(1);
    });

    test("Length of atoms[0]", () => {
        expect(atoms[0].length).toBe(127);
    });

    test("Elem of first 126 atom is C", () => {
        for(let i = 0; i < atoms[0].length-1; i++){
            expect(atoms[0][i].elem).toBe('C');
        }
    });

    test("Elem the last atom is Ge", () => {
        expect(atoms[0][126].elem).toBe('Ge');
    });
    
    test("x of atoms[0][0]", () => {    
        expect(atoms[0][0].x).toBeCloseTo(-0.0022);

    });

    test("y of atoms[0][0]", () => {
        expect(atoms[0][0].y).toBeCloseTo(-0.0022);
    });

    test("z of atoms[0][0]", () => {
        expect(atoms[0][0].z).toBeCloseTo(-0.0022);
    });

    test("Bonds of atoms[0] are empty", () => {
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].bonds).toEqual([]);        
        }
    });

});


describe('Function CUBE', ()=>{ 
    // CUBE returns 'atoms'
    const data = fs.readFileSync('tests/auto/data/1fas.cube', 'utf-8');
    let atoms = $3Dmol.Parsers.CUBE(data, {});
    let len = atoms[0].length; // 913

    test("Atom is not empty", ()=>{
        expect(atoms).not.toEqual([[]]); //if lines.length < 6, atoms is empty
    });

    test("Length of 'atoms' is equal to 1", ()=>{
        expect(atoms.length).toBe(1); //if lines.length < 6, atoms is empty
    });

    test("Length of element inside 'atoms' is equal to 913", ()=>{
        // start = 0; lines.length = 913 based on the file 1fas.cube (from line 7 to line 919)
        expect(atoms[0].length).toBe(913);
    });

    test("All atom hetflag are equal to true", ()=>{
        for(let i = 0; i < len; i++){
            expect(atoms[0][0].hetflag).toBe(true);
        }
    });

    test("Atom serial is equal to index", ()=>{
        for(let i = 0; i < len; i++){
            expect(atoms[0][i].serial).toBe(i);
        }
    });

    test("Elem of atom 0 is H", ()=>{
        expect(atoms[0][0].elem).toBe('H');
    });

    test("Atom 0 is bonded to 8,9,11", ()=>{
        expect(atoms[0][0].bonds).toEqual([8,9,11]);
    });
    
    test("Bond order of atom 0 is ", ()=>{
        expect(atoms[0][0].bondOrder).toEqual([1,1,1]);
    });
    
    test("Properties of all atoms are empty", ()=>{
        for(let i = 0; i < len; i++){
            expect(atoms[0][i].properties).toEqual({});
        }
    });

    test("x of atom [0][0]", ()=>{
        expect(atoms[0][0].x).toBeCloseTo(46.148);
        
    });

    test("y of atom [0][0]", ()=>{
        expect(atoms[0][0].y).toBeCloseTo(16.581);
    });

    test("z of atom [0][0]", ()=>{
        expect(atoms[0][0].z).toBeCloseTo(2.104);
    });

});

describe('Function CUBE\nassignBonds: false', ()=>{ 
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

});

//Comment out for now for time purpose 
/*
describe('Function XYZ\ninput: md.xyz, options:{}', ()=>{
    const data = fs.readFileSync('tests/auto/data/md.xyz', 'utf-8')
    let atomCount = 17411;
    let atoms = $3Dmol.Parsers.XYZ(data, {}); //assignbonds

    test("Atom is empy when input data has less than or equal to 3 lines.", ()=>{
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

    test("Length of first element in atoms is equal to 17411 (atom count)", ()=>{
        expect(atoms[0].length).toBe(atomCount); // atom count 
    });
    
    test("All atom hetflag are equal to true", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][0].hetflag).toBe(true);
        }
    });

    test("Atom serial is equal to index", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].serial).toBe(i);
        }
    });

    test("Atom with odd index has elem 'H', otherwise 'O'", ()=>{
        for(let i = 0; i < atomCount; i++){
            if(i % 2 == 1){ // odd
                expect(atoms[0][i].elem).toBe('H');
            }else{
                expect(atoms[0][i].elem).toBe('O');
            }
        }
    });

    test("Atom with odd index has atom 'H', otherwise 'O'", ()=>{
        // same as elem
        for(let i = 0; i < atomCount; i++){
            if(i % 2 == 1){ // odd
                expect(atoms[0][i].atom).toBe('H');
            }else{
                expect(atoms[0][i].atom).toBe('O');
            }
        }
    });
    
    test("X of atom is 0", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].x).toBe(0);
        }
    });

    test("Y of atom 1", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].y).toBe(0);
        }
    });

    test("Z of atom 1", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].z).toBe(0);
        }
    });

    test("Atom bonds are all empty", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].bonds).toEqual([]);
        }
    });

    test("Atom bondOrder are all empty", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].bondOrder).toEqual([]);
        }
    });

    test("Atom properties are all empty", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].properties).toEqual({});
        }
    });
   
});
 */
/* onemol gives no change */

describe('Function XYZ\ninput: h-bn.xyz, options: assignBonds', ()=>{
    const data = fs.readFileSync('tests/auto/data/h-bn.xyz', 'utf-8');
    // assignBonds affects bonds and bondOrder
    test("Bonds is empty when assignBonds is false", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {assignBonds:false}); //assignbonds
        expect(atoms[0][0].bonds).toEqual([]);
    });

    test("BondOrder is empty when assignBonds is false", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {assignBonds:false}); //assignbonds
        expect(atoms[0][0].bonds).toEqual([]);
    });

    test("BondOrder is not empty when assignBonds is true", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {assignBonds:true}); //assignbonds
        expect(atoms[0][0].bondOrder).not.toEqual([]);
    });

    test("BondOrder is not empty when assignBonds is true", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {assignBonds:true}); //assignbonds
        expect(atoms[0][0].bondOrder).not.toEqual([]);
    });
});

describe('Function XYZ\ninput: h-bn.xyz options: multimodel', ()=>{
    const data = fs.readFileSync('tests/auto/data/h-bn.xyz', 'utf-8');
    test("Length of atoms is 1 when multimodel is undefined", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {}); //assignbonds
        expect(atoms.length).toBe(1);
    });

    test("Length of atoms is greater than 1 when multimodel is true", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {multimodel:true}); //assignbonds
        expect(atoms.length).toBe(2);
    });

    test("Length of atoms is 1 when onemol is true", ()=>{
        let atoms = $3Dmol.Parsers.XYZ(data, {multimodel:true, onemol:true}); 
        expect(atoms.length).toBe(1);
    });

});


describe('Function SDF\nparseV2000', ()=>{
    const data = fs.readFileSync('tests/test_structs/aromaticsdf.sdf', 'utf-8')
    // undefined keepH
    let atoms = $3Dmol.Parsers.SDF(data, {});
    let atomCount = atoms[0].length;
    /*
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
    */

    test("Atoms is not empty", () => {
        expect(atoms).not.toEqual([[]]);
    });

    test("Length of atoms is 1", ()=>{
        expect(atoms.length).toBe(1); 
    });

    test("Length of atoms[0] is equal to atomCount (19)", ()=>{
        expect(atoms[0].length).toBe(atomCount); // atom count 
    });
    
    test("Atom serial is equal to index", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].serial).toBe(i); 
        }
    });

    test("X of atom 1", ()=>{
        expect(atoms[0][1].x).toBeCloseTo(-12.153);
    });

    test("Y of atom 1", ()=>{
        expect(atoms[0][1].y).toBeCloseTo(2.39);
    });

    test("Z of atom 1", ()=>{
        expect(atoms[0][1].z).toBeCloseTo(5.307); 
    });

    test("All atom hetflag is true", ()=>{
        for(let i = 0; i < atomCount; i++){
            expect(atoms[0][i].hetflag).toBeTruthy(); 
        }
    });

    //Elem = H; noH = false
    test("Elem of atoms[0][8] is 'H'", ()=>{
        expect(atoms[0][8].elem).toBe('H');
    });

    test("Atom of atoms[0][8] is equal to elem", ()=>{
        expect(atoms[0][8].atom).toBe(atoms[0][8].elem);
    });

    test("Test bonds of atoms[0][8]", ()=>{
        expect(atoms[0][8].bonds).toEqual([6]);
    });

    test("Test bonds of atoms[0][8]", ()=>{
        expect(atoms[0][8].bondOrder).toEqual([1]);
    });


});

// keepH = false; noH = true;
describe('Function SDF\nparseV2000 options: keepH', ()=>{
    const data = fs.readFileSync('tests/test_structs/aromaticsdf.sdf', 'utf-8')
    let atoms = $3Dmol.Parsers.SDF(data, {keepH:false});
    let atoms2 = $3Dmol.Parsers.SDF(data, {});

    //Elem = 'H'
    test("Serial of atoms[0][8] when keepH is false should be different from atoms[0][8] when keepH is undefined", ()=>{
        expect(atoms[0][8].serial).not.toBe(atoms2[0][8].serial);
    });

    test("Bonds of atoms[0][8] when keepH is false is different from atoms[0][8] when keepH is undefined", ()=>{
        expect(atoms[0][8].bonds).not.toEqual(atoms2[0][8].bonds);
    });

    test("BondOrder of atoms[0][8] when keepH is false is different from atoms[0][8] when keepH is undefined", ()=>{
        expect(atoms[0][8].bonds).not.toEqual(atoms2[0][8].bonds);
    });

    test("x of atoms[0][8] when keepH is false is different from atoms[0][8] when keepH is undefined", ()=>{
        expect(atoms[0][8].x).not.toEqual(atoms2[0][8].x);
    });

    test("y of atoms[0][8] when keepH is false is different from atoms[0][8] when keepH is undefined", ()=>{
        expect(atoms[0][8].x).not.toEqual(atoms2[0][8].x);
    });

    test("z of atoms[0][8] when keepH is false is different from atoms[0][8] when keepH is undefined", ()=>{
        expect(atoms[0][8].x).not.toEqual(atoms2[0][8].x);
    });
    
    // Elem = C; keepH = false; noH = true
    test("Serial of atoms[0][1] when keepH is false should be the same as atoms[0][1] when keepH is undefined", ()=>{
        expect(atoms[0][1].serial).toBe(atoms2[0][1].serial);
    });

    test("Bonds of atoms[0][1] when keepH is false should be the same as atoms[0][1] when keepH is undefined", ()=>{
        expect(atoms[0][1].bonds).toEqual(atoms2[0][1].bonds); 
    });

    test("BondOrder of atoms[0][1] when keepH is false should be the same as atoms[0][1] when keepH is undefined", ()=>{
        expect(atoms[0][1].bondOrder).toEqual(atoms2[0][1].bondOrder);
    });

    test("Properties of atoms[0][1] when keepH is false should be the same as atoms[0][1] when keepH is undefined", ()=>{
        expect(atoms[0][1].properties).toEqual(atoms2[0][1].properties);
    });

    test("x of atoms[0][1] when keepH is false is the same as atoms[0][1] when keepH is undefined", ()=>{
        expect(atoms[0][1].x).toEqual(atoms2[0][1].x);
    });

    test("y of atoms[0][1] when keepH is false is the same as atoms[0][1] when keepH is undefined", ()=>{
        expect(atoms[0][1].y).toEqual(atoms2[0][1].y);
    });

    test("z of atoms[0][1] when keepH is false is the same as atoms[0][1] when keepH is undefined", ()=>{
        expect(atoms[0][1].z).toEqual(atoms2[0][1].z);
    });
    
});

//multimodel
describe('Function SDF\nparseV2000 options: multimodel', ()=>{
    const data = fs.readFileSync('tests/test_structs/aromaticsdf.sdf', 'utf-8')
    let atoms = $3Dmol.Parsers.SDF(data, {multimodel:true});

    test("Length of atoms is 2 when multimodel is true", ()=>{
        expect(atoms.length).toBe(2);
    });

    test("Atoms[1] is empty", ()=>{
        expect(atoms[1]).toEqual([]);
    });

    test("Length of atoms is 1 when onemol is true", ()=>{
        let one_atoms = $3Dmol.Parsers.SDF(data, {multimodel:true, onemol:true});
        expect(one_atoms.length).toBe(1);
    });

});


describe('Function json', ()=>{
    const data = '{"m":[{"a":[{"x":85,"y":144},{"x":102.32050807568878,"y":134},{"x":67.67949192431124,"y":134},{"x":85,"y":164}],"b":[{"b":0,"e":1},{"b":0,"e":2},{"b":0,"e":3}]},{"a":[{"x":213,"y":131},{"x":230.32050807568876,"y":121},{"x":247.64101615137753,"y":131},{"x":264.9615242270663,"y":121}],"b":[{"b":0,"e":1},{"b":1,"e":2},{"b":2,"e":3}]}]}'
    let atoms = $3Dmol.Parsers.json(data, {});

     test("Atom is not empty ", ()=>{
        expect(atoms).not.toBe([]);
    });

});


describe('Function cif\nInput: multiple.cif', ()=>{
    const data = fs.readFileSync('tests/test_structs/multiple.cif', 'utf-8')
    let atoms = $3Dmol.Parsers.cif(data, {});
    let atoms2 = $3Dmol.Parsers.cif(data, {assignBonds:false});
    
    test("Atoms is not empty", ()=>{
        expect(atoms).not.toBe([]);
    });
    
    test("Test atoms length", ()=>{
        expect(atoms.length).toBe(98);
    });

    test("Every element in atoms has the same length", ()=>{
        for(let i = 0; i < atoms.length; i++){
            expect(atoms[i].length).toBe(43);
        }
    });

    test("Test atoms[0][0] x", ()=>{
        expect(atoms[0][0].x).toBeCloseTo(5.9063);
    });

    test("Test atoms[0][0] y", ()=>{
        expect(atoms[0][0].y).toBeCloseTo(24.9584);
    });

    test("Test atoms[0][0] z", ()=>{
        expect(atoms[0][0].z).toBeCloseTo(25.9808);
    });

    test("All atom chain are undefined", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].chain).not.toBeDefined();
            }
        }
    });

    test("All atom resi are undefined", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].resi).not.toBeDefined();
            }
        }
    });

    test("All atom resn are undefined", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].resn).not.toBeDefined();
            }
        }
    });

    test("All atom hetflag are true", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].hetflag).toBe(true);
            }
        }
    });

    test("All atom atom are undefined", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].atom).not.toBeDefined();
            }
        }
    });

    test("Test atoms[0][0] elem", ()=>{
        expect(atoms[0][0].elem).toBe('C');
    });

    test("Test atoms[0][0] bonds", ()=>{
        expect(atoms[0][0].bonds).toEqual([1,28,5]);
    });

    test("Test atoms[0][0] bondOrder", ()=>{
        expect(atoms[0][0].bondOrder).toEqual([1,1,1]);
    });

    test("All atom ss are 'c'", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].ss).toBe('c');
            }
        }
    });

    test("Test atoms serial", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].serial).toBe(j);
            }
        }
    });
    /*
    test("Test atoms bondOrder", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].bondOrder).toEqual(Array(atoms[i][j].bonds.length).fill(1));
            }
        }
    });
    */
    test("Test atoms properties", ()=>{
        for(let i = 0; i < atoms.length; i++){
            for(let j = 0; j < atoms[0].length; j++){
                expect(atoms[i][j].properties).toEqual({});
            }
        }
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

});



describe('Function MOL2\nInput: multiple.mol2', ()=>{
    const data = fs.readFileSync('tests/test_structs/multiple.mol2', 'utf-8')
    let atoms = $3Dmol.Parsers.MOL2(data, {});// cannot use multimodel
    
    test("Atoms is not empty", ()=>{
        expect(atoms).not.toBe([[]]);
    });

    test("Test atoms length", ()=>{
        expect(atoms.length).toBe(1);
    });

    test("natoms is 43 (atoms[0].length)", ()=>{
        expect(atoms[0].length).toBe(43);

    });

    test("Elem of atoms[0][6] is N", ()=>{
        expect(atoms[0][6].elem).toBe('N');
    });

    test("Atom of atoms[0][6] is N.am", ()=>{
        expect(atoms[0][6].atom).toBe('N.am');
    });

    test("Atom serial is the serial id in mol2 file", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].serial).toBe(i+1);
        }
    });
    // pick 1 example to test
    test("Test x of atoms[0][6]", ()=>{
        expect(atoms[0][6].x).toBeCloseTo(3.5113);
    });

    test("Test y of atoms[0][6]", ()=>{
        expect(atoms[0][6].y).toBeCloseTo(23.9905);
    });

    test("Test z of atoms[0][6]", ()=>{
        expect(atoms[0][6].z).toBeCloseTo(25.9254);
    });

    test("Atom index is the index in 'atoms'", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].index).toBe(i);
        }
    });

    test("Test bonds of atoms[0][21]", ()=>{
        expect(atoms[0][21].bonds).toEqual([20,22,23]);
    });

    test("Test bondOrder of atoms[0][21]", ()=>{
        expect(atoms[0][21].bondOrder).toEqual([1,2,1]);
    });
    
    test("Test properties of atoms[0][21]", ()=>{
        let charge = 0.242;
        expect(atoms[0][21].properties).toEqual({'charge' : charge, 
                                                    'partialCharge' : charge});
    });

});

//noH
describe('Function MOL2\nInput: benzene.mol2 options: keepH', ()=>{
    const data = fs.readFileSync('tests/test_structs/benzene.mol2', 'utf-8')

    test("Atom is defined when elem is 'H' and keepH is true or undefined", ()=>{
        let atoms = $3Dmol.Parsers.MOL2(data, {}); // noH = true
        expect(atoms[0][6]).toBeDefined();
    });

    test("Atom is undefined when elem is 'H' and keepH is false", ()=>{
        let atoms = $3Dmol.Parsers.MOL2(data, {keepH:false}); // noH = true
        expect(atoms[0][6]).toBeUndefined();
    });

});


describe('function PDBQT', ()=>{
    const data = fs.readFileSync('tests/test_structs/1B5S.pdb', 'utf-8')
    let atoms = $3Dmol.Parsers.PDBQT(data, {});
    
    test("Atoms is not empty", ()=>{
        expect(atoms).not.toBe([[]]);
    });

    test("PDBQT input", ()=>{
        expect(atoms).toBeDefined();
    });
});

/*
describe('function PQR', ()=>{
    const data = fs.readFileSync('tests/auto/data/1fas.pqr', 'utf-8')
    let atoms = $3Dmol.Parsers.PQR(data, {});
    
    test("PQR input", ()=>{
        expect(atoms).toBeDefined();
    });
});
*/
/*
describe('function MMTF', ()=>{
    const data = fs.readFileSync('tests/auto/data/', 'utf-8')
    let atoms = $3Dmol.Parsers.parseV3000(data);
    
    test("mmtf input", ()=>{
        expect(atoms).toBeDefined();
    });
});

describe('function PRMTOP', ()=>{
    const data = fs.readFileSync('tests/auto/data/model1.prmtop', 'utf-8')
    let atoms = $3Dmol.Parsers.PRMTOP(data);
    
    test("PRMTOP input", ()=>{
        expect(atoms).toBeDefined();
    });
});
*/