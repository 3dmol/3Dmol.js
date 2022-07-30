/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');

describe('Function XYZ | Input: C111tiny.xyz | options:{} |', ()=>{
    const data = fs.readFileSync('../auto/data/C111tiny.xyz', 'utf-8')
    let atomCount = 17411;
    let atoms = $3Dmol.Parsers.XYZ(data, {}); 

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

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });   

});

describe('Function XYZ | Input: C111tiny.xyz | options: assignBonds |', ()=>{
    const data = fs.readFileSync('../auto/data/C111tiny.xyz', 'utf-8');
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
    const data = fs.readFileSync('../auto/data/C111tiny.xyz', 'utf-8');

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