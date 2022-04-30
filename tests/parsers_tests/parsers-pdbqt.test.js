/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');

describe('Function PDBQT | Input:1B5S.pdb |', ()=>{
    const data = fs.readFileSync('tests/test_structs/1B5S.pdb', 'utf-8')
    let atoms = $3Dmol.Parsers.PDBQT(data, {});

    test("Length of atoms should be 1", ()=>{
        expect(atoms.length).toBe(1);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function PDBQT | Input:1B5S.pdb | option: multimodel, onemol |', ()=>{
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