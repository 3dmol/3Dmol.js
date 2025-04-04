
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require('path');

describe('Function PDBQT | Input:1SMT.pdb |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname,'../test_structs/1SMT.pdb'), 'utf-8')
    let atoms = $3Dmol.Parsers.PDBQT(data, {});

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function PDBQT | Input:1SMT.pdb | option: multimodel, onemol |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname,'../test_structs/1SMT.pdb'), 'utf-8')
    
    test("Atoms should match the snapshot (multimodel:true, onemol:false/undefined)", ()=>{
        let atoms = $3Dmol.Parsers.PDBQT(data, {multimodel:true, onemol:false});
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot (multimodel:true, onemol:true)", ()=>{
        let atoms = $3Dmol.Parsers.PDBQT(data, {multimodel:true, onemol:true});
        expect(atoms).toMatchSnapshot();
    });

});
