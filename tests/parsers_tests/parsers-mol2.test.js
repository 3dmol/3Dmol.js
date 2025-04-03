
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require('path');

describe('Function MOL2 | Input: multiple.mol2 |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname,'../test_structs/multiple.mol2'), 'utf-8')
    let atoms = $3Dmol.Parsers.MOL2(data, {});// cannot use multimodel

    test("Atoms length is 1", ()=>{
        expect(atoms.length).toBe(1);
    });
    
    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    
});

// noH
describe('Function MOL2 | Input: multiple.mol2 | options: keepH |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname,'../test_structs/multiple.mol2'), 'utf-8')
    let atoms = $3Dmol.Parsers.MOL2(data, {keepH:false}); // noH = true
    
    test("Atoms should match the snapshot (keepH: false)", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    
});
