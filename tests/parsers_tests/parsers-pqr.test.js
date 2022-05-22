/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');

describe('Function PQR | Input:1fas.pqr |', ()=>{
    const data = fs.readFileSync('tests/auto/data/1fas.pqr', 'utf-8')
    let atoms = $3Dmol.Parsers.PQR(data, {});

    test("Length of atoms should be 1", ()=>{
        expect(atoms.length).toBe(1);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function PQR | Input:1fas.pqr | options:multimodel, onemol |', ()=>{
    const data = fs.readFileSync('tests/auto/data/1fas.pqr', 'utf-8')

    test("Atoms should match the snapshot (multimodel:true, onemol:false)", ()=>{
        let atoms = $3Dmol.Parsers.PQR(data, {multimodel:true, onemol:false});
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot (multimodel:true, onemol:true)", ()=>{
        let atoms = $3Dmol.Parsers.PQR(data, {multimodel:true, onemol:true});
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function PQR | Input:1fas.pqr | options:noSecondaryStructure |', ()=>{
    const data = fs.readFileSync('tests/auto/data/1fas.pqr', 'utf-8')

    test("Atoms should match the snapshot (noSecondaryStructure:true)", ()=>{
        let atoms = $3Dmol.Parsers.PQR(data, {noSecondaryStructure:true}); // computeStruct = false
        expect(atoms).toMatchSnapshot();
    });

});