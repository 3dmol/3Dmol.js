
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require('path');

describe('Function LAMMPSTRJ | Input:Al.lammpstrj |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../auto/data/Al.lammpstrj'), 'utf-8')
    let atoms = $3Dmol.Parsers.LAMMPSTRJ(data, {});

    test("Atoms should match the snapshot (assignBonds is undefined/false)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function LAMMPSTRJ | Input:Al.lammpstrj | option:assignBonds |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../auto/data/Al.lammpstrj'), 'utf-8')
    let atoms = $3Dmol.Parsers.LAMMPSTRJ(data, {assignBonds:true});

    test("Atoms should match the snapshot (assignBonds:true)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});
