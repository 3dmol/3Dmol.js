/**
*@jest-environment jsdom
*/
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');

const mmtfstr = fs.readFileSync(__dirname + '/testfile.mmtf', 'utf8');

describe('Function MMTF | options:{} |', ()=>{
    let atoms = $3Dmol.Parsers.MMTF(mmtfstr, {});

    test("Atom should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function MMTF | options:doAssembly |', ()=>{
    let atoms = $3Dmol.Parsers.MMTF(mmtfstr, {doAssembly:true});  // !noAssembly = true

    test("Atom should match the snapshot (doAssembly:true)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function MMTF | options:keepH |', ()=>{
    let atoms = $3Dmol.Parsers.MMTF(mmtfstr, {keepH:true});  // noH = false

    test("Atom should match the snapshot (keepH:true)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function MMTF | options:noComputeSecondaryStructure |', ()=>{
    let atoms = $3Dmol.Parsers.MMTF(mmtfstr, {noComputeSecondaryStructure:true});  // computeStruct = false -> won't do computeSecondaryStructure()

    test("Atom should match the snapshot (noComputeSecondaryStructure:true)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function MMTF | options:altLoc |', ()=>{
    let atoms = $3Dmol.Parsers.MMTF(mmtfstr, {altLoc:'B'}); // default: A

    test("Atom should match the snapshot (altLoc:B)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function MMTF | options:assemblyIndex |', ()=>{
    let atoms = $3Dmol.Parsers.MMTF(mmtfstr, {assemblyIndex:1}); // default: 0

    test("Atom should match the snapshot (assemblyIndex:1)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});