/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function () { };
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require('path');

describe('Function cif | Input: multiple.cif | option: assignBonds |', () => {
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/multiple.cif'), 'utf-8')
    let atoms = $3Dmol.Parsers.cif(data, {});
    let atoms2 = $3Dmol.Parsers.cif(data, { assignBonds: false });

    test("Atoms length is 98", () => {
        expect(atoms.length).toBe(98);
    });

    test("Atoms should match the snapshot", () => {
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms bonds are empty when assignBonds is false", () => {
        for (let i = 0; i < atoms2.length; i++) {
            for (let j = 0; j < atoms2[0].length; j++) {
                expect(atoms2[i][j].bonds).toEqual([]);
            }
        }
    });

    test("Atoms bondOrder are empty when assignBonds is false", () => {
        for (let i = 0; i < atoms2.length; i++) {
            for (let j = 0; j < atoms2[0].length; j++) {
                expect(atoms2[i][j].bondOrder).toEqual([]);
            }
        }
    });

    test("Atoms should match the snapshot (assignBonds: false)", () => {
        expect(atoms2).toMatchSnapshot();
    });

});

describe('Function cif | Input: multiple.cif | option: duplicateAssemblyAtoms, dontConnectDuplicatedAtoms |', () => {
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/multiple.cif'), 'utf-8')
    let atoms = $3Dmol.Parsers.cif(data, { duplicateAssemblyAtoms: true, dontConnectDuplicatedAtoms: false });

    test("Atoms should match the snapshot (duplicateAssemblyAtoms:true, dontConnectDuplicatedAtoms:false)", () => {
        expect(atoms).toMatchSnapshot();
    });

});