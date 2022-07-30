/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');

describe('Function VASP | Input: CONTCAR |', () => {
    const data = fs.readFileSync('../test_structs/CONTCAR', 'utf-8');
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