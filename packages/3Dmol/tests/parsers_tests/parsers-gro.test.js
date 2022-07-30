/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');

describe('Function GRO | Input:2water.gro |', ()=>{
    const data = fs.readFileSync('../auto/data/2water.gro', 'utf-8')
    let atoms = $3Dmol.Parsers.GRO(data);
    let atomCount = atoms[0].length; // 6

    test("Atoms is empty when input data is not greater than 3 lines",()=>{
        let test_data1 = "line1\nline2\nline3";
        let test_atoms1 = $3Dmol.Parsers.GRO(test_data1, {});
        expect(test_atoms1).toEqual([]);
    });

    test("Atoms is empty when atomCount is not a number", ()=>{
        let test_data3 = "line1\natomCount\nline3\nline4"; // atomCount = 6
        let test_atoms3 = $3Dmol.Parsers.GRO(test_data3, {});
        expect(test_atoms3).toEqual([]);
    });

    test("Atoms is empty when atomCount is not greater than 0", ()=>{
        let test_data4 = "line1\n-1\nline3\nline4";
        let test_atoms4 = $3Dmol.Parsers.GRO(test_data4, {});
        expect(test_atoms4).toEqual([]);
    });

    test("Atoms is empty when the number of lines of input is less than atomCount + 3", ()=>{
        let test_data5 = "line1\n2\nline3\nline4"; // 4 < 2 + 3
        let test_atoms5 = $3Dmol.Parsers.GRO(test_data5, {});
        expect(test_atoms5).toEqual([]);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});