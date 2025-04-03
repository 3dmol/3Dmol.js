
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require("path");

describe('Function SDF | parseV2000 |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/aromaticsdf.sdf'), 'utf-8')
    // undefined keepH
    let atoms = $3Dmol.Parsers.SDF(data, {});
    let atomCount = atoms[0].length;
    
    test("Atoms is empty when input data is not greater than 3 lines",()=>{
        let test_data1 = "line1\nline2\nline3";
        let test_atoms1 = $3Dmol.Parsers.SDF(test_data1, {});
        expect(test_atoms1).toEqual([[]]);
    });

    test("Atoms is empty when length of line 4 of input data is less than 38",()=>{
        let test_data2 = "line1\nline2\nline3\nlessthan38";
        let test_atoms2 = $3Dmol.Parsers.SDF(test_data2, {});
        expect(test_atoms2).toEqual([[]]);
    });

    test("Atoms is empty when atomCount is not a number", ()=>{
        let test_data3 = "line1\nline2\nline3\natomCount";
        let test_atoms3 = $3Dmol.Parsers.SDF(test_data3, {});
        expect(test_atoms3).toEqual([[]]);
    });

    test("Atoms is empty when atomCount is not greater than 0", ()=>{
        let test_data4 = "line1\nline2\nline3\n -1";
        let test_atoms4 = $3Dmol.Parsers.SDF(test_data4, {});
        expect(test_atoms4).toEqual([[]]);
    });

    test("Atoms is empty when the number of lines of input is less than 4 + atomCount + bondCount", ()=>{
        let test_data5 = "line1\nline2\nline3\n  1 0 ";// atomCount=0, bondCount=1
        let test_atoms5 = $3Dmol.Parsers.SDF(test_data5, {});
        expect(test_atoms5).toEqual([[]]);
    });

    test("Length of atoms is 1", ()=>{
        expect(atoms.length).toBe(1); 
    });
    
    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

// keepH = false; noH = true;
describe('Function SDF | parseV2000 | options: keepH |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/aromaticsdf.sdf'), 'utf-8')
    let atoms = $3Dmol.Parsers.SDF(data, {keepH:false});

    test("Atoms should match the snapshot when keepH is false", ()=>{
        expect(atoms).toMatchSnapshot();
    });
    
});

// multimodel
describe('Function SDF | parseV2000 | options: multimodel, onemol |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname, '../test_structs/aromaticsdf.sdf'), 'utf-8')

    test("Atoms should match the snapshot when multimodel is true and onemol is false/undefinedd", ()=>{
        let atoms = $3Dmol.Parsers.SDF(data, {multimodel:true, onemol:false}); 
        expect(atoms).toMatchSnapshot();
    });

    test("Atoms should match the snapshot when multimodel is true and onemol is true", ()=>{
        let atoms = $3Dmol.Parsers.SDF(data, {multimodel:true, onemol:true}); 
        expect(atoms).toMatchSnapshot();
    });

});
