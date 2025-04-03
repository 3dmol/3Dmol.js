
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require("path");

describe('Function PRMTOP | Input:model1.prmtop |', ()=>{
    const data = fs.readFileSync(path.resolve(__dirname,'../auto/data/model1.prmtop'), 'utf-8')
    let atoms = $3Dmol.Parsers.PRMTOP(data);
    
    test("Atoms is empty if input data has less than 1 line", ()=>{
        let str = '';
        let result = $3Dmol.Parsers.PRMTOP(str);
        expect(result).toEqual([]);
    });

    test("Atoms is empty if first line of input does not include 'VERSION'", ()=>{
        let str = 'line1\n';
        let result = $3Dmol.Parsers.PRMTOP(str);
        expect(result).toEqual([]);
    });

    test("Length of atoms should be 1", ()=>{
        expect(atoms.length).toBe(1);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});
