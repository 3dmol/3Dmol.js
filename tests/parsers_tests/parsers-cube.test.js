
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');
const path = require('path');

describe('Function CUBE | Input:1fas.cube |', ()=>{
    // CUBE returns 'atoms'
    const data = fs.readFileSync(path.resolve(__dirname,'../auto/data/1fas.cube'), 'utf-8');
    let atoms = $3Dmol.Parsers.CUBE(data, {}); // assignBonds is undefined (same as true)

    test("Atom is empty if input has less than 6 lines", () => {
        let str = "line1\nline2\nline3";
        let result = $3Dmol.Parsers.CUBE(str);
        expect(result).toEqual([[]]); // if lines.length < 6, atoms is empty
    });

    test("Length of atoms is equal to 1", ()=>{
        expect(atoms.length).toBe(1);
    });

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});

describe('Function CUBE | Input:1fas.cube | assignBonds: false |', ()=>{ 
    const data = fs.readFileSync(path.resolve(__dirname,'../auto/data/1fas.cube'), 'utf-8');
    let atoms = $3Dmol.Parsers.CUBE(data, {assignBonds:false});

    test("All bonds are empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].bonds).toEqual([]);
        }
    });
    
    test("All bondOrder are empty when assignBonds is false", ()=>{
        for(let i = 0; i < atoms[0].length; i++){
            expect(atoms[0][i].bondOrder).toEqual([]);
        }
    });

    test("Atoms should match the snapshot (assignBonds:false)", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});
