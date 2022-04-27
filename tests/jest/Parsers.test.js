/**
 *@jest-environment jsdom
 */

global.$ = require('jquery');
global.URL.createObjectURL = function () {};
let $3Dmol = require('../../build/3Dmol.js');
const fs = require('fs');
/*
test("CUBE exists", () => {
    let result = $3Dmol.Parsers.CUBE('test', {multimodel:true, frames:true});
    expect(result).toBeDefined();
});
*/
describe('function VASP', () => {
  const data = fs.readFileSync('tests/test_structs/CONTCAR', 'utf-8');
  let atoms = $3Dmol.Parsers.VASP(data);
  test('test VASP input', () => {
    expect(atoms).toBeDefined();
  });

  /*
    test("Input str has less than three lines", () => {
        let str = "";
        let result = $3Dmol.Parsers.VASP(str);
        expect(result).toEqual([[]]);
    });

    test("Second line of vasp input should be a number", () => {
        let str = "";
        let result = $3Dmol.Parsers.VASP(str);
        expect(result).toEqual([[]]);
    });
    */
});

describe('function CUBE', () => {
  // CUBE returns 'atoms'
  const data = fs.readFileSync('tests/auto/data/1fas.cube', 'utf-8');
  let atoms = $3Dmol.Parsers.CUBE(data, { });
  let len = atoms[0].length; // 913

  test('Input data has more than or equal to 6 lines', () => {
    expect(atoms).not.toEqual([[]]); //if lines.length < 6, atoms is empty
  });

  test("Length of 'atoms' is equal to 1", () => {
    expect(atoms.length).toBe(1); //if lines.length < 6, atoms is empty
  });

  test("Length of element inside 'atoms' is equal to 913", () => {
    // start = 0; lines.length = 913 based on the file 1fas.cube (from line 7 to line 919)
    expect(atoms[0].length).toBe(913);
  });

  test('All atom hetflag are equal to true', () => {
    for (let i = 0; i < len; i++) {
      expect(atoms[0][0].hetflag).toBe(true);
    }
  });

  test('Atom serial is equal to atom index', () => {
    for (let i = 0; i < len; i++) {
      expect(atoms[0][i].serial).toBe(i);
    }
  });

  test('Elem of atom 1 is H', () => {
    expect(atoms[0][0].elem).toBe('H');
  });

  //   test('Atom 1 is bonded to 8,9,11', () => {
  //     //for(let i = 0; i < len; i++){
  //     expect(atoms[0][0].bonds).toEqual([8, 9, 11]);
  //     //}
  //   });

  test('Bond order of atom 1 is ', () => {
    //for(let i = 0; i < len; i++){
    expect(atoms[0][0].bondOrder).toEqual([1, 1, 1]);
    //}
  });

  test('Properties of all atoms are empty', () => {
    for (let i = 0; i < len; i++) {
      expect(atoms[0][i].properties).toEqual({});
    }
  });

  test('x of atom 1', () => {
    //for(let i = 0; i < len; i++){
    expect(atoms[0][0].x).toBeCloseTo(46.148);
    //}
  });

  test('y of atom 1', () => {
    //for(let i = 0; i < len; i++){
    expect(atoms[0][0].y).toBeCloseTo(16.581);
    //}
  });

  test('z of atom 1', () => {
    //for(let i = 0; i < len; i++){
    expect(atoms[0][0].z).toBeCloseTo(2.104);
    //}
  });
});

describe('function XYZ', () => {
  const data = fs.readFileSync('tests/auto/data/md.xyz', 'utf-8');
  let atoms = $3Dmol.Parsers.XYZ(data, {}); //assignbonds

  test('XYZ input', () => {
    expect(atoms).toBeDefined();
  });
  //may need to test onemol
});
/*
describe('function parseV2000', ()=>{
    const data = fs.readFileSync('tests/auto/data/v2000.mol', 'utf-8')
    let atoms = $3Dmol.Parsers.parseV2000(data);
    
    test("parseV2000 input", ()=>{
        expect(atoms).toBeDefined();
    });
});

describe('function parseV3000', ()=>{
    const data = fs.readFileSync('tests/auto/data/v3000.mol', 'utf-8')
    let atoms = $3Dmol.Parsers.parseV3000(data);
    
    test("parseV3000 input", ()=>{
        expect(atoms).toBeDefined();
    });
});
*/

describe('function SDF', () => {
  const data = fs.readFileSync('tests/test_structs/aromaticsdf.sdf', 'utf-8');
  let atoms = $3Dmol.Parsers.SDF(data, {});

  test('SDF input', () => {
    expect(atoms).toBeDefined();
  });
});

//describe('function json', () => {});
/*
describe('function cif', ()=>{
    const data = fs.readFileSync('tests/auto/data/1jpy.cif', 'utf-8')
    let atoms = $3Dmol.Parsers.cif(data, {});
    
    test("cif input", ()=>{
        expect(atoms).toBeDefined();
    });
});
*/
describe('function MOL2', () => {
  const data = fs.readFileSync('tests/test_structs/multiple.mol2', 'utf-8');
  let atoms = $3Dmol.Parsers.MOL2(data, {});

  test('MOL2 input', () => {
    expect(atoms).toBeDefined();
  });
});

describe('function PDBQT', () => {
  const data = fs.readFileSync('tests/test_structs/1B5S.pdb', 'utf-8');
  let atoms = $3Dmol.Parsers.PDBQT(data, {});

  test('PDBQT input', () => {
    expect(atoms).toBeDefined();
  });
});

describe('function PQR', () => {
  const data = fs.readFileSync('tests/auto/data/1fas.pqr', 'utf-8');
  let atoms = $3Dmol.Parsers.PQR(data, {});

  test('PQR input', () => {
    expect(atoms).toBeDefined();
  });
});
/*
describe('function MMTF', ()=>{
    const data = fs.readFileSync('tests/auto/data/', 'utf-8')
    let atoms = $3Dmol.Parsers.parseV3000(data);
    
    test("mmtf input", ()=>{
        expect(atoms).toBeDefined();
    });
});
*/
describe('function PRMTOP', () => {
  const data = fs.readFileSync('tests/auto/data/model1.prmtop', 'utf-8');
  let atoms = $3Dmol.Parsers.PRMTOP(data);

  test('PRMTOP input', () => {
    expect(atoms).toBeDefined();
  });
});
