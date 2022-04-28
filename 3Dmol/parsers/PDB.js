/* eslint-disable radix */
import {Matrix4} from '../WebGL/math';
import assignPDBBonds from './util/assignPDBBonds';
import atomNameToElem from './util/atomNameToElem';
import bondTable from './util/bondTable';
import computeSecondaryStructure from './util/computeSecondaryStructure';
import processSymmetries from './util/processSymmeties';

const isEmpty = obj => {
  let name;
  for (name in obj) {
    return false;
  }
  return true;
};

/**
 * // make sure bonds are actually two way
 * @param {import('../specs').AtomSpec[]} atomsarray
 * @param {number[]} serialToIndex
 */
const validateBonds = (atomsarray, serialToIndex) => {
  for (let i = 0, n = atomsarray.length; i < n; i++) {
    const atom = atomsarray[i];
    for (let b = 0; b < atom.bonds.length; b++) {
      const a2i = atom.bonds[b];
      const atom2 = atomsarray[a2i];
      const atomi = serialToIndex[atom.serial || 0];
      if (atom2 && atomi) {
        const a1i = atom2.bonds.indexOf(atomi);
        if (a1i < 0) {
          atom2.bonds.push(atomi);
          atom2.bondOrder.push(atom.bondOrder[b]);
        }
      }
    }
  }
};

// return one model worth of pdb, returns atoms, modelData, and remaining lines
/**
 *
 * @param {string[]} lines
 * @param {import("../specs").ParserOptionsSpec} options
 * @param {any} sslookup
* @returns {[import('../specs').AtomSpec[], {symmetries: Matrix4[];cryst?: Record<string, number> | undefined; matrix?: import("../WebGL/math").Matrix3}, any[]]}
 */
const getSinglePDB = (lines, options, sslookup) => {
  /** @type {Array<Partial<import("../specs").AtomSpec>>} */
  const atoms = [];
  const assignbonds = options.assignBonds === undefined ? true : options.assignBonds;
  const noH = !options.keepH; // suppress hydrogens by default
  const ignoreStruct = !!options.noSecondaryStructure;
  const computeStruct = !options.noComputeSecondaryStructure;
  const noAssembly = !options.doAssembly; // don't assemble by default
  const selAltLoc = options.altLoc ? options.altLoc : 'A'; // default alternate location to select if present
  /** @type {{symmetries: Matrix4[], cryst?: Record<string,number>}} */
  const modelData = {symmetries: []};
  let atom;
  let remainingLines = [];
  const serialToIndex = []; // map from pdb serial to index in atoms
  const seenbonds = {}; // sometimes connect records are duplicated as an unofficial means of relaying bond orders

  for (let i = 0; i < lines.length; i++) {
    let line = lines[i].replace(/^\s*/, ''); // remove indent
    const recordName = line.substring(0, 6);
    let startChain;
    let startResi;
    let endResi;

    if (recordName.indexOf('END') === 0) {
      remainingLines = lines.slice(i + 1);
      if (recordName === 'END') {
        // as opposed to ENDMDL
        // reset secondary structure
        for (const prop in sslookup) {
          if (sslookup[prop]) {
            delete sslookup[prop];
          }
        }
      }
      break;
    } else if (recordName === 'ATOM  ' || recordName === 'HETATM') {
      let hetflag;
      let elem;
      const altLoc = line.substring(16, 1);
      if (altLoc !== ' ' && altLoc !== selAltLoc && selAltLoc !== '*') continue;
      
      const serial = parseInt(line.substring(6, 5));
      atom = line.substring(12, 4).replace(/ /g, '');
      const resn = line.substring(17, 3).replace(/ /g, '');
      const chain = line.substring(21, 1);
      
      const resi = parseInt(line.substring(22, 4));
      const icode = line.substring(26, 1);
      const x = parseFloat(line.substring(30, 8));
      const y = parseFloat(line.substring(38, 8));
      const z = parseFloat(line.substring(46, 8));
      const b = parseFloat(line.substring(60, 8));
      elem = line.substring(76, 2).replace(/ /g, '');
      if (elem === '' || typeof bondTable[elem] === 'undefined') {
        // for some incorrect PDB files
        elem = atomNameToElem(line.substring(12, 2), line[0] === 'A');
      } else {
        elem = elem[0].toUpperCase() + elem.substring(1).toLowerCase();
      }

      if (elem === 'H' && noH) continue;
      if (recordName[0] === 'H') hetflag = true;
      else hetflag = false;
      serialToIndex[serial] = atoms.length;
      atoms.push({
        resn,
        x,
        y,
        z,
        elem,
        hetflag,
        altLoc,
        chain,
        resi,
        icode,
        rescode: resi + (icode !== ' ' ? `^${icode}` : ''),

        // resi
        // and
        // icode
        serial,
        atom,
        bonds: [],
        ss: 'c',
        bondOrder: [],
        properties: {},
        b,
        pdbline: line,
      });
    } else if (recordName === 'SHEET ') {
      startChain = line.substring(21, 1);
      startResi = parseInt(line.substring(22, 4));
      endResi = parseInt(line.substring(33, 4));
      if (!(startChain in sslookup)) {
        sslookup[startChain] = {};
      }
      // mark start and end with additional character
      sslookup[startChain][startResi] = 's1';
      for (let res = startResi + 1; res < endResi; res++) {
        sslookup[startChain][res] = 's';
      }
      sslookup[startChain][endResi] = 's2';
    } else if (recordName === 'CONECT') {
      // MEMO: We don't have to parse SSBOND, LINK because both are
      // also
      // described in CONECT. But what about 2JYT???
      
      const from = parseInt(line.substring(6, 5));
      const fromindex = serialToIndex[from];
      const fromAtom = atoms[fromindex];
      const coffsets = [11, 16, 21, 26];
      for (let j = 0; j < 4; j++) {
        
        const to = parseInt(line.substring(coffsets[j], 5));
        const toindex = serialToIndex[to];
        const toAtom = atoms[toindex] || {};
        if (!fromAtom.bonds) fromAtom.bonds = [];
        if (!fromAtom.bondOrder) fromAtom.bondOrder = [];
        if (fromAtom !== undefined && toAtom !== undefined) {
          // duplicated conect records indicate bond order
          if (!seenbonds[[fromindex, toindex]]) {
            seenbonds[[fromindex, toindex]] = 1;
            if (
              fromAtom.bonds.length === 0 ||
              fromAtom.bonds[fromAtom.bonds.length - 1] !== toindex
            ) {
              fromAtom.bonds.push(toindex);
              fromAtom.bondOrder.push(1);
            }
          } else {
            // update bond order
            seenbonds[[fromindex, toindex]] += 1;

            for (let bi = 0; bi < fromAtom.bonds.length; bi++) {
              if (fromAtom.bonds[bi] === toindex) {
                const newbo = seenbonds[[fromindex, toindex]];
                if (newbo >= 4) {
                  // aromatic
                  fromAtom.bondOrder[bi] = 1;
                } else {
                  fromAtom.bondOrder[bi] = newbo;
                }
              }
            }
          }
        }
      }
    } else if (recordName === 'HELIX ') {
      startChain = line.substring(19, 1);
      startResi = parseInt(line.substring(21, 4));
      endResi = parseInt(line.substring(33, 4));

      if (!(startChain in sslookup)) {
        sslookup[startChain] = {};
      }

      sslookup[startChain][startResi] = 'h1';

      for (let res = startResi + 1; res < endResi; res++) {
        sslookup[startChain][res] = 'h';
      }

      sslookup[startChain][endResi] = 'h2';
    } else if (!noAssembly && recordName === 'REMARK' && line.substring(13, 5) === 'BIOMT') {
      let n;
      const matrix = new Matrix4();
      for (n = 1; n <= 3; n++) {
        line = lines[i].replace(/^\s*/, '');
        if (parseInt(line.substring(18, 1)) === n) {
          // check for all
          // three lines
          // by matching #
          // @ end of
          // "BIOMT" to n
          matrix.elements[n - 1] = parseFloat(line.substring(23, 10));
          matrix.elements[n - 1 + 4] = parseFloat(line.substring(33, 10));
          matrix.elements[n - 1 + 8] = parseFloat(line.substring(43, 10));
          matrix.elements[n - 1 + 12] = parseFloat(line.substring(53));
          i++;
        } else {
          while (line.substring(13, 5) === 'BIOMT') {
            i++;
            line = lines[i].replace(/^\s*/, '');
          }
        }
      }
      matrix.elements[3] = 0;
      matrix.elements[7] = 0;
      matrix.elements[11] = 0;
      matrix.elements[15] = 1;
      modelData.symmetries.push(matrix);
      i--; // set i back
    } else if (recordName === 'CRYST1') {
      const a = parseFloat(line.substring(7, 8));
      const b = parseFloat(line.substring(16, 8));
      const c = parseFloat(line.substring(25, 8));
      const alpha = parseFloat(line.substring(34, 6));
      const beta = parseFloat(line.substring(41, 6));
      const gamma = parseFloat(line.substring(48, 6));
      modelData.cryst = {a, b, c, alpha, beta, gamma};
    } else if (recordName === 'ANISOU') {
      const serial = parseInt(line.substring(6, 5));
      const anisouAtomIndex = serialToIndex[serial];
      const anisouAtom = atoms[anisouAtomIndex];
      if (anisouAtom) {
        const vals = line.substring(30).trim().split(/\s+/);
        const uMat = {
          u11: parseInt(vals[0]),
          u22: parseInt(vals[1]),
          u33: parseInt(vals[2]),
          u12: parseInt(vals[3]),
          u13: parseInt(vals[4]),
          u23: parseInt(vals[5]),
        };

        anisouAtom.uMat = uMat;
      }
    }
  }

  // fix any "one-way" bonds in CONECT records
  validateBonds(/** @type {import('../specs').AtomSpec[]} */ (atoms), serialToIndex);
  // assign bonds - yuck, can't count on connect records
  if (assignbonds) assignPDBBonds(/** @type {import('../specs').AtomSpec[]} */ (atoms));

  if (!noAssembly) processSymmetries(modelData.symmetries, atoms, options, modelData.cryst);

  if (computeStruct && !ignoreStruct) {
    computeSecondaryStructure(atoms);
  }

  // Assign secondary structures from pdb file
  if (!isEmpty(sslookup)) {
    for (let i = 0; i < atoms.length; i++) {
      atom = atoms[i];
      if (atom === undefined) continue;
      if (atom.chain && atom.chain in sslookup && atom.resi && atom.resi in sslookup[atom.chain]) {
        const code = sslookup[atom.chain][atom.resi];
        atom.ss = code[0];
        if (code.length > 1) {
          if (code[1] === '1') atom.ssbegin = true;
          else if (code[1] === '2') atom.ssend = true;
        }
      }
    }
  }
  // console.log("assign structure " + ((new Date()).getTime() - starttime));
  return [/** @type {import('../specs').AtomSpec[]} */(atoms), modelData, remainingLines];
};

// parse pdb file from str and create atoms
// if computeStruct is true will always perform secondary structure
// analysis,
// otherwise only do analysis of SHEET/HELIX comments are missing
/**
 * @param {string} str
 * @param {import("../specs").ParserOptionsSpec} options - keepH (do not strip hydrogens), noSecondaryStructure,
 *            assignbonds (default true, calculate implicit bonds)
 *            (do not compute ss), altLoc (which alternate location to select, if present; '*' to load all)
 */
export default function parsePDB(str, options) {
  options = options || {};
  /** @type {import("../specs").ParserResult} */
  const atoms = []; // a separate list for each model
  const sslookup = {}; // stores SHEET and HELIX info, which is shared across models
  atoms.modelData = [];
  let lines = str.split(/\r?\n|\r/);
  while (lines.length > 0) {
    const pdbinfo = getSinglePDB(lines, options, sslookup);
    const modelatoms = pdbinfo[0];
    const modelData = pdbinfo[1];
    lines = pdbinfo[2];

    if (modelatoms.length === 0) {
      continue; // happens when there are blank lines
    }
    if (options.multimodel && options.onemol && atoms.length > 0) {
      // merge into existing atoms
      const inc = atoms[0].length;
      for (let i = 0; i < modelatoms.length; i++) {
        // renumber
        const atom = modelatoms[i];
        atom.index = i;
        for (let b = 0; b < atom.bonds.length; b++) {
          atom.bonds[b] += inc;
        }
        atoms[0].push(atom);
      }
    } else {
      atoms.modelData.push(modelData);
      atoms.push(modelatoms);
    }

    if (!options.multimodel) {
      break;
    }
  }

  return atoms;
}
