/* eslint-disable radix */

import assignPDBBonds from "./util/assignPDBBonds";
import computeSecondaryStructure from "./util/computeSecondaryStructure";

/**
 * Parse a pqr file from str and create atoms. A pqr file is assumed to be a
 * whitespace delimited PDB with charge and radius fields.
 *
 * @param {string} str
 * @param {import("../specs").ParserOptionsSpec} options - noSecondaryStructure (do not compute ss)
 */
export default function parsePRQ(str, options) {
  /** @type {import("../specs").ParserResult} */
  const atoms = [[]];
  const computeStruct = !options.noSecondaryStructure;
  atoms.modelData = [{symmetries: []}];
  const serialToIndex = []; // map from pdb serial to index in atoms
  const lines = str.split(/\r?\n|\r/);
  let line;
  for (let i = 0; i < lines.length; i++) {
    line = lines[i].replace(/^\s*/, ''); // remove indent
    const recordName = line.substring(0, 6);

    if (recordName.indexOf('END') === 0) {
      if (options.multimodel) {
        if (!options.onemol) atoms.push([]);
        continue;
      } else {
        break;
      }
    } else if (recordName === 'ATOM  ' || recordName === 'HETATM') {
      // I would have liked to split based solely on whitespace, but
      // it seems that there is no guarantee that all the fields will
      // be filled out (e.g. the chain) so this doesn't work
      let hetflag;
      const serial = parseInt(line.substring(6, 5));
      const atom = line.substring(12, 4).replace(/ /g, '');
      const resn = line.substring(17, 3).trim();
      const chain = line.substring(21, 1);
      const resi = parseInt(line.substring(22, 4));
      // however let's split the coordinates, charge and radius by
      // whitespace
      // to support extra precision
      const vals = line.substring(30).trim().split(/\s+/);
      const x = parseFloat(vals[0]);
      const y = parseFloat(vals[1]);
      const z = parseFloat(vals[2]);
      const charge = parseFloat(vals[3]);
      const radius = parseFloat(vals[4]);

      let elem = atom[0];
      if (atom.length > 1 && atom[1].toUpperCase() !== atom[1]) {
        // slight hack - identify two character elements by the
        // second character in the atom name being lowercase
        elem = atom.substring(0, 2);
      }

      if (line[0] === 'H') hetflag = true;
      else hetflag = false;
      serialToIndex[serial] = atoms[atoms.length - 1].length;
      atoms[atoms.length - 1].push({
        resn,
        x,
        y,
        z,
        elem,
        hetflag,
        chain,
        resi,
        serial,
        atom,
        bonds: [],
        ss: 'c',
        bondOrder: [],
        properties: {
          charge,
          partialCharge: charge,
          radius,
        },
        pdbline: line,
      });
    } else if (recordName === 'CONECT') {
      // MEMO: We don't have to parse SSBOND, LINK because both are
      // also
      // described in CONECT. But what about 2JYT???
      const from = parseInt(line.substring(6, 5));
      const fromAtom = atoms[atoms.length - 1][serialToIndex[from]];
      for (let j = 0; j < 4; j++) {
        const to = parseInt(line.substring([11, 16, 21, 26][j], 5));
        const toAtom = atoms[atoms.length - 1][serialToIndex[to]];
        if (fromAtom && toAtom) {
          if (!fromAtom.bonds) fromAtom.bonds = [];
          if (!fromAtom.bondOrder) fromAtom.bondOrder = [];
          fromAtom.bonds.push(serialToIndex[to]);
          fromAtom.bondOrder.push(1);
        }
      }
    }
  }

  // assign bonds - yuck, can't count on connect records
  for (let i = 0; i < atoms.length; i++) {
    assignPDBBonds(/** @type {import("../specs").AtomSpec[]} */(atoms[i]));
    if (computeStruct) computeSecondaryStructure(atoms[i]);
  }

  return atoms;
};
