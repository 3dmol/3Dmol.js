  // read an XYZ file from str and return result

import { Matrix3 } from "../WebGL/math";
import assignBonds from "./util/assignBonds";

  /**
   * @param {string} str
   * @param {import("../specs").ParserOptionsSpec} options
   */
export default function parseXYZ(str, options) {
    options = options || {};
      /** @type {import("../specs").ParserResult} */

    let atoms = [[]];
    const assignbonds = options.assignBonds === undefined ? true : options.assignBonds;
    const lines = str.split(/\r?\n|\r/);
    while (lines.length > 0) {
      if (lines.length < 3)
        break;
      const atomCount = parseInt(lines[0]);
      if (isNaN(atomCount) || atomCount <= 0)
        break;
      if (lines.length < atomCount + 2)
        break;

      const latticeRe = /Lattice\s*=\s*["\{\}]([^"\{\}]+)["\{\}]\s*/gi;
      const latticeMatch = latticeRe.exec(lines[1]);
      if (latticeMatch != null && latticeMatch.length > 1) {
        const lattice = new Float32Array(latticeMatch[1].split(/\s+/).map(Number.parseFloat));
        const matrix = new Matrix3(
          lattice[0],
          lattice[3],
          lattice[6],
          lattice[1],
          lattice[4],
          lattice[7],
          lattice[2],
          lattice[5],
          lattice[8]
        );
        atoms.modelData = [{ cryst: { matrix } }];
      }

      let offset = 2;
      const start = atoms[atoms.length - 1].length;
      const end = start + atomCount;
      for (let i = start; i < end; i++) {
        const line = lines[offset++];
        const tokens = line.replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
        const atom = {};
        atom.serial = i;
        const elem = tokens[0];
        atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1, 1).toLowerCase();
        atom.x = parseFloat(tokens[1]);
        atom.y = parseFloat(tokens[2]);
        atom.z = parseFloat(tokens[3]);
        atom.hetflag = true;
        atom.bonds = [];
        atom.bondOrder = [];
        atom.properties = {};
        atoms[atoms.length - 1][i] = atom;
        if (tokens.length >= 7) {
          atom.dx = parseFloat(tokens[4]);
          atom.dy = parseFloat(tokens[5]);
          atom.dz = parseFloat(tokens[6]);
        }
      }

      if (options.multimodel) {
        atoms.push([]);
        lines.splice(0, offset);
      } else {
        break;
      }
    }

    if (assignbonds) {
      for (let i = 0; i < atoms.length; i++) {
        assignBonds(atoms[i]);
      }
    }

    if (options.onemol) {
      const temp = atoms;
      atoms = [];
      atoms.push(temp[0]);
      for (let i = 1; i < temp.length; i++) {
        const offset = atoms[0].length;
        for (let j = 0; j < temp[i].length; j++) {
          const a = temp[i][j];
          if (!a.bonds) a.bonds = [];
          for (let k = 0; k < a.bonds.length; k++) {
            a.bonds[k] += offset;
          }
          a.index = atoms[0].length;
          a.serial = atoms[0].length;
          atoms[0].push(a);
        }
      }
    }

    return atoms;
  };