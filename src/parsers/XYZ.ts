import { ParserOptionsSpec } from './ParserOptionsSpec';
// read an XYZ file from str and return result

import { Matrix3 } from "../WebGL";
import { assignBonds } from "./utils/assignBonds";

/**
 * @param {string} str
 * @param {ParserOptionsSpec} options
 * @category Parsers
 */
export function XYZ(str: string, options: ParserOptionsSpec) {
  options = options || {};
  var atoms: any[][] & Record<string, any> = [[]];
  var assignbonds =
    options.assignBonds === undefined ? true : options.assignBonds;
  var lines = str.split(/\r?\n|\r/);
  while (lines.length > 0) {
    if (lines.length < 3) break;
    var atomCount = parseInt(lines[0]);
    if (isNaN(atomCount) || atomCount <= 0) break;
    if (lines.length < atomCount + 2) break;

    var lattice_re = /Lattice\s*=\s*["\{\}]([^"\{\}]+)["\{\}]\s*/gi;
    var lattice_match = lattice_re.exec(lines[1]);
    if (lattice_match != null && lattice_match.length > 1) {
      var lattice = new Float32Array(lattice_match[1].split(/\s+/) as any);
      var matrix = new Matrix3(
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
      atoms.modelData = [{ cryst: { matrix: matrix } }];
    }

    var offset = 2;
    var start = atoms[atoms.length - 1].length;
    var end = start + atomCount;
    for (var i = start; i < end; i++) {
      var line = lines[offset++];
      var tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
      var atom: Record<string, any> = {};
      atom.serial = i;
      var elem = tokens[0];
      atom.atom = atom.elem =
        elem[0].toUpperCase() + elem.substr(1, 1).toLowerCase();
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
    var temp = atoms;
    atoms = [];
    atoms.push(temp[0]);
    for (let i = 1; i < temp.length; i++) {
      let offset = atoms[0].length;
      for (let j = 0; j < temp[i].length; j++) {
        let a = temp[i][j];
        for (let k = 0; k < a.bonds.length; k++) {
          a.bonds[k] = a.bonds[k] + offset;
        }
        a.index = atoms[0].length;
        a.serial = atoms[0].length;
        atoms[0].push(a);
      }
    }
  }

  return atoms;
}
