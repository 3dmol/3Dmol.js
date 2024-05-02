import { AtomSpec } from "specs";
import { ParserOptionsSpec } from "./ParserOptionsSpec";
import { assignPDBBonds } from "./utils/assignPDBBonds";
import { atomNameToElem } from "./utils/atomNameToElem";

/**
 * Parse a gro file from str and create atoms
 *
 * @param {string} str
 * @param {ParserOptionsSpec} options
 * @category Parsers
 * @returns {Array<AtomSpec[]>} - Returns a 2D array of type AtomSpec
 */

export function GRO(str: string, options: ParserOptionsSpec) {
  const allatoms: AtomSpec[][] & { box?: string[] } = [];
  const lines = str.split(/\r?\n|\r/);
  while (lines.length > 0) {
    const atomCount = parseInt(lines[1]);
    const breakCondition =
      lines.length < 3 ||
      isNaN(atomCount) ||
      atomCount <= 0 ||
      lines.length < atomCount + 3;
    if (breakCondition) break;
    const atoms: AtomSpec[] = [];
    allatoms.push(atoms);
    let offset = 2;
    const start = atoms.length;
    const end = start + atomCount;
    for (let i = start; i < end; i++) {
      const line = lines[offset++];
      const atom: AtomSpec = {};
      atom.serial = i;
      atom.atom = line.slice(10, 15).trim();
      atom.elem = atomNameToElem(atom.atom, true);
      //coordinates are in nM, convert to A
      atom.x = 10.0 * parseFloat(line.slice(20, 28));
      atom.y = 10.0 * parseFloat(line.slice(28, 36));
      atom.z = 10.0 * parseFloat(line.slice(36, 44));
      atom.resi = parseInt(line.slice(0, 5));
      atom.resn = line.slice(5, 10).trim();
      atom.bonds = [];
      atom.bondOrder = [];
      atom.properties = {};
      if (line.length > 44) {
        atom.dx = 10.0 * parseFloat(line.slice(44, 52));
        atom.dy = 10.0 * parseFloat(line.slice(52, 60));
        atom.dz = 10.0 * parseFloat(line.slice(60, 68));
      }
      atoms[i] = atom;
    } //for all atoms

    if (lines.length <= offset + 3) {
      //single line left, assume it is the box
      const last = lines[offset++];
      const box = last.trim().split(/\s+/);
      if (box.length === 3) {
        for (let b = 0; b < 3; b++) {
          box[b] = (parseFloat(box[b]) * 10.0).toString();
        }
        allatoms.box = box;
      }
    }
    lines.splice(0, ++offset);
  }

  for (let i = 0; i < allatoms.length; i++) {
    assignPDBBonds(allatoms[i], options);
  }
  return allatoms;
}
