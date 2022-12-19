import { assignPDBBonds } from "./utils/assignPDBBonds";
import { atomNameToElem } from "./utils/atomNameToElem";

/**
 * Parse a gro file from str and create atoms
  * 
  * @param {string}
  *            str
  * @param {ParserOptionsSpec}
  *            options* 
  * @category Parsers

 */
export function GRO(str /*, options*/) {
  var allatoms: any[][] & Record<string, any> = [];
  var lines = str.split(/\r?\n|\r/);
  while (lines.length > 0) {
    if (lines.length < 3) break;
    var atomCount = parseInt(lines[1]);
    if (isNaN(atomCount) || atomCount <= 0) break;
    if (lines.length < atomCount + 3) break;
    var atoms: any[] = [];
    allatoms.push(atoms);
    var offset = 2;
    var start = atoms.length;
    var end = start + atomCount;
    for (var i = start; i < end; i++) {
      var line = lines[offset++];
      var atom: Record<string, any> = {};
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
      var last = lines[offset++];
      var box = last.trim().split(/\s+/);
      if (box.length == 3) {
        for (var b = 0; b < 3; b++) {
          box[b] = parseFloat(box[b]) * 10.0;
        }
        allatoms.box = box;
      }
    }
    lines.splice(0, ++offset);
  }

  for (let i = 0; i < allatoms.length; i++) {
    assignPDBBonds(allatoms[i]);
  }
  return allatoms;
}
