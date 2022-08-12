// parse SYBYL mol2 file from string - assumed to only contain one molecule
// tag
/**
 * @param {string}
 *            str
 * @param {ParserOptionsSpec}
 *            options
 */
export function MOL2(str, options) {
  var atoms: any[][] & Record<string,any> = [[]];
  var noH = false;
  if (typeof options.keepH !== "undefined") noH = !options.keepH;

  // Note: these regex's work, though they don't match '<TRIPOS>'
  // correctly - something to do with angle brackets
  var mol_pos = str.search(/@<TRIPOS>MOLECULE/);
  var atom_pos = str.search(/@<TRIPOS>ATOM/);

  // Assuming both Molecule and Atom sections exist
  if (mol_pos == -1 || atom_pos == -1) return atoms;

  var lines = str.substr(mol_pos, str.length).split(/\r?\n|\r/);
  while (lines.length > 0) {
    // serial is atom's index in file; index is atoms index in 'atoms'
    var serialToIndex: number[] = [];
    var tokens = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
    var natoms = parseInt(tokens[0]);
    var nbonds = 0;

    if (tokens.length > 1) nbonds = parseInt(tokens[1]);

    var offset = 4;
    var i;
    // Continue until 'Atom' section
    for (i = 3; i < lines.length; i++) {
      if (lines[i] == "@<TRIPOS>ATOM") {
        offset = i + 1;
        break;
      }
    }

    var start = atoms[atoms.length - 1].length;
    var end = start + natoms;
    var line;
    // Process ATOMS
    for (i = start; i < end; i++) {
      line = lines[offset++];
      tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
      var atom: Record<string, any> = {};
      // get element
      var elem = tokens[5].split(".")[0];
      atom.atom = atom.elem =
        elem[0].toUpperCase() + elem.substr(1).toLowerCase();
      if (atom.elem == "H" && noH) {
        // ignore
      } else {
        // 'index' is this atom's index in 'atoms'; 'serial' is this
        // atom's
        // serial id in mol2 file
        var index = atoms[atoms.length - 1].length;
        var serial = parseInt(tokens[0]);
        atom.serial = serial;
        // atom.serial = i;

        atom.x = parseFloat(tokens[2]);
        atom.y = parseFloat(tokens[3]);
        atom.z = parseFloat(tokens[4]);
        atom.atom = tokens[5];
        var charge = parseFloat(tokens[8]);

        atom.index = index;
        atom.bonds = [];
        atom.bondOrder = [];
        atom.properties = {
          charge: charge,
          partialCharge: charge,
        };
        serialToIndex[serial] = index;

        atoms[atoms.length - 1].push(atom);
      }
    }

    // Process BONDS
    var bonds_found = false;
    while (offset < lines.length) {
      if (lines[offset++] == "@<TRIPOS>BOND") {
        bonds_found = true;
        break;
      }
    }

    if (bonds_found && nbonds) {
      for (i = 0; i < nbonds; i++) {
        line = lines[offset++];

        tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        var from = parseInt(tokens[1]);
        var fromAtom = atoms[atoms.length - 1][serialToIndex[from]];
        var to = parseInt(tokens[2]);
        var toAtom = atoms[atoms.length - 1][serialToIndex[to]];

        // Won't be able to read aromatic bonds correctly...
        var order = parseInt(tokens[3]);
        if (isNaN(order)) order = 1;

        if (fromAtom !== undefined && toAtom !== undefined) {
          fromAtom.bonds.push(serialToIndex[to]);
          fromAtom.bondOrder.push(order);
          toAtom.bonds.push(serialToIndex[from]);
          toAtom.bondOrder.push(order);
        }
      }
    }
    if (options.multimodel) {
      if (!options.onemol) atoms.push([]);
      lines.splice(0, offset);
      str = lines.join("\n"); //update for str.search
      continue;
    } else {
      break;
    }
  }
  return atoms;
}
