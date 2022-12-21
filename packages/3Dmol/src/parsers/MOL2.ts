
let SYBYLtoElem = {
  'C.1': 'C',
  'C1': 'C',
  'C.2': 'C',
  'C2': 'C',
  'C.3': 'C',
  'C3': 'C',  
  'C.ar': 'C',
  'Car': 'C',
  'C.cat': 'C',
  'Ccat': 'C',
  'H.spc' :'H',
  'Hspc':'H',
  'H.t3p':'H',
  'Ht3p': 'H',
  'N.1':'N',
  'N1':'N',
  'N.2':'N',
  'N2':'N',
  'N.3':'N',
  'N3':'N',
  'N.4':'N',
  'N4':'N',
  'N.am':'N',
  'Nam':'N',
  'N.ar':'N',
  'Nar':'N',
  'N.p13':'N',
  'Np13':'N',    
  'O.2':'O',
  'O2':'O',
  'O.3':'O',
  'O3':'O',
  'O.co2':'O',
  'Oco2':'O',
  'O.spc':'O',
  'Ospc':'O',    
  'O.t3p':'O',
  'Ot3p':'O',  
  'P.3':'P',
  'P3':'P',
  'S.2':'S',
  'S2':'S',  
  'S.3':'S',
  'S3':'S',  
  'S.o':'S',
  'So':'S',  
  'S.o2':'S',
  'So2':'S'
};

// parse SYBYL mol2 file from string - assumed to only contain one molecule
// tag
/**
 * @param {string}
 *            str
 * @param {ParserOptionsSpec}
 *            options
 * @category Parsers
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
    var i: number;
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
      var elem = tokens[5];
      if(SYBYLtoElem[elem] !== undefined) {        
        elem = SYBYLtoElem[elem];
      } else {
        elem = elem.split(".")[0];
        elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();
      }

      atom.atom = tokens[1];
      atom.elem = elem;
        
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
