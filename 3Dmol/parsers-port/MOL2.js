// parse SYBYL mol2 file from string - assumed to only contain one molecule
// tag
/**
 * @param {string} str
 * @param {import("../specs").ParserOptionsSpec} options
 */
export default function parseMOL2(str, options) {
  /** @type {import("../specs").ParserResult} */
  const atoms = [[]];
  let noH = false;
  if (typeof options.keepH !== 'undefined') noH = !options.keepH;

  // Note: these regex's work, though they don't match '<TRIPOS>'
  // correctly - something to do with angle brackets
  const molPos = str.search(/@<TRIPOS>MOLECULE/);
  const atomPos = str.search(/@<TRIPOS>ATOM/);

  // Assuming both Molecule and Atom sections exist
  if (molPos === -1 || atomPos === -1) return atoms;

  const lines = str.substr(molPos, str.length).split(/\r?\n|\r/);
  while (lines.length > 0) {
    // serial is atom's index in file; index is atoms index in 'atoms'
    const serialToIndex = [];
    let tokens = lines[2].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
    // eslint-disable-next-line radix
    const natoms = parseInt(tokens[0]);
    let nbonds = 0;

    // eslint-disable-next-line radix
    if (tokens.length > 1) nbonds = parseInt(tokens[1]);

    let offset = 4;
    let i;
    // Continue until 'Atom' section
    for (i = 3; i < lines.length; i++) {
      if (lines[i] === '@<TRIPOS>ATOM') {
        offset = i + 1;
        break;
      }
    }

    const start = atoms[atoms.length - 1].length;
    const end = start + natoms;
    let line;
    // Process ATOMS
    for (i = start; i < end; i++) {
      line = lines[offset++];
      tokens = line.replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
      const atom = {};
      // get element
      const elem = tokens[5].split('.')[0];
      atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();
      if (atom.elem === 'H' && noH) {
        // ignore
      } else {
        // 'index' is this atom's index in 'atoms'; 'serial' is this
        // atom's
        // serial id in mol2 file
        const index = atoms[atoms.length - 1].length;
        // eslint-disable-next-line radix
        const serial = parseInt(tokens[0]);
        atom.serial = serial;
        // atom.serial = i;
        atom.x = parseFloat(tokens[2]);
        atom.y = parseFloat(tokens[3]);
        atom.z = parseFloat(tokens[4]);
        atom.atom = tokens[5];
        const charge = parseFloat(tokens[8]);

        atom.index = index;
        atom.bonds = [];
        atom.bondOrder = [];
        atom.properties = {
          charge,
          partialCharge: charge,
        };
        serialToIndex[serial] = index;

        atoms[atoms.length - 1].push(atom);
      }
    }

    // Process BONDS
    let bondsFound = false;
    while (offset < lines.length) {
      if (lines[offset++] === '@<TRIPOS>BOND') {
        bondsFound = true;
        break;
      }
    }

    if (bondsFound && nbonds) {
      for (i = 0; i < nbonds; i++) {
        line = lines[offset++];

        tokens = line.replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
        // eslint-disable-next-line radix
        const from = parseInt(tokens[1]);
        const fromAtom = atoms[atoms.length - 1][serialToIndex[from]];
        // eslint-disable-next-line radix
        const to = parseInt(tokens[2]);
        const toAtom = atoms[atoms.length - 1][serialToIndex[to]];

        // Won't be able to read aromatic bonds correctly...
        // eslint-disable-next-line radix
        let order = parseInt(tokens[3]);
        if (Number.isNaN(order)) order = 1;

        if (fromAtom !== undefined && toAtom !== undefined) {
          if (!fromAtom.bonds) fromAtom.bonds = [];
          fromAtom.bonds.push(serialToIndex[to]);
          if (!fromAtom.bondOrder) fromAtom.bondOrder = [];
          fromAtom.bondOrder.push(order);
          if (!toAtom.bonds) toAtom.bonds = [];
          toAtom.bonds.push(serialToIndex[from]);
          if (!toAtom.bondOrder) toAtom.bondOrder = [];
          toAtom.bondOrder.push(order);
        }
      }
    }
    if (options.multimodel) {
      if (!options.onemol) atoms.push([]);
      lines.splice(0, offset);
      str = lines.join('\n'); // update for str.search
      continue;
    } else {
      break;
    }
  }
  return atoms;
}
