/**
 * @param {!Array.<string>} lines
 * @param {import("../specs").ParserOptionsSpec} options
 * @returns {!Array.<Array<Object>>}
 */
const parseV2000 = (lines, options) => {
  /** @type {import("../specs").ParserResult} */
  const atoms = [[]];
  let noH = false;
  if (typeof options.keepH !== 'undefined') noH = !options.keepH;

  while (lines.length > 0) {
    if (lines.length < 4) break;
    // eslint-disable-next-line radix
    const atomCount = parseInt(lines[3].substring(0, 3));
    if (isNaN(atomCount) || atomCount <= 0) break;
    // eslint-disable-next-line radix
    const bondCount = parseInt(lines[3].substring(3, 3));
    let offset = 4;
    if (lines.length < 4 + atomCount + bondCount) break;

    // serial is atom's index in file; index is atoms index in 'atoms'
    const serialToIndex = [];
    const start = atoms[atoms.length - 1].length;
    const end = start + atomCount;
    let i;
    let line;
    for (i = start; i < end; i++, offset++) {
      line = lines[offset];
      const atom = {};
      const elem = line.substr(31, 3).replace(/ /g, '');
      atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();

      if (atom.elem !== 'H' || !noH) {
        atom.serial = i;
        serialToIndex[i] = atoms[atoms.length - 1].length;
        atom.x = parseFloat(line.substr(0, 10));
        atom.y = parseFloat(line.substr(10, 10));
        atom.z = parseFloat(line.substr(20, 10));
        atom.hetflag = true;
        atom.bonds = [];
        atom.bondOrder = [];
        atom.properties = {};
        atom.index = atoms[atoms.length - 1].length;
        atoms[atoms.length - 1].push(atom);
      }
    }

    for (i = 0; i < bondCount; i++, offset++) {
      line = lines[offset];
      // eslint-disable-next-line radix
      const from = serialToIndex[parseInt(line.substring(0, 3)) - 1 + start];
      // eslint-disable-next-line radix
      const to = serialToIndex[parseInt(line.substring(3, 3)) - 1 + start];
      // eslint-disable-next-line radix
      const order = parseInt(line.substring(6, 3));
      if (typeof from != 'undefined' && typeof to != 'undefined') {
        // @ts-ignore
        atoms[atoms.length - 1][from].bonds.push(to);
        // @ts-ignore
        atoms[atoms.length - 1][from].bondOrder.push(order);
        // @ts-ignore
        atoms[atoms.length - 1][to].bonds.push(from);
        // @ts-ignore
        atoms[atoms.length - 1][to].bondOrder.push(order);
      }
    }
    if (options.multimodel) {
      if (!options.onemol) atoms.push([]);
      while (lines[offset] !== '$$$$') offset++;
      lines.splice(0, ++offset);
    } else {
      break;
    }
  }
  return atoms;
};

/**
 * @param {!Array.<string>} lines
 * @param {import("../specs").ParserOptionsSpec} options
 * @returns {!Array.<!Array<!Object>>}
 */
const parseV3000 = (lines, options) => {
  /** @type {import("../specs").ParserResult} */
  const atoms = [[]];
  let noH = false;
  if (typeof options.keepH !== 'undefined') noH = !options.keepH;

  while (lines.length > 0) {
    if (lines.length < 8) break;

    if (!lines[4].startsWith('M  V30 BEGIN CTAB')) break;
    if (!lines[5].startsWith('M  V30 COUNTS') || lines[5].length < 14) break;

    const counts = lines[5].substring(13).match(/\S+/g);

    if (!counts || counts.length < 2) break;

    // eslint-disable-next-line radix
    const atomCount = parseInt(counts[0]);
    if (Number.isNaN(atomCount) || atomCount <= 0) break;
    // eslint-disable-next-line radix
    const bondCount = parseInt(counts[1]);
    let offset = 7;

    if (lines.length < 8 + atomCount + bondCount)
      // header, bgn+end CTAB, counts, END
      break;

    // serial is atom's index in file; index is atoms index in 'atoms'
    const serialToIndex = [];
    const start = atoms[atoms.length - 1].length;
    const end = start + atomCount;
    let i;
    let line;
    for (i = start; i < end; i++, offset++) {
      line = lines[offset];
      const atomParts = line.substr(6).match(/\S+/g);
      if (atomParts && atomParts.length > 4) {
        const atom = {};
        const elem = atomParts[1].replace(/ /g, '');
        atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();

        if (atom.elem !== 'H' || !noH) {
          atom.serial = i;
          serialToIndex[i] = atoms[atoms.length - 1].length;
          atom.x = parseFloat(atomParts[2]);
          atom.y = parseFloat(atomParts[3]);
          atom.z = parseFloat(atomParts[4]);
          atom.hetflag = true;
          atom.bonds = [];
          atom.bondOrder = [];
          atom.properties = {};
          atom.index = atoms[atoms.length - 1].length;
          atoms[atoms.length - 1].push(atom);
        }
      }
    }

    if (lines[offset] === 'M  V30 END ATOM') offset++;
    else break;

    if (bondCount !== 0 && lines[offset] === 'M  V30 BEGIN BOND') offset++;
    else break;

    for (i = 0; i < bondCount; i++, offset++) {
      line = lines[offset];
      const bondParts = line.substr(6).match(/\S+/g);
      if (bondParts && bondParts.length > 3) {
        const from = serialToIndex[parseInt(bondParts[2]) - 1 + start];
        const to = serialToIndex[parseInt(bondParts[3]) - 1 + start];
        const order = parseInt(bondParts[1]);
        if (typeof from != 'undefined' && typeof to != 'undefined') {
          // @ts-ignore
          atoms[atoms.length - 1][from].bonds.push(to);
          // @ts-ignore
          atoms[atoms.length - 1][from].bondOrder.push(order);
          // @ts-ignore
          atoms[atoms.length - 1][to].bonds.push(from);
          // @ts-ignore
          atoms[atoms.length - 1][to].bondOrder.push(order);
        }
      }
    }
    if (options.multimodel) {
      if (!options.onemol) {
        atoms.push([]);
      }
      while (lines[offset] !== '$$$$') {
        offset++;
      }
      lines.splice(0, ++offset);
    } else {
      break;
    }
  }
  return atoms;
};

// put atoms specified in sdf fromat in str into atoms
// adds to atoms, does not replace
/**
 * @param {string} str
 * @param {import("../specs").ParserOptionsSpec} options
 */
export default function parseSDF(str, options) {
  let molformat = 'V2000';
  const lines = str.split(/\r?\n|\r/);
  if (lines.length > 3 && lines[3].length > 38) {
    molformat = lines[3].substr(34, 5);
  }
  if (molformat === 'V2000') {
    return parseV2000(lines, options);
  }
  if (molformat === 'V3000') {
    return parseV3000(lines, options);
  }
  return [[]];
}
