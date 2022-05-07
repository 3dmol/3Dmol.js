// This parses the ChemDoodle json file format. Although this is registered
// for the json file extension, other chemical json file formats exist that
// this can not parse. Check which one you have and do not assume that
// .json can be parsed
/**
 * 
 * @param {Record<string, any>} str 
 * @param {import("../specs").ParserOptionsSpec} options 
 * @returns {import("../specs").ParserResult}
 */
export default function parseCDJSON(str, options) {
  /** @type {import("../specs").ParserResult} */
  const atoms = [[]];
  if (typeof str === 'string') {
    // Str is usually automatically parsed by JQuery
    str = JSON.parse(str);
  }
  const molecules = str.m;
  const atomsInFile = molecules[0].a; // Assumes there is at least one
  const bondsInFile = molecules[0].b; // molecule and ignores any more

  // Ignores any shapes
  const styles = molecules[0].s;
  const parseStyle =
    options !== undefined && options.parseStyle !== undefined
      ? options.parseStyle
      : styles !== undefined;

  const offset = atoms[atoms.length - 1].length; // When adding atoms their index will be

  // Offset by the number of existing atoms
  for (let i = 0; i < atomsInFile.length; i++) {
    const currentAtom = atomsInFile[i];
    const atom = {};
    atom.id = currentAtom.i; // Probably won't exist. Doesn't seem to

    // break anything.
    atom.x = currentAtom.x;
    atom.y = currentAtom.y;
    atom.z = currentAtom.z || 0; // Default value if file is 2D

    atom.bonds = [];
    atom.bondOrder = [];

    const elem = currentAtom.l || 'C';
    atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();

    atom.serial = atoms[atoms.length - 1].length;
    if (parseStyle) {
      atom.style = styles[currentAtom.s || 0];
    }
    atoms[atoms.length - 1].push(atom);
  }
  for (let i = 0; i < bondsInFile.length; i++) {
    const currentBond = bondsInFile[i];
    const beginIndex = currentBond.b + offset;
    const endIndex = currentBond.e + offset;
    const bondOrder = currentBond.o || 1;

    const firstAtom = atoms[atoms.length - 1][beginIndex];
    const secondAtom = atoms[atoms.length - 1][endIndex];

    if (!firstAtom.bonds) firstAtom.bonds = [];
    firstAtom.bonds.push(endIndex);
    if (!firstAtom.bondOrder) firstAtom.bondOrder = [];
    firstAtom.bondOrder.push(bondOrder);
    if (!secondAtom.bonds) secondAtom.bonds = [];
    secondAtom.bonds.push(beginIndex);
    if (!secondAtom.bondOrder) secondAtom.bondOrder = [];
    secondAtom.bondOrder.push(bondOrder);
  }
  return atoms;
}
