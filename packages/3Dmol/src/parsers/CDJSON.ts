/** 
 * This parses the ChemDoodle json file format. Although this is registered
 * for the json file extension, other chemical json file formats exist that
 * this can not parse. Check which one you have and do not assume that
 * .json can be parsed
 * 
 * @param {string} str
 * @param {ParserOptionsSpec} options
 * @category Parsers
 */
export function CDJSON(str, options) {
  var atoms: any[][] & Record<string, any> = [[]];
  if (typeof str === "string") {
    // Str is usually automatically parsed by JQuery
    str = JSON.parse(str);
  }
  var molecules = str.m;
  var atomsInFile = molecules[0].a; // Assumes there is at least one
  var bondsInFile = molecules[0].b; // molecule and ignores any more
  // Ignores any shapes
  var styles = molecules[0].s;
  var parseStyle =
    options !== undefined && options.parseStyle !== undefined
      ? options.parseStyle
      : styles !== undefined;

  var offset = atoms[atoms.length - 1].length; // When adding atoms their index will be
  // Offset by the number of existing atoms

  for (var i = 0; i < atomsInFile.length; i++) {
    var currentAtom = atomsInFile[i];
    var atom: Record<string, any> = {};
    atom.id = currentAtom.i; // Probably won't exist. Doesn't seem to
    // break anything.
    atom.x = currentAtom.x;
    atom.y = currentAtom.y;
    atom.z = currentAtom.z || 0; // Default value if file is 2D

    atom.bonds = [];
    atom.bondOrder = [];

    var elem = currentAtom.l || "C";
    atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();

    atom.serial = atoms[atoms.length - 1].length;
    if (parseStyle) {
      atom.style = styles[currentAtom.s || 0];
    }
    atoms[atoms.length - 1].push(atom);
  }
  for (let i = 0; i < bondsInFile.length; i++) {
    let currentBond = bondsInFile[i];
    let beginIndex = currentBond.b + offset;
    let endIndex = currentBond.e + offset;
    let bondOrder = currentBond.o || 1;

    let firstAtom = atoms[atoms.length - 1][beginIndex];
    let secondAtom = atoms[atoms.length - 1][endIndex];

    firstAtom.bonds.push(endIndex);
    firstAtom.bondOrder.push(bondOrder);
    secondAtom.bonds.push(beginIndex);
    secondAtom.bondOrder.push(bondOrder);
  }
  return atoms;
}
