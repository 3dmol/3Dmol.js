import { assignBonds } from "./utils/assignBonds";

/**
 * Parse a lammps trajectory file from str and create atoms
 * @category Parsers

 */
export function LAMMPSTRJ(str, options) {
  var atoms: any[] = [];
  var dic = {
    id: "serial",
    type: "atom",
    element: "elem",
    q: "charge",
    radius: "radius",
    x: "x",
    xu: "x",
    xs: "x",
    xsu: "x",
    y: "y",
    yu: "y",
    ys: "y",
    ysu: "y",
    z: "z",
    zu: "z",
    zs: "z",
    zsu: "z",
  };
  var lines = str.split(/\r?\n|\r/);
  var offset = 0;
  var atomCount = 0;
  var start = 0;
  while (start < lines.length - 9) {
    for (var j = start; j < lines.length; j++) {
      if (lines[j].match(/ITEM: NUMBER OF ATOMS/))
        atomCount = parseInt(lines[j + 1]);
      if (lines[j].match(/ITEM: ATOMS/)) {
        offset = j + 1;
        break;
      }
    }
    var types = lines[offset - 1].replace("ITEM: ATOMS ", "").split(" ");
    atoms.push([]);
    for (let j = offset; j < offset + atomCount; j++) {
      var atom: Record<string, any> = {};
      var properties = {};
      var tokens = lines[j].split(" ");
      for (var k = 0; k < tokens.length; k++) {
        var prop = dic[types[k]];
        if (prop != undefined) {
          if (prop == "serial") atom[prop] = parseInt(tokens[k]);
          else if (prop == "x" || prop == "y" || prop === "z")
            atom[prop] = parseFloat(tokens[k]);
          else if (prop == "charge" || prop == "radius")
            properties[prop] = parseFloat(tokens[k]);
          else atom[prop] = tokens[k];
        }
        atom.properties = properties;
        atom.bonds = [];
        atom.bondOrder = [];
      }
      atoms[atoms.length - 1][j - offset] = atom;
    }
    start = offset + atomCount - 1;
  }
  if (options.assignBonds) {
    for (var i = 0; i < atoms.length; i++) assignBonds(atoms[i]);
  }
  return atoms;
}
