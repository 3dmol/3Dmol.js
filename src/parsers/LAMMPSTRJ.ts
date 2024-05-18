import { AtomSpec } from "specs";
import { ParserOptionsSpec } from "./ParserOptionsSpec";
import { assignBonds } from "./utils/assignBonds";

/**
 * Parse a lammps trajectory file from str and create atoms
 * 
 * @category Parsers
*/
const dic: Record<string,string> = {
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

export function LAMMPSTRJ(str: string, options: ParserOptionsSpec) {
  const atoms: AtomSpec[][] = [];
  const lines = str.split(/\r?\n|\r/);
  let offset = 0;
  let atomCount = 0;
  let start = 0;
  while (start < lines.length - 9) {
    for (var j = start; j < lines.length; j++) {
      if (lines[j].match(/ITEM: NUMBER OF ATOMS/))
        atomCount = parseInt(lines[j + 1]);
      if (lines[j].match(/ITEM: ATOMS/)) {
        offset = j + 1;
        break;
      }
    }
    const types = lines[offset - 1].replace("ITEM: ATOMS ", "").split(" ");
    atoms.push([]);
    for (let j = offset; j < offset + atomCount; j++) {
      const atom: AtomSpec = {};
      const properties = {};
      const tokens = lines[j].split(" ");
      for (let k = 0; k < tokens.length; k++) {
        const prop = dic[types[k]];
        if (prop !== undefined) {
          if (prop === "serial") atom[prop] = parseInt(tokens[k]);
          else if (prop === "x" || prop === "y" || prop === "z")
            atom[prop] = parseFloat(tokens[k]);
          else if (prop === "charge" || prop === "radius")
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
    for (let i = 0; i < atoms.length; i++) assignBonds(atoms[i], options);
  }
  return atoms;
}
