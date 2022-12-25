// parse pdb file from str and create atoms
// if computeStruct is true will always perform secondary structure
// analysis,
// otherwise only do analysis of SHEET/HELIX comments are missing

import { getSinglePDB } from "./utils/getSinglePDB";

/**
 * @param {string} str
 * @param {ParserOptionsSpec} options - keepH (do not strip hydrogens), noSecondaryStructure,
 *            assignbonds (default true, calculate implicit bonds)
 *            (do not compute ss), altLoc (which alternate location to select, if present; '*' to load all)
 * @category Parsers
 * 
 */
export function PDB(str, options) {
  options = options || {};
  var atoms: any[] & Record<string, any> = []; //a separate list for each model
  var sslookup = {}; //stores SHEET and HELIX info, which is shared across models
  atoms.modelData = [];
  var lines = str.split(/\r?\n|\r/);
  while (lines.length > 0) {
    var pdbinfo = getSinglePDB(lines, options, sslookup);
    var modelatoms = pdbinfo[0];
    var modelData = pdbinfo[1];
    lines = pdbinfo[2];

    if (modelatoms.length == 0) {
      continue; //happens when there are blank lines
    }
    if (options.multimodel && options.onemol && atoms.length > 0) {
      //merge into existing atoms
      var inc = atoms[0].length;
      for (var i = 0; i < modelatoms.length; i++) {
        //renumber
        var atom = modelatoms[i];
        atom.index = i;
        for (var b = 0; b < atom.bonds.length; b++) {
          atom.bonds[b] += inc;
        }
        atoms[0].push(atom);
      }
    } else {
      atoms.modelData.push(modelData);
      atoms.push(modelatoms);
    }

    if (!options.multimodel) {
      break;
    }
  }

  return atoms;
}
