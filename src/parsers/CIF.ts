import { ParserOptionsSpec } from './ParserOptionsSpec';
// puts atoms specified in mmCIF fromat in str into atoms

import { assignBonds } from "./utils/assignBonds";
import { computeSecondaryStructure } from "./utils/computeSecondaryStructure";
import { processSymmetries } from "./utils/processSymmetries";
import { conversionMatrix3, Matrix4, Vector3,  } from "../WebGL"

/**
 * @param {string} str
 * @param {ParserOptionsSpec} options
 * @category Parsers

 */
export function CIF(str: string, options: ParserOptionsSpec = {}) {
  var atoms: any[] & Record<string, any> = [];
  var noAssembly = !options.doAssembly; // don't assemble by default
  var modelData: any[] & Record<string, any> = (atoms.modelData = []);
  var assignbonds =
    options.assignBonds === undefined ? true : options.assignBonds;

  //coordinate conversion
  var fractionalToCartesian = function (cmat, x, y, z) {
    return new Vector3(x, y, z).applyMatrix3(cmat);
  };

  // Used to handle quotes correctly
  function splitRespectingQuotes(string, separator) {
    var sections: any[] = [];
    var sectionStart = 0;
    var sectionEnd = 0;
    while (sectionEnd < string.length) {
      while (
        string.substr(sectionEnd, separator.length) !== separator &&
        sectionEnd < string.length
      ) {
        // currently does not support escaping quotes
        if (string[sectionEnd] === "'") {
          sectionEnd++;
          while (sectionEnd < string.length && string[sectionEnd] !== "'") {
            sectionEnd++;
          }
        } else if (string[sectionEnd] === '"') {
          sectionEnd++;
          while (sectionEnd < string.length && string[sectionEnd] !== '"') {
            sectionEnd++;
          }
        }
        sectionEnd++;
      }
      sections.push(string.substr(sectionStart, sectionEnd - sectionStart));
      sectionStart = sectionEnd = sectionEnd + separator.length;
    }
    return sections;
  }

  var lines = str.split(/\r?\n|\r/);
  // Filter text to remove comments, trailing spaces, and empty lines
  var linesFiltered: string[] = [];
  var trimDisabled = false;
  for (let lineNum = 0; lineNum < lines.length; lineNum++) {
    // first remove comments
    // incorrect if #'s are allowed in strings
    // comments might only be allowed at beginning of line, not sure
    var line = lines[lineNum].split("#")[0];

    // inside data blocks, the string must be left verbatim
    // datablocks are started with a ';' at the beginning of a line
    // and ended with a ';' on its own line.
    if (trimDisabled) {
      if (line[0] === ";") {
        trimDisabled = false;
      }
    } else {
      if (line[0] === ";") {
        trimDisabled = true;
      }
    }

    if (trimDisabled || line !== "") {
      if (!trimDisabled) {
        line = line.trim();
        if (line[0] === "_") {
          // Replace dot separating category from data item with underscore. Dots aren't guarenteed, to makes
          // files consistent.
          var dot = line.split(/\s/)[0].indexOf(".");
          if (dot > -1) {
            let lineArr = line.split('')
            lineArr[dot] = "_";
            line = lineArr.join('');
            line = line.substr(0, dot) + "_" + line.substr(dot + 1);
          }
        }
      }
      linesFiltered.push(line);
    }
  }

  var lineNum = 0;
  while (lineNum < linesFiltered.length) {
    while (
      !linesFiltered[lineNum].startsWith("data_") ||
      linesFiltered[lineNum] === "data_global"
    ) {
      lineNum++;
    }
    lineNum++;

    // Process the lines and puts all of the data into an object.
    var mmCIF: Record<string, any> = {};
    while (
      lineNum < linesFiltered.length &&
      !linesFiltered[lineNum].startsWith("data_")
    ) {
      if (linesFiltered[lineNum][0] === undefined) {
        lineNum++;
      } else if (linesFiltered[lineNum][0] === "_") {
        var dataItemName = linesFiltered[lineNum].split(/\s/)[0].toLowerCase();
        var dataItem = (mmCIF[dataItemName] = mmCIF[dataItemName] || []);

        // if nothing left on the line go to the next one
        var restOfLine = linesFiltered[lineNum].substr(
          linesFiltered[lineNum].indexOf(dataItemName) + dataItemName.length
        );
        if (restOfLine === "") {
          lineNum++;
          if (linesFiltered[lineNum][0] === ";") {
            var dataBlock = linesFiltered[lineNum].substr(1);
            lineNum++;
            while (linesFiltered[lineNum] !== ";") {
              dataBlock = dataBlock + "\n" + linesFiltered[lineNum];
              lineNum++;
            }
            dataItem.push(dataBlock);
          } else {
            dataItem.push(linesFiltered[lineNum]);
          }
        } else {
          dataItem.push(restOfLine.trim());
        }
        lineNum++;
      } else if (linesFiltered[lineNum].substr(0, 5) === "loop_") {
        lineNum++;
        var dataItems: any[] = [];
        while (
          linesFiltered[lineNum] === "" ||
          linesFiltered[lineNum][0] === "_"
        ) {
          if (linesFiltered[lineNum] !== "") {
            let dataItemName = linesFiltered[lineNum]
              .split(/\s/)[0]
              .toLowerCase();
            let dataItem = (mmCIF[dataItemName] = mmCIF[dataItemName] || []);
            dataItems.push(dataItem);
          }
          lineNum++;
        }

        var currentDataItem = 0;
        while (
          lineNum < linesFiltered.length &&
          linesFiltered[lineNum][0] !== "_" &&
          !linesFiltered[lineNum].startsWith("loop_") &&
          !linesFiltered[lineNum].startsWith("data_")
        ) {
          let line = splitRespectingQuotes(linesFiltered[lineNum], " ");
          for (var field = 0; field < line.length; field++) {
            if (line[field] !== "") {
              dataItems[currentDataItem].push(line[field]);
              currentDataItem = (currentDataItem + 1) % dataItems.length;
            }
          }
          lineNum++;
        }
      } else {
        lineNum++;
      }
    }

    modelData.push({ symmetries: [] });

    // Pulls atom information out of the data
    atoms.push([]);
    var atomCount =
      mmCIF._atom_site_id !== undefined
        ? mmCIF._atom_site_id.length
        : mmCIF._atom_site_label.length;

    var conversionMatrix;
    if (mmCIF._cell_length_a !== undefined) {
      var a = parseFloat(mmCIF._cell_length_a);
      var b = parseFloat(mmCIF._cell_length_b);
      var c = parseFloat(mmCIF._cell_length_c);
      var alpha_deg = parseFloat(mmCIF._cell_angle_alpha) || 90;
      var beta_deg = parseFloat(mmCIF._cell_angle_beta) || 90;
      var gamma_deg = parseFloat(mmCIF._cell_angle_gamma) || 90;

      conversionMatrix = conversionMatrix3(
        a,
        b,
        c,
        alpha_deg,
        beta_deg,
        gamma_deg
      );
      modelData[modelData.length - 1].cryst = {
        a: a,
        b: b,
        c: c,
        alpha: alpha_deg,
        beta: beta_deg,
        gamma: gamma_deg,
      };
    }

    for (var i = 0; i < atomCount; i++) {
      if (
        mmCIF._atom_site_group_pdb !== undefined &&
        mmCIF._atom_site_group_pdb[i] === "TER"
      )
        continue;
      var atom: Record<string, any> = {};
      if (mmCIF._atom_site_cartn_x !== undefined) {
        atom.x = parseFloat(mmCIF._atom_site_cartn_x[i]);
        atom.y = parseFloat(mmCIF._atom_site_cartn_y[i]);
        atom.z = parseFloat(mmCIF._atom_site_cartn_z[i]);
      } else {
        var coords = fractionalToCartesian(
          conversionMatrix,
          parseFloat(mmCIF._atom_site_fract_x[i]),
          parseFloat(mmCIF._atom_site_fract_y[i]),
          parseFloat(mmCIF._atom_site_fract_z[i])
        );
        atom.x = coords.x;
        atom.y = coords.y;
        atom.z = coords.z;
      }
      atom.chain = mmCIF._atom_site_auth_asym_id
        ? mmCIF._atom_site_auth_asym_id[i]
        : undefined;
      atom.resi = mmCIF._atom_site_auth_seq_id
        ? parseInt(mmCIF._atom_site_auth_seq_id[i])
        : undefined;
      atom.resn = mmCIF._atom_site_auth_comp_id
        ? mmCIF._atom_site_auth_comp_id[i].trim()
        : undefined;
      atom.atom = mmCIF._atom_site_auth_atom_id
        ? mmCIF._atom_site_auth_atom_id[i].replace(/"/gm, "")
        : undefined; //"primed" names are in quotes
      atom.hetflag =
        !mmCIF._atom_site_group_pdb ||
        mmCIF._atom_site_group_pdb[i] === "HETA" ||
        mmCIF._atom_site_group_pdb[i] === "HETATM";
      var elem = "X";
      if (mmCIF._atom_site_type_symbol) {
        elem = mmCIF._atom_site_type_symbol[i].replace(/\(?\+?\d+.*/, "");
      } else if (mmCIF._atom_site_label) {
        //first two components are concatenated, then separated by underscore
        //best I can do is assume second component, if present, starts with a number
        elem = mmCIF._atom_site_label[i].split("_")[0].replace(/\(?\d+.*/, "");
      }
      atom.elem = elem[0].toUpperCase() + elem.substr(1, 1).toLowerCase();
      atom.bonds = [];
      atom.ss = "c";
      atom.serial = i;
      atom.bondOrder = [];
      atom.properties = {};
      atoms[atoms.length - 1].push(atom);
    }

    if (mmCIF._pdbx_struct_oper_list_id !== undefined && !noAssembly) {
      for (let i = 0; i < mmCIF._pdbx_struct_oper_list_id.length; i++) {
        var matrix11 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[1][1]"][i]
        );
        var matrix12 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[1][2]"][i]
        );
        var matrix13 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[1][3]"][i]
        );
        var vector1 = parseFloat(mmCIF["_pdbx_struct_oper_list_vector[1]"][i]);
        var matrix21 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[2][1]"][i]
        );
        var matrix22 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[2][2]"][i]
        );
        var matrix23 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[2][3]"][i]
        );
        var vector2 = parseFloat(mmCIF["_pdbx_struct_oper_list_vector[2]"][i]);
        var matrix31 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[3][1]"][i]
        );
        var matrix32 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[3][2]"][i]
        );
        var matrix33 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[3][3]"][i]
        );
        var vector3 = parseFloat(mmCIF["_pdbx_struct_oper_list_vector[3]"][i]);

        var matrix = new Matrix4(
          matrix11,
          matrix12,
          matrix13,
          vector1,
          matrix21,
          matrix22,
          matrix23,
          vector2,
          matrix31,
          matrix32,
          matrix33,
          vector3
        );
        modelData[modelData.length - 1].symmetries.push(matrix);
      }
    }
    var parseTerm = function (term) {
      var negative = term.match("-");
      term = term.replace(/[-xyz]/g, "");
      var fractionParts = term.split("/");

      var numerator, denominator;
      if (fractionParts[1] === undefined) {
        denominator = 1;
      } else {
        denominator = parseInt(fractionParts[1]);
      }
      if (fractionParts[0] === "") {
        numerator = 1;
      } else {
        numerator = parseInt(fractionParts[0]);
      }
      return (numerator / denominator) * (negative ? -1 : 1);
    };
    if (mmCIF._symmetry_equiv_pos_as_xyz !== undefined && !noAssembly) {
      for (var sym = 0; sym < mmCIF._symmetry_equiv_pos_as_xyz.length; sym++) {
        var transform = mmCIF._symmetry_equiv_pos_as_xyz[sym].replace(
          /["' ]/g,
          ""
        );
        var componentStrings = transform.split(",").map(function (val) {
          return val.replace(/-/g, "+-");
        });
        let matrix = new Matrix4(
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          1
        );
        for (let coord = 0; coord < 3; coord++) {
          var terms = componentStrings[coord].split("+");
          for (let t = 0; t < terms.length; t++) {
            var term = terms[t];
            if (term === "") continue;
            var coefficient = parseTerm(term);
            if (term.match("x")) {
              matrix.elements[coord + 0] = coefficient;
            } else if (term.match("y")) {
              matrix.elements[coord + 4] = coefficient;
            } else if (term.match("z")) {
              matrix.elements[coord + 8] = coefficient;
            } else {
              matrix.elements[coord + 12] = coefficient;
            }
          }
        }
        var conversionMatrix4 = conversionMatrix.getMatrix4();
        var conversionInverse = new Matrix4().getInverse(
          conversionMatrix4,
          true
        );
        matrix = new Matrix4().multiplyMatrices(
          matrix,
          conversionInverse
        );
        matrix = new Matrix4().multiplyMatrices(
          conversionMatrix4,
          matrix
        );
        modelData[modelData.length - 1].symmetries.push(matrix);
      }
    }
  }
  for (let i = 0; i < atoms.length; i++) {
    if (assignbonds) assignBonds(atoms[i]);
    computeSecondaryStructure(atoms[i],options.hbondCutoff);
    processSymmetries(
      modelData[i].symmetries,
      atoms[i],
      options,
      modelData[i].cryst
    );
    if (
      options.duplicateAssemblyAtoms &&
      !options.dontConnectDuplicatedAtoms &&
      assignbonds
    )
      assignBonds(atoms[i]);
  }

  return atoms;
}
