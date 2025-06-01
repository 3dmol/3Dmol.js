import { ParserOptionsSpec } from "./ParserOptionsSpec";

import { computeSecondaryStructure } from "./utils/computeSecondaryStructure";
import { processSymmetries } from "./utils/processSymmetries";
import { conversionMatrix3, Matrix3, Matrix4, Vector3 } from "../WebGL";
import { assignPDBBonds } from "./utils/assignPDBBonds";
import { AtomSpec } from "specs";

//coordinate conversion
const fractionalToCartesian = function (
  cmat: Matrix3,
  x: number,
  y: number,
  z: number
) {
  return new Vector3(x, y, z).applyMatrix3(cmat);
};

/**
 * Puts atoms specified in mmCIF fromat in str into atoms
 *
 * @param {string} str
 * @param {ParserOptionsSpec} options
 * @category Parsers
 */
export function CIF(str: string, options: ParserOptionsSpec = {}) {
  const atoms: Array<AtomSpec[]> & { modelData?: unknown } = [];
  const noAssembly = !options.doAssembly; // don't assemble by default
  const modelData = (atoms.modelData = []);
  const assignbonds =
    options.assignBonds === undefined ? true : options.assignBonds;

  // Used to handle quotes correctly
  function splitRespectingQuotes(string: string, separator: string) {
    const sections: string[] = [];
    let sectionStart = 0;
    let sectionEnd = 0;
    while (sectionEnd < string.length) {
      while (
        string.substring(sectionEnd, sectionEnd + separator.length) !==
          separator &&
        sectionEnd < string.length
      ) {
        // currently does not support escaping quotes
        if (string[sectionEnd] === "'") {
          sectionEnd++;
          while (sectionEnd < string.length && string[sectionEnd] !== "'") {
            sectionEnd++;
          }
          //biopython apparently generates invalid string literals so if we think we are done but aren't at a separator keep going
          while (
            string.substring(sectionEnd, sectionEnd + separator.length) !==
              separator &&
            sectionEnd < string.length
          ) {
            sectionEnd++;
          }
        } else if (string[sectionEnd] === '"') {
          sectionEnd++;
          while (sectionEnd < string.length && string[sectionEnd] !== '"') {
            sectionEnd++;
          }
          sectionEnd++;
        } else {
          sectionEnd++;
        }
      }
      sections.push(string.substring(sectionStart, sectionEnd));
      sectionStart = sectionEnd = sectionEnd + separator.length;
    }
    return sections;
  }

  const lines = str.split(/\r?\n|\r/);
  // Filter text to remove comments, trailing spaces, and empty lines
  const linesFiltered: string[] = [];
  let trimDisabled = false;
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
          const dot = line.split(/\s/)[0].indexOf(".");
          if (dot > -1) {
            let lineArr = line.split("");
            lineArr[dot] = "_";
            line = lineArr.join("");
            line = line.substring(0, dot) + "_" + line.substring(dot + 1);
          }
        }
      }
      linesFiltered.push(line);
    }
  }

  let lineNum = 0;
  while (lineNum < linesFiltered.length) {
    while (
      !linesFiltered[lineNum].startsWith("data_") ||
      linesFiltered[lineNum] === "data_global"
    ) {
      lineNum++;
    }
    lineNum++;

    // Process the lines and puts all of the data into an object.
    const mmCIF: Record<string, any> = {};
    while (
      lineNum < linesFiltered.length &&
      !linesFiltered[lineNum].startsWith("data_")
    ) {
      if (linesFiltered[lineNum][0] === undefined) {
        lineNum++;
      } else if (linesFiltered[lineNum][0] === "_") {
        const dataItemName = linesFiltered[lineNum]
          .split(/\s/)[0]
          .toLowerCase();
        const dataItem = (mmCIF[dataItemName] = mmCIF[dataItemName] || []);

        // if nothing left on the line go to the next one
        const restOfLine = linesFiltered[lineNum].substring(
          linesFiltered[lineNum].indexOf(dataItemName) + dataItemName.length
        );
        if (restOfLine === "") {
          lineNum++;
          if (linesFiltered[lineNum][0] === ";") {
            let dataBlock = linesFiltered[lineNum].substring(1);
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
      } else if (linesFiltered[lineNum].substring(0, 5) === "loop_") {
        lineNum++;
        const dataItems = [];
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

        let currentDataItem = 0;
        while (
          lineNum < linesFiltered.length &&
          linesFiltered[lineNum][0] !== "_" &&
          !linesFiltered[lineNum].startsWith("loop_") &&
          !linesFiltered[lineNum].startsWith("data_")
        ) {
          let line = splitRespectingQuotes(linesFiltered[lineNum], " ");
          for (let field = 0; field < line.length; field++) {
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
    const atomCount =
      mmCIF._atom_site_id !== undefined
        ? mmCIF._atom_site_id.length
        : mmCIF._atom_site_label.length;

    let conversionMatrix: Matrix3;
    if (mmCIF._cell_length_a !== undefined) {
      const a = parseFloat(mmCIF._cell_length_a);
      const b = parseFloat(mmCIF._cell_length_b);
      const c = parseFloat(mmCIF._cell_length_c);
      const alpha_deg = parseFloat(mmCIF._cell_angle_alpha) || 90;
      const beta_deg = parseFloat(mmCIF._cell_angle_beta) || 90;
      const gamma_deg = parseFloat(mmCIF._cell_angle_gamma) || 90;

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

    for (let i = 0; i < atomCount; i++) {
      if (
        mmCIF._atom_site_group_pdb !== undefined &&
        mmCIF._atom_site_group_pdb[i] === "TER"
      )
        continue;
      const atom: AtomSpec = {};
      if (mmCIF._atom_site_cartn_x !== undefined) {
        atom.x = parseFloat(mmCIF._atom_site_cartn_x[i]);
        atom.y = parseFloat(mmCIF._atom_site_cartn_y[i]);
        atom.z = parseFloat(mmCIF._atom_site_cartn_z[i]);
      } else {
        const coords = fractionalToCartesian(
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
        : mmCIF._atom_site_label_asym_id
        ? mmCIF._atom_site_label_asym_id[i]
        : undefined;
      atom.lchain = mmCIF._atom_site_label_asym_id
        ? mmCIF._atom_site_label_asym_id[i]
        : undefined;
      atom.resi = mmCIF._atom_site_auth_seq_id
        ? parseInt(mmCIF._atom_site_auth_seq_id[i])
        : mmCIF._atom_site_label_seq_id
        ? mmCIF._atom_site_label_seq_id[i]
        : undefined;
      atom.resn = mmCIF._atom_site_auth_comp_id
        ? mmCIF._atom_site_auth_comp_id[i].trim()
        : mmCIF._atom_site_label_comp_id
        ? mmCIF._atom_site_label_comp_id[i].trim()
        : undefined;
      atom.atom = mmCIF._atom_site_auth_atom_id
        ? mmCIF._atom_site_auth_atom_id[i].replace(/"/gm, "")
        : mmCIF._atom_site_label_atom_id
        ? mmCIF._atom_site_label_atom_id[i].replace(/"/gm, "")
        : undefined; //"primed" names are in quotes
      atom.hetflag =
        !mmCIF._atom_site_group_pdb ||
        mmCIF._atom_site_group_pdb[i] === "HETA" ||
        mmCIF._atom_site_group_pdb[i] === "HETATM";
      if(mmCIF._atom_site_b_iso_or_equiv ) {
        atom.b = parseFloat(mmCIF._atom_site_b_iso_or_equiv[i]);
      }
      let elem = "X";
      if (mmCIF._atom_site_type_symbol) {
        elem = mmCIF._atom_site_type_symbol[i].replace(/\(?\+?\d+.*/, "");
      } else if (mmCIF._atom_site_label) {
        //first two components are concatenated, then separated by underscore
        //best I can do is assume second component, if present, starts with a number
        elem = mmCIF._atom_site_label[i].split("_")[0].replace(/\(?\d+.*/, "");
      }
      atom.elem = elem[0].toUpperCase() + elem.substring(1, 2).toLowerCase();
      atom.bonds = [];
      atom.ss = "c";
      atom.serial = i;
      atom.bondOrder = [];
      atom.properties = {};
      atoms[atoms.length - 1].push(atom);
    }

    if (mmCIF._pdbx_struct_oper_list_id !== undefined && !noAssembly) {
      for (let i = 0; i < mmCIF._pdbx_struct_oper_list_id.length; i++) {
        const matrix11 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[1][1]"][i]
        );
        const matrix12 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[1][2]"][i]
        );
        const matrix13 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[1][3]"][i]
        );
        const vector1 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_vector[1]"][i]
        );
        const matrix21 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[2][1]"][i]
        );
        const matrix22 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[2][2]"][i]
        );
        const matrix23 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[2][3]"][i]
        );
        const vector2 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_vector[2]"][i]
        );
        const matrix31 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[3][1]"][i]
        );
        const matrix32 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[3][2]"][i]
        );
        const matrix33 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_matrix[3][3]"][i]
        );
        const vector3 = parseFloat(
          mmCIF["_pdbx_struct_oper_list_vector[3]"][i]
        );

        const matrix = new Matrix4(
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
    const parseTerm = function (term: string) {
      const negative = term.match("-");
      term = term.replace(/[-xyz]/g, "");
      const fractionParts = term.split("/");

      let numerator: number, denominator: number;
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
      for (let sym = 0; sym < mmCIF._symmetry_equiv_pos_as_xyz.length; sym++) {
        const transform = mmCIF._symmetry_equiv_pos_as_xyz[sym].replace(
          /["' ]/g,
          ""
        );
        const componentStrings = transform
          .split(",")
          .map(function (val: string) {
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
          const terms = componentStrings[coord].split("+");
          for (let t = 0; t < terms.length; t++) {
            const term = terms[t];
            if (term === "") continue;
            const coefficient = parseTerm(term);
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
        const conversionMatrix4 = conversionMatrix.getMatrix4();
        const conversionInverse = new Matrix4().getInverse(
          conversionMatrix4,
          true
        );
        matrix = new Matrix4().multiplyMatrices(matrix, conversionInverse);
        matrix = new Matrix4().multiplyMatrices(conversionMatrix4, matrix);
        modelData[modelData.length - 1].symmetries.push(matrix);
      }
    }
  }
  for (let i = 0; i < atoms.length; i++) {
    if (
      assignbonds &&
      !(options.duplicateAssemblyAtoms && !options.dontConnectDuplicatedAtoms)
    ) {
      assignPDBBonds(atoms[i], options);
    }
    computeSecondaryStructure(atoms[i], options.hbondCutoff);
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
      assignPDBBonds(atoms[i], options);
  }

  return atoms;
}
