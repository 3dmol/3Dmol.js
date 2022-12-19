// @ts-nocheck
/**
 * Parse a prmtop file from str and create atoms
 * 
  * @param {string}
  *            str
  * @param {ParserOptionsSpec}
  *            options - noSecondaryStructure (do not compute ss)
  *  
 * @category Parsers

 */
export function PRMTOP(str /*, options*/) {
  var atoms = [];
  var atomIndex;
  var count = 0;
  var lines = str.split(/\r?\n|\r/);
  if (lines.length > 0 && lines[0].includes("VERSION")) {
    var sectionList = lines.filter(function (line) {
      //store the relevant section lists
      return (
        line.includes("POINTERS") ||
        line.includes("ATOM_NAME") ||
        line.includes("CHARGE") ||
        line.includes("RADII") ||
        line.includes("BONDS_INC_HYDROGEN") ||
        line.includes("BONDS_WITHOUT_HYDROGEN")
      );
    });
    var index = getIndex("POINTERS");
    if (index == -1) return [];
    var col = getColEleSize(index);
    var atomCount = parseInt(lines[index + 1].slice(0, col[1]));
    if (isNaN(atomCount) || atomCount <= 0) return [];
    index = getIndex("ATOM_NAME");
    if (index == -1) return [];
    col = getColEleSize(index);
    var noOfCol = col[0];
    for (let i = 0; i < atomCount / col[0]; i++) {
      if (i == parseInt(atomCount / col[0])) noOfCol = atomCount % col[0];
      for (let j = 0; j < noOfCol; j++) {
        let atom = {};
        let properties = { charge: "", radii: "" };
        atom.serial = count;
        atom.x = 0;
        atom.y = 0;
        atom.z = 0;
        atom.atom = lines[index + 1].slice(col[1] * j, col[1] * (j + 1));
        atom.elem = lines[index + 1].slice(col[1] * j, col[1] * j + 1);
        atom.properties = properties;
        atom.bonds = [];
        atom.bondOrder = [];
        atoms.push(atom);
        count++;
      }
      index++;
    }
    index = getIndex("CHARGE");
    if (index != -1) {
      col = getColEleSize(index);
      count = 0;
      noOfCol = col[0];
      for (let i = 0; i < atomCount / col[0]; i++) {
        if (i == parseInt(atomCount / col[0])) noOfCol = atomCount % col[0];
        for (let j = 0; j < noOfCol; j++) {
          atoms[count].properties.charge = parseFloat(
            lines[index + 1].slice(col[1] * j, col[1] * (j + 1))
          );
          count++;
        }
        index++;
      }
    }
    index = getIndex("RADII");
    if (index != -1) {
      col = getColEleSize(index);
      count = 0;
      noOfCol = col[0];
      for (let i = 0; i < atomCount / col[0]; i++) {
        if (i == parseInt(atomCount / col[0])) noOfCol = atomCount % col[0];
        for (let j = 0; j < noOfCol; j++) {
          atoms[count].properties.radii = parseFloat(
            lines[index + 1].slice(col[1] * j, col[1] * (j + 1))
          );
          count++;
        }
        index++;
      }
    }
    index = getIndex("BONDS_WITHOUT_HYDROGEN");
    if (index != -1) {
      col = getColEleSize(index);
      count = 0;
      noOfCol = col[0];
      index = index + 1;
      while (!lines[index].match(/^%FLAG/)) {
        if (lines[index + 1].match(/^%FLAG/))
          //its the last line
          noOfCol = atomCount % col[0];
        for (let j = 0; j < noOfCol; j++) {
          if (count % 3 == 0) {
            atomIndex = parseInt(
              lines[index].slice(col[1] * j, col[1] * (j + 1)) / 3
            );
          } else if (count % 3 == 1) {
            atoms[atomIndex].bonds.push(
              parseInt(lines[index].slice(col[1] * j, col[1] * (j + 1)) / 3)
            );
          }
          count++;
        }
        index++;
      }
    }
    index = getIndex("BONDS_INC_HYDROGEN");
    if (index != -1) {
      col = getColEleSize(index);
      count = 0;
      noOfCol = col[0];
      index = index + 1;
      while (!lines[index].match(/^%FLAG/)) {
        if (lines[index + 1].match(/^%FLAG/))
          //its the last line
          noOfCol = atomCount % col[0];
        for (let j = 0; j < noOfCol; j++) {
          if (count % 3 == 0) {
            atomIndex = parseInt(
              lines[index].slice(col[1] * j, col[1] * (j + 1)) / 3
            );
          } else if (count % 3 == 1) {
            atoms[atomIndex].bonds.push(
              parseInt(lines[index].slice(col[1] * j, col[1] * (j + 1)) / 3)
            );
          }
          count++;
        }
        index++;
      }
    }
  } else {
    return [];
  }

  function getIndex(section) {
    var index = lines.indexOf(
      sectionList.filter(function (line) {
        return line.includes(section);
      })[0]
    ); //returns the index of the line containing FLAG POINTERS
    if (Number.isInteger(index) && index > 0) {
      while (!lines[index].includes("FORMAT"))
        //doing this so as to take comments into consideration
        index++;
      return index;
    } else {
      return -1;
    }
  }
  function getColEleSize(i) {
    var numberOfCol = lines[i].match(/\((\d*)\S*/); // stores the number of columns
    var elementSize = lines[i].match(/[a-zA-Z](\d*)\)\s*/);
    if (elementSize == null) {
      elementSize = lines[i].match(/[a-zA-Z](\d*)\.\d*\)\s*/); //stores the element size
    }
    return [numberOfCol[1], elementSize[1]];
  }
  return [atoms];
}
