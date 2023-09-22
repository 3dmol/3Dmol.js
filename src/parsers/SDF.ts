import { ParserOptionsSpec } from "./ParserOptionsSpec";

var parseV2000 = function (lines: any, options: ParserOptionsSpec) {
  var atoms: any & Record<string, any> = [[]];
  var noH = false;
  if (typeof options.keepH !== "undefined") noH = !options.keepH;

  while (lines.length > 0) {
    if (lines.length < 4) break;
    var atomCount = parseInt(lines[3].substring(0, 3));
    if (isNaN(atomCount) || atomCount <= 0) break;
    var bondCount = parseInt(lines[3].substring(3, 6));
    var offset = 4;
    if (lines.length < 4 + atomCount + bondCount) break;

    // Serial is atom's index in file; index is atoms index in 'atoms'
    var serialToIndex: number[] = [];
    var start = atoms[atoms.length - 1].length;
    var end = start + atomCount;
    var i: number, line: string;
    for (i = start; i < end; i++, offset++) {
      line = lines[offset];
      var atom: Record<string, any> = {};
      var elem = line.substring(31, 34).replace(/ /g, "");
      atom.atom = atom.elem =
        elem[0].toUpperCase() + elem.substring(1).toLowerCase();

      if (atom.elem !== "H" || !noH) {
        atom.serial = i;
        serialToIndex[i] = atoms[atoms.length - 1].length;
        atom.x = parseFloat(line.substring(0, 10));
        atom.y = parseFloat(line.substring(10, 20));
        atom.z = parseFloat(line.substring(20, 30));
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
      var from = serialToIndex[parseInt(line.substring(0, 3)) - 1 + start];
      var to = serialToIndex[parseInt(line.substring(3, 6)) - 1 + start];
      var order = parseFloat(line.substring(6));
      if (typeof from != "undefined" && typeof to != "undefined") {
        atoms[atoms.length - 1][from].bonds.push(to);
        atoms[atoms.length - 1][from].bondOrder.push(order);
        atoms[atoms.length - 1][to].bonds.push(from);
        atoms[atoms.length - 1][to].bondOrder.push(order);
      }
    }
    if (options.multimodel) {
      if (!options.onemol) atoms.push([]);
      while (lines[offset] !== "$$$$" && offset < lines.length) offset++;
      lines.splice(0, ++offset);
    } else {
      break;
    }
  }
  return atoms;
};

/**
 * @param {!Array.<string>} lines
 * @param {ParserOptionsSpec} options
 * @returns {!Array.<!Array<!Object>>}
*/

var parseV3000 = function (lines: any, options: ParserOptionsSpec) {
  var atoms: any[][] & Record<string, any> = [[]];
  var noH = false;
  if (typeof options.keepH !== "undefined") noH = !options.keepH;

  while (lines.length > 0) {
    if (lines.length < 8) break;

    if (!lines[4].startsWith("M  V30 BEGIN CTAB")) break;
    if (!lines[5].startsWith("M  V30 COUNTS") || lines[5].length < 14) break;

    var counts = lines[5].substring(13).match(/\S+/g);

    if (counts.length < 2) break;

    var atomCount = parseInt(counts[0]);
    if (isNaN(atomCount) || atomCount <= 0) break;
    var bondCount = parseInt(counts[1]);
    var offset = 7;

    if (lines.length < 8 + atomCount + bondCount)
      //header, bgn+end CTAB, counts, END
      break;

    // serial is atom's index in file; index is atoms index in 'atoms'
    var serialToIndex: number[] = [];
    var start = atoms[atoms.length - 1].length;
    var end = start + atomCount;
    var i: number, line: string;
    for (i = start; i < end; i++, offset++) {
      line = lines[offset];
      var atomParts = line.substring(6).match(/\S+/g);
      if (atomParts!.length > 4) {
        var atom: Record<string, any> = {};
        var elem = atomParts![1].replace(/ /g, "");
        atom.atom = atom.elem =
          elem[0].toUpperCase() + elem.substring(1).toLowerCase();

        if (atom.elem !== "H" || !noH) {
          atom.serial = i;
          serialToIndex[i] = atoms[atoms.length - 1].length;
          atom.x = parseFloat(atomParts![2]);
          atom.y = parseFloat(atomParts![3]);
          atom.z = parseFloat(atomParts![4]);
          atom.hetflag = true;
          atom.bonds = [];
          atom.bondOrder = [];
          atom.properties = {};
          atom.index = atoms[atoms.length - 1].length;
          atoms[atoms.length - 1].push(atom);
        }
      }
    }

    if (lines[offset] === "M  V30 END ATOM") offset++;
    else break;

    if (bondCount !== 0 && lines[offset] === "M  V30 BEGIN BOND") offset++;
    else break;

    for (i = 0; i < bondCount; i++, offset++) {
      line = lines[offset];
      var bondParts = line.substring(6).match(/\S+/g);
      if (bondParts!.length > 3) {
        var from = serialToIndex[parseInt(bondParts![2]) - 1 + start];
        var to = serialToIndex[parseInt(bondParts![3]) - 1 + start];
        var order = parseFloat(bondParts![1]);
        if (typeof from != "undefined" && typeof to != "undefined") {
          atoms[atoms.length - 1][from].bonds.push(to);
          atoms[atoms.length - 1][from].bondOrder.push(order);
          atoms[atoms.length - 1][to].bonds.push(from);
          atoms[atoms.length - 1][to].bondOrder.push(order);
        }
      }
    }
    if (options.multimodel) {
      if (!options.onemol) {
        atoms.push([]);
      }
      while (lines[offset] !== "$$$$" && offset < lines.length) {
        offset++;
      }
      lines.splice(0, ++offset);
    } else {
      break;
    }
  }
  return atoms;
};

/**
 * @param {string}
 *            str
 * @param {ParserOptionsSpec}
 *            options
 * @category Parsers
*/

export function SDF(str: string, options: ParserOptionsSpec) {
  var molformat = "V2000";
  var lines = str.split(/\r?\n|\r/);
  if (lines.length > 3 && lines[3].length > 38) {
    molformat = lines[3].substring(34, 39);
  }
  if (molformat === "V2000") {
    return parseV2000(lines, options);
  } else if (molformat === "V3000") {
    return parseV3000(lines, options);
  }
  return [['']];
}
