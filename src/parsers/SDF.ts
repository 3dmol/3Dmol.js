import { ParserOptionsSpec } from "./ParserOptionsSpec";

/**
 * Apply a 3DMOL_STYLE JSON block to atom bond styles and model data.
 *
 * JSON format:
 *   {"stick":{"radius":0.15},"bonds":{"1-5":{"color1":"red","color2":"blue"}}}
 *
 * - "bonds" maps "serial1-serial2" (1-based) to BondStyle properties,
 *   set on atom.bondStyles[] directly. Supported properties:
 *     - color1, color2: overall bond half-colors
 *     - dashedBondConfig: { solidColor, dashedColor } for fractional bond color overrides
 *     - dashedBondFlip: override which side gets the dashed line (boolean)
 *     - singleBond, radius, iswire
 *   When serials are given in reverse order, per-half color properties are swapped.
 * - All other keys are AtomStyleSpec, stored on atoms.modelData for the frame.
 */
const apply3DmolStyle = function (
  jsonStr: string,
  atoms: any,
  serialToIndex: number[],
  start: number
) {
  let styleObj;
  try {
    styleObj = JSON.parse(jsonStr);
  } catch (e) {
    console.warn("3DMOL_STYLE: invalid JSON: " + jsonStr);
    return;
  }

  const frameIdx = atoms.length - 1;
  const frameAtoms = atoms[frameIdx];

  if (styleObj.bonds && typeof styleObj.bonds === "object") {
    const bonds = styleObj.bonds;
    for (const bondKey in bonds) {
      if (!bonds.hasOwnProperty(bondKey)) continue;
      const parts = bondKey.split("-");
      if (parts.length !== 2) continue;

      const serial1 = parseInt(parts[0]) - 1 + start;
      const serial2 = parseInt(parts[1]) - 1 + start;
      const idx1 = serialToIndex[serial1];
      const idx2 = serialToIndex[serial2];
      if (idx1 == null || idx2 == null) continue;

      const bondProps = bonds[bondKey];

      // Rendering only checks bondStyles on the lower-indexed atom
      const lowerIdx = idx1 < idx2 ? idx1 : idx2;
      const higherIdx = idx1 < idx2 ? idx2 : idx1;
      const lowerAtom = frameAtoms[lowerIdx];

      const bondPos = lowerAtom.bonds.indexOf(higherIdx);
      if (bondPos < 0) continue;

      if (!lowerAtom.bondStyles) lowerAtom.bondStyles = [];
      if (!lowerAtom.bondStyles[bondPos]) lowerAtom.bondStyles[bondPos] = {};

      // If serials were given in reverse order, swap per-half color properties
      if (idx1 > idx2) {
        if (bondProps.color1 || bondProps.color2) {
          [bondProps.color1, bondProps.color2] = [bondProps.color2, bondProps.color1];
        }
      }

      Object.assign(lowerAtom.bondStyles[bondPos], bondProps);
    }
    delete styleObj.bonds;
  }

  if (Object.keys(styleObj).length > 0) {
    atoms.modelData[frameIdx].style = styleObj;
  }
};

/**
 * Parse SDF data items (property blocks) that follow M  END.
 *
 * Per the CTfile spec, each data item consists of:
 *   - A data header line starting with ">" containing a field name in < >
 *   - One or more data lines (up to 200 chars each per spec)
 *   - A blank line terminator
 *
 * The record ends at the $$$$ delimiter.
 *
 * @returns offset past $$$$ (or end of lines if missing)
 */
const parseSDProperties = function (
  lines: string[],
  offset: number,
  atoms: any,
  serialToIndex: number[],
  start: number
): number {
  // Advance past M  END (V3000 callers may still be inside structure blocks)
  while (offset < lines.length && !lines[offset].startsWith("M  END")) offset++;
  if (offset < lines.length) offset++; // past M  END

  while (offset < lines.length) {
    const line = lines[offset];
    let trimmed = line.trim();
    offset++;

    if (trimmed === "$$$$") return offset;

    // Data header: starts with >, field name in angle brackets
    const match = line.match(/^>\s*<([^>]+)>/);
    if (!match) continue;

    const propName = match[1].trim();

    // Collect data lines until blank line, next header, or $$$$
    const data: string[] = [];
    while (offset < lines.length) {
      trimmed = lines[offset].trim();
      if (trimmed === "" || trimmed === "$$$$" || lines[offset].charAt(0) === '>') break;
      data.push(lines[offset]);
      offset++;
    }

    if (propName === "3DMOL_STYLE" && data.length > 0) {
      apply3DmolStyle(data.join(" ").trim(), atoms, serialToIndex, start);
    }
  }

  return offset;
};

/** Add a bond between two atoms (bidirectional). */
const addBond = function (curFrame: any[], from: number, to: number, order: number) {
  curFrame[from].bonds.push(to);
  curFrame[from].bondOrder.push(order);
  curFrame[to].bonds.push(from);
  curFrame[to].bondOrder.push(order);
};

const parseV2000 = function (lines: any, options: ParserOptionsSpec) {
  const atoms: any & Record<string, any> = [[]];
  const modelData = atoms.modelData = [] as any[];
  modelData.push({});  // entry for frame 0
  const noH = typeof options.keepH !== "undefined" ? !options.keepH : false;
  let lineIndex = 0;

  while (lineIndex < lines.length) {
    if (lines.length - lineIndex < 4) break;
    const atomCount = parseInt(lines[lineIndex + 3].substring(0, 3));
    if (isNaN(atomCount) || atomCount <= 0) break;
    const bondCount = parseInt(lines[lineIndex + 3].substring(3, 6));
    let offset = 4;
    if (lines.length - lineIndex < 4 + atomCount + bondCount) break;

    // Serial is atom's index in file; index is atoms index in 'atoms'
    const serialToIndex: number[] = [];
    const curFrame = atoms[atoms.length - 1];
    const start = curFrame.length;
    const end = start + atomCount;
    for (let i = start; i < end; i++, offset++) {
      const line = lines[lineIndex + offset];
      const atom: Record<string, any> = {};
      const elem = line.substring(31, 34).replace(/ /g, "");
      atom.atom = atom.elem =
        elem[0].toUpperCase() + elem.substring(1).toLowerCase();

      if (atom.elem !== "H" || !noH) {
        atom.serial = i;
        serialToIndex[i] = curFrame.length;
        atom.x = parseFloat(line.substring(0, 10));
        atom.y = parseFloat(line.substring(10, 20));
        atom.z = parseFloat(line.substring(20, 30));
        atom.hetflag = true;
        atom.bonds = [];
        atom.bondOrder = [];
        atom.properties = {};
        atom.index = curFrame.length;
        curFrame.push(atom);
      }
    }

    for (let i = 0; i < bondCount; i++, offset++) {
      const line = lines[lineIndex + offset];
      const from = serialToIndex[parseInt(line.substring(0, 3)) - 1 + start];
      const to = serialToIndex[parseInt(line.substring(3, 6)) - 1 + start];
      const order = parseFloat(line.substring(6));
      if (typeof from != "undefined" && typeof to != "undefined") {
        addBond(curFrame, from, to, order);
      }
    }

    // Parse M  END and property blocks; returns absolute line index past $$$$
    lineIndex = parseSDProperties(lines, lineIndex + offset, atoms, serialToIndex, start);

    if (!options.multimodel) break;
    if (!options.onemol) {
      atoms.push([]);
      modelData.push({});  // 1:1 with frames — prevents undefined on frame switch
    }
  }
  return atoms;
};

/**
 * @param {!Array.<string>} lines
 * @param {ParserOptionsSpec} options
 * @returns {!Array.<!Array<!Object>>}
*/

const parseV3000 = function (lines: any, options: ParserOptionsSpec) {
  const atoms: any[][] & Record<string, any> = [[]];
  const modelData = atoms.modelData = [] as any[];
  modelData.push({});  // entry for frame 0
  const noH = typeof options.keepH !== "undefined" ? !options.keepH : false;
  let lineIndex = 0;

  while (lineIndex < lines.length) {
    if (lines.length - lineIndex < 8) break;

    if (!lines[lineIndex + 4].startsWith("M  V30 BEGIN CTAB")) break;
    if (!lines[lineIndex + 5].startsWith("M  V30 COUNTS") || lines[lineIndex + 5].length < 14) break;

    const counts = lines[lineIndex + 5].substring(13).match(/\S+/g);

    if (counts.length < 2) break;

    const atomCount = parseInt(counts[0]);
    if (isNaN(atomCount) || atomCount <= 0) break;
    const bondCount = parseInt(counts[1]);
    let offset = 7;

    if (lines.length - lineIndex < 8 + atomCount + bondCount)
      //header, bgn+end CTAB, counts, END
      break;

    // serial is atom's index in file; index is atoms index in 'atoms'
    const serialToIndex: number[] = [];
    const curFrame = atoms[atoms.length - 1];
    const start = curFrame.length;
    const end = start + atomCount;
    for (let i = start; i < end; i++, offset++) {
      const line = lines[lineIndex + offset];
      const atomParts = line.substring(6).match(/\S+/g);
      if (atomParts!.length > 4) {
        const atom: Record<string, any> = {};
        const elem = atomParts![1].replace(/ /g, "");
        atom.atom = atom.elem =
          elem[0].toUpperCase() + elem.substring(1).toLowerCase();

        if (atom.elem !== "H" || !noH) {
          atom.serial = i;
          serialToIndex[i] = curFrame.length;
          atom.x = parseFloat(atomParts![2]);
          atom.y = parseFloat(atomParts![3]);
          atom.z = parseFloat(atomParts![4]);
          atom.hetflag = true;
          atom.bonds = [];
          atom.bondOrder = [];
          atom.properties = {};
          atom.index = curFrame.length;
          curFrame.push(atom);
        }
      }
    }

    if (lines[lineIndex + offset] === "M  V30 END ATOM") offset++;
    else break;

    if (bondCount !== 0 && lines[lineIndex + offset] === "M  V30 BEGIN BOND") offset++;
    else break;

    for (let i = 0; i < bondCount; i++, offset++) {
      const line = lines[lineIndex + offset];
      const bondParts = line.substring(6).match(/\S+/g);
      if (bondParts!.length > 3) {
        const from = serialToIndex[parseInt(bondParts![2]) - 1 + start];
        const to = serialToIndex[parseInt(bondParts![3]) - 1 + start];
        const order = parseFloat(bondParts![1]);
        if (typeof from != "undefined" && typeof to != "undefined") {
          addBond(curFrame, from, to, order);
        }
      }
    }

    // Parse remaining V3000 structure blocks and property blocks
    lineIndex = parseSDProperties(lines, lineIndex + offset, atoms, serialToIndex, start);

    if (!options.multimodel) break;
    if (!options.onemol) {
      atoms.push([]);
      modelData.push({});  // 1:1 with frames
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
  let molformat = "V2000";
  const lines = str.split(/\r?\n|\r/);
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
