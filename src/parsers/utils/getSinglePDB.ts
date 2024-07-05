import { Matrix4 } from "../../WebGL";
import { atomNameToElem } from "./atomNameToElem";
import { bondTable } from "./bondLength";
import { computeSecondaryStructure } from "./computeSecondaryStructure";
import { isEmpty } from "./isEmpty";
import { processSymmetries } from "./processSymmetries";
import { assignPDBBonds } from "./assignPDBBonds";
import { validateBonds } from "./validateBonds";
import { ParserOptionsSpec } from "../ParserOptionsSpec";
import { AtomSpec, Cryst } from "specs";
// Return one model worth of pdb, returns atoms, modelData, and remaining lines
export function getSinglePDB(
  lines: string[],
  options: ParserOptionsSpec,
  sslookup: { [x: string]: { [x: string]: string }; hasOwnProperty?: any }
): [
    AtomSpec[],
    {
      symmetries: Matrix4[];
      cryst: Omit<Cryst, "origin" | "size" | "unit" | "matrix4" | "matrix">;
    },
    string[]
  ] {
  const atoms: AtomSpec[] = [];
  const assignbonds =
    options.assignBonds === undefined ? true : options.assignBonds;
  const noH = !options.keepH; // suppress hydrogens by default
  const ignoreStruct = !!options.noSecondaryStructure;
  const computeStruct = !options.noComputeSecondaryStructure;
  const noAssembly = !options.doAssembly; // don't assemble by default
  const selAltLoc = options.altLoc ? options.altLoc : "A"; //default alternate location to select if present
  const modelData: {
    symmetries: Matrix4[];
    cryst: Omit<Cryst, "origin" | "size" | "unit" | "matrix4" | "matrix">;
  } = { symmetries: [], cryst: undefined };
  //atom name
  let atom: string;
  let remainingLines = [];

  const serialToIndex: number[] = []; // map from pdb serial to index in atoms
  let line: string | string[];
  const seenbonds: Record<string, number> = {}; //sometimes connect records are duplicated as an unofficial means of relaying bond orders

  for (let i = 0; i < lines.length; i++) {
    line = lines[i].replace(/^\s*/, ""); // remove indent
    const recordName = line.substring(0, 6);

    let startChain: string, startResi: number, endResi: number;

    if (recordName.indexOf("END") === 0) {
      remainingLines = lines.slice(i + 1);
      if (recordName === "END") {
        //as opposed to ENDMDL
        //reset secondary structure
        for (const prop in sslookup) {
          if (sslookup.hasOwnProperty(prop)) {
            delete sslookup[prop];
          }
        }
      }
      break;
    } else if (recordName === "ATOM  " || recordName === "HETATM") {
      let resn: string,
        chain: any,
        resi: string | number,
        icode: string,
        x: number,
        y: number,
        z: number,
        hetflag: boolean,
        elem: string | string[],
        serial: number,
        altLoc: string,
        b: number;
      altLoc = line.substring(16, 17);
      if (altLoc !== " " && altLoc !== selAltLoc && selAltLoc !== "*") continue;
      serial = parseInt(line.substring(6, 11));
      atom = line.substring(12, 16).replace(/ /g, "");
      resn = line.substring(17, 20).replace(/ /g, "");
      chain = line.substring(21, 22);
      resi = parseInt(line.substring(22, 26));
      icode = line.substring(26, 27);
      x = parseFloat(line.substring(30, 38));
      y = parseFloat(line.substring(38, 46));
      z = parseFloat(line.substring(46, 54));
      b = parseFloat(line.substring(60, 68));
      elem = line.substring(76, 78).replace(/ /g, "");
      if (elem === "" || bondTable[elem] === undefined) {
        // for some incorrect PDB files
        elem = atomNameToElem(line.substring(12, 14), line[0] == "A");
      } else {
        elem = elem[0].toUpperCase() + elem.substring(1).toLowerCase();
      }

      if (elem === "H" && noH) continue;
      if (recordName[0] == "H") hetflag = true;
      else hetflag = false;
      serialToIndex[serial] = atoms.length;
      atoms.push({
        resn,
        x,
        y,
        z,
        elem,
        hetflag,
        altLoc,
        chain,
        resi,
        icode: icode,
        rescode: resi + (icode !== " " ? "^" + icode : ""), // combo
        // resi
        // and
        // icode
        serial,
        atom,
        bonds: [],
        ss: "c",
        bondOrder: [],
        properties: {},
        b,
        pdbline: line,
      });
    } else if (recordName === "SHEET ") {
      startChain = line.substring(21, 22);
      startResi = parseInt(line.substring(22, 26));
      endResi = parseInt(line.substring(33, 37));
      if (!(startChain in sslookup)) {
        sslookup[startChain] = {};
      }
      //mark start and end with additional character
      sslookup[startChain][startResi] = "s1";
      for (let res = startResi + 1; res < endResi; res++) {
        sslookup[startChain][res] = "s";
      }
      sslookup[startChain][endResi] = "s2";
    } else if (recordName === "CONECT") {
      // MEMO: We don't have to parse SSBOND, LINK because both are
      // also
      // described in CONECT. But what about 2JYT???
      const from = parseInt(line.substring(6, 11));
      const fromindex = serialToIndex[from];
      const fromAtom = atoms[fromindex];
      const coffsets = [11, 16, 21, 26];
      for (let j = 0; j < 4; j++) {
        const to = parseInt(line.substring(coffsets[j], coffsets[j] + 5));
        const toindex = serialToIndex[to];
        let from_to = fromindex + ":" + toindex;
        const toAtom = atoms[toindex];
        if (fromAtom !== undefined && toAtom !== undefined) {
          // duplicated conect records indicate bond order
          if (!seenbonds[from_to]) {
            seenbonds[from_to] = 1;
            if (
              fromAtom.bonds.length == 0 ||
              fromAtom.bonds[fromAtom.bonds.length - 1] !== toindex
            ) {
              fromAtom.bonds.push(toindex);
              fromAtom.bondOrder.push(1);
            }
          } else {
            //update bond order
            seenbonds[from_to] += 1;

            for (let bi = 0; bi < fromAtom.bonds.length; bi++) {
              if (fromAtom.bonds[bi] == toindex) {
                const newbo = seenbonds[from_to];
                if (newbo >= 4) {
                  //aromatic
                  fromAtom.bondOrder[bi] = 1;
                } else {
                  fromAtom.bondOrder[bi] = newbo;
                }
              }
            }
          }
        }
      }
    } else if (recordName === "HELIX ") {
      startChain = line.substring(19, 20);
      startResi = parseInt(line.substring(21, 25));
      endResi = parseInt(line.substring(33, 37));
      if (!(startChain in sslookup)) {
        sslookup[startChain] = {};
      }
      sslookup[startChain][startResi] = "h1";
      for (let res = startResi + 1; res < endResi; res++) {
        sslookup[startChain][res] = "h";
      }
      sslookup[startChain][endResi] = "h2";
    } else if (
      !noAssembly &&
      recordName === "REMARK" &&
      line.substring(13, 18) === "BIOMT"
    ) {
      let n: number;
      let matrix = new Matrix4();
      for (n = 1; n <= 3; n++) {
        line = lines[i].replace(/^\s*/, "");
        if (parseInt(line.substring(18, 19)) == n) {
          // check for all
          // three lines
          // by matching #
          // @ end of
          // "BIOMT" to n
          matrix.elements[n - 1] = parseFloat(line.substring(23, 33));
          matrix.elements[n - 1 + 4] = parseFloat(line.substring(33, 43));
          matrix.elements[n - 1 + 8] = parseFloat(line.substring(43, 53));
          matrix.elements[n - 1 + 12] = parseFloat(line.substring(53));
          i++;
        } else {
          while (line.substring(13, 18) === "BIOMT") {
            i++;
            line = lines[i].replace(/^\s*/, "");
          }
        }
      }
      matrix.elements[3] = 0;
      matrix.elements[7] = 0;
      matrix.elements[11] = 0;
      matrix.elements[15] = 1;
      modelData.symmetries.push(matrix);
      i--; // set i back
    } else if (recordName === "CRYST1") {
      let a: number,
        b: number,
        c: number,
        alpha: number,
        beta: number,
        gamma: number;
      a = parseFloat(line.substring(7, 15));
      b = parseFloat(line.substring(16, 24));
      c = parseFloat(line.substring(25, 33));
      alpha = parseFloat(line.substring(34, 40));
      beta = parseFloat(line.substring(41, 47));
      gamma = parseFloat(line.substring(48, 54));
      modelData.cryst = {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
      };
    } else if (recordName === "ANISOU") {
      const serial = parseInt(line.substring(6, 11));
      const anisouAtomIndex = serialToIndex[serial];
      const anisouAtom = atoms[anisouAtomIndex];
      if (anisouAtom) {
        const vals = line.substring(30).trim().split(/\s+/);
        const uMat = {
          u11: parseInt(vals[0]),
          u22: parseInt(vals[1]),
          u33: parseInt(vals[2]),
          u12: parseInt(vals[3]),
          u13: parseInt(vals[4]),
          u23: parseInt(vals[5]),
        };

        anisouAtom.uMat = uMat;
      }
    }
  }

  //fix any "one-way" bonds in CONECT records
  validateBonds(atoms, serialToIndex);
  // assign bonds - yuck, can't count on connect records
  if (assignbonds) assignPDBBonds(atoms, options);

  if (!noAssembly)
    processSymmetries(modelData.symmetries, atoms, options, modelData.cryst);

  if (computeStruct && !ignoreStruct) {
    computeSecondaryStructure(atoms, options.hbondCutoff);
  }

  // Assign secondary structures from pdb file
  if (!isEmpty(sslookup)) {
    for (let i = 0; i < atoms.length; i++) {
      const atom = atoms[i];
      if (atom === undefined) continue;
      if (atom.chain in sslookup && atom.resi in sslookup[atom.chain]) {
        const code = sslookup[atom.chain][atom.resi];
        atom.ss = code[0];
        if (code.length > 1) {
          if (code[1] == "1") atom.ssbegin = true;
          else if (code[1] == "2") atom.ssend = true;
        }
      }
    }
  }
  //console.log("assign structure " + ((new Date()).getTime() - starttime));

  return [atoms, modelData, remainingLines];
}
