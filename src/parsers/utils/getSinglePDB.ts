import { Matrix4 } from "../../WebGL";
import { atomNameToElem } from "./atomNameToElem";
import { bondTable } from "./bondLength";
import { computeSecondaryStructure } from "./computeSecondaryStructure";
import { isEmpty } from "./isEmpty";
import { processSymmetries } from "./processSymmetries";
import { assignPDBBonds } from "./assignPDBBonds";
import { validateBonds } from "./validateBonds";

//return one model worth of pdb, returns atoms, modelData, and remaining lines
export function getSinglePDB(lines, options, sslookup) {
  var atoms: any[] = [];
  var assignbonds =
    options.assignBonds === undefined ? true : options.assignBonds;
  var noH = !options.keepH; // suppress hydrogens by default
  var ignoreStruct = !!options.noSecondaryStructure;
  var computeStruct = !options.noComputeSecondaryStructure;
  var noAssembly = !options.doAssembly; // don't assemble by default
  var selAltLoc = options.altLoc ? options.altLoc : "A"; //default alternate location to select if present
  var modelData: Record<string, any> = { symmetries: [] };
  var atom;
  var remainingLines = [];

  var hasStruct = false;
  var serialToIndex: number[] = []; // map from pdb serial to index in atoms
  var line;
  var seenbonds: Record<any, any> = {}; //sometimes connect records are duplicated as an unofficial means of relaying bond orders

  for (let i = 0; i < lines.length; i++) {
    line = lines[i].replace(/^\s*/, ""); // remove indent
    var recordName = line.substr(0, 6);
    var startChain, startResi, endChain, endResi;

    if (recordName.indexOf("END") == 0) {
      remainingLines = lines.slice(i + 1);
      if (recordName == "END") {
        //as opposed to ENDMDL
        //reset secondary structure
        for (var prop in sslookup) {
          if (sslookup.hasOwnProperty(prop)) {
            delete sslookup[prop];
          }
        }
      }
      break;
    } else if (recordName == "ATOM  " || recordName == "HETATM") {
      var resn, chain, resi, icode, x, y, z, hetflag, elem, serial, altLoc, b;
      altLoc = line.substr(16, 1);
      if (altLoc != " " && altLoc != selAltLoc && selAltLoc != "*") continue;
      serial = parseInt(line.substr(6, 5));
      atom = line.substr(12, 4).replace(/ /g, "");
      resn = line.substr(17, 3).replace(/ /g, "");
      chain = line.substr(21, 1);
      resi = parseInt(line.substr(22, 4));
      icode = line.substr(26, 1);
      x = parseFloat(line.substr(30, 8));
      y = parseFloat(line.substr(38, 8));
      z = parseFloat(line.substr(46, 8));
      b = parseFloat(line.substr(60, 8));
      elem = line.substr(76, 2).replace(/ /g, "");
      if (elem === "" || typeof bondTable[elem] === "undefined") {
        // for some incorrect PDB files
        elem = atomNameToElem(line.substr(12, 2), line[0] == "A");
      } else {
        elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();
      }

      if (elem == "H" && noH) continue;
      if (recordName[0] == "H") hetflag = true;
      else hetflag = false;
      serialToIndex[serial] = atoms.length;
      atoms.push({
        resn: resn,
        x: x,
        y: y,
        z: z,
        elem: elem,
        hetflag: hetflag,
        altLoc: altLoc,
        chain: chain,
        resi: resi,
        icode: icode,
        rescode: resi + (icode != " " ? "^" + icode : ""), // combo
        // resi
        // and
        // icode
        serial: serial,
        atom: atom,
        bonds: [],
        ss: "c",
        bondOrder: [],
        properties: {},
        b: b,
        pdbline: line,
      });
    } else if (recordName == "SHEET ") {
      hasStruct = true;
      startChain = line.substr(21, 1);
      startResi = parseInt(line.substr(22, 4));
      endChain = line.substr(32, 1);
      endResi = parseInt(line.substr(33, 4));
      if (!(startChain in sslookup)) {
        sslookup[startChain] = {};
      }
      //mark start and end with additional character
      sslookup[startChain][startResi] = "s1";
      for (var res = startResi + 1; res < endResi; res++) {
        sslookup[startChain][res] = "s";
      }
      sslookup[startChain][endResi] = "s2";
    } else if (recordName == "CONECT") {
      // MEMO: We don't have to parse SSBOND, LINK because both are
      // also
      // described in CONECT. But what about 2JYT???
      var from = parseInt(line.substr(6, 5));
      var fromindex = serialToIndex[from];
      var fromAtom = atoms[fromindex];
      var coffsets = [11, 16, 21, 26];
      for (let j = 0; j < 4; j++) {
        var to = parseInt(line.substr(coffsets[j], 5));
        var toindex = serialToIndex[to];
        let from_to = fromindex+":"+toindex;
        var toAtom = atoms[toindex];
        if (fromAtom !== undefined && toAtom !== undefined) {
          // duplicated conect records indicate bond order
          if (!seenbonds[from_to]) {
            seenbonds[from_to] = 1;
            if (
              fromAtom.bonds.length == 0 ||
              fromAtom.bonds[fromAtom.bonds.length - 1] != toindex
            ) {
              fromAtom.bonds.push(toindex);
              fromAtom.bondOrder.push(1);
            }
          } else {
            //update bond order
            seenbonds[from_to] += 1;

            for (let bi = 0; bi < fromAtom.bonds.length; bi++) {
              if (fromAtom.bonds[bi] == toindex) {
                var newbo = seenbonds[from_to];
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
    } else if (recordName == "HELIX ") {
      hasStruct = true;
      startChain = line.substr(19, 1);
      startResi = parseInt(line.substr(21, 4));
      endChain = line.substr(31, 1);
      endResi = parseInt(line.substr(33, 4));
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
      recordName == "REMARK" &&
      line.substr(13, 5) == "BIOMT"
    ) {
      var n;
      var matrix = new Matrix4();
      for (n = 1; n <= 3; n++) {
        line = lines[i].replace(/^\s*/, "");
        if (parseInt(line.substr(18, 1)) == n) {
          // check for all
          // three lines
          // by matching #
          // @ end of
          // "BIOMT" to n
          matrix.elements[n - 1] = parseFloat(line.substr(23, 10));
          matrix.elements[n - 1 + 4] = parseFloat(line.substr(33, 10));
          matrix.elements[n - 1 + 8] = parseFloat(line.substr(43, 10));
          matrix.elements[n - 1 + 12] = parseFloat(line.substr(53));
          i++;
        } else {
          while (line.substr(13, 5) == "BIOMT") {
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
    } else if (recordName == "CRYST1") {
      let a, b, c, alpha, beta, gamma;
      a = parseFloat(line.substr(7, 8));
      b = parseFloat(line.substr(16, 8));
      c = parseFloat(line.substr(25, 8));
      alpha = parseFloat(line.substr(34, 6));
      beta = parseFloat(line.substr(41, 6));
      gamma = parseFloat(line.substr(48, 6));
      modelData.cryst = {
        a: a,
        b: b,
        c: c,
        alpha: alpha,
        beta: beta,
        gamma: gamma,
      };
    } else if (recordName == "ANISOU") {
      let serial = parseInt(line.substr(6, 5));
      var anisouAtomIndex = serialToIndex[serial];
      var anisouAtom = atoms[anisouAtomIndex];
      if (anisouAtom) {
        var vals = line.substr(30).trim().split(/\s+/);
        var uMat = {
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
  if (assignbonds) assignPDBBonds(atoms);

  if (!noAssembly)
    processSymmetries(modelData.symmetries, atoms, options, modelData.cryst);

  if (computeStruct && !ignoreStruct) {
    computeSecondaryStructure(atoms, options.hbondCutoff);
  }

  // Assign secondary structures from pdb file
  if (!isEmpty(sslookup)) {
    for (let i = 0; i < atoms.length; i++) {
      atom = atoms[i];
      if (atom === undefined) continue;
      if (atom.chain in sslookup && atom.resi in sslookup[atom.chain]) {
        var code = sslookup[atom.chain][atom.resi];
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
