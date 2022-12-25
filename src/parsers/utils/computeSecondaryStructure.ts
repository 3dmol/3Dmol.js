import { assignBackboneHBonds } from "./assignBackboneHBonds";

export function computeSecondaryStructure(atomsarray, hbondCutoff) {
  assignBackboneHBonds(atomsarray, hbondCutoff);

  // compute, per residue, what the secondary structure is
  var chres = {}; // lookup by chain and resid
  var i, il, c, r; // i: used in for loop, il: length of atomsarray
  var atom, val;

  //identify helices first
  for (i = 0, il = atomsarray.length; i < il; i++) {
    atom = atomsarray[i];

    if (typeof chres[atom.chain] === "undefined") chres[atom.chain] = [];

    if (isFinite(atom.hbondDistanceSq)) {
      var other = atom.hbondOther;
      if (typeof chres[other.chain] === "undefined") chres[other.chain] = [];

      if (Math.abs(other.resi - atom.resi) === 4) {
        // helix
        chres[atom.chain][atom.resi] = "h";
      }
    }
  }

  // plug gaps in helices
  for (c in chres) {
    for (r = 1; r < chres[c].length - 1; r++) {
      var valbefore = chres[c][r - 1];
      var valafter = chres[c][r + 1];
      val = chres[c][r];
      if (valbefore == "h" && valbefore == valafter && val != valbefore) {
        chres[c][r] = valbefore;
      }
    }
  }

  //now potential sheets - but only if mate not part of helix
  for (i = 0, il = atomsarray.length; i < il; i++) {
    atom = atomsarray[i];

    if (
      isFinite(atom.hbondDistanceSq) &&
      chres[atom.chain][atom.resi] != "h" &&
      atom.ss != "h"
    ) {
      chres[atom.chain][atom.resi] = "maybesheet";
    }
  }

  //sheets must bond to other sheets
  for (let i = 0, il = atomsarray.length; i < il; i++) {
    atom = atomsarray[i];

    if (
      isFinite(atom.hbondDistanceSq) &&
      chres[atom.chain][atom.resi] == "maybesheet"
    ) {
      let other = atom.hbondOther;
      let otherval = chres[other.chain][other.resi];
      if (otherval == "maybesheet" || otherval == "s") {
        // true sheet
        chres[atom.chain][atom.resi] = "s";
        chres[other.chain][other.resi] = "s";
      }
    }
  }

  // plug gaps in sheets and remove singletons
  for (let c in chres) {
    for (let r = 1; r < chres[c].length - 1; r++) {
      let valbefore = chres[c][r - 1];
      let valafter = chres[c][r + 1];
      val = chres[c][r];
      if (valbefore == "s" && valbefore == valafter && val != valbefore) {
        chres[c][r] = valbefore;
      }
    }
    for (let r = 0; r < chres[c].length; r++) {
      let val = chres[c][r];
      if (val == "h" || val == "s") {
        if (chres[c][r - 1] != val && chres[c][r + 1] != val)
          delete chres[c][r];
      }
    }
  }

  // assign to all atoms in residue, keep track of start
  for (i = 0, il = atomsarray.length; i < il; i++) {
    atom = atomsarray[i];
    val = chres[atom.chain][atom.resi];

    //clear hbondOther to eliminate circular references that prohibit serialization
    delete atom.hbondOther;
    delete atom.hbondDistanceSq;
    if (typeof val == "undefined" || val == "maybesheet") continue;
    atom.ss = val;
    if (chres[atom.chain][atom.resi - 1] != val) atom.ssbegin = true;
    if (chres[atom.chain][atom.resi + 1] != val) atom.ssend = true;
  }
}
