// this is optimized for proteins where it is assumed connected
// atoms are on the same or next residue

import { areConnected } from "./areConnected";
import { assignBonds } from "./assignBonds";
import { standardResidues } from "./standardResidues";

/*
 * @param {AtomSpec[]}
 *            atomsarray
 */
export function assignPDBBonds(atomsarray) {
  // assign bonds - yuck, can't count on connect records
  var protatoms: any[] = [];
  var hetatoms: any[] = [];
  var i, n;
  for (i = 0, n = atomsarray.length; i < n; i++) {
    var atom: any = atomsarray[i];
    atom.index = i;
    if (atom.hetflag || !standardResidues.has(atom.resn)) hetatoms.push(atom);
    else protatoms.push(atom);
  }

  assignBonds(hetatoms);

  // sort by resid
  protatoms.sort(function (a: any, b: any) {
    if (a.chain != b.chain) return a.chain < b.chain ? -1 : 1;
    return a.resi - b.resi;
  });

  // for identifying connected residues
  var currentResi = -1;
  var reschain = -1;
  var lastResConnected;

  for (i = 0, n = protatoms.length; i < n; i++) {
    var ai = protatoms[i];

    if (ai.resi !== currentResi) {
      currentResi = ai.resi;
      if (!lastResConnected) reschain++;

      lastResConnected = false;
    }

    ai.reschain = reschain;

    for (var j = i + 1; j < protatoms.length; j++) {
      var aj = protatoms[j];
      if (aj.chain != ai.chain) break;
      if (aj.resi - ai.resi > 1)
        // can't be connected
        break;
      if (areConnected(ai, aj)) {
        if (ai.bonds.indexOf(aj.index) === -1) {
          // only add if not already there
          ai.bonds.push(aj.index);
          ai.bondOrder.push(1);
          aj.bonds.push(ai.index);
          aj.bondOrder.push(1);
        }

        if (ai.resi !== aj.resi) lastResConnected = true;
      }
    }
  }
}
