// This is optimized for proteins where it is assumed connected atoms are on the same or next residue

import { AtomSpec } from "specs";
import { areConnected } from "./areConnected";
import { assignBonds } from "./assignBonds";
import { standardResidues } from "./standardResidues";
import { ParserOptionsSpec } from "parsers/ParserOptionsSpec";


/**
 * @param {AtomSpec[]}
 *            atomsarray
 */

export function assignPDBBonds(atomsarray: AtomSpec[], options: ParserOptionsSpec) {
  // assign bonds - yuck, can't count on connect records
  const protatoms: Array<AtomSpec> = [];
  const hetatoms: Array<AtomSpec> = [];
  for (let i = 0, n = atomsarray.length; i < n; i++) {
    const atom = atomsarray[i];
    atom.index = i;
    if (atom.hetflag || !standardResidues.has(atom.resn)) hetatoms.push(atom);
    else protatoms.push(atom);
  }

  assignBonds(hetatoms, options);

  // sort by resid
  protatoms.sort(function (a, b) {
    if (a.chain !== b.chain) return a.chain < b.chain ? -1 : 1;
    return a.resi - b.resi;
  });

  // for identifying connected residues
  let currentResi = -1;
  let reschain = -1;
  let lastResConnected: boolean;

  for (let i = 0, n = protatoms.length; i < n; i++) {
    const ai = protatoms[i];

    if (ai.resi !== currentResi) {
      currentResi = ai.resi;
      if (!lastResConnected) reschain++;

      lastResConnected = false;
    }

    ai.reschain = reschain;

    for (let j = i + 1; j < protatoms.length; j++) {
      const aj = protatoms[j];
      if (aj.chain !== ai.chain || aj.resi - ai.resi > 1) break;

      if (areConnected(ai, aj, options)) {
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
