// This will identify all hydrogen bonds between backbone
// atoms; assume atom names are correct, only identifies
// single closest hbond

import { AtomSpec } from "specs";
// interface Atoms {index: number; atom: string; hbondDistanceSq: number; hbondOther: any; hetflag:any}
export function assignBackboneHBonds(
  atomsarray: Array<AtomSpec>,
  hbondCutoff: number
) {
  const maxlength = hbondCutoff || 3.2;
  const maxlengthSq = maxlength * maxlength;
  const atoms = [];

  for (let i = 0, n = atomsarray.length; i < n; i++) {
    atomsarray[i].index = i;
    // only consider 'N' and 'O'
    const atom = atomsarray[i];
    if (!atom.hetflag && (atom.atom === "N" || atom.atom === "O")) {
      atoms.push(atom);
      atom.hbondOther = null;
      atom.hbondDistanceSq = Number.POSITIVE_INFINITY;
    }
  }

  atoms.sort(function (a, b) {
    return a.z - b.z;
  });
  for (let i = 0, n = atoms.length; i < n; i++) {
    const ai = atoms[i];

    for (let j = i + 1; j < n; j++) {
      const aj = atoms[j];
      const zdiff = aj.z - ai.z;
      if (zdiff > maxlength)
        // can't be connected
        break;
      if (aj.atom == ai.atom) continue; // can't be connected, but later might be
      const ydiff = Math.abs(aj.y - ai.y);
      if (ydiff > maxlength) continue;
      const xdiff = Math.abs(aj.x - ai.x);
      if (xdiff > maxlength) continue;
      const dist = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
      if (dist > maxlengthSq) continue;

      if (aj.chain == ai.chain && Math.abs(aj.resi - ai.resi) < 4) continue; // ignore bonds between too close residues
      // select closest hbond
      if (dist < ai.hbondDistanceSq) {
        ai.hbondOther = aj;
        ai.hbondDistanceSq = dist;
      }
      if (dist < aj.hbondDistanceSq) {
        aj.hbondOther = ai;
        aj.hbondDistanceSq = dist;
      }
    }
  }
}
