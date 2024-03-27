import { AtomSpec } from "specs";

// Make sure bonds are actually two way
export function validateBonds(atomsarray: AtomSpec[], serialToIndex: number[]) {
  for (let i = 0, n = atomsarray.length; i < n; i++) {
    const atom = atomsarray[i];
    for (let b = 0; b < atom.bonds.length; b++) {
      const a2i = atom.bonds[b];
      const atom2 = atomsarray[a2i];
      const atomi = serialToIndex[atom.serial];
      if (atom2 && atomi) {
        const a1i = atom2.bonds.indexOf(atomi);
        if (a1i < 0) {
          atom2.bonds.push(atomi);
          atom2.bondOrder.push(atom.bondOrder[b]);
        }
      }
    }
  }
}
