// this will identify all hydrogen bonds between backbone
// atoms; assume atom names are correct, only identifies
// single closest hbond
export function assignBackboneHBonds(atomsarray, hbondCutoff) {
  let maxlength = hbondCutoff || 3.2;
  let maxlengthSq = maxlength*maxlength;
  let atoms = [];

  for (let i = 0, n = atomsarray.length; i < n; i++) {
    atomsarray[i].index = i;
    // only consider 'N' and 'O'
    var atom = atomsarray[i];
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
    var ai = atoms[i];

    for (let j = i + 1; j < n; j++) {
      var aj = atoms[j];
      var zdiff = aj.z - ai.z;
      if (zdiff > maxlength)
        // can't be connected
        break;
      if (aj.atom == ai.atom) continue; // can't be connected, but later might be
      var ydiff = Math.abs(aj.y - ai.y);
      if (ydiff > maxlength) continue;
      var xdiff = Math.abs(aj.x - ai.x);
      if (xdiff > maxlength) continue;
      var dist = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
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
