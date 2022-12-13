//make sure bonds are actually two way
export function validateBonds (atomsarray, serialToIndex) {
  for (var i = 0, n = atomsarray.length; i < n; i++) {
      var atom = atomsarray[i];
      for(var b = 0; b < atom.bonds.length; b++) {
          var a2i = atom.bonds[b];
          var atom2 = atomsarray[a2i];
          var atomi = serialToIndex[atom.serial];
          if(atom2 && atomi) {
              var a1i = atom2.bonds.indexOf(atomi);
              if(a1i < 0) {
                  atom2.bonds.push(atomi);
                  atom2.bondOrder.push(atom.bondOrder[b]);
              }
          }
      }
  }
};