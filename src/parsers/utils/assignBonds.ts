import { areConnected } from "./areConnected";

/*
 * @param {AtomSpec[]} atoms
 */
export function assignBonds(atoms) {
  // assign bonds - yuck, can't count on connect records

  for (var i = 0, n = atoms.length; i < n; i++) {
    // Don't reindex if atoms are already indexed
    if (!atoms[i].index) atoms[i].index = i;
  }

  var grid = {};
  var MAX_BOND_LENGTH = 4.95; // (largest bond length, Cs) 2.25 * 2 * 1.1 (fudge factor)

  for (var index = 0; index < atoms.length; index++) {
    var atom = atoms[index];
    var x = Math.floor(atom.x / MAX_BOND_LENGTH);
    var y = Math.floor(atom.y / MAX_BOND_LENGTH);
    var z = Math.floor(atom.z / MAX_BOND_LENGTH);
    if (!grid[x]) {
      grid[x] = {};
    }
    if (!grid[x][y]) {
      grid[x][y] = {};
    }
    if (!grid[x][y][z]) {
      grid[x][y][z] = [];
    }

    grid[x][y][z].push(atom);
  }

  var findConnections = function (points, otherPoints) {
    for (var i = 0; i < points.length; i++) {
      var atom1 = points[i];
      for (var j = 0; j < otherPoints.length; j++) {
        var atom2 = otherPoints[j];

        if (areConnected(atom1, atom2)) {
          //gracefully handle one-sided bonds
          var a2i = atom1.bonds.indexOf(atom2.index);
          var a1i = atom2.bonds.indexOf(atom1.index);
          if (a2i == -1 && a1i == -1) {
            atom1.bonds.push(atom2.index);
            atom1.bondOrder.push(1);
            atom2.bonds.push(atom1.index);
            atom2.bondOrder.push(1);
          } else if (a2i == -1) {
            atom1.bonds.push(atom2.index);
            atom1.bondOrder.push(atom2.bondOrder[a1i]);
          } else if (a1i == -1) {
            atom2.bonds.push(atom1.index);
            atom2.bondOrder.push(atom1.bondOrder[a2i]);
          }
        }
      }
    }
  };

  /*const*/ var OFFSETS = [
    { x: 0, y: 0, z: 1 },
    { x: 0, y: 1, z: -1 },
    { x: 0, y: 1, z: 0 },
    { x: 0, y: 1, z: 1 },
    { x: 1, y: -1, z: -1 },
    { x: 1, y: -1, z: 0 },
    { x: 1, y: -1, z: 1 },
    { x: 1, y: 0, z: -1 },
    { x: 1, y: 0, z: 0 },
    { x: 1, y: 0, z: 1 },
    { x: 1, y: 1, z: -1 },
    { x: 1, y: 1, z: 0 },
    { x: 1, y: 1, z: 1 },
  ];
  for (let xg in grid) {
    let x = parseInt(xg);
    for (let yg in grid[x]) {
      let y = parseInt(yg);
      for (let zg in grid[x][y]) {
        let z = parseInt(zg);
        let points = grid[x][y][z];

        for (let i = 0; i < points.length; i++) {
          let atom1 = points[i];
          for (let j = i + 1; j < points.length; j++) {
            let atom2 = points[j];
            if (areConnected(atom1, atom2)) {
              if (atom1.bonds.indexOf(atom2.index) == -1) {
                atom1.bonds.push(atom2.index);
                atom1.bondOrder.push(1);
                atom2.bonds.push(atom1.index);
                atom2.bondOrder.push(1);
              }
            }
          }
        }

        for (let o = 0; o < OFFSETS.length; o++) {
          let offset = OFFSETS[o];
          if (
            !grid[x + offset.x] ||
            !grid[x + offset.x][y + offset.y] ||
            !grid[x + offset.x][y + offset.y][z + offset.z]
          )
            continue;

          let otherPoints = grid[x + offset.x][y + offset.y][z + offset.z];
          findConnections(points, otherPoints);
        }
      }
    }
  }
}
