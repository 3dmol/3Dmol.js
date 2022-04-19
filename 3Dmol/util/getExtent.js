/**
 * computes the bounding box around the provided atoms
 * @param {import("../specs").AtomSpec[]} atomlist
 * @param {boolean} [ignoreSymmetries]
 * @return {Array}
 */
export default function getExtent(atomlist, ignoreSymmetries) {
  let xmin;
  let ymin;
  let zmin;
  let xmax;
  let ymax;
  let zmax;
  let xsum;
  let ysum;
  let zsum;
  let cnt;
  const includeSym = !ignoreSymmetries;

  xmin = ymin = zmin = 9999;
  xmax = ymax = zmax = -9999;
  xsum = ysum = zsum = cnt = 0;

  if (atomlist.length === 0)
    return [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
  for (let i = 0; i < atomlist.length; i++) {
    const atom = atomlist[i];
    if (
      typeof atom == 'undefined' ||
      !Number.isFinite(atom.x) ||
      !Number.isFinite(atom.y) ||
      !Number.isFinite(atom.z)
    )
      continue;
    cnt += 1;
    xsum += atom.x;
    ysum += atom.y;
    zsum += atom.z;

    xmin = xmin < atom.x ? xmin : atom.x;
    ymin = ymin < atom.y ? ymin : atom.y;
    zmin = zmin < atom.z ? zmin : atom.z;
    xmax = xmax > atom.x ? xmax : atom.x;
    ymax = ymax > atom.y ? ymax : atom.y;
    zmax = zmax > atom.z ? zmax : atom.z;

    if (atom.symmetries && includeSym) {
      for (let n = 0; n < atom.symmetries.length; n++) {
        cnt += 1;
        xsum += atom.symmetries[n].x;
        ysum += atom.symmetries[n].y;
        zsum += atom.symmetries[n].z;
        xmin = xmin < atom.symmetries[n].x ? xmin : atom.symmetries[n].x;
        ymin = ymin < atom.symmetries[n].y ? ymin : atom.symmetries[n].y;
        zmin = zmin < atom.symmetries[n].z ? zmin : atom.symmetries[n].z;
        xmax = xmax > atom.symmetries[n].x ? xmax : atom.symmetries[n].x;
        ymax = ymax > atom.symmetries[n].y ? ymax : atom.symmetries[n].y;
        zmax = zmax > atom.symmetries[n].z ? zmax : atom.symmetries[n].z;
      }
    }
  }

  return [
    [xmin, ymin, zmin],
    [xmax, ymax, zmax],
    [xsum / cnt, ysum / cnt, zsum / cnt],
  ];
}
