/* eslint-disable max-classes-per-file */
import MarchingCube from "./MarchingCube";
/* 
  There are many // @ts-ignore's in this file because the data arrays are nullable
*/
/*  ProteinSurface.js by biochem_fan

Ported and modified for Javascript based on EDTSurf,
  whose license is as follows.

Permission to use, copy, modify, and distribute this program for any
purpose, with or without fee, is hereby granted, provided that this
copyright notice and the reference information appear in all copies or
substantial portions of the Software. It is provided "as is" without
express or implied warranty. 

Reference:
http://zhanglab.ccmb.med.umich.edu/EDTSurf/
D. Xu, Y. Zhang (2009) Generating Triangulated Macromolecular Surfaces
by Euclidean Distance Transform. PLoS ONE 4(12): e8140.

=======

TODO: Improved performance on Firefox
      Reduce memory consumption
      Refactor!
 */

// dkoes
// Surface calculations.  This must be safe to use within a web worker.
if (typeof console === 'undefined') {
  // this should only be true inside of a webworker
  // @ts-ignore
  console = {
    log() {},
  };
}

// constants for vpbits bitmasks
/** @const */
const INOUT = 1;
/** @const */
const ISDONE = 2;
/** @const */
const ISBOUND = 4;
const probeRadius = 1.4;
const defaultScaleFactor = 2;
const vdwRadii = {
  H: 1.2,
  Li: 1.82,
  Na: 2.27,
  K: 2.75,
  C: 1.7,
  N: 1.55,
  O: 1.52,
  F: 1.47,
  P: 1.8,
  S: 1.8,
  CL: 1.75,
  BR: 1.85,
  SE: 1.9,
  ZN: 1.39,
  CU: 1.4,
  NI: 1.63,
  X: 2,
};

// a little class for 3d array, should really generalize this and
// use throughout...
class PointGrid {
  constructor(length, width, height) {
    // the standard says this is zero initialized
    this.data = new Int32Array(length * width * height * 3);
    this.width = width;
    this.height = height;
  }

  // set position x,y,z to pt, which has ix,iy,and iz
  set(x, y, z, pt) {
    const index = ((x * this.width + y) * this.height + z) * 3;
    this.data[index] = pt.ix;
    this.data[index + 1] = pt.iy;
    this.data[index + 2] = pt.iz;
  }

  get(x, y, z) {
    const index = ((x * this.width + y) * this.height + z) * 3;
    return {
      ix: this.data[index],
      iy: this.data[index + 1],
      iz: this.data[index + 2],
    };
  }
}

export default class ProteinSurface {
  ptranx = 0;
  ptrany = 0;
  ptranz = 0;
  // also have to adjust offset used to find non-shown
  // atoms
  pHeight = 0;
  pWidth = 0;
  pLength = 0;
  cutRadius = 0;
  /** @type {Uint8Array|null} */
  vpBits = null; // uint8 array of bitmasks
  /** @type {Float64Array|null} */
  vpDistance = null; // floatarray of _squared_ distances
  /** @type {Int32Array|null} */
  vpAtomID = null; // intarray
  pminx = 0;
  pminy = 0;
  pminz = 0;
  pmaxx = 0;
  pmaxy = 0;
  pmaxz = 0;

  depty = {};
  widxz = {};
  faces;
  verts;
  nb = [
    new Int32Array([1, 0, 0]),
    new Int32Array([-1, 0, 0]),
    new Int32Array([0, 1, 0]),
    new Int32Array([0, -1, 0]),
    new Int32Array([0, 0, 1]),
    new Int32Array([0, 0, -1]),
    new Int32Array([1, 1, 0]),
    new Int32Array([1, -1, 0]),
    new Int32Array([-1, 1, 0]),
    new Int32Array([-1, -1, 0]),
    new Int32Array([1, 0, 1]),
    new Int32Array([1, 0, -1]),
    new Int32Array([-1, 0, 1]),
    new Int32Array([-1, 0, -1]),
    new Int32Array([0, 1, 1]),
    new Int32Array([0, 1, -1]),
    new Int32Array([0, -1, 1]),
    new Int32Array([0, -1, -1]),
    new Int32Array([1, 1, 1]),
    new Int32Array([1, 1, -1]),
    new Int32Array([1, -1, 1]),
    new Int32Array([-1, 1, 1]),
    new Int32Array([1, -1, -1]),
    new Int32Array([-1, -1, 1]),
    new Int32Array([-1, 1, -1]),
    new Int32Array([-1, -1, -1]),
  ];

  scaleFactor = defaultScaleFactor; // 2 is .5A grid; if this is made user configurable,


  boundPoint;

  /** @param {import("./specs").AtomSpec} atom */
  static getVDWIndex(atom) {
    if (!atom.elem || typeof vdwRadii[atom.elem] == 'undefined') {
      return 'X';
    }
    return atom.elem;
  }

  getFacesAndVertices(atomlist) {
    const atomsToShow = {};
    for (let i = 0, il = atomlist.length; i < il; i++) atomsToShow[atomlist[i]] = true;
    const vertices = this.verts;
    for (let i = 0, il = vertices.length; i < il; i++) {
      vertices[i].x = vertices[i].x / this.scaleFactor - this.ptranx;
      vertices[i].y = vertices[i].y / this.scaleFactor - this.ptrany;
      vertices[i].z = vertices[i].z / this.scaleFactor - this.ptranz;
    }

    const finalfaces = [];
    for (let i = 0, il = this.faces.length; i < il; i += 3) {
      // var f = faces[i];
      const fa = this.faces[i];
      const fb = this.faces[i + 1];
      const fc = this.faces[i + 2];
      const a = vertices[fa].atomid;
      const b = vertices[fb].atomid;
      const c = vertices[fc].atomid;

      // must be a unique face for each atom
      let which = a;
      if (b < which) which = b;
      if (c < which) which = c;
      if (!atomsToShow[which]) {
        continue;
      }
      //            var av = vertices[faces[i]];
      //            var bv = vertices[faces[i+1]];
      //            var cv = vertices[faces[i+2]];
      if (fa !== fb && fb !== fc && fa !== fc) {
        finalfaces.push(fa);
        finalfaces.push(fb);
        finalfaces.push(fc);
      }
    }

    // try to help the garbage collector
    this.vpBits = null; // uint8 array of bitmasks
    this.vpDistance = null; // floatarray
    this.vpAtomID = null; // intarray

    return {
      vertices,
      faces: finalfaces,
    };
  }

  initparm(extent, btype, volume) {
    if (volume > 1000000)
      // heuristical decrease resolution to avoid large memory consumption
      this.scaleFactor = defaultScaleFactor / 2;

    const margin = (1 / this.scaleFactor) * 5.5; // need margin to avoid

    // boundary/round off effects
    this.pminx = extent[0][0];
    this.pmaxx = extent[1][0];
    this.pminy = extent[0][1];
    this.pmaxy = extent[1][1];
    this.pminz = extent[0][2];
    this.pmaxz = extent[1][2];

    if (!btype) {
      this.pminx -= margin;
      this.pminy -= margin;
      this.pminz -= margin;
      this.pmaxx += margin;
      this.pmaxy += margin;
      this.pmaxz += margin;
    } else {
      this.pminx -= probeRadius + margin;
      this.pminy -= probeRadius + margin;
      this.pminz -= probeRadius + margin;
      this.pmaxx += probeRadius + margin;
      this.pmaxy += probeRadius + margin;
      this.pmaxz += probeRadius + margin;
    }

    this.pminx = Math.floor(this.pminx * this.scaleFactor) / this.scaleFactor;
    this.pminy = Math.floor(this.pminy * this.scaleFactor) / this.scaleFactor;
    this.pminz = Math.floor(this.pminz * this.scaleFactor) / this.scaleFactor;
    this.pmaxx = Math.ceil(this.pmaxx * this.scaleFactor) / this.scaleFactor;
    this.pmaxy = Math.ceil(this.pmaxy * this.scaleFactor) / this.scaleFactor;
    this.pmaxz = Math.ceil(this.pmaxz * this.scaleFactor) / this.scaleFactor;

    this.ptranx = -this.pminx;
    this.ptrany = -this.pminy;
    this.ptranz = -this.pminz;

    this.pLength = Math.ceil(this.scaleFactor * (this.pmaxx - this.pminx)) + 1;
    this.pWidth = Math.ceil(this.scaleFactor * (this.pmaxy - this.pminy)) + 1;
    this.pHeight = Math.ceil(this.scaleFactor * (this.pmaxz - this.pminz)) + 1;

    this.boundingatom(btype);
    this.cutRadius = probeRadius * this.scaleFactor;

    this.vpBits = new Uint8Array(this.pLength * this.pWidth * this.pHeight);
    this.vpDistance = new Float64Array(this.pLength * this.pWidth * this.pHeight); // float 32

    // doesn't
    // play
    // nicely
    // with
    // native
    // floats
    this.vpAtomID = new Int32Array(this.pLength * this.pWidth * this.pHeight);
    // console.log("Box size: ", pLength, pWidth, pHeight, vpBits.length);
  }

  boundingatom(btype) {
    const tradius = [];
    let txz;
    let tdept;
    let sradius;

    for (const i in vdwRadii) {
      if (!vdwRadii.hasOwnProperty(i)) continue;
      const r = vdwRadii[i];
      if (!btype) tradius[i] = r * this.scaleFactor + 0.5;
      else tradius[i] = (r + probeRadius) * this.scaleFactor + 0.5;

      sradius = tradius[i] * tradius[i];
      this.widxz[i] = Math.floor(tradius[i]) + 1;
      this.depty[i] = new Int32Array(this.widxz[i] * this.widxz[i]);
      let indx = 0;
      for (let j = 0; j < this.widxz[i]; j++) {
        for (let k = 0; k < this.widxz[i]; k++) {
          txz = j * j + k * k;
          if (txz > sradius) this.depty[i][indx] = -1; // outside
          else {
            tdept = Math.sqrt(sradius - txz);
            this.depty[i][indx] = Math.floor(tdept);
          }
          indx += 1;
        }
      }
    }
  }

  fillvoxels(atoms, atomlist) {
    // seqterm,bool
    // atomtype,atom*
    // proseq,bool bcolor)
    // @ts-ignore
    for (let i = 0, il = this.vpBits.length; i < il; i++) {
      // @ts-ignore
      this.vpBits[i] = 0;
      // @ts-ignore
      this.vpDistance[i] = -1.0;
      // @ts-ignore
      this.vpAtomID[i] = -1;
    }

    for (const i in atomlist) {
      const atom = atoms[atomlist[i]];
      if (atom === undefined) continue;
      this.fillAtom(atom, atoms);
    }

    // @ts-ignore
    for (let i = 0, il = this.vpBits.length; i < il; i++) if (this.vpBits[i] & INOUT) this.vpBits[i] |= ISDONE;
  }

  fillAtom(atom, atoms) {
    let ox;
    let oy;
    let oz;
    let mi;
    let mj;
    let mk;
    let i;
    let j;
    let k;
    let si;
    let sj;
    let sk;
    let ii;
    let jj;
    let kk;
    let n;
    const cx = Math.floor(0.5 + this.scaleFactor * (atom.x + this.ptranx));
    const cy = Math.floor(0.5 + this.scaleFactor * (atom.y + this.ptrany));
    const cz = Math.floor(0.5 + this.scaleFactor * (atom.z + this.ptranz));

    const at = ProteinSurface.getVDWIndex(atom);
    let nind = 0;
    const pWH = this.pWidth * this.pHeight;

    for (i = 0, n = this.widxz[at]; i < n; i++) {
      for (j = 0; j < n; j++) {
        if (this.depty[at][nind] !== -1) {
          for (ii = -1; ii < 2; ii++) {
            for (jj = -1; jj < 2; jj++) {
              for (kk = -1; kk < 2; kk++) {
                if (ii !== 0 && jj !== 0 && kk !== 0) {
                  mi = ii * i;
                  mk = kk * j;
                  for (k = 0; k <= this.depty[at][nind]; k++) {
                    mj = k * jj;
                    si = cx + mi;
                    sj = cy + mj;
                    sk = cz + mk;
                    if (
                      si < 0 ||
                      sj < 0 ||
                      sk < 0 ||
                      si >= this.pLength ||
                      sj >= this.pWidth ||
                      sk >= this.pHeight
                    )
                      continue;
                    const index = si * pWH + sj * this.pHeight + sk;

                    // @ts-ignore
                    if (!(this.vpBits[index] & INOUT)) {
                      // @ts-ignore
                      this.vpBits[index] |= INOUT;
                      // @ts-ignore
                      this.vpAtomID[index] = atom.serial;
                    } else {
                      // @ts-ignore
                      const atom2 = atoms[this.vpAtomID[index]];
                      if (atom2.serial !== atom.serial) {
                        ox = cx + mi - Math.floor(0.5 + this.scaleFactor * (atom2.x + this.ptranx));
                        oy = cy + mj - Math.floor(0.5 + this.scaleFactor * (atom2.y + this.ptrany));
                        oz = cz + mk - Math.floor(0.5 + this.scaleFactor * (atom2.z + this.ptranz));
                        if (mi * mi + mj * mj + mk * mk < ox * ox + oy * oy + oz * oz)
                          // @ts-ignore
                          this.vpAtomID[index] = atom.serial;
                      }
                    }
                  } // k
                } // if
              } // kk
            } // jj
          } // ii
        } // if
        nind += 1;
      } // j
    } // i
  }

  fillvoxelswaals(atoms, atomlist) {
    // @ts-ignore
    for (let i = 0, il = this.vpBits.length; i < il; i++) this.vpBits[i] &= ~ISDONE; // not isdone

    for (const i in atomlist) {
      const atom = atoms[atomlist[i]];
      if (atom === undefined) continue;

      this.fillAtomWaals(atom, atoms);
    }
  }

  fillAtomWaals(atom, atoms) {
    let ox;
    let oy;
    let oz;
    let nind = 0;
    let mi;
    let mj;
    let mk;
    let si;
    let sj;
    let sk;
    let i;
    let j;
    let k;
    let ii;
    let jj;
    let kk;
    let n;
    const cx = Math.floor(0.5 + this.scaleFactor * (atom.x + this.ptranx));
    const cy = Math.floor(0.5 + this.scaleFactor * (atom.y + this.ptrany));
    const cz = Math.floor(0.5 + this.scaleFactor * (atom.z + this.ptranz));

    const at = ProteinSurface.getVDWIndex(atom);
    const pWH = this.pWidth * this.pHeight;
    for (i = 0, n = this.widxz[at]; i < n; i++) {
      for (j = 0; j < n; j++) {
        if (this.depty[at][nind] !== -1) {
          for (ii = -1; ii < 2; ii++) {
            for (jj = -1; jj < 2; jj++) {
              for (kk = -1; kk < 2; kk++) {
                if (ii !== 0 && jj !== 0 && kk !== 0) {
                  mi = ii * i;
                  mk = kk * j;
                  for (k = 0; k <= this.depty[at][nind]; k++) {
                    mj = k * jj;
                    si = cx + mi;
                    sj = cy + mj;
                    sk = cz + mk;
                    if (
                      si < 0 ||
                      sj < 0 ||
                      sk < 0 ||
                      si >= this.pLength ||
                      sj >= this.pWidth ||
                      sk >= this.pHeight
                    )
                      continue;
                    const index = si * pWH + sj * this.pHeight + sk;
                    // @ts-ignore
                    if (!(this.vpBits[index] & ISDONE)) {
                      // @ts-ignore
                      this.vpBits[index] |= ISDONE;
                      // @ts-ignore
                      this.vpAtomID[index] = atom.serial;
                    } else {
                      // @ts-ignore
                      const atom2 = atoms[this.vpAtomID[index]];
                      if (atom2.serial !== atom.serial) {
                        ox = cx + mi - Math.floor(0.5 + this.scaleFactor * (atom2.x + this.ptranx));
                        oy = cy + mj - Math.floor(0.5 + this.scaleFactor * (atom2.y + this.ptrany));
                        oz = cz + mk - Math.floor(0.5 + this.scaleFactor * (atom2.z + this.ptranz));
                        if (mi * mi + mj * mj + mk * mk < ox * ox + oy * oy + oz * oz)
                          // @ts-ignore
                          this.vpAtomID[index] = atom.serial;
                      }
                    }
                  } // k
                } // if
              } // kk
            } // jj
          } // ii
        } // if
        nind += 1;
      } // j
    } // i
  }

  buildboundary() {
    const pWH = this.pWidth * this.pHeight;
    for (let i = 0; i < this.pLength; i++) {
      for (let j = 0; j < this.pHeight; j++) {
        for (let k = 0; k < this.pWidth; k++) {
          const index = i * pWH + k * this.pHeight + j;
          // @ts-ignore
          if (this.vpBits[index] & INOUT) {
            let ii = 0;
            while (ii < 26) {
              const ti = i + this.nb[ii][0];
              const tj = j + this.nb[ii][2];
              const tk = k + this.nb[ii][1];
              if (
                ti > -1 &&
                ti < this.pLength &&
                tk > -1 &&
                tk < this.pWidth &&
                tj > -1 &&
                tj < this.pHeight &&
                // @ts-ignore
                !(this.vpBits[ti * pWH + tk * this.pHeight + tj] & INOUT)
              ) {
                // @ts-ignore
                this.vpBits[index] |= ISBOUND;
                break;
              } else ii += 1;
            }
          }
        }
      }
    }
  }

  fastdistancemap() {
    let i;
    let j;
    let k;
    let n;

    /** @type {PointGrid|null} */
    let boundPoint = new PointGrid(this.pLength, this.pWidth, this.pHeight);
    const pWH = this.pWidth * this.pHeight;
    const cutRSq = this.cutRadius * this.cutRadius;

    let inarray = [];
    let outarray = [];

    let index;

    for (i = 0; i < this.pLength; i++) {
      for (j = 0; j < this.pWidth; j++) {
        for (k = 0; k < this.pHeight; k++) {
          index = i * pWH + j * this.pHeight + k;
          // @ts-ignore
          this.vpBits[index] &= ~ISDONE; // isdone = false
          // @ts-ignore
          if (this.vpBits[index] & INOUT) {
            // @ts-ignore
            if (this.vpBits[index] & ISBOUND) {
              const triple = {
                ix: i,
                iy: j,
                iz: k,
              };
              boundPoint.set(i, j, k, triple);
              inarray.push(triple);
              // @ts-ignore
              this.vpDistance[index] = 0;
              // @ts-ignore
              this.vpBits[index] |= ISDONE;
              // @ts-ignore
              this.vpBits[index] &= ~ISBOUND;
            }
          }
        }
      }
    }

    do {
      outarray = this.fastoneshell(inarray, boundPoint);
      inarray = [];
      for (i = 0, n = outarray.length; i < n; i++) {
        index = pWH * outarray[i].ix + this.pHeight * outarray[i].iy + outarray[i].iz;
        // @ts-ignore
        this.vpBits[index] &= ~ISBOUND;
        // @ts-ignore
        if (this.vpDistance[index] <= 1.0404 * cutRSq) {
          inarray.push({
            ix: outarray[i].ix,
            iy: outarray[i].iy,
            iz: outarray[i].iz,
          });
        }
      }
    } while (inarray.length !== 0);

    inarray = [];
    outarray = [];
    boundPoint = null;

    let cutsf = this.scaleFactor - 0.5;
    if (cutsf < 0) cutsf = 0;
    const cutoff = cutRSq - 0.5 / (0.1 + cutsf);
    for (i = 0; i < this.pLength; i++) {
      for (j = 0; j < this.pWidth; j++) {
        for (k = 0; k < this.pHeight; k++) {
          index = i * pWH + j * this.pHeight + k;
          // @ts-ignore
          this.vpBits[index] &= ~ISBOUND;
          // ses solid
          // @ts-ignore
          if (this.vpBits[index] & INOUT) {
            if (
              // @ts-ignore
              !(this.vpBits[index] & ISDONE) ||
              // @ts-ignore
              (this.vpBits[index] & ISDONE && this.vpDistance[index] >= cutoff)
            ) {
              // @ts-ignore
              this.vpBits[index] |= ISBOUND;
            }
          }
        }
      }
    }
  };

  fastoneshell(inarray, boundPoint) {
    // *allocout,voxel2
    // ***boundPoint, int*
    // outnum, int *elimi)
    let tx;
    let ty;
    let tz;
    let dx;
    let dy;
    let dz;
    let i;
    let j;
    let n;
    let square;
    let bp;
    let index;
    const outarray = [];
    if (inarray.length === 0) return outarray;

    const tnv = {
      ix: -1,
      iy: -1,
      iz: -1,
    };
    const pWH = this.pWidth * this.pHeight;
    for (i = 0, n = inarray.length; i < n; i++) {
      tx = inarray[i].ix;
      ty = inarray[i].iy;
      tz = inarray[i].iz;
      bp = boundPoint.get(tx, ty, tz);

      for (j = 0; j < 6; j++) {
        tnv.ix = tx + this.nb[j][0];
        tnv.iy = ty + this.nb[j][1];
        tnv.iz = tz + this.nb[j][2];

        if (
          tnv.ix < this.pLength &&
          tnv.ix > -1 &&
          tnv.iy < this.pWidth &&
          tnv.iy > -1 &&
          tnv.iz < this.pHeight &&
          tnv.iz > -1
        ) {
          index = tnv.ix * pWH + this.pHeight * tnv.iy + tnv.iz;

          // @ts-ignore
          if (this.vpBits[index] & INOUT && !(this.vpBits[index] & ISDONE)) {
            boundPoint.set(tnv.ix, tnv.iy, tz + this.nb[j][2], bp);
            dx = tnv.ix - bp.ix;
            dy = tnv.iy - bp.iy;
            dz = tnv.iz - bp.iz;
            square = dx * dx + dy * dy + dz * dz;
            // @ts-ignore
            this.vpDistance[index] = square;
            // @ts-ignore
            this.vpBits[index] |= ISDONE;
            // @ts-ignore
            this.vpBits[index] |= ISBOUND;

            outarray.push({
              ix: tnv.ix,
              iy: tnv.iy,
              iz: tnv.iz,
            });
          // @ts-ignore
          } else if (this.vpBits[index] & INOUT && this.vpBits[index] & ISDONE) {
            dx = tnv.ix - bp.ix;
            dy = tnv.iy - bp.iy;
            dz = tnv.iz - bp.iz;
            square = dx * dx + dy * dy + dz * dz;
            // @ts-ignore
            if (square < this.vpDistance[index]) {
              boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

              // @ts-ignore
              this.vpDistance[index] = square;
              // @ts-ignore
              if (!(this.vpBits[index] & ISBOUND)) {
                // @ts-ignore
                this.vpBits[index] |= ISBOUND;
                outarray.push({
                  ix: tnv.ix,
                  iy: tnv.iy,
                  iz: tnv.iz,
                });
              }
            }
          }
        }
      }
    }

    // console.log("part1", positout);
    for (i = 0, n = inarray.length; i < n; i++) {
      tx = inarray[i].ix;
      ty = inarray[i].iy;
      tz = inarray[i].iz;
      bp = boundPoint.get(tx, ty, tz);

      for (j = 6; j < 18; j++) {
        tnv.ix = tx + this.nb[j][0];
        tnv.iy = ty + this.nb[j][1];
        tnv.iz = tz + this.nb[j][2];

        if (
          tnv.ix < this.pLength &&
          tnv.ix > -1 &&
          tnv.iy < this.pWidth &&
          tnv.iy > -1 &&
          tnv.iz < this.pHeight &&
          tnv.iz > -1
        ) {
          index = tnv.ix * pWH + this.pHeight * tnv.iy + tnv.iz;

          // @ts-ignore
          if (this.vpBits[index] & INOUT && !(this.vpBits[index] & ISDONE)) {
            boundPoint.set(tnv.ix, tnv.iy, tz + this.nb[j][2], bp);

            dx = tnv.ix - bp.ix;
            dy = tnv.iy - bp.iy;
            dz = tnv.iz - bp.iz;
            square = dx * dx + dy * dy + dz * dz;
            // @ts-ignore
            this.vpDistance[index] = square;
            // @ts-ignore
            this.vpBits[index] |= ISDONE;
            // @ts-ignore
            this.vpBits[index] |= ISBOUND;

            outarray.push({
              ix: tnv.ix,
              iy: tnv.iy,
              iz: tnv.iz,
            });
          // @ts-ignore
          } else if (this.vpBits[index] & INOUT && this.vpBits[index] & ISDONE) {
            dx = tnv.ix - bp.ix;
            dy = tnv.iy - bp.iy;
            dz = tnv.iz - bp.iz;
            square = dx * dx + dy * dy + dz * dz;
            // @ts-ignore
            if (square < this.vpDistance[index]) {
              boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
              // @ts-ignore
              this.vpDistance[index] = square;
              // @ts-ignore
              if (!(this.vpBits[index] & ISBOUND)) {
                // @ts-ignore
                this.vpBits[index] |= ISBOUND;
                outarray.push({
                  ix: tnv.ix,
                  iy: tnv.iy,
                  iz: tnv.iz,
                });
              }
            }
          }
        }
      }
    }

    // console.log("part2", positout);
    for (i = 0, n = inarray.length; i < n; i++) {
      tx = inarray[i].ix;
      ty = inarray[i].iy;
      tz = inarray[i].iz;
      bp = boundPoint.get(tx, ty, tz);

      for (j = 18; j < 26; j++) {
        tnv.ix = tx + this.nb[j][0];
        tnv.iy = ty + this.nb[j][1];
        tnv.iz = tz + this.nb[j][2];

        if (
          tnv.ix < this.pLength &&
          tnv.ix > -1 &&
          tnv.iy < this.pWidth &&
          tnv.iy > -1 &&
          tnv.iz < this.pHeight &&
          tnv.iz > -1
        ) {
          index = tnv.ix * pWH + this.pHeight * tnv.iy + tnv.iz;

          // @ts-ignore
          if (this.vpBits[index] & INOUT && !(this.vpBits[index] & ISDONE)) {
            boundPoint.set(tnv.ix, tnv.iy, tz + this.nb[j][2], bp);

            dx = tnv.ix - bp.ix;
            dy = tnv.iy - bp.iy;
            dz = tnv.iz - bp.iz;
            square = dx * dx + dy * dy + dz * dz;
            // @ts-ignore
            this.vpDistance[index] = square;
            // @ts-ignore
            this.vpBits[index] |= ISDONE;
            // @ts-ignore
            this.vpBits[index] |= ISBOUND;

            outarray.push({
              ix: tnv.ix,
              iy: tnv.iy,
              iz: tnv.iz,
            });
          // @ts-ignore
          } else if (this.vpBits[index] & INOUT && this.vpBits[index] & ISDONE) {
            dx = tnv.ix - bp.ix;
            dy = tnv.iy - bp.iy;
            dz = tnv.iz - bp.iz;
            square = dx * dx + dy * dy + dz * dz;
            // @ts-ignore
            if (square < this.vpDistance[index]) {
              boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

              // @ts-ignore
              this.vpDistance[index] = square;
              // @ts-ignore
              if (!(this.vpBits[index] & ISBOUND)) {
                // @ts-ignore
                this.vpBits[index] |= ISBOUND;
                outarray.push({
                  ix: tnv.ix,
                  iy: tnv.iy,
                  iz: tnv.iz,
                });
              }
            }
          }
        }
      }
    }

    // console.log("part3", positout);
    return outarray;
  };

  marchingcubeinit(stype) {
    // @ts-ignore
    for (let i = 0, lim = this.vpBits.length; i < lim; i++) {
      if (stype === 1) {
        // vdw
        // @ts-ignore
        this.vpBits[i] &= ~ISBOUND;
      } else if (stype === 4) {
        // ses
        // @ts-ignore
        this.vpBits[i] &= ~ISDONE;
        // @ts-ignore
        if (this.vpBits[i] & ISBOUND) this.vpBits[i] |= ISDONE;
        // @ts-ignore
        this.vpBits[i] &= ~ISBOUND;
      } else if (stype === 2) {
        // after vdw
        // @ts-ignore
        if (this.vpBits[i] & ISBOUND && this.vpBits[i] & ISDONE) this.vpBits[i] &= ~ISBOUND;
        // @ts-ignore
        else if (this.vpBits[i] & ISBOUND && !(this.vpBits[i] & ISDONE)) this.vpBits[i] |= ISDONE;
      } else if (stype === 3) {
        // sas
        // @ts-ignore
        this.vpBits[i] &= ~ISBOUND;
      }
    }
  };

  marchingcube(stype) {
    this.marchingcubeinit(stype);
    this.verts = [];
    this.faces = [];
    // @ts-ignore
    MarchingCube.march(this.vpBits, this.verts, this.faces, {
      smooth: 1,
      nX: this.pLength,
      nY: this.pWidth,
      nZ: this.pHeight,
    });

    const pWH = this.pWidth * this.pHeight;
    for (let i = 0, vlen = this.verts.length; i < vlen; i++) {
      // @ts-ignore
      this.verts[i].atomid = this.vpAtomID[this.verts[i].x * pWH + this.pHeight * this.verts[i].y + this.verts[i].z];
    }

    // @ts-ignore
    MarchingCube.laplacianSmooth(1, this.verts, this.faces);
  };
}
