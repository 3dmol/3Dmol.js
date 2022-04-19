/* eslint-disable max-classes-per-file */
/* eslint-disable no-cond-assign */
import pako from 'pako';
import {CC} from './colors';
import GLShape from './GLShape';
import base64ToArray from './util/base64ToArray';
import {Object3D} from './WebGL/core';
import {Texture, VolumetricMaterial} from './WebGL/materials';
import {Matrix4, Vector3} from './WebGL/math';
import {Mesh} from './WebGL/objects';
import {Sphere} from './WebGL/shapes';
import Parsers from "./Parsers"

/**
 * VolumeData stores volumetric data. This includes file parsing
 * functionality.
 *
 * @class
 * @param {string} str - volumetric data
 * @param {string} format - format of supplied data (cube, dx, vasp); append .gz if compressed
 * @param {Object} options - normalize (zero mean, unit variance), negate
 */
export class VolumeData {
  /**
   * @param {string | Uint8Array} str
   * @param {string} format
   * @param {{ negate: any; normalize: any; }} [options]
   */
  constructor(str, format, options) {
    this.unit = {
      x: 1,
      y: 1,
      z: 1,
    }; // scale of each voxel
    this.origin = {
      x: 0,
      y: 0,
      z: 0,
    }; // origin (bottom "left", not center)
    this.size = {
      x: 0,
      y: 0,
      z: 0,
    }; // number of voxels in each direction
    this.data = new Float32Array([]); // actual floating point data, arranged

    // x->y->z
    this.matrix = null; // if set must transform data
    format = format.toLowerCase();

    if (/\.gz$/.test(format)) {
      // unzip gzipped files
      format = format.replace(/\.gz$/, '');
      try {
        if (this[format] && this[format].isbinary) {
          if (typeof str == 'string') {
            // assume base64 encoded
            str = base64ToArray(str);
          }
          str = pako.inflate(str);
        } else {
          str = new TextDecoder('utf-8').decode(pako.inflate(str));
        }
      } catch (err) {
        console.log(err);
      }
    }

    if (this[format]) {
      if (this[format].isbinary && typeof str == 'string') {
        str = base64ToArray(str);
      }
      this[format](str);
    }

    if (options) {
      if (options.negate) {
        for (let i = 0, n = this.data.length; i < n; i++) {
          this.data[i] = -this.data[i];
        }
      }
      if (options.normalize) {
        let total = 0.0;
        for (let i = 0, n = this.data.length; i < n; i++) {
          total += this.data[i];
        }
        const mean = total / this.data.length;
        console.log(`computed mean: ${mean}`);
        total = 0;
        for (let i = 0, n = this.data.length; i < n; i++) {
          const diff = this.data[i] - mean;
          total += diff * diff; // variance is ave of squared difference with mean
        }
        const variance = total / this.data.length;
        // console.log("Computed variance: "+variance);
        // now normalize
        for (let i = 0, n = this.data.length; i < n; i++) {
          this.data[i] = (this.data[i] - mean) / variance;
        }
      }
    }
  }

  /**
   * @function VolumeData.getIndex
   * @param {number} x ,y,z - the coordinates
   * @returns - index into flat array closest to provided coordinate; -1 if invalid
   * @param {number} y
   * @param {number} z
   */
  getIndex(x, y, z) {
    if (this.matrix) {
      // all transformation is done through matrix multiply
      if (!this.inversematrix) {
        this.inversematrix = new Matrix4().getInverse(this.matrix);
      }
      let pt = new Vector3(x, y, z);
      pt = pt.applyMatrix4(this.inversematrix);
      x = pt.x;
      y = pt.y;
      z = pt.z;
    } else {
      // use simple origin/unit transform
      x -= this.origin.x;
      y -= this.origin.y;
      z -= this.origin.z;

      x /= this.unit.x;
      y /= this.unit.y;
      z /= this.unit.z;
    }
    x = Math.round(x);
    y = Math.round(y);
    z = Math.round(z);

    if (x < 0 || x >= this.size.x) return -1;
    if (y < 0 || y >= this.size.y) return -1;
    if (z < 0 || z >= this.size.z) return -1;

    return x * this.size.y * this.size.z + y * this.size.z + z;
  }

  /**
   * @function VolumeData.getVal
   * @param {number} x ,y,z - the coordinates
   * @returns - value closest to provided coordinate; zero if coordinate invalid
   * @param {any} y
   * @param {any} z
   */
  getVal(x, y, z) {
    const i = this.getIndex(x, y, z);
    if (i < 0) return 0;
    return this.data[i];
  }

  /**
   * @param {number} index
   */
  getCoordinates(index) {
    let x = index / (this.size.y * this.size.z);
    let y = index % (this.size.y * this.size.z);
    let z = index % this.size.z;

    x *= this.unit.x;
    y *= this.unit.y;
    z *= this.unit.z;

    x += this.origin.x;
    y += this.origin.y;
    z += this.origin.z;

    return {x, y, z};
  }
  
  /*
   * parse vasp data
   * Essentially this parser converts the CHGCAR data into
   * cube data. It has been adapted from 'chg2cube.pl' found in
   * http://theory.cm.utexas.edu/vtsttools/
   */
  /**
   * @param {string} str
   */
  vasp(str) {
    const lines = str.replace(/^\s+/, '').split(/[\n\r]/);

    const atomicData = Parsers.vasp(str)[0];
    const natoms = atomicData.length;

    if (natoms === 0) {
      console.log(
        'No good formating of CHG or CHGCAR file, not atomic information provided in the file.'
      );
      // @ts-ignore
      this.data = [];
      return;
    }

    // Assume atomic units
    //    var unittype = "bohr/hartree";
    const lUnits = 1.889725992;
    const eUnits = 0.036749309;

    // copied from Parsers.vasp
    const convFactor = parseFloat(lines[1]);
    // This is how Vasp reads in the basis We need the lUnits in order to
    // compute the volume of the cell. Afterwards to obtain the axis for the
    // voxels we have to remove this unit and divide by the number of voxels in
    // each dimension
    let v;
    v = lines[2].replace(/^\s+/, '').split(/\s+/);
    let xVec = new Vector3(parseFloat(v[0]), parseFloat(v[1]), parseFloat(v[2])).multiplyScalar(
      convFactor * lUnits
    );
    v = lines[3].replace(/^\s+/, '').split(/\s+/);
    let yVec = new Vector3(parseFloat(v[0]), parseFloat(v[1]), parseFloat(v[2])).multiplyScalar(
      convFactor * lUnits
    );
    v = lines[4].replace(/^\s+/, '').split(/\s+/);
    let zVec = new Vector3(parseFloat(v[0]), parseFloat(v[1]), parseFloat(v[2])).multiplyScalar(
      convFactor * lUnits
    );

    // correct volume for non-orthognal box (expansion by minors)
    let vol =
      xVec.x * (yVec.y * zVec.z - zVec.y * yVec.z) -
      yVec.x * (xVec.y * zVec.z - zVec.y * xVec.z) +
      zVec.x * (xVec.y * yVec.z - yVec.y * xVec.z);

    vol = Math.abs(vol) / lUnits ** 3;
    const volScale = 1.0 / vol; // This Only for CHGCAR files

    // We splice the structure information
    // 2 (header) + 3 (vectors) + 2 (atoms) + 1 (vaspMode) + natoms (coords) + 1 (blank line)
    lines.splice(0, 2 + 3 + 2 + 1 + natoms + 1);

    const lineArr = lines[0].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');

    // @ts-ignore
    const nX = Math.abs(lineArr[0]);
    // @ts-ignore
    const nY = Math.abs(lineArr[1]);
    // @ts-ignore
    const nZ = Math.abs(lineArr[2]);

    const origin = new Vector3(0, 0, 0);
    this.origin = new Vector3(0, 0, 0);

    this.size = {x: nX, y: nY, z: nZ};
    this.unit = new Vector3(xVec.x, yVec.y, zVec.z);

    // resize the vectors accordingly
    xVec = xVec.multiplyScalar(1 / (lUnits * nX));
    yVec = yVec.multiplyScalar(1 / (lUnits * nY));
    zVec = zVec.multiplyScalar(1 / (lUnits * nZ));

    if (
      xVec.y !== 0 ||
      xVec.z !== 0 ||
      yVec.x !== 0 ||
      yVec.z !== 0 ||
      zVec.x !== 0 ||
      zVec.y !== 0
    ) {
      // need a transformation matrix
      this.matrix = new Matrix4(
        xVec.x,
        yVec.x,
        zVec.x,
        0,
        xVec.y,
        yVec.y,
        zVec.y,
        0,
        xVec.z,
        yVec.z,
        zVec.z,
        0,
        0,
        0,
        0,
        1
      );
      // include translation in matrix
      this.matrix = this.matrix.multiplyMatrices(
        this.matrix,
        new Matrix4().makeTranslation(origin.x, origin.y, origin.z)
      );
      // all translation and scaling done by matrix, so reset origin and unit
      this.origin = new Vector3(0, 0, 0);
      this.unit = new Vector3(1, 1, 1);
    }

    lines.splice(0, 1); // Remove the dimension line
    let raw = lines.join(' ');

    raw = raw.replace(/^\s+/, '');
    // @ts-ignore
    raw = raw.split(/[\s\r]+/);
    // @ts-ignore
    raw.splice(nX * nY * nZ + 1);

    // @ts-ignore
    const preConvertedData = new Float32Array(raw); // We still have to format it to get the density

    for (let i = 0; i < preConvertedData.length; i++) {
      preConvertedData[i] = preConvertedData[i] * volScale * eUnits;
    }

    this.data = preConvertedData;

    // console.log(xVec);
    // console.log(yVec);
    // console.log(zVec);
    // console.log(this.unit);
    // console.log(this.origin);
    // console.log(this.matrix);
    // console.log(this.data);
  }

  // parse dx data - does not support all features of the file format
  /**
   * @param {string} str
   */
  dx(str) {
    const lines = str.split(/[\n\r]+/);
    let i;
    let m;
    const recounts = /gridpositions\s+counts\s+(\d+)\s+(\d+)\s+(\d+)/;
    const reorig = /^origin\s+(\S+)\s+(\S+)\s+(\S+)/;
    const redelta = /^delta\s+(\S+)\s+(\S+)\s+(\S+)/;
    const follows = /data follows/;

    for (i = 0; i < lines.length; i++) {
      let line = lines[i];
      if ((m = recounts.exec(line))) {
        const nX = parseInt(m[1]);
        const nY = parseInt(m[2]);
        const nZ = parseInt(m[3]);
        this.size = {x: nX, y: nY, z: nZ};
      } else if ((m = redelta.exec(line))) {
        const xunit = parseFloat(m[1]);
        if (parseFloat(m[2]) !== 0 || parseFloat(m[3]) !== 0) {
          console.log('Non-orthogonal delta matrix not currently supported in dx format');
        }
        i += 1;
        line = lines[i];
        m = redelta.exec(line);
        if (m == null) {
          console.log('Parse error in dx delta matrix');
          return;
        }

        const yunit = parseFloat(m[2]);
        if (parseFloat(m[1]) !== 0 || parseFloat(m[3]) !== 0) {
          console.log('Non-orthogonal delta matrix not currently supported in dx format');
        }

        i += 1;
        line = lines[i];
        m = redelta.exec(line);
        if (m == null) {
          console.log('Parse error in dx delta matrix');
          return;
        }

        const zunit = parseFloat(m[3]);
        if (parseFloat(m[1]) !== 0 || parseFloat(m[2]) !== 0) {
          console.log('Non-orthogonal delta matrix not currently supported in dx format');
        }
        this.unit = new Vector3(xunit, yunit, zunit);
      } else if ((m = reorig.exec(line))) {
        const xorig = parseFloat(m[1]);
        const yorig = parseFloat(m[2]);
        const zorig = parseFloat(m[3]);
        this.origin = new Vector3(xorig, yorig, zorig);
      } else if ((m = follows.exec(line))) {
        break;
      }
    }
    i += 1;
    if (!this.size || !this.origin || !this.unit || !this.size) {
      console.log('Error parsing dx format');
      return;
    }
    let raw = lines.splice(i).join(' ');
    // @ts-ignore
    raw = raw.split(/[\s\r]+/);
    // @ts-ignore
    this.data = new Float32Array(raw);
  }

  // parse cube data
  /**
   * @param {string} str
   */
  cube(str) {
    const lines = str.split(/\r?\n/);

    if (lines.length < 6) return;

    const parsed = (/** @type {{modelData: { cryst: any}[]}} */(Parsers.cube(str)))
    const {cryst} = parsed.modelData[0];

    const lineArr = lines[2].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');

    const atomsnum = parseFloat(lineArr[0]); // includes sign, which indicates presence of oribital line in header
    const natoms = Math.abs(atomsnum);

    this.origin = cryst.origin;
    this.size = cryst.size;
    this.unit = cryst.unit;
    this.matrix = cryst.matrix4;

    let headerlines = 6;
    if (atomsnum < 0) headerlines += 1; // see: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
    let raw = lines.splice(natoms + headerlines).join(' ');
    raw = raw.replace(/^\s+/, '');
    // @ts-ignore
    raw = raw.split(/[\s\r]+/);
    // @ts-ignore
    this.data = new Float32Array(raw);
  }

  // parse cp4 files
  /**
   * @param {Iterable<number>} binData
   */
  ccp4(binData) {
    // http://www.ccp4.ac.uk/html/maplib.html#description
    // code from ngl: https://github.com/arose/ngl/blob/master/js/ngl/parser.js
    const header = {};
    const bin = new Int8Array(binData);
    const intView = new Int32Array(bin.buffer, 0, 56);
    const floatView = new Float32Array(bin.buffer, 0, 56);
    const dv = new DataView(bin.buffer);

    // 53  MAP         Character string 'MAP ' to identify file type
    header.MAP = String.fromCharCode(
      dv.getUint8(52 * 4),
      dv.getUint8(52 * 4 + 1),
      dv.getUint8(52 * 4 + 2),
      dv.getUint8(52 * 4 + 3)
    );

    // 54  MACHST      Machine stamp indicating machine type which wrote file
    //                 17 and 17 for big-endian or 68 and 65 for little-endian
    header.MACHST = [dv.getUint8(53 * 4), dv.getUint8(53 * 4 + 1)];

    // swap byte order when big endian
    if (header.MACHST[0] === 17 && header.MACHST[1] === 17) {
      const n = bin.byteLength;
      for (let i = 0; i < n; i += 4) {
        dv.setFloat32(i, dv.getFloat32(i), true);
      }
    }

    header.NX = intView[0]; // NC - columns (fastest changing)
    header.NY = intView[1]; // NR - rows
    header.NZ = intView[2]; // NS - sections (slowest changing)

    // mode
    //  0 image : signed 8-bit bytes range -128 to 127
    //  1 image : 16-bit halfwords
    //  2 image : 32-bit reals
    //  3 transform : complex 16-bit integers
    //  4 transform : complex 32-bit reals
    //  6 image : unsigned 16-bit range 0 to 65535
    // 16 image: unsigned char * 3 (for rgb data, non-standard)
    //
    // Note: Mode 2 is the normal mode used in the CCP4 programs.
    //       Other modes than 2 and 0 may NOT WORK
    header.MODE = intView[3];

    // start
    header.NXSTART = intView[4]; // NCSTART - first column
    header.NYSTART = intView[5]; // NRSTART - first row
    header.NZSTART = intView[6]; // NSSTART - first section

    // intervals
    header.MX = intView[7]; // intervals along x
    header.MY = intView[8]; // intervals along y
    header.MZ = intView[9]; // intervals along z

    // cell length (Angstroms in CCP4)
    header.xlen = floatView[10];
    header.ylen = floatView[11];
    header.zlen = floatView[12];

    // cell angle (Degrees)
    header.alpha = floatView[13];
    header.beta = floatView[14];
    header.gamma = floatView[15];

    // axis correspondence (1,2,3 for X,Y,Z)
    header.MAPC = intView[16]; // column
    header.MAPR = intView[17]; // row
    header.MAPS = intView[18]; // section

    // density statistics
    header.DMIN = floatView[19];
    header.DMAX = floatView[20];
    header.DMEAN = floatView[21];

    // space group number 0 or 1 (default=0)
    header.ISPG = intView[22];

    // number of bytes used for symmetry data (0 or 80)
    header.NSYMBT = intView[23];

    // Flag for skew transformation, =0 none, =1 if foll
    header.LSKFLG = intView[24];

    // 26-34  SKWMAT  Skew matrix S (in order S11, S12, S13, S21 etc) if
    //                LSKFLG .ne. 0.
    // 35-37  SKWTRN  Skew translation t if LSKFLG != 0.
    //                Skew transformation is from standard orthogonal
    //                coordinate frame (as used for atoms) to orthogonal
    //                map frame, as Xo(map) = S * (Xo(atoms) - t)
    // 38      future use       (some of these are used by the MSUBSX routines
    //  .          "              in MAPBRICK, MAPCONT and FRODO)
    //  .          "   (all set to zero by default)
    //  .          "
    // 52          "
    // 50-52 origin in X,Y,Z used for transforms
    header.originX = floatView[49];
    header.originY = floatView[50];
    header.originZ = floatView[51];

    // 53  MAP         Character string 'MAP ' to identify file type
    // => see top of this parser
    // 54  MACHST      Machine stamp indicating machine type which wrote file
    // => see top of this parser
    // Rms deviation of map from mean density
    header.ARMS = floatView[54];

    // 56      NLABL           Number of labels being used
    // 57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)
    // console.log("Map has min,mean,average,rmsddv: "+header.DMIN+","+header.DMAX+","+header.DMEAN+","+header.ARMS);
    // create transformation matrix, code mostly copied from ngl
    const h = header;
    const basisX = [h.xlen, 0, 0];

    const basisY = [
      h.ylen * Math.cos((Math.PI / 180.0) * h.gamma),
      h.ylen * Math.sin((Math.PI / 180.0) * h.gamma),
      0,
    ];

    const basisZ = [
      h.zlen * Math.cos((Math.PI / 180.0) * h.beta),
      (h.zlen *
        (Math.cos((Math.PI / 180.0) * h.alpha) -
          Math.cos((Math.PI / 180.0) * h.gamma) * Math.cos((Math.PI / 180.0) * h.beta))) /
        Math.sin((Math.PI / 180.0) * h.gamma),
      0,
    ];
    basisZ[2] = Math.sqrt(
      h.zlen *
        h.zlen *
        Math.sin((Math.PI / 180.0) * h.beta) *
        Math.sin((Math.PI / 180.0) * h.beta) -
        basisZ[1] * basisZ[1]
    );

    const basis = [0, basisX, basisY, basisZ];
    const nxyz = [0, h.MX, h.MY, h.MZ];
    const mapcrs = [0, h.MAPC, h.MAPR, h.MAPS];

    this.matrix = new Matrix4();

    this.matrix.set(
      basis[mapcrs[1]][0] / nxyz[mapcrs[1]],
      basis[mapcrs[2]][0] / nxyz[mapcrs[2]],
      basis[mapcrs[3]][0] / nxyz[mapcrs[3]],
      0,

      basis[mapcrs[1]][1] / nxyz[mapcrs[1]],
      basis[mapcrs[2]][1] / nxyz[mapcrs[2]],
      basis[mapcrs[3]][1] / nxyz[mapcrs[3]],
      0,

      basis[mapcrs[1]][2] / nxyz[mapcrs[1]],
      basis[mapcrs[2]][2] / nxyz[mapcrs[2]],
      basis[mapcrs[3]][2] / nxyz[mapcrs[3]],
      0,

      0,
      0,
      0,
      1
    );
    // include translation in matrix, NXSTART etc are an offset in grid space
    this.matrix = this.matrix.multiplyMatrices(
      this.matrix,
      new Matrix4().makeTranslation(
        h.NXSTART + h.originX,
        h.NYSTART + h.originY,
        h.NZSTART + h.originZ
      )
    );
    // all translation and scaling done by matrix, so reset origin and unit
    this.origin = new Vector3(0, 0, 0);
    this.unit = new Vector3(1, 1, 1);
    this.size = {x: header.NX, y: header.NY, z: header.NZ};
    this.dimensionorder = [header.MAPC, header.MAPR, header.MAPS];
    const data = new Float32Array(bin.buffer, 1024 + header.NSYMBT);
    // data must by (slowest changing) x,y,z (fastest changing)
    const {NX} = header;
    const {NY} = header;
    const {NZ} = header;
    this.data = new Float32Array(NX * NY * NZ);
    for (let i = 0; i < NX; i++) {
      for (let j = 0; j < NY; j++) {
        for (let k = 0; k < NZ; k++) {
          // should I be concerned that I'm not using mapc?
          this.data[(i * NY + j) * NZ + k] = data[(k * NY + j) * NX + i];
        }
      }
    }
  }
}

// @ts-ignore
VolumeData.prototype.ccp4.isbinary = true;

// interpolation function used from http://hevi.info/do-it-yourself/interpolating-and-array-to-fit-another-size/
/**
 * @param {string | any[]} data
 * @param {number} fitCount
 */
function interpolateArray(data, fitCount) {
  /**
   * @param {number} before
   * @param {number} after
   * @param {number} atPoint
   */
  function linearInterpolate(before, after, atPoint) {
    return before + (after - before) * atPoint;
  }
  const newData = [];
  const springFactor = (data.length - 1) / (fitCount - 1);
  newData[0] = data[0]; // for new allocation
  for (let i = 1; i < fitCount - 1; i++) {
    const tmp = i * springFactor;
    const before = Math.floor(tmp).toFixed();
    const after = Math.ceil(tmp).toFixed();
    // @ts-ignore
    const atPoint = tmp - before;
    newData[i] = linearInterpolate(data[before], data[after], atPoint);
  }
  newData[fitCount - 1] = data[data.length - 1]; // for new allocation
  return newData;
}

/**
 * A GLVolumetricRender is a "shape" for representing volumetric data as a density distribution.
 */
export class GLVolumetricRender {
  /** @type {Object3D | null} */
  shapeObj = null;
  /** @type {Object3D | null} */
  renderedShapeObj = null;
  // need to create transformation matrix that maps model points into
  // texture space
  // need extent (bounding box dimensions), maxdepth (box diagonal),
  // texmatrix (conversion from model to texture coords), minunit,
  texmatrix;
  extent;
  maxdepth;
  minunit;

  hidden = false;
  /**
   * @type {any}
   */
  data

  /** @type {number|undefined} */
  shapePosition;

  /**
   * @param {VolumeData} data
   * @param {import('./specs').VolumetricRendererSpec} spec
   */
  constructor(data, spec) {
    spec = spec || {};
    const transferfn = spec.transferfn || [];

    this.subsamples = spec.subsamples || 5.0;

    const TRANSFER_BUFFER_SIZE = 256;
    const transferfunctionbuffer = [];
    // arrange points based on position property
    transferfn.forEach((/** @type {{ value: string | number; }} */ a) => {
      // @ts-ignore
      a.value = parseFloat(a.value);
    });
    transferfn.sort((/** @type {{ value: number; }} */ a, /** @type {{ value: number; }} */ b) => a.value - b.value);
    this.min = transferfn[0].value;
    if (transferfn.length === 0) transferfn.push(transferfn[0]); // need at least two
    this.max = transferfn[transferfn.length - 1].value;

    // create and fill an array of interpolated values per 2 colors
    let pos1;
    let pos2;
    let color1;
    let color2;
    let R;
    let G;
    let B;
    let A;
    let alpha1;
    let alpha2;
    for (let i = 0; i < transferfn.length - 1; i++) {
      color1 = CC.color(transferfn[i].color);
      color2 = CC.color(transferfn[i + 1].color);
      alpha1 = transferfn[i].opacity;
      alpha2 = transferfn[i + 1].opacity;
      pos1 = Math.floor(((transferfn[i].value - this.min) * TRANSFER_BUFFER_SIZE) / (this.max - this.min));
      pos2 = Math.floor(((transferfn[i + 1].value - this.min) * TRANSFER_BUFFER_SIZE) / (this.max - this.min));
      if (pos1 === pos2) continue;
      R = interpolateArray([color1.r * 255, color2.r * 255], pos2 - pos1);
      G = interpolateArray([color1.g * 255, color2.g * 255], pos2 - pos1);
      B = interpolateArray([color1.b * 255, color2.b * 255], pos2 - pos1);
      A = interpolateArray([alpha1 * 255, alpha2 * 255], pos2 - pos1);

      for (let j = 0; j < R.length; j++) {
        transferfunctionbuffer.push(R[j]);
        transferfunctionbuffer.push(G[j]);
        transferfunctionbuffer.push(B[j]);
        transferfunctionbuffer.push(A[j]); // opacity will be added later
      }
    }

    this.transferfunctionbuffer = new Uint8ClampedArray(transferfunctionbuffer);

    // possibly non-orthnormal basis if matrix
    if (data.matrix) {
      // figure out bounding box of transformed grid
      const start = new Vector3(0, 0, 0);
      const end = new Vector3(data.size.x, data.size.y, data.size.z);
      const unit = new Vector3(1, 1, 1);

      start.applyMatrix4(data.matrix);
      end.applyMatrix4(data.matrix);
      unit.applyMatrix4(data.matrix).sub(start);

      this.extent = [
        [start.x, start.y, start.z],
        [end.x, end.y, end.z],
      ];

      // check all corners, these may not be the farthest apart
      for (let i = 1; i < 7; i++) {
        end.x = i & 1 ? data.size.x : 0;
        end.y = i & 2 ? data.size.y : 0;
        end.z = i & 4 ? data.size.z : 0;
        end.applyMatrix4(data.matrix);
        this.extent[0][0] = Math.min(this.extent[0][0], end.x);
        this.extent[0][1] = Math.min(this.extent[0][1], end.y);
        this.extent[0][2] = Math.min(this.extent[0][2], end.z);
        this.extent[1][0] = Math.max(this.extent[1][0], end.x);
        this.extent[1][1] = Math.max(this.extent[1][1], end.y);
        this.extent[1][2] = Math.max(this.extent[1][2], end.z);
      }

      const xoff = end.x - start.x;
      const yoff = end.y - start.y;
      const zoff = end.z - start.z;
      this.maxdepth = Math.sqrt(xoff * xoff + yoff * yoff + zoff * zoff);

      this.minunit = Math.min(Math.min(unit.x, unit.y), unit.z);

      // invert onto grid, then scale by grid dimensions to get
      // normalized texture coordinates
      this.texmatrix = new Matrix4()
        .identity()
        .scale({x: data.size.x, y: data.size.y, z: data.size.z});
      this.texmatrix = this.texmatrix.multiplyMatrices(data.matrix, this.texmatrix);

      this.texmatrix = this.texmatrix.getInverse(this.texmatrix);
    } else {
      this.texmatrix = new Matrix4().identity();
      const xoff = data.unit.x * data.size.x;
      const yoff = data.unit.y * data.size.y;
      const zoff = data.unit.z * data.size.z;
      // scale doesn't apply to the translation vector, so preapply it
      this.texmatrix.makeTranslation(
        -data.origin.x / xoff,
        -data.origin.y / yoff,
        -data.origin.z / zoff
      );
      this.texmatrix.scale({x: 1.0 / xoff, y: 1.0 / yoff, z: 1.0 / zoff});
      this.minunit = Math.min(Math.min(data.unit.x, data.unit.y), data.unit.z);

      // need the bounding box so we can intersect with rays
      this.extent = [
        [data.origin.x, data.origin.y, data.origin.z],
        [data.origin.x + xoff, data.origin.y + yoff, data.origin.z + zoff],
      ];

      this.maxdepth = Math.sqrt(xoff * xoff + yoff * yoff + zoff * zoff);
    }

    // use GLShape to construct box
    const shape = new GLShape();
    shape.addBox({
      corner: {x: this.extent[0][0], y: this.extent[0][1], z: this.extent[0][2]},
      dimensions: {
        w: this.extent[1][0] - this.extent[0][0],
        h: this.extent[1][1] - this.extent[0][1],
        d: this.extent[1][2] - this.extent[0][2],
      },
    });

    this.geo = shape.finalize();
    this.boundingSphere = new Sphere();
    this.boundingSphere.center = {
      x: (this.extent[0][0] + this.extent[1][0]) / 2.0,
      y: (this.extent[0][1] + this.extent[1][1]) / 2.0,
      z: (this.extent[0][2] + this.extent[1][2]) / 2.0,
    };
    this.boundingSphere.radius = this.maxdepth / 2;

    // volume selectivity based on given coords and distance
    if (spec.coords !== undefined && spec.seldist !== undefined) {
      const mask = new Uint8Array(data.data.length);
      // for each coordinate
      const d = spec.seldist;
      const d2 = d * d;
      for (let i = 0, n = spec.coords.length; i < n; i++) {
        const c = spec.coords[i];
        const minx = c.x - d;
        const miny = c.y - d;
        const minz = c.z - d;
        const maxx = c.x + d;
        const maxy = c.y + d;
        const maxz = c.z + d;
        if (data.getIndex(minx, miny, minz) >= 0 || data.getIndex(maxx, maxy, maxz) >= 0) {
          // bounding box overlaps grid
          // iterate over the grid points in the seldist bounding box
          // minunit may be inefficient if axes have very different units. oh well.
          for (let x = minx; x < maxx; x += this.minunit) {
            for (let y = miny; y < maxy; y += this.minunit) {
              for (let z = minz; z < maxz; z += this.minunit) {
                const idx = data.getIndex(x, y, z);
                if (idx >= 0 && !mask[idx]) {
                  // if not already masked, check distance
                  const distsq =
                    (x - c.x) * (x - c.x) + (y - c.y) * (y - c.y) + (z - c.z) * (z - c.z);
                  if (distsq < d2) {
                    mask[idx] = 1;
                  }
                }
              }
            }
          }
        }
      }
      // any place mask is zero, make infinite in data
      for (let i = 0, n = data.data.length; i < n; i++) {
        if (mask[i] === 0) data.data[i] = Infinity;
      }
    }
  }

  /**
   * Initialize webgl objects for rendering
   * @param {Object3D} group
   *
   */
  globj(group) {
    if (this.renderedShapeObj) {
      group.remove(this.renderedShapeObj);
      this.renderedShapeObj = null;
    }

    if (this.hidden) return;

    this.shapeObj = new Object3D();
    let material = null;

    const texture = new Texture(this.data, true);
    const transfertexture = new Texture(this.transferfunctionbuffer, false);
    texture.needsUpdate = true;
    transfertexture.needsUpdate = true;
    transfertexture.flipY = false;

    material = new VolumetricMaterial({
      transferfn: transfertexture,
      transfermin: this.min,
      transfermax: this.max,
      map: texture,
      extent: this.extent,
      maxdepth: this.maxdepth,
      texmatrix: this.texmatrix,
      unit: this.minunit,
      subsamples: this.subsamples,
    });

    const mesh = new Mesh(this.geo, material);
    this.shapeObj.add(mesh);

    this.renderedShapeObj = this.shapeObj.clone();
    group.add(this.renderedShapeObj);
  }

  /**
   * @param {{ remove: (arg0: Object3D) => void; }} group
   */
  removegl(group) {
    if (this.renderedShapeObj) {
      // dispose of geos and materials
      if (this.renderedShapeObj.geometry !== undefined) this.renderedShapeObj.geometry.dispose();
      if (this.renderedShapeObj.material !== undefined) this.renderedShapeObj.material.dispose();
      group.remove(this.renderedShapeObj);
      this.renderedShapeObj = null;
    }
    this.shapeObj = null;
  }

  get position() {
    return this.boundingSphere.center;
  }

  get x() {
    return this.boundingSphere.center.x;
  }

  get y() {
    return this.boundingSphere.center.y;
  }

  get z() {
    return this.boundingSphere.center.z;
  }
}
