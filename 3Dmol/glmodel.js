/* eslint-disable no-underscore-dangle */
/* eslint-disable eqeqeq */
// A model is a collection of related atoms.  Bonds are only allowed between
// atoms in the same model.  An atom is uniquely specified by its model id and
// its serial number.
// A glmodel knows how to apply the styles on each atom to create a gl object
import netcdfjs from 'netcdfjs';
import pako from 'pako';
import $ from 'jquery';
import {CC, getColorFromStyle, elementColors} from './colors';
import GLDraw from './GLDraw';
import Gradient from './Gradient';
import deepCopy from './util/deepCopy';
import extend from './util/extend';
import getbin from './util/getbin';
import getExtent from './util/getExtent';
import makeFunction from './util/makeFunction';
import {Object3D, Geometry} from './WebGL/core';
import Color from './WebGL/core/Color';
import Parsers from './parsers';
import {
  LineBasicMaterial,
  MeshLambertMaterial,
  SphereImposterMaterial,
  StickImposterMaterial,
  InstancedMaterial,
} from './WebGL/materials';
import {Vector3, Matrix3, conversionMatrix3, Matrix4} from './WebGL/math';
import {Mesh, LinePieces, Line} from './WebGL/objects';
import {Cylinder, Sphere} from './WebGL/shapes';
import getAtomProperty from './util/getAtomProperty';
import specStringToObject from './util/specStringToObject';
import getPropertyRange from './util/getPropertyRange';
import drawCartoon from './drawCartoon';
/**
 * GLModel represents a group of related atoms
 * @constructor
 * @param {number=} mid
 * @param {Object=} defaultcolors Object defining default atom colors as atom => color property value pairs
 * @see download
 */
// class variables go here
const defaultAtomStyle = {
  line: {},
};

const defaultlineWidth = 1.0;

// Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
// https://en.wikipedia.org/wiki/Van_der_Waals_radius
const vdwRadii = {
  H: 1.2,
  He: 1.4,
  Li: 1.82,
  Be: 1.53,
  B: 1.92,
  C: 1.7,
  N: 1.55,
  O: 1.52,
  F: 1.47,
  Ne: 1.54,
  Na: 2.27,
  Mg: 1.73,
  Al: 1.84,
  Si: 2.1,
  P: 1.8,
  S: 1.8,
  Cl: 1.75,
  Ar: 1.88,
  K: 2.75,
  Ca: 2.31,
  Ni: 1.63,
  Cu: 1.4,
  Zn: 1.39,
  Ga: 1.87,
  Ge: 2.11,
  As: 1.85,
  Se: 1.9,
  Br: 1.85,
  Kr: 2.02,
  Rb: 3.03,
  Sr: 2.49,
  Pd: 1.63,
  Ag: 1.72,
  Cd: 1.58,
  In: 1.93,
  Sn: 2.17,
  Sb: 2.06,
  Te: 2.06,
  I: 1.98,
  Xe: 2.16,
  Cs: 3.43,
  Ba: 2.68,
  Pt: 1.75,
  Au: 1.66,
  Hg: 1.55,
  Tl: 1.96,
  Pb: 2.02,
  Bi: 2.07,
  Po: 1.97,
  At: 2.02,
  Rn: 2.2,
  Fr: 3.48,
  Ra: 2.83,
  U: 1.86,
};

// prop : It is used to add the option for property in context menu in the 3dmol ui
// the code for prop can be found under /ui/ui.js -> UI -> ContextMenu -> setProperties -> submit.ui.on
// gui : It is used to generate forms for different features in the 3dmol ui
// the code for gui can be found under /ui/form.js -> Form (Form defination)
// floatType : separates integer from float since these are used in
// input validation of the 3dmol ui

// type is irrelivent here becuase htey are are invalid
const validExtras = {
  // valid atom specs are ok too
  model: {type: 'string', valid: false}, // a single model or list of models from which atoms should be selected
  bonds: {type: 'number', valid: false, gui: true}, // overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms
  predicate: {type: 'string', valid: false}, // user supplied function that gets passed an {AtomSpec} and should return true if the atom should be selected
  invert: {type: 'boolean', valid: false, gui: true}, // if set, inverts the meaning of the selection
  byres: {type: 'boolean', valid: false, gui: true}, // if set, expands the selection to include all atoms of any residue that has any atom selected
  expand: {type: 'number', valid: false, gui: false}, // expands the selection to include all atoms within a given distance from the selection
  within: {type: 'string', valid: false}, // intersects the selection with the set of atoms within a given distance from another selection
  and: {type: 'string', valid: false}, // and boolean logic
  or: {type: 'string', valid: false}, // or boolean logic
  not: {type: 'string', valid: false}, // not boolean logic
};

const validLineSpec = {
  hidden: {type: 'boolean', gui: true},
  linewidth: {type: 'number', floatType: true, gui: true, step: 0.1, default: defaultlineWidth},
  colorscheme: {type: 'colorscheme', gui: true},
  color: {type: 'color', gui: true},
  opacity: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0, max: 1},
};

const validCrossSpec = {
  hidden: {type: 'boolean', gui: true},
  linewidth: {
    type: 'number',
    floatType: true,
    gui: false,
    step: 0.1,
    default: defaultlineWidth,
    min: 0,
  }, // deprecated
  colorscheme: {type: 'colorscheme', gui: true},
  color: {type: 'color', gui: true},
  radius: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0.1},
  scale: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0},
  opacity: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0, max: 1},
};

const validStickSpec = {
  hidden: {type: 'boolean', gui: true},
  colorscheme: {type: 'colorscheme', gui: true},
  color: {type: 'color', gui: true},
  radius: {type: 'number', floatType: true, gui: true, step: 0.1, default: 0.25, min: 0.1},
  singleBonds: {type: 'boolean', gui: true},
  opacity: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0, max: 1},
};

const validSphereSpec = {
  hidden: {type: 'boolean', gui: false}, // needed in the new gui it has separate function to hide the spheres
  singleBonds: {type: 'boolean', gui: true},
  colorscheme: {type: 'colorscheme', gui: true},
  color: {type: 'color', gui: true},
  radius: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1.5, min: 0},
  scale: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1.0, min: 0.1},
  opacity: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0, max: 1},
};

const validCartoonSpec = {
  style: {validItems: ['trace', 'oval', 'rectangle', 'parabola', 'edged'], gui: true},
  color: {type: 'color', gui: true, spectrum: true},
  arrows: {type: 'boolean', gui: true},
  ribbon: {type: 'boolean', gui: true},
  hidden: {type: 'boolean', gui: true},
  tubes: {type: 'boolean', gui: true},
  thickness: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0},
  width: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0},
  opacity: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0, max: 1},
};

// class functions
// return true if a and b represent the same style
function sameObj(a, b) {
  if (a && b) return JSON.stringify(a) == JSON.stringify(b);
  return a == b;
}

/** @typedef SphereStyleSpec
 * @prop {boolean} hidden - do not show atom
 * @prop {number} radius - override van der waals radius
 * @prop {number} scale - scale radius by specified amount
 * @prop {import('./specs').ColorSchemeSpec} colorscheme - element based coloring
 * @prop {import('./specs').ColorSpec} color - fixed coloring, overrides colorscheme
 * @prop {number} opacity - opacity, must be the same for all atoms in the model
 */

function addLine(vertexArray, colorArray, offset, p1, p2, c1) {
  // make line from p1 to p2, does not incremeant counts
  vertexArray[offset] = p1.x;
  vertexArray[offset + 1] = p1.y;
  vertexArray[offset + 2] = p1.z;
  colorArray[offset] = c1.r;
  colorArray[offset + 1] = c1.g;
  colorArray[offset + 2] = c1.b;
  vertexArray[offset + 3] = p2.x;
  vertexArray[offset + 4] = p2.y;
  vertexArray[offset + 5] = p2.z;
  colorArray[offset + 3] = c1.r;
  colorArray[offset + 4] = c1.g;
  colorArray[offset + 5] = c1.b;
}

// bonds as cylinders
const defaultStickRadius = 0.25;

function drawSphereImposter(geo, center, radius, C) {
  // create flat square
  const geoGroup = geo.updateGeoGroup(4);
  let i;
  const startv = geoGroup.vertices;
  const start = startv * 3;
  const {vertexArray} = geoGroup;
  const {colorArray} = geoGroup;

  // use center point for each vertex
  for (i = 0; i < 4; i++) {
    vertexArray[start + 3 * i] = center.x;
    vertexArray[start + 3 * i + 1] = center.y;
    vertexArray[start + 3 * i + 2] = center.z;
  }

  // same colors for all 4 vertices
  const {normalArray} = geoGroup;
  for (i = 0; i < 4; i++) {
    colorArray[start + 3 * i] = C.r;
    colorArray[start + 3 * i + 1] = C.g;
    colorArray[start + 3 * i + 2] = C.b;
  }

  normalArray[start + 0] = -radius;
  normalArray[start + 1] = radius;
  normalArray[start + 2] = 0;

  normalArray[start + 3] = -radius;
  normalArray[start + 4] = -radius;
  normalArray[start + 5] = 0;

  normalArray[start + 6] = radius;
  normalArray[start + 7] = -radius;
  normalArray[start + 8] = 0;

  normalArray[start + 9] = radius;
  normalArray[start + 10] = radius;
  normalArray[start + 11] = 0;

  geoGroup.vertices += 4;

  // two faces
  const {faceArray} = geoGroup;
  const faceoffset = geoGroup.faceidx; // not number faces, but index
  faceArray[faceoffset + 0] = startv;
  faceArray[faceoffset + 1] = startv + 1;
  faceArray[faceoffset + 2] = startv + 2;
  faceArray[faceoffset + 3] = startv + 2;
  faceArray[faceoffset + 4] = startv + 3;
  faceArray[faceoffset + 5] = startv;
  geoGroup.faceidx += 6;
}

function negateColor(c) {
  // set sign bit
  let n = -c;
  if (n == 0) n = -0.0001;
  return n;
}

function drawStickImposter(geo, from, to, radius, color) {
  // we need the four corners - two have from coord, two have to coord, the normal
  // is the opposing point, from which we can get the normal and length
  // also need the radius
  const geoGroup = geo.updateGeoGroup(4);
  const startv = geoGroup.vertices;
  const start = startv * 3;
  const {vertexArray} = geoGroup;
  const {colorArray} = geoGroup;
  const {radiusArray} = geoGroup;
  const {normalArray} = geoGroup;
  // encode extra bits of information in the color
  const {r} = color;
  const {g} = color;
  const {b} = color;

  /* for sticks, always draw caps, but we could in theory set caps in color */
  // 4 vertices, distinguish between p1 and p2 with neg blue
  let pos = start;
  for (let i = 0; i < 4; i++) {
    vertexArray[pos] = from.x;
    normalArray[pos] = to.x;
    colorArray[pos] = r;
    pos++;
    vertexArray[pos] = from.y;
    normalArray[pos] = to.y;
    colorArray[pos] = g;
    pos++;
    vertexArray[pos] = from.z;
    normalArray[pos] = to.z;
    if (i < 2) colorArray[pos] = b;
    else colorArray[pos] = negateColor(b);
    pos++;
  }

  geoGroup.vertices += 4;

  radiusArray[startv] = -radius;
  radiusArray[startv + 1] = radius;
  radiusArray[startv + 2] = -radius;
  radiusArray[startv + 3] = radius;

  // two faces
  const {faceArray} = geoGroup;
  const faceoffset = geoGroup.faceidx; // not number faces, but index
  faceArray[faceoffset + 0] = startv;
  faceArray[faceoffset + 1] = startv + 1;
  faceArray[faceoffset + 2] = startv + 2;
  faceArray[faceoffset + 3] = startv + 2;
  faceArray[faceoffset + 4] = startv + 3;
  faceArray[faceoffset + 5] = startv;
  geoGroup.faceidx += 6;
}

function adjustCoord(x1, x2, margin, adjust) {
  // return new value of x2 that isn't more than margin away
  const dist = x2 - x1;
  if (dist < -margin) {
    return x2 + adjust;
  }
  if (dist > margin) {
    return x2 - adjust;
  }
  return x2;
}

// return true if atom value matches property val
function propertyMatches(atomval, val) {
  if (atomval == val) {
    return true;
  }
  if (typeof val == 'string' && typeof atomval == 'number') {
    // support numerical integer ranges, e.g. resi: 3-7
    const match = val.match(/(-?\d+)\s*-\s*(-?\d+)/);
    if (match) {
      // eslint-disable-next-line radix
      const lo = Number.parseInt(match[1]);
      // eslint-disable-next-line radix
      const hi = Number.parseInt(match[2]);
      if (match && atomval >= lo && atomval <= hi) {
        return true;
      }
    }
  }
  return false;
}

// make a deep copy of a selection object and create caches of expensive
// selections.  We create a copy so caches are not attached to user
// supplied objects where the user might change them invalidating the cache.
// This does not support arbitrary
// javascript objects, but support enough for eveything that is
// used in selections: number, string, boolean, functions; as well
// as arrays and nested objects with values of the aformentioned
// types.
function deepCopyAndCache(selobject, model) {
  if (typeof selobject != 'object' || selobject === null) return selobject;
  if (selobject.__cache_created) return selobject; // already done
  const copy = {};
  for (const key in selobject) {
    const item = selobject[key];
    if (Array.isArray(item)) {
      // handle array separatly from other typeof == "object"
      // elements
      copy[key] = [];
      for (let i = 0; i < item.length; i++) {
        copy[key].push(deepCopyAndCache(item[i], model));
      }
    } else if (typeof item === 'object' && key != 'properties') {
      copy[key] = deepCopyAndCache(item, model);
    } else {
      copy[key] = item;
    }

    // create caches of expensive selection types - the cache
    // stores the atoms matching the selection type
    if (key == 'and' || key == 'or') {
      // create a list of sets of matching atoms indexes for
      // each sub-selection
      const results = [];
      for (const subSelection of copy[key]) {
        const set = new Set();
        for (const match of model.selectedAtoms(subSelection)) {
          set.add(match.index);
        }
        results.push(set);
      }

      if (key == 'and') {
        // get the intersection of two sets
        const intersect = function intersect(first, other) {
          const result = new Set();
          for (const elem of other) {
            if (first.has(elem)) {
              result.add(elem);
            }
          }
          return result;
        };

        let intersection = new Set(results[0]);
        for (const set of results.splice(1)) {
          intersection = intersect(intersection, set);
        }
        copy[key].__cached_results = intersection;
      } else if (key == 'or') {
        const union = new Set();
        for (const set of results) {
          for (const elem of Array.from(set)) {
            union.add(elem);
          }
        }

        copy[key].__cached_results = union;
      }
    }
  }
  copy.__cache_created = true;
  return copy;
}

function squaredDistance(atom1, atom2) {
  const xd = atom2.x - atom1.x;
  const yd = atom2.y - atom1.y;
  const zd = atom2.z - atom1.z;
  return xd * xd + yd * yd + zd * zd;
}

export default class GLModel {
  /**
   * The dictionaries are for dropdown menus and validation of the viewer
   */
  static validElements = [
    'H',
    'Li',
    'LI',
    'Na',
    'NA',
    'K',
    'C',
    'N',
    'O',
    'F',
    'P',
    'S',
    'CL',
    'Cl',
    'BR',
    'Br',
    'SE',
    'Se',
    'ZN',
    'Zn',
    'CU',
    'Cu',
    'NI',
    'Ni',
  ];

  static validAtomSpecs = {
    resn: {type: 'string', valid: true, prop: true, gui: true}, // Parent residue name
    x: {type: 'number', floatType: true, valid: false, step: 0.1, prop: true}, // Atom's x coordinate
    y: {type: 'number', floatType: true, valid: false, step: 0.1, prop: true}, // Atom's y coordinate
    z: {type: 'number', floatType: true, valid: false, step: 0.1, prop: true}, // Atom's z coordinate
    color: {type: 'color', gui: false}, // Atom's color, as hex code
    surfaceColor: {type: 'color', gui: false}, // Hex code for color to be used for surface patch over this atom
    elem: {type: 'element', gui: true, prop: true}, // Element abbreviation (e.g. 'H', 'Ca', etc)
    hetflag: {type: 'boolean', valid: false, gui: true}, // Set to true if atom is a heteroatom
    chain: {type: 'string', gui: true, prop: true}, // Chain this atom belongs to, if specified in input file (e.g 'A' for chain A)
    resi: {type: 'array_range', gui: true}, // Residue number
    icode: {type: 'number', valid: false, step: 0.1},
    rescode: {type: 'number', valid: false, step: 0.1, prop: true},
    serial: {type: 'number', valid: false, step: 0.1}, // Atom's serial id numbermodels
    atom: {type: 'string', valid: false, gui: true, prop: true}, // Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon)
    bonds: {type: 'array', valid: false}, // Array of atom ids this atom is bonded to
    ss: {type: 'string', valid: false}, // Secondary structure identifier (for cartoon render; e.g. 'h' for helix)
    singleBonds: {type: 'boolean', valid: false}, // true if this atom forms only single bonds or no bonds at all
    bondOrder: {type: 'array', valid: false}, // Array of this atom's bond orders, corresponding to bonds identfied by 'bonds'
    properties: {type: 'properties', valid: false}, // Optional mapping of additional properties
    b: {type: 'number', floatType: true, valid: false, step: 0.1, prop: true}, // Atom b factor data
    pdbline: {type: 'string', valid: false}, // If applicable, this atom's record entry from the input PDB file (used to output new PDB from models)
    clickable: {type: 'boolean', valid: false, gui: false}, // Set this flag to true to enable click selection handling for this atom
    contextMenuEnabled: {type: 'boolean', valid: false, gui: false}, // Set this flag to true to enable click selection handling for this atom
    callback: {type: 'function', valid: false}, // Callback click handler function to be executed on this atom and its parent viewer
    invert: {type: 'boolean', valid: false}, // for selection, inverts the meaning of the selection
    // unsure about this
    reflectivity: {type: 'number', floatType: true, gui: false, step: 0.1}, // for describing the reflectivity of a model
    altLoc: {type: 'invalid', valid: false}, // alternative location, e.g. in PDB
    sym: {type: 'number', gui: false}, // which symmetry
  };

  static validSurfaceSpecs = {
    opacity: {type: 'number', floatType: true, gui: true, step: 0.01, default: 1, min: 0, max: 1},
    colorscheme: {type: 'colorscheme', gui: true},
    color: {type: 'color', gui: true},
    voldata: {type: 'number', floatType: true, gui: false},
    volscheme: {type: 'number', floatType: true, gui: false},
    map: {type: 'number', gui: false},
  };

  static validAtomStyleSpecs = {
    line: {validItems: validLineSpec, type: 'form', gui: true}, // draw bonds as lines
    cross: {validItems: validCrossSpec, type: 'form', gui: true}, // draw atoms as crossed lines (aka stars)
    stick: {validItems: validStickSpec, type: 'form', gui: true}, // draw bonds as capped cylinders
    sphere: {validItems: validSphereSpec, type: 'form', gui: true}, // draw atoms as spheres
    cartoon: {validItems: validCartoonSpec, type: 'form', gui: true}, // draw cartoon representation of secondary structure
    colorfunc: {validItems: null, type: 'js', valid: false},
    clicksphere: {validItems: validSphereSpec, type: 'form'}, // invisible style for click handling
  };

  static validAtomSelectionSpecs = extend(GLModel.validAtomSpecs, validExtras);

  static validLabelResSpecs = {
    font: {type: 'string', gui: true},
    fontSize: {type: 'number', floatType: true, gui: true, step: 1, default: 12, min: 1},
    fontColor: {type: 'color', gui: true},
    fontOpacity: {
      type: 'number',
      floatType: true,
      gui: true,
      step: 0.1,
      default: 1,
      min: 0,
      max: 1,
    },
    borderThickness: {type: 'number', floatType: true, gui: true, step: 0.1, default: 1, min: 0},
    borderColor: {type: 'color', gui: true},
    borderOpacity: {
      type: 'number',
      floatType: true,
      gui: true,
      step: 0.1,
      default: 1,
      min: 0,
      max: 1,
    },
    backgroundColor: {type: 'color', gui: true},
    backgroundOpacity: {
      type: 'number',
      floatType: true,
      gui: true,
      step: 0.1,
      default: 1,
      min: 0,
      max: 1,
    },
    position: {type: 'array', valid: false},
    inFront: {type: 'boolean', gui: true},
    showBackground: {type: 'boolean', gui: true},
    fixed: {type: 'boolean', gui: true},
    alignment: {
      validItems: [
        'topLeft',
        'topCenter',
        'topRight',
        'centerLeft',
        'center',
        'centerRight',
        'bottomLeft',
        'bottomCenter',
        'bottomRight',
      ],
      gui: true,
    },
    scale: {type: 'boolean', gui: false},
  };

  // private variables
  atoms = [];
  /** @type {any[] & Partial<{numFrames: any; url:any; origIndex:any; path:any}>} */
  frames = [];
  /** @type {number[]|null} */
  box = null;
  /** @type {Array<any> | null} */
  atomdfs = null; // depth first search over connected components
  hidden = false;
  /** @type{Object3D|null} */
  molObj = null;
  /** @type {Object3D|null} */
  renderedMolObj = null;
  /** @type {Color[]|null} */
  lastColors = null;
  modelData = {};
  /** @type {Array<any> | null} */
  modelDatas = null; // if there is different modelData per frame
  idMatrix = new Matrix4();
  dontDuplicateAtoms = true;
  /** @type {boolean|undefined|{shapes: any[], labels: any[]}} */
  unitCellObjects = undefined;

  /**
   * @param {number} mid
   * @param {Record<string, any>} [options]
   */
  constructor(mid, options) {
    this.id = mid;
    const {defaultColor} = elementColors;
    this.defaultColor = defaultColor;

    options = options || {};
    this.ElementColors = options.defaultcolors
      ? options.defaultcolors
      : elementColors.defaultColors;

    // drawing functions must be associated with model object since
    // geometries can't span multiple canvases
    // sphere drawing
    this.defaultSphereRadius = options.defaultSphereRadius ? options.defaultSphereRadius : 1.5;
    this.defaultCartoonQuality = options.cartoonQuality ? options.cartoonQuality : 5;
  }

  /**
   * Return object representing internal state of
   * the model appropriate for passing to setInternalState
   *
   * @function GLModel#getInternalState
   */
  getInternalState() {
    return {atoms: this.atoms, frames: this.frames};
  }

  /**
   * Overwrite the internal model state with the passed state.
   *
   * @function GLModel#setInternalState
   */
  setInternalState(state) {
    this.atoms = state.atoms;
    this.frames = state.frames;
    this.molObj = null;
  }

  /**
   * Returns crystallographic information if present.
   *
   * @function GLModel#getCrystData
   *
   */
  getCrystData() {
    if (this.modelData.cryst) {
      // add the matrix if it is missing
      if (!this.modelData.cryst.matrix) {
        const {cryst} = this.modelData;
        this.modelData.cryst.matrix = conversionMatrix3(
          cryst.a,
          cryst.b,
          cryst.c,
          cryst.alpha,
          cryst.beta,
          cryst.gamma
        );
      }
      return this.modelData.cryst;
    }
    return null;
  }

  /**
         * Set crystallographic information using three angles and three lengths
         *
         * @function GLModel#setCrystData
         * @param {number} a - length of unit cell side
         * @param {number} b - length of unit cell side
         * @param {number} c - length of unit cell side
         * @param {number} alpha - unit cell angle in degrees (default 90)
         * @param {number} beta - unit cell angle in degrees (default 90)
         * @param {number} gamma - unit cell angle in degrees (default 90)
         
         */
  setCrystData(a, b, c, alpha, beta, gamma) {
    // I am assuming these
    a = a || 1.0;
    b = b || 1.0;
    c = c || 1.0;
    alpha = alpha || 90;
    beta = beta || 90;
    gamma = gamma || 90;

    const matrix = conversionMatrix3(a, b, c, alpha, beta, gamma);
    this.modelData.cryst = {
      a,
      b,
      c,
      alpha,
      beta,
      gamma,
      matrix,
    };
  }

  /**
   * Set the crystallographic matrix to the given matrix.
   *
   * This function removes `a`, `b`, `c`, `alpha`, `beta`, `gamma` from
   * the crystal data.
   *
   * @function GLModel#setCrystMatrix
   * @param {Matrix3} matrix - unit cell matrix
   */
  setCrystMatrix(matrix) {
    matrix = matrix || new Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);

    this.modelData.cryst = {
      matrix,
    };
  }

  /**
   * Returns list of rotational/translational matrices if there is BIOMT data
   * Otherwise returns a list of just the ID matrix
   *
   * @function GLModel#getSymmetries
   * @return {Array<Matrix4>}
   *
   */
  getSymmetries() {
    if (typeof this.modelData.symmetries == 'undefined') {
      this.modelData.symmetries = [this.idMatrix];
    }
    return this.modelData.symmetries;
  }

  /**
   * Sets symmetries based on specified matrices in list
   *
   * @function GLModel#setSymmetries
   * @param {Array<Matrix4>} list
   *
   */
  setSymmetries(list) {
    if (typeof list == 'undefined') {
      // delete sym data
      this.modelData.symmetries = [this.idMatrix];
    } else {
      this.modelData.symmetries = list;
    }
  }

  /**
   * Returns model id number
   *
   * @function GLModel#getID
   * @return {number} Model ID
   */
  getID() {
    return this.id;
  }

  /**
   * Returns model's frames property, a list of atom lists
   *
   * @function GLModel#getNumFrames
   * @return {number}
   */
  getNumFrames() {
    return this.frames.numFrames != undefined ? this.frames.numFrames : this.frames.length;
  }

  /**
   * Sets model's atomlist to specified frame
   * Sets to last frame if framenum out of range
   *
   * @function GLModel#setFrame
   * @param {number} framenum - model's atoms are set to this index in frames list
   * @param {import("./GLViewer").default} [viewer]
   * @return {Promise<void>}
   */
  setFrame(framenum, viewer) {
    // viewer only passed internally for unit cell
    const numFrames = this.getNumFrames();
    const model = this;
    return new Promise((resolve, reject) => {
      if (numFrames == 0) {
        // return;
        resolve();
      }
      if (framenum < 0 || framenum >= numFrames) {
        framenum = numFrames - 1;
      }
      if (model.frames.url != undefined) {
        const {url} = model.frames;
        getbin(`${url}/traj/frame/${framenum}/${model.frames.path}`, undefined, 'POST')
          .then(buffer => {
            const values = new Float32Array(buffer, 44);
            let count = 0;
            for (let i = 0; i < model.atoms.length; i++) {
              model.atoms[i].x = values[count++];
              model.atoms[i].y = values[count++];
              model.atoms[i].z = values[count++];
            }
            // if a box was provided, check to see if we need to wrap connected components
            if (model.box && model.atomdfs) {
              this.adjustCoordinatesToBox();
            }
            resolve();
          })
          .catch(reject);
      } else {
        model.atoms = this.frames[framenum];
        resolve();
      }
      model.molObj = null;
      if (model.modelDatas && framenum < model.modelDatas.length) {
        model.modelData = model.modelDatas[framenum];
        if (model.unitCellObjects && viewer) {
          viewer.removeUnitCell(model);
          viewer.addUnitCell(model);
        }
      }
    });
  }

  /**
   * Add atoms as frames of model
   *
   * @function GLModel#addFrame
   * @param {import('./specs').AtomSpec} atom - atoms to be added
   */
  addFrame(atom) {
    this.frames.push(atom);
  }

  /**
   * If model atoms have dx, dy, dz properties (in some xyz files), vibrate populates the model's frame property based on parameters.
   * Model can then be animated
   *
   * @function GLModel#vibrate
   * @param {number} numFrames - number of frames to be created, default to 10
   * @param {number} amplitude - amplitude of distortion, default to 1 (full)
   * @param {boolean} [bothWays] - if true, extend both in positive and negative directions by numFrames
   * @param {import("./GLViewer").default} [viewer] - required if arrowSpec is provided
   * @param {import('./specs').ArrowSpec} [arrowSpec] - specification for drawing animated arrows. If color isn't specified, atom color (sphere, stick, line preference) is used.
   * @example
    $3Dmol.download("pdb:4UAA",viewer,{},function(){
      viewer.setStyle({},{stick:{}});
      viewer.vibrate(10, 1);
      viewer.animate({loop: "forward",reps: 1});
      viewer.zoomTo();
            viewer.render();
        });
   */
  vibrate(numFrames, amplitude, bothWays, viewer, arrowSpec) {
    amplitude = amplitude || 1;
    numFrames = numFrames || 10;
    let start = 0;
    let end = numFrames;
    if (bothWays) {
      start = -numFrames;
      end = numFrames;
    }

    // to enable multiple setting of vibrate with bothWays, must record original position
    if (this.frames !== undefined && this.frames.origIndex !== undefined) {
      this.setFrame(this.frames.origIndex);
    } else {
      this.setFrame(0);
    }

    if (start < end) this.frames = []; // clear
    if (bothWays) this.frames.origIndex = numFrames;

    for (let i = start; i < end; i++) {
      const newAtoms = [];
      const currframe = this.frames.length;
      if (i == 0 && !arrowSpec) {
        // still need to calculate if drawing arrows
        this.frames.push(this.atoms);
        continue;
      }
      for (let j = 0; j < this.atoms.length; j++) {
        const dx = getAtomProperty(this.atoms[j], 'dx');
        const dy = getAtomProperty(this.atoms[j], 'dy');
        const dz = getAtomProperty(this.atoms[j], 'dz');
        const newVector = new Vector3(dx, dy, dz);
        const starting = new Vector3(this.atoms[j].x, this.atoms[j].y, this.atoms[j].z);
        const mult = (i * amplitude) / numFrames;
        newVector.multiplyScalar(mult);
        starting.add(newVector);
        /** @type {import('./specs').AtomSpec} */
        const newAtom = {};
        for (const k in this.atoms[j]) {
          newAtom[k] = this.atoms[j][k];
        }
        newAtom.x = starting.x;
        newAtom.y = starting.y;
        newAtom.z = starting.z;
        newAtoms.push(newAtom);
        if (viewer && arrowSpec) {
          const spec = extend({}, arrowSpec);
          const arrowend = new Vector3(dx, dy, dz);
          arrowend.multiplyScalar(amplitude);
          arrowend.add(starting);

          spec.start = starting;
          spec.end = arrowend;
          spec.frame = currframe;
          if (!spec.color) {
            let s = newAtom.style.sphere;
            if (!s) s = newAtom.style.stick;
            if (!s) s = newAtom.style.line;
            spec.color = getColorFromStyle(newAtom, s);
          }
          viewer.addArrow(spec);
        }
      }
      this.frames.push(newAtoms);
    }
  }

  // set default style and colors for atoms
  setAtomDefaults(atoms) {
    for (let i = 0; i < atoms.length; i++) {
      const atom = atoms[i];
      if (atom) {
        atom.style = atom.style || deepCopy(defaultAtomStyle);
        atom.color = atom.color || this.ElementColors[atom.elem] || this.defaultColor;
        atom.model = this.id;
        if (atom.clickable || atom.hoverable)
          atom.intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
      }
    }
  }

  /** add atoms to this model from molecular data string
   *
   * @function GLModel#addMolData
   * @param {string|Uint8Array} [data] - atom structure file input data string, for gzipped input use ArrayBuffer
   * @param {string} [format] - input file string format (e.g 'pdb', 'sdf', 'sdf.gz', etc.)
   * @param {import('./specs').ParserOptionsSpec} [options] - format dependent options. Attributes depend on the input format
   */
  addMolData(data, format, options) {
    options = options || {};
    if (!data) throw new Error("No data provided ot GLModel.addMolData");
    const parsedAtoms = GLModel.parseMolData(data, format, options);
    this.dontDuplicateAtoms = !options.duplicateAssemblyAtoms;
    const mData = parsedAtoms.modelData;
    if (mData) {
      if (Array.isArray(mData)) {
        this.modelData = mData[0];
        if (options.frames) {
          this.modelDatas = mData;
        }
      } else {
        this.modelData = mData;
      }
    }

    if (parsedAtoms.box) {
      this.box = parsedAtoms.box;
    } else {
      this.box = null;
    }

    if (this.frames.length == 0) {
      // first call
      for (let i = 0; i < parsedAtoms.length; i++) {
        if (parsedAtoms[i].length != 0) this.frames.push(parsedAtoms[i]);
      }
      if (this.frames[0]) this.atoms = this.frames[0];
    } else {
      // subsequent calls
      if (options.frames) {
        // add to new frame
        for (let i = 0; i < parsedAtoms.length; i++) {
          this.frames.push(parsedAtoms[i]);
        }
      } else {
        // add atoms to current frame
        for (let i = 0; i < parsedAtoms.length; i++) {
          this.addAtoms(/** @type {import('./specs').AtomSpec[]} */(parsedAtoms[i]));
        }
      }
    }

    for (let i = 0; i < this.frames.length; i++) {
      this.setAtomDefaults(this.frames[i]);
    }

    if (options.vibrate && options.vibrate.frames && options.vibrate.amplitude) {
      // fill in vibrational modes
      this.vibrate(options.vibrate.frames, options.vibrate.amplitude);
    }

    if (options.style) {
      this.setStyle({}, options.style);
    }
  }

  setDontDuplicateAtoms(dup) {
    this.dontDuplicateAtoms = dup;
  }

  setModelData(mData) {
    this.modelData = mData;
  }

  /** return list of atoms selected by sel, this is specific to glmodel
   *
   * @function GLModel#selectedAtoms
   * @param {import('./specs').AtomSelectionSpec} sel
   * @param {any} [from]
   * @return {Array<import('./specs').AtomSpec>}
   * @example
      $3Dmol.download("pdb:4wwy",viewer,{},function(){
            var atoms = viewer.selectedAtoms({chain:'A'});
            for(var i = 0, n = atoms.length; i < n; i++) {
               atoms[i].b = 0.0;
            }
            viewer.setStyle({cartoon:{colorscheme:{prop:'b',gradient: 'roygb',min:0,max:30}}});
            viewer.render();
        });
   */
  selectedAtoms(sel, from) {
    let ret = [];

    // make a copy of the selection to allow caching results without
    // the possibility for the user to change the selection and this
    // code not noticing the changes
    sel = deepCopyAndCache(sel || {}, this);

    if (!from) from = this.atoms;
    const aLength = from.length;
    for (let i = 0; i < aLength; i++) {
      const atom = from[i];
      if (atom) {
        if (this.atomIsSelected(atom, sel)) ret.push(atom);
      }
    }

    // expand selection by some distance
    if (sel.expand) {
      // get atoms in expanded bounding box
      const exdist = parseFloat(/** @type {string} */ (sel.expand));
      const expand = this.expandAtomList(ret, exdist);
      const retlen = ret.length;
      const thresh = exdist * exdist;
      for (let i = 0; i < expand.length; i++) {
        for (let j = 0; j < retlen; j++) {
          const dist = squaredDistance(expand[i], ret[j]);
          if (dist < thresh && dist > 0) {
            ret.push(expand[i]);
          }
        }
      }
    }

    // selection within distance of sub-selection
    if (sel.within && sel.within.sel && sel.within.distance) {
      // get atoms in second selection
      const sel2 = this.selectedAtoms(sel.within.sel, this.atoms);
      const within = {};
      let dist = parseFloat(/** @type {string} */ (sel.within.distance));
      const thresh = dist * dist;
      for (let i = 0; i < sel2.length; i++) {
        for (let j = 0; j < ret.length; j++) {
          dist = squaredDistance(sel2[i], ret[j]);
          if (dist < thresh && dist > 0) {
            within[j] = 1;
          }
        }
      }
      const newret = [];
      if (sel.within.invert) {
        for (let j = 0; j < ret.length; j++) {
          if (!within[j]) newret.push(ret[j]);
        }
      } else {
        for (const j in within) {
          newret.push(ret[j]);
        }
      }
      ret = newret;
    }

    // byres selection flag
    if (sel.byres) {
      // Keep track of visited residues, visited atoms, and atom stack
      const vResis = {};
      const vAtoms = [];
      const stack = [];

      for (let i = 0; i < ret.length; i++) {
        // Check if atom is part of a residue, and that the residue hasn't been traversed yet
        let atom = ret[i];
        let c = atom.chain;
        let r = atom.resi;
        if (vResis[c] === undefined) vResis[c] = {};
        if (atom.resi && vResis[c][r] === undefined) {
          // Perform a depth-first search of atoms with the same resi
          vResis[c][r] = true;
          stack.push(atom);
          while (stack.length > 0) {
            atom = stack.pop();
            c = atom.chain;
            r = atom.resi;
            if (vAtoms[atom.index] === undefined) {
              vAtoms[atom.index] = true;
              for (let j = 0; j < atom.bonds.length; j++) {
                const atom2 = this.atoms[atom.bonds[j]];
                if (
                  vAtoms[atom2.index] === undefined &&
                  atom2.resi &&
                  atom2.chain == c &&
                  atom2.resi == r
                ) {
                  stack.push(atom2);
                  ret.push(atom2);
                }
              }
            }
          }
        }
      }
    }

    return ret;
  }

  /** Add list of new atoms to model.  Adjusts bonds appropriately.
   *
   * @function GLModel#addAtoms
   * @param {import('./specs').AtomSpec[]} newatoms
   * @example
   * var atoms = [{elem: 'C', x: 0, y: 0, z: 0, bonds: [1,2], bondOrder: [1,2]}, {elem: 'O', x: -1.5, y: 0, z: 0, bonds: [0]},{elem: 'O', x: 1.5, y: 0, z: 0, bonds: [0], bondOrder: [2]}];
     
      viewer.setBackgroundColor(0xffffffff);
      var m = viewer.addModel();
      m.addAtoms(atoms);
      m.setStyle({},{stick:{}});
      viewer.zoomTo();
      viewer.render();
   */
  addAtoms(newatoms) {
    this.molObj = null;
    const start = this.atoms.length;
    const indexmap = [];
    // mapping from old index to new index
    let i;
    for (i = 0; i < newatoms.length; i++) {
      if (typeof newatoms[i].index == 'undefined') newatoms[i].index = i;
      if (typeof newatoms[i].serial == 'undefined') newatoms[i].serial = i;
      indexmap[newatoms[i].index] = start + i;
    }

    // copy and push newatoms onto atoms
    for (i = 0; i < newatoms.length; i++) {
      const olda = newatoms[i];
      const nindex = indexmap[olda.index];
      const a = extend({}, olda);
      a.index = nindex;
      a.bonds = [];
      a.bondOrder = [];
      a.model = this.id;
      a.style = a.style || deepCopy(defaultAtomStyle);
      if (typeof a.color == 'undefined') a.color = this.ElementColors[a.elem] || this.defaultColor;
      // copy over all bonds contained in selection,
      // updating indices appropriately
      const nbonds = olda.bonds ? olda.bonds.length : 0;
      for (let j = 0; j < nbonds; j++) {
        const neigh = indexmap[olda.bonds[j]];
        if (typeof neigh != 'undefined') {
          a.bonds.push(neigh);
          a.bondOrder.push(olda.bondOrder ? olda.bondOrder[j] : 1);
        }
      }
      this.atoms.push(a);
    }
  }

  /** Remove specified atoms from model
   *
   * @function GLModel#removeAtoms
   * @param {Array<any>} badatoms - list of atoms
   */
  removeAtoms(badatoms) {
    this.molObj = null;
    // make map of all baddies
    const baddies = [];
    let i;
    for (i = 0; i < badatoms.length; i++) {
      baddies[badatoms[i].index] = true;
    }

    // create list of good atoms
    const newatoms = [];
    for (i = 0; i < this.atoms.length; i++) {
      const a = this.atoms[i];
      if (!baddies[a.index]) newatoms.push(a);
    }

    // clear it all out
    this.atoms = [];
    // and add back in to get updated bonds
    this.addAtoms(newatoms);
  }

  /** Set atom style of selected atoms
   *
   * @function GLModel#setStyle
   * @param {import('./specs').AtomSelectionSpec | import('./specs').AtomStyleSpec} sel
   * @param {import('./specs').AtomStyleSpec} [style]
   * @param {boolean} [add] - if true, add to current style, don't replace
   @example
    $3Dmol.download("pdb:4UB9",viewer,{},function(){
            viewer.setBackgroundColor(0xffffffff);
            viewer.setStyle({chain:'A'},{line:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.setStyle({chain:'B'},{line:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.setStyle({chain:'C'},{cross:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.setStyle({chain:'D'},{cross:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.setStyle({chain:'E'},{cross:{radius:2.0,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.setStyle({chain:'F'},{stick:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.setStyle({chain:'G'},{stick:{radius:0.8,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.setStyle({chain:'H'},{stick:{singleBonds:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
            viewer.render();
        });
   */
  setStyle(sel, style, add) {
    if (typeof style === 'undefined' && typeof add == 'undefined') {
      // if a single argument is provided, assume it is a style and select all
      style = /** @type {import('./specs').AtomStyleSpec} */ (sel);
      sel = {};
    }

    // if type is just a string, promote it to an object
    if (typeof style === 'string') {
      style = specStringToObject(style);
    }
    // report to console if this is not a valid selector
    for (const s in sel) {
      if (!GLModel.validAtomSelectionSpecs[s]) {
        console.log(`Unknown selector ${s}`);
      }
    }
    // report to console if this is not a valid style
    for (const s in style) {
      if (!GLModel.validAtomStyleSpecs[s]) {
        console.log(`Unknown style ${s}`);
      }
    }

    let changedAtoms = false;
    // somethings we only calculate if there is a change in a certain
    // style, although these checks will only catch cases where both
    // are either null or undefined
    const setStyleHelper = atomArr => {
      const selected = this.selectedAtoms(
        /** @type {import('./specs').AtomSelectionSpec} */ (sel),
        atomArr
      );
      for (let i = 0; i < atomArr.length; i++) {
        if (atomArr[i]) atomArr[i].capDrawn = false; // reset for proper stick render
      }

      for (let i = 0; i < selected.length; i++) {
        changedAtoms = true;
        if (selected[i].clickable || selected[i].hoverable)
          selected[i].intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};

        if (!add) selected[i].style = {};
        for (const s in style) {
          if (style[s]) {
            selected[i].style[s] = selected[i].style[s] || {}; // create distinct object for each atom
            Object.assign(selected[i].style[s], style[s]);
          }
        }
      }
    };

    setStyleHelper(this.atoms);
    for (let i = 0; i < this.frames.length; i++) {
      if (this.frames[i] !== this.atoms) setStyleHelper(this.frames[i]);
    }

    if (changedAtoms) this.molObj = null; // force rebuild
  }

  /** Set clickable and callback of selected atoms
   *
   * @function GLModel#setClickable
   * @param {import('./specs').AtomSelectionSpec} sel - atom selection to apply clickable settings to
   * @param {boolean} clickable - whether click-handling is enabled for the selection
   * @param {import('./specs').AnyFunc|string} callbackSrc - function called when an atom in the selection is clicked
   */
  setClickable(sel, clickable, callbackSrc) {
    // report to console if this is not a valid selector
    let s;
    for (s in sel) {
      if (!GLModel.validAtomSelectionSpecs[s]) {
        console.log(`Unknown selector ${s}`);
      }
    }

    // make sure clickable is a boolean
    clickable = !!clickable;
    const callback = makeFunction(callbackSrc);
    if (callback === null) {
      console.log('Callback is not a function');
      return;
    }

    let i;
    const selected = this.selectedAtoms(sel, this.atoms);
    const len = selected.length;
    for (i = 0; i < len; i++) {
      selected[i].intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
      selected[i].clickable = clickable;
      if (callback) selected[i].callback = callback;
    }

    if (len > 0) this.molObj = null; // force rebuild to get correct intersection shapes
  }

  /**
   * Set hoverable and callback of selected atoms
   * @function GLModel#setHoverable
   * @param {import('./specs').AtomSelectionSpec} sel - atom selection to apply hoverable settings to
   * @param {boolean} hoverable - whether hover-handling is enabled for the selection
   * @param {import('./specs').AnyFunc | string} hoverCallbackSrc - function called when an atom in the selection is hovered over
   * @param {import('./specs').AnyFunc | string} unhoverCallbackSrc - function called when the mouse moves out of the hover area
   */
  setHoverable(sel, hoverable, hoverCallbackSrc, unhoverCallbackSrc) {
    let s;
    for (s in sel) {
      if (!GLModel.validAtomSelectionSpecs[s]) {
        console.log(`Unknown selector ${s}`);
      }
    }

    // make sure hoverable is a boolean
    hoverable = !!hoverable;
    const hoverCallback = makeFunction(hoverCallbackSrc);
    const unhoverCallback = makeFunction(unhoverCallbackSrc);

    // report to console if hoverCallback is not a valid function
    if (hoverCallback === null) {
      console.log('Hover_callback is not a function');
      return;
    }
    // report to console if unhoverCallback is not a valid function
    if (unhoverCallback === null) {
      console.log('Unhover_callback is not a function');
      return;
    }

    let i;
    const selected = this.selectedAtoms(sel, this.atoms);
    const len = selected.length;
    for (i = 0; i < len; i++) {
      selected[i].intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
      selected[i].hoverable = hoverable;
      if (hoverCallback) selected[i].hoverCallback = hoverCallback;
      if (unhoverCallback) selected[i].unhoverCallback = unhoverCallback;
    }

    if (len > 0) this.molObj = null; // force rebuild to get correct intersection shapes
  }

  /** enable context menu of selected atoms
   *
   * @function GLModel#enableContextMenu
   * @param {import('./specs').AtomSelectionSpec} sel - atom selection to apply hoverable settings to
   * @param {boolean} contextMenuEnabled - whether contextMenu-handling is enabled for the selection
   */
  enableContextMenu(sel, contextMenuEnabled) {
    let s;
    for (s in sel) {
      if (!GLModel.validAtomSelectionSpecs[s]) {
        console.log(`Unknown selector ${s}`);
      }
    }

    // make sure contextMenuEnabled is a boolean
    contextMenuEnabled = !!contextMenuEnabled;

    let i;
    const selected = this.selectedAtoms(sel, this.atoms);
    const len = selected.length;
    for (i = 0; i < len; i++) {
      selected[i].intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
      selected[i].contextMenuEnabled = contextMenuEnabled;
    }

    if (len > 0) this.molObj = null; // force rebuild to get correct intersection shapes
  }

  /** given a mapping from element to color, set atom colors
   *
   * @function GLModel#setColorByElement
   * @param {import("./specs").AtomSelectionSpec} sel
   * @param {Color[]} colors
   */
  setColorByElement(sel, colors) {
    const atoms = this.selectedAtoms(sel, this.atoms);
    if (this.molObj !== null && sameObj(colors, this.lastColors)) return; // don't recompute
    this.lastColors = colors;

    if (atoms.length > 0) this.molObj = null; // force rebuild
    for (let i = 0; i < atoms.length; i++) {
      const a = atoms[i];
      if (typeof colors[a.elem] !== 'undefined') {
        a.color = colors[a.elem];
      }
    }
  }

  /**
   * @function GLModel.setColorByProperty
   * @param {import("./specs").AtomSelectionSpec} sel
   * @param {string} prop
   * @param {import("./Gradient").default} scheme
   */
  setColorByProperty(sel, prop, scheme, range) {
    let i;
    let a;
    const atoms = this.selectedAtoms(sel, this.atoms);
    this.lastColors = null; // don't bother memoizing
    if (atoms.length > 0) this.molObj = null; // force rebuild

    if (typeof Gradient.builtinGradients[scheme] != 'undefined') {
      scheme = new Gradient.builtinGradients[scheme]();
    }

    if (!range) {
      // no explicit range, get from scheme
      range = scheme.range();
    }

    if (!range) {
      // no range in scheme, compute the range for this model
      range = getPropertyRange(atoms, prop);
    }
    // now apply colors using scheme
    for (i = 0; i < atoms.length; i++) {
      a = atoms[i];
      const val = getAtomProperty(a, prop);
      if (val != null) {
        a.color = scheme.valueToHex(parseFloat(a.properties[prop]), range);
      }
    }
  }

  /**
   * @function GLModel#setColorByFunction
   * @deprecated use setStyle and colorfunc attribute
   * @param {import('./specs').AtomSelectionSpec} sel - selection object
   * @param {unknown} colorfun - function to be used to set the color
   @example
    $3Dmol.download("pdb:4UAA",viewer,{},function(){
            viewer.setBackgroundColor(0xffffffff);
            var colorAsSnake = function(atom) {
              return atom.resi % 2 ? 'white': 'green'
            };
            viewer.setStyle( {}, { cartoon: {colorfunc: colorAsSnake }});
            viewer.render();
        });
  
   */
  setColorByFunction(sel, colorfun) {
    // eslint-disable-next-line no-var
    var atoms = this.selectedAtoms(sel, this.atoms);
    if (typeof colorfun !== 'function') return;
    this.lastColors = null; // don't bother memoizing
    if (atoms.length > 0) this.molObj = null; // force rebuild

    // now apply colorfun
    for (let i = 0; i < atoms.length; i++) {
      const a = atoms[i];
      a.color = colorfun(a);
    }
  }

  /** Convert the model into an object in the format of a ChemDoodle JSON model.
   *
   * @function GLModel#toCDObject
   * @param {boolean} [includeStyles] - whether or not to include style information. Defaults to false.
   * @return {Object}
   */
  toCDObject(includeStyles) {
    /** @type {any} */
    const out = {a: [], b: []};
    if (includeStyles) {
      out.s = [];
    }
    for (let i = 0; i < this.atoms.length; i++) {
      const atomJSON = {};
      const atom = this.atoms[i];
      atomJSON.x = atom.x;
      atomJSON.y = atom.y;
      atomJSON.z = atom.z;
      if (atom.elem != 'C') {
        atomJSON.l = atom.elem;
      }
      if (includeStyles) {
        let s = 0;
        while (s < out.s.length && JSON.stringify(atom.style) !== JSON.stringify(out.s[s])) {
          s++;
        }
        if (s === out.s.length) {
          out.s.push(atom.style);
        }
        if (s !== 0) {
          atomJSON.s = s;
        }
      }

      out.a.push(atomJSON);

      for (let b = 0; b < atom.bonds.length; b++) {
        const firstAtom = i;
        const secondAtom = atom.bonds[b];
        if (firstAtom >= secondAtom) continue;
        const bond = {
          b: firstAtom,
          e: secondAtom,
        };
        const bondOrder = atom.bondOrder[b];
        if (bondOrder != 1) {
          bond.o = bondOrder;
        }
        out.b.push(bond);
      }
    }
    return out;
  }

  /** manage the globj for this model in the possed modelGroup - if it has to be regenerated, remove and add
   *
   * @function GLModel#globj
   * @param {Object3D} group
   * @param {Object} options
   */
  globj(group, options) {
    if (this.molObj === null || options.regen) {
      // have to regenerate
      this.molObj = this.createMolObj(this.atoms, options);
      if (this.renderedMolObj) {
        // previously rendered, remove
        group.remove(this.renderedMolObj);
        this.renderedMolObj = null;
      }
      this.renderedMolObj = this.molObj.clone();
      if (this.hidden && this.renderedMolObj) {
        this.renderedMolObj.setVisible(false);
        this.molObj.setVisible(false);
      }
      group.add(this.renderedMolObj);
    }
  }

  /** return a VRML string representation of the model.  Does not include VRML header information
   * @function GLModel#exportVRML
   * @return VRML
   */
  exportVRML() {
    // todo: export spheres and cylinder objects instead of all mesh
    const tmpobj = this.createMolObj(this.atoms, {supportsImposters: false, supportsAIA: false});
    return tmpobj.vrml();
  }

  /** Remove any renderable mol object from scene
   *
   * @function GLModel#removegl
   * @param {Object3D} group
   */
  removegl(group) {
    if (this.renderedMolObj) {
      // dispose of geos and materials
      if (this.renderedMolObj.geometry !== undefined) this.renderedMolObj.geometry.dispose();
      if (this.renderedMolObj.material !== undefined) this.renderedMolObj.material.dispose();
      group.remove(this.renderedMolObj);
      this.renderedMolObj = null;
    }
    this.molObj = null;
  }

  /** @function hide
   * Don't show this model in future renderings. Keep all styles and state
   * so it can be efficiencly shown again.
   * @example
      $3Dmol.download("pdb:3ucr",viewer,{},function(){
      viewer.setStyle({},{stick:{}});
      viewer.getModel().hide();
      viewer.render();
      });
   * @function GLModel#hide
   */
  hide() {
    this.hidden = true;
    if (this.renderedMolObj) this.renderedMolObj.setVisible(false);
    if (this.molObj) this.molObj.setVisible(false);
  }

  /** @function show
   * Unhide a hidden model (see GLModel#hide)
   * @example
      $3Dmol.download("pdb:3ucr",viewer,{},function(){
      viewer.setStyle({},{stick:{}});
      viewer.getModel().hide();
      viewer.render(  )
      viewer.getModel().show()
      viewer.render();
      });
   * @function GLModel#show
   */
  show() {
    this.hidden = false;
    if (this.renderedMolObj) this.renderedMolObj.setVisible(true);
    if (this.molObj) this.molObj.setVisible(true);
  }

  /** Create labels for atoms that show the value of the passed property.
   * @function GLModel#addPropertyLabels
   *
   * @param {string} prop - property name
   * @param {import('./specs').AtomSelectionSpec} sel
   * @param {import("./GLViewer").default} viewer
   * @param {import('./specs').LabelSpec} style
   */
  addPropertyLabels(prop, sel, viewer, style) {
    // eslint-disable-next-line no-var
    var selectedAtoms = this.selectedAtoms(sel, this.atoms);
    const mystyle = deepCopy(style);
    for (let i = 0; i < selectedAtoms.length; i++) {
      const a = selectedAtoms[i];
      // eslint-disable-next-line no-nested-ternary
      const label = a[prop] ? `${a[prop]}` : a.properties[prop] ? `${a.properties[prop]}` : '';

      // this was here but typeof(<a boolean expresion>) is always the string "boolean"
      // } else if (typeof (a.properties[prop] != 'undefined')) {
      //   label = String(a.properties[prop]);
      // }

      if (label != null) {
        mystyle.position = a;
        viewer.addLabel(label, mystyle);
      }
    }
  }

  /** Create labels for residues of selected atoms.
   * Will create a single label at the center of mass of all atoms
   * with the same chain,resn, and resi.
   * @function GLModel#addResLabels
   *
   * @param {import('./specs').AtomSelectionSpec} sel
   * @param {import("./GLViewer").default} viewer
   * @param {import('./specs').LabelSpec} style
   * @param {boolean} byframe - if true, create labels for every individual frame, not just current; frames must be loaded already
   */
  addResLabels(sel, viewer, style, byframe) {
    /** @type {Array<import("./Label").default>} */
    const createdLabels = [];

    /**
     *
     * @param {GLModel} model
     * @param {number} [framenum]
     */
    const helper = (model, framenum) => {
      const atoms = model.selectedAtoms(sel, this.atoms);
      /** @type {Record<string, Record<string, Array<import('./specs').AtomSpec>>>} */
      const bylabel = {};
      // collect by chain:resn:resi
      for (const a of atoms) {
        // added default to avoid undefined
        const c = a.chain || '';
        const {resn} = a;
        const {resi} = a;
        const label = `${resn}${resi}`;
        if (!bylabel[c]) bylabel[c] = {};
        if (!bylabel[c][label]) bylabel[c][label] = [];
        bylabel[c][label].push(a);
      }

      const mystyle = deepCopy(style);
      // now compute centers of mass
      for (const labels of Object.values(bylabel)) {
        for (const [label, labeledAtoms] of Object.entries(labels)) {
          const sum = new Vector3(0, 0, 0);
          for (let i = 0; i < labeledAtoms.length; i++) {
            const a = labeledAtoms[i];
            sum.x += a.x;
            sum.y += a.y;
            sum.z += a.z;
          }
          sum.divideScalar(labeledAtoms.length);
          mystyle.position = sum;
          mystyle.frame = framenum;
          const l = viewer.addLabel(label, mystyle, undefined, true);
          createdLabels.push(l);
        }
      }
    };

    if (byframe) {
      const n = this.getNumFrames();
      const savedatoms = this.atoms;
      for (let i = 0; i < n; i++) {
        if (this.frames[i]) {
          this.atoms = this.frames[i];
          helper(this, i);
        }
      }
      this.atoms = savedatoms;
    } else {
      helper(this);
    }
    return createdLabels;
  }

  /**
   * Set coordinates from remote trajectory file.
   * @function GLModel#setCoordinatesFromURL
   * @param {string} url - contains the url where mdsrv has been hosted
   * @param {string} path - contains the path of the file (<root>/filename)
   * @return {Promise}
   */
  setCoordinatesFromURL(url, path) {
    // @ts-ignore
    this.frames = [];
    const self = this;
    if (this.box) this.setupDFS();

    return new Promise((resolve, reject) => {
      if (!url.startsWith('http')) url = `http://${url}`;
      $.get(`${url}/traj/numframes/${path}`, numFrames => {
        // not sure if this needs to parse 0x as hex
        // eslint-disable-next-line radix
        if (!Number.isNaN(Number.parseInt(numFrames))) {
          self.frames.push(this.atoms);
          self.frames.numFrames = numFrames;
          self.frames.url = url;
          self.frames.path = path;
          self
            .setFrame(0)
            .then(() => {
              resolve(undefined);
            })
            .catch(reject);
        }
      });
    });
  }

  /**
   * Set coordinates for the atoms from provided trajectory file.
   * @function GLModel#setCoordinates
   * @param {string} str - contains the data of the file
   * @param {string} format - contains the format of the file (mdcrd, inpcrd, pdb, netcdf, or array).  Arrays should be TxNx3 where T is the number of timesteps and N the number of atoms.
     @example
        let m = viewer.addModel()  //create an empty model
        m.addAtoms([{x:0,y:0,z:0,elem:'C'},{x:2,y:0,z:0,elem:'C'}]) //provide a list of dictionaries representing the atoms
        viewer.setStyle({'sphere':{}})
        m.setCoordinates([[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [2.8888888359069824, 0.0, 0.0]], [[0.0, 0.0, 0.0], [3.777777671813965, 0.0, 0.0]], [[0.0, 0.0, 0.0], [4.666666507720947, 0.0, 0.0]], [[0.0, 0.0, 0.0], [5.55555534362793, 0.0, 0.0]], [[0.0, 0.0, 0.0], [6.44444465637207, 0.0, 0.0]], [[0.0, 0.0, 0.0], [7.333333492279053, 0.0, 0.0]], [[0.0, 0.0, 0.0], [8.222222328186035, 0.0, 0.0]], [[0.0, 0.0, 0.0], [9.11111068725586, 0.0, 0.0]], [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]],'array');
        viewer.animate({loop: "forward",reps: 1});
        viewer.zoomTo();
        viewer.zoom(0.5);
        viewer.render();
   */
  setCoordinates(str, format) {
    format = format || '';
    if (!str) return []; // leave an empty model

    if (/\.gz$/.test(format)) {
      // unzip gzipped files
      format = format.replace(/\.gz$/, '');
      try {
        str = pako.inflate(str, {
          to: 'string',
        });
      } catch (err) {
        console.log(err);
      }
    }
    const supportedFormats = {mdcrd: '', inpcrd: '', pdb: '', netcdf: '', array: ''};
    if (supportedFormats[format]) {
      this.frames = [];
      const atomCount = this.atoms.length;
      const values = GLModel.parseCrd(str, format);
      let count = 0;
      while (count < values.length) {
        const temp = [];
        for (let i = 0; i < atomCount; i++) {
          const newAtom = {};
          for (const k in this.atoms[i]) {
            newAtom[k] = this.atoms[i][k];
          }
          temp[i] = newAtom;
          temp[i].x = values[count++];
          temp[i].y = values[count++];
          temp[i].z = values[count++];
        }

        this.frames.push(temp);
      }
      this.atoms = this.frames[0];
      return this.frames;
    }
    return [];
  }

  /** given a selection specification, return true if atom is selected.
   * Does not support context-aware selectors like expand/within/byres.
   *
   * @function GLModel#atomIsSelected
   * @param {import('./specs').AtomSpec} atom
   * @param {import('./specs').AtomSelectionSpec} sel
   * @return {boolean}
   */
  atomIsSelected(atom, sel) {
    if (typeof sel === 'undefined') return true; // undef gets all
    const invert = !!sel.invert;
    let ret = true;
    for (const key in sel) {
      if (key == 'and' || key == 'or' || key == 'not') {
        // boolean operators
        if (key == 'not') {
          // @ts-ignore
          if (this.atomIsSelected(atom, sel[key])) {
            ret = false;
            break;
          }
        } else {
          // "or" and "and"
          // these selections are expensive so when called via
          // selectedAtoms shoudl be cached - but if atomIsSelected
          // is called directly create the cache
          // @ts-ignore
          if (sel[key].__cached_results === undefined) {
            sel = deepCopyAndCache(sel, this);
          }

          // @ts-ignore
          ret = sel[key].__cached_results.has(atom.index);
          if (!ret) {
            break;
          }
        }
      } else if (key === 'predicate') {
        // a user supplied function for evaluating atoms
        if (sel.predicate && !sel.predicate(atom)) {
          ret = false;
          break;
        }
      } else if (key == 'properties' && atom[key]) {
        // @ts-ignore
        for (const propkey in sel.properties) {
          if (propkey.startsWith('__cache')) continue;
          if (typeof atom.properties[propkey] === 'undefined') {
            ret = false;
            break;
          }
          // @ts-ignore
          if (atom.properties[propkey] != sel.properties[propkey]) {
            ret = false;
            break;
          }
        }
      } else if (
        sel[key] &&
        key != 'props' &&
        key != 'invert' &&
        key != 'model' &&
        key != 'byres' &&
        key != 'expand' &&
        key != 'within' &&
        key != 'and' &&
        key != 'or' &&
        key != 'not' &&
        !key.startsWith('__cache')
      ) {
        // if something is in sel, atom must have it
        if (typeof atom[key] === 'undefined') {
          ret = false;
          break;
        }
        let isokay = false;
        if (key === 'bonds') {
          // special case counting number of bonds, for selecting nonbonded mostly
          const val = sel[key];
          if (val != atom.bonds.length) {
            ret = false;
            break;
          }
        } else if (Array.isArray(sel[key])) {
          // can be any of the listed values
          const valarr = sel[key];
          const atomval = atom[key];
          for (let i = 0; i < valarr.length; i++) {
            if (propertyMatches(atomval, valarr[i])) {
              isokay = true;
              break;
            }
          }
          if (!isokay) {
            ret = false;
            break;
          }
        } else {
          // single match
          const val = sel[key];
          if (!propertyMatches(atom[key], val)) {
            ret = false;
            break;
          }
        }
      }
    }

    return invert ? !ret : ret;
  }

  // private functions
  // go through all the atoms and regenerate their geometries
  // we try to have one geometry for each style since this is much much
  // faster
  // at some point we should optimize this to avoid unnecessary
  // recalculation
  /** param {AtomSpec[]} atoms */
  createMolObj(atoms, options) {
    options = options || {};

    const ret = new Object3D();
    const cartoonAtoms = [];
    const lineGeometries = {};
    const crossGeometries = {};

    let drawSphereFunc = this.drawAtomSphere;
    let sphereGeometry = null;
    let stickGeometry = null;
    if (options.supportsImposters) {
      drawSphereFunc = this.drawAtomImposter.bind(this);
      sphereGeometry = new Geometry(true);
      sphereGeometry.imposter = true;
      stickGeometry = new Geometry(true, true);
      stickGeometry.imposter = true;
      stickGeometry.sphereGeometry = new Geometry(true); // for caps
      stickGeometry.sphereGeometry.imposter = true;
      stickGeometry.drawnCaps = {};
    } else if (options.supportsAIA) {
      drawSphereFunc = this.drawAtomInstanced.bind(this);
      sphereGeometry = new Geometry(false, true, true);
      sphereGeometry.instanced = true;
      stickGeometry = new Geometry(true); // don't actually have instanced sticks
    } else {
      sphereGeometry = new Geometry(true);
      stickGeometry = new Geometry(true);
    }

    let i;
    let j;
    let n;
    let testOpacities;
    const opacities = {};
    const range = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
    for (i = 0, n = atoms.length; i < n; i++) {
      const atom = atoms[i];
      // recreate gl info for each atom as necessary
      // set up appropriate intersection spheres for clickable atoms
      if (atom && atom.style) {
        if ((atom.clickable || atom.hoverable) && atom.intersectionShape === undefined)
          atom.intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};

        testOpacities = {line: undefined, cross: undefined, stick: undefined, sphere: undefined};
        for (j in testOpacities) {
          if (atom.style[j]) {
            if (atom.style[j].opacity) testOpacities[j] = parseFloat(atom.style[j].opacity);
            else testOpacities[j] = 1;
          } else testOpacities[j] = undefined;

          if (opacities[j]) {
            if (testOpacities[j] != undefined && opacities[j] != testOpacities[j]) {
              console.log(`Warning: ${j} opacity is ambiguous`);
              opacities[j] = 1;
            }
          } else opacities[j] = testOpacities[j];
        }

        drawSphereFunc(atom, sphereGeometry);
        this.drawAtomClickSphere(atom);
        // @ts-ignore
        this.drawAtomCross(atom, crossGeometries);
        // @ts-ignore
        this.drawBondLines(atom, atoms, lineGeometries);
        this.drawBondSticks(atom, atoms, stickGeometry);

        if (typeof atom.style.cartoon !== 'undefined' && !atom.style.cartoon.hidden) {
          // gradient color scheme range
          if (
            atom.style.cartoon.color === 'spectrum' &&
            typeof atom.resi === 'number' &&
            !atom.hetflag
          ) {
            if (atom.resi < range[0]) range[0] = atom.resi;
            if (atom.resi > range[1]) range[1] = atom.resi;
          }

          cartoonAtoms.push(atom);
        }
      }
    }
    // create cartoon if needed - this is a whole model analysis
    if (cartoonAtoms.length > 0) {
      drawCartoon(ret, cartoonAtoms, range, this.defaultCartoonQuality);
    }

    // add sphere geometry
    if (sphereGeometry && sphereGeometry.vertices > 0) {
      // Initialize buffers in geometry
      sphereGeometry.initTypedArrays();
      let sphereMaterial = null;
      let sphere = null;

      // create appropriate material
      if (sphereGeometry.imposter) {
        sphereMaterial = new SphereImposterMaterial({
          ambient: 0x000000,
          vertexColors: true,
          reflectivity: 0,
        });
      } else if (sphereGeometry.instanced) {
        sphere = new Geometry(true);
        GLDraw.drawSphere(sphere, {x: 0, y: 0, z: 0}, 1, new Color(0.5, 0.5, 0.5));
        sphere.initTypedArrays();
        sphereMaterial = new InstancedMaterial({
          sphereMaterial: new MeshLambertMaterial({
            ambient: 0x000000,
            vertexColors: true,
            reflectivity: 0,
          }),
          sphere,
        });
      } else {
        // regular mesh
        sphereMaterial = new MeshLambertMaterial({
          ambient: 0x000000,
          vertexColors: true,
          reflectivity: 0,
        });
      }
      if (opacities.sphere < 1 && opacities.sphere >= 0) {
        sphereMaterial.transparent = true;
        sphereMaterial.opacity = opacities.sphere;
      }

      sphere = new Mesh(sphereGeometry, sphereMaterial);
      ret.add(sphere);
    }

    // add stick geometry
    if (stickGeometry.vertices > 0) {
      let stickMaterial = null;
      let ballMaterial = null;
      let balls = stickGeometry.sphereGeometry;
      if (!balls || typeof balls.vertices === 'undefined' || balls.vertices == 0) balls = null; // no balls

      // Initialize buffers in geometry
      stickGeometry.initTypedArrays();
      if (balls) balls.initTypedArrays();

      // create material
      const matvals = {ambient: 0x000000, vertexColors: true, reflectivity: 0};

      if (stickGeometry.imposter) {
        stickMaterial = new StickImposterMaterial(matvals);
        ballMaterial = new SphereImposterMaterial(matvals);
      } else {
        stickMaterial = new MeshLambertMaterial(matvals);
        ballMaterial = new MeshLambertMaterial(matvals);

        if (stickMaterial.wireframe) {
          stickGeometry.setUpWireframe();
          if (balls) balls.setUpWireframe();
        }
      }

      if (opacities.stick < 1 && opacities.stick >= 0) {
        stickMaterial.transparent = true;
        stickMaterial.opacity = opacities.stick;
        ballMaterial.transparent = true;
        ballMaterial.opacity = opacities.stick;
      }
      const sticks = new Mesh(stickGeometry, stickMaterial);
      ret.add(sticks);

      if (balls) {
        const stickspheres = new Mesh(balls, ballMaterial);
        ret.add(stickspheres);
      }
    }

    // var linewidth;
    // add any line geometries, distinguished by line width
    let linewidth;
    for (i in lineGeometries) {
      if (lineGeometries[i]) {
        linewidth = i;
        const lineMaterial = new LineBasicMaterial({
          linewidth,
          vertexColors: true,
        });
        if (opacities.line < 1 && opacities.line >= 0) {
          lineMaterial.transparent = true;
          lineMaterial.opacity = opacities.line;
        }

        lineGeometries[i].initTypedArrays();

        const line = new Line(lineGeometries[i], lineMaterial, LinePieces);

        ret.add(line);
      }
    }

    // add any cross geometries
    for (i in crossGeometries) {
      if (crossGeometries[i]) {
        linewidth = i;
        const crossMaterial = new LineBasicMaterial({
          linewidth,
          vertexColors: true,
        });
        if (opacities.cross < 1 && opacities.cross >= 0) {
          crossMaterial.transparent = true;
          crossMaterial.opacity = opacities.cross;
        }

        crossGeometries[i].initTypedArrays();

        const cross = new Line(crossGeometries[i], crossMaterial, LinePieces);

        ret.add(cross);
      }
    }

    // for BIOMT assembly
    if (
      this.dontDuplicateAtoms &&
      this.modelData.symmetries &&
      this.modelData.symmetries.length > 0
    ) {
      const finalRet = new Object3D();
      let t;
      for (t = 0; t < this.modelData.symmetries.length; t++) {
        let transformedRet = new Object3D();
        transformedRet = ret.clone();
        transformedRet.matrix.copy(this.modelData.symmetries[t]);
        transformedRet.matrixAutoUpdate = false;
        finalRet.add(transformedRet);
      }
      return finalRet;
    }

    return ret;
  }

  // from atom, return a normalized vector v that is orthogonal and along which
  // it is appropraite to draw multiple bonds
  getSideBondV(atom, atom2, i) {
    let i2;
    let j2;
    let atom3;
    let p3;
    let dir2;
    const self = this;
    const p1 = new Vector3(atom.x, atom.y, atom.z);
    const p2 = new Vector3(atom2.x, atom2.y, atom2.z);
    const dir = p2.clone();
    let v = null;
    dir.sub(p1);

    /**
     * 
     * @param {import('./specs').AtomSpec} crossAtom1 
     * @param {import('./specs').AtomSpec} crossAtom2 
     * @returns 
     */
    function getGoodCross(crossAtom1, crossAtom2) {
      // get vector 2 different neighboring atom
      // find most divergent neighbor
      let bestv = null;
      let bestlen = -1;
      for (let j = 0, n = crossAtom1.bonds.length; j < n; j++) {
        if (crossAtom1.bonds[j] != crossAtom2.index) {
          j2 = crossAtom1.bonds[j];
          atom3 = self.atoms[j2];
          p3 = new Vector3(atom3.x, atom3.y, atom3.z);

          dir2 = p3.clone();
          dir2.sub(p1);

          v = dir2.clone();
          v.cross(dir);
          const l = v.lengthSq();
          if (l > bestlen) {
            bestlen = l;
            bestv = v;
            if (bestlen > 0.1) {
              return bestv;
            }
          }
        }
      }
      return bestv;
    }

    if (atom.bonds.length === 1) {
      if (atom2.bonds.length === 1) {
        v = dir.clone();
        if (Math.abs(v.x) > 0.0001) v.y += 1;
        else v.x += 1;
      } else {
        i2 = (i + 1) % atom2.bonds.length;
        j2 = atom2.bonds[i2];
        atom3 = this.atoms[j2];
        p3 = new Vector3(atom3.x, atom3.y, atom3.z);

        dir2 = p3.clone();
        dir2.sub(p1);

        v = dir2.clone();
        v.cross(dir);
      }
    } else {
      v = getGoodCross(atom, atom2);
      if (!v) return null;
      if (v.lengthSq() < 0.01) {
        const v2 = getGoodCross(atom2, atom);
        if (v2 != null) v = v2; // can be null if no other neighbors
      }
    }

    // especially for C#C (triple bond) dir and dir2
    // may be opposites resulting in a zero v
    if (v.lengthSq() < 0.01) {
      v = dir.clone();
      if (Math.abs(v.x) > 0.0001) v.y += 1;
      else v.x += 1;
    }

    v.cross(dir);
    v.normalize();

    return v;

    // v.multiplyScalar(r * 1.5);
  }

  // recurse over the current atoms to establish a depth first order
  setupDFS() {
    this.atomdfs = [];
    const self = this;

    const visited = new Int8Array(this.atoms.length);
    visited.fill(0);

    function search(i, prev, component) {
      // add i to component and recursive explore connected atoms
      component.push([i, prev]);
      const atom = self.atoms[i];
      visited[i] = 1;
      for (let b = 0; b < atom.bonds.length; b++) {
        const nexti = atom.bonds[b];
        if (self.atoms[nexti] && !visited[nexti]) {
          search(nexti, i, component);
        }
      }
    };

    for (let i = 0; i < this.atoms.length; i++) {
      const atom = this.atoms[i];
      if (atom && !visited[i]) {
        const component = [];
        search(i, -1, component);
        this.atomdfs.push(component);
      }
    }
  }

  // go over current atoms in depth first order and ensure that connected
  // attoms aren't split across the box
  adjustCoordinatesToBox() {
    if (!this.box) return;
    if (!this.atomdfs) return;
    const bx = this.box[0];
    const by = this.box[1];
    const bz = this.box[2];
    const mx = bx * 0.9;
    const my = by * 0.9;
    const mz = bz * 0.9;

    for (let c = 0; c < this.atomdfs.length; c++) {
      // for each connected component
      const component = this.atomdfs[c];
      for (let i = 1; i < component.length; i++) {
        // compare each atom to its previous and prevent coordinates from wrapping
        const atom = this.atoms[component[i][0]];
        const prev = this.atoms[component[i][1]];
        atom.x = adjustCoord(prev.x, atom.x, mx, bx);
        atom.y = adjustCoord(prev.y, atom.y, my, by);
        atom.z = adjustCoord(prev.z, atom.z, mz, bz);
      }
    }
  }

  // return proper radius for atom given style
  /**
   *
   * @param {import('./specs').AtomSpec} atom
   * @param {import('./specs').AtomStyleSpec} style
   * @return {number}
   *
   */
  getRadiusFromStyle(atom, style) {
    let r = this.defaultSphereRadius;
    if (typeof style.radius != 'undefined') r = style.radius;
    else if (vdwRadii[atom.elem]) r = vdwRadii[atom.elem];
    else if (atom.elem.length > 1) {
      // see if adjusting case helps
      let e = atom.elem;
      e = `${e[1].toLowerCase()}${e.substring(1)}`;
      if (vdwRadii[e]) r = vdwRadii[e];
    }

    if (typeof style.scale != 'undefined') r *= style.scale;
    return r;
  }

  /** Register atom shaped click handlers */
  drawAtomClickSphere(atom) {
    if (!atom.style.clicksphere) return;
    const style = atom.style.clicksphere;
    if (style.hidden) return;

    const radius = this.getRadiusFromStyle(atom, style);

    if ((atom.clickable === true || atom.hoverable) && atom.intersectionShape !== undefined) {
      const center = new Vector3(atom.x, atom.y, atom.z);
      atom.intersectionShape.sphere.push(new Sphere(center, radius));
    }
  }

  // sphere drawing
  // See also: drawCylinder
  /**
   *
   * @param {import('./specs').AtomSpec} atom
   * @param {Geometry} geo
   */
  drawAtomSphere(atom, geo) {
    if (!atom.style.sphere) return;
    const style = atom.style.sphere;
    if (style.hidden) return;

    const C = getColorFromStyle(atom, style);

    const radius = this.getRadiusFromStyle(atom, style);

    if ((atom.clickable === true || atom.hoverable) && atom.intersectionShape !== undefined) {
      const center = new Vector3(atom.x, atom.y, atom.z);
      atom.intersectionShape.sphere.push(new Sphere(center, radius));
    }

    GLDraw.drawSphere(geo, atom, radius, C);
  }

  // cross drawing
  /** @typedef CrossStyleSpec
   * @prop {boolean} hidden - do not show
   * @prop {number} linewidth *deprecated due to vanishing browser support*
   * @prop {number} radius
   * @prop {number} scale - scale radius by specified amount
   * @prop {import('./specs').ColorSchemeSpec} colorscheme - element based coloring
   * @prop {import('./specs').ColorSpec} color - fixed coloring, overrides colorscheme
   * @prop {number} opacity - opacity, must be the same for all atoms in the model
   */
  /**
   *
   * @param {import('./specs').AtomSpec} atom
   * @param {Geometry[]} geos
   */
  drawAtomCross(atom, geos) {
    if (!atom.style.cross) return;
    const style = atom.style.cross;
    if (style.hidden) return;
    const linewidth = style.linewidth || defaultlineWidth;
    if (!geos[linewidth]) geos[linewidth] = new Geometry();

    const geoGroup = geos[linewidth].updateGeoGroup(6);

    const delta = this.getRadiusFromStyle(atom, style);

    const points = [
      [delta, 0, 0],
      [-delta, 0, 0],
      [0, delta, 0],
      [0, -delta, 0],
      [0, 0, delta],
      [0, 0, -delta],
    ];

    const clickable = atom.clickable || atom.hoverable;
    if (clickable && atom.intersectionShape === undefined)
      atom.intersectionShape = {sphere: [], cylinder: [], line: []};

    const c = getColorFromStyle(atom, style);

    const {vertexArray} = geoGroup;
    const {colorArray} = geoGroup;

    for (let j = 0; j < 6; j++) {
      const offset = geoGroup.vertices * 3;

      geoGroup.vertices++;
      vertexArray[offset] = atom.x + points[j][0];
      vertexArray[offset + 1] = atom.y + points[j][1];
      vertexArray[offset + 2] = atom.z + points[j][2];
      colorArray[offset] = c.r;
      colorArray[offset + 1] = c.g;
      colorArray[offset + 2] = c.b;

      if (clickable) {
        const point = new Vector3(points[j][0], points[j][1], points[j][2]);

        // decrease cross size for selection to prevent misselection from atom overlap
        point.multiplyScalar(0.1);
        point.set(point.x + atom.x, point.y + atom.y, point.z + atom.z);
        atom.intersectionShape.line.push(point);
      }
    }
  }

  /** @typedef LineStyleSpec
   * @prop {boolean} hidden - do not show line
   * @prop {number} linewidth *deprecated due to vanishing browser support*
   * @prop {import('./specs').ColorSchemeSpec} colorscheme - element based coloring
   * @prop {import('./specs').ColorSpec} color - fixed coloring, overrides colorscheme
   * @prop {number} opacity - opacity, must be the same for all atoms in the model
   */
  // bonds - both atoms must match bond style
  // standardize on only drawing for lowest to highest
  /**
   *
   * @param {import('./specs').AtomSpec}
   *            atom
   * @param {import('./specs').AtomSpec[]} atoms
   * @param {Geometry[]} geos
   */
  drawBondLines(atom, atoms, geos) {
    if (!atom.style.line) return;
    const style = atom.style.line;
    if (style.hidden) return;
    let p1a;
    let p1b;
    let p2a;
    let p2b;
    // have a separate geometry for each linewidth
    const linewidth = style.linewidth || defaultlineWidth;

    if (!geos[linewidth]) geos[linewidth] = new Geometry();
    /** @type {import("./WebGL/core").geometryGroup} */
    const geoGroup = geos[linewidth].updateGeoGroup(6 * atom.bonds.length); // reserve enough space even for triple bonds

    const {vertexArray} = geoGroup;
    const {colorArray} = geoGroup;

    for (let i = 0; i < atom.bonds.length; i++) {
      const j = atom.bonds[i]; // our neighbor

      const atom2 = atoms[j];
      if (!atom2.style.line) continue; // don't sweat the details

      if (atom.index >= atom2.index)
        // only draw if less, this way we can do multi bonds correctly
        continue;
      const p1 = new Vector3(atom.x, atom.y, atom.z);
      const p2 = new Vector3(atom2.x, atom2.y, atom2.z);
      const mp = p1.clone().add(p2).multiplyScalar(0.5);
      let singleBond = false;

      const atomneedsi = atom.clickable || atom.hoverable;
      const atom2needsi = atom2.clickable || atom2.hoverable;

      if (atomneedsi || atom2needsi) {
        if (atomneedsi) {
          if (atom.intersectionShape === undefined)
            atom.intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
          atom.intersectionShape.line.push(p1);
          atom.intersectionShape.line.push(mp);
        }
        if (atom2needsi) {
          if (atom2.intersectionShape === undefined)
            atom2.intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
          atom2.intersectionShape.line.push(mp);
          atom2.intersectionShape.line.push(p2);
        }
      }
      let c1 = getColorFromStyle(atom, atom.style.line);
      let c2 = getColorFromStyle(atom2, atom2.style.line);

      if (atom.bondStyles && atom.bondStyles[i]) {
        const bstyle = atom.bondStyles[i];
        if (!bstyle.iswire) {
          continue;
        }
        if (bstyle.singleBond) singleBond = true;
        if (typeof bstyle.color1 != 'undefined') {
          c1 = CC.color(bstyle.color1);
        }
        if (typeof bstyle.color2 != 'undefined') {
          c2 = CC.color(bstyle.color2);
        }
      }

      const offset = geoGroup.vertices * 3;
      let mpa;
      let mpb;

      if (atom.bondOrder[i] > 1 && atom.bondOrder[i] < 4 && !singleBond) {
        const v = this.getSideBondV(atom, atom2, i);
        const dir = p2.clone();
        dir.sub(p1);

        if (atom.bondOrder[i] == 2 && v) {
          // double
          v.multiplyScalar(0.1);
          p1a = p1.clone();
          p1a.add(v);
          p1b = p1.clone();
          p1b.sub(v);

          p2a = p1a.clone();
          p2a.add(dir);
          p2b = p1b.clone();
          p2b.add(dir);

          if (c1 == c2) {
            geoGroup.vertices += 4;
            addLine(vertexArray, colorArray, offset, p1a, p2a, c1);
            addLine(vertexArray, colorArray, offset + 6, p1b, p2b, c1);
          } else {
            geoGroup.vertices += 8;
            dir.multiplyScalar(0.5);
            mpa = p1a.clone();
            mpa.add(dir);
            mpb = p1b.clone();
            mpb.add(dir);

            addLine(vertexArray, colorArray, offset, p1a, mpa, c1);
            addLine(vertexArray, colorArray, offset + 6, mpa, p2a, c2);
            addLine(vertexArray, colorArray, offset + 12, p1b, mpb, c1);
            addLine(vertexArray, colorArray, offset + 18, mpb, p2b, c2);
          }
        } else if (atom.bondOrder[i] == 3 && v) {
          // triple
          v.multiplyScalar(0.1);
          p1a = p1.clone();
          p1a.add(v);
          p1b = p1.clone();
          p1b.sub(v);

          p2a = p1a.clone();
          p2a.add(dir);
          p2b = p1b.clone();
          p2b.add(dir);

          if (c1 == c2) {
            geoGroup.vertices += 6;
            addLine(vertexArray, colorArray, offset, p1, p2, c1);
            addLine(vertexArray, colorArray, offset + 6, p1a, p2a, c1);
            addLine(vertexArray, colorArray, offset + 12, p1b, p2b, c1);
          } else {
            geoGroup.vertices += 12;
            dir.multiplyScalar(0.5);
            mpa = p1a.clone();
            mpa.add(dir);
            mpb = p1b.clone();
            mpb.add(dir);

            addLine(vertexArray, colorArray, offset, p1, mp, c1);
            addLine(vertexArray, colorArray, offset + 6, mp, p2, c2);
            addLine(vertexArray, colorArray, offset + 12, p1a, mpa, c1);
            addLine(vertexArray, colorArray, offset + 18, mpa, p2a, c2);
            addLine(vertexArray, colorArray, offset + 24, p1b, mpb, c1);
            addLine(vertexArray, colorArray, offset + 30, mpb, p2b, c2);
          }
        }
      } else {
        // single bond
        if (c1 == c2) {
          geoGroup.vertices += 2;
          addLine(vertexArray, colorArray, offset, p1, p2, c1);
        } else {
          geoGroup.vertices += 4;
          addLine(vertexArray, colorArray, offset, p1, mp, c1);
          addLine(vertexArray, colorArray, offset + 6, mp, p2, c2);
        }
      }
    }
  }

  /** @typedef StickStyleSpec
   * @prop {boolean} hidden - do not show
   * @prop {number} radius
   * @prop {boolean} singleBonds - draw all bonds as single bonds if set
   * @prop {import('./specs').ColorSchemeSpec} colorscheme - element based coloring
   * @prop {import('./specs').ColorSpec} color - fixed coloring, overrides colorscheme
   * @prop {number} opacity - opacity, must be the same for all atoms in the model
   */
  /**
   * draws cylinders and small spheres (at bond radius)
   * @param {any} atom
   * @param {any} atoms
   * @param {any} geo
   * @returns
   */
  drawBondSticks(atom, atoms, geo) {
    if (!atom.style.stick) return;
    const style = atom.style.stick;
    if (style.hidden) return;

    const atomBondR = style.radius || defaultStickRadius;
    let bondR = atomBondR;
    const atomSingleBond = style.singleBonds || false;
    let fromCap = 0;
    let toCap = 0;
    let atomneedsi;
    let atom2needsi;
    let i;
    let singleBond;
    let bstyle;
    let cylinder1a;
    let cylinder1b;
    let cylinder1c;
    let cylinder2a;
    let cylinder2b;
    let cylinder2c;

    let C1 = getColorFromStyle(atom, style);

    let mp;
    let mp2;
    let mp3;

    if (!atom.capDrawn && atom.bonds.length < 4) fromCap = 2;

    let drawCyl = GLDraw.drawCylinder; // mesh cylinder
    if (geo.imposter) drawCyl = drawStickImposter;

    for (i = 0; i < atom.bonds.length; i++) {
      const j = atom.bonds[i]; // our neighbor
      const atom2 = atoms[j]; // parsePDB, etc should only add defined bonds
      mp = mp2 = mp3 = null;
      if (atom.index < atom2.index) {
        // only draw if less, this
        // lets us combine
        // cylinders of the same
        // color
        const style2 = atom2.style;
        if (!style2.stick || style2.stick.hidden) continue; // don't sweat the details

        let C2 = getColorFromStyle(atom2, style2.stick);

        // support bond specific styles
        bondR = atomBondR;
        singleBond = atomSingleBond;
        if (atom.bondStyles && atom.bondStyles[i]) {
          bstyle = atom.bondStyles[i];
          if (bstyle.iswire) {
            continue;
          }
          if (bstyle.radius) bondR = bstyle.radius;
          if (bstyle.singleBond) singleBond = true;
          if (typeof bstyle.color1 != 'undefined') {
            C1 = CC.color(bstyle.color1);
          }
          if (typeof bstyle.color2 != 'undefined') {
            C2 = CC.color(bstyle.color2);
          }
        }
        const p1 = new Vector3(atom.x, atom.y, atom.z);
        const p2 = new Vector3(atom2.x, atom2.y, atom2.z);

        // draw cylinders
        if (atom.bondOrder[i] === 1 || singleBond || atom.bondOrder[i] > 3) {
          // TODO: aromatics at 4
          if (!atom2.capDrawn && atom2.bonds.length < 4) toCap = 2;

          if (C1 != C2) {
            mp = new Vector3().addVectors(p1, p2).multiplyScalar(0.5);
            drawCyl(geo, p1, mp, bondR, C1, fromCap, 0);
            drawCyl(geo, mp, p2, bondR, C2, 0, toCap);
          } else {
            drawCyl(geo, p1, p2, bondR, C1, fromCap, toCap);
          }

          atomneedsi = atom.clickable || atom.hoverable;
          atom2needsi = atom2.clickable || atom2.hoverable;

          if (atomneedsi || atom2needsi) {
            if (!mp) mp = new Vector3().addVectors(p1, p2).multiplyScalar(0.5);
            if (atomneedsi) {
              const cylinder1 = new Cylinder(p1, mp, bondR);
              const sphere1 = new Sphere(p1, bondR);
              atom.intersectionShape.cylinder.push(cylinder1);
              atom.intersectionShape.sphere.push(sphere1);
            }
            if (atom2needsi) {
              const cylinder2 = new Cylinder(p2, mp, bondR);
              const sphere2 = new Sphere(p2, bondR);
              atom2.intersectionShape.cylinder.push(cylinder2);
              atom2.intersectionShape.sphere.push(sphere2);
            }
          }
        } else if (atom.bondOrder[i] > 1) {
          // multi bond caps
          let mfromCap = 0;
          let mtoCap = 0;

          if (bondR != atomBondR) {
            // assume jmol style multiple bonds - the radius doesn't fit within atom sphere
            mfromCap = 2;
            mtoCap = 2;
          }

          const dir = p2.clone();
          let v = null;
          dir.sub(p1);

          let r;
          let p1a;
          let p1b;
          let p2a;
          let p2b;
          v = this.getSideBondV(atom, atom2, i);

          if (atom.bondOrder[i] == 2 && v) {
            r = bondR / 2.5;

            v.multiplyScalar(r * 1.5);
            p1a = p1.clone();
            p1a.add(v);
            p1b = p1.clone();
            p1b.sub(v);

            p2a = p1a.clone();
            p2a.add(dir);
            p2b = p1b.clone();
            p2b.add(dir);

            if (C1 != C2) {
              mp = new Vector3().addVectors(p1a, p2a).multiplyScalar(0.5);
              mp2 = new Vector3().addVectors(p1b, p2b).multiplyScalar(0.5);
              drawCyl(geo, p1a, mp, r, C1, mfromCap, 0);
              drawCyl(geo, mp, p2a, r, C2, 0, mtoCap);
              drawCyl(geo, p1b, mp2, r, C1, mfromCap, 0);
              drawCyl(geo, mp2, p2b, r, C2, 0, mtoCap);
            } else {
              drawCyl(geo, p1a, p2a, r, C1, mfromCap, mtoCap);
              drawCyl(geo, p1b, p2b, r, C1, mfromCap, mtoCap);
            }

            atomneedsi = atom.clickable || atom.hoverable;
            atom2needsi = atom2.clickable || atom2.hoverable;

            if (atomneedsi || atom2needsi) {
              if (!mp) mp = new Vector3().addVectors(p1a, p2a).multiplyScalar(0.5);
              if (!mp2) mp2 = new Vector3().addVectors(p1b, p2b).multiplyScalar(0.5);
              if (atomneedsi) {
                cylinder1a = new Cylinder(p1a, mp, r);
                cylinder1b = new Cylinder(p1b, mp2, r);
                atom.intersectionShape.cylinder.push(cylinder1a);
                atom.intersectionShape.cylinder.push(cylinder1b);
              }
              if (atom2needsi) {
                cylinder2a = new Cylinder(p2a, mp, r);
                cylinder2b = new Cylinder(p2b, mp2, r);
                atom2.intersectionShape.cylinder.push(cylinder2a);
                atom2.intersectionShape.cylinder.push(cylinder2b);
              }
            }
          } else if (atom.bondOrder[i] == 3 && v) {
            r = bondR / 4;
            v.cross(dir);
            v.normalize();
            v.multiplyScalar(r * 3);

            p1a = p1.clone();
            p1a.add(v);
            p1b = p1.clone();
            p1b.sub(v);

            p2a = p1a.clone();
            p2a.add(dir);
            p2b = p1b.clone();
            p2b.add(dir);

            if (C1 != C2) {
              mp = new Vector3().addVectors(p1a, p2a).multiplyScalar(0.5);
              mp2 = new Vector3().addVectors(p1b, p2b).multiplyScalar(0.5);
              mp3 = new Vector3().addVectors(p1, p2).multiplyScalar(0.5);
              drawCyl(geo, p1a, mp, r, C1, mfromCap, 0);
              drawCyl(geo, mp, p2a, r, C2, 0, mtoCap);
              drawCyl(geo, p1, mp3, r, C1, fromCap, 0);
              drawCyl(geo, mp3, p2, r, C2, 0, toCap);
              drawCyl(geo, p1b, mp2, r, C1, mfromCap, 0);
              drawCyl(geo, mp2, p2b, r, C2, 0, mtoCap);
            } else {
              drawCyl(geo, p1a, p2a, r, C1, mfromCap, mtoCap);
              drawCyl(geo, p1, p2, r, C1, fromCap, toCap);
              drawCyl(geo, p1b, p2b, r, C1, mfromCap, mtoCap);
            }

            atomneedsi = atom.clickable || atom.hoverable;
            atom2needsi = atom2.clickable || atom2.hoverable;

            if (atomneedsi || atom2needsi) {
              if (!mp) mp = new Vector3().addVectors(p1a, p2a).multiplyScalar(0.5);
              if (!mp2) mp2 = new Vector3().addVectors(p1b, p2b).multiplyScalar(0.5);
              if (!mp3) mp3 = new Vector3().addVectors(p1, p2).multiplyScalar(0.5);

              if (atomneedsi) {
                cylinder1a = new Cylinder(p1a.clone(), mp.clone(), r);
                cylinder1b = new Cylinder(p1b.clone(), mp2.clone(), r);
                cylinder1c = new Cylinder(p1.clone(), mp3.clone(), r);
                atom.intersectionShape.cylinder.push(cylinder1a);
                atom.intersectionShape.cylinder.push(cylinder1b);
                atom.intersectionShape.cylinder.push(cylinder1c);
              }
              if (atom2needsi) {
                cylinder2a = new Cylinder(p2a.clone(), mp.clone(), r);
                cylinder2b = new Cylinder(p2b.clone(), mp2.clone(), r);
                cylinder2c = new Cylinder(p2.clone(), mp3.clone(), r);
                atom2.intersectionShape.cylinder.push(cylinder2a);
                atom2.intersectionShape.cylinder.push(cylinder2b);
                atom2.intersectionShape.cylinder.push(cylinder2c);
              }
            }
          }
        }
      }
    }

    // draw non bonded heteroatoms as spheres
    let drawSphere = false;
    let numsinglebonds = 0;
    let differentradii = false;
    // also, if any bonds were drawn as multiples, need sphere
    for (i = 0; i < atom.bonds.length; i++) {
      singleBond = atomSingleBond;
      if (atom.bondStyles && atom.bondStyles[i]) {
        bstyle = atom.bondStyles[i];
        if (bstyle.singleBond) singleBond = true;
        if (bstyle.radius && bstyle.radius != atomBondR) {
          differentradii = true;
        }
      }
      if (singleBond || atom.bondOrder[i] == 1) {
        numsinglebonds++;
      }
    }

    if (differentradii) {
      // jmol style double/triple bonds - no sphere
      if (numsinglebonds > 0) drawSphere = true; // unless needed as a cap
    } else if (numsinglebonds == 0 && atom.bonds.length > 0) {
      drawSphere = true;
    }

    if (drawSphere) {
      bondR = atomBondR;
      // do not use bond style as this can be variable, particularly
      // with jmol export of double/triple bonds
      if (geo.imposter) {
        drawSphereImposter(geo.sphereGeometry, atom, bondR, C1);
      } else {
        GLDraw.drawSphere(geo, atom, bondR, C1);
      }
    }
  }

  drawAtomImposter(atom, geo) {
    if (!atom.style.sphere) return;
    const style = atom.style.sphere;
    if (style.hidden) return;

    const radius = this.getRadiusFromStyle(atom, style);
    const C = getColorFromStyle(atom, style);

    if ((atom.clickable === true || atom.hoverable) && atom.intersectionShape !== undefined) {
      const center = new Vector3(atom.x, atom.y, atom.z);
      atom.intersectionShape.sphere.push(new Sphere(center, radius));
    }

    drawSphereImposter(geo, atom, radius, C);
  }

  drawAtomInstanced(atom, geo) {
    if (!atom.style.sphere) return;
    const style = atom.style.sphere;
    if (style.hidden) return;

    const radius = this.getRadiusFromStyle(atom, style);
    const C = getColorFromStyle(atom, style);

    const geoGroup = geo.updateGeoGroup(1);
    const startv = geoGroup.vertices;
    const start = startv * 3;
    const {vertexArray} = geoGroup;
    const {colorArray} = geoGroup;
    const {radiusArray} = geoGroup;

    vertexArray[start] = atom.x;
    vertexArray[start + 1] = atom.y;
    vertexArray[start + 2] = atom.z;

    colorArray[start] = C.r;
    colorArray[start + 1] = C.g;
    colorArray[start + 2] = C.b;

    radiusArray[startv] = radius;

    if ((atom.clickable === true || atom.hoverable) && atom.intersectionShape !== undefined) {
      const center = new Vector3(atom.x, atom.y, atom.z);
      atom.intersectionShape.sphere.push(new Sphere(center, radius));
    }

    geoGroup.vertices += 1;
  }

  /** returns a list of atoms in the expanded bounding box, but not in the current one
   *
   *  Bounding box:
   *
   *    [ [ xmin, ymin, zmin ],
   *      [ xmax, ymax, zmax ],
   *      [ xctr, yctr, zctr ] ]
   *
   * */
  expandAtomList(atomList, amt) {
    if (amt <= 0) return atomList;

    const pb = getExtent(atomList); // previous bounding box
    /** @type {Array<Array<number>>} */
    const nb = [[], [], []]; // expanded bounding box

    for (let i = 0; i < 3; i++) {
      nb[0][i] = pb[0][i] - amt;
      nb[1][i] = pb[1][i] + amt;
      nb[2][i] = pb[2][i];
    }

    // look in added box "shell" for new atoms
    const expand = [];
    for (let i = 0; i < this.atoms.length; i++) {
      const {x} = this.atoms[i];
      const {y} = this.atoms[i];
      const {z} = this.atoms[i];

      if (
        x >= nb[0][0] &&
        x <= nb[1][0] &&
        y >= nb[0][1] &&
        y <= nb[1][1] &&
        z >= nb[0][2] &&
        z <= nb[1][2]
      ) {
        if (
          !(
            x >= pb[0][0] &&
            x <= pb[1][0] &&
            y >= pb[0][1] &&
            y <= pb[1][1] &&
            z >= pb[0][2] &&
            z <= pb[1][2]
          )
        ) {
          expand.push(this.atoms[i]);
        }
      }
    }
    return expand;
  }

  // static functions
  /**
   * add atomSpecs to validAtomSelectionSpecs
   * @function GLModel#addAtomSpecs
   * @param {Array} customAtomSpecs - array of strings that can be used as atomSelectionSpecs
   * this is to prevent the 'Unknown Selector x' message on the console for the strings passed
   *
   */
  static addAtomSpecs(customAtomSpecs) {
    for (let i = 0; i < customAtomSpecs.length; i++) {
      if (!GLModel.validAtomSelectionSpecs[customAtomSpecs[i]]) {
        GLModel.validAtomSelectionSpecs[customAtomSpecs[i]] = {};
      }
    }
  }

  static parseCrd(data, format) {
    let values = []; // this will contain the all the float values in the

    // file.
    let counter = 0;
    if (format == 'pdb') {
      let index = data.indexOf('\nATOM');
      while (index != -1) {
        while (
          data.slice(index, index + 5) == '\nATOM' ||
          data.slice(index, index + 7) == '\nHETATM'
        ) {
          values[counter++] = parseFloat(data.slice(index + 31, index + 39));
          values[counter++] = parseFloat(data.slice(index + 39, index + 47));
          values[counter++] = parseFloat(data.slice(index + 47, index + 55));
          index = data.indexOf('\n', index + 54);
          if (data.slice(index, index + 4) == '\nTER') index = data.indexOf('\n', index + 5);
        }
        index = data.indexOf('\nATOM', index);
      }
    } else if (format == 'netcdf') {
      const reader = new netcdfjs(data);
      values = [...reader.getDataVariable('coordinates')];
    } else if (format == 'array' || Array.isArray(data)) {
      return data.flat(2);
    } else {
      let index = data.indexOf('\n'); // remove the first line containing title
      if (format == 'inpcrd') {
        index = data.indexOf('\n', index + 1); // remove second line w/#atoms
      }

      data = data.slice(index + 1);
      values = data.match(/\S+/g).map(parseFloat);
    }
    return values;
  }

  /**
   * 
   * @param {string|Uint8Array} data 
   * @param {import('./specs').FileFomats|string|undefined} format 
   * @param {import('./specs').ParserOptionsSpec} options 
   * @returns {import('./specs').ParserResult}
   */
  static parseMolData(data, format, options) {
    format = format || '';
    if (!data) return []; // leave an empty model

    if (/\.gz$/.test(format)) {
      // unzip gzipped files
      format = format.replace(/\.gz$/, '');
      try {
        data = pako.inflate(data, {to: 'string'});
      } catch (err) {
        console.log(err);
      }
    }

    if (typeof Parsers[format] == 'undefined') {
      // let someone provide a file name and get format from extension
      format = format.split('.').pop();
      if (typeof Parsers[format] == 'undefined') {
        console.log(`Unknown format: ${format}`);
        // try to guess correct format from data contents
        if (data instanceof Uint8Array) {
          format = 'mmtf'; // currently only supported binary format?
        } else if (data.match(/^@<TRIPOS>MOLECULE/gm)) {
          format = 'mol2';
        } else if (data.match(/^data_/gm) && data.match(/^loop_/gm)) {
          format = 'cif';
        } else if (data.match(/^HETATM/gm) || data.match(/^ATOM/gm)) {
          format = 'pdb';
        } else if (data.match(/ITEM: TIMESTEP/gm)) {
          format = 'lammpstrj';
        } else if (data.match(/^.*\n.*\n.\s*(\d+)\s+(\d+)/gm)) {
          format = 'sdf'; // could look at line 3
        } else if (data.match(/^%VERSION\s+VERSION_STAMP/gm)) {
          format = 'prmtop';
        } else {
          format = 'xyz';
        }
        console.log(`Best guess: ${format}`);
      }
    }
    const parse = Parsers[format];
    const parsedAtoms = parse(data, options);

    return parsedAtoms;
  }
}
