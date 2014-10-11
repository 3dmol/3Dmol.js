// Specifications for various object types used in 3Dmol.js
// This is primarily for documentation 

/**
 * GLViewer input specification
 * @typedef ViewerSpec
 */
var ViewerSpec = {};
ViewerSpec.order;
ViewerSpec.defaultcolors;
/** 
 * @type {function($3Dmol.GLViewer)} */
ViewerSpec.callback;

/**
 * Object literal Atom representation.  Can be used as a selection specification to 
 * select all atoms with matching properties
 * @typedef AtomSpec
 * @struct
 * @prop {string} resn - Parent residue name
 * @prop {number} x - Atom's x coordinate
 * @prop {number} y - Atom's y coordinate
 * @prop {number} z - Atom's z coordinate
 * @prop {number} color - Atom's color, as hex code
 * @prop {number} surfaceColor - Hex code for color to be used for surface patch over this atom
 * @prop {string} elem - Element abbreviation (e.g. 'H', 'Ca', etc)
 * @prop {boolean} hetflag - Set to true if atom is a heteroatom
 * @prop {string} chain - Chain this atom belongs to, if specified in input file (e.g 'A' for chain A)
 * @prop {number} resi - Residue number 
 * @prop {number} icode
 * @prop {number} rescode
 * @prop {number} serial - Atom's serial id number
 * @prop {string} atom - Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon)
 * @prop {Array.<number>} bonds - Array of atom ids this atom is bonded to
 * @prop {string} ss - Secondary structure identifier (for cartoon render; e.g. 'h' for helix)
 * @prop {boolean} singleBonds - true if this atom forms only single bonds or no bonds at all
 * @prop {Array.<number>} bondOrder - Array of this atom's bond orders, corresponding to bonds identfied by 'bonds'
 * @prop {Object} properties - Optional mapping of additional properties
 * @prop {number} b - Atom b factor data
 * @prop {string} pdbline - If applicable, this atom's record entry from the input PDB file (used to output new PDB from models)
 * @prop {boolean} clickable - Set this flag to true to enable click selection handling for this atom
 * @prop {function(this, $3Dmol.GLViewer)} callback - Callback click handler function to be executed on this atom and its parent viewer
 * @prop {AtomStyleSpec} style - Atom style specification
 */
var AtomSpec = {};
AtomSpec.resn;
AtomSpec.x;
AtomSpec.y;
AtomSpec.z;
AtomSpec.color;
AtomSpec.surfaceColor;
AtomSpec.elem;
AtomSpec.hetflag;
AtomSpec.chain;
AtomSpec.resi;
AtomSpec.icode;
AtomSpec.rescode;
AtomSpec.serial;
AtomSpec.atom;
AtomSpec.bonds;
AtomSpec.ss;
AtomSpec.singleBonds;
AtomSpec.bondOrder;
AtomSpec.properties;
AtomSpec.b;
AtomSpec.pdbline;
/** @type {IntersectionShapes} */
AtomSpec.intersectionShape;
AtomSpec.clickable;
/** @type {function(AtomSpec, $3Dmol.GLViewer)} */
AtomSpec.callback;

/** 
 * @typedef AtomStyleSpec
 */
var AtomStyleSpec = {};

AtomSpec.style = {};
/** @type {atomstyle} */
AtomSpec.style.line;
/** @type {atomstyle} */
AtomSpec.style.cross;
/** @type {atomstyle} */
AtomSpec.style.sphere;
/** @type {atomstyle} */
AtomSpec.style.stick;
/** @type {atomstyle} */
AtomSpec.style.cartoon;
AtomSpec.style.cartoon.gradient;


/**
 * Label type specification
 * @typedef
 */
var LabelSpec = {};

/** Label text font style
 * @type {string} */
LabelSpec.font;

/** Label text font pt size
 * @type {number} */
LabelSpec.fontSize;

/** Label font color - specify with an object with r, g, b, and a (alpha) values
 * @type {colorlike | $3Dmol.Color} */
LabelSpec.fontColor;

LabelSpec.borderThickness;
/** @type {colorlike} */
LabelSpec.borderColor;
/** @type {colorlike} */
LabelSpec.backgroundColor;
/**
 * Label position
 * @type {$3Dmol.Vector3}
 */
LabelSpec.position;

/** labels always rendered in front of model(s) if true
 * 
 * @type {boolean}
 */
LabelSpec.inFront;


/** 
 * GLShape style specification
 * @typedef
 */
var ShapeSpec = {};
/** @type {$3Dmol.Color} */
ShapeSpec.color;
ShapeSpec.wireframe;
ShapeSpec.alpha;
ShapeSpec.side;
ShapeSpec.clickable;
/** @type {function($3Dmol.GLShape, $3Dmol.GLViewer)} */
ShapeSpec.callback;

/**
 * Specification for adding custom shape
 * @typedef
 */
var CustomShapeSpec = {};
CustomShapeSpec.vertexArr;
CustomShapeSpec.faceArr;
CustomShapeSpec.normalArr;
CustomShapeSpec.lineArr;

/**
 * Sphere shape specification
 * @typedef
 */
var SphereSpec = {};
SphereSpec.radius;
/** @type {$3Dmol.Vector3} */
SphereSpec.center;

/**
 * Arrow shape specification
 * @typedef
 */
var ArrowSpec = {};
/** @var {$3Dmol.Vector3} ArrowSpec.start - Arrow start point*/
ArrowSpec.start;
/** @property {$3Dmol.Vector3} */
ArrowSpec.end;
ArrowSpec.radius;
ArrowSpec.radiusRatio;
ArrowSpec.mid;


/**
 * Volumetric data specification
 * @typedef
 */
var VolSpec = {};
VolSpec.isoval;
VolSpec.voxel;
