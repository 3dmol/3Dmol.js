/** 
 * All of the functionality of $3Dmol.js is contained within the
 * $3Dmol global namespace
 * @namespace */
var $3Dmol = {};

var colorlike = {};
colorlike.r;
colorlike.g;
colorlike.b;
colorlike.a;

/**
 * Create and initialize an appropriate viewer at supplied HTML element using specification in config
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {ViewerSpec} config Viewer specification
 * @return {$3Dmol.GLViewer} GLViewer
 * 
 * @example
 * // Assume there exists an HTML div with id "gldiv"
 * var element = $("#gldiv");
 * 
 * // Viewer config - properties 'defaultcolors' and 'callback'
 * var config = {defaultcolors: $3Dmol.rasmolElementColors,
 *               callback : function(viewer) {
 *                            //'data' is a string containing molecule data in pdb format  
 *                            viewer.addModel(data, "pdb");
 *                            viewer.zoomTo();
 *                            viewer.render();
 *                          }  
 *                        
 *               };
 * 
 * // Create GLViewer within 'gldiv' and execute callback
 * var myviewer = $3Dmol.createViewer(element, config);
 *      
 */
$3Dmol.createViewer = function(element, config) {};

/**
  * Contains a dictionary of embedded viewers created from HTML elements
  * with a the viewer_3Dmoljs css class indexed by their id (or numerically
  * if they do not have an id).
 */
$3Dmol.viewers = {}

/**
 * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
 * 
 * @function $3Dmol.download
 * @param {string} query String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
 * @param {$3Dmol.GLViewer} viewer - Add new model to existing viewer
 * @example
 * var myviewer = $3Dmol.createViewer(gldiv);
 * 
 * // GLModel 'm' created and loaded into glviewer for PDB id 2POR
 * var m = $3Dmol.download('pdb: 2POR', myviewer);
 * 
 * @return {$3Dmol.GLModel} GLModel
 */ 
$3Dmol.download = function(query, viewer) {};

/**
 * @ignore
 * @param {$3Dmol.Object3D} group
 * @param {AtomSpec} atomlist
 * @param {$3Dmol.ColorScheme} gradientscheme
 */
$3Dmol.drawCartoon = function(group, atomlist, gradientscheme) {};

/** Preset element coloring - from individual element colors to entire mappings (e.g. '$3Dmol.elementColors.Jmol' colors atoms with Jmol stylings)
 * @struct
 */
$3Dmol.elementColors = {};

$3Dmol.elementColors.defaultColor;

/** @property Jmol-like element colors*/
$3Dmol.elementColors.Jmol = {};

/** @property rasmol-like element colors */
$3Dmol.elementColors.rasmol = {};

$3Dmol.elementColors.defaultColors = $3Dmol.elementColors.rasmol;

$3Dmol.elementColors.greenCarbon = $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.greenCarbon['C'] = 0x00ff00;

$3Dmol.elementColors.cyanCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.cyanCarbon['C'] = 0x00ffff;

$3Dmol.elementColors.magentaCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.magentaCarbon['C'] = 0xff00ff;

$3Dmol.elementColors.yellowCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.yellowCarbon['C'] = 0xffff00;

$3Dmol.elementColors.whiteCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.whiteCarbon['C'] = 0xffffff;

$3Dmol.elementColors.orangeCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.orangeCarbon['C'] = 0xff6600;

$3Dmol.elementColors.purpleCarbon =  $.extend({},$3Dmol.elementColors.defaultColors);
$3Dmol.elementColors.purpleCarbon['C'] = 0x800080;

//Specification arguments 
//TODO: flesh out the annotations

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

//Viewer
// The constructor
/**
 * WebGL $3Dmol viewer
 * Note: The preferred method of instantiating a GLViewer is through {@link $3Dmol.createViewer} 
 * 
 * @constructor 
 * @param {Object} element HTML element within which to create viewer
 * @param {function} callback - Callback function to be immediately executed on this viewer
 * @param {Object} defaultcolors - Object defining default atom colors as atom => color property value pairs for all models within this viewer
 */
$3Dmol.GLViewer = function(element, callback, defaultcolors) {};

/**
 * Set the background color (default white)
 * 
 * @function $3Dmol.GLViewer#setBackgroundColor
 * @param {number} hex Hexcode specified background color
 * @param {number} a Alpha level (default 1.0)
 * 
 * @example
 * 
 * //Set 'myviewer' background color to white
 * myviewer.setBackgroundColor(0xffffff)
 * 
 */
$3Dmol.GLViewer.setBackgroundColor = function(hex, a) {};

/**
 * Set viewer width
 * 
 * @function $3Dmol.GLViewer#setWidth
 * @param {number} w - Width in pixels
 */
$3Dmol.GLViewer.setWidth = function(w) {};

/**
 * Set viewer height
 * 
 * @function $3Dmol.GLViewer#setHeight
 * @param {number} h - Height in pixels
 */
$3Dmol.GLViewer.setHeight = function(h) {};

/**
 * Resize viewer according to containing HTML element's dimensions
 * 
 * @function $3Dmol.GLViewer#resize
*/
$3Dmol.GLViewer.resize = function() {};

/**
 * Return specified model
 * 
 * @function $3Dmol.GLViewer#getModel
 * @param {number=} [id=last model id] - Retrieve model with specified id
 * @default Returns last model added to viewer
 * @return {$3Dmol.GLModel}
 * 
 * @example
 * // Retrieve reference to first GLModel added
 * var m = glviewer.getModel(0);
 */
$3Dmol.GLViewer.getModel = function(id) {};

$3Dmol.GLViewer.getView = function() {};

/** @param {Array.<number>} arg */
$3Dmol.GLViewer.setView = function(arg) {};

// apply styles, models, etc in viewer
/**
 * Render current state of viewer, after 
 * adding/removing models, applying styles, etc.
 * 
 * @function $3Dmol.GLViewer#render
 */
$3Dmol.GLViewer.render = function() {};

/**
 * Return pdb output of selected atoms (if atoms from pdb input)
 * 
 * @function $3Dmol.GLViewer#pdbData  
 * @param {Object=} [sel] - Selection specification specifying model and atom properties to select.  Default: all atoms in viewer
 * @return {string} PDB string of selected atoms
 */
$3Dmol.GLViewer.pdbData = function(sel) {};

/**
 * Zoom to center of atom selection
 * 
 * @function $3Dmol.GLViewer#zoomTo
 * @param {Object=} [sel] - Selection specification specifying model and atom properties to select. Default: all atoms in viewer
 * 
 * @example
 * // Assuming we have created a model of a protein with multiple chains (e.g. from a PDB file), focus on atoms in chain B
 * glviewer.zoomTo({chain: 'B'});
 * 
 * // Focus on centroid of all atoms of all models in this viewer
 * glviewer.zoomTo();  // (equivalent to glviewer.zoomTo({}) )
 */
$3Dmol.GLViewer.zoomTo = function(sel) {};

/**
 * Add label to viewer
 * 
 * @function $3Dmol.GLViewer#addLabel
 * @param {string} text - Label text
 * @param {LabelSpec} data - Label style specification
 * @return {$3Dmol.Label}
 * 
 * @example
 * 
 * // Assuming glviewer contains a model representing a protein, label all alpha carbons with their residue name
 * 
 * // Select all alpha carbons (have property atom : "CA") from last model added
 * var atoms = glviewer.getModel().selectedAtoms({atom:"CA"});
 * var labels = [];
 * 
 * for (var a in atoms) {
 *     var atom = atoms[a];
 * 
 *     // Create label at alpha carbon's position displaying atom's residue and residue number
 *     var labelText = atom.resname + " " + atom.resi;
 *      
 *     var l = glviewer.createLabel(labelText, {fontSize: 12, position: {x: atom.x, y: atom.y, z: atom.z});
 * 
 *     labels.push(l);
 * }
 * 
 * // Render labels
 * glviewer.render();
 */
$3Dmol.GLViewer.addLabel = function(text, data) {};

/**
 * Remove label from viewer
 * 
 * @function $3Dmol.GLViewer#removeLabel
 * @param {$3Dmol.Label} label - $3Dmol label
 * 
 * @example
 * // Remove labels created in [addLabel example]{@link $3Dmol.GLViewer#addLabel}
 * 
 * for (var i = 0; i < labels.length; i++) {
 *     glviewer.removeLabel(label);
 * }
 * 
 * glviewer.render();
 */
$3Dmol.GLViewer.removeLabel = function(label) {};

//Modify label style
/**
 * Modify existing label's style
 * 
 * @function $3Dmol.GLViewer#setLabelStyle
 * @param {$3Dmol.Label} label - $3Dmol label
 * @param {LabelSpec} stylespec - Label style specification
 * @return {$3Dmol.Label}
 */
$3Dmol.GLViewer.setLabelStyle = function(label, stylespec) {};

//Change label text
/**
 * Modify existing label's text
 * 
 * @function $3Dmol.GLViewer#setLabelText
 * @param {$3Dmol.Label} label - $3Dmol label
 * @param {String} text - Label text
 * @return {$3Dmol.Label}
 */
$3Dmol.GLViewer.setLabelText = function(label, text) {};

/**
 * Add shape object to viewer 
 * @see {@link $3Dmol.GLShape}
 * 
 * @function $3Dmol.GLViewer#addShape
 * @param {ShapeSpec} shapeSpec - style specification for label
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLViewer.addShape = function(shapeSpec) {};

/**
 * Remove shape object from viewe
 *
 * @function $3Dmol.GLViewer#removeShape
 * @param {$3Dmol.GLShape} shape - Reference to shape object to remove
 */
$3Dmol.GLViewer.removeShape = function(shape) {};

/**
 * Create and add sphere shape. This method provides a shorthand 
 * way to create a spherical shape object
 * 
 * @function $3Dmol.GLViewer#addSphere
 * @param {SphereSpec} spec - Sphere shape style specification
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLViewer.addSphere = function(spec) {};

/**
 * Create and add arrow shape
 * 
 * @function $3Dmol.GLViewer#addArrow
 * @param {ArrowSpec} spec - Style specification
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLViewer.addArrow = function(spec) {};

/**
 * Create and add cylinder shape
 * 
 * @function $3Dmol.GLViewer#addArrow
 * @param {CylinderSpec} spec - Style specification
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLViewer.addCylinder = function(spec) {};

/**
 * Add custom shape component from user supplied function
 * 
 * @function $3Dmol.GLViewer#addCustom
 * @param {CustomSpec} spec - Style specification
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLViewer.addCustom = function(spec) {};

/**
 * Construct isosurface from volumetric data in gaussian cube format
 * 
 * @function $3Dmol.GLViewer#addVolumetricData
 * @param {String} data - Input file contents 
 * @param {String} format - Input file format (currently only supports "cube")
 * @param {VolSpec} spec - Shape style specification
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLViewer.addVolumetricData = function(data, format, spec) {};

/**
 * Create and add model to viewer, given molecular data and its format 
 * (pdb, sdf, xyz, or mol2)
 * 
 * @function $3Dmol.GLViewer#addModel
 * @param {string} data - Input data
 * @param {string} format - Input format ('pdb', 'sdf', 'xyz', or 'mol2')
 * @return {$3Dmol.GLModel}
 */
$3Dmol.GLViewer.addModel = function(data, format) {};

/**
 * Delete specified model from viewer
 * 
 * @function $3Dmol.GLViewer#removeModel
 * @param {$3Dmol.GLModel} model
 */
$3Dmol.GLViewer.removeModel = function(model) {};

/** 
 * Delete all existing models
 * @function $3Dmol.GLViewer#removeAllModels
 */
$3Dmol.GLViewer.removeAllModels = function() {};

/**
 * Create a new model from atoms specified by sel.
 * If extract, removes selected atoms from existing models 
 * 
 * @function $3Dmol.GLViewer#createModelFrom
 * @param {Object} sel - Atom selection specification
 * @param {boolean=} extract - If true, remove selected atoms from existing models
 * @return {$3Dmol.GLModel}
 */
$3Dmol.GLViewer.createModelFrom = function(sel, extract) {};

/**
 * Add surface representation to atoms
 * 
 * @param {$3Dmol.SurfaceType} type - Surface type
 * @param {Object} style - optional style specification for surface material (e.g. for different coloring scheme, etc)
 * @param {AtomSpec} atomsel - Show surface for atoms in this selection
 * @param {AtomSpec} allsel - Use atoms in this selection to calculate surface; may be larger group than 'atomsel' 
 * @param {AtomSpec} focus - Optionally begin rendering surface specified atoms
 * 
 * @return {number} surfid - Identifying number for this surface
 */
$3Dmol.GLViewer.addSurface = function(type, style, atomsel, allsel, focus) {};

/**
 * Adds an explicit mesh as a surface object.
 * 
 * @param {$3Dmol.Mesh}
 *            mesh
 * @param {Object}
 *            style
 * @returns {Number} surfid
 */
$3Dmol.GLViewer.addMesh = function(mesh) {};

/**
 * $3Dmol surface types
 * @enum {number}
 */
$3Dmol.SurfaceType = {};
$3Dmol.SurfaceType.VDW;
$3Dmol.SurfaceType.MS;
$3Dmol.SurfaceType.SAS;
$3Dmol.SurfaceType.SES;

/** 
 * Render surface synchronously if true
 * @param {boolean} [$3Dmol.SyncSurface=false]
 * @type {boolean} */
$3Dmol.syncSurface;

/**
 * Set the surface material to something else, must render change
 * 
 * @param {number} surf - Surface ID to apply changes to
 * @param {matSpec} style - new material style specification
 */ 
$3Dmol.GLViewer.setSurfaceMaterialStyle = function(surf, style) {};

/**
 * Remove surface with given ID
 * 
 * @param {number} surf - surface id
 */
$3Dmol.GLViewer.removeSurface = function(surf) {};

/**
 * Set style properties to all selected atoms
 * 
 * @function $3Dmol.GLViewer#setStyle
 * @param {AtomSpec} sel - Atom selection specification
 * @param {AtomSpec.style} style - Style spec to apply to specified atoms
 */
$3Dmol.GLViewer.setStyle = function(sel, style) {};

/**
 * Add style properties to all selected atoms
 * 
 * @function $3Dmol.GLViewer#addStyle
 * @param {Object} sel - Atom selection specification
 * @param {Object} style - style spec to add to specified atoms
 */
$3Dmol.GLViewer.addStyle = function(sel, style) {};

/**
 * @function $3Dmol.GLViewer#setColorByProperty
 * @param {type} sel
 * @param {type} prop
 * @param {type} scheme
 */
$3Dmol.GLViewer.setColorByProperty = function(sel, prop, scheme) {};

/**
 * @function $3Dmol.GLViewer#setColorByElement
 * @param {type} sel
 * @param {type} colors
 */
$3Dmol.GLViewer.setColorByElement = function(sel, colors) {};

/** Clear all models, surfaces, and shapes from viewer */
$3Dmol.GLViewer.clear = function() {};

// props is a list of objects that select certain atoms and enumerate
// properties for those atom
/**
 * Add specified properties to all atoms matching input argument
 * @param {AtomSpec} props
 */
$3Dmol.GLViewer.mapAtomProperties = function(props) {};

/**
 * GLModel represents a group of related atoms
 * @constructor 
 * @param {number=} mid 
 * @param {Object=} defaultcolors Object defining default atom colors as atom => color property value pairs
 * @see $3Dmol.download
 */
$3Dmol.GLModel = function(mid, defaultcolors) {};

/**
 * Returns model id number
 * 
 * @function $3Dmol.GLMode#getID
 * @return {number} Model ID
 */
$3Dmol.GLModel.getID = function() {};

/** add atoms to this model from molecular data string
 * 
 * @function $3Dmol.GLModel#addMolData
 * @param {string} data - atom structure file input data string
 * @param {string} format - input file string format (e.g 'pdb', 'sdf', etc.)
 */
$3Dmol.GLModel.addMolData = function(data, format) {};

/** given a selection specification, return true if atom is selected
 * 
 * @function $3Dmol.GLModel#atomIsSelected
 * @param {type} atom
 * @param {type} sel
 * @return {boolean}
 */
$3Dmol.GLModel.atomIsSelected = function(atom, sel) {};

/** return list of atoms selected by sel, this is specific to glmodel
 * 
 * @function $3Dmol.GLModel#selectedAtoms
 * @param {type} sel
 * @return {Array.<Object>}
 */
$3Dmol.GLModel.selectedAtoms = function(sel) {};

/** Add list of new atoms to model
 * 
 * @function $3Dmol.GLModel#addAtoms
 * @param {type} newatoms
 */
$3Dmol.GLModel.addAtoms = function(newatoms) {};

/** Remove specified atoms from model
 * 
 * @function $3Dmol.GLModel#removeAtoms
 * @param {type} badatoms
 * @return {removeAtoms}
 */
$3Dmol.GLModel.removeAtoms = function(badatoms) {};

/** Set atom style of selected atoms
 * 
 * @function $3Dmol.GLModel#setStyle
 * @param {type} sel
 * @param {type} style
 * @param {type} add
 */
$3Dmol.GLModel.setStyle = function(sel, style, add) {};

/** given a mapping from element to color, set atom colors
 * 
 * @function $3Dmol.GLModel#setColorByElement
 * @param {type} sel
 * @param {type} colors
 */
$3Dmol.GLModel.setColorByElement = function(sel, colors) {};

/**
 * @function $3Dmol.GLModelSetColorByProperty
 * @param {type} sel
 * @param {type} prop
 * @param {type} scheme
 */
$3Dmol.GLModel.setColorByProperty = function(sel, prop, scheme) {};

/** manage the globj for this model in the possed modelGroup - if it has to be regenerated, remove and add
 * 
 * @function $3Dmol.GLModel#globj
 * @param {$3Dmol.Object3D} group
 */
$3Dmol.GLModel.globj = function(group) {};

/** Remove any renderable mol object from scene
 * 
 * @function $3Dmol.GLModel#removegl
 * @param {$3Dmol.Object3D} group
 */
$3Dmol.GLModel.removegl = function(group) {};


$3Dmol.LabelCount;

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
 * Renderable labels
 * @constructor $3Dmol.Label
 * @extends {LabelSpec}
 * @param {string} tag - Label text
 * @param {Object} parameters Label style and font specifications
 */
$3Dmol.Label = function(text, parameters) {};
        
$3Dmol.Label.id;    
/** @type {LabelSpec} */
$3Dmol.Label.stylespec;
$3Dmol.Label.canvas;
$3Dmol.Label.context;
/** @type {$3Dmol.Sprite} */
$3Dmol.Label.sprite;
$3Dmol.Label.text;

$3Dmol.Label.prototype.setContext = function() {};
$3Dmol.Label.prototype.dispose = function() {};

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

/**
 * A GLShape is a collection of user specified shapes.
 * 
 * @constructor $3Dmol.GLShape
 * @extends {ShapeSpec}
 * @param {number} sid - Unique identifier
 * @param {ShapeSpec} stylespec - shape style specification
 */
$3Dmol.GLShape = function(sid, stylespec) {};

$3Dmol.GLShape.id;
$3Dmol.GLShape.boundingSphere;
/** @type {IntersectionShapes} */
$3Dmol.GLShape.intersectionShape;
/** @type {$3Dmol.Vector3} */
$3Dmol.GLShape.position;
$3Dmol.GLShape.x;
$3Dmol.GLShape.y;
$3Dmol.GLShape.z;

/** Update shape with new style specification
 * @param {ShapeSpec} newspec
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLShape.updateStyle = function(newspec) {};

/**
 * Creates a custom shape from supplied vertex and face arrays
 * @param {CustomSpec} customSpec
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLShape.addCustom = function(customSpec) {};
        
       
/**
 * Creates a sphere shape
 * @param {SphereSpec} sphereSpec
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLShape.addSphere = function(sphereSpec) {};    

/**
 * Creates an arrow shape
 * @param {ArrowSpec} arrowSpec
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLShape.addArrow = function(arrowSpec) {};

/**
 * Creates a Cylinder shape
 * @param {CylinderSpec} cylinderSpec
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLShape.addCylinder = function(cylinderSpec) {};

/** 
 * Creates custom shape from volumetric data 
 * @param {string} data - Volumetric input data 
 * @param {string} fmt - Input data format (e.g. 'cube' for cube file format)
 * @param {VolSpec} volSpec - Volumetric data shape specification
 * @return {$3Dmol.GLShape}
 */
$3Dmol.GLShape.addVolumetricData = function(data, fmt, volSpec) {};

/**
 * Initialize webgl objects for rendering
 * @param {$3Dmol.Object3D} group
 * 
 */  
$3Dmol.GLShape.globj = function(group) {};

$3Dmol.ShapeIDCount;


//color schemes
/** Color mapping scheme
 * @interface
 * @param {number} min
 * @param {number} max
 */
$3Dmol.ColorScheme = function(min, max) {};

/**
 * Map value to hex color
 * @param {number} val
 * @param {number} range
 * @returns {number}
 */
$3Dmol.ColorScheme.valueToHex = function(val, range) {};
$3Dmol.ColorScheme.jmolID = function() {};
//return range used for color mapping, null if none set
$3Dmol.ColorScheme.range = function() {};

/**
 * Color scheme red to white to blue, for charges
 * @constructor
 * @implements {$3Dmol.ColorScheme}
 */
$3Dmol.RWB = function(min, max) {};

/**
 * rainbow gradient, but without purple to match jmol
 * @constructor
 * @implements {$3Dmol.ColorScheme}
 */
$3Dmol.ROYGB = function(min, max) {};

/**
 * rainbow gradient with constant saturation, all the way to purple!
 * @constructor
 * @implements {$3Dmol.ColorScheme}
 */
$3Dmol.Sinebow = function(min, max) {};


