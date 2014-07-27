/** 
 * All of the functionality of WebMol.js is contained within the
 * WebMol global namespace
 * @namespace */
var WebMol = {};

var colorlike = {};
colorlike.r;
colorlike.g;
colorlike.b;
colorlike.a;

/**
 * Create and initialize an appropriate viewer at supplied HTML element using specification in config
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {ViewerSpec} config Viewer specification
 * @return {WebMol.GLViewer} GLViewer
 * 
 * @example
 * // Assume there exists an HTML div with id "gldiv"
 * var element = $("#gldiv");
 * 
 * // Viewer config - properties 'defaultcolors' and 'callback'
 * var config = {defaultcolors: WebMol.rasmolElementColors,
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
 * var myviewer = WebMol.createViewer(element, config);
 *      
 */
WebMol.createViewer = function(element, config) {};

/**
 * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
 * 
 * @function WebMol.download
 * @param {string} query String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
 * @param {WebMol.GLViewer} viewer - Add new model to existing viewer
 * @example
 * var myviewer = WebMol.createViewer(gldiv);
 * 
 * // GLModel 'm' created and loaded into glviewer for PDB id 2POR
 * var m = WebMol.download('pdb: 2POR', myviewer);
 * 
 * @return {WebMol.GLModel} GLModel
 */ 
WebMol.download = function(query, viewer) {};

/**
 * @ignore
 * @param {WebMol.Object3D} group
 * @param {AtomSpec} atomlist
 * @param {WebMol.ColorScheme} gradientscheme
 */
WebMol.drawCartoon = function(group, atomlist, gradientscheme) {};

/** @property JmolElementColors - Jmol style atom colors (default color scheme) */
WebMol.JmolElementColors = {};
/** @property rasmolElementColors - Rasmol style atom colors */
WebMol.rasmolElementColors = {};

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
 * @type {function(WebMol.GLViewer)} */
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
 * @prop {function(this, WebMol.GLViewer)} callback - Callback click handler function to be executed on this atom and its parent viewer
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
/** @type {function(AtomSpec, WebMol.GLViewer)} */
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
 * WebGL WebMol viewer
 * Note: The preferred method of instantiating a GLViewer is through {@link WebMol.createViewer} 
 * 
 * @constructor 
 * @param {Object} element HTML element within which to create viewer
 * @param {function} callback - Callback function to be immediately executed on this viewer
 * @param {Object} defaultcolors - Object defining default atom colors as atom => color property value pairs for all models within this viewer
 */
WebMol.GLViewer = function(element, callback, defaultcolors) {};

/**
 * Set the background color (default white)
 * 
 * @function WebMol.GLViewer#setBackgroundColor
 * @param {number} hex Hexcode specified background color
 * @param {number} a Alpha level (default 1.0)
 * 
 * @example
 * 
 * //Set 'myviewer' background color to white
 * myviewer.setBackgroundColor(0xffffff)
 * 
 */
WebMol.GLViewer.setBackgroundColor = function(hex, a) {};

/**
 * Set viewer width
 * 
 * @function WebMol.GLViewer#setWidth
 * @param {number} w - Width in pixels
 */
WebMol.GLViewer.setWidth = function(w) {};

/**
 * Set viewer height
 * 
 * @function WebMol.GLViewer#setHeight
 * @param {number} h - Height in pixels
 */
WebMol.GLViewer.setHeight = function(h) {};

/**
 * Resize viewer according to containing HTML element's dimensions
 * 
 * @function WebMol.GLViewer#resize
*/
WebMol.GLViewer.resize = function() {};

/**
 * Return specified model
 * 
 * @function WebMol.GLViewer#getModel
 * @param {number=} [id=last model id] - Retrieve model with specified id
 * @default Returns last model added to viewer
 * @return {WebMol.GLModel}
 * 
 * @example
 * // Retrieve reference to first GLModel added
 * var m = glviewer.getModel(0);
 */
WebMol.GLViewer.getModel = function(id) {};

WebMol.GLViewer.getView = function() {};

/** @param {Array.<number>} arg */
WebMol.GLViewer.setView = function(arg) {};

// apply styles, models, etc in viewer
/**
 * Render current state of viewer, after 
 * adding/removing models, applying styles, etc.
 * 
 * @function WebMol.GLViewer#render
 */
WebMol.GLViewer.render = function() {};

/**
 * Return pdb output of selected atoms (if atoms from pdb input)
 * 
 * @function WebMol.GLViewer#pdbData  
 * @param {Object=} [sel] - Selection specification specifying model and atom properties to select.  Default: all atoms in viewer
 * @return {string} PDB string of selected atoms
 */
WebMol.GLViewer.pdbData = function(sel) {};

/**
 * Zoom to center of atom selection
 * 
 * @function WebMol.GLViewer#zoomTo
 * @param {Object=} [sel] - Selection specification specifying model and atom properties to select. Default: all atoms in viewer
 * 
 * @example
 * // Assuming we have created a model of a protein with multiple chains (e.g. from a PDB file), focus on atoms in chain B
 * glviewer.zoomTo({chain: 'B'});
 * 
 * // Focus on centroid of all atoms of all models in this viewer
 * glviewer.zoomTo();  // (equivalent to glviewer.zoomTo({}) )
 */
WebMol.GLViewer.zoomTo = function(sel) {};

/**
 * Add label to viewer
 * 
 * @function WebMol.GLViewer#addLabel
 * @param {string} text - Label text
 * @param {LabelSpec} data - Label style specification
 * @return {WebMol.Label}
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
WebMol.GLViewer.addLabel = function(text, data) {};

/**
 * Remove label from viewer
 * 
 * @function WebMol.GLViewer#removeLabel
 * @param {WebMol.Label} label - WebMol label
 * 
 * @example
 * // Remove labels created in [addLabel example]{@link WebMol.GLViewer#addLabel}
 * 
 * for (var i = 0; i < labels.length; i++) {
 *     glviewer.removeLabel(label);
 * }
 * 
 * glviewer.render();
 */
WebMol.GLViewer.removeLabel = function(label) {};

//Modify label style
/**
 * Modify existing label's style
 * 
 * @function WebMol.GLViewer#setLabelStyle
 * @param {WebMol.Label} label - WebMol label
 * @param {LabelSpec} stylespec - Label style specification
 * @return {WebMol.Label}
 */
WebMol.GLViewer.setLabelStyle = function(label, stylespec) {};

//Change label text
/**
 * Modify existing label's text
 * 
 * @function WebMol.GLViewer#setLabelText
 * @param {WebMol.Label} label - WebMol label
 * @param {String} text - Label text
 * @return {WebMol.Label}
 */
WebMol.GLViewer.setLabelText = function(label, text) {};

/**
 * Add shape object to viewer 
 * @see {@link WebMol.GLShape}
 * 
 * @function WebMol.GLViewer#addShape
 * @param {ShapeSpec} shapeSpec - style specification for label
 * @return {WebMol.GLShape}
 */
WebMol.GLViewer.addShape = function(shapeSpec) {};

/**
 * Remove shape object from viewe
 *
 * @function WebMol.GLViewer#removeShape
 * @param {WebMol.GLShape} shape - Reference to shape object to remove
 */
WebMol.GLViewer.removeShape = function(shape) {};

/**
 * Create and add sphere shape. This method provides a shorthand 
 * way to create a spherical shape object
 * 
 * @function WebMol.GLViewer#addSphere
 * @param {SphereSpec} spec - Sphere shape style specification
 * @return {WebMol.GLShape}
 */
WebMol.GLViewer.addSphere = function(spec) {};

/**
 * Create and add arrow shape
 * 
 * @function WebMol.GLViewer#addArrow
 * @param {ArrowSpec} spec - Style specification
 * @return {WebMol.GLShape}
 */
WebMol.GLViewer.addArrow = function(spec) {};

/**
 * Add custom shape component from user supplied function
 * 
 * @function WebMol.GLViewer#addCustom
 * @param {CustomSpec} spec - Style specification
 * @return {WebMol.GLShape}
 */
WebMol.GLViewer.addCustom = function(spec) {};

/**
 * Construct isosurface from volumetric data in gaussian cube format
 * 
 * @function WebMol.GLViewer#addVolumetricData
 * @param {String} data - Input file contents 
 * @param {String} format - Input file format (currently only supports "cube")
 * @param {VolSpec} spec - Shape style specification
 * @return {WebMol.GLShape}
 */
WebMol.GLViewer.addVolumetricData = function(data, format, spec) {};

/**
 * Create and add model to viewer, given molecular data and its format 
 * (pdb, sdf, xyz, or mol2)
 * 
 * @function WebMol.GLViewer#addModel
 * @param {string} data - Input data
 * @param {string} format - Input format ('pdb', 'sdf', 'xyz', or 'mol2')
 * @return {WebMol.GLModel}
 */
WebMol.GLViewer.addModel = function(data, format) {};

/**
 * Delete specified model from viewer
 * 
 * @function WebMol.GLViewer#removeModel
 * @param {WebMol.GLModel} model
 */
WebMol.GLViewer.removeModel = function(model) {};

/** 
 * Delete all existing models
 * @function WebMol.GLViewer#removeAllModels
 */
WebMol.GLViewer.removeAllModels = function() {};

/**
 * Create a new model from atoms specified by sel.
 * If extract, removes selected atoms from existing models 
 * 
 * @function WebMol.GLViewer#createModelFrom
 * @param {Object} sel - Atom selection specification
 * @param {boolean=} extract - If true, remove selected atoms from existing models
 * @return {WebMol.GLModel}
 */
WebMol.GLViewer.createModelFrom = function(sel, extract) {};

/**
 * Add surface representation to atoms
 * 
 * @param {WebMol.SurfaceType} type - Surface type
 * @param {Object} style - optional style specification for surface material (e.g. for different coloring scheme, etc)
 * @param {AtomSpec} atomsel - Show surface for atoms in this selection
 * @param {AtomSpec} allsel - Use atoms in this selection to calculate surface; may be larger group than 'atomsel' 
 * @param {AtomSpec} focus - Optionally begin rendering surface specified atoms
 * 
 * @return {number} surfid - Identifying number for this surface
 */
WebMol.GLViewer.addSurface = function(type, style, atomsel, allsel, focus) {};

/**
 * WebMol surface types
 * @enum {number}
 */
WebMol.SurfaceType = {};
WebMol.SurfaceType.VDW;
WebMol.SurfaceType.MS;
WebMol.SurfaceType.SAS;
WebMol.SurfaceType.SES;

/** 
 * Render surface synchronously if true
 * @param {boolean} [WebMol.SyncSurface=false]
 * @type {boolean} */
WebMol.syncSurface;

/**
 * Set the surface material to something else, must render change
 * 
 * @param {number} surf - Surface ID to apply changes to
 * @param {matSpec} style - new material style specification
 */ 
WebMol.GLViewer.setSurfaceMaterialStyle = function(surf, style) {};

/**
 * Remove surface with given ID
 * 
 * @param {number} surf - surface id
 */
WebMol.GLViewer.removeSurface = function(surf) {};

/**
 * Set style properties to all selected atoms
 * 
 * @function WebMol.GLViewer#setStyle
 * @param {AtomSpec} sel - Atom selection specification
 * @param {AtomSpec.style} style - Style spec to apply to specified atoms
 */
WebMol.GLViewer.setStyle = function(sel, style) {};

/**
 * Add style properties to all selected atoms
 * 
 * @function WebMol.GLViewer#addStyle
 * @param {Object} sel - Atom selection specification
 * @param {Object} style - style spec to add to specified atoms
 */
WebMol.GLViewer.addStyle = function(sel, style) {};

/**
 * @function WebMol.GLViewer#setColorByProperty
 * @param {type} sel
 * @param {type} prop
 * @param {type} scheme
 */
WebMol.GLViewer.setColorByProperty = function(sel, prop, scheme) {};

/**
 * @function WebMol.GLViewer#setColorByElement
 * @param {type} sel
 * @param {type} colors
 */
WebMol.GLViewer.setColorByElement = function(sel, colors) {};

/** Clear all models, surfaces, and shapes from viewer */
WebMol.GLViewer.clear = function() {};

// props is a list of objects that select certain atoms and enumerate
// properties for those atom
/**
 * Add specified properties to all atoms matching input argument
 * @param {AtomSpec} props
 */
WebMol.GLViewer.mapAtomProperties = function(props) {};

/**
 * GLModel represents a group of related atoms
 * @constructor 
 * @param {number=} mid 
 * @param {Object=} defaultcolors Object defining default atom colors as atom => color property value pairs
 * @see WebMol.download
 */
WebMol.GLModel = function(mid, defaultcolors) {};

/**
 * Returns model id number
 * 
 * @function WebMol.GLMode#getID
 * @return {number} Model ID
 */
WebMol.GLModel.getID = function() {};

/** add atoms to this model from molecular data string
 * 
 * @function WebMol.GLModel#addMolData
 * @param {string} data - atom structure file input data string
 * @param {string} format - input file string format (e.g 'pdb', 'sdf', etc.)
 */
WebMol.GLModel.addMolData = function(data, format) {};

/** given a selection specification, return true if atom is selected
 * 
 * @function WebMol.GLModel#atomIsSelected
 * @param {type} atom
 * @param {type} sel
 * @return {boolean}
 */
WebMol.GLModel.atomIsSelected = function(atom, sel) {};

/** return list of atoms selected by sel, this is specific to glmodel
 * 
 * @function WebMol.GLModel#selectedAtoms
 * @param {type} sel
 * @return {Array.<Object>}
 */
WebMol.GLModel.selectedAtoms = function(sel) {};

/** Add list of new atoms to model
 * 
 * @function WebMol.GLModel#addAtoms
 * @param {type} newatoms
 */
WebMol.GLModel.addAtoms = function(newatoms) {};

/** Remove specified atoms from model
 * 
 * @function WebMol.GLModel#removeAtoms
 * @param {type} badatoms
 * @return {removeAtoms}
 */
WebMol.GLModel.removeAtoms = function(badatoms) {};

/** Set atom style of selected atoms
 * 
 * @function WebMol.GLModel#setStyle
 * @param {type} sel
 * @param {type} style
 * @param {type} add
 */
WebMol.GLModel.setStyle = function(sel, style, add) {};

/** given a mapping from element to color, set atom colors
 * 
 * @function WebMol.GLModel#setColorByElement
 * @param {type} sel
 * @param {type} colors
 */
WebMol.GLModel.setColorByElement = function(sel, colors) {};

/**
 * @function WebMol.GLModelSetColorByProperty
 * @param {type} sel
 * @param {type} prop
 * @param {type} scheme
 */
WebMol.GLModel.setColorByProperty = function(sel, prop, scheme) {};

/** manage the globj for this model in the possed modelGroup - if it has to be regenerated, remove and add
 * 
 * @function WebMol.GLModel#globj
 * @param {WebMol.Object3D} group
 */
WebMol.GLModel.globj = function(group) {};

/** Remove any renderable mol object from scene
 * 
 * @function WebMol.GLModel#removegl
 * @param {WebMol.Object3D} group
 */
WebMol.GLModel.removegl = function(group) {};


WebMol.LabelCount;

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
 * @type {colorlike | WebMol.Color} */
LabelSpec.fontColor;

LabelSpec.borderThickness;
/** @type {colorlike} */
LabelSpec.borderColor;
/** @type {colorlike} */
LabelSpec.backgroundColor;
/**
 * Label position
 * @type {WebMol.Vector3}
 */
LabelSpec.position;

/** labels always rendered in front of model(s) if true
 * 
 * @type {boolean}
 */
LabelSpec.inFront;

/**
 * Renderable labels
 * @constructor WebMol.Label
 * @extends {LabelSpec}
 * @param {string} tag - Label text
 * @param {Object} parameters Label style and font specifications
 */
WebMol.Label = function(text, parameters) {};
        
WebMol.Label.id;    
/** @type {LabelSpec} */
WebMol.Label.stylespec;
WebMol.Label.canvas;
WebMol.Label.context;
/** @type {WebMol.Sprite} */
WebMol.Label.sprite;
WebMol.Label.text;

WebMol.Label.prototype.setContext = function() {};
WebMol.Label.prototype.dispose = function() {};

/** 
 * GLShape style specification
 * @typedef
 */
var ShapeSpec = {};
/** @type {WebMol.Color} */
ShapeSpec.color;
ShapeSpec.wireframe;
ShapeSpec.alpha;
ShapeSpec.side;
ShapeSpec.clickable;
/** @type {function(WebMol.GLShape, WebMol.GLViewer)} */
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
/** @type {WebMol.Vector3} */
SphereSpec.center;

/**
 * Arrow shape specification
 * @typedef
 */
var ArrowSpec = {};
/** @var {WebMol.Vector3} ArrowSpec.start - Arrow start point*/
ArrowSpec.start;
/** @property {WebMol.Vector3} */
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
 * @constructor WebMol.GLShape
 * @extends {ShapeSpec}
 * @param {number} sid - Unique identifier
 * @param {ShapeSpec} stylespec - shape style specification
 */
WebMol.GLShape = function(sid, stylespec) {};

WebMol.GLShape.id;
WebMol.GLShape.boundingSphere;
/** @type {IntersectionShapes} */
WebMol.GLShape.intersectionShape;
/** @type {WebMol.Vector3} */
WebMol.GLShape.position;
WebMol.GLShape.x;
WebMol.GLShape.y;
WebMol.GLShape.z;

/** Update shape with new style specification
 * @param {ShapeSpec} newspec
 * @return {WebMol.GLShape}
 */
WebMol.GLShape.updateStyle = function(newspec) {};

/**
 * Creates a custom shape from supplied vertex and face arrays
 * @param {CustomSpec} customSpec
 * @return {WebMol.GLShape}
 */
WebMol.GLShape.addCustom = function(customSpec) {};
        
       
/**
 * Creates a sphere shape
 * @param {SphereSpec} sphereSpec
 * @return {WebMol.GLShape}
 */
WebMol.GLShape.addSphere = function(sphereSpec) {};    

/**
 * Creates an arrow shape
 * @param {ArrowSpec} arrowSpec
 * @return {WebMol.GLShape}
 */
WebMol.GLShape.addArrow = function(arrowSpec) {};

/** 
 * Creates custom shape from volumetric data 
 * @param {string} data - Volumetric input data 
 * @param {string} fmt - Input data format (e.g. 'cube' for cube file format)
 * @param {VolSpec} volSpec - Volumetric data shape specification
 * @return {WebMol.GLShape}
 */
WebMol.GLShape.addVolumetricData = function(data, fmt, volSpec) {};

/**
 * Initialize webgl objects for rendering
 * @param {WebMol.Object3D} group
 * 
 */  
WebMol.GLShape.globj = function(group) {};

WebMol.ShapeIDCount;


//color schemes
/** Color mapping scheme
 * @interface
 * @param {number} min
 * @param {number} max
 */
WebMol.ColorScheme = function(min, max) {};

/**
 * Map value to hex color
 * @param {number} val
 * @param {number} range
 * @returns {number}
 */
WebMol.ColorScheme.valueToHex = function(val, range) {};
WebMol.ColorScheme.jmolID = function() {};
//return range used for color mapping, null if none set
WebMol.ColorScheme.range = function() {};

/**
 * Color scheme red to white to blue, for charges
 * @constructor
 * @implements {WebMol.ColorScheme}
 */
WebMol.RWB = function(min, max) {};

/**
 * rainbow gradient, but without purple to match jmol
 * @constructor
 * @implements {WebMol.ColorScheme}
 */
WebMol.ROYGB = function(min, max) {};

/**
 * rainbow gradient with constant saturation, all the way to purple!
 * @constructor
 * @implements {WebMol.ColorScheme}
 */
WebMol.Sinebow = function(min, max) {};


