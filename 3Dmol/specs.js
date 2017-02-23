// Specifications for various object types used in 3Dmol.js
// This is primarily for documentation 
(function() {
/**
 * GLViewer input specification
 * @typedef ViewerSpec
 * @prop {Object} defaultcolors - map of elements to colors
 * @prop {boolean} nomouse - if true, disable handling of mouse events
 * @prop {ColorSpec} backgroundColor - color of background
 */

/**
 * Atom representation. Depending on the input file format, not all fields may be defined.
 * @typedef AtomSpec
 * @prop {string} resn - Parent residue name
 * @prop {number} x - Atom's x coordinate
 * @prop {number} y - Atom's y coordinate
 * @prop {number} z - Atom's z coordinate
 * @prop {ColorSpec} color - Atom's color, as hex code or built-in color string
 * @prop {ColorSpec} surfaceColor - Hex code for color to be used for surface patch over this atom
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
 * @prop {boolean} invert - for selection, inverts the meaning of the selection
 */

 /**
  * Parser options specification. Used to specify the options of a GLModel.  Depending on the input file format, not all fields may be defined.
  * @typedef ParserOptionsSpec
  * @prop {boolean} duplicateAssemblyAtoms- Set to true if you wish to diplicate assembly atoms otherwise false ; supported by all
  * @prop {boolean} frames - true if you want to add to a new frame and false otherwise ; supported by all
  * @prop {object} vibrate - object specifying the vibration behavior ; supported by all
  * @prop {number} vibrate.frames - number of frames to be created, default to 10 ; supported by all
  * @prop {number} vibrate.amplitude -amplitude of distortion, default to 1 (full) ; supported by all
  * @prop {boolean} multimodel - specifies weather or not multiple models are being defined ; supported by xyz,sdf, or mol2
  * @prop {boolean} onemol -specifies weather or not the model is of one molecule ; Supported by xyz , sdf , mol2
  * @prop {boolean} keepH - do not strip hydrogens ; supported by sdf,mol2
  * @prop {object} parseStyle - used to define ChemDoodle styles ; upported by cdjson
  * @prop {boolean} doAssembly - boolean dictating weather or not to do assembly ; supported by mcif
  * @prop {boolean} noSecondaryStructure - boolean dictating the presence of a secondary structure ; supported by pdb
  * @prob {boolean} noComputeSecondaryStructure - do not compute ss ; supported by pdb
  * @prob {string} altLoc -which alternate location to select, if present; '*' to load all ; supported by pdb
  * @prob {number} assemblyIndex - index of the assembly in symmetry ; supported by mmtf
  */

/**
*3 dimensional vector 
*@typedef Vector3
*@prop {number} x - x coordinate
*@prop {number} y - y coordinate
*@prop {number} z - z coordinate


*/

/**
 * Atom selection object. Used to specify what atoms should be selected.  Can include
 * any field from {@link AtomSpec} in which case atoms must equal the specified value.  
 * All fields must match for the selection to hold. If values
 * are provided as a list, then only one value of the list must match.
 * 
 * @typedef AtomSelectionSpec
 * @prop {AtomSpec} ... - any field from {@link AtomSpec}, values may be singletons or lists
 * @prop {GLModel} model - a single model or list of models from which atoms should be selected
 * @prop {number} bonds - overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms
 * @prop {function} predicate - user supplied function that gets passed an {AtomSpec} and should return true if the atom should be selected
 * @prop {boolean} invert - if set, inverts the meaning of the selection
 * @prop {boolean} byres - if set, expands the selection to include all atoms of any residue that has any atom selected
 * @prop {number} expand - expands the selection to include all atoms within a given distance from the selection
 * @prop {WithinSelectionSpec} within - intersects the selection with the set of atoms within a given distance from another selection
 
 * @example
 * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
                  
                  viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
                  viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                           borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                           position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
                  viewer.setStyle({chain:'A'},{cross:{hidden:true}});
                  viewer.setStyle({chain:'B'},{cross:{hidden:false,
                                                      colorscheme:'greenCarbon'}});
                  viewer.setStyle({chain:'C'},{cross:{hidden:false,
                                                      radius:0.5}});
                  viewer.setStyle({chain:'D'},{cross:{hidden:false}});
                  viewer.setStyle({chain:'E'},{cross:{hidden:false,
                                                      color:'black'}});
                  
                  viewer.render();

                  
                });
 */

/**
 * Within selection object. Used to find the subset of an atom selection that is within
 * some distance from another atom selection. When added as a field of an {@link AtomSelectionSpec},
 * intersects the set of atoms in that selection with the set of atoms within a given
 * distance from the given {@link AtomSelectionSpec}.
 
 * @typedef WithinSelectionSpec
 * @example
 $3Dmol.download("pdb:2EJ0",viewer,{},function(){
                  
                  viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
                  viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                           borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                           position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
            
 * viewer.setStyle({chain: 'A', within:{distance: 10, sel:{chain: 'B'}}}, {sphere:{}}); 
                  
                  viewer.render();

                  
                });// stylizes atoms in chain A that are within 10 angstroms of an atom in chain B
 *
 * @prop {number} distance - the distance in angstroms away from the atom selection to include atoms in the parent selection
*  @prop {boolean} invert - if set, selects atoms not within distance range for intersection
 * @prop {AtomSelectionSpec} sel - the selection of atoms against which to measure the distance from the parent atom selection
 */



/** 
 * @typedef AtomStyleSpec
 * @prop {LineStyleSpec} line - draw bonds as lines
 * @prop {CrossStyleSpec} cross - draw atoms as crossed lines (aka stars)
 * @prop {StickStyleSpec} stick  - draw bonds as capped cylinders
 * @prop {SphereStyleSpec} sphere - draw atoms as spheres
 * @prop {CartoonStyleSpec} cartoon - draw cartoon representation of secondary structure
 */

/** 
 * @typedef SurfaceStyleSpec
 * @prop {number} opacity - sets the transparency: 0 to hide, 1 for fully opaque
 * @prop {ColorschemeSpec} colorscheme - element based coloring
 * @prop {ColorSpec} color - fixed coloring, overrides colorscheme
 * @prop {$3Dmol.VolumeData} voldata - volumetric data for vertex coloring
 * @prop {$3Dmol.Gradient} volscheme - coloring scheme for mapping volumetric data to vertex color
 * @prop {Object} map - specifies a numeric atom property (prop) and color mapping (scheme) such as {@link $3Dmol.Gradient.RWB}.  Deprecated, use colorscheme instead.
 * 
 * @example
 * var setStyles = function(volumedata){
                    var data = new $3Dmol.VolumeData(volumedata, "cube");
                    viewer.addSurface("VDW", {opacity:0.85, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)},{chain:'A'});
                    viewer.mapAtomProperties($3Dmol.applyPartialCharges);
                    viewer.addSurface($3Dmol.SurfaceType.SAS, {map:{prop:'partialCharge',scheme:new $3Dmol.Gradient.RWB(-.05,.05)}, opacity:1.0},{chain:'B'});
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: data, color:'red'},{chain:'C'});
                    viewer.addSurface($3Dmol.SurfaceType.SAS, {opacity:0.85,voldata: data, colorscheme:'greenCarbon'},{chain:'D'});
                 
              viewer.render();
              };
              $3Dmol.download("pdb:4DLN",viewer,{},function(){
                  
                  $.get("volData/1fas.cube",setStyles);
                });

 */

/** 
 * Isosurface style specification
 * @typedef IsoSurfaceSpec
 * @prop {number} isoval - specifies the isovalue to draw surface at
 * @propr {boolean} voxel - if true uses voxel style rendering
 * @prop {ColorSpec} color - solid color
 * @prop {number} opacity - transparency, between 0 and 1
 * @prop {boolean} wireframe - draw as wireframe, not surface
 * @prop {number} linewidth - width of line for wireframe rendering **No longer supported by most browsers**
 * @prop {number} smoothness - amount to smooth surface (default 1)
 * @prop {AtomSelectionSpec} sel - selection around which to show data
 * @prop {list} coords - coordinates around which to include data
 * @prop {number} seldist - distance around selection/coords to include data [default = 2.0]
 * @prop {boolean} clickable - if true, user can click on object to trigger callback
 * @prop {function} callback - function to call on click 
 */

/** 
 * GLShape style specification
 * @typedef ShapeSpec
 * @prop {ColorSpec} color - solid color
 * @prop {number} alpha - transparency
 * @prop {boolean} wireframe - draw as wireframe, not surface
 * @prop {number} linewidth - width of line for wireframe rendering **No longer supported by most browsers**
 * @prop {boolean} clickable - if true, user can click on object to trigger callback
 * @prop {function} callback - function to call on click 
 */


/**
 * Specification for adding custom shape. Extends {@link ShapeSpec}.
 * @typedef CustomShapeSpec
 * @augments ShapeSpec
 * @prop {Array.<$3Dmol.Vector3>} vertexArr - List of vertex positions
 * @prop {Array.<$3Dmol.Vector3>} normalArr - List of normal vectors for each vertex
 * @prop {Array.<number>} faceArr - List of triangles to build the custom shape. Each triangle is defined by the indices of 3 vertices in vertexArr, so the array length should be 3 times the number of faces.
 * @prop {ColorSpec | Array.<ColorSpec>} color - Either a single color for the whole object or an array specifying the color at each vertex.
 */

/**
 * Sphere shape specification. Extends {@link ShapeSpec}  
 * 
 * @typedef SphereSpec   
 * @prop {$3Dmol.Vector3} center
 * @prop {number} radius
 * 
 */


/**
 * Arrow shape specification.  Extends {@link ShapeSpec}  
 * @typedef ArrowSpec
 * @prop {$3Dmol.Vector3} start
 * @prop {$3Dmol.Vector3} end
 * @prop {number} radius
 * @prop {number} radiusRatio - ratio of arrow base to cylinder (1.618034 default)
 * @prop {number} mid - relative position of arrow base (0.618034 default)
 */


/**
 * Cylinder shape specification.  Extends {@link ShapeSpec}  
 * @typedef CylinderSpec
 * @prop {$3Dmol.Vector3} start
 * @prop {$3Dmol.Vector3} end
 * @prop {number} radius
 * @prop {boolean} fromCap
 * @prop {boolean} toCap
 */

/**
 * Line shape specification.  Extends {@link ShapeSpec}  (but defaults to wireframe)
 * @typedef LineSpec
 * @prop {$3Dmol.Vector3} start
 * @prop {$3Dmol.Vector3} end
 */


/**
* File formats supported by 3Dmol.js
* @typedef FileFormats
* @prop cdjson,json  Chemical JSON format
* @prop cube Gaussian cube format
* @prop gro  Gromacs topology format, need to add coordinates to resulting model.
* @prop mcif,cif Crystallographic Information File, the successor to PDB that makes you miss the PDB file format
* @prop mmtf Macromolecular Transmission Format, the successor to PDB that is totally awesome
* @prop mol2 Sybyl Mol2 format 
* @prop pdb The venerable Protein Data Bank format
* @prop pqr Like PDB but with partial charges which are read into the partialcharge atom property
* @prop prmtop Amber topology file, must add coordinates
* @prop sdf MDL MOL format, supports muliple models and meta data
* @prop vasp VASP format (CONTCAR, POSCAR)
* @prop xyz XYZ cartesian coordinates format
*/
	
});
