// Specifications for various object types used in 3Dmol.js
// This is primarily for documentation
(function() {
/**
 * GLViewer input specification
 * @typedef ViewerSpec
 * @prop {function} callback - Callback function to be immediately executed on this viewer
 * @prop {Object} defaultcolors - Object defining default atom colors as atom => color property value pairs for all models within this viewer
 * @prop {boolean} nomouse - Whether to disable disable handling of mouse events.  
                If you want to use your own mouse handlers, set this then bind your handlers to the canvas object.  
                The default 3Dmol.js handlers are available for use: 
                'mousedown touchstart': viewer._handleMouseDown,
                'DOMMouseScroll mousewheel': viewer._handleMouseScroll
                'mousemove touchmove': viewer._handleMouseMove                
 * @prop {string} backgroundColor - Color of the canvas background
 * @prop {number} backgroundAlpha - Alpha transparency of canvas background
 * @prop {number} camerax
 * @prop {number} hoverDuration
 * @prop {string} id - id of the canvas
 * @prop {number} cartoonQuality - default 5
 * @prop {number} row
 * @prop {number} col
 * @prop {number} rows
 * @prop {number} cols
 * @prop canvas
 * @prop viewers
 * @prop minimumZoomToDistance
 * @prop lowerZoomLimit
 * @prop upperZoomLimit
 * @prop {boolean} antialias
 * @prop {boolean} control_all
 * @prop {boolean} orthographic
 * @prop {boolean} disableFog - Disable fog, default to false
 */

/**
 * Grid GLViewer input specification
 * @typedef ViewerGridSpec
 * @prop {number} rows - number of rows in grid
 * @prop {number} cols - number of columns in grid
 * @prop {boolean} control_all - if true, mouse events are linked
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
 * @prop {AtomSpec} ... - any field from {@link AtomSpec}, values may be singletons or lists. Integer numerical ranges are supported as strings.
 * @prop {GLModel} model - a single model or list of models from which atoms should be selected.  Can also specify by numerical creation order.  Reverse indexing is allowed (-1 specifies last added model).
 * @prop {number} bonds - overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms
 * @prop {function} predicate - user supplied function that gets passed an {AtomSpec} and should return true if the atom should be selected
 * @prop {boolean} invert - if set, inverts the meaning of the selection
 * @prop {boolean} byres - if set, expands the selection to include all atoms of any residue that has any atom selected
 * @prop {number} expand - expands the selection to include all atoms within a given distance from the selection
 * @prop {WithinSelectionSpec} within - intersects the selection with the set of atoms within a given distance from another selection
 * @prop {AtomSelectionSpec[]} and - take the intersection of the provided lists of {AtomSelectionSpec}s
 * @prop {AtomSelectionSpec[]} or - take the union of the provided lists of {AtomSelectionSpec}s
 * @prop {AtomSelectionSpec} not - take the inverse of the provided {AtomSelectionSpec}

 * @example
 * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
                  viewer.setStyle({chain:'B'},{cartoon:{color:'spectrum'}});
                  viewer.setStyle({chain:'B',invert:true},{cartoon:{}});
                  viewer.setStyle({bonds: 0},{sphere:{radius:0.5}}); //water molecules
                  viewer.setStyle({resn:'PMP',byres:true,expand:5},{stick:{colorscheme:"greenCarbon"}});
                  viewer.setStyle({resi:["91-95","42-50"]},{cartoon:{color:"green",thickness:1.0}});
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

                  viewer.setStyle({chain: 'A', within:{distance: 10, sel:{chain: 'B'}}}, {sphere:{}});
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
 * @prop {ClickSphereStyleSpec} clicksphere - invisible style for click handling only
 */

/**
 * @typedef SurfaceStyleSpec
 * @prop {number} opacity - sets the transparency: 0 to hide, 1 for fully opaque
 * @prop {ColorschemeSpec} colorscheme - element based coloring
 * @prop {ColorSpec} color - fixed coloring, overrides colorscheme
 * @prop {$3Dmol.VolumeData} voldata - volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified
 * @prop {$3Dmol.Gradient} volscheme - coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields.
 * @prop {string} volformat - format of voldata if not a $3Dmol.VolumeData object
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

                  $.get("data/1fas.cube",setStyles);
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
 * @prop {list} coords - coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates
 * @prop {number} seldist - distance around coords to include data [default = 2.0]
 * @prop {$3Dmol.VolumeData} voldata - volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified
 * @prop {$3Dmol.Gradient} volscheme - coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields.
 * @prop {string} volformat - format of voldata if not a $3Dmol.VolumeData object 
 * @prop {boolean} clickable - if true, user can click on object to trigger callback
 * @prop {function} callback - function to call on click
 */

/**
 * VolumetricRenderer style specification
 * @typedef VolumetricRendererSpec
 * @prop {list} transferfn - list of objects containing @color, @opacity and @value properties to specify color per voxel data value
 * @prop {number} subsamples - number of times to sample each voxel approximately (default 5)
 */
    //TODO: implement pruning of data, should start with box
    // @prop {list} coords - coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates
    // @prop {number} seldist - distance around coords to include data [default = 2.0]

/**
 * GLShape style specification
 * @typedef ShapeSpec
 * @prop {ColorSpec} color - solid color
 * @prop {number} alpha - transparency
 * @prop {boolean} wireframe - draw as wireframe, not surface
 * @prop {boolean} hidden - if true, do not display object
 * @prop {number} linewidth - width of line for wireframe rendering **No longer supported by most browsers**
 * @prop {boolean} clickable - if true, user can click on object to trigger callback
 * @prop {function} callback - function to call on click
 * @prop {number} frame - if set, only display in this frame of an animation
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
 * @typedef SphereShapeSpec
 * @prop {$3Dmol.Vector3} center
 * @prop {number} radius
 *
 */

/**
 * Box shape specification. Extends {@link ShapeSpec}
 *
 * @typedef BoxSpec
 * @prop {$3Dmol.Vector3} corner - bottom corner of box
 * @prop {$3Dmol.Vector3} center - alternative to corner: center of box
 * @prop {Object} dimensions - {w:width, h:height, d:depth}; can be either scalars or vectors (for not-axis aligned boxes)
 *
 */


/**
 * Arrow shape specification.  Extends {@link ShapeSpec}
 * @typedef ArrowSpec
 * @prop {$3Dmol.Vector3} start
 * @prop {$3Dmol.Vector3} end
 * @prop {number} radius
 * @prop {ColorSpec} color
 * @prop {boolean} hidden
 * @prop {number} radiusRatio - ratio of arrow base to cylinder (1.618034 default)
 * @prop {number} mid - relative position of arrow base (0.618034 default)
 * @prop {number} midpos - position of arrow base in length units, if negative positioned from end instead of start.  Overrides mid.
 */


/**
 * Cylinder shape specification.  Extends {@link ShapeSpec}
 * @typedef CylinderSpec
 * @prop {$3Dmol.Vector3} start
 * @prop {$3Dmol.Vector3} end
 * @prop {number} radius
 * @prop {$3Dmol.CAP} fromCap - 0 for none, 1 for flat, 2 for round
 * @prop {$3Dmol.CAP} toCap - 0 for none, 1 for flat, 2 for round
 * @prop {boolean} dashed
 */

/**
 * Curve shape specification.  Extends {@link ShapeSpec}
 * @typedef CurveSpec
 * @prop {$3Dmol.Vector3} points - list of (x,y,z) points to interpolate between to make curve
 * @prop {number} smooth - amount of interpolation
 * @prop {number} radius
 * @prop {boolean} fromArrow - if an arrow should be drawn at the start
 * @prop {boolean} toArrow - if an arrow should be drawn at the end
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
* @prop sdf MDL MOL format, supports multiple models and meta data
* @prop vasp VASP format (CONTCAR, POSCAR)
* @prop xyz XYZ cartesian coordinates format
*/

})();
