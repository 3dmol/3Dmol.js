/**
 * All of the functionality of $3Dmol.js is contained within the
 * $3Dmol global namespace
 */
declare namespace $3Dmol {

  type Matrix3 = any;
  type Matrix4 = any;
  type Object3D = {};

  /**
* Renderable labels
* @param tag - Label text
* @param parameters - Label style and font specifications
*/
  class Label {
    constructor(tag: string, parameters: LabelSpec);
  }

  /**
* GLModel represents a group of related atoms
* @param [defaultcolors] - Object defining default atom colors as atom => color property value pairs
*/
  class GLModel {
    constructor(mid?: number, defaultcolors?: any);
    /**
     * Return object representing internal state of
     * the model appropriate for passing to setInternalState
     */
    getInternalState(): void;
    /**
     * Overwrite the internal model state with the passed state.
     */
    setInternalState(): void;
    /**
     * Returns crystallographic information if present.
     */
    getCrystData(): void;
    /**
     * Set crystallographic information using three angles and three lengths
     * @param a - length of unit cell side
     * @param b - length of unit cell side
     * @param c - length of unit cell side
     * @param alpha - unit cell angle in degrees (default 90)
     * @param beta - unit cell angle in degrees (default 90)
     * @param gamma - unit cell angle in degrees (default 90)
     */
    setCrystData(a: number, b: number, c: number, alpha: number, beta: number, gamma: number): void;
    /**
     * Set the crystallographic matrix to the given matrix.
     *
     * This function removes `a`, `b`, `c`, `alpha`, `beta`, `gamma` from
     * the crystal data.
     * @param matrix - unit cell matrix
     */
    setCrystMatrix(matrix: $3Dmol.Matrix3): void;
    /**
     * Returns list of rotational/translational matrices if there is BIOMT data
     * Otherwise returns a list of just the ID matrix
     */
    getSymmetries(): $3Dmol.Matrix4[];
    /**
     * Sets symmetries based on specified matrices in list
     */
    setSymmetries(list: $3Dmol.Matrix4[]): void;
    /**
     * Returns model id number
     * @returns Model ID
     */
    getID(): number;
    /**
     * Returns model's frames property, a list of atom lists
     */
    getNumFrames(): number;
    /**
     * Sets model's atomlist to specified frame
     * Sets to last frame if framenum out of range
     * @param framenum - model's atoms are set to this index in frames list
     */
    setFrame(framenum: number): Promise<void>;
    /**
     * Add atoms as frames of model
     * @param atom - atoms to be added
     */
    addFrame(atom: AtomSpec): void;
    /**
     * If model atoms have dx, dy, dz properties (in some xyz files), vibrate populates the model's frame property based on parameters.
     * Model can then be animated
     * @example
     * $3Dmol.download("pdb:4UAA",viewer,{},function(){
     *             viewer.setStyle({},{stick:{}});
     *             viewer.vibrate(10, 1);
     *             viewer.animate({loop: "forward",reps: 1});
     *
     *             viewer.zoomTo();
     *                   viewer.render();
     *               });
     * @param numFrames - number of frames to be created, default to 10
     * @param amplitude - amplitude of distortion, default to 1 (full)
     * @param bothWays - if true, extend both in positive and negative directions by numFrames
     * @param viewer - required if arrowSpec is provided
     * @param arrowSpec - specification for drawing animated arrows. If color isn't specified, atom color (sphere, stick, line preference) is used.
     */
    vibrate(numFrames: number, amplitude: number, bothWays: boolean, viewer: GLViewer, arrowSpec: ArrowSpec): void;
    /**
     * add atoms to this model from molecular data string
     * @param data - atom structure file input data string, for gzipped input use ArrayBuffer
     * @param format - input file string format (e.g 'pdb', 'sdf', 'sdf.gz', etc.)
     * @param options - format dependent options. Attributes depend on the input format
     */
    addMolData(data: string | ArrayBuffer, format: string, options: ParserOptionsSpec): void;
    /**
     * given a selection specification, return true if atom is selected.
     * Does not support context-aware selectors like expand/within/byres.
     */
    atomIsSelected(atom: AtomSpec, sel: AtomSelectionSpec): boolean;
    /**
     * return list of atoms selected by sel, this is specific to glmodel
     * @example
     * $3Dmol.download("pdb:4wwy",viewer,{},function(){
     *                   var atoms = viewer.selectedAtoms({chain:'A'});
     *                   for(var i = 0, n = atoms.length; i < n; i++) {
     *                      atoms[i].b = 0.0;
     *                   }
     *                   viewer.setStyle({cartoon:{colorscheme:{prop:'b',gradient: 'roygb',min:0,max:30}}});
     *                   viewer.render();
     *               });
     */
    selectedAtoms(sel: AtomSelectionSpec): object[];
    /**
     * Add list of new atoms to model.  Adjusts bonds appropriately.
     * @example
     * var atoms = [{elem: 'C', x: 0, y: 0, z: 0, bonds: [1,2], bondOrder: [1,2]}, {elem: 'O', x: -1.5, y: 0, z: 0, bonds: [0]},{elem: 'O', x: 1.5, y: 0, z: 0, bonds: [0], bondOrder: [2]}];
     *
     *             viewer.setBackgroundColor(0xffffffff);
     *             var m = viewer.addModel();
     *             m.addAtoms(atoms);
     *             m.setStyle({},{stick:{}});
     *             viewer.zoomTo();
     *             viewer.render();
     */
    addAtoms(newatoms: any): void;
    /**
     * Remove specified atoms from model
     * @param badatoms - list of atoms
     */
    removeAtoms(badatoms: any): void;
    /**
     * Set atom style of selected atoms
     * @example
     * $3Dmol.download("pdb:4UB9",viewer,{},function(){
     *                   viewer.setBackgroundColor(0xffffffff);
     *
     *                   viewer.setStyle({chain:'A'},{line:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.setStyle({chain:'B'},{line:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.setStyle({chain:'C'},{cross:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.setStyle({chain:'D'},{cross:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.setStyle({chain:'E'},{cross:{radius:2.0,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.setStyle({chain:'F'},{stick:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.setStyle({chain:'G'},{stick:{radius:0.8,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.setStyle({chain:'H'},{stick:{singleBonds:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
     *                   viewer.render();
     *               });
     * @param add - if true, add to current style, don't replace
     */
    setStyle(sel: AtomSelectionSpec, style: AtomStyleSpec, add: boolean): void;
    /**
     * Set clickable and callback of selected atoms
     * @param sel - atom selection to apply clickable settings to
     * @param clickable - whether click-handling is enabled for the selection
     * @param callback - function called when an atom in the selection is clicked
     */
    setClickable(sel: AtomSelectionSpec, clickable: boolean, callback: (...params: any[]) => any): void;
    /**
     * Set hoverable and callback of selected atoms
     * @param sel - atom selection to apply hoverable settings to
     * @param hoverable - whether hover-handling is enabled for the selection
     * @param hover_callback - function called when an atom in the selection is hovered over
     * @param unhover_callback - function called when the mouse moves out of the hover area
     */
    setHoverable(sel: AtomSelectionSpec, hoverable: boolean, hover_callback: (...params: any[]) => any, unhover_callback: (...params: any[]) => any): void;
    /**
     * enable context menu of selected atoms
     * @param sel - atom selection to apply hoverable settings to
     * @param contextMenuEnabled - whether contextMenu-handling is enabled for the selection
     */
    enableContextMenu(sel: AtomSelectionSpec, contextMenuEnabled: boolean): void;
    /**
     * given a mapping from element to color, set atom colors
     */
    setColorByElement(sel: any, colors: any): void;
    static setColorByProperty(sel: any, prop: any, gradient: any): void;
    /**
     * @example
     * $3Dmol.download("pdb:4UAA",viewer,{},function(){
     *                   viewer.setBackgroundColor(0xffffffff);
     *                   var colorAsSnake = function(atom) {
     *                     return atom.resi % 2 ? 'white': 'green'
     *                   };
     *
     *                   viewer.setStyle( {}, { cartoon: {colorfunc: colorAsSnake }});
     *
     *                   viewer.render();
     *               });
     * @param sel - selection object
     * @param func - function to be used to set the color
     */
    setColorByFunction(sel: any, func: any): void;
    /**
     * Convert the model into an object in the format of a ChemDoodle JSON model.
     * @param whether - or not to include style information. Defaults to false.
     */
    toCDObject(whether: boolean): any;
    /**
     * manage the globj for this model in the possed modelGroup - if it has to be regenerated, remove and add
     * @param Object - options
     */
    globj(group: $3Dmol.Object3D, Object: any): void;
    /**
     * return a VRML string representation of the model.  Does not include VRML header information
     * @returns VRML
     */
    exportVRML(): any;
    /**
     * Remove any renderable mol object from scene
     */
    removegl(group: $3Dmol.Object3D): void;
    /**
     * @example
     * $3Dmol.download("pdb:3ucr",viewer,{},function(){
     *             viewer.setStyle({},{stick:{}});
     *             viewer.getModel().hide();
     *             viewer.render();
     *             });
     */
    hide(): void;
    /**
     * @example
     * $3Dmol.download("pdb:3ucr",viewer,{},function(){
     *             viewer.setStyle({},{stick:{}});
     *             viewer.getModel().hide();
     *             viewer.render(  )
     *             viewer.getModel().show()
     *             viewer.render();
     *             });
     */
    show(): void;
    /**
     * Create labels for atoms that show the value of the passed property.
     * @param prop - property name
     */
    addPropertyLabels(prop: string, sel: AtomSelectionSpec, viewer: $3Dmol.GLViewer, options: LabelSpec): void;
    /**
     * Create labels for residues of selected atoms.
     * Will create a single label at the center of mass of all atoms
     * with the same chain,resn, and resi.
     * @param byframe - if true, create labels for every individual frame, not just current; frames must be loaded already
     */
    addResLabels(sel: AtomSelectionSpec, viewer: $3Dmol.GLViewer, options: LabelSpec, byframe: boolean): void;
    /**
     * Set coordinates from remote trajectory file.
     * @param url - contains the url where mdsrv has been hosted
     * @param path - contains the path of the file (<root>/filename)
     */
    setCoordinatesFromURL(url: string, path: string): Promise<void>;
    /**
     * Set coordinates for the atoms from provided trajectory file.
     * @example
     * let m = viewer.addModel()  //create an empty model
     *          m.addAtoms([{x:0,y:0,z:0,elem:'C'},{x:2,y:0,z:0,elem:'C'}]) //provide a list of dictionaries representing the atoms
     *          viewer.setStyle({'sphere':{}})
     *          m.setCoordinates([[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [2.8888888359069824, 0.0, 0.0]], [[0.0, 0.0, 0.0], [3.777777671813965, 0.0, 0.0]], [[0.0, 0.0, 0.0], [4.666666507720947, 0.0, 0.0]], [[0.0, 0.0, 0.0], [5.55555534362793, 0.0, 0.0]], [[0.0, 0.0, 0.0], [6.44444465637207, 0.0, 0.0]], [[0.0, 0.0, 0.0], [7.333333492279053, 0.0, 0.0]], [[0.0, 0.0, 0.0], [8.222222328186035, 0.0, 0.0]], [[0.0, 0.0, 0.0], [9.11111068725586, 0.0, 0.0]], [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]],'array');
     *          viewer.animate({loop: "forward",reps: 1});
     *          viewer.zoomTo();
     *          viewer.zoom(0.5);
     *          viewer.render();
     * @param str - contains the data of the file
     * @param format - contains the format of the file (mdcrd, inpcrd, pdb, netcdf, or array).  Arrays should be TxNx3 where T is the number of timesteps and N the number of atoms.
     */
    setCoordinates(str: string, format: string): void;
    /**
     * add atomSpecs to validAtomSelectionSpecs
     * @param customAtomSpecs - array of strings that can be used as atomSelectionSpecs
     * this is to prevent the 'Unknown Selector x' message on the console for the strings passed
     */
    addAtomSpecs(customAtomSpecs: any[]): void;
  }




  /**
* WebGL-based 3Dmol.js viewer
* Note: The preferred method of instantiating a GLViewer is through {@link $3Dmol.createViewer}
* @param element - HTML element within which to create viewer
* @param config - Object containing optional configuration for the viewer
*/
  class GLViewer {
    constructor(element: any, config: ViewerSpec);
    /**
     * Return a list of objects that intersect that at the specified viewer position.
     * @param undefined - x position in screen coordinates
     * @param undefined - y position in screen coordinates
     * @param undefined - list of objects or selection object specifying what object to check for targeting
     */
    targetedObjects(): void;
    /**
     * Convert model coordinates to screen coordinates.
     * @param undefined - an object or list of objects with x,y,z attributes (e.g. an atom)
     * @returns - and object or list of {x: screenX, y: screenY}
     */
    modelToScreen(): any | Array<any>;
    /**
     * For a given screen (x,y) displacement return model displacement
     */
    screenOffsetToModel(): void;
    /**
     * Distance from screen coordinate to model coordinate assuming screen point
     * is projected to the same depth as model coordinate
     */
    screenToModelDistance(): void;
    /**
     * Set a callback to call when the view has potentially changed.
     */
    setViewChangeCallback(): void;
    /**
     * Set a callback to call when the view has potentially changed.
     */
    setStateChangeCallback(): void;
    /**
     * Return object representing internal state of
     * the viewer appropriate for passing to setInternalState
     */
    getInternalState(): void;
    /**
     * Overwrite internal state of the viewer with passed  object
     * which should come from getInternalState.
     */
    setInternalState(): void;
    /**
     * Set lower and upper limit stops for zoom.
     * @example
     * $.get("data/set1_122_complex.mol2", function(moldata) {
     *                 var m = viewer.addModel(moldata);
     *                 viewer.setStyle({stick:{colorscheme:"Jmol"}});
     *                 viewer.setZoomLimits(100,200);
     *                 viewer.zoomTo();
     *                 viewer.zoom(10); //will not zoom all the way
     *                 viewer.render();
     *             });
     * @param undefined - limit on zoom in (positive number).  Default 0.
     * @param undefined - limit on zoom out (positive number).  Default infinite.
     */
    setZoomLimits(): void;
    /**
     * Set camera parameters (distance to the origin and field of view)
     * @example
     * $.get("data/set1_122_complex.mol2", function(data) {
     *                 var m = viewer.addModel(data);
     *                 viewer.setStyle({stick:{}});
     *                 viewer.zoomTo();
     *                 viewer.setCameraParameters({ fov: 10 , z: 300 });
     *                 viewer.render();
     *             });
     * @param undefined - new camera parameters, with possible fields
     *                       being fov for the field of view, z for the
     *                       distance to the origin, and orthographic (boolean)
     *                       for kind of projection (default false).
     */
    setCameraParameters(): void;
    /**
     * Return image URI of viewer contents (base64 encoded).
     */
    pngURI(): void;
    /**
     * Return a promise that resolves to an animated PNG image URI of
     *          viewer contents (base64 encoded) for nframes of viewer changes.
     */
    apngURI(): Promise<string>;
    /**
     * Set the duration of the hover delay
     * @param [hoverDuration] - an optional parameter that denotes
     *            the duration of the hover delay (in milliseconds) before the hover action is called
     */
    setHoverDuration(hoverDuration?: number): void;
    /**
     * Change the viewer's container element
     * Also useful if the original container element was removed from the DOM.
     * @param element - Either HTML element or string identifier. Defaults to the element used to initialize the viewer.
     */
    setContainer(element: any | string): void;
    /**
     * Set the background color (default white)
     * @example
     * viewer.setBackgroundColor(0x000000);
     * @param hex - Hexcode specified background color, or standard color spec
     * @param a - Alpha level (default 1.0)
     */
    setBackgroundColor(hex: number, a: number): void;
    /**
     * Set view projection scheme.  Either orthographic or perspective.
     * Default is perspective.  Orthographic can also be enabled on viewer creation
     * by setting orthographic to true in the config object.
     * @example
     * viewer.setViewStyle({style:"outline"});
     *               $.get('data/1fas.pqr', function(data){
     *                   viewer.addModel(data, "pqr");
     *                   $.get("data/1fas.cube",function(volumedata){
     *                       viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
     *                   });
     *                   viewer.zoomTo();
     *
     *                   viewer.setProjection("orthographic");
     *                   viewer.render(callback);
     *               });
     */
    setProjection(): void;
    /**
     * Set global view styles.
     * @example
     * viewer.setViewStyle({style:"outline"});
     *               $.get('data/1fas.pqr', function(data){
     *                   viewer.addModel(data, "pqr");
     *                   $.get("data/1fas.cube",function(volumedata){
     *                       viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
     *                   });
     *                   viewer.zoomTo();
     *                   viewer.render(callback);
     *               });
     */
    setViewStyle(): void;
    /**
     * Set viewer width
     * @param w - Width in pixels
     */
    setWidth(w: number): void;
    /**
     * Set viewer height
     * @param h - Height in pixels
     */
    setHeight(h: number): void;
    /**
     * Resize viewer according to containing HTML element's dimensions
     */
    resize(): void;
    /**
     * Return specified model
     * @example
     * // Retrieve reference to first GLModel added var m =
     *    $3Dmol.download("pdb:1UBQ",viewer,{},function(m1){
     *                   $3Dmol.download("pdb:1UBI", viewer,{}, function(m2) {
     *                     viewer.zoomTo();
     *                     m1.setStyle({cartoon: {color:'green'}});
     *                     //could use m2 here as well
     *                     viewer.getModel().setStyle({cartoon: {color:'blue'}});
     *                     viewer.render();
     *                 })
     *               });
     * @param [id = last model id] - Retrieve model with specified id
     */
    getModel(id?: number): $3Dmol.GLModel;
    /**
     * Continuously rotate a scene around the specified axis.
     *
     * Call `$3Dmol.GLViewer.spin(false)` to stop spinning.
     * @param [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around.
     *            Default "y".  View relative (rather than model relative) axes are prefixed with v.
     * @param [speed] - Speed multiplier for spinning the viewer. 1 is default and a negative
     *             value reverses the direction of the spin.
     */
    spin(axis?: string, speed?: number): void;
    /**
     * Rotate scene by angle degrees around axis
     * @example
     * $3Dmol.download('cid:4000', viewer, {}, function() {
     *       viewer.setStyle({stick:{}});
     *       viewer.zoomTo();
     *       viewer.rotate(90,'y',1);
     *       viewer.render(callback);
     *     });
     * @param [angle] - Angle, in degrees, to rotate by.
     * @param [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around.
     *            Default "y".  View relative (rather than model relative) axes are prefixed with v.
     *            Axis can also be specified as a vector.
     * @param [animationDuration] - an optional parameter that denotes
     *            the duration of the rotation animation. Default 0 (no animation)
     * @param [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     */
    rotate(angle?: number, axis?: string, animationDuration?: number, fixedPath?: boolean): void;
    /**
     * Returns an array representing the current viewpoint.
     * Translation, zoom, and rotation quaternion.
     * @returns [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y, q.z, q.w ]
     */
    getView(): number[];
    /**
     * Sets the view to the specified translation, zoom, and rotation.
     * @param arg - Array formatted identically to the return value of getView
     */
    setView(arg: number[]): void;
    /**
     * Render current state of viewer, after
     * adding/removing models, applying styles, etc.
     */
    render(): void;
    /**
     * return list of atoms selected by sel
     */
    selectedAtoms(sel: AtomSelectionSpec): object[];
    /**
     * Return pdb output of selected atoms (if atoms from pdb input)
     * @param [sel] - Selection specification specifying model and atom properties to select.  Default: all atoms in viewer
     * @returns PDB string of selected atoms
     */
    pdbData(sel?: any): string;
    /**
     * Zoom current view by a constant factor
     * @example
     * $.get('data/4csv.pdb', function(data) {
     *       viewer.addModel(data,'pdb');
     *       viewer.setStyle({cartoon:{},stick:{}});
     *       viewer.zoomTo()
     *       viewer.zoom(2,1000);
     *       viewer.render();
     *     });
     * @param [factor] - Magnification factor. Values greater than 1
     *            will zoom in, less than one will zoom out. Default 2.
     * @param [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation
     */
    zoom(factor?: number, animationDuration?: number, fixedPath?: boolean): void;
    /**
     * Translate current view by x,y screen coordinates
     * This pans the camera rather than translating the model.
     * @example
     * $.get('data/4csv.pdb', function(data) {
     *       viewer.addModel(data,'pdb');
     *       viewer.setStyle({cartoon:{},stick:{}});
     *       viewer.zoomTo();
     *       viewer.translate(200,50);
     *       viewer.rotate(90,'z');
     *       viewer.render(callback);
     *     });
     * @param x - Relative change in view coordinates of camera
     * @param y - Relative change in view coordinates of camera
     * @param [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     */
    translate(x: number, y: number, animationDuration?: number, fixedPath?: boolean): void;
    /**
     * Translate current models by x,y screen coordinates
     * This translates the models relative to the current view. It does
     * not change the center of rotation.
     * @example
     * $.get('data/4csv.pdb', function(data) {
     *       viewer.addModel(data,'pdb');
     *       viewer.setStyle({cartoon:{},stick:{}});
     *       viewer.zoomTo();
     *       viewer.translateScene(200,50);
     *       viewer.rotate(90,'z'); // will no longer be around model center
     *       viewer.render(callback);
     *     });
     * @param x - Relative change in x screen coordinate
     * @param y - Relative change in y screen coordinate
     * @param [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     */
    translateScene(x: number, y: number, animationDuration?: number, fixedPath?: boolean): void;
    /**
     * Adjust slab to fully enclose selection (default everything).
     * @param [sel] - Selection specification specifying model and atom
     *            properties to select. Default: all atoms in viewer
     */
    fitSlab(sel?: any): void;
    /**
     * Re-center the viewer around the provided selection (unlike zoomTo, does not zoom).
     * @example
     * // if the user were to pass the animationDuration value to
     *           // the function like so viewer.zoomTo({resn:'STI'},1000);
     *         //   the program would center on resn 'STI' over the course
     *         //   of 1 second(1000 milleseconds).
     *  // Reposition to centroid of all atoms of all models in this
     * //viewer glviewer.center();
     *     $.get('data/4csv.pdb', function(data) {
     *       viewer.addModel(data,'pdb');
     *       viewer.setStyle({cartoon:{},stick:{}});
     *       viewer.center();
     *       viewer.render(callback);
     *     });
     * @param [sel] - Selection specification specifying model and atom
     *            properties to select. Default: all atoms in viewer
     * @param [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     */
    center(sel?: any, animationDuration?: number, fixedPath?: boolean): void;
    /**
     * Zoom to center of atom selection.  The slab will be set appropriately for
     * the selection, unless an empty selection is provided, in which case there will be no slab.
     * @example
     * $.get('data/1fas.pqr', function(data){
     *                   viewer.addModel(data, "pqr");
     *                   $.get("data/1fas.cube",function(volumedata){
     *                       viewer.addSurface($3Dmol.SurfaceType.VDW, {
     *                           opacity:0.85,
     *                           voldata: new $3Dmol.VolumeData(volumedata, "cube"),
     *                           volscheme: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'charge'))
     *                       },{});
     *
     *                   viewer.render();
     *                   });
     *                   viewer.zoomTo();
     *                 });
     * @param [sel] - Selection specification specifying model and atom
     *            properties to select. Default: all atoms in viewer
     * @param [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     */
    zoomTo(sel?: any, animationDuration?: number, fixedPath?: boolean): void;
    /**
     * Set slab of view (contents outside of slab are clipped).
     * Must call render to update.
     * @param near - near clipping plane distance
     * @param far - far clipping plane distance
     */
    setSlab(near: number, far: number): void;
    /**
     * Get slab of view (contents outside of slab are clipped).
     * @property near - near clipping plane distance
     * @property far - far clipping plane distance
     */
    getSlab(): {
      near: number;
      far: number;
    };
    /**
     * Add label to viewer
     * @example
     * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
     *
     *                   viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
     *                   viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
     *                                            borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
     *                                            position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
     *                   viewer.setStyle({chain:'A'},{cross:{hidden:true}});
     *                   viewer.setStyle({chain:'B'},{cross:{hidden:false,
     *                                                       linewidth:1.0,
     *                                                       colorscheme:'greenCarbon'}});
     *                   viewer.setStyle({chain:'C'},{cross:{hidden:false,
     *                                                       linewidth:1.0,
     *                                                       radius:0.5}});
     *                   viewer.setStyle({chain:'D'},{cross:{hidden:false,
     *                                                       linewidth:10.0}});
     *                   viewer.setStyle({chain:'E'},{cross:{hidden:false,
     *                                                       linewidth:1.0,
     *                                                       color:'black'}});
     *
     *                   viewer.render();
     *
     *
     *                 });
     * @param text - Label text
     * @param options - Label style specification
     * @param sel - Set position of label to center of this selection
     * @param noshow - if true, do not immediately display label - when adding multiple labels this is more efficient
     */
    addLabel(text: string, options: LabelSpec, sel: AtomSelection, noshow: boolean): $3Dmol.Label;
    /**
     * Add residue labels.  This will generate one label per a
     * residue within the selected atoms.  The label will be at the
     * centroid of the atoms and styled according to the passed style.
     * The label text will be [resn][resi]
     * @param byframe - if true, create labels for every individual frame, not just current
     *
     * * @example
     *              $3Dmol.download("mmtf:2ll5",viewer,{},function(){
     *                   viewer.setStyle({stick:{radius:0.15},cartoon:{}});
     *                   viewer.addResLabels({hetflag:false}, {font: 'Arial', fontColor:'black',showBackground:false, screenOffset: {x:0,y:0}});
     *                   viewer.zoomTo();
     *                   viewer.render();
     *                 });
     */
    addResLabels(sel: any, style: any, byframe: boolean): void;
    /**
     * Add property labels.  This will generate one label per a selected
     * atom at the atom's coordinates with the property value as the label text.
     * @param prop - property name
     * @param style - * @example
     *              $3Dmol.download("cid:5291",viewer,{},function(){
     *                   viewer.setStyle({stick: {radius:.2}});
     *                   viewer.addPropertyLabels("index",{not:{elem:'H'}}, {fontColor:'black',font: 'sans-serif', fontSize: 28, showBackground:false,alignment:'center'});
     *                   viewer.zoomTo();
     *                   viewer.render();
     *                 });
     */
    addPropertyLabels(prop: string, sel: any, style: any): void;
    /**
     * Remove label from viewer
     * @example
     * // Remove labels created in
     *          $3Dmol.download("pdb:2EJ0",viewer,{},function(){
     *                   var toremove = viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
     *                   viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
     *                                            borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
     *                                            position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
     *                   viewer.removeLabel(toremove);
     *                   viewer.render();
     *
     *
     *                 });
     * @param label - $3Dmol label
     */
    removeLabel(label: $3Dmol.Label): void;
    /**
     * Remove all labels from viewer
     * @example
     * $3Dmol.download("pdb:1ubq",viewer,{},function(){
     *
     *                viewer.addResLabels();
     *                viewer.setStyle({},{stick:{}});
     *                viewer.render( ); //show labels
     *
     *                viewer.removeAllLabels();
     *                viewer.render(); //hide labels
     *         });
     */
    removeAllLabels(): void;
    /**
     * Modify existing label's style
     * @param label - $3Dmol label
     * @param stylespec - Label style specification
     */
    setLabelStyle(label: $3Dmol.Label, stylespec: any): $3Dmol.Label;
    /**
     * Modify existing label's text
     * @param label - $3Dmol label
     * @param text - Label text
     */
    setLabelText(label: $3Dmol.Label, text: string): $3Dmol.Label;
    /**
     * Add shape object to viewer
     * @param shapeSpec - style specification for label
     */
    addShape(shapeSpec: ShapeSpec): $3Dmol.GLShape;
    /**
     * Remove shape object from viewer
     * @param shape - Reference to shape object to remove
     */
    removeShape(shape: $3Dmol.GLShape): void;
    /**
     * Remove all shape objects from viewer
     */
    removeAllShapes(): void;
    /**
     * Create and add sphere shape. This method provides a shorthand
     * way to create a spherical shape object
     * @example
     * viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
     *
     *          viewer.render();
     * @param spec - Sphere shape style specification
     */
    addSphere(spec: SphereShapeSpec): $3Dmol.GLShape;
    /**
     * Create and add box shape. This method provides a shorthand
     * way to create a box shape object
     * @example
     * viewer.addLine({color:'red',start:{x:0,y:0,z:0},end:{x:5,y:0,z:0}});
     *          viewer.addLine({color:'blue',start:{x:0,y:0,z:0},end:{x:0,y:5,z:0}});
     *          viewer.addLine({color:'green',start:{x:0,y:0,z:0},end:{x:0,y:0,z:5}});
     *
     *          viewer.addBox({center:{x:0,y:0,z:0},dimensions: {w:3,h:4,d:2},color:'magenta'});
     *          viewer.zoomTo();
     *          viewer.rotate(45, {x:1,y:1,z:1});
     *          viewer.render();
     * @param spec - Box shape style specification
     */
    addBox(spec: BoxSpec): $3Dmol.GLShape;
    /**
     * Create and add arrow shape
     * @example
     * $3Dmol.download("pdb:4DM7",viewer,{},function(){
     *
     *                   viewer.setBackgroundColor(0xffffffff);
     *                   viewer.addArrow({
     *                       start: {x:-10.0, y:0.0, z:0.0},
     *                       end: {x:0.0, y:-10.0, z:0.0},
     *                       radius: 1.0,
     *                       radiusRadio:1.0,
     *                       mid:1.0,
     *                       clickable:true,
     *                       callback:function(){
     *                           this.color.setHex(0xFF0000FF);
     *                           viewer.render( );
     *                       }
     *                   });
     *                   viewer.render();
     *                 });
     * @param spec - Style specification
     */
    addArrow(spec: ArrowSpec): $3Dmol.GLShape;
    /**
     * Create and add cylinder shape
     * @example
     * viewer.setBackgroundColor(0xffffffff);
     *               viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
     *                                   end:{x:10.0,y:0.0,z:0.0},
     *                                   radius:1.0,
     *                                   fromCap:1,
     *                                   toCap:2,
     *                                   color:'red',
     *                                   hoverable:true,
     *                                   clickable:true,
     *                                   callback:function(){ this.color.setHex(0x00FFFF00);viewer.render( );},
     *                                   hover_callback: function(){ viewer.render( );},
     *                                   unhover_callback: function(){ this.color.setHex(0xFF000000);viewer.render( );}
     *                                  });
     *               viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
     *                                   end:{x:0.0,y:10.0,z:0.0},
     *                                   radius:0.5,
     *                                   fromCap:false,
     *                                   toCap:true,
     *                                   color:'teal'});
     *               viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
     *                                   end:{x:20.0,y:0.0,z:0.0},
     *                                   radius:1.0,
     *                                   color:'black',
     *                                   fromCap:false,
     *                                   toCap:false});
     *               viewer.render();
     * @param spec - Style specification
     */
    addCylinder(spec: CylinderSpec): $3Dmol.GLShape;
    /**
     * Create and add Curve shape
     * @example
     * viewer.addCurve({points: [{x:0.0,y:0.0,z:0.0}, {x:5.0,y:3.0,z:0.0}, {x:5.0,y:7.0,z:0.0}, {x:0.0,y:10.0,z:0.0}],
     *                                   radius:0.5,
     *                                   smooth: 10,
     *                                   fromArrow:false,
     *                                   toArrow: true,
     *                                   color:'orange',
     *                                   });
     *               viewer.addCurve({points: [{x:-1,y:0.0,z:0.0}, {x:-5.0,y:5.0,z:0.0}, {x:-2,y:10.0,z:0.0}],
     *                                   radius:1,
     *                                   fromArrow:true,
     *                                   toArrow: false,
     *                                   color:'purple',
     *                                   });
     *               viewer.zoomTo();
     *               viewer.render();
     * @param spec - Style specification
     */
    addCurve(spec: CurveSpec): $3Dmol.GLShape;
    /**
     * Create and add line shape
     * @example
     * $3Dmol.download("pdb:2ABJ",viewer,{},function(){
     *
     *                   viewer.setViewStyle({style:"outline"});
     *                   viewer.setStyle({chain:'A'},{sphere:{hidden:true}});
     *                   viewer.setStyle({chain:'D'},{sphere:{radius:3.0}});
     *                   viewer.setStyle({chain:'G'},{sphere:{colorscheme:'greenCarbon'}});
     *                   viewer.setStyle({chain:'J'},{sphere:{color:'blue'}});
     *                   viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
     *                   viewer.render();
     *               });
     * @param spec - Style specification, can specify dashed, dashLength, and gapLength
     */
    addLine(spec: LineSpec): $3Dmol.GLShape;
    /**
     * Create and add unit cell visualization.
     * @example
     * $.get('data/1jpy.cif', function(data) {
     *                   let m = viewer.addModel(data);
     *                   viewer.addUnitCell(m, {box:{color:'purple'},alabel:'X',blabel:'Y',clabel:'Z',alabelstyle: {fontColor: 'black',backgroundColor:'white',inFront:true,fontSize:40},astyle:{color:'darkred', radius:5,midpos: -10}});
     *                   viewer.zoomTo();
     *                   viewer.render();
     *     });
     * @param model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
     * @param spec - visualization style
     */
    addUnitCell(model: GLModel, spec: UnitCellStyleSpec): void;
    /**
     * Remove unit cell visualization from model.
     * @example
     * $.get('data/icsd_200866.cif', function(data) {
     *                   let m = viewer.addModel(data);
     *                   viewer.setStyle({sphere:{}})
     *                   viewer.addUnitCell();
     *                   viewer.zoomTo();
     *                   viewer.removeUnitCell();
     *                   viewer.render();
     *             });
     * @param model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
     */
    removeUnitCell(model: GLModel): void;
    /**
     * Replicate atoms in model to form a super cell of the specified dimensions.
     * Original cell will be centered as much as possible.
     * @example
     * $.get('data/icsd_200866.cif', function(data) {
     *                   let m = viewer.addModel(data);
     *                   viewer.setStyle({sphere:{scale:.25}})
     *                   viewer.addUnitCell();
     *                   viewer.zoomTo();
     *                   viewer.replicateUnitCell(3,2,1,m);
     *                   viewer.render();
     *             });
     * @param A - number of times to replicate cell in X dimension.
     * @param B - number of times to replicate cell in Y dimension.  If absent, X value is used.
     * @param C - number of times to replicate cell in Z dimension.  If absent, Y value is used.
     * @param model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
     */
    replicateUnitCell(A: integer, B: integer, C: integer, model: GLModel): void;
    /**
     * Add custom shape component from user supplied function
     * @example
     * function triangle(viewer) {
     *     var vertices = [];
     *     var normals = [];
     *     var colors = [];
     *     var r = 20;
     *     //triangle
     *     vertices.push(new $3Dmol.Vector3(0,0,0));
     *     vertices.push(new $3Dmol.Vector3(r,0,0));
     *     vertices.push(new $3Dmol.Vector3(0,r,0));
     *
     *     normals.push(new $3Dmol.Vector3(0,0,1));
     *     normals.push(new $3Dmol.Vector3(0,0,1));
     *     normals.push(new $3Dmol.Vector3(0,0,1));
     *
     *     colors.push({r:1,g:0,b:0});
     *     colors.push({r:0,g:1,b:0});
     *     colors.push({r:0,g:0,b:1});
     *
     *     var faces = [ 0,1,2 ];
     *
     *     var spec = {vertexArr:vertices, normalArr: normals, faceArr:faces,color:colors};
     *     viewer.addCustom(spec);
     * }
     *             triangle(viewer);
     *             viewer.render();
     * @param spec - Style specification
     */
    addCustom(spec: CustomSpec): $3Dmol.GLShape;
    /**
     * Construct isosurface from volumetric data in gaussian cube format
     * @example
     * $.get('data/bohr.cube', function(data) {
     *
     *       viewer.addVolumetricData(data, "cube", {isoval: -0.01, color: "red", opacity: 0.95});
     *       viewer.setStyle({cartoon:{},stick:{}});
     *       viewer.zoomTo();
     *       viewer.render();
     *     });
     * @param data - Input file contents
     * @param format - Input file format
     * @param or - {VolumetricRenderSpec} spec - Shape style specification
     */
    addVolumetricData(data: string, format: string, or: IsoSurfaceSpec): $3Dmol.GLShape;
    /**
     * Construct isosurface from volumetric data.  This is more flexible
     * than addVolumetricData, but can not be used with py3Dmol.
     * @example
     * $.get('../test_structs/benzene-homo.cube', function(data){
     *                   var voldata = new $3Dmol.VolumeData(data, "cube");
     *                   viewer.addIsosurface(voldata, {isoval: 0.01,
     *                                                  color: "blue"});
     *                   viewer.addIsosurface(voldata, {isoval: -0.01,
     *                                                  color: "red"});
     *                   viewer.zoomTo();
     *                   viewer.render();
     *                 });
     * @param data - volumetric data
     * @param spec - Shape style specification
     */
    addIsosurface(data: $3Dmol.VolumeData, spec: IsoSurfaceSpec): $3Dmol.GLShape;
    /**
     * Create volumetric renderer for volumetricData
     * @param data - volumetric data
     * @param spec - specification of volumetric render
     */
    addVolumetricRender(data: $3Dmol.VolumeData, spec: VolumetricRenderSpec): $3Dmol.GLShape;
    /**
     * Return true if volumetric rendering is supported (WebGL 2.0 required)
     */
    hasVolumetricRender(): boolean;
    /**
     * Enable/disable fog for content far from the camera
     * @param fog - whether to enable or disable the fog
     */
    enableFog(fog: boolean): void;
    /**
     * Sets the atomlists of all models in the viewer to specified frame.
     * Shapes and labels can also be displayed by frame.
     * Sets to last frame if framenum out of range
     * @param framenum - fame index to use, starts at zero
     */
    setFrame(framenum: number): Promise;
    /**
     * Gets the current viewer frame.
     */
    getFrame(): void;
    /**
     * Returns the number of frames that the model with the most frames in the viewer has
     */
    getNumFrames(): number;
    /**
     * Animate all models in viewer from their respective frames
     * @param options - can specify interval (speed of animation), loop (direction
     * of looping, 'backward', 'forward' or 'backAndForth'), step interval between frames ('step'), and reps (numer of repetitions, 0 indicates infinite loop)
     */
    animate(options: any): void;
    /**
     * Stop animation of all models in viewer
     */
    stopAnimate(): void;
    /**
     * Return true if viewer is currently being animated, false otherwise
     */
    isAnimated(): boolean;
    /**
     * Create and add model to viewer, given molecular data and its format
     * @example
     * viewer.setViewStyle({style:"outline"});
     *               $.get('data/1fas.pqr', function(data){
     *                   viewer.addModel(data, "pqr");
     *                   $.get("data/1fas.cube",function(volumedata){
     *                       viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
     *
     *                   viewer.render();
     *                   });
     *                   viewer.zoomTo();
     *               });
     * @param data - Input data
     * @param format - Input format ('pdb', 'sdf', 'xyz', 'pqr', or 'mol2')
     * @param options - format dependent options. Attributes depend on the input file format.
     */
    addModel(data: string, format: string, options: ParserOptionsSpec): $3Dmol.GLModel;
    /**
     * Given multimodel file and its format, add atom data to the viewer as separate models
     * and return list of these models
     * @param data - Input data
     * @param format - Input format (see {@link FileFormats})
     */
    addModels(data: string, format: string): $3Dmol.GLModel[];
    /**
     * Create and add model to viewer. Given multimodel file and its format,
     * different atomlists are stored in model's frame
     * property and model's atoms are set to the 0th frame
     * @example
     * $.get('../test_structs/multiple2.xyz', function(data){
     *                   viewer.addModelsAsFrames(data, "xyz");
     *                   viewer.animate({loop: "forward",reps: 1});
     *                   viewer.setStyle({stick:{colorscheme:'magentaCarbon'}});
     *                   viewer.zoomTo();
     *                   viewer.render();
     *               });
     * @param data - Input data
     * @param format - Input format (see {@link FileFormats})
     */
    addModelsAsFrames(data: string, format: string): $3Dmol.GLModel;
    /**
     * Create and add model to viewer. Given multimodel file and its format,
     * all atoms are added to one model
     * @example
     * $.get('../test_structs/multiple.sdf', function(data){
     *                   viewer.addAsOneMolecule(data, "sdf");
     *                   viewer.zoomTo();
     *                   viewer.render();
     *               });
     * @param data - Input data
     * @param format - Input format (see {@link FileFormats})
     */
    addAsOneMolecule(data: string, format: string): $3Dmol.GLModel;
    /**
     * Delete specified model from viewer
     */
    removeModel(model: $3Dmol.GLModel): void;
    /**
     * Delete all existing models
     */
    removeAllModels(): void;
    /**
     * Export one or all of the loaded models into ChemDoodle compatible JSON.
     * @param includeStyles - Whether or not to include style information.
     * @param modelID - Optional parameter for which model to export. If left out, export all of them.
     */
    exportJSON(includeStyles: boolean, modelID: number): string;
    /**
     * return a VRML string representation of the scene.  Include VRML header information
     * @returns VRML
     */
    exportVRML(): any;
    /**
     * Create a new model from atoms specified by sel.
     * If extract, removes selected atoms from existing models
     * @param sel - Atom selection specification
     * @param [extract] - If true, remove selected atoms from existing models
     */
    createModelFrom(sel: any, extract?: boolean): $3Dmol.GLModel;
    /**
     * Set style properties to all selected atoms
     * @example
     * viewer.setBackgroundColor(0xffffffff);
     *        $3Dmol.download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
     *         m.setStyle({chain:'A'},{'cartoon':{color:'spectrum'}});
     *         m.setStyle({chain:'C'},{'cartoon':{style:'trace',color:'blue'}});
     *         m.setStyle({chain:'E'},{'cartoon':{tubes:true,arrows:true,color:'green',opacity:0.75}});
     *         m.setStyle({chain:'B'},{'cartoon':{color:'red',opacity:0.5}});
     *         m.setStyle({chain:'D'},{'cartoon':{style:'trace',color:'grey',opacity:0.75}});
     *         m.setStyle({chain:'F'},{'cartoon':{arrows:true,color:'white'}});
     *        // viewer.addStyle({chain:'B'},{line:{}});
     *        viewer.zoomTo();
     *        viewer.render();
     *     });
     * @param sel - Atom selection specification
     * @param style - Style spec to apply to specified atoms
     */
    setStyle(sel: AtomSelectionSpec, style: AtomStyleSpec): void;
    /**
     * Add style properties to all selected atoms
     * @example
     * $3Dmol.download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
     *        viewer.setStyle({cartoon:{}});
     *        //keep cartoon style, but show thick sticks for chain A
     *        viewer.addStyle({chain:'A'},{stick:{radius:.5,colorscheme:"magentaCarbon"}});
     *        viewer.zoomTo();
     *        viewer.render();
     *        });
     * @param sel - Atom selection specification
     * @param style - style spec to add to specified atoms
     */
    addStyle(sel: AtomSelectionSpec, style: AtomStyleSpec): void;
    /**
     * Set click-handling properties to all selected atomsthis.
     * @example
     * $3Dmol.download("cid:307900",viewer,{},function(){
     *
     *                    viewer.setStyle({},{sphere:{}});
     *                    viewer.setClickable({},true,function(atom,viewer,event,container) {
     *                        viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8});
     *                    });
     *                    viewer.render();
     *         });
     * @param sel - atom selection to apply clickable settings to
     * @param clickable - whether click-handling is enabled for the selection
     * @param callback - function called when an atom in the selection is clicked
     */
    setClickable(sel: AtomSelectionSpec, clickable: boolean, callback: (...params: any[]) => any): void;
    /**
     * Set hoverable and callback of selected atoms
     * @example
     * $3Dmol.download("pdb:1ubq",viewer,{},function(){
     *
     *                viewer.setHoverable({},true,function(atom,viewer,event,container) {
     *                    if(!atom.label) {
     *                     atom.label = viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
     *                    }
     *                },
     *                function(atom) {
     *                    if(atom.label) {
     *                     viewer.removeLabel(atom.label);
     *                     delete atom.label;
     *                    }
     *                 }
     *                );
     *                viewer.setStyle({},{stick:{}});
     *                viewer.render();
     *         });
     * @param sel - atom selection to apply hoverable settings to
     * @param hoverable - whether hover-handling is enabled for the selection
     * @param hover_callback - function called when an atom in the selection is hovered over
     * @param unhover_callback - function called when the mouse moves out of the hover area
     */
    setHoverable(sel: AtomSelectionSpec, hoverable: boolean, hover_callback: (...params: any[]) => any, unhover_callback: (...params: any[]) => any): void;
    /**
     * enable context menu and callback of selected atoms
     * @param sel - atom selection to apply hoverable settings to
     * @param contextMenuEnabled - whether contextMenu-handling is enabled for the selection
     */
    enableContextMenu(sel: AtomSelectionSpec, contextMenuEnabled: boolean): void;
    /**
     * If  atoms have dx, dy, dz properties (in some xyz files), vibrate populates each model's frame property based on parameters.
     * Models can then be animated
     * @param numFrames - number of frames to be created, default to 10
     * @param amplitude - amplitude of distortion, default to 1 (full)
     * @param bothWays - if true, extend both in positive and negative directions by numFrames
     * @param arrowSpec - specification for drawing animated arrows. If color isn't specified, atom color (sphere, stick, line preference) is used.
     */
    vibrate(numFrames: number, amplitude: number, bothWays: boolean, arrowSpec: ArrowSpec): void;
    setColorByProperty(sel: AtomSelectionSpec, prop: any, scheme: any): void;
    setColorByElement(sel: AtomSelectionSpec, colors: any): void;
    /**
     * Adds an explicit mesh as a surface object.
     * @returns surfid
     */
    addMesh(mesh: $3Dmol.Mesh, style: any): number;
    /**
     * Add surface representation to atoms
     * @param type - Surface type (VDW, MS, SAS, or SES)
     * @param style - optional style specification for surface material (e.g. for different coloring scheme, etc)
     * @param atomsel - Show surface for atoms in this selection
     * @param allsel - Use atoms in this selection to calculate surface; may be larger group than 'atomsel'
     * @param focus - Optionally begin rendering surface specified atoms
     * @param surfacecallback - function to be called after setting the surface
     * @returns promise - Returns a promise that ultimately resovles to the surfid.  Returns surfid immediately if surfacecallback is specified.  Returned promise has a [surfid, GLViewer, style, atomsel, allsel, focus] fields for immediate access.
     */
    addSurface(type: $3Dmol.SurfaceType | string, style: SurfaceStyleSpec, atomsel: AtomSelectionSpec, allsel: AtomSelectionSpec, focus: AtomSelectionSpec, surfacecallback: (...params: any[]) => any): Promise;
    /**
     * Set the surface material to something else, must render change
     * @example
     * $.get("data/9002806.cif",function(data){
     *             viewer.addModel(data);
     *             viewer.setStyle({stick:{}});
     *             let surf = viewer.addSurface("SAS");
     *             surf.then(function() {
     *                 viewer.setSurfaceMaterialStyle(surf.surfid, {color:'blue',opacity:0.5});
     *                 viewer.render();
     *                 });
     *            });
     * @param surf - Surface ID to apply changes to
     * @param style - new material style specification
     */
    setSurfaceMaterialStyle(surf: number, style: SurfaceStyleSpec): void;
    /**
     * Return surface object
     * @param surf - surface id
     */
    getSurface(surf: number): void;
    /**
     * Remove surface with given ID
     * @param surf - surface id
     */
    removeSurface(surf: number): void;
    /**
     * Remove all surfaces.
     */
    removeAllSurfaces(): void;
    /**
     * Clear scene of all objects
     */
    clear(): void;
    /**
     * Add specified properties to all atoms matching input argument
     * @example
     * $.get('../test_structs/b.sdf', function(data){
     *                       viewer.addModel(data,'sdf');
     *                       let props = [];
     *                       //make the atom index a property x
     *                       for(let i = 0; i < 8; i++) {
     *                         props.push({index:i,props:{'x':i}});
     *                       }
     *                       viewer.mapAtomProperties(props);
     *                       viewer.setStyle({sphere:{colorscheme:{gradient:'roygb',prop:'x',min:0,max:8}}});
     *                       viewer.zoomTo();
     *                       viewer.render();
     *                     });
     * @param props, - either array of atom selectors with associated props, or function that takes atom and sets its properties
     * @param sel - subset of atoms to work on - model selection must be specified here
     */
    mapAtomProperties(props: any, sel: AtomSelectionSpec): void;
    /**
     * Synchronize this view matrix of this viewer to the passed viewer.
     * When the viewpoint of this viewer changes, the other viewer will
     * be set to this viewer's view.
     */
    linkViewer(otherview: $3Dmol.GLViewer): void;
    /**
     * Return the z distance between the model and the camera
     * @returns distance
     */
    getPerceivedDistance(): number;
    /**
     * Set the distance between the model and the camera
     * Essentially zooming. Useful while stereo rendering.
     */
    setPerceivedDistance(): void;
    /**
     * Used for setting an approx value of eyeSeparation. Created for calling by StereoViewer object
     * @returns camera x position
     */
    setAutoEyeSeparation(): number;
  }


  /**
* A GLShape is a collection of user specified shapes.
* @param sid - Unique identifier
* @param stylespec - shape style specification
*/
  class GLShape extends ShapeSpec {
    constructor(sid: number, stylespec: ShapeSpec);
    /**
     * Update shape with new style specification
     * @example
     * let sphere = viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
     *             sphere.updateStyle({color:'yellow',opacity:0.5});
     *             viewer.render();
     */
    updateStyle(newspec: ShapeSpec): $3Dmol.GLShape;
    /**
     * Creates a custom shape from supplied vertex and face arrays
     */
    addCustom(customSpec: CustomShapeSpec): $3Dmol.GLShape;
    /**
     * Creates a sphere shape
     * @example
     * viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
     *
     *          viewer.render();
     */
    addSphere(sphereSpec: SphereSpec): $3Dmol.GLShape;
    /**
     * Creates a box
     * @example
     * var shape = viewer.addShape({color:'red'});
     *          shape.addBox({corner: {x:1,y:2,z:0}, dimensions: {w: 4, h: 2, d: 6}});
     *          shape.addBox({corner: {x:-5,y:-3,z:0},
     *                        dimensions: { w: {x:1,y:1,z:0},
     *                                      h: {x:-1,y:1,z:0},
     *                                      d: {x:0,y:0,z:1} }});
     *          viewer.zoomTo();
     *          viewer.rotate(30);
     *          viewer.render();
     */
    addBox(boxSpec: BoxSpec): $3Dmol.GLShape;
    /**
     * Creates a cylinder shape
     * @example
     * viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
     *                                   end:{x:10.0,y:0.0,z:0.0},
     *                                   radius:1.0,
     *                                   fromCap:1,
     *                                   toCap:2,
     *                                   color:'red',
     *                                   hoverable:true,
     *                                   clickable:true,
     *                                   callback:function(){ this.color.setHex(0x00FFFF00);viewer.render( );},
     *                                   hover_callback: function(){ viewer.render( );},
     *                                   unhover_callback: function(){ this.color.setHex(0xFF000000);viewer.render( );}
     *                                  });
     *               viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
     *                                   end:{x:0.0,y:10.0,z:0.0},
     *                                   radius:0.5,
     *                                   fromCap:false,
     *                                   toCap:true,
     *                                   color:'teal'});
     *               viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
     *                                   end:{x:20.0,y:0.0,z:0.0},
     *                                   radius:1.0,
     *                                   color:'black',
     *                                   fromCap:false,
     *                                   toCap:false});
     *               viewer.render();
     */
    addCylinder(cylinderSpec: CylinderSpec): $3Dmol.GLShape;
    /**
     * Creates a dashed cylinder shape
     */
    addDashedCylinder(cylinderSpec: CylinderSpec): $3Dmol.GLShape;
    /**
     * Creates a curved shape
     */
    addCurve(curveSpec: CurveSpec): $3Dmol.GLShape;
    /**
     * Creates a line shape
     * @example
     * $3Dmol.download("pdb:2ABJ",viewer,{},function(){
     *                   viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
     *                   viewer.render(callback);
     *               });
     */
    addLine(lineSpec: LineSpec): $3Dmol.GLShape;
    /**
     * Creates an arrow shape
     * @example
     * $3Dmol.download("pdb:4DM7",viewer,{},function(){
     *                   viewer.setBackgroundColor(0xffffffff);
     *                   viewer.addArrow({
     *                       start: {x:-10.0, y:0.0, z:0.0},
     *                       end: {x:0.0, y:-10.0, z:0.0},
     *                       radius: 1.0,
     *                       radiusRadio:1.0,
     *                       mid:1.0,
     *                       clickable:true,
     *                       callback:function(){
     *                           this.color.setHex(0xFF0000FF);
     *                           viewer.render( );
     *                       }
     *                   });
     *                   viewer.render();
     *                 });
     */
    addArrow(arrowSpec: ArrowSpec): $3Dmol.GLShape;
    /**
     * Create isosurface from voluemetric data.
     * @example
     * //the user can specify a selected region for the isosurface
     *          $.get('../test_structs/benzene-homo.cube', function(data){
     *                   var voldata = new $3Dmol.VolumeData(data, "cube");
     *                   viewer.addIsosurface(voldata, {isoval: 0.01,
     *                                                  color: "blue",
     *                                                  alpha: 0.5,
     *                                                  smoothness: 10});
     *                   viewer.addIsosurface(voldata, {isoval: -0.01,
     *                                                  color: "red",
     *                                                  smoothness: 5,
     *                                                  opacity:0.5,
     *                                                  wireframe:true,
     *                                                  clickable:true,
     *                                                  callback:
     *                                                  function() {
     *                                                      this.opacity = 0.0;
     *                                                      viewer.render( );
     *                                                  }});
     *                   viewer.setStyle({}, {stick:{}});
     *                   viewer.zoomTo();
     *                   viewer.render();
     *                 });
     * @param data - volumetric input data
     * @param isoSpec - volumetric data shape specification
     */
    addIsosurface(data: $3Dmol.VolumeData, isoSpec: IsoSurfaceSpec): void;
  }

  /**
 * $3Dmol.VolumeData stores volumetric data. This includes file parsing
 * functionality.
 * @param str - volumetric data
 * @param format - format of supplied data (cube, dx, vasp); append .gz if compressed
 * @param options - normalize (zero mean, unit variance), negate
 */
  class VolumeData {
    constructor(str: string, format: string, options: any);
    /**
     * @param x,y,z - the coordinates
     * @returns - index into flat array closest to provided coordinate; -1 if invalid
     */
    static getIndex(x, y, z: number): any;
    /**
     * @param x,y,z - the coordinates
     * @returns - value closest to provided coordinate; zero if coordinate invalid
     */
    static getVal(x, y, z: number): any;
  }

  /*TODO: find docstring*/
  class Mesh {
  }

  namespace Gradient {
    /**
     * Map value to hex color
     */
    function valueToHex(val: number, range: number): number;
    interface RWB extends $3Dmol.Gradient {
    }
    /**
     * Color scheme red to white to blue, for charges
     * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
     */
    class RWB implements $3Dmol.Gradient {
    }
    interface ROYGB extends $3Dmol.Gradient {
    }
    /**
     * rainbow gradient, but without purple to match jmol
     * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
     */
    class ROYGB implements $3Dmol.Gradient {
    }
    interface Sinebow extends $3Dmol.Gradient {
    }
    /**
     * rainbow gradient with constant saturation, all the way to purple!
     * Reverse gradients are supported when min>max so that the colors are displayed in reverse order.
     */
    class Sinebow implements $3Dmol.Gradient {
    }
  }

  /**
   * Create and initialize an appropriate viewer at supplied HTML element using specification in config
   * @example
   * var viewer = $3Dmol.createViewer(
   *      'gldiv', //id of div to create canvas in
   *      {
   *        defaultcolors: $3Dmol.elementColors.rasmol,
   *        backgroundColor: 'black'
   *      }
   *    );
   *
   * @param element - Either HTML element or string identifier
   * @param config - Viewer configuration
   * @param shared_viewer_resources - shared resources between viewers' renderers
   * @returns GLViewer, null if unable to instantiate WebGL
   */
  function createViewer(element: any | string, config: ViewerSpec, shared_viewer_resources: any): $3Dmol.GLViewer;

  declare type GridSpec = any;

  /**
   * Create and initialize an appropriate a grid of viewers that share a WebGL canvas
   * @example
   * var viewers = $3Dmol.createViewerGrid(
   *      'gldiv', //id of div to create canvas in
   *      {
   *        rows: 2,
   *        cols: 2,
   *        control_all: true  //mouse controls all viewers
   *      },
   *      { backgroundColor: 'lightgrey' }
   *    );
   *    $.get('data/1jpy.cif', function(data) {
   *      var viewer = viewers[0][0];
   *      viewer.addModel(data,'cif');
   *      viewer.setStyle({sphere:{}});
   *      viewer.zoomTo();
   *      viewer.render( );
   *
   *      viewer = viewers[0][1];
   *      viewer.addModel(data,'cif');
   *      viewer.setStyle({stick:{}});
   *      viewer.zoomTo();
   *      viewer.render( );
   *
   *      viewer = viewers[1][0];
   *      viewer.addModel(data,'cif');
   *      viewer.setStyle({cartoon:{color:'spectrum'}});
   *      viewer.zoomTo();
   *      viewer.render( );
   *
   *      viewer = viewers[1][1];
   *      viewer.addModel(data,'cif');
   *      viewer.setStyle({cartoon:{colorscheme:'chain'}});
   *      viewer.zoomTo();
   *      viewer.render();
   *
   *
   *    });
   * @param element - Either HTML element or string identifier
   * @param config - grid configuration
   * @param viewer_config - Viewer specification to apply to all subviewers
   * @returns [[$3Dmol.GLViewer]] 2D array of GLViewers
   */
  function createViewerGrid(element: any | string, config: GridSpec, viewer_config: ViewerGridSpec): any;
  /**
   * Contains a dictionary of embedded viewers created from HTML elements
   * with a the viewer_3Dmoljs css class indexed by their id (or numerically
   * if they do not have an id).
   */
  var viewers: any;
  /**
   * Download binary data (e.g. a gzipped file) into an array buffer and provide
   * arraybuffer to callback.
   * @param uri - location of data
   * @param callback - Function to call with arraybuffer as argument.
   * @param request - type of request
   */
  function getbin(uri: string, callback: (...params: any[]) => any, request: string): Promise<any>;
  /**
   * Convert a base64 encoded string to a Uint8Array
   * @param base64 - encoded string
   */
  function base64ToArray(base64: string): void;
  /**
   * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
   * @example
   * viewer.setBackgroundColor(0xffffffff);
   *        $3Dmol.download('pdb:2nbd',viewer,{onemol: true,multimodel: true},function(m) {
   *         m.setStyle({'cartoon':{colorscheme:{prop:'ss',map:$3Dmol.ssColors.Jmol}}});
   *        viewer.zoomTo();
   *        viewer.render(callback);
   *     });
   * @param query - String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
   * @param viewer - Add new model to existing viewer
   * @param options - Specify additional options
   *                           format: file format to download, if multiple are available, default format is pdb
   *                           pdbUri: URI to retrieve PDB files, default URI is http://www.rcsb.org/pdb/files/
   * @param callback - Function to call with model as argument after data is loaded.
   * @returns GLModel, Promise if callback is not provided
   */
  function download(query: string, viewer: $3Dmol.GLViewer, options: any, callback: (...params: any[]) => any): $3Dmol.GLModel;
  /**
   * $3Dmol surface types
   */
  enum SurfaceType {
    VDW = 1,
    MS = 2,
    SAS = 3,
    SES = 4
  }
  function mergeGeos(geometry: $3Dmol.Geometry, mesh: $3Dmol.Mesh): undefined;
  /**
   * Render surface synchronously if true
   */
  var syncSurface: boolean;
  /**
   * Parse a string that represents a style or atom selection and convert it
   * into an object.  The goal is to make it easier to write out these specifications
   * without resorting to json. Objects cannot be defined recursively.
   * ; - delineates fields of the object
   * : - if the field has a value other than an empty object, it comes after a colon
   * , - delineates key/value pairs of a value object
   *     If the value object consists of ONLY keys (no = present) the keys are
   *     converted to a list.  Otherwise a object of key/value pairs is created with
   *     any missing values set to null
   * = OR ~ - separates key/value pairs of a value object, if not provided value is null
   *     twiddle is supported since = has special meaning in URLs
   * @param (String) - str
   */
  function specStringToObject(String: string): any;
  /**
   * computes the bounding box around the provided atoms
   */
  function getExtent(atomlist: AtomSpec[]): any[];
  class Camera {
  }
  class EventDispatcher { }

  class Raycaster {
  }
  class Projector {
  }
  /**
   * Sprite render plugin
   */
  class SpritePlugin {
  }
  /**
   * Line and Mesh material types
   */
  class Material {
  }
  class LineBasicMaterial {
  }
  class MeshLambertMaterial {
  }
  class MeshDoubleLambertMaterial {
  }
  class MeshOutlineMaterial {
  }
  class ImposterMaterial {
  }
  class VolumetricMaterial {
  }
  class SpriteMaterial {
  }
  class Texture {
  }
  class Quaternion {
  }
  class Vector2 {
  }
  class Ray {
  }
  class Line {
  }
  class Sprite {
  }
  /**
   * Simplified webGL renderer
   */
  function Renderer(): void;
  class Scene {
  }
  class Fog {
  }
  class Sphere {
  }
  class Cylinder {
  }
  class Triangle {
  }
  /**
   * Preset secondary structure color scheme
   */
  var ssColors: any;
  /**
   * Preset element coloring - from individual element colors to entire mappings (e.g. '$3Dmol.elementColors.Jmol' colors atoms with Jmol stylings)
   */
  var elementColors: any;
  /**
   * @property built - in color schemes
   * The user can pass all of these values directly as the colorscheme and they will use the respective colorscheme
   */
  var builtinColorSchemes: {
    built: any;
  };


  declare type Color = number;
  declare type AtomStyle = any;

  /**
   * Return proper color for atom given style
   */
  function getColorFromStyle(atom: AtomSpec, style: AtomStyle): $3Dmol.Color;
  /**
   * Enum for cylinder cap styles.
   */
  enum CAP {
    NONE = 0,
    FLAT = 1,
    ROUND = 2
  }

  declare type Point = any;
  declare type geometry = any;

  /**
   * Lower level utilities for creating WebGL shape geometries.
   * These are not intended for general consumption.
   */
  namespace GLDraw {
    /**
     * Create a cylinder
     * @param fromCap - 0 for none, 1 for flat, 2 for round
     * @param toCap - = 0 for none, 1 for flat, 2 for round
     */
    function drawCylinder(geo: geometry, from: $3Dmol.Point, to: Point, radius: number, color: $3Dmol.Color, fromCap: $3Dmol.CAP, toCap: $3Dmol.CAP): void;
    /**
     * Create a cone
     */
    function drawCone(geo: geometry, from: Point, to: Point, radius: number, color: $3Dmol.Color): void;
    /**
     * Create a sphere.
     */
    function drawSphere(geo: geometry, pos: Point, radius: number, color: $3Dmol.Color): void;
  }



  /**
   * Color mapping gradients
   */
  interface Gradient {
  }

  /**
   * $3Dmol.StateManager - StateManager creates the space to preserve the state of the ui and sync it with the GLViewer
   * @param glviewer - StateManager is required to have interaction between glviewer and the ui.
   * @param config - Loads the user defined parameters to generate the ui and handle state
   */
  class StateManager {
    constructor(glviewer: $3Dmol.GLViewer, config: any);
    /**
     * Add Selection from the ui to glviewer
     * @param spec - Object that contains the output from the form
     * @param sid - If surface id being edited then sid is set to some string
     * @returns String
     */
    addSelection(spec: any, sid: string): any;
    /**
     * Return true if the selections contain at least one atom
     * @param sel - Atom selection spec
     * @returns Boolean
     */
    checkAtoms(sel: AtomSelectionSpec): any;
    /**
     * Toggle the hidden property of the selection
     * @param sid - Selection id
     */
    toggleHide(sid: string): void;
    /**
     * Add style and renders it into the viewport
     * @param spec - Output object of style form
     * @param sid - Selection Id
     * @param stid - Style Id
     * @returns String
     */
    addStyle(spec: string, sid: string, stid: string): any;
    /**
     * Removes the style specified by stid
     * @param sid - Selection id
     * @param stid - Style Id
     */
    removeStyle(sid: string, stid: string): void;
    /**
     * Toggle hidden property of a style
     * @param sid - Selection Id
     * @param stid - Style Id
     */
    toggleHideStyle(sid: string, stid: string): void;
    /**
     * Adds surface to the viewport
     * @param property - Surface output object
     * @param callback - callback
     * @returns String
     */
    addSurface(property: any, callback: (...params: any[]) => any): any;
    /**
     * Removes surface from the viewport
     * @param id - Surface Id
     */
    removeSurface(id: string): void;
    /**
     * Edit the exisiting surface in the viewport
     * @param surfaceProperty - Surface Style
     */
    editSurface(surfaceProperty: any): void;
    /**
     * Returns the list of ids of selections that are created so far
     * @returns <Array of selection ids>
     */
    getSelectionList(): any;
    /**
     * Opens context menu when called from glviewer
     * @param atom - Atom spec obtained from context menu event
     * @param x - x coordinate of mouse on viewport
     * @param y - y coordinate of mouse on the viewport
     */
    openContextMenu(atom: AtomSpec, x: number, y: number): void;
    /**
     * Adds Label to the viewport specific to the selection
     * @param labelValue - Output object from label form of Context Menu
     */
    addLabel(labelValue: any): void;
    /**
     * Adds atom label to the viewport
     * @param labelValue - Output object from propertyMenu form of Context Menu
     * @param atom - Atom spec that are to be added in the label
     */
    addAtomLabel(labelValue: any, atom: AtomSpec): void;
    /**
     * Executes hide context menu and process the label if needed
     * @param processContextMenu - Specify the need to process the values in the context menu
     */
    exitContextMenu(processContextMenu: boolean): void;
    /**
     * Removes the atom label from the viewpoer
     * @param atom - Atom spec
     */
    removeAtomLabel(atom: AtomSpec): void;
    /**
     * Add model to the viewport
     * @param modelDesc - Model Toolbar output
     */
    addModel(modelDesc: any): void;
    /**
     * Updates the state variable for selections and styles and trigger ui to show the
     * ui elements for these selections and styles.
     * @param selSpec - Atom Selection Spec
     * @param styleSpec - Atom Style Spec
     */
    createSelectionAndStyle(selSpec: AtomSelectionSpec, styleSpec: AtomStyleSpec): void;
    /**
     * Creates selection and add surface with reference to that selection
     * and triggers updates in the ui
     * @param surfaceType - Type of surface to be created
     * @param sel - Atom selection spec
     * @param style - Atom style spec
     * @param sid - selection id
     */
    createSurface(surfaceType: string, sel: AtomSelectionSpec, style: AtomStyleSpec, sid: string): void;
    /**
     * Sets the value of title in ModelToolBar
     * @param title - Model title
     */
    setModelTitle(title: string): void;
    /**
     * Updates the UI on viewport change
     */
    updateUI(): void;
  }


  /** @todo add type */
  declare type VolumetricRenderSpec = any;

  /**
   * A GLVolumetricRender is a "shape" for representing volumetric data as a density distribution.
   * @param data - volumetric data
   * @param spec - specification of volumetric render
   */
  class GLVolumetricRender {
    constructor(data: $3Dmol.VolumeData, spec: VolumetricRenderSpec);
  }

  /**
   * Color representation.
   * @property 0xAF10AB - any hex number
   * @property <html - color name>
   */
  declare type ColorSpec = {
    html: string;
  } | number;

  /**
   * @example
   * //Using a function in order to define the colors.
   *   $3Dmol.download("pdb:4UAA",viewer,{},function(){
   *                   viewer.setBackgroundColor(0xffffffff);
   *                   var colorAsSnake = function(atom) {
   *                     return atom.resi % 2 ? 'white': 'green'
   *                   };
   *
   *                   viewer.setStyle( {chain:'A'}, { cartoon: {colorfunc: colorAsSnake }});
   *                   viewer.setStyle( {chain:'B'}, { stick: {colorscheme: 'yellowCarbon'}});
   *
   *                   viewer.render();
   *               });
   * @property <html - color>Carbon   - use default element colors but with carbon set to specify html color string
   * @property ssPyMOL - PyMol secondary colorscheme
   * @property ssJmol - Jmol secondary colorscheme
   * @property Jmol - Jmol primary colorscheme
   * @property default - default colorscheme
   * @property amino - amino acid colorscheme
   * @property shapely - shapely protien colorscheme
   * @property nucleic - nucleic acid colorscheme
   * @property chain - standard chain colorscheme
   * @property chainHetatm - chain Hetatm colorscheme
   * @property prop - atomSpec property. Example 'b'. See AtomSpec.
   * @property gradient - Allows the user to provide a gradient to the colorscheme.  Is either a $3Dmol.Gradient object or the name of a built-in gradient (rwb, roygb, sinebow)
   * @property undefined - min value for gradient
   * @property undefined - max value for gradient
   * @property undefined - mid point value for gradient (for rwb)
   * @property map - map of a certain AtomSpec property to a color of the form `{'prop': 'elem', map:$3Dmol.elementColors.greenCarbon}` Allows the user to provide a mapping of elements to colors to the colorscheme.  This can be done with any properties, and not just 'elem'.
   * @property colorfunc - Allows the user to provide a function for setting the colorschemes.
   */
  declare type ColorschemeSpec = {
    html: string;
    ssPyMOL: string;
    ssJmol: string;
    Jmol: string;
    default: string;
    amino: string;
    shapely: string;
    nucleic: string;
    chain: string;
    chainHetatm: string;
    prop: string;
    gradient: Gradient;
    map: any;
    colorfunc: (...params: any[]) => any;
  };

  /**
   * A visualization of protein or nucleic acid secondary structure.  Applying this to other molecules will not show anything.
   * @example
   * $3Dmol.download("pdb:4ZD3",viewer,{},function(){
   *                   viewer.setBackgroundColor(0xffffffff);
   *                   viewer.setViewStyle({style:"outline"});
   *                   viewer.setStyle({},{cartoon:{}});
   *                   viewer.render();
   *               });
   * @property color - strand color, may specify as 'spectrum' which will apply reversed gradient based on residue number
   * @property style - style of cartoon rendering (trace, oval, rectangle
   *       (default), parabola, edged)
   * @property ribbon - whether to use constant strand width, disregarding
   *       secondary structure; use thickness to adjust radius
   * @property arrows - whether to add arrows showing beta-sheet
   *       directionality; does not apply to trace or ribbon
   * @property tubes - whether to display alpha helices as simple cylinders;
   *       does not apply to trace
   * @property thickness - cartoon strand thickness, default is 0.4
   * @property width - cartoon strand width, default is secondary
   *       structure-dependent; does not apply to trace or ribbon
   * @property opacity - set opacity from 0-1; transparency is set per-chain
   *       with a warning outputted in the event of ambiguity
   * @property In - nucleic acids, the base cylinders obtain their color from the
   *       atom to which the cylinder is drawn, which is 'N1' for purines (resn:
   *       'A', 'G', 'DA', 'DG') and 'N3' for pyrimidines (resn: 'C', 'U', 'DC',
   *       'DT'). The different nucleobases can therefore be distinguished as
   *       follows:
   */
  declare type CartoonStyleSpec = {
    color: ColorSpec;
    style: string;
    ribbon: boolean;
    arrows: boolean;
    tubes: boolean;
    thickness: number;
    width: number;
    opacity: number;
    In: any;
  };

  /**
   * @property hidden - do not show
   * @property linewidth - *deprecated due to vanishing browser support*
   * @property scale - scale radius by specified amount
   * @property colorscheme - element based coloring
   * @property color - fixed coloring, overrides colorscheme
   * @property opacity - opacity, must be the same for all atoms in the model
   */
  declare type CrossStyleSpec = {
    hidden: boolean;
    linewidth: number;
    radius: number;
    scale: number;
    colorscheme: ColorschemeSpec;
    color: ColorSpec;
    opacity: number;
  };

  /**
   * @property hidden - do not show line
   * @property linewidth - *deprecated due to vanishing browser support*
   * @property colorscheme - element based coloring
   * @property color - fixed coloring, overrides colorscheme
   * @property opacity - opacity, must be the same for all atoms in the model
   */
  declare type LineStyleSpec = {
    hidden: boolean;
    linewidth: number;
    colorscheme: ColorschemeSpec;
    color: ColorSpec;
    opacity: number;
  };

  /**
   * @property hidden - do not show atom
   * @property radius - override van der waals radius
   * @property scale - scale radius by specified amount
   * @property colorscheme - element based coloring
   * @property color - fixed coloring, overrides colorscheme
   * @property opacity - opacity, must be the same for all atoms in the model
   */
  declare type SphereStyleSpec = {
    hidden: boolean;
    radius: number;
    scale: number;
    colorscheme: ColorschemeSpec;
    color: ColorSpec;
    opacity: number;
  };

  /**
   * @property hidden - do not show
   * @property singleBonds - draw all bonds as single bonds if set
   * @property colorscheme - element based coloring
   * @property color - fixed coloring, overrides colorscheme
   * @property opacity - opacity, must be the same for all atoms in the model
   */
  declare type StickStyleSpec = {
    hidden: boolean;
    radius: number;
    singleBonds: boolean;
    colorscheme: ColorschemeSpec;
    color: ColorSpec;
    opacity: number;
  };

  /**
   * Unit Cell shape specification.
   * @property box - line style used to draw box
   * @property astyle - arrow specification of "a" axis
   * @property bstyle - arrow specification of "b" axis
   * @property cstyle - arrow specification of "c" axis
   * @property alabel - label for a axis
   * @property alabelstyle - label style for a axis
   * @property blabel - label for b axis
   * @property blabelstyle - label style for b axis
   * @property clabel - label for c axis
   * @property clabelstyle - label style for c axis
   */
  declare type UnitCellStyleSpec = {
    box: LineStyleSpec;
    astyle: ArrowSpec;
    bstyle: ArrowSpec;
    cstyle: ArrowSpec;
    alabel: string;
    alabelstyle: LabelSpec;
    blabel: string;
    blabelstyle: LabelSpec;
    clabel: string;
    clabelstyle: LabelSpec;
  };

  /**
   * Label type specification
   * @property font - font name, default sans-serif
   * @property fontSize - height of text, default 18
   * @property fontColor - font color, default white
   * @property fontOpacity - font opacity, default 1
   * @property borderThickness - line width of border around label, default 0
   * @property borderColor - color of border, default backgroundColor
   * @property borderOpacity - color of border
   * @property backgroundColor - color of background, default black
   * @property backgroundOpacity - opacity of background, default 1
   * @property position - x,y,z coordinates for label
   * @property screenOffset - x,y _pixel_ offset of label from position
   * @property inFront - always put labels in from of model
   * @property showBackground - show background rounded rectangle, default true
   * @property fixed - sets the label to change with the model when zooming
   * @property useScreen - position is in screen (not model) coordinates which are pixel offsets from upper left corner.
   * @property backgroundImage - An element to draw into the label.  Any CanvasImageSource is allowed.
   * @property alignment - how to orient the label w/respect to position: topLeft (default), topCenter, topRight, centerLeft, center, centerRight, bottomLeft, bottomCenter, bottomRight
   * @property frame - if set, only display in this frame of an animation
   */
  declare type LabelSpec = {
    font: string;
    fontSize: number;
    fontColor: ColorSpec;
    fontOpacity: number;
    borderThickness: number;
    borderColor: ColorSpec;
    borderOpacity: string;
    backgroundColor: ColorSpec;
    backgroundOpacity: string;
    position: $3Dmol.Vector3;
    screenOffset: $3Dmol.Vector2;
    inFront: boolean;
    showBackground: boolean;
    fixed: boolean;
    useScreen: boolean;
    backgroundImage: any;
    alignment: string;
    frame: number;
  };

  /**
   * Parser options specification. Used to specify the options of a GLModel.  Depending on the input file format, not all fields may be defined.
   * @property frames - true if you want to add to a new frame and false otherwise ; supported by all
   * @property vibrate - object specifying the vibration behavior ; supported by all
   * @property vibrate.frames - number of frames to be created, default to 10 ; supported by all
   * @property vibrate.amplitude - amplitude of distortion, default to 1 (full) ; supported by all
   * @property multimodel - specifies weather or not multiple models are being defined ; supported by xyz,sdf, or mol2
   * @property onemol - specifies weather or not the model is of one molecule ; Supported by xyz , sdf , mol2
   * @property keepH - do not strip hydrogens ; supported by sdf,mol2
   * @property parseStyle - used to define ChemDoodle styles ; supported by cdjson
   * @property doAssembly - boolean dictating weather or not to do assembly ; supported by mcif
   * @property duplicateAssemblyAtoms- - Set to true if you wish to duplicate assembly atoms otherwise false ; supported by all formats with symmetries.  Not duplicating will result in faster rendering but it will not be possible to individually style symmetries.
   * @property normalizeAssembly - shift symmetry mates so their centroid is in the unit cell
   * @property dontConnectDuplicatedAtoms - do not detect bonds between symmetries generated with duplicateAssemblyAtoms (cif only - other formats never make bonds between symmetries)
   * @property noSecondaryStructure - boolean dictating the presence of a secondary structure ; supported by pdb
   * @property noComputeSecondaryStructure - do not compute ss ; supported by pdb
   * @property altLoc - which alternate location to select, if present; '*' to load all ; supported by pdb
   * @property assemblyIndex - index of the assembly in symmetry ; supported by mmtf
   * @property assignBonds - for formats without explicit bonds (e.g. PDB, xyz) infer bonding (default true).
   */
  declare type ParserOptionsSpec = {
    frames: boolean;
    vibrate: {
      frames: number;
      amplitude: number;
    };
    multimodel: boolean;
    onemol: boolean;
    keepH: boolean;
    parseStyle: any;
    doAssembly: boolean;
    duplicateAssemblyAtoms: boolean;
    normalizeAssembly: boolean;
    dontConnectDuplicatedAtoms: boolean;
    noSecondaryStructure: boolean;
    noComputeSecondaryStructure: boolean;
    altLoc: string;
    assemblyIndex: number;
    assignBonds: boolean;
  };

  /**
   * GLViewer input specification
   * @property callback - Callback function to be immediately executed on this viewer
   * @property defaultcolors - Object defining default atom colors as atom => color property value pairs for all models within this viewer
   * @property nomouse - Whether to disable disable handling of mouse events.
   *                 If you want to use your own mouse handlers, set this then bind your handlers to the canvas object.
   *                 The default 3Dmol.js handlers are available for use:
   *                 'mousedown touchstart': viewer._handleMouseDown,
   *                 'DOMMouseScroll mousewheel': viewer._handleMouseScroll
   *                 'mousemove touchmove': viewer._handleMouseMove
   * @property backgroundColor - Color of the canvas background
   * @property backgroundAlpha - Alpha transparency of canvas background
   * @property id - id of the canvas
   * @property cartoonQuality - default 5
   * @property disableFog - Disable fog, default to false
   */
  declare type ViewerSpec = {
    callback: (...params: any[]) => any;
    defaultcolors: any;
    nomouse: boolean;
    backgroundColor: string;
    backgroundAlpha: number;
    camerax: number;
    hoverDuration: number;
    id: string;
    cartoonQuality: number;
    row: number;
    col: number;
    rows: number;
    cols: number;
    canvas: any;
    viewers: any;
    minimumZoomToDistance: any;
    lowerZoomLimit: any;
    upperZoomLimit: any;
    antialias: boolean;
    control_all: boolean;
    orthographic: boolean;
    disableFog: boolean;
  };

  /**
   * Grid GLViewer input specification
   * @property rows - number of rows in grid
   * @property cols - number of columns in grid
   * @property control_all - if true, mouse events are linked
   */
  declare type ViewerGridSpec = {
    rows: number;
    cols: number;
    control_all: boolean;
  };

  /**
   * Atom representation. Depending on the input file format, not all fields may be defined.
   * @property resn - Parent residue name
   * @property x - Atom's x coordinate
   * @property y - Atom's y coordinate
   * @property z - Atom's z coordinate
   * @property color - Atom's color, as hex code or built-in color string
   * @property surfaceColor - Hex code for color to be used for surface patch over this atom
   * @property elem - Element abbreviation (e.g. 'H', 'Ca', etc)
   * @property hetflag - Set to true if atom is a heteroatom
   * @property chain - Chain this atom belongs to, if specified in input file (e.g 'A' for chain A)
   * @property resi - Residue number
   * @property serial - Atom's serial id number
   * @property atom - Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon)
   * @property bonds - Array of atom ids this atom is bonded to
   * @property ss - Secondary structure identifier (for cartoon render; e.g. 'h' for helix)
   * @property singleBonds - true if this atom forms only single bonds or no bonds at all
   * @property bondOrder - Array of this atom's bond orders, corresponding to bonds identfied by 'bonds'
   * @property properties - Optional mapping of additional properties
   * @property b - Atom b factor data
   * @property pdbline - If applicable, this atom's record entry from the input PDB file (used to output new PDB from models)
   * @property clickable - Set this flag to true to enable click selection handling for this atom
   * @property callback - Callback click handler function to be executed on this atom and its parent viewer
   * @property invert - for selection, inverts the meaning of the selection
   */
  declare type AtomSpec = {
    resn: string;
    x: number;
    y: number;
    z: number;
    color: ColorSpec;
    surfaceColor: ColorSpec;
    elem: string;
    hetflag: boolean;
    chain: string;
    resi: number;
    icode: number;
    rescode: number;
    serial: number;
    atom: string;
    bonds: number[];
    ss: string;
    singleBonds: boolean;
    bondOrder: number[];
    properties: any;
    b: number;
    pdbline: string;
    clickable: boolean;
    callback: (...params: any[]) => any;
    invert: boolean;
  };

  /**
   * 3 dimensional vector
   * @property x - x coordinate
   * @property y - y coordinate
   * @property z - z coordinate
   */
  declare type Vector3 = {
    x: number;
    y: number;
    z: number;
  };

  /**
   * Atom selection object. Used to specify what atoms should be selected.  Can include
   * any field from {@link AtomSpec} in which case atoms must equal the specified value.
   * All fields must match for the selection to hold. If values
   * are provided as a list, then only one value of the list must match.
   * @example
   * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
   *                   viewer.setStyle({chain:'B'},{cartoon:{color:'spectrum'}});
   *                   viewer.setStyle({chain:'B',invert:true},{cartoon:{}});
   *                   viewer.setStyle({bonds: 0},{sphere:{radius:0.5}}); //water molecules
   *                   viewer.setStyle({resn:'PMP',byres:true,expand:5},{stick:{colorscheme:"greenCarbon"}});
   *                   viewer.setStyle({resi:["91-95","42-50"]},{cartoon:{color:"green",thickness:1.0}});
   *                   viewer.render();
   *
   *
   *                 });
   * @property ... - any field from {@link AtomSpec}, values may be singletons or lists. Integer numerical ranges are supported as strings.
   * @property model - a single model or list of models from which atoms should be selected.  Can also specify by numerical creation order.  Reverse indexing is allowed (-1 specifies last added model).
   * @property bonds - overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms
   * @property predicate - user supplied function that gets passed an {AtomSpec} and should return true if the atom should be selected
   * @property invert - if set, inverts the meaning of the selection
   * @property byres - if set, expands the selection to include all atoms of any residue that has any atom selected
   * @property expand - expands the selection to include all atoms within a given distance from the selection
   * @property within - intersects the selection with the set of atoms within a given distance from another selection
   * @property and - take the intersection of the provided lists of {AtomSelectionSpec}s
   * @property or - take the union of the provided lists of {AtomSelectionSpec}s
   * @property not - take the inverse of the provided {AtomSelectionSpec}
   */
  declare type AtomSelectionSpec = {
    model: $3Dmol.GLModel;
    bonds: number;
    predicate: (...params: any[]) => any;
    invert: boolean;
    byres: boolean;
    expand: number;
    within: WithinSelectionSpec;
    and: AtomSelectionSpec[];
    or: AtomSelectionSpec[];
    not: AtomSelectionSpec;
  };

  /**
   * Within selection object. Used to find the subset of an atom selection that is within
   * some distance from another atom selection. When added as a field of an {@link AtomSelectionSpec},
   * intersects the set of atoms in that selection with the set of atoms within a given
   * distance from the given {@link AtomSelectionSpec}.
   * @example
   * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
   *
   *                   viewer.setStyle({chain: 'A', within:{distance: 10, sel:{chain: 'B'}}}, {sphere:{}});
   *                   viewer.render();
   *                 });// stylizes atoms in chain A that are within 10 angstroms of an atom in chain B
   * @property distance - the distance in angstroms away from the atom selection to include atoms in the parent selection
   * @property invert - if set, selects atoms not within distance range for intersection
   * @property sel - the selection of atoms against which to measure the distance from the parent atom selection
   */
  declare type WithinSelectionSpec = {
    distance: number;
    invert: boolean;
    sel: AtomSelectionSpec;
  };

  /**
   * @property line - draw bonds as lines
   * @property cross - draw atoms as crossed lines (aka stars)
   * @property stick - draw bonds as capped cylinders
   * @property sphere - draw atoms as spheres
   * @property cartoon - draw cartoon representation of secondary structure
   * @property clicksphere - invisible style for click handling only
   */
  declare type AtomStyleSpec = {
    line: LineStyleSpec;
    cross: CrossStyleSpec;
    stick: StickStyleSpec;
    sphere: SphereStyleSpec;
    cartoon: CartoonStyleSpec;
    clicksphere: ClickSphereStyleSpec;
  };

  /** TODO: find docs */
  declare type ClickSphereStyleSpec = any;

  /**
   * @example
   * var setStyles = function(volumedata){
   *                     var data = new $3Dmol.VolumeData(volumedata, "cube");
   *                     viewer.addSurface("VDW", {opacity:0.85, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)},{chain:'A'});
   *                     viewer.mapAtomProperties($3Dmol.applyPartialCharges);
   *                     viewer.addSurface($3Dmol.SurfaceType.SAS, {map:{prop:'partialCharge',scheme:new $3Dmol.Gradient.RWB(-.05,.05)}, opacity:1.0},{chain:'B'});
   *                     viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: data, color:'red'},{chain:'C'});
   *                     viewer.addSurface($3Dmol.SurfaceType.SAS, {opacity:0.85,voldata: data, colorscheme:'greenCarbon'},{chain:'D'});
   *
   *               viewer.render();
   *               };
   *               $3Dmol.download("pdb:4DLN",viewer,{},function(){
   *
   *                   $.get("data/1fas.cube",setStyles);
   *                 });
   * @property opacity - sets the transparency: 0 to hide, 1 for fully opaque
   * @property colorscheme - element based coloring
   * @property color - fixed coloring, overrides colorscheme
   * @property voldata - volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified
   * @property volscheme - coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields.
   * @property volformat - format of voldata if not a $3Dmol.VolumeData object
   * @property map - specifies a numeric atom property (prop) and color mapping (scheme) such as {@link $3Dmol.Gradient.RWB}.  Deprecated, use colorscheme instead.
   */
  declare type SurfaceStyleSpec = {
    opacity: number;
    colorscheme: ColorschemeSpec;
    color: ColorSpec;
    voldata: $3Dmol.VolumeData;
    volscheme: $3Dmol.Gradient;
    volformat: string;
    map: any;
  };

  /**
   * Isosurface style specification
   * @property isoval - specifies the isovalue to draw surface at
   * @property color - solid color
   * @property opacity - transparency, between 0 and 1
   * @property wireframe - draw as wireframe, not surface
   * @property linewidth - width of line for wireframe rendering **No longer supported by most browsers**
   * @property smoothness - amount to smooth surface (default 1)
   * @property coords - coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates
   * @property seldist - distance around coords to include data [default = 2.0]
   * @property voldata - volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified
   * @property volscheme - coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields.
   * @property volformat - format of voldata if not a $3Dmol.VolumeData object
   * @property clickable - if true, user can click on object to trigger callback
   * @property callback - function to call on click
   */
  declare type IsoSurfaceSpec = {
    isoval: number;
    color: ColorSpec;
    opacity: number;
    wireframe: boolean;
    linewidth: number;
    smoothness: number;
    coords: Array<any>;
    seldist: number;
    voldata: $3Dmol.VolumeData;
    volscheme: $3Dmol.Gradient;
    volformat: string;
    clickable: boolean;
    callback: (...params: any[]) => any;
  };

  /**
   * VolumetricRenderer style specification
   * @property transferfn - list of objects containing @color, @opacity and @value properties to specify color per voxel data value
   * @property subsamples - number of times to sample each voxel approximately (default 5)
   */
  declare type VolumetricRendererSpec = {
    transferfn: Array<any>
    subsamples: number;
  };

  /**
   * GLShape style specification
   * @property color - solid color
   * @property alpha - transparency
   * @property wireframe - draw as wireframe, not surface
   * @property hidden - if true, do not display object
   * @property linewidth - width of line for wireframe rendering **No longer supported by most browsers**
   * @property clickable - if true, user can click on object to trigger callback
   * @property callback - function to call on click
   * @property frame - if set, only display in this frame of an animation
   */
  declare abstract class ShapeSpec {
    color: ColorSpec;
    alpha: number;
    wireframe: boolean;
    hidden: boolean;
    linewidth: number;
    clickable: boolean;
    callback: (...params: any[]) => any;
    frame: number;
  };

  /**
   * Specification for adding custom shape. Extends {@link ShapeSpec}.
   * @property vertexArr - List of vertex positions
   * @property normalArr - List of normal vectors for each vertex
   * @property faceArr - List of triangles to build the custom shape. Each triangle is defined by the indices of 3 vertices in vertexArr, so the array length should be 3 times the number of faces.
   * @property color - Either a single color for the whole object or an array specifying the color at each vertex.
   */
  declare type CustomShapeSpec = {
    vertexArr: $3Dmol.Vector3[];
    normalArr: $3Dmol.Vector3[];
    faceArr: number[];
    color: ColorSpec | ColorSpec[];
  };

  /**
   * Sphere shape specification. Extends {@link ShapeSpec}
   */
  declare type SphereShapeSpec = {
    center: $3Dmol.Vector3;
    radius: number;
  };

  /**
   * Box shape specification. Extends {@link ShapeSpec}
   * @property corner - bottom corner of box
   * @property center - alternative to corner: center of box
   * @property dimensions - {w:width, h:height, d:depth}; can be either scalars or vectors (for not-axis aligned boxes)
   */
  declare type BoxSpec = {
    corner: $3Dmol.Vector3;
    center: $3Dmol.Vector3;
    dimensions: any;
  };

  /**
   * Arrow shape specification.  Extends {@link ShapeSpec}
   * @property radiusRatio - ratio of arrow base to cylinder (1.618034 default)
   * @property mid - relative position of arrow base (0.618034 default)
   * @property midpos - position of arrow base in length units, if negative positioned from end instead of start.  Overrides mid.
   */
  declare type ArrowSpec = {
    start: $3Dmol.Vector3;
    end: $3Dmol.Vector3;
    radius: number;
    color: ColorSpec;
    hidden: boolean;
    radiusRatio: number;
    mid: number;
    midpos: number;
  };

  /**
   * Cylinder shape specification.  Extends {@link ShapeSpec}
   * @property fromCap - 0 for none, 1 for flat, 2 for round
   * @property toCap - 0 for none, 1 for flat, 2 for round
   */
  declare type CylinderSpec = {
    start: $3Dmol.Vector3;
    end: $3Dmol.Vector3;
    radius: number;
    fromCap: $3Dmol.CAP;
    toCap: $3Dmol.CAP;
    dashed: boolean;
  };

  /**
   * Curve shape specification.  Extends {@link ShapeSpec}
   * @property points - list of (x,y,z) points to interpolate between to make curve
   * @property smooth - amount of interpolation
   * @property fromArrow - if an arrow should be drawn at the start
   * @property toArrow - if an arrow should be drawn at the end
   */
  declare type CurveSpec = {
    points: $3Dmol.Vector3;
    smooth: number;
    radius: number;
    fromArrow: boolean;
    toArrow: boolean;
  };

  /**
   * Line shape specification.  Extends {@link ShapeSpec}  (but defaults to wireframe)
   */
  declare type LineSpec = {
    start: $3Dmol.Vector3;
    end: $3Dmol.Vector3;
  };

  /**
   * File formats supported by 3Dmol.js
   * @property cdjson,json - Chemical JSON format
   * @property cube - Gaussian cube format
   * @property gro - Gromacs topology format, need to add coordinates to resulting model.
   * @property mcif,cif - Crystallographic Information File, the successor to PDB that makes you miss the PDB file format
   * @property mmtf - Macromolecular Transmission Format, the successor to PDB that is totally awesome
   * @property mol2 - Sybyl Mol2 format
   * @property pdb - The venerable Protein Data Bank format
   * @property pqr - Like PDB but with partial charges which are read into the partialcharge atom property
   * @property prmtop - Amber topology file, must add coordinates
   * @property sdf - MDL MOL format, supports multiple models and meta data
   * @property vasp - VASP format (CONTCAR, POSCAR)
   * @property xyz - XYZ cartesian coordinates format
   */
  declare type FileFormats =
    "cdjson " |
    "json" |
    "cube" |
    "gro" |
    "mcif, cif" |
    "mmtf" |
    "mol2" |
    "pdb" |
    "pqr" |
    "prmtop" |
    "sdf" |
    "vasp" |
    "xyz";


  declare namespace UI {
    declare type Icons = any;
  }

  /**
 * $3Dmol.UI - UI creates panels in the viewer to assist control of the viewport
 * @param stateManager - StateManager is required to have interaction between glviewer and the ui.
 * @param config - Loads the user defined parameters to generate the ui
 * @param parentElement - Refers the parent division used to hold the canvas for 3Dmol.js
 */
  class UI {
    constructor(stateManager: $3Dmol.StateManager, config: any, parentElement: any);
    /**
     * This is a colection of contructor to make different input element
     */
    Form(): void;
    /**
     * This is a colection of contructor to make different input element
     */
    Form(): void;
    /**
     * Generates the object to hold different icons present Icons : move, rotate, pencil, listArrow, option, minus, plus, painbrush, select, movie.play, move.pause, movie.stop, movie.next, move.previous, tick, cross, edit, remove, list, style, visible, invisible, mouse, nomouse, label, surface, molecule, change
     */
    Icons(): void;
    /**
     * Resize the panel with respect to the new viewport
     */
    resize(): void;


    /**
     * ModelToolbar is part of $3Dmol.UI to edit or change the model loaded into the viewer
     */
    static ModelToolbar(): void;

    /**
     * Selection box creates the UI panel to manipulate selections and style that are drawn
     * on the viewport
     * @param icon - takes the svg code for the icon that is to be used to display
     * selection box
     * @returns Jquery element of div
     */
    static SelectionBox(icon: $3Dmol.UI.Icons): any;

    /**
     * Card for manipulation of a selection form and related styles
     */
    static Selection(): void;

    /**
     * Creates StyleBox for listing out different styles inside the selection
     * @param selId - Id of the selection for which the style box is created
     * @param side - Alignment of text inside the box
     */
    static StyleBox(selId: string, side: string): void;

    /**
     * Add alert messages to different panels
     * @param config - Configuraiton for alert box display
     */
    static AlertBox(config: any): void;

    /**
     * Creates the panel for manipulation of labels on the viewport
     */
    static ContextMenu(): void;

    /**
     * Property object used in property menu
     * @param key - Name of the atom property
     * @param value - Value of the property
     */
    static Property(key: string, value: any): void;

    /**
     * Creates UI panel for surface manipulations
     */
    static SurfaceMenu(): void;

    /**
     * Creates cards for manipulation of surface
     */
    static Surface(): void;

    /**
     * Sets the css position property left and top for the element
     * @param jquery - html element
     * @param left - : css left property
     * @param top - : css top peroperty
     */
    static setPosition(jquery: any, left: number, top: number): void;

    /**
     * Sets the location of the element relative to the parseInt
     * as per position types
     * @param parent - jquery object
     * @param child - jquery object
     * @param x_type - 'left|right'
     * @param y_type - 'top|bottom'
     * @param x_offset - Offset x values in pixels
     * @param y_offset - Offset y values in pixels
     */
    static setLocation(parent: any, child: any, x_type: string, y_type: string, x_offset: number, y_offset: number): void;
  }
}