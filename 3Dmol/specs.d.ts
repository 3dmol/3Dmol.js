import { elementColors } from "./colors";
import CAP from "./enum/CAP";
import GLModel from "./GLModel";
import GLViewer from "./GLViewer";
import Gradient from "./Gradient";
import VolumeData from "./VolumeData";
import Color from "./WebGL/core/Color";
import { Matrix3, Matrix4, Vector2, Vector3 } from "./WebGL/math";

export type AnyFunc = (...args: any[] | undefined[]) => any

export type ColorSpec = Color | string | number;

/**
 * @typedef ColorschemeSpec
 * Built in colorschemes
 *
 * @example //Using a function in order to define the colors. 
  download("pdb:4UAA",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  var colorAsSnake = function(atom) {
                    return atom.resi % 2 ? 'white': 'green'
                  };

                  viewer.setStyle( {chain:'A'}, { cartoon: {colorfunc: colorAsSnake }});
                  viewer.setStyle( {chain:'B'}, { stick: {colorscheme: 'yellowCarbon'}});

                  viewer.render();
              });
 * @prop {string} <html color>Carbon   - use default element colors but with carbon set to specify html color string
 * @prop {string} ssPyMOL - PyMol secondary colorscheme
 * @prop {string} ssJmol - Jmol secondary colorscheme
 * @prop {string} Jmol - Jmol primary colorscheme
 * @prop {string} default - default colorscheme
 * @prop {string} amino - amino acid colorscheme
 * @prop {string} shapely - shapely protien colorscheme
 * @prop {string} nucleic - nucleic acid colorscheme
 * @prop {string} chain - standard chain colorscheme
 * @prop {string} chainHetatm - chain Hetatm colorscheme
 * @prop {string} prop - atomSpec property. Example 'b'. See AtomSpec.
 * @prop {Gradient} gradient - Allows the user to provide a gradient to the colorscheme.  Is either a Gradient object or the name of a built-in gradient (rwb, roygb, sinebow)
 * @prop {min} - min value for gradient
 * @prop {max} - max value for gradient
 * @prop {mid} - mid point value for gradient (for rwb)
 * @prop {object} map - map of a certain AtomSpec property to a color of the form `{'prop': 'elem', map:elementColors.greenCarbon}` Allows the user to provide a mapping of elements to colors to the colorscheme.  This can be done with any properties, and not just 'elem'.
 * @prop {function} colorfunc - Allows the user to provide a function for setting the colorschemes.
 */
export type ColorSchemeSpec = Partial<{
  ssPyMOL: string;
  ssJmol: string;
  Jmol: string;
  default: string;
  amino: string;
  shapely: string;
  nucleic: string;
  chain: string;
  chainHetatm: string;
  prop: unknown;
  Carbon: ColorSpec;
  gradient: Gradient;
  min: ColorSpec;
  mid: ColorSpec;
  max: ColorSpec;
  map: any;
  colorfunc: AnyFunc;
  colorscheme: string | ColorSchemeSpec;
  color: ColorSpec;
}>;

/**  GLViewer input specification */
export type ViewerSpec = Partial<{
  callback: () => any;
  defaultcolors: Record<string, number>;
  nomouse: boolean;
  backgroundColor: string | number;
  backgroundAlpha: number;
  camerax: string | number;
  hoverDuration: number;
  id: string;
  cartoonQuality: number;
  row: number;
  col: number;
  rows: number;
  cols: number;
  canvas: HTMLCanvasElement;
  viewers: Array<GLViewer>;
  minimumZoomToDistance: number;
  lowerZoomLimit: number;
  upperZoomLimit: number;
  antialias: boolean;
  controlAll: boolean;
  orthographic: boolean;
  disableFog: boolean;
  defaultColors: Record<string, number>;
  style: unknown;
}>;

/** Grid GLViewer input specification */
export type ViewerGridSpec = Partial<{
  rows: number;
  cols: number;
  row: number;
  col: number;
  viewers: Array<Array<GLViewer | null>>;
  canvas: HTMLCanvasElement;
  controlAll: boolean;
}>

export type GridSpec = Partial<{
  rows: number;
  cols: number;
  controlAll: boolean;
}>;

/** Atom representation. Depending on the input file format, not all fields may be defined. */
export type AtomSpec = {
  index: number;
  bonds: Array<number>;
  bondOrder: Array<number>;
  properties: Record<string, any>;
  elem: string;
  x: number;
  y: number;
  z: number;
  style: AtomStyleSpec;
} & Partial<{
  resn: string;
  surfaceColor: ColorSpec;
  hetflag: boolean;
  chain: string;
  resi: number;
  icode: number | string;
  rescode: number | string;
  serial: number;
  atom: string;
  ss: string;
  singleBonds: boolean;
  b: number | string;
  pdbline: string;
  clickable: boolean;
  hoverable: boolean;
  callback: (e: GLViewer) => any;
  invert: boolean;
  intersectionShape: any;
  bondStyles: any[];
  symmetries: Array<Vector3Like>;
  color: ColorSpec;
  hoverCallback: AnyFunc;
  unhoverCallback: AnyFunc;
  contextMenuEnabled: boolean;
  hbondOther: AtomSpec;
  hbondDistanceSq: number;
  altLoc: string;
  uMat: {
    u11: number;
    u22: number;
    u33: number;
    u12: number;
    u13: number;
    u23: number;
  }
  ssbegin:boolean;
  ssend:boolean;
  reschain:number;
}>;

/** 3 dimensional vector */
export type Vector3Like = Vector3 | { x: number, y: number, z: number };

/**
 * Atom selection object. Used to specify what atoms should be selected.  Can include
 * any field from {@link AtomSpec} in which case atoms must equal the specified value.
 * All fields must match for the selection to hold. If values
 * are provided as a list, then only one value of the list must match.
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
export type AtomSelectionSpec = Partial<{
  model: GLModel;
  bonds: number;
  predicate: (a: AtomSpec) => boolean;
  invert: boolean;
  byres: boolean;
  expand: number | string;
  within: WithinSelectionSpec;
  and: Array<AtomSelectionSpec>;
  or: Array<AtomSelectionSpec>;
  not: AtomSelectionSpec;
  clickable: boolean;
  hoverable: boolean;
  contextMenuEnabled: boolean;
}>

/**
 * Within selection object. Used to find the subset of an atom selection that is within
 * some distance from another atom selection. When added as a field of an {@link AtomSelectionSpec},
 * intersects the set of atoms in that selection with the set of atoms within a given
 * distance from the given {@link AtomSelectionSpec}.
 * @example
  $3Dmol.download("pdb:2EJ0",viewer,{},function(){
    viewer.setStyle({chain: 'A', within:{distance: 10, sel:{chain: 'B'}}}, {sphere:{}});
    viewer.render();
  });// stylizes atoms in chain A that are within 10 angstroms of an atom in chain B
 *
 */
export type WithinSelectionSpec = Partial<{
  distance: number | string;
  invert: boolean;
  sel: AtomSelectionSpec;
}>;


type LineStyleSpec = any;
type StickStyleSpec = any;
type SphereStyleSpec = any;
type ClickSphereStyleSpec = any;
export type AtomStyleSpec = Partial<{
  line: LineStyleSpec;
  cross: CrossStyleSpec;
  stick: StickStyleSpec;
  sphere: SphereStyleSpec;
  cartoon: CartoonStyleSpec;
  clicksphere: ClickSphereStyleSpec;
  radius: number;
  scale: number;
}>;


/**
 * @example
  var setStyles = function(volumedata){
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
export type SurfaceStyleSpec = Partial<{
  opacity: number;
  colorscheme: ColorSchemeSpec;
  color: ColorSpec;
  voldata: VolumeData;
  volscheme: Gradient;
  volformat: string;
  /** @deprecated */
  map: any;
}>;

/** Isosurface style specification */
export type IsoSurfaceSpec = Partial<{
  isoval: number;
  voxel: boolean;
  color: ColorSpec;
  opacity: number;
  wireframe: boolean;
  linewidth: number;
  smoothness: number;
  cords: Array<Vector3Like>;
  seldist: number;
  voldata: VolumeData;
  volscheme: Gradient;
  volformat: string;
  clickable: boolean;
  callback: AnyFunc;
  coords: Array<Vector3Like>;
  selectedRegion: Array<Vector3Like>;
  hoverable: boolean;
  radius: number;
  selectedOffset: number;
  transferfn: unknown;
}>;

/** VolumetricRenderer style specification */
export type VolumetricRendererSpec = Partial<{
  transferfn: Array<{ opacity: number; value: number; color: ColorSpec }>;
  subsamples: number;
  seldist: number;
  coords: Array<Vector3Like>;
}>;

/** GLShape style specification */
export type ShapeSpec = Partial<{
  color: ColorSpec;
  alpha: number;
  wireframe: boolean;
  hidden: boolean;
  linewidth: number;
  clickable: boolean;
  callback: AnyFunc | string
  frame: number;
  voldata: VolumeData;
  volscheme: Gradient;
  hoverable: boolean;
  opacity: number;
  side: number;
  hover_callback: AnyFunc | string
  unhover_callback: () => AnyFunc | string
}>;

/** Specification for adding custom shape. Extends {@link ShapeSpec}. */
export type CustomShapeSpec = {
  vertexArr: Array<Vector3Like>;
  normalArr: Array<Vector3Like>;
  faceArr: Array<number>;
} & Partial<{
  color: ColorSpec | Array<ColorSpec>;
}> & ShapeSpec;

/** Sphere shape specification. Extends {@link ShapeSpec} */
export type SphereShapeSpec = Partial<{
  center: Vector3Like;
  radius: number;
}> & ShapeSpec;

/** Box shape specification. Extends {@link ShapeSpec} */
export type BoxSpec = Partial<{
  corner: Vector3Like;
  center: Vector3Like;
  dimensions: { w: number | Vector3Like, h: number | Vector3Like, d: number | Vector3Like };
}> & ShapeSpec;

/** Arrow shape specification.  Extends {@link ShapeSpec} */
export type ArrowSpec = Partial<{
  start: Vector3Like;
  end: Vector3Like;
  radius: number;
  color: ColorSpec;
  hidden: boolean;
  radiusRatio: number;
  mid: number;
  midpos: number;
  dir: Vector3Like;
  length: number;
}> & ShapeSpec;

/** Cylinder shape specification.  Extends {@link ShapeSpec} */
export type CylinderSpec = Partial<{
  start: Vector3Like;
  end: Vector3Like;
  radius: number;
  fromCap: number;
  toCap: number;
  dashed: boolean;
  dashLength: number;
  gapLength: number;
}> & ShapeSpec;

/** Curve shape specification.  Extends {@link ShapeSpec} */
export type CurveSpec = Partial<{
  points: Array<Vector3Like>;
  smooth: number;
  radius: number;
  fromArrow: boolean;
  toArrow: boolean;
  fromCap: number;
  toCap: number;
}> & ShapeSpec;

/** Line shape specification.  Extends {@link ShapeSpec}  (but defaults to wireframe) */
export type LineSpec = Partial<{
  start: Vector3Like;
  end: Vector3Like;
  dashed: boolean;
}> & ShapeSpec;

/**
* File formats supported by 3Dmol.js
* - cdjson,json  Chemical JSON format
* - cube Gaussian cube format
* - gro  Gromacs topology format, need to add coordinates to resulting model.
* - mcif,cif Crystallographic Information File, the successor to PDB that makes you miss the PDB file format
* - mmtf Macromolecular Transmission Format, the successor to PDB that is totally awesome
* - mol2 Sybyl Mol2 format
* - pdb The venerable Protein Data Bank format
* - pqr Like PDB but with partial charges which are read into the partialcharge atom property
* - prmtop Amber topology file, must add coordinates
* - sdf MDL MOL format, supports multiple models and meta data
* - vasp VASP format (CONTCAR, POSCAR)
* - xyz XYZ cartesian coordinates format
*/
export type FileFomats =
  "cdjson" |
  "cube" |
  "gro" |
  "mcif" |
  "mmtf" |
  "mol2" |
  "pdb" |
  "pqr" |
  "prmtop" |
  "sdf" |
  "vasp" |
  "xyz" |
  "json" |
  "cif"

/**
* Parser options specification. Used to specify the options of a GLModel.  Depending on the input file format, not all fields may be defined.
*/
export type ParserOptionsSpec = Partial<{
  frames: boolean;
  vibrate: Partial<{
    frames: number;
    amplitude: number;
  }>;
  multimodel: boolean;
  onemol: boolean;
  keepH: boolean;
  parseStyle: any;
  doAssembly: boolean;
  duplicateAssemblyAtoms: boolean;
  dontConnectDuplicateAtoms: boolean;
  normalizeAssembly: boolean;
  noSecondaryStructure: boolean;
  noComputeSecondaryStructure: boolean;
  altLoc: string;
  assemblyIndex: number;
  assignBonds: boolean;
  style: AtomStyleSpec;
  cartoonQuality: number;
  defaultcolors: Record<string, number>;
  dontConnectDuplicatedAtoms: boolean;
}>;

/**
 * A visualization of protein or nucleic acid secondary structure.  Applying this to other molecules will not show anything.
 * @example download("pdb:4ZD3",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  viewer.setViewStyle({style:"outline"});
                  viewer.setStyle({},{cartoon:{}});
                  viewer.render();
              });
 */
export type CartoonStyleSpec = Partial<{
  /** strand color, may specify as 'spectrum' which will apply reversed gradient based on residue number */
  color: ColorSpec;
  /** style of cartoon rendering (trace, oval, rectangle (default), parabola, edged) */
  style: string;
  /** whether to use constant strand width, disregarding secondary structure; use thickness to adjust radius */
  ribbon: boolean;
  /** whether to add arrows showing beta-sheet directionality; does not apply to trace or ribbon */
  arrows: boolean;
  /** whether to display alpha helices as simple cylinders; does not apply to trace */
  tubes: boolean;
  /** cartoon strand thickness, default is 0.4 */
  thickness: number;
  /** cartoon strand width, default is secondary structure-dependent; does not apply to trace or ribbon */
  width: number;
  /** set opacity from 0-1; transparency is set per-chain with a warning outputted in the event of ambiguity */
  opacity: number;
  hidden: boolean;
}>;

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
  export type CrossStyleSpec = Partial<{
    hidden: boolean;
    linewidth: number;
    radius: number;
    scale: number;
    colorscheme: ColorSchemeSpec;
    color: ColorSpec;
    opacity: number;
  }>

/**
  * Label type specification
  * @prop {string} font - font name, default sans-serif
  * @prop {number} fontSize - height of text, default 18
  * @prop {import('./colors').ColorSpec} fontColor - font color, default white
  * @prop {number} fontOpacity - font opacity, default 1
  * @prop {number} borderThickness - line width of border around label, default 0
  * @prop {import('./colors').ColorSpec} borderColor - color of border, default backgroundColor
  * @prop {string} borderOpacity - color of border
  * @prop {import('./colors').ColorSpec} backgroundColor - color of background, default black
  * @prop {string} backgroundOpacity - opacity of background, default 1
  * @prop {import("./WebGL/math").Vector3} position - x,y,z coordinates for label
  * @prop {import("./WebGL/math").Vector2} screenOffset - x,y _pixel_ offset of label from position
  * @prop {boolean} inFront - always put labels in from of model
  * @prop {boolean} showBackground - show background rounded rectangle, default true
  * @prop {boolean} fixed - sets the label to change with the model when zooming
  * @prop {boolean} useScreen - position is in screen (not model) coordinates which are pixel offsets from upper left corner.
  * @prop {Object} backgroundImage - An element to draw into the label.  Any CanvasImageSource is allowed.
  * @prop {string} alignment - how to orient the label w/respect to position: topLeft (default), topCenter, topRight, centerLeft, center, centerRight, bottomLeft, bottomCenter, bottomRight
  * @prop {number} frame - if set, only display in this frame of an animation
  */
export type LabelSpec = Partial<{
  font: string;
  fontSize: number | string;
  fontColor: ColorSpec;
  fontOpacity: number;
  borderThickness: number;
  borderColor: ColorSpec;
  borderOpacity: string;
  backgroundColor: ColorSpec;
  backgroundOpacity: string;
  position: Vector3Like;
  screenOffset: Vector2;
  inFront: boolean | 'false' | '0';
  showBackground: boolean | string | number;
  fixed: boolean;
  useScreen: boolean;
  backgroundImage: CanvasImageSource;
  alignment: string;
  frame: number;
  padding: number;
  bold: boolean;
  backgroundWidth: number;
  backgroundHeight: number;
  backgroundGradient: any;
}>;

/**
 * Unit Cell shape specification.
 * @typedef UnitCellStyleSpec
 * @prop {import('./specs').LineStyleSpec} box - line style used to draw box
 * @prop {import('./specs').ArrowSpec} astyle - arrow specification of "a" axis
 * @prop {import('./specs').ArrowSpec} bstyle - arrow specification of "b" axis
 * @prop {import('./specs').ArrowSpec} cstyle - arrow specification of "c" axis
 * @prop {string} alabel - label for a axis
 * @prop {import('./specs').LabelSpec} alabelstyle - label style for a axis
 * @prop {string} blabel - label for b axis
 * @prop {import('./specs').LabelSpec} blabelstyle - label style for b axis
 * @prop {string} clabel - label for c axis
 * @prop {import('./specs').LabelSpec} clabelstyle - label style for c axis
 */
export type UnitCellStyleSpec = Partial<{
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
}>;

export type ModelDataType = {
  symmetries?: Matrix4[];
  cryst?: {
    matrix?: Matrix3
  }
}


export type ParserResult = Array<Array<Partial<AtomSpec>>>
  & {
    modelData?: Array<ModelDataType>; 
    origin?: Vector3;
    box?: [number, number, number];
  };

