// Specifications for various object types used in 3Dmol.js
// This is primarily for documentation

import { Matrix4, Vector3 } from "WebGL";
import { AtomStyleSpec, BondStyle, GLModel } from "./GLModel";
import { GLViewer } from "./GLViewer";
import type { GLShape } from "./GLShape";
import { ColorSpec } from "./colors";

/** Ray intersection with a clickable or hoverable atom or shape. */
export interface PickingIntersection {
  /** Atom or shape intersected by the picking ray. */
  clickable: AtomSpec | GLShape;
  /** Distance along the picking ray. */
  distance: number;
}

/**
 * Atom representation. Depending on the input file format, not all fields may be defined.
 */
export interface AtomSpec {
  /** Parent residue name */
  resn?: string;
  /** Residue label name */
  lresn?: string;
  /**  Atom's x coordinate  */
  x?: number;
  /**  Atom's y coordinate  */
  y?: number;
  /**  Atom's z coordinate  */
  z?: number;
  /**  Atom's color, as hex code or built-in color string */
  color?: ColorSpec;
  /**  Hex code for color to be used for surface patch over this atom */
  surfaceColor?: ColorSpec;
  /**  Element abbreviation (e.g. 'H', 'Ca', etc) */
  elem?: string;
  /**  Set to true if atom is a heteroatom */
  hetflag?: boolean;
  /**  Chain this atom belongs to, if specified in input file (e.g 'A' for chain A) */
  chain?: string;
  /** Label chain (not author) */
  lchain?: string;
  /**  Residue number */
  resi?: number;
  /** Label residues number (not author) */
  lresi?: number; 
  icode?: string;
  rescode?: string;
  /** Atom's serial id number */
  serial?: number;
  /** Index of atom in molecule */
  index?: number;
  /** Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon) */
  atom?: string;
  /** Indices of bonded atoms in this model's atom array (not input-file serial numbers). */
  bonds?: number[];
  /** Secondary structure identifier (for cartoon render; e.g. 'h' for helix) */
  ss?: string;
  ssbegin?: boolean;
  ssend?: boolean;
  /** true if this atom forms only single bonds or no bonds at all */
  singleBonds?: boolean;
  /** Array of this atom's bond orders, corresponding to the atom indices in 'bonds' */
  bondOrder?: number[];
  /** Optional mapping of additional properties */
  properties?: Record<string, any>;
  /** Atom b factor data */
  b?: number;
  /** If applicable, this atom's record entry from the input PDB file (used to output new PDB from models) */
  pdbline?: string;
  /** Set this flag to true to enable click selection handling for this atom */
  clickable?: boolean;
  /** Click handler receiving the atom, viewer, mouse or touch event, viewer container, and intersections sorted nearest first. */
  callback?: (atom: AtomSpec, viewer: GLViewer, event: MouseEvent | TouchEvent,
    container: HTMLElement | null, intersects: PickingIntersection[]) => void;
  /** Set this flag to true to enable hover selection handling for this atom */
  hoverable?: boolean;
  /** Hover handler receiving the atom, viewer, mouse or touch event, viewer container, and intersections sorted nearest first. */
  hover_callback?: (atom: AtomSpec, viewer: GLViewer, event: MouseEvent | TouchEvent,
    container: HTMLElement | null, intersects: PickingIntersection[]) => void;
  /** Unhover handler with the same arguments as the hover handler.
   * The event and intersections can be undefined when hover is cleared. */
  unhover_callback?: (atom: AtomSpec, viewer: GLViewer, event: MouseEvent | TouchEvent | undefined,
    container: HTMLElement | null, intersects: PickingIntersection[] | undefined) => void;
  /** for selection, inverts the meaning of the selection */
  invert?: boolean;
  /** style of atom */
  style?: AtomStyleSpec;
  /** custom bond styling by position in bonds */
  bondStyles?: BondStyle[];
  intersectionShape?: any;
  capDrawn?: boolean;
  model?: number;
  contextMenuEnabled?: boolean;
  hbondDistanceSq?: number;
  hbondOther?: any;
  altLoc?: string;
  reschain?: number;
  uMat?: Record<string, number>;
  symmetries?: Vector3[];
  sym?: any;
  dx?: number;
  dy?: number;
  dz?: number;
}

// Ordinary scalar atom fields can match either one value or a list of values.
// Preserve structured and array-valued atom fields, which use different semantics.
type AtomSelectionValue<T> = [Extract<T, string | number | boolean>] extends [never]
  ? T
  : T | Extract<T, string | number | boolean>[];

type AtomSelectionFields = {
  [K in keyof AtomSpec]: AtomSelectionValue<AtomSpec[K]>;
};

/**
 * Atom selection object. Used to specify what atoms should be selected.  Can include
 * any field from {@link AtomSpec} in which case atoms must equal the specified value.
 * All fields must match for the selection to hold. Ordinary string, number, and
 * boolean atom fields also accept lists; only one value of each list must match.
 * Selection controls such as invert and byres keep their individual types.
 * @extends AtomSpec
 * @example
 * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
 *  viewer.setStyle({chain:'B'},{cartoon:{color:'spectrum'}});
 *  viewer.setStyle({chain:'B',invert:true},{cartoon:{}});
 *  viewer.setStyle({bonds: 0},{sphere:{radius:0.5}}); //water molecules
 *  viewer.setStyle({resn:'PMP',byres:true,expand:5},{stick:{colorscheme:"greenCarbon"}});
 *  viewer.setStyle({resi:["91-95","42-50"]},{cartoon:{color:"green",thickness:1.0}});
 *  viewer.render();
 * });
 */
export interface AtomSelectionSpec
  extends Omit<AtomSelectionFields, "bonds" | "model" | "index" | "resi" | "invert"> {
  /** a single model or list of models from which atoms should be selected.  Can also specify by numerical creation order.  Reverse indexing is allowed (-1 specifies last added model). */
  model?: GLModel | number | GLModel[] | number[];
  /** frame index of individual frame to style; will apply to all frames if not set */
  frame?: number;
  /** index of the atom or atoms to select */
  index?: number | number[];
  /** overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms */
  bonds?: number;
  /** overloaded to allow ranges and lists of residues, e.g. {resi: "5-10"} or {resi: [5,10,32]} */
  resi?: number | SelectionRange | (number | SelectionRange)[];
  /** user supplied function that gets passed an {@link AtomSpec} and should return true if the atom should be selected */
  predicate?: (atom: AtomSpec) => boolean;
  /** if set, inverts the meaning of the selection */
  invert?: boolean;
  /** if set, expands the selection to include all atoms of any residue that has any atom selected */
  byres?: boolean;
  /** expands the selection to include all atoms within a given distance from the selection */
  expand?: number | string;
  /** intersects the selection with the set of atoms within a given distance from another selection */
  within?: WithinSelectionSpec;
  /** take the intersection of the provided lists of {@link AtomSelectionSpec}s */
  and?: AtomSelectionSpec[] & { __cached_results?: any };
  /** take the union of the provided lists of {@link AtomSelectionSpec}s */
  or?: AtomSelectionSpec[] & { __cached_results?: any };
  /** take the inverse of the provided {@link AtomSelectionSpec} */
  not?: AtomSelectionSpec;
}

/**
 * Within selection object. Used to find the subset of an atom selection that is within
 * some distance from another atom selection. When added as a field of an {@link AtomSelectionSpec},
 * intersects the set of atoms in that selection with the set of atoms within a given
 * distance from the given {@link AtomSelectionSpec}.
 * @example
 * $3Dmol.download("pdb:2EJ0",viewer,{},function(){
 *  viewer.setStyle({chain: 'A', within:{distance: 10, sel:{chain: 'B'}}}, {sphere:{}});
 *  viewer.render();
 * });// stylizes atoms in chain A that are within 10 angstroms of an atom in chain B
 *
 */
export interface WithinSelectionSpec {
  /** the distance in angstroms away from the atom selection to include atoms in the parent selection */
  distance?: number;
  /** if set, selects atoms not within distance range for intersection */
  invert?: boolean;
  /** the selection of atoms against which to measure the distance from the parent atom selection */
  sel?: AtomSelectionSpec;
}

export type Cryst = {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
  origin: Vector3;
  size: { x: number; y: number; z: number };
  unit: Vector3;
  matrix4: Matrix4;
  matrix: unknown;
};
export type SelectionRange = `${number}-${number}`;
