import { ColorSpec } from "../colors";

/**
 * Atom representation. Depending on the input file format, not all fields may be defined.
 */
export type AtomSpec = Partial<{
  /** Parent residue name */
  resn: string;
  /**  Atom's x coordinate  */
  x: number;
  /**  Atom's y coordinate  */
  y: number;
  /**  Atom's z coordinate  */
  z: number;
  /**  Atom's color, as hex code or built-in color string */
  color: ColorSpec;
  /**  Hex code for color to be used for surface patch over this atom */
  surfaceColor: ColorSpec;
  /**  Element abbreviation (e.g. 'H', 'Ca', etc) */
  elem: string;
  /**  Set to true if atom is a heteroatom */
  hetflag: boolean;
  /**  Chain this atom belongs to, if specified in input file (e.g 'A' for chain A) */
  chain: string;
  /**  Residue number */
  resi: number;
  icode: number;
  rescode: number;
  /** Atom's serial id number */
  serial: number;
  /** Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon) */
  atom: string;
  /** Array of atom ids this atom is bonded to */
  bonds: number[];
  /** Secondary structure identifier (for cartoon render; e.g. 'h' for helix) */
  ss: string;
  /** true if this atom forms only single bonds or no bonds at all */
  singleBonds: boolean;
  /** Array of this atom's bond orders, corresponding to bonds identfied by 'bonds' */
  bondOrder: number[];
  /** Optional mapping of additional properties */
  properties: Record<string, any>;
  /** Atom b factor data */
  b: number;
  /** If applicable, this atom's record entry from the input PDB file (used to output new PDB from models) */
  pdbline: string;
  /** Set this flag to true to enable click selection handling for this atom */
  clickable: boolean;
  /** Callback click handler function to be executed on this atom and its parent viewer */
  callback: (
    atom: AtomSpec,
    viewer: unknown /** glviewer has not been ported yet */
  ) => void;
  /** for selection, inverts the meaning of the selection */
  invert: boolean;
}>;
