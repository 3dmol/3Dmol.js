import { AtomStyleSpec } from "GLModel";


/**
  * Options for reading molecular structures with GLViewer.addModel or GLModel.addMolData.
  * Options are format-dependent; unsupported options have no effect. Use multimodel
  * to request multiple structures, frames to append them as trajectory frames to an
  * existing model, or onemol to combine them when the parser supports it. See {@link FileFormats}.
  * @category Parsers
  */
 export interface ParserOptionsSpec {
  /** On subsequent addMolData calls, append parsed structures as new frames instead of adding their atoms to the current frame. Does not itself enable reading multiple structures; use multimodel for that. */
  frames?: boolean;
  /** Generate vibration frames from atom displacement vectors when both frames and amplitude are supplied and nonzero. */
  vibrate?: {
    /** Number of vibration frames to create. */
    frames?: number;
    /** Scale applied to the atom displacement vectors. */
    amplitude?: number
  };
  /** Read multiple models when supported by the input format. */
  multimodel?: boolean;
  /** Combine parsed models into one atom list when supported by the format, typically together with multimodel. */
  onemol?: boolean;
  /** Retain hydrogens in parsers that can strip them. Defaults to false for pdb/pdbqt, mmtf, and bcif, and true for sdf and mol2. Other parsers may always retain hydrogens. */
  keepH?: boolean;
  /** used to define ChemDoodle styles ; supported by cdjson */
  parseStyle?: unknown;
  /** Apply assembly or symmetry transformations from the input; supported by cif/mcif, pdb, mmtf, and bcif. */
  doAssembly?: boolean;
  /**  Set to true if you wish to duplicate assembly atoms otherwise false ; supported by all formats with symmetries.  Not duplicating will result in faster rendering but it will not be possible to individually style symmetries. */
  duplicateAssemblyAtoms?: boolean;
  /** Set to true with duplicateAssemblyAtoms to individually wrap atoms (from symmetries) into unit cell  */
  wrapAtoms?: boolean;
  /** shift symmetry mates so their centroid is in the unit cell */
  normalizeAssembly?: boolean;
  /** do not detect bonds between symmetries generated with duplicateAssemblyAtoms (cif only - other formats never make bonds between symmetries) */
  dontConnectDuplicatedAtoms?: boolean;
  /** Do not compute secondary structure; supported by pdb, mmtf, and pqr. */
  noSecondaryStructure?: boolean;
  /** Do not compute secondary structure; supported by pdb, mmtf, cif/mcif, and bcif.
   * CIF does not read file annotations, so disabling computation leaves residues as coil. */
  noComputeSecondaryStructure?: boolean;
  /** Maximum hydrogen-bond distance in angstroms when computing secondary structure; supported by pdb, pqr, mmtf, cif, and bcif. */
  hbondCutoff?: number;
  /** Alternate location to select (default A), in addition to atoms without an alternate location. Use '*' to load all; supported by pdb/pdbqt, mmtf, and bcif. */
  altLoc?: string;
  /** index of the assembly in symmetry ; supported by mmtf */
  assemblyIndex?: number;
  /** Infer bonds from coordinates where supported. Defaults to true for pdb/pdbqt, xyz, cif/mcif, and bcif; LAMMPSTRJ requires true explicitly. Disabling inference does not remove bonds supplied by the file. */
  assignBonds?: boolean;
  /** When inferring bonds, do not assign bonds to common cations. */
  unboundCations?: boolean;
  /** Apply this style after parsing, replacing the default or file-supplied atom styles. */
  style?: AtomStyleSpec;
 };