/**
* File formats supported by 3Dmol.js
* @typedef FileFormats
* @category Parsers
* @prop bcif Binary Crystallographic Information File format
* @prop cdjson,json  Chemical JSON format
* @prop cube Gaussian cube format
* @prop gro Gromacs structure format, including atom coordinates.
* @prop lammpstrj LAMMPS trajectory format
* @prop mcif,cif Crystallographic Information File, the successor to PDB that makes you miss the PDB file format
* @prop mmtf Macromolecular Transmission Format, the successor to PDB that is totally awesome
* @prop mol2 Sybyl Mol2 format
* @prop pdb The venerable Protein Data Bank format
* @prop pdbqt AutoDock PDBQT format
* @prop pqr PDB-like format with charge and radius, stored in atom properties. Charges are available as properties.charge and properties.partialCharge.
* @prop prmtop Amber topology file, must add coordinates
* @prop sdf MDL MOL format, supports multiple models and meta data
* @prop vasp VASP format (CONTCAR, POSCAR)
* @prop xyz XYZ cartesian coordinates format
*/

import { AtomStyleSpec } from "GLModel";


/**
  * Parser options specification. Used to specify the options of a GLModel.  Depending on the input file format, not all fields may be defined.
  * @category Parsers
  */
 export interface ParserOptionsSpec {
  /** true if you want to add to a new frame and false otherwise ; supported by all */
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
  /** Retain hydrogens in parsers that can strip them, including pdb, pdbqt, sdf, mol2, mmtf, and bcif. */
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
  /** Alternate location to select, if present; '*' loads all; supported by pdb, mmtf, and bcif. */
  altLoc?: string;
  /** index of the assembly in symmetry ; supported by mmtf */
  assemblyIndex?: number;
  /** for formats without explicit bonds (e.g. PDB, xyz) infer bonding (default true).  */
  assignBonds?: boolean;
  /** for formats without explicit bonds, if assigning bonds to not assign bonds to common cations */
  unboundCations?: boolean;
  /** set model to this style after parsing */
  style?: AtomStyleSpec;
 };