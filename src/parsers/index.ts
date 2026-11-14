import { ParserOptionsSpec } from './ParserOptionsSpec';
import { VASP } from "./VASP";
import { CUBE } from "./CUBE";
import { XYZ } from "./XYZ";
import { SDF } from "./SDF";
import { CDJSON } from "./CDJSON";
import { CIF } from "./CIF";
import { MOL2 } from "./MOL2";
import { PDB } from "./PDB";
import { PQR } from "./PQR";
import { MMTFparser } from "./MMTF";
import { PRMTOP } from "./PRMTOP";
import { GRO } from "./GRO";
import { LAMMPSTRJ } from "./LAMMPSTRJ";
import {BCIF} from "./BCIF";

export { bondLength, setBondLength } from './utils/bondLength';

// Keep this standalone typedef beside a runtime declaration so the documentation
// transform preserves it when TypeScript removes type-only imports and interfaces.
/**
* Molecular structure formats accepted by {@link GLViewer#addModel} and {@link GLModel#addMolData}.
* These are distinct from the density-map formats accepted by {@link VolumeData} and the
* trajectory formats accepted by {@link GLModel#setCoordinates}. Append .gz for gzip input.
* @typedef FileFormats
* @category Parsers
* @prop bcif Binary Crystallographic Information File format
* @prop cdjson,json ChemDoodle JSON format
* @prop cube Gaussian cube format
* @prop gro Gromacs structure format, including atom coordinates.
* @prop lammpstrj LAMMPS trajectory format
* @prop mcif,cif Crystallographic Information File, including macromolecular mmCIF
* @prop mmtf Macromolecular Transmission Format (binary data)
* @prop mol2 Sybyl Mol2 format
* @prop pdb The venerable Protein Data Bank format
* @prop pdbqt AutoDock PDBQT format
* @prop pqr PDB-like format with charge and radius, stored in atom properties. Charges are available as properties.charge and properties.partialCharge.
* @prop prmtop Amber topology file, must add coordinates
* @prop sdf MDL MOL/SDF format, including V2000, V3000, multiple models, and metadata
* @prop vasp VASP format (CONTCAR, POSCAR)
* @prop xyz XYZ cartesian coordinates format
*/

/** Structure parser registry, keyed by format names from {@link FileFormats}.
 * Most applications should use GLViewer.addModel or GLModel.addMolData, which
 * initialize model defaults and handle compressed input. Direct parsers return
 * arrays of atom arrays, with format-specific metadata attached to the result.
 * See {@link ParserOptionsSpec} for format-dependent options.
 * @category Parsers
 */
export const Parsers = {
  vasp: VASP,
  VASP,
  cube: CUBE,
  CUBE,
  xyz: XYZ,
  XYZ,
  sdf: SDF,
  SDF,
  json: CDJSON,
  cdjson: CDJSON,
  CDJSON,
  mcif: CIF,
  cif: CIF,
  CIF,
  mol2: MOL2,
  MOL2,
  pdb: PDB,
  PDB,
  pdbqt: PDB,
  PDBQT: PDB,
  pqr: PQR,
  PQR,
  mmtf:MMTFparser, //need to avoid name collision
  MMTF:MMTFparser,
  prmtop: PRMTOP,
  PRMTOP,
  gro: GRO,
  GRO,
  lammpstrj: LAMMPSTRJ,
  LAMMPSTRJ,
  bcif: BCIF,
  BCIF
} as Record<string, (str: string, options: ParserOptionsSpec) => any>;
