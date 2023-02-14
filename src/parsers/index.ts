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

export { bondLength, setBondLength } from './utils/bondLength';

// @category Parsers
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
} as Record<string, (str: string, options: ParserOptionsSpec) => any>;
