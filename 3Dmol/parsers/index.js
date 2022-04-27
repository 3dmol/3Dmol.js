// @ts-nocheck
/* eslint-disable vars-on-top */
/* eslint-disable eqeqeq */
import parseCDJSON from './CDJSON';
import parseCIF from './CIF';
import parseCUBE from './CUBE';
import parseGRO from './GRO';
import parseLAMMPSTRJ from './LAMMPSTRJ';
import parseMMTF from './MMTF';
import parseMOL2 from './MOL2';
import parsePDB from './PDB';
import parsePRQ from './PQR';
import parsePRMTOP from './PRMTOP';
import parseSDF from './SDF';
import ParseVASP from './VASP';
import parseXYZ from './XYZ';

/**
 * Parsers stores functions for parsing molecular data. They all take a string of molecular data
 * and options. The default behavior is to only read the first model in the case of multimodel files, and
 * all parsers return a list of atom list(s)
 *
 * Parsers.<ext> corresponds to the parsers for files with extension ext
 */
const Parsers = {
  vasp: ParseVASP,
  VASP: ParseVASP,
  cube: parseCUBE,
  CUBE: parseCUBE,
  xyz: parseXYZ,
  XYZ: parseXYZ,
  sdf: parseSDF,
  SDF: parseSDF,
  cdjson: parseCDJSON,
  json: parseCDJSON,
  mcif: parseCIF,
  cif: parseCIF,
  mol2: parseMOL2,
  MOL2: parseMOL2,
  pdb: parsePDB,
  PDB: parsePDB,
  pdbqt: parsePDB,
  PDBQT: parsePDB,
  prq: parsePRQ,
  PQR: parsePRQ,
  mmtf: parseMMTF,
  MMTF: parseMMTF,
  prmtop: parsePRMTOP,
  PRMTOP: parsePRMTOP,
  gro: parseGRO,
  GRO: parseGRO,
  lammpstrj: parseLAMMPSTRJ,
  LAMMPSTRJ: parseLAMMPSTRJ,
};

export default Parsers;