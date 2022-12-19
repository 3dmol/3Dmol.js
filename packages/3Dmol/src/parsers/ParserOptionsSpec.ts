/**
* File formats supported by 3Dmol.js
* @typedef FileFormats
* @category Parsers
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


/**
  * Parser options specification. Used to specify the options of a GLModel.  Depending on the input file format, not all fields may be defined.
  * @typedef ParserOptionsSpec
  * @category Parsers
  * @prop {boolean} frames - true if you want to add to a new frame and false otherwise ; supported by all
  * @prop {object} vibrate - object specifying the vibration behavior ; supported by all
  * @prop {number} vibrate.frames - number of frames to be created, default to 10 ; supported by all
  * @prop {number} vibrate.amplitude -amplitude of distortion, default to 1 (full) ; supported by all
  * @prop {boolean} multimodel - specifies weather or not multiple models are being defined ; supported by xyz,sdf, or mol2
  * @prop {boolean} onemol -specifies weather or not the model is of one molecule ; Supported by xyz , sdf , mol2
  * @prop {boolean} keepH - do not strip hydrogens ; supported by sdf,mol2
  * @prop {object} parseStyle - used to define ChemDoodle styles ; supported by cdjson
  * @prop {boolean} doAssembly - boolean dictating weather or not to do assembly ; supported by mcif
  * @prop {boolean} duplicateAssemblyAtoms - Set to true if you wish to duplicate assembly atoms otherwise false ; supported by all formats with symmetries.  Not duplicating will result in faster rendering but it will not be possible to individually style symmetries.
  * @prop {boolean} normalizeAssembly - shift symmetry mates so their centroid is in the unit cell
  * @prop {boolean} dontConnectDuplicatedAtoms - do not detect bonds between symmetries generated with duplicateAssemblyAtoms (cif only - other formats never make bonds between symmetries)
  * @prop {boolean} noSecondaryStructure - boolean dictating the presence of a secondary structure ; supported by pdb
  * @prop {boolean} noComputeSecondaryStructure - do not compute ss ; supported by pdb, mmtf, cif
  * @prop {number} hbondCutoff - maximum distance used for identifying hydrogen bonds when computing secondary structure; supported by pdb, mmtf, cif
  * @prop {string} altLoc -which alternate location to select, if present; '*' to load all ; supported by pdb
  * @prop {number} assemblyIndex - index of the assembly in symmetry ; supported by mmtf
  * @prop {boolean} assignBonds - for formats without explicit bonds (e.g. PDB, xyz) infer bonding (default true). 
  */
 export type ParserOptionsSpec = Partial<{
  frames: boolean,
  vibrate: {
    frames: number,
    amplitude: number
  },
  multimodel: boolean,
  onemol: boolean,
  keepH: boolean,
  parseStyle: unknown,
  doAssembly: boolean,
  duplicateAssemblyAtoms: boolean,
  normalizeAssembly: boolean,
  dontConnectDuplicatedAtoms: boolean,
  noSecondaryStructure: boolean,
  noComputeSecondaryStructure: boolean,
  hbondCutoff: number,
  altLoc: string,
  assemblyIndex: number,
  assignBonds: boolean
 }>