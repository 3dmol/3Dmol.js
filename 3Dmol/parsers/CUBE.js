import { Vector3, Matrix4 } from "../WebGL/math";
import anumToSymbol from "./util/anumToSymbol";
import assignBonds from "./util/assignBonds";

/**
 * @param {string} str
 * @param {import("../specs").ParserOptionsSpec} options
 */
export default function parseCUBE(str, options) {
  options = options || {};
  const atoms = [[]];
  let lines = str.split(/\r?\n/);
  const assignbonds = options.assignBonds === undefined ? true : options.assignBonds;

  if (lines.length < 6) return atoms;

  let lineArr = lines[2].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');

  const natoms = Math.abs(parseFloat(lineArr[0]));

  const cryst = {};
  const origin = new Vector3(
    parseFloat(lineArr[1]),
    parseFloat(lineArr[2]),
    parseFloat(lineArr[3])
  );
  cryst.origin = origin;

  lineArr = lines[3].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
  lineArr = lines[3].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');

  // might have to convert from bohr units to angstroms
  // there is a great deal of confusion here:
  // n>0 means angstroms: http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm
  // n<0 means angstroms: http://paulbourke.net/dataformats/cube/
  // always assume bohr: openbabel source code
  // always assume angstrom: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
  // we are going to go with n<0 means angstrom - note this is just the first n
  // @ts-ignore
  const convFactor = lineArr[0] > 0 ? 0.529177 : 1;
  origin.multiplyScalar(convFactor);

  // @ts-ignore
  const nX = Math.abs(lineArr[0]);
  const xVec = new Vector3(
    parseFloat(lineArr[1]),
    parseFloat(lineArr[2]),
    parseFloat(lineArr[3])
  ).multiplyScalar(convFactor);

  lineArr = lines[4].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
  // @ts-ignore
  const nY = Math.abs(lineArr[0]);
  const yVec = new Vector3(
    parseFloat(lineArr[1]),
    parseFloat(lineArr[2]),
    parseFloat(lineArr[3])
  ).multiplyScalar(convFactor);

  lineArr = lines[5].replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
  // @ts-ignore
  const nZ = Math.abs(lineArr[0]);
  const zVec = new Vector3(
    parseFloat(lineArr[1]),
    parseFloat(lineArr[2]),
    parseFloat(lineArr[3])
  ).multiplyScalar(convFactor);

  cryst.size = {x: nX, y: nY, z: nZ};
  cryst.unit = new Vector3(xVec.x, yVec.y, zVec.z);

  // eslint-disable-next-line eqeqeq
  if (xVec.y != 0 || xVec.z != 0 || yVec.x != 0 || yVec.z != 0 || zVec.x != 0 || zVec.y != 0) {
    // need a transformation matrix
    cryst.matrix4 = new Matrix4(
      xVec.x,
      yVec.x,
      zVec.x,
      0,
      xVec.y,
      yVec.y,
      zVec.y,
      0,
      xVec.z,
      yVec.z,
      zVec.z,
      0,
      0,
      0,
      0,
      1
    );
    // include translation in matrix
    const t = new Matrix4().makeTranslation(origin.x, origin.y, origin.z);
    cryst.matrix4 = cryst.matrix4.multiplyMatrices(t, cryst.matrix4);
    cryst.matrix = cryst.matrix4.matrix3FromTopLeft();
    // all translation and scaling done by matrix, so reset origin and unit
    cryst.origin = new Vector3(0, 0, 0);
    cryst.unit = new Vector3(1, 1, 1);
  }

  // @ts-ignore
  atoms.modelData = [{cryst}];

  // Extract atom portion; send to new GLModel...
  lines = lines.splice(6, natoms);

  const start = atoms[atoms.length - 1].length;
  const end = start + lines.length;

  for (let i = start; i < end; ++i) {
    const atom = {};
    atom.serial = i;
    const line = lines[i - start];
    const tokens = line.replace(/^\s+/, '').replace(/\s+/g, ' ').split(' ');
    atom.elem = anumToSymbol[tokens[0]];
    atom.x = parseFloat(tokens[2]) * convFactor;
    atom.y = parseFloat(tokens[3]) * convFactor;
    atom.z = parseFloat(tokens[4]) * convFactor;

    atom.hetflag = true;
    atom.bonds = [];
    atom.bondOrder = [];
    atom.properties = {};
    // @ts-ignore
    atoms[atoms.length - 1].push(atom);
  }

  if (assignbonds) {
    for (let i = 0; i < atoms.length; i++) assignBonds(atoms[i]);
  }
  return atoms;
}
