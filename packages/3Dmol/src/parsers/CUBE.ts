import { Vector3, Matrix4 } from "../WebGL";
import { assignBonds } from "./utils/assignBonds";
import { anumToSymbol } from "./utils/anumToSymbol";

/**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options
 * @category Parsers
     
*/
export function CUBE(str, options) {
  options = options || {};
  var atoms: any[][] & Record<string, any> = [[]];
  var lines = str.split(/\r?\n/);
  var assignbonds = options.assignBonds === undefined ? true : options.assignBonds;

  if (lines.length < 6)
      return atoms;

  var lineArr = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");

  var natoms = Math.abs(parseFloat(lineArr[0]));

  let cryst: Record<string, any> = {};
  var origin = cryst.origin = new Vector3(parseFloat(lineArr[1]),
          parseFloat(lineArr[2]), parseFloat(lineArr[3]));

  lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
  lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
 
  // might have to convert from bohr units to angstroms
  // there is a great deal of confusion here:
  // n>0 means angstroms: http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm
  // n<0 means angstroms: http://paulbourke.net/dataformats/cube/
  // always assume bohr: openbabel source code
  // always assume angstrom: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
  // we are going to go with n<0 means angstrom - note this is just the first n
  var convFactor = (lineArr[0] > 0) ? 0.529177 : 1;
  origin.multiplyScalar(convFactor);

  var nX = Math.abs(lineArr[0]);
  var xVec = new Vector3(parseFloat(lineArr[1]),
          parseFloat(lineArr[2]), parseFloat(lineArr[3]))
          .multiplyScalar(convFactor);

  lineArr = lines[4].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
  var nY = Math.abs(lineArr[0]);
  var yVec = new Vector3(parseFloat(lineArr[1]),
          parseFloat(lineArr[2]), parseFloat(lineArr[3]))
          .multiplyScalar(convFactor);

  lineArr = lines[5].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
  var nZ = Math.abs(lineArr[0]);
  var zVec = new Vector3(parseFloat(lineArr[1]),
          parseFloat(lineArr[2]), parseFloat(lineArr[3]))
          .multiplyScalar(convFactor);

  cryst.size = {x:nX, y:nY, z:nZ};
  cryst.unit = new Vector3(xVec.x, yVec.y, zVec.z);

  if (xVec.y != 0 || xVec.z != 0 || yVec.x != 0 || yVec.z != 0 || zVec.x != 0
          || zVec.y != 0) {
      //need a transformation matrix
      cryst.matrix4 =  new Matrix4(xVec.x, yVec.x, zVec.x, 0, xVec.y, yVec.y, zVec.y, 0, xVec.z, yVec.z, zVec.z, 0, 0,0,0,1);
      // include translation in matrix
      let t = new Matrix4().makeTranslation(origin.x, origin.y, origin.z);
      cryst.matrix4 = cryst.matrix4.multiplyMatrices(t,cryst.matrix4);
      cryst.matrix = cryst.matrix4.matrix3FromTopLeft();
      // all translation and scaling done by matrix, so reset origin and unit
      cryst.origin = new Vector3(0,0,0);
      cryst.unit = new Vector3(1,1,1);
  }

  atoms.modelData = [{cryst:cryst}];


  // Extract atom portion; send to new GLModel...
  lines = lines.splice(6, natoms);

  var start = atoms[atoms.length-1].length;
  var end = start + lines.length;

  for (var i = start; i < end; ++i) {
      var atom: Record<string, any> = {};
      atom.serial = i;
      var line = lines[i - start];
      var tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
              " ");
      atom.elem = anumToSymbol[tokens[0]];
      atom.x = parseFloat(tokens[2]) * convFactor;
      atom.y = parseFloat(tokens[3]) * convFactor;
      atom.z = parseFloat(tokens[4]) * convFactor;

      atom.hetflag = true;
      atom.bonds = [];
      atom.bondOrder = [];
      atom.properties = {};
      atoms[atoms.length-1].push(atom);

  }
  
  if(assignbonds) {
      for (let i = 0; i < atoms.length; i++)
          assignBonds(atoms[i]);
  }
  return atoms;
};