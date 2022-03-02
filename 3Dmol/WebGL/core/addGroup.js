// @ts-check

import { geometryGroup } from "./GeometryGroup";

const BUFFERSIZE = 65535; //limited to 16bit indices

export function addGroup(geo) {
  var ret = new geometryGroup(geo.geometryGroups.length);
  geo.geometryGroups.push(ret);
  geo.groups = geo.geometryGroups.length;

  ret.vertexArray = new Float32Array(BUFFERSIZE * 3);
  ret.colorArray = new Float32Array(BUFFERSIZE * 3);

  //TODO: instantiating uint arrays according to max number of vertices
  // is dangerous, since there exists the possibility that there will be 
  // more face or line indices than vertex points - but so far that doesn't
  // seem to be the case for any of the renders 
  if (geo.mesh) {
    ret.normalArray = new Float32Array(BUFFERSIZE * 3);
    ret.faceArray = new Uint16Array(BUFFERSIZE * 6);
    ret.lineArray = new Uint16Array(BUFFERSIZE * 6);
  }
  if (geo.radii) {
    ret.radiusArray = new Float32Array(BUFFERSIZE);
  }
  ret.useOffset = geo.offset;


  return ret;
};