// @ts-check
import { Material } from "../materials/Material";
import { LineBasicMaterial } from "../materials/LineBasicMaterial";
import { Vector3 } from "../math/Vector3";
import { Color } from "./Color";
import { Geometry } from "./Geometry";

const BUFFERSIZE = 65535; //limited to 16bit indices

/**
 * @param {Geometry} geo
 */
export function addGroup(geo) {
  var ret = new GeometryGroup(geo.geometryGroups.length);
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
}

export class GeometryGroup {
  /**
   * @type{WebGLBuffer | null}
   */
  __webglFaceBuffer = null;
  /**
   * @type{WebGLBuffer | null}
   */
  __webglLineBuffer = null;
  /**
   * @type{WebGLBuffer | null}
   */
  __webglRadiusBuffer = null;
  /**
   * @type{WebGLBuffer | null}
   */
  __webglOffsetBuffer = null;
  /**
   * @type{WebGLBuffer | null}
   */
  __webglVertexBuffer = null;
  /**
   * @type{WebGLBuffer | null}
   */
  __webglNormalBuffer = null;
  /**
   * @type{WebGLBuffer | null}
   */
  __webglColorBuffer = null;
  /**
   * @type {any}
   */
  useOffset;
  /**
   * @param {number} id
   */
  constructor(id) {
    this.id = id || 0;
    //for performance reasons, callers must directly modify these
    /**
     *
     */
    this.vertexArray = null;
    this.colorArray = null;
    this.normalArray = null;
    this.faceArray = null;
    this.radiusArray = null;
    //this.adjFaceArray=null;
    this.lineArray = null;
    this.vertices = 0;
    this.faceidx = 0;
    this.lineidx = 0;
  }

  /**
   * @param {(arg0: any, arg1: any, arg2: any) => any} setcolor
   */
  setColors(setcolor) {
    //apply a function that takes the vertex coordinate and returns a color
    var v = this.vertexArray;
    var c = this.colorArray;
    if (v.length != c.length) {
      console.log("Cannot re-color geometry group due to mismatched lengths.");
      return;
    }
    for (var i = 0; i < v.length; i += 3) {
      var col = setcolor(v[i], v[i + 1], v[i + 2]);
      if (!(col instanceof Color)) {
        col = globalThis.$3Dmol.CC.color(col);
      }
      c[i] = col.r;
      c[i + 1] = col.g;
      c[i + 2] = col.b;
    }
  }

  getNumVertices() {
    return this.vertices;
  }

  getVertices() {
    return this.vertexArray;
  }

  getCentroid() {
    var centroid = new Vector3();
    var offset, x, y, z;

    for (var i = 0; i < this.vertices; ++i) {
      offset = i * 3;

      x = this.vertexArray[offset];
      y = this.vertexArray[offset + 1];
      z = this.vertexArray[offset + 2];

      centroid.x += x;
      centroid.y += y;
      centroid.z += z;
    }

    //divideScalar checks for 0 denom
    centroid.divideScalar(this.vertices);

    return centroid;
  }

  setNormals() {
    var faces = this.faceArray;
    var verts = this.vertexArray;
    var norms = this.normalArray;

    if (!this.vertices || !this.faceidx) return;

    //vertex indices
    var a,
      b,
      c,
      //and actual vertices
      vA,
      vB,
      vC,
      norm;

    for (var i = 0; i < faces.length / 3; ++i) {
      a = faces[i * 3] * 3;
      b = faces[i * 3 + 1] * 3;
      c = faces[i * 3 + 2] * 3;

      vA = new Vector3(verts[a], verts[a + 1], verts[a + 2]);
      vB = new Vector3(verts[b], verts[b + 1], verts[b + 2]);
      vC = new Vector3(verts[c], verts[c + 1], verts[c + 2]);

      vA.subVectors(vA, vB);
      vC.subVectors(vC, vB);
      vC.cross(vA);

      //face normal
      norm = vC;
      norm.normalize();

      norms[a] += norm.x;
      norms[b] += norm.x;
      norms[c] += norm.x;
      norms[a + 1] += norm.y;
      norms[b + 1] += norm.y;
      norms[c + 1] += norm.y;
      norms[a + 2] += norm.z;
      norms[b + 2] += norm.z;
      norms[c + 2] += norm.z;
    }
  }

  setLineIndices() {
    if (!this.faceidx) return;

    if (
      this.lineArray &&
      this.lineArray.length == this.faceidx * 2 &&
      this.lineidx == this.faceidx * 2
    )
      return; //assume already computed

    var faceArr = this.faceArray,
      lineArr = (this.lineArray = new Uint16Array(this.faceidx * 2));
    this.lineidx = this.faceidx * 2;

    for (var i = 0; i < this.faceidx / 3; ++i) {
      var faceoffset = i * 3;
      var lineoffset = faceoffset * 2;
      var a = faceArr[faceoffset],
        b = faceArr[faceoffset + 1],
        c = faceArr[faceoffset + 2];

      lineArr[lineoffset] = a;
      lineArr[lineoffset + 1] = b;
      lineArr[lineoffset + 2] = a;
      lineArr[lineoffset + 3] = c;
      lineArr[lineoffset + 4] = b;
      lineArr[lineoffset + 5] = c;
    }
  }

  /**
   * @param {string} indent
   * @param {Material} material
   */
  vrml(indent, material) {
    var ret = "";
    ret +=
      indent +
      "Shape {\n" +
      indent +
      " appearance Appearance {\n" +
      indent +
      "  material Material {\n" +
      indent +
      "   diffuseColor " +
      // @ts-ignore
      material.color.r +
      " " +
      // @ts-ignore
      material.color.g +
      " " +
      // @ts-ignore
      material.color.b +
      "\n";
    if (material.transparent) {
      ret += indent + "   transparency " + (1.0 - material.opacity) + "\n";
    }
    ret += indent + "  }\n"; //material
    ret += indent + " }\n"; //appearance

    var oldindent = indent;
    indent += " "; //inshape
    if (material instanceof LineBasicMaterial) {
      ret +=
        indent +
        "geometry IndexedLineSet {\n" +
        indent +
        " colorPerVertex TRUE\n" +
        indent +
        " coord Coordinate {\n" +
        indent +
        "  point [\n";
      let x, y, z;
      for (let i = 0; i < this.vertices; ++i) {
        let offset = i * 3;
        x = this.vertexArray[offset];
        y = this.vertexArray[offset + 1];
        z = this.vertexArray[offset + 2];
        ret += indent + "   " + x + " " + y + " " + z + ",\n";
      }
      ret += indent + "  ]\n";
      ret += indent + " }\n"; //end coordinate

      if (this.colorArray) {
        ret += indent + " color Color {\n" + indent + "  color [\n";
        for (let i = 0; i < this.vertices; ++i) {
          let offset = i * 3;
          x = this.colorArray[offset];
          y = this.colorArray[offset + 1];
          z = this.colorArray[offset + 2];
          ret += indent + "   " + x + " " + y + " " + z + ",\n";
        }
        ret += indent + "  ]\n";
        ret += indent + " }\n"; //end color
      }

      ret += indent + " coordIndex [\n";
      for (let i = 0; i < this.vertices; i += 2) {
        ret += indent + "  " + i + ", " + (i + 1) + ", -1,\n";
      }
      ret += indent + " ]\n";
      ret += indent + "}\n"; //geometry
    } else {
      //faces
      ret +=
        indent +
        "geometry IndexedFaceSet {\n" +
        indent +
        " colorPerVertex TRUE\n" +
        indent +
        " normalPerVertex TRUE\n" +
        indent +
        " solid FALSE\n";

      //vertices
      ret += indent + " coord Coordinate {\n" + indent + "  point [\n";
      let x, y, z;
      for (let i = 0; i < this.vertices; ++i) {
        let offset = i * 3;
        x = this.vertexArray[offset];
        y = this.vertexArray[offset + 1];
        z = this.vertexArray[offset + 2];
        ret += indent + "   " + x + " " + y + " " + z + ",\n";
      }
      ret += indent + "  ]\n";
      ret += indent + " }\n"; //end coordinate

      //normals
      ret += indent + " normal Normal {\n" + indent + "  vector [\n";
      for (let i = 0; i < this.vertices; ++i) {
        let offset = i * 3;
        x = this.normalArray[offset];
        y = this.normalArray[offset + 1];
        z = this.normalArray[offset + 2];
        ret += indent + "   " + x + " " + y + " " + z + ",\n";
      }
      ret += indent + "  ]\n";
      ret += indent + " }\n"; //end normal

      //colors
      if (this.colorArray) {
        ret += indent + " color Color {\n" + indent + "  color [\n";
        for (let i = 0; i < this.vertices; ++i) {
          let offset = i * 3;
          x = this.colorArray[offset];
          y = this.colorArray[offset + 1];
          z = this.colorArray[offset + 2];
          ret += indent + "   " + x + " " + y + " " + z + ",\n";
        }
        ret += indent + "  ]\n";
        ret += indent + " }\n"; //end color
      }

      //faces
      ret += indent + " coordIndex [\n";
      for (var i = 0; i < this.faceidx; i += 3) {
        x = this.faceArray[i];
        y = this.faceArray[i + 1];
        z = this.faceArray[i + 2];
        ret += indent + "  " + x + ", " + y + ", " + z + ", -1,\n";
      }
      ret += indent + " ]\n"; //end faces
      ret += indent + "}\n"; //geometry
    }

    ret += oldindent + "}"; //shape
    return ret;
  }

  /**
   * @param {boolean} mesh
   * @param {boolean} reallocatemem
   */
  truncateArrayBuffers(mesh, reallocatemem) {
    mesh = mesh === true ? true : false;

    let vertexArr = this.vertexArray;
    let colorArr = this.colorArray;
    let normalArr = this.normalArray;
    let faceArr = this.faceArray;
    let lineArr = this.lineArray;
    let radiusArr = this.radiusArray;

    //subarray to avoid copying and reallocating memory
    this.vertexArray = vertexArr.subarray(0, this.vertices * 3);
    this.colorArray = colorArr.subarray(0, this.vertices * 3);

    if (mesh) {
      this.normalArray = normalArr.subarray(0, this.vertices * 3);
      this.faceArray = faceArr.subarray(0, this.faceidx);

      if (this.lineidx > 0)
        //not always set so reclaim memory
        this.lineArray = lineArr.subarray(0, this.lineidx);
      else this.lineArray = new Uint16Array(0);
    } else {
      this.normalArray = new Float32Array(0);
      this.faceArray = new Uint16Array(0);
      this.lineArray = new Uint16Array(0);
    }
    if (radiusArr) {
      this.radiusArray = radiusArr.subarray(0, this.vertices);
    }

    if (reallocatemem) {
      //actually copy smaller arrays to save memory
      if (this.normalArray)
        this.normalArray = new Float32Array(this.normalArray);
      if (this.faceArray) this.faceArray = new Uint16Array(this.faceArray);
      if (this.lineArray) this.lineArray = new Uint16Array(this.lineArray);
      if (this.vertexArray)
        this.vertexArray = new Float32Array(this.vertexArray);
      if (this.colorArray) this.colorArray = new Float32Array(this.colorArray);
      if (this.radiusArray)
        this.radiusArray = new Float32Array(this.radiusArray);
    }
    this.__inittedArrays = true;
  }
}
