import type { Material } from './../materials/Material';
import { LineBasicMaterial } from '../materials/LineBasicMaterial';
import { EventDispatcher } from "./EventDispatcher";
import { Vector3 } from "../math";
import { CC, Color, Colored } from "../../colors";
const BUFFERSIZE = 65535; //limited to 16bit indices
export class GeometryGroup {
  id: number;
  vertexArray: Float32Array | null = null;
  colorArray: Float32Array | null = null;
  normalArray: Float32Array | null = null;
  radiusArray: Float32Array | null = null;
  faceArray: Uint16Array | null = null;
  lineArray: Uint16Array | null = null;
  vertices: number = 0;
  faceidx: number = 0;
  lineidx: number = 0;
  __inittedArrays = false;
  useOffset: unknown;

  constructor(id = 0) {
    this.id = id;
  }

  setColors(setcolor: (r: number, g: number, b: number) => Color | number): void {
    //apply a function that takes the vertex coordinate and returns a color
    var v = this.vertexArray;
    var c = this.colorArray;
    if (!v) throw new Error("vertex array not initialized");
    if (!c) throw new Error("color array not initialized");
    if (v.length != c.length) {
      console.log("Cannot re-color geometry group due to mismatched lengths.");
      return;
    }
    for (var i = 0; i < v.length; i += 3) {
      var col = setcolor(v[i], v[i + 1], v[i + 2]);
      if (!(col instanceof Color)) {
        col = CC.color(col);
      }
      c[i] = col.r;
      c[i + 1] = col.g;
      c[i + 2] = col.b;
    }
  }

  getNumVertices(): number {
    return this.vertices;
  }

  getVertices() {
    return this.vertexArray;
  }

  getCentroid() {
    if (!this.vertexArray) throw new Error("vertex array not initialized");
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

  //setup normals - vertex and face array must exist
  setNormals(): void {
    var faces = this.faceArray;
    var verts = this.vertexArray;
    var norms = this.normalArray;

    if (!this.vertices || !this.faceidx) return;
    if (!faces) throw new Error("face array not initialized");
    if (!verts) throw new Error("vertex array not initialized");
    if (!norms) throw new Error("normal array not initialized");

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

  /* sets line index array from face arr
  Note - assumes all faces are triangles (i.e. there will
  be an extra diagonal for four-sided faces - user should
  specify linearr for custom shape generation to show wireframe squares
  as rectangles rather than two triangles) */
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

    if (!faceArr) throw new Error("face array not initialized");

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

  vrml(indent: string, material?: Material) {
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
      material?.color?.r +
      " " +
      material?.color?.g +
      " " +
      material?.color?.b +
      "\n";
    if (material?.transparent) {
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
        x = this.vertexArray?.[offset];
        y = this.vertexArray?.[offset + 1];
        z = this.vertexArray?.[offset + 2];
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
        x = this.vertexArray?.[offset];
        y = this.vertexArray?.[offset + 1];
        z = this.vertexArray?.[offset + 2];
        ret += indent + "   " + x + " " + y + " " + z + ",\n";
      }
      ret += indent + "  ]\n";
      ret += indent + " }\n"; //end coordinate

      //normals
      ret += indent + " normal Normal {\n" + indent + "  vector [\n";
      for (let i = 0; i < this.vertices; ++i) {
        let offset = i * 3;
        x = this.normalArray?.[offset];
        y = this.normalArray?.[offset + 1];
        z = this.normalArray?.[offset + 2];
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
        x = this.faceArray?.[i];
        y = this.faceArray?.[i + 1];
        z = this.faceArray?.[i + 2];
        ret += indent + "  " + x + ", " + y + ", " + z + ", -1,\n";
      }
      ret += indent + " ]\n"; //end faces
      ret += indent + "}\n"; //geometry
    }

    ret += oldindent + "}"; //shape
    return ret;
  }

  truncateArrayBuffers(mesh = true, reallocatemem = false) {

    var vertexArr = this.vertexArray,
      colorArr = this.colorArray,
      normalArr = this.normalArray,
      faceArr = this.faceArray,
      lineArr = this.lineArray,
      radiusArr = this.radiusArray;

    //subarray to avoid copying and reallocating memory
    this.vertexArray = vertexArr?.subarray(0, this.vertices * 3) || null;
    this.colorArray = colorArr?.subarray(0, this.vertices * 3) || null;

    if (mesh) {
      this.normalArray = normalArr?.subarray(0, this.vertices * 3) || null;
      this.faceArray = faceArr?.subarray(0, this.faceidx) || null;

      if (this.lineidx > 0)
        //not always set so reclaim memory
        this.lineArray = lineArr?.subarray(0, this.lineidx) || null;
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
 
export class Geometry extends EventDispatcher {
  id: number;
  name: string = "";
  hasTangents: boolean =  false;
  dynamic: boolean = true; // the intermediate typed arrays will be deleted when set to false;
  radii: boolean;
  mesh: boolean;
  offset: boolean;
  verticesNeedUpdate: boolean = false;
  elementsNeedUpdate: boolean = false;
  normalsNeedUpdate: boolean = false;
  colorsNeedUpdate: boolean = false;
  buffersNeedUpdate: boolean = false;
  imposter: boolean = false;
  instanced: boolean = false;
  geometryGroups: GeometryGroup[] = [];
  groups: number = 0;
  sphereGeometry?: Geometry;
  drawnCaps?: any;
  
  constructor(mesh = false, radii = false, offset = false) {
    super();
    this.id = GeometryIDCount++;
    this.mesh = mesh; // Does this geometry represent a mesh (i.e. do we need Face/Line index buffers?)
    this.radii = radii;
    this.offset = offset; //offset buffer used for instancing
  }

  //Get geometry group to accomodate addVertices new vertices - create
  // new group if necessary
  updateGeoGroup(addVertices = 0): GeometryGroup {
    var retGroup =
      this.groups > 0 ? this.geometryGroups[this.groups - 1] : null;

    if (
      !retGroup ||
      retGroup.vertices + addVertices > (retGroup?.vertexArray?.length || 0) / 3
    )
      retGroup = this.addGeoGroup();

    return retGroup;
  }

  //return comma separated list of IndexedFace (or Line) sets from geometry groups
  vrml(indent: string, material?: Material): string {
    var ret = "";
    var len = this.geometryGroups.length;
    for (var g = 0; g < len; g++) {
      var geoGroup = this.geometryGroups[g];
      ret += geoGroup.vrml(indent, material) + ",\n";
    }
    return ret;
  }

  addGeoGroup() {
    var ret = new GeometryGroup(this.geometryGroups.length);
    this.geometryGroups.push(ret);
    this.groups = this.geometryGroups.length;
  
    ret.vertexArray = new Float32Array(BUFFERSIZE * 3);
    ret.colorArray = new Float32Array(BUFFERSIZE * 3);
  
    //TODO: instantiating uint arrays according to max number of vertices
    // is dangerous, since there exists the possibility that there will be
    // more face or line indices than vertex points - but so far that doesn't
    // seem to be the case for any of the renders
    if (this.mesh) {
      ret.normalArray = new Float32Array(BUFFERSIZE * 3);
      ret.faceArray = new Uint16Array(BUFFERSIZE * 6);
      ret.lineArray = new Uint16Array(BUFFERSIZE * 6);
    }
    if (this.radii) {
      ret.radiusArray = new Float32Array(BUFFERSIZE);
    }
    ret.useOffset = this.offset;
  
    return ret;
  }

  setUpNormals(...args: Parameters<GeometryGroup["setNormals"]>) {
    for (var g = 0; g < this.groups; g++) {
      var geoGroup = this.geometryGroups[g];

      geoGroup.setNormals(...args);
    }
  }

  setColors(...setcolor: Parameters<GeometryGroup["setColors"]>): void {
    var len = this.geometryGroups.length;
    for (var g = 0; g < len; g++) {
      var geoGroup = this.geometryGroups[g];
      geoGroup.setColors(...setcolor);
    }
  }


  setUpWireframe(...lineIndexArgs: Parameters<GeometryGroup["setLineIndices"]>) {
    for (var g = 0; g < this.groups; g++) {
      var geoGroup = this.geometryGroups[g];

      geoGroup.setLineIndices(...lineIndexArgs);
    }
  }

  //After vertices, colors, etc are collected in regular or typed arrays,
  //  create typed arrays from regular arrays if they don't already exist,
  initTypedArrays() {
    for (var g = 0; g < this.groups; g++) {
      var group = this.geometryGroups[g];

      if (group.__inittedArrays === true) continue;

      //do not actually reallocate smaller memory here because
      //of the performance hit - if you know your geometry is small,
      //truncate manually with the second parameter true
      group.truncateArrayBuffers(this.mesh, false);
    }
  }

  dispose() {
    this.dispatchEvent({ type: "dispose" });
  }

  get vertices (): number {
    var vertices = 0;
    for (var g = 0; g < this.groups; g++)
      vertices += this.geometryGroups[g].vertices;

    return vertices;
  }
}

export let GeometryIDCount = 0;