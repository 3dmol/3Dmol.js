// @ts-check
import { Material } from "../materials/Material";
import { EventDispatcher } from "./EventDispatcher";
import { addGroup, GeometryGroup } from "./GeometryGroup";

export var GeometryIDCount = 0;

export class Geometry extends EventDispatcher {
  __webglVertexBuffer;
  __webglColorBuffer;
  __webglInit;
  buffersNeedUpdate = false;
  colorsNeedUpdate = false;
  dynamic = true; // the intermediate typed arrays will be deleted when set to false
  elementsNeedUpdate = false;
  /**
   * @type {GeometryGroup[]}
   */
  geometryGroups = [];
  groups = 0;
  hasTangents = false;
  name = "";
  normalsNeedUpdate = false;
  verticesNeedUpdate = false;

  /**
   * @param {boolean} mesh
   * @param {boolean} radii
   * @param {boolean} offset
   */
  constructor(mesh, radii, offset) {
    super();
    this.id = GeometryIDCount++;
    this.mesh = mesh === true ? true : false; // Does this geometry represent a mesh (i.e. do we need Face/Line index buffers?)
    this.radii = radii || false;
    this.offset = offset || false; //offset buffer used for instancing
  }

  get vertices() {
    var vertices = 0;
    for (var g = 0; g < this.groups; g++)
      vertices += this.geometryGroups[g].vertices;

    return vertices;
  }

  /**
   * @param {number} addVertices
   */
  updateGeoGroup(addVertices) {
    addVertices = addVertices || 0;

    var retGroup =
      this.groups > 0 ? this.geometryGroups[this.groups - 1] : null;

    if (
      !retGroup ||
      retGroup.vertices + addVertices > retGroup.vertexArray.length / 3
    )
      retGroup = addGroup(this);

    return retGroup;
  }

  /**
   * @param {string} indent
   * @param {Material} material
   */
  vrml(indent, material) {
    var ret = "";
    var len = this.geometryGroups.length;
    for (var g = 0; g < len; g++) {
      var geoGroup = this.geometryGroups[g];
      ret += geoGroup.vrml(indent, material) + ",\n";
    }
    return ret;
  }

  addGeoGroup() {
    return addGroup(this);
  }

  /**
   * @param {boolean} three
   */
  setUpNormals(three) {
    three = three || false;

    for (var g = 0; g < this.groups; g++) {
      var geoGroup = this.geometryGroups[g];

      // @ts-ignore
      geoGroup.setNormals(three);
    }
  }

  /**
   * @param {any} setcolor
   */
  setColors(setcolor) {
    var len = this.geometryGroups.length;
    for (var g = 0; g < len; g++) {
      var geoGroup = this.geometryGroups[g];
      geoGroup.setColors(setcolor);
    }
  }

  setUpWireframe() {
    for (var g = 0; g < this.groups; g++) {
      var geoGroup = this.geometryGroups[g];

      geoGroup.setLineIndices();
    }
  }

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
}
