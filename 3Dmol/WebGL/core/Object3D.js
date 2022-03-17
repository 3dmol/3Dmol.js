// @ts-check

import { Geometry } from "./Geometry";
import { Material } from "../materials/Material";
import { Matrix4 } from "../math/matrix4";
import { Quaternion } from "../math/Quaternion";
import { Vector3 } from "../math/Vector3";
import { Scene } from "./Scene";

export var Object3DIDCount = 0;

//Object3D base constructor function
/** @this {Object3D} */
export class Object3D {
  /**
   * @type {Geometry | undefined}
   */
  geometry;
  /**
   * @type {Material | undefined}
   */
  material;
  constructor() {
    this.id = Object3DIDCount++;

    this.name = "";

    /**
     * @type {Object3D | undefined}
     */
    this.parent = undefined;
    /**
     * @type {Object3D[]}
     */
    this.children = [];

    this.position = new Vector3();
    this.rotation = new Vector3();
    this.matrix = new Matrix4();
    this.matrixWorld = new Matrix4();
    this.quaternion = new Quaternion();
    this.eulerOrder = "XYZ";

    this.up = new Vector3(0, 1, 0);
    this.scale = new Vector3(1, 1, 1);

    this.matrixAutoUpdate = true;
    this.matrixWorldNeedsUpdate = true;
    this.rotationAutoUpdate = true;
    this.useQuaternion = false;

    this.visible = true;
  }

  /**
   * @param {any} vector
   */
  lookAt(vector) {
    this.matrix.lookAt(vector, this.position, this.up);

    if (this.rotationAutoUpdate) {
      if (this.useQuaternion === true)
        this.quaternion.copy(this.matrix.decompose()[1]);
      else
        this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
    }
  }

  /**
   * @param {Object3D} object
   */
  add(object) {
    if (object === this) {
      console.error("Can't add Object3D to itself");
      return;
    }

    object.parent = this;
    this.children.push(object);

    //add to the scene (i.e. follow up this instance's parents until reach the top)

    /**
     * @type {Object3D}
     */
    var scene = this;

    while (scene.parent !== undefined) scene = scene.parent;

    if (scene !== undefined && scene instanceof Scene)
      scene.__addObject(object);
  }

  /**
   * @param {Object3D} object
   */
  remove(object) {
    var index = this.children.indexOf(object);

    if (index !== -1) {
      object.parent = undefined;
      this.children.splice(index, 1);

      //Remove from scene
      /**
       * @type {Object3D}
       */
      var scene = this;

      while (scene.parent !== undefined) scene = scene.parent;

      if (scene !== undefined && scene instanceof Scene)
        scene.__removeObject(object);
    }
  }

  /**
   * @param {string} indent
   */
  vrml(indent) {
    //attempt to pretty print
    if (!indent) indent = " ";
    //all objects have a transformation (usually identity)
    //not quite sure if getting rotation right here..
    var theta = 2 * Math.atan2(this.quaternion.lengthxyz(), this.quaternion.w);
    var x = 0,
      y = 0,
      z = 0;
    if (theta != 0) {
      let st = Math.sin(theta / 2);
      x = this.quaternion.x / st;
      y = this.quaternion.y / st;
      z = this.quaternion.z / st;
    }
    var ret =
      indent +
      "Transform {\n" +
      indent +
      " center " +
      this.position.x +
      " " +
      this.position.y +
      " " +
      this.position.z +
      "\n" +
      indent +
      " rotation " +
      x +
      " " +
      y +
      " " +
      z +
      " " +
      theta +
      "\n" +
      indent +
      " children [\n";

    // added a material null check here
    if (this.geometry && this.material) {
      ret += this.geometry.vrml(indent, this.material);
    }
    for (var i = 0; i < this.children.length; i++) {
      ret += this.children[i].vrml(indent + " ") + ",\n";
    }
    ret += " ]\n";
    ret += "}";
    return ret;
  }

  updateMatrix() {
    this.matrix.setPosition(this.position);

    if (this.useQuaternion === false)
      this.matrix.setRotationFromEuler(this.rotation, this.eulerOrder);
    else this.matrix.setRotationFromQuaternion(this.quaternion);

    //TODO: Do I need this??
    if (this.scale.x !== 1 || this.scale.y !== 1 || this.scale.z !== 1)
      this.matrix.scale(this.scale);

    this.matrixWorldNeedsUpdate = true;
  }

  /**
   * @param {boolean} [force]
   */
  updateMatrixWorld(force) {
    if (this.matrixAutoUpdate === true) this.updateMatrix();

    if (this.matrixWorldNeedsUpdate === true || force === true) {
      if (this.parent === undefined) this.matrixWorld.copy(this.matrix);
      else
        this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
    }

    this.matrixWorldNeedsUpdate = false;

    //Update matrices of all children
    for (var i = 0; i < this.children.length; i++) {
      this.children[i].updateMatrixWorld(true);
    }
  }

  /**
   * @param {Object3D} [object]
   */
  clone(object) {
    if (object === undefined) object = new Object3D();

    object.name = this.name;

    object.up.copy(this.up);
    object.position.copy(this.position);
    object.rotation.copy(this.rotation);
    object.eulerOrder = this.eulerOrder;
    object.scale.copy(this.scale);

    object.rotationAutoUpdate = this.rotationAutoUpdate;
    object.matrix.copy(this.matrix);
    object.matrixWorld.copy(this.matrixWorld);
    object.quaternion.copy(this.quaternion);
    object.matrixAutoUpdate = this.matrixAutoUpdate;
    object.matrixWorldNeedsUpdate = this.matrixWorldNeedsUpdate;

    object.useQuaternion = this.useQuaternion;

    object.visible = this.visible;

    for (var i = 0; i < this.children.length; i++) {
      var child = this.children[i];
      object.add(child.clone());
    }

    return object;
  }

  /**
   * @param {boolean} val
   */
  setVisible(val) {
    //recursively set visibility
    this.visible = val;
    for (var i = 0; i < this.children.length; i++) {
      var child = this.children[i];
      child.setVisible(val);
    }
  }
}
