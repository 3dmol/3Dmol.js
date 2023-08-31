import type { Material } from './../materials/Material';
import { Matrix4, Quaternion, Vector3 } from "../math";
import type { Geometry } from './Geometry';
import type { Fog } from '../Fog';
import { Color, ColorConstructorArg } from "../../colors";
import { Sprite } from 'WebGL/objects';

export let Object3DIDCount = 0;
// Object3D base constructor function
export class Object3D {
  id = Object3DIDCount++;
  name = "";
  parent?: Object3D;
  children: Array<Object3D> = [];
  position = new Vector3();
  rotation: Vector3 | number = new Vector3();
  matrix = new Matrix4();
  matrixWorld = new Matrix4();
  quaternion = new Quaternion();
  eulerOrder = "XYZ";
  up = new Vector3(0, 1, 0);
  scale = new Vector3(1, 1, 1);
  matrixAutoUpdate = true;
  matrixWorldNeedsUpdate = true;
  rotationAutoUpdate = true;
  useQuaternion = false;
  visible = true;
  geometry?: Geometry;
  material?: Material;

  lookAt(vector: Vector3) {
    this.matrix.lookAt(vector, this.position, this.up);
    if (this.rotationAutoUpdate) {
      if (this.useQuaternion === true)
        console.error("Unimplemented math operation.");
      else if (this.rotation instanceof Vector3)
        this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
    }
  }

  //add child object
  add<T extends Object3D>(object: T): void {
    if (object === (this as Object3D)) {
      console.error("Can't add $3Dmol.Object3D to itself");
      return;
    }

    object.parent = this;
    this.children.push(object);

    //add to the scene (i.e. follow up this instance's parents until reach the top)

    var scene = this as Object3D;

    while (scene.parent !== undefined) scene = scene.parent;

    if (scene !== undefined && scene instanceof Scene)
      scene.__addObject(object);
  }

  remove<T extends Object3D>(object: T): void {
    var index = this.children.indexOf(object);

    if (index !== -1) {
      object.parent = undefined;
      this.children.splice(index, 1);

      //Remove from scene

      var scene = this as Object3D;

      while (scene.parent !== undefined) scene = scene.parent;

      if (scene !== undefined && scene instanceof Scene)
        scene.__removeObject(object);
    }
  }

  //convert to vrml
  vrml(indent?: string) {
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

    if (this.geometry) {
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

    if (this.useQuaternion === false && this.rotation instanceof Vector3) {
      this.matrix.setRotationFromEuler(this.rotation, this.eulerOrder);
    } else {
      this.matrix.setRotationFromQuaternion(this.quaternion);
    }

    //TODO: Do I need this??
    if (this.scale.x !== 1 || this.scale.y !== 1 || this.scale.z !== 1)
      this.matrix.scale(this.scale);

    this.matrixWorldNeedsUpdate = true;
  }

  updateMatrixWorld(force?: boolean) {
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

  clone(object?: Object3D): Object3D {
    if (object === undefined) object = new Object3D();

    object.name = this.name;

    object.up.copy(this.up);
    object.position.copy(this.position);
    if (
      object.rotation instanceof Vector3 &&
      this.rotation instanceof Vector3
    ) {
      object.rotation.copy(this.rotation);
    } else {
      object.rotation = this.rotation;
    }
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

  setVisible(val: boolean): void {
    //recursively set visibility
    this.visible = val;
    for (var i = 0; i < this.children.length; i++) {
      var child = this.children[i];
      child.setVisible(val);
    }
  }
}


/*
 * Scene class
 */
/* @constructor */
export class Scene extends Object3D {
  fog: Fog | null = null;
  //may not need...
  overrideMaterial: Material | null = null;
  matrixAutoUpdate = false;
  __objects = [] as Object3D[];
  __lights = [] as Light[];
  __objectsAdded = [] as Object3D[];
  __objectsRemoved = [] as Object3D[];
  __webglSprites: Sprite[];

  __addObject<T extends Object3D>(object: T) {
    //Directional Lighting
    if (object instanceof Light) {
      if (this.__lights.indexOf(object as unknown as Light) === -1) this.__lights.push(object as unknown as Light);

      //TODO: Do I need this??
      if ((object as unknown as Light).target && (object as unknown as Light).target.parent === undefined)
        this.add((object as unknown as Light).target);
    }

    //Rotation group
    else {
      if (this.__objects.indexOf(object) === -1) {
        this.__objects.push(object);
        this.__objectsAdded.push(object);

        //Check if previously removed

        var idx = this.__objectsRemoved.indexOf(object);

        if (idx !== -1) this.__objectsRemoved.splice(idx, 1);
      }
    }

    //Add object's children

    for (var i = 0; i < object.children.length; i++)
      this.__addObject(object.children[i]);
  }

  __removeObject<T extends Object3D>(object: T) {
    var idx;
    if (object instanceof Light) {
      idx = this.__lights.indexOf(object as unknown as Light);

      if (idx !== -1) this.__lights.splice(idx, 1);
    }

    //Object3D
    else {
      idx = this.__objects.indexOf(object);

      if (idx !== -1) {
        this.__objects.splice(idx, 1);
        this.__objectsRemoved.push(object);

        //Check if previously added

        var ai = this.__objectsAdded.indexOf(object);

        if (ai !== -1) this.__objectsAdded.splice(idx, 1);
      }
    }

    //Remove object's children
    for (var i = 0; i < object.children.length; i++)
      this.__removeObject(object.children[i]);
  }
}


export class Light extends Object3D {
  color: Color;
  intensity: any;
  position = new Vector3(0, 1, 0);
  target = new Object3D();
  castShadow = false;
  onlyShadow = false;
  constructor(hex?: ColorConstructorArg, intensity: number = 1) {
    super();
    this.color = new Color(hex);
    this.intensity = intensity;
  }
}
