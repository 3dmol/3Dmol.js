import { Matrix4, Quaternion, Vector3 } from "../math";
import { Scene } from "../Scene";

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
  geometry;
  material;

  lookAt(vector: Vector3) {
    this.matrix.lookAt(vector, this.position, this.up);
    if (this.rotationAutoUpdate) {
      if (this.useQuaternion === true)
        this.quaternion.copy(this.matrix.decompose()[1] as Quaternion);
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

  setVisible(val): void {
    //recursively set visibility
    this.visible = val;
    for (var i = 0; i < this.children.length; i++) {
      var child = this.children[i];
      child.setVisible(val);
    }
  }
}

