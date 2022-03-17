// @ts-check
/*
 * Scene class
 */

import { Object3D } from "./Object3D";
import { Light } from "../objects/Light";

export class Scene extends Object3D {
  constructor() {
    super();

    this.fog = null;

    //May not need...
    this.overrideMaterial = null;

    this.matrixAutoUpdate = false;

    /**
     * @type {any[]}
     */
    this.__objects = [];
    /**
     * @type {Light[]}
     */
    this.__lights = [];

    /**
     * @type {any[]}
     */
    this.__objectsAdded = [];
    /**
     * @type {any[]}
     */
    this.__objectsRemoved = [];
    this.__webglSprites = undefined;
    this.__webglObjects = undefined;
    this.__webglObjectsImmediate = undefined;
    this.__webglFlares = undefined;
  }

  /**
   * @param {Object3D} object
   */
  __addObject(object) {
    //Directional Lighting
    if (object instanceof Light) {
      if (this.__lights.indexOf(object) === -1) this.__lights.push(object);

      //TODO: Do I need this??
      if (object.target && object.target.parent === undefined)
        this.add(object.target);
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

  /**
   * @param {Object3D} object
   */
  __removeObject(object) {
    var idx;
    if (object instanceof Light) {
      idx = this.__lights.indexOf(object);

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
