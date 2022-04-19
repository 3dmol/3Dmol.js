/* eslint-disable no-param-reassign */
/* eslint-disable max-classes-per-file */
import {Object3D} from './core';
// @ts-ignore
import {LineBasicMaterial, SpriteMaterial, MeshLambertMaterial as MeshBasicMaterial} from './materials';

/*
 * $3Dmol Mesh and Line objects
 */

export const LineStrip = 0;
export const LinePieces = 1;

// Line Object
export class Line extends Object3D {
  constructor(geometry, material, type) {
    super();

    this.geometry = geometry;
    // TODO: update material and type to webgl
    this.material =
      material !== undefined ? material : new LineBasicMaterial({color: Math.random() * 0xffffff});
    this.type = type !== undefined ? type : LineStrip;
  }

  clone(object) {
    if (object === undefined) object = new Line(this.geometry, this.material, this.type);

    Object3D.prototype.clone.call(this, object);

    return object;
  }
}

// Mesh Object
export class Mesh extends Object3D {
  doubleSided;

  constructor(geometry, material) {
    super();

    this.geometry = geometry;
    this.material =
      material !== undefined
        ? material
        : new MeshBasicMaterial({color: Math.random() * 0xffffff, wireframe: true});
  }

  clone(object) {
    if (object === undefined) object = new Mesh(this.geometry, this.material);

    Object3D.prototype.clone.call(this, object);

    return object;
  }
}

// Sprite object
// @ts-ignore
export class Sprite extends Object3D {
  constructor(material) {
    super();

    this.material = material !== undefined ? material : new SpriteMaterial();

    this.rotation3d = this.rotation;
    this.rotation = 0;
  }

  updateMatrix() {
    this.matrix.setPosition(this.position);

    this.rotation3d.set(0, 0, this.rotation);
    this.matrix.setRotationFromEuler(this.rotation3d);

    if (this.scale.x !== 1 || this.scale.y !== 1) this.matrix.scale(this.scale);

    this.matrixWorldNeedsUpdate = true;
  }

  clone(object) {
    if (object === undefined) object = new Sprite(this.material);

    Object3D.prototype.clone.call(this, object);

    return object;
  }
}
