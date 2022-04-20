/* eslint-disable no-underscore-dangle */
// Intersection sphere and box shapes.
import {Vector3} from './math';

// Intersection sphere for sphere, stick render
export class Sphere {
  /** @type {{min: import('../specs').Vector3Like, max: import('../specs').Vector3Like}|undefined} */
  box;
  /**
   * 
   * @param {Vector3} [center] 
   * @param {number} [radius] 
   */
  constructor(center, radius) {
    this.center = center !== undefined ? center : new Vector3();
    this.radius = radius !== undefined ? radius : 0;
  }

  set(center, radius) {
    this.center.copy(center);
    this.radius = radius;
    return this;
  }

  copy(sphere) {
    this.center.copy(sphere.center);
    this.radius = sphere.radius;
    return this;
  }

  applyMatrix4(matrix) {
    this.center.applyMatrix4(matrix);
    this.radius *= matrix.getMaxScaleOnAxis();
    return this;
  }

  translate(offset) {
    this.center.add(offset);
    return this;
  }

  equals(sphere) {
    return sphere.center.equals(this.center) && sphere.radius === this.radius;
  }

  clone() {
    return new Sphere().copy(this);
  }
}

// Bounding cylinder for stick render
/** @constructor */
export function Cylinder(c1, c2, radius) {
  this.c1 = c1 !== undefined ? c1 : new Vector3();

  this.c2 = c2 !== undefined ? c2 : new Vector3();

  this.direction = new Vector3().subVectors(this.c2, this.c1).normalize();

  this.radius = radius !== undefined ? radius : 0;
}

const vector = new Vector3();

Cylinder.prototype = {
  constructor: Cylinder,

  copy(cylinder) {
    this.c1.copy(cylinder.c1);
    this.c2.copy(cylinder.c2);
    this.direction.copy(cylinder.direction);
    this.radius = cylinder.radius;

    return this;
  },

  lengthSq() {
    return vector.subVectors(this.c2, this.c1).lengthSq();
  },

  applyMatrix4(matrix) {
    this.direction.add(this.c1).applyMatrix4(matrix);
    this.c1.applyMatrix4(matrix);
    this.c2.applyMatrix4(matrix);
    this.direction.sub(this.c1).normalize();
    this.radius *= matrix.getMaxScaleOnAxis();

    return this;
  },
};

// plane specified by three points
/** @constructor */
export function Triangle(a, b, c) {
  this.a = a !== undefined ? a : new Vector3();

  this.b = b !== undefined ? b : new Vector3();

  this.c = c !== undefined ? c : new Vector3();
}

const v1 = new Vector3();

Triangle.prototype = {
  constructor: Triangle,

  copy(triangle) {
    this.a.copy(triangle.a);
    this.b.copy(triangle.b);
    this.c.copy(triangle.c);

    return this;
  },

  applyMatrix4(matrix) {
    this.a.applyMatrix4(matrix);
    this.b.applyMatrix4(matrix);
    this.c.applyMatrix4(matrix);

    return this;
  },

  getNormal() {
    const norm = this.a.clone();
    norm.sub(this.b);
    v1.subVectors(this.c, this.b);

    norm.cross(v1);
    norm.normalize();

    return norm;
  },
};
