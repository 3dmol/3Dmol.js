import { clamp } from ".";

// A 3 Vector
export class Vector3 {
  constructor(x, y, z) {
    this.x = x || 0.0;
    this.y = y || 0.0;
    this.z = z || 0.0;
  }

  set(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;

    return this;
  }

  copy(v) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;

    return this;
  }

  add(v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;

    return this;
  }

  addVectors(a, b) {
    this.x = a.x + b.x;
    this.y = a.y + b.y;
    this.z = a.z + b.z;

    return this;
  }

  multiplyVectors(a, b) {
    //elementwise
    this.x = a.x * b.x;
    this.y = a.y * b.y;
    this.z = a.z * b.z;

    return this;
  }

  sub(v) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;

    return this;
  }

  subVectors(a, b) {
    this.x = a.x - b.x;
    this.y = a.y - b.y;
    this.z = a.z - b.z;

    return this;
  }

  multiplyScalar(s) {
    this.x *= s;
    this.y *= s;
    this.z *= s;

    return this;
  }

  divideScalar(s) {
    if (s !== 0) {
      this.x /= s;
      this.y /= s;
      this.z /= s;
    } else {
      this.x = 0;
      this.y = 0;
      this.z = 0;
    }

    return this;
  }

  max(s) {
    this.x = Math.max(this.x, s.x);
    this.y = Math.max(this.y, s.y);
    this.z = Math.max(this.z, s.z);

    return this;
  }

  min(s) {
    this.x = Math.min(this.x, s.x);
    this.y = Math.min(this.y, s.y);
    this.z = Math.min(this.z, s.z);

    return this;
  }

  distanceTo(v) {
    return Math.sqrt(this.distanceToSquared(v));
  }

  distanceToSquared(v) {
    var dx = this.x - v.x;
    var dy = this.y - v.y;
    var dz = this.z - v.z;

    return dx * dx + dy * dy + dz * dz;
  }

  applyMatrix3(m) {
    var x = this.x,
      y = this.y,
      z = this.z;

    var e = m.elements;
    //column major ordering
    this.x = e[0] * x + e[3] * y + e[6] * z;
    this.y = e[1] * x + e[4] * y + e[7] * z;
    this.z = e[2] * x + e[5] * y + e[8] * z;

    return this;
  }

  applyMatrix4(m) {
    var x = this.x,
      y = this.y,
      z = this.z;

    var e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
    this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
    this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

    return this;
  }

  applyProjection(m) {
    // input: Matrix4 projection matrix

    var x = this.x,
      y = this.y,
      z = this.z;

    var e = m.elements;
    var d = e[3] * x + e[7] * y + e[11] * z + e[15];

    this.x = (e[0] * x + e[4] * y + e[8] * z + e[12]) / d;
    this.y = (e[1] * x + e[5] * y + e[9] * z + e[13]) / d;
    this.z = (e[2] * x + e[6] * y + e[10] * z + e[14]) / d;

    return this;
  }

  applyQuaternion(q) {
    // save values
    var x = this.x;
    var y = this.y;
    var z = this.z;

    var qx = q.x;
    var qy = q.y;
    var qz = q.z;
    var qw = q.w;

    // compute this as
    // t = 2 * cross(q.xyz, v)
    // newv = v + q.w * t + cross(q.xyz, t)
    // this from molecularmusings
    // http://molecularmusings.wordpress.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    var t = {};
    t.x = 2 * (y * qz - z * qy);
    t.y = 2 * (z * qx - x * qz);
    t.z = 2 * (x * qy - y * qx);

    // cross t with q
    var t2 = {};
    t2.x = t.y * qz - t.z * qy;
    t2.y = t.z * qx - t.x * qz;
    t2.z = t.x * qy - t.y * qx;

    this.x = x + qw * t.x + t2.x;
    this.y = y + qw * t.y + t2.y;
    this.z = z + qw * t.z + t2.z;

    return this;
  }

  negate() {
    return this.multiplyScalar(-1);
  }

  dot(v) {
    return this.x * v.x + this.y * v.y + this.z * v.z;
  }

  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }

  lengthSq() {
    return this.x * this.x + this.y * this.y + this.z * this.z;
  }

  normalize() {
    return this.divideScalar(this.length());
  }

  cross(v) {
    var x = this.x,
      y = this.y,
      z = this.z;

    this.x = y * v.z - z * v.y;
    this.y = z * v.x - x * v.z;
    this.z = x * v.y - y * v.x;

    return this;
  }

  crossVectors(a, b) {
    this.x = a.y * b.z - a.z * b.y;
    this.y = a.z * b.x - a.x * b.z;
    this.z = a.x * b.y - a.y * b.x;

    return this;
  }

  getPositionFromMatrix(m) {
    this.x = m.elements[12];
    this.y = m.elements[13];
    this.z = m.elements[14];

    return this;
  }

  setEulerFromRotationMatrix(m, order) {
    // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

    var te = m.elements;
    var m11 = te[0],
      m12 = te[4],
      m13 = te[8];
    // var m21 = te[1];
    var m22 = te[5],
      m23 = te[9];
    // var m31 = te[2];
    var m32 = te[6],
      m33 = te[10];

    if (order === undefined || order === "XYZ") {
      this.y = Math.asin(clamp(m13, -1, 1));

      if (Math.abs(m13) < 0.99999) {
        this.x = Math.atan2(-m23, m33);
        this.z = Math.atan2(-m12, m11);
      } else {
        this.x = Math.atan2(m32, m22);
        this.z = 0;
      }
    } else {
      console.error(
        "Error with vector's setEulerFromRotationMatrix: Unknown order: " +
          order
      );
    }

    return this;
  }

  rotateAboutVector(axis, ang) {
    axis.normalize();
    var cosang = Math.cos(ang);
    var sinang = Math.sin(ang);
    // Rodrigues' rotation formula, from wikipedia

    var term1 = this.clone().multiplyScalar(cosang);
    var term2 = axis.clone().cross(this).multiplyScalar(sinang);
    var term3 = axis
      .clone()
      .multiplyScalar(axis.clone().dot(this))
      .multiplyScalar(1 - cosang);

    var rot = term1.add(term2).add(term3);

    this.x = rot.x;
    this.y = rot.y;
    this.z = rot.z;

    return this;
  }

  setFromMatrixPosition(m) {
    var e = m.elements;

    this.x = e[12];
    this.y = e[13];
    this.z = e[14];

    return this;
  }

  transformDirection(m) {
    // input: THREE.Matrix4 affine matrix
    // vector interpreted as a direction

    var x = this.x,
      y = this.y,
      z = this.z;
    var e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z;
    this.y = e[1] * x + e[5] * y + e[9] * z;
    this.z = e[2] * x + e[6] * y + e[10] * z;

    return this.normalize();
  }

  clone() {
    return new Vector3(this.x, this.y, this.z);
  }
}
