/* eslint-disable prefer-destructuring */
/* eslint-disable no-param-reassign */
/* eslint-disable max-classes-per-file */
/*
 * math-like functionality
 * quaternion, vector, matrix
 */
export function clamp(x, min, max) {
  return Math.min(Math.max(x, min), max);
}

const degreeToRadiansFactor = Math.PI / 180;

export function degToRad(deg) {
  return deg * degreeToRadiansFactor;
}

// Quaternion
/** @constructor */
// Quaternion
/** @constructor */
export class Quaternion {
  constructor(x, y, z, w) {
    this.x = x || 0;
    this.y = y || 0;
    this.z = z || 0;
    this.w = w !== undefined ? w : 1;
  }

  set(x, y, z, w) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;

    return this;
  }

  copy(q) {
    this.x = q.x;
    this.y = q.y;
    this.z = q.z;
    this.w = q.w;

    return this;
  }

  conjugate() {
    this.x *= -1;
    this.y *= -1;
    this.z *= -1;

    return this;
  }

  inverse() {
    return this.conjugate().normalize();
  }

  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w);
  }

  lengthxyz() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }

  normalize() {
    let l = this.length();

    if (l === 0) {
      this.x = 0;
      this.y = 0;
      this.z = 0;
      this.w = 1;
    } else {
      l = 1 / l;

      this.x *= l;
      this.y *= l;
      this.z *= l;
      this.w *= l;
    }

    return this;
  }

  multiply(q) {
    return this.multiplyQuaternions(this, q);
  }

  multiplyScalar(s) {
    this.x *= s;
    this.y *= s;
    this.z *= s;
    this.w *= s;
    return this;
  }

  multiplyQuaternions(a, b) {
    const qax = a.x;
    const qay = a.y;
    const qaz = a.z;
    const qaw = a.w;
    const qbx = b.x;
    const qby = b.y;
    const qbz = b.z;
    const qbw = b.w;

    this.x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
    this.y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
    this.z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
    this.w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;
    return this;
  }

  sub(q) {
    this.x -= q.x;
    this.y -= q.y;
    this.z -= q.z;
    this.w -= q.w;
    return this;
  }

  clone() {
    return new Quaternion(this.x, this.y, this.z, this.w);
  }

  setFromEuler(e) {
    const c1 = Math.cos(e.x / 2);
    const c2 = Math.cos(e.y / 2);
    const c3 = Math.cos(e.z / 2);
    const s1 = Math.sin(e.x / 2);
    const s2 = Math.sin(e.y / 2);
    const s3 = Math.sin(e.z / 2);

    this.x = s1 * c2 * c3 + c1 * s2 * s3;
    this.y = c1 * s2 * c3 - s1 * c2 * s3;
    this.z = c1 * c2 * s3 + s1 * s2 * c3;
    this.w = c1 * c2 * c3 - s1 * s2 * s3;

    return this;
  }
}

// A 2 Vector
export class Vector2 {
  constructor(x, y) {
    this.x = x || 0.0;
    this.y = y || 0.0;
  }

  set(x, y) {
    this.x = x;
    this.y = y;

    return this;
  }

  subVectors(a, b) {
    this.x = a.x - b.x;
    this.y = a.y - b.y;

    return this;
  }

  copy(v) {
    this.x = v.x;
    this.y = v.y;

    return this;
  }

  clone() {
    return new Vector2(this.x, this.y);
  }
}

// A 3 Vector
export class Vector3 {
  // unaccounted for assignents to vector3 properties in other parts of the code
  // look in glcartoon.js for example
  color;
  resi;
  style;
  smoothen;
  atom;
  skip;


  /** @type {number} */
  x;
  /** @type {number} */
  y;
  /** @type {number} */
  z;
  
  constructor(x, y, z) {
    this.x = x || 0.0;
    this.y = y || 0.0;
    this.z = z || 0.0;
    this.atomid = undefined;
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
    // elementwise
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

  // accumulate maximum
  max(s) {
    this.x = Math.max(this.x, s.x);
    this.y = Math.max(this.y, s.y);
    this.z = Math.max(this.z, s.z);

    return this;
  }

  // accumulate min
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
    const dx = this.x - v.x;
    const dy = this.y - v.y;
    const dz = this.z - v.z;

    return dx * dx + dy * dy + dz * dz;
  }

  applyMatrix3(m) {
    const {x} = this;
    const {y} = this;
    const {z} = this;

    const e = m.elements;
    // column major ordering
    this.x = e[0] * x + e[3] * y + e[6] * z;
    this.y = e[1] * x + e[4] * y + e[7] * z;
    this.z = e[2] * x + e[5] * y + e[8] * z;

    return this;
  }

  applyMatrix4(m) {
    const {x} = this;
    const {y} = this;
    const {z} = this;

    const e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
    this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
    this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

    return this;
  }

  applyProjection(m) {
    // input: Matrix4 projection matrix
    const {x} = this;
    const {y} = this;
    const {z} = this;

    const e = m.elements;
    const d = e[3] * x + e[7] * y + e[11] * z + e[15];

    this.x = (e[0] * x + e[4] * y + e[8] * z + e[12]) / d;
    this.y = (e[1] * x + e[5] * y + e[9] * z + e[13]) / d;
    this.z = (e[2] * x + e[6] * y + e[10] * z + e[14]) / d;

    return this;
  }

  applyQuaternion(q) {
    // save values
    const {x} = this;
    const {y} = this;
    const {z} = this;

    const qx = q.x;
    const qy = q.y;
    const qz = q.z;
    const qw = q.w;

    // compute this as
    // t = 2 * cross(q.xyz, v)
    // newv = v + q.w * t + cross(q.xyz, t)
    // this from molecularmusings
    // http://molecularmusings.wordpress.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    const t = {};
    t.x = 2 * (y * qz - z * qy);
    t.y = 2 * (z * qx - x * qz);
    t.z = 2 * (x * qy - y * qx);

    // cross t with q
    const t2 = {};
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
    const {x} = this;
    const {y} = this;
    const {z} = this;

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
    const te = m.elements;
    const m11 = te[0];
    const m12 = te[4];
    const m13 = te[8];
    // var m21 = te[1];
    const m22 = te[5];
    const m23 = te[9];
    // var m31 = te[2];
    const m32 = te[6];
    const m33 = te[10];

    if (order === undefined || order === 'XYZ') {
      this.y = Math.asin(clamp(m13, -1, 1));

      if (Math.abs(m13) < 0.99999) {
        this.x = Math.atan2(-m23, m33);
        this.z = Math.atan2(-m12, m11);
      } else {
        this.x = Math.atan2(m32, m22);
        this.z = 0;
      }
    } else {
      console.error(`Error with vector's setEulerFromRotationMatrix: Unknown order: ${order}`);
    }

    return this;
  }

  rotateAboutVector(axis, ang) {
    axis.normalize();
    const cosang = Math.cos(ang);
    const sinang = Math.sin(ang);
    // Rodrigues' rotation formula, from wikipedia
    const term1 = this.clone().multiplyScalar(cosang);
    const term2 = axis.clone().cross(this).multiplyScalar(sinang);
    const term3 = axis
      .clone()
      .multiplyScalar(axis.clone().dot(this))
      .multiplyScalar(1 - cosang);

    const rot = term1.add(term2).add(term3);

    this.x = rot.x;
    this.y = rot.y;
    this.z = rot.z;

    return this;
  }

  setFromMatrixPosition(m) {
    const e = m.elements;

    this.x = e[12];
    this.y = e[13];
    this.z = e[14];

    return this;
  }

  // unproject is defined after Matrix4
  transformDirection(m) {
    // input: THREE.Matrix4 affine matrix
    // vector interpreted as a direction
    const {x} = this;
    const {y} = this;
    const {z} = this;
    const e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z;
    this.y = e[1] * x + e[5] * y + e[9] * z;
    this.z = e[2] * x + e[6] * y + e[10] * z;

    return this.normalize();
  }

  clone() {
    return new Vector3(this.x, this.y, this.z);
  }
}

// Matrices

// Matrix3
/** @constructor */
// Matrices
// Matrix3
/** @constructor */
export class Matrix3 {
  constructor(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
    this.elements = new Float32Array(9);

    this.set(
      n11 !== undefined ? n11 : 1,
      n12 || 0,
      n13 || 0,
      n21 || 0,
      n22 !== undefined ? n22 : 1,
      n23 || 0,
      n31 || 0,
      n32 || 0,
      n33 !== undefined ? n33 : 1
    );
  }

  set(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
    const te = this.elements;

    te[0] = n11;
    te[3] = n12;
    te[6] = n13;
    te[1] = n21;
    te[4] = n22;
    te[7] = n23;
    te[2] = n31;
    te[5] = n32;
    te[8] = n33;

    return this;
  }

  identity() {
    this.set(1, 0, 0, 0, 1, 0, 0, 0, 1);

    return this;
  }

  copy(m) {
    const me = m.elements;

    this.set(me[0], me[3], me[6], me[1], me[4], me[7], me[2], me[5], me[8]);
  }

  multiplyScalar(s) {
    const te = this.elements;

    te[0] *= s;
    te[3] *= s;
    te[6] *= s;
    te[1] *= s;
    te[4] *= s;
    te[7] *= s;
    te[2] *= s;
    te[5] *= s;
    te[8] *= s;

    return this;
  }

  getInverse3(matrix) {
    // input: Matrix3
    const me = matrix.elements;
    const te = this.elements;

    te[0] = me[4] * me[8] - me[5] * me[7];
    te[3] = me[6] * me[5] - me[3] * me[8];
    te[6] = me[3] * me[7] - me[6] * me[4];
    te[1] = me[7] * me[2] - me[1] * me[8];
    te[4] = me[0] * me[8] - me[6] * me[2];
    te[7] = me[1] * me[6] - me[0] * me[7];
    te[2] = me[1] * me[5] - me[2] * me[4];
    te[5] = me[2] * me[3] - me[0] * me[5];
    te[8] = me[0] * me[4] - me[1] * me[3];

    const det = me[0] * te[0] + me[3] * te[1] + me[6] * te[2];
    this.multiplyScalar(1.0 / det);

    return this;
  }

  getInverse(matrix, throwOnInvertible) {
    // input: Matrix4
    const me = matrix.elements;
    const te = this.elements;

    te[0] = me[10] * me[5] - me[6] * me[9];
    te[1] = -me[10] * me[1] + me[2] * me[9];
    te[2] = me[6] * me[1] - me[2] * me[5];
    te[3] = -me[10] * me[4] + me[6] * me[8];
    te[4] = me[10] * me[0] - me[2] * me[8];
    te[5] = -me[6] * me[0] + me[2] * me[4];
    te[6] = me[9] * me[4] - me[5] * me[8];
    te[7] = -me[9] * me[0] + me[1] * me[8];
    te[8] = me[5] * me[0] - me[1] * me[4];

    const det = me[0] * te[0] + me[1] * te[3] + me[2] * te[6];

    // no inverse
    if (det === 0) {
      const msg = "Matrix3.getInverse(): can't invert matrix, determinant is 0";

      if (throwOnInvertible || false) {
        throw new Error(msg);
      } else {
        console.warn(msg);
      }

      this.identity();

      return this;
    }

    this.multiplyScalar(1.0 / det);

    return this;
  }

  // https://en.wikipedia.org/wiki/Determinant
  getDeterminant() {
    const m = this.elements;

    /*
     * |a b c| |d e f| |g h i|
     */
    const determinant =
      m[0] * m[4] * m[8] + // +aei
      m[1] * m[5] * m[6] + // +bfg
      m[2] * m[3] * m[7] - // +cdh
      m[2] * m[4] * m[6] - // -ceg
      m[1] * m[3] * m[8] - // -bdi
      m[0] * m[5] * m[7]; // -afh
    return determinant;
  }

  transpose() {
    let tmp;
    const m = this.elements;

    tmp = m[1];
    m[1] = m[3];
    m[3] = tmp;
    tmp = m[2];
    m[2] = m[6];
    m[6] = tmp;
    tmp = m[5];
    m[5] = m[7];
    m[7] = tmp;

    return this;
  }

  clone() {
    const te = this.elements;

    return new Matrix3(te[0], te[3], te[6], te[1], te[4], te[7], te[2], te[5], te[8]);
  }
}

export function square(n) {
  return n * n;
}

// return conversion matrix given crystal unit cell parameters
export function conversionMatrix3(a, b, c, alpha, beta, gamma) {
  // convert to radians
  alpha = (alpha * Math.PI) / 180;
  beta = (beta * Math.PI) / 180;
  gamma = (gamma * Math.PI) / 180;
  const sqr = square;
  const cosAlpha = Math.cos(alpha);
  const cosBeta = Math.cos(beta);
  const cosGamma = Math.cos(gamma);
  const sinGamma = Math.sin(gamma);
  const conversionMatrix = new Matrix3(
    a,
    b * cosGamma,
    c * cosBeta,
    0,
    b * sinGamma,
    (c * (cosAlpha - cosBeta * cosGamma)) / sinGamma,
    0,
    0,
    (c *
      Math.sqrt(
        1 - sqr(cosAlpha) - sqr(cosBeta) - sqr(cosGamma) + 2 * cosAlpha * cosBeta * cosGamma
      )) /
      sinGamma
  );
  return conversionMatrix;
}

// Matrix 4
export class Matrix4 {
  constructor(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44) {
    if (typeof n12 === 'undefined' && typeof n11 !== 'undefined') {
      // passing list like initialization
      this.elements = new Float32Array(n11);
    } else {
      this.elements = new Float32Array(16);

      this.elements[0] = n11 !== undefined ? n11 : 1;
      this.elements[4] = n12 || 0;
      this.elements[8] = n13 || 0;
      this.elements[12] = n14 || 0;
      this.elements[1] = n21 || 0;
      this.elements[5] = n22 !== undefined ? n22 : 1;
      this.elements[9] = n23 || 0;
      this.elements[13] = n24 || 0;
      this.elements[2] = n31 || 0;
      this.elements[6] = n32 || 0;
      this.elements[10] = n33 !== undefined ? n33 : 1;
      this.elements[14] = n34 || 0;
      this.elements[3] = n41 || 0;
      this.elements[7] = n42 || 0;
      this.elements[11] = n43 || 0;
      this.elements[15] = n44 !== undefined ? n44 : 1;
    }
  }

  // eslint-disable-next-line no-unused-vars, class-methods-use-this
  makeScale(x, y, z) {
    throw new Error("Method not implemented.");
  }

  set(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44) {
    const te = this.elements;

    te[0] = n11;
    te[4] = n12;
    te[8] = n13;
    te[12] = n14;
    te[1] = n21;
    te[5] = n22;
    te[9] = n23;
    te[13] = n24;
    te[2] = n31;
    te[6] = n32;
    te[10] = n33;
    te[14] = n34;
    te[3] = n41;
    te[7] = n42;
    te[11] = n43;
    te[15] = n44;

    return this;
  }

  identity() {
    this.set(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);

    return this;
  }

  copy(m) {
    const me = m.elements;

    this.set(
      me[0],
      me[4],
      me[8],
      me[12],
      me[1],
      me[5],
      me[9],
      me[13],
      me[2],
      me[6],
      me[10],
      me[14],
      me[3],
      me[7],
      me[11],
      me[15]
    );

    return this;
  }

  matrix3FromTopLeft() {
    const te = this.elements;
    return new Matrix3(te[0], te[4], te[8], te[1], te[5], te[9], te[2], te[6], te[10]);
  }

  setRotationFromEuler(v, order) {
    const te = this.elements;

    const {x} = v;
    const {y} = v;
    const {z} = v;
    const a = Math.cos(x);
    const b = Math.sin(x);
    const c = Math.cos(y);
    const d = Math.sin(y);
    const e = Math.cos(z);
    const f = Math.sin(z);

    if (order === undefined || order === 'XYZ') {
      const ae = a * e;
      const af = a * f;
      const be = b * e;
      const bf = b * f;

      te[0] = c * e;
      te[4] = -c * f;
      te[8] = d;

      te[1] = af + be * d;
      te[5] = ae - bf * d;
      te[9] = -b * c;

      te[2] = bf - ae * d;
      te[6] = be + af * d;
      te[10] = a * c;
    } else console.error(`Error with matrix4 setRotationFromEuler. Order: ${order}`);

    return this;
  }

  setRotationFromQuaternion(q) {
    const te = this.elements;

    const {x} = q;
    const {y} = q;
    const {z} = q;
    const {w} = q;
    const x2 = x + x;
    const y2 = y + y;
    const z2 = z + z;
    const xx = x * x2;
    const xy = x * y2;
    const xz = x * z2;
    const yy = y * y2;
    const yz = y * z2;
    const zz = z * z2;
    const wx = w * x2;
    const wy = w * y2;
    const wz = w * z2;

    te[0] = 1 - (yy + zz);
    te[4] = xy - wz;
    te[8] = xz + wy;

    te[1] = xy + wz;
    te[5] = 1 - (xx + zz);
    te[9] = yz - wx;

    te[2] = xz - wy;
    te[6] = yz + wx;
    te[10] = 1 - (xx + yy);

    return this;
  }

  multiplyMatrices(a, b) {
    const ae = a.elements;
    const be = b.elements;
    const te = this.elements;

    const a11 = ae[0];
    const a12 = ae[4];
    const a13 = ae[8];
    const a14 = ae[12];
    const a21 = ae[1];
    const a22 = ae[5];
    const a23 = ae[9];
    const a24 = ae[13];
    const a31 = ae[2];
    const a32 = ae[6];
    const a33 = ae[10];
    const a34 = ae[14];
    const a41 = ae[3];
    const a42 = ae[7];
    const a43 = ae[11];
    const a44 = ae[15];

    const b11 = be[0];
    const b12 = be[4];
    const b13 = be[8];
    const b14 = be[12];
    const b21 = be[1];
    const b22 = be[5];
    const b23 = be[9];
    const b24 = be[13];
    const b31 = be[2];
    const b32 = be[6];
    const b33 = be[10];
    const b34 = be[14];
    const b41 = be[3];
    const b42 = be[7];
    const b43 = be[11];
    const b44 = be[15];

    te[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
    te[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
    te[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
    te[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

    te[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
    te[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
    te[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
    te[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

    te[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
    te[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
    te[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
    te[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

    te[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
    te[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
    te[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
    te[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

    return this;
  }

  multiplyScalar(s) {
    const te = this.elements;

    te[0] *= s;
    te[4] *= s;
    te[8] *= s;
    te[12] *= s;
    te[1] *= s;
    te[5] *= s;
    te[9] *= s;
    te[13] *= s;
    te[2] *= s;
    te[6] *= s;
    te[10] *= s;
    te[14] *= s;
    te[3] *= s;
    te[7] *= s;
    te[11] *= s;
    te[15] *= s;

    return this;
  }

  makeTranslation(x, y, z) {
    this.set(1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1);

    return this;
  }

  // snap values close to integers to their integer value
  // useful and identifying identity matrices
  snap(digits) {
    if (!digits) digits = 4;
    const mult = 10 ** 4;
    const te = this.elements;
    for (let i = 0; i < 16; i++) {
      const rounded = Math.round(te[i]);
      if (rounded === Math.round(te[i] * mult) / mult) {
        te[i] = rounded;
      }
    }
    return this;
  }

  transpose() {
    const te = this.elements;
    let tmp;

    tmp = te[1];
    te[1] = te[4];
    te[4] = tmp;
    tmp = te[2];
    te[2] = te[8];
    te[8] = tmp;
    tmp = te[6];
    te[6] = te[9];
    te[9] = tmp;

    tmp = te[3];
    te[3] = te[12];
    te[12] = tmp;
    tmp = te[7];
    te[7] = te[13];
    te[13] = tmp;
    tmp = te[11];
    te[11] = te[14];
    te[14] = tmp;

    return this;
  }

  setPosition(v) {
    const te = this.elements;

    te[12] = v.x;
    te[13] = v.y;
    te[14] = v.z;

    return this;
  }

  translate(v) {
    const te = this.elements;

    te[12] += v.x;
    te[13] += v.y;
    te[14] += v.z;

    return this;
  }

  getInverse(m, throwOnInvertible) {
    // based on
    // http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
    const te = this.elements;
    const me = m.elements;

    const n11 = me[0];
    const n12 = me[4];
    const n13 = me[8];
    const n14 = me[12];
    const n21 = me[1];
    const n22 = me[5];
    const n23 = me[9];
    const n24 = me[13];
    const n31 = me[2];
    const n32 = me[6];
    const n33 = me[10];
    const n34 = me[14];
    const n41 = me[3];
    const n42 = me[7];
    const n43 = me[11];
    const n44 = me[15];

    te[0] =
      n23 * n34 * n42 -
      n24 * n33 * n42 +
      n24 * n32 * n43 -
      n22 * n34 * n43 -
      n23 * n32 * n44 +
      n22 * n33 * n44;
    te[4] =
      n14 * n33 * n42 -
      n13 * n34 * n42 -
      n14 * n32 * n43 +
      n12 * n34 * n43 +
      n13 * n32 * n44 -
      n12 * n33 * n44;
    te[8] =
      n13 * n24 * n42 -
      n14 * n23 * n42 +
      n14 * n22 * n43 -
      n12 * n24 * n43 -
      n13 * n22 * n44 +
      n12 * n23 * n44;
    te[12] =
      n14 * n23 * n32 -
      n13 * n24 * n32 -
      n14 * n22 * n33 +
      n12 * n24 * n33 +
      n13 * n22 * n34 -
      n12 * n23 * n34;
    te[1] =
      n24 * n33 * n41 -
      n23 * n34 * n41 -
      n24 * n31 * n43 +
      n21 * n34 * n43 +
      n23 * n31 * n44 -
      n21 * n33 * n44;
    te[5] =
      n13 * n34 * n41 -
      n14 * n33 * n41 +
      n14 * n31 * n43 -
      n11 * n34 * n43 -
      n13 * n31 * n44 +
      n11 * n33 * n44;
    te[9] =
      n14 * n23 * n41 -
      n13 * n24 * n41 -
      n14 * n21 * n43 +
      n11 * n24 * n43 +
      n13 * n21 * n44 -
      n11 * n23 * n44;
    te[13] =
      n13 * n24 * n31 -
      n14 * n23 * n31 +
      n14 * n21 * n33 -
      n11 * n24 * n33 -
      n13 * n21 * n34 +
      n11 * n23 * n34;
    te[2] =
      n22 * n34 * n41 -
      n24 * n32 * n41 +
      n24 * n31 * n42 -
      n21 * n34 * n42 -
      n22 * n31 * n44 +
      n21 * n32 * n44;
    te[6] =
      n14 * n32 * n41 -
      n12 * n34 * n41 -
      n14 * n31 * n42 +
      n11 * n34 * n42 +
      n12 * n31 * n44 -
      n11 * n32 * n44;
    te[10] =
      n12 * n24 * n41 -
      n14 * n22 * n41 +
      n14 * n21 * n42 -
      n11 * n24 * n42 -
      n12 * n21 * n44 +
      n11 * n22 * n44;
    te[14] =
      n14 * n22 * n31 -
      n12 * n24 * n31 -
      n14 * n21 * n32 +
      n11 * n24 * n32 +
      n12 * n21 * n34 -
      n11 * n22 * n34;
    te[3] =
      n23 * n32 * n41 -
      n22 * n33 * n41 -
      n23 * n31 * n42 +
      n21 * n33 * n42 +
      n22 * n31 * n43 -
      n21 * n32 * n43;
    te[7] =
      n12 * n33 * n41 -
      n13 * n32 * n41 +
      n13 * n31 * n42 -
      n11 * n33 * n42 -
      n12 * n31 * n43 +
      n11 * n32 * n43;
    te[11] =
      n13 * n22 * n41 -
      n12 * n23 * n41 -
      n13 * n21 * n42 +
      n11 * n23 * n42 +
      n12 * n21 * n43 -
      n11 * n22 * n43;
    te[15] =
      n12 * n23 * n31 -
      n13 * n22 * n31 +
      n13 * n21 * n32 -
      n11 * n23 * n32 -
      n12 * n21 * n33 +
      n11 * n22 * n33;

    const det = n11 * te[0] + n21 * te[4] + n31 * te[8] + n41 * te[12];

    if (det === 0) {
      const msg = "Matrix4.getInverse(): can't invert matrix, determinant is 0";

      if (throwOnInvertible || false) {
        throw new Error(msg);
      } else {
        console.warn(msg);
      }

      this.identity();

      return this;
    }

    this.multiplyScalar(1 / det);

    return this;
  }

  isReflected() {
    const te = this.elements;

    const m0 = te[0];
    const m3 = te[4];
    const m6 = te[8];
    const m1 = te[1];
    const m4 = te[5];
    const m7 = te[9];
    const m2 = te[2];
    const m5 = te[6];
    const m8 = te[10];

    const determinant =
      m0 * m4 * m8 + // +aei
      m1 * m5 * m6 + // +bfg
      m2 * m3 * m7 - // +cdh
      m2 * m4 * m6 - // -ceg
      m1 * m3 * m8 - // -bdi
      m0 * m5 * m7; // -afh

    return determinant < 0;
  }

  scale(v) {
    const te = this.elements;
    const {x} = v;
    const {y} = v;
    const {z} = v;

    te[0] *= x;
    te[4] *= y;
    te[8] *= z;
    te[1] *= x;
    te[5] *= y;
    te[9] *= z;
    te[2] *= x;
    te[6] *= y;
    te[10] *= z;
    te[3] *= x;
    te[7] *= y;
    te[11] *= z;

    return this;
  }

  getMaxScaleOnAxis() {
    const te = this.elements;

    const scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
    const scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
    const scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

    return Math.sqrt(Math.max(scaleXSq, Math.max(scaleYSq, scaleZSq)));
  }

  makeFrustum(left, right, bottom, top, near, far) {
    const te = this.elements;

    const x = (2 * near) / (right - left);
    const y = (2 * near) / (top - bottom);

    const a = (right + left) / (right - left);
    const b = (top + bottom) / (top - bottom);
    const c = -(far + near) / (far - near);
    const d = (-2 * far * near) / (far - near);

    te[0] = x;
    te[4] = 0;
    te[8] = a;
    te[12] = 0;
    te[1] = 0;
    te[5] = y;
    te[9] = b;
    te[13] = 0;
    te[2] = 0;
    te[6] = 0;
    te[10] = c;
    te[14] = d;
    te[3] = 0;
    te[7] = 0;
    te[11] = -1;
    te[15] = 0;

    return this;
  }

  makePerspective(fov, aspect, near, far) {
    const ymax = near * Math.tan(degToRad(fov * 0.5));
    const ymin = -ymax;
    const xmin = ymin * aspect;
    const xmax = ymax * aspect;

    return this.makeFrustum(xmin, xmax, ymin, ymax, near, far);
  }

  makeOrthographic(left, right, top, bottom, near, far) {
    const te = this.elements;
    const w = 1.0 / (right - left);
    const h = 1.0 / (top - bottom);
    const p = 1.0 / (far - near);

    const x = (right + left) * w;
    const y = (top + bottom) * h;
    const z = (far + near) * p;

    te[0] = 2 * w;
    te[4] = 0;
    te[8] = 0;
    te[12] = -x;
    te[1] = 0;
    te[5] = 2 * h;
    te[9] = 0;
    te[13] = -y;
    te[2] = 0;
    te[6] = 0;
    te[10] = -2 * p;
    te[14] = -z;
    te[3] = 0;
    te[7] = 0;
    te[11] = 0;
    te[15] = 1;

    return this;
  }

  isEqual(m) {
    const me = m.elements;
    const te = this.elements;

    if (
      te[0] === me[0] &&
      te[4] === me[4] &&
      te[8] === me[8] &&
      te[12] === me[12] &&
      te[1] === me[1] &&
      te[5] === me[5] &&
      te[9] === me[9] &&
      te[13] === me[13] &&
      te[2] === me[2] &&
      te[6] === me[6] &&
      te[10] === me[10] &&
      te[14] === me[14] &&
      te[3] === me[3] &&
      te[7] === me[7] &&
      te[11] === me[11] &&
      te[15] === me[15]
    ) {
      return true;
    }
    return false;
  }

  clone() {
    const te = this.elements;

    return new Matrix4(
      te[0],
      te[4],
      te[8],
      te[12],
      te[1],
      te[5],
      te[9],
      te[13],
      te[2],
      te[6],
      te[10],
      te[14],
      te[3],
      te[7],
      te[11],
      te[15]
    );
  }

  isIdentity() {
    const te = this.elements;

    if (
      te[0] === 1 &&
      te[4] === 0 &&
      te[8] === 0 &&
      te[12] === 0 &&
      te[1] === 0 &&
      te[5] === 1 &&
      te[9] === 0 &&
      te[13] === 0 &&
      te[2] === 0 &&
      te[6] === 0 &&
      te[10] === 1 &&
      te[14] === 0 &&
      te[3] === 0 &&
      te[7] === 0 &&
      te[11] === 0 &&
      te[15] === 1
    ) {
      return true;
    }
    return false;
  }

  // return true if elements are with digits of identity
  isNearlyIdentity(digits) {
    const snapped = this.clone().snap(digits);
    return snapped.isIdentity();
  }
}

const mRotation = new Matrix4();
const mScale = new Matrix4();
const v1 = new Vector3();
const matrix = new Matrix4();
const zAlloc = new Vector3();
const xAlloc = new Vector3();
const yAlloc = new Vector3();

// Out of order assignment to avoid use before defign
Vector3.prototype.unproject = function unproject(camera) {
  matrix.multiplyMatrices(camera.matrixWorld, matrix.getInverse(camera.projectionMatrix));
  return this.applyMatrix4(matrix);
}

Matrix3.prototype.getMatrix4 = function getMatrix4() {
  const m = this.elements;
  return new Matrix4(m[0], m[3], m[6], 0, m[1], m[4], m[7], 0, m[2], m[5], m[8], 0);
}

Matrix4.prototype.decompose = function decompose(translation, rotation, scale) {
  const te = this.elements;

  // grab the axis vectors
  xAlloc.set(te[0], te[1], te[2]);
  yAlloc.set(te[4], te[5], te[6]);
  zAlloc.set(te[8], te[9], te[10]);

  translation = translation instanceof Vector3 ? translation : new Vector3();
  rotation = rotation instanceof Quaternion ? rotation : new Quaternion();
  scale = scale instanceof Vector3 ? scale : new Vector3();

  scale.x = xAlloc.length();
  scale.y = yAlloc.length();
  scale.z = zAlloc.length();

  translation.x = te[12];
  translation.y = te[13];
  translation.z = te[14];

  // scale the rotation part
  matrix.copy(this);

  matrix.elements[0] /= scale.x;
  matrix.elements[1] /= scale.x;
  matrix.elements[2] /= scale.x;

  matrix.elements[4] /= scale.y;
  matrix.elements[5] /= scale.y;
  matrix.elements[6] /= scale.y;

  matrix.elements[8] /= scale.z;
  matrix.elements[9] /= scale.z;
  matrix.elements[10] /= scale.z;

  rotation.setFromRotationMatrix(matrix);

  return [translation, rotation, scale];
}

Matrix4.prototype.lookAt = function lookAt(eye, target, up) {
  const te = this.elements;

  zAlloc.subVectors(eye, target).normalize();

  if (zAlloc.length() === 0) {
    zAlloc.z = 1;
  }

  xAlloc.crossVectors(up, zAlloc).normalize();

  if (xAlloc.length() === 0) {
    zAlloc.x += 0.0001;
    xAlloc.crossVectors(up, zAlloc).normalize();
  }

  yAlloc.crossVectors(zAlloc, xAlloc);

  te[0] = xAlloc.x;
  te[4] = yAlloc.x;
  te[8] = zAlloc.x;
  te[1] = xAlloc.y;
  te[5] = yAlloc.y;
  te[9] = zAlloc.y;
  te[2] = xAlloc.z;
  te[6] = yAlloc.z;
  te[10] = zAlloc.z;

  return this;
}

/// Return scale factor present in trnsformation matrix
Matrix4.prototype.getScale = function getScale(scale) {
  const te = this.elements;
  scale = scale instanceof Vector3 ? scale : new Vector3();
  // grab the axis vectors
  xAlloc.set(te[0], te[1], te[2]);
  yAlloc.set(te[4], te[5], te[6]);
  zAlloc.set(te[8], te[9], te[10]);

  scale.x = xAlloc.length();
  scale.y = yAlloc.length();
  scale.z = zAlloc.length();

  return scale;
}

Matrix4.prototype.compose = function compose(translation, rotation, scale) {
  const te = this.elements;

  mRotation.identity();
  mRotation.setRotationFromQuaternion(rotation);

  mScale.makeScale(scale.x, scale.y, scale.z);

  this.multiplyMatrices(mRotation, mScale);

  te[12] = translation.x;
  te[13] = translation.y;
  te[14] = translation.z;

  return this;
}

Matrix4.prototype.getPosition = function getPosition() {
  console.warn(
    "DEPRECATED: Matrix4's .getPosition() has been removed. Use Vector3.getPositionFromMatrix( matrix ) instead."
  );

  const te = this.elements;
  return v1.set(te[12], te[13], te[14]);
}

/** @constructor */
/** @constructor */
export class Ray {
  constructor(origin, direction) {
    this.origin = origin !== undefined ? origin : new Vector3();

    this.direction = direction !== undefined ? direction : new Vector3();
  }

  set(origin, direction) {
    this.origin.copy(origin);
    this.direction.copy(direction);

    return this;
  }

  copy(ray) {
    this.origin.copy(ray.origin);
    this.direction.copy(ray.direction);

    return this;
  }

  at(t, optionalTarget) {
    const result = optionalTarget || new Vector3();

    return result.copy(this.direction).multiplyScalar(t).add(this.origin);
  }

  recast(t) {
    this.origin.copy(this.at(t, v1));

    return this;
  }

  closestPointToPoint(point, optionalTarget) {
    const result = optionalTarget || new Vector3();
    result.subVectors(point, this.origin);
    const directionDistance = result.dot(this.direction);

    // returns a point on this ray
    return result.copy(this.direction).multiplyScalar(directionDistance).add(this.origin);
  }

  distanceToPoint(point) {
    const directionDistance = v1.subVectors(point, this.origin).dot(this.direction);
    v1.copy(this.direction).multiplyScalar(directionDistance).add(this.origin);
    return v1.distanceTo(point);
  }

  // eslint-disable-next-line class-methods-use-this
  isIntersectionCylinder() {}

  isIntersectionSphere(sphere) {
    return this.distanceToPoint(sphere.center) <= sphere.radius;
  }

  isIntersectionPlane(plane) {
    const denominator = plane.normal.dot(this.direction);

    // plane and ray are not perpendicular
    if (denominator !== 0) return true;

    if (plane.distanceToPoint(this.origin) === 0) return true;

    return false;
  }

  distanceToPlane(plane) {
    const denominator = plane.normal.dot(this.direction);
    if (denominator === 0) {
      // line is coplanar
      if (plane.distanceToPoint(this.origin) === 0) return 0;

      // ray is parallel
      return undefined;
    }

    const t = -(this.origin.dot(plane.normal) + plane.constant) / denominator;

    return t;
  }

  intersectPlane(plane, optionalTarget) {
    const t = this.distanceToPlane(plane);

    if (t === undefined) return undefined;

    return this.at(t, optionalTarget);
  }

  applyMatrix4(matrix4) {
    this.direction.add(this.origin).applyMatrix4(matrix4);
    this.origin.applyMatrix4(matrix4);
    this.direction.sub(this.origin);

    return this;
  }

  equals(ray) {
    return ray.origin.equals(this.origin) && ray.direction.equals(this.direction);
  }

  clone() {
    return new Ray().copy(this);
  }
}
