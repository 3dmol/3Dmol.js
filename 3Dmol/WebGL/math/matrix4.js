// @ts-check

import { Matrix3 } from "./Matrix3";
import { Quaternion } from "./Quaternion";
import { Vector3 } from "./Vector3";




// Matrix 4
export class Matrix4 {
  constructor(
    n11,
    n12,
    n13,
    n14,
    n21,
    n22,
    n23,
    n24,
    n31,
    n32,
    n33,
    n34,
    n41,
    n42,
    n43,
    n44
  ) {

    if (typeof (n12) === 'undefined' && typeof (n11) !== 'undefined') {
      // passing list like initialization
      this.elements = new Float32Array(n11);
    } else {
      var te = this.elements = new Float32Array(16);

      te[0] = (n11 !== undefined) ? n11 : 1;
      te[4] = n12 || 0;
      te[8] = n13 || 0;
      te[12] = n14 || 0;
      te[1] = n21 || 0;
      te[5] = (n22 !== undefined) ? n22 : 1;
      te[9] = n23 || 0;
      te[13] = n24 || 0;
      te[2] = n31 || 0;
      te[6] = n32 || 0;
      te[10] = (n33 !== undefined) ? n33 : 1;
      te[14] = n34 || 0;
      te[3] = n41 || 0;
      te[7] = n42 || 0;
      te[11] = n43 || 0;
      te[15] = (n44 !== undefined) ? n44 : 1;
    }
  }

  set(
    n11,
    n12,
    n13,
    n14,
    n21,
    n22,
    n23,
    n24,
    n31,
    n32,
    n33,
    n34,
    n41,
    n42,
    n43,
    n44
  ) {
    var te = this.elements;

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
    this.set(

      1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1

    );

    return this;
  }

  copy(m) {
    var me = m.elements;

    this.set(

      me[0], me[4], me[8], me[12], me[1], me[5], me[9], me[13], me[2], me[6],
      me[10], me[14], me[3], me[7], me[11], me[15]

    );

    return this;
  }

  matrix3FromTopLeft() {
    var te = this.elements;
    return new Matrix3(te[0], te[4], te[8], te[1], te[5], te[9],
      te[2], te[6], te[10]);
  }

  setRotationFromEuler(v, order) {

    var te = this.elements;

    var x = v.x, y = v.y, z = v.z;
    var a = Math.cos(x), b = Math.sin(x);
    var c = Math.cos(y), d = Math.sin(y);
    var e = Math.cos(z), f = Math.sin(z);

    if (order === undefined || order === 'XYZ') {

      var ae = a * e, af = a * f, be = b * e, bf = b * f;

      te[0] = c * e;
      te[4] = -c * f;
      te[8] = d;

      te[1] = af + be * d;
      te[5] = ae - bf * d;
      te[9] = -b * c;

      te[2] = bf - ae * d;
      te[6] = be + af * d;
      te[10] = a * c;

    }

    else
      console.error("Error with matrix4 setRotationFromEuler. Order: "
        + order);

    return this;

  }

  setRotationFromQuaternion(q) {
    var te = this.elements;

    var x = q.x, y = q.y, z = q.z, w = q.w;
    var x2 = x + x, y2 = y + y, z2 = z + z;
    var xx = x * x2, xy = x * y2, xz = x * z2;
    var yy = y * y2, yz = y * z2, zz = z * z2;
    var wx = w * x2, wy = w * y2, wz = w * z2;

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
    var ae = a.elements;
    var be = b.elements;
    var te = this.elements;

    var a11 = ae[0], a12 = ae[4], a13 = ae[8], a14 = ae[12];
    var a21 = ae[1], a22 = ae[5], a23 = ae[9], a24 = ae[13];
    var a31 = ae[2], a32 = ae[6], a33 = ae[10], a34 = ae[14];
    var a41 = ae[3], a42 = ae[7], a43 = ae[11], a44 = ae[15];

    var b11 = be[0], b12 = be[4], b13 = be[8], b14 = be[12];
    var b21 = be[1], b22 = be[5], b23 = be[9], b24 = be[13];
    var b31 = be[2], b32 = be[6], b33 = be[10], b34 = be[14];
    var b41 = be[3], b42 = be[7], b43 = be[11], b44 = be[15];

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
    var te = this.elements;

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

    this.set(

      1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1

    );

    return this;

  }

  snap(digits) {
    if (!digits) digits = 4;
    let mult = Math.pow(10, 4);
    let te = this.elements;
    for (let i = 0; i < 16; i++) {
      let rounded = Math.round(te[i]);
      if (rounded == Math.round(te[i] * mult) / mult) {
        te[i] = rounded;
      }
    }
    return this;
  }

  transpose() {
    var te = this.elements;
    var tmp;

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
    var te = this.elements;

    te[12] = v.x;
    te[13] = v.y;
    te[14] = v.z;

    return this;
  }

  translate(v) {
    var te = this.elements;

    te[12] += v.x;
    te[13] += v.y;
    te[14] += v.z;

    return this;
  }

  getInverse(m, throwOnInvertible) {
    // based on
    // http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
    var te = this.elements;
    var me = m.elements;

    var n11 = me[0], n12 = me[4], n13 = me[8], n14 = me[12];
    var n21 = me[1], n22 = me[5], n23 = me[9], n24 = me[13];
    var n31 = me[2], n32 = me[6], n33 = me[10], n34 = me[14];
    var n41 = me[3], n42 = me[7], n43 = me[11], n44 = me[15];

    te[0] = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34
      * n43 - n23 * n32 * n44 + n22 * n33 * n44;
    te[4] = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34
      * n43 + n13 * n32 * n44 - n12 * n33 * n44;
    te[8] = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24
      * n43 - n13 * n22 * n44 + n12 * n23 * n44;
    te[12] = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12
      * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;
    te[1] = n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34
      * n43 + n23 * n31 * n44 - n21 * n33 * n44;
    te[5] = n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34
      * n43 - n13 * n31 * n44 + n11 * n33 * n44;
    te[9] = n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24
      * n43 + n13 * n21 * n44 - n11 * n23 * n44;
    te[13] = n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11
      * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34;
    te[2] = n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34
      * n42 - n22 * n31 * n44 + n21 * n32 * n44;
    te[6] = n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34
      * n42 + n12 * n31 * n44 - n11 * n32 * n44;
    te[10] = n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11
      * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44;
    te[14] = n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11
      * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34;
    te[3] = n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33
      * n42 + n22 * n31 * n43 - n21 * n32 * n43;
    te[7] = n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33
      * n42 - n12 * n31 * n43 + n11 * n32 * n43;
    te[11] = n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11
      * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43;
    te[15] = n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11
      * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33;

    var det = n11 * te[0] + n21 * te[4] + n31 * te[8] + n41 * te[12];

    if (det === 0) {

      var msg = "Matrix4.getInverse(): can't invert matrix, determinant is 0";

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
    let te = this.elements;

    let m0 = te[0],
      m3 = te[4],
      m6 = te[8],
      m1 = te[1],
      m4 = te[5],
      m7 = te[9],
      m2 = te[2],
      m5 = te[6],
      m8 = te[10];

    let determinant = m0 * m4 * m8 // +aei
      + m1 * m5 * m6 // +bfg
      + m2 * m3 * m7 // +cdh
      - m2 * m4 * m6 // -ceg
      - m1 * m3 * m8 // -bdi
      - m0 * m5 * m7;// -afh

    return determinant < 0;
  }

  scale(v) {
    var te = this.elements;
    var x = v.x, y = v.y, z = v.z;

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

    var te = this.elements;

    var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
    var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
    var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

    return Math.sqrt(Math.max(scaleXSq, Math.max(scaleYSq, scaleZSq)));

  }

  makeFrustum(left, right, bottom, top, near, far) {
    var te = this.elements;

    var x = 2 * near / (right - left);
    var y = 2 * near / (top - bottom);

    var a = (right + left) / (right - left);
    var b = (top + bottom) / (top - bottom);
    var c = -(far + near) / (far - near);
    var d = -2 * far * near / (far - near);

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
    // @ts-ignore
    var ymax = near * Math.tan(Math.degToRad(fov * 0.5));
    var ymin = -ymax;
    var xmin = ymin * aspect;
    var xmax = ymax * aspect;

    return this.makeFrustum(xmin, xmax, ymin, ymax, near, far);
  }

  makeOrthographic(left, right, top, bottom, near, far) {

    var te = this.elements;
    var w = 1.0 / (right - left);
    var h = 1.0 / (top - bottom);
    var p = 1.0 / (far - near);

    var x = (right + left) * w;
    var y = (top + bottom) * h;
    var z = (far + near) * p;

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
    var me = m.elements;
    var te = this.elements;

    if (te[0] == me[0] && te[4] == me[4] && te[8] == me[8]
      && te[12] == me[12] && te[1] == me[1] && te[5] == me[5]
      && te[9] == me[9] && te[13] == me[13] && te[2] == me[2]
      && te[6] == me[6] && te[10] == me[10] && te[14] == me[14]
      && te[3] == me[3] && te[7] == me[7] && te[11] == me[11]
      && te[15] == me[15]) {
      return true;
    } else {
      return false;
    }
  }

  clone() {
    var te = this.elements;

    return new Matrix4(

      te[0], te[4], te[8], te[12], te[1], te[5], te[9], te[13], te[2], te[6],
      te[10], te[14], te[3], te[7], te[11], te[15]

    );
  }

  isIdentity() {
    var te = this.elements;

    if (te[0] == 1 && te[4] == 0 && te[8] == 0 && te[12] == 0 && te[1] == 0
      && te[5] == 1 && te[9] == 0 && te[13] == 0 && te[2] == 0
      && te[6] == 0 && te[10] == 1 && te[14] == 0 && te[3] == 0
      && te[7] == 0 && te[11] == 0 && te[15] == 1) {
      return true;
    } else {
      return false;
    }
  }

  isNearlyIdentity(digits) {
    let snapped = this.clone().snap(digits);
    return snapped.isIdentity();
  }

  // module level allocations used 
  lookAt(eye, target, up) {

    var te = this.elements;

    z.subVectors(eye, target).normalize();

    if (z.length() === 0) {

      z.z = 1;

    }

    x.crossVectors(up, z).normalize();

    if (x.length() === 0) {

      z.x += 0.0001;
      x.crossVectors(up, z).normalize();

    }

    y.crossVectors(z, x);

    te[0] = x.x;
    te[4] = y.x;
    te[8] = z.x;
    te[1] = x.y;
    te[5] = y.y;
    te[9] = z.y;
    te[2] = x.z;
    te[6] = y.z;
    te[10] = z.z;

    return this;
  };

  // module level allocations used 
  /**
   * @deprecated
   */
  getPosition() {
    console
      .warn('DEPRECATED: Matrix4\'s .getPosition() has been removed. Use Vector3.getPositionFromMatrix( matrix ) instead.');

    var te = this.elements;
    return v1.set(te[12], te[13], te[14]);
  };

  // module level allocations used 
  compose(translation, rotation, scale) {
    var te = this.elements;

    mRotation.identity();
    mRotation.setRotationFromQuaternion(rotation);

    // TODO: investigate this fn
    // 
    mScale.makeScale(scale.x, scale.y, scale.z);

    this.multiplyMatrices(mRotation, mScale);

    te[12] = translation.x;
    te[13] = translation.y;
    te[14] = translation.z;

    return this;

  };

  // module level allocations used 
  /// Return scale factor present in trnsformation matrix
  getScale(scale) {
    var te = this.elements;
    scale = (scale instanceof Vector3) ? scale
      : new Vector3();
    // grab the axis vectors
    x.set(te[0], te[1], te[2]);
    y.set(te[4], te[5], te[6]);
    z.set(te[8], te[9], te[10]);

    scale.x = x.length();
    scale.y = y.length();
    scale.z = z.length();

    return scale;
  };


  decompose() {
    return function (translation, rotation, scale) {

      var te = this.elements;

      // grab the axis vectors
      x.set(te[0], te[1], te[2]);
      y.set(te[4], te[5], te[6]);
      z.set(te[8], te[9], te[10]);

      translation = (translation instanceof Vector3) ? translation
        : new Vector3();
      rotation = (rotation instanceof Quaternion) ? rotation
        : new Quaternion();
      scale = (scale instanceof Vector3) ? scale
        : new Vector3();

      scale.x = x.length();
      scale.y = y.length();
      scale.z = z.length();

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

    };
  }
}

// tmp allocation for lookAt
let x = new Vector3();
let y = new Vector3();
let z = new Vector3();
// tmp allocation for getPosition
let v1 = new Vector3();
// tmp allocation for decompose
// note uses the same private temp vectors as Matrix4.lookAt()
let matrix = new Matrix4();
// tmp matrix allocation for compose
let mRotation = new Matrix4()
let mScale = new Matrix4(); 