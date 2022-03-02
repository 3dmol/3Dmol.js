export class Quaternion {
  constructor(x, y, z, w) {
    this.x = x || 0;
    this.y = y || 0;
    this.z = z || 0;
    this.w = (w !== undefined) ? w : 1;
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

    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z
      + this.w * this.w);
  }

  lengthxyz() {

    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }

  normalize() {

    var l = this.length();

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

    var qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
    var qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;

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
    var c1 = Math.cos(e.x / 2);
    var c2 = Math.cos(e.y / 2);
    var c3 = Math.cos(e.z / 2);
    var s1 = Math.sin(e.x / 2);
    var s2 = Math.sin(e.y / 2);
    var s3 = Math.sin(e.z / 2);

    this.x = s1 * c2 * c3 + c1 * s2 * s3;
    this.y = c1 * s2 * c3 - s1 * c2 * s3;
    this.z = c1 * c2 * s3 + s1 * s2 * c3;
    this.w = c1 * c2 * c3 - s1 * s2 * s3;

    return this;
  }
};

