/*
 * math-like functionality
 * quaternion, vector, matrix
 */

var $3Dmol = $3Dmol || {};
$3Dmol.Math = {

    clamp : function(x, min, max) {
        return Math.min(Math.max(x, min), max);
    },

    degToRad : function() {
        var degreeToRadiansFactor = Math.PI / 180;

        return function(deg) {
            return deg * degreeToRadiansFactor;
        };

    }()

};

// Quaternion
/** @constructor */
$3Dmol.Quaternion = function(x, y, z, w) {

    this.x = x || 0;
    this.y = y || 0;
    this.z = z || 0;
    this.w = (w !== undefined) ? w : 1;

};

$3Dmol.Quaternion.prototype = {

    constructor : $3Dmol.Quaternion,

    set : function(x, y, z, w) {

        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;

        return this;
    },

    copy : function(q) {

        this.x = q.x;
        this.y = q.y;
        this.z = q.z;
        this.w = q.w;

        return this;
    },

    conjugate : function() {

        this.x *= -1;
        this.y *= -1;
        this.z *= -1;

        return this;
    },

    inverse : function() {

        return this.conjugate().normalize();
    },

    length : function() {

        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z
                + this.w * this.w);
    },

    lengthxyz : function() {

        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    },

    normalize : function() {

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

    },

    multiply : function(q) {

        return this.multiplyQuaternions(this, q);
    },

    multiplyScalar : function(s) {
        this.x *= s;
        this.y *= s;
        this.z *= s;
        this.w *= s;
        return this;
    },

    multiplyQuaternions : function(a, b) {

        var qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
        var qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;

        this.x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
        this.y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
        this.z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
        this.w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;
        return this;
    },

    sub : function(q) {
        this.x -= q.x;
        this.y -= q.y;
        this.z -= q.z;
        this.w -= q.w;
        return this;
    },

    clone : function() {
        return new $3Dmol.Quaternion(this.x, this.y, this.z, this.w);
    },
    setFromEuler : function(e) {
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

// A 2 Vector
/** @constructor */
$3Dmol.Vector2 = function(x, y) {

    this.x = x || 0.0;
    this.y = y || 0.0;
};

$3Dmol.Vector2.prototype = {

    constructor : $3Dmol.Vector2,

    set : function(x, y) {

        this.x = x;
        this.y = y;

        return this;
    },

    subVectors : function(a, b) {

        this.x = a.x - b.x;
        this.y = a.y - b.y;

        return this;
    },

    copy : function(v) {

        this.x = v.x;
        this.y = v.y;

        return this;
    },

    clone : function() {

        return new $3Dmol.Vector2(this.x, this.y);
    }

};

// A 3 Vector

$3Dmol.Vector3 = function(x, y, z) {
    this.x = x || 0.0;
    this.y = y || 0.0;
    this.z = z || 0.0;
};
/** @this {$3Dmol.Vector3} */
$3Dmol.Vector3.prototype = {

    constructor : $3Dmol.Vector3,

    set : function(x, y, z) {

        this.x = x;
        this.y = y;
        this.z = z;

        return this;
    },

    copy : function(v) {

        this.x = v.x;
        this.y = v.y;
        this.z = v.z;

        return this;
    },

    add : function(v) {

        this.x += v.x;
        this.y += v.y;
        this.z += v.z;

        return this;
    },

    addVectors : function(a, b) {

        this.x = a.x + b.x;
        this.y = a.y + b.y;
        this.z = a.z + b.z;

        return this;
    },

    multiplyVectors: function(a, b) { //elementwise
        this.x = a.x * b.x;
        this.y = a.y * b.y;
        this.z = a.z * b.z;

        return this;
    },
    sub : function(v) {

        this.x -= v.x;
        this.y -= v.y;
        this.z -= v.z;

        return this;
    },

    subVectors : function(a, b) {

        this.x = a.x - b.x;
        this.y = a.y - b.y;
        this.z = a.z - b.z;

        return this;
    },

    multiplyScalar : function(s) {

        this.x *= s;
        this.y *= s;
        this.z *= s;

        return this;
    },

    divideScalar : function(s) {

        if (s !== 0) {
            this.x /= s;
            this.y /= s;
            this.z /= s;
        }

        else {
            this.x = 0;
            this.y = 0;
            this.z = 0;
        }

        return this;
    },

    // accumulate maximum
    max : function(s) {

        this.x = Math.max(this.x, s.x);
        this.y = Math.max(this.y, s.y);
        this.z = Math.max(this.z, s.z);

        return this;
    },

    // accumulate min
    min : function(s) {

        this.x = Math.min(this.x, s.x);
        this.y = Math.min(this.y, s.y);
        this.z = Math.min(this.z, s.z);

        return this;
    },
    distanceTo : function(v) {
        return Math.sqrt(this.distanceToSquared(v));
    },

    distanceToSquared : function(v) {
        var dx = this.x - v.x;
        var dy = this.y - v.y;
        var dz = this.z - v.z;

        return dx * dx + dy * dy + dz * dz;
    },

    applyMatrix3 : function(m) {
        
       var x = this.x, y = this.y, z = this.z;

        var e = m.elements;
        //column major ordering
        this.x = e[0] * x + e[3] * y + e[6] * z;
        this.y = e[1] * x + e[4] * y + e[7] * z;
        this.z = e[2] * x + e[5] * y + e[8] * z;

        return this;       
    },
    
    applyMatrix4 : function(m) {

        var x = this.x, y = this.y, z = this.z;

        var e = m.elements;

        this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
        this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
        this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

        return this;
    },

    applyProjection : function(m) {

        // input: $3Dmol.Matrix4 projection matrix

        var x = this.x, y = this.y, z = this.z;

        var e = m.elements;
        var d = (e[3] * x + e[7] * y + e[11] * z + e[15]);

        this.x = (e[0] * x + e[4] * y + e[8] * z + e[12]) / d;
        this.y = (e[1] * x + e[5] * y + e[9] * z + e[13]) / d;
        this.z = (e[2] * x + e[6] * y + e[10] * z + e[14]) / d;

        return this;
    },

    applyQuaternion : function(q) {

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
    },

    negate : function() {

        return this.multiplyScalar(-1);
    },

    dot : function(v) {

        return this.x * v.x + this.y * v.y + this.z * v.z;
    },

    length : function() {

        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    },

    lengthSq : function() {

        return (this.x * this.x + this.y * this.y + this.z * this.z);
    },

    normalize : function() {

        return this.divideScalar(this.length());
    },

    cross : function(v) {

        var x = this.x, y = this.y, z = this.z;

        this.x = y * v.z - z * v.y;
        this.y = z * v.x - x * v.z;
        this.z = x * v.y - y * v.x;

        return this;
    },

    crossVectors : function(a, b) {

        this.x = a.y * b.z - a.z * b.y;
        this.y = a.z * b.x - a.x * b.z;
        this.z = a.x * b.y - a.y * b.x;

        return this;
    },

    getPositionFromMatrix : function(m) {

        this.x = m.elements[12];
        this.y = m.elements[13];
        this.z = m.elements[14];

        return this;
    },

    setEulerFromRotationMatrix : function(m, order) {

        // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

        var te = m.elements;
        var m11 = te[0], m12 = te[4], m13 = te[8];
        // var m21 = te[1];
        var m22 = te[5], m23 = te[9];
        // var m31 = te[2];
        var m32 = te[6], m33 = te[10];

        if (order === undefined || order === 'XYZ') {

            this.y = Math.asin($3Dmol.Math.clamp(m13, -1, 1));

            if (Math.abs(m13) < 0.99999) {

                this.x = Math.atan2(-m23, m33);
                this.z = Math.atan2(-m12, m11);

            } else {

                this.x = Math.atan2(m32, m22);
                this.z = 0;

            }
        }

        else {
            console
                    .error("Error with vector's setEulerFromRotationMatrix: Unknown order: "
                            + order);
        }

        return this;

    },

    rotateAboutVector : function(axis, ang) {

        axis.normalize();
        var cosang = Math.cos(ang);
        var sinang = Math.sin(ang);
        // Rodrigues' rotation formula, from wikipedia

        var term1 = this.clone().multiplyScalar(cosang);
        var term2 = (axis.clone().cross(this)).multiplyScalar(sinang);
        var term3 = axis.clone().multiplyScalar(axis.clone().dot(this))
                .multiplyScalar(1 - cosang);

        var rot = term1.add(term2).add(term3);

        this.x = rot.x;
        this.y = rot.y;
        this.z = rot.z;

        return this;
    },

    setFromMatrixPosition : function(m) {

        var e = m.elements;

        this.x = e[12];
        this.y = e[13];
        this.z = e[14];

        return this;

    },
    // unproject is defined after Matrix4

    transformDirection : function(m) {

        // input: THREE.Matrix4 affine matrix
        // vector interpreted as a direction

        var x = this.x, y = this.y, z = this.z;
        var e = m.elements;

        this.x = e[0] * x + e[4] * y + e[8] * z;
        this.y = e[1] * x + e[5] * y + e[9] * z;
        this.z = e[2] * x + e[6] * y + e[10] * z;

        return this.normalize();
    },

    clone : function() {
        return new $3Dmol.Vector3(this.x, this.y, this.z);
    }

};

// Matrices

// Matrix3
/** @constructor */
$3Dmol.Matrix3 = function(n11, n12, n13, n21, n22, n23, n31, n32, n33) {

    this.elements = new Float32Array(9);

    this.set((n11 !== undefined) ? n11 : 1, n12 || 0, n13 || 0, n21 || 0,
            (n22 !== undefined) ? n22 : 1, n23 || 0, n31 || 0, n32 || 0,
            (n33 !== undefined) ? n33 : 1);

};

$3Dmol.square = function(n) {
        return n*n;
};
    
//return conversion matrix given crystal unit cell parameters
$3Dmol.conversionMatrix3 = function(a, b, c, alpha, beta, gamma) {
    //convert to radians
    alpha = alpha * Math.PI / 180;
    beta = beta * Math.PI / 180;
    gamma = gamma * Math.PI / 180;
    let sqr = $3Dmol.square;
    let cos_alpha = Math.cos(alpha);
    let cos_beta = Math.cos(beta);
    let cos_gamma = Math.cos(gamma);
    let sin_gamma = Math.sin(gamma);
    let conversionMatrix = new $3Dmol.Matrix3(
            a, b*cos_gamma, c*cos_beta,
            0, b*sin_gamma, c*(cos_alpha-cos_beta*cos_gamma)/sin_gamma,
            0, 0, c*Math.sqrt(1-sqr(cos_alpha)-sqr(cos_beta)-sqr(cos_gamma)+2*cos_alpha*cos_beta*cos_gamma)/sin_gamma);
    return conversionMatrix;
};

$3Dmol.Matrix3.prototype = {

    constructor : $3Dmol.Matrix3,

    set : function(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
        var te = this.elements;

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
    },

    identity : function() {
        this.set(1, 0, 0, 0, 1, 0, 0, 0, 1);

        return this;
    },

    copy : function(m) {
        var me = m.elements;

        this.set(me[0], me[3], me[6], me[1], me[4], me[7], me[2], me[5], me[8]);
    },

    multiplyScalar : function(s) {
        var te = this.elements;

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
    },
    

    getInverse3 : function(matrix) {
        // input: Matrix3

        let me = matrix.elements;
        let te = this.elements;
        
        te[0] = me[4] * me[8] - me[5] * me[7];
        te[3] = me[6] * me[5] - me[3] * me[8];
        te[6] = me[3] * me[7] - me[6] * me[4];
        te[1] = me[7] * me[2] - me[1] * me[8];
        te[4] = me[0] * me[8] - me[6] * me[2];
        te[7] = me[1] * me[6] - me[0] * me[7];
        te[2] = me[1] * me[5] - me[2] * me[4];
        te[5] = me[2] * me[3] - me[0] * me[5];
        te[8] = me[0] * me[4] - me[1] * me[3];

        let det = me[0] * te[0] + me[3]* te[1] + me[6]*te[2];
        this.multiplyScalar(1.0 / det);

        return this;
    },
    
    getInverse : function(matrix, throwOnInvertible) {
        // input: Matrix4

        var me = matrix.elements;
        var te = this.elements;

        te[0] = me[10] * me[5] - me[6] * me[9];
        te[1] = -me[10] * me[1] + me[2] * me[9];
        te[2] = me[6] * me[1] - me[2] * me[5];
        te[3] = -me[10] * me[4] + me[6] * me[8];
        te[4] = me[10] * me[0] - me[2] * me[8];
        te[5] = -me[6] * me[0] + me[2] * me[4];
        te[6] = me[9] * me[4] - me[5] * me[8];
        te[7] = -me[9] * me[0] + me[1] * me[8];
        te[8] = me[5] * me[0] - me[1] * me[4];

        var det = me[0] * te[0] + me[1] * te[3] + me[2] * te[6];

        // no inverse

        if (det === 0) {

            var msg = "Matrix3.getInverse(): can't invert matrix, determinant is 0";

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
    },

    // https://en.wikipedia.org/wiki/Determinant
    getDeterminant : function() {
        var m = this.elements;

        /*
         * |a b c| |d e f| |g h i|
         */

        var determinant = m[0] * m[4] * m[8] // +aei
                + m[1] * m[5] * m[6] // +bfg
                + m[2] * m[3] * m[7] // +cdh
                - m[2] * m[4] * m[6] // -ceg
                - m[1] * m[3] * m[8] // -bdi
                - m[0] * m[5] * m[7];// -afh
        return determinant;
    },

    getMatrix4 : function() {
      var m = this.elements;
      return new $3Dmol.Matrix4(m[0],m[3],m[6],0,m[1],m[4],m[7],0,m[2],m[5],m[8],0);
    },
    
    transpose : function() {
        var tmp, m = this.elements;

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
    },

    clone : function() {
        var te = this.elements;

        return new $3Dmol.Matrix3(

        te[0], te[3], te[6], te[1], te[4], te[7], te[2], te[5], te[8]

        );
    }

};

// Matrix 4
/** @constructor */
$3Dmol.Matrix4 = function(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32,
        n33, n34, n41, n42, n43, n44) {

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
};

$3Dmol.Matrix4.prototype = {

    constructor : $3Dmol.Matrix4,

    set : function(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34,
            n41, n42, n43, n44) {
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
    },

    identity : function() {
        this.set(

        1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1

        );

        return this;
    },

    copy : function(m) {
        var me = m.elements;

        this.set(

        me[0], me[4], me[8], me[12], me[1], me[5], me[9], me[13], me[2], me[6],
                me[10], me[14], me[3], me[7], me[11], me[15]

        );

        return this;
    },

    matrix3FromTopLeft : function() {
        var te = this.elements;
        return new $3Dmol.Matrix3(te[0], te[4], te[8], te[1], te[5], te[9],
                te[2], te[6], te[10]);
    },

    setRotationFromEuler : function(v, order) {

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

    },

    setRotationFromQuaternion : function(q) {
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
    },

    lookAt : function() {
        var x = new $3Dmol.Vector3();
        var y = new $3Dmol.Vector3();
        var z = new $3Dmol.Vector3();

        return function(eye, target, up) {

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

    }(),

    multiplyMatrices : function(a, b) {
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
    },

    multiplyScalar : function(s) {
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
    },

    makeTranslation : function(x, y, z) {

        this.set(

        1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1

        );

        return this;

    },
    
    //snap values close to integers to their integer value
    //useful and identifying identity matrices
    snap : function(digits) {
        if(!digits) digits = 4;
        let mult = Math.pow(10,4);
        let te = this.elements;
        for(let i = 0; i < 16; i++) {
            let rounded = Math.round(te[i]);
            if(rounded == Math.round(te[i]*mult)/mult) {
                te[i] = rounded;
            }
        }
        return this;
    },

    transpose : function() {
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
    },

    getPosition : function() {
        var v1 = new $3Dmol.Vector3();

        return function() {

            console
                    .warn('DEPRECATED: Matrix4\'s .getPosition() has been removed. Use Vector3.getPositionFromMatrix( matrix ) instead.');

            var te = this.elements;
            return v1.set(te[12], te[13], te[14]);
        };

    }(),

    setPosition : function(v) {
        var te = this.elements;

        te[12] = v.x;
        te[13] = v.y;
        te[14] = v.z;

        return this;
    },
    
    translate : function(v) {
        var te = this.elements;

        te[12] += v.x;
        te[13] += v.y;
        te[14] += v.z;

        return this;
    },

    getInverse : function(m, throwOnInvertible) {
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
    },

    isReflected : function() {
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
    },

    compose : function() {
        var mRotation = new $3Dmol.Matrix4(), mScale = new $3Dmol.Matrix4();

        return function(translation, rotation, scale) {

            var te = this.elements;

            mRotation.identity();
            mRotation.setRotationFromQuaternion(rotation);

            mScale.makeScale(scale.x, scale.y, scale.z);

            this.multiplyMatrices(mRotation, mScale);

            te[12] = translation.x;
            te[13] = translation.y;
            te[14] = translation.z;

            return this;

        };
    }(),
    /// Return scale factor present in trnsformation matrix
    getScale : function() {
        var x = new $3Dmol.Vector3(), y = new $3Dmol.Vector3(), z = new $3Dmol.Vector3();

        return function(scale) {

            var te = this.elements;
            scale = (scale instanceof $3Dmol.Vector3) ? scale
                    : new $3Dmol.Vector3();
            // grab the axis vectors
            x.set(te[0], te[1], te[2]);
            y.set(te[4], te[5], te[6]);
            z.set(te[8], te[9], te[10]);

            scale.x = x.length();
            scale.y = y.length();
            scale.z = z.length();

            return scale;

        };
    }(),
    decompose : function() {
        var x = new $3Dmol.Vector3(), y = new $3Dmol.Vector3(), z = new $3Dmol.Vector3(), matrix = new $3Dmol.Matrix4();

        return function(translation, rotation, scale) {

            var te = this.elements;

            // grab the axis vectors
            x.set(te[0], te[1], te[2]);
            y.set(te[4], te[5], te[6]);
            z.set(te[8], te[9], te[10]);

            translation = (translation instanceof $3Dmol.Vector3) ? translation
                    : new $3Dmol.Vector3();
            rotation = (rotation instanceof $3Dmol.Quaternion) ? rotation
                    : new $3Dmol.Quaternion();
            scale = (scale instanceof $3Dmol.Vector3) ? scale
                    : new $3Dmol.Vector3();

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

            return [ translation, rotation, scale ];

        };
    }(),

    scale : function(v) {
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
    },

    getMaxScaleOnAxis : function() {

        var te = this.elements;

        var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
        var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
        var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

        return Math.sqrt(Math.max(scaleXSq, Math.max(scaleYSq, scaleZSq)));

    },

    makeFrustum : function(left, right, bottom, top, near, far) {
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
    },

    makePerspective : function(fov, aspect, near, far) {
        var ymax = near * Math.tan($3Dmol.Math.degToRad(fov * 0.5));
        var ymin = -ymax;
        var xmin = ymin * aspect;
        var xmax = ymax * aspect;

        return this.makeFrustum(xmin, xmax, ymin, ymax, near, far);
    },

    makeOrthographic : function(left, right, top, bottom, near, far) {

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

    },

    isEqual : function(m) {
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
    },

    clone : function() {
        var te = this.elements;

        return new $3Dmol.Matrix4(

        te[0], te[4], te[8], te[12], te[1], te[5], te[9], te[13], te[2], te[6],
                te[10], te[14], te[3], te[7], te[11], te[15]

        );
    },

    isIdentity : function() {
        var te = this.elements;

        if (te[0] == 1 && te[4] == 0 && te[8] == 0 && te[12] == 0 && te[1] == 0
                && te[5] == 1 && te[9] == 0 && te[13] == 0 && te[2] == 0
                && te[6] == 0 && te[10] == 1 && te[14] == 0 && te[3] == 0
                && te[7] == 0 && te[11] == 0 && te[15] == 1) {
            return true;
        } else {
            return false;
        }
    },
    
    //return true if elements are with digits of identity
    isNearlyIdentity : function(digits) {
        let snapped = this.clone().snap(digits);
        return snapped.isIdentity();
    }

};

$3Dmol.Vector3.prototype.unproject = function() {

    var matrix = new $3Dmol.Matrix4();

    return function unproject(camera) {

        matrix.multiplyMatrices(camera.matrixWorld, matrix
                .getInverse(camera.projectionMatrix));
        return this.applyMatrix4(matrix);

    };

}();

/** @constructor */
$3Dmol.Ray = function(origin, direction) {

    this.origin = (origin !== undefined) ? origin : new $3Dmol.Vector3();

    this.direction = (direction !== undefined) ? direction
            : new $3Dmol.Vector3();

};

// TODO: Remove methods we don't need (intersectPlane ??)
$3Dmol.Ray.prototype = {

    constructor : $3Dmol.Ray,

    set : function(origin, direction) {

        this.origin.copy(origin);
        this.direction.copy(direction);

        return this;

    },

    copy : function(ray) {

        this.origin.copy(ray.origin);
        this.direction.copy(ray.direction);

        return this;

    },

    at : function(t, optionalTarget) {

        var result = optionalTarget || new $3Dmol.Vector3();

        return result.copy(this.direction).multiplyScalar(t).add(this.origin);

    },

    recast : function() {

        var v1 = new $3Dmol.Vector3();

        return function(t) {
            this.origin.copy(this.at(t, v1));

            return this;
        };

    }(),

    closestPointToPoint : function(point, optionalTarget) {

        var result = optionalTarget || new $3Dmol.Vector3();
        result.subVectors(point, this.origin);
        var directionDistance = result.dot(this.direction);

        // returns a point on this ray
        return result.copy(this.direction).multiplyScalar(directionDistance)
                .add(this.origin);

    },

    distanceToPoint : function() {

        var v1 = new $3Dmol.Vector3();

        return function(point) {
            var directionDistance = v1.subVectors(point, this.origin).dot(
                    this.direction);
            v1.copy(this.direction).multiplyScalar(directionDistance).add(
                    this.origin);
            return v1.distanceTo(point);
        };

    }(),

    isIntersectionCylinder : function() {

    },

    isIntersectionSphere : function(sphere) {
        return (this.distanceToPoint(sphere.center) <= sphere.radius);

    },

    isIntersectionPlane : function(plane) {

        var denominator = plane.normal.dot(this.direction);

        // plane and ray are not perpendicular
        if (denominator !== 0)
            return true;

        if (plane.distanceToPoint(this.origin) === 0)
            return true;

        return false;

    },

    distanceToPlane : function(plane) {

        var denominator = plane.normal.dot(this.direction);
        if (denominator === 0) {

            // line is coplanar
            if (plane.distanceToPoint(this.origin) === 0)
                return 0;

            // ray is parallel
            return undefined;
        }

        var t = -(this.origin.dot(plane.normal) + plane.constant) / denominator;

        return t;

    },

    intersectPlane : function(plane, optionalTarget) {

        var t = this.distanceToPlane(plane);

        if (t === undefined)
            return undefined;

        return this.at(t, optionalTarget);

    },

    applyMatrix4 : function(matrix4) {

        this.direction.add(this.origin).applyMatrix4(matrix4);
        this.origin.applyMatrix4(matrix4);
        this.direction.sub(this.origin);

        return this;

    },

    equals : function(ray) {

        return ray.origin.equals(this.origin)
                && ray.direction.equals(this.direction);

    },

    clone : function() {

        return new $3Dmol.Ray().copy(this);

    }

};
