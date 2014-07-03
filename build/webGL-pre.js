/*
* math-like functionality
* quaternion, vector, matrix
*/

//Math functions
var WebMol = WebMol || {};

WebMol.Math = {

    clamp : function(x, min, max) {
        return Math.min( Math.max( x, min ), max );
    },

    degToRad : function() {
       var degreeToRadiansFactor = Math.PI / 180;
       
       return function(deg) {
           return deg * degreeToRadiansFactor;
       };
    
    }()
    
};


// Quaternion

WebMol.Quaternion = function(x, y, z, w) {

    this.x = x || 0;
    this.y = y || 0;
    this.z = z || 0;
    this.w = (w !== undefined) ? w : 1;

};

WebMol.Quaternion.prototype = {

    constructor : WebMol.Quaternion,

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
        
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w);
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

    multiplyQuaternions : function(a, b) {

        var qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
        var qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;

        this.x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
        this.y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
        this.z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
        this.w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;

    }
};

//A 2 Vector
WebMol.Vector2 = function(x, y) {
    
    this.x = x || 0.0;
    this.y = y || 0.0;
};

WebMol.Vector2.prototype = {
    
    constructor : WebMol.Vector2,
   
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
        
        return new WebMol.Vector2(this.x, this.y);
    }    
   
};

//A 3 Vector

WebMol.Vector3 = function(x, y, z) {
    this.x = x || 0.0;
    this.y = y || 0.0;
    this.z = z || 0.0;
};

WebMol.Vector3.prototype =  {
    
    constructor : WebMol.Vector3,
    
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
    

    distanceTo: function(v) {
        return Math.sqrt(this.distanceToSquared(v));
    },

    distanceToSquared: function(v) {
        var dx = this.x - v.x;
        var dy = this.y - v.y;
        var dz = this.z - v.z;

        return dx * dx + dy * dy + dz * dz;
    },
    
    applyMatrix4 : function(m) {
    
        var x = this.x, y = this.y, z = this.z;
        
        var e = m.elements;
        
        this.x = e[0]*x + e[4]*y + e[8]*z + e[12];
        this.y = e[1]*x + e[5]*y + e[9]*z + e[13];
        this.z = e[2]*x + e[6]*y + e[10]*z + e[14];
        
        return this;
    },
    
    applyProjection : function(m) {
        
        //input: WebMol.Matrix4 projection matrix
        
        var x = this.x, y = this.y, z = this.z;
        
        var e = m.elements;
        var d = ( e[3]*x + e[7]*y + e[11]*z + e[15]);
        
        this.x = (e[0]*x + e[4]*y + e[8]*z + e[12]) / d;
        this.y = (e[1]*x + e[5]*y + e[9]*z + e[13]) / d;
        this.z = (e[2]*x + e[6]*y + e[10]*z + e[14]) / d;
        
        return this;
    },
    
    applyQuaternion : function(q) { 
        
        var x = this.x;
        var y = this.y;
        var z = this.z;
        
        var qx = q.x;
        var qy = q.y;
        var qz = q.z;
        var qw = q.w;
        
        // calculate quaternion * vector
        
        var ix = qw * x + qy * z - qz * y;
        var iy = qw * y + qz * x - qx * z;
        var iz = qw * z + qx * y - qy * x;
        var iw = -qw * x - qy * y - qz * z;
        
        // calculate result * inverse quaternion
        
        this.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
        this.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
        this.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;
        
        return this;
    },
    
    negate : function() {
        
        return this.multiplyScalar(-1);
    },
    
    dot : function(v) {
        
        return this.x * v.x + this.y * v.y + this.z * v.z;
    },
    
    length : function() {
        
        return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    },
    
    lengthSq : function() {
    
        return (this.x*this.x + this.y*this.y + this.z*this.z);
    },
    
    normalize : function() {
        
        return this.divideScalar( this.length() );
    },
    
    cross : function (v) {
        
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

    setEulerFromRotationMatrix : function (m, order) {

        // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

        var te = m.elements;
        var m11 = te[0], m12 = te[4], m13 = te[8];
        var m21 = te[1], m22 = te[5], m23 = te[9];
        var m31 = te[2], m32 = te[6], m33 = te[10];

        if ( order === undefined || order === 'XYZ' ) {

            this.y = Math.asin( WebMol.Math.clamp( m13, -1, 1 ) );

            if ( Math.abs( m13 ) < 0.99999 ) {

                this.x = Math.atan2( - m23, m33 );
                this.z = Math.atan2( - m12, m11 );

            } else {

                this.x = Math.atan2( m32, m22 );
                this.z = 0;

            }
        }
        
        else {
            console.error("Error with vector's setEulerFromRotationMatrix: Unknown order: " + order);
        }
        
        return this;

    },
    
    clone : function() {
        return new WebMol.Vector3(this.x, this.y, this.z);
    }
    
};

//Matrices

//Matrix3

WebMol.Matrix3 = function(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
    
    this.elements = new Float32Array(9);
    
    this.set(
        (n11 !== undefined) ? n11 : 1, n12 || 0, n13 || 0,
        n21 || 0, (n22 !== undefined) ? n22 : 1, n23 || 0,
        n31 || 0, n32 || 0, (n33 !== undefined) ? n33 : 1
    );
    
};

WebMol.Matrix3.prototype = {
    
    constructor : WebMol.Matrix3,    
   
    set : function(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
        var te = this.elements;
        
        te[0] = n11; te[3] = n12; te[6] = n13;
        te[1] = n21; te[4] = n22; te[7] = n23;
        te[2] = n31; te[5] = n32; te[8] = n33;
        
        return this;
    },
    
    identity : function() {   
        this.set(
            1,0,0,
            0,1,0,
            0,0,1
        );
        
        return this;
    },
   
    copy : function(m) {
        var me = m.elements;
       
        this.set(
            me[0], me[3], me[6],
            me[1], me[4], me[7],
            me[2], me[5], me[8]
        );
    },
    
    multiplyScalar: function ( s ) {
        var te = this.elements;

        te[0] *= s; te[3] *= s; te[6] *= s;
        te[1] *= s; te[4] *= s; te[7] *= s;
        te[2] *= s; te[5] *= s; te[8] *= s;

        return this;
    },

    getInverse: function ( matrix, throwOnInvertible ) {
        // input: Matrix4

        var me = matrix.elements;
        var te = this.elements;

        te[ 0 ] =   me[10] * me[5] - me[6] * me[9];
        te[ 1 ] = - me[10] * me[1] + me[2] * me[9];
        te[ 2 ] =   me[6] * me[1] - me[2] * me[5];
        te[ 3 ] = - me[10] * me[4] + me[6] * me[8];
        te[ 4 ] =   me[10] * me[0] - me[2] * me[8];
        te[ 5 ] = - me[6] * me[0] + me[2] * me[4];
        te[ 6 ] =   me[9] * me[4] - me[5] * me[8];
        te[ 7 ] = - me[9] * me[0] + me[1] * me[8];
        te[ 8 ] =   me[5] * me[0] - me[1] * me[4];

        var det = me[ 0 ] * te[ 0 ] + me[ 1 ] * te[ 3 ] + me[ 2 ] * te[ 6 ];

        // no inverse

        if ( det === 0 ) {

            var msg = "Matrix3.getInverse(): can't invert matrix, determinant is 0";

            if ( throwOnInvertible || false ) {

                throw new Error( msg ); 

            } else {

                console.warn( msg );

            }

            this.identity();

            return this;

        }

        this.multiplyScalar( 1.0 / det );

        return this;
    },
    
    transpose: function () {
        var tmp, m = this.elements;

        tmp = m[1]; m[1] = m[3]; m[3] = tmp;
        tmp = m[2]; m[2] = m[6]; m[6] = tmp;
        tmp = m[5]; m[5] = m[7]; m[7] = tmp;

        return this;
    },
    
    clone: function () {
        var te = this.elements;

        return new WebMol.Matrix3(

            te[0], te[3], te[6],
            te[1], te[4], te[7],
            te[2], te[5], te[8]

        );
    }
   
};

//Matrix 4

WebMol.Matrix4 = function(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44) {

    var te = this.elements = new Float32Array( 16 );
    
    te[0] = ( n11 !== undefined ) ? n11 : 1; te[4] = n12 || 0; te[8] = n13 || 0; te[12] = n14 || 0;
    te[1] = n21 || 0; te[5] = ( n22 !== undefined ) ? n22 : 1; te[9] = n23 || 0; te[13] = n24 || 0;
    te[2] = n31 || 0; te[6] = n32 || 0; te[10] = ( n33 !== undefined ) ? n33 : 1; te[14] = n34 || 0;
    te[3] = n41 || 0; te[7] = n42 || 0; te[11] = n43 || 0; te[15] = ( n44 !== undefined ) ? n44 : 1;

};

WebMol.Matrix4.prototype = {

    constructor : WebMol.Matrix4,

    set: function ( n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44 ) {
        var te = this.elements;

        te[0] = n11; te[4] = n12; te[8] = n13; te[12] = n14;
        te[1] = n21; te[5] = n22; te[9] = n23; te[13] = n24;
        te[2] = n31; te[6] = n32; te[10] = n33; te[14] = n34;
        te[3] = n41; te[7] = n42; te[11] = n43; te[15] = n44;

        return this;
    },

    identity: function () {
        this.set(

            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1

        );

        return this;
    },

    copy: function ( m ) {
        var me = m.elements;

        this.set(

            me[0], me[4], me[8], me[12],
            me[1], me[5], me[9], me[13],
            me[2], me[6], me[10], me[14],
            me[3], me[7], me[11], me[15]

        );

        return this;
    },

    setRotationFromEuler: function ( v, order ) {

        var te = this.elements;

        var x = v.x, y = v.y, z = v.z;
        var a = Math.cos( x ), b = Math.sin( x );
        var c = Math.cos( y ), d = Math.sin( y );
        var e = Math.cos( z ), f = Math.sin( z );

        if ( order === undefined || order === 'XYZ' ) {

            var ae = a * e, af = a * f, be = b * e, bf = b * f;

            te[0] = c * e;
            te[4] = - c * f;
            te[8] = d;

            te[1] = af + be * d;
            te[5] = ae - bf * d;
            te[9] = - b * c;

            te[2] = bf - ae * d;
            te[6] = be + af * d;
            te[10] = a * c;

        } 
        
        else
            console.error("Error with matrix4 setRotationFromEuler. Order: " + order);

        return this;

    },

    setRotationFromQuaternion: function ( q ) {
        var te = this.elements;

        var x = q.x, y = q.y, z = q.z, w = q.w;
        var x2 = x + x, y2 = y + y, z2 = z + z;
        var xx = x * x2, xy = x * y2, xz = x * z2;
        var yy = y * y2, yz = y * z2, zz = z * z2;
        var wx = w * x2, wy = w * y2, wz = w * z2;

        te[0] = 1 - ( yy + zz );
        te[4] = xy - wz;
        te[8] = xz + wy;

        te[1] = xy + wz;
        te[5] = 1 - ( xx + zz );
        te[9] = yz - wx;

        te[2] = xz - wy;
        te[6] = yz + wx;
        te[10] = 1 - ( xx + yy );

        return this;
    },

    lookAt: function() {
        var x = new WebMol.Vector3();
        var y = new WebMol.Vector3();
        var z = new WebMol.Vector3();

        return function ( eye, target, up ) {

            var te = this.elements;

            z.subVectors( eye, target ).normalize();

            if ( z.length() === 0 ) {

                z.z = 1;

            }

            x.crossVectors( up, z ).normalize();

            if ( x.length() === 0 ) {

                z.x += 0.0001;
                x.crossVectors( up, z ).normalize();

            }

            y.crossVectors( z, x );


            te[0] = x.x; te[4] = y.x; te[8] = z.x;
            te[1] = x.y; te[5] = y.y; te[9] = z.y;
            te[2] = x.z; te[6] = y.z; te[10] = z.z;

            return this;
        };

    }(),

    multiplyMatrices: function ( a, b ) {
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
    
    multiplyScalar: function ( s ) {
        var te = this.elements;
    
        te[0] *= s; te[4] *= s; te[8] *= s; te[12] *= s;
        te[1] *= s; te[5] *= s; te[9] *= s; te[13] *= s;
        te[2] *= s; te[6] *= s; te[10] *= s; te[14] *= s;
        te[3] *= s; te[7] *= s; te[11] *= s; te[15] *= s;
    
        return this;
    },
    
    transpose: function () {
        var te = this.elements;
        var tmp;

        tmp = te[1]; te[1] = te[4]; te[4] = tmp;
        tmp = te[2]; te[2] = te[8]; te[8] = tmp;
        tmp = te[6]; te[6] = te[9]; te[9] = tmp;

        tmp = te[3]; te[3] = te[12]; te[12] = tmp;
        tmp = te[7]; te[7] = te[13]; te[13] = tmp;
        tmp = te[11]; te[11] = te[14]; te[14] = tmp;

        return this;
    },

    getPosition: function() {
        var v1 = new WebMol.Vector3();

        return function () {

            console.warn( 'DEPRECATED: Matrix4\'s .getPosition() has been removed. Use Vector3.getPositionFromMatrix( matrix ) instead.' );

            var te = this.elements;
            return v1.set( te[12], te[13], te[14] );
        };

    }(),

    setPosition: function ( v ) {
        var te = this.elements;

        te[12] = v.x;
        te[13] = v.y;
        te[14] = v.z;

        return this;
    },

    getInverse: function ( m, throwOnInvertible ) {
        // based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
        var te = this.elements;
        var me = m.elements;

        var n11 = me[0], n12 = me[4], n13 = me[8], n14 = me[12];
        var n21 = me[1], n22 = me[5], n23 = me[9], n24 = me[13];
        var n31 = me[2], n32 = me[6], n33 = me[10], n34 = me[14];
        var n41 = me[3], n42 = me[7], n43 = me[11], n44 = me[15];

        te[0] = n23*n34*n42 - n24*n33*n42 + n24*n32*n43 - n22*n34*n43 - n23*n32*n44 + n22*n33*n44;
        te[4] = n14*n33*n42 - n13*n34*n42 - n14*n32*n43 + n12*n34*n43 + n13*n32*n44 - n12*n33*n44;
        te[8] = n13*n24*n42 - n14*n23*n42 + n14*n22*n43 - n12*n24*n43 - n13*n22*n44 + n12*n23*n44;
        te[12] = n14*n23*n32 - n13*n24*n32 - n14*n22*n33 + n12*n24*n33 + n13*n22*n34 - n12*n23*n34;
        te[1] = n24*n33*n41 - n23*n34*n41 - n24*n31*n43 + n21*n34*n43 + n23*n31*n44 - n21*n33*n44;
        te[5] = n13*n34*n41 - n14*n33*n41 + n14*n31*n43 - n11*n34*n43 - n13*n31*n44 + n11*n33*n44;
        te[9] = n14*n23*n41 - n13*n24*n41 - n14*n21*n43 + n11*n24*n43 + n13*n21*n44 - n11*n23*n44;
        te[13] = n13*n24*n31 - n14*n23*n31 + n14*n21*n33 - n11*n24*n33 - n13*n21*n34 + n11*n23*n34;
        te[2] = n22*n34*n41 - n24*n32*n41 + n24*n31*n42 - n21*n34*n42 - n22*n31*n44 + n21*n32*n44;
        te[6] = n14*n32*n41 - n12*n34*n41 - n14*n31*n42 + n11*n34*n42 + n12*n31*n44 - n11*n32*n44;
        te[10] = n12*n24*n41 - n14*n22*n41 + n14*n21*n42 - n11*n24*n42 - n12*n21*n44 + n11*n22*n44;
        te[14] = n14*n22*n31 - n12*n24*n31 - n14*n21*n32 + n11*n24*n32 + n12*n21*n34 - n11*n22*n34;
        te[3] = n23*n32*n41 - n22*n33*n41 - n23*n31*n42 + n21*n33*n42 + n22*n31*n43 - n21*n32*n43;
        te[7] = n12*n33*n41 - n13*n32*n41 + n13*n31*n42 - n11*n33*n42 - n12*n31*n43 + n11*n32*n43;
        te[11] = n13*n22*n41 - n12*n23*n41 - n13*n21*n42 + n11*n23*n42 + n12*n21*n43 - n11*n22*n43;
        te[15] = n12*n23*n31 - n13*n22*n31 + n13*n21*n32 - n11*n23*n32 - n12*n21*n33 + n11*n22*n33;

        var det = me[ 0 ] * te[ 0 ] + me[ 1 ] * te[ 4 ] + me[ 2 ] * te[ 8 ] + me[ 3 ] * te[ 12 ];

        if ( det === 0 ) {

            var msg = "Matrix4.getInverse(): can't invert matrix, determinant is 0";

            if ( throwOnInvertible || false ) {

                throw new Error( msg ); 

            } else {

                console.warn( msg );

            }

            this.identity();

            return this;
        }

        this.multiplyScalar( 1 / det );

        return this;
    },

    compose: function() {
        var mRotation = new WebMol.Matrix4(),
            mScale = new WebMol.Matrix4();
        
        return function ( translation, rotation, scale ) {

            var te = this.elements;

            mRotation.identity();
            mRotation.setRotationFromQuaternion( rotation );

            mScale.makeScale( scale.x, scale.y, scale.z );

            this.multiplyMatrices( mRotation, mScale );

            te[12] = translation.x;
            te[13] = translation.y;
            te[14] = translation.z;

            return this;

        };
    }(),

    decompose: function() {
        var x = new WebMol.Vector3(),
            y = new WebMol.Vector3(),
            z = new WebMol.Vector3(),
            matrix = new WebMol.Matrix4();

        return function ( translation, rotation, scale ) {

            var te = this.elements;

            // grab the axis vectors
            x.set( te[0], te[1], te[2] );
            y.set( te[4], te[5], te[6] );
            z.set( te[8], te[9], te[10] );

            translation = ( translation instanceof WebMol.Vector3 ) ? translation : new WebMol.Vector3();
            rotation = ( rotation instanceof WebMol.Quaternion ) ? rotation : new WebMol.Quaternion();
            scale = ( scale instanceof Webmol.Vector3 ) ? scale : new WebMol.Vector3();

            scale.x = x.length();
            scale.y = y.length();
            scale.z = z.length();

            translation.x = te[12];
            translation.y = te[13];
            translation.z = te[14];

            // scale the rotation part

            matrix.copy( this );

            matrix.elements[0] /= scale.x;
            matrix.elements[1] /= scale.x;
            matrix.elements[2] /= scale.x;

            matrix.elements[4] /= scale.y;
            matrix.elements[5] /= scale.y;
            matrix.elements[6] /= scale.y;

            matrix.elements[8] /= scale.z;
            matrix.elements[9] /= scale.z;
            matrix.elements[10] /= scale.z;

            rotation.setFromRotationMatrix( matrix );

            return [ translation, rotation, scale ];

        };
    }(),

    scale: function ( v ) {
        var te = this.elements;
        var x = v.x, y = v.y, z = v.z;

        te[0] *= x; te[4] *= y; te[8] *= z;
        te[1] *= x; te[5] *= y; te[9] *= z;
        te[2] *= x; te[6] *= y; te[10] *= z;
        te[3] *= x; te[7] *= y; te[11] *= z;

        return this;
    },
    
    getMaxScaleOnAxis : function() {
        
        var te = this.elements;
        
        var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
        var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
        var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];
        
        return Math.sqrt(Math.max(scaleXSq, Math.max(scaleYSq, scaleZSq)));
        
    },

    makeFrustum: function ( left, right, bottom, top, near, far ) {
        var te = this.elements;
        var x = 2 * near / ( right - left );
        var y = 2 * near / ( top - bottom );

        var a = ( right + left ) / ( right - left );
        var b = ( top + bottom ) / ( top - bottom );
        var c = - ( far + near ) / ( far - near );
        var d = - 2 * far * near / ( far - near );

        te[0] = x;  te[4] = 0;  te[8] = a;  te[12] = 0;
        te[1] = 0;  te[5] = y;  te[9] = b;  te[13] = 0;
        te[2] = 0;  te[6] = 0;  te[10] = c; te[14] = d;
        te[3] = 0;  te[7] = 0;  te[11] = - 1;   te[15] = 0;

        return this;
    },

    makePerspective: function ( fov, aspect, near, far ) {
        var ymax = near * Math.tan( WebMol.Math.degToRad( fov * 0.5 ) );
        var ymin = - ymax;
        var xmin = ymin * aspect;
        var xmax = ymax * aspect;

        return this.makeFrustum( xmin, xmax, ymin, ymax, near, far );
    },
    

    clone: function () {
        var te = this.elements;

        return new WebMol.Matrix4(

            te[0], te[4], te[8], te[12],
            te[1], te[5], te[9], te[13],
            te[2], te[6], te[10], te[14],
            te[3], te[7], te[11], te[15]

        );
    }
    
};

WebMol.Ray = function(origin, direction) {
    
    this.origin = (origin !== undefined) ? 
        origin : new WebMol.Vector3();
        
    this.direction = (direction !== undefined) ?
        direction : new WebMol.Vector3();
      
};

//TODO: Remove methods we don't need (intersectPlane ??)
WebMol.Ray.prototype = {
    
    constructor : WebMol.Ray,
     
    set : function(origin, direction){
        
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
        
        var result = optionalTarget || new WebMol.Vector3();
        
        return result.copy(this.direction).multiplyScalar(t).add(this.origin);
        
    },
    
    recast : function() {
        
        var v1 = new WebMol.Vector3();
        
        return function(t) {
            this.origin.copy(this.at(t, v1));
            
            return this;
        };
        
    }(),
    
    closestPointToPoint : function(point, optionalTarget) {
        
        var result = optionalTarget || new WebMol.Vector3();
        result.subVectors(point, this.origin);
        var directionDistance = result.dot(this.direction);
        
        //returns a point on this ray
        return result.copy(this.direction).multiplyScalar(directionDistance).add(this.origin);
        
    },
    
    distanceToPoint : function() {
        
        var v1 = new WebMol.Vector3();
        
        return function(point) {
            var directionDistance = v1.subVectors(point, this.origin).dot(this.direction);
            v1.copy(this.direction).multiplyScalar(directionDistance).add(this.origin);
            
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
        
        //plane and ray are not perpendicular
        if (denominator !== 0) 
            return true;
        
        if (plane.distanceToPoint(this.origin) === 0) 
            return true;
        
        return false;
        
    },
    
    distanceToPlane : function(plane) {
       
       var denominator = plane.normal.dot(this.direction);
       if (denominator === 0) {
           
           //line is coplanar
       if (plane.distanceToPoint(this.origin) === 0)
           return 0;
       
       //ray is parallel
           return undefined;
       }
       
       var t = - (this.origin.dot(plane.normal) + plane.constant) / denominator;
       
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
       
       return ray.origin.equals(this.origin) && ray.direction.equals(this.direction);
       
    },
    
    clone : function() {
    
       return new WebMol.Ray().copy(this);
    
    }
 
     
};

//Intersection sphere and box shapes.  


//Intersection sphere for sphere, stick render
WebMol.Sphere = function(center, radius) {

    this.center = (center !== undefined) ? 
        center : new WebMol.Vector3();
        
    this.radius = (radius !== undefined) ?
        radius : 0;
        
};

WebMol.Sphere.prototype = {
    
    constructor : WebMol.Sphere,
    
    set : function(center, radius) {
        
        this.center.copy(center);
        this.radius = radius;
        
        return this;
        
    },
    
    copy : function(sphere) {
        
        this.center.copy(sphere.center);
        this.radius = sphere.radius;
        
        return this;
        
    },
    
    applyMatrix4 : function(matrix) {
        
        this.center.applyMatrix4(matrix);
        this.radius = this.radius * matrix.getMaxScaleOnAxis();
        
        return this;
        
    },
    
    translate : function(offset) {
        
        this.center.add(offset);
        
        return this;
        
    },
    
    equals : function(sphere) {
        
        return sphere.center.equals(this.center) && (sphere.radius === this.radius);
        
    },
       
    clone : function() {
        
        return new WebMol.Sphere().copy(this);
        
    }

};


//Bounding cylinder for stick render  
WebMol.Cylinder = function(c1, c2, radius) {

    this.c1 = (c1 !== undefined) ?
        c1 : new WebMol.Vector3();

    this.c2 = (c2 !== undefined) ?
        c2 : new WebMol.Vector3();
        
    this.direction = new WebMol.Vector3().subVectors(this.c2, this.c1).normalize();

    this.radius = (radius !== undefined) ?
        radius : 0;
    
};

WebMol.Cylinder.prototype = {

    constructor : WebMol.Cylinder,

    copy : function(cylinder) {

        this.c1.copy(cylinder.c1);
        this.c2.copy(cylinder.c2);
        this.direction.copy(cylinder.direction);
        this.radius = cylinder.radius;

        return this;

    },
    
    lengthSq : function() {
    
        var vector = new WebMol.Vector3();
        
        return function(){
            return vector.subVectors(this.c2, this.c1).lengthSq();
        };
        
    }(),

    applyMatrix4 : function(matrix) {
        
        this.direction.add(this.c1).applyMatrix4(matrix);
        this.c1.applyMatrix4(matrix);
        this.c2.applyMatrix4(matrix);
        this.direction.sub(this.c1).normalize();
        this.radius = this.radius * matrix.getMaxScaleOnAxis();

        return this;

    }

};


//plane specified by three points
WebMol.Triangle = function(a, b, c){
   
    this.a = (a !== undefined) ?
        a : new WebMol.Vector3();

    this.b = (b !== undefined) ?
        b : new WebMol.Vector3();
    
    this.c = (c !== undefined) ?
        c : new WebMol.Vector3();   
  
};

WebMol.Triangle.prototype = {

    constructor : WebMol.Triangle,
    
    copy : function(triangle) {
        
        this.a.copy(triangle.a);
        this.b.copy(triangle.b);
        this.c.copy(triangle.c);
        
        return this;
        
    },
    
    applyMatrix4 : function(matrix) {
        
        this.a.applyMatrix4(matrix);
        this.b.applyMatrix4(matrix);
        this.c.applyMatrix4(matrix);
        
        return this;
        
    },
    
    getNormal : function() {
        
        var v1 = new WebMol.Vector3();
        
        return function() {
            
            var norm = this.a.clone();
            norm.sub(this.b);
            v1.subVectors(this.c, this.b);
            
            norm.cross(v1);
            norm.normalize();
            
            return norm;
            
        };
        
    }()

};


/* core Object3D
 * Base class for Scene, Camera, Geometry
 * Geometry class
 */

var WebMol = WebMol || {};

//Event Handling
WebMol.EventDispatcher = function() {
  
    var listeners = {};
    
    this.addEventListener = function(type, listener) {
        if (listeners[type] === undefined)
            listeners[type] = [];
        
        if (listeners[type].indexOf(listener) === -1)
            listeners[type].push(listener);
    };  
    
    this.removeEventListener = function(type, listener) {
        
        var index = listeners[type].indexOf(listener);
        
        if (index !== -1)
            listeners[type].splice(index, 1);
              
    };
    
    this.dispatchEvent = function(event) {
        
        var listenerArray = listeners[event.type];
        
        if (listenerArray !== undefined) {
            event.target = this;
            
            for (var i = 0, l = listenerArray.length; i < l; i++)
                listenerArray[i].call(this, event);
                
        }
            
    };
    
};


//Object3D base constructor function
WebMol.Object3D = function() {
    
    this.id = WebMol.Object3DIDCount++;
    
    this.name = "";
    
    this.parent = undefined;
    this.children = [];
    
    this.position = new WebMol.Vector3();
    this.rotation = new WebMol.Vector3();
    this.matrix = new WebMol.Matrix4();
    this.matrixWorld = new WebMol.Matrix4();
    this.quaternion = new WebMol.Quaternion();
    this.eulerOrder = 'XYZ';
    
    this.up = new WebMol.Vector3(0, 1, 0);
    this.scale = new WebMol.Vector3(1, 1, 1);
    
    this.matrixAutoUpdate = true;
    this.matrixWorldNeedsUpdate = true;
    this.rotationAutoUpdate = true;
    this.useQuaternion = false;
    
    this.visible = true;
    
};

WebMol.Object3D.prototype = {
    
    constructor : WebMol.Object3D,
    
    lookAt : function(vector) {
        
        this.matrix.lookAt(vector, this.position, this.up);
        
        if (this.rotationAutoUpdate) {
            
            if (this.useQuaternion === true) 
                this.quaternion.copy(this.matrix.decompose()[1]);
            else
                this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
        }
    },
    
    //add child object
    add : function(object) {
        if (object === this){
            console.error("Can't add WebMol.Object3D to itself");
            return;
        }
        
        object.parent = this;
        this.children.push(object);
        
        //add to the scene (i.e. follow up this instance's parents until reach the top)
        
        var scene = this;
        
        while (scene.parent !== undefined)
            scene = scene.parent;
            
        if (scene !== undefined && scene instanceof WebMol.Scene) 
            scene.__addObject(object);
        
    },
    
    remove : function(object) {
        
        var index = this.children.indexOf(object);
        
        if (index !== -1) {
            
            object.parent = undefined;
            this.children.splice(index, 1);
            
            //Remove from scene
            
            var scene = this;
            
            while (scene.parent !== undefined)
                scene = scene.parent;
                
            if (scene !== undefined && scene instanceof WebMol.Scene)
                scene.__removeObject(object);
                
        }
    },
    
    updateMatrix : function() {
        
        this.matrix.setPosition(this.position);
        
        if (this.useQuaternion === false) 
            this.matrix.setRotationFromEuler(this.rotation, this.eulerOrder);
        else
            this.matrix.setRotationFromQuaternion(this.quaternion);
        
        //TODO: Do I need this??
        if (this.scale.x !== 1 || this.scale.y !== 1 || this.scale.z !== 1)
            this.matrix.scale(this.scale);
            
        this.matrixWorldNeedsUpdate = true;
        
    },
    
    updateMatrixWorld : function(force) {
        
        if (this.matrixAutoUpdate === true) 
            this.updateMatrix();
        
        if (this.matrixWorldNeedsUpdate === true || force === true) {
            
            if (this.parent === undefined)
                this.matrixWorld.copy(this.matrix);
            else
                this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
                
        }
        
        this.matrixWorldNeedsUpdate = false;
        
        //Update matrices of all children
        for (var i in this.children) {
            this.children[i].updateMatrixWorld(true);
        }
    },
    
    clone : function(object) {
        
        if (object === undefined)
            object = new WebMol.Object3D();
            
        object.name = this.name;
        
        object.up.copy(this.up);
        object.position.copy(this.position);
        object.rotation.copy(this.rotation);
        object.eulerOrder = this.eulerOrder;
        object.scale.copy(this.scale);
        
        object.rotationAutoUpdate = this.rotationAutoUpdate;
        object.matrix.copy(this.matrix);
        object.matrixWorld.copy(this.matrixWorld);
        object.quaternion.copy(this.quaternion);
        
        object.useQuaternion = this.useQuaternion;
        
        object.visible = this.visible;
        
        for (var i in this.children) {
            var child = this.children[i];
            object.add(child.clone());
        }
        
        return object;
        
    }
    
};

WebMol.Object3DIDCount = 0;

//Geometry class
//TODO: What can I remove - how can I optimize ?
WebMol.Geometry = (function() {

    //return truncated typed array, including its buffer
    // type == 0 => Uint16Array; type == 1 => Float32Array
    //TODO: Should integrate this directly into geometryGroup's truncateArrayBuffers method
    var truncateArrayBuffer = function(arr, type, end) {
        
        if (arr === null || arr === undefined) {
            return (type === 0) ? new Uint16Array() : new Float32Array();
        }
        
        if (type === 0)
            return new Uint16Array(arr.buffer.slice(arr.byteOffset, end*2));
        else if (type === 1) 
            return new Float32Array(arr.buffer.slice(arr.byteOffset, end*4));
    };
    
    
    var geometryGroup = function(id) {
        this.id = id || 0;
        this.__vertexArray = null;
        this.__colorArray = null;
        this.__normalArray = null;
        this.__faceArray = null;
        this.__lineArray = null;
        this.vertices = 0;
        this.faceidx = 0;
        this.lineidx = 0;
    };
    
    geometryGroup.prototype.getCentroid = function() {
        
        var centroid = new WebMol.Vector3();
        var offset, x, y, z;
        
        for (var i = 0; i < this.vertices; ++i) {
            offset = i*3;
            
            x = this.__vertexArray[offset]; y = this.__vertexArray[offset+1]; z = this.__vertexArray[offset+2];
            
            centroid.x += x; centroid.y += y; centroid.z += z;
        }
        
        //divideScalar checks for 0 denom
        centroid.divideScalar(this.vertices);
        
        return centroid;
    };
    
    //setup normals - vertex and face array must exist
    geometryGroup.prototype.setNormals = function() {        
        
        var faces = this.__faceArray;
        var verts = this.__vertexArray;
        var norms = this.__normalArray;
        
        if (! this.vertices || ! this.faceidx) 
            return;
        
        //vertex indices
        var a, b, c, d,
        //and actual vertices
        vA, vB, vC, norm;
            
        for (var i = 0; i < faces.length / 3; ++i) {
            
            a = faces[i * 3] * 3;
            b = faces[i * 3 + 1] * 3;
            c = faces[i * 3 + 2] * 3;
            
            vA = new WebMol.Vector3(verts[a], verts[a+1], verts[a+2]);
            vB = new WebMol.Vector3(verts[b], verts[b+1], verts[b+2]);
            vC = new WebMol.Vector3(verts[c], verts[c+1], verts[c+2]);
            
            vA.subVectors(vA, vB);
            vC.subVectors(vC, vB);
            vC.cross(vA);
            
            //face normal
            norm = vC;
            norm.normalize();
            
            norms[a] += norm.x; norms[b] += norm.x; norms[c] += norm.x;
            norms[a + 1] += norm.y; norms[b + 1] += norm.y; norms[c + 1] += norm.y;
            norms[a + 2] += norm.z; norms[b + 2] += norm.z; norms[c + 2] += norm.z;
            
        }             
                
    };
    
    //sets line index array from face arr
    //Note - assumes all faces are triangles (i.e. there will
    //be an extra diagonal for four-sided faces - user should 
    //specify linearr for custom shape generation to show wireframe squares
    //as rectangles rather than two triangles)
    geometryGroup.prototype.setLineIndices = function() {
        
        if (! this.faceidx)
            return;
                    
        var faceArr = this.__faceArray, lineArr = this.__lineArray = new Uint16Array(this.faceidx*2);      
        this.lineidx = this.faceidx*2;         
        var faceoffset;
            
        for (var i = 0; i < this.faceidx / 3; ++i) {
            
            faceoffset = i*3; lineoffset = faceoffset*2;          
            var a = faceArr[faceoffset], b = faceArr[faceoffset+1], c = faceArr[faceoffset+2];
            
            lineArr[lineoffset] = a; lineArr[lineoffset+1] = b;
            lineArr[lineoffset+2] = a; lineArr[lineoffset+3] = c;
            lineArr[lineoffset+4] = b; lineArr[lineoffset+5] = c;
            
        }
    };
    
    geometryGroup.prototype.truncateArrayBuffers = function(mesh) {
        
        mesh = (mesh === true) ? true : false;
        
        var vertexArr = this.__vertexArray,
            colorArr = this.__colorArray,
            normalArr = this.__normalArray,
            faceArr = this.__faceArray,
            lineArr = this.__lineArray;
                       
        this.__vertexArray = truncateArrayBuffer(vertexArr, 1, this.vertices*3);
        this.__colorArray = truncateArrayBuffer(colorArr, 1, this.vertices*3);
        
        if (mesh) {
            this.__normalArray = truncateArrayBuffer(normalArr, 1, this.vertices*3);
            this.__faceArray = truncateArrayBuffer(faceArr, 0, this.faceidx);
            this.__lineArray = truncateArrayBuffer(lineArr, 0, this.lineidx);
        }
        else {
            this.__normalArray = truncateArrayBuffer(normalArr, 1, 0);
            this.__faceArray = truncateArrayBuffer(faceArr, 0, 0);
            this.__lineArray = truncateArrayBuffer(lineArr, 0, 0);            
        }
        
        this.__inittedArrays = true;        
        
    };
    
    var addGroup = function(geo) {
        var ret = new geometryGroup(geo.geometryGroups.length);
        geo.geometryGroups.push(ret);
        geo.groups = geo.geometryGroups.length;
        
        ret.__vertexArray = new Float32Array(65535*3);
        ret.__colorArray = new Float32Array(65535*3);
        
        //TODO: instantiating uint arrays according to max number of vertices
        // is dangerous, since there exists the possibility that there will be 
        // more face or line indices than vertex points - but so far that doesn't
        // seem to be the case for any of the renders 
        if (geo.mesh) {
            ret.__normalArray = new Float32Array(65535*3);
            ret.__faceArray = new Uint16Array(65535*6);
            ret.__lineArray = new Uint16Array(65535*6);
        }
        
        
        return ret;
    };
        
    var Geometry = function(mesh) {
        
        WebMol.EventDispatcher.call(this);
        
        this.id = WebMol.GeometryIDCount++;
    
        this.name = '';
    
        this.hasTangents = false;
    
        this.dynamic = true; // the intermediate typed arrays will be deleted when set to false
        this.mesh = (mesh === true) ? true : false; // Does this geometry represent a mesh (i.e. do we need Face/Line index buffers?)
        // update flags
    
        this.verticesNeedUpdate = false;
        this.elementsNeedUpdate = false;
        this.normalsNeedUpdate = false;
        this.colorsNeedUpdate = false;
    
        this.buffersNeedUpdate = false;
        
        this.geometryGroups = [];
        this.groups = 0;
        
    };
    
    
    Geometry.prototype = {
        
        constructor : Geometry,

        //Get geometry group to accomodate addVertices new vertices - create 
        // new group if necessary       
        updateGeoGroup : function(addVertices) {
        
            addVertices = addVertices || 0;
            
            var retGroup = this.groups > 0 ? this.geometryGroups[ this.groups - 1 ] : null;
            
            if (!retGroup || retGroup.vertices + addVertices > 65535) 
                retGroup = addGroup(this);
                
            return retGroup;
            
        },
        
        addGeoGroup : function() {
            return addGroup(this);  
        },
        
        setUpNormals : function(three) {
            
            three = three || false;
            
            for ( var g in this.geometryGroups ) {
            
                var geoGroup = this.geometryGroups[g];            
                
                geoGroup.setNormals(three);
                
            }  
                      
        },
        
        setUpWireframe : function() {
            for (var g in this.geometryGroups ) {
                var geoGroup = this.geometryGroups[g];
                
                geoGroup.setLineIndices();
            }
        },
        
        //After vertices, colors, etc are collected in regular or typed arrays,
        // either create typed arrays from regular arrays if they don't already exist,
        // or shorten last typed array
        initTypedArrays : function() {
                
            for (var g in this.geometryGroups) {
                
                var group = this.geometryGroups[g];
                
                if (group.__inittedArrays === true)
                    continue;
                
                group.truncateArrayBuffers(this.mesh);
            }
            
        
        },
        
        dispose : function() {
            this.dispatchEvent( {type : 'dispose'} );
        }
    };

    
    return Geometry;
    
})();

Object.defineProperty(WebMol.Geometry.prototype, "vertices", {
    
    get : function() {
        var vertices = 0;
        for (var g in this.geometryGroups)
            vertices += this.geometryGroups[g].vertices;
            
        return vertices;
    } 
        
});

WebMol.GeometryIDCount = 0;


//Raycaster

WebMol.Raycaster = (function() {
    
    var Raycaster = function(origin, direction, far, near) {
        
        this.ray = new WebMol.Ray(origin, direction);
        
        if (this.ray.direction.lengthSq() > 0) 
            this.ray.direction.normalize();
        
        this.near = near || 0;
        this.far = far || Infinity;
    
    };
    
    var sphere = new WebMol.Sphere();
    var cylinder = new WebMol.Cylinder();
    var triangle = new WebMol.Triangle();
    var w_0 = new WebMol.Vector3(); // for cylinders, cylinder.c1 - ray.origin
    var v1 = new WebMol.Vector3(); // all purpose local vector
    var v2 = new WebMol.Vector3();
    var v3 = new WebMol.Vector3();
    //var facePlane = new WebMol.Plane();
    var localRay = new WebMol.Ray();
    var intersectPoint = new WebMol.Vector3();
    var matrixPosition = new WebMol.Vector3();
    
    var inverseMatrix = new WebMol.Matrix4();
        
    var descSort = function(a, b) {
        return a.distance - b.distance;
    };

    // [-1, 1]
    var clamp = function(x) {
        return Math.min(Math.max(x, -1), 1);
    };
    
    //object is a Sphere or (Bounding) Box
    var intersectObject = function(group, clickable, raycaster, intersects) {
        
        matrixPosition.getPositionFromMatrix(group.matrixWorld);
        
        if ((clickable.clickable !== true) || (clickable.intersectionShape === undefined))
            return intersects;
        
        var intersectionShape = clickable.intersectionShape;
        var precision = raycaster.linePrecision;
        precision *= group.matrixWorld.getMaxScaleOnAxis();
        var precisionSq = precision*precision;

        //Check for intersection with clickable's bounding sphere, if it exists
        if (clickable.boundingSphere !== undefined && clickable.boundingSphere instanceof WebMol.Sphere) {
            sphere.copy(clickable.boundingSphere);
            sphere.applyMatrix4(group.matrixWorld);
            
            if (!raycaster.ray.isIntersectionSphere(sphere)) {
                return intersects;
            }
        }
        
        //Iterate through intersection objects
        var i, il,
            norm, normProj, cylProj, rayProj,
            distance, closestDistSq, denom, discriminant,
            s, t, s_c, t_c;
        //triangle faces
        for (i = 0, il = intersectionShape.triangle.length; i < il; i++) {
            
            if (intersectionShape.triangle[i] instanceof WebMol.Triangle) {
                
                triangle.copy(intersectionShape.triangle[i]);
                triangle.applyMatrix4(group.matrixWorld);
                
                norm = triangle.getNormal();
                
                normProj = raycaster.ray.direction.dot(norm);
                
                //face culling
                if (normProj >= 0)
                    continue;
                
                w_0.subVectors(triangle.a, raycaster.ray.origin);
                
                distance = (norm.dot(w_0)) / normProj;
                
                if (distance < 0)
                    continue;
                    
                //intersects with plane, check if P inside triangle
                v1.copy(raycaster.ray.direction).multiplyScalar(distance).add(raycaster.ray.origin);
                v1.sub(triangle.a); // from pt a to intersection point P
                
                v2.copy(triangle.b).sub(triangle.a); // from pt a to b
                v3.copy(triangle.c).sub(triangle.a); // from pt a to c
                var b_dot_c = v2.dot(v3);
                var b_sq = v2.lengthSq();
                var c_sq = v3.lengthSq();
                
                // P = A + s(v2) + t(v3), inside trianle if 0 <= s, t <=1  and (s + t) <=0
                
                t = ( b_sq*v1.dot(v3) - b_dot_c*v1.dot(v2) ) / ( b_sq*c_sq - b_dot_c*b_dot_c );
                
                if (t < 0 || t > 1)
                    continue;
                
                s = ( v1.dot(v2) - t*b_dot_c ) / b_sq;
                
                if ( (s < 0 || s > 1) || s + t > 1)
                    continue;
                    
                else
                    intersects.push({clickable : clickable,
                                     distance : distance});  
            }
        }
        
        //cylinders
        for (i = 0, il = intersectionShape.cylinder.length; i < il; i++) {
            
            if (intersectionShape.cylinder[i] instanceof WebMol.Cylinder){
                
                cylinder.copy(intersectionShape.cylinder[i]);
                cylinder.applyMatrix4(group.matrixWorld);
                
                w_0.subVectors(cylinder.c1, raycaster.ray.origin); 
                
                cylProj = w_0.dot(cylinder.direction); // Dela
                rayProj = w_0.dot(raycaster.ray.direction); // Epsilon
                
                normProj = clamp(raycaster.ray.direction.dot(cylinder.direction)); // Beta
                
                denom = 1 - normProj*normProj;
                
                if (denom === 0.0)
                    continue;
                
                s_c = (normProj*rayProj - cylProj) / denom;
                t_c = (rayProj - normProj*cylProj) / denom;
                
                v1.copy(cylinder.direction).multiplyScalar(s_c).add(cylinder.c1);  // Q_c
                v2.copy(raycaster.ray.direction).multiplyScalar(t_c).add(raycaster.ray.origin); // P_c
                
                closestDistSq = v3.subVectors(v1, v2).lengthSq();
                var radiusSq = cylinder.radius*cylinder.radius;
                
                //Smoothing?
                //if (closestDistSq > radiusSq) radiusSq += precisionSq;
                
                // closest distance between ray and cylinder axis not greater than cylinder radius;
                // might intersect this cylinder between atom and bond midpoint
                if (closestDistSq <= radiusSq){

                    //Find points where ray intersects sides of cylinder
                    discriminant = (normProj*cylProj - rayProj)*(normProj*cylProj - rayProj) - 
                            denom*(w_0.lengthSq() - cylProj*cylProj - radiusSq);
                    
                    // ray tangent to cylinder?
                    if (discriminant <= 0)
                        t = distance = Math.sqrt(closestDistSq);
                    else
                        t = distance = ( (rayProj - normProj*cylProj) - Math.sqrt(discriminant) ) / denom; 
                    
                    //find closest intersection point; make sure it's between atom's position and cylinder midpoint
                    
                    s = normProj*t - cylProj;
                    
                    //does not intersect cylinder between atom and midpoint,
                    // or intersects cylinder behind camera
                    if (s < 0 || s*s > cylinder.lengthSq() || t < 0)
                        continue;
                    
                    else
                        intersects.push({clickable : clickable,
                                         distance : distance});
                    
                }
                    
                
            }
            
        }
         
        //lines
        for (i = 0, il = intersectionShape.line.length; i < il; i += 2) {
            
            v1.copy(intersectionShape.line[i]);
            v1.applyMatrix4(group.matrixWorld);
            v2.copy(intersectionShape.line[i+1]);
            v2.applyMatrix4(group.matrixWorld);
            
            v3.subVectors(v2, v1);
            var bondLengthSq = v3.lengthSq();
            v3.normalize();
            
            w_0.subVectors(v1, raycaster.ray.origin);
            
            lineProj = w_0.dot(v3);
            rayProj = w_0.dot(raycaster.ray.direction);
            
            normProj = clamp(raycaster.ray.direction.dot(v3));
            
            denom = 1 - normProj*normProj;
            
            if (denom === 0.0)
                continue;
            
            s_c = (normProj*rayProj - lineProj) / denom;
            t_c = (rayProj - normProj*lineProj) / denom;
            
            v1.add(v3.multiplyScalar(s_c)); // Q_c
            v2.copy(raycaster.ray.direction).multiplyScalar(t_c).add(raycaster.ray.origin); // P_c
            
            closestDistSq = v3.subVectors(v2, v1).lengthSq();
            
            if (closestDistSq < precisionSq && s_c*s_c < bondLengthSq)
                intersects.push({clickable : clickable,
                                 distance : t_c
                                });
            
        }

        for (i = 0, il = intersectionShape.sphere.length; i < il; i++) {
            //sphere
            if (intersectionShape.sphere[i] instanceof WebMol.Sphere) {
                
                sphere.copy(intersectionShape.sphere[i]);
                sphere.applyMatrix4(group.matrixWorld);
                
                if (raycaster.ray.isIntersectionSphere(sphere)) {
                    
                    v1.subVectors(sphere.center, raycaster.ray.origin);
                    
                    //distance from ray origin to point on the ray normal to sphere's center
                    //must be less than sphere's radius (since ray intersects sphere)
                    var distanceToCenter = v1.dot(raycaster.ray.direction);
                    
                    discriminant = distanceToCenter*distanceToCenter - (v1.lengthSq() - sphere.radius*sphere.radius);
                    
                    //Don't select if sphere center behind camera
                    if (distanceToCenter < 0) 
                        return intersects;
                    
                    //ray tangent to sphere?
                    if (discriminant <= 0)
                        distance = distanceToCenter;
                    
                    //This is reversed if sphere is closer than ray origin.  Do we have 
                    //to worry about handling that case?
                    else 
                        distance = distanceToCenter - Math.sqrt(discriminant);
    
                    intersects.push({clickable : clickable, 
                                     distance : distance});
                    return intersects;
                }
            }        
       }
        
    };   
       
    Raycaster.prototype.precision = 0.0001;
    Raycaster.prototype.linePrecision = 0.2;
    
    Raycaster.prototype.set = function(origin, direction) {
        
        this.ray.set(origin, direction);
          
    };
    
    Raycaster.prototype.intersectObjects = function(group, objects) {
        
        var intersects = [];
        
        for (var i = 0, l = objects.length; i < l; i++)            
            intersectObject(group, objects[i], this, intersects);
            
        intersects.sort(descSort);
        
        return intersects;
        
    };
    
    return Raycaster;
    
})();


//WebMol Projecion 
//TODO: can probably strip this down a lot (only used for selection handling)
WebMol.Projector = function () {

    var _viewMatrix = new WebMol.Matrix4(),
    _viewProjectionMatrix = new WebMol.Matrix4();

    this.projectVector = function ( vector, camera ) {

        camera.matrixWorldInverse.getInverse( camera.matrixWorld );

        _viewProjectionMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

        return vector.applyProjection( _viewProjectionMatrix );

    };

    this.unprojectVector = function ( vector, camera ) {

        camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);

        _viewProjectionMatrix.multiplyMatrices(camera.matrixWorld, camera.projectionMatrixInverse);

        return vector.applyProjection( _viewProjectionMatrix );

    };

};/*
 * Simplified Perspective Camera
 */


WebMol.Camera = function(fov, aspect, near, far) {
    
    WebMol.Object3D.call(this);
    
    this.fov = fov !== undefined ? fov : 50;
    this.aspect = aspect !== undefined ? aspect : 1;
    this.near = near !== undefined ? near : 0.1;
    this.far = far !== undefined ? far : 2000;

    this.projectionMatrix = new WebMol.Matrix4();
    this.projectionMatrixInverse = new WebMol.Matrix4();
    this.matrixWorldInverse = new WebMol.Matrix4();
    
    this.updateProjectionMatrix();
        
};

//Inherit Object3D's prototyped methods
WebMol.Camera.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Camera.prototype.lookAt = function(vector){
    
    //Why is the parameter order switched (compared to Object3D)?
    this.matrix.lookAt(this.position, vector, this.up);
    
    if (this.rotationAutoUpdate) {    
        
        if (this.useQuaternion === false) 
            this.rotation.setEulerFromRotationMatrix( this.matrix, this.eulerOrder );
        else
            this.quaternion.copy( this.matrix.decompose()[ 1 ] );    
            
    }
    
};

WebMol.Camera.prototype.updateProjectionMatrix = function () {

    this.projectionMatrix.makePerspective( this.fov, this.aspect, this.near, this.far );

};


//Render plugins go here

/**
 * Sprite render plugin
 */

WebMol.SpritePlugin = function () {

    var _gl, _renderer, _precision, _sprite = {};

    this.init = function ( renderer ) {

        _gl = renderer.context;
        _renderer = renderer;

        _precision = renderer.getPrecision();

        _sprite.vertices = new Float32Array( 8 + 8 );
        _sprite.faces    = new Uint16Array( 6 );

        var i = 0;

        _sprite.vertices[ i++ ] = -1; _sprite.vertices[ i++ ] = -1; // vertex 0
        _sprite.vertices[ i++ ] = 0;  _sprite.vertices[ i++ ] = 0;  // uv 0

        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = -1; // vertex 1
        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = 0;  // uv 1

        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = 1;  // vertex 2
        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = 1;  // uv 2

        _sprite.vertices[ i++ ] = -1; _sprite.vertices[ i++ ] = 1;  // vertex 3
        _sprite.vertices[ i++ ] = 0;  _sprite.vertices[ i++ ] = 1;  // uv 3

        i = 0;

        _sprite.faces[ i++ ] = 0; _sprite.faces[ i++ ] = 1; _sprite.faces[ i++ ] = 2;
        _sprite.faces[ i++ ] = 0; _sprite.faces[ i++ ] = 2; _sprite.faces[ i++ ] = 3;

        _sprite.vertexBuffer  = _gl.createBuffer();
        _sprite.elementBuffer = _gl.createBuffer();

        _gl.bindBuffer( _gl.ARRAY_BUFFER, _sprite.vertexBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, _sprite.vertices, _gl.STATIC_DRAW );

        _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, _sprite.elementBuffer );
        _gl.bufferData( _gl.ELEMENT_ARRAY_BUFFER, _sprite.faces, _gl.STATIC_DRAW );

        _sprite.program = createProgram( WebMol.ShaderLib.sprite, _precision );

        _sprite.attributes = {};
        _sprite.uniforms = {};

        _sprite.attributes.position           = _gl.getAttribLocation ( _sprite.program, "position" );
        _sprite.attributes.uv                 = _gl.getAttribLocation ( _sprite.program, "uv" );

        _sprite.uniforms.uvOffset             = _gl.getUniformLocation( _sprite.program, "uvOffset" );
        _sprite.uniforms.uvScale              = _gl.getUniformLocation( _sprite.program, "uvScale" );

        _sprite.uniforms.rotation             = _gl.getUniformLocation( _sprite.program, "rotation" );
        _sprite.uniforms.scale                = _gl.getUniformLocation( _sprite.program, "scale" );
        _sprite.uniforms.alignment            = _gl.getUniformLocation( _sprite.program, "alignment" );

        _sprite.uniforms.color                = _gl.getUniformLocation( _sprite.program, "color" );
        _sprite.uniforms.map                  = _gl.getUniformLocation( _sprite.program, "map" );
        _sprite.uniforms.opacity              = _gl.getUniformLocation( _sprite.program, "opacity" );

        _sprite.uniforms.useScreenCoordinates = _gl.getUniformLocation( _sprite.program, "useScreenCoordinates" );
        _sprite.uniforms.sizeAttenuation      = _gl.getUniformLocation( _sprite.program, "sizeAttenuation" );
        _sprite.uniforms.screenPosition       = _gl.getUniformLocation( _sprite.program, "screenPosition" );
        _sprite.uniforms.modelViewMatrix      = _gl.getUniformLocation( _sprite.program, "modelViewMatrix" );
        _sprite.uniforms.projectionMatrix     = _gl.getUniformLocation( _sprite.program, "projectionMatrix" );

        _sprite.uniforms.fogType              = _gl.getUniformLocation( _sprite.program, "fogType" );
        _sprite.uniforms.fogDensity           = _gl.getUniformLocation( _sprite.program, "fogDensity" );
        _sprite.uniforms.fogNear              = _gl.getUniformLocation( _sprite.program, "fogNear" );
        _sprite.uniforms.fogFar               = _gl.getUniformLocation( _sprite.program, "fogFar" );
        _sprite.uniforms.fogColor             = _gl.getUniformLocation( _sprite.program, "fogColor" );

        _sprite.uniforms.alphaTest            = _gl.getUniformLocation( _sprite.program, "alphaTest" );

    };

    this.render = function ( scene, camera, viewportWidth, viewportHeight ) {

        var sprites = scene.__webglSprites,
            nSprites = sprites.length;

        if ( ! nSprites ) return;

        var attributes = _sprite.attributes,
            uniforms = _sprite.uniforms;

        var invAspect = viewportHeight / viewportWidth;

        var halfViewportWidth = viewportWidth * 0.5,
            halfViewportHeight = viewportHeight * 0.5;

        // setup gl

        _gl.useProgram( _sprite.program );

        _gl.enableVertexAttribArray( attributes.position );
        _gl.enableVertexAttribArray( attributes.uv );

        _gl.disable( _gl.CULL_FACE );
        _gl.enable( _gl.BLEND );

        _gl.bindBuffer( _gl.ARRAY_BUFFER, _sprite.vertexBuffer );
        _gl.vertexAttribPointer( attributes.position, 2, _gl.FLOAT, false, 2 * 8, 0 );
        _gl.vertexAttribPointer( attributes.uv, 2, _gl.FLOAT, false, 2 * 8, 8 );

        _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, _sprite.elementBuffer );

        _gl.uniformMatrix4fv( uniforms.projectionMatrix, false, camera.projectionMatrix.elements );

        _gl.activeTexture( _gl.TEXTURE0 );
        _gl.uniform1i( uniforms.map, 0 );

        var oldFogType = 0;
        var sceneFogType = 0;
        var fog = scene.fog;

        if ( fog ) {

            _gl.uniform3f( uniforms.fogColor, fog.color.r, fog.color.g, fog.color.b );

            _gl.uniform1f( uniforms.fogNear, fog.near );
            _gl.uniform1f( uniforms.fogFar, fog.far );

            _gl.uniform1i( uniforms.fogType, 1 );
            oldFogType = 1;
            sceneFogType = 1;


        } 
        
        else {

            _gl.uniform1i( uniforms.fogType, 0 );
            oldFogType = 0;
            sceneFogType = 0;

        }


        // update positions and sort

        var i, sprite, material, screenPosition, size, fogType, scale = [];

        for( i = 0; i < nSprites; i ++ ) {

            sprite = sprites[ i ];
            material = sprite.material;

            if ( ! sprite.visible || material.opacity === 0 ) continue;

            if ( ! material.useScreenCoordinates ) {

                sprite._modelViewMatrix.multiplyMatrices( camera.matrixWorldInverse, sprite.matrixWorld );
                sprite.z = - sprite._modelViewMatrix.elements[ 14 ];

            } else {

                sprite.z = - sprite.position.z;

            }

        }

        sprites.sort( painterSortStable );

        // render all sprites

        for( i = 0; i < nSprites; i ++ ) {

            sprite = sprites[ i ];
            material = sprite.material;

            if ( ! sprite.visible || material.opacity === 0 ) continue;

            if ( material.map && material.map.image && material.map.image.width ) {

                _gl.uniform1f( uniforms.alphaTest, material.alphaTest );

                if ( material.useScreenCoordinates === true ) {

                    _gl.uniform1i( uniforms.useScreenCoordinates, 1 );
                    _gl.uniform3f(
                        uniforms.screenPosition,
                        ( ( sprite.position.x * _renderer.devicePixelRatio ) - halfViewportWidth  ) / halfViewportWidth,
                        ( halfViewportHeight - ( sprite.position.y * _renderer.devicePixelRatio ) ) / halfViewportHeight,
                        Math.max( 0, Math.min( 1, sprite.position.z ) )
                    );

                    scale[ 0 ] = _renderer.devicePixelRatio;
                    scale[ 1 ] = _renderer.devicePixelRatio;

                } else {

                    _gl.uniform1i( uniforms.useScreenCoordinates, 0 );
                    _gl.uniform1i( uniforms.sizeAttenuation, material.sizeAttenuation ? 1 : 0 );
                    _gl.uniformMatrix4fv( uniforms.modelViewMatrix, false, sprite._modelViewMatrix.elements );

                    scale[ 0 ] = 1;
                    scale[ 1 ] = 1;

                }

                if ( scene.fog && material.fog ) {

                    fogType = sceneFogType;

                } else {

                    fogType = 0;

                }

                if ( oldFogType !== fogType ) {

                    _gl.uniform1i( uniforms.fogType, fogType );
                    oldFogType = fogType;

                }

                size = 1 / ( material.scaleByViewport ? viewportHeight : 1 );

                scale[ 0 ] *= size * invAspect * sprite.scale.x;
                scale[ 1 ] *= size * sprite.scale.y;

                _gl.uniform2f( uniforms.uvScale, material.uvScale.x, material.uvScale.y );
                _gl.uniform2f( uniforms.uvOffset, material.uvOffset.x, material.uvOffset.y );
                _gl.uniform2f( uniforms.alignment, material.alignment.x, material.alignment.y );

                _gl.uniform1f( uniforms.opacity, material.opacity );
                _gl.uniform3f( uniforms.color, material.color.r, material.color.g, material.color.b );

                _gl.uniform1f( uniforms.rotation, sprite.rotation );
                _gl.uniform2fv( uniforms.scale, scale );

                //_renderer.setBlending( material.blending, material.blendEquation, material.blendSrc, material.blendDst );
                _renderer.setDepthTest( material.depthTest );
                _renderer.setDepthWrite( material.depthWrite );
                _renderer.setTexture( material.map, 0 );

                _gl.drawElements( _gl.TRIANGLES, 6, _gl.UNSIGNED_SHORT, 0 );

            }

        }

        // restore gl

        _gl.enable( _gl.CULL_FACE );

    };

    function createProgram ( shader, precision ) {

        var program = _gl.createProgram();

        var fragmentShader = _gl.createShader( _gl.FRAGMENT_SHADER );
        var vertexShader = _gl.createShader( _gl.VERTEX_SHADER );

        var prefix = "precision " + precision + " float;\n";

        _gl.shaderSource( fragmentShader, prefix + shader.fragmentShader );
        _gl.shaderSource( vertexShader, prefix + shader.vertexShader );

        _gl.compileShader( fragmentShader );
        _gl.compileShader( vertexShader );
        
        if ( ! _gl.getShaderParameter(fragmentShader, _gl.COMPILE_STATUS) || ! _gl.getShaderParameter(vertexShader,_gl.COMPILE_STATUS) ) {

                console.error(_gl.getShaderInfoLog(fragmentShader));
                console.error("could not initialize shader");
                return null;
        }

        _gl.attachShader( program, fragmentShader );
        _gl.attachShader( program, vertexShader );

        _gl.linkProgram( program );

        if (! _gl.getProgramParameter(program, _gl.LINK_STATUS) )
                console.error("Could not initialize shader");

        return program;

    }

    function painterSortStable ( a, b ) {

        if ( a.z !== b.z ) {

            return b.z - a.z;

        } else {

            return b.id - a.id;

        }

    }

};
/* 
 * WebMol Lighting
 */

//TODO: Strip down this class - do I really use all of these instance variables?
WebMol.Light = function(hex, intensity) {
    
    WebMol.Object3D.call(this);
    
    this.color = new WebMol.Color(hex);
    this.position = new WebMol.Vector3( 0, 1, 0 );
    this.target = new WebMol.Object3D();

    this.intensity = ( intensity !== undefined ) ? intensity : 1;

    this.castShadow = false;
    this.onlyShadow = false;
    
};

WebMol.Light.prototype = Object.create(WebMol.Object3D.prototype);
/* 
 * Line and Mesh material types
 */

//Material base class

WebMol.Material = function () {

    WebMol.EventDispatcher.call( this );

    this.id = WebMol.MaterialIdCount ++;

    this.name = '';
    
    //TODO: Which of these instance variables can I remove??
    this.side = WebMol.FrontSide;

    this.opacity = 1;
    this.transparent = false;

    this.blending = WebMol.NormalBlending;

    this.depthTest = true;
    this.depthWrite = true;

    this.polygonOffset = false;
    this.polygonOffsetFactor = 0;
    this.polygonOffsetUnits = 0;

    this.alphaTest = 0;

    this.visible = true;

    this.needsUpdate = true;

};


WebMol.Material.prototype.setValues = function ( values ) {

    if ( values === undefined ) return;

    for ( var key in values ) {

        var newValue = values[ key ];

        if ( newValue === undefined ) {

            console.warn( 'WebMol.Material: \'' + key + '\' parameter is undefined.' );
            continue;

        }

        if ( key in this ) {

            var currentValue = this[ key ];

            if ( currentValue instanceof WebMol.Color && newValue instanceof WebMol.Color ) {

                currentValue.copy( newValue );

            } else if ( currentValue instanceof WebMol.Color ) {

                currentValue.set( newValue );

            } else if ( currentValue instanceof WebMol.Vector3 && newValue instanceof WebMol.Vector3 ) {

                currentValue.copy( newValue );

            } else {

                this[ key ] = newValue;

            }

        }

    }

};
//TODO: might want to look into blending equations
WebMol.Material.prototype.clone = function ( material ) {

    if ( material === undefined ) material = new WebMol.Material();

    material.name = this.name;

    material.side = this.side;

    material.opacity = this.opacity;
    material.transparent = this.transparent;

    material.blending = this.blending;

    material.depthTest = this.depthTest;
    material.depthWrite = this.depthWrite;

    material.polygonOffset = this.polygonOffset;
    material.polygonOffsetFactor = this.polygonOffsetFactor;
    material.polygonOffsetUnits = this.polygonOffsetUnits;

    material.alphaTest = this.alphaTest;

    material.overdraw = this.overdraw;

    material.visible = this.visible;

    return material;

};

WebMol.Material.prototype.dispose = function () {

    this.dispatchEvent( { type: 'dispose' } );

};

WebMol.MaterialIdCount = 0;

//Line basic material

WebMol.LineBasicMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    
    this.linewidth = 1;
    this.linecap = 'round';
    this.linejoin = 'round';
    
    this.vertexColors = false;
    
    this.fog = true;
    
    this.setValues(parameters);
    
};

WebMol.LineBasicMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.LineBasicMaterial.prototype.clone = function() {
  
    var material = new WebMol.LineBasicMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy();
    
};

//Mesh Lambert material

WebMol.MeshLambertMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    this.ambient = new WebMol.Color(0xfffff);
    this.emissive = new WebMol.Color(0x000000);
    
    //TODO: Which of these instance variables do I really need?
    this.wrapAround = false;
    this.wrapRGB = new WebMol.Vector3(1,1,1);
    
    this.map = null;
    
    this.lightMap = null;
    
    this.specularMap = null;
    
    this.envMap = null;
    this.reflectivity = 1;
    this.refractionRatio = 0.98;
    
    this.fog = true;
    
    this.wireframe = false;
    this.wireframeLinewidth = 1;
    this.wireframeLinecap = 'round';
    this.wireframeLinejoin = 'round';
    
    this.shading = WebMol.SmoothShading;
    
    this.vertexColors = WebMol.NoColors;
    
    this.skinning = false;
    
    this.setValues(parameters);
    
};

WebMol.MeshLambertMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.MeshLambertMaterial.prototype.clone = function() {
  
    var material = new WebMol.MeshLambertMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy(this.color);
    material.ambient.copy(this.ambient);
    material.emissive.copy(this.emissive);
    
    material.wrapAround = this.wrapAround;
    material.wrapRGB.copy(this.wrapRGB);
    
    material.map = this.map;
    
    material.lightMap = this.lightMap;
    
    material.specularMap = this.specularMap;
    
    material.envMap = this.envMap;
    material.combine = this.combine;
    material.reflectivity = this.reflectivity;
    material.refractionRatio = this.refractionRatio;
    
    material.fog = this.fog;
    
    material.shading = this.shading;
    
    material.vertexColors = this.vertexColors;
    
    material.skinning = this.skinning;
    material.morphTargets = this.morphTargets;
    material.morphNormals = this.morphNormals;
    
    return material;
    
};


//Sprite material
WebMol.SpriteMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    this.map = new WebMol.Texture();
    
    this.useScreenCoordinates = true;
    this.depthTest = !this.useScreenCoordinates;
    this.sizeAttenuation = !this.useScreenCoordinates;
    this.scaleByViewPort = !this.sizeAttenuation;
    this.alignment = WebMol.SpriteAlignment.center.clone();
    
    this.fog = false; // use scene fog
    
    this.uvOffset = new WebMol.Vector2(0, 0);
    this.uvScale = new WebMol.Vector2(1, 1);
    
    this.setValues(parameters);
    
    parameters = parameters || {};
    
    if (parameters.depthTest === undefined)
        this.depthTest = !this.useScreenCoordinates;
    if (parameters.sizeAttenuation === undefined)
        this.sizeAttenuation = !this.useScreenCoordinates;
    if (parameters.scaleByViewPort === undefined)
        this.scaleByViewPort = !this.sizeAttenuation;
    
};

WebMol.SpriteMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.SpriteMaterial.prototype.clone = function() {
    
    var material = new WebMol.SpriteMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy(this.color);
    material.map = this.map;
    
    material.useScreenCoordinates = useScreenCoordinates;
    material.sizeAttenuation = this.sizeAttenuation;
    material.scaleByViewport = this.scaleByViewPort;
    material.alignment.copy(this.alignment);
    
    material.uvOffset.copy(this.uvOffset);
    
    return material;
    
};

//Alignment for Sprites

WebMol.SpriteAlignment = {};
WebMol.SpriteAlignment.topLeft = new WebMol.Vector2(1, -1);
WebMol.SpriteAlignment.topCenter = new WebMol.Vector2(0, -1);
WebMol.SpriteAlignment.topRight = new WebMol.Vector2(-1, -1);
WebMol.SpriteAlignment.centerLeft = new WebMol.Vector2(1, 0);
WebMol.SpriteAlignment.center = new WebMol.Vector2(0, 0);
WebMol.SpriteAlignment.centerRight = new WebMol.Vector2(-1, 0);
WebMol.SpriteAlignment.bottomLeft = new WebMol.Vector2(1, 1);
WebMol.SpriteAlignment.bottomCenter = new WebMol.Vector2(0, 1);
WebMol.SpriteAlignment.bottomRight = new WebMol.Vector2(-1, 1);


//Texture
//We really only create textures from 2d rendering contexts (to display text labels)

WebMol.Texture = function(image) {

    WebMol.EventDispatcher.call(this);
    
    this.id = WebMol.TextureIdCount++;
    
    this.name = "";
    
    this.image = image;
    this.mipmaps = [];
    
    this.mapping = new WebMol.UVMapping();
    
    this.wrapS = WebMol.ClampToEdgeWrapping;
    this.wrapT = WebMol.ClampToEdgeWrapping;
    
    this.magFilter = WebMol.LinearFilter;
    this.minFilter = WebMol.LinearMipMapLinearFilter;
    
    this.anisotropy = 1;
    
    this.format = WebMol.RGBAFormat;
    this.type = WebMol.UnsignedByteType;
    
    this.offset = new WebMol.Vector2(0, 0);
    this.repeat = new WebMol.Vector2(1, 1);
    
    this.generateMipmaps = true;
    this.premultiplyAlpha = false;
    this.flipY = true;
    this.unpackAlignment = 4;
    
    this.needsUpdate = false;
    this.onUpdate = null;
    
};

WebMol.Texture.prototype = {

    constructor : WebMol.Texture,
    
    clone : function(texture) {
        
        if (texture === undefined)
            texture = new WebMol.Texture();
        
        texture.image = this.image;
        texture.mipmaps = this.mipmaps.slice(0);
        
        texture.mapping = this.mapping;
        
        texture.wrapS = this.wrapS;
        texture.wrapT = this.wrapT;
        
        texture.magFilter = this.magFilter;
        texture.minFilter = this.minFilter;
        
        texture.anisotropy = this.anisotropy;
        
        texture.format = this.format;
        texture.type = this.type;
        
        texture.offset.copy(this.offset);
        texture.repeat.copy(this.repeat);
        
        texture.generateMipmaps = this.generateMipmaps;
        texture.premultiplyAlpha = this.premultiplyAlpha;
        texture.flipY = this.flipY;
        texture.unpackAlignment = this.unpackAlignment;
        
        return texture;
        
    },
    
    dispose : function() {
        
        this.dispatchEvent( {type: 'dispose'});
        
    }    
    
};

WebMol.TextureIdCount = 0;


/* 
 * WebMol Mesh and Line objects
 */


//Line Object

WebMol.Line = function (geometry, material, type) {

    WebMol.Object3D.call(this);

    this.geometry = geometry;
        //TODO: update material and type to webgl
    this.material = (material !== undefined) ? material : new WebMol.LineBasicMaterial( { color: Math.random() * 0xffffff } );
    this.type = (type !== undefined) ? type : WebMol.LineStrip;

};

WebMol.LineStrip = 0;
WebMol.LinePieces = 1;

WebMol.Line.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Line.prototype.clone = function (object) {

    if (object === undefined) object = new WebMol.Line(this.geometry, this.material, this.type);

    WebMol.Object3D.prototype.clone.call(this, object);

    return object;

};


//Mesh Object

WebMol.Mesh = function(geometry, material) {

    WebMol.Object3D.call(this);

    this.geometry = geometry;
    this.material = (material !== undefined) ? material : new WebMol.MeshBasicMaterial( { color: Math.random() * 0xffffff, wireframe: true } );

};

WebMol.Mesh.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Mesh.prototype.clone = function (object) {

    if (object === undefined) object = new WebMol.Mesh(this.geometry, this.material);

    WebMol.Object3D.prototype.clone.call(this, object);

    return object;

};


//Sprite object

WebMol.Sprite = function(material) {
    
    WebMol.Object3D.call(this);
    
    this.material = (material !== undefined) ? material : new WebMol.SpriteMaterial();

    this.rotation3d = this.rotation;
    this.rotation = 0;
    
};

WebMol.Sprite.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Sprite.prototype.updateMatrix = function() {
    
    this.matrix.setPosition(this.position);
    
    this.rotation3d.set(0, 0, this.rotation);
    this.matrix.setRotationFromEuler(this.rotation3d);
    
    if (this.scale.x !== 1 || this.scale.y !== 1)
        this.matrix.scale(this.scale);
    
    this.matrixWorldNeedsUpdate = true;
    
};

WebMol.Sprite.prototype.clone = function(object) {
    
    if (object === undefined)
        object = new WebMol.Sprite(this.material);
    
    WebMol.Object3D.prototype.clone.call(this, object);
    
    return object;
    
};
/**
Simplified webGL renderer 
 */

WebMol.Renderer = function ( parameters ) {
    
    parameters = parameters || {};
    
    var _canvas = parameters.canvas !== undefined ? parameters.canvas : document.createElement( 'canvas' ),

    _precision = parameters.precision !== undefined ? parameters.precision : 'highp',

    _alpha = parameters.alpha !== undefined ? parameters.alpha : true,
    _premultipliedAlpha = parameters.premultipliedAlpha !== undefined ? parameters.premultipliedAlpha : true,
    _antialias = parameters.antialias !== undefined ? parameters.antialias : false,
    _stencil = parameters.stencil !== undefined ? parameters.stencil : true,
    _preserveDrawingBuffer = parameters.preserveDrawingBuffer !== undefined ? parameters.preserveDrawingBuffer : false,

    _clearColor = parameters.clearColor !== undefined ? new WebMol.Color( parameters.clearColor ) : new WebMol.Color( 0x000000 ),
    _clearAlpha = parameters.clearAlpha !== undefined ? parameters.clearAlpha : 0;
    
    this.domElement = _canvas;
    this.context = null;
    this.devicePixelRatio = parameters.devicePixelRatio !== undefined ? 
    parameters.devicePixelRatio : (self.devicePixelRatio !== undefined) ? 
                                   self.devicePixelRatio : 1;

    // clearing

    this.autoClear = true;
    this.autoClearColor = true;
    this.autoClearDepth = true;
    this.autoClearStencil = true;

    // scene graph

    this.sortObjects = true;

    this.autoUpdateObjects = true;
    this.autoUpdateScene = true;
    
    this.renderPluginsPost = [];
    
    // info

    this.info = {

        memory: {
    
        programs: 0,
        geometries: 0,
        textures: 0
    
        },
    
        render: {
    
        calls: 0,
        vertices: 0,
        faces: 0,
        points: 0
    
        }

    };

    // internal properties

    var _this = this,

    _programs = [],
    _programs_counter = 0,

    // internal state cache

    _currentProgram = null,
    _currentFramebuffer = null,
    _currentMaterialId = -1,
    _currentGeometryGroupHash = null,
    _currentCamera = null,
    _geometryGroupCounter = 0,

    _usedTextureUnits = 0,

    // GL state cache

    _oldDoubleSided = -1,
    _oldFlipSided = -1,

    _oldBlending = -1,

    _oldBlendEquation = -1,
    _oldBlendSrc = -1,
    _oldBlendDst = -1,

    _oldDepthTest = -1,
    _oldDepthWrite = -1,

    _oldPolygonOffset = null,
    _oldPolygonOffsetFactor = null,
    _oldPolygonOffsetUnits = null,

    _oldLineWidth = null,

    _viewportX = 0,
    _viewportY = 0,
    _viewportWidth = 0,
    _viewportHeight = 0,
    _currentWidth = 0,
    _currentHeight = 0,

    _enabledAttributes = {},

     // camera matrices cache

    _projScreenMatrix = new WebMol.Matrix4(),

    _vector3 = new WebMol.Vector3(),

    // light arrays cache

    _direction = new WebMol.Vector3(),

    _lightsNeedUpdate = true,

    _lights = {

            ambient: [0,0,0],
            directional: { length: 0, colors: [], positions: [] },
            point: { length: 0, colors: [], positions: [], distances: [] },
            spot: { length: 0, colors: [], positions: [], distances: [], directions: [], anglesCos: [], exponents: [] },
            hemi: { length: 0, skyColors: [], groundColors: [], positions: [] }

    };

    // initialize

    var _gl;

    initGL();

    setDefaultGLState();

    this.context = _gl;    

    // API

    this.getContext = function () {

            return _gl;

    };

    this.getPrecision = function () {

            return _precision;

    };
    
    this.setClearColorHex = function ( hex, alpha ) {

            _clearColor.setHex( hex );
            _clearAlpha = alpha;

            _gl.clearColor( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha );

    };

    this.setSize = function ( width, height ) {

            _canvas.width = width * this.devicePixelRatio;
            _canvas.height = height * this.devicePixelRatio;

            _canvas.style.width = width + 'px';
            _canvas.style.height = height + 'px';

            this.setViewport( 0, 0, _canvas.width, _canvas.height );

    };

    this.setViewport = function ( x, y, width, height ) {

            _viewportX = x !== undefined ? x : 0;
            _viewportY = y !== undefined ? y : 0;

            _viewportWidth = width !== undefined ? width : _canvas.width;
            _viewportHeight = height !== undefined ? height : _canvas.height;

            _gl.viewport( _viewportX, _viewportY, _viewportWidth, _viewportHeight );

    };

    this.clear = function ( color, depth, stencil ) {

            var bits = 0;

            if ( color === undefined || color ) bits |= _gl.COLOR_BUFFER_BIT;
            if ( depth === undefined || depth ) bits |= _gl.DEPTH_BUFFER_BIT;
            if ( stencil === undefined || stencil ) bits |= _gl.STENCIL_BUFFER_BIT;

            _gl.clear( bits );

    };

    this.clearTarget = function ( renderTarget, color, depth, stencil ) {

            this.setRenderTarget( renderTarget );
            this.clear( color, depth, stencil );

    };

    this.setMaterialFaces = function ( material ) {

            var doubleSided = material.side === WebMol.DoubleSide;
            var flipSided = material.side === WebMol.BackSide;

            if ( _oldDoubleSided !== doubleSided ) {

                if ( doubleSided ) {

                    _gl.disable( _gl.CULL_FACE );

                } else {

                    _gl.enable( _gl.CULL_FACE );

                }

                _oldDoubleSided = doubleSided;

            }

            if ( _oldFlipSided !== flipSided ) {

                if ( flipSided ) {

                    _gl.frontFace( _gl.CW );

                } else {

                    _gl.frontFace( _gl.CCW );

                }

                _oldFlipSided = flipSided;

            }    

    };
    
    this.setDepthTest = function ( depthTest ) {

            if ( _oldDepthTest !== depthTest ) {

                if ( depthTest ) {

                    _gl.enable( _gl.DEPTH_TEST );

                } else {

                    _gl.disable( _gl.DEPTH_TEST );

                }

                _oldDepthTest = depthTest;

            }

    };

    this.setDepthWrite = function ( depthWrite ) {

            if ( _oldDepthWrite !== depthWrite ) {

                    _gl.depthMask( depthWrite );
                    _oldDepthWrite = depthWrite;

            }

    };

    this.setBlending = function( blending ) {

            if (blending === WebMol.NoBlending) 
                    _gl.disable( _gl.BLEND );

            else {
                    _gl.enable( _gl.BLEND );
                    _gl.blendEquationSeparate( _gl.FUNC_ADD, _gl.FUNC_ADD );
                    _gl.blendFuncSeparate( _gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA, _gl.ONE, _gl.ONE_MINUS_SRC_ALPHA );
            }

            _oldBlending = blending;
    };
    
    // Plugins
    
    this.addPostPlugin = function(plugin) {

        plugin.init(this);
        this.renderPluginsPost.push(plugin);

    };

    // Sorting

    function numericalSort ( a, b ) {

            return b[ 0 ] - a[ 0 ];

    }

    function enableAttribute( attribute ) {

        if ( ! _enabledAttributes[ attribute ] ) {

            _gl.enableVertexAttribArray( attribute );
            _enabledAttributes[ attribute ] = true;

        }

    }

    function disableAttributes() {

        for ( var attribute in _enabledAttributes ) {

            if ( _enabledAttributes[ attribute ] ) {

                _gl.disableVertexAttribArray( attribute );
                _enabledAttributes[ attribute ] = false;

            }

        }

    } 

    function setPolygonOffset ( polygonOffset, factor, units) {

        if ( _oldPolygonOffset !== polygonOffset ) {

            if (polygonOffset)
                _gl.enable( _gl.POLYGON_OFFSET_FILL );
            else
                _gl.disable( _gl.POLYGON_OFFSET_FILL );
        }
    }

    function setLineWidth ( width ) {

        if ( width !== _oldLineWidth ) {
            _gl.lineWidth(width);
            _oldLineWidth = width;
        }

    }
    
    var onGeometryDispose = function(event) {
        
        var geometry = event.target;
        geometry.removeEventListener('dispose', onGeometryDispose);
        
        deallocateGeometry(geometry);
        
        _this.info.memory.geometries--;
        
    };
    
    var onTextureDispose = function(event) {

        var texture = event.target;

        texture.removeEventListener('dispose', onTextureDispose);

        deallocateTexture(texture);

        _this.info.memory.textures--;


    };
    
    var onMaterialDispose = function(event) {
        
        var material = event.target;
        material.removeEventListener('dispose', onMaterialDispose);
        
        deallocateMaterial(material);
        
    };
    
    var deallocateGeometry = function(geometry) {
        
        geometry.__webglInit = undefined;
        
        if (geometry.__webglVertexBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglVertexBuffer);
        
        if (geometry.__webglColorBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglColorBuffer);
        
        if (geometry.geometryGroups !== undefined) {
            
            for (var g = 0, gl = geometry.groups; g < gl; g++) {  
                
                var geometryGroup = geometry.geometryGroups[g];

                if (geometryGroup.__webglVertexBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglVertexBuffer);

                if (geometryGroup.__webglColorBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglColorBuffer);
                
                if (geometryGroup.__webglNormalBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglNormalBuffer);  
                
                if (geometryGroup.__webglFaceBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglFaceBuffer);
                    
                if (geometryGroup.__webglLineBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglLineBuffer);
                    
            }
        }
    };
    
    var deallocateMaterial = function (material) {

        var program = material.program;

        if ( program === undefined ) return;

        material.program = undefined;

        // only deallocate GL program if this was the last use of shared program
        // assumed there is only single copy of any program in the _programs list
        // (that's how it's constructed)

        var i, il, programInfo;
        var deleteProgram = false;

        for ( i = 0, il = _programs.length; i < il; i ++ ) {

            programInfo = _programs[ i ];

            if ( programInfo.program === program ) {

                programInfo.usedTimes --;

                if ( programInfo.usedTimes === 0 ) {

                    deleteProgram = true;

                }

                break;

            }

        }

        if ( deleteProgram === true ) {

            // avoid using array.splice, this is costlier than creating new array from scratch

            var newPrograms = [];

            for ( i = 0, il = _programs.length; i < il; i ++ ) {

                programInfo = _programs[ i ];

                if ( programInfo.program !== program ) {

                    newPrograms.push( programInfo );

                }

            }

            _programs = newPrograms;

            _gl.deleteProgram( program );

            _this.info.memory.programs --;

        }

    };
    
    var deallocateTexture = function(texture) {

        if (texture.image && texture.image.__webglTextureCube) {

            // cube texture

            _gl.deleteTexture(texture.image.__webglTextureCube);

        } 
        
        else {

            // 2D texture

            if ( ! texture.__webglInit ) return;

            texture.__webglInit = false;
            _gl.deleteTexture( texture.__webglTexture );

        }

    };

    //Compile and return shader
    function getShader (type, str) {

        var shader;

        if (type === "fragment")
            shader = _gl.createShader( _gl.FRAGMENT_SHADER );
        else if (type === "vertex")
            shader = _gl.createShader( _gl.VERTEX_SHADER );

        _gl.shaderSource(shader, str);
        _gl.compileShader(shader);

        if ( ! _gl.getShaderParameter(shader, _gl.COMPILE_STATUS) ) {

            console.error(_gl.getShaderInfoLog(shader));
            console.error("could not initialize shader");
            return null;

        }

        return shader;

    } 


    //Compile appropriate shaders (if necessary) from source code and attach to gl program.
    function buildProgram(fragmentShader, vertexShader, uniforms, parameters) {

        var p, pl, d, program, code;
        var chunks = [];

        chunks.push(fragmentShader);
        chunks.push(vertexShader);
        
        for (p in parameters) {
            chunks.push(p);
            chunks.push(parameters[p]);
        }
        
        code = chunks.join();

        //check if program has already been compiled

        for (p = 0, pl = _programs.length; p < pl; p++) {

            var programInfo = _programs[p];

            if (programInfo.code === code) {

                programInfo.usedTimes++;

                return programInfo.program;
            }
        }

        //Set up new program and compile shaders

        program = _gl.createProgram();
        
        //set up precision
        var precision = _precision;
        var prefix = "precision " + precision + " float;";
        
        var prefix_vertex = [
                             prefix
                            ].join("\n");
                            
        var prefix_fragment = [                               
                               parameters.wireframe ? "#define WIREFRAME 1" : "",
                               prefix
                              ].join("\n");
        
        var glFragmentShader = getShader("fragment", prefix_fragment + fragmentShader);
        var glVertexShader = getShader("vertex", prefix_vertex + vertexShader);

        _gl.attachShader(program, glVertexShader);
        _gl.attachShader(program, glFragmentShader);

        _gl.linkProgram(program);

        if (! _gl.getProgramParameter(program, _gl.LINK_STATUS) )
            console.error("Could not initialize shader");

        //gather and cache uniform variables and attributes

        program.uniforms = {};
        program.attributes = {};

        var identifiers, u, a, i;

        //uniform vars
        identifiers = 
            [ 'viewMatrix', 'modelViewMatrix', 'projectionMatrix', 'normalMatrix', 'modelMatrix', 'cameraPosition' ];

        //custom uniform vars
        for (u in uniforms) 
            identifiers.push(u);

        for (i = 0; i < identifiers.length; i++) {

            var uniformVar = identifiers[i];
            program.uniforms[uniformVar] = _gl.getUniformLocation(program, uniformVar);

        }

        //attributes
        identifiers = 
            [ 'position', 'normal', 'color', 'lineDistance' ];

        /*
        for (a in attributes)
                identifiers.push(a);
        */

        for (i = 0; i < identifiers.length; i++) {

            var attributeVar = identifiers[i];
            program.attributes[attributeVar] = _gl.getAttribLocation(program, attributeVar);
        }

        program.id = _programs_counter++;
        _programs.push( {program: program, code: code, usedTimes: 1} );
        _this.info.memory.programs = _programs.length;

        return program;
    }

    //TODO: need to set up shader attributes and uniforms as attributes on material object after attaching prgm
    //We need to attach appropriate uniform variables to material after shaders have been chosen
    this.initMaterial = function ( material, lights, fog, object ) {

        material.addEventListener('dispose', onMaterialDispose);

        var u, a, identifiers, i, parameters, maxLightCount, maxBones, maxShadows, shaderID;

        if (material instanceof WebMol.LineBasicMaterial)
            shaderID = "basic";
        else if (material instanceof WebMol.MeshLambertMaterial)
            shaderID = "lambert";

        if (shaderID) {

            var shader = WebMol.ShaderLib[shaderID];
            material.shaderType = shaderID;
            material.vertexShader = shader.vertexShader;
            material.fragmentShader = shader.fragmentShader;
            material.uniforms = WebMol.ShaderUtils.clone(shader.uniforms);
            //TODO: set material uniforms to shader uniform variables

        }
        
        parameters = {
            wireframe: material.wireframe
        };

        material.program = buildProgram(material.fragmentShader, material.vertexShader, material.uniforms, parameters);

    };

    function setProgram( camera, lights, fog, material, object ) {

        if ( material.needsUpdate ) {

            if (material.program)
                deallocateMaterial(material);

                _this.initMaterial( material, lights, fog, object );
                material.needsUpdate = false;
        }

        var refreshMaterial = false;

        //p_uniforms: uniformVarName => uniformLocation
        //m_uniforms: uniformVarName => uniformJsVal
        var program = material.program,
            p_uniforms = program.uniforms,
            m_uniforms = material.uniforms;

        if (program != _currentProgram) {        
            _gl.useProgram(program);
            _currentProgram = program;

            refreshMaterial = true;
        }

        if (material.id != _currentMaterialId) {
            _currentMaterialId = material.id;
            refreshMaterial = true;
        }

        if (camera != _currentCamera) {    
            _currentCamera = camera;
            refreshMaterial = true;
        }

        //Send projection matrix to uniform variable in shader
        if (refreshMaterial) {

            //Load projection, model-view matrices for perspective
            _gl.uniformMatrix4fv(p_uniforms.projectionMatrix, false, camera.projectionMatrix.elements);
            _gl.uniformMatrix4fv(p_uniforms.modelViewMatrix, false, object._modelViewMatrix.elements);

            //Set up correct fog uniform vals
            m_uniforms.fogColor.value = fog.color;
            m_uniforms.fogNear.value = fog.near;
            m_uniforms.fogFar.value = fog.far;

            //Set up lights for lambert shader
            if (material.shaderType === "lambert") {

                //load view and normal matrices for directional and object lighting
                _gl.uniformMatrix4fv(p_uniforms.viewMatrix, false, camera.matrixWorldInverse.elements);
                _gl.uniformMatrix3fv(p_uniforms.normalMatrix, false, object._normalMatrix.elements);
                //_gl.uniformMatrix4fv(p_uniforms.modelMatrix, false, object.matrixWorld.elements);

                if (_lightsNeedUpdate) {
                    setupLights(program, lights);
                    _lightsNeedUpdate = false;
                }

                //Set up correct light uniform var vals
                m_uniforms.ambientLightColor.value = _lights.ambient;
                m_uniforms.directionalLightColor.value = _lights.directional.colors;
                m_uniforms.directionalLightDirection.value = _lights.directional.positions;
                m_uniforms.ambient.value = material.ambient;
                m_uniforms.emissive.value = material.emissive;

            }

            //opacity, diffuse, emissive, etc
            m_uniforms.opacity.value = material.opacity;
            m_uniforms.diffuse.value = material.color;

            //Load any other material specific uniform variables to gl shaders
            loadMaterialUniforms(p_uniforms, m_uniforms);

        }

        return program;

    }

    function loadMaterialUniforms(p_uniforms, m_uniforms) {
        var uniformVar, type, uniformVal, uniformLoc;

        for (uniformVar in m_uniforms) {
            if (! p_uniforms[uniformVar])
                continue;

            type = m_uniforms[uniformVar].type;
            uniformVal = m_uniforms[uniformVar].value;
            uniformLoc = p_uniforms[uniformVar];

            //single float
            if (type === 'f')
                _gl.uniform1f(uniformLoc, uniformVal);
            //array of floats
            else if (type === 'fv')
                _gl.uniform3fv(uniformLoc, uniformVal);
            //color - r,g,b floats
            else if (type === 'c')
                _gl.uniform3f(uniformLoc, uniformVal.r, uniformVal.g, uniformVal.b);

        }

    }

    this.renderBuffer = function ( camera, lights, fog, material, geometryGroup, object ) {

        if ( ! material.visible )
            return;

        var program, attributes, linewidth, primitives, a, attribute, i, il;

        //Sets up proper vertex and fragment shaders and attaches them to webGL program
        //Also sets appropriate uniform variables 
        program = setProgram(camera, lights, fog, material, object);

        attributes = program.attributes;

        var updateBuffers = false,
            wireframeBit = material.wireframe ? 1 : 0,
            geometryGroupHash = (geometryGroup.id * 0xffffff) + (program.id * 2) + wireframeBit;

        if (geometryGroupHash !== _currentGeometryGroupHash) {
            _currentGeometryGroupHash = geometryGroupHash;
            updateBuffers = true;
        }

        //rebind shader attributes to appropriate (and already initialized) gl buffers
        if (updateBuffers) {

            disableAttributes();

            // Vertices
            if (attributes.position >= 0) {            
                _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglVertexBuffer );
                enableAttribute( attributes.position );
                _gl.vertexAttribPointer( attributes.position, 3, _gl.FLOAT, false, 0, 0 );    
            }

            // Colors
            if (attributes.color >= 0) {
                _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglColorBuffer);
                enableAttribute( attributes.color );
                _gl.vertexAttribPointer( attributes.color, 3, _gl.FLOAT, false, 0, 0 );
            }

            // Normals (lambert shader only)
            if (attributes.normal >=0) {
                _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglNormalBuffer );
                enableAttribute( attributes.normal );
                _gl.vertexAttribPointer( attributes.normal, 3, _gl.FLOAT, false, 0, 0 );
            }

        }

        //Render
        var faceCount, lineCount;
        //lambert shaders - draw triangles
        //TODO: make sure geometryGroup's face count is setup correctly
        if (object instanceof WebMol.Mesh) {
            
            if (material.wireframe) {
                lineCount = geometryGroup.lineidx;
                setLineWidth(material.wireframeLinewidth);
                
                if (updateBuffers)
                    _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglLineBuffer );
                
                _gl.drawElements( _gl.LINES, lineCount, _gl.UNSIGNED_SHORT, 0 );
            }
            
            else {
                faceCount = geometryGroup.faceidx;

                if (updateBuffers)
                    _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglFaceBuffer );
                
                _gl.drawElements( _gl.TRIANGLES, faceCount, _gl.UNSIGNED_SHORT, 0 );
                
            }


            _this.info.render.calls++;
            _this.info.render.vertices += faceCount;
            _this.info.render.faces += faceCount / 3;
        }

        //basic shaders - draw lines
        else if (object instanceof WebMol.Line) {
            lineCount = geometryGroup.vertices;

            setLineWidth(material.linewidth);
            _gl.drawArrays( _gl.LINES, 0, lineCount );

            _this.info.render.calls++;
        }

    };

    //rendering
    function renderObjects ( renderList, reverse, materialType, camera, lights, fog, useBlending, overrideMaterial)  {

        var webglObject, object, buffer, material, start, end, delta;

        //Forward or backward render

        if (reverse) {
            start = renderList.length - 1;
            end = -1;
            delta = -1;
        }

        else {
            start = 0;
            end = renderList.length;
            delta = 1;
        }

        for (var i = start; i !== end; i += delta) {

            webglObject = renderList[i];

            if (webglObject.render) {

                object = webglObject.object;
                buffer = webglObject.buffer;
                material = webglObject[materialType];

                if ( ! material )
                    continue;

                if (useBlending)
                    _this.setBlending(material.blending);

                _this.setDepthTest(material.depthTest);
                _this.setDepthWrite(material.depthWrite);
                setPolygonOffset(material.polygonOffset, material.polygonOffsetFactor, material.polygonOffsetUnits);

                _this.setMaterialFaces(material);

                _this.renderBuffer(camera, lights, fog, material, buffer, object);
            }
        }

    }
    
    this.render = function ( scene, camera, renderTarget, forceClear ) {

        if ( camera instanceof WebMol.Camera === false )  {

            console.error( 'WebMol.Renderer.render: camera is not an instance of WebMol.Camera.' );
            return;

        }

        var i, il,

        webglObject, object,
        renderList,

        lights = scene.__lights,
        fog = scene.fog;

        // reset caching for this frame

        _currentMaterialId = -1;
        _lightsNeedUpdate = true;

        // update scene graph

        if ( this.autoUpdateScene ) scene.updateMatrixWorld();

        // update camera matrices
        //Pretty sure camera's parent is always going to be undefined for our purposes...
        if ( camera.parent === undefined ) camera.updateMatrixWorld();

        camera.matrixWorldInverse.getInverse( camera.matrixWorld );

        _projScreenMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

        // update WebGL objects

        if ( this.autoUpdateObjects ) this.initWebGLObjects( scene );


        _this.info.render.calls = 0;
        _this.info.render.vertices = 0;
        _this.info.render.faces = 0;
        _this.info.render.points = 0;

        _currentWidth = _viewportWidth;
        _currentHeight = _viewportHeight;

        if ( this.autoClear || forceClear ) {

            this.clear( this.autoClearColor, this.autoClearDepth, this.autoClearStencil );

        }

        // set matrices for regular objects (frustum culled)

        renderList = scene.__webglObjects;

        for ( i = 0, il = renderList.length; i < il; i ++ ) {

            webglObject = renderList[ i ];
            object = webglObject.object;

            webglObject.render = false;

            if ( object.visible ) {        
                setupMatrices( object, camera );
                unrollBufferMaterial( webglObject );
                webglObject.render = true;
            }
        }

        // set matrices for immediate objects

        var material = null;

        // opaque pass (front-to-back order)

        this.setBlending( WebMol.NoBlending );

        renderObjects( scene.__webglObjects, true, "opaque", camera, lights, fog, false, material );
        
        //prime depth buffer
        renderObjects( scene.__webglObjects, true, "blank", camera, lights, fog, true, material );

        // transparent pass (back-to-front order)

        renderObjects( scene.__webglObjects, false, "transparent", camera, lights, fog, true, material );

        // Render plugins (e.g. sprites), and reset state
        
        renderPlugins(this.renderPluginsPost, scene, camera);

        // Ensure depth buffer writing is enabled so it can be cleared on next render

        this.setDepthTest( true );
        this.setDepthWrite( true );

        //_gl.finish();

    };
    
    function renderPlugins(plugins, scene, camera) {
        
        //Reset state once regardless
        //This should also fix cartoon render bug (after transparent surface render)
        
        _currentGeometryGroupHash = -1;
        _currentProgram = null;
        _currentCamera = null;
        _oldBlending = -1;
        _oldDepthWrite = -1;
        _oldDepthTest = -1;
        _oldDoubleSided = -1;
        _currentMaterialId = -1;
        _oldFlipSided = -1;
        
        
        if (!plugins.length)
            return;
        
        for (var i = 0, il = plugins.length; i < il; i++) {
            
            _lightsNeedUpdate = true;
            
            plugins[i].render(scene, camera, _currentWidth, _currentHeight);
            
            //Reset state after plugin render
            _currentGeometryGroupHash = -1;
            _currentProgram = null;
            _currentCamera = null;
            _oldBlending = -1;
            _oldDepthWrite = -1;
            _oldDepthTest = -1;
            _oldDoubleSided = -1;
            _currentMaterialId = -1;
            _oldFlipSided = -1;       
                
        }  
        
    }

    this.initWebGLObjects = function ( scene ) {

        if ( !scene.__webglObjects ) {

            scene.__webglObjects = [];
            scene.__webglObjectsImmediate = [];
            scene.__webglSprites = [];
            scene.__webglFlares = [];

        }

        //Add objects; this sets up buffers for each geometryGroup
        if (scene.__objectsAdded.length) {
            
            while(scene.__objectsAdded.length){
                addObject(scene.__objectsAdded[0], scene);
                scene.__objectsAdded.splice(0, 1);
            }
            
            //Force buffer update during render
            //Hackish fix for initial cartoon-render-then-transparent-surface bug
            _currentGeometryGroupHash = -1;
            
        }

        while (scene.__objectsRemoved.length) {

            removeObject(scene.__objectsRemoved[ 0 ], scene);
            scene.__objectsRemoved.splice(0, 1);

        }

        // update must be called after objects adding / removal
        //This sends typed arrays to GL buffers for each geometryGroup
        for ( var o = 0, ol = scene.__webglObjects.length; o < ol; o ++ ) {

            updateObject(scene.__webglObjects[ o ].object);

        }

    };
    
    // Objects adding

    function addObject (object, scene) {

        var g, gl, geometry, material, geometryGroup;

        if ( !object.__webglInit ) {

            object.__webglInit = true;

            object._modelViewMatrix = new WebMol.Matrix4();
            object._normalMatrix = new WebMol.Matrix3();

            if (object.geometry !== undefined && object.geometry.__webglInit === undefined) {

                object.geometry.__webglInit = true;
                object.geometry.addEventListener('dispose', onGeometryDispose);

            }
            
            if (object instanceof WebMol.Mesh || object instanceof WebMol.Line) {
                geometry = object.geometry;
                material = object.material;           
    
                for (g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
    
                    geometryGroup = geometry.geometryGroups[ g ];
                    
                    geometryGroup.id = _geometryGroupCounter++;

                    // initialise VBO on the first access

                    if ( !geometryGroup.__webglVertexBuffer ) {
                            
                        if (object instanceof WebMol.Mesh) {
                            createMeshBuffers(geometryGroup);
                            geometry.elementsNeedUpdate = true;
                            geometry.normalsNeedUpdate = true;
                        }
                            
                        else if (object instanceof WebMol.Line)
                            createLineBuffers(geometryGroup);

                        geometry.verticesNeedUpdate = true;
                        geometry.colorsNeedUpdate = true;

                    }
                        
                }
                
            }
        
        }
        
        if ( ! object.__webglActive ) {
            
            if (object instanceof WebMol.Mesh || object instanceof WebMol.Line) {
                
                geometry = object.geometry;

                for ( g = 0, gl = geometry.geometryGroups.length; g < gl; g++ ) {
                    geometryGroup = geometry.geometryGroups[g];

                    addBuffer(scene.__webglObjects, geometryGroup, object);
                }
                
            }
            
            //Sprite
            else if (object instanceof WebMol.Sprite) 
                scene.__webglSprites.push(object);
         
                     
            object.__webglActive = true;
            
        }

    }

    function updateObject ( object ) {

        var geometry = object.geometry, material = object.material,
                geometryGroup, customAttributesDirty;
        
        if ( object instanceof WebMol.Mesh || object instanceof WebMol.Line ) {
            
            for (var g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
                
                geometryGroup = geometry.geometryGroups[ g ];

                if ( geometry.verticesNeedUpdate || geometry.elementsNeedUpdate || geometry.colorsNeedUpdate || geometry.normalsNeedUpdate) {
                    setBuffers( geometryGroup, _gl.DYNAMIC_DRAW );
                }
            }
            
            geometry.verticesNeedUpdate = false;
            geometry.elementsNeedUpdate = false;
            geometry.normalsNeedUpdate = false;
            geometry.colorsNeedUpdate = false;

            geometry.buffersNeedUpdate = false;

        }

    }
    
    function removeObject( object, scene ) {

        if (object instanceof WebMol.Mesh || object instanceof WebMol.Line )
            removeInstances(scene.__webglObjects, object);

        else if (object instanceof WebMol.Sprite)
            removeInstancesDirect(scene.__webglSprites, object);
            
        object.__webglActive = false;

    }

    function removeInstances( objList, object ) {

        for (var o = objList.length - 1; o >= 0; --o) {

            if (objList[o].object === object) 
                objList.splice(o, 1);

        }
    }

    function removeInstancesDirect( objList, object ) {

        for (var o = objList.length - 1; o >= 0; --o) {

            if (objList[o] === object) 
                objList.splice(o, 1);

        }
    }

    function unrollBufferMaterial( globject ) {

        var object = globject.object;
        var material = object.material;

        if ( material.transparent) {                    
            globject.opaque = null;
            globject.transparent = material;
            var blankMaterial = material.clone();
            blankMaterial.opacity = 0;
            globject.blank = blankMaterial;
        }

        else {
            globject.opaque = material;
            globject.transparent = null;
        }

    }

    function setBuffers( geometryGroup, hint, line ) {

        var vertexArray = geometryGroup.__vertexArray;
        var colorArray = geometryGroup.__colorArray;
         
        //vertex buffers
        _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglVertexBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, vertexArray, hint );        

        //color buffers
        _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglColorBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, colorArray, hint );    
              
        
        //normal buffers
        if (geometryGroup.__normalArray !== undefined && geometryGroup.__webglNormalBuffer !== undefined) {
            var normalArray = geometryGroup.__normalArray;
            _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglNormalBuffer );
            _gl.bufferData( _gl.ARRAY_BUFFER, normalArray, hint );       
             
        }
        
        //face (index) buffers
        if (geometryGroup.__faceArray !== undefined && geometryGroup.__webglFaceBuffer !== undefined) {
            var faceArray = geometryGroup.__faceArray;
            _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglFaceBuffer );
            _gl.bufferData( _gl.ELEMENT_ARRAY_BUFFER, faceArray, hint );  
                      
        }
        
        //line (index) buffers (for wireframe)
        if (geometryGroup.__lineArray !== undefined && geometryGroup.__webglLineBuffer !== undefined) {
            var lineArray = geometryGroup.__lineArray;            
            _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglLineBuffer );
            _gl.bufferData( _gl.ELEMENT_ARRAY_BUFFER, lineArray, hint );
        }

    }
    
    //Creates appropriate gl buffers for geometry chunk
    //TODO: do we need line buffer for mesh objects?
    //Also, can we integrate this with createLineBuffers?
    function createMeshBuffers ( geometryGroup ) {

        geometryGroup.__webglVertexBuffer = _gl.createBuffer();
        geometryGroup.__webglNormalBuffer = _gl.createBuffer();
        geometryGroup.__webglColorBuffer = _gl.createBuffer();

        geometryGroup.__webglFaceBuffer = _gl.createBuffer();
        geometryGroup.__webglLineBuffer = _gl.createBuffer();

        _this.info.memory.geometries++;
    }
    
    function createLineBuffers ( geometry ) {
        
        geometry.__webglVertexBuffer = _gl.createBuffer();
        geometry.__webglColorBuffer = _gl.createBuffer();
        
        _this.info.memory.geometries++;
    }

    function addBuffer (objlist, buffer, object) {

        objlist.push(
            {
                buffer: buffer,
                object: object,
                opaque: null,
                transparent: null
            }
        );

    }

    function setupMatrices (object, camera) {

        object._modelViewMatrix.multiplyMatrices( camera.matrixWorldInverse, object.matrixWorld );

        object._normalMatrix.getInverse( object._modelViewMatrix );
        object._normalMatrix.transpose();

    }

    function isPowerOfTwo ( value ) {

        return ( value & ( value - 1 ) ) === 0;

    }
    
    // Fallback filters for non-power-of-2 textures

    function filterFallback ( f ) {

        return _gl.LINEAR;

    }

    function setTextureParameters ( textureType, texture, isImagePowerOfTwo ) {

        if ( isImagePowerOfTwo ) {

            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, paramToGL( texture.wrapS ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, paramToGL( texture.wrapT ) );

            _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, paramToGL( texture.magFilter ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, paramToGL( texture.minFilter ) );

        } else {

            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE );
            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE );

            _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, filterFallback( texture.magFilter ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, filterFallback( texture.minFilter ) );

        }

    }
    
    this.setTexture = function (texture, slot) {

        if (texture.needsUpdate) {

            if ( !texture.__webglInit ) {

                texture.__webglInit = true;

                texture.addEventListener('dispose', onTextureDispose);

                texture.__webglTexture = _gl.createTexture();

                _this.info.memory.textures++;

            }

            _gl.activeTexture(_gl.TEXTURE0 + slot);
            _gl.bindTexture(_gl.TEXTURE_2D, texture.__webglTexture);

            _gl.pixelStorei(_gl.UNPACK_FLIP_Y_WEBGL, texture.flipY);
            _gl.pixelStorei(_gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, texture.premultiplyAlpha);
            _gl.pixelStorei(_gl.UNPACK_ALIGNMENT, texture.unpackAlignment);

            var image = texture.image,
            isImagePowerOfTwo = isPowerOfTwo(image.width) && isPowerOfTwo(image.height),
            glFormat = paramToGL(texture.format),
            glType = paramToGL(texture.type);

            setTextureParameters(_gl.TEXTURE_2D, texture, isImagePowerOfTwo);

            var mipmap, mipmaps = texture.mipmaps;

            // regular Texture (image, video, canvas)

            // use manually created mipmaps if available
            // if there are no manual mipmaps
            // set 0 level mipmap and then use GL to generate other mipmap levels

            if ( mipmaps.length > 0 && isImagePowerOfTwo ) {

                for ( var i = 0, il = mipmaps.length; i < il; i ++ ) {
                    mipmap = mipmaps[ i ];
                    _gl.texImage2D( _gl.TEXTURE_2D, i, glFormat, glFormat, glType, mipmap );
                }
                
                texture.generateMipmaps = false;
            } 
            
            else 
                _gl.texImage2D( _gl.TEXTURE_2D, 0, glFormat, glFormat, glType, texture.image );

            
            if ( texture.generateMipmaps && isImagePowerOfTwo ) _gl.generateMipmap( _gl.TEXTURE_2D );

            texture.needsUpdate = false;

            if ( texture.onUpdate ) texture.onUpdate();

        } else {

            _gl.activeTexture( _gl.TEXTURE0 + slot );
            _gl.bindTexture( _gl.TEXTURE_2D, texture.__webglTexture );

        }

    };
    
    // Map constants to WebGL constants

    function paramToGL ( p ) {

        if ( p === WebMol.UnsignedByteType ) return _gl.UNSIGNED_BYTE;
        if ( p === WebMol.RGBAFormat ) return _gl.RGBA;

        return 0;

    }
    
    function setupLights ( program, lights ) {
        var l, ll, light, n,
        r = 0, g = 0, b = 0,
        color,
        position,
        intensity,
        distance,
        
        zlights = _lights,
        
        dirColors = zlights.directional.colors,
        dirPositions = zlights.directional.positions,
        
        dirCount = 0,
        dirLength = 0,
        dirOffset = 0;
        
        for ( l = 0, ll = lights.length; l < ll; l++) {
            
            light = lights[l];
            
            color = light.color;
            intensity = light.intensity;
            distance = light.distance;
            
            if (light instanceof WebMol.Light) {
                
                dirCount++;
                
                _direction.getPositionFromMatrix(light.matrixWorld);
                _vector3.getPositionFromMatrix(light.target.matrixWorld);
                _direction.sub(_vector3);
                _direction.normalize();
                
                if (_direction.x === 0 && _direction.y === 0 && _direction.z === 0)
                    continue;
                
                dirPositions[dirOffset] = _direction.x;
                dirPositions[dirOffset + 1] = _direction.y;
                dirPositions[dirOffset + 2] = _direction.z;

                dirColors[dirOffset] = color.r * intensity;
                dirColors[dirOffset + 1] = color.g * intensity;
                dirColors[dirOffset + 2] = color.b * intensity;
                
                dirOffset += 3;
                
                dirLength++;    
            }
        
        }

        zlights.ambient[0] = r;
        zlights.ambient[1] = g;
        zlights.ambient[2] = b;
        zlights.directional.length = dirLength;
    }

    function initGL () {

        try {

            if ( ! ( _gl = _canvas.getContext( 'experimental-webgl', { alpha: _alpha, premultipliedAlpha: _premultipliedAlpha, antialias: _antialias, stencil: _stencil, preserveDrawingBuffer: _preserveDrawingBuffer } ) ) ) {

                throw 'Error creating WebGL context.';

            }

        } catch ( error ) {

            console.error( error );

        }

    }

    function setDefaultGLState () {

        _gl.clearColor( 0, 0, 0, 1 );
        _gl.clearDepth( 1 );
        _gl.clearStencil( 0 );

        _gl.enable( _gl.DEPTH_TEST );
        _gl.depthFunc( _gl.LEQUAL );

        _gl.frontFace( _gl.CCW );
        _gl.cullFace( _gl.BACK );
        _gl.enable( _gl.CULL_FACE );

        _gl.enable( _gl.BLEND );
        _gl.blendEquation( _gl.FUNC_ADD );
        _gl.blendFunc( _gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA );

        _gl.clearColor( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha );

    }
    
    this.addPostPlugin(new WebMol.SpritePlugin());
        
};

/*
 * Scene class
 */

WebMol.Scene = function() {
    
    WebMol.Object3D.call(this);
    
    this.fog = null;
    
    //May not need...
    this.overrideMaterial = null;
    
    this.matrixAutoUpdate = false;
    
    this.__objects = [];
    this.__lights = [];
    
    this.__objectsAdded = [];
    this.__objectsRemoved = [];
    
};

WebMol.Scene.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Scene.prototype.__addObject = function(object) {
    
    //Directional Lighting
    if (object instanceof WebMol.Light) {
        
        if (this.__lights.indexOf(object) === -1)
            this.__lights.push(object);
        
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
            
            if (idx !== -1)
                this.__objectsRemoved.splice(i, 1);
                
        }
    }
    
    //Add object's children
    
    for (var i in object.children) 
        this.__addObject(object.children[i]);
    
};

WebMol.Scene.prototype.__removeObject = function(object) {
    
    var idx;
    if (object instanceof WebMol.Light) {
        
        idx = this.__lights.indexOf(object);
        
        if (idx !== -1)
            this.__lights.splice(idx, 1);
            
    }
    
    //Object3D
    else {
        
        idx = this.__objects.indexOf(object);
        
        if (idx !== -1) {
            
            this.__objects.splice(idx, 1);
            this.__objectsRemoved.push(object);
            
            //Check if previously added
            
            var ai = this.__objectsAdded.indexOf(object);
            
            if (ai !== -1) 
                this.__objectsAdded.splice(idx, 1);
                
        }
    
    }
    
    //Remove object's children
    for (var i in object.children)
        this.__removeObject(object.children[i]);
    
};


/*
 * Fog Class
 */


WebMol.Fog = function ( hex, near, far ) {

    this.name = '';

    this.color = new WebMol.Color( hex );

    this.near = ( near !== undefined ) ? near : 1;
    this.far = ( far !== undefined ) ? far : 1000;

};

WebMol.Fog.prototype.clone = function () {

    return new WebMol.Fog( this.color.getHex(), this.near, this.far );

};/* 
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

WebMol.ShaderUtils = {
    
    clone: function ( uniforms_src ) {
        
        var u, p, parameter, parameter_src, uniforms_clone = {};
        
        for (u in uniforms_src) {
            uniforms_clone[u] = {};
            uniforms_clone[u].type = uniforms_src[u].type;
            
            var srcValue = uniforms_src[u].value;
            
            if (srcValue instanceof WebMol.Color)
                uniforms_clone[u].value = srcValue.clone();
            else if (typeof srcValue === "number")
                uniforms_clone[u].value = srcValue;
            else if (srcValue instanceof Array) 
                uniforms_clone[u].value = [];
            else
                console.error("Error copying shader uniforms from ShaderLib: unknown type for uniform");
            
        }
        
        return uniforms_clone;
    }
};

WebMol.ShaderLib = { 
    basic : {
        fragmentShader : [
"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 diffuse;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vColor;",

"void main() {",
    
"    gl_FragColor = vec4( diffuse, opacity );",
"    gl_FragColor = gl_FragColor * vec4( vColor, opacity );",
    
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",    
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
    
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"
                                                     
].join("\n"),
        
        vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",

"attribute vec3 position;",
"attribute vec3 color;",

"varying vec3 vColor;",

"void main() {",

"    vColor = color;",
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
"    gl_Position = projectionMatrix * mvPosition;",

"}"
        
].join("\n"),
    
        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            diffuse: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogColor: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000}
        }

    },
    
    
    lambert : { 
        fragmentShader : [

"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vLightFront;",
"varying vec3 vColor;",

"void main() {",
    
"    gl_FragColor = vec4( vec3 ( 1.0 ), opacity );",
    
"    #ifndef WIREFRAME",
"    gl_FragColor.xyz *= vLightFront;",
"    #endif",
    
"    gl_FragColor = gl_FragColor * vec4( vColor, opacity );",
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",
    
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
    
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"


].join("\n"),
       
       vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 ambient;",
"uniform vec3 diffuse;",
"uniform vec3 emissive;",
"uniform vec3 ambientLightColor;",
"uniform vec3 directionalLightColor[ 1 ];",
"uniform vec3 directionalLightDirection[ 1 ];",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"varying vec3 vColor;",
"varying vec3 vLightFront;",

"void main() {",
    
"    vColor = color;",
    
"    vec3 objectNormal = normal;",  
"    vec3 transformedNormal = normalMatrix * objectNormal;",    
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
    
"    vLightFront = vec3( 0.0 );",
    
"    transformedNormal = normalize( transformedNormal );",
    
"    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );",
"    vec3 dirVector = normalize( lDirection.xyz );",
"    float dotProduct = dot( transformedNormal, dirVector );",
"    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );",
    
"    vLightFront += directionalLightColor[ 0 ] * directionalLightWeighting;",
"    vLightFront = vLightFront * diffuse + ambient * ambientLightColor + emissive;",
    
"    gl_Position = projectionMatrix * mvPosition;",
"}"
           
].join("\n"),

        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            diffuse: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogColor: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000},           
            ambient: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            emissive: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            ambientLightColor: { type: 'fv', value: [] },
            directionalLightColor: { type: 'fv', value: [] },
            directionalLightDirection: { type: 'fv', value: [] }
        }

    },
    
    
    sprite: {
        
        fragmentShader: [
                                                         
"uniform vec3 color;",
"uniform sampler2D map;",
"uniform float opacity;",

"uniform int fogType;",
"uniform vec3 fogColor;",
"uniform float fogDensity;",
"uniform float fogNear;",
"uniform float fogFar;",
"uniform float alphaTest;",

"varying vec2 vUV;",

"void main() {",
    
"    vec4 texture = texture2D(map, vUV);",
    
"    if (texture.a < alphaTest) discard;",
    
"    gl_FragColor = vec4(color * texture.xyz, texture.a * opacity);",
    
"    if (fogType > 0) {",
        
"        float depth = gl_FragCoord.z / gl_FragCoord.w;",
"        float fogFactor = 0.0;",
        
"        if (fogType == 1) {",
"            fogFactor = smoothstep(fogNear, fogFar, depth);",
"        }",
        
"        else {",
"            const float LOG2 = 1.442695;",
"            float fogFactor = exp2(- fogDensity * fogDensity * depth * depth * LOG2);",
"            fogFactor = 1.0 - clamp(fogFactor, 0.0, 1.0);",
"        }",
        
"        gl_FragColor = mix(gl_FragColor, vec4(fogColor, gl_FragColor.w), fogFactor);",
        
"    }",
"}"                                              
            
].join("\n"),
        
        vertexShader: [

"uniform int useScreenCoordinates;",
"uniform int sizeAttenuation;",     
"uniform vec3 screenPosition;",
"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform float rotation;",
"uniform vec2 scale;",
"uniform vec2 alignment;",
"uniform vec2 uvOffset;",
"uniform vec2 uvScale;",

"attribute vec2 position;",
"attribute vec2 uv;",

"varying vec2 vUV;",

"void main() {",
    
"    vUV = uvOffset + uv * uvScale;",
    
"    vec2 alignedPosition = position + alignment;",
    
"    vec2 rotatedPosition;",
"    rotatedPosition.x = ( cos(rotation) * alignedPosition.x - sin(rotation) * alignedPosition.y ) * scale.x;",
"    rotatedPosition.y = ( sin(rotation) * alignedPosition.x + cos(rotation) * alignedPosition.y ) * scale.y;",
    
"    vec4 finalPosition;",
    
"    if(useScreenCoordinates != 0) {",
"        finalPosition = vec4(screenPosition.xy + rotatedPosition, screenPosition.z, 1.0);",
"    }",
    
"    else {",
"        finalPosition = projectionMatrix * modelViewMatrix * vec4(0.0, 0.0, 0.0, 1.0);",
"        finalPosition.xy += rotatedPosition * (sizeAttenuation == 1 ? 1.0 : finalPosition.z);",
"    }",
    
"    gl_Position = finalPosition;",
    
"}"
       
].join("\n"),

        uniforms : {
            
        }
        
    }
    
};