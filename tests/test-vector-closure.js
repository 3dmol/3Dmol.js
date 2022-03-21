const $3Dmol = {};

/** @constructor */
$3Dmol.Vector3 = function(x, y, z) {
    this.x = x || 0.0;
    this.y = y || 0.0;
    this.z = z || 0.0;
};

$3Dmol.Vector3.prototype =  {
    
    constructor : Vector3,
    
    set(x, y, z) {
        
        this.x = x;
        this.y = y;
        this.z = z;
        
        return this;
    },
    
    copy(v) {
        
        this.x = v.x;
        this.y = v.y;
        this.z = v.z;
        
        return this;  
    },
    
    add(v) {
        
        this.x += v.x;
        this.y += v.y;
        this.z += v.z;  
        
        return this;
    },
    
    addVectors(a, b) {
        
        this.x = a.x + b.x;
        this.y = a.y + b.y;
        this.z = a.z + b.z;
        
        return this;
    },
    
    sub(v) {
        
        this.x -= v.x;
        this.y -= v.y;
        this.z -= v.z;
        
        return this;
    },
    
    subVectors(a, b) {
        
        this.x = a.x - b.x;
        this.y = a.y - b.y;
        this.z = a.z - b.z;
        
        return this;
    },
    
    multiplyScalar(s) {
        
        this.x *= s;
        this.y *= s;
        this.z *= s;
        
        return this;
    },
    
    divideScalar(s) {
        
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
    

    distanceTo(v) {
        return Math.sqrt(this.distanceToSquared(v));
    },

    distanceToSquared(v) {
        const dx = this.x - v.x;
        const dy = this.y - v.y;
        const dz = this.z - v.z;

        return dx * dx + dy * dy + dz * dz;
    },
    
    applyMatrix4(m) {
    
        const {x} = this; const {y} = this; const {z} = this;
        
        const e = m.elements;
        
        this.x = e[0]*x + e[4]*y + e[8]*z + e[12];
        this.y = e[1]*x + e[5]*y + e[9]*z + e[13];
        this.z = e[2]*x + e[6]*y + e[10]*z + e[14];
        
        return this;
    },
    
    applyProjection(m) {
        
        // input: $3Dmol.Matrix4 projection matrix
        
        const {x} = this; const {y} = this; const {z} = this;
        
        const e = m.elements;
        const d = ( e[3]*x + e[7]*y + e[11]*z + e[15]);
        
        this.x = (e[0]*x + e[4]*y + e[8]*z + e[12]) / d;
        this.y = (e[1]*x + e[5]*y + e[9]*z + e[13]) / d;
        this.z = (e[2]*x + e[6]*y + e[10]*z + e[14]) / d;
        
        return this;
    },
    
    applyQuaternion(q) { 
        
        const {x} = this;
        const {y} = this;
        const {z} = this;
        
        const qx = q.x;
        const qy = q.y;
        const qz = q.z;
        const qw = q.w;
        
        // calculate quaternion * vector
        
        const ix = qw * x + qy * z - qz * y;
        const iy = qw * y + qz * x - qx * z;
        const iz = qw * z + qx * y - qy * x;
        const iw = -qw * x - qy * y - qz * z;
        
        // calculate result * inverse quaternion
        
        this.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
        this.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
        this.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;
        
        return this;
    },
    
    negate() {
        
        return this.multiplyScalar(-1);
    },
    
    dot(v) {
        
        return this.x * v.x + this.y * v.y + this.z * v.z;
    },
    
    length() {
        
        return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    },
    
    lengthSq() {
    
        return (this.x*this.x + this.y*this.y + this.z*this.z);
    },
    
    normalize() {
        
        return this.divideScalar( this.length() );
    },
    
    cross (v) {
        
        const {x} = this; const {y} = this; const {z} = this;
        
        this.x = y * v.z - z * v.y;
        this.y = z * v.x - x * v.z;
        this.z = x * v.y - y * v.x;
        
        return this;
    },
    
    crossVectors(a, b) {
        
        this.x = a.y * b.z - a.z * b.y;
        this.y = a.z * b.x - a.x * b.z;
        this.z = a.x * b.y - a.y * b.x;
        
        return this;
    },
    
    getPositionFromMatrix(m) {
        
        this.x = m.elements[12];
        this.y = m.elements[13];
        this.z = m.elements[14];
        
        return this;
    },

    setEulerFromRotationMatrix (m, order) {

        // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

        const te = m.elements;
        const m11 = te[0]; const m12 = te[4]; const m13 = te[8];
        const m21 = te[1]; const m22 = te[5]; const m23 = te[9];
        const m31 = te[2]; const m32 = te[6]; const m33 = te[10];

        if ( order === undefined || order === 'XYZ' ) {

            this.y = Math.asin( $3Dmol.Math.clamp( m13, -1, 1 ) );

            if ( Math.abs( m13 ) < 0.99999 ) {

                this.x = Math.atan2( - m23, m33 );
                this.z = Math.atan2( - m12, m11 );

            } else {

                this.x = Math.atan2( m32, m22 );
                this.z = 0;

            }
        }
        
        else {
            console.error(`Error with vector's setEulerFromRotationMatrix: Unknown order: ${  order}`);
        }
        
        return this;

    },
    
    clone() {
        return new $3Dmol.Vector3(this.x, this.y, this.z);
    }
    
};

