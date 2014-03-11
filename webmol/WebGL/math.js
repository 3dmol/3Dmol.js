/* 
 * math-like functionality 
 * quaternion, vector, matrix
 */

// Quaternion

WebMol.Quaternion = function(x, y, z, w) {

    this.x = x || 0;
    this.y = y || 0;
    this.z = z || 0;
    this.w = (w !== undefined) ? w : 1;
    
};

WebMol.extend(WebMol.Quaternion.prototype, {
    
    set: function(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
        
        return this;
    },
    
    copy: function(q) {
        this.x = q.x;
        this.y = q.y;
        this.z = q.z
        this.w = q.w;
        
        return this;
    },
    
    conjugate: function() {
        this.x *= -1;
        this.y *= -1;
        this.z *= -1;
        
        return this;
    },
    
    inverse: function() {
        return this.conjugate.normalize();  
    },
    
    length: function() {
        return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z + this.w*this.w);
    },
    
    normalize: function() {
        var l = this.length();
        
        if (l === 0) {
            this.x = 0;
            this.y = 0;
            this.z = 0;
            this.w = 1;
        }
        
        else {
            l = 1 / l;
            
            this.x *= l;
            this.y *= l;
            this.z *= l;
            this.w *= l;
        }
        
        return this;
        
    },
    
    multiply: function(q) {
        return this.multiplyQuaternions(this, q);
    },
    
    multiplyQuaternions: function(a, b) {
        
        var qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
        var qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;
        
        this.x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
        this.y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
        this.z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
        this.w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;
        
    }
    
});

WebMol.Vertex = function(x, y, z) {
	this.x = x || 0.0;
	this.y = y || 0.0;
	this.z = z || 0.0;
};

WebMol.Vertex.prototype = {
	
	constructor : WebMol.Vertex,
	
	add : function( x, y, z ) {
		this.x += x;
		this.y += y;
		this.z += z;
	},
	
	sub : function( x, y, z ) {
		this.x -= x;
		this.y -= y;
		this.z -= z;
	},
	
	normalize : function() {
		var normFactor = 1 / Math.sqrt( (this.x * this.x) + (this.y * this.y) + (this.z * this.z) );
		
		this.multiplyByScalar(normFactor);
	},
	
	multiplyByScalar : function(s) {
		this.x *= s;
		this.y *= s;
		this.z *= s;
	},
	
	clone : function() {
		return new WebMol.Vertex(this.x, this.y, this.z);
	}
	
};
