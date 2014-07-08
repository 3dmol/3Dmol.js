//Intersection sphere and box shapes.  


//Intersection sphere for sphere, stick render
/** @constructor */
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
/** @constructor */
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
/** @constructor */
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


