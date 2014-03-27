//Intersection sphere and box shapes.  
//May want to extend this to other uses (e.g. drawing pharmacophores?


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

WebMol.Box3 = function(min, max) {
    
    this.min = (min !== undefined) ?
        min : new WebMol.Vector3(Infinity, Infinity, Infinity);
        
    this.max = (max !== undefined) ?
        max : new WebMol.Vector3(-Infinity, -Infinity, -Infinity);
        
};

//Axis aligned bounding box
WebMol.Box3.prototype = {

    constructor : WebMol.Box3,

    set : function(min, max) {
        
        this.min.copy(min);
        this.max.copy(max);
        
        return this;
        
    }




    
};


//Intersection bounding box for line, stick (cylinder), ribbon render
