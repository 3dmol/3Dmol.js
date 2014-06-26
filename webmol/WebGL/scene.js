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
    
    var i;
    
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
            
            i = this.__objectsRemoved.indexOf(object);
            
            if (i !== -1)
                this.__objectsRemoved.splice(i, 1);
                
        }
    }
    
    //Add object's children
    
    for (i= 0, il = object.children.length; i < il; i++) 
        this.__addObject(object.children[i]);
    
};

WebMol.Scene.prototype.__removeObject = function(object) {
    
    var i, il;
    
    if (object instanceof WebMol.Light) {
        
        i = this.__lights.indexOf(object);
        
        if (i !== -1)
            this.__lights.splice(i, 1);
            
    }
    
    //Object3D
    else {
        
        i = this.__objects.indexOf(object);
        
        if (i !== -1) {
            
            this.__objects.splice(i, 1);
            this.__objectsRemoved.push(object);
            
            //Check if previously added
            
            var ai = this.__objectsAdded.indexOf(object);
            
            if (ai !== -1) 
                this.__objectsAdded.splice(i, 1);
                
        }
    
    }
    
    //Remove object's children
    for (i = 0, il = object.children.length; i < il; i++)
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

};