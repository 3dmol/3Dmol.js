/*
 * Scene class
 */
/** @constructor */
$3Dmol.Scene = function() {
    
    $3Dmol.Object3D.call(this);
    
    this.fog = null;
    
    //May not need...
    this.overrideMaterial = null;
    
    this.matrixAutoUpdate = false;
    
    this.__objects = [];
    this.__lights = [];
    
    this.__objectsAdded = [];
    this.__objectsRemoved = [];
    
};

$3Dmol.Scene.prototype = Object.create($3Dmol.Object3D.prototype);

$3Dmol.Scene.prototype.__addObject = function(object) {
    
    //Directional Lighting
    if (object instanceof $3Dmol.Light) {
        
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
                this.__objectsRemoved.splice(idx, 1);
                
        }
    }
    
    //Add object's children
    
    for (var i = 0; i < object.children.length; i++) 
        this.__addObject(object.children[i]);
    
};

$3Dmol.Scene.prototype.__removeObject = function(object) {
    
    var idx;
    if (object instanceof $3Dmol.Light) {
        
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
    for (var i = 0; i < object.children.length; i++)
        this.__removeObject(object.children[i]);
    
};


/*
 * Fog Class
 */

/** @constructor */
$3Dmol.Fog = function ( hex, near, far ) {

    this.name = '';

    this.color = new $3Dmol.Color( hex );

    this.near = ( near !== undefined ) ? near : 1;
    this.far = ( far !== undefined ) ? far : 1000;

};

$3Dmol.Fog.prototype.clone = function () {

    return new $3Dmol.Fog( this.color.getHex(), this.near, this.far );

};