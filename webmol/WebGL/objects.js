/* 
 * WebMol Mesh and Line objects
 */


//Line Object
/** @constructor */
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
/** @constructor */
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
/** @constructor */
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
