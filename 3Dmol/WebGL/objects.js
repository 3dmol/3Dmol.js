/* 
 * $3Dmol Mesh and Line objects
 */


//Line Object
/** @constructor */
$3Dmol.Line = function (geometry, material, type) {

    $3Dmol.Object3D.call(this);

    this.geometry = geometry;
        //TODO: update material and type to webgl
    this.material = (material !== undefined) ? material : new $3Dmol.LineBasicMaterial( { color: Math.random() * 0xffffff } );
    this.type = (type !== undefined) ? type : $3Dmol.LineStrip;

};

$3Dmol.LineStrip = 0;
$3Dmol.LinePieces = 1;

$3Dmol.Line.prototype = Object.create($3Dmol.Object3D.prototype);

$3Dmol.Line.prototype.clone = function (object) {

    if (object === undefined) object = new $3Dmol.Line(this.geometry, this.material, this.type);

    $3Dmol.Object3D.prototype.clone.call(this, object);

    return object;

};


//Mesh Object
/** @constructor */
$3Dmol.Mesh = function(geometry, material) {

    $3Dmol.Object3D.call(this);

    this.geometry = geometry;
    this.material = (material !== undefined) ? material : new $3Dmol.MeshBasicMaterial( { color: Math.random() * 0xffffff, wireframe: true } );

};

$3Dmol.Mesh.prototype = Object.create($3Dmol.Object3D.prototype);

$3Dmol.Mesh.prototype.clone = function (object) {

    if (object === undefined) object = new $3Dmol.Mesh(this.geometry, this.material);

    $3Dmol.Object3D.prototype.clone.call(this, object);

    return object;

};


//Sprite object
/** @constructor */
$3Dmol.Sprite = function(material) {
    
    $3Dmol.Object3D.call(this);
    
    this.material = (material !== undefined) ? material : new $3Dmol.SpriteMaterial();

    this.rotation3d = this.rotation;
    this.rotation = 0;
    
};

$3Dmol.Sprite.prototype = Object.create($3Dmol.Object3D.prototype);

$3Dmol.Sprite.prototype.updateMatrix = function() {
    
    this.matrix.setPosition(this.position);
    
    this.rotation3d.set(0, 0, this.rotation);
    this.matrix.setRotationFromEuler(this.rotation3d);
    
    if (this.scale.x !== 1 || this.scale.y !== 1)
        this.matrix.scale(this.scale);
    
    this.matrixWorldNeedsUpdate = true;
    
};

$3Dmol.Sprite.prototype.clone = function(object) {
    
    if (object === undefined)
        object = new $3Dmol.Sprite(this.material);
    
    $3Dmol.Object3D.prototype.clone.call(this, object);
    
    return object;
    
};
