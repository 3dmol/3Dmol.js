/* 
 * $3Dmol Mesh and Line objects
 */


//Line Object
/** @constructor */
export function Line(geometry, material, type) {

    Object3D.call(this);

    this.geometry = geometry;
        //TODO: update material and type to webgl
    this.material = (material !== undefined) ? material : new LineBasicMaterial( { color: Math.random() * 0xffffff } );
    this.type = (type !== undefined) ? type : LineStrip;

};

export var LineStrip = 0;
export var LinePieces = 1;

Line.prototype = Object.create(Object3D.prototype);

Line.prototype.clone = function (object) {

    if (object === undefined) object = new Line(this.geometry, this.material, this.type);

    Object3D.prototype.clone.call(this, object);

    return object;

};


//Mesh Object
/** @constructor */
export function Mesh(geometry, material) {

    Object3D.call(this);

    this.geometry = geometry;
    this.material = (material !== undefined) ? material : new MeshBasicMaterial( { color: Math.random() * 0xffffff, wireframe: true } );

};

Mesh.prototype = Object.create(Object3D.prototype);

Mesh.prototype.clone = function (object) {

    if (object === undefined) object = new Mesh(this.geometry, this.material);

    Object3D.prototype.clone.call(this, object);

    return object;

};


//Sprite object
/** @constructor */
export function Sprite(material) {
    
    Object3D.call(this);
    
    this.material = (material !== undefined) ? material : new SpriteMaterial();

    this.rotation3d = this.rotation;
    this.rotation = 0;
    
};

Sprite.prototype = Object.create(Object3D.prototype);

Sprite.prototype.updateMatrix = function() {
    
    this.matrix.setPosition(this.position);
    
    this.rotation3d.set(0, 0, this.rotation);
    this.matrix.setRotationFromEuler(this.rotation3d);
    
    if (this.scale.x !== 1 || this.scale.y !== 1)
        this.matrix.scale(this.scale);
    
    this.matrixWorldNeedsUpdate = true;
    
};

Sprite.prototype.clone = function(object) {
    
    if (object === undefined)
        object = new Sprite(this.material);
    
    Object3D.prototype.clone.call(this, object);
    
    return object;
    
};
