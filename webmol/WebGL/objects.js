/* 
 * WebMol Mesh and Line objects
 */

//Line Object
WebMol.Line = function ( geometry, material, type ) {

	WebMol.Object3D.call( this );

	this.geometry = geometry;
        //TODO: update material and type to webgl
	this.material = ( material !== undefined ) ? material : new THREE.LineBasicMaterial( { color: Math.random() * 0xffffff } );
	this.type = ( type !== undefined ) ? type : THREE.LineStrip;

};

WebMol.LineStrip = 0;
WebMol.LinePieces = 1;

WebMol.Line.prototype = Object.create( WebMol.Object3D.prototype );

WebMol.Line.prototype.clone = function ( object ) {

	if ( object === undefined ) object = new WebMol.Line( this.geometry, this.material, this.type );

	WebMol.Object3D.prototype.clone.call( this, object );

	return object;

};

//Mesh Object
WebMol.Mesh = function ( geometry, material ) {

	WebMol.Object3D.call( this );

	this.geometry = geometry;
	this.material = ( material !== undefined ) ? material : new THREE.MeshBasicMaterial( { color: Math.random() * 0xffffff, wireframe: true } );

};

WebMol.Mesh.prototype = Object.create( WebMol.Object3D.prototype );

WebMol.Mesh.prototype.clone = function ( object ) {

	if ( object === undefined ) object = new WebMol.Mesh( this.geometry, this.material );

	WebMol.Object3D.prototype.clone.call( this, object );

	return object;

};
