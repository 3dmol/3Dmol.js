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

	if ( this.geometry ) {

		if ( ! this.geometry.boundingSphere ) {

			this.geometry.computeBoundingSphere();

		}

	}

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

	if ( this.geometry !== undefined ) {

		if ( this.geometry.boundingSphere === null ) {

			this.geometry.computeBoundingSphere();

		}

		this.updateMorphTargets();

	}

};

WebMol.Mesh.prototype = Object.create( WebMol.Object3D.prototype );

WebMol.Mesh.prototype.updateMorphTargets = function () {

	if ( this.geometry.morphTargets.length > 0 ) {

		this.morphTargetBase = -1;
		this.morphTargetForcedOrder = [];
		this.morphTargetInfluences = [];
		this.morphTargetDictionary = {};

		for ( var m = 0, ml = this.geometry.morphTargets.length; m < ml; m ++ ) {

			this.morphTargetInfluences.push( 0 );
			this.morphTargetDictionary[ this.geometry.morphTargets[ m ].name ] = m;

		}

	}

};

WebMol.Mesh.prototype.getMorphTargetIndexByName = function ( name ) {

	if ( this.morphTargetDictionary[ name ] !== undefined ) {

		return this.morphTargetDictionary[ name ];

	}

	console.log( "THREE.Mesh.getMorphTargetIndexByName: morph target " + name + " does not exist. Returning 0." );

	return 0;

};

WebMol.Mesh.prototype.clone = function ( object ) {

	if ( object === undefined ) object = new WebMol.Mesh( this.geometry, this.material );

	WebMol.Object3D.prototype.clone.call( this, object );

	return object;

};
