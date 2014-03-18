/* 
 * WebMol Lighting
 */

//TODO: Strip down this class - do I really use all of these instance variables?
WebMol.Light = function(hex, intensity) {
    
    WebMol.Object3D.call(this);
    
    this.color = new WebMol.Color(hex);
    this.position = new WebMol.Vector3( 0, 1, 0 );
    this.target = new WebMol.Object3D();

    this.intensity = ( intensity !== undefined ) ? intensity : 1;

    this.castShadow = false;
    this.onlyShadow = false;
    
};

WebMol.Light.prototype = Object.create(WebMol.Object3D.prototype);
