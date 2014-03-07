/* 
 * WebMol Lighting
 */

WebMol.Light = function(hex) {
	
	WebMol.Object3D.call(this);
	
	this.color = new WebMol.CC.cache[hex];
	
};

WebMol.Light.prototype = Object.create(WebMol.Object3D.prototype);
