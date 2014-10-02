
$3Dmol.Light = function(hex, intensity) {
    
    $3Dmol.Object3D.call(this);
    
    this.color = new $3Dmol.Color(hex);
    this.position = new $3Dmol.Vector3( 0, 1, 0 );
    this.target = new $3Dmol.Object3D();

    this.intensity = ( intensity !== undefined ) ? intensity : 1;

    this.castShadow = false;
    this.onlyShadow = false;
    
};

$3Dmol.Light.prototype = Object.create($3Dmol.Object3D.prototype);
