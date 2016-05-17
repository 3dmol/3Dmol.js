/*
 * Simplified Perspective Camera
 */

/** @constructor */
$3Dmol.Camera = function(fov, aspect, near, far, ortho) {
    
    $3Dmol.Object3D.call(this);
    
    this.fov = fov !== undefined ? fov : 50;
    this.aspect = aspect !== undefined ? aspect : 1;
    this.near = near !== undefined ? near : 0.1;
    this.far = far !== undefined ? far : 2000;

    this.projectionMatrix = new $3Dmol.Matrix4();
    this.projectionMatrixInverse = new $3Dmol.Matrix4();
    this.matrixWorldInverse = new $3Dmol.Matrix4();
    
    var center = this.position.z;
    this.right = center * Math.tan(Math.PI / 180 * fov);
    this.left = -this.right;
    this.top = this.right / this.aspect;
    this.bottom = -this.top;
    
    this.ortho = !!ortho;
    
    this.updateProjectionMatrix();
        
};

//Inherit Object3D's prototyped methods
$3Dmol.Camera.prototype = Object.create($3Dmol.Object3D.prototype);

$3Dmol.Camera.prototype.lookAt = function(vector){
    
    //Why is the parameter order switched (compared to Object3D)?
    this.matrix.lookAt(this.position, vector, this.up);
    
    if (this.rotationAutoUpdate) {    
        
        if (this.useQuaternion === false) 
            this.rotation.setEulerFromRotationMatrix( this.matrix, this.eulerOrder );
        else
            this.quaternion.copy( this.matrix.decompose()[ 1 ] );    
            
    }
    
};

$3Dmol.Camera.prototype.updateProjectionMatrix = function () {

    if(this.ortho) {
        this.projectionMatrix.makeOrthographic( this.left, this.right, this.top, this.bottom, this.near, this.far );
    } else {
        this.projectionMatrix.makePerspective( this.fov, this.aspect, this.near, this.far );
    }
    
    this.projectionMatrixInverse.getInverse(this.projectionMatrix);
};


