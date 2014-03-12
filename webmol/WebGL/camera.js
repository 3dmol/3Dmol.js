/*
 * Simplified Perspective Camera
 */


WebMol.Camera = function(fov, aspect, near, far) {
    
    //Invoke Object3D constructor function (i.e. initialize the same instance variables)
    WebMol.Object3D.call(this);
    
    this.fov = fov !== undefined ? fov : 50;
    this.aspect = aspect !== undefined ? aspect : 1;
    this.near = near !== undefined ? near : 0.1;
    this.far = far !== undefined ? far : 2000;
    
    //TODO: Replace these with own matrix class
    this.projectionMatrix = new THREE.Matrix4();
    this.matrixWorldInverse = new THREE.Matrix4();
    
    this.updateProjectionMatrix();
        
};

//Inherit Object3D's prototyped methods
WebMol.Camera.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Camera.prototype.lookAt = function(vector){
    
    //Why is the parameter order switched (compared to Object3D)?
    this.matrix.lookAt(this.position, vector, this.up);
    
    if (this.rotationAutoUpdate) {    
        
        if (this.useQuaternion === false) 
            this.rotation.setEulerFromRotationMatrix( this.matrix, this.eulerOrder );
        else
            this.quaternion.copy( this.matrix.decompose()[ 1 ] );    
            
    }
    
};

WebMol.Camera.prototype.updateProjectionMatrix = function () {

    this.projectionMatrix.makePerspective( this.fov, this.aspect, this.near, this.far );

};


