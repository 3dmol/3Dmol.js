/* core Object3D
 * Base class for Scene, Camera, Geometry
 * Geometry class
 */

var WebMol = WebMol || {};

//Event Handling
WebMol.EventDispatcher = function() {
  
    var listeners = {};
    
    this.addEventListener = function(type, listener) {
        if (listeners[type] === undefined)
            listeners[type] = [];
        
        if (listeners[type].indexOf(listener) === -1)
            listeners[type].push(listener);
    };  
    
    this.removeEventListener = function(type, listener) {
        
        var index = listeners[type].indexOf(listener);
        
        if (index !== -1)
            listeners[type].splice(index, 1);
              
    };
    
    this.dispatchEvent = function(event) {
        
        var listenerArray = listeners[event.type];
        
        if (listenerArray !== undefined) {
            event.target = this;
            
            for (var i = 0, l = listenerArray.length; i < l; i++)
                listenerArray[i].call(this, event);
                
        }
            
    };
    
};


//Object3D base constructor function
WebMol.Object3D = function() {
    
    this.id = WebMol.Object3DIDCount++;
    
    this.name = "";
    
    this.parent = undefined;
    this.children = [];
    
    this.position = new WebMol.Vector3();
    this.rotation = new WebMol.Vector3();
    this.matrix = new WebMol.Matrix4();
    this.matrixWorld = new WebMol.Matrix4();
    this.quaternion = new WebMol.Quaternion();
    this.eulerOrder = 'XYZ';
    
    this.up = new WebMol.Vector3(0, 1, 0);
    this.scale = new WebMol.Vector3(1, 1, 1);
    
    this.matrixAutoUpdate = true;
    this.matrixWorldNeedsUpdate = true;
    this.rotationAutoUpdate = true;
    this.useQuaternion = false;
    
    this.visible = true;
    
};

WebMol.Object3D.prototype = {
    
    constructor : WebMol.Object3D,
    
    lookAt : function(vector) {
        
        this.matrix.lookAt(vector, this.position, this.up);
        
        if (this.rotationAutoUpdate) {
            
            if (this.useQuaternion === true) 
                this.quaternion.copy(this.matrix.decompose()[1]);
            else
                this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
        }
    },
    
    //add child object
    add : function(object) {
        if (object === this){
            console.error("Can't add WebMol.Object3D to itself");
            return;
        }
        
        object.parent = this;
        this.children.push(object);
        
        //add to the scene (i.e. follow up this instance's parents until reach the top)
        
        var scene = this;
        
        while (scene.parent !== undefined)
            scene = scene.parent;
            
        if (scene !== undefined && scene instanceof WebMol.Scene) 
            scene.__addObject(object);
        
    },
    
    remove : function(object) {
        
        var index = this.children.indexOf(object);
        
        if (index !== -1) {
            
            object.parent = undefined;
            this.children.splice(index, 1);
            
            //Remove from scene
            
            var scene = this;
            
            while (scene.parent !== undefined)
                scene = scene.parent;
                
            if (scene !== undefined && scene instanceof WebMol.Scene)
                scene.__removeObject(object);
                
        }
    },
    
    updateMatrix : function() {
        
        this.matrix.setPosition(this.position);
        
        if (this.useQuaternion === false) 
            this.matrix.setRotationFromEuler(this.rotation, this.eulerOrder);
        else
            this.matrix.setRotationFromQuaternion(this.quaternion);
        
        //TODO: Do I need this??
        if (this.scale.x !== 1 || this.scale.y !== 1 || this.scale.z !== 1)
            this.matrix.scale(this.scale);
            
        this.matrixWorldNeedsUpdate = true;
        
    },
    
    updateMatrixWorld : function(force) {
        
        if (this.matrixAutoUpdate === true) 
            this.updateMatrix();
        
        if (this.matrixWorldNeedsUpdate === true || force === true) {
            
            if (this.parent === undefined)
                this.matrixWorld.copy(this.matrix);
            else
                this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
                
        }
        
        this.matrixWorldNeedsUpdate = false;
        
        //Update matrices of all children
        for (var i in this.children) {
            this.children[i].updateMatrixWorld(true);
        }
    },
    
    clone : function(object) {
        
        if (object === undefined)
            object = new WebMol.Object3D();
            
        object.name = this.name;
        
        object.up.copy(this.up);
        object.position.copy(this.position);
        object.rotation.copy(this.rotation);
        object.eulerOrder = this.eulerOrder;
        object.scale.copy(this.scale);
        
        object.rotationAutoUpdate = this.rotationAutoUpdate;
        object.matrix.copy(this.matrix);
        object.matrixWorld.copy(this.matrixWorld);
        object.quaternion.copy(this.quaternion);
        
        object.useQuaternion = this.useQuaternion;
        
        object.visible = this.visible;
        
        for (var i in this.children) {
            var child = this.children[i];
            object.add(child.clone());
        }
        
        return object;
        
    }
    
};

WebMol.Object3DIDCount = 0;

//Geometry class
//TODO: What can I remove - how can I optimize ?
WebMol.Geometry = (function() {

    //return truncated typed array, including its buffer
    // type == 0 => Uint16Array; type == 1 => Float32Array
    //TODO: Should integrate this directly into geometryGroup's truncateArrayBuffers method
    var truncateArrayBuffer = function(arr, type, end) {
        
        if (arr === null || arr === undefined) {
            return (type === 0) ? new Uint16Array() : new Float32Array();
        }
        
        if (type === 0)
            return new Uint16Array(arr.buffer.slice(arr.byteOffset, end*2));
        else if (type === 1) 
            return new Float32Array(arr.buffer.slice(arr.byteOffset, end*4));
    };
    
    
    var geometryGroup = function(id) {
        this.id = id || 0;
        this.__vertexArray = null;
        this.__colorArray = null;
        this.__normalArray = null;
        this.__faceArray = null;
        this.__lineArray = null;
        this.vertices = 0;
        this.faceidx = 0;
        this.lineidx = 0;
    };
    
    geometryGroup.prototype.getCentroid = function() {
        
        var centroid = new WebMol.Vector3();
        var offset, x, y, z;
        
        for (var i = 0; i < this.vertices; ++i) {
            offset = i*3;
            
            x = this.__vertexArray[offset], y = this.__vertexArray[offset+1], z = this.__vertexArray[offset+2];
            
            centroid.x += x, centroid.y += y, centroid.z += z;
        }
        
        //divideScalar checks for 0 denom
        centroid.divideScalar(this.vertices);
        
        return centroid;
    };
    
    //setup normals - vertex and face array must exist
    geometryGroup.prototype.setNormals = function() {        
        
        var faces = this.__faceArray;
        var verts = this.__vertexArray;
        var norms = this.__normalArray;
        
        if (! this.vertices || ! this.faceidx) 
            return;
        
        //vertex indices
        var a, b, c, d,
        //and actual vertices
        vA, vB, vC, norm;
            
        for (var i = 0; i < faces.length / 3; ++i) {
            
            a = faces[i * 3] * 3;
            b = faces[i * 3 + 1] * 3;
            c = faces[i * 3 + 2] * 3;
            
            vA = new WebMol.Vector3(verts[a], verts[a+1], verts[a+2]);
            vB = new WebMol.Vector3(verts[b], verts[b+1], verts[b+2]);
            vC = new WebMol.Vector3(verts[c], verts[c+1], verts[c+2]);
            
            vA.subVectors(vA, vB);
            vC.subVectors(vC, vB);
            vA.cross(vC);
            
            //face normal
            norm = vA;
            norm.normalize();
            
            norms[a] += norm.x, norms[b] += norm.x, norms[c] += norm.x;
            norms[a + 1] += norm.y, norms[b + 1] += norm.y, norms[c + 1] += norm.y;
            norms[a + 2] += norm.z, norms[b + 2] += norm.z, norms[c + 2] += norm.z;
            
        }             
                
    };
    
    //sets line index array from face arr
    //Note - assumes all faces are triangles (i.e. there will
    //be an extra diagonal for four-sided faces - user should 
    //specify linearr for custom shape generation to show wireframe squares
    //as rectangles rather than two triangles)
    geometryGroup.prototype.setLineIndices = function() {
        
        if (! this.faceidx)
            return;
                    
        var faceArr = this.__faceArray, lineArr = this.__lineArray = new Uint16Array(this.faceidx*2);      
        this.lineidx = this.faceidx*2;         
        var faceoffset;
            
        for (var i = 0; i < this.faceidx; ++i) {
            
            faceoffset = i*3, lineoffset = faceoffset*2;          
            var a = faceArr[faceoffset], b = faceArr[faceoffset+1], c = faceArr[faceoffset+2];
            
            lineArr[lineoffset] = a, lineArr[lineoffset+1] = b;
            lineArr[lineoffset+2] = a, lineArr[lineoffset+3] = c;
            lineArr[lineoffset+4] = b, lineArr[lineoffset+5] = c;
            
        }
    };
    
    geometryGroup.prototype.truncateArrayBuffers = function(mesh) {
        
        var mesh = (mesh === true) ? true : false;
        
        var vertexArr = this.__vertexArray,
            colorArr = this.__colorArray,
            normalArr = this.__normalArray,
            faceArr = this.__faceArray,
            lineArr = this.__lineArray;
                       
        this.__vertexArray = truncateArrayBuffer(vertexArr, 1, this.vertices*3);
        this.__colorArray = truncateArrayBuffer(colorArr, 1, this.vertices*3);
        
        if (mesh) {
            this.__normalArray = truncateArrayBuffer(normalArr, 1, this.vertices*3);
            this.__faceArray = truncateArrayBuffer(faceArr, 0, this.faceidx);
            this.__lineArray = truncateArrayBuffer(lineArr, 0, this.lineidx);
        }
        else {
            this.__normalArray = truncateArrayBuffer(normalArr, 1, 0);
            this.__faceArray = truncateArrayBuffer(faceArr, 0, 0);
            this.__lineArray = truncateArrayBuffer(lineArr, 0, 0);            
        }
        
        this.__inittedArrays = true;        
        
    };
    
    var addGroup = function(geo) {
        var ret = new geometryGroup(geo.geometryGroups.length);
        geo.geometryGroups.push(ret);
        geo.groups = geo.geometryGroups.length;
        
        ret.__vertexArray = new Float32Array(65535*3);
        ret.__colorArray = new Float32Array(65535*3);
        
        //TODO: instantiating uint arrays according to max number of vertices
        // is dangerous, since there exists the possibility that there will be 
        // more face or line indices than vertex points - but so far that doesn't
        // seem to be the case for any of the renders 
        if (geo.mesh) {
            ret.__normalArray = new Float32Array(65535*3);
            ret.__faceArray = new Uint16Array(65535*6);
            ret.__lineArray = new Uint16Array(65535*6);
        }
        
        
        return ret;
    };
        
    Geometry = function(mesh) {
        
        WebMol.EventDispatcher.call(this);
        
        this.id = WebMol.GeometryIDCount++;
    
        this.name = '';
    
        this.hasTangents = false;
    
        this.dynamic = true; // the intermediate typed arrays will be deleted when set to false
        this.mesh = (mesh === true) ? true : false; // Does this geometry represent a mesh (i.e. do we need Face/Line index buffers?)
        // update flags
    
        this.verticesNeedUpdate = false;
        this.elementsNeedUpdate = false;
        this.normalsNeedUpdate = false;
        this.colorsNeedUpdate = false;
    
        this.buffersNeedUpdate = false;
        
        this.geometryGroups = [];
        this.groups = 0;
        
    };
    
    
    Geometry.prototype = {
        
        constructor : Geometry,

        //Get geometry group to accomodate addVertices new vertices - create 
        // new group if necessary       
        updateGeoGroup : function(addVertices) {
        
            addVertices = addVertices || 0;
            
            var retGroup = this.groups > 0 ? this.geometryGroups[ this.groups - 1 ] : null;
            
            if (!retGroup || retGroup.vertices + addVertices > 65535) 
                retGroup = addGroup(this);
                
            return retGroup;
            
        },
        
        addGeoGroup : function() {
            return addGroup(this);  
        },
        
        setUpNormals : function(three) {
            
            var three = three || false;
            
            for ( var g in this.geometryGroups ) {
            
                var geoGroup = this.geometryGroups[g];            
                
                geoGroup.setNormals(three);
                
            }  
                      
        },
        
        //After vertices, colors, etc are collected in regular or typed arrays,
        // either create typed arrays from regular arrays if they don't already exist,
        // or shorten last typed array
        initTypedArrays : function() {
                
            for (var g in this.geometryGroups) {
                
                var group = this.geometryGroups[g];
                
                if (group.__inittedArrays === true)
                    continue;
                
                group.truncateArrayBuffers(this.mesh);
            }
            
        
        },
        
        dispose : function() {
            this.dispatchEvent( {type : 'dispose'} );
        }
    };

    
    return Geometry;
    
})();

Object.defineProperty(WebMol.Geometry.prototype, "vertices", {
    
    get : function() {
        var vertices = 0;
        for (g in this.geometryGroups)
            vertices += this.geometryGroups[g].vertices;
            
        return vertices;
    } 
        
});

WebMol.GeometryIDCount = 0;


//Raycaster

WebMol.Raycaster = (function() {
    
    function Raycaster(origin, direction, far, near) {
        
        this.ray = new WebMol.Ray(origin, direction);
        
        if (this.ray.direction.lengthSq() > 0) 
            this.ray.direction.normalize();
        
        this.near = near || 0;
        this.far = far || Infinity;
    
    };
    
    var sphere = new WebMol.Sphere();
    var cylinder = new WebMol.Cylinder();
    var triangle = new WebMol.Triangle();
    var w_0 = new WebMol.Vector3(); // for cylinders, cylinder.c1 - ray.origin
    var v1 = new WebMol.Vector3(); // all purpose local vector
    var v2 = new WebMol.Vector3();
    var v3 = new WebMol.Vector3();
    //var facePlane = new WebMol.Plane();
    var localRay = new WebMol.Ray();
    var intersectPoint = new WebMol.Vector3();
    var matrixPosition = new WebMol.Vector3();
    
    var inverseMatrix = new WebMol.Matrix4();
        
    var descSort = function(a, b) {
        return a.distance - b.distance;
    };

    // [-1, 1]
    var clamp = function(x) {
        return Math.min(Math.max(x, -1), 1);
    };
    
    //object is a Sphere or (Bounding) Box
    var intersectObject = function(group, clickable, raycaster, intersects) {
        
        matrixPosition.getPositionFromMatrix(group.matrixWorld);
        
        if ((clickable.clickable !== true) || (clickable.intersectionShape === undefined))
            return intersects;
        
        var intersectionShape = clickable.intersectionShape;
        var precision = raycaster.linePrecision;
        precision *= group.matrixWorld.getMaxScaleOnAxis();
        var precisionSq = precision*precision;

        //Check for intersection with clickable's bounding sphere, if it exists
        if (clickable.boundingSphere !== undefined && clickable.boundingSphere instanceof WebMol.Sphere) {
            sphere.copy(clickable.boundingSphere);
            sphere.applyMatrix4(group.matrixWorld);
            
            if (!raycaster.ray.isIntersectionSphere(sphere)) {
                return intersects;
            }
        }
        
        //triangle faces
        for (var i in intersectionShape.triangle) {
            
            if (intersectionShape.triangle[i] instanceof WebMol.Triangle) {
                
                triangle.copy(intersectionShape.triangle[i]);
                triangle.applyMatrix4(group.matrixWorld);
                
                var norm = triangle.getNormal();
                
                var normProj = raycaster.ray.direction.dot(norm);
                
                //face culling
                if (normProj >= 0)
                    continue;
                
                w_0.subVectors(triangle.a, raycaster.ray.origin);
                
                var distance = (norm.dot(w_0)) / normProj;
                
                if (distance < 0)
                    continue;
                    
                //intersects with plane, check if P inside triangle
                v1.copy(raycaster.ray.direction).multiplyScalar(distance).add(raycaster.ray.origin);
                v1.sub(triangle.a); // from pt a to intersection point P
                
                v2.copy(triangle.b).sub(triangle.a); // from pt a to b
                v3.copy(triangle.c).sub(triangle.a); // from pt a to c
                var b_dot_c = v2.dot(v3);
                var b_sq = v2.lengthSq();
                var c_sq = v3.lengthSq();
                
                // P = A + s(v2) + t(v3), inside trianle if 0 <= s, t <=1  and (s + t) <=0
                var s, t;
                
                t = ( b_sq*v1.dot(v3) - b_dot_c*v1.dot(v2) ) / ( b_sq*c_sq - b_dot_c*b_dot_c );
                
                if (t < 0 || t > 1)
                    continue;
                
                s = ( v1.dot(v2) - t*b_dot_c ) / b_sq;
                
                if ( (s < 0 || s > 1) || s + t > 1)
                    continue;
                    
                else
                    intersects.push({clickable : clickable,
                                     distance : distance});  
            }
        }
        
        //cylinders
        for (var i in intersectionShape.cylinder) {
            
            if (intersectionShape.cylinder[i] instanceof WebMol.Cylinder){
                
                cylinder.copy(intersectionShape.cylinder[i]);
                cylinder.applyMatrix4(group.matrixWorld);
                
                w_0.subVectors(cylinder.c1, raycaster.ray.origin); 
                
                var cylProj = w_0.dot(cylinder.direction); // Dela
                var rayProj = w_0.dot(raycaster.ray.direction); // Epsilon
                
                var normProj = clamp(raycaster.ray.direction.dot(cylinder.direction)); // Beta
                
                var denom = 1 - normProj*normProj;
                
                if (denom === 0.0)
                    continue;
                
                var s_c = (normProj*rayProj - cylProj) / denom;
                var t_c = (rayProj - normProj*cylProj) / denom;
                
                v1.copy(cylinder.direction).multiplyScalar(s_c).add(cylinder.c1);  // Q_c
                v2.copy(raycaster.ray.direction).multiplyScalar(t_c).add(raycaster.ray.origin); // P_c
                
                var closestDistSq = v3.subVectors(v1, v2).lengthSq();
                var radiusSq = cylinder.radius*cylinder.radius;
                
                //Smoothing?
                //if (closestDistSq > radiusSq) radiusSq += precisionSq;
                
                // closest distance between ray and cylinder axis not greater than cylinder radius;
                // might intersect this cylinder between atom and bond midpoint
                if (closestDistSq <= radiusSq){
                    var distance;
                    
                    //Find points where ray intersects sides of cylinder
                    var discriminant = (normProj*cylProj - rayProj)*(normProj*cylProj - rayProj) - 
                            denom*(w_0.lengthSq() - cylProj*cylProj - radiusSq);
                    
                    var t;
                    // ray tangent to cylinder?
                    if (discriminant <= 0)
                        t = distance = Math.sqrt(closestDistSq);
                    else
                        t = distance = ( (rayProj - normProj*cylProj) - Math.sqrt(discriminant) ) / denom; 
                    
                    //find closest intersection point; make sure it's between atom's position and cylinder midpoint
                    
                                          
                    var s = normProj*t - cylProj;
                    
                    //does not intersect cylinder between atom and midpoint,
                    // or intersects cylinder behind camera
                    if (s < 0 || s*s > cylinder.lengthSq() || t < 0)
                        continue;
                    
                    else
                        intersects.push({clickable : clickable,
                                         distance : distance});
                    
                }
                    
                
            }
            
        }
         
        //lines
        for (var i = 0, il = intersectionShape.line.length; i < il; i += 2) {
            
            v1.copy(intersectionShape.line[i]);
            v1.applyMatrix4(group.matrixWorld);
            v2.copy(intersectionShape.line[i+1]);
            v2.applyMatrix4(group.matrixWorld);
            
            v3.subVectors(v2, v1);
            var bondLengthSq = v3.lengthSq();
            v3.normalize();
            
            w_0.subVectors(v1, raycaster.ray.origin);
            
            var lineProj = w_0.dot(v3);
            var rayProj = w_0.dot(raycaster.ray.direction);
            
            var normProj = clamp(raycaster.ray.direction.dot(v3));
            
            var denom = 1 - normProj*normProj;
            
            if (denom === 0.0)
                continue;
            
            var s_c = (normProj*rayProj - lineProj) / denom;
            var t_c = (rayProj - normProj*lineProj) / denom;
            
            v1.add(v3.multiplyScalar(s_c)); // Q_c
            v2.copy(raycaster.ray.direction).multiplyScalar(t_c).add(raycaster.ray.origin); // P_c
            
            var closestDistSq = v3.subVectors(v2, v1).lengthSq();
            
            if (closestDistSq < precisionSq && s_c*s_c < bondLengthSq)
                intersects.push({clickable : clickable,
                                 distance : t_c
                                });
            
        }

        for (var i = 0; i < intersectionShape.sphere.length; i++) {
            //sphere
            if (intersectionShape.sphere[i] instanceof WebMol.Sphere) {
                
                sphere.copy(intersectionShape.sphere[i]);
                sphere.applyMatrix4(group.matrixWorld);
                
                if (raycaster.ray.isIntersectionSphere(sphere)) {
                    
                    var distance;
                    
                    v1.subVectors(sphere.center, raycaster.ray.origin);
                    
                    //distance from ray origin to point on the ray normal to sphere's center
                    //must be less than sphere's radius (since ray intersects sphere)
                    var distanceToCenter = v1.dot(raycaster.ray.direction);
                    
                    var discriminant = distanceToCenter*distanceToCenter - (v1.lengthSq() - sphere.radius*sphere.radius);
                    
                    //Don't select if sphere center behind camera
                    if (distanceToCenter < 0) 
                        return intersects;
                    
                    //ray tangent to sphere?
                    if (discriminant <= 0)
                        distance = distanceToCenter;
                    
                    //This is reversed if sphere is closer than ray origin.  Do we have 
                    //to worry about handling that case?
                    else 
                        distance = distanceToCenter - Math.sqrt(discriminant);
    
                    intersects.push({clickable : clickable, 
                                     distance : distance});
                    return intersects;
                }
            }        
       }
        
    };   
       
    Raycaster.prototype.precision = 0.0001;
    Raycaster.prototype.linePrecision = 0.2;
    
    Raycaster.prototype.set = function(origin, direction) {
        
        this.ray.set(origin, direction);
          
    };
    
    Raycaster.prototype.intersectObjects = function(group, objects) {
        
        var intersects = [];
        
        for (var i = 0, l = objects.length; i < l; i++)            
            intersectObject(group, objects[i], this, intersects);
            
        intersects.sort(descSort);
        
        return intersects;
        
    };
    
    return Raycaster;
    
})();


//WebMol Projecion 
//TODO: can probably strip this down a lot (only used for selection handling)
WebMol.Projector = function () {

    var _viewMatrix = new WebMol.Matrix4(),
    _viewProjectionMatrix = new WebMol.Matrix4();

    this.projectVector = function ( vector, camera ) {

        camera.matrixWorldInverse.getInverse( camera.matrixWorld );

        _viewProjectionMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

        return vector.applyProjection( _viewProjectionMatrix );

    };

    this.unprojectVector = function ( vector, camera ) {

        camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);

        _viewProjectionMatrix.multiplyMatrices(camera.matrixWorld, camera.projectionMatrixInverse);

        return vector.applyProjection( _viewProjectionMatrix );

    };

};