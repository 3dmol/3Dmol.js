/* core Object3D
 * Base class for Scene, Camera, Geometry
 * Geometry class
 */

// Event Handling
/** @this {$3Dmol.EventDispatcher} */
$3Dmol.EventDispatcher = function() {
  
    const listeners = {};
    
    this.addEventListener = function(type, listener) {
        if (listeners[type] === undefined)
            listeners[type] = [];
        
        if (listeners[type].indexOf(listener) === -1)
            listeners[type].push(listener);
    };  
    
    this.removeEventListener = function(type, listener) {
        
        const index = listeners[type].indexOf(listener);
        
        if (index !== -1)
            listeners[type].splice(index, 1);
              
    };
    
    this.dispatchEvent = function(event) {
        
        const listenerArray = listeners[event.type];
        
        if (listenerArray !== undefined) {
            event.target = this;
            
            for (let i = 0, l = listenerArray.length; i < l; i++)
                listenerArray[i].call(this, event);
                
        }
            
    };
    
};

$3Dmol.Color = function( color ){
    
    if ( arguments.length > 1) {
            this.r = arguments[0] || 0.0;
            this.g = arguments[1] || 0.0;
            this.b = arguments[2] || 0.0;

            return this;
    }
    
    return this.set(color);
                
};

$3Dmol.Color.prototype = {
    
    constructor: $3Dmol.Color,
    
    r: 0.0, g: 0.0, b: 0.0,
    
    set(val) {
        
            if (val instanceof $3Dmol.Color) 
                return val.clone();

            if (typeof val === 'number')
                this.setHex(val);
            
            else if (typeof val === 'object' && "r" in val && "g" in val && "b" in val) {
                this.r = val.r;
                this.g = val.g;
                this.b = val.b;
            }
    },
    
    setHex(hex) {
        
            hex = Math.floor(hex);

            this.r = (hex >> 16 & 255) / 255;
            this.g = (hex >> 8 & 255) / 255;
            this.b = (hex & 255) / 255;                                                                                     
        
            return this;
    },
    
    getHex() {
        const R = Math.round(this.r*255);
        const G = Math.round(this.g*255);
        const B = Math.round(this.b*255);
        return R<<16 | G << 8 | B;
    },
    
    clone() {
            return new $3Dmol.Color(this.r, this.g, this.b);
    },
        
    copy(color) {
        this.r = color.r;
        this.g = color.g;
        this.b = color.b;
        
        return this;
    },
    
    // return object that represents color components from 0 to 255
    scaled() {
        const ret = {};
        ret.r = Math.round(this.r*255);
        ret.g = Math.round(this.g*255);
        ret.b = Math.round(this.b*255);
        ret.a = 1.0;
        return ret;
    }
    
};

// Object3D base constructor function
/** @this {$3Dmol.Object3D} */
$3Dmol.Object3D = function() {
    
    this.id = $3Dmol.Object3DIDCount+=1;
    
    this.name = "";
    
    this.parent = undefined;
    this.children = [];
    
    this.position = new $3Dmol.Vector3();
    this.rotation = new $3Dmol.Vector3();
    this.matrix = new $3Dmol.Matrix4();
    this.matrixWorld = new $3Dmol.Matrix4();
    this.quaternion = new $3Dmol.Quaternion();
    this.eulerOrder = 'XYZ';
    
    this.up = new $3Dmol.Vector3(0, 1, 0);
    this.scale = new $3Dmol.Vector3(1, 1, 1);
    
    this.matrixAutoUpdate = true;
    this.matrixWorldNeedsUpdate = true;
    this.rotationAutoUpdate = true;
    this.useQuaternion = false;
    
    this.visible = true;
    
};

$3Dmol.Object3D.prototype = {
    
    constructor : $3Dmol.Object3D,
    
    lookAt(vector) {
        
        this.matrix.lookAt(vector, this.position, this.up);
        
        if (this.rotationAutoUpdate) {
            
            if (this.useQuaternion === true) 
                this.quaternion.copy(this.matrix.decompose()[1]);
            else
                this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
        }
    },
    
    // add child object
    add(object) {
        if (object === this){
            console.error("Can't add $3Dmol.Object3D to itself");
            return;
        }
        
        object.parent = this;
        this.children.push(object);
        
        // add to the scene (i.e. follow up this instance's parents until reach the top)
        
        let scene = this;
        
        while (scene.parent !== undefined)
            scene = scene.parent;
            
        if (scene !== undefined && scene instanceof $3Dmol.Scene) 
            scene.__addObject(object);
        
    },
    
    remove(object) {
        
        const index = this.children.indexOf(object);
        
        if (index !== -1) {
            
            object.parent = undefined;
            this.children.splice(index, 1);
            
            // Remove from scene
            
            let scene = this;
            
            while (scene.parent !== undefined)
                scene = scene.parent;
                
            if (scene !== undefined && scene instanceof $3Dmol.Scene)
                scene.__removeObject(object);
                
        }
    },
    
    // convert to vrml
    vrml(indent) {
        // attempt to pretty print
        if(!indent) indent = ' ';
        // all objects have a transformation (usually identity)
        // not quite sure if getting rotation right here..
        const theta = 2*Math.atan2(this.quaternion.lengthxyz(),this.quaternion.w);
        let x = 0; let y = 0; let z = 0;
        if(theta !== 0) {
            const st = Math.sin(theta/2);
            x = this.quaternion.x/st;
            y = this.quaternion.y/st;
            z = this.quaternion.z/st;
        } 
        let ret = `${indent}Transform {\n${ 
        indent} center ${this.position.x} ${this.position.y} ${this.position.z}\n${
        indent} rotation ${x} ${y} ${z} ${theta}\n${
        indent} children [\n`;
        
        if(this.geometry) {
            ret += this.geometry.vrml(indent, this.material);
        }
        for(let i = 0; i < this.children.length; i++) {
            ret += `${this.children[i].vrml(`${indent} `)},\n`;
        }
        ret += ' ]\n';
        ret += '}';
        return ret;
    },
    
    updateMatrix() {
        
        this.matrix.setPosition(this.position);
        
        if (this.useQuaternion === false) 
            this.matrix.setRotationFromEuler(this.rotation, this.eulerOrder);
        else
            this.matrix.setRotationFromQuaternion(this.quaternion);
        
        // TODO: Do I need this??
        if (this.scale.x !== 1 || this.scale.y !== 1 || this.scale.z !== 1)
            this.matrix.scale(this.scale);
            
        this.matrixWorldNeedsUpdate = true;
        
    },
    
    updateMatrixWorld(force) {
        
        if (this.matrixAutoUpdate === true) 
            this.updateMatrix();
        
        if (this.matrixWorldNeedsUpdate === true || force === true) {
            
            if (this.parent === undefined)
                this.matrixWorld.copy(this.matrix);
            else
                this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
                
        }
        
        this.matrixWorldNeedsUpdate = false;
        
        // Update matrices of all children
        for (let i = 0; i < this.children.length; i++) {
            this.children[i].updateMatrixWorld(true);
        }
    },
    
    clone(object) {
        
        if (object === undefined)
            object = new $3Dmol.Object3D();
            
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
        object.matrixAutoUpdate = this.matrixAutoUpdate;
        object.matrixWorldNeedsUpdate = this.matrixWorldNeedsUpdate;
        
        object.useQuaternion = this.useQuaternion;
        
        object.visible = this.visible;
        
        for (let i = 0; i < this.children.length; i++) {
            const child = this.children[i];
            object.add(child.clone());
        }
        
        return object;
        
    },
    
    setVisible(val) { // recursively set visibility
        this.visible = val;
        for (let i = 0; i < this.children.length; i++) {
            const child = this.children[i];
            child.setVisible(val);
        }
    }
};

$3Dmol.Object3DIDCount = 0;

// Geometry class
// TODO: What can I remove - how can I optimize ?
$3Dmol.Geometry = (function() {
   
    const BUFFERSIZE = 65535; // limited to 16bit indices
    
    
    /** @constructor */
    const geometryGroup = function(id) {
        this.id = id || 0;
        // for performance reasons, callers must directly modify these
        this.vertexArray = null;
        this.colorArray = null;
        this.normalArray = null;
        this.faceArray = null;
        this.radiusArray = null;
        // this.adjFaceArray=null;
        this.lineArray = null;
        this.vertices = 0;
        this.faceidx = 0;
        this.lineidx = 0;
        
    };
    
    geometryGroup.prototype.setColors = function(setcolor) {
        // apply a function that takes the vertex coordinate and returns a color
        const v = this.vertexArray;
        const c = this.colorArray;
        if(v.length !== c.length) {
            console.log("Cannot re-color geometry group due to mismatched lengths.");
            return;
        }
        for(let i = 0; i < v.length; i+= 3) {
            let col = setcolor(v[i],v[i+1],v[i+2]);
            if(!(col instanceof $3Dmol.Color)) {
                col = $3Dmol.CC.color(col);
            }
            c[i] = col.r;
            c[i+1] = col.g;
            c[i+2] = col.b;
        }
    };
    geometryGroup.prototype.getNumVertices = function() {
        return this.vertices;
    };
    
    geometryGroup.prototype.getVertices = function() {
        return this.vertexArray;
    };
    
    
    geometryGroup.prototype.getCentroid = function() {
        
        const centroid = new $3Dmol.Vector3();
        let offset; let x; let y; let z;
        
        for (let i = 0; i < this.vertices; ++i) {
            offset = i*3;
            
            x = this.vertexArray[offset]; y = this.vertexArray[offset+1]; z = this.vertexArray[offset+2];
            
            centroid.x += x; centroid.y += y; centroid.z += z;
        }
        
        // divideScalar checks for 0 denom
        centroid.divideScalar(this.vertices);
        
        return centroid;
    };
    
    // setup normals - vertex and face array must exist
    geometryGroup.prototype.setNormals = function() {        
        
        const faces = this.faceArray;
        const verts = this.vertexArray;
        const norms = this.normalArray;
        
        if (! this.vertices || ! this.faceidx) 
            return;
        
        // vertex indices
        let a; let b; let c; 
        // and actual vertices
        let vA; let vB; let vC; let norm;
            
        for (let i = 0; i < faces.length / 3; ++i) {
            
            a = faces[i * 3] * 3;
            b = faces[i * 3 + 1] * 3;
            c = faces[i * 3 + 2] * 3;
            
            vA = new $3Dmol.Vector3(verts[a], verts[a+1], verts[a+2]);
            vB = new $3Dmol.Vector3(verts[b], verts[b+1], verts[b+2]);
            vC = new $3Dmol.Vector3(verts[c], verts[c+1], verts[c+2]);
            
            vA.subVectors(vA, vB);
            vC.subVectors(vC, vB);
            vC.cross(vA);
            
            // face normal
            norm = vC;
            norm.normalize();
            
            norms[a] += norm.x; norms[b] += norm.x; norms[c] += norm.x;
            norms[a + 1] += norm.y; norms[b + 1] += norm.y; norms[c + 1] += norm.y;
            norms[a + 2] += norm.z; norms[b + 2] += norm.z; norms[c + 2] += norm.z;
            
        }             
                
    };
    
    // sets line index array from face arr
    // Note - assumes all faces are triangles (i.e. there will
    // be an extra diagonal for four-sided faces - user should 
    // specify linearr for custom shape generation to show wireframe squares
    // as rectangles rather than two triangles)
    geometryGroup.prototype.setLineIndices = function() {
        
        if (! this.faceidx)
            return;
                    
        if(this.lineArray && this.lineArray.length === this.faceidx*2 && this.lineidx === this.faceidx*2) 
            return; // assume already computed
        
        const faceArr = this.faceArray; const lineArr = this.lineArray = new Uint16Array(this.faceidx*2);      
        this.lineidx = this.faceidx*2;         
            
        for (let i = 0; i < this.faceidx / 3; ++i) {
            
            const faceoffset = i*3;
            const lineoffset = faceoffset*2;
            const a = faceArr[faceoffset]; const b = faceArr[faceoffset+1]; const c = faceArr[faceoffset+2];
            
            lineArr[lineoffset] = a; lineArr[lineoffset+1] = b;
            lineArr[lineoffset+2] = a; lineArr[lineoffset+3] = c;
            lineArr[lineoffset+4] = b; lineArr[lineoffset+5] = c;
            
        }
    };
    
    geometryGroup.prototype.vrml = function(indent, material) {
        let ret = '';
        ret += `${indent}Shape {\n${
               indent} appearance Appearance {\n${
               indent}  material Material {\n${ 
               indent}   diffuseColor ${material.color.r} ${material.color.g} ${material.color.b}\n`;
        if(material.transparent) {
            ret += `${indent}   transparency ${1.0-material.opacity}\n`;            
        }
        ret += `${indent }  }\n`; // material
        ret += `${indent} }\n`; // appearance
        
        const oldindent = indent;
        indent += ' '; // inshape
        if(material instanceof $3Dmol.LineBasicMaterial) {
            ret += `${indent}geometry IndexedLineSet {\n${ 
            indent} colorPerVertex TRUE\n${
            indent} coord Coordinate {\n${
            indent}  point [\n`;
            let x; let y; let z;
            for (let i = 0; i < this.vertices; ++i) {
                const offset = i*3;                
                x = this.vertexArray[offset]; y = this.vertexArray[offset+1]; z = this.vertexArray[offset+2];
                ret += `${indent}   ${x} ${y} ${z},\n`;
            }
            ret += `${indent}  ]\n`;
            ret += `${indent} }\n`; // end coordinate
            
            if(this.colorArray) {
                ret += `${indent} color Color {\n${
                    indent}  color [\n`;
                for (let i = 0; i < this.vertices; ++i) {
                    const offset = i*3;                
                    x = this.colorArray[offset]; y = this.colorArray[offset+1]; z = this.colorArray[offset+2];
                    ret += `${indent}   ${x} ${y} ${z},\n`;
                }
                ret += `${indent}  ]\n`;
                ret += `${indent} }\n`; // end color
            }
            
            ret += `${indent} coordIndex [\n`;
            for(let i = 0; i < this.vertices; i += 2) {
                ret += `${indent}  ${i}, ${i+1}, -1,\n`;
            }
            ret += `${indent} ]\n`;
            ret += `${indent}}\n`; // geometry
        } else {
            // faces
            ret += `${indent}geometry IndexedFaceSet {\n${ 
            indent} colorPerVertex TRUE\n${
            indent} normalPerVertex TRUE\n${
            indent} solid FALSE\n`;
            
            // vertices
            ret += `${indent} coord Coordinate {\n${
            indent}  point [\n`;
            let x; let y; let z;
            for (let i = 0; i < this.vertices; ++i) {
                const offset = i*3;                
                x = this.vertexArray[offset]; y = this.vertexArray[offset+1]; z = this.vertexArray[offset+2];
                ret += `${indent}   ${x} ${y} ${z},\n`;
            }
            ret += `${indent}  ]\n`;
            ret += `${indent} }\n`; // end coordinate
            
            // normals
            ret += `${indent} normal Normal {\n${
                   indent}  vector [\n`;
            for (let i = 0; i < this.vertices; ++i) {
                const offset = i*3;                
                x = this.normalArray[offset]; y = this.normalArray[offset+1]; z = this.normalArray[offset+2];
                ret += `${indent}   ${x} ${y} ${z},\n`;
            }
            ret += `${indent}  ]\n`;
            ret += `${indent} }\n`; // end normal
            
            // colors
            if(this.colorArray) {
                ret += `${indent} color Color {\n${
                    indent}  color [\n`;
                for (let i = 0; i < this.vertices; ++i) {
                    const offset = i*3;                
                    x = this.colorArray[offset]; y = this.colorArray[offset+1]; z = this.colorArray[offset+2];
                    ret += `${indent}   ${x} ${y} ${z},\n`;
                }
                ret += `${indent}  ]\n`;
                ret += `${indent} }\n`; // end color
            }
            
            // faces
            ret += `${indent} coordIndex [\n`;
            for(let i = 0; i < this.faceidx; i += 3) {
                x = this.faceArray[i]; y = this.faceArray[i+1]; z = this.faceArray[i+2];                
                ret += `${indent}  ${x}, ${y}, ${z}, -1,\n`;
            }
            ret += `${indent} ]\n`; // end faces 
            ret += `${indent}}\n`; // geometry            
        }
        
        ret += `${oldindent}}`; // shape
        return ret;
    };
    
    geometryGroup.prototype.truncateArrayBuffers = function(mesh, reallocatemem) {
        
        mesh = (mesh === true);
        
        const vertexArr = this.vertexArray;
            const colorArr = this.colorArray;
            const normalArr = this.normalArray;
            const faceArr = this.faceArray;
            const lineArr = this.lineArray;
            const radiusArr = this.radiusArray;

        // subarray to avoid copying and reallocating memory
        this.vertexArray = vertexArr.subarray(0,this.vertices*3);
        this.colorArray = colorArr.subarray(0,this.vertices*3);
        
        if (mesh) {
            this.normalArray = normalArr.subarray(0,this.vertices*3);
            this.faceArray = faceArr.subarray(0,this.faceidx); 
            
            if(this.lineidx > 0) // not always set so reclaim memory
                this.lineArray = lineArr.subarray(0,this.lineidx); 
            else
                this.lineArray = new Uint16Array(0);
                        
        }        
        else {
            this.normalArray = new Float32Array(0);
            this.faceArray = new Uint16Array(0);
            this.lineArray = new Uint16Array(0);
        }
        if (radiusArr) {
            this.radiusArray = radiusArr.subarray(0, this.vertices);
        }
        
        if(reallocatemem) { 
            // actually copy smaller arrays to save memory
            if(this.normalArray) this.normalArray = new Float32Array(this.normalArray);
            if(this.faceArray) this.faceArray = new Uint16Array(this.faceArray);
            if(this.lineArray) this.lineArray = new Uint16Array(this.lineArray);
            if(this.vertexArray) this.vertexArray = new Float32Array(this.vertexArray);
            if(this.colorArray) this.colorArray = new Float32Array(this.colorArray);
            if(this.radiusArray) this.radiusArray = new Float32Array(this.radiusArray);
        }
        this.__inittedArrays = true;                
    };
    
    const addGroup = function(geo) {
        const ret = new geometryGroup(geo.geometryGroups.length);
        geo.geometryGroups.push(ret);
        geo.groups = geo.geometryGroups.length;
        
        ret.vertexArray = new Float32Array(BUFFERSIZE*3);
        ret.colorArray = new Float32Array(BUFFERSIZE*3);
        
        // TODO: instantiating uint arrays according to max number of vertices
        // is dangerous, since there exists the possibility that there will be 
        // more face or line indices than vertex points - but so far that doesn't
        // seem to be the case for any of the renders 
        if (geo.mesh) {
            ret.normalArray = new Float32Array(BUFFERSIZE*3);
            ret.faceArray = new Uint16Array(BUFFERSIZE*6);
            ret.lineArray = new Uint16Array(BUFFERSIZE*6);
        }
        if (geo.radii) {
            ret.radiusArray = new Float32Array(BUFFERSIZE);
        }
        ret.useOffset = geo.offset;
        
        
        return ret;
    };
    /** @constructor */
    const Geometry = function(mesh, radii,offset) {
        
        $3Dmol.EventDispatcher.call(this);
        
        this.id = $3Dmol.GeometryIDCount+=1;
    
        this.name = '';
    
        this.hasTangents = false;
    
        this.dynamic = true; // the intermediate typed arrays will be deleted when set to false
        this.mesh = (mesh === true); // Does this geometry represent a mesh (i.e. do we need Face/Line index buffers?)
        this.radii = radii || false;
        this.offset = offset || false; // offset buffer used for instancing
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

        // Get geometry group to accomodate addVertices new vertices - create 
        // new group if necessary       
        updateGeoGroup(addVertices) {
        
            addVertices = addVertices || 0;
            
            let retGroup = this.groups > 0 ? this.geometryGroups[ this.groups - 1 ] : null;
            
            if (!retGroup || retGroup.vertices + addVertices > retGroup.vertexArray.length/3) 
                retGroup = addGroup(this);
                
            return retGroup;
            
        },
        
        // return comma separated list of IndexedFace (or Line) sets from geometry groups
        vrml(indent, material) {
            let ret = '';
            const len = this.geometryGroups.length;
            for (let g = 0; g < len; g++) {
                const geoGroup = this.geometryGroups[g];
                ret += `${geoGroup.vrml(indent, material) },\n`;                
            }
            return ret;
        },
        
        addGeoGroup() {
            return addGroup(this);  
        },
        
        setUpNormals(three) {
            
            three = three || false;
            
            for (let g = 0; g < this.groups; g++) {
            
                const geoGroup = this.geometryGroups[g];            
                
                geoGroup.setNormals(three);
                
            }  
                      
        },
        
        setColors(setcolor) {
            const len = this.geometryGroups.length;
            for (let g = 0; g < len; g++) {
                
                const geoGroup = this.geometryGroups[g];                            
                geoGroup.setColors(setcolor);
                
            }  
        },
        
        setUpWireframe() {
            for (let g = 0; g < this.groups; g++) {
                const geoGroup = this.geometryGroups[g];
                
                geoGroup.setLineIndices();
            }
        },
        
        // After vertices, colors, etc are collected in regular or typed arrays,
        //  create typed arrays from regular arrays if they don't already exist,
        initTypedArrays() {
                
            for (let g = 0; g < this.groups; g++) {
                
                const group = this.geometryGroups[g];
                
                if (group.__inittedArrays === true)
                    continue;
                
                // do not actually reallocate smaller memory here because
                // of the performance hit - if you know your geometry is small,
                // truncate manually with the second parameter true
                group.truncateArrayBuffers(this.mesh, false);
            }
            
        
        },
        
        dispose() {
            this.dispatchEvent( {type : 'dispose'} );
        }
    };

    
    return Geometry;
    
})();

Object.defineProperty($3Dmol.Geometry.prototype, "vertices", {
    
    /** @this {$3Dmol.Geometry} */
    get() {
        let vertices = 0;
        for (let g = 0; g < this.groups; g++)
            vertices += this.geometryGroups[g].vertices;
            
        return vertices;
    } 
        
});

$3Dmol.GeometryIDCount = 0;


// Raycaster
/** @constructor */
$3Dmol.Raycaster = (function() {
    
    const Raycaster = function(origin, direction, far, near) {
        
        this.ray = new $3Dmol.Ray(origin, direction);
        
        if (this.ray.direction.lengthSq() > 0) 
            this.ray.direction.normalize();
        
        this.near = near || 0;
        this.far = far || Infinity;
    
    };
    
    const sphere = new $3Dmol.Sphere();
    const cylinder = new $3Dmol.Cylinder();
    const triangle = new $3Dmol.Triangle();
    const w0 = new $3Dmol.Vector3(); // for cylinders, cylinder.c1 - ray.origin
    const v1 = new $3Dmol.Vector3(); // all purpose local vector
    const v2 = new $3Dmol.Vector3();
    const v3 = new $3Dmol.Vector3();
    // var facePlane = new $3Dmol.Plane();
    const matrixPosition = new $3Dmol.Vector3();
            
    const descSort = function(a, b) {
        return a.distance - b.distance;
    };

    // [-1, 1]
    const clamp = function(x) {
        return Math.min(Math.max(x, -1), 1);
    };
    
    // object is a Sphere or (Bounding) Box
    const intersectObject = function(group, clickable, raycaster, intersects) {
        
        matrixPosition.getPositionFromMatrix(group.matrixWorld);
        
        if (clickable.intersectionShape === undefined)
            return intersects;       
        const {intersectionShape} = clickable;
        let precision = raycaster.linePrecision;
        precision *= group.matrixWorld.getMaxScaleOnAxis();
        const precisionSq = precision*precision;

        // Check for intersection with clickable's bounding sphere, if it exists
        if (clickable.boundingSphere !== undefined && clickable.boundingSphere instanceof $3Dmol.Sphere) {
            sphere.copy(clickable.boundingSphere);
            sphere.applyMatrix4(group.matrixWorld);          
            if (!raycaster.ray.isIntersectionSphere(sphere)) {             
              return intersects;
            }
        }     


        // Iterate through intersection objects
        let i; let il;
            let norm; let normProj; let cylProj; let rayProj;
            let distance; let closestDistSq; let denom; let discriminant;
            let s; let t; let sC; let tC;
        // triangle faces
        for (i = 0, il = intersectionShape.triangle.length; i < il; i++) {
            
            if (intersectionShape.triangle[i] instanceof $3Dmol.Triangle) {
                
                triangle.copy(intersectionShape.triangle[i]);
                triangle.applyMatrix4(group.matrixWorld);
                
                norm = triangle.getNormal();
                
                normProj = raycaster.ray.direction.dot(norm);
                
                // face culling
                if (normProj >= 0)
                    continue;
                
                w0.subVectors(triangle.a, raycaster.ray.origin);
                
                distance = (norm.dot(w0)) / normProj;
                
                if (distance < 0)
                    continue;
                    
                // intersects with plane, check if P inside triangle
                v1.copy(raycaster.ray.direction).multiplyScalar(distance).add(raycaster.ray.origin);
                v1.sub(triangle.a); // from pt a to intersection point P
                
                v2.copy(triangle.b).sub(triangle.a); // from pt a to b
                v3.copy(triangle.c).sub(triangle.a); // from pt a to c
                const bDotC = v2.dot(v3);
                const bSq = v2.lengthSq();
                const cSq = v3.lengthSq();
                
                // P = A + s(v2) + t(v3), inside trianle if 0 <= s, t <=1  and (s + t) <=0
                
                t = ( bSq*v1.dot(v3) - bDotC*v1.dot(v2) ) / ( bSq*cSq - bDotC*bDotC );
                
                if (t < 0 || t > 1)
                    continue;
                
                s = ( v1.dot(v2) - t*bDotC ) / bSq;
                
                if ( (s < 0 || s > 1) || s + t > 1)
                    continue;
                    
                else{
                    intersects.push({clickable,
                                     distance});
				}  
            }
        }    
        // cylinders
        for (i = 0, il = intersectionShape.cylinder.length; i < il; i++) {
            
            if (intersectionShape.cylinder[i] instanceof $3Dmol.Cylinder){
                
                cylinder.copy(intersectionShape.cylinder[i]);
                cylinder.applyMatrix4(group.matrixWorld);
                
                w0.subVectors(cylinder.c1, raycaster.ray.origin); 
                
                cylProj = w0.dot(cylinder.direction); // Dela
                rayProj = w0.dot(raycaster.ray.direction); // Epsilon
                
                normProj = clamp(raycaster.ray.direction.dot(cylinder.direction)); // Beta
                
                denom = 1 - normProj*normProj;
                
                if (denom === 0.0)
                    continue;
                
                sC = (normProj*rayProj - cylProj) / denom;
                tC = (rayProj - normProj*cylProj) / denom;
                
                v1.copy(cylinder.direction).multiplyScalar(sC).add(cylinder.c1);  // Q_c
                v2.copy(raycaster.ray.direction).multiplyScalar(tC).add(raycaster.ray.origin); // P_c
                
                closestDistSq = v3.subVectors(v1, v2).lengthSq();
                const radiusSq = cylinder.radius*cylinder.radius;
                
                // Smoothing?
                // if (closestDistSq > radiusSq) radiusSq += precisionSq;
                
                // closest distance between ray and cylinder axis not greater than cylinder radius;
                // might intersect this cylinder between atom and bond midpoint
                if (closestDistSq <= radiusSq){

                    // Find points where ray intersects sides of cylinder
                    discriminant = (normProj*cylProj - rayProj)*(normProj*cylProj - rayProj) - 
                            denom*(w0.lengthSq() - cylProj*cylProj - radiusSq);
                    
                    // ray tangent to cylinder?
                    if (discriminant <= 0)
                        t = distance = Math.sqrt(closestDistSq);
                    else
                        t = distance = ( (rayProj - normProj*cylProj) - Math.sqrt(discriminant) ) / denom; 
                    
                    // find closest intersection point; make sure it's between atom's position and cylinder midpoint
                    
                    s = normProj*t - cylProj;
                    
                    // does not intersect cylinder between atom and midpoint,
                    // or intersects cylinder behind camera
                    if (s < 0 || s*s > cylinder.lengthSq() || t < 0)
                        continue;
                    
                    else
                        intersects.push({clickable,
                                         distance});
                    
                }
                    
                
            }
            
        }       
        // lines
        for (i = 0, il = intersectionShape.line.length; i < il; i += 2) {
            
            v1.copy(intersectionShape.line[i]);
            v1.applyMatrix4(group.matrixWorld);
            v2.copy(intersectionShape.line[i+1]);
            v2.applyMatrix4(group.matrixWorld);
            
            v3.subVectors(v2, v1);
            const bondLengthSq = v3.lengthSq();
            v3.normalize();
            
            w0.subVectors(v1, raycaster.ray.origin);
            
            const lineProj = w0.dot(v3);
            rayProj = w0.dot(raycaster.ray.direction);
            
            normProj = clamp(raycaster.ray.direction.dot(v3));
            
            denom = 1 - normProj*normProj;
            
            if (denom === 0.0)
                continue;
            
            sC = (normProj*rayProj - lineProj) / denom;
            tC = (rayProj - normProj*lineProj) / denom;
            
            v1.add(v3.multiplyScalar(sC)); // Q_c
            v2.copy(raycaster.ray.direction).multiplyScalar(tC).add(raycaster.ray.origin); // P_c
            
            closestDistSq = v3.subVectors(v2, v1).lengthSq();
            
            if (closestDistSq < precisionSq && sC*sC < bondLengthSq)
                intersects.push({clickable,
                                 distance : tC
                                });
            
        }
        for (i = 0, il = intersectionShape.sphere.length; i < il; i++) {
            // sphere
            if (intersectionShape.sphere[i] instanceof $3Dmol.Sphere) {
                
                sphere.copy(intersectionShape.sphere[i]);
                sphere.applyMatrix4(group.matrixWorld);
                
                if (raycaster.ray.isIntersectionSphere(sphere)) {
                    
                    v1.subVectors(sphere.center, raycaster.ray.origin);
                    
                    // distance from ray origin to point on the ray normal to sphere's center
                    // must be less than sphere's radius (since ray intersects sphere)
                    const distanceToCenter = v1.dot(raycaster.ray.direction);
                    
                    discriminant = distanceToCenter*distanceToCenter - (v1.lengthSq() - sphere.radius*sphere.radius);
                    
                    // Don't select if sphere center behind camera
                    if (distanceToCenter < 0) 
                        return intersects;
                    
                    // ray tangent to sphere?
                    if (discriminant <= 0)
                        distance = distanceToCenter;
                    
                    // This is reversed if sphere is closer than ray origin.  Do we have 
                    // to worry about handling that case?
                    else 
                        distance = distanceToCenter - Math.sqrt(discriminant);
    
                    intersects.push({clickable, 
                                     distance});                    
                }
            }        
       }
       return intersects; 
    };   
       
    Raycaster.prototype.precision = 0.0001;
    Raycaster.prototype.linePrecision = 0.2;
    
    Raycaster.prototype.set = function(origin, direction) {
        
        this.ray.set(origin, direction);
          
    };
    
    Raycaster.prototype.setFromCamera =  function ( ) { 
        const _viewProjectionMatrix = new $3Dmol.Matrix4();
        return function(coords, camera ) {    

            if ( !camera.ortho ) {            
                this.ray.origin.setFromMatrixPosition( camera.matrixWorld );
                this.ray.direction.set( coords.x, coords.y, coords.z);

                camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);
                _viewProjectionMatrix.multiplyMatrices(camera.matrixWorld, camera.projectionMatrixInverse);
                this.ray.direction.applyProjection( _viewProjectionMatrix );                                
                this.ray.direction.sub( this.ray.origin ).normalize();
    
            } else {
                this.ray.origin.set( coords.x, coords.y, ( camera.near + camera.far ) / ( camera.near - camera.far ) ).unproject( camera ); 
                this.ray.direction.set( 0, 0, - 1 ).transformDirection( camera.matrixWorld );
    
            } 
        };
    }();
    
    Raycaster.prototype.intersectObjects = function(group, objects) {     
        const intersects = [];
        
        for (let i = 0, l = objects.length; i < l; i++)            
            intersectObject(group, objects[i], this, intersects);
            
        intersects.sort(descSort);
        
        return intersects;
        
    };
    
    return Raycaster;
    
})();


// $3Dmol Projection 
/** @constructor */
$3Dmol.Projector = function () {

    const _viewProjectionMatrix = new $3Dmol.Matrix4();

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
