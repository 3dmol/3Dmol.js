//A GLShapeCollection is a collection of user specified shapes. Includes
// build in sphere and arrow shapes, as well as generic user specified shapes



//Wrapper around Geometry allowing for quick toggle between mesh and wireframe 
//after initial vertex typed array contsruction
WebMol.Shape = function(type, color) {
        
    this.type = type || 'wireframe';
    this.color = color || {r: 0, g: 0, b: 0, a: 1.0};
    this.vertices = 0; 
    this.__vertexArray = null;
    this.__colorArray = null;
    this.__normalArray = null;
    this.__faceArray = null;
    this.__lineArray = null;
    this.initialized = false;
    //TODO: Can/should vertices be updated (e.g. scaling a sphere) ?
    
    this.geo = new WebMol.Geometry();
    this.geo.geometryChunks = [];
    this.geo.geometryChunks.push( new geometryChunk() );
    
};

WebMol.Shape.prototype = {
    
    constructor : WebMol.Shape,
    
    updateGeometry : function() {
        
        if (this.initialized) {
            
            var geoGroup = this.geo.geometryChunks[0];
            geoGroup.__vertexArray = this.__vertexArray;
            geoGroup.__colorArray = this.__colorArray;
            geoGroup.__normalArray = this.__normalArray;
            geoGroup.__faceArray = this.__faceArray;
            geoGroup.__lineArray = this.__lineArray;
            
            geoGroup.vertices = this.vertices;
            
            geoGroup.__inittedArrays = true;
            
        }

    }

};


WebMol.GLShapeCollection = (function() {

    
    // construct vertices around origin for given radius, memoize results
    var sphereVertexCache = {
        cache : {},
        getVerticesForRadius : function(radius) {

            if (typeof (this.cache[radius]) != "undefined")
                return this.cache[radius];

            var obj = {
                vertices : [],
                verticesRows : [],
                normals : []
            };
            // scale quality with radius heuristically
            var widthSegments = 12;
            var heightSegments = 10;
            if (radius < 1) {
                widthSegments = 8;
                heightSegments = 6;
            }

            var phiStart = 0;
            var phiLength = Math.PI * 2;

            var thetaStart = 0;
            var thetaLength = Math.PI;

            var x, y, vertices = [], uvs = [];

            for (y = 0; y <= heightSegments; y++) {

                var verticesRow = [];
                for (x = 0; x <= widthSegments; x++) {

                    var u = x / widthSegments;
                    var v = y / heightSegments;

                    var vertex = {};
                    vertex.x = -radius * Math.cos(phiStart + u * phiLength)
                            * Math.sin(thetaStart + v * thetaLength);
                    vertex.y = radius
                            * Math.cos(thetaStart + v * thetaLength);
                    vertex.z = radius * Math.sin(phiStart + u * phiLength)
                            * Math.sin(thetaStart + v * thetaLength);

                    var n = new WebMol.Vector3(vertex.x, vertex.y, vertex.z);
                    n.normalize();

                    obj.vertices.push(vertex);
                    obj.normals.push(n);

                    verticesRow.push(obj.vertices.length - 1);

                }

                obj.verticesRows.push(verticesRow);

            }

            this.cache[radius] = obj;
            return obj;
        }
    }; 
    
    var drawSphere = function(shape) {
        
        var center = shape.center;
        var radius = shape.radius;
        var shapeObj = shape.shape;
        var geo = shapeObj.geo;
        
        var color = shape.color;
        var C = WebMol.CC.color(color);
        
        if ((shape.clickable === true) && (shape.intersectionShape !== undefined)) {
            var intersectionSphere = new WebMol.Sphere(new WebMol.Vector(center.x, center.y, center.z), radius);
            shape.intersectionShape.sphere = intersectionSphere;
        }
        
        var vobj = sphereVertexCache.getVerticesForRadius(radius);
        
        var vertices = vobj.vertices;
        var normals = vobj.normals;
        
        var vertexArr = new Float32Array(vertices.length * 3);
        var colorArr = new Float32Array(vertices.length * 3);
        var normalArr = new Float32Array(vertices.length * 3);
        
        //Element array indices for drawing faces or lines
        var faceArr = [];
        var lineArr = [];
        
        var start = shapeObj.vertices;
        
        for (var i in vertices) {
            var v = vertices[i];
            var index = i * 3;
            
            vertexArr[index] = v.x, vertexArr[index+1] = v.y, vertexArr[index+2] = v.z;
            colorArr[index] = C.r, colorArr[index+1] = C.g, colorArr[index+2] = C.b;
        }
        
        shapeObj.vertices = vertices.length;

        var verticesRows = vobj.verticesRows;
        var h = verticesRows.length - 1;
       
        //Calculate normals
        for (var y = 0; y < h; y++) {
            var w = verticesRows[y].length - 1;
            for (var x = 0; x < w; x++) {

                var v1 = verticesRows[y][x + 1] + start, v1offset = v1 * 3;
                var v2 = verticesRows[y][x] + start, v2offset = v2 * 3;
                var v3 = verticesRows[y + 1][x] + start, v3offset = v3 * 3;
                var v4 = verticesRows[y + 1][x + 1] + start, v4offset = v4 * 3;

                var n1 = normals[v1 - start];
                var n2 = normals[v2 - start];
                var n3 = normals[v3 - start];
                var n4 = normals[v4 - start];
                var face, norm;
                if (Math.abs(vertices[v1 - start].y) === radius) {
                    //face = [v1, v3, v4];
                    //norm = [n1, n3, n4];
                    
                    normalArr[v1offset] = n1.x, normalArr[v3offset] = n3.x, normalArr[v4offset] = n4.x;
                    normalArr[v1offset+1] = n1.y, normalArr[v3offset+1] = n3.y, normalArr[v4offset+1] = n4.y;
                    normalArr[v1offset+2] = n1.z, normalArr[v3offset+2] = n3.z, normalArr[v4offset+2] = n4.z;
                    
                    faceArr.push(v1), faceArr.push(v3), faceArr.push(v4);
                    
                    lineArr.push(v1), lineArr.push(v3), lineArr.push(v1),
                    lineArr.push(v4), lineArr.push(v2), lineArr.push(v4);
                    
                } else if (Math.abs(vertices[v3 - start].y) === radius) {
                    //face = [v1, v2, v3];            
                    //norm = [n1, n2, n3];
                    
                    normalArr[v1offset] = n1.x, normalArr[v2offset] = n2.x, normalArr[v3offset] = n3.x;
                    normalArr[v1offset+1] = n1.y, normalArr[v2offset+1] = n2.y, normalArr[v3offset+1] = n3.y;
                    normalArr[v1offset+2] = n1.z, normalArr[v2offset+2] = n2.z, normalArr[v3offset+2] = n3.z;
                    
                    faceArr.push(v1), faceArr.push(v2), faceArr.push(v3);
                    
                    lineArr.push(v1), lineArr.push(v2), lineArr.push(v1),
                    lineArr.push(v3), lineArr.push(v2), lineArr.push(v3);
                    
                } else {
                    //face = [v1, v2, v3, v4];
                    //norm = [n1, n2, n3, n4];
                    
                    normalArr[v1offset] = n1.x, normalArr[v2offset] = n2.x, normalArr[v4offset] = n4.x;
                    normalArr[v1offset+1] = n1.y, normalArr[v2offset+1] = n2.y, normalArr[v4offset+1] = n4.y;
                    normalArr[v1offset+2] = n1.z, normalArr[v2offset+2] = n2.z, normalArr[v4offset+2] = n4.z;
                    
                    normalArr[v2offset] = n2.x, normalArr[v3offset] = n3.x, normalArr[v4offset] = n4.x;
                    normalArr[v2offset+1] = n2.y, normalArr[v3offset+1] = n3.y, normalArr[v4offset+1] = n4.y;
                    normalArr[v2offset+2] = n2.z, normalArr[v3offset+2] = n3.z, normalArr[v4offset+2] = n4.z;
                    
                    faceArr.push(v1), faceArr.push(v2), faceArr.push(v4);
                    faceArr.push(v2), faceArr.push(v3), faceArr.push(v4);
                    
                    lineArr.push(v1), lineArr.push(v2), lineArr.push(v1), lineArr.push(v4),
                    lineArr.push(v2), lineArr.push(v3), lineArr.push(v3), lineArr.push(v4);
                    
                }

            }
        }
        
        
        shapeObj.__vertexArray = vertexArr;
        shapeObj.__colorArray = colorArr;
        shapeObj.__normalArray = normalArr;
        shapeObj.__faceArray = new Uint16Array(faceArr);
        shapeObj.__lineArray = new Uint16Array(lineArr);
        
    };
    
    var drawArrow = function(shape) {
        
    };
    
    var drawGeneric = function(shape) {
        
    };
    
    var GLShapeCollection = function(sid) {
        
        var id = sid;
        
        //shapes, like atoms, have a style property (either sphere, arrow, or mesh),
        //and an optional 'solid' property (if true, draw as mesh, otherwise wireframe)
        var shapes = [];
        var shapeColObj = null;
        var renderedShapeColObj = null;
        
        //Construct a sphere and add to the collection
        this.addSphere = function(center, radius, type, color) {
            
            var shapeIndex = shapes.length;
            //TODO: Almost certainly will want to change this - how can clients easily access/modify spheres
            // they've added (e.g. between solid and wireframe) ?
            var sphereShape = {
                index : shapeIndex,
                style : "sphere",
                center : center || {x: 0, y: 0, z: 0},
                radius : radius || 1.5,
                type : (type === "wireframe" || type === "mesh") ? type : "wireframe",
                color : color || {r: 0, g: 0, b: 0, a: 1.0},
                shape : new WebMol.Shape(type, color)  
            };
            
            drawSphere(sphereShape);
            sphereShape.shape.initialized = true;
            
            shapes.push(sphereShape);
            
            return sphereShape;

        };
        
        this.removeShapes = function(badshapes) {
            
            shapeColObj = null;
            
            var badies = [];
            for (var i in badshapes){
                badies[badshapes.index] = true;
            }
            
            var newShapes = [];
            for (var i in shapes) {
                var shape = shapes[i];
                if (!badies[shape.index])
                    newShapes.push(shape);
            }
            
            shapes = newShapes;
        };

        //Go through shapes and regenerate geometries (i.e. reconstruct typed arrays)
        var createShapeObj = function() {
            
            console.log("creating for "+id);
            var ret = new WebMol.Object3D();
            
            var wireGeo = new WebMol.Geometry();
            wireGeo.geometryChunks = [];
            var meshGeo = new WebMol.Geometry();
            meshGeo.geometryChunks = [];
            
            for ( var i = 0; i < shapes.length; i++) {
                var shape = shapes[i];
                var shapeObj = shape.shape;
                
                
                // recreate gl info for each shape as necessary
                // set up appropriate intersection spheres for clickable shapes
                if (shape && shapeObj.initialized) {
                    
                    //TODO: How to represent arrows? Cylinder + faces ?
                    if (shape.clickable && shape.intersectionShape === undefined)
                        shape.intersectionShape = {sphere: null, cylinder: [], triangle : []};                    
                    
                    var wireframe = shape.type === 'wireframe';
                    
                    var geo = wireframe ? wireGeo : meshGeo;
                    
                    //merge shape's geometry with appropriate geo
                    shapeObj.updateGeometry();
                    geo.geometryChunks.push(shapeObj.geo.geometryChunks[0]);
                    geo.vertices += shapeObj.vertices;
    
                }
            }
            
            if (wireGeo.vertices > 0) {
                
                var wireframeMaterial = new WebMol.MeshLambertMaterial({
                    wireframe : true,
                    vertexColors : true, 
                    ambient : 0x000000,
                    reflectivity : 0 
                });
                
                var wireframeShape = new WebMol.Mesh(wireGeo, wireframeMaterial);
                ret.add(wireframeShape);
            }
 
            if (meshGeo.vertices > 0) {
                
                var meshMaterial = new WebMol.MeshLambertMaterial({
                    wireframe : false,
                    vertexColors : true, 
                    ambient : 0x000000,
                    reflectivity : 0 
                });
                
                var meshShape = new WebMol.Mesh(meshGeo, meshMaterial);
                ret.add(meshShape);
            }
            
            return ret;
               
        };   
        
        this.globj = function(group) {
            var time1 = new Date();
            if (shapeColObj === null) {
                shapeColObj = createShapeObj();
                var time2 = new Date();
                console.log("shape object creation time: " + (time2 - time1) );
                if (renderedShapeColObj) {
                    group.remove(renderedShapeColObj);
                    renderedShapeColObj = null;
                }
                renderedShapeColObj = shapeColObj.clone();
                group.add(renderedShapeColObj);
            }
        };
    
    };
    

    
    return GLShapeCollection;
    
}());
