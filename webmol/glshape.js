//A GLShape is a collection of user specified shapes. Includes
// build in sphere and arrow shapes, as well as custom user specified shapes


WebMol.GLShape = (function() {
    
    //Preset component builders
    
    //Sphere component 
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
    
    var drawSphere = function(shape, geo, spec) {
        
        var pos = spec.position, radius = spec.radius;
        
        if ((shape.clickable === true) && (shape.intersectionShape !== undefined)) {
            var center = new WebMol.Vector3(pos.x, pos.y, pos.z);
            shape.intersectionShape.sphere = new WebMol.Sphere(center, radius);
        }
                                                             
        var color = shape.color;
        var C = WebMol.CC.color(color);

        var x, y;
        var vobj = sphereVertexCache.getVerticesForRadius(radius);                
                    
        var vertices = vobj.vertices;
        var normals = vobj.normals;
        
        geoGroup = geo.updateGeoGroup(65536);
        var start = geoGroup.vertices;
        
        for (var i = 0, il = vertices.length; i < il; ++i) {
            var offset = 3*(start + i);   
            var v = vertices[i];
            
            geoGroup.__vertexArray[offset] = (v.x + pos.x);
            geoGroup.__vertexArray[offset+1] = (v.y + pos.y);
            geoGroup.__vertexArray[offset+2] = (v.z + pos.z);
            
            geoGroup.__colorArray[offset] = C.r;
            geoGroup.__colorArray[offset+1] = C.g;
            geoGroup.__colorArray[offset+2] = C.b;
           
        }
        
        geoGroup.vertices += vertices.length;
        
        var verticesRows = vobj.verticesRows;
        var h = verticesRows.length - 1;
        
        //var color = [C.r, C.g, C.b];
        for (y = 0; y < h; y++) {
            var w = verticesRows[y].length - 1;
            for (x = 0; x < w; x++) {
                
                var faceoffset = geoGroup.faceidx;
                
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
                    
                    geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v3offset] = n3.x, geoGroup.__normalArray[v4offset] = n4.x;
                    geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v3offset+1] = n3.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                    geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v3offset+2] = n3.z, geoGroup.__normalArray[v4offset+2] = n4.z;

                    geoGroup.__faceArray[faceoffset] = v1; 
                    geoGroup.__faceArray[faceoffset+1] = v3;
                    geoGroup.__faceArray[faceoffset+2] = v4;
                    
                    geoGroup.faceidx += 3;
                    
                } else if (Math.abs(vertices[v3 - start].y) === radius) {
                    //face = [v1, v2, v3];            
                    //norm = [n1, n2, n3];
                    
                    geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v3offset] = n3.x;
                    geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v3offset+1] = n3.y;
                    geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v3offset+2] = n3.z;

                    geoGroup.__faceArray[faceoffset] = v1;
                    geoGroup.__faceArray[faceoffset+1] = v2;
                    geoGroup.__faceArray[faceoffset+2] = v3;
                    
                    geoGroup.faceidx += 3;
                    
                } else {
                    //face = [v1, v2, v3, v4];
                    //norm = [n1, n2, n3, n4];
                    
                    geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v4offset] = n4.x;
                    geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                    geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v4offset+2] = n4.z;
                    
                    geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v3offset] = n3.x, geoGroup.__normalArray[v4offset] = n4.x;
                    geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v3offset+1] = n3.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                    geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v3offset+2] = n3.z, geoGroup.__normalArray[v4offset+2] = n4.z;
                    
                    geoGroup.__faceArray[faceoffset] = v1;
                    geoGroup.__faceArray[faceoffset+1] = v2;
                    geoGroup.__faceArray[faceoffset+2] = v4;
                    
                    geoGroup.__faceArray[faceoffset+3] = v2;
                    geoGroup.__faceArray[faceoffset+4] = v3;
                    geoGroup.__faceArray[faceoffset+5] = v4;
                    
                    geoGroup.faceidx += 6;
                }

            }
        }

    };
    
    var drawCylinder = function() {
    
    };
    
    var drawArrow = function(shape) {
        
    };
    

    var GLShape = function(stylespec) {
    
        this.id = WebMol.ShapeIDCount++;
        this.color = stylespec.color || new WebMol.Color();
        this.wireframe = stylespec.wireframe ? true : false;
        this.alpha = stylespec.alpha ? WebMol.Math.clamp(stylespec.alpha, 0.0, 1.0) : 1.0;
        
        var shapeObj = null;
        var renderedShapeObj = null;
        var components = [];
        
        var geo = new WebMol.Geometry(true);
        
        var createShapeComponent = function(fn) {
            
            //Force creation of new geometryGroup for each added component
            var geoGroup = geo.updateGeoGroup(65536);
            var vertexArr = [], normalArr = [], faceArr = [], lineArr = [];
            
            fn(vertexArr, normalArr, faceArr, lineArr);
            
            geoGroup.__vertexArray = new Float32Array(vertexArr);
            geoGroup.__normalArray = new Float32Array(normalArr);
            geoGroup.__colorArray = new Float32Array(vertexArr.length);

            for (var i = 0; i < geoGroup.__colorArray.length / 3; ++i) {
                geoGroup.__colorArray[i*3] = this.color.r;
                geoGroup.__colorArray[i*3 + 1] = this.color.g;
                geoGroup.__colorArray[i*3 + 2] = this.color.b;      
            }
            
            geoGroup.__faceArray = new Uint16Array(faceArr);
            geoGroup.__lineArray = new Uint16Array(lineArr);
            
        };
        
        this.addSphere = function(sphereSpec) {
          
            sphereSpec.position = sphereSpec.position || {x: 0, y: 0, z: 0};
            sphereSpec.radius = sphereSpec.radius ? WebMol.Math.clamp(sphereSpec.radius, 0, Infinity) : 1.5;
            
            drawSphere(shape, geo, sphereSpec);
                
        };        
    
        this.globj = function(group) {
            
            shapeObj = new WebMol.Object3D();
            var material = new WebMol.MeshLambertMaterial({
                wireframe : this.wireframe,
                vertexColors : true,
                ambient : 0x000000,
                reflectivity : 0    
            });
            
            var mesh = new WebMol.Mesh(geo, material);
            shapeObj.add(mesh);
                 
            if (renderedShapeObj) {
                group.remove(renderedShapeObj);
                renderedShapeObj = null;
            }
            renderedShapeObj = shapeObj.clone();
            group.add(renderedShapeObj);
            
        };
    
    };
    
       
    return GLShape;
    
}());

WebMol.ShapeIDCount = 0;