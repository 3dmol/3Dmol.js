//A GLShape is a collection of user specified shapes. Includes
// build in sphere and arrow shapes, as well as custom user specified shapes


WebMol.GLShape = (function() {

    var updateColor = function(geo, color) {
        
        var C = color || WebMol.CC.color(color);
        geo.colorsNeedUpdate = true;
        
        for (var g in geo.geometryGroups) {
            
            var geoGroup = geo.geometryGroups[g];
            var colorArr = geoGroup.__colorArray;
            
            for (var i = 0, il = geoGroup.vertices; i < il; ++i) {
                colorArr[i*3] = C.r, colorArr[i*3+1] = C.g, colorArr[i*3+2] = C.b;                       
            }
        }
        
    };
    
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
    
    var drawSphere = function(shape, geoGroup, spec) {
        
        var pos = spec.position, radius = spec.radius;        
        
        var center = new WebMol.Vector3(pos.x, pos.y, pos.z);
        shape.intersectionShape.sphere.push( new WebMol.Sphere(center, radius) );                                                                  

        var x, y;
        var vobj = sphereVertexCache.getVerticesForRadius(radius);                
                    
        var vertices = vobj.vertices;
        var normals = vobj.normals;
        
        var start = geoGroup.vertices;
        
        for (var i = 0, il = vertices.length; i < il; ++i) {
            var offset = 3*(start + i);   
            var v = vertices[i];
            
            geoGroup.__vertexArray[offset] = (v.x + pos.x);
            geoGroup.__vertexArray[offset+1] = (v.y + pos.y);
            geoGroup.__vertexArray[offset+2] = (v.z + pos.z);            
           
        }
        
        geoGroup.vertices += vertices.length;
        
        var verticesRows = vobj.verticesRows;
        var h = verticesRows.length - 1;
        
        //var color = [C.r, C.g, C.b];
        for (y = 0; y < h; y++) {
            var w = verticesRows[y].length - 1;
            for (x = 0; x < w; x++) {
                
                var faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;
                
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
                    
                    geoGroup.__lineArray[lineoffset] = v1, geoGroup.__lineArray[lineoffset+1] = v3, geoGroup.__lineArray[lineoffset+2] = v1;
                    geoGroup.__lineArray[lineoffset+3] = v4, geoGroup.__lineArray[lineoffset+4] = v3, geoGroup.__lineArray[lineoffset+5] = v4;
                    
                    geoGroup.faceidx += 3;
                    geoGroup.lineidx += 6;
                    
                } else if (Math.abs(vertices[v3 - start].y) === radius) {
                    //face = [v1, v2, v3];            
                    //norm = [n1, n2, n3];
                    
                    geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v3offset] = n3.x;
                    geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v3offset+1] = n3.y;
                    geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v3offset+2] = n3.z;

                    geoGroup.__faceArray[faceoffset] = v1;
                    geoGroup.__faceArray[faceoffset+1] = v2;
                    geoGroup.__faceArray[faceoffset+2] = v3;
                    
                    geoGroup.__lineArray[lineoffset] = v1, geoGroup.__lineArray[lineoffset+1] = v2, geoGroup.__lineArray[lineoffset+2] = v1;
                    geoGroup.__lineArray[lineoffset+3] = v3, geoGroup.__lineArray[lineoffset+4] = v2, geoGroup.__lineArray[lineoffset+5] = v3;
                    
                    geoGroup.faceidx += 3;
                    geoGroup.lineidx += 6;
                    
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
                    
                    geoGroup.__lineArray[lineoffset] = v1, geoGroup.__lineArray[lineoffset+1] = v2,
                    geoGroup.__lineArray[lineoffset+2] = v1, geoGroup.__lineArray[lineoffset+3] = v4;
                    
                    geoGroup.__lineArray[lineoffset+4] = v2, geoGroup.__lineArray[lineoffset+5] = v3,
                    geoGroup.__lineArray[lineoffset+6] = v3, geoGroup.__lineArray[lineoffset+7] = v4;
                    
                    geoGroup.faceidx += 6;
                    geoGroup.lineidx += 8;
                    
                }

            }
        }

    };
    
    var drawCylinder = function(shape, geoGroup, spec) {
        
        var from = spec.from, to = spec.to, radius = spec.radius;

        if (!from || !from)
            return;

        // vertices
        
        var dir = to.clone();
        dir.sub(from);

        // get orthonormal vector
        var nvecs = [];
        nvecs[0] = dir.clone();
        if (Math.abs(nvecs[0].x) > .0001)
            nvecs[0].y += 1;
        else
            nvecs[0].x += 1;
        nvecs[0].cross(dir);
        nvecs[0].normalize();

        nvecs[0] = nvecs[0];
        // another orth vector
        nvecs[4] = nvecs[0].clone();
        nvecs[4].crossVectors(nvecs[0], dir);
        nvecs[4].normalize();
        nvecs[8] = nvecs[0].clone().negate();
        nvecs[12] = nvecs[4].clone().negate();

        // now quarter positions
        nvecs[2] = nvecs[0].clone().add(nvecs[4]).normalize();
        nvecs[6] = nvecs[4].clone().add(nvecs[8]).normalize();
        nvecs[10] = nvecs[8].clone().add(nvecs[12]).normalize();
        nvecs[14] = nvecs[12].clone().add(nvecs[0]).normalize();

        // eights
        nvecs[1] = nvecs[0].clone().add(nvecs[2]).normalize();
        nvecs[3] = nvecs[2].clone().add(nvecs[4]).normalize();
        nvecs[5] = nvecs[4].clone().add(nvecs[6]).normalize();
        nvecs[7] = nvecs[6].clone().add(nvecs[8]).normalize();
        nvecs[9] = nvecs[8].clone().add(nvecs[10]).normalize();
        nvecs[11] = nvecs[10].clone().add(nvecs[12]).normalize();
        nvecs[13] = nvecs[12].clone().add(nvecs[14]).normalize();
        nvecs[15] = nvecs[14].clone().add(nvecs[0]).normalize();

        //var start = geo.vertices.length;
        var start = geoGroup.vertices;
        // add vertices, opposing vertices paired together
        for ( var i = 0, n = nvecs.length; i < n; ++i) {
            var offset = 3*(start + 2*i);
            var bottom = nvecs[i].clone().multiplyScalar(radius).add(from);
            var top = nvecs[i].clone().multiplyScalar(radius).add(to);

            geoGroup.__vertexArray[offset] = bottom.x;
            geoGroup.__vertexArray[offset+1] = bottom.y;
            geoGroup.__vertexArray[offset+2] = bottom.z;             
            
            geoGroup.__vertexArray[offset+3] = top.x;
            geoGroup.__vertexArray[offset+4] = top.y;
            geoGroup.__vertexArray[offset+5] = top.z;                    
            
        }
        
        geoGroup.vertices += 32;
        
        // now faces
        var face, norm, offset, faceoffset, lineoffset;
        var n_vertices = 0;
        for ( var i = 0, n = nvecs.length - 1; i < n; ++i) {
        
            var ti = start + 2 * i, offset = ti * 3;
            faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;
            
            var t1 = ti, t1offset = t1 * 3;
            var t2 = ti + 1, t2offset = t2 * 3;
            var t3 = ti + 3, t3offset = t3 * 3;
            var t4 = ti + 2, t4offset = t4 * 3;
            
            //face = [t1, t2, t4], [t2, t3, t4];    
            //face = [t1, t2, t3, t4];
                
            norm = [ nvecs[i], nvecs[i], nvecs[i + 1], nvecs[i + 1]];
            var n1, n2, n3, n4;
            n1 = n2 = nvecs[i];
            n3 = n4 = nvecs[i + 1];
            
            geoGroup.__normalArray[t1offset] = n1.x, geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t4offset] = n4.x;
            geoGroup.__normalArray[t1offset+1] = n1.y, geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t4offset+1] = n4.y;
            geoGroup.__normalArray[t1offset+2] = n1.z, geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t4offset+2] = n4.z;
            
            geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t3offset] = n3.x, geoGroup.__normalArray[t4offset] = n4.x;
            geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t3offset+1] = n3.y, geoGroup.__normalArray[t4offset+1] = n4.y;
            geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t3offset+2] = n3.z, geoGroup.__normalArray[t4offset+2] = n4.z;
            
            geoGroup.__faceArray[faceoffset] = t1, geoGroup.__faceArray[faceoffset+1] = t2, geoGroup.__faceArray[faceoffset+2] = t4;
            geoGroup.__faceArray[faceoffset+3] = t2, geoGroup.__faceArray[faceoffset+4] = t3, geoGroup.__faceArray[faceoffset+5] = t4;
            
            geoGroup.__lineArray[lineoffset] = t1, geoGroup.__lineArray[lineoffset+1] = t2;
            geoGroup.__lineArray[lineoffset+2] = t1, geoGroup.__lineArray[lineoffset+3] = t4;
            
            geoGroup.__lineArray[lineoffset+4] = t2, geoGroup.__lineArray[lineoffset+5] = t3;
            geoGroup.__lineArray[lineoffset+6] = t3, geoGroup.__lineArray[lineoffset+7] = t4;
            
            geoGroup.faceidx += 6;
            geoGroup.lineidx += 8;
            
        }
        // final face

        face = [start + 30, start + 31, start + 1, start];
        norm = [ nvecs[15], nvecs[15], nvecs[0], nvecs[0] ];
        
        faceoffset = geoGroup.faceidx, lineidx = geoGroup.lineidx;
        
        var t1 = face[0], t1offset = t1 * 3;
        var t2 = face[1], t2offset = t2 * 3;
        var t3 = face[2], t3offset = t3 * 3;
        var t4 = face[3], t4offset = t4 * 3;
        var n1, n2, n3, n4;
        
        n1 = n2 = nvecs[15];
        n3 = n4 = nvecs[0];

        geoGroup.__normalArray[t1offset] = n1.x, geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t4offset] = n4.x;
        geoGroup.__normalArray[t1offset+1] = n1.y, geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t4offset+1] = n4.y;
        geoGroup.__normalArray[t1offset+2] = n1.z, geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t4offset+2] = n4.z;
        
        geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t3offset] = n3.x, geoGroup.__normalArray[t4offset] = n4.x;
        geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t3offset+1] = n3.y, geoGroup.__normalArray[t4offset+1] = n4.y;
        geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t3offset+2] = n3.z, geoGroup.__normalArray[t4offset+2] = n4.z;

        geoGroup.__faceArray[faceoffset] = t1, geoGroup.__faceArray[faceoffset+1] = t2, geoGroup.__faceArray[faceoffset+2] = t4;
        geoGroup.__faceArray[faceoffset+3] = t2, geoGroup.__faceArray[faceoffset+4] = t3, geoGroup.__faceArray[faceoffset+5] = t4;
        
        geoGroup.__lineArray[lineoffset] = t1, geoGroup.__lineArray[lineoffset+1] = t2;
        geoGroup.__lineArray[lineoffset+2] = t1, geoGroup.__lineArray[lineoffset+3] = t4;
        
        geoGroup.__lineArray[lineoffset+4] = t2, geoGroup.__lineArray[lineoffset+5] = t3;
        geoGroup.__lineArray[lineoffset+6] = t3, geoGroup.__lineArray[lineoffset+7] = t4;
        
        geoGroup.faceidx += 6;        
        geoGroup.lineidx += 8;
        
    };
    
    var drawArrow = function(shape) {
        
    };
    
    //Update a bounding sphere's position and radius
    //from list of centroids and new points
    var updateBoundingFromPoints = function(sphere, components, points) {       
           
        sphere.center.set(0,0,0);
            
        if (components.length > 0) {
                            
            for (var i = 0, il = components.length; i < il; ++i) {
                var centroid = components[i].centroid;
                sphere.center.add(centroid);            
            }                
            
            sphere.center.divideScalar(components.length);           
        }
       
        var maxRadiusSq = sphere.radius*sphere.radius;
        
        for (var i = 0, il = points.length / 3; i < il; i++) {              
            var x = points[i*3], y = points[i*3 + 1], z = points[i*3 + 2];
            var radiusSq = sphere.center.distanceToSquared({x:x, y:y, z:z});                
            maxRadiusSq = Math.max(maxRadiusSq, radiusSq);
        }
        
        sphere.radius = Math.sqrt(maxRadiusSq);                        

    };
    

    var GLShape = function(stylespec) {
    
        this.id = WebMol.ShapeIDCount++;
        this.color = stylespec.color || new WebMol.Color();
        this.wireframe = stylespec.wireframe ? true : false;
        this.alpha = stylespec.alpha ? WebMol.Math.clamp(stylespec.alpha, 0.0, 1.0) : 1.0;
        
        this.boundingSphere = new WebMol.Sphere();
        
        //Click handling
        this.clickable = stylespec.clickable ? true : false;
        this.callback = typeof(stylespec.callback) == "function" ? stylespec.callback : null;
        this.intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
        
        //Keep track of shape components and their centroids
        var components = [];
        var shapeObj = null;
        var renderedShapeObj = null;
        
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
        
        //TODO: Refactor so 'drawSphere' method automatically updates bounding sphere as vertices are added
        this.addSphere = function(sphereSpec) {
          
            sphereSpec.position = sphereSpec.position || {x: 0, y: 0, z: 0};
            sphereSpec.radius = sphereSpec.radius ? WebMol.Math.clamp(sphereSpec.radius, 0, Infinity) : 1.5;
            
            var geoGroup = geo.addGeoGroup();
            drawSphere(shape, geoGroup, sphereSpec);
            geo.initTypedArrays();
            
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup, //has to be last group added
                centroid : new WebMol.Vector3(sphereSpec.position.x, sphereSpec.position.y, sphereSpec.position.z)
            });
            
            updateBoundingFromPoints( this.boundingSphere, components, geo.updateGeoGroup().__vertexArray );
        };        
        
        this.addCylinder = function(cylSpec) {
            
            cylSpec.from = cylSpec.from || {};
            cylSpec.to = cylSpec.to || {};
            
            cylSpec.from = new WebMol.Vector3(cylSpec.from.x || 0, cylSpec.from.y || 0, cylSpec.from.z || 0);
            cylSpec.to = new WebMol.Vector3(cylSpec.to.x || 1, cylSpec.to.y || 0, cylSpec.to.z || 0);
            cylSpec.radius = cylSpec.radius || 1.0;
            
            var geoGroup = geo.addGeoGroup();
            
            drawCylinder(shape, geoGroup, cylSpec);
                
        };
    
        //TODO: Adding multiple overlapping shapes in wireframe mode should obscure overlapped meshes
        this.globj = function(group) {
            
            geo.initTypedArrays();
            
            updateColor(geo, this.color);
            
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

    Object.defineProperty(GLShape.prototype, "position", {
        
        get : function() {
            return this.boundingSphere.center;
        } 
            
    });    
    
    Object.defineProperty(GLShape.prototype, "x", {
        
        get : function() {
            return this.boundingSphere.center.x;
        } 
            
    });  
 
    Object.defineProperty(GLShape.prototype, "y", {
        
        get : function() {
            return this.boundingSphere.center.y;
        } 
            
    });   

    Object.defineProperty(GLShape.prototype, "z", {
        
        get : function() {
            return this.boundingSphere.center.z;
        } 
            
    });  
  
    return GLShape;
    
}());

WebMol.ShapeIDCount = 0;