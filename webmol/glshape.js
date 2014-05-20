//A GLShape is a collection of user specified shapes. Includes
// build in sphere and arrow shapes, as well as custom user specified shapes


WebMol.GLShape = (function() {
    
    //Marching cube, to match with protein surface generation
    var ISDONE = 2;

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

            if (typeof (this.cache[radius]) !== "undefined")
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
        
        var pos = spec.center, radius = spec.radius;        
        
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
    
    var drawArrow = function(shape, geoGroup, spec) {
        
        var from = spec.start, end = spec.end, radius = spec.radius, radiusRatio = spec.radiusRatio, mid = spec.mid;

        if (!from || !end)
            return;
        
        // vertices
        
        var dir = end.clone();
        dir.sub(from).multiplyScalar(mid);
        var to = from.clone().add(dir);
        var negDir = dir.clone().negate();
        
        shape.intersectionShape.cylinder.push( new WebMol.Cylinder(from.clone(), to.clone(), radius) );
        shape.intersectionShape.sphere.push( new WebMol.Sphere(from.clone(), radius) );
        
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
            var offset = 3*(start + 3*i);
            var bottom = nvecs[i].clone().multiplyScalar(radius).add(from);
            var top = nvecs[i].clone().multiplyScalar(radius).add(to);
            var conebase = nvecs[i].clone().multiplyScalar(radius*radiusRatio).add(to);

            geoGroup.__vertexArray[offset] = bottom.x;
            geoGroup.__vertexArray[offset+1] = bottom.y;
            geoGroup.__vertexArray[offset+2] = bottom.z;             
            
            geoGroup.__vertexArray[offset+3] = top.x;
            geoGroup.__vertexArray[offset+4] = top.y;
            geoGroup.__vertexArray[offset+5] = top.z; 
            
            geoGroup.__vertexArray[offset+6] = conebase.x;
            geoGroup.__vertexArray[offset+7] = conebase.y;
            geoGroup.__vertexArray[offset+8] = conebase.z;
            
            if (i > 0) {
                var prev_x = geoGroup.__vertexArray[offset-3];
                var prev_y = geoGroup.__vertexArray[offset-2];
                var prev_z = geoGroup.__vertexArray[offset-1];
                
                var c = new WebMol.Vector3(prev_x, prev_y, prev_z);
                var b = end.clone(), b2 = to.clone();
                var a = new WebMol.Vector3(conebase.x, conebase.y, conebase.z);
                
                shape.intersectionShape.triangle.push( new WebMol.Triangle(a, b, c) );
                shape.intersectionShape.triangle.push( new WebMol.Triangle(c.clone(), b2, a.clone()) );
            }
        }
        
        geoGroup.vertices += 48;
        offset = geoGroup.vertices*3;
        
        //caps
        geoGroup.__vertexArray[offset] = from.x;
        geoGroup.__vertexArray[offset+1] = from.y;
        geoGroup.__vertexArray[offset+2] = from.z;
        
        geoGroup.__vertexArray[offset+3] = to.x;
        geoGroup.__vertexArray[offset+4] = to.y;
        geoGroup.__vertexArray[offset+5] = to.z;
        
        geoGroup.__vertexArray[offset+6] = end.x;
        geoGroup.__vertexArray[offset+7] = end.y;
        geoGroup.__vertexArray[offset+8] = end.z;
        
        geoGroup.vertices += 3;
        
        // now faces
        var face, norm, offset, faceoffset, lineoffset;
        var n_vertices = 0;
        var fromi = geoGroup.vertices - 3, toi = geoGroup.vertices - 2, endi = geoGroup.vertices - 1;
        var fromoffset = fromi*3, tooffset = toi*3, endoffset = endi*3;
        for (var i = 0, n = nvecs.length - 1; i < n; ++i) {
        
            var ti = start + 3 * i, offset = ti * 3;
            faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;
            
            var t1 = ti, t1offset = t1 * 3;
            var t2 = ti + 1, t2offset = t2 * 3;
            var t2b = ti + 2, t2boffset = t2b * 3;
            var t3 = ti + 4, t3offset = t3 * 3;
            var t3b = ti + 5, t3boffset = t3b * 3;
            var t4 = ti + 3, t4offset = t4 * 3;
            
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
            
            geoGroup.__normalArray[t2boffset] = n2.x, geoGroup.__normalArray[t3boffset] = n3.x;
            geoGroup.__normalArray[t2boffset+1] = n2.y, geoGroup.__normalArray[t3boffset+1] = n3.y;
            geoGroup.__normalArray[t2boffset+2] = n2.z, geoGroup.__normalArray[t3boffset+2] = n3.z;
            
            //sides
            geoGroup.__faceArray[faceoffset] = t1, geoGroup.__faceArray[faceoffset+1] = t2, geoGroup.__faceArray[faceoffset+2] = t4;
            geoGroup.__faceArray[faceoffset+3] = t2, geoGroup.__faceArray[faceoffset+4] = t3, geoGroup.__faceArray[faceoffset+5] = t4;
            //caps
            geoGroup.__faceArray[faceoffset+6] = t1, geoGroup.__faceArray[faceoffset+7] = t4, geoGroup.__faceArray[faceoffset+8] = fromi;
            geoGroup.__faceArray[faceoffset+9] = t2b, geoGroup.__faceArray[faceoffset+10] = toi, geoGroup.__faceArray[faceoffset+11] = t3b;
            //arrowhead
            geoGroup.__faceArray[faceoffset+12] = t2b, geoGroup.__faceArray[faceoffset+13] = endi, geoGroup.__faceArray[faceoffset+14] = t3b;
            
            //sides
            geoGroup.__lineArray[lineoffset] = t1, geoGroup.__lineArray[lineoffset+1] = t2;
            geoGroup.__lineArray[lineoffset+2] = t1, geoGroup.__lineArray[lineoffset+3] = t4;           
            //geoGroup.__lineArray[lineoffset+4] = t2, geoGroup.__lineArray[lineoffset+5] = t3;
            geoGroup.__lineArray[lineoffset+4] = t3, geoGroup.__lineArray[lineoffset+5] = t4;
            //caps
            geoGroup.__lineArray[lineoffset+6] = t1, geoGroup.__lineArray[lineoffset+7] = t4;
            //geoGroup.__lineArray[lineoffset+10] = t1, geoGroup.__lineArray[lineoffset+11] = fromi;           
            //geoGroup.__lineArray[lineoffset+12] = t4, geoGroup.__lineArray[lineoffset+13] = fromi;
            
            geoGroup.__lineArray[lineoffset+8] = t2b, geoGroup.__lineArray[lineoffset+9] = t2; //toi   
            geoGroup.__lineArray[lineoffset+10] = t2b, geoGroup.__lineArray[lineoffset+11] = t3b;
            geoGroup.__lineArray[lineoffset+12] = t3, geoGroup.__lineArray[lineoffset+13] = t3b; //toi
            //arrowhead
            geoGroup.__lineArray[lineoffset+14] = t2b, geoGroup.__lineArray[lineoffset+15] = endi;
            geoGroup.__lineArray[lineoffset+16] = t2b, geoGroup.__lineArray[lineoffset+17] = t3b;
            geoGroup.__lineArray[lineoffset+18] = endi, geoGroup.__lineArray[lineoffset+19] = t3b;
                             
            geoGroup.faceidx += 15;
            geoGroup.lineidx += 20;
            
        }
        // final face

        face = [start + 45, start + 46, start + 1, start, start + 47, start + 2];
        norm = [ nvecs[15], nvecs[15], nvecs[0], nvecs[0] ];
        
        faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;
        
        var t1 = face[0], t1offset = t1 * 3;
        var t2 = face[1], t2offset = t2 * 3;
        var t2b = face[4], t2boffset = t2b * 3;
        var t3 = face[2], t3offset = t3 * 3;
        var t3b = face[5], t3boffset = t3b * 3;
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
        
        geoGroup.__normalArray[t2boffset] = n2.x, geoGroup.__normalArray[t3boffset] = n3.x;
        geoGroup.__normalArray[t2boffset+1] = n2.y, geoGroup.__normalArray[t3boffset+1] = n3.y;
        geoGroup.__normalArray[t2boffset+2] = n2.z, geoGroup.__normalArray[t3boffset+2] = n3.z;
        
        //Cap normals
        dir.normalize(), negDir.normalize();
        geoGroup.__normalArray[fromoffset] = negDir.x, geoGroup.__normalArray[tooffset] = geoGroup.__normalArray[endoffset] = dir.x;
        geoGroup.__normalArray[fromoffset+1] = negDir.y, geoGroup.__normalArray[tooffset+1] = geoGroup.__normalArray[endoffset+1] = dir.y;
        geoGroup.__normalArray[fromoffset+2] = negDir.z, geoGroup.__normalArray[tooffset+2] = geoGroup.__normalArray[endoffset+2] = dir.z;
        
        //Final side
        geoGroup.__faceArray[faceoffset] = t1, geoGroup.__faceArray[faceoffset+1] = t2, geoGroup.__faceArray[faceoffset+2] = t4;
        geoGroup.__faceArray[faceoffset+3] = t2, geoGroup.__faceArray[faceoffset+4] = t3, geoGroup.__faceArray[faceoffset+5] = t4;
        //final caps
        geoGroup.__faceArray[faceoffset+6] = t1, geoGroup.__faceArray[faceoffset+7] = t4, geoGroup.__faceArray[faceoffset+8] = fromi;
        geoGroup.__faceArray[faceoffset+9] = t2b, geoGroup.__faceArray[faceoffset+10] = toi, geoGroup.__faceArray[faceoffset+11] = t3b;
        //final arrowhead
        geoGroup.__faceArray[faceoffset+12] = t2b, geoGroup.__faceArray[faceoffset+13] = endi, geoGroup.__faceArray[faceoffset+14] = t3b;
        
        //sides
        geoGroup.__lineArray[lineoffset] = t1, geoGroup.__lineArray[lineoffset+1] = t2;
        geoGroup.__lineArray[lineoffset+2] = t1, geoGroup.__lineArray[lineoffset+3] = t4;           
        //geoGroup.__lineArray[lineoffset+4] = t2, geoGroup.__lineArray[lineoffset+5] = t3;
        geoGroup.__lineArray[lineoffset+4] = t3, geoGroup.__lineArray[lineoffset+5] = t4;
        //caps
        geoGroup.__lineArray[lineoffset+6] = t1, geoGroup.__lineArray[lineoffset+7] = t4;
        //geoGroup.__lineArray[lineoffset+10] = t1, geoGroup.__lineArray[lineoffset+11] = fromi;           
        //geoGroup.__lineArray[lineoffset+12] = t4, geoGroup.__lineArray[lineoffset+13] = fromi;

        geoGroup.__lineArray[lineoffset+8] = t2b, geoGroup.__lineArray[lineoffset+9] = t2; //toi        
        geoGroup.__lineArray[lineoffset+10] = t2b, geoGroup.__lineArray[lineoffset+11] = t3b;
        geoGroup.__lineArray[lineoffset+12] = t3, geoGroup.__lineArray[lineoffset+13] = t3b; //toi
        //arrowhead
        geoGroup.__lineArray[lineoffset+14] = t2b, geoGroup.__lineArray[lineoffset+15] = endi;
        geoGroup.__lineArray[lineoffset+16] = t2b, geoGroup.__lineArray[lineoffset+17] = t3b;
        geoGroup.__lineArray[lineoffset+18] = endi, geoGroup.__lineArray[lineoffset+19] = t3b; 
        
        geoGroup.faceidx += 15;        
        geoGroup.lineidx += 20;    

        
    };
    
    //handles custom shape generation from user supplied arrays
    //May need to generate normal and/or line indices
    var drawCustom = function(shape, geoGroup, customSpec) {
        
        var vertexArr = customSpec.vertexArr, normalArr = customSpec.normalArr, faceArr = customSpec.faceArr, lineArr = customSpec.lineArr;        
        
        if (vertexArr.length === 0 || faceArr.length === 0) {
            console.error("Error adding custom shape component: No vertices and/or face indices supplied!");
        }
        
        geoGroup.vertices = vertexArr.length, geoGroup.faceidx = faceArr.length;
        
        var offset, v, a, b, c;
        
        for (var i = 0; i < geoGroup.vertices; ++i) {            
            offset = i*3;
            v = vertexArr[i];    
            geoGroup.__vertexArray[offset] = v.x, geoGroup.__vertexArray[offset+1] = v.y, geoGroup.__vertexArray[offset+2] = v.z;           
        }
        
        for (var i = 0; i < geoGroup.faceidx / 3; ++i) {
            offset = i*3;
            a = faceArr[offset], b = faceArr[offset+1], c = faceArr[offset+2];
            var vA = new WebMol.Vector3(), vB = new WebMol.Vector3(), vC = new WebMol.Vector3();
            shape.intersectionShape.triangle.push( new WebMol.Triangle( vA.copy(vertexArr[a]), vB.copy(vertexArr[b]), vC.copy(vertexArr[c]) ) );
        }
                  
        geoGroup.__faceArray = new Uint16Array(faceArr);
     
        geoGroup.truncateArrayBuffers(true);
        
        if (normalArr.length < geoGroup.vertices)
            geoGroup.setNormals();
        else {
            
            geoGroup.__normalArray = new Float32Array(geoGroup.vertices*3);
            var n;
            for (var i = 0; i < geoGroup.vertices; ++i) {
                offset = i*3;
                n = normalArr[i];
                geoGroup.__normalArray[offset] = n.x, geoGroup.__normalArray[offset+1] = n.y, geoGroup.__normalArray[offset+2] = n.z;
            }
        }
            
        
        if (! lineArr.length)
            geoGroup.setLineIndices(); 
        else
            geoGroup.__lineArray = new Uint16Array(lineArr);
        
        geoGroup.lineidx = geoGroup.__lineArray.length;
             
    };
    
    //Read a cube file - generate model and possibly shape(s)
    var parseCube = function(shape, str, isoval, voxel) {
        
        var smooth = !(voxel);
        
        var lines = str.replace(/^\s+/, "").split(/[\n\r]+/);
        
        if (lines.length < 6)
            return;
            
        var lineArr = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");       
          
        var natoms = Math.abs(parseFloat(lineArr[0]));        
        var origin = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3]));
        
        lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        
        //might have to convert from bohr units to angstroms
        var convFactor = (parseFloat(lineArr[0]) > 0) ? 0.529177 : 1;
        
        origin.multiplyScalar(convFactor);
        
        var nX = Math.abs(lineArr[0]);
        var xVec = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3])).multiplyScalar(convFactor);
        
        lineArr = lines[4].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        
        var nY = Math.abs(lineArr[0]);
        var yVec = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3])).multiplyScalar(convFactor);
        
        lineArr = lines[5].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        
        var nZ = Math.abs(lineArr[0]);
        var zVec = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3])).multiplyScalar(convFactor);
        
        //lines.splice(6, natoms).join("\n");
        
        lines = new Float32Array(lines.splice(natoms+7).join(" ").replace(/^\s+/, "").split(/[\s\r]+/));        
        
        var vertnums = new Int16Array(nX*nY*nZ);        
        for (var i = 0; i < vertnums.length; ++i)
            vertnums[i] = -1;

        var bitdata = new Uint8Array(nX*nY*nZ);
        
        for (var i = 0; i < lines.length; ++i) {
            var val = (isoval >= 0) ? lines[i] - isoval : isoval - lines[i];
            
            if (val > 0)
                bitdata[i] |= ISDONE;
            
        }
        
        var verts = [], faces = [];
        
        WebMol.MarchingCube.march(bitdata, verts, faces, {
            fulltable : true,
            scale : xVec.length(),
            origin : origin,
            nX : nX,
            nY : nY,
            nZ : nZ        
        });
        
        if (smooth)
            WebMol.MarchingCube.laplacianSmooth(10, verts, faces);        
        
        shape.addCustom({vertexArr:verts, 
                         faceArr:faces});
                
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
    

    var GLShape = function(sid, stylespec) {
        
        stylespec = stylespec || {};
        WebMol.ShapeIDCount++;
        this.id = sid || WebMol.ShapeIDCount;
        this.color = stylespec.color || new WebMol.Color();
        this.wireframe = stylespec.wireframe ? true : false;
        this.alpha = stylespec.alpha ? WebMol.Math.clamp(stylespec.alpha, 0.0, 1.0) : 1.0;
        this.side = (stylespec.side !== undefined) ? stylespec.side : WebMol.DoubleSide;
        
        this.boundingSphere = new WebMol.Sphere();
        
        //Click handling
        this.clickable = stylespec.clickable ? true : false;
        this.callback = typeof(stylespec.callback) === "function" ? stylespec.callback : null;
        this.intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
        
        //Keep track of shape components and their centroids
        var components = [];
        var shapeObj = null;
        var renderedShapeObj = null;
        
        var geo = new WebMol.Geometry(true);
        
        this.addCustom = function(customSpec) {
            
            customSpec.vertexArr = customSpec.vertexArr || [];
            customSpec.faceArr = customSpec.faceArr || [];
            customSpec.normalArr = customSpec.normalArr || [];
            customSpec.lineArr = customSpec.lineArr || [];
            
            //Force creation of new geometryGroup for each added component
            var geoGroup = geo.addGeoGroup();
            drawCustom(this, geoGroup, customSpec);            
            geoGroup.truncateArrayBuffers(true);
    
            for (var i = 0; i < geoGroup.__colorArray.length / 3; ++i) {
                geoGroup.__colorArray[i*3] = this.color.r;
                geoGroup.__colorArray[i*3 + 1] = this.color.g;
                geoGroup.__colorArray[i*3 + 2] = this.color.b;      
            }            
            
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup,
                centroid : geoGroup.getCentroid() 
            });
            
            updateBoundingFromPoints( this.boundingSphere, components, geoGroup.__vertexArray );            
        };
        
        //TODO: Refactor so 'drawSphere' method automatically updates bounding sphere as vertices are added
        this.addSphere = function(sphereSpec) {
          
            sphereSpec.center = sphereSpec.center || {x: 0, y: 0, z: 0};
            sphereSpec.radius = sphereSpec.radius ? WebMol.Math.clamp(sphereSpec.radius, 0, Infinity) : 1.5;
            
            var geoGroup = geo.addGeoGroup();
            drawSphere(this, geoGroup, sphereSpec);
            geoGroup.truncateArrayBuffers(true);            
            
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup, //has to be last group added
                centroid : new WebMol.Vector3(sphereSpec.center.x, sphereSpec.center.y, sphereSpec.center.z)
            });
            
            updateBoundingFromPoints( this.boundingSphere, components, geoGroup.__vertexArray );
        };        
        
        this.addArrow = function(arrowSpec) {
            
            arrowSpec.start = arrowSpec.start || {};
            arrowSpec.end = arrowSpec.end || {};
            
            arrowSpec.start = new WebMol.Vector3(arrowSpec.start.x || 0, arrowSpec.start.y || 0, arrowSpec.start.z || 0);
            arrowSpec.end = new WebMol.Vector3(arrowSpec.end.x || 3, arrowSpec.end.y || 0, arrowSpec.end.z || 0);
            arrowSpec.radius = arrowSpec.radius || 0.25;
            
            arrowSpec.radiusRatio = arrowSpec.radiusRatio || 1.618034;
            arrowSpec.mid = ( 0 < arrowSpec.mid && arrowSpec.mid < 1) ? arrowSpec.mid : 0.618034;
            
            var geoGroup = geo.addGeoGroup();
            
            drawArrow(this, geoGroup, arrowSpec);
            geoGroup.truncateArrayBuffers(true);
            
            var centroid = new WebMol.Vector3();
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup,
                centroid : centroid.addVectors(arrowSpec.start, arrowSpec.end).multiplyScalar(0.5)
            });
            
            updateBoundingFromPoints( this.boundingSphere, components, geoGroup.__vertexArray );
                            
        };
        
        this.addVolumetricData = function(str, fmt, isoval, vxl) {
            isoval = (typeof(isoval) === "number") ? isoval : 0.0;
            vxl = (vxl) ? true : false;
            
            switch(fmt) {
                case "cube":
                    parseCube(this, str, isoval, vxl);
                    break;
            }              
            
        };
    
        //TODO: Adding multiple overlapping shapes in wireframe mode should obscure overlapped meshes ...
        //         (i.e. generate a new isosurface ... ?)
        this.globj = function(group) {
            
            geo.initTypedArrays();
            
            updateColor(geo, this.color);
            
            shapeObj = new WebMol.Object3D();
            var material = new WebMol.MeshLambertMaterial({
                wireframe : this.wireframe,
                vertexColors : true,
                ambient : 0x000000,
                reflectivity : 0,
                side : this.side,
                transparent : (this.alpha < 1) ? true : false,
                opacity : this.alpha
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