/**
 * A GLShape is a collection of user specified shapes.
 * 
 * @constructor $3Dmol.GLShape
 * @extends {ShapeSpec}
 * @param {number} sid - Unique identifier
 * @param {ShapeSpec} stylespec - shape style specification
 */
$3Dmol.GLShape = (function() {

    // Marching cube, to match with protein surface generation
    var ISDONE = 2;

    
    var finalizeGeo = function(geo) {
        //to avoid creating a bunch of geometries, we leave geoGroup untruncated
        //until render is called, at which point we truncate; 
        //successive called up updateGeo will return a new geometry
        var geoGroup = geo.updateGeoGroup(0);
        if(geoGroup.vertices > 0) {
            geoGroup.truncateArrayBuffers(true, true);
        }
    };
    
    /**
     * 
     * @param {$3Dmol.Geometry}
     *            geo
     * @param {$3Dmol.Color |
     *            colorlike} color
     */
    var updateColor = function(geo, color) {

        var C = color || $3Dmol.CC.color(color);
        geo.colorsNeedUpdate = true;
        
        var r,g,b;
        if(! (color.constructor === Array)) {
            r = color.r;
            g = color.g;
            b = color.b;
        }


        for ( var gg in geo.geometryGroups) {

            var geoGroup = geo.geometryGroups[gg];
            var colorArr = geoGroup.colorArray;

            for (var i = 0, il = geoGroup.vertices; i < il; ++i) {
            
                if( color.constructor === Array) {
                    var c = color[i];
                    r = c.r;
                    g = c.g;
                    b = c.b;
                }

                colorArr[i * 3] = r;
                colorArr[i * 3 + 1] = g;
                colorArr[i * 3 + 2] = b;
            }
        }

    };


    /**
     * @param {$3Dmol.GLShape}
     *            shape
     * @param {geometryGroup}
     *            geoGroup
     * @param {ArrowSpec}
     *            spec
     */
    var drawArrow = function(shape, geo, spec) {

        var from = spec.start, end = spec.end, radius = spec.radius, radiusRatio = spec.radiusRatio, mid = spec.mid;

        if (!(from && end))
            return;

        var geoGroup = geo.updateGeoGroup(51);

        // vertices

        var dir = end.clone();
        dir.sub(from).multiplyScalar(mid);
        var to = from.clone().add(dir);
        var negDir = dir.clone().negate();

        shape.intersectionShape.cylinder.push(new $3Dmol.Cylinder(from.clone(),
                to.clone(), radius));
        shape.intersectionShape.sphere.push(new $3Dmol.Sphere(from.clone(),
                radius));

        // get orthonormal vector
        var nvecs = [];
        nvecs[0] = dir.clone();
        if (Math.abs(nvecs[0].x) > 0.0001)
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

        var start = geoGroup.vertices;
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;
        var faceArray = geoGroup.faceArray;
        var normalArray = geoGroup.normalArray;
        var lineArray = geoGroup.lineArray;

        var offset, i, n;
        // add vertices, opposing vertices paired together
        for (i = 0, n = nvecs.length; i < n; ++i) {
            offset = 3 * (start + 3 * i);
            var bottom = nvecs[i].clone().multiplyScalar(radius).add(from);
            var top = nvecs[i].clone().multiplyScalar(radius).add(to);
            var conebase = nvecs[i].clone()
                    .multiplyScalar(radius * radiusRatio).add(to);

            vertexArray[offset] = bottom.x;
            vertexArray[offset + 1] = bottom.y;
            vertexArray[offset + 2] = bottom.z;

            vertexArray[offset + 3] = top.x;
            vertexArray[offset + 4] = top.y;
            vertexArray[offset + 5] = top.z;

            vertexArray[offset + 6] = conebase.x;
            vertexArray[offset + 7] = conebase.y;
            vertexArray[offset + 8] = conebase.z;

            if (i > 0) {
                var prev_x = vertexArray[offset - 3];
                var prev_y = vertexArray[offset - 2];
                var prev_z = vertexArray[offset - 1];

                var c = new $3Dmol.Vector3(prev_x, prev_y, prev_z);
                var b = end.clone(), b2 = to.clone();
                var a = new $3Dmol.Vector3(conebase.x, conebase.y, conebase.z);

                shape.intersectionShape.triangle.push(new $3Dmol.Triangle(a, b,
                        c));
                shape.intersectionShape.triangle.push(new $3Dmol.Triangle(c
                        .clone(), b2, a.clone()));
            }
        }

        geoGroup.vertices += 48;
        offset = geoGroup.vertices * 3;

        // caps
        vertexArray[offset] = from.x;
        vertexArray[offset + 1] = from.y;
        vertexArray[offset + 2] = from.z;

        vertexArray[offset + 3] = to.x;
        vertexArray[offset + 4] = to.y;
        vertexArray[offset + 5] = to.z;

        vertexArray[offset + 6] = end.x;
        vertexArray[offset + 7] = end.y;
        vertexArray[offset + 8] = end.z;

        geoGroup.vertices += 3;

        // now faces
        var face, norm, faceoffset, lineoffset;
        var t1, t2, t2b, t3, t3b, t4, t1offset, t2offset, t2boffset, t3offset, t3boffset, t4offset;
        var n1, n2, n3, n4;
        var n_vertices = 0;
        var fromi = geoGroup.vertices - 3, toi = geoGroup.vertices - 2, endi = geoGroup.vertices - 1;
        var fromoffset = fromi * 3, tooffset = toi * 3, endoffset = endi * 3;
        for (i = 0, n = nvecs.length - 1; i < n; ++i) {

            var ti = start + 3 * i;
            offset = ti * 3;
            faceoffset = geoGroup.faceidx;
            lineoffset = geoGroup.lineidx;

            t1 = ti;
            t1offset = t1 * 3;
            t2 = ti + 1;
            t2offset = t2 * 3;
            t2b = ti + 2;
            t2boffset = t2b * 3;
            t3 = ti + 4;
            t3offset = t3 * 3;
            t3b = ti + 5;
            t3boffset = t3b * 3;
            t4 = ti + 3;
            t4offset = t4 * 3;

            // face = [t1, t2, t4], [t2, t3, t4];
            // face = [t1, t2, t3, t4];

            norm = [ nvecs[i], nvecs[i], nvecs[i + 1], nvecs[i + 1] ];

            n1 = n2 = nvecs[i];
            n3 = n4 = nvecs[i + 1];

            normalArray[t1offset] = n1.x;
            normalArray[t2offset] = n2.x;
            normalArray[t4offset] = n4.x;
            normalArray[t1offset + 1] = n1.y;
            normalArray[t2offset + 1] = n2.y;
            normalArray[t4offset + 1] = n4.y;
            normalArray[t1offset + 2] = n1.z;
            normalArray[t2offset + 2] = n2.z;
            normalArray[t4offset + 2] = n4.z;

            normalArray[t2offset] = n2.x;
            normalArray[t3offset] = n3.x;
            normalArray[t4offset] = n4.x;
            normalArray[t2offset + 1] = n2.y;
            normalArray[t3offset + 1] = n3.y;
            normalArray[t4offset + 1] = n4.y;
            normalArray[t2offset + 2] = n2.z;
            normalArray[t3offset + 2] = n3.z;
            normalArray[t4offset + 2] = n4.z;

            normalArray[t2boffset] = n2.x;
            normalArray[t3boffset] = n3.x;
            normalArray[t2boffset + 1] = n2.y;
            normalArray[t3boffset + 1] = n3.y;
            normalArray[t2boffset + 2] = n2.z;
            normalArray[t3boffset + 2] = n3.z;

            // sides
            faceArray[faceoffset] = t1;
            faceArray[faceoffset + 1] = t2;
            faceArray[faceoffset + 2] = t4;
            faceArray[faceoffset + 3] = t2;
            faceArray[faceoffset + 4] = t3;
            faceArray[faceoffset + 5] = t4;
            // caps
            faceArray[faceoffset + 6] = t1;
            faceArray[faceoffset + 7] = t4;
            faceArray[faceoffset + 8] = fromi;
            faceArray[faceoffset + 9] = t2b;
            faceArray[faceoffset + 10] = toi;
            faceArray[faceoffset + 11] = t3b;
            // arrowhead
            faceArray[faceoffset + 12] = t2b;
            faceArray[faceoffset + 13] = endi;
            faceArray[faceoffset + 14] = t3b;

            // sides
            lineArray[lineoffset] = t1;
            lineArray[lineoffset + 1] = t2;
            lineArray[lineoffset + 2] = t1;
            lineArray[lineoffset + 3] = t4;
            // lineArray[lineoffset+4] = t2, lineArray[lineoffset+5] = t3;
            lineArray[lineoffset + 4] = t3;
            lineArray[lineoffset + 5] = t4;
            // caps
            lineArray[lineoffset + 6] = t1;
            lineArray[lineoffset + 7] = t4;
            // lineArray[lineoffset+10] = t1, lineArray[lineoffset+11] = fromi;
            // lineArray[lineoffset+12] = t4, lineArray[lineoffset+13] = fromi;

            lineArray[lineoffset + 8] = t2b;
            lineArray[lineoffset + 9] = t2; // toi
            lineArray[lineoffset + 10] = t2b;
            lineArray[lineoffset + 11] = t3b;
            lineArray[lineoffset + 12] = t3;
            lineArray[lineoffset + 13] = t3b; // toi
            // arrowhead
            lineArray[lineoffset + 14] = t2b;
            lineArray[lineoffset + 15] = endi;
            lineArray[lineoffset + 16] = t2b;
            lineArray[lineoffset + 17] = t3b;
            lineArray[lineoffset + 18] = endi;
            lineArray[lineoffset + 19] = t3b;

            geoGroup.faceidx += 15;
            geoGroup.lineidx += 20;

        }
        // final face

        face = [ start + 45, start + 46, start + 1, start, start + 47,
                start + 2 ];
        norm = [ nvecs[15], nvecs[15], nvecs[0], nvecs[0] ];

        faceoffset = geoGroup.faceidx;
        lineoffset = geoGroup.lineidx;

        t1 = face[0];
        t1offset = t1 * 3;
        t2 = face[1];
        t2offset = t2 * 3;
        t2b = face[4];
        t2boffset = t2b * 3;
        t3 = face[2];
        t3offset = t3 * 3;
        t3b = face[5];
        t3boffset = t3b * 3;
        t4 = face[3];
        t4offset = t4 * 3;

        n1 = n2 = nvecs[15];
        n3 = n4 = nvecs[0];

        normalArray[t1offset] = n1.x;
        normalArray[t2offset] = n2.x;
        normalArray[t4offset] = n4.x;
        normalArray[t1offset + 1] = n1.y;
        normalArray[t2offset + 1] = n2.y;
        normalArray[t4offset + 1] = n4.y;
        normalArray[t1offset + 2] = n1.z;
        normalArray[t2offset + 2] = n2.z;
        normalArray[t4offset + 2] = n4.z;

        normalArray[t2offset] = n2.x;
        normalArray[t3offset] = n3.x;
        normalArray[t4offset] = n4.x;
        normalArray[t2offset + 1] = n2.y;
        normalArray[t3offset + 1] = n3.y;
        normalArray[t4offset + 1] = n4.y;
        normalArray[t2offset + 2] = n2.z;
        normalArray[t3offset + 2] = n3.z;
        normalArray[t4offset + 2] = n4.z;

        normalArray[t2boffset] = n2.x;
        normalArray[t3boffset] = n3.x;
        normalArray[t2boffset + 1] = n2.y;
        normalArray[t3boffset + 1] = n3.y;
        normalArray[t2boffset + 2] = n2.z;
        normalArray[t3boffset + 2] = n3.z;

        // Cap normals
        dir.normalize();
        negDir.normalize();
        normalArray[fromoffset] = negDir.x;
        normalArray[tooffset] = normalArray[endoffset] = dir.x;
        normalArray[fromoffset + 1] = negDir.y;
        normalArray[tooffset + 1] = normalArray[endoffset + 1] = dir.y;
        normalArray[fromoffset + 2] = negDir.z;
        normalArray[tooffset + 2] = normalArray[endoffset + 2] = dir.z;

        // Final side
        faceArray[faceoffset] = t1;
        faceArray[faceoffset + 1] = t2;
        faceArray[faceoffset + 2] = t4;
        faceArray[faceoffset + 3] = t2;
        faceArray[faceoffset + 4] = t3;
        faceArray[faceoffset + 5] = t4;
        // final caps
        faceArray[faceoffset + 6] = t1;
        faceArray[faceoffset + 7] = t4;
        faceArray[faceoffset + 8] = fromi;
        faceArray[faceoffset + 9] = t2b;
        faceArray[faceoffset + 10] = toi;
        faceArray[faceoffset + 11] = t3b;
        // final arrowhead
        faceArray[faceoffset + 12] = t2b;
        faceArray[faceoffset + 13] = endi;
        faceArray[faceoffset + 14] = t3b;

        // sides
        lineArray[lineoffset] = t1;
        lineArray[lineoffset + 1] = t2;
        lineArray[lineoffset + 2] = t1;
        lineArray[lineoffset + 3] = t4;
        // lineArray[lineoffset+4] = t2, lineArray[lineoffset+5] = t3;
        lineArray[lineoffset + 4] = t3;
        lineArray[lineoffset + 5] = t4;
        // caps
        lineArray[lineoffset + 6] = t1;
        lineArray[lineoffset + 7] = t4;
        // lineArray[lineoffset+10] = t1, lineArray[lineoffset+11] = fromi;
        // lineArray[lineoffset+12] = t4, lineArray[lineoffset+13] = fromi;

        lineArray[lineoffset + 8] = t2b;
        lineArray[lineoffset + 9] = t2; // toi
        lineArray[lineoffset + 10] = t2b;
        lineArray[lineoffset + 11] = t3b;
        lineArray[lineoffset + 12] = t3;
        lineArray[lineoffset + 13] = t3b; // toi
        // arrowhead
        lineArray[lineoffset + 14] = t2b;
        lineArray[lineoffset + 15] = endi;
        lineArray[lineoffset + 16] = t2b;
        lineArray[lineoffset + 17] = t3b;
        lineArray[lineoffset + 18] = endi;
        lineArray[lineoffset + 19] = t3b;

        geoGroup.faceidx += 15;
        geoGroup.lineidx += 20;

    };

    //helper function for adding an appropriately sized mesh
    var addCustomGeo = function(shape, geo, mesh, color, clickable) {
        var geoGroup = geo.addGeoGroup();
        var vertexArr = mesh.vertexArr, normalArr = mesh.normalArr, 
            faceArr = mesh.faceArr;

        geoGroup.vertices = vertexArr.length;
        geoGroup.faceidx = faceArr.length;

        var offset, v, a, b, c, i, il;
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;
        
        if(! (color.constructor === Array)) {
            var r = color.r;
            var g = color.g;
            var b = color.b;
        }
        for (i = 0, il = geoGroup.vertices; i < il; ++i) {
            offset = i * 3;
            v = vertexArr[i];
            vertexArray[offset] = v.x;
            vertexArray[offset + 1] = v.y;
            vertexArray[offset + 2] = v.z;
            
            if( color.constructor === Array) {
                var c = color[i];
                var r = c.r;
                var g = c.g;
                var b = c.b;
            }
            
            colorArray[offset] = r;
            colorArray[offset + 1] = g;
            colorArray[offset + 2] = b;
        }
        
        if(clickable) {
            for (i = 0, il = geoGroup.faceidx / 3; i < il; ++i) {
                offset = i * 3;
                a = faceArr[offset];
                b = faceArr[offset + 1];
                c = faceArr[offset + 2];
                var vA = new $3Dmol.Vector3(), vB = new $3Dmol.Vector3(), vC = new $3Dmol.Vector3();
                shape.intersectionShape.triangle.push(new $3Dmol.Triangle(vA
                        .copy(vertexArr[a]), vB.copy(vertexArr[b]), vC
                        .copy(vertexArr[c])));
            }
        }
        
        if(clickable) {
            
            var center = new $3Dmol.Vector3(0,0,0);
            var cnt = 0;
            for(var g = 0; g < geo.geometryGroups.length; g++) {
                center.add(geo.geometryGroups[g].getCentroid());
                cnt++;
            }
            center.divideScalar(cnt);
            
            
            updateBoundingFromPoints(shape.boundingSphere, {centroid: center}, vertexArray);
        }

        geoGroup.faceArray = new Uint16Array(faceArr);

        geoGroup.truncateArrayBuffers(true, true);

        if (normalArr.length < geoGroup.vertices)
            geoGroup.setNormals();
        else {

            var normalArray = geoGroup.normalArray = new Float32Array(geoGroup.vertices * 3);
            var n;
            for (i = 0, il = geoGroup.vertices; i < il; ++i) {
                offset = i * 3;
                n = normalArr[i];
                normalArray[offset] = n.x;
                normalArray[offset + 1] = n.y;
                normalArray[offset + 2] = n.z;
            }
        }
        
        geoGroup.setLineIndices();
        geoGroup.lineidx = geoGroup.lineArray.length;
    };
    

    
    // handles custom shape generation from user supplied arrays
    // May need to generate normal and/or line indices
    /**
     * @param {$3Dmol.GLShape}
     *            shape
     * @param {geometry}
     *            geo
     * @param {CustomSpec}
     *            customSpec
     */
    var drawCustom = function(shape, geo, customSpec) {
        var mesh = customSpec;
        var vertexArr = mesh.vertexArr, normalArr = mesh.normalArr, 
        faceArr = mesh.faceArr;
        if (vertexArr.length === 0 || faceArr.length === 0) {
            console
                    .warn("Error adding custom shape component: No vertices and/or face indices supplied!");
        }

        var color = customSpec.color;
        if(typeof(color) == 'undefined') {
            color = shape.color;
        }
        color =  $3Dmol.CC.color(color);

        //var firstgeo = geo.geometryGroups.length;
        var splits = $3Dmol.splitMesh(mesh);
        for(var i = 0, n = splits.length; i < n; i++) {
            addCustomGeo(shape, geo, splits[i], splits[i].colorArr ? splits[i].colorArr : color, customSpec.clickable);
        } 
    }; 

    // Update a bounding sphere's position and radius
    // from list of centroids and new points
    /**
     * @param {$3Dmol.Sphere}
     *            sphere
     * @param {Object}
     *            components
     * @param {Array}
     *            points
     */
    var updateBoundingFromPoints = function(sphere, components, points) {

        sphere.center.set(0, 0, 0);

        var i, il;

        if (components.length > 0) {

            for (i = 0, il = components.length; i < il; ++i) {
                var centroid = components[i].centroid;
                sphere.center.add(centroid);
            }

            sphere.center.divideScalar(components.length);
        }

        var maxRadiusSq = sphere.radius * sphere.radius;

        for (i = 0, il = points.length / 3; i < il; i++) {
            var x = points[i * 3], y = points[i * 3 + 1], z = points[i * 3 + 2];
            var radiusSq = sphere.center.distanceToSquared({
                x : x,
                y : y,
                z : z
            });
            maxRadiusSq = Math.max(maxRadiusSq, radiusSq);
        }

        sphere.radius = Math.sqrt(maxRadiusSq);

    };

    /**
     * 
     * @param {$3Dmol.GLShape}
     *            shape
     * @param {ShapeSpec}
     *            stylespec
     * @returns {undefined}
     */
    var updateFromStyle = function(shape, stylespec) {
        if(typeof(stylespec.color) != 'undefined') {
            shape.color = stylespec.color || new $3Dmol.Color();
            if(! (stylespec.color instanceof $3Dmol.Color))
                shape.color = $3Dmol.CC.color(stylespec.color);
        }
        shape.wireframe = stylespec.wireframe ? true : false;
        //opacity is the preferred nomenclature, support alpha for backwards compat
        shape.opacity = stylespec.alpha ? $3Dmol.Math.clamp(stylespec.alpha, 0.0,
                1.0) : 1.0;
        if(typeof(stylespec.opacity) != 'undefined') {
            shape.opacity = $3Dmol.Math.clamp(stylespec.opacity, 0.0, 1.0);
        }
        shape.side = (stylespec.side !== undefined) ? stylespec.side
                : $3Dmol.DoubleSide;

        shape.linewidth = typeof(stylespec.linewidth) == 'undefined' ? 1 : stylespec.linewidth;
        // Click handling
        shape.clickable = stylespec.clickable ? true : false;
        shape.callback = typeof (stylespec.callback) === "function" ? stylespec.callback
                : null;
        shape.hidden = stylespec.hidden;
    };

    /**
     * Custom renderable shape
     * 
     * @constructor $3Dmol.GLShape
     * 
     * @param {Object}
     *            stylespec
     * @returns {$3Dmol.GLShape}
     */
    function GLShape(stylespec) {

        stylespec = stylespec || {};
        $3Dmol.ShapeIDCount++;

        this.boundingSphere = new $3Dmol.Sphere();
        /** @type {IntersectionShapes} */
        this.intersectionShape = {
            sphere : [],
            cylinder : [],
            line : [],
            triangle : []
        };

        updateFromStyle(this, stylespec);

        // Keep track of shape components and their centroids
        var components = [];
        var shapeObj = null;
        var renderedShapeObj = null;

        var geo = new $3Dmol.Geometry(true);
        var linegeo = new $3Dmol.Geometry(true);

        /** Update shape with new style specification
         * @param {ShapeSpec} newspec
         * @return {$3Dmol.GLShape}
         */
        this.updateStyle = function(newspec) {

            for ( var prop in newspec) {
                stylespec[prop] = newspec[prop];
            }

            updateFromStyle(this, stylespec);
        };

        /**
         * Creates a custom shape from supplied vertex and face arrays
         * @function $3Dmol.GLShape#addCustom
         * @param {CustomSpec} customSpec
         * @return {$3Dmol.GLShape}
         */
        this.addCustom = function(customSpec) {

            customSpec.vertexArr = customSpec.vertexArr || [];
            customSpec.faceArr = customSpec.faceArr || [];
            customSpec.normalArr = customSpec.normalArr || [];

            var firstgeo = geo.geometryGroups.length;
            // will split mesh as needed
            drawCustom(this, geo, customSpec);
        };

        /**
         * Creates a sphere shape
         * @function $3Dmol.GLShape#addSphere
         * @param {SphereSpec} sphereSpec
         * @return {$3Dmol.GLShape}
         */
        this.addSphere = function(sphereSpec) {

            sphereSpec.center = sphereSpec.center || {
                x : 0,
                y : 0,
                z : 0
            };
            sphereSpec.radius = sphereSpec.radius ? $3Dmol.Math.clamp(
                    sphereSpec.radius, 0, Infinity) : 1.5;
            sphereSpec.color = $3Dmol.CC.color(sphereSpec.color);
            
            this.intersectionShape.sphere.push(new $3Dmol.Sphere(
                    sphereSpec.center, sphereSpec.radius));

            $3Dmol.GLDraw.drawSphere(geo, sphereSpec.center,
                    sphereSpec.radius, sphereSpec.color);

            components.push({
                centroid : new $3Dmol.Vector3(sphereSpec.center.x,
                        sphereSpec.center.y, sphereSpec.center.z)
            });
            var geoGroup = geo.updateGeoGroup(0);
            
            updateBoundingFromPoints(this.boundingSphere, components,
                    geoGroup.vertexArray);
        };

        /**
         * Creates a cylinder shape
         * @function $3Dmol.GLShape#addCylinder
         * @param {CylinderSpec} cylinderSpec
         * @return {$3Dmol.GLShape}
         */
        this.addCylinder = function(cylinderSpec) {

            cylinderSpec.start = cylinderSpec.start || {};
            cylinderSpec.end = cylinderSpec.end || {};

            var start = new $3Dmol.Vector3(cylinderSpec.start.x || 0,
                    cylinderSpec.start.y || 0, cylinderSpec.start.z || 0);
            var end = new $3Dmol.Vector3(cylinderSpec.end.x,
                    cylinderSpec.end.y || 0, cylinderSpec.end.z || 0);
        	if(typeof(end.x) == 'undefined') end.x = 3; //show something even if undefined

            var radius = cylinderSpec.radius || 0.1;
            var color = $3Dmol.CC.color(cylinderSpec.color);
            
            this.intersectionShape.cylinder.push(new $3Dmol.Cylinder(start, end, radius));

            $3Dmol.GLDraw.drawCylinder(geo, start, end, radius, color, cylinderSpec.fromCap, cylinderSpec.toCap);            
            
            var centroid = new $3Dmol.Vector3();
            components.push({
                centroid : centroid.addVectors(start,end).multiplyScalar(0.5)
            });
            var geoGroup = geo.updateGeoGroup(0);
            updateBoundingFromPoints(this.boundingSphere, components,
                    geoGroup.vertexArray);

        };

        /**
         * Creates a line shape
         * @function $3Dmol.GLShape#addLine         
         * @param {LineSpec} lineSpec
         * @return {$3Dmol.GLShape}
         */
        this.addLine = function(lineSpec) {
            lineSpec.start = lineSpec.start || {};
            lineSpec.end = lineSpec.end || {};

            var start = new $3Dmol.Vector3(lineSpec.start.x || 0,
                    lineSpec.start.y || 0, lineSpec.start.z || 0);
            var end = new $3Dmol.Vector3(lineSpec.end.x,
                    lineSpec.end.y || 0, lineSpec.end.z || 0);            
            if(typeof(end.x) == 'undefined') end.x = 3; //show something even if undefined

            var geoGroup = geo.updateGeoGroup(2);

            //make line from start to end
            //for consistency with rest of shapes, uses vertices and lines rather
            //than a separate line geometry
            var vstart = geoGroup.vertices;
            var i = vstart*3;
            var vertexArray = geoGroup.vertexArray;
            vertexArray[i] = start.x;
            vertexArray[i+1] = start.y;
            vertexArray[i+2] = start.z;
            vertexArray[i+3] = end.x;
            vertexArray[i+4] = end.y;
            vertexArray[i+5] = end.z;
            geoGroup.vertices += 2;
            
            var lineArray = geoGroup.lineArray;
            var li =  geoGroup.lineidx;
            lineArray[li] = vstart;
            lineArray[li+1] = vstart+1;
            geoGroup.lineidx += 2;
            
        }
        /**
         * Creates an arrow shape
         * @function $3Dmol.GLShape#addArrow        
         * @param {ArrowSpec} arrowSpec
         * @return {$3Dmol.GLShape}
         */
        this.addArrow = function(arrowSpec) {

            arrowSpec.start = arrowSpec.start || {};
            arrowSpec.end = arrowSpec.end || {};

            arrowSpec.start = new $3Dmol.Vector3(arrowSpec.start.x || 0,
                    arrowSpec.start.y || 0, arrowSpec.start.z || 0);

            if (arrowSpec.dir instanceof $3Dmol.Vector3
                    && arrowSpec.length instanceof number) {
                var end = arrowSpec.dir.clone().multiplyScalar(arrowSpec.length).add(
                        start);
                arrowSpec.end = end;
            }

            else {
                arrowSpec.end = new $3Dmol.Vector3(arrowSpec.end.x,
                        arrowSpec.end.y || 0, arrowSpec.end.z || 0);
            	if(typeof(arrowSpec.end.x) == 'undefined') arrowSpec.end.x = 3; //show something even if undefined
            }

            arrowSpec.radius = arrowSpec.radius || 0.1;

            arrowSpec.radiusRatio = arrowSpec.radiusRatio || 1.618034;
            arrowSpec.mid = (0 < arrowSpec.mid && arrowSpec.mid < 1) ? arrowSpec.mid
                    : 0.618034;


            drawArrow(this, geo, arrowSpec);

            var centroid = new $3Dmol.Vector3();
            components.push({
                centroid : centroid.addVectors(arrowSpec.start, arrowSpec.end)
                        .multiplyScalar(0.5)
            });
            var geoGroup = geo.updateGeoGroup(0);
            updateBoundingFromPoints(this.boundingSphere, components,
                    geoGroup.vertexArray);

        };
        
        /**
         * Create isosurface from voluemetric data.
         * @function $3Dmol.GLShape#addIsosurface         
         * @param {$3Dmol.VolumeData} data - volumetric input data
         * @param {IsoSurfaceSpec} isoSpec - volumetric data shape specification
         */
        this.addIsosurface = function(data, volSpec) {
            var isoval = (volSpec.isoval !== undefined && typeof (volSpec.isoval) === "number") ? volSpec.isoval
                    : 0.0;
            var voxel = (volSpec.voxel) ? true : false;
            var smoothness = (volSpec.smoothness === undefined) ? 1 : volSpec.smoothness;

            var nX = data.size.x;
            var nY = data.size.y;
            var nZ = data.size.z;
            var vertnums = new Int16Array(nX * nY * nZ);
            var vals = data.data;
            var i, il;

            for (i = 0, il = vertnums.length; i < il; ++i)
                vertnums[i] = -1;

            //mark locations partitioned by isoval
            var bitdata = new Uint8Array(nX * nY * nZ);

            for (i = 0, il = vals.length; i < il; ++i) {
                var val = (isoval >= 0) ? vals[i] - isoval : isoval - vals[i];

                if (val > 0)
                    bitdata[i] |= ISDONE;

            }

            var verts = [], faces = [];

            $3Dmol.MarchingCube.march(bitdata, verts, faces, {
                fulltable : true,
                voxel : voxel,
                unitCube : data.unit,
                origin : data.origin,
                matrix: data.matrix,
                nX : nX,
                nY : nY,
                nZ : nZ
            });

            if (!voxel && smoothness > 0)
                $3Dmol.MarchingCube.laplacianSmooth(smoothness, verts, faces);

            drawCustom(this, geo, {
                vertexArr : verts,
                faceArr : faces,
                normalArr : [],
                clickable : volSpec.clickable
            });
           
            this.updateStyle(volSpec);
            
            //computing bounding sphere from vertices
            var origin = new $3Dmol.Vector3(data.origin.x, data.origin.y, data.origin.z);
            var size = new $3Dmol.Vector3(data.size.x*data.unit.x, data.size.y*data.unit.y, data.size.z*data.unit.z);            

            var total = new $3Dmol.Vector3(0,0,0);
            var maxv = origin.clone();
            var minv = origin.clone().add(size);
            for(var i = 0; i < verts.length; i++) {
                total.add(verts[i]);
                maxv.max(verts[i]);
                minv.min(verts[i]);
            }
            total.divideScalar(verts.length);
            var len1 = total.distanceTo(minv);
            var len2 = total.distanceTo(maxv);
            this.boundingSphere.center = total;
            this.boundingSphere.radius = Math.max(len1,len2);
           
        };
        
        /** 
         * @deprecated Use addIsosurface instead
         * Creates custom shape from volumetric data 
         * @param {string} data - Volumetric input data 
         * @param {string} fmt - Input data format (e.g. 'cube' for cube file format)
         * @param {IsoSurfaceSpec} isoSpec - Volumetric data shape specification
         * @return {$3Dmol.GLShape}
         */
        this.addVolumetricData = function(data, fmt, volSpec) {
            var data = new $3Dmol.VolumeData(data, fmt);
            this.addIsosurface(data, volSpec);
        };

        /**
         * Initialize webgl objects for rendering
         * @param {$3Dmol.Object3D} group
         * 
         */  
        this.globj = function(group) {

            if (renderedShapeObj) {
                group.remove(renderedShapeObj);
                renderedShapeObj = null;
            }
            
            if(this.hidden)
                return;
            finalizeGeo(geo);
            geo.initTypedArrays();

            if(typeof(this.color) != 'undefined')
                updateColor(geo, this.color);

            shapeObj = new $3Dmol.Object3D();
            var material = null;
            if(this.side == $3Dmol.DoubleSide) {
                var material = new $3Dmol.MeshDoubleLambertMaterial({
                    wireframe : this.wireframe,
                    side : this.side,
                    transparent : (this.opacity < 1) ? true : false,
                    opacity : this.opacity,
                    wireframeLinewidth: this.linewidth
                });
            } else {
                var material = new $3Dmol.MeshLambertMaterial({
                    wireframe : this.wireframe,
                    side : this.side,
                    transparent : (this.opacity < 1) ? true : false,
                    opacity : this.opacity,
                    wireframeLinewidth: this.linewidth
                });
            }
            
            var mesh = new $3Dmol.Mesh(geo, material);

            shapeObj.add(mesh);
            
            var lineMaterial = new $3Dmol.LineBasicMaterial({
                linewidth : this.linewidth,
                color: this.color
            });
            var line = new $3Dmol.Line(linegeo, lineMaterial,
                    $3Dmol.LinePieces);
            shapeObj.add(line);

            renderedShapeObj = shapeObj.clone();
            group.add(renderedShapeObj);

        };

        this.removegl = function(group) {
            if (renderedShapeObj) {
                // dispose of geos and materials
                if (renderedShapeObj.geometry !== undefined)
                    renderedShapeObj.geometry.dispose();
                if (renderedShapeObj.material !== undefined)
                    renderedShapeObj.material.dispose();
                group.remove(renderedShapeObj);
                renderedShapeObj = null;
            }
            shapeObj = null;
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

$3Dmol.ShapeIDCount = 0;


$3Dmol.splitMesh = function(mesh) {
	    var MAXVERT = 64000; //webgl only supports 2^16 elements, leave a little breathing room (require at least 2)
    //peel off 64k vertices rsvh into their own mesh
    //duplicating vertices and normals as necessary to preserve faces and lines
	
        if(mesh.vertexArr.length < MAXVERT) return [mesh]; //typical case
        
        var nverts = mesh.vertexArr.length;
        var slices = [{vertexArr: [], normalArr: [], faceArr: []}];
        if(mesh.colorArr) slices.colorArr = [];
        var vertSlice = []; //indexed by original vertex to get current slice
        var vertIndex =[]; //indexed by original vertex to get index within slice
        var currentSlice = 0;
        
        //for each face, make sure all three vertices (or copies) are in the same slice
        var faces = mesh.faceArr;
        var vs = [0,0,0];
        for(var i = 0, nf = faces.length; i < nf; i += 3) {
            var slice = slices[currentSlice];
            for(var j = 0; j < 3; j++) {
                //process each vertex to make sure it is assigned a slice
                //all vertices of a face must belong to the same slice
                var v = faces[i+j];
                if(vertSlice[v] !== currentSlice) { //true if undefined
                    vertSlice[v] = currentSlice;
                    vertIndex[v] = slice.vertexArr.length;
                    slice.vertexArr.push(mesh.vertexArr[v]);
                    if(mesh.normalArr && mesh.normalArr[v]) slice.normalArr.push(mesh.normalArr[v]);
                    if(mesh.colorArr && mesh.colorArr[v]) slice.colorArr.push(mesh.colorArr[v]);
                }
                slice.faceArr.push(vertIndex[v]);
            }
            
            if(slice.vertexArr.length >= MAXVERT) {
                //new slice
                slices.push({vertexArr: [], normalArr: [], faceArr: []});
                if(mesh.colorArr) slices.colorArr = [];
                currentSlice++;
            }
        }
        return slices;
    }
