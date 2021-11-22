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
        if (geoGroup.vertices > 0) {
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

        color = color || $3Dmol.CC.color(color);
        geo.colorsNeedUpdate = true;

        var r, g, b;
        if (color.constructor !== Array) {
            r = color.r;
            g = color.g;
            b = color.b;
        }


        for (var gg in geo.geometryGroups) {

            var geoGroup = geo.geometryGroups[gg];
            var colorArr = geoGroup.colorArray;

            for (var i = 0, il = geoGroup.vertices; i < il; ++i) {

                if (color.constructor === Array) {
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

        var from = spec.start, end = spec.end, radius = spec.radius,
            radiusRatio = spec.radiusRatio, mid = spec.mid, midoffset = spec.midpos;

        if (!(from && end))
            return;

        var geoGroup = geo.updateGeoGroup(51);

        // vertices

        var dir = end.clone().sub(from);
        if (midoffset) { //absolute offset, convert to relative
            let length = dir.length();
            if (midoffset > 0) mid = midoffset / length;
            else mid = (length + midoffset) / length;
        }

        dir.multiplyScalar(mid);

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

            norm = [nvecs[i], nvecs[i], nvecs[i + 1], nvecs[i + 1]];

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

        face = [start + 45, start + 46, start + 1, start, start + 47,
        start + 2];
        norm = [nvecs[15], nvecs[15], nvecs[0], nvecs[0]];

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

    // Update a bounding sphere's position and radius
    // from list of centroids and new points
    /**
     * @param {$3Dmol.Sphere}
     *            sphere
     * @param {Object}
     *            components, centroid of all objects in shape
     * @param {Array}
     *            points, flat array of all points in shape
     * @param {int} numPoints, number of valid poitns in points
     */
    var updateBoundingFromPoints = function(sphere, components, points, numPoints) {

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
        if (points.length / 3 < numPoints)
            numPoints = points.length / 3;

        for (i = 0, il = numPoints; i < il; i++) {
            var x = points[i * 3], y = points[i * 3 + 1], z = points[i * 3 + 2];
            var radiusSq = sphere.center.distanceToSquared({
                x: x,
                y: y,
                z: z
            });
            maxRadiusSq = Math.max(maxRadiusSq, radiusSq);
        }

        sphere.radius = Math.sqrt(maxRadiusSq);

    };

    //helper function for adding an appropriately sized mesh
    var addCustomGeo = function(shape, geo, mesh, color, clickable) {
        var geoGroup = geo.addGeoGroup();
        var vertexArr = mesh.vertexArr, normalArr = mesh.normalArr,
            faceArr = mesh.faceArr;

        geoGroup.vertices = vertexArr.length;
        geoGroup.faceidx = faceArr.length;

        var offset, v, a, b, c, i, il, r, g;
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;

        if (color.constructor !== Array) {
            r = color.r;
            g = color.g;
            b = color.b;
        }
        for (i = 0, il = geoGroup.vertices; i < il; ++i) {
            offset = i * 3;
            v = vertexArr[i];
            vertexArray[offset] = v.x;
            vertexArray[offset + 1] = v.y;
            vertexArray[offset + 2] = v.z;

            if (color.constructor === Array) {
                c = color[i];
                r = c.r;
                g = c.g;
                b = c.b;
            }

            colorArray[offset] = r;
            colorArray[offset + 1] = g;
            colorArray[offset + 2] = b;
        }

        if (clickable) {
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

        if (clickable) {

            var center = new $3Dmol.Vector3(0, 0, 0);
            var cnt = 0;
            for (let g = 0; g < geo.geometryGroups.length; g++) {
                center.add(geo.geometryGroups[g].getCentroid());
                cnt++;
            }
            center.divideScalar(cnt);


            updateBoundingFromPoints(shape.boundingSphere, { centroid: center }, vertexArray, geoGroup.vertices);
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
        var vertexArr = mesh.vertexArr;
        var faceArr = mesh.faceArr;
        if (vertexArr.length === 0 || faceArr.length === 0) {
            console
                .warn("Error adding custom shape component: No vertices and/or face indices supplied!");
        }

        var color = customSpec.color;
        if (typeof (color) == 'undefined') {
            color = shape.color;
        }
        color = $3Dmol.CC.color(color);

        //var firstgeo = geo.geometryGroups.length;
        var splits = $3Dmol.splitMesh(mesh);
        for (var i = 0, n = splits.length; i < n; i++) {
            addCustomGeo(shape, geo, splits[i], splits[i].colorArr ? splits[i].colorArr : color, customSpec.clickable);
        }
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
        if (typeof (stylespec.color) != 'undefined') {
            shape.color = stylespec.color || new $3Dmol.Color();
            if (!(stylespec.color instanceof $3Dmol.Color))
                shape.color = $3Dmol.CC.color(stylespec.color);
        } else {
            shape.color = $3Dmol.CC.color(0);
        }
        shape.wireframe = stylespec.wireframe ? true : false;
        //opacity is the preferred nomenclature, support alpha for backwards compat
        shape.opacity = stylespec.alpha ? $3Dmol.Math.clamp(stylespec.alpha, 0.0,
            1.0) : 1.0;
        if (typeof (stylespec.opacity) != 'undefined') {
            shape.opacity = $3Dmol.Math.clamp(stylespec.opacity, 0.0, 1.0);
        }
        shape.side = (stylespec.side !== undefined) ? stylespec.side
            : $3Dmol.DoubleSide;

        shape.linewidth = typeof (stylespec.linewidth) == 'undefined' ? 1 : stylespec.linewidth;
        // Click handling
        shape.clickable = stylespec.clickable ? true : false;
        shape.callback = $3Dmol.makeFunction(stylespec.callback);
        shape.hoverable = stylespec.hoverable ? true : false;
        shape.hover_callback = $3Dmol.makeFunction(stylespec.hover_callback);
        shape.unhover_callback = $3Dmol.makeFunction(stylespec.unhover_callback);

        shape.hidden = stylespec.hidden;
        shape.frame = stylespec.frame;
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
            sphere: [],
            cylinder: [],
            line: [],
            triangle: []
        };

        updateFromStyle(this, stylespec);

        // Keep track of shape components and their centroids
        var components = [];
        var shapeObj = null;
        var renderedShapeObj = null;

        var geo = new $3Dmol.Geometry(true);
        var linegeo = new $3Dmol.Geometry(true);

        /** Update shape with new style specification
         * @function $3Dmol.GLShape#updateStyle
         * @param {ShapeSpec} newspec
         * @return {$3Dmol.GLShape}
         */
        this.updateStyle = function(newspec) {

            for (var prop in newspec) {
                stylespec[prop] = newspec[prop];
            }

            updateFromStyle(this, stylespec);
            
            if (newspec.voldata && newspec.volscheme) {
                $3Dmol.adjustVolumeStyle(newspec);
    
                //convert volumetric data into colors
                const scheme = newspec.volscheme;
                const voldata = newspec.voldata;
                const cc = $3Dmol.CC;
                const range = scheme.range() || [-1, 1];
                geo.setColors(function(x, y, z) {
                    let val = voldata.getVal(x, y, z);
                    let col = cc.color(scheme.valueToHex(val, range));
                    return col;
                });
                delete this.color;
            }    
        };

        /**
         * Creates a custom shape from supplied vertex and face arrays
         * @function $3Dmol.GLShape#addCustom
         * @param {CustomShapeSpec} customSpec
         * @return {$3Dmol.GLShape}
         
         */
        this.addCustom = function(customSpec) {

            customSpec.vertexArr = customSpec.vertexArr || [];
            customSpec.faceArr = customSpec.faceArr || [];
            customSpec.normalArr = customSpec.normalArr || [];

            // will split mesh as needed
            drawCustom(this, geo, customSpec);
        };

        /**
         * Creates a sphere shape
         * @function $3Dmol.GLShape#addSphere
         * @param {SphereSpec} sphereSpec
         * @return {$3Dmol.GLShape}
         @example 
         viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
         
         viewer.render();
         */
        this.addSphere = function(sphereSpec) {

            sphereSpec.center = sphereSpec.center || {
                x: 0,
                y: 0,
                z: 0
            };
            sphereSpec.radius = sphereSpec.radius ? $3Dmol.Math.clamp(
                sphereSpec.radius, 0, Infinity) : 1.5;
            sphereSpec.color = $3Dmol.CC.color(sphereSpec.color);

            this.intersectionShape.sphere.push(new $3Dmol.Sphere(
                sphereSpec.center, sphereSpec.radius));

            $3Dmol.GLDraw.drawSphere(geo, sphereSpec.center,
                sphereSpec.radius, sphereSpec.color);

            components.push({
                centroid: new $3Dmol.Vector3(sphereSpec.center.x,
                    sphereSpec.center.y, sphereSpec.center.z)
            });
            var geoGroup = geo.updateGeoGroup(0);

            updateBoundingFromPoints(this.boundingSphere, components,
                geoGroup.vertexArray, geoGroup.vertices);
        };


        /**
         * Creates a box
         * @function $3Dmol.GLShape#addBox
         * @param {BoxSpec} boxSpec
         * @return {$3Dmol.GLShape}
         @example 
         var shape = viewer.addShape({color:'red'});
         shape.addBox({corner: {x:1,y:2,z:0}, dimensions: {w: 4, h: 2, d: 6}});
         shape.addBox({corner: {x:-5,y:-3,z:0},
                       dimensions: { w: {x:1,y:1,z:0},
                                     h: {x:-1,y:1,z:0},
                                     d: {x:0,y:0,z:1} }});
         viewer.zoomTo();
         viewer.rotate(30);
         viewer.render();
         */
        this.addBox = function(boxSpec) {

            var dim = boxSpec.dimensions || { w: 1, h: 1, d: 1 };

            //dimensions may be scalar or vector quantities
            var w = dim.w;
            if (typeof (dim.w) == "number") {
                w = { x: dim.w, y: 0, z: 0 };
            }
            var h = dim.h;
            if (typeof (dim.h) == "number") {
                h = { x: 0, y: dim.h, z: 0 };
            }
            var d = dim.d;
            if (typeof (dim.d) == "number") {
                d = { x: 0, y: 0, z: dim.d };
            }

            //can position using corner OR center
            var c = boxSpec.corner;
            if (c == undefined) {
                if (boxSpec.center !== undefined) {

                    c = {
                        x: boxSpec.center.x - 0.5 * (w.x + h.x + d.x),
                        y: boxSpec.center.y - 0.5 * (w.y + h.y + d.y),
                        z: boxSpec.center.z - 0.5 * (w.z + h.z + d.z)
                    };
                } else { // default to origin
                    c = { x: 0, y: 0, z: 0 };
                }
            }

            //8 vertices
            var uv =
                [{ x: c.x, y: c.y, z: c.z },
                { x: c.x + w.x, y: c.y + w.y, z: c.z + w.z },
                { x: c.x + h.x, y: c.y + h.y, z: c.z + h.z },
                { x: c.x + w.x + h.x, y: c.y + w.y + h.y, z: c.z + w.z + h.z },
                { x: c.x + d.x, y: c.y + d.y, z: c.z + d.z },
                { x: c.x + w.x + d.x, y: c.y + w.y + d.y, z: c.z + w.z + d.z },
                { x: c.x + h.x + d.x, y: c.y + h.y + d.y, z: c.z + h.z + d.z },
                { x: c.x + w.x + h.x + d.x, y: c.y + w.y + h.y + d.y, z: c.z + w.z + h.z + d.z }];

            //but.. so that we can have sharp issues, we want a unique normal
            //for each face - since normals are associated with vertices, need to duplicate 

            //bottom
            // 0 1
            // 2 3
            //top
            // 4 5
            // 6 7
            var verts = [];
            var faces = [];
            //bottom
            verts.splice(verts.length, 0, uv[0], uv[1], uv[2], uv[3]);
            faces.splice(faces.length, 0, 0, 2, 1, 1, 2, 3);
            var foff = 4;
            //front
            verts.splice(verts.length, 0, uv[2], uv[3], uv[6], uv[7]);
            faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
            foff += 4;
            //back
            verts.splice(verts.length, 0, uv[4], uv[5], uv[0], uv[1]);
            faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
            foff += 4;
            //top
            verts.splice(verts.length, 0, uv[6], uv[7], uv[4], uv[5]);
            faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
            foff += 4;
            //right
            verts.splice(verts.length, 0, uv[3], uv[1], uv[7], uv[5]);
            faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
            foff += 4;
            //left
            verts.splice(verts.length, 0, uv[2], uv[6], uv[0], uv[4]); // fix: was 2 0 6 4 , was flipped! will this ruin anything?
            // and is this the reason for having double sided lambert shading? the box had a flipped face
            faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
            foff += 4;

            var spec = $3Dmol.extend({}, boxSpec);
            spec.vertexArr = verts;
            spec.faceArr = faces;
            spec.normalArr = [];
            drawCustom(this, geo, spec);

            var centroid = new $3Dmol.Vector3();
            components.push({
                centroid: centroid.addVectors(uv[0], uv[7]).multiplyScalar(0.5)
            });
            var geoGroup = geo.updateGeoGroup(0);
            updateBoundingFromPoints(this.boundingSphere, components, geoGroup.vertexArray, geoGroup.vertices);
        };

        /**
         * Creates a cylinder shape
         * @function $3Dmol.GLShape#addCylinder
         * @param {CylinderSpec} cylinderSpec
         * @return {$3Dmol.GLShape}
         @example
              viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
                                  end:{x:10.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  fromCap:1,
                                  toCap:2,
                                  color:'red',
                                  hoverable:true,
                                  clickable:true,
                                  callback:function(){ this.color.setHex(0x00FFFF00);viewer.render( );},
                                  hover_callback: function(){ viewer.render( );},
                                  unhover_callback: function(){ this.color.setHex(0xFF000000);viewer.render( );}
                                 });
              viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
                                  end:{x:0.0,y:10.0,z:0.0},
                                  radius:0.5,
                                  fromCap:false,
                                  toCap:true,
                                  color:'teal'});
              viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
                                  end:{x:20.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  color:'black',
                                  fromCap:false,
                                  toCap:false});
              viewer.render();
         */
        this.addCylinder = function(cylinderSpec) {

            cylinderSpec.start = cylinderSpec.start || {};
            cylinderSpec.end = cylinderSpec.end || {};


            var start = new $3Dmol.Vector3(cylinderSpec.start.x || 0,
                cylinderSpec.start.y || 0, cylinderSpec.start.z || 0);
            var end = new $3Dmol.Vector3(cylinderSpec.end.x,
                cylinderSpec.end.y || 0, cylinderSpec.end.z || 0);
            if (typeof (end.x) == 'undefined') end.x = 3; //show something even if undefined

            var radius = cylinderSpec.radius || 0.1;
            var color = $3Dmol.CC.color(cylinderSpec.color);

            this.intersectionShape.cylinder.push(new $3Dmol.Cylinder(start, end, radius));

            $3Dmol.GLDraw.drawCylinder(geo, start, end, radius, color, cylinderSpec.fromCap, cylinderSpec.toCap);

            var centroid = new $3Dmol.Vector3();
            components.push({
                centroid: centroid.addVectors(start, end).multiplyScalar(0.5)
            });
            var geoGroup = geo.updateGeoGroup(0);
            updateBoundingFromPoints(this.boundingSphere, components,
                geoGroup.vertexArray, geoGroup.vertices);

        };

        /**
         * Creates a dashed cylinder shape
         * @function $3Dmol.GLShape#addDashedCylinder
         * @param {CylinderSpec} cylinderSpec
         * @return {$3Dmol.GLShape}       
         */
        this.addDashedCylinder = function(cylinderSpec) {
            cylinderSpec.start = cylinderSpec.start || {};
            cylinderSpec.end = cylinderSpec.end || {};
            cylinderSpec.dashLength = cylinderSpec.dashLength || 0.25;
            cylinderSpec.gapLength = cylinderSpec.gapLength || 0.25;

            var start = new $3Dmol.Vector3(cylinderSpec.start.x || 0,
                cylinderSpec.start.y || 0, cylinderSpec.start.z || 0);
            var end = new $3Dmol.Vector3(cylinderSpec.end.x,
                cylinderSpec.end.y || 0, cylinderSpec.end.z || 0);
            if (typeof (end.x) == 'undefined') end.x = 3; //show something even if undefined

            var radius = cylinderSpec.radius || 0.1;
            var color = $3Dmol.CC.color(cylinderSpec.color);

            var cylinderLength = Math.sqrt(Math.pow((start.x - end.x), 2) + Math.pow((start.y - end.y), 2) + Math.pow((start.z - end.z), 2));

            var count = cylinderLength / (cylinderSpec.gapLength + cylinderSpec.dashLength);

            var new_start = new $3Dmol.Vector3(cylinderSpec.start.x || 0,
                cylinderSpec.start.y || 0, cylinderSpec.start.z || 0);
            var new_end = new $3Dmol.Vector3(cylinderSpec.end.x,
                cylinderSpec.end.y || 0, cylinderSpec.end.z || 0);

            var gapVector = new $3Dmol.Vector3((end.x - start.x) / (cylinderLength / cylinderSpec.gapLength), (end.y - start.y) / (cylinderLength / cylinderSpec.gapLength), (end.z - start.z) / (cylinderLength / cylinderSpec.gapLength));
            var dashVector = new $3Dmol.Vector3((end.x - start.x) / (cylinderLength / cylinderSpec.dashLength), (end.y - start.y) / (cylinderLength / cylinderSpec.dashLength), (end.z - start.z) / (cylinderLength / cylinderSpec.dashLength));

            for (var place = 0; place < count; place++) {
                new_end = new $3Dmol.Vector3(new_start.x + dashVector.x, new_start.y + dashVector.y, new_start.z + dashVector.z);

                this.intersectionShape.cylinder.push(new $3Dmol.Cylinder(new_start, new_end, radius));

                $3Dmol.GLDraw.drawCylinder(geo, new_start, new_end, radius, color, cylinderSpec.fromCap, cylinderSpec.toCap);

                new_start = new $3Dmol.Vector3(new_end.x + gapVector.x, new_end.y + gapVector.y, new_end.z + gapVector.z);

            }
            var centroid = new $3Dmol.Vector3();
            components.push({
                centroid: centroid.addVectors(start, end).multiplyScalar(0.5)
            });
            var geoGroup = geo.updateGeoGroup(0);
            updateBoundingFromPoints(this.boundingSphere, components,
                geoGroup.vertexArray, geoGroup.vertices);
        };

        /**
         * Creates a curved shape
         * @function $3Dmol.GLShape#addCurve
         * @param {CurveSpec} curveSpec
         * @return {$3Dmol.GLShape}
         */
        this.addCurve = function(curveSpec) {

            curveSpec.points = curveSpec.points || [];
            curveSpec.smooth = curveSpec.smooth || 10;
            if (typeof (curveSpec.fromCap) == "undefined") curveSpec.fromCap = 2;
            if (typeof (curveSpec.toCap) == "undefined") curveSpec.toCap = 2;

            //subdivide into smoothed spline points
            var points = $3Dmol.subdivide_spline(curveSpec.points, curveSpec.smooth);

            if (points.length < 3) {
                console.log("Too few points in addCurve");
                return;
            }

            var radius = curveSpec.radius || 0.1;
            var color = $3Dmol.CC.color(curveSpec.color);
            //TODO TODO - this is very inefficient, should create our
            //own water tight model with proper normals...


            //if arrows are requested, peel off enough points to fit
            //at least 2*r of arrowness
            var start = 0;
            var end = points.length - 1;
            var segmentlen = points[0].distanceTo(points[1]);
            var npts = Math.ceil(2 * radius / segmentlen);
            if (curveSpec.toArrow) {
                end -= npts;
                let arrowspec = {
                    start: points[end],
                    end: points[points.length - 1],
                    radius: radius,
                    color: color,
                    mid: 0.0001
                };
                this.addArrow(arrowspec);
            }
            if (curveSpec.fromArrow) {
                start += npts;
                let arrowspec = {
                    start: points[start],
                    end: points[0],
                    radius: radius,
                    color: color,
                    mid: 0.0001
                };
                this.addArrow(arrowspec);
            }

            var midway = Math.ceil(points.length / 2);
            var middleSpec = { radius: radius, color: color, fromCap: 2, toCap: 2 };
            for (var i = start; i < end; i++) {
                middleSpec.start = points[i];
                middleSpec.end = points[i + 1];
                middleSpec.fromCap = 2;
                middleSpec.toCap = 2;
                if (i < midway) {
                    middleSpec.fromCap = 2;
                    middleSpec.toCap = 0;
                } else if (i > midway) {
                    middleSpec.fromCap = 0;
                    middleSpec.toCap = 2;
                } else {
                    middleSpec.fromCap = 2;
                    middleSpec.toCap = 2;
                }

                this.addCylinder(middleSpec);
            }


        };

        //*****BFM Add cone, copied from cylinder
        this.addCone = function(coneSpec) {
            coneSpec.start = coneSpec.start || {};
            coneSpec.end = coneSpec.end || {};

            var start = new $3Dmol.Vector3(coneSpec.start.x || 0,
                        coneSpec.start.y || 0, coneSpec.start.z || 0);
            var end = new $3Dmol.Vector3(coneSpec.end.x,
                        coneSpec.end.y || 0, coneSpec.end.z || 0);
            if(typeof(end.x) == 'undefined') end.x = 3; //show something even if undefined

            var radius = coneSpec.radius || 0.1;
            var color = $3Dmol.CC.color(coneSpec.color);

            // No intersection support
            // this.intersectionShape.cylinder.push(new $3Dmol.Cylinder(start, end, radius));

            $3Dmol.GLDraw.drawCone(geo, start, end, radius, color);

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
         @example
         $3Dmol.download("pdb:2ABJ",viewer,{},function(){
                  viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
                  viewer.render(callback);
              });

         */
        this.addLine = function(lineSpec) {
            lineSpec.start = lineSpec.start || {};
            lineSpec.end = lineSpec.end || {};

            var start = new $3Dmol.Vector3(lineSpec.start.x || 0,
                lineSpec.start.y || 0, lineSpec.start.z || 0);
            var end = new $3Dmol.Vector3(lineSpec.end.x,
                lineSpec.end.y || 0, lineSpec.end.z || 0);
            if (typeof (end.x) == 'undefined') end.x = 3; //show something even if undefined

            var geoGroup = geo.updateGeoGroup(2);

            //make line from start to end
            //for consistency with rest of shapes, uses vertices and lines rather
            //than a separate line geometry
            var vstart = geoGroup.vertices;
            var i = vstart * 3;
            var vertexArray = geoGroup.vertexArray;
            vertexArray[i] = start.x;
            vertexArray[i + 1] = start.y;
            vertexArray[i + 2] = start.z;
            vertexArray[i + 3] = end.x;
            vertexArray[i + 4] = end.y;
            vertexArray[i + 5] = end.z;
            geoGroup.vertices += 2;

            var lineArray = geoGroup.lineArray;
            var li = geoGroup.lineidx;
            lineArray[li] = vstart;
            lineArray[li + 1] = vstart + 1;
            geoGroup.lineidx += 2;

            var centroid = new $3Dmol.Vector3();
            components.push({
                centroid: centroid.addVectors(start, end).multiplyScalar(0.5)
            });
            geoGroup = geo.updateGeoGroup(0);
            updateBoundingFromPoints(this.boundingSphere, components,
                geoGroup.vertexArray, geoGroup.vertices);
        };

        /**
         * Creates an arrow shape
         * @function $3Dmol.GLShape#addArrow        
         * @param {ArrowSpec} arrowSpec
         * @return {$3Dmol.GLShape}
         @example
          $3Dmol.download("pdb:4DM7",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  viewer.addArrow({
                      start: {x:-10.0, y:0.0, z:0.0},
                      end: {x:0.0, y:-10.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      clickable:true,
                      callback:function(){
                          this.color.setHex(0xFF0000FF);
                          viewer.render( );
                      }
                  });
                  viewer.render();
                });
         */
        this.addArrow = function(arrowSpec) {

            arrowSpec.start = arrowSpec.start || {};
            arrowSpec.end = arrowSpec.end || {};

            arrowSpec.start = new $3Dmol.Vector3(arrowSpec.start.x || 0,
                arrowSpec.start.y || 0, arrowSpec.start.z || 0);

            if (arrowSpec.dir instanceof $3Dmol.Vector3 && typeof (arrowSpec.length) === 'number') {
                var end = arrowSpec.dir.clone().multiplyScalar(arrowSpec.length).add(
                    arrowSpec.start);
                arrowSpec.end = end;
            }

            else {
                arrowSpec.end = new $3Dmol.Vector3(arrowSpec.end.x,
                    arrowSpec.end.y || 0, arrowSpec.end.z || 0);
                if (typeof (arrowSpec.end.x) == 'undefined') arrowSpec.end.x = 3; //show something even if undefined
            }

            arrowSpec.radius = arrowSpec.radius || 0.1;

            arrowSpec.radiusRatio = arrowSpec.radiusRatio || 1.618034;

            arrowSpec.mid = (0 < arrowSpec.mid && arrowSpec.mid < 1) ? arrowSpec.mid
                : 0.618034;

            drawArrow(this, geo, arrowSpec);

            var centroid = new $3Dmol.Vector3();
            components.push({
                centroid: centroid.addVectors(arrowSpec.start, arrowSpec.end)
                    .multiplyScalar(0.5)
            });
            var geoGroup = geo.updateGeoGroup(0);
            updateBoundingFromPoints(this.boundingSphere, components,
                geoGroup.vertexArray, geoGroup.vertices);

        };


        var distance_from = function(c1, c2) {
            return Math.sqrt(Math.pow((c1.x - c2.x), 2) + Math.pow((c1.y - c2.y), 2) + Math.pow((c1.z - c2.z), 2));
        };

        var inSelectedRegion = function(coordinate, selectedRegion, radius) {

            for (var i = 0; i < selectedRegion.length; i++) {
                if (distance_from(selectedRegion[i], coordinate) <= radius)
                    return true;
            }
            return false;
        };

        /**
         * Create isosurface from voluemetric data.
         * @function $3Dmol.GLShape#addIsosurface         
         * @param {$3Dmol.VolumeData} data - volumetric input data
         * @param {IsoSurfaceSpec} isoSpec - volumetric data shape specification
         * @example //the user can specify a selected region for the isosurface 
         $.get('../test_structs/benzene-homo.cube', function(data){
                  var voldata = new $3Dmol.VolumeData(data, "cube");
                  viewer.addIsosurface(voldata, {isoval: 0.01,
                                                 color: "blue",
                                                 alpha: 0.5,
                                                 smoothness: 10});
                  viewer.addIsosurface(voldata, {isoval: -0.01,
                                                 color: "red",
                                                 smoothness: 5,
                                                 opacity:0.5,
                                                 wireframe:true,
                                                 clickable:true,
                                                 callback:
                                                 function() {
                                                     this.opacity = 0.0;
                                                     viewer.render( );
                                                 }});
                  viewer.setStyle({}, {stick:{}});
                  viewer.zoomTo();
                  viewer.render();
                });
         */
        this.addIsosurface = function(data, volSpec, callback) {//may want to cache the arrays geneerated when selectedRegion ==true

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

            var bitdata = new Uint8Array(nX * nY * nZ);

            //mark locations partitioned by isoval
            for (i = 0, il = vals.length; i < il; ++i) {
                var val = (isoval >= 0) ? vals[i] - isoval : isoval - vals[i];
                if (val > 0)
                    bitdata[i] |= ISDONE;

            }

            var verts = [], faces = [];

            $3Dmol.MarchingCube.march(bitdata, verts, faces, {
                fulltable: true,
                voxel: voxel,
                unitCube: data.unit,
                origin: data.origin,
                matrix: data.matrix,
                nX: nX,
                nY: nY,
                nZ: nZ
            });

            if (!voxel && smoothness > 0)
                $3Dmol.MarchingCube.laplacianSmooth(smoothness, verts, faces);
            var vertexmapping = [];
            var newvertices = [];
            var newfaces = [];

            if (volSpec.selectedRegion && volSpec.coords === undefined) {
                volSpec.coords = volSpec.selectedRegion; //backwards compat for incorrectly documented feature
            }
            if (volSpec.coords !== undefined) {

                var xmax = volSpec.coords[0].x,
                    ymax = volSpec.coords[0].y,
                    zmax = volSpec.coords[0].z,
                    xmin = volSpec.coords[0].x,
                    ymin = volSpec.coords[0].y,
                    zmin = volSpec.coords[0].z;

                for (let i = 0; i < volSpec.coords.length; i++) {
                    if (volSpec.coords[i].x > xmax)
                        xmax = volSpec.coords[i].x;
                    else if (volSpec.coords[i].x < xmin)
                        xmin = volSpec.coords[i].x;
                    if (volSpec.coords[i].y > ymax)
                        ymax = volSpec.coords[i].y;
                    else if (volSpec.coords[i].y < ymin)
                        ymin = volSpec.coords[i].y;
                    if (volSpec.coords[i].z > zmax)
                        zmax = volSpec.coords[i].z;
                    else if (volSpec.coords[i].z < zmin)
                        zmin = volSpec.coords[i].z;
                }

                var rad = 2;
                if (volSpec.radius !== undefined) {
                    rad = volSpec.radius; //backwards compat
                }
                if (volSpec.selectedOffset !== undefined) { //backwards compat
                    rad = volSpec.selectedOffset;
                }
                if (volSpec.seldist !== undefined) {
                    rad = volSpec.seldist;
                }

                xmin -= rad;
                xmax += rad;
                ymin -= rad;
                ymax += rad;
                zmin -= rad;
                zmax += rad;

                // accounts for radius
                for (let i = 0; i < verts.length; i++) {
                    if (verts[i].x > xmin &&
                        verts[i].x < xmax &&
                        verts[i].y > ymin &&
                        verts[i].y < ymax &&
                        verts[i].z > zmin &&
                        verts[i].z < zmax &&
                        inSelectedRegion(verts[i],
                            volSpec.coords, rad)) {
                        vertexmapping.push(newvertices.length);
                        newvertices.push(verts[i]);

                    } else {
                        vertexmapping.push(-1);
                    }

                }
                for (let i = 0; i + 2 < faces.length; i += 3) {
                    if (vertexmapping[faces[i]] !== -1 &&
                        vertexmapping[faces[i + 1]] !== -1 &&
                        vertexmapping[faces[i + 2]] !== -1) {
                        newfaces.push(faces[i] - (faces[i] - vertexmapping[faces[i]]));
                        newfaces.push(faces[i + 1] - (faces[i + 1] - vertexmapping[faces[i + 1]]));
                        newfaces.push(faces[i + 2] - (faces[i + 2] - vertexmapping[faces[i + 2]]));
                    }
                }
                verts = newvertices;
                faces = newfaces;
            }

            drawCustom(this, geo, {
                vertexArr: verts,
                faceArr: faces,
                normalArr: [],
                clickable: volSpec.clickable,
                hoverable: volSpec.hoverable
            });

            this.updateStyle(volSpec);

            //computing bounding sphere from vertices
            var origin = new $3Dmol.Vector3(data.origin.x, data.origin.y, data.origin.z);
            var size = new $3Dmol.Vector3(data.size.x * data.unit.x, data.size.y * data.unit.y, data.size.z * data.unit.z);

            var total = new $3Dmol.Vector3(0, 0, 0);
            var maxv = origin.clone();
            var minv = origin.clone().add(size);
            for (let i = 0; i < verts.length; i++) {
                total.add(verts[i]);
                maxv.max(verts[i]);
                minv.min(verts[i]);
            }
            total.divideScalar(verts.length);
            var len1 = total.distanceTo(minv);
            var len2 = total.distanceTo(maxv);
            this.boundingSphere.center = total;
            this.boundingSphere.radius = Math.max(len1, len2);
            if (typeof callback == "function")
                callback();
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
            data = new $3Dmol.VolumeData(data, fmt);
            this.addIsosurface(data, volSpec);
        };

        //for internal use, truncate buffers to save memory
        this.finalize = function() {
            finalizeGeo(geo);
            geo.initTypedArrays();
            return geo;
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

            if (this.hidden)
                return;
            finalizeGeo(geo);
            geo.initTypedArrays();

            if (this.wireframe) {
                geo.setUpWireframe();
            }

            if (typeof (this.color) != 'undefined')
                updateColor(geo, this.color);

            shapeObj = new $3Dmol.Object3D();
            var material = null;

            if (this.side == $3Dmol.DoubleSide) {
                material = new $3Dmol.MeshDoubleLambertMaterial({
                    wireframe: this.wireframe,
                    side: this.side,
                    transparent: (this.opacity < 1) ? true : false,
                    opacity: this.opacity,
                    wireframeLinewidth: this.linewidth,
                    vertexColors: $3Dmol.VertexColors
                });
            } else {
                material = new $3Dmol.MeshLambertMaterial({
                    wireframe: this.wireframe,
                    side: this.side,
                    transparent: (this.opacity < 1) ? true : false,
                    opacity: this.opacity,
                    wireframeLinewidth: this.linewidth,
                    vertexColors: $3Dmol.VertexColors
                });
            }

            var mesh = new $3Dmol.Mesh(geo, material);

            shapeObj.add(mesh);

            var lineMaterial = new $3Dmol.LineBasicMaterial({
                linewidth: this.linewidth,
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

    }

    Object.defineProperty(GLShape.prototype, "position", {

        get: function() {
            return this.boundingSphere.center;
        }

    });

    Object.defineProperty(GLShape.prototype, "x", {

        get: function() {
            return this.boundingSphere.center.x;
        }

    });

    Object.defineProperty(GLShape.prototype, "y", {

        get: function() {
            return this.boundingSphere.center.y;
        }

    });

    Object.defineProperty(GLShape.prototype, "z", {

        get: function() {
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

    if (mesh.vertexArr.length < MAXVERT) return [mesh]; //typical case

    var slices = [{ vertexArr: [], normalArr: [], faceArr: [] }];
    if (mesh.colorArr) slices.colorArr = [];
    var vertSlice = []; //indexed by original vertex to get current slice
    var vertIndex = []; //indexed by original vertex to get index within slice
    var currentSlice = 0;

    //for each face, make sure all three vertices (or copies) are in the same slice
    var faces = mesh.faceArr;
    for (let i = 0, nf = faces.length; i < nf; i += 3) {
        let slice = slices[currentSlice];
        for (let j = 0; j < 3; j++) {
            //process each vertex to make sure it is assigned a slice
            //all vertices of a face must belong to the same slice
            var v = faces[i + j];
            if (vertSlice[v] !== currentSlice) { //true if undefined
                vertSlice[v] = currentSlice;
                vertIndex[v] = slice.vertexArr.length;
                slice.vertexArr.push(mesh.vertexArr[v]);
                if (mesh.normalArr && mesh.normalArr[v]) slice.normalArr.push(mesh.normalArr[v]);
                if (mesh.colorArr && mesh.colorArr[v]) slice.colorArr.push(mesh.colorArr[v]);
            }
            slice.faceArr.push(vertIndex[v]);
        }

        if (slice.vertexArr.length >= MAXVERT) {
            //new slice
            slices.push({ vertexArr: [], normalArr: [], faceArr: [] });
            if (mesh.colorArr) slices.colorArr = [];
            currentSlice++;
        }
    }
    return slices;
};

