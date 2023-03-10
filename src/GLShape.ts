import { Geometry, Material } from "./WebGL";
import { Sphere, Cylinder, Triangle } from "./WebGL/shapes";
import { Vector3, XYZ } from "./WebGL/math";
import { clamp } from "./WebGL/math";
import { DoubleSide } from "./WebGL";
import { Color, CC, ColorSpec } from "./colors";
import { MarchingCube } from "./ProteinSurface4";
import { VolumeData } from "./VolumeData";
import { MeshDoubleLambertMaterial, MeshLambertMaterial, Object3D, Coloring, Mesh, LineBasicMaterial, Line, LineStyle } from "./WebGL";
import { CAP, GLDraw } from "./GLDraw"
import { subdivide_spline } from "./glcartoon";
import { adjustVolumeStyle, extend, Func, makeFunction } from "./utilities";
import { Gradient, GradientType } from "./Gradient";
import { LineStyleSpec } from "GLModel";
import { VolumetricRendererSpec } from "VolumetricRender";


/**
 * A GLShape is a collection of user specified shapes.
 * 
 * @class
 * @extends {ShapeSpec}
 * @param {number} sid - Unique identifier
 * @param {ShapeSpec} stylespec - shape style specification
 */
export class GLShape {

    // Marching cube, to match with protein surface generation
    private static ISDONE = 2;

    private static finalizeGeo(geo) {
        //to avoid creating a bunch of geometries, we leave geoGroup untruncated
        //until render is called, at which point we truncate; 
        //successive called up updateGeo will return a new geometry
        var geoGroup = geo.updateGeoGroup(0);
        if (geoGroup.vertices > 0) {
            geoGroup.truncateArrayBuffers(true, true);
        }
    };

    /* 
     * 
     * @param {Geometry}
     *            geo
     * @param {Color | colorlike} color
     */
    static updateColor(geo: Geometry, color) {

        color = color || CC.color(color);
        geo.colorsNeedUpdate = true;

        var r, g, b;
        if (color.constructor !== Array) {
            r = color.r;
            g = color.g;
            b = color.b;
        }


        for (let gg in geo.geometryGroups) {

            let geoGroup = geo.geometryGroups[gg];
            let colorArr = geoGroup.colorArray;

            for (let i = 0, il = geoGroup.vertices; i < il; ++i) {

                if (color.constructor === Array) {
                    let c = color[i];
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


    /*
     * @param {GLShape}
     *            shape
     * @param {geometryGroup}
     *            geoGroup
     * @param {ArrowSpec}
     *            spec
     */
    static drawArrow(shape: GLShape, geo: Geometry, spec: ArrowSpec) {

        var from = spec.start, end = spec.end, radius = spec.radius,
            radiusRatio = spec.radiusRatio, mid = spec.mid, midoffset = spec.midpos;

        if (!(from && end))
            return;

        var geoGroup = geo.updateGeoGroup(51);

        // vertices

        var dir = new Vector3(end.x, end.y, end.z).sub(from);
        if (midoffset) { //absolute offset, convert to relative
            let length = dir.length();
            if (midoffset > 0) mid = midoffset / length;
            else mid = (length + midoffset) / length;
        }

        dir.multiplyScalar(mid);

        var to = new Vector3(from.x, from.y, from.z).add(dir);
        var negDir = dir.clone().negate();

        let fromv = new Vector3(from.x, from.y, from.z);
        shape.intersectionShape.cylinder.push(new Cylinder(fromv, to.clone(), radius));
        shape.intersectionShape.sphere.push(new Sphere(fromv, radius));

        // get orthonormal vector
        var nvecs = [];
        nvecs[0] = dir.clone();
        if (Math.abs(nvecs[0].x) > 0.0001)
            nvecs[0].y += 1;
        else
            nvecs[0].x += 1;
        nvecs[0].cross(dir);
        nvecs[0].normalize();

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

                var c = new Vector3(prev_x, prev_y, prev_z);
                var b = new Vector3(end.x, end.y, end.z), b2 = to.clone();
                var a = new Vector3(conebase.x, conebase.y, conebase.z);

                shape.intersectionShape.triangle.push(new Triangle(a, b, c));
                shape.intersectionShape.triangle.push(new Triangle(c.clone(), b2, a.clone()));
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
        var face, faceoffset, lineoffset;
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

            //   norm = [nvecs[i], nvecs[i], nvecs[i + 1], nvecs[i + 1]];

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
        //norm = [nvecs[15], nvecs[15], nvecs[0], nvecs[0]];

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
    /*
     * @param {Sphere}
     *            sphere
     * @param {Object}
     *            components, centroid of all objects in shape
     * @param {Array}
     *            points, flat array of all points in shape
     * @param {int} numPoints, number of valid poitns in points
     */
    static updateBoundingFromPoints(sphere: Sphere, components, points, numPoints: number) {

        sphere.center.set(0, 0, 0);

        //previously I weighted each component's center equally, but I think
        //it is better to use all points
        let xmin = Infinity, ymin = Infinity, zmin = Infinity;
        let xmax = -Infinity, ymax = -Infinity, zmax = -Infinity;
        if (sphere.box) {
            xmin = sphere.box.min.x;
            xmax = sphere.box.max.x;
            ymin = sphere.box.min.y;
            ymax = sphere.box.max.y;
            zmin = sphere.box.min.z;
            zmax = sphere.box.max.z;
        }

        for (let i = 0, il = numPoints; i < il; i++) {
            var x = points[i * 3], y = points[i * 3 + 1], z = points[i * 3 + 2];
            if (x < xmin) xmin = x;
            if (y < ymin) ymin = y;
            if (z < zmin) zmin = z;
            if (x > xmax) xmax = x;
            if (y > ymax) ymax = y;
            if (z > zmax) zmax = z;
        }

        sphere.center.set((xmax + xmin) / 2, (ymax + ymin) / 2, (zmax + zmin) / 2);
        sphere.radius = sphere.center.distanceTo({ x: xmax, y: ymax, z: zmax });
        sphere.box = { min: { x: xmin, y: ymin, z: zmin }, max: { x: xmax, y: ymax, z: zmax } };
    };

    //helper function for adding an appropriately sized mesh
    private static addCustomGeo(shape: GLShape, geo: Geometry, mesh, color, clickable) {
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
                var vA = new Vector3(), vB = new Vector3(), vC = new Vector3();
                shape.intersectionShape.triangle.push(new Triangle(vA.copy(vertexArr[a]),
                    vB.copy(vertexArr[b]), vC.copy(vertexArr[c])));
            }
        }

        if (clickable) {

            var center = new Vector3(0, 0, 0);
            var cnt = 0;
            for (let g = 0; g < geo.geometryGroups.length; g++) {
                center.add(geo.geometryGroups[g].getCentroid());
                cnt++;
            }
            center.divideScalar(cnt);


            GLShape.updateBoundingFromPoints(shape.boundingSphere, { centroid: center }, vertexArray, geoGroup.vertices);
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
    /*
     * @param {$3Dmol.GLShape}
     *            shape
     * @param {geometry}
     *            geo
     * @param {CustomShapeSpec}
     *            customSpec
     */
    static drawCustom = function (shape: GLShape, geo: Geometry, customSpec: CustomShapeSpec) {
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
        color = CC.color(color);

        //var firstgeo = geo.geometryGroups.length;
        var splits = splitMesh(mesh);
        for (var i = 0, n = splits.length; i < n; i++) {
            GLShape.addCustomGeo(shape, geo, splits[i], splits[i].colorArr ? splits[i].colorArr : color, customSpec.clickable);
        }
    };


    /*
     * 
     * @param {$3Dmol.GLShape}
     *            shape
     * @param {ShapeSpec}
     *            stylespec
     * @returns {undefined}
     */
    static updateFromStyle(shape: GLShape, stylespec: ShapeSpec) {
        if (typeof (stylespec.color) != 'undefined') {
            shape.color = stylespec.color || new Color();
            if (!(stylespec.color instanceof Color))
                shape.color = CC.color(stylespec.color);
        } else {
            shape.color = CC.color(0);
        }
        shape.wireframe = stylespec.wireframe ? true : false;
        //opacity is the preferred nomenclature, support alpha for backwards compat
        shape.opacity = stylespec.alpha ? clamp(stylespec.alpha, 0.0,
            1.0) : 1.0;
        if (typeof (stylespec.opacity) != 'undefined') {
            shape.opacity = clamp(stylespec.opacity, 0.0, 1.0);
        }
        shape.side = (stylespec.side !== undefined) ? stylespec.side : DoubleSide;
        shape.linewidth = typeof (stylespec.linewidth) == 'undefined' ? 1 : stylespec.linewidth;
        // Click handling
        shape.clickable = stylespec.clickable ? true : false;
        shape.callback = makeFunction(stylespec.callback);
        shape.hoverable = stylespec.hoverable ? true : false;
        shape.hover_callback = makeFunction(stylespec.hover_callback);
        shape.unhover_callback = makeFunction(stylespec.unhover_callback);
        shape.contextMenuEnabled = !!stylespec.contextMenuEnabled;

        shape.hidden = stylespec.hidden;
        shape.frame = stylespec.frame;
    };

    boundingSphere: Sphere;
    intersectionShape: any;
    color: any;
    hidden = false;
    wireframe = false;
    opacity = 1;
    linewidth = 1;
    clickable = false;
    callback: Func;
    hoverable = false;
    hover_callback: Func;
    unhover_callback: Func;
    contextMenuEnabled: boolean = false;
    frame: any;
    side = DoubleSide;
    shapePosition: any;

    private geo: Geometry;
    private linegeo: Geometry;
    private stylespec: any;
    private components: any;
    private shapeObj: any;
    private renderedShapeObj: any;
    /**
     * Custom renderable shape
     * 
     * @constructor 
     * 
     * @param {ShapeSpec} stylespec
     */
    constructor(stylespec: ShapeSpec) {

        this.stylespec = stylespec || {};

        this.boundingSphere = new Sphere();
        /** @type {IntersectionShapes} */
        this.intersectionShape = {
            sphere: [],
            cylinder: [],
            line: [],
            triangle: []
        };

        GLShape.updateFromStyle(this, this.stylespec);

        // Keep track of shape components and their centroids
        this.components = [];
        this.shapeObj = null;
        this.renderedShapeObj = null;

        this.geo = new Geometry(true);
        this.linegeo = new Geometry(true);
    };

    /** Update shape with new style specification
     * @param {ShapeSpec} newspec
       @example 
        let sphere = viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
        sphere.updateStyle({color:'yellow',opacity:0.5});
        viewer.render();
     */
    updateStyle(newspec: ShapeSpec) {

        for (var prop in newspec) {
            this.stylespec[prop] = newspec[prop];
        }

        GLShape.updateFromStyle(this, this.stylespec);

        if (newspec.voldata && newspec.volscheme) {
            adjustVolumeStyle(newspec);

            //convert volumetric data into colors
            const scheme = newspec.volscheme;
            const voldata = newspec.voldata;
            const cc = CC;
            const range = scheme.range() || [-1, 1];
            this.geo.setColors(function (x, y, z) {
                let val = voldata.getVal(x, y, z);
                let col = cc.color(scheme.valueToHex(val, range));
                return col;
            });
            delete this.color;
        }
    };

    /**
     * Creates a custom shape from supplied vertex and face arrays
     * @param {CustomShapeSpec} customSpec     
     */
    public addCustom(customSpec: CustomShapeSpec) {

        customSpec.vertexArr = customSpec.vertexArr || [];
        customSpec.faceArr = customSpec.faceArr || [];
        customSpec.normalArr = customSpec.normalArr || [];

        // will split mesh as needed
        GLShape.drawCustom(this, this.geo, customSpec);
    };

    /**
     * Creates a sphere shape
     * @param {SphereSpec} sphereSpec
     @example 
     viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
     
     viewer.render();
     */
    public addSphere(sphereSpec: SphereSpec) {

        if (!sphereSpec.center) {
            sphereSpec.center = new Vector3(0, 0, 0);
        }

        sphereSpec.radius = sphereSpec.radius ? clamp(sphereSpec.radius, 0, Infinity) : 1.5;
        sphereSpec.color = CC.color(sphereSpec.color);

        this.intersectionShape.sphere.push(new Sphere(sphereSpec.center, sphereSpec.radius));

        GLDraw.drawSphere(this.geo, sphereSpec.center,
            sphereSpec.radius, sphereSpec.color, sphereSpec.quality);

        this.components.push({
            centroid: new Vector3(sphereSpec.center.x,
                sphereSpec.center.y, sphereSpec.center.z)
        });
        var geoGroup = this.geo.updateGeoGroup(0);

        GLShape.updateBoundingFromPoints(this.boundingSphere, this.components,
            geoGroup.vertexArray, geoGroup.vertices);
    };


    /**
     * Creates a box
     * @param {BoxSpec} boxSpec
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
    public addBox(boxSpec: BoxSpec) {

        var dim = boxSpec.dimensions || { w: 1, h: 1, d: 1 };

        //dimensions may be scalar or vector quantities
        var w: XYZ;
        if (typeof (dim.w) == "number") {
            w = { x: dim.w, y: 0, z: 0 };
        } else {
            w = dim.w;
        }
        var h: XYZ;
        if (typeof (dim.h) == "number") {
            h = { x: 0, y: dim.h, z: 0 };
        } else {
            h = dim.h;
        }
        var d: XYZ;
        if (typeof (dim.d) == "number") {
            d = { x: 0, y: 0, z: dim.d };
        } else {
            d = dim.d;
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

        var spec = extend({}, boxSpec);
        spec.vertexArr = verts;
        spec.faceArr = faces;
        spec.normalArr = [];
        GLShape.drawCustom(this, this.geo, spec);

        var centroid = new Vector3();
        this.components.push({
            centroid: centroid.addVectors(uv[0], uv[7]).multiplyScalar(0.5)
        });
        var geoGroup = this.geo.updateGeoGroup(0);
        GLShape.updateBoundingFromPoints(this.boundingSphere, this.components, geoGroup.vertexArray, geoGroup.vertices);
    };

    /**
     * Creates a cylinder shape
     * @param {CylinderSpec} cylinderSpec
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
    public addCylinder(cylinderSpec: CylinderSpec) {

        var start: Vector3;
        var end: Vector3;
        if (!cylinderSpec.start) {
            start = new Vector3(0, 0, 0);
        } else {
            start = new Vector3(cylinderSpec.start.x || 0,
                cylinderSpec.start.y || 0, cylinderSpec.start.z || 0);
        }

        if (!cylinderSpec.end) {
            end = new Vector3(0, 0, 0);
        } else {
            end = new Vector3(cylinderSpec.end.x,
                cylinderSpec.end.y || 0, cylinderSpec.end.z || 0);
            if (typeof (end.x) == 'undefined') end.x = 3; //show something even if undefined
        }
        var radius = cylinderSpec.radius || 0.1;
        var color = CC.color(cylinderSpec.color);

        this.intersectionShape.cylinder.push(new Cylinder(start, end, radius));

        GLDraw.drawCylinder(this.geo, start, end, radius, color, cylinderSpec.fromCap, cylinderSpec.toCap);

        var centroid = new Vector3();
        this.components.push({
            centroid: centroid.addVectors(start, end).multiplyScalar(0.5)
        });
        var geoGroup = this.geo.updateGeoGroup(0);
        GLShape.updateBoundingFromPoints(this.boundingSphere, this.components,
            geoGroup.vertexArray, geoGroup.vertices);

    };

    /**
     * Creates a dashed cylinder shape
     * @param {CylinderSpec} cylinderSpec
     */
    public addDashedCylinder(cylinderSpec: CylinderSpec) {
  
        cylinderSpec.dashLength = cylinderSpec.dashLength || 0.25;
        cylinderSpec.gapLength = cylinderSpec.gapLength || 0.25;

        var start: Vector3;
        if (!cylinderSpec.start) start = new Vector3(0, 0, 0);
        else {
            start = new Vector3(cylinderSpec.start.x || 0,
                cylinderSpec.start.y || 0, cylinderSpec.start.z || 0);
        }

        var end: Vector3;
        if (!cylinderSpec.end) end = new Vector3(3, 0, 0);
        else {
            end = new Vector3(cylinderSpec.end.x,
                cylinderSpec.end.y || 0, cylinderSpec.end.z || 0);
            if (typeof (end.x) == 'undefined') end.x = 3; //show something even if undefined
        }

        var radius = cylinderSpec.radius || 0.1;
        var color = CC.color(cylinderSpec.color);

        var cylinderLength = Math.sqrt(Math.pow((start.x - end.x), 2) + Math.pow((start.y - end.y), 2) + Math.pow((start.z - end.z), 2));

        var count = cylinderLength / (cylinderSpec.gapLength + cylinderSpec.dashLength);

        var new_start = new Vector3(cylinderSpec.start.x || 0,
            cylinderSpec.start.y || 0, cylinderSpec.start.z || 0);
        var new_end = new Vector3(cylinderSpec.end.x,
            cylinderSpec.end.y || 0, cylinderSpec.end.z || 0);

        var gapVector = new Vector3((end.x - start.x) / (cylinderLength / cylinderSpec.gapLength), (end.y - start.y) / (cylinderLength / cylinderSpec.gapLength), (end.z - start.z) / (cylinderLength / cylinderSpec.gapLength));
        var dashVector = new Vector3((end.x - start.x) / (cylinderLength / cylinderSpec.dashLength), (end.y - start.y) / (cylinderLength / cylinderSpec.dashLength), (end.z - start.z) / (cylinderLength / cylinderSpec.dashLength));

        for (var place = 0; place < count; place++) {
            new_end = new Vector3(new_start.x + dashVector.x, new_start.y + dashVector.y, new_start.z + dashVector.z);

            this.intersectionShape.cylinder.push(new Cylinder(new_start, new_end, radius));

            GLDraw.drawCylinder(this.geo, new_start, new_end, radius, color, cylinderSpec.fromCap, cylinderSpec.toCap);

            new_start = new Vector3(new_end.x + gapVector.x, new_end.y + gapVector.y, new_end.z + gapVector.z);

        }
        var centroid = new Vector3();
        this.components.push({
            centroid: centroid.addVectors(start, end).multiplyScalar(0.5)
        });
        var geoGroup = this.geo.updateGeoGroup(0);
        GLShape.updateBoundingFromPoints(this.boundingSphere, this.components,
            geoGroup.vertexArray, geoGroup.vertices);
    };

    /**
     * Creates a curved shape
     * @param {CurveSpec} curveSpec
     */
    public addCurve(curveSpec: CurveSpec) {

        curveSpec.points = curveSpec.points || [];
        curveSpec.smooth = curveSpec.smooth || 10;
        if (typeof (curveSpec.fromCap) == "undefined") curveSpec.fromCap = 2;
        if (typeof (curveSpec.toCap) == "undefined") curveSpec.toCap = 2;

        //subdivide into smoothed spline points
        var points = subdivide_spline(curveSpec.points, curveSpec.smooth);

        if (points.length < 3) {
            console.log("Too few points in addCurve");
            return;
        }

        var radius = curveSpec.radius || 0.1;
        var color = CC.color(curveSpec.color);
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
                color: color as ColorSpec,
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
                color: color as ColorSpec,
                mid: 0.0001
            };
            this.addArrow(arrowspec);
        }

        var midway = Math.ceil(points.length / 2);
        var middleSpec: any = { radius: radius, color: color, fromCap: 2, toCap: 2 };
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

    /**
     * Creates a line shape
     * @param {LineSpec} lineSpec
     @example
     $3Dmol.download("pdb:2ABJ",viewer,{},function(){
              viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
              viewer.render(callback);
          });
    
     */
    public addLine(lineSpec: LineSpec) {

        var start: Vector3;
        var end: Vector3;
        if (!lineSpec.start) {
            start = new Vector3(0, 0, 0);
        } else {
            start = new Vector3(lineSpec.start.x || 0,
                lineSpec.start.y || 0, lineSpec.start.z || 0);
        }
        if (!lineSpec.end) {
            end = new Vector3(3, 0, 0);
        } else {
            end = new Vector3(lineSpec.end.x,
                lineSpec.end.y || 0, lineSpec.end.z || 0);
            if (typeof (end.x) == 'undefined') end.x = 3; //show something even if undefined
        }

        var geoGroup = this.geo.updateGeoGroup(2);

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

        var centroid = new Vector3();
        this.components.push({
            centroid: centroid.addVectors(start, end).multiplyScalar(0.5)
        });
        geoGroup = this.geo.updateGeoGroup(0);
        GLShape.updateBoundingFromPoints(this.boundingSphere, this.components,
            geoGroup.vertexArray, geoGroup.vertices);
    };


    /**
     * Creates an arrow shape
     * @param {ArrowSpec} arrowSpec
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
    public addArrow(arrowSpec: ArrowSpec) {

        if (!arrowSpec.start) {
            arrowSpec.start = new Vector3(0, 0, 0);
        } else {
            arrowSpec.start = new Vector3(arrowSpec.start.x || 0,
                arrowSpec.start.y || 0, arrowSpec.start.z || 0);
        }

        if (arrowSpec.dir instanceof Vector3 && typeof (arrowSpec.length) === 'number') {
            var end = arrowSpec.dir.clone().multiplyScalar(arrowSpec.length).add(
                arrowSpec.start);
            arrowSpec.end = end;
        }
        else if (!arrowSpec.end) {
            arrowSpec.end = new Vector3(3, 0, 0);
        } else {
            arrowSpec.end = new Vector3(arrowSpec.end.x,
                arrowSpec.end.y || 0, arrowSpec.end.z || 0);
            if (typeof (arrowSpec.end.x) == 'undefined') arrowSpec.end.x = 3; //show something even if undefined
        }

        arrowSpec.radius = arrowSpec.radius || 0.1;

        arrowSpec.radiusRatio = arrowSpec.radiusRatio || 1.618034;

        arrowSpec.mid = (0 < arrowSpec.mid && arrowSpec.mid < 1) ? arrowSpec.mid
            : 0.618034;

        GLShape.drawArrow(this, this.geo, arrowSpec);

        var centroid = new Vector3();
        this.components.push({
            centroid: centroid.addVectors(arrowSpec.start, arrowSpec.end)
                .multiplyScalar(0.5)
        });
        var geoGroup = this.geo.updateGeoGroup(0);
        GLShape.updateBoundingFromPoints(this.boundingSphere, this.components,
            geoGroup.vertexArray, geoGroup.vertices);
    };


    static distance_from(c1: XYZ, c2: XYZ) {
        return Math.sqrt(Math.pow((c1.x - c2.x), 2) + Math.pow((c1.y - c2.y), 2) + Math.pow((c1.z - c2.z), 2));
    };

    static inSelectedRegion(coordinate: XYZ, selectedRegion, radius: number) {

        for (var i = 0; i < selectedRegion.length; i++) {
            if (GLShape.distance_from(selectedRegion[i], coordinate) <= radius)
                return true;
        }
        return false;
    };

    /**
     * Create isosurface from voluemetric data.
     * @param {VolumeData} data - volumetric input data
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
    addIsosurface(data, volSpec:IsoSurfaceSpec, callback?) {//may want to cache the arrays geneerated when selectedRegion ==true

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
                bitdata[i] |= GLShape.ISDONE;

        }

        var verts = [], faces = [];
        MarchingCube.march(bitdata, verts, faces, {
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
            MarchingCube.laplacianSmooth(smoothness, verts, faces);
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
                    GLShape.inSelectedRegion(verts[i],
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

        GLShape.drawCustom(this, this.geo, {
            vertexArr: verts,
            faceArr: faces,
            normalArr: [],
            clickable: volSpec.clickable,
            hoverable: volSpec.hoverable
        });

        this.updateStyle(volSpec);

        //computing bounding sphere from vertices
        var origin = new Vector3(data.origin.x, data.origin.y, data.origin.z);
        var size = new Vector3(data.size.x * data.unit.x, data.size.y * data.unit.y, data.size.z * data.unit.z);

        var total = new Vector3(0, 0, 0);
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
     */
    public addVolumetricData(data, fmt, volSpec: IsoSurfaceSpec) {
        data = new VolumeData(data, fmt);
        this.addIsosurface(data, volSpec);
    };

    //for internal use, truncate buffers to save memory
    finalize() {
        GLShape.finalizeGeo(this.geo);
        this.geo.initTypedArrays();
        return this.geo;
    };

    /*
     * Initialize webgl objects for rendering
     * @param {$3Dmol.Object3D} group
     * 
     */
    globj(group) {

        if (this.renderedShapeObj) {
            group.remove(this.renderedShapeObj);
            this.renderedShapeObj = null;
        }

        if (this.hidden)
            return;
        GLShape.finalizeGeo(this.geo);
        this.geo.initTypedArrays();

        if (this.wireframe) {
            this.geo.setUpWireframe();
        }

        if (typeof (this.color) != 'undefined')
            GLShape.updateColor(this.geo, this.color);

        this.shapeObj = new Object3D();
        var material = null;

        if (this.side == DoubleSide) {
            material = new MeshDoubleLambertMaterial({
                wireframe: this.wireframe,
                side: this.side,
                transparent: (this.opacity < 1) ? true : false,
                opacity: this.opacity,
                wireframeLinewidth: this.linewidth,
                vertexColors: Coloring.VertexColors
            });
        } else {
            material = new MeshLambertMaterial({
                wireframe: this.wireframe,
                side: this.side,
                transparent: (this.opacity < 1) ? true : false,
                opacity: this.opacity,
                wireframeLinewidth: this.linewidth,
                vertexColors: Coloring.VertexColors
            });
        }

        var mesh = new Mesh(this.geo, material);

        this.shapeObj.add(mesh);

        var lineMaterial = new LineBasicMaterial({
            linewidth: this.linewidth,
            color: this.color
        });
        var line = new Line(this.linegeo, lineMaterial as Material, LineStyle.LinePieces);
        this.shapeObj.add(line);

        this.renderedShapeObj = this.shapeObj.clone();
        group.add(this.renderedShapeObj);

    };

    removegl(group) {
        if (this.renderedShapeObj) {
            // dispose of geos and materials
            if (this.renderedShapeObj.geometry !== undefined)
                this.renderedShapeObj.geometry.dispose();
            if (this.renderedShapeObj.material !== undefined)
                this.renderedShapeObj.material.dispose();
            group.remove(this.renderedShapeObj);
            this.renderedShapeObj = null;
        }
        this.shapeObj = null;
    };

    get position() {
        return this.boundingSphere.center;
    }

    get x() {
        return this.boundingSphere.center.x;
    }
    get y() {
        return this.boundingSphere.center.y;
    }
    get z() {
        return this.boundingSphere.center.z;
    }
};



export function splitMesh(mesh) {
    var MAXVERT = 64000; //webgl only supports 2^16 elements, leave a little breathing room (require at least 2)
    //peel off 64k vertices rsvh into their own mesh
    //duplicating vertices and normals as necessary to preserve faces and lines

    if (mesh.vertexArr.length < MAXVERT) return [mesh]; //typical case

    var slices: any = [{ vertexArr: [], normalArr: [], faceArr: [] }];
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

/**
 * GLShape style specification
 */
export interface ShapeSpec {
    /** Either a single color for the whole object or an array specifying the color at each vertex ({@link CustomShapeSpec}). */
    color?: ColorSpec | ColorSpec[];
    alpha?: number; //prefer opacity
    /** transparency, between 0 (invisible) and 1 (opaque) */
    opacity?: number;
    /** draw as wireframe, not surface */
    wireframe?: boolean;
    /** if true, do not display object */
    hidden?: boolean;
    /** width of line for wireframe rendering **No longer supported by most browsers** */
    linewidth?: number;
    /** if true, user can click on object to trigger callback */
    clickable?: boolean;
    /** function to call on click */
    callback?: Func;
    /** if true, user can hover on object to trigger callback */
    hoverable?: boolean;
    /** hover callback */
    hover_callback?: Func;
    /** unhover callback */
    unhover_callback?: Func;
    /** if true, user can right-click or long press to trigger callback */
    contextMenuEnabled?: boolean;
    /** if set, only display in this frame of an animation */
    frame?: number;
    side?: number;
    voldata?: VolumeData;
    volscheme?: GradientType
};


/**
 * Isosurface style specification
 * @extends ShapeSpec
 */
export interface IsoSurfaceSpec extends ShapeSpec {
    /** specifies the isovalue to draw surface at */
    isoval?: number;
    /** if true uses voxel style rendering */
    voxel?: boolean;
    /** amount to smooth surface (default 1) */
    smoothness?: number;
    /** coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates */
    coords?: XYZ[];
    /** distance around coords to include data [default = 2.0] */
    seldist?: number;
    /** volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified */
    voldata?: VolumeData;
    /** coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields. */
    volscheme?: GradientType;
    /** format of voldata if not a $3Dmol.VolumeData object */
    volformat?: string;

    selectedRegion?: XYZ[]; //deprecated
    selectedOffset?: number; //deprecated
    radius?: number; //also deprecated
};



/**
 * Arrow shape specification.  
  * @extends ShapeSpec
 */
export interface ArrowSpec extends ShapeSpec {
    /** starting position */
    start?: XYZ;
    /** ending position */
    end?: XYZ;
    /** direction to extend from start (instead of specifying end) */
    dir?: XYZ;
    /** length to extend in dir direction from start (instead of specifying end) */
    length?: number;
    /** radius (default 0.1A) */
    radius?: number;
    /** color */
    color?: ColorSpec;
    /** hidden */
    hidden?: boolean;
    /** ratio of arrow base to cylinder (1.618034 default) */
    radiusRatio?: number;
    /** relative position of arrow base (0.618034 default) */
    mid?: number;
    /** position of arrow base in length units, if negative positioned from end instead of start.  Overrides mid. */
    midpos?: number;
};

/**
 * Cylinder shape specification.  
 * @extends ShapeSpec
 * 
 */
export interface CylinderSpec extends ShapeSpec {
    /** starting vector */
    start?: XYZ;
    /** ending position */
    end?: XYZ;
    /** radius */
    radius?: number;
    /** Place a cap at the start (none, flat or round) */
    fromCap?: CAP;
    /** Place a cap at the end (none, flat or round) */
    toCap?: CAP;
    /** Make the cylinder dashed. */
    dashed?: boolean;
    /** Length of dashes (default 0.25) */
    dashLength?: number;
    /** Length of gaps (default 0.25) */
    gapLength?: number;
};

/**
 * Curve shape specification.  
 * @extends ShapeSpec
 */
export interface CurveSpec extends ShapeSpec {
    /** Sequence of points to draw curve through */
    points?: XYZ[];
    /** amount of interpolation */
    smooth?: number;
    /** radius of curve */
    radius?: number;
    /** if an arrow should be drawn at the start */
    fromArrow?: boolean;
    /** if an arrow should be drawn at the end */
    toArrow?: boolean;
    /** have cap at start */
    fromCap?: CAP;
    /** have cap at end */
    toCap?: CAP;
};

/**
 * Line shape specification.  Default to wireframe.
 * @extends ShapeSpec
 */
export interface LineSpec extends ShapeSpec {
    /** Starting position */
    start?: XYZ;
    /** Ending position */
    end?: XYZ;
    /** make dashed */
    dashed?: boolean;
};

/**
 * Box shape specification. 
 * @extends ShapeSpec
 */
export interface BoxSpec extends ShapeSpec {
    /** bottom corner of box */
    corner?: XYZ;
    /** center of box */
    center?: XYZ;
    /** width, height, depth of box */
    dimensions?: {
        w: number | XYZ;
        h: number | XYZ;
        d: number | XYZ;
    };
};


/**
 * Specification for adding custom shape. 
 * @extends ShapeSpec
 */
export interface CustomShapeSpec extends ShapeSpec {
    /** List of vertex positions */
    vertexArr?: XYZ[];
    /** List of normal vectors for each vertex */
    normalArr?: XYZ[];
    /** List of triangles to build the custom shape. Each triangle is defined by the indices of 3 vertices in vertexArr, so the array length should be 3 times the number of faces. */
    faceArr?: number[];
};

/**
 * Sphere shape specification. Extends {@link ShapeSpec}.
 */
export interface SphereSpec extends ShapeSpec {
    /** center of sphere */
    center?: XYZ;
    /** radius of sphere */
    radius?: number;
    /** quality metric, higher uses more triangles (default 2) */
    quality?: number;
};

