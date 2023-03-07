//glcartoon.js
//This contains all the routines for rendering a cartoon given a set
//of atoms with assigned secondary structure

import { Vector3 } from "./WebGL/math";
import { Triangle, Sphere } from "./WebGL/shapes";
import { MeshDoubleLambertMaterial, Mesh, Geometry, Material, Coloring } from "./WebGL";
import { Gradient } from "./Gradient";
import { CC, ColorSpec } from "./colors";
import { GLDraw } from "./GLDraw";
import { isNumeric, getColorFromStyle } from "./utilities";

/**
 * A visualization of protein or nucleic acid secondary structure.  Applying this to other molecules will not show anything.
  In nucleic acids, the base cylinders obtain their color from the  atom to which the cylinder is drawn, which is 'N1' for purines (resn: 'A', 'G', 'DA', 'DG') and 
  'N3' for pyrimidines (resn: 'C', 'U', 'DC', 'DT'). 
  The different nucleobases can therefore be distinguished as by setting the colors 
  of each of these atoms. The backbone color is set from the 'P' atoms ('O5' for the 5' terminus).

 * 
 * @example $3Dmol.download("pdb:4ZD3",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  viewer.setViewStyle({style:"outline"});
                  viewer.setStyle({},{cartoon:{}});
                  viewer.render();
              });
 */
export interface CartoonStyleSpec {
    /** do not show  */
    hidden?: boolean;
    /** strand color, may specify as 'spectrum' which will apply reversed gradient based on residue number */
    color?: ColorSpec;
    /**  style of cartoon rendering (trace, oval, rectangle (default), parabola, edged) */
    style?: string;
    /**  whether to use constant strand width, disregarding
    *       secondary structure; use thickness to adjust radius */
    ribbon?: boolean;
    /** whether to add arrows showing beta-sheet
    *       directionality; does not apply to trace or ribbon */
    arrows?: boolean;
    /** whether to display alpha helices as simple cylinders;
    *       does not apply to trace */
    tubes?: boolean;
    /** cartoon strand thickness, default is 0.4 */
    thickness?: number;
    /** cartoon strand width, default is secondary
    *       structure-dependent; does not apply to trace or ribbon */
    width?: number;
    /** set opacity from 0-1; transparency is set per-chain
    *       with a warning outputted in the event of ambiguity */
    opacity?: number
};

// helper functions
// Catmull-Rom subdivision
export function subdivide_spline(_points, DIV) { // points as Vector3
    var ret = [];
    var points = _points;
    points = []; // Smoothing test
    points.push(_points[0]);

    var i, lim, size;
    var p0, p1, p2, p3, v0, v1;

    for (i = 1, lim = _points.length - 1; i < lim; i++) {
        p1 = _points[i];
        p2 = _points[i + 1];
        if (p1.smoothen) {
            var np = new Vector3((p1.x + p2.x) / 2,
                (p1.y + p2.y) / 2, (p1.z + p2.z) / 2);
            np.atom = p1.atom;
            points.push(np);
        }
        else
            points.push(p1);
    }
    points.push(_points[_points.length - 1]);

    for (i = -1, size = points.length; i <= size - 3; i++) {
        p0 = points[(i === -1) ? 0 : i];
        p1 = points[i + 1];
        p2 = points[i + 2];
        p3 = points[(i === size - 3) ? size - 1 : i + 3];
        v0 = new Vector3().subVectors(p2, p0).multiplyScalar(0.5);
        v1 = new Vector3().subVectors(p3, p1).multiplyScalar(0.5);
        if (p2.skip)
            continue;

        for (var j = 0; j < DIV; j++) {
            var t = 1.0 / DIV * j;
            var x = p1.x + t * v0.x + t * t *
                (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x) + t * t * t *
                (2 * p1.x - 2 * p2.x + v0.x + v1.x);
            var y = p1.y + t * v0.y + t * t *
                (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y) + t * t * t *
                (2 * p1.y - 2 * p2.y + v0.y + v1.y);
            var z = p1.z + t * v0.z + t * t *
                (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z) + t * t * t *
                (2 * p1.z - 2 * p2.z + v0.z + v1.z);

            var pt = new Vector3(x, y, z);
            if (j < DIV / 2) {
                pt.atom = p1.atom;
            } else {
                pt.atom = p2.atom;
            }

            ret.push(pt);
        }
    }
    ret.push(points[points.length - 1]);
    return ret;
};


const coilWidth = 0.5;
const helixSheetWidth = 1.3;
const nucleicAcidWidth = 0.8;
const defaultThickness = 0.4;
const baseThickness    = 0.4;

function drawThinStrip(geo: Geometry, p1, p2, colors) {

    var offset, vertoffset;
    var color, colori;

    for (var i = 0, lim = p1.length; i < lim; i++) {

        colori = Math.round(i * (colors.length - 1) / lim);
        color = CC.color(colors[colori]);

        var geoGroup = geo.updateGeoGroup(2);
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;
        var faceArray = geoGroup.faceArray;
        offset = geoGroup.vertices;
        vertoffset = offset * 3;

        vertexArray[vertoffset] = p1[i].x;
        vertexArray[vertoffset + 1] = p1[i].y;
        vertexArray[vertoffset + 2] = p1[i].z;

        vertexArray[vertoffset + 3] = p2[i].x;
        vertexArray[vertoffset + 4] = p2[i].y;
        vertexArray[vertoffset + 5] = p2[i].z;

        for (var j = 0; j < 6; ++j) {
            colorArray[vertoffset + 3 * j] = color.r;
            colorArray[vertoffset + 1 + 3 * j] = color.g;
            colorArray[vertoffset + 2 + 3 * j] = color.b;
        }

        if (i > 0) {
            var faces = [offset, offset + 1, offset - 1, offset - 2];
            var faceoffset = geoGroup.faceidx;

            faceArray[faceoffset] = faces[0];
            faceArray[faceoffset + 1] = faces[1];
            faceArray[faceoffset + 2] = faces[3];
            faceArray[faceoffset + 3] = faces[1];
            faceArray[faceoffset + 4] = faces[2];
            faceArray[faceoffset + 5] = faces[3];

            geoGroup.faceidx += 6;
        }

        geoGroup.vertices += 2;
    }

};

function drawShapeStrip(geo: Geometry, points, colors, div, thickness, opacity, shape) {

    // points is a 2D array, dimensionality given by [num = cross-sectional
    // resolution][len = length of strip]
    var i, j, num, len;
    num = points.length;
    if (num < 2 || points[0].length < 2)
        return;

    for (i = 0; i < num; i++) { // spline to generate greater length-wise
        // resolution
        points[i] = subdivide_spline(points[i], div);
    }
    len = points[0].length;

    if (!thickness) // if thickness is 0, we can use a smaller geometry than
        // this function generates
        return drawThinStrip(geo, points[0], points[num - 1], colors);

    var axis, cs_shape, cs_bottom, cs_top, last_cs_bottom, last_cs_top;

    // cache the available cross-sectional shapes
    var cs_ellipse = [], cs_rectangle = [], cs_parabola = [];
    for (j = 0; j < num; j++) {
        cs_ellipse.push(0.25 + 1.5 *
            Math.sqrt((num - 1) * j - Math.pow(j, 2)) / (num - 1));
        cs_rectangle.push(0.5);
        cs_parabola.push(2 * (Math.pow(j / num, 2) - j / num) + 0.6);
    }

    /*
     * face_refs array is used to generate faces from vertexArray
     * iteratively. As we move through each cross-sectional segment of
     * points, we draw lateral faces backwards to the previous
     * cross-sectional segment.
     * 
     * To correctly identify the points needed to make each face we use this
     * array as a lookup table for the relative indices of each needed point
     * in the vertices array.
     * 
     * 4 points are used to create 2 faces.
     */

    var face_refs = [];
    for (j = 0; j < num * 2 - 1; j++) {
        /*
         * [curr vertex in curr cross-section, next vertex in curr
         * cross-section, next vertex in prev cross-section, curr vertex in
         * prev cross-section]
         */
        face_refs[j] = [j, j + 1, j + 1 - 2 * num, j - 2 * num];
    }
    // last face is different. easier to conceptualize this by drawing a
    // diagram
    face_refs[num * 2 - 1] = [j, j + 1 - 2 * num, j + 1 - 4 * num,
        j - 2 * num];

    var v_offset, va_offset, f_offset;
    var currentAtom;
    var color, colori;
    var vertexArray, colorArray, faceArray, face;
    let geoGroup = geo.updateGeoGroup();

    for (i = 0; i < len; i++) {
        let gnum = geo.groups;
        let replicating = false;
        geoGroup = geo.updateGeoGroup(2 * num); // ensure vertex capacity
        if (gnum != geo.groups && i > 0) {
            //we created a new geo - need to replicate vertices at edge 
            //(but not faces)
            i = i - 1;
            replicating = true;
        }
        colori = Math.round(i * (colors.length - 1) / len);
        color = CC.color(colors[colori]);

        last_cs_bottom = cs_bottom;
        last_cs_top = cs_top;
        cs_bottom = [];
        cs_top = [];
        axis = [];

        if (points[0][i].atom !== undefined) // TODO better edge case
        // handling
        {
            currentAtom = points[0][i].atom;
            if (shape === "oval")
                cs_shape = cs_ellipse;
            else if (shape === "rectangle")
                cs_shape = cs_rectangle;
            else if (shape === "parabola")
                cs_shape = cs_parabola;
        }
        if (!cs_shape)
            cs_shape = cs_rectangle;

        // calculate thickness at each width point, from cross-sectional
        // shape
        var toNext, toSide;
        for (j = 0; j < num; j++) {
            if (i < len - 1)
                toNext = points[j][i + 1].clone().sub(points[j][i]);
            else
                toNext = points[j][i - 1].clone().sub(points[j][i])
                    .negate();

            if (j < num - 1)
                toSide = points[j + 1][i].clone().sub(points[j][i]);
            else
                toSide = points[j - 1][i].clone().sub(points[j][i])
                    .negate();

            axis[j] = toSide.cross(toNext).normalize().multiplyScalar(
                thickness * cs_shape[j]);
        }

        // generate vertices by applying cross-sectional shape thickness to
        // input points
        for (j = 0; j < num; j++)
            cs_bottom[j] = points[j][i].clone().add(
                axis[j].clone().negate());
        for (j = 0; j < num; j++)
            cs_top[j] = points[j][i].clone().add(axis[j]);

        /*
         * Until this point the vertices have been dealt with as
         * Vector3() objects, but we need to serialize them into the
         * geoGroup.vertexArray, where every three indices represents the
         * next vertex. The colorArray is analogous.
         * 
         * In the following for-loops, j iterates through VERTICES so we
         * need to index them in vertexArray by 3*j + either 0, 1, or 2 for
         * xyz or rgb component.
         */

        vertexArray = geoGroup.vertexArray;
        colorArray = geoGroup.colorArray;
        faceArray = geoGroup.faceArray;
        v_offset = geoGroup.vertices;
        va_offset = v_offset * 3; // in case geoGroup already contains
        // vertices

        // bottom edge of cross-section, vertices [0, num)
        for (j = 0; j < num; j++) {
            vertexArray[va_offset + 3 * j + 0] = cs_bottom[j].x;
            vertexArray[va_offset + 3 * j + 1] = cs_bottom[j].y;
            vertexArray[va_offset + 3 * j + 2] = cs_bottom[j].z;
        }

        // top edge of cross-section, vertices [num, 2*num)
        // add these backwards, so that each cross-section's vertices are
        // added sequentially to vertexArray
        for (j = 0; j < num; j++) {
            vertexArray[va_offset + 3 * j + 0 + 3 * num] = cs_top[num - 1 - j].x;
            vertexArray[va_offset + 3 * j + 1 + 3 * num] = cs_top[num - 1 - j].y;
            vertexArray[va_offset + 3 * j + 2 + 3 * num] = cs_top[num - 1 - j].z;
        }

        for (j = 0; j < 2 * num; ++j) {
            colorArray[va_offset + 3 * j + 0] = color.r;
            colorArray[va_offset + 3 * j + 1] = color.g;
            colorArray[va_offset + 3 * j + 2] = color.b;
        }

        if (i > 0 && !replicating) {

            for (j = 0; j < num * 2; j++) {

                // get VERTEX indices of the 4 points of a rectangular face
                // (as opposed to literal vertexArray indices)
                face = [v_offset + face_refs[j][0],
                v_offset + face_refs[j][1],
                v_offset + face_refs[j][2],
                v_offset + face_refs[j][3]];

                f_offset = geoGroup.faceidx;

                // need 2 triangles to draw a face between 4 points
                faceArray[f_offset] = face[0];
                faceArray[f_offset + 1] = face[1];
                faceArray[f_offset + 2] = face[3];

                faceArray[f_offset + 3] = face[1];
                faceArray[f_offset + 4] = face[2];
                faceArray[f_offset + 5] = face[3];

                geoGroup.faceidx += 6;

                // TODO implement clickable the right way. midpoints of
                // strand between consecutive atoms
            }

            if (currentAtom.clickable || currentAtom.hoverable) {
                var faces = [];

                faces.push(new Triangle(last_cs_bottom[0],
                    cs_bottom[0], cs_bottom[num - 1]));
                faces.push(new Triangle(last_cs_bottom[0],
                    cs_bottom[num - 1], last_cs_bottom[num - 1]));

                faces.push(new Triangle(last_cs_bottom[num - 1],
                    cs_bottom[num - 1], cs_top[num - 1]));
                faces.push(new Triangle(last_cs_bottom[num - 1],
                    cs_top[num - 1], last_cs_top[num - 1]));

                faces.push(new Triangle(cs_top[0], last_cs_top[0],
                    last_cs_top[num - 1]));
                faces.push(new Triangle(cs_top[num - 1], cs_top[0],
                    last_cs_top[num - 1]));

                faces.push(new Triangle(cs_bottom[0],
                    last_cs_bottom[0], last_cs_top[0]));
                faces.push(new Triangle(cs_top[0], cs_bottom[0],
                    last_cs_top[0]));

                for (j in faces) {
                    currentAtom.intersectionShape.triangle.push(faces[j]);
                }
            }
        }

        geoGroup.vertices += 2 * num;
    }

    // for terminal faces
    vertexArray = geoGroup.vertexArray;
    colorArray = geoGroup.colorArray;
    faceArray = geoGroup.faceArray;
    v_offset = geoGroup.vertices;
    va_offset = v_offset * 3;
    f_offset = geoGroup.faceidx;

    for (i = 0; i < num - 1; i++) // "rear" face
    {
        face = [i, i + 1, 2 * num - 2 - i, 2 * num - 1 - i];

        f_offset = geoGroup.faceidx;

        faceArray[f_offset] = face[0];
        faceArray[f_offset + 1] = face[1];
        faceArray[f_offset + 2] = face[3];

        faceArray[f_offset + 3] = face[1];
        faceArray[f_offset + 4] = face[2];
        faceArray[f_offset + 5] = face[3];

        geoGroup.faceidx += 6;
    }

    for (i = 0; i < num - 1; i++) // "front" face
    {
        face = [v_offset - 1 - i, v_offset - 2 - i,
        v_offset - 2 * num + i + 1, v_offset - 2 * num + i];

        f_offset = geoGroup.faceidx;

        faceArray[f_offset] = face[0];
        faceArray[f_offset + 1] = face[1];
        faceArray[f_offset + 2] = face[3];

        faceArray[f_offset + 3] = face[1];
        faceArray[f_offset + 4] = face[2];
        faceArray[f_offset + 5] = face[3];

        geoGroup.faceidx += 6;
    }

};

function drawPlainStrip(geo, points, colors, div, thickness, opacity) {
    if ((points.length) < 2)
        return;

    var p1, p2;
    p1 = points[0];
    p2 = points[points.length - 1];

    p1 = subdivide_spline(p1, div);
    p2 = subdivide_spline(p2, div);
    if (!thickness)
        return drawThinStrip(geo, p1, p2, colors);

    // var vs = geo.vertices, fs = geo.faces;
    var vs = [];
    var axis, p1v, p2v, a1v, a2v;

    var faces = [[0, 2, -6, -8], [-4, -2, 6, 4], [7, -1, -5, 3],
    [-3, 5, 1, -7]];

    var offset, vertoffset, faceoffset;
    var color, colori;
    var currentAtom, lastAtom;
    var i, lim, j;
    var face1, face2, face3;
    var geoGroup, vertexArray, colorArray, faceArray;

    for (i = 0, lim = p1.length; i < lim; i++) {

        colori = Math.round(i * (colors.length - 1) / lim);
        color = CC.color(colors[colori]);

        vs.push(p1v = p1[i]); // 0
        vs.push(p1v); // 1
        vs.push(p2v = p2[i]); // 2
        vs.push(p2v); // 3
        if (i < lim - 1) {
            var toNext = p1[i + 1].clone().sub(p1[i]);
            var toSide = p2[i].clone().sub(p1[i]);
            axis = toSide.cross(toNext).normalize().multiplyScalar(
                thickness);
        }
        vs.push(a1v = p1[i].clone().add(axis)); // 4
        vs.push(a1v); // 5
        vs.push(a2v = p2[i].clone().add(axis)); // 6
        vs.push(a2v); // 7

        if (p1v.atom !== undefined)
            currentAtom = p1v.atom;

        geoGroup = geo.updateGeoGroup(8);
        vertexArray = geoGroup.vertexArray;
        colorArray = geoGroup.colorArray;
        faceArray = geoGroup.faceArray;
        offset = geoGroup.vertices;
        vertoffset = offset * 3;

        vertexArray[vertoffset] = p1v.x;
        vertexArray[vertoffset + 1] = p1v.y;
        vertexArray[vertoffset + 2] = p1v.z;
        vertexArray[vertoffset + 3] = p1v.x;
        vertexArray[vertoffset + 4] = p1v.y;
        vertexArray[vertoffset + 5] = p1v.z;
        vertexArray[vertoffset + 6] = p2v.x;
        vertexArray[vertoffset + 7] = p2v.y;
        vertexArray[vertoffset + 8] = p2v.z;
        vertexArray[vertoffset + 9] = p2v.x;
        vertexArray[vertoffset + 10] = p2v.y;
        vertexArray[vertoffset + 11] = p2v.z;
        vertexArray[vertoffset + 12] = a1v.x;
        vertexArray[vertoffset + 13] = a1v.y;
        vertexArray[vertoffset + 14] = a1v.z;
        vertexArray[vertoffset + 15] = a1v.x;
        vertexArray[vertoffset + 16] = a1v.y;
        vertexArray[vertoffset + 17] = a1v.z;
        vertexArray[vertoffset + 18] = a2v.x;
        vertexArray[vertoffset + 19] = a2v.y;
        vertexArray[vertoffset + 20] = a2v.z;
        vertexArray[vertoffset + 21] = a2v.x;
        vertexArray[vertoffset + 22] = a2v.y;
        vertexArray[vertoffset + 23] = a2v.z;

        for (j = 0; j < 8; ++j) {
            colorArray[vertoffset + 3 * j] = color.r;
            colorArray[vertoffset + 1 + 3 * j] = color.g;
            colorArray[vertoffset + 2 + 3 * j] = color.b;
        }

        if (i > 0) {

            // both points have distinct atoms
            var diffAtoms = ((lastAtom !== undefined && currentAtom !== undefined) && lastAtom.serial !== currentAtom.serial);

            for (j = 0; j < 4; j++) {

                var face = [offset + faces[j][0], offset + faces[j][1],
                offset + faces[j][2], offset + faces[j][3]];

                faceoffset = geoGroup.faceidx;

                faceArray[faceoffset] = face[0];
                faceArray[faceoffset + 1] = face[1];
                faceArray[faceoffset + 2] = face[3];
                faceArray[faceoffset + 3] = face[1];
                faceArray[faceoffset + 4] = face[2];
                faceArray[faceoffset + 5] = face[3];

                geoGroup.faceidx += 6;

                if (currentAtom.clickable || lastAtom.clickable || currentAtom.hoverable || lastAtom.hoverable) {

                    var p1a = vs[face[3]].clone(), p1b = vs[face[0]]
                        .clone(), p2a = vs[face[2]].clone(), p2b = vs[face[1]]
                            .clone();

                    p1a.atom = vs[face[3]].atom || null; // should be
                    // same
                    p2a.atom = vs[face[2]].atom || null;

                    p1b.atom = vs[face[0]].atom || null; // should be
                    // same
                    p2b.atom = vs[face[1]].atom || null;

                    if (diffAtoms) {
                        var m1 = p1a.clone().add(p1b).multiplyScalar(0.5);
                        var m2 = p2a.clone().add(p2b).multiplyScalar(0.5);
                        var m = p1a.clone().add(p2b).multiplyScalar(0.5);

                        if (j % 2 === 0) {
                            if (lastAtom.clickable || lastAtom.hoverable) {
                                face1 = new Triangle(m1, m, p1a);
                                face2 = new Triangle(m2, p2a, m);
                                face3 = new Triangle(m, p2a, p1a);
                                lastAtom.intersectionShape.triangle
                                    .push(face1);
                                lastAtom.intersectionShape.triangle
                                    .push(face2);
                                lastAtom.intersectionShape.triangle
                                    .push(face3);
                            }

                            if (currentAtom.clickable || currentAtom.hoverable) {
                                face1 = new Triangle(p1b, p2b, m);
                                face2 = new Triangle(p2b, m2, m);
                                face3 = new Triangle(p1b, m, m1);
                                currentAtom.intersectionShape.triangle
                                    .push(face1);
                                currentAtom.intersectionShape.triangle
                                    .push(face2);
                                currentAtom.intersectionShape.triangle
                                    .push(face3);
                            }
                        } else {
                            if (currentAtom.clickable || currentAtom.hoverable) {
                                face1 = new Triangle(m1, m, p1a);
                                face2 = new Triangle(m2, p2a, m);
                                face3 = new Triangle(m, p2a, p1a);
                                currentAtom.intersectionShape.triangle
                                    .push(face1);
                                currentAtom.intersectionShape.triangle
                                    .push(face2);
                                currentAtom.intersectionShape.triangle
                                    .push(face3);
                            }

                            if (lastAtom.clickable || lastAtom.hoverable) {
                                face1 = new Triangle(p1b, p2b, m);
                                face2 = new Triangle(p2b, m2, m);
                                face3 = new Triangle(p1b, m, m1);
                                lastAtom.intersectionShape.triangle
                                    .push(face1);
                                lastAtom.intersectionShape.triangle
                                    .push(face2);
                                lastAtom.intersectionShape.triangle
                                    .push(face3);
                            }
                        }

                    }

                    // face for single atom
                    else if (currentAtom.clickable || currentAtom.hoverable) {
                        face1 = new Triangle(p1b, p2b, p1a);
                        face2 = new Triangle(p2b, p2a, p1a);
                        currentAtom.intersectionShape.triangle.push(face1);
                        currentAtom.intersectionShape.triangle.push(face2);
                    }

                }

            }
        }

        geoGroup.vertices += 8;
        lastAtom = currentAtom;
    }

    var vsize = vs.length - 8; // Cap

    geoGroup = geo.updateGeoGroup(8);
    vertexArray = geoGroup.vertexArray;
    colorArray = geoGroup.colorArray;
    faceArray = geoGroup.faceArray;
    offset = geoGroup.vertices;
    vertoffset = offset * 3;
    faceoffset = geoGroup.faceidx;

    for (i = 0; i < 4; i++) {
        vs.push(vs[i * 2]);
        vs.push(vs[vsize + i * 2]);

        var v1 = vs[i * 2], v2 = vs[vsize + i * 2];

        vertexArray[vertoffset + 6 * i] = v1.x;
        vertexArray[vertoffset + 1 + 6 * i] = v1.y;
        vertexArray[vertoffset + 2 + 6 * i] = v1.z;
        vertexArray[vertoffset + 3 + 6 * i] = v2.x;
        vertexArray[vertoffset + 4 + 6 * i] = v2.y;
        vertexArray[vertoffset + 5 + 6 * i] = v2.z;

        colorArray[vertoffset + 6 * i] = color.r;
        colorArray[vertoffset + 1 + 6 * i] = color.g;
        colorArray[vertoffset + 2 + 6 * i] = color.b;
        colorArray[vertoffset + 3 + 6 * i] = color.r;
        colorArray[vertoffset + 4 + 6 * i] = color.g;
        colorArray[vertoffset + 5 + 6 * i] = color.b;

    }

    vsize += 8;

    face1 = [offset, offset + 2, offset + 6, offset + 4];
    face2 = [offset + 1, offset + 5, offset + 7, offset + 3];

    faceArray[faceoffset] = face1[0];
    faceArray[faceoffset + 1] = face1[1];
    faceArray[faceoffset + 2] = face1[3];
    faceArray[faceoffset + 3] = face1[1];
    faceArray[faceoffset + 4] = face1[2];
    faceArray[faceoffset + 5] = face1[3];
    faceArray[faceoffset + 6] = face2[0];
    faceArray[faceoffset + 7] = face2[1];
    faceArray[faceoffset + 8] = face2[3];
    faceArray[faceoffset + 9] = face2[1];
    faceArray[faceoffset + 10] = face2[2];
    faceArray[faceoffset + 11] = face2[3];

    geoGroup.faceidx += 12;
    geoGroup.vertices += 8;

    // TODO: Add intersection planes for caps
};


function drawStrip(geo, points, colors, div, thickness, opacity, shape) {
    if (!shape || shape === "default")
        shape = "rectangle";
    if (shape === 'edged')
        drawPlainStrip(geo, points, colors, div, thickness, opacity);
    else if (shape === "rectangle" || shape === "oval" || shape === "parabola")
        drawShapeStrip(geo, points, colors, div, thickness, opacity, shape);
};

// check if given atom is an alpha carbon
function isAlphaCarbon(atom) {
    return atom && atom.elem === "C" && atom.atom === "CA"; // note that
    // calcium is
    // also CA
};

// check whether two atoms are members of the same residue or subsequent,
// connected residues (a before b)
function inConnectedResidues(a, b) {
    if (a && b && a.chain === b.chain) {
        if (!a.hetflag && !b.hetflag && (a.reschain === b.reschain) &&
            (a.resi === b.resi || a.resi === b.resi - 1))
            return true;
        if (a.resi < b.resi) {
            // some PDBs have gaps in the numbering but the residues are
            // still connected
            // assume if within 4A they are connected
            var dx = a.x - b.x;
            var dy = a.y - b.y;
            var dz = a.z - b.z;
            var dist = dx * dx + dy * dy + dz * dz;
            if (a.atom == "CA" && b.atom == "CA" && dist < 16.0) //protein residues not connected
                return true; // calpha dist
            else if ((a.atom == "P" || b.atom == "P") && dist < 64.0) //dna
                return true;
        }
    }

    return false;
};

// add geo to the group
function setGeo(group, geo, opacity, outline, setNormals) {

    if (geo == null || geo.vertices == 0) return;
    if (setNormals) {
        geo.initTypedArrays();
        geo.setUpNormals();
    }

    var cartoonMaterial = new MeshDoubleLambertMaterial();
    cartoonMaterial.vertexColors = Coloring.FaceColors;
    if (typeof (opacity) === "number" && opacity >= 0 && opacity < 1) {
        cartoonMaterial.transparent = true;
        cartoonMaterial.opacity = opacity;
    }
    cartoonMaterial.outline = outline;
    var cartoonMesh = new Mesh(geo, cartoonMaterial as Material);
    group.add(cartoonMesh);
};

function addBackbonePoints(points, num, smoothen, backbonePt,
    orientPt, prevOrientPt, backboneAtom, atoms, atomi) {
    var widthScalar, i, delta, v, addArrowPoints, testStyle;

    if (!backbonePt || !orientPt || !backboneAtom)
        return;

    // the side vector points along the axis from backbone atom to
    // orientation atom (eg. CA to O, in peptides)
    var sideVec = orientPt.sub(backbonePt);
    sideVec.normalize();

    //find next atom like this one
    var forwardVec = atoms[atomi];
    for (i = atomi + 1; i < atoms.length; i++) {
        forwardVec = atoms[i];
        if (forwardVec.atom == backboneAtom.atom)
            break;
    }
    // the forward vector points along the axis from backbone atom to next
    // backbone atom
    forwardVec = forwardVec ? new Vector3(forwardVec.x,
        forwardVec.y, forwardVec.z) : new Vector3(0, 0, 0);
    forwardVec.sub(backbonePt);

    // adjustments for proper beta arrow appearance
    if (backboneAtom.ss === "arrow start") {
        var adjustment = forwardVec.clone().multiplyScalar(0.3).cross(
            orientPt); // adjust perpendicularly to strand face
        backbonePt.add(adjustment);

        var upVec = forwardVec.clone().cross(sideVec).normalize();
        sideVec.rotateAboutVector(upVec, 0.43);
    }

    // determine from cartoon style or secondary structure how wide the
    // strand should be here
    // ribbon shape should have same width as thickness
    if (backboneAtom.style.cartoon.ribbon) {
        widthScalar = backboneAtom.style.cartoon.thickness || defaultThickness;

    } else // depending on secondary structure, multiply the orientation
    // vector by some scalar
    {
        if (!backboneAtom.style.cartoon.width) {
            if (backboneAtom.ss === "c") {
                if (backboneAtom.atom === "P")
                    widthScalar = nucleicAcidWidth;
                else
                    widthScalar = coilWidth;
            } else if (backboneAtom.ss === "arrow start") {
                widthScalar = helixSheetWidth;
                addArrowPoints = true;

            } else if (backboneAtom.ss === "arrow end")
                widthScalar = coilWidth;

            else if (backboneAtom.ss === "h" &&
                backboneAtom.style.cartoon.tubes ||
                backboneAtom.ss === "tube start")
                widthScalar = coilWidth;

            else
                widthScalar = helixSheetWidth;
        } else
            widthScalar = backboneAtom.style.cartoon.width;
    }

    // make sure the strand orientation doesn't twist more than 90 degrees
    if (prevOrientPt != null && sideVec.dot(prevOrientPt) < 0)
        sideVec.negate();

    sideVec.multiplyScalar(widthScalar);
    for (i = 0; i < num; i++) {
        // produces NUM incremental points from backbone atom minus
        // orientation vector
        // to backbone atom plus orientation vector
        delta = -1 + i * 2 / (num - 1); // -1 to 1 incrementing by num
        v = new Vector3(backbonePt.x + delta * sideVec.x,
            backbonePt.y + delta * sideVec.y, backbonePt.z + delta * sideVec.z);
        v.atom = backboneAtom;
        if (smoothen && backboneAtom.ss === "s")
            v.smoothen = true;
        points[i].push(v); // a num-length array of arrays, where each
        // inner array contains length-wise points
        // along the backbone offset by some constant pertaining to its cell
        // in the outer array
    }

    if (addArrowPoints) {

        sideVec.multiplyScalar(2);
        for (i = 0; i < num; i++) {
            delta = -1 + i * 2 / (num - 1); // -1 to 1 incrementing by num
            v = new Vector3(backbonePt.x + delta * sideVec.x,
                backbonePt.y + delta * sideVec.y, backbonePt.z + delta * sideVec.z);
            v.atom = backboneAtom;
            v.smoothen = false;
            v.skip = true;
            points[i].push(v);
        }
    }

    // make sure the strand is all the same style
    testStyle = backboneAtom.style.cartoon.style || 'default';
    if (points.style) {
        if (points.style != testStyle) {
            console
                .log("Warning: a cartoon chain's strand-style is ambiguous");
            points.style = 'default';
        }

    } else
        points.style = testStyle;

    // revert ss keywords used for arrow rendering back to original value
    if (backboneAtom.ss === "arrow start" || backboneAtom.ss === "arrow end")
        backboneAtom.ss = "s";

    return addArrowPoints;
};

// proteins na backbone na terminus nucleobases
const cartoonAtoms = {
    "C": true, "CA": true, "O": true, "P": true, "OP2": true,
    "O2P": true, "O5'": true, "O3'": true, "C5'": true,
    "C2'": true, "O5*": true, "O3*": true, "C5*": true,
    "C2*": true, "N1": true, "N3": true
};
const purResns = { "DA": true, "DG": true, "A": true, "G": true };
const pyrResns = { "DT": true, "DC": true, "U": true, "C": true, "T": true };
const naResns = { "DA": true, "DG": true, "A": true, "G": true, "DT": true, "DC": true, "U": true, "C": true, "T": true };


export function drawCartoon(group, atomList, gradientrange, quality = 10) {

    let num = quality;
    let div = quality;

    var cartoon, prev, curr, next, currColor, nextColor, thickness, i;
    var backbonePt, orientPt, prevOrientPt, terminalPt, termOrientPt, baseStartPt, baseEndPt;
    var tubeStart, tubeEnd, drawingTube;
    var shapeGeo = new Geometry(true); // for shapes that don't need normals computed
    var geo = new Geometry(true);
    var colors = [];
    var points: any = [];
    var opacity = 1;
    var outline = false;

    var gradients: any = {};
    for (var g in Gradient.builtinGradients) {
        if (Gradient.builtinGradients.hasOwnProperty(g)) {
            //COUNTER INTUITIVE - spectrum reverses direction to gradient to match other tools
            gradients[g] = new Gradient.builtinGradients[g](gradientrange[1], gradientrange[0]);
        }
    }

    var cartoonColor = function (next, cartoon) { //atom and cartoon style object
        if (gradientrange && cartoon.color === 'spectrum') {
            if (cartoon.colorscheme in gradients) {
                return gradients[cartoon.colorscheme].valueToHex(next.resi);
            } else {
                return gradients.sinebow.valueToHex(next.resi);
            }
        }
        else {
            return getColorFromStyle(next, cartoon).getHex();
        }
    };

    for (i = 0; i < num; i++)
        points[i] = [];

    // first determine where beta sheet arrows and alpha helix tubes belong
    var inSheet = false;
    var inHelix = false; //only considering tube styled helices
    var atoms = [];
    for (i in atomList) {
        next = atomList[i];
        if (next.elem === 'C' && next.atom === 'CA') {
            var connected = inConnectedResidues(curr, next);

            // last two residues in a beta sheet become arrowhead
            if (connected && next.ss === "s") {
                inSheet = true;
            } else if (inSheet) {
                if (curr && prev && curr.style.cartoon.arrows && prev.style.cartoon.arrows) {
                    curr.ss = "arrow end";
                    prev.ss = "arrow start";
                }
                inSheet = false;
            }

            // first and last residues in a helix are used to draw tube
            if (connected && (curr.ss === "h" || curr.ss == "tube start") && curr.style.cartoon.tubes) {
                if (!inHelix && curr.ss != "tube start" && next.style.cartoon.tubes) {
                    next.ss = "tube start";
                    inHelix = true;
                }
            } else if (inHelix) {
                if (curr.ss === "tube start") {
                    curr.ss = "tube end"; //only one residue
                } else if (prev && prev.style.cartoon.tubes) {
                    prev.ss = "tube end";
                }
                inHelix = false;
            }
            prev = curr;
            curr = next;
        }
        if (next && next.atom in cartoonAtoms) {
            atoms.push(next);
        }
    }
    if (inHelix && curr.style.cartoon.tubes) {
        curr.ss = "tube end";
        inHelix = false;
    }

    var flushGeom = function (connect) {
        //write out points, update geom,etc
        if (points[0].length > 0) {
            drawStrip(geo, points, colors, div, thickness, opacity, points.style);
        }

        var saved = [], savedc = null;
        if (connect) {
            //recycle last point to first point of next points array
            for (i = 0; i < num; i++) {
                saved[i] = points[i][points[i].length - 1];
            }
            savedc = colors[colors.length - 1];
        }
        points = [];
        for (i = 0; i < num; i++)
            points[i] = [];
        colors = [];

        if (connect) {
            for (i = 0; i < num; i++) {
                points[i].push(saved[i]);
            }
            colors.push(savedc);
        }

        setGeo(group, geo, opacity, outline, true);
        setGeo(group, shapeGeo, opacity, outline, false);
        geo = new Geometry(true);
        shapeGeo = new Geometry(true);
    };

    // then accumulate points
    curr = undefined;
    for (var a = 0; a < atoms.length; a++) {
        next = atoms[a];

        var nextresn = next.resn.trim();
        var inNucleicAcid = nextresn in naResns;
        opacity = 1;
        // determine cartoon style
        cartoon = next.style.cartoon;
        if (curr && curr.style.cartoon)
            opacity = curr.style.cartoon.opacity;
        if (curr && curr.style.cartoon && curr.style.cartoon.outline)
            outline = curr.style.cartoon.outline;

        // create a new geometry when opacity changes
        //this should work fine if opacity is set by chain, but will
        //break if it changes within the chain
        if (curr && curr.style.cartoon && (!next.style.cartoon ||
            curr.style.cartoon.opacity != next.style.cartoon.opacity)) {
            flushGeom(curr.chain == next.chain);
        }

        if (cartoon.style === "trace") // draw cylinders connecting
        // consecutive 'backbone' atoms
        {
            /*
             * "trace" style just draws cylinders between consecutive
             * 'backbone' atoms, such as alpha carbon for polypeptides and
             * phosphorus for DNA.
             */
            if (next.hetflag) {
                ; //ignore non-protein atoms
            } else if (next.elem === 'C' && next.atom === 'CA' ||
                inNucleicAcid && next.atom === "P" ||
                next.atom === 'BB') {
                // determine cylinder color
                nextColor = cartoonColor(next, cartoon);

                // determine cylinder thickness
                if (isNumeric(cartoon.thickness))
                    thickness = cartoon.thickness;
                else
                    thickness = defaultThickness;

                if (inConnectedResidues(curr, next)) {
                    // if both atoms are same color, draw single cylinder
                    if (nextColor == currColor) {
                        var color = CC.color(nextColor);
                        GLDraw.drawCylinder(shapeGeo, curr, next,
                            thickness, color, 2, 2);
                    }

                    else // otherwise draw cylinders for each color
                    // (split down the middle)
                    {
                        var midpoint = new Vector3().addVectors(
                            curr, next).multiplyScalar(0.5);
                        var color1 = CC.color(currColor);
                        var color2 = CC.color(nextColor);
                        GLDraw.drawCylinder(shapeGeo, curr,
                            midpoint, thickness, color1, 2, 0);
                        GLDraw.drawCylinder(shapeGeo, midpoint,
                            next, thickness, color2, 0, 2);
                    } // note that an atom object can be duck-typed as a
                    // Vector3 in this case
                }

                if ((next.clickable === true || next.hoverable) && (next.intersectionShape !== undefined)) {
                    //can click on joints to get alpha carbons
                    var center = new Vector3(next.x, next.y, next.z);
                    next.intersectionShape.sphere.push(new Sphere(center, thickness));
                }

                curr = next;
                currColor = nextColor;
            }
        } else // draw default-style cartoons based on secondary structure
        {
            // draw backbone through these atoms
            if (isAlphaCarbon(next) || inNucleicAcid && (next.atom === "P" || next.atom.indexOf('O5') == 0)) {
                if (drawingTube) {
                    if (next.ss === "tube end") {
                        drawingTube = false;
                        tubeEnd = new Vector3(next.x, next.y, next.z);
                        GLDraw.drawCylinder(shapeGeo, tubeStart,
                            tubeEnd, 2, CC.color(currColor), 1,
                            1);
                        next.ss = "h";

                    }
                    else if (curr.chain != next.chain || curr.ss === "tube end") { //don't span chains no matter what, check for short tubes (less than ideal)
                        drawingTube = false;
                        curr.ss = "h";
                        tubeEnd = new Vector3(curr.x, curr.y, curr.z);
                        GLDraw.drawCylinder(shapeGeo, tubeStart,
                            tubeEnd, 2, CC.color(currColor), 1,
                            1);
                    }
                    else
                        continue; // don't accumulate strand points while
                    // in the middle of drawing a tube
                }

                // end of a chain of connected residues (of same style)
                if (curr && (!inConnectedResidues(curr, next) || curr.ss === "tube start")) {
                    if (curr.ss === "tube start") {
                        drawingTube = true;
                        tubeStart = new Vector3(curr.x, curr.y,
                            curr.z);
                        curr.ss = "h";
                    }

                    if (baseEndPt) // draw the last base if it's a NA chain
                    {
                        if (terminalPt)
                            baseStartPt = new Vector3().addVectors(
                                curr, terminalPt).multiplyScalar(0.5);
                        else
                            baseStartPt = new Vector3(curr.x,
                                curr.y, curr.z);

                        GLDraw.drawCylinder(shapeGeo, baseStartPt,
                            baseEndPt, baseThickness, CC
                                .color(baseEndPt.color), 0, 2);
                        addBackbonePoints(points, num,
                            true, terminalPt, termOrientPt,
                            prevOrientPt, curr, atoms, a);
                        colors.push(nextColor);

                        baseStartPt = null;
                        baseEndPt = null;
                    }

                    // draw accumulated strand points
                    if (points[0].length > 0)
                        drawStrip(geo, points, colors, div, thickness,
                            opacity, points.style);

                    // clear arrays for points and colors
                    points = [];
                    for (i = 0; i < num; i++)
                        points[i] = [];
                    colors = [];
                }

                // reached next residue (potentially the first residue)
                if (curr === undefined || curr.rescode != next.rescode || curr.resi != next.resi) {
                    if (baseEndPt) // draw last NA residue's base
                    {
                        // start the cylinder at the midpoint between
                        // consecutive backbone atoms
                        baseStartPt = new Vector3().addVectors(curr,
                            next).multiplyScalar(0.5);
                        var startFix = baseStartPt.clone().sub(baseEndPt)
                            .multiplyScalar(0.02); // TODO: apply this
                        // as function of
                        // thickness
                        baseStartPt.add(startFix);

                        GLDraw.drawCylinder(shapeGeo, baseStartPt,
                            baseEndPt, baseThickness, CC
                                .color(baseEndPt.color), 0, 2);
                        baseStartPt = null;
                        baseEndPt = null;
                    }

                    // determine color and thickness of the next strand
                    // segment
                    nextColor = cartoonColor(next, cartoon);
                    colors.push(nextColor);
                    if (isNumeric(cartoon.thickness))
                        thickness = cartoon.thickness;
                    else
                        thickness = defaultThickness;

                    curr = next; // advance backbone
                    backbonePt = new Vector3(curr.x, curr.y, curr.z);
                    backbonePt.resi = curr.resi;
                    currColor = nextColor;
                }

                // click handling
                if ((next.clickable === true || next.hoverable === true) &&
                    (next.intersectionShape === undefined || next.intersectionShape.triangle === undefined))
                    next.intersectionShape = {
                        sphere: null,
                        cylinder: [],
                        line: [],
                        triangle: []
                    };
            }
            // atoms used to orient the backbone strand
            else if (isAlphaCarbon(curr) && next.atom === "O" ||
                inNucleicAcid && curr.atom === "P" &&
                (next.atom === "OP2" || next.atom === "O2P") ||
                inNucleicAcid && curr.atom.indexOf("O5") == 0 &&
                next.atom.indexOf("C5") == 0) {
                orientPt = new Vector3(next.x, next.y, next.z);
                orientPt.resi = next.resi;
                if (next.atom === "OP2" || next.atom === "O2P") // for NA 3'
                    // terminus
                    termOrientPt = new Vector3(next.x, next.y,
                        next.z);
            }

            // NA 3' terminus is an edge case, need a vector for most recent
            // O3'
            else if (inNucleicAcid && next.atom.indexOf("O3") == 0) {
                terminalPt = new Vector3(next.x, next.y, next.z);
            }

            // atoms used for drawing the NA base cylinders (diff for
            // purines and pyramidines)
            else if ((next.atom === "N1" && (nextresn in purResns)) ||
                (next.atom === "N3" && (nextresn in pyrResns))) {
                baseEndPt = new Vector3(next.x, next.y, next.z);
                baseEndPt.color = getColorFromStyle(next, cartoon)
                    .getHex();
            }

            // when we have a backbone point and orientation point in the
            // same residue, accumulate strand points
            if (orientPt && backbonePt && orientPt.resi === backbonePt.resi) {
                addBackbonePoints(points, num, true,
                    backbonePt, orientPt, prevOrientPt, curr, atoms,
                    a);
                prevOrientPt = orientPt;
                backbonePt = null;
                orientPt = null;
                colors.push(nextColor);
            }
        }

    }

    if (baseEndPt) // draw last NA base if needed
    {
        if (terminalPt)
            baseStartPt = new Vector3().addVectors(curr, terminalPt)
                .multiplyScalar(0.5);
        else
            baseStartPt = new Vector3(curr.x, curr.y, curr.z);

        GLDraw.drawCylinder(shapeGeo, baseStartPt, baseEndPt, baseThickness,
            CC.color(baseEndPt.color), 0, 2);
        addBackbonePoints(points, num, true, terminalPt,
            termOrientPt, prevOrientPt, curr, atoms, a);
        colors.push(nextColor);
    }

    // for default style, draw the last strand
    flushGeom(false);
};



