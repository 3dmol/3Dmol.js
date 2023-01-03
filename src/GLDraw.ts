import { Vector3 } from "./WebGL/math";
import { Geometry } from "./WebGL";


//define enum values
/**
 * Enum for cylinder cap styles.
 * @readonly
 * @enum 
 * @property NONE
 * @property FLAT
 * @property ROUND
 */
export enum CAP {
    NONE = 0,
    FLAT = 1,
    ROUND = 2
};


/**
 * Lower level utilities for creating WebGL shape geometries.
 * These are not intended for general consumption.
 * @namespace 
  */
export namespace GLDraw {

    // Rotation matrix around z and x axis -
    // according to y basis vector
    // TODO: Try to optimize this (square roots?)
    function getRotationMatrix(dx: number, dy: number, dz: number) {
        var dxy = Math.sqrt(dx * dx + dy * dy);
        var dyz;

        var sinA, cosA, sinB, cosB;

        // about z axis - Phi
        if (dxy < 0.0001) {
            sinA = 0;
            cosA = 1;
        }

        else {
            sinA = -dx / dxy;
            cosA = dy / dxy;
        }

        // recast dy in terms of new axes - z is the same

        dy = -sinA * dx + cosA * dy;
        dyz = Math.sqrt(dy * dy + dz * dz);

        // about new x axis - Theta

        if (dyz < 0.0001) {
            sinB = 0;
            cosB = 1;
        }

        else {
            sinB = dz / dyz;
            cosB = dy / dyz;
        }

        var rot = new Float32Array(9);
        rot[0] = cosA;
        rot[1] = sinA;
        rot[2] = 0;
        rot[3] = -sinA * cosB;
        rot[4] = cosA * cosB;
        rot[5] = sinB;
        rot[6] = sinA * sinB;
        rot[7] = -cosA * sinB;
        rot[8] = cosB;

        return rot;

    };


    // memoize capped cylinder for given radius cylVertexCache
    class CylVertexCache {

        // memoize both rounded and flat caps (hemisphere and circle)
        cache: any = {};


        // Ortho normal vectors for cylinder radius/ sphere cap equator and cones
        // Direction is j basis (0,1,0)
        basisVectors: any;

        constructor() {

            //initialize basisVectors
            let nvecs = [];

            let subdivisions = 4; // including the initial 2, eg. 4 => 16 subintervals
            let N = Math.pow(2, subdivisions);  // eg. 2**4 = 16 subintervals in total
            let i = 2;  // start with 2 subdivisions already done
            let M = Math.pow(2, i); // 4
            let spacing = N / M;  // 16/4 = 4; if there were 5 subdivs, then 32/4 = 8.
            let j: number;

            nvecs[0] = new Vector3(-1, 0, 0);
            nvecs[spacing] = new Vector3(0, 0, 1);
            nvecs[spacing * 2] = new Vector3(1, 0, 0);
            nvecs[spacing * 3] = new Vector3(0, 0, -1);

            for (i = 3; i <= subdivisions; i++) {
                // eg. i=3, we need to add 2**(3-1) = 4 new vecs. Call it M.
                // their spacing is N/M, eg. N=16, M=4, N/M=4; M=8, N/M=2.
                // they start off at half this spacing
                // and are equal to the average of the two vectors on either side
                M = Math.pow(2, (i - 1));
                spacing = N / M;
                for (j = 0; j < (M - 1); j++) {
                    nvecs[spacing / 2 + j * spacing] = nvecs[j * spacing].clone().add(nvecs[(j + 1) * spacing]).normalize();
                }
                // treat the last one specially so it wraps around to zero
                j = M - 1;
                nvecs[spacing / 2 + j * spacing] = nvecs[j * spacing].clone().add(nvecs[0]).normalize();
            }

            this.basisVectors = nvecs;
        };

        getVerticesForRadius(radius, cap, capType) {
            if (typeof (this.cache) !== "undefined" && this.cache[radius] !== undefined)
                if (this.cache[radius][cap + capType] !== undefined)
                    return this.cache[radius][cap + capType];

            var w = this.basisVectors.length;
            var nvecs = [], norms = [];
            var n;

            for (var i = 0; i < w; i++) {
                // bottom
                nvecs.push(this.basisVectors[i].clone().multiplyScalar(radius));
                // top
                nvecs.push(this.basisVectors[i].clone().multiplyScalar(radius));

                // NOTE: this normal is used for constructing sphere caps -
                // cylinder normals taken care of in drawCylinder
                n = this.basisVectors[i].clone().normalize();
                norms.push(n);
                norms.push(n);
            }

            // norms[0]

            var verticesRows = [];

            // Require that heightSegments is even and >= 2
            // Equator points at h/2 (theta = pi/2)
            // (repeated) polar points at 0 and h (theta = 0 and pi)
            var heightSegments = 10, widthSegments = w; // 16 or however many
            // basis vectors for
            // cylinder

            if (heightSegments % 2 !== 0 || !heightSegments) {
                console.error("heightSegments must be even");
                return null;
            }

            var phiStart = 0;
            var phiLength = Math.PI * 2;

            var thetaStart = 0;
            var thetaLength = Math.PI;

            var x: number, y:number;
            var polar = false, equator = false;

            for (y = 0; y <= heightSegments; y++) {

                polar = (y === 0 || y === heightSegments) ? true : false;
                equator = (y === heightSegments / 2) ? true : false;

                var verticesRow = [], toRow = [];

                for (x = 0; x <= widthSegments; x++) {

                    // Two vertices rows for equator pointing to previously
                    // constructed cyl points
                    if (equator) {
                        var xi = (x < widthSegments) ? 2 * x : 0;
                        toRow.push(xi + 1);
                        verticesRow.push(xi);

                        continue;
                    }

                    var u = x / widthSegments;
                    var v = y / heightSegments;

                    // Only push first polar point

                    if (!polar || x === 0) {

                        if (x < widthSegments) {
                            var vertex = new Vector3();
                            vertex.x = -radius *
                                Math.cos(phiStart + u * phiLength) *
                                Math.sin(thetaStart + v * thetaLength);
                            if (cap == 1)
                                vertex.y = 0;
                            else
                                vertex.y = radius * Math.cos(thetaStart + v * thetaLength);

                            vertex.z = radius *
                                Math.sin(phiStart + u * phiLength) *
                                Math.sin(thetaStart + v * thetaLength);

                            if (Math.abs(vertex.x) < 1e-5)
                                vertex.x = 0;
                            if (Math.abs(vertex.y) < 1e-5)
                                vertex.y = 0;
                            if (Math.abs(vertex.z) < 1e-5)
                                vertex.z = 0;

                            if (cap == CAP.FLAT) {
                                n = new Vector3(0, Math.cos(thetaStart + v * thetaLength), 0);
                                n.normalize();
                            }
                            else {
                                n = new Vector3(vertex.x, vertex.y, vertex.z);
                                n.normalize();
                            }

                            nvecs.push(vertex);
                            norms.push(n);

                            verticesRow.push(nvecs.length - 1);
                        }

                        // last point is just the first point for this row
                        else {
                            verticesRow.push(nvecs.length - widthSegments);
                        }

                    }

                    // x > 0; index to already added point
                    else if (polar)
                        verticesRow.push(nvecs.length - 1);

                }

                // extra equator row
                if (equator)
                    verticesRows.push(toRow);

                verticesRows.push(verticesRow);

            }

            var obj = {
                vertices: nvecs,
                normals: norms,
                verticesRows: verticesRows,
                w: widthSegments,
                h: heightSegments
            };

            if (!(radius in this.cache)) this.cache[radius] = {};
            this.cache[radius][cap + capType] = obj;

            return obj;

        }
    };

    var cylVertexCache = new CylVertexCache();


    /** 
     * Create a cylinder 
     * @memberof GLDraw 
     * @param {geometry}
     *            geo
     * @param {Point}
     *            from
     * @param {Point}
     *            to
     * @param {float}
     *            radius
     * @param {Color}
     *            color
     * @param {CAP} fromCap - 0 for none, 1 for flat, 2 for round
     * @param {CAP} toCap = 0 for none, 1 for flat, 2 for round
     *            
     * */
    export function drawCylinder(geo: Geometry, from, to, radius: number, color, fromCap:CAP = 0, toCap:CAP = 0) {
        if (!from || !to)
            return;

        // vertices
        var drawcaps = toCap || fromCap;
        color = color || { r: 0, g: 0, b: 0 };

        var e = getRotationMatrix(to.x-from.x, to.y-from.y, to.z-from.z);
        // get orthonormal vectors from cache
        // TODO: Will have orient with model view matrix according to direction

        var vobj = cylVertexCache.getVerticesForRadius(radius, toCap, "to");
        // w (n) corresponds to the number of orthonormal vectors for cylinder
        // (default 16)
        var n = vobj.w, h = vobj.h;

        // get orthonormal vector
        var n_verts = (drawcaps) ? h * n + 2 : 2 * n;

        var geoGroup = geo.updateGeoGroup(n_verts);

        var vertices = vobj.vertices, normals = vobj.normals, verticesRows = vobj.verticesRows;
        var toRow = verticesRows[h / 2], fromRow = verticesRows[h / 2 + 1];

        var start = geoGroup.vertices;
        var offset, faceoffset;
        var i, x, y, z;

        var vertexArray = geoGroup.vertexArray;
        var normalArray = geoGroup.normalArray;
        var colorArray = geoGroup.colorArray;
        var faceArray = geoGroup.faceArray;
        // add vertices, opposing vertices paired together
        for (i = 0; i < n; ++i) {

            var vi = 2 * i;

            x = e[0] * vertices[vi].x + e[3] * vertices[vi].y + e[6] * vertices[vi].z;
            y = e[1] * vertices[vi].x + e[4] * vertices[vi].y + e[7] * vertices[vi].z;
            z = e[5] * vertices[vi].y + e[8] * vertices[vi].z;

            // var xn = x/radius, yn = y/radius, zn = z/radius;

            offset = 3 * (start + vi);
            faceoffset = geoGroup.faceidx;

            // from
            vertexArray[offset] = x + from.x;
            vertexArray[offset + 1] = y + from.y;
            vertexArray[offset + 2] = z + from.z;
            // to
            vertexArray[offset + 3] = x + to.x;
            vertexArray[offset + 4] = y + to.y;
            vertexArray[offset + 5] = z + to.z;

            // normals
            normalArray[offset] = x;
            normalArray[offset + 3] = x;
            normalArray[offset + 1] = y;
            normalArray[offset + 4] = y;
            normalArray[offset + 2] = z;
            normalArray[offset + 5] = z;

            // colors
            colorArray[offset] = color.r;
            colorArray[offset + 3] = color.r;
            colorArray[offset + 1] = color.g;
            colorArray[offset + 4] = color.g;
            colorArray[offset + 2] = color.b;
            colorArray[offset + 5] = color.b;

            // faces
            // 0 - 2 - 1
            faceArray[faceoffset] = fromRow[i] + start;
            faceArray[faceoffset + 1] = fromRow[i + 1] + start;
            faceArray[faceoffset + 2] = toRow[i] + start;
            // 1 - 2 - 3
            faceArray[faceoffset + 3] = toRow[i] + start;
            faceArray[faceoffset + 4] = fromRow[i + 1] + start;
            faceArray[faceoffset + 5] = toRow[i + 1] + start;

            geoGroup.faceidx += 6;

        }

        // SPHERE CAPS
        if (drawcaps) {
            // h - sphere rows, verticesRows.length - 2

            var ystart = (toCap) ? 0 : h / 2;
            var yend = (fromCap) ? h + 1 : h / 2 + 1;
            var v1, v2, v3, v4, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4, nz1, nz2, nz3, nz4, v1offset, v2offset, v3offset, v4offset;

            for (y = ystart; y < yend; y++) {
                if (y === h / 2)
                    continue;
                // n number of points for each level (verticesRows[i].length -
                // 1)
                var cap = (y <= h / 2) ? to : from;
                var toObj = cylVertexCache.getVerticesForRadius(radius, toCap, "to");
                var fromObj = cylVertexCache.getVerticesForRadius(radius, fromCap, "from");
                if (cap === to) {
                    vertices = toObj.vertices;
                    normals = toObj.normals;
                    verticesRows = toObj.verticesRows;
                } else if (cap == from) {
                    vertices = fromObj.vertices;
                    normals = fromObj.normals;
                    verticesRows = fromObj.verticesRows;
                }
                for (x = 0; x < n; x++) {

                    faceoffset = geoGroup.faceidx;

                    v1 = verticesRows[y][x + 1];
                    v1offset = (v1 + start) * 3;
                    v2 = verticesRows[y][x];
                    v2offset = (v2 + start) * 3;
                    v3 = verticesRows[y + 1][x];
                    v3offset = (v3 + start) * 3;
                    v4 = verticesRows[y + 1][x + 1];
                    v4offset = (v4 + start) * 3;

                    // rotate sphere vectors
                    x1 = e[0] * vertices[v1].x + e[3] * vertices[v1].y + e[6] * vertices[v1].z;
                    x2 = e[0] * vertices[v2].x + e[3] * vertices[v2].y + e[6] * vertices[v2].z;
                    x3 = e[0] * vertices[v3].x + e[3] * vertices[v3].y + e[6] * vertices[v3].z;
                    x4 = e[0] * vertices[v4].x + e[3] * vertices[v4].y + e[6] * vertices[v4].z;

                    y1 = e[1] * vertices[v1].x + e[4] * vertices[v1].y + e[7] * vertices[v1].z;
                    y2 = e[1] * vertices[v2].x + e[4] * vertices[v2].y + e[7] * vertices[v2].z;
                    y3 = e[1] * vertices[v3].x + e[4] * vertices[v3].y + e[7] * vertices[v3].z;
                    y4 = e[1] * vertices[v4].x + e[4] * vertices[v4].y + e[7] * vertices[v4].z;

                    z1 = e[5] * vertices[v1].y + e[8] * vertices[v1].z;
                    z2 = e[5] * vertices[v2].y + e[8] * vertices[v2].z;
                    z3 = e[5] * vertices[v3].y + e[8] * vertices[v3].z;
                    z4 = e[5] * vertices[v4].y + e[8] * vertices[v4].z;

                    vertexArray[v1offset] = x1 + cap.x;
                    vertexArray[v2offset] = x2 + cap.x;
                    vertexArray[v3offset] = x3 + cap.x;
                    vertexArray[v4offset] = x4 + cap.x;

                    vertexArray[v1offset + 1] = y1 + cap.y;
                    vertexArray[v2offset + 1] = y2 + cap.y;
                    vertexArray[v3offset + 1] = y3 + cap.y;
                    vertexArray[v4offset + 1] = y4 + cap.y;

                    vertexArray[v1offset + 2] = z1 + cap.z;
                    vertexArray[v2offset + 2] = z2 + cap.z;
                    vertexArray[v3offset + 2] = z3 + cap.z;
                    vertexArray[v4offset + 2] = z4 + cap.z;

                    colorArray[v1offset] = color.r;
                    colorArray[v2offset] = color.r;
                    colorArray[v3offset] = color.r;
                    colorArray[v4offset] = color.r;

                    colorArray[v1offset + 1] = color.g;
                    colorArray[v2offset + 1] = color.g;
                    colorArray[v3offset + 1] = color.g;
                    colorArray[v4offset + 1] = color.g;

                    colorArray[v1offset + 2] = color.b;
                    colorArray[v2offset + 2] = color.b;
                    colorArray[v3offset + 2] = color.b;
                    colorArray[v4offset + 2] = color.b;

                    nx1 = e[0] * normals[v1].x + e[3] * normals[v1].y + e[6] * normals[v1].z;
                    nx2 = e[0] * normals[v2].x + e[3] * normals[v2].y + e[6] * normals[v2].z;
                    nx3 = e[0] * normals[v3].x + e[3] * normals[v3].y + e[6] * normals[v3].z;
                    nx4 = e[0] * normals[v4].x + e[3] * normals[v4].y + e[6] * normals[v4].z;

                    ny1 = e[1] * normals[v1].x + e[4] * normals[v1].y + e[7] * normals[v1].z;
                    ny2 = e[1] * normals[v2].x + e[4] * normals[v2].y + e[7] * normals[v2].z;
                    ny3 = e[1] * normals[v3].x + e[4] * normals[v3].y + e[7] * normals[v3].z;
                    ny4 = e[1] * normals[v4].x + e[4] * normals[v4].y + e[7] * normals[v4].z;

                    nz1 = e[5] * normals[v1].y + e[8] * normals[v1].z;
                    nz2 = e[5] * normals[v2].y + e[8] * normals[v2].z;
                    nz3 = e[5] * normals[v3].y + e[8] * normals[v3].z;
                    nz4 = e[5] * normals[v4].y + e[8] * normals[v4].z;

                    // if (Math.abs(vobj.sphereVertices[v1].y) === radius) {

                    if (y === 0) {//to center circle
                        // face = [v1, v3, v4];
                        // norm = [n1, n3, n4];

                        normalArray[v1offset] = nx1;
                        normalArray[v3offset] = nx3;
                        normalArray[v4offset] = nx4;
                        normalArray[v1offset + 1] = ny1;
                        normalArray[v3offset + 1] = ny3;
                        normalArray[v4offset + 1] = ny4;
                        normalArray[v1offset + 2] = nz1;
                        normalArray[v3offset + 2] = nz3;
                        normalArray[v4offset + 2] = nz4;

                        faceArray[faceoffset] = v1 + start;
                        faceArray[faceoffset + 1] = v3 + start;
                        faceArray[faceoffset + 2] = v4 + start;

                        geoGroup.faceidx += 3;

                    }

                    // else if (Math.abs(vobj.sphereVertices[v3].y) === radius)
                    // {
                    else if (y === yend - 1) {//from end center circle
                        // face = [v1, v2, v3];
                        // norm = [n1, n2, n3];

                        normalArray[v1offset] = nx1;
                        normalArray[v2offset] = nx2;
                        normalArray[v3offset] = nx3;
                        normalArray[v1offset + 1] = ny1;
                        normalArray[v2offset + 1] = ny2;
                        normalArray[v3offset + 1] = ny3;
                        normalArray[v1offset + 2] = nz1;
                        normalArray[v2offset + 2] = nz2;
                        normalArray[v3offset + 2] = nz3;

                        faceArray[faceoffset] = v1 + start;
                        faceArray[faceoffset + 1] = v2 + start;
                        faceArray[faceoffset + 2] = v3 + start;

                        geoGroup.faceidx += 3;

                    }

                    else { // the rest of the circles
                        // face = [v1, v2, v3, v4];
                        // norm = [n1, n2, n3, n4];

                        normalArray[v1offset] = nx1;
                        normalArray[v2offset] = nx2;
                        normalArray[v4offset] = nx4;
                        normalArray[v1offset + 1] = ny1;
                        normalArray[v2offset + 1] = ny2;
                        normalArray[v4offset + 1] = ny4;
                        normalArray[v1offset + 2] = nz1;
                        normalArray[v2offset + 2] = nz2;
                        normalArray[v4offset + 2] = nz4;

                        normalArray[v2offset] = nx2;
                        normalArray[v3offset] = nx3;
                        normalArray[v4offset] = nx4;
                        normalArray[v2offset + 1] = ny2;
                        normalArray[v3offset + 1] = ny3;
                        normalArray[v4offset + 1] = ny4;
                        normalArray[v2offset + 2] = nz2;
                        normalArray[v3offset + 2] = nz3;
                        normalArray[v4offset + 2] = nz4;

                        faceArray[faceoffset] = v1 + start;
                        faceArray[faceoffset + 1] = v2 + start;
                        faceArray[faceoffset + 2] = v4 + start;

                        faceArray[faceoffset + 3] = v2 + start;
                        faceArray[faceoffset + 4] = v3 + start;
                        faceArray[faceoffset + 5] = v4 + start;

                        geoGroup.faceidx += 6;
                    }

                }
            }

        }

        geoGroup.vertices += n_verts;
    };


    /** Create a cone 
     * @memberof GLDraw
     * @param {geometry}
     *            geo
     * @param {Point}
     *            from
     * @param {Point}
     *            to
     * @param {float}
     *            radius
     * @param {Color}
     *            color
     *            */
    export function drawCone (geo: Geometry, from, to, radius: number, color?) {
        if (!from || !to)
            return;

        //TODO: check if from and to do not contain x,y,z and if  so generate a center based on the passed selections

        color = color || { r: 0, g: 0, b: 0 };

        let ndir = new Vector3(to.x-from.x, to.y-from.y, to.z-from.z);
        var e = getRotationMatrix(ndir.x, ndir.y, ndir.z);
        ndir = ndir.normalize();

        // n vertices around bottom plust the two points
        var n = cylVertexCache.basisVectors.length;
        var basis = cylVertexCache.basisVectors;
        var n_verts = n + 2;

        //setup geo structures
        var geoGroup = geo.updateGeoGroup(n_verts);
        var start = geoGroup.vertices;
        var offset, faceoffset;
        var i, x, y, z;
        var vertexArray = geoGroup.vertexArray;
        var normalArray = geoGroup.normalArray;
        var colorArray = geoGroup.colorArray;
        var faceArray = geoGroup.faceArray;

        offset = start * 3;
        //base point first vertex
        vertexArray[offset] = from.x;
        vertexArray[offset + 1] = from.y;
        vertexArray[offset + 2] = from.z;
        normalArray[offset] = -ndir.x;
        normalArray[offset + 1] = -ndir.y;
        normalArray[offset + 2] = -ndir.z;
        colorArray[offset] = color.r;
        colorArray[offset + 1] = color.g;
        colorArray[offset + 2] = color.b;

        //second vertex top
        vertexArray[offset + 3] = to.x;
        vertexArray[offset + 4] = to.y;
        vertexArray[offset + 5] = to.z;

        normalArray[offset + 3] = ndir.x;
        normalArray[offset + 4] = ndir.y;
        normalArray[offset + 5] = ndir.z;
        colorArray[offset + 3] = color.r;
        colorArray[offset + 4] = color.g;
        colorArray[offset + 5] = color.b;

        offset += 6;

        // add circle vertices
        for (i = 0; i < n; ++i) {
            var vec = basis[i].clone();
            vec.multiplyScalar(radius);
            x = e[0] * vec.x + e[3] * vec.y + e[6] * vec.z;
            y = e[1] * vec.x + e[4] * vec.y + e[7] * vec.z;
            z = e[5] * vec.y + e[8] * vec.z;

            // from
            vertexArray[offset] = x + from.x;
            vertexArray[offset + 1] = y + from.y;
            vertexArray[offset + 2] = z + from.z;

            // normals
            normalArray[offset] = x;
            normalArray[offset + 1] = y;
            normalArray[offset + 2] = z;

            // colors
            colorArray[offset] = color.r;
            colorArray[offset + 1] = color.g;
            colorArray[offset + 2] = color.b;

            offset += 3;

        }
        geoGroup.vertices += (n + 2);
        //faces
        faceoffset = geoGroup.faceidx;
        for (i = 0; i < n; i++) {
            //two neighboring circle vertices
            var v1 = start + 2 + i;
            var v2 = start + 2 + ((i + 1) % n);

            faceArray[faceoffset] = v1;
            faceArray[faceoffset + 1] = v2;
            faceArray[faceoffset + 2] = start;
            faceoffset += 3;
            faceArray[faceoffset] = v1;
            faceArray[faceoffset + 1] = v2;
            faceArray[faceoffset + 2] = start + 1;
            faceoffset += 3;
        }
        geoGroup.faceidx += 6 * n;
    };


    // Sphere component sphereVertexCache
    class  SphereVertexCache {
        private cache = new Map<number, Map<number, any>>(); //sphereQuality then radius
        constructor() {}

        getVerticesForRadius(radius, sphereQuality) {
            sphereQuality = sphereQuality || 2;

            if (!this.cache.has(sphereQuality))  {
                this.cache.set(sphereQuality, new Map<number,any>());
            }
            let radiusCache = this.cache.get(sphereQuality);
            if (radiusCache.has(radius))
                return radiusCache.get(radius);

            var obj = {
                vertices: [],
                verticesRows: [],
                normals: []
            };
            // scale quality with radius heuristically
            var widthSegments = 16 * sphereQuality;
            var heightSegments = 10 * sphereQuality;
            if (radius < 1) {
                widthSegments = 10 * sphereQuality;
                heightSegments = 8 * sphereQuality;
            }

            var phiStart = 0;
            var phiLength = Math.PI * 2;

            var thetaStart = 0;
            var thetaLength = Math.PI;

            var x, y;

            for (y = 0; y <= heightSegments; y++) {

                let verticesRow = [];
                for (x = 0; x <= widthSegments; x++) {

                    let u = x / widthSegments;
                    let v = y / heightSegments;

                    let vx = -radius * Math.cos(phiStart + u * phiLength) *
                        Math.sin(thetaStart + v * thetaLength);
                    let vy = radius * Math.cos(thetaStart + v * thetaLength);
                    let vz = radius * Math.sin(phiStart + u * phiLength) *
                        Math.sin(thetaStart + v * thetaLength);

                    var n = new Vector3(vx, vy, vz);
                    n.normalize();

                    obj.vertices.push({x: vx, y: vy, z: vz});
                    obj.normals.push(n);

                    verticesRow.push(obj.vertices.length - 1);

                }

                obj.verticesRows.push(verticesRow);

            }

            radiusCache.set(radius, obj);
            return obj;
        }

    };
    var sphereVertexCache = new SphereVertexCache();

    /** Create a sphere.
     * @memberof GLDraw
     * @param {geometry}
     *            geo
     * @param {Point}
     *            pos
     * @param {float}
     *            radius
     * @param {Color}
     *            color
     * @param {number} quality of sphere (default 2, higher increases number of triangles)
     */
    export function drawSphere(geo:Geometry, pos, radius, color, sphereQuality?) {

        var vobj = sphereVertexCache.getVerticesForRadius(radius, sphereQuality);

        var vertices = vobj.vertices;
        var normals = vobj.normals;

        var geoGroup = geo.updateGeoGroup(vertices.length);

        var start = geoGroup.vertices;
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;
        var faceArray = geoGroup.faceArray;
        var lineArray = geoGroup.lineArray;
        var normalArray = geoGroup.normalArray;

        for (let i = 0, il = vertices.length; i < il; ++i) {
            let offset = 3 * (start + i);
            let v = vertices[i];

            vertexArray[offset] = (v.x + pos.x);
            vertexArray[offset + 1] = (v.y + pos.y);
            vertexArray[offset + 2] = (v.z + pos.z);

            colorArray[offset] = color.r;
            colorArray[offset + 1] = color.g;
            colorArray[offset + 2] = color.b;

        }

        geoGroup.vertices += vertices.length;

        let verticesRows = vobj.verticesRows;
        let h = verticesRows.length - 1;

        for (let y = 0; y < h; y++) {
            let w = verticesRows[y].length - 1;
            for (let x = 0; x < w; x++) {

                let faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;

                let v1 = verticesRows[y][x + 1] + start, v1offset = v1 * 3;
                let v2 = verticesRows[y][x] + start, v2offset = v2 * 3;
                let v3 = verticesRows[y + 1][x] + start, v3offset = v3 * 3;
                let v4 = verticesRows[y + 1][x + 1] + start, v4offset = v4 * 3;

                let n1 = normals[v1 - start];
                let n2 = normals[v2 - start];
                let n3 = normals[v3 - start];
                let n4 = normals[v4 - start];

                if (Math.abs(vertices[v1 - start].y) === radius) {
                    // face = [v1, v3, v4];
                    // norm = [n1, n3, n4];

                    normalArray[v1offset] = n1.x;
                    normalArray[v3offset] = n3.x;
                    normalArray[v4offset] = n4.x;
                    normalArray[v1offset + 1] = n1.y;
                    normalArray[v3offset + 1] = n3.y;
                    normalArray[v4offset + 1] = n4.y;
                    normalArray[v1offset + 2] = n1.z;
                    normalArray[v3offset + 2] = n3.z;
                    normalArray[v4offset + 2] = n4.z;

                    faceArray[faceoffset] = v1;
                    faceArray[faceoffset + 1] = v3;
                    faceArray[faceoffset + 2] = v4;

                    lineArray[lineoffset] = v1;
                    lineArray[lineoffset + 1] = v3;
                    lineArray[lineoffset + 2] = v1;
                    lineArray[lineoffset + 3] = v4;
                    lineArray[lineoffset + 4] = v3;
                    lineArray[lineoffset + 5] = v4;

                    geoGroup.faceidx += 3;
                    geoGroup.lineidx += 6;

                } else if (Math.abs(vertices[v3 - start].y) === radius) {
                    // face = [v1, v2, v3];
                    // norm = [n1, n2, n3];

                    normalArray[v1offset] = n1.x;
                    normalArray[v2offset] = n2.x;
                    normalArray[v3offset] = n3.x;
                    normalArray[v1offset + 1] = n1.y;
                    normalArray[v2offset + 1] = n2.y;
                    normalArray[v3offset + 1] = n3.y;
                    normalArray[v1offset + 2] = n1.z;
                    normalArray[v2offset + 2] = n2.z;
                    normalArray[v3offset + 2] = n3.z;

                    faceArray[faceoffset] = v1;
                    faceArray[faceoffset + 1] = v2;
                    faceArray[faceoffset + 2] = v3;

                    lineArray[lineoffset] = v1;
                    lineArray[lineoffset + 1] = v2;
                    lineArray[lineoffset + 2] = v1;
                    lineArray[lineoffset + 3] = v3;
                    lineArray[lineoffset + 4] = v2;
                    lineArray[lineoffset + 5] = v3;

                    geoGroup.faceidx += 3;
                    geoGroup.lineidx += 6;

                } else {
                    // face = [v1, v2, v3, v4];
                    // norm = [n1, n2, n3, n4];

                    normalArray[v1offset] = n1.x;
                    normalArray[v2offset] = n2.x;
                    normalArray[v4offset] = n4.x;
                    normalArray[v1offset + 1] = n1.y;
                    normalArray[v2offset + 1] = n2.y;
                    normalArray[v4offset + 1] = n4.y;
                    normalArray[v1offset + 2] = n1.z;
                    normalArray[v2offset + 2] = n2.z;
                    normalArray[v4offset + 2] = n4.z;

                    normalArray[v2offset] = n2.x;
                    normalArray[v3offset] = n3.x;
                    normalArray[v4offset] = n4.x;
                    normalArray[v2offset + 1] = n2.y;
                    normalArray[v3offset + 1] = n3.y;
                    normalArray[v4offset + 1] = n4.y;
                    normalArray[v2offset + 2] = n2.z;
                    normalArray[v3offset + 2] = n3.z;
                    normalArray[v4offset + 2] = n4.z;

                    faceArray[faceoffset] = v1;
                    faceArray[faceoffset + 1] = v2;
                    faceArray[faceoffset + 2] = v4;

                    faceArray[faceoffset + 3] = v2;
                    faceArray[faceoffset + 4] = v3;
                    faceArray[faceoffset + 5] = v4;

                    lineArray[lineoffset] = v1;
                    lineArray[lineoffset + 1] = v2;
                    lineArray[lineoffset + 2] = v1;
                    lineArray[lineoffset + 3] = v4;

                    lineArray[lineoffset + 4] = v2;
                    lineArray[lineoffset + 5] = v3;
                    lineArray[lineoffset + 6] = v3;
                    lineArray[lineoffset + 7] = v4;

                    geoGroup.faceidx += 6;
                    geoGroup.lineidx += 8;

                }

            }
        }

    };

}
