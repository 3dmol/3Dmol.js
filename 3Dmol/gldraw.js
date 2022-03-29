//

// const $3Dmol = $3Dmol || {};

// define enum values
/**
 * Enum for cylinder cap styles.
 * @readonly
 * @enum {number} $3Dmol.CAP
 */
$3Dmol.CAP = {
    NONE : 0,
    FLAT : 1,
    ROUND : 2
};


/**
 * Lower level utilities for creating WebGL shape geometries.
 * These are not intended for general consumption.
 * @namespace $3Dmol.GLDraw
  */
$3Dmol.GLDraw = (function() {

    const draw = {}; // object for exporting functions

    // Rotation matrix around z and x axis -
    // according to y basis vector
    // TODO: Try to optimize this (square roots?)
    const getRotationMatrix = function() {

        const d = new $3Dmol.Vector3();
        // var rot = new (9);

        return function(dir) {

            d.set(dir[0], dir[1], dir[2]);

            const dx = d.x; let dy = d.y; const dz = d.z;

            const dxy = Math.sqrt(dx * dx + dy * dy);
            let dyz;

            let sinA; let cosA; let sinB; let cosB;

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

            const rot = new Float32Array(9);
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

    }();
    
    // Ortho normal vectors for cylinder radius/ sphere cap equator and cones
    // Direction is j basis (0,1,0)
    const basisVectors = function() {

        const nvecs = [];

        const subdivisions = 4; // including the initial 2, eg. 4 => 16 subintervals
        const N = 2**subdivisions;  // eg. 2**4 = 16 subintervals in total
        let i = 2;  // start with 2 subdivisions already done
        let M = 2**i; // 4
        let spacing = N/M;  // 16/4 = 4; if there were 5 subdivs, then 32/4 = 8.
        let j;

        nvecs[0] = new $3Dmol.Vector3(-1, 0, 0);
        nvecs[spacing] = new $3Dmol.Vector3(0, 0, 1);
        nvecs[spacing*2] = new $3Dmol.Vector3(1, 0, 0);
        nvecs[spacing*3] = new $3Dmol.Vector3(0, 0, -1);

        for ( i = 3; i <= subdivisions; i ++ ) {
            // eg. i=3, we need to add 2**(3-1) = 4 new vecs. Call it M.
            // their spacing is N/M, eg. N=16, M=4, N/M=4; M=8, N/M=2.
            // they start off at half this spacing
            // and are equal to the average of the two vectors on either side
            M = 2**(i-1);
            spacing = N/M;
            for ( j = 0; j < (M-1); j ++ ) {
                nvecs[spacing/2 + j*spacing] = nvecs[j*spacing].clone().add(nvecs[(j+1)*spacing]).normalize();
            }
            // treat the last one specially so it wraps around to zero
            j = M - 1;
            nvecs[spacing/2 + j*spacing] = nvecs[j*spacing].clone().add(nvecs[0]).normalize();
        }

        /*
         * nvecs[0] = new $3Dmol.Vector3(-1,0,0); nvecs[1] = new
         * $3Dmol.Vector3(0,0,1); nvecs[2] = new $3Dmol.Vector3(1,0,0);
         * nvecs[3] = new $3Dmol.Vector3(0,0,-1);
         */
        return nvecs;

    }();

    // memoize capped cylinder for given radius
    const cylVertexCache = {

        // memoize both rounded and flat caps (hemisphere and circle)
        cache :{} ,

        getVerticesForRadius(radius, cap, capType) {
            if(typeof(this.chache) !== "undefined" && this.cache[radius] !== undefined)
                if(this.cache[radius][cap+capType] !== undefined)                                                      
                    return this.cache[radius][cap+capType];

            const w = basisVectors.length;
            const nvecs = []; const norms = [];
            let n;

            for (let i = 0; i < w; i++) {
                // bottom
                nvecs.push(basisVectors[i].clone().multiplyScalar(radius));
                // top
                nvecs.push(basisVectors[i].clone().multiplyScalar(radius));

                // NOTE: this normal is used for constructing sphere caps -
                // cylinder normals taken care of in drawCylinder
                n = basisVectors[i].clone().normalize();
                norms.push(n);
                norms.push(n);
            }

            // norms[0]

            const verticesRows = [];

            // Require that heightSegments is even and >= 2
            // Equator points at h/2 (theta = pi/2)
            // (repeated) polar points at 0 and h (theta = 0 and pi)
            const heightSegments = 10; const widthSegments = w; // 16 or however many
                                                        // basis vectors for
                                                        // cylinder

            if (heightSegments % 2 !== 0 || !heightSegments) {
                console.error("heightSegments must be even");

                return null;
            }

            const phiStart = 0;
            const phiLength = Math.PI * 2;

            const thetaStart = 0;
            const thetaLength = Math.PI;

            let x; let y;
            let polar = false; let equator = false;

            for (y = 0; y <= heightSegments; y++) {

                polar = !!((y === 0 || y === heightSegments));
                equator = (y === heightSegments / 2);

                const verticesRow = []; const toRow = [];

                for (x = 0; x <= widthSegments; x++) {

                    // Two vertices rows for equator pointing to previously
                    // constructed cyl points
                    if (equator) {
                        const xi = (x < widthSegments) ? 2 * x : 0;
                        toRow.push(xi + 1);
                        verticesRow.push(xi);

                        continue;
                    }

                    const u = x / widthSegments;
                    const v = y / heightSegments;

                    // Only push first polar point

                    if (!polar || x === 0) {

                        if (x < widthSegments) {
                            const vertex = new $3Dmol.Vector3();
                            vertex.x = -radius *
                                    Math.cos(phiStart + u * phiLength) *
                                    Math.sin(thetaStart + v * thetaLength);
                            if(cap===1)
                                vertex.y=0;
                            else
                                vertex.y=radius * Math.cos(thetaStart + v * thetaLength);

                            vertex.z = radius *
                                    Math.sin(phiStart + u * phiLength) *
                                    Math.sin(thetaStart + v * thetaLength);

                            if (Math.abs(vertex.x) < 1e-5)
                                vertex.x = 0;
                            if (Math.abs(vertex.y) < 1e-5)
                                vertex.y = 0;
                            if (Math.abs(vertex.z) < 1e-5)
                                vertex.z = 0;

                            if (cap === $3Dmol.CAP.FLAT) {
                                n = new $3Dmol.Vector3(0, Math.cos(thetaStart + v * thetaLength), 0);
                                n.normalize();
                            }
                            else {
                                n = new $3Dmol.Vector3(vertex.x, vertex.y, vertex.z);
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

            const obj = {
                vertices : nvecs,
                normals : norms,
                verticesRows,
                w : widthSegments,
                h : heightSegments
            };
            
            if(!(radius in this.cache)) this.cache[radius]={};
            this.cache[radius][cap+capType] = obj;

            return obj;

        }
    };
    

    // creates a cylinder
    let drawnC = 0;
    
    /** Create a cylinder 
     * @function $3Dmol.GLDraw.drawCylinder
     * @param {geometry}
     *            geo
     * @param {Point}
     *            from
     * @param {Point}
     *            to
     * @param {float}
     *            radius
     * @param {$3Dmol.Color}
     *            color
     * @param {$3Dmol.CAP} fromCap - 0 for none, 1 for flat, 2 for round
     * @param {$3Dmol.CAP} toCap = 0 for none, 1 for flat, 2 for round
     *            
     * */
    draw.drawCylinder = function(geo, from, to, radius, color, fromCap, toCap) {
        if (!from || !to)
            return;
        drawnC+=1;
        // vertices
        const drawcaps = toCap || fromCap;
        color = color || {r:0, g:0, b:0};

        /** @type {Array.<number>} */
        const dir = [ to.x, to.y, to.z ];
        dir[0] -= from.x;
        dir[1] -= from.y;
        dir[2] -= from.z;

        const e = getRotationMatrix(dir);
        // get orthonormal vectors from cache
        // TODO: Will have orient with model view matrix according to direction

        const vobj = cylVertexCache.getVerticesForRadius(radius, toCap, "to");
        // w (n) corresponds to the number of orthonormal vectors for cylinder
        // (default 16)
        const n = vobj.w; const {h} = vobj;
        
        // get orthonormal vector
        const nVerts = (drawcaps) ? h * n + 2 : 2 * n;

        const geoGroup = geo.updateGeoGroup(nVerts);

        let {vertices} = vobj; let {normals} = vobj; let {verticesRows} = vobj;
        const toRow = verticesRows[h / 2]; const fromRow = verticesRows[h / 2 + 1];

        const start = geoGroup.vertices;
        let offset; let faceoffset;
        let i; let x; let y; let z;

        const {vertexArray} = geoGroup;
        const {normalArray} = geoGroup;
        const {colorArray} = geoGroup;
        const {faceArray} = geoGroup;
        // add vertices, opposing vertices paired together
        for (i = 0; i < n; ++i) {

            const vi = 2 * i;

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
            
            const ystart = (toCap) ? 0 : h / 2;
            const yend = (fromCap) ? h + 1 : h / 2 + 1;
            let v1; let v2; let v3; let v4; let x1; let x2; let x3; let x4; let y1; let y2; let y3; let y4; let z1; let z2; let z3; let z4; let nx1; let nx2; let nx3; let nx4; let ny1; let ny2; let ny3; let ny4; let nz1; let nz2; let nz3; let nz4; let v1offset; let v2offset; let v3offset; let v4offset;

            for (y = ystart; y < yend; y++) {
                if (y === h / 2)
                    continue;
                // n number of points for each level (verticesRows[i].length -
                // 1)
                const cap = (y <= h / 2) ? to : from;
                const toObj = cylVertexCache.getVerticesForRadius(radius, toCap, "to");
                const fromObj = cylVertexCache.getVerticesForRadius(radius, fromCap, "from");
                if(cap===to){
                    vertices = toObj.vertices;
                    normals = toObj.normals;
                    verticesRows = toObj.verticesRows;
                }else if(cap===from){
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

                    if (y === 0) {// to center circle
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
                    else if (y === yend - 1) {// from end center circle
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

        geoGroup.vertices += nVerts;
    };


    /** Create a cone 
     * @function $3Dmol.GLDraw.drawCone
     * @param {geometry}
     *            geo
     * @param {Point}
     *            from
     * @param {Point}
     *            to
     * @param {float}
     *            radius
     * @param {$3Dmol.Color}
     *            color
     *            */
    draw.drawCone = function(geo, from, to, radius, color) {
        if (!from || !to)
            return;

        // check if from and to do not contain x,y,z and if  so generate a center based on the passed selections

        color = color || {r:0, g:0, b:0};

        const dir =[to.x, to.y, to.z ];        
        dir.x -= from.x;
        dir.y -= from.y;
        dir.z -= from.z;

        const e = getRotationMatrix(dir);


        // n vertices around bottom plust the two points
        const n = basisVectors.length;
        const basis = basisVectors;
        const nVerts =  n + 2;

        
        // setup geo structures
        const geoGroup = geo.updateGeoGroup(nVerts);
        const start = geoGroup.vertices;    
        let offset; let faceoffset;
        let i; let x; let y; let z;
        const {vertexArray} = geoGroup;
        const {normalArray} = geoGroup;
        const {colorArray} = geoGroup;
        const {faceArray} = geoGroup;
        
        offset = start*3;
        const ndir = new $3Dmol.Vector3(dir[0],dir[1],dir[2]).normalize();
        // base point first vertex
        vertexArray[offset] = from.x;
        vertexArray[offset+1] = from.y;
        vertexArray[offset+2] = from.z;
        normalArray[offset] = -ndir.x;
        normalArray[offset + 1] = -ndir.y;
        normalArray[offset + 2] = -ndir.z;
        colorArray[offset] = color.r;
        colorArray[offset + 1] = color.g;
        colorArray[offset + 2] = color.b;
        
        // second vertex top
        vertexArray[offset+3] = to.x;
        vertexArray[offset+4] = to.y;
        vertexArray[offset+5] = to.z;
        
        normalArray[offset+3] = ndir.x;
        normalArray[offset+4] = ndir.y;
        normalArray[offset+5] = ndir.z;
        colorArray[offset+3] = color.r;
        colorArray[offset + 4] = color.g;
        colorArray[offset + 5] = color.b;
        
        offset += 6;
        
        // add circle vertices
        for (i = 0; i < n; ++i) {
            const vec = basis[i].clone();
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
        geoGroup.vertices += (n+2);
        // faces
        faceoffset = geoGroup.faceidx;
        for( i = 0; i < n; i++) {
            // two neighboring circle vertices
            const v1 = start+2+i;
            const v2 = start+2+ ((i+1)%n);
            
            faceArray[faceoffset] = v1;
            faceArray[faceoffset+1] = v2;
            faceArray[faceoffset+2] = start;
            faceoffset += 3;
            faceArray[faceoffset] = v1;
            faceArray[faceoffset+1] = v2;
            faceArray[faceoffset+2] = start+1;
            faceoffset += 3;
        }
        geoGroup.faceidx += 6*n;
    };

    
    // Sphere component
    const sphereVertexCache = {
        cache : {},
        getVerticesForRadius(radius) {

            if (typeof (this.cache[radius]) !== "undefined")
                return this.cache[radius];

            const obj = {
                vertices : [],
                verticesRows : [],
                normals : []
            };
            // scale quality with radius heuristically
            const sphereQuality = 1;
            let widthSegments = 16 * sphereQuality;
            let heightSegments = 10 * sphereQuality;
            if (radius < 1) {
                widthSegments = 10 * sphereQuality;
                heightSegments = 8 * sphereQuality;
            }

            const phiStart = 0;
            const phiLength = Math.PI * 2;

            const thetaStart = 0;
            const thetaLength = Math.PI;

            let x; let y;

            for (y = 0; y <= heightSegments; y++) {

                const verticesRow = [];
                for (x = 0; x <= widthSegments; x++) {

                    const u = x / widthSegments;
                    const v = y / heightSegments;

                    const vertex = {};
                    vertex.x = -radius * Math.cos(phiStart + u * phiLength) *
                            Math.sin(thetaStart + v * thetaLength);
                    vertex.y = radius * Math.cos(thetaStart + v * thetaLength);
                    vertex.z = radius * Math.sin(phiStart + u * phiLength) *
                            Math.sin(thetaStart + v * thetaLength);

                    const n = new $3Dmol.Vector3(vertex.x, vertex.y, vertex.z);
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

    /** Create a sphere.
     * @function $3Dmol.GLDraw.drawSphere
     * @param {geometry}
     *            geo
     * @param {Point}
     *            pos
     * @param {float}
     *            radius
     * @param {$3Dmol.Color}
     *            color
     */
    draw.drawSphere = function(geo, pos, radius, color) {

        let x; let y;
        const vobj = sphereVertexCache.getVerticesForRadius(radius);

        const {vertices} = vobj;
        const {normals} = vobj;

        const geoGroup = geo.updateGeoGroup(vertices.length);

        const start = geoGroup.vertices;
        const {vertexArray} = geoGroup;
        const {colorArray} = geoGroup;
        const {faceArray} = geoGroup;
        const {lineArray} = geoGroup;
        const {normalArray} = geoGroup;

        for (let i = 0, il = vertices.length; i < il; ++i) {
            const offset = 3 * (start + i);
            const v = vertices[i];

            vertexArray[offset] = (v.x + pos.x);
            vertexArray[offset + 1] = (v.y + pos.y);
            vertexArray[offset + 2] = (v.z + pos.z);

            colorArray[offset] = color.r;
            colorArray[offset + 1] = color.g;
            colorArray[offset + 2] = color.b;

        }

        geoGroup.vertices += vertices.length;

        const {verticesRows} = vobj;
        const h = verticesRows.length - 1;

        for (y = 0; y < h; y++) {
            const w = verticesRows[y].length - 1;
            for (x = 0; x < w; x++) {

                const faceoffset = geoGroup.faceidx; const lineoffset = geoGroup.lineidx;

                const v1 = verticesRows[y][x + 1] + start; const v1offset = v1 * 3;
                const v2 = verticesRows[y][x] + start; const v2offset = v2 * 3;
                const v3 = verticesRows[y + 1][x] + start; const v3offset = v3 * 3;
                const v4 = verticesRows[y + 1][x + 1] + start; const v4offset = v4 * 3;

                const n1 = normals[v1 - start];
                const n2 = normals[v2 - start];
                const n3 = normals[v3 - start];
                const n4 = normals[v4 - start];

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

    return draw;

})();
