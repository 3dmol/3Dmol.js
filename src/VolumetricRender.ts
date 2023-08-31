import { Sphere } from "./WebGL/shapes";
import { Vector3, Matrix4 } from "./WebGL/math";
import { VolumetricMaterial, Mesh, Texture, Object3D } from "./WebGL";
import { CC } from "./colors";
import { GLShape } from "./GLShape";



/**
 * VolumetricRenderer style specification
*/
export interface VolumetricRendererSpec {
    /** list of objects containing @color, @opacity and @value properties to specify color per voxel data value */
    transferfn?: { color: unknown; opacity: unknown; value: unknown }[];
    /** number of times to sample each voxel approximately (default 5) */
    subsamples?: number;
    /**  coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates */
    coords?: Array<Vector3>;
    /** distance around coords to include data [default = 2.0] */
    seldist?: number; 
};

/**
 * A GLVolumetricRender is a "shape" for representing volumetric data as a density distribution.
 *
 * @class
 *
 * @param {VolumeData} data - volumetric data
 * @param {VolumetricRenderSpec} spec - specification of volumetric render
 * @returns {$3Dmol.GLShape}
 */
export class GLVolumetricRender {

    static interpolateArray(data, fitCount) {
        function linearInterpolate(before, after, atPoint) {
            return before + (after - before) * atPoint;
        }
        var newData = [];
        var springFactor = (data.length - 1) / (fitCount - 1);
        newData[0] = data[0]; // for new allocation
        for (var i = 1; i < fitCount - 1; i++) {
            var tmp = i * springFactor;
            var before = Math.floor(tmp);
            var after = Math.ceil(tmp);
            var atPoint = tmp - before;
            newData[i] = linearInterpolate(data[before], data[after], atPoint);
        }
        newData[fitCount - 1] = data[data.length - 1]; // for new allocation
        return newData;
    }

    hidden = false;
    boundingSphere = new Sphere();
    shapePosition: any;
    renderedShapeObj: any = null;
    shapeObj: any = null;
    geo: any;
    subsamples = 5.0;
    data: any = null;
    transferfunctionbuffer: any = [];
    min: number = 0;
    max: number = 0;
    extent: any;
    maxdepth: number;
    texmatrix: any;
    minunit: any;

    constructor(data, spec) {
        spec = spec || {};
        var transferfn = Object.assign([], spec.transferfn);
        this.subsamples = spec.subsamples || 5.0;

        let TRANSFER_BUFFER_SIZE = 256;

        // arrange points based on position property
        transferfn.forEach(function (a) { a.value = parseFloat(a.value); });
        transferfn.sort(function (a, b) { return a.value - b.value; });
        this.min = transferfn[0].value;
        if (transferfn.length == 0) transferfn.push(transferfn[0]); //need at least two
        this.max = transferfn[transferfn.length - 1].value;

        // create and fill an array of interpolated values per 2 colors
        var pos1, pos2, color1, color2, R, G, B, A, alpha1, alpha2;
        for (let i = 0; i < transferfn.length - 1; i++) {
            color1 = CC.color(transferfn[i].color);
            color2 = CC.color(transferfn[i + 1].color);
            alpha1 = transferfn[i].opacity;
            alpha2 = transferfn[i + 1].opacity;
            pos1 = Math.floor((transferfn[i].value - this.min) * TRANSFER_BUFFER_SIZE / (this.max - this.min));
            pos2 = Math.floor((transferfn[i + 1].value - this.min) * TRANSFER_BUFFER_SIZE / (this.max - this.min));
            if (pos1 == pos2)
                continue;
            R = GLVolumetricRender.interpolateArray([color1.r * 255, color2.r * 255], pos2 - pos1);
            G = GLVolumetricRender.interpolateArray([color1.g * 255, color2.g * 255], pos2 - pos1);
            B = GLVolumetricRender.interpolateArray([color1.b * 255, color2.b * 255], pos2 - pos1);
            A = GLVolumetricRender.interpolateArray([alpha1 * 255, alpha2 * 255], pos2 - pos1);

            for (let j = 0; j < R.length; j++) {
                this.transferfunctionbuffer.push(R[j]);
                this.transferfunctionbuffer.push(G[j]);
                this.transferfunctionbuffer.push(B[j]);
                this.transferfunctionbuffer.push(A[j]); // opacity will be added later
            }
        }

        this.transferfunctionbuffer = new Uint8ClampedArray(this.transferfunctionbuffer);

        //need to create transformation matrix that maps model points into
        //texture space
        // need extent (bounding box dimensions), maxdepth (box diagonal), 
        // texmatrix (conversion from model to texture coords), minunit,
        // possibly non-orthnormal basis if matrix
        if (data.matrix) {
            //figure out bounding box of transformed grid
            let start = new Vector3(0, 0, 0);
            let end = new Vector3(data.size.x, data.size.y, data.size.z);
            let unit = new Vector3(1, 1, 1);

            start.applyMatrix4(data.matrix);
            end.applyMatrix4(data.matrix);
            unit.applyMatrix4(data.matrix).sub(start);

            this.extent = [[start.x, start.y, start.z], [end.x, end.y, end.z]];

            //check all corners, these may not be the farthest apart
            for (let i = 1; i < 7; i++) {
                end.x = (i & 1) ? data.size.x : 0;
                end.y = (i & 2) ? data.size.y : 0;
                end.z = (i & 4) ? data.size.z : 0;
                end.applyMatrix4(data.matrix);
                this.extent[0][0] = Math.min(this.extent[0][0], end.x);
                this.extent[0][1] = Math.min(this.extent[0][1], end.y);
                this.extent[0][2] = Math.min(this.extent[0][2], end.z);
                this.extent[1][0] = Math.max(this.extent[1][0], end.x);
                this.extent[1][1] = Math.max(this.extent[1][1], end.y);
                this.extent[1][2] = Math.max(this.extent[1][2], end.z);
            }

            let xoff = end.x - start.x;
            let yoff = end.y - start.y;
            let zoff = end.z - start.z;
            this.maxdepth = Math.sqrt(xoff * xoff + yoff * yoff + zoff * zoff);

            this.minunit = Math.min(Math.min(unit.x, unit.y), unit.z);

            //invert onto grid, then scale by grid dimensions to get
            //normalized texture coordinates
            this.texmatrix = new Matrix4().identity().scale({ x: data.size.x, y: data.size.y, z: data.size.z });
            this.texmatrix = this.texmatrix.multiplyMatrices(data.matrix, this.texmatrix);

            this.texmatrix = this.texmatrix.getInverse(this.texmatrix);

        } else {
            this.texmatrix = new Matrix4().identity();
            let xoff = data.unit.x * data.size.x;
            let yoff = data.unit.y * data.size.y;
            let zoff = data.unit.z * data.size.z;
            //scale doesn't apply to the translation vector, so preapply it
            this.texmatrix.makeTranslation(-data.origin.x / xoff, -data.origin.y / yoff, -data.origin.z / zoff);
            this.texmatrix.scale({ x: 1.0 / xoff, y: 1.0 / yoff, z: 1.0 / zoff });
            this.minunit = Math.min(Math.min(data.unit.x, data.unit.y), data.unit.z);

            //need the bounding box so we can intersect with rays
            this.extent = [[data.origin.x, data.origin.y, data.origin.z],
            [data.origin.x + xoff, data.origin.y + yoff, data.origin.z + zoff]];

            this.maxdepth = Math.sqrt(xoff * xoff + yoff * yoff + zoff * zoff);
        }

        //use GLShape to construct box
        var shape = new GLShape({});
        shape.addBox({
            corner: { x: this.extent[0][0], y: this.extent[0][1], z: this.extent[0][2] },
            dimensions: {
                w: this.extent[1][0] - this.extent[0][0],
                h: this.extent[1][1] - this.extent[0][1],
                d: this.extent[1][2] - this.extent[0][2]
            }
        });

        this.geo = shape.finalize();
        this.boundingSphere.center = new Vector3(
            (this.extent[0][0] + this.extent[1][0]) / 2.0,
            (this.extent[0][1] + this.extent[1][1]) / 2.0,
            (this.extent[0][2] + this.extent[1][2]) / 2.0
        );
        this.boundingSphere.radius = this.maxdepth / 2;

        // volume selectivity based on given coords and distance
        if (spec.coords !== undefined && spec.seldist !== undefined) {
            let mask = new Uint8Array(data.data.length);
            //for each coordinate
            let d = spec.seldist;
            let d2 = d * d;
            for (let i = 0, n = spec.coords.length; i < n; i++) {
                let c = spec.coords[i];
                let minx = c.x - d, miny = c.y - d, minz = c.z - d;
                let maxx = c.x + d, maxy = c.y + d, maxz = c.z + d;
                if (data.getIndex(minx, miny, minz) >= 0 || data.getIndex(maxx, maxy, maxz) >= 0) {
                    //bounding box overlaps grid
                    //iterate over the grid points in the seldist bounding box
                    //minunit may be inefficient if axes have very different units. oh well.
                    for (let x = minx; x < maxx; x += this.minunit) {
                        for (let y = miny; y < maxy; y += this.minunit) {
                            for (let z = minz; z < maxz; z += this.minunit) {
                                let idx = data.getIndex(x, y, z);
                                if (idx >= 0 && !mask[idx]) {
                                    //if not already masked, check distance
                                    let distsq = (x - c.x) * (x - c.x) + (y - c.y) * (y - c.y) + (z - c.z) * (z - c.z);
                                    if (distsq < d2) {
                                        mask[idx] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            //any place mask is zero, make infinite in data
            for (let i = 0, n = data.data.length; i < n; i++) {
                if (mask[i] == 0) data.data[i] = Infinity;
            }
        }
        this.data = data;
    }

    /**
     * Initialize webgl objects for rendering
     * @param {Object3D} group
     *
     */
    globj(group) {

        if (this.renderedShapeObj) {
            group.remove(this.renderedShapeObj);
            this.renderedShapeObj = null;
        }

        if (this.hidden)
            return;

        this.shapeObj = new Object3D();
        var material = null;

        var texture = new Texture(this.data, true);
        var transfertexture = new Texture(this.transferfunctionbuffer, false);
        texture.needsUpdate = true;
        transfertexture.needsUpdate = true;
        transfertexture.flipY = false;

        material = new VolumetricMaterial({
            transferfn: transfertexture,
            transfermin: this.min,
            transfermax: this.max,
            map: texture,
            extent: this.extent,
            maxdepth: this.maxdepth,
            texmatrix: this.texmatrix,
            unit: this.minunit,
            subsamples: this.subsamples,
        });

        var mesh = new Mesh(this.geo, material);
        this.shapeObj.add(mesh);

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
}

