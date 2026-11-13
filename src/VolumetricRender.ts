import { Sphere } from "./WebGL/shapes";
import { Vector3, Matrix4, XYZ } from "./WebGL/math";
import { VolumetricMaterial, Mesh, Texture, Object3D, Material } from "./WebGL";
import { disposeObjectResources } from "./WebGL/disposeObjectResources";
import { CC } from "./colors";
import type { ColorSpec } from "./colors";
import { GLShape } from "./GLShape";
import { AtomSelectionSpec } from "specs";
import { GLViewer } from "GLViewer";



/**
 * VolumetricRenderer style specification
*/
export interface VolumetricRendererSpec {
    /** list of objects containing @color, @opacity and @value properties to specify color per voxel data value.
     * Requires at least two points with finite numeric values (numbers or numeric strings) spanning a nonzero range that can be sampled without numeric overflow.
     * Opacity ranges from 0 (transparent) to 1 (opaque); numeric strings are also accepted.
     * The input array and point objects are not modified.
     * Missing or invalid transfer functions log a warning and skip rendering. */
    transferfn?: { color: ColorSpec; opacity: number | string; value: number | string }[];
    /** number of times to sample each voxel approximately (default 5) */
    subsamples?: number;
    /**  coordinates around which to include data; use viewer.selectedAtoms() to convert an AtomSelectionSpec to coordinates */
    coords?: XYZ[];
    /** Selection around which to include data, resolved using the supplied viewer when coords are omitted. */
    selection?: AtomSelectionSpec;
    /** distance around coords or selection to include data [default = 2.0] */
    seldist?: number; 
};

/**
 * A GLVolumetricRender is a "shape" for representing volumetric data as a density distribution.
 *
 * @class
 *
 * @param {VolumeData} data - volumetric data
 * @param {VolumetricRendererSpec} spec - specification of volumetric render
 * @param {GLViewer} [viewer] - Viewer used to resolve spec.selection when coords are not supplied.
 */
export class GLVolumetricRender {

    static interpolateArray(data: string | any[], fitCount: number) {
        function linearInterpolate(before: number, after: number, atPoint: number) {
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

    constructor(data: { matrix: { elements: any; }; size: XYZ; 
                unit: XYZ; origin: XYZ; data: number[]; getIndex: (arg0: number, arg1: number, arg2: number) => number; }, 
        spec: VolumetricRendererSpec, viewer?: GLViewer) {
        spec = spec || {};
        // Array.from also exposes sparse entries so they are validated below.
        var transferfn: any[] = Array.isArray(spec.transferfn) ? Array.from(spec.transferfn) : [];
        if (transferfn.length < 2 || transferfn.some(point => !point ||
            (typeof point.value !== "number" && typeof point.value !== "string") ||
            !Number.isFinite(parseFloat(String(point.value))))) {
            console.warn("Volumetric rendering requires at least two transfer function points with finite values.");
            this.hidden = true;
            return;
        }
        this.subsamples = spec.subsamples || 5.0;

        let TRANSFER_BUFFER_SIZE = 256;

        // Normalize private point records, leaving the caller's values and order intact.
        for (let i = 0; i < transferfn.length; i++) {
            const point = transferfn[i];
            transferfn[i] = {
                value: parseFloat(point.value), color: point.color, opacity: point.opacity
            };
        }
        transferfn.sort(function (a: { value: number; }, b: { value: number; }) { return a.value - b.value; });
        this.min = transferfn[0].value;
        this.max = transferfn[transferfn.length - 1].value;
        if (!Number.isFinite((this.max - this.min) * TRANSFER_BUFFER_SIZE) || this.min === this.max) {
            console.warn("Volumetric rendering requires a nonzero transfer function range that can be sampled without numeric overflow.");
            this.hidden = true;
            return;
        }

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

        // Sampling distances are lengths, independent of grid rotation or reflection.
        let spacing: number[];
        if (data.matrix) {
            const e = data.matrix.elements;
            spacing = [
                Math.hypot(e[0], e[1], e[2]),
                Math.hypot(e[4], e[5], e[6]),
                Math.hypot(e[8], e[9], e[10])
            ];
        } else {
            spacing = [Math.abs(data.unit.x), Math.abs(data.unit.y), Math.abs(data.unit.z)];
        }
        if (spacing.some(step => !Number.isFinite(step) || step <= 0)) {
            console.warn("Volumetric rendering requires finite, nonzero voxel spacing.");
            this.hidden = true;
            return;
        }
        this.minunit = Math.min(...spacing);

        //need to create transformation matrix that maps model points into
        //texture space
        // need extent (bounding box dimensions), maxdepth (box diagonal), 
        // texmatrix (conversion from model to texture coords), minunit,
        // possibly non-orthnormal basis if matrix
        if (data.matrix) {
            //figure out bounding box of transformed grid
            let start = new Vector3(0, 0, 0);
            let end = new Vector3(data.size.x, data.size.y, data.size.z);

            start.applyMatrix4(data.matrix);
            end.applyMatrix4(data.matrix);

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

        if (spec.coords === undefined && spec.selection !== undefined) {
            if(viewer) {
                spec.coords = viewer.selectedAtoms(spec.selection) as XYZ[];
            } else {
                console.log('Need to provide viewer to volumetric renderer if selection specified.');
            }
        }
        this.data = data;
        // Select voxel centers without modifying the caller's volume.
        if (spec.coords !== undefined) {
            const d = spec.seldist === undefined ? 2.0 : spec.seldist;
            if (!Number.isFinite(d) || d < 0) {
                console.warn("Volumetric selection requires a finite, nonnegative distance.");
                this.hidden = true;
                return;
            }
            const d2 = d * d;
            const size = data.size;
            // Keep axis-aligned origins and spacing at their original precision.
            const world = data.matrix ? data.matrix.elements : [
                data.unit.x, 0, 0, 0,
                0, data.unit.y, 0, 0,
                0, 0, data.unit.z, 0,
                data.origin.x, data.origin.y, data.origin.z, 1
            ];
            const determinant =
                world[0] * (world[5] * world[10] - world[9] * world[6]) -
                world[4] * (world[1] * world[10] - world[9] * world[2]) +
                world[8] * (world[1] * world[6] - world[5] * world[2]);
            if (data.matrix && (!Number.isFinite(determinant) || determinant === 0)) {
                console.warn("Volumetric selection requires an invertible grid transform.");
                this.hidden = true;
                return;
            }
            const inverse = data.matrix ? new Matrix4().getInverse(data.matrix as Matrix4) : null;
            const inv = inverse && inverse.elements;
            // A world-space sphere becomes an ellipsoid in grid coordinates.
            const rx = d * (inv ? Math.hypot(inv[0], inv[4], inv[8]) : 1 / Math.abs(data.unit.x));
            const ry = d * (inv ? Math.hypot(inv[1], inv[5], inv[9]) : 1 / Math.abs(data.unit.y));
            const rz = d * (inv ? Math.hypot(inv[2], inv[6], inv[10]) : 1 / Math.abs(data.unit.z));
            const masked = new Float32Array(data.data.length).fill(Infinity);
            for (const c of spec.coords) {
                const center = inverse ? new Vector3(c.x, c.y, c.z).applyMatrix4(inverse) : new Vector3(
                    (c.x - data.origin.x) / data.unit.x,
                    (c.y - data.origin.y) / data.unit.y,
                    (c.z - data.origin.z) / data.unit.z
                );
                // Conservative bounds include edge candidates; the distance test is exact.
                const minx = Math.max(0, Math.floor(center.x - rx));
                const maxx = Math.min(size.x - 1, Math.ceil(center.x + rx));
                const miny = Math.max(0, Math.floor(center.y - ry));
                const maxy = Math.min(size.y - 1, Math.ceil(center.y + ry));
                const minz = Math.max(0, Math.floor(center.z - rz));
                const maxz = Math.min(size.z - 1, Math.ceil(center.z + rz));
                for (let x = minx; x <= maxx; x++) {
                    for (let y = miny; y <= maxy; y++) {
                        const px = world[0] * x + world[4] * y + world[12] - c.x;
                        const py = world[1] * x + world[5] * y + world[13] - c.y;
                        const pz = world[2] * x + world[6] * y + world[14] - c.z;
                        let idx = (x * size.y + y) * size.z + minz;
                        for (let z = minz; z <= maxz; z++, idx++) {
                            if (masked[idx] !== Infinity) continue;
                            const dx = px + world[8] * z;
                            const dy = py + world[9] * z;
                            const dz = pz + world[10] * z;
                            if (dx * dx + dy * dy + dz * dz < d2) masked[idx] = data.data[idx];
                        }
                    }
                }
            }
            this.data = Object.assign(Object.create(Object.getPrototypeOf(data)), data, { data: masked });
        }
    }

    /**
     * Initialize webgl objects for rendering
     * @param {Object3D} group
     *
     */
    globj(group: { remove: (arg0: any) => void; add: (arg0: any) => void; }) {

        if (this.renderedShapeObj) {
            group.remove(this.renderedShapeObj);
            disposeObjectResources(this.renderedShapeObj, {
                geometries: false, volumeTextures: true
            });
            this.renderedShapeObj = null;
        }
        this.shapeObj = null;

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

        var mesh = new Mesh(this.geo, (material as Material));
        this.shapeObj.add(mesh);

        this.renderedShapeObj = this.shapeObj.clone();
        group.add(this.renderedShapeObj);
    };

    removegl(group: { remove: (arg0: any) => void; }) {
        if (this.renderedShapeObj) {
            group.remove(this.renderedShapeObj);
            disposeObjectResources(this.renderedShapeObj, { geometries: false, volumeTextures: true });
            this.renderedShapeObj = null;
        }
        // Invalid volume input may return from the constructor before creating geometry.
        this.geo?.dispose();
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

