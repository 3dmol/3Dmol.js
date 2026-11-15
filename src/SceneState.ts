// Serializable scene values contain data only: never revive executable callbacks.
import { Color } from './colors';
import { GradientType, getGradient } from './Gradient';
import { Geometry, GeometryGroup } from './WebGL/core/Geometry';
import { Matrix4 } from './WebGL/math';
import { VolumeData } from './VolumeData';

const omitted = new Set(['__proto__', 'constructor', 'prototype', 'callback', 'hover_callback',
    'unhover_callback', 'intersectionShape', 'capDrawn']);
const arrays = { Float32Array, Float64Array, Int8Array, Uint8Array, Uint8ClampedArray,
    Int16Array, Uint16Array, Int32Array, Uint32Array };

/** Make a detached, JSON-safe copy of supported scene data. */
export function encodeSceneValue(value: any, ancestors = new Set<object>()): any {
    if (value === null || value === undefined || typeof value === 'string' || typeof value === 'boolean') return value;
    if (typeof value === 'number') return Number.isFinite(value) ? value : { __3dmol: 'Number', value: String(value) };
    if (typeof value !== 'object') return undefined;
    // Selections may refer to model objects; the portable form is the model ID.
    if (typeof value.getID === 'function') return value.getID();
    if (value instanceof Color) return value.getHex();
    if (ArrayBuffer.isView(value)) return { __3dmol: value.constructor.name,
        values: Array.from(value as any, (number: number) => Number.isFinite(number) ? number : String(number)) };
    if (value instanceof ArrayBuffer) return { __3dmol: 'ArrayBuffer', values: Array.from(new Uint8Array(value)) };
    if (value instanceof Matrix4) return { __3dmol: 'Matrix4', elements: Array.from(value.elements) };
    if (value instanceof VolumeData) return { __3dmol: 'VolumeData',
        data: encodeSceneValue(value.data), size: encodeSceneValue(value.size),
        origin: encodeSceneValue(value.origin), unit: encodeSceneValue(value.unit), matrix: encodeSceneValue(value.matrix) };
    if (value instanceof GradientType) return { __3dmol: 'Gradient', gradient: value.gradient,
        min: value.range()?.[0], max: value.range()?.[1], mid: (value as any).mid,
        colors: encodeSceneValue((value as any).colors) };
    if (ancestors.has(value)) return undefined;
    ancestors.add(value);
    let result: any;
    if (Array.isArray(value)) result = value.map(item => encodeSceneValue(item, ancestors) ?? null);
    else {
        result = {};
        for (const key of Object.keys(value)) {
            if (omitted.has(key)) continue;
            const item = encodeSceneValue(value[key], ancestors);
            if (item !== undefined) result[key] = item;
        }
    }
    ancestors.delete(value);
    return result;
}

/** Restore data classes from a snapshot using an explicit type allowlist. */
export function decodeSceneValue(value: any): any {
    if (!value || typeof value !== 'object') return value;
    if (Array.isArray(value)) return value.map(decodeSceneValue);
    if (value.__3dmol === 'Number') return Number(value.value);
    if (Object.prototype.hasOwnProperty.call(arrays, value.__3dmol)) return new arrays[value.__3dmol](value.values.map(Number));
    if (value.__3dmol === 'ArrayBuffer') return new Uint8Array(value.values).buffer;
    if (value.__3dmol === 'Matrix4') {
        const matrix = new Matrix4();
        matrix.elements.set(value.elements);
        return matrix;
    }
    if (value.__3dmol === 'VolumeData') {
        const volume = Object.create(VolumeData.prototype);
        for (const key of ['data', 'size', 'origin', 'unit', 'matrix']) volume[key] = decodeSceneValue(value[key]);
        return volume;
    }
    if (value.__3dmol === 'Gradient') return getGradient({ gradient: value.gradient, min: value.min,
        max: value.max, mid: value.mid, colors: decodeSceneValue(value.colors) });
    const result = {};
    for (const key of Object.keys(value)) if (!omitted.has(key)) result[key] = decodeSceneValue(value[key]);
    return result;
}

const floatBuffers = ['vertexArray', 'colorArray', 'normalArray', 'radiusArray', 'alphaArray'];

/** Snapshot only populated geometry buffers, without WebGL resources or spare capacity. */
export function snapshotGeometry(geometry: Geometry): any {
    if (!geometry) return null;
    return { mesh: geometry.mesh, radii: geometry.radii, offset: geometry.offset, alpha: geometry.alpha,
        groups: geometry.geometryGroups.filter(group => group.vertices > 0).map(group => {
            const result: any = { vertices: group.vertices, faceidx: group.faceidx, lineidx: group.lineidx };
            for (const key of floatBuffers) if (group[key]) {
                const length = group.vertices * (key === 'radiusArray' || key === 'alphaArray' ? 1 : 3);
                result[key] = Array.from(group[key].subarray(0, length));
            }
            if (group.faceArray) result.faceArray = Array.from(group.faceArray.subarray(0, group.faceidx));
            if (group.lineArray) result.lineArray = Array.from(group.lineArray.subarray(0, group.lineidx));
            if (group.atomArray?.length) {
                const indices = new Map(), atoms = [];
                result.atomIndices = group.atomArray.slice(0, group.vertices).map(atom => {
                    if (!indices.has(atom)) { indices.set(atom, atoms.length); atoms.push(encodeSceneValue(atom)); }
                    return indices.get(atom);
                });
                result.atoms = atoms;
            }
            return result;
        }) };
}

export function restoreGeometry(snapshot: any): Geometry {
    const geometry = new Geometry(snapshot?.mesh ?? true, !!snapshot?.radii, !!snapshot?.offset, !!snapshot?.alpha);
    for (const record of snapshot?.groups || []) {
        const group = new GeometryGroup(geometry.geometryGroups.length);
        group.vertices = record.vertices; group.faceidx = record.faceidx; group.lineidx = record.lineidx;
        for (const key of floatBuffers) if (record[key]) group[key] = new Float32Array(record[key]);
        if (record.faceArray) group.faceArray = new Uint16Array(record.faceArray);
        if (record.lineArray) group.lineArray = new Uint16Array(record.lineArray);
        const atoms = (record.atoms || []).map(decodeSceneValue);
        group.atomArray = (record.atomIndices || []).map(index => atoms[index]);
        group.__inittedArrays = true;
        geometry.geometryGroups.push(group);
    }
    geometry.groups = geometry.geometryGroups.length;
    return geometry;
}
