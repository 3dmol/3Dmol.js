global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const faces = [0, 1, 2, 0, 2, 3];
const triangleEdges = [0, 1, 0, 2, 1, 2, 0, 2, 0, 3, 2, 3];
const vertices = [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 1, y: 1, z: 0 }, { x: 0, y: 1, z: 0 }];

// Obtain the internal Geometry class through a real public shape.
const sample = new $3Dmol.GLShape({});
sample.addCustom({ vertexArr: vertices, faceArr: faces });
let rendered;
sample.globj({ add(object) { rendered = object; }, remove() {} });
const Geometry = rendered.children[0].geometry.constructor;

function quad(imposter = true) {
    const geometry = new Geometry(true);
    geometry.imposter = imposter;
    const group = geometry.addGeoGroup();
    group.vertices = 4;
    group.faceidx = faces.length;
    group.faceArray.set(faces);
    return { geometry, group };
}

describe("imposter line-index storage", () => {
    test("skips the full-capacity line allocation while retaining the face allocation", () => {
        const Native = global.Uint16Array, lengths = [];
        global.Uint16Array = new Proxy(Native, {
            construct(target, args) { lengths.push(args[0]); return Reflect.construct(target, args); }
        });
        let imposter, mesh;
        try {
            imposter = quad(true);
            mesh = quad(false);
        } finally { global.Uint16Array = Native; }
        expect(imposter.group.lineArray).toBeNull();
        expect(mesh.group.lineArray.length).toBe(mesh.group.faceArray.length);
        // One face array for each group, plus the ordinary mesh's line array.
        expect(lengths).toEqual([
            imposter.group.faceArray.length, mesh.group.faceArray.length, mesh.group.lineArray.length
        ]);
    });

    test("initialization retains the established empty line-array result", () => {
        const { geometry, group } = quad();
        expect(group.lineArray).toBeNull();
        geometry.initTypedArrays();
        expect(group.lineArray).toBeInstanceOf(Uint16Array);
        expect(group.lineArray).toHaveLength(0);
        expect(group.lineArray.buffer.byteLength).toBe(0);
        expect(Array.from(group.faceArray)).toEqual(faces);
    });

    test.each([false, true])("wireframe indices can be requested after initialization: %s", initialized => {
        const { geometry, group } = quad();
        if (initialized) geometry.initTypedArrays();
        geometry.setUpWireframe();
        expect(Array.from(group.lineArray)).toEqual(triangleEdges);
        expect(group.lineidx).toBe(triangleEdges.length);
        expect(group.lineArray.buffer.byteLength).toBe(triangleEdges.length * Uint16Array.BYTES_PER_ELEMENT);
        geometry.initTypedArrays();
        expect(Array.from(group.lineArray)).toEqual(triangleEdges);
    });

    test("ordinary meshes retain explicit perimeter edges written before initialization", () => {
        const { geometry, group } = quad(false);
        const perimeter = [0, 1, 1, 2, 2, 3, 3, 0];
        const storage = group.lineArray.buffer;
        group.lineArray.set(perimeter);
        group.lineidx = perimeter.length;
        geometry.initTypedArrays();
        expect(Array.from(group.lineArray)).toEqual(perimeter);
        expect(group.lineArray.buffer).toBe(storage);
    });

    test("public custom shapes still generate both triangles' wireframe edges", () => {
        const shape = new $3Dmol.GLShape({ wireframe: true });
        shape.addCustom({ vertexArr: vertices, faceArr: faces });
        let object;
        shape.globj({ add(value) { object = value; }, remove() {} });
        const geometry = object.children[0].geometry, group = geometry.geometryGroups[0];
        expect(geometry.imposter).toBe(false);
        expect(Array.from(group.faceArray)).toEqual(faces);
        expect(Array.from(group.lineArray)).toEqual(triangleEdges);
        expect(group.lineidx).toBe(triangleEdges.length);
    });

    test("group rollover preserves radius, alpha, normal and offset data", () => {
        const geometry = new Geometry(true, true, true, true);
        geometry.imposter = true;
        const first = geometry.updateGeoGroup(4);
        first.vertices = first.vertexArray.length / 3;
        const second = geometry.updateGeoGroup(4);
        expect(second).not.toBe(first);
        expect(geometry.groups).toBe(2);
        for (const [i, group] of geometry.geometryGroups.entries()) {
            expect(group.lineArray).toBeNull();
            expect(group.useOffset).toBe(true);
            group.vertices = 4;
            group.faceidx = faces.length;
            group.faceArray.set(faces);
            group.normalArray[0] = i + 1;
            group.radiusArray[0] = i + 2;
            group.alphaArray[0] = 0.5;
        }
        geometry.initTypedArrays();
        for (const [i, group] of geometry.geometryGroups.entries()) {
            expect(group.normalArray).toHaveLength(12);
            expect(group.normalArray[0]).toBe(i + 1);
            expect(group.radiusArray).toHaveLength(4);
            expect(group.radiusArray[0]).toBe(i + 2);
            expect(group.alphaArray).toHaveLength(4);
            expect(group.alphaArray[0]).toBe(0.5);
            expect(group.lineArray).toHaveLength(0);
            expect(Array.from(group.faceArray)).toEqual(faces);
        }
    });
});
