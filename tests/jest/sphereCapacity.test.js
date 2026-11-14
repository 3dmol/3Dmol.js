global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const FULL = 65535;

function atom(x, bonds = [], order = 1) {
    return { elem: "C", x, y: 0, z: 0, bonds, bondOrder: bonds.map(() => order) };
}
function geometries(atoms, style) {
    const m = new $3Dmol.GLModel(0);
    m.addAtoms(atoms);
    if (style) m.setStyle({}, style);
    let root;
    m.globj({ add(object) { root = object; }, remove() {} }, { supportsImposters: true });
    return root.children.map(child => child.geometry);
}
function expectBackingCapacity(group, vertices) {
    for (const name of ["vertexArray", "colorArray", "normalArray"]) {
        expect(group[name].buffer.byteLength).toBe(vertices * 3 * Float32Array.BYTES_PER_ELEMENT);
    }
    expect(group.faceArray.buffer.byteLength).toBe(Math.floor(vertices / 4) * 6 * Uint16Array.BYTES_PER_ELEMENT);
}
const Geometry = geometries([atom(0)], { sphere: {} })[0].constructor;

describe("sphere imposter allocation hints", () => {
    test("small sphere models allocate exactly four vertices and six indices per atom", () => {
        const atoms = [atom(0), atom(3), atom(6)];
        const geometry = geometries(atoms, { sphere: { radius: 0.5 } })[0];
        expect(geometry.groups).toBe(1);
        const group = geometry.geometryGroups[0];
        expectBackingCapacity(group, atoms.length * 4);
        expect(group.vertices).toBe(12);
        expect(group.faceArray.length).toBe(18);
        expect(group.lineArray.length).toBe(0);
        expect([0, 12, 24].map(offset => group.vertexArray[offset])).toEqual([0, 3, 6]);
        expect(Array.from(group.faceArray.slice(-6))).toEqual([8, 9, 10, 10, 11, 8]);
    });

    test("opaque and translucent sphere groups each use the model bound and retain per-vertex alpha", () => {
        const atoms = [1, 0.25, 1, 0.25].map((opacity, i) => ({
            ...atom(i * 3), style: { sphere: { opacity } }
        }));
        const meshes = geometries(atoms);
        expect(meshes).toHaveLength(2);
        const opaque = meshes.find(g => !g.alpha).geometryGroups[0];
        const translucent = meshes.find(g => g.alpha).geometryGroups[0];
        for (const group of [opaque, translucent]) {
            expectBackingCapacity(group, atoms.length * 4);
            expect(group.vertices).toBe(8);
            expect(group.faceArray.length).toBe(12);
        }
        expect(opaque.alphaArray).toBeNull();
        expect(translucent.alphaArray.buffer.byteLength).toBe(atoms.length * 4 * Float32Array.BYTES_PER_ELEMENT);
        expect(Array.from(translucent.alphaArray)).toEqual(Array(8).fill(0.25));
    });

    test("joint caps and nonbonded stick spheres fit one quad per atom", () => {
        const atoms = [atom(0, [1], 2), atom(1.5, [0], 2), atom(5)];
        const meshes = geometries(atoms, { stick: { showNonBonded: true } });
        const sticks = meshes.find(g => g.sphereGeometry);
        const caps = sticks.sphereGeometry;
        expect(meshes).toContain(caps);
        expect(caps.geometryGroups[0].vertices).toBe(atoms.length * 4);
        expectBackingCapacity(caps.geometryGroups[0], atoms.length * 4);
        expectBackingCapacity(sticks.geometryGroups[0], FULL);
    });

    test("aromatic ring and stick buffers keep their full capacity", () => {
        const atoms = Array.from({ length: 6 }, (_, i) => ({
            ...atom(2 * Math.cos(i * Math.PI / 3), [(i + 5) % 6, (i + 1) % 6], 4),
            y: 2 * Math.sin(i * Math.PI / 3)
        }));
        const meshes = geometries(atoms, { stick: { aromaticStyle: "circle" } });
        const sticks = meshes.find(g => g.sphereGeometry);
        const ring = meshes.find(g => g.radii && !g.sphereGeometry);
        expect(ring.geometryGroups[0].vertices).toBe(4);
        expectBackingCapacity(ring.geometryGroups[0], FULL);
        expectBackingCapacity(sticks.geometryGroups[0], FULL);
        expectBackingCapacity(sticks.sphereGeometry.geometryGroups[0], atoms.length * 4);
    });

    test("a reservation larger than the hint receives sufficient space and later groups retain the hint", () => {
        const geometry = new Geometry(true, true, false, true, 4);
        geometry.imposter = true;
        const first = geometry.updateGeoGroup(12);
        expectBackingCapacity(first, 12);
        expect(first.radiusArray.length).toBe(12);
        expect(first.alphaArray.length).toBe(12);
        first.vertices = 12;
        const next = geometry.updateGeoGroup(4);
        expect(next).not.toBe(first);
        expectBackingCapacity(next, 4);
        expect(next.radiusArray.length).toBe(4);
        expect(next.alphaArray.length).toBe(4);
    });

    test("hints clamp to valid capacities, round up and default on nonfinite input", () => {
        for (const [hint, capacity] of [[0, 4], [-8, 4], [6.2, 7], [100000, FULL], [NaN, FULL], [Infinity, FULL], [-Infinity, FULL]]) {
            const geometry = new Geometry(true, false, false, false, hint);
            geometry.imposter = true;
            expectBackingCapacity(geometry.addGeoGroup(), capacity);
        }
    });

    test("default imposters and ordinary meshes retain their previous capacity", () => {
        const imposter = new Geometry(true), ordinary = new Geometry(true, false, false, false, 4);
        imposter.imposter = true;
        expectBackingCapacity(imposter.addGeoGroup(), FULL);
        const mesh = ordinary.addGeoGroup();
        expect(mesh.vertexArray.length).toBe(FULL * 3);
        expect(mesh.normalArray.length).toBe(FULL * 3);
        expect(mesh.colorArray.length).toBe(FULL * 3);
        expect(mesh.faceArray.length).toBe(FULL * 6);
        expect(mesh.lineArray.length).toBe(FULL * 6);
    });
});
