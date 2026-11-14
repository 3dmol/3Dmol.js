global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const QUADS = 16383;
const FACE_CAPACITY = QUADS * 6;
const pattern = [0, 1, 2, 2, 3, 0];

function modelGeometry(atoms, style) {
    const m = new $3Dmol.GLModel(0);
    m.addAtoms(atoms);
    m.setStyle({}, style);
    let root;
    m.globj({ add(object) { root = object; }, remove() {} }, { supportsImposters: true });
    return root.children[0].geometry;
}
function atom(x, elem = "C", bonds = []) {
    return { elem, x, y: 0, z: 0, bonds, bondOrder: bonds.map(() => 1) };
}
function expectQuads(geometry, counts) {
    expect(geometry.geometryGroups.map(g => g.vertices)).toEqual(counts.map(n => n * 4));
    expect(geometry.geometryGroups.map(g => g.faceidx)).toEqual(counts.map(n => n * 6));
    for (const group of geometry.geometryGroups) {
        expect(group.faceArray.length).toBe(group.faceidx);
        expect(group.faceArray.buffer.byteLength).toBe(FACE_CAPACITY * Uint16Array.BYTES_PER_ELEMENT);
        // Check every index without creating thousands of assertion objects.
        expect(group.faceArray.every((index, i) => index === Math.floor(i / 6) * 4 + pattern[i % 6])).toBe(true);
        expect(Array.from(group.faceArray.slice(-6))).toEqual(pattern.map(offset => group.vertices - 4 + offset));
    }
}
const Geometry = modelGeometry([atom(0)], { sphere: {} }).constructor;

describe("imposter quad face storage", () => {
    test.each([QUADS, QUADS + 1])("%i sphere quads fill and roll over without truncating triangles", count => {
        const atoms = Array.from({ length: count }, (_, i) => atom(i * 3));
        const geometry = modelGeometry(atoms, { sphere: {} });
        expectQuads(geometry, count === QUADS ? [QUADS] : [QUADS, 1]);
        const first = geometry.geometryGroups[0];
        expect(first.vertexArray[(first.vertices - 4) * 3]).toBe((QUADS - 1) * 3);
        if (count > QUADS) expect(geometry.geometryGroups[1].vertexArray[0]).toBe(QUADS * 3);
    });

    test("color-split sticks remain continuous when a bond's halves straddle groups", () => {
        const bonds = 8192, atoms = [];
        for (let i = 0; i < bonds; i++) {
            atoms.push(atom(i * 3, "C", [2 * i + 1]));
            atoms.push(atom(i * 3 + 1.5, "O", [2 * i]));
        }
        const geometry = modelGeometry(atoms, { stick: {} });
        expectQuads(geometry, [QUADS, 1]);
        const [first, second] = geometry.geometryGroups;
        const from = (bonds - 1) * 3, midpoint = from + 0.75, to = from + 1.5;
        const last = (first.vertices - 4) * 3;
        for (let corner = 0; corner < 4; corner++) {
            expect(first.vertexArray[last + corner * 3]).toBe(from);
            expect(first.normalArray[last + corner * 3]).toBe(midpoint);
            expect(second.vertexArray[corner * 3]).toBe(midpoint);
            expect(second.normalArray[corner * 3]).toBe(to);
        }
    });

    test("the aromatic ring writer retains the last quad on each side of a group boundary", () => {
        const m = new $3Dmol.GLModel(0), geometry = new Geometry(true, true);
        geometry.imposter = true;
        const normal = { x: 0, y: 0, z: 1 }, color = new $3Dmol.Color(0xff0000);
        for (let i = 0; i <= QUADS; i++) {
            m.drawRingImposter(geometry, { x: i * 3, y: 0, z: 0 }, normal, 2, 0.5, color);
        }
        geometry.initTypedArrays();
        expectQuads(geometry, [QUADS, 1]);
        const [first, second] = geometry.geometryGroups;
        expect(first.vertexArray[(first.vertices - 4) * 3]).toBe((QUADS - 1) * 3);
        expect(second.vertexArray[0]).toBe(QUADS * 3);
        expect(Array.from(first.radiusArray.slice(-4))).toEqual([-2.5, -2.5, 2.5, 2.5]);
        expect(Array.from(second.radiusArray)).toEqual([-2.5, -2.5, 2.5, 2.5]);
    });

    test("only imposter face capacity is reduced; ordinary mesh storage is unchanged", () => {
        const imposter = new Geometry(true), ordinary = new Geometry(true);
        imposter.imposter = true;
        const compact = imposter.addGeoGroup(), mesh = ordinary.addGeoGroup();
        expect(compact.faceArray).toHaveLength(FACE_CAPACITY);
        expect(compact.lineArray).toBeNull();
        expect(mesh.faceArray).toHaveLength(65535 * 6);
        expect(mesh.lineArray).toHaveLength(65535 * 6);
        expect(compact.vertexArray.length).toBe(mesh.vertexArray.length);
        expect(compact.normalArray.length).toBe(mesh.normalArray.length);
        expect(compact.colorArray.length).toBe(mesh.colorArray.length);
    });
});
