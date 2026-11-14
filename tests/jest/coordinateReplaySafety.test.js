global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const imposters = { supportsImposters: true, supportsAIA: false };
const group = { add() {}, remove() {} };

function bondedAtoms(length = 1.5, order = 1) {
    return [
        { elem: "C", x: 0, y: 0, z: 0, bonds: [1], bondOrder: [order] },
        { elem: "O", x: length, y: 0, z: 0, bonds: [0], bondOrder: [order] }
    ];
}
function model(atoms, style) {
    const m = new $3Dmol.GLModel(0);
    m.addAtoms(atoms);
    m.addFrame(m.getInternalState().atoms);
    m.setStyle({}, style);
    return m;
}
function draw(m, options = imposters) {
    m.globj(group, options);
    return m.renderedMolObj;
}
function geometries(object) {
    const result = [];
    function visit(node) {
        if (node.geometry) result.push(node.geometry);
        node.children.forEach(visit);
    }
    visit(object);
    return result;
}
function buffers(object) {
    return geometries(object).map(geometry => geometry.geometryGroups.map(g => ({
        vertices: Array.from(g.vertexArray || []), normals: Array.from(g.normalArray || []),
        radii: Array.from(g.radiusArray || []), colors: Array.from(g.colorArray || []),
        faces: Array.from(g.faceArray || [])
    })));
}
function move(atoms) {
    atoms.forEach((a, i) => { a.x += 2 + i; a.y += 3; a.z += 1; });
}
function expectRefusal(m, previous) {
    const before = buffers(previous);
    const geos = geometries(previous);
    for (const g of geos) g.verticesNeedUpdate = g.normalsNeedUpdate = false;
    expect(m.syncAtomPositions()).toBe(false);
    expect(buffers(previous)).toEqual(before);
    expect(geos.every(g => !g.verticesNeedUpdate && !g.normalsNeedUpdate)).toBe(true);
}

const unsupported = [
    ["cross mixed with spheres", () => bondedAtoms(), { sphere: {}, cross: {} }, imposters],
    ["cartoon mixed with spheres", () => bondedAtoms().map((a, i) => ({
        ...a, elem: "C", atom: "CA", resi: i + 1, resn: "ALA", chain: "A", ss: "c"
    })), { sphere: {}, cartoon: { style: "trace" } }, imposters],
    ["cartoon with no initial backbone geometry", () => bondedAtoms(), { sphere: {}, cartoon: {} }, imposters],
    ["tessellated spheres mixed with lines", () => bondedAtoms(), { sphere: {}, line: {} }, {}],
    ["instanced spheres mixed with lines", () => bondedAtoms(), { sphere: {}, line: {} }, { supportsAIA: true }],
    ["tessellated sticks mixed with lines", () => bondedAtoms(), { stick: {}, line: {} }, {}],
    ["dashed stick style", () => bondedAtoms(4), { sphere: {}, stick: { dashedBonds: true } }, imposters],
    ["fractional single stick", () => bondedAtoms(4, 0.5), { sphere: {}, stick: {} }, imposters],
    ["tiny multiple-bond offsets", () => bondedAtoms(1.5, 2), { sphere: {}, stick: { radius: 0.000001 } }, imposters],
    ["coincident stick endpoints", () => bondedAtoms(0), { sphere: {}, stick: {} }, imposters]
];

describe("coordinate replay safety", () => {
    test.each(unsupported)("%s refuses before writes and fast setFrame falls back to a complete rebuild", async (_name, atoms, style, options) => {
        const m = model(atoms(), style), previous = draw(m, options);
        move(m.getInternalState().atoms);
        expectRefusal(m, previous);
        await m.setFrame(0, { fast: true });
        expect(draw(m, options)).not.toBe(previous);
        const expected = model(atoms(), style);
        move(expected.getInternalState().atoms);
        expect(buffers(draw(m, options))).toEqual(buffers(draw(expected, options)));
    });

    test("a degenerate aromatic circle rebuilds when the ring gains an area", async () => {
        const atoms = Array.from({ length: 6 }, (_, i) => ({
            elem: "C", x: i, y: 0, z: 0,
            bonds: [(i + 5) % 6, (i + 1) % 6], bondOrder: [4, 4]
        }));
        const style = { stick: { aromaticStyle: "circle" } };
        const m = model(atoms, style), previous = draw(m);
        const initialMeshes = previous.children.length;
        const spread = list => list.forEach((a, i) => {
            a.x = 2 + 2 * Math.cos(i * Math.PI / 3);
            a.y = 2 * Math.sin(i * Math.PI / 3);
        });
        spread(m.getInternalState().atoms);
        expectRefusal(m, previous);
        await m.setFrame(0, { fast: true });
        const rebuilt = draw(m);
        expect(rebuilt).not.toBe(previous);
        expect(rebuilt.children.length).toBeGreaterThan(initialMeshes);
        const expected = model(atoms, style);
        spread(expected.getInternalState().atoms);
        expect(buffers(rebuilt)).toEqual(buffers(draw(expected)));
    });

    test.each(["clickable", "hoverable", "contextMenuEnabled"])("current-frame %s atoms refuse replay before writes", async flag => {
        const style = { sphere: {}, stick: {}, line: {}, clicksphere: { radius: 2 } };
        const m = model(bondedAtoms(), style), previous = draw(m);
        const atoms = m.getInternalState().atoms;
        // The flag appears only after the geometry was built, on the last atom.
        atoms[1][flag] = true;
        move(atoms);
        expectRefusal(m, previous);
        await m.setFrame(0, { fast: true });
        expect(draw(m)).not.toBe(previous);
        const expectedAtoms = bondedAtoms();
        expectedAtoms[1][flag] = true;
        move(expectedAtoms);
        const expected = model(expectedAtoms, style);
        draw(expected);
        expect(buffers(draw(m))).toEqual(buffers(draw(expected)));
        expect(atoms[1].intersectionShape).toEqual(expected.getInternalState().atoms[1].intersectionShape);
        if (flag !== "contextMenuEnabled") {
            expect(atoms[1].intersectionShape.sphere.length).toBeGreaterThan(0);
            expect(atoms[1].intersectionShape.sphere[0].center.x).toBe(atoms[1].x);
        }
    });

    test("a missing atom refuses without partially updating earlier atoms", () => {
        const m = model(bondedAtoms(), { sphere: {}, stick: {} }), previous = draw(m);
        const atoms = m.getInternalState().atoms;
        atoms[0].x += 5;
        atoms[1] = undefined;
        expectRefusal(m, previous);
    });

    test.each([
        ["ordinary imposters and lines", { sphere: {}, stick: {}, line: {} }, imposters, 1.5],
        ["hidden cross and cartoon", { sphere: {}, cross: { hidden: true }, cartoon: { hidden: true } }, imposters, 1.5],
        ["hidden tessellated sphere and stick", { line: {}, sphere: { hidden: true }, stick: { hidden: true } }, {}, 1.5],
        ["coincident single-bond lines", { line: {} }, imposters, 0]
    ])("%s still replays and matches a full rebuild", async (_name, style, options, length) => {
        const m = model(bondedAtoms(length), style), previous = draw(m, options);
        const before = buffers(previous);
        move(m.getInternalState().atoms);
        await m.setFrame(0, { fast: true });
        expect(draw(m, options)).toBe(previous);
        expect(buffers(previous)).not.toEqual(before);
        const expected = model(bondedAtoms(length), style);
        move(expected.getInternalState().atoms);
        expect(buffers(previous)).toEqual(buffers(draw(expected, options)));
    });
});
