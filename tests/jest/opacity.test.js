
global.$ = require("jquery");
global.URL.createObjectURL = function () { };
let $3Dmol = require("../../build/3Dmol.js");

// Per-atom sphere opacity (issue #166) -- structural tests.
//
// Opacity is a visual trait, but pixel comparison is a poor way to test one: a reference
// image is hostage to whichever renderer produced it, so it fails for reasons unrelated to
// the code. What IS deterministic is the geometry the model builds, so these assert on
// structure instead of appearance:
//
//   1. no opacity styled          -> nothing new is allocated; the pre-feature shape
//   2. one opacity for all atoms  -> the legacy material path, exactly as before
//   3. atoms disagreeing          -> the translucent twin appears, carrying per-vertex alpha
//
// Case 2 is the regression guard: sphere opacity is a pre-existing whole-model style, so a
// uniform value must keep behaving as it always has. Case 3 is the feature itself.
//
// The model is driven directly rather than through a viewer -- no canvas, no WebGL context,
// no mocking -- and globj() is handed supportsImposters explicitly, the same way exportVRML
// pins its options, so the imposter path is exercised deterministically.

const PDB = [
    "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00 90.00           N",
    "ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00 90.00           C",
    "ATOM      3  N   ALA A   2      8.000   0.000   0.000  1.00 40.00           N",
    "ATOM      4  CA  ALA A   2      9.500   0.000   0.000  1.00 40.00           C",
    "END",
].join("\n");

// globj() only ever calls add/remove on the group it is given.
const stubGroup = () => ({ add() { }, remove() { } });

function build(styleFn) {
    const model = new $3Dmol.GLModel(0, {});
    model.addMolData(PDB, "pdb", {});
    styleFn(model);
    model.globj(stubGroup(), { supportsImposters: true, supportsAIA: false, regen: true });
    return model.molObj;
}

const groupsOf = (geo) => (geo.geometryGroups || []).filter((g) => g && g.vertices > 0);
// Meshes whose geometry carries per-vertex alpha: the translucent twin, and only it.
const alphaMeshes = (obj) => obj.children.filter((c) => c.geometry && c.geometry.alpha);
const plainMeshes = (obj) => obj.children.filter((c) => c.geometry && !c.geometry.alpha);
const anyAlphaArray = (obj) =>
    plainMeshes(obj).some((m) => groupsOf(m.geometry).some((g) => g.alphaArray));

describe("per-atom sphere opacity (#166)", () => {

    test("no opacity styled: no twin geometry and no alpha buffer", () => {
        const obj = build((m) => m.setStyle({}, { sphere: {} }));

        expect(alphaMeshes(obj)).toHaveLength(0);
        expect(anyAlphaArray(obj)).toBe(false);
    });

    // Swept rather than spot-checked: the legacy path has to hold across the range, and a
    // single value could pass by luck. Each case also asserts the molecule renders at ONE
    // opacity -- if the gate ever leaked, some atoms would divert to the twin and the
    // remaining material value would disagree with what was asked for.
    test.each([0.1, 0.25, 0.4, 0.5, 0.75, 0.9])(
        "one opacity for every atom (%p): legacy material path, consistent across the molecule",
        (value) => {
            const obj = build((m) => m.setStyle({}, { sphere: { opacity: value } }));

            // Every atom agrees, so this is indistinguishable from the pre-existing
            // whole-model API and must stay on it: the value rides the material, not a
            // per-vertex buffer.
            expect(alphaMeshes(obj)).toHaveLength(0);
            expect(anyAlphaArray(obj)).toBe(false);

            const translucent = plainMeshes(obj).filter((m) => m.material && m.material.transparent);
            expect(translucent.length).toBeGreaterThan(0);

            // One value, molecule-wide -- not merely "the first mesh happens to be right".
            const distinct = [...new Set(translucent.map((m) => m.material.opacity))];
            expect(distinct).toHaveLength(1);
            expect(distinct[0]).toBeCloseTo(value);
        }
    );

    test("opacity 1.0 for every atom: fully opaque, nothing marked transparent", () => {
        const obj = build((m) => m.setStyle({}, { sphere: { opacity: 1.0 } }));

        expect(alphaMeshes(obj)).toHaveLength(0);
        expect(anyAlphaArray(obj)).toBe(false);
        // 1.0 is not translucency: the material must not be flagged transparent, or the
        // model would be routed through the blending path for nothing.
        expect(plainMeshes(obj).some((m) => m.material && m.material.transparent)).toBe(false);
    });

    test("atoms disagreeing: twin geometry carries per-vertex alpha", () => {
        const obj = build((m) => {
            m.setStyle({}, { sphere: { opacity: 1.0 } });
            m.setStyle({ resi: 2 }, { sphere: { opacity: 0.25 } });
        });

        const twins = alphaMeshes(obj);
        expect(twins).toHaveLength(1);
        const twin = twins[0];

        // material.opacity stays 1: the vertices own the alpha, and setting both would fade
        // every atom twice.
        expect(twin.material.transparent).toBe(true);
        expect(twin.material.opacity).toBeCloseTo(1.0);
        // No depth writes -- back-to-front ordering owns ghost-vs-ghost, and a near ghost
        // stamping the depth buffer would discard the far ghosts behind it.
        expect(twin.material.depthWrite).toBe(false);

        const values = [];
        for (const g of groupsOf(twin.geometry)) {
            expect(g.alphaArray).toBeTruthy();
            values.push(...Array.from(g.alphaArray.slice(0, g.vertices)));
        }
        // Only the residue styled 0.25 was rerouted, so every vertex in the twin carries it.
        expect(values.length).toBeGreaterThan(0);
        for (const v of values) expect(v).toBeCloseTo(0.25);

        // ...and the opaque atoms stayed behind in the ordinary sphere geometry.
        expect(plainMeshes(obj).some((m) => groupsOf(m.geometry).length > 0)).toBe(true);
    });
});
