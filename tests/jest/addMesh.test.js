global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function fixture({ wireframe = false } = {}) {
    // Obtain real mesh/geometry/material instances from the public shape API.
    const shape = new $3Dmol.GLShape({});
    shape.addCustom({
        vertexArr: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }],
        faceArr: [0, 1, 2]
    });
    let root;
    shape.globj({ add(object) { root = object; } });
    const mesh = root.children[0];
    const Mesh = mesh.constructor;
    // External meshes may supply face indices without precomputed wireframe edges.
    mesh.geometry.geometryGroups[0].lineidx = 0;
    mesh.geometry.geometryGroups[0].lineArray = new Uint16Array(0);
    mesh.material.wireframe = wireframe;

    // Initialize the real viewer and surface render path without a WebGL context.
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"),
            getAspect: () => 1, setSize() {}, setClearColorHex() {},
            setViewport() {}, supportedExtensions: () => ({})
        };
    });
    let viewer;
    try {
        viewer = new $3Dmol.GLViewer(document.createElement("div"), { nomouse: true });
    } finally {
        setup.mockRestore();
    }
    viewer.show = () => {};
    return { viewer, mesh, Mesh };
}

describe("explicit mesh surfaces", () => {
    test("first render creates one mesh and later renders preserve its resources", () => {
        const { viewer, mesh, Mesh } = fixture();
        const id = viewer.addMesh(mesh);
        viewer.render();
        const first = viewer.getSurface(id)[0].lastGL;
        expect(first).toBeInstanceOf(Mesh);
        expect(viewer.modelGroup.children).toEqual([first]);
        expect(first.geometry).toBe(mesh.geometry);
        expect(first.material).toBe(mesh.material);
        expect(first.geometry.geometryGroups[0].lineidx).toBe(0);

        viewer.render();
        expect(viewer.getSurface(id)[0].lastGL).toBe(first);
        viewer.render(undefined, { regen: true });
        const regenerated = viewer.getSurface(id)[0].lastGL;
        expect(regenerated).not.toBe(first);
        expect(viewer.modelGroup.children).toEqual([regenerated]);
        expect(regenerated.geometry).toBe(mesh.geometry);
        expect(regenerated.material).toBe(mesh.material);
    });

    test("a supplied wireframe material builds the triangle's edge indices", () => {
        const { viewer, mesh } = fixture({ wireframe: true });
        const group = mesh.geometry.geometryGroups[0];
        expect(group.lineidx).toBe(0);
        const id = viewer.addMesh(mesh);
        viewer.render();
        expect(viewer.getSurface(id)[0].lastGL.material).toBe(mesh.material);
        expect(mesh.material.wireframe).toBe(true);
        expect(group.lineidx).toBe(6);
        expect(Array.from(group.lineArray)).toEqual([0, 1, 0, 2, 1, 2]);
    });

    test("restyling replaces the material while retaining the supplied geometry", () => {
        const { viewer, mesh } = fixture();
        const id = viewer.addMesh(mesh);
        viewer.render();
        viewer.setSurfaceMaterialStyle(id, { color: "red", opacity: 0.5, wireframe: true });
        viewer.render();
        const drawn = viewer.getSurface(id)[0].lastGL;
        expect(viewer.modelGroup.children).toEqual([drawn]);
        expect(drawn.geometry).toBe(mesh.geometry);
        expect(drawn.material).not.toBe(mesh.material);
        expect(drawn.material.opacity).toBe(0.5);
        expect(drawn.material.wireframe).toBe(true);
    });

    test("distinct surface IDs can be removed independently", () => {
        const { viewer, mesh } = fixture();
        const first = viewer.addMesh(mesh);
        const second = viewer.addMesh(fixture().mesh);
        expect(second).not.toBe(first);
        viewer.render();
        const remaining = viewer.getSurface(second)[0].lastGL;
        viewer.removeSurface(first);
        expect(viewer.getSurface(first)).toBeUndefined();
        expect(viewer.modelGroup.children).toEqual([remaining]);
        viewer.removeAllSurfaces();
        expect(viewer.modelGroup.children).toEqual([]);
    });

    test("removing an unrendered mesh leaves no scene object", () => {
        const { viewer, mesh } = fixture();
        const id = viewer.addMesh(mesh);
        viewer.removeSurface(id);
        viewer.render();
        expect(viewer.getSurface(id)).toBeUndefined();
        expect(viewer.modelGroup.children).toEqual([]);
    });
});
