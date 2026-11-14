global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const flush = async () => { await Promise.resolve(); await Promise.resolve(); };

function viewer() {
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"), getAspect: () => 1,
            setSize() {}, setClearColorHex() {}, setViewport() {}, supportedExtensions: () => ({})
        };
    });
    let result;
    try {
        result = new $3Dmol.GLViewer(document.createElement("div"), { nomouse: true });
    } finally {
        setup.mockRestore();
    }
    result.show = jest.fn();
    return result;
}

function mesh(material) {
    const shape = new $3Dmol.GLShape({});
    shape.addCustom({
        vertexArr: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }],
        faceArr: [0, 1, 2]
    });
    let root;
    shape.globj({ add(object) { root = object; } });
    const result = root.children[0];
    if (material) result.material = material;
    return result;
}

function watch(resource) {
    const disposed = jest.fn();
    resource.addEventListener("dispose", disposed);
    return disposed;
}

describe("surface removal", () => {
    test("unknown IDs return the viewer without redrawing or scheduling cleanup", async () => {
        const v = viewer();
        const retire = jest.spyOn(v, "disposeUnusedSurfaceMaterials");
        for (const id of [0, 42, -1, NaN]) expect(v.removeSurface(id)).toBe(v);
        await flush();
        expect(v.show).not.toHaveBeenCalled();
        expect(retire).not.toHaveBeenCalled();
        expect(v.modelGroup.children).toEqual([]);
        expect(Object.keys(v.surfaces)).toEqual([]);
    });

    test("repeated removal before and after deferred disposal retires resources once", async () => {
        const v = viewer(), source = mesh(), id = v.addMesh(source);
        const geometryDisposed = watch(source.geometry), materialDisposed = watch(source.material);
        v.render();
        v.show.mockClear();

        expect(v.removeSurface(id)).toBe(v);
        expect(v.removeSurface(id)).toBe(v);
        expect(v.getSurface(id)).toBeUndefined();
        expect(v.modelGroup.children).toEqual([]);
        expect(geometryDisposed).toHaveBeenCalledTimes(1);
        expect(materialDisposed).not.toHaveBeenCalled();
        await flush();
        expect(materialDisposed).toHaveBeenCalledTimes(1);

        expect(v.removeSurface(id)).toBe(v);
        await flush();
        expect(v.show).toHaveBeenCalledTimes(1);
        expect(geometryDisposed).toHaveBeenCalledTimes(1);
        expect(materialDisposed).toHaveBeenCalledTimes(1);
    });

    test.each(["removeAllSurfaces", "clear"])("removing an ID after %s is harmless", async method => {
        const v = viewer(), source = mesh(), id = v.addMesh(source);
        const disposed = watch(source.material);
        v.render();
        v[method]();
        v.show.mockClear();
        expect(v.removeSurface(id)).toBe(v);
        await flush();
        expect(disposed).toHaveBeenCalledTimes(1);
        expect(v.show).not.toHaveBeenCalled();
        expect(v.modelGroup.children).toEqual([]);
    });

    test("repeated removal preserves another surface and its shared material", async () => {
        const v = viewer(), first = mesh(), second = mesh(first.material);
        const a = v.addMesh(first), b = v.addMesh(second);
        const sharedDisposed = watch(first.material), secondGeometryDisposed = watch(second.geometry);
        v.render();
        const remaining = v.getSurface(b)[0].lastGL;
        v.show.mockClear();
        v.removeSurface(a);
        expect(v.removeSurface(a)).toBe(v);
        await flush();
        expect(v.getSurface(b)[0].lastGL).toBe(remaining);
        expect(v.modelGroup.children).toEqual([remaining]);
        expect(sharedDisposed).not.toHaveBeenCalled();
        expect(secondGeometryDisposed).not.toHaveBeenCalled();
        expect(v.show).toHaveBeenCalledTimes(1);

        v.removeSurface(b);
        await flush();
        expect(sharedDisposed).toHaveBeenCalledTimes(1);
        expect(secondGeometryDisposed).toHaveBeenCalledTimes(1);
    });

    test("an existing empty surface is removed and redraws once", async () => {
        const v = viewer();
        v.surfaces[0] = [];
        expect(v.removeSurface(0)).toBe(v);
        expect(Object.prototype.hasOwnProperty.call(v.surfaces, 0)).toBe(false);
        expect(v.show).toHaveBeenCalledTimes(1);
        expect(v.removeSurface(0)).toBe(v);
        await flush();
        expect(v.show).toHaveBeenCalledTimes(1);
    });
});
