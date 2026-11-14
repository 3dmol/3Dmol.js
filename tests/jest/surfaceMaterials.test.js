global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const flush = async () => { await Promise.resolve(); await Promise.resolve(); };

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

function viewer() {
    // Initialize all viewer state and exercise real surface methods without a GPU.
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"),
            getAspect: () => 1, setSize() {}, setClearColorHex() {},
            setViewport() {}, supportedExtensions: () => ({})
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

function watch(resource) {
    const listener = jest.fn();
    resource.addEventListener("dispose", listener);
    return listener;
}

describe("surface material ownership", () => {
    test("restyling retains the displayed material until render replaces it", async () => {
        const v = viewer(), source = mesh(), id = v.addMesh(source);
        const disposed = watch(source.material), geometryDisposed = watch(source.geometry);
        v.render();
        const previous = v.getSurface(id)[0].lastGL;
        v.setSurfaceMaterialStyle(id, { color: "red" });
        await flush();
        v.show();
        expect(v.modelGroup.children).toEqual([previous]);
        expect(previous.material).toBe(source.material);
        expect(disposed).not.toHaveBeenCalled();

        v.render();
        const current = v.getSurface(id)[0];
        expect(current.lastGL.material).toBe(current.mat);
        expect(current.lastGL.material).not.toBe(source.material);
        expect(disposed).not.toHaveBeenCalled();
        await flush();
        expect(disposed).toHaveBeenCalledTimes(1);
        expect(geometryDisposed).not.toHaveBeenCalled();
        v.removeAllSurfaces();
        await flush();
    });

    test("several restyles retire unused intermediate materials exactly once", async () => {
        const v = viewer(), id = v.addMesh(mesh()), retired = [];
        v.render();
        for (let i = 0; i < 4; i++) {
            retired.push(watch(v.getSurface(id)[0].mat));
            v.setSurfaceMaterialStyle(id, { color: "red" });
        }
        const current = watch(v.getSurface(id)[0].mat);
        v.render();
        retired.forEach(listener => expect(listener).not.toHaveBeenCalled());
        await flush();
        retired.forEach(listener => expect(listener).toHaveBeenCalledTimes(1));
        expect(current).not.toHaveBeenCalled();
        v.removeAllSurfaces();
        await flush();
        expect(current).toHaveBeenCalledTimes(1);
    });

    test("removal before redraw retires both the displayed and replacement materials", async () => {
        const v = viewer(), id = v.addMesh(mesh());
        v.render();
        const previous = watch(v.getSurface(id)[0].mat);
        v.setSurfaceMaterialStyle(id, { color: "blue" });
        const replacement = watch(v.getSurface(id)[0].mat);
        v.removeSurface(id);
        expect(v.modelGroup.children).toEqual([]);
        expect(previous).not.toHaveBeenCalled();
        expect(replacement).not.toHaveBeenCalled();
        await flush();
        expect(previous).toHaveBeenCalledTimes(1);
        expect(replacement).toHaveBeenCalledTimes(1);
    });

    test("a material shared by surfaces survives until the final owner is removed", async () => {
        const v = viewer(), first = mesh(), shared = first.material;
        const a = v.addMesh(first), b = v.addMesh(mesh(shared)), disposed = watch(shared);
        v.render();
        v.setSurfaceMaterialStyle(a, { color: "red" });
        v.render();
        await flush();
        expect(disposed).not.toHaveBeenCalled();
        expect(v.getSurface(b)[0].lastGL.material).toBe(shared);
        v.removeSurface(a);
        await flush();
        expect(disposed).not.toHaveBeenCalled();
        v.removeSurface(b);
        expect(disposed).not.toHaveBeenCalled();
        await flush();
        expect(disposed).toHaveBeenCalledTimes(1);
    });

    test("an unrendered owner retains a shared material until it is restyled", async () => {
        const v = viewer(), first = mesh(), shared = first.material, disposed = watch(shared);
        const a = v.addMesh(first);
        v.render();
        const b = v.addMesh(mesh(shared));
        v.removeSurface(a);
        await flush();
        expect(disposed).not.toHaveBeenCalled();
        expect(v.getSurface(b)[0].lastGL).toBeUndefined();
        v.setSurfaceMaterialStyle(b, { color: "blue" });
        await flush();
        expect(disposed).toHaveBeenCalledTimes(1);
        v.removeAllSurfaces();
        await flush();
    });

    test("reattaching a material before cleanup preserves its new owner", async () => {
        const v = viewer(), first = mesh(), material = first.material, disposed = watch(material);
        const oldId = v.addMesh(first);
        v.render();
        v.removeSurface(oldId);
        const newId = v.addMesh(mesh(material));
        await flush();
        expect(disposed).not.toHaveBeenCalled();
        v.render();
        expect(v.getSurface(newId)[0].lastGL.material).toBe(material);
        v.removeAllSurfaces();
        await flush();
        expect(disposed).toHaveBeenCalledTimes(1);
    });

    test("a disposal listener can reattach another material awaiting retirement", async () => {
        const v = viewer(), first = mesh(), second = mesh();
        const disposedFirst = watch(first.material), disposedSecond = watch(second.material);
        v.addMesh(first);
        v.addMesh(second);
        v.render();
        let reattached;
        first.material.addEventListener("dispose", () => { reattached = v.addMesh(mesh(second.material)); });
        v.removeAllSurfaces();
        await flush();
        expect(disposedFirst).toHaveBeenCalledTimes(1);
        expect(disposedSecond).not.toHaveBeenCalled();
        expect(v.getSurface(reattached)[0].mat).toBe(second.material);
        v.render();
        expect(v.getSurface(reattached)[0].lastGL.material).toBe(second.material);
        v.removeSurface(reattached);
        await flush();
        expect(disposedSecond).toHaveBeenCalledTimes(1);
    });

    test("partial, finished, and regenerated draws retain their current material", async () => {
        const v = viewer(), id = v.addMesh(mesh()), surface = v.getSurface(id)[0];
        const disposed = watch(surface.mat);
        surface.done = false;
        for (let i = 0; i < 3; i++) {
            v.render();
            await flush();
            expect(disposed).not.toHaveBeenCalled();
        }
        surface.done = true;
        v.render();
        v.render();
        v.render(undefined, { regen: true });
        await flush();
        expect(disposed).not.toHaveBeenCalled();
        v.removeAllSurfaces();
        await flush();
        expect(disposed).toHaveBeenCalledTimes(1);
    });
});
