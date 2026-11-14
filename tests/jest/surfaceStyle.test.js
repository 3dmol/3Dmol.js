global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const flush = async () => { await Promise.resolve(); await Promise.resolve(); };
let viewers, workers, jobs;

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
    viewers.push(result);
    return result;
}

function mesh(wireframe = false) {
    const shape = new $3Dmol.GLShape({});
    shape.addCustom({
        vertexArr: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }],
        faceArr: [0, 1, 2]
    });
    let root;
    shape.globj({ add(object) { root = object; } });
    const result = root.children[0];
    result.material.wireframe = wireframe;
    for (const group of result.geometry.geometryGroups) {
        group.lineidx = 0;
        group.lineArray = new Uint16Array(0);
    }
    return result;
}

function pendingSurface(wireframe) {
    const v = viewer(), model = new $3Dmol.GLModel(0);
    model.addAtoms([{ x: 0, y: 0, z: 0, elem: "C", bonds: [] }]);
    v.models.push(model);
    // Control delivery timing while retaining real surface generation and merging.
    v.carveUpExtent = () => [0, 1].map(() => ({
        extent: [[0, 0, 0], [1, 1, 1]], atoms: [0], toshow: [0]
    }));
    const promise = v.addSurface("VDW", { color: "red", wireframe });
    expect(jobs).toHaveLength(2);
    return { v, promise, part: v.getSurface(promise.surfid)[0] };
}

function deliver(index) {
    const x = index * 2;
    jobs[index].worker.onmessage({ data: {
        vertices: [{ x, y: 0, z: 0, atomid: 0 }, { x: x + 1, y: 0, z: 0, atomid: 0 }, { x, y: 1, z: 0, atomid: 0 }],
        faces: [0, 1, 2]
    } });
}

describe("surface style changes", () => {
    let previousAPI, previousWorker, previousSync;
    beforeEach(() => {
        viewers = []; workers = []; jobs = [];
        previousAPI = global.$3Dmol;
        previousWorker = global.Worker;
        previousSync = $3Dmol.syncSurface;
        global.$3Dmol = { ...$3Dmol, SurfaceWorker: "blob:controlled-style-worker" };
        global.Worker = class {
            constructor() { this.terminate = jest.fn(); workers.push(this); }
            postMessage(message) { if (message.type !== -1) jobs.push({ worker: this, message }); }
        };
        $3Dmol.setSyncSurface(false);
    });
    afterEach(async () => {
        try {
            for (const v of viewers) v.removeAllSurfaces();
            await flush();
        } finally {
            for (const worker of workers) worker.terminate();
            $3Dmol.setSyncSurface(previousSync);
            global.Worker = previousWorker;
            global.$3Dmol = previousAPI;
            jest.restoreAllMocks();
        }
    });

    test("restyling replaces stored metadata without merging omitted settings", () => {
        const v = viewer(), id = v.addMesh(mesh(true)), surface = v.getSurface(id);
        const previous = surface.style;
        previous.opacity = 0.25;
        const replacement = { color: "blue" };
        v.setSurfaceMaterialStyle(id, replacement);
        expect(surface.style).toBe(replacement);
        expect(surface.style.opacity).toBeUndefined();
        expect(surface.style.wireframe).toBeUndefined();
        expect(surface[0].mat.wireframe).toBe(false);
        expect(previous).toEqual({ wireframe: true, opacity: 0.25 });
    });

    test("enabling wireframe on completed geometry preserves cached edge arrays", () => {
        const v = viewer(), id = v.addMesh(mesh()), part = v.getSurface(id)[0];
        v.render();
        const group = part.geo.geometryGroups[0];
        expect(group.lineidx).toBe(0);
        v.setSurfaceMaterialStyle(id, { color: "red", wireframe: true });
        v.render();
        expect(Array.from(group.lineArray)).toEqual([0, 1, 0, 2, 1, 2]);
        const edges = group.lineArray;
        v.render(undefined, { regen: true });
        expect(group.lineArray).toBe(edges);
    });

    test("disabling wireframe stops obsolete setup on subsequent draws", () => {
        const v = viewer(), id = v.addMesh(mesh(true)), part = v.getSurface(id)[0];
        v.render();
        const setup = jest.spyOn(part.geo, "setUpWireframe");
        v.setSurfaceMaterialStyle(id, { color: "blue", wireframe: false });
        v.render();
        v.render(undefined, { regen: true });
        expect(setup).not.toHaveBeenCalled();
        expect(part.lastGL.material.wireframe).toBe(false);
    });

    test("each component uses its material flag rather than shared style metadata", () => {
        const v = viewer(), first = mesh(true), second = mesh(false), id = v.addMesh(first);
        const surface = v.getSurface(id);
        surface.push({ ...surface[0], geo: second.geometry, mat: second.material });
        surface.style = { wireframe: false };
        v.render();
        expect(first.geometry.geometryGroups[0].lineidx).toBe(6);
        expect(second.geometry.geometryGroups[0].lineidx).toBe(0);
        surface.style.wireframe = true;
        v.render(undefined, { regen: true });
        expect(second.geometry.geometryGroups[0].lineidx).toBe(0);
    });

    test.each([false, true])("enabling wireframe covers later chunks (first chunk already received: %s)", async received => {
        const { v, promise, part } = pendingSurface(false);
        if (received) {
            deliver(0);
            expect(part.geo.geometryGroups[0].lineidx).toBe(0);
        }
        v.setSurfaceMaterialStyle(promise.surfid, { color: "red", wireframe: true });
        if (!received) deliver(0);
        deliver(1);
        await promise;
        expect(part.geo.geometryGroups.map(group => group.lineidx)).toEqual([6, 6]);
        expect(part.lastGL.material.wireframe).toBe(true);
        for (const worker of workers) expect(worker.terminate).toHaveBeenCalledTimes(1);
    });

    test("disabling wireframe during generation avoids edges on later chunks", async () => {
        const { v, promise, part } = pendingSurface(true);
        deliver(0);
        expect(part.geo.geometryGroups[0].lineidx).toBe(6);
        v.setSurfaceMaterialStyle(promise.surfid, { color: "blue", wireframe: false });
        deliver(1);
        await promise;
        expect(part.geo.geometryGroups[1].lineidx).toBe(0);
        expect(part.lastGL.material.wireframe).toBe(false);
    });
});
