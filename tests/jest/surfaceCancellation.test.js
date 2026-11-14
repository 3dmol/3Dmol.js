global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

let viewers, workers, jobs;
const flush = async () => { for (let i = 0; i < 5; i++) await Promise.resolve(); };
const chunk = () => ({
    vertices: [
        { x: 0, y: 0, z: 0, atomid: 0 },
        { x: 1, y: 0, z: 0, atomid: 0 },
        { x: 0, y: 1, z: 0, atomid: 0 }
    ], faces: [0, 1, 2]
});

function viewer(symmetry = false) {
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"), getAspect: () => 1,
            setSize() {}, setClearColorHex() {}, setViewport() {}, supportedExtensions: () => ({}),
            dispose: jest.fn()
        };
    });
    let v;
    try {
        v = new $3Dmol.GLViewer(document.createElement("div"), { nomouse: true });
    } finally {
        setup.mockRestore();
    }
    v.show = jest.fn();
    viewers.push(v);
    const model = new $3Dmol.GLModel(0);
    model.addAtoms([{ x: 0, y: 0, z: 0, elem: "C", color: 0x22aaee, bonds: [] }]);
    if (symmetry) model.setSymmetries([new $3Dmol.Matrix4(), new $3Dmol.Matrix4()]);
    v.models.push(model);
    v.carveUpExtent = () => [0, 1, 2].map(() => ({
        extent: [[0, 0, 0], [2, 2, 2]], atoms: [0], toshow: [0]
    }));
    return v;
}

function pending(v = viewer(), callback) {
    const firstWorker = workers.length, firstJob = jobs.length;
    const result = v.addSurface("VDW", { color: "red" }, {}, undefined, undefined, callback);
    const id = callback ? result : result.surfid;
    return {
        v, result, id, surface: v.getSurface(id),
        workers: workers.slice(firstWorker), jobs: jobs.slice(firstJob)
    };
}

function deliver(job) {
    job.handler({ data: chunk() });
}

function expectStopped(surface) {
    for (const worker of surface.workers) {
        expect(worker.terminate).toHaveBeenCalledTimes(1);
        expect(worker.onmessage).toBeNull();
        expect(worker.onerror).toBeNull();
    }
}

describe("canceling surface generation", () => {
    let previousAPI, previousWorker, previousSync;
    beforeEach(() => {
        jest.useFakeTimers();
        viewers = []; workers = []; jobs = [];
        previousAPI = global.$3Dmol;
        previousWorker = global.Worker;
        previousSync = $3Dmol.syncSurface;
        global.$3Dmol = { ...$3Dmol, SurfaceWorker: "blob:controlled-cancellation-worker" };
        global.Worker = class {
            constructor() { this.terminate = jest.fn(); workers.push(this); }
            postMessage(message) {
                if (message.type !== -1) jobs.push({ worker: this, handler: this.onmessage, error: this.onerror });
            }
        };
        $3Dmol.setSyncSurface(false);
    });
    afterEach(async () => {
        try {
            for (const v of viewers) if (!v.isDestroyed()) v.removeAllSurfaces();
            await flush();
        } finally {
            $3Dmol.setSyncSurface(previousSync);
            global.Worker = previousWorker;
            global.$3Dmol = previousAPI;
            jest.useRealTimers();
            jest.restoreAllMocks();
        }
    });

    test.each(["removeSurface", "removeAllSurfaces", "clear", "destroy"])(
        "%s stops workers, settles the Promise, and ignores queued delivery", async method => {
            const surface = pending(), part = surface.surface[0];
            const disposed = jest.fn();
            part.geo.addEventListener("dispose", disposed);
            const render = jest.spyOn(surface.v, "render");
            setTimeout(() => deliver(surface.jobs[0]), 10);
            setTimeout(() => surface.jobs[1].error({ message: "obsolete failure" }), 20);
            expect(surface.v[method](surface.id)).toBe(surface.v);
            await expect(surface.result).resolves.toBe(surface.id);
            expectStopped(surface);
            expect(disposed).toHaveBeenCalledTimes(1);
            render.mockClear();
            jest.runOnlyPendingTimers();
            await flush();
            expect(render).not.toHaveBeenCalled();
            expect(surface.v.getSurface(surface.id)).toBeUndefined();
            expect(part.geo.geometryGroups).toHaveLength(0);
        }
    );

    test("old deliveries cannot modify a replacement with the same numeric ID", async () => {
        const old = pending();
        old.v.removeSurface(old.id);
        const next = pending(old.v);
        expect(next.id).toBe(old.id);
        const render = jest.spyOn(old.v, "render");
        for (const job of old.jobs) deliver(job);
        expect(render).not.toHaveBeenCalled();
        expect(next.surface[0].geo.geometryGroups).toHaveLength(0);
        for (const worker of next.workers) expect(worker.terminate).not.toHaveBeenCalled();
        for (const job of next.jobs) deliver(job);
        await expect(next.result).resolves.toBe(next.id);
        expect(next.surface[0].done).toBe(true);
        expect(next.surface[0].geo.geometryGroups).toHaveLength(3);
        expectStopped(old);
        expectStopped(next);
    });

    test("removing all surfaces preserves the registry and replacements created during disposal", async () => {
        const old = pending(), registry = old.v.surfaces;
        let next;
        old.surface[0].geo.addEventListener("dispose", () => { next = pending(old.v); });
        old.v.removeAllSurfaces();
        expect(old.v.surfaces).toBe(registry);
        expect(next.id).toBe(old.id);
        expect(old.v.getSurface(next.id)).toBe(next.surface);
        await expect(old.result).resolves.toBe(old.id);
        for (const job of old.jobs) deliver(job);
        expect(next.surface[0].geo.geometryGroups).toHaveLength(0);
        for (const worker of next.workers) expect(worker.terminate).not.toHaveBeenCalled();
        for (const job of next.jobs) deliver(job);
        await expect(next.result).resolves.toBe(next.id);
        expect(next.surface[0].done).toBe(true);
        expectStopped(old);
        expectStopped(next);
    });

    test.each([false, true])("ordinary cancellation preserves callback shape (symmetry=%s)", async symmetry => {
        const callback = jest.fn(), surface = pending(viewer(symmetry), callback);
        expect(surface.result).toBe(surface.id);
        surface.v.removeSurface(surface.id);
        await flush();
        expect(callback).toHaveBeenCalledTimes(1);
        expect(callback).toHaveBeenCalledWith(symmetry ? [surface.id] : surface.id);
        for (const job of surface.jobs) deliver(job);
        await flush();
        expect(callback).toHaveBeenCalledTimes(1);
    });

    test("destruction suppresses a queued surface callback", async () => {
        const callback = jest.fn(), surface = pending(viewer(), callback);
        surface.v.destroy();
        await flush();
        expect(callback).not.toHaveBeenCalled();
        expectStopped(surface);
    });

    test("removal from a render callback stops remaining worker chunks", async () => {
        const surface = pending();
        surface.v.show.mockImplementationOnce(() => surface.v.removeSurface(surface.id));
        deliver(surface.jobs[0]);
        await expect(surface.result).resolves.toBe(surface.id);
        expectStopped(surface);
        const render = jest.spyOn(surface.v, "render");
        deliver(surface.jobs[1]);
        expect(render).not.toHaveBeenCalled();
    });

    test("removal from a coloring callback disposes the unmerged temporary geometry", async () => {
        const surface = pending(), part = surface.surface[0];
        surface.v.setSurfaceMaterialStyle(surface.id, {
            colorfunc: () => { surface.v.removeSurface(surface.id); return 0xff0000; }
        });
        const dispose = jest.spyOn(Object.getPrototypeOf(part.geo), "dispose");
        deliver(surface.jobs[0]);
        await expect(surface.result).resolves.toBe(surface.id);
        expectStopped(surface);
        expect(part.geo.geometryGroups).toHaveLength(0);
        expect(dispose.mock.contexts).toContain(part.geo);
        expect(dispose.mock.contexts.some(geo => geo !== part.geo)).toBe(true);
    });

    test("initial coloring can cancel before any worker is created", async () => {
        const v = viewer();
        const promise = v.addSurface("VDW", { colorfunc: () => { v.clear(); return 0xff0000; } });
        await expect(promise).resolves.toBe(promise.surfid);
        expect(workers).toHaveLength(0);
        expect(v.getSurface(promise.surfid)).toBeUndefined();
    });

    test("synchronous geometry is still immediate and removal cancels its completion redraw", async () => {
        $3Dmol.setSyncSurface(true);
        jest.spyOn($3Dmol.GLViewer, "generateMeshSyncHelper").mockImplementation(chunk);
        const v = viewer(), promise = v.addSurface("VDW", { color: "red" });
        expect(v.getSurface(promise.surfid)[0].geo.geometryGroups).toHaveLength(3);
        v.removeSurface(promise.surfid);
        const render = jest.spyOn(v, "render");
        await expect(promise).resolves.toBe(promise.surfid);
        expect(render).not.toHaveBeenCalled();
    });

    test("synchronous chunks stop when a render callback clears the scene", async () => {
        $3Dmol.setSyncSurface(true);
        const generate = jest.spyOn($3Dmol.GLViewer, "generateMeshSyncHelper").mockImplementation(chunk);
        const v = viewer();
        v.show.mockImplementationOnce(() => v.clear());
        const promise = v.addSurface("VDW", { color: "red" });
        await expect(promise).resolves.toBe(promise.surfid);
        expect(generate).toHaveBeenCalledTimes(1);
        expect(v.getSurface(promise.surfid)).toBeUndefined();
    });

    test("a normal worker error still rejects once and releases all workers", async () => {
        const surface = pending(), error = { message: "failed", filename: "surface.js", lineno: 1 };
        jest.spyOn(console, "log").mockImplementation(() => {});
        const result = expect(surface.result).rejects.toBe(error);
        surface.jobs[0].error(error);
        await result;
        expectStopped(surface);
        surface.v.removeSurface(surface.id);
        expectStopped(surface);
    });
});
