global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

let viewers, workers, jobs;
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
    viewers.push(result);
    const model = new $3Dmol.GLModel(0);
    model.addAtoms([
        { x: 0, y: 0, z: 0, elem: "C", color: 0x22aaee, bonds: [], chain: "A" },
        { x: 2, y: 0, z: 0, elem: "O", color: 0xaacc33, bonds: [], chain: "B" }
    ]);
    result.models.push(model);
    return result;
}

function pending(style, count = 2, v = viewer()) {
    const firstJob = jobs.length, firstWorker = workers.length;
    // Keep real geometry merging and coloring; control only worker delivery timing.
    v.carveUpExtent = () => Array.from({ length: count }, () => ({
        extent: [[0, 0, 0], [3, 1, 1]], atoms: [0, 1], toshow: [0, 1]
    }));
    const promise = v.addSurface("VDW", style);
    expect(jobs.length - firstJob).toBe(count);
    return {
        v, promise, part: v.getSurface(promise.surfid)[0],
        jobs: jobs.slice(firstJob), workers: workers.slice(firstWorker)
    };
}

function deliver(surface, index, atom = index % 2) {
    const x = index * 2;
    surface.jobs[index].worker.onmessage({ data: {
        vertices: [
            { x, y: 0, z: 0, atomid: atom },
            { x: x + 1, y: 0, z: 0, atomid: atom },
            { x, y: 1, z: 0, atomid: atom }
        ],
        faces: [0, 1, 2]
    } });
}

function expectColors(surface, colorAt) {
    for (const group of surface.part.geo.geometryGroups) {
        for (let i = 0; i < group.vertices; i++) {
            const offset = i * 3;
            const point = Array.from(group.vertexArray.slice(offset, offset + 3));
            const color = $3Dmol.CC.color(typeof colorAt === "function" ? colorAt(group.atomArray[i], point) : colorAt);
            expect(Array.from(group.colorArray.slice(offset, offset + 3)))
                .toEqual(Array.from(new Float32Array([color.r, color.g, color.b])));
        }
    }
    expect(surface.part.lastGL.material).toBe(surface.part.mat);
}

async function finish(surface) {
    await surface.promise;
    for (const worker of surface.workers) expect(worker.terminate).toHaveBeenCalledTimes(1);
}

async function restyle(initial, next, colorAt) {
    const surface = pending(initial);
    deliver(surface, 0);
    surface.v.setSurfaceMaterialStyle(surface.promise.surfid, next);
    deliver(surface, 1);
    await finish(surface);
    expectColors(surface, colorAt);
    return surface;
}

function volume(getVal, colorAt) {
    return {
        voldata: { getVal },
        volscheme: { range: () => [-2, 4], valueToHex: typeof colorAt === "function" ? colorAt : () => colorAt }
    };
}

describe("surface coloring while worker chunks are pending", () => {
    let previousAPI, previousWorker, previousSync, hadAPI, hadWorker;
    beforeEach(() => {
        viewers = []; workers = []; jobs = [];
        hadAPI = Object.prototype.hasOwnProperty.call(global, "$3Dmol");
        hadWorker = Object.prototype.hasOwnProperty.call(global, "Worker");
        previousAPI = global.$3Dmol;
        previousWorker = global.Worker;
        previousSync = $3Dmol.syncSurface;
        global.$3Dmol = { ...$3Dmol, SurfaceWorker: "blob:controlled-color-worker" };
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
            for (const worker of workers) if (!worker.terminate.mock.calls.length) worker.terminate();
            $3Dmol.setSyncSurface(previousSync);
            if (hadWorker) global.Worker = previousWorker;
            else delete global.Worker;
            if (hadAPI) global.$3Dmol = previousAPI;
            else delete global.$3Dmol;
            jest.restoreAllMocks();
        }
    });

    test("unchanged surfaces retain cached atom colors without repeating callbacks", async () => {
        const colorfunc = jest.fn(atom => atom.chain === "A" ? 0xff0000 : 0x00ff00);
        const surface = pending({ colorfunc });
        expect(colorfunc).toHaveBeenCalledTimes(2);
        deliver(surface, 0);
        deliver(surface, 1);
        await finish(surface);
        expectColors(surface, atom => atom.chain === "A" ? 0xff0000 : 0x00ff00);
        expect(colorfunc).toHaveBeenCalledTimes(2);
    });

    test("late chunks follow a uniform red to blue change", async () => {
        await restyle({ color: "red" }, { color: "blue" }, 0x0000ff);
    });

    test("nonprimary colors and the current material agree across chunks", async () => {
        const surface = await restyle({ color: "red" }, { color: 0x6496c8, opacity: 0.4 }, 0x6496c8);
        expect(surface.part.mat.opacity).toBe(0.4);
        expect(surface.part.mat.transparent).toBe(true);
        expect(surface.part.mat.color.getHex()).toBe(0x6496c8);
    });

    test("clearing a fixed color restores atom colors on every chunk", async () => {
        await restyle({ color: "red" }, {}, atom => atom.elem === "C" ? 0x22aaee : 0xaacc33);
    });

    test("element mappings apply to atoms first seen in late chunks", async () => {
        await restyle({ color: "red" }, { colorscheme: { C: 0x2244aa, O: 0x44aa22 } },
            atom => atom.elem === "C" ? 0x2244aa : 0x44aa22);
    });

    test("custom coloring receives late atom metadata once per actual vertex", async () => {
        const colorfunc = jest.fn(atom => atom.chain === "A" ? 0xffcc00 : 0x00ccff);
        await restyle({ color: "red" }, { colorfunc }, atom => atom.chain === "A" ? 0xffcc00 : 0x00ccff);
        expect(colorfunc).toHaveBeenCalledTimes(6);
    });

    test("switching to volumetric coloring samples current coordinates only", async () => {
        const getVal = jest.fn(x => x);
        await restyle({ color: "red" }, volume(getVal, x => x < 2 ? 0x123456 : 0xabcdef),
            (_atom, point) => point[0] < 2 ? 0x123456 : 0xabcdef);
        expect(getVal).toHaveBeenCalledTimes(6);
    });

    test("switching from volumetric coloring restores atom colors", async () => {
        await restyle(volume(() => 1, 0xff3300), {}, atom => atom.elem === "C" ? 0x22aaee : 0xaacc33);
    });

    test("replacing both volume and gradient updates late geometry", async () => {
        await restyle(volume(() => -1, 0xff3300), volume(x => x, x => x < 2 ? 0x336699 : 0x663399),
            (_atom, point) => point[0] < 2 ? 0x336699 : 0x663399);
    });

    test("fixed colors retain precedence over volumetric coloring when restyled", async () => {
        const getVal = jest.fn(() => { throw new Error("the fixed color should take precedence"); });
        await restyle({ color: "red" }, { ...volume(getVal, 0x00ff00), color: 0x6496c8 }, 0x6496c8);
        expect(getVal).not.toHaveBeenCalled();
    });

    test("obsolete volume data is not accessed after a style replacement", async () => {
        let stale = false;
        const old = volume(() => {
            if (stale) throw new Error("obsolete volume sampled");
            return 1;
        }, 0xff3300);
        const surface = pending(old);
        deliver(surface, 0);
        surface.v.setSurfaceMaterialStyle(surface.promise.surfid, { color: 0x336699 });
        stale = true;
        deliver(surface, 1);
        await finish(surface);
        expectColors(surface, 0x336699);
    });

    test("restyling before the first chunk applies to every chunk", async () => {
        const surface = pending({ color: "red" });
        surface.v.setSurfaceMaterialStyle(surface.promise.surfid, { color: "blue" });
        deliver(surface, 0);
        deliver(surface, 1);
        await finish(surface);
        expectColors(surface, 0x0000ff);
    });

    test("successive restyles replace pending colors", async () => {
        const surface = pending({ color: "red" }, 3);
        surface.v.setSurfaceMaterialStyle(surface.promise.surfid, { color: "green" });
        deliver(surface, 0);
        surface.v.setSurfaceMaterialStyle(surface.promise.surfid, { color: "blue" });
        deliver(surface, 1);
        surface.v.setSurfaceMaterialStyle(surface.promise.surfid, {
            colorfunc: atom => atom.elem === "C" ? 0x112233 : 0xaabbcc
        });
        deliver(surface, 2);
        await finish(surface);
        expectColors(surface, atom => atom.elem === "C" ? 0x112233 : 0xaabbcc);
    });

    test("simultaneous surfaces keep independent colors", async () => {
        const v = viewer(), first = pending({ color: "red" }, 2, v), second = pending({ color: "yellow" }, 2, v);
        deliver(first, 0);
        deliver(second, 0);
        v.setSurfaceMaterialStyle(first.promise.surfid, { color: "blue" });
        deliver(second, 1);
        deliver(first, 1);
        await Promise.all([finish(first), finish(second)]);
        expectColors(first, 0x0000ff);
        expectColors(second, 0xffff00);
    });
});
