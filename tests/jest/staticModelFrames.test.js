global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const fs = require("fs");
const path = require("path");
const ts = require("typescript");
const trajectoryModule = {};
new Function("exports", ts.transpileModule(
    fs.readFileSync(path.join(__dirname, "../../src/ui/trajectory.ts"), "utf8"),
    { compilerOptions: { module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2019 } }
).outputText)(trajectoryModule);
const { TrajectoryControls } = trajectoryModule;
const originalFetch = global.fetch;
const group = { add() {}, remove() {} };
const extensions = { supportsImposters: true, supportsAIA: false };
let viewers = [];

function model(coordinates) {
    const m = new $3Dmol.GLModel(0);
    m.addMolData("1\ninitial\nC 0 0 0\n", "xyz");
    if (coordinates) m.setCoordinates(coordinates, "array");
    m.setStyle({}, { sphere: { radius: 0.5 } });
    return m;
}
function draw(m) {
    m.globj(group, extensions);
    return m.renderedMolObj;
}
function viewer(models) {
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        const canvas = document.createElement("canvas");
        this.renderer = {
            domElement: canvas, getCanvas: () => canvas, getAspect: () => 1,
            getXRatio: () => 1, getYRatio: () => 1, isLost: () => false,
            setSize() {}, setClearColorHex() {}, setViewport() {}, dispose() {},
            supportedExtensions: () => ({})
        };
    });
    const container = document.createElement("div");
    document.body.appendChild(container);
    let v;
    try { v = new $3Dmol.GLViewer(container, { nomouse: true }); }
    finally { setup.mockRestore(); }
    v.models = models;
    v.render = jest.fn(() => models.forEach(draw));
    viewers.push(v);
    return v;
}
function response(x) {
    const buffer = new ArrayBuffer(56);
    new Float32Array(buffer, 44).set([x, 0, 0]);
    return { ok: true, status: 200, arrayBuffer: async () => buffer };
}

afterEach(() => {
    for (const v of viewers) {
        const container = v.container;
        v.destroy();
        container.remove();
    }
    viewers = [];
    global.fetch = originalFetch;
    jest.useRealTimers();
    jest.restoreAllMocks();
});

describe("single-frame models during trajectory updates", () => {
    test("explicit model updates refresh direct edits and picking geometry", async () => {
        const m = model();
        m.setClickable({}, true, () => {});
        const previous = draw(m), atom = m.selectedAtoms({})[0];
        atom.x = 7;
        await m.setFrame(0);
        expect(draw(m)).not.toBe(previous);
        expect(atom.intersectionShape.sphere[0].center.x).toBe(7);
    });

    test("repeating an explicit viewer frame reevaluates external color-function state", async () => {
        const m = model(), v = viewer([m]);
        let color = 0xff0000;
        m.setStyle({}, { sphere: { colorfunc: () => color } });
        await v.setFrame(0);
        const previous = draw(m);
        expect(Array.from(previous.children[0].geometry.geometryGroups[0].colorArray.slice(0, 3))).toEqual([1, 0, 0]);
        color = 0x0000ff;
        await v.setFrame(0);
        const current = draw(m);
        expect(current).not.toBe(previous);
        expect(Array.from(current.children[0].geometry.geometryGroups[0].colorArray.slice(0, 3))).toEqual([0, 0, 1]);
    });

    test("a replacement local frame is selected and rebuilt", async () => {
        const m = model(), previous = draw(m);
        const state = m.getInternalState();
        const replacement = state.atoms.map(atom => ({ ...atom, x: 9 }));
        state.frames[0] = replacement;
        await m.setFrame(0);
        expect(m.getInternalState().atoms).toBe(replacement);
        expect(draw(m)).not.toBe(previous);
        expect(m.selectedAtoms({})[0].x).toBe(9);
    });

    test("style and picking setters still rebuild the single-frame model", async () => {
        const m = model();
        let previous = draw(m);
        m.setStyle({}, { sphere: { radius: 2 } });
        await m.setFrame(0);
        expect(draw(m)).not.toBe(previous);
        previous = draw(m);
        m.setClickable({}, true, () => {});
        await m.setFrame(0);
        expect(draw(m)).not.toBe(previous);
        expect(m.selectedAtoms({})[0].intersectionShape.sphere[0].radius).toBe(2);
    });

    test.each([false, true])("refreshes frame metadata and unit cells (metadata changed: %s)", async changed => {
        const m = model();
        const metadata = { cryst: { a: 10, b: 20, c: 30 } };
        m.modelDatas = [metadata];
        m.setModelData(changed ? {} : metadata);
        const previous = draw(m);
        m.unitCellObjects = {};
        m.viewer = { removeUnitCell: jest.fn(), addUnitCell: jest.fn() };
        await m.setFrame(3);
        expect(m.modelData).toBe(metadata);
        expect(m.viewer.removeUnitCell).toHaveBeenCalledWith(m);
        expect(m.viewer.addUnitCell).toHaveBeenCalledWith(m);
        expect(draw(m)).not.toBe(previous);
    });

    test.each([false, true])("one-frame remote models reload and ignore superseded responses (failure: %s)", async failOld => {
        const m = model(), previous = draw(m);
        const frames = m.getInternalState().frames;
        Object.assign(frames, { numFrames: 1, url: "https://trajectory.test", path: "sample" });
        const requests = [];
        global.fetch = jest.fn(() => new Promise((resolve, reject) => requests.push({ resolve, reject })));
        const older = m.setFrame(0);
        const newer = m.setFrame(0);
        expect(global.fetch).toHaveBeenCalledTimes(2);
        expect(global.fetch).toHaveBeenCalledWith("https://trajectory.test/traj/frame/0/sample", {
            method: "POST", body: undefined
        });
        requests[1].resolve(response(12));
        await newer;
        const current = draw(m);
        expect(current).not.toBe(previous);
        expect(m.selectedAtoms({})[0].x).toBe(12);
        if (failOld) requests[0].reject(new Error("superseded failure"));
        else requests[0].resolve(response(3));
        await older;
        expect(m.selectedAtoms({})[0].x).toBe(12);
        expect(draw(m)).toBe(current);
    });

    test.each(["coordinates", "style", "color-function state"])("animation refreshes direct changes to %s between frames", async change => {
        jest.useFakeTimers();
        const fixed = model(), moving = model([[[0, 0, 0]], [[5, 0, 0]]]);
        let color = 0xff0000;
        fixed.setStyle({}, { sphere: { radius: 0.5, colorfunc: () => color } });
        const v = viewer([fixed, moving]);
        v.animate({ interval: 20, reps: 1 });
        await jest.advanceTimersByTimeAsync(0);
        const previous = draw(fixed), atom = fixed.selectedAtoms({})[0];
        if (change === "coordinates") atom.x = 7;
        else if (change === "style") atom.style.sphere.radius = 2;
        else color = 0x0000ff;
        // These edits do not themselves invalidate the existing model geometry.
        expect(fixed.molObj).not.toBeNull();
        await jest.advanceTimersByTimeAsync(20);
        const current = draw(fixed), geometry = current.children[0].geometry.geometryGroups[0];
        expect(current).not.toBe(previous);
        if (change === "coordinates") expect(geometry.vertexArray[0]).toBe(7);
        else if (change === "style") expect(geometry.normalArray[0]).toBe(-2);
        else expect(Array.from(geometry.colorArray.slice(0, 3))).toEqual([0, 0, 1]);
        expect(moving.selectedAtoms({})[0].x).toBe(5);
        expect(v.render).toHaveBeenCalledTimes(2);
        expect(v.isAnimated()).toBe(false);
    });

    test("trajectory controls refresh direct edits to a single-frame model while seeking", async () => {
        jest.useFakeTimers();
        const fixed = model(), moving = model([[[0, 0, 0]], [[5, 0, 0]]]);
        let color = 0xff0000;
        fixed.setStyle({}, { sphere: { radius: 0.5, colorfunc: () => color } });
        const v = viewer([fixed, moving]), controls = new TrajectoryControls(v);
        try {
            await controls.setSettings({ frame: 1 });
            const previous = draw(fixed), atom = fixed.selectedAtoms({})[0];
            atom.x = 9;
            atom.style.sphere.radius = 2;
            color = 0x0000ff;
            await controls.setSettings({ frame: 0 });
            const current = draw(fixed), geometry = current.children[0].geometry.geometryGroups[0];
            expect(current).not.toBe(previous);
            expect(geometry.vertexArray[0]).toBe(9);
            expect(geometry.normalArray[0]).toBe(-2);
            expect(Array.from(geometry.colorArray.slice(0, 3))).toEqual([0, 0, 1]);
        } finally { controls.destroy(); }
    });

    test("explicit fast updates still reuse geometry when only coordinates change", async () => {
        const m = model([[[0, 0, 0]], [[5, 0, 0]]]), previous = draw(m);
        const geometry = previous.children[0].geometry.geometryGroups[0];
        const colors = Array.from(geometry.colorArray), faces = Array.from(geometry.faceArray);
        await m.setFrame(1, { fast: true });
        expect(draw(m)).toBe(previous);
        expect(geometry.vertexArray[0]).toBe(5);
        expect(Array.from(geometry.colorArray)).toEqual(colors);
        expect(Array.from(geometry.faceArray)).toEqual(faces);
    });
});
