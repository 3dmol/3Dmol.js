global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const xyz = "1\nfixture\nC 0 0 0\n";
const response = () => ({ ok: true, text: async () => xyz });
const flush = async () => { for (let i = 0; i < 12; i++) await Promise.resolve(); };
const originalFetch = global.fetch;
let containers;

beforeEach(() => {
    containers = [];
    jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        const canvas = this.config.canvas || document.createElement("canvas");
        this.renderer = {
            domElement: canvas, getCanvas: () => canvas, getAspect: () => 1,
            setSize() {}, setClearColorHex() {}, setViewport() {},
            supportedExtensions: () => ({}), render: jest.fn(), dispose: jest.fn()
        };
    });
});

afterEach(() => {
    for (const container of containers) {
        $3Dmol.viewers[container.id]?.destroy();
        container.remove();
    }
    global.fetch = originalFetch;
    jest.restoreAllMocks();
});

function container(autoload = false) {
    const element = document.createElement("div");
    element.id = `destroy-loading-${containers.length}`;
    if (autoload) {
        element.className = "viewer_3Dmoljs";
        element.dataset.href = "https://example.invalid/first.xyz";
        element.dataset.href2 = "https://example.invalid/second.xyz";
    }
    document.body.appendChild(element);
    containers.push(element);
    return element;
}

describe("loading completion after viewer destruction", () => {
    test("autoload releases its registry entry and ignores an in-flight response", async () => {
        let finish;
        global.fetch = jest.fn(() => new Promise(resolve => { finish = resolve; }));
        const element = container(true), callback = jest.fn();
        $3Dmol.autoload(undefined, callback);
        const viewer = $3Dmol.viewers[element.id];
        const addModel = jest.spyOn(viewer, "addModel");
        viewer.destroy();
        expect($3Dmol.viewers[element.id]).toBeUndefined();
        finish(response());
        await flush();
        expect(global.fetch).toHaveBeenCalledTimes(1);
        expect(addModel).not.toHaveBeenCalled();
        expect(callback).not.toHaveBeenCalled();
        expect($3Dmol.processing_autoinit).toBe(false);
    });

    test("destroying during the autoload render prevents later zoom and callbacks", async () => {
        global.fetch = jest.fn(async () => response());
        const element = container(true), callback = jest.fn();
        delete element.dataset.href2;
        $3Dmol.autoload(undefined, callback);
        const viewer = $3Dmol.viewers[element.id];
        const zoomTo = jest.spyOn(viewer, "zoomTo");
        viewer.setViewChangeCallback(() => viewer.destroy());
        await flush();
        expect(viewer.isDestroyed()).toBe(true);
        expect(zoomTo).not.toHaveBeenCalled();
        expect(callback).not.toHaveBeenCalled();
    });

    test("live autoload still loads every model and invokes its callback", async () => {
        global.fetch = jest.fn(async () => response());
        const element = container(true), callback = jest.fn();
        $3Dmol.autoload(undefined, callback);
        await flush();
        const viewer = $3Dmol.viewers[element.id];
        expect(viewer.isDestroyed()).toBe(false);
        expect(viewer.getModel(0).selectedAtoms({})).toHaveLength(1);
        expect(viewer.getModel(1).selectedAtoms({})).toHaveLength(1);
        expect(global.fetch).toHaveBeenCalledTimes(2);
        expect(callback).toHaveBeenCalledTimes(1);
        expect(callback).toHaveBeenCalledWith(viewer);
    });

    test.each([false, true])("download preserves completion without redrawing a destroyed viewer (callback: %s)", async callbackMode => {
        let finish;
        global.fetch = jest.fn(() => new Promise(resolve => { finish = resolve; }));
        const viewer = new $3Dmol.GLViewer(container(), { nomouse: true });
        const callback = jest.fn();
        let immediate;
        const loaded = callbackMode ? new Promise(resolve => {
            immediate = $3Dmol.download("url:https://example.invalid/model.xyz", viewer, {}, model => {
                callback(model);
                resolve(model);
            });
        }) : $3Dmol.download("url:https://example.invalid/model.xyz", viewer, {});
        await flush();
        viewer.destroy();
        const zoomTo = jest.spyOn(viewer, "zoomTo");
        const render = jest.spyOn(viewer, "render");
        finish(response());
        const model = await loaded;
        expect(model.selectedAtoms({})).toHaveLength(1);
        expect(zoomTo).not.toHaveBeenCalled();
        expect(render).not.toHaveBeenCalled();
        if (callbackMode) {
            expect(model).toBe(immediate);
            expect(callback).toHaveBeenCalledTimes(1);
            expect(callback).toHaveBeenCalledWith(model);
        }
    });
});
