global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const realSetupRenderer = $3Dmol.GLViewer.prototype.setupRenderer;

let viewers, containers, observers, registrations, previousObservers;

function viewer(config = {}) {
    const container = document.createElement("div");
    document.body.appendChild(container);
    containers.push(container);
    const v = new $3Dmol.GLViewer(container, config);
    v.WIDTH = v.HEIGHT = 400;
    viewers.push(v);
    return v;
}
const flush = async () => { await Promise.resolve(); await Promise.resolve(); };

beforeEach(() => {
    jest.useFakeTimers();
    viewers = []; containers = []; observers = [];
    previousObservers = [window.ResizeObserver, window.IntersectionObserver];
    class Observer {
        constructor(callback) {
            this.callback = callback;
            this.observe = jest.fn(); this.unobserve = jest.fn(); this.disconnect = jest.fn();
            observers.push(this);
        }
    }
    window.ResizeObserver = window.IntersectionObserver = Observer;
    registrations = [
        jest.spyOn(document.body, "addEventListener"),
        jest.spyOn(window, "addEventListener"),
        jest.spyOn(HTMLCanvasElement.prototype, "addEventListener")
    ];
    jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        const canvas = this.config.canvas || document.createElement("canvas");
        this.renderer = {
            domElement: canvas, getCanvas: () => canvas, getAspect: () => 1,
            getXRatio: () => 1, getYRatio: () => 1, isLost: () => false,
            setSize: jest.fn(), setClearColorHex: jest.fn(), setViewport: jest.fn(),
            supportedExtensions: () => ({}), render: jest.fn(), dispose: jest.fn()
        };
    });
});

afterEach(async () => {
    for (const v of viewers) if (typeof v.destroy === "function") v.destroy();
    await flush();
    // Also clean up baseline runs, where destroy() does not exist yet.
    for (const spy of registrations) {
        spy.mock.calls.forEach((args, index) => spy.mock.contexts[index].removeEventListener(...args));
    }
    for (const container of containers) container.remove();
    [window.ResizeObserver, window.IntersectionObserver] = previousObservers;
    jest.clearAllTimers();
    jest.useRealTimers();
    jest.restoreAllMocks();
});

describe("viewer destruction", () => {
    test("disposes scene resources and controls exactly once without redrawing", async () => {
        const v = viewer(), renderer = v.renderer, canvas = v.getCanvas();
        const model = new $3Dmol.GLModel(0);
        model.addAtoms([{ x: 0, y: 0, z: 0, elem: "C", bonds: [] }]);
        v.models.push(model);
        v.setStyle({ sphere: {} });
        const shape = v.addSphere({ center: { x: 2, y: 0, z: 0 }, radius: 1 });
        const label = v.addLabel("label");
        v.render();
        const modelRemove = jest.spyOn(model, "removegl"), shapeRemove = jest.spyOn(shape, "removegl");
        const labelDispose = jest.spyOn(label, "dispose"), ui = { destroy: jest.fn() };
        v.ui = ui;
        renderer.render.mockClear();
        expect(v.destroy()).toBe(v);
        expect(v.destroy()).toBe(v);
        await flush();
        expect(v.isDestroyed()).toBe(true);
        expect(renderer.dispose).toHaveBeenCalledTimes(1);
        expect(modelRemove).toHaveBeenCalledTimes(1);
        expect(shapeRemove).toHaveBeenCalledTimes(1);
        expect(labelDispose).toHaveBeenCalledTimes(1);
        expect(ui.destroy).toHaveBeenCalledTimes(1);
        expect(v.models).toEqual([]);
        expect(v.labels).toEqual([]);
        expect(v.shapes).toEqual([]);
        expect(v.ui).toBeUndefined();
        expect(canvas.isConnected).toBe(false);
        expect(canvas._3dmol_viewer).toBeUndefined();
        expect(renderer.render).not.toHaveBeenCalled();
    });

    test("detaches global and canvas events and ignores observer deliveries queued before disconnect", () => {
        const mouseup = jest.spyOn($3Dmol.GLViewer.prototype, "_handleMouseUp");
        const mousedown = jest.spyOn($3Dmol.GLViewer.prototype, "_handleMouseDown");
        const v = viewer(), canvas = v.getCanvas(), container = v.container;
        const resize = jest.spyOn(v, "resize");
        expect(observers).toHaveLength(2);
        v.destroy();
        window.dispatchEvent(new Event("resize"));
        document.body.dispatchEvent(new Event("mouseup"));
        document.body.dispatchEvent(new Event("touchend"));
        canvas.dispatchEvent(new Event("mousedown"));
        canvas.dispatchEvent(new Event("webglcontextlost"));
        for (const observer of observers) {
            expect(observer.disconnect).toHaveBeenCalledTimes(1);
            observer.callback([{ target: container, isIntersecting: true }]);
        }
        expect(resize).not.toHaveBeenCalled();
        expect(mouseup).not.toHaveBeenCalled();
        expect(mousedown).not.toHaveBeenCalled();
    });

    test.each([false, true])("cancels motion, model, spin, hover and touch timers (fixed path: %s)", async fixed => {
        const v = viewer(), renderer = v.renderer;
        v.getNumFrames = () => 4;
        v.setFrame = jest.fn().mockResolvedValue(v);
        v.zoom(2, 80, fixed);
        v.animate({ interval: 20 });
        v.spin("y");
        v.hoverables = [{}];
        const hover = jest.spyOn(v, "handleHoverSelection").mockImplementation(() => {});
        const event = { preventDefault() {}, pageX: 10, pageY: 10 };
        v._handleMouseMove(event);
        v._handleMouseDown({ ...event, targetTouches: [{ pageX: 10, pageY: 10 }] });
        expect(jest.getTimerCount()).toBeGreaterThanOrEqual(5);
        renderer.render.mockClear();
        v.destroy();
        await jest.advanceTimersByTimeAsync(2000);
        expect(jest.getTimerCount()).toBe(0);
        expect(v.isAnimated()).toBe(false);
        expect(v.setFrame).not.toHaveBeenCalled();
        expect(hover).not.toHaveBeenCalled();
        expect(renderer.render).not.toHaveBeenCalled();
    });

    test("a pending frame cannot render or restart animation after destruction", async () => {
        const v = viewer();
        let finish;
        v.getNumFrames = () => 4;
        v.setFrame = jest.fn(() => new Promise(resolve => { finish = resolve; }));
        const render = jest.spyOn(v, "render");
        v.animate({ interval: 20 });
        await jest.advanceTimersByTimeAsync(0);
        v.destroy();
        finish();
        await flush();
        await jest.advanceTimersByTimeAsync(100);
        expect(render).not.toHaveBeenCalled();
        expect(v.isAnimated()).toBe(false);
        expect(jest.getTimerCount()).toBe(0);
    });

    test.each(["state", "view"])("destruction inside a %s callback stops the active render", kind => {
        const v = viewer(), renderer = v.renderer, complete = jest.fn();
        if (kind === "state") v.setStateChangeCallback(() => v.destroy());
        else v.setViewChangeCallback(() => v.destroy());
        expect(() => v.render(complete)).not.toThrow();
        expect(v.isDestroyed()).toBe(true);
        expect(complete).not.toHaveBeenCalled();
        expect(renderer.dispose).toHaveBeenCalledTimes(1);
        expect(v.render(complete)).toBe(v);
        expect(v.resize()).toBe(v);
        expect(complete).not.toHaveBeenCalled();
    });

    test("destruction inside unhover does not start the next hover callback", () => {
        const v = viewer(), next = { hover_callback: jest.fn() };
        v.current_hover = { unhover_callback: () => v.destroy() };
        v.setHover(next);
        expect(next.hover_callback).not.toHaveBeenCalled();
        expect(v.current_hover).toBeNull();
    });

    test("unlinks both directions and drops embedded-viewer registry entries", () => {
        const a = viewer(), b = viewer(), c = viewer();
        a.linkViewer(b); c.linkViewer(a); c.linkViewer(b);
        $3Dmol.viewers.__destroyTest = a;
        a.destroy();
        expect(a.linkedViewers).toEqual([]);
        expect(c.linkedViewers).toEqual([b]);
        expect(b.linkedFrom.has(a)).toBe(false);
        expect($3Dmol.viewers.__destroyTest).toBeUndefined();
        expect(c.linkViewer(a)).toBe(c);
        expect(c.linkedViewers).toEqual([b]);
        expect(() => c.render()).not.toThrow();
    });

    test("preserves a shared grid canvas until the final viewer is destroyed", () => {
        const container = document.createElement("div");
        document.body.appendChild(container); containers.push(container);
        const grid = $3Dmol.createViewerGrid(container, { rows: 1, cols: 2 });
        const [a, b] = grid[0]; viewers.push(a, b);
        const canvas = a.getCanvas(), bRenderer = b.renderer;
        b.destroy();
        expect(grid[0]).toEqual([a, null]);
        expect(canvas.isConnected).toBe(true);
        expect(canvas._3dmol_viewer).toBe(a);
        expect(bRenderer.dispose).toHaveBeenCalledTimes(1);
        expect(() => a.resize()).not.toThrow();
        expect(() => a.render()).not.toThrow();
        a.destroy();
        expect(grid[0]).toEqual([null, null]);
        expect(canvas.isConnected).toBe(false);
    });

    test("retains a standalone caller-supplied canvas", () => {
        const canvas = document.createElement("canvas"), v = viewer({ canvas });
        const container = v.container;
        v.destroy();
        expect(canvas.parentNode).toBe(container);
        expect(canvas._3dmol_viewer).toBeUndefined();
    });

    test("context recovery retires the old renderer before final destruction retires its replacement", () => {
        const v = viewer(), oldRenderer = v.renderer, oldCanvas = v.getCanvas();
        oldRenderer.isLost = () => true;
        v.getWidth = v.getHeight = () => 400;
        v.setupRenderer = function () {
            // Exercise real replacement ownership without requiring a browser GL context.
            this.WIDTH = this.HEIGHT = 0;
            realSetupRenderer.call(this);
        };
        v.resize();
        const replacement = v.renderer, dispose = jest.spyOn(replacement, "dispose");
        expect(replacement).not.toBe(oldRenderer);
        expect(oldRenderer.dispose).toHaveBeenCalledTimes(1);
        expect(oldCanvas._3dmol_viewer).toBeUndefined();
        expect(v.getCanvas().isConnected).toBe(true);
        v.destroy();
        expect(oldRenderer.dispose).toHaveBeenCalledTimes(1);
        expect(dispose).toHaveBeenCalledTimes(1);
    });

    test("settles an unfinished APNG capture without restoring old callbacks", async () => {
        const v = viewer(), oldCallback = jest.fn();
        v.setViewChangeCallback(oldCallback);
        const capture = v.apngURI(3);
        v.destroy();
        await expect(capture).resolves.toBeUndefined();
        expect(v.viewChangeCallback).toBeNull();
        expect(oldCallback).not.toHaveBeenCalled();
        expect(v.captureCancellations.size).toBe(0);
    });
});
