global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

let viewers, containers, frames, nextFrame;
const event = (extra = {}) => ({ preventDefault() {}, pageX: 200, pageY: 200, which: 1, ...extra });
const wheel = v => v._handleMouseScroll(event({ wheelDelta: 30 }));
function tick() {
    const pending = [...frames.values()];
    frames.clear();
    pending.forEach(callback => callback(16));
}
function viewer(config = {}) {
    const container = document.createElement("div");
    document.body.append(container);
    containers.push(container);
    const v = new $3Dmol.GLViewer(container, { nomouse: true, autoDestroy: false, ...config });
    viewers.push(v);
    v.getWidth = v.getHeight = () => 400;
    v.WIDTH = v.HEIGHT = 400;
    v.show();
    v.renderer.render.mockClear();
    return v;
}
function pair(config) {
    const deferred = viewer(config), immediate = viewer(config);
    immediate.requestInteractiveRender = () => immediate.show();
    return [deferred, immediate];
}
function drag(v, moves) {
    v._handleMouseDown(event());
    moves.forEach(move => v._handleMouseMove(event(move)));
}

beforeEach(() => {
    jest.useFakeTimers();
    viewers = []; containers = []; frames = new Map(); nextFrame = 0;
    jest.spyOn(window, "requestAnimationFrame").mockImplementation(callback => {
        const id = nextFrame++;
        frames.set(id, callback);
        return id;
    });
    jest.spyOn(window, "cancelAnimationFrame").mockImplementation(id => frames.delete(id));
    jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        const canvas = this.config.canvas || document.createElement("canvas");
        this.renderer = {
            domElement: canvas, getCanvas: () => canvas, getAspect: () => 1,
            getXRatio: () => 1, getYRatio: () => 1, isLost: () => false,
            setSize: jest.fn(), setClearColorHex: jest.fn(), setViewport: jest.fn(),
            supportedExtensions: () => ({}), autoUpdateScene: true, dispose: jest.fn(),
            render: jest.fn((scene, camera) => {
                scene.updateMatrixWorld();
                camera.updateMatrixWorld();
                camera.matrixWorldInverse.getInverse(camera.matrixWorld);
            })
        };
    });
});
afterEach(() => {
    viewers.forEach(v => v.destroy());
    containers.forEach(container => container.remove());
    jest.clearAllTimers();
    jest.useRealTimers();
    jest.restoreAllMocks();
});

test("wheel bursts preserve every zoom step but draw once, including RAF id zero", () => {
    const [v, immediate] = pair();
    const changed = jest.fn();
    v.setViewChangeCallback(changed);
    for (let i = 0; i < 12; i++) { wheel(v); wheel(immediate); }
    expect(v.getView()).toEqual(immediate.getView());
    expect(v.renderer.render).not.toHaveBeenCalled();
    expect(changed).not.toHaveBeenCalled();
    expect(frames.size).toBe(1);
    tick();
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
    expect(changed).toHaveBeenCalledTimes(1);
    expect(frames.size).toBe(0);
});

test.each([false, true])("rotation, zoom and pan within a frame match immediate input (orthographic: %s)", orthographic => {
    const [v, immediate] = pair({ orthographic });
    const moves = [
        { pageX: 215, pageY: 212 },
        { pageX: 225, pageY: 224, shiftKey: true },
        { pageX: 230, pageY: 229, shiftKey: true },
        { pageX: 240, pageY: 235, ctrlKey: true }
    ];
    for (const target of [v, immediate]) { wheel(target); drag(target, moves); }
    expect(v.getView()).toEqual(immediate.getView());
    expect(v.renderer.render).not.toHaveBeenCalled();
    tick();
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
    expect(v.camera.projectionMatrix.elements).toEqual(immediate.camera.projectionMatrix.elements);
    expect(v.modelGroup.matrixWorld.elements).toEqual(immediate.modelGroup.matrixWorld.elements);
});

test("touch pinch accumulates zoom steps before the single redraw", () => {
    const [v, immediate] = pair();
    for (const target of [v, immediate]) {
        target._handleMouseDown(event({ targetTouches: [{ pageX: 150, pageY: 200 }, { pageX: 250, pageY: 200 }] }));
        for (const x of [260, 275, 290]) {
            target._handleMouseMove(event({ targetTouches: [{ pageX: 150, pageY: 200 }, { pageX: x, pageY: 200 }] }));
        }
    }
    expect(v.getView()).toEqual(immediate.getView());
    expect(v.renderer.render).not.toHaveBeenCalled();
    tick();
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
});

test("picking uses pending transforms without drawing", () => {
    const [v, immediate] = pair({ orthographic: true });
    const shapes = [v, immediate].map(target => {
        const shape = target.addSphere({ center: { x: 15, y: 10, z: 0 }, radius: 5, clickable: true });
        target.render(); target.renderer.render.mockClear();
        wheel(target);
        drag(target, [{ pageX: 260, pageY: 225 }]);
        return shape;
    });
    const screen = immediate.modelToScreen({ x: 15, y: 10, z: 0 });
    const hits = v.targetedObjects(screen.x / 200 - 1, 1 - screen.y / 200, [shapes[0]]);
    expect(hits[0].clickable).toBe(shapes[0]);
    expect(v.modelGroup.matrixWorld.elements).toEqual(immediate.modelGroup.matrixWorld.elements);
    expect(v.renderer.render).not.toHaveBeenCalled();
});

test("coordinate queries synchronize once per input and reuse the result", () => {
    const [v, immediate] = pair({ orthographic: true });
    const update = jest.spyOn(v.scene, "updateMatrixWorld");
    for (let i = 0; i < 2; i++) {
        wheel(v); wheel(immediate);
        const point = { x: 10, y: 4, z: 0 };
        expect(v.modelToScreen(point)).toEqual(immediate.modelToScreen(point));
        expect(v.screenOffsetToModel(20, 30)).toEqual(immediate.screenOffsetToModel(20, 30));
        expect(v.screenToModelDistance({ x: 190, y: 210 }, point)).toEqual(immediate.screenToModelDistance({ x: 190, y: 210 }, point));
        expect(update).toHaveBeenCalledTimes(i + 1);
    }
    expect(v.renderer.render).not.toHaveBeenCalled();
});

test("programmatic camera translation following pending zoom matches immediate projection", () => {
    const [v, immediate] = pair({ orthographic: true });
    wheel(v); wheel(immediate);
    v.translate(15, 25); immediate.translate(15, 25);
    expect(v.getView()).toEqual(immediate.getView());
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
    expect(frames.size).toBe(0);
});

test("mouseup flushes before selection and cancels the queued duplicate", () => {
    const v = viewer();
    v._handleMouseDown(event());
    wheel(v);
    const click = jest.spyOn(v, "handleClickSelection").mockImplementation(() => {
        expect(v.renderer.render).toHaveBeenCalledTimes(1);
    });
    v._handleMouseUp(event());
    expect(click).toHaveBeenCalledTimes(1);
    expect(v.isDragging).toBe(false);
    expect(window.cancelAnimationFrame).toHaveBeenCalledWith(0);
    tick();
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
});

test("PNG capture flushes the latest view synchronously", () => {
    const v = viewer();
    const toDataURL = jest.spyOn(v.getCanvas(), "toDataURL").mockImplementation(() => {
        expect(v.renderer.render).toHaveBeenCalledTimes(1);
        return "data:image/png;base64,latest";
    });
    wheel(v);
    expect(v.pngURI()).toBe("data:image/png;base64,latest");
    expect(toDataURL).toHaveBeenCalledWith("image/png");
    expect(frames.size).toBe(0);
});

test("PNG capture flushes every pending viewer that shares the captured canvas", () => {
    const canvas = document.createElement("canvas");
    const a = viewer({ canvas }), b = viewer({ canvas });
    jest.spyOn(canvas, "toDataURL").mockImplementation(() => {
        expect(a.renderer.render).toHaveBeenCalledTimes(1);
        expect(b.renderer.render).toHaveBeenCalledTimes(1);
        return "data:image/png;base64,both";
    });
    wheel(a); wheel(b);
    expect(a.pngURI()).toBe("data:image/png;base64,both");
    expect(frames.size).toBe(0);
});

test.each(["render", "setView"])("explicit %s draws synchronously and consumes the queued redraw", method => {
    const v = viewer();
    wheel(v);
    if (method === "render") v.render();
    else v.setView(v.getView());
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
    tick();
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
});

test("a view callback may queue the next frame without losing it", () => {
    const v = viewer();
    const changed = jest.fn(() => { if (changed.mock.calls.length === 1) wheel(v); });
    v.setViewChangeCallback(changed);
    wheel(v); tick();
    expect(frames.size).toBe(1);
    tick();
    expect(v.renderer.render).toHaveBeenCalledTimes(2);
    expect(changed).toHaveBeenCalledTimes(2);
    expect(frames.size).toBe(0);
});

test("destroy cancels pending drawing and a callback can destroy during a flush", () => {
    const v = viewer(), render = v.renderer.render;
    wheel(v); v.destroy(); tick();
    expect(render).not.toHaveBeenCalled();
    const other = viewer(), draw = other.renderer.render;
    other.setViewChangeCallback(() => other.destroy());
    drag(other, [{ pageX: 201, pageY: 201 }]);
    expect(() => other._handleMouseUp(event())).not.toThrow();
    expect(draw).toHaveBeenCalledTimes(1);
    expect(other.isDestroyed()).toBe(true);
    expect(frames.size).toBe(0);
});

test("PNG capture retains its canvas if a draw callback destroys the viewer", () => {
    const v = viewer();
    jest.spyOn(v.getCanvas(), "toDataURL").mockReturnValue("data:image/png;base64,last");
    v.setViewChangeCallback(() => v.destroy());
    wheel(v);
    expect(v.pngURI()).toBe("data:image/png;base64,last");
    expect(v.isDestroyed()).toBe(true);
});

test("hover continuation remains immediate while a redraw is pending", () => {
    const v = viewer();
    const unhover = jest.fn();
    v.current_hover = { unhover_callback: unhover };
    v.hoverables = [];
    wheel(v);
    v._handleMouseMove(event({ pageX: 450 }));
    expect(unhover).toHaveBeenCalledTimes(1);
    expect(v.current_hover).toBeNull();
    expect(v.renderer.render).not.toHaveBeenCalled();
});

test("linking flushes pending views and shared stereo input stays synchronous", () => {
    const canvas = document.createElement("canvas");
    const a = viewer({ canvas, control_all: true }), b = viewer({ canvas, control_all: true });
    wheel(a); wheel(b);
    a.linkViewer(b); b.linkViewer(a);
    expect(frames.size).toBe(0);
    a.renderer.render.mockClear(); b.renderer.render.mockClear();
    const start = a.rotationGroup.position.z;
    wheel(a); wheel(b);
    expect(a.getView()).toEqual(b.getView());
    const wd = 30 * 600 / 630;
    const step = z => z - (150 - z) * 0.85 * wd / 400;
    expect(a.rotationGroup.position.z).toBeCloseTo(step(step(start)));
    expect(a.renderer.render).toHaveBeenCalledTimes(2);
    expect(b.renderer.render).toHaveBeenCalledTimes(2);
    expect(frames.size).toBe(0);
});

test("moving documents cancels on the old window and queues future work on the new one", () => {
    const v = viewer();
    const iframe = document.createElement("iframe");
    document.body.append(iframe); containers.push(iframe);
    const nextWindow = iframe.contentWindow;
    nextWindow.requestAnimationFrame = jest.fn(() => 42);
    nextWindow.cancelAnimationFrame = jest.fn();
    const host = nextWindow.document.createElement("div");
    nextWindow.document.body.append(host);
    wheel(v);
    v.setContainer(host);
    expect(window.cancelAnimationFrame).toHaveBeenCalledWith(0);
    expect(frames.size).toBe(0);
    expect(v.renderer.render).toHaveBeenCalledTimes(1);
    wheel(v);
    expect(nextWindow.requestAnimationFrame).toHaveBeenCalledTimes(1);
    v.destroy();
    expect(nextWindow.cancelAnimationFrame).toHaveBeenCalledWith(42);
});
