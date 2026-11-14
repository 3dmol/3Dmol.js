global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

let containers;

function viewer() {
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"), getAspect: () => 1,
            setSize() {}, setClearColorHex() {}, setViewport() {}, supportedExtensions: () => ({})
        };
    });
    const container = document.createElement("div");
    document.body.appendChild(container);
    containers.push(container);
    let result;
    try {
        result = new $3Dmol.GLViewer(container, { nomouse: true });
    } finally {
        setup.mockRestore();
    }
    result.WIDTH = result.HEIGHT = 400;
    result.camera.updateMatrixWorld(true);
    result.show = jest.fn();
    result.render = jest.fn();
    return result;
}

beforeEach(() => {
    containers = [];
    jest.useFakeTimers();
});

afterEach(() => {
    jest.clearAllTimers();
    jest.useRealTimers();
    for (const container of containers) container.remove();
    jest.restoreAllMocks();
});

function trajectory(v) {
    v.getNumFrames = () => 4;
    v.setFrame = jest.fn().mockResolvedValue(v);
    return v;
}

describe("camera motion and model animation accounting", () => {
    test.each([false, true])("stopping models leaves camera motion active (fixed path: %s)", async fixed => {
        const v = trajectory(viewer());
        v.zoom(2, 60, fixed);
        v.animate({ interval: 10 });
        v.animate({ interval: 10 });
        await jest.advanceTimersByTimeAsync(0);
        expect(v.setFrame).toHaveBeenCalledTimes(2);

        expect(v.stopAnimate()).toBe(v);
        expect(v.isAnimated()).toBe(true);
        v.stopAnimate();
        expect(v.isAnimated()).toBe(true);
        await jest.advanceTimersByTimeAsync(60);
        expect(v.isAnimated()).toBe(false);
        expect(v.setFrame).toHaveBeenCalledTimes(2);
        expect(v.rotationGroup.position.z).toBeCloseTo(75);
    });

    test.each([false, true])("old camera completion preserves a restarted trajectory (fixed path: %s)", async fixed => {
        const v = trajectory(viewer());
        let finishOldFrame;
        v.setFrame.mockImplementationOnce(() => new Promise(resolve => { finishOldFrame = resolve; }));
        v.zoom(2, 20, fixed);
        v.animate({ interval: 40 });
        await jest.advanceTimersByTimeAsync(0);
        v.stopAnimate();
        v.animate({ interval: 40 });
        await jest.advanceTimersByTimeAsync(0);
        finishOldFrame(v);
        await jest.advanceTimersByTimeAsync(20);
        expect(v.isAnimated()).toBe(true);
        expect(v.render).toHaveBeenCalledTimes(1);
        await jest.advanceTimersByTimeAsync(80);
        expect(v.setFrame).toHaveBeenCalledTimes(4);
        expect(v.isAnimated()).toBe(true);
        v.stopAnimate();
        expect(v.isAnimated()).toBe(false);
    });

    test("stop and restart inside render cancels all old model timers without losing camera ownership", async () => {
        const v = trajectory(viewer());
        v.render.mockImplementationOnce(() => {
            v.stopAnimate();
            v.animate({ interval: 40 });
        });
        v.zoom(2, 20);
        v.animate({ interval: 40 });
        v.animate({ interval: 40 });
        await jest.advanceTimersByTimeAsync(0);
        await jest.advanceTimersByTimeAsync(20);
        expect(v.isAnimated()).toBe(true);
        const frames = v.setFrame.mock.calls.length;
        await jest.advanceTimersByTimeAsync(80);
        expect(v.setFrame.mock.calls.length).toBeGreaterThan(frames);
        v.stopAnimate();
        expect(v.isAnimated()).toBe(false);
    });
});

describe("setView input validation", () => {
    test("invalid input leaves every view component unchanged and does not redraw", () => {
        const v = viewer();
        const initial = [1, 2, 3, -40, 0, 0, 0, 1, 5, 6];
        v.setView(initial);
        v.show.mockClear();
        const short = [1], sparse = new Array(8);
        const invalid = [undefined, null, {}, "12345678", short, sparse,
            initial.slice(0, 7), initial.slice(0, 9), [...initial, 0],
            [0, 0, 0, 0, 0, 0, 0, NaN], [0, 0, 0, 0, 0, 0, 0, Infinity],
            [0, 0, 0, 0, 0, 0, 0, "1"], [0, 0, 0, 0, 0, 0, 0, null]];
        for (const value of invalid) {
            expect(() => v.setView(value)).not.toThrow();
            expect(v.getView()).toEqual(initial.slice(0, 8));
            expect(v.rotationGroup.position.x).toBe(5);
            expect(v.rotationGroup.position.y).toBe(6);
        }
        expect(v.show).not.toHaveBeenCalled();
    });

    test("eight- and ten-component views preserve the legacy pan fields and nolink flag", () => {
        const v = viewer(), full = [1, 2, 3, -40, 0, 0, 0, 1, 5, 6];
        expect(v.setView(full, true)).toBe(v);
        expect(v.getView()).toEqual(full.slice(0, 8));
        expect(v.rotationGroup.position.x).toBe(5);
        expect(v.rotationGroup.position.y).toBe(6);
        expect(v.show).toHaveBeenLastCalledWith(true);
        const compact = [3, 2, 1, -20, 0, 0, 1, 0];
        v.setView(compact);
        expect(v.getView()).toEqual(compact);
        expect(v.rotationGroup.position.x).toBe(5);
        expect(v.rotationGroup.position.y).toBe(6);
    });
});

describe("camera and scene positioning", () => {
    test.each([false, true])("camera z updates position, projection and clipping together (orthographic: %s)", ortho => {
        const v = viewer();
        const before = Array.from(v.camera.projectionMatrix.elements);
        v.setCameraParameters({ z: 300, fov: 30, orthographic: ortho });
        expect(v.camera.position.z).toBe(300);
        expect(v.camera.position.x).toBe(0);
        expect(v.camera.fov).toBe(30);
        expect(v.camera.ortho).toBe(ortho);
        expect(v.getPerceivedDistance()).toBe(300);
        expect(v.camera.near).toBe(250);
        expect(v.camera.far).toBe(350);
        expect(v.camera.right).toBeCloseTo(300 * Math.tan(Math.PI / 6));
        expect(v.scene.fog.near).toBe(290);
        expect(v.scene.fog.far).toBe(350);
        expect(Array.from(v.camera.projectionMatrix.elements)).not.toEqual(before);
        v.setCameraParameters({ fov: 20 });
        expect(v.camera.position.z).toBe(300);
        expect(v.camera.near).toBe(250);
    });

    test.each([false, true])("animated scene translation reaches the immediate destination (fixed path: %s)", async fixed => {
        const immediate = viewer(), animated = viewer();
        for (const v of [immediate, animated]) {
            v.modelGroup.position.set(1, 2, 3);
            v.rotationGroup.quaternion.set(0, 0, Math.SQRT1_2, Math.SQRT1_2);
        }
        const initial = animated.modelGroup.position.clone();
        const rotation = animated.rotationGroup.quaternion.clone();
        const cameraTarget = animated.lookingAt.clone();
        immediate.translateScene(40, -20);
        expect(immediate.modelGroup.position.distanceTo(initial)).toBeGreaterThan(0);
        animated.translateScene(40, -20, 60, fixed);
        await jest.advanceTimersByTimeAsync(20);
        const target = immediate.modelGroup.position;
        expect(animated.modelGroup.position.distanceTo(initial)).toBeGreaterThan(0);
        expect(animated.modelGroup.position.distanceTo(target)).toBeGreaterThan(0);
        await jest.advanceTimersByTimeAsync(40);
        expect(animated.modelGroup.position.x).toBeCloseTo(target.x);
        expect(animated.modelGroup.position.y).toBeCloseTo(target.y);
        expect(animated.modelGroup.position.z).toBeCloseTo(target.z);
        for (const key of ["x", "y", "z", "w"]) {
            expect(animated.rotationGroup.quaternion[key]).toBeCloseTo(rotation[key]);
        }
        expect(animated.lookingAt).toEqual(cameraTarget);
        expect(animated.isAnimated()).toBe(false);
    });
});
