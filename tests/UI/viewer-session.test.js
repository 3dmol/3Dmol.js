/** @jest-environment node */
const { createApp, serialize } = require('../../viewer.js');

const flushMicrotasks = async () => {
    await Promise.resolve();
    await Promise.resolve();
    await Promise.resolve();
};

function sessionState(x) {
    return { frame: 0, models: [{ atoms: [{ index: 0, x, y: 0, z: 0, elem: 'C', color: 0 }], frames: [] }] };
}

// Run the real URL/session coordinator with controllable frame loading. No DOM,
// canvas, molecular parser, or simulated Socket.IO transport is needed here.
async function fixture() {
    const controls = [], frames = [], applied = [];
    let models = [], camera = [0, 0, 0, -30, 0, 0, 0, 1], currentFrame = 0;
    const message = { hidden: true, dataset: {}, textContent: '' };
    const win = {
        location: { href: 'https://example.test/viewer.html' },
        document: { hidden: false, getElementById: id => id === 'viewer-message' ? message : null },
        AbortController,
        history: { pushState: jest.fn(), replaceState: jest.fn() },
        setTimeout, clearTimeout, setInterval, clearInterval,
        addEventListener: jest.fn(), removeEventListener: jest.fn(), leaveViewerSession: jest.fn()
    };
    const viewer = {
        getNextModelId: () => models.length,
        getModel: id => models[id],
        getView: () => camera.slice(),
        setView: value => { camera = value.slice(); },
        getSlab: () => ({ near: -50, far: 50 }),
        getFrame: () => currentFrame,
        stopAnimate: jest.fn(), zoomTo: jest.fn(), render: jest.fn(),
        setHoverable: jest.fn(), setClickable: jest.fn(),
        clear: () => { models = []; },
        setInternalState: internal => {
            applied.push(internal.models[0].atoms[0].x);
            models = internal.models.map((model, id) => ({
                getID: () => id, getInternalState: () => model,
                selectedAtoms: () => model.atoms, isHidden: () => false
            }));
        },
        setFrame: frame => {
            let resolve, reject;
            const promise = new Promise((yes, no) => { resolve = yes; reject = no; });
            frames.push({ resolve, reject });
            return promise.then(() => { currentFrame = frame; });
        }
    };
    class StateManager {
        constructor() {
            this.selections = new Map(); this.surfaces = new Map(); this.modelTitles = new Map();
            this.ui = {
                element: { isConnected: true, hidden: false }, isOpen: false,
                getTrajectorySettings: () => ({ frame: currentFrame, playing: false }),
                setTrajectorySettings: async settings => { currentFrame = settings.frame; }
            };
            controls.push(this.ui.element);
        }
        subscribe() { return () => {}; }
        batch(callback) { callback(); }
        refresh() {}
        getViewSettings() { return {}; }
        destroy() { this.ui.element.isConnected = false; }
    }
    const app = createApp({ window: win, viewer, library: { StateManager } });
    await app.ready;
    app.sessionChanged('Classroom', true);
    return {
        app, viewer, win, frames, applied, message,
        connectedControls: () => controls.filter(element => element.isConnected).length,
        x: () => models[0]?.getInternalState().atoms[0].x
    };
}

describe('classroom state updates', () => {
    let current;
    beforeEach(() => { jest.useFakeTimers(); });
    afterEach(() => {
        current?.app.destroy(); current = undefined;
        jest.useRealTimers();
    });

    test('waits for the newest state queued while an earlier frame is unresolved', async () => {
        current = await fixture();
        const { app, frames, applied, x, connectedControls } = current;
        const resolvedAt = [];
        const track = promise => promise.then(() => { resolvedAt.push(x()); });
        const first = track(app.acceptSessionState(sessionState(1)));
        await flushMicrotasks();
        expect(frames).toHaveLength(1);
        expect(applied).toEqual([1]);

        const second = track(app.acceptSessionState(sessionState(2)));
        await flushMicrotasks();
        const third = track(app.acceptSessionState(sessionState(3)));
        const finished = Promise.all([first, second, third]);
        await flushMicrotasks();
        expect(applied).toEqual([1]);
        expect(resolvedAt).toEqual([]);
        expect(connectedControls()).toBe(0);

        frames[0].resolve();
        await flushMicrotasks();
        expect(applied).toEqual([1, 3]);
        expect(frames).toHaveLength(2);
        expect(resolvedAt).toEqual([]);
        expect(connectedControls()).toBe(0);

        frames[1].resolve();
        await finished;
        expect(resolvedAt).toEqual([3, 3, 3]);
        expect(connectedControls()).toBe(1);
        expect(app.state.ui.element.hidden).toBe(true);
    });

    test.each(['resolve', 'reject'])('URL navigation supersedes in-flight and queued states when the old frame later %ss', async completion => {
        current = await fixture();
        const { app, viewer, win, frames, applied, message, connectedControls } = current;
        const first = app.acceptSessionState(sessionState(1));
        await flushMicrotasks();
        const second = app.acceptSessionState(sessionState(2));
        const finished = Promise.all([first, second]).then(() => 'resolved', error => {
            app.showError(error); return error;
        });
        const camera = [5, 6, 7, -80, 0, 0, 0, 1];
        const navigation = app.loadURL(serialize({ models: [], selections: [], surfaces: [], view: { camera } }, win.location.href));
        // Packets received during URL restoration must not restart the session.
        await app.acceptSessionState(sessionState(3));
        await navigation;
        const navigationState = app.state;
        expect(current.x()).toBeUndefined();
        expect(connectedControls()).toBe(1);

        frames[0][completion](completion === 'reject' ? new Error('Old frame request failed') : undefined);
        expect(await finished).toBe('resolved');
        expect(applied).toEqual([1]);
        expect(app.state).toBe(navigationState);
        expect(viewer.getView()).toEqual(camera);
        expect(connectedControls()).toBe(1);
        expect(message.hidden).toBe(true);
    });

    test('applies the newest queued state after an older frame fails', async () => {
        current = await fixture();
        const { app, frames, applied, connectedControls } = current;
        const first = app.acceptSessionState(sessionState(1));
        await flushMicrotasks();
        const second = app.acceptSessionState(sessionState(2));
        await flushMicrotasks();
        const third = app.acceptSessionState(sessionState(3));
        const finished = Promise.all([first, second, third]).then(() => 'resolved', error => error);

        frames[0].reject(new Error('Earlier frame unavailable'));
        await flushMicrotasks();
        expect(applied).toEqual([1, 3]);
        expect(frames).toHaveLength(2);
        expect(connectedControls()).toBe(0);
        frames[1].resolve();
        expect(await finished).toBe('resolved');
        expect(current.x()).toBe(3);
        expect(connectedControls()).toBe(1);
    });
});
