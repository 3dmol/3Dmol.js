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
async function fixture(options = {}) {
    const controls = [], frames = [], applied = [];
    let models = [], camera = [0, 0, 0, -30, 0, 0, 0, 1], currentFrame = 0;
    const message = { hidden: true, dataset: {}, textContent: '' };
    const win = {
        location: { href: options.url || 'https://example.test/viewer.html' },
        document: { hidden: false, getElementById: id => id === 'viewer-message' ? message : null },
        AbortController,
        history: { pushState: jest.fn(), replaceState: jest.fn() },
        setTimeout, clearTimeout, setInterval, clearInterval,
        addEventListener: jest.fn(), removeEventListener: jest.fn(), leaveViewerSession: jest.fn(),
        jQuery: {}, io: {}, initSessions: jest.fn(), joinSession: jest.fn()
    };
    const viewer = {
        getNextModelId: () => models.length,
        getModel: id => models[id],
        getView: () => camera.slice(),
        setView: value => { camera = value.slice(); },
        getSlab: () => ({ near: -50, far: 50 }),
        getFrame: () => currentFrame,
        stopAnimate: jest.fn(), zoomTo: jest.fn(), render: jest.fn(),
        setHoverable: jest.fn(), setClickable: jest.fn(), setSlab: jest.fn(),
        clear: () => { models = []; },
        setInternalState: jest.fn((internal, options) => {
            if (options?.partial && !internal.models) return;
            applied.push(internal.models[0].atoms[0].x);
            models = internal.models.map((model, id) => ({
                getID: () => id, getInternalState: () => model,
                selectedAtoms: () => model.atoms, isHidden: () => false
            }));
        }),
        setFrame: frame => {
            let resolve, reject;
            const promise = new Promise((yes, no) => { resolve = yes; reject = no; });
            frames.push({ resolve, reject });
            return promise.then(() => { currentFrame = frame; });
        }
    };
    class StateManager {
        constructor(_viewer, config) {
            this.selections = new Map(); this.surfaces = new Map(); this.modelTitles = new Map();
            this.controls = { edit: config.edit ?? true, trajectory: config.trajectory ?? true };
            this.view = {};
            this.ui = {
                element: { isConnected: true, hidden: false }, isOpen: false,
                getTrajectorySettings: () => ({ frame: currentFrame, playing: false }),
                setTrajectorySettings: async settings => { if (settings.frame !== undefined) currentFrame = settings.frame; }
            };
            controls.push(this.ui.element);
        }
        subscribe() { return () => {}; }
        batch(callback) { callback(); }
        refresh() {}
        getViewSettings() { return this.view; }
        setViewSettings(settings) { this.view = settings; }
        getControls() { return { ...this.controls }; }
        setControls(settings) { Object.assign(this.controls, settings); if (!this.controls.edit) this.ui.isOpen = false; }
        destroy() { this.ui.element.isConnected = false; }
    }
    class LearningSession {
        ready = Promise.resolve();
        join = name => win.joinSession(name);
        leave = () => win.leaveViewerSession();
        destroy = () => this.leave();
    }
    const app = createApp({ window: win, viewer, library: { StateManager, LearningSession } });
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


    test('captures compact scene state and omits transient picking labels', async () => {
        current = await fixture();
        current.viewer.getInternalState = jest.fn(() => ({ models: [], shapes: [], surfaces: [], labels: [
            { text: 'Residue', style: {} }, { text: 'Query answer', style: { _sessionTransient: true } }
        ] }));
        const snapshot = current.app.getSessionState();
        expect(current.viewer.getInternalState).toHaveBeenCalledWith({ compact: true });
        expect(snapshot.labels).toEqual([{ text: 'Residue', style: {} }]);
        expect(snapshot.frame).toBe(0);
        expect(snapshot.viewSettings.slab).toEqual({ near: -50, far: 50 });
    });

    test('camera-only polling does not request a molecular state snapshot', async () => {
        current = await fixture(); current.app.sessionChanged(null, false);
        const changed = jest.fn(); current.app.subscribeSessionChanges(changed);
        jest.advanceTimersByTime(200); changed.mockClear();
        current.viewer.setView([2, 3, 4, -50, 0, 0, 0, 1]);
        jest.advanceTimersByTime(200);
        expect(changed).not.toHaveBeenCalled();
    });

    test('disables editing immediately and restores controls when leaving', async () => {
        current = await fixture();
        expect(current.app.state.getControls()).toEqual({ edit: false, trajectory: false });
        expect(current.app.state.ui.isOpen).toBe(false);
        expect(current.app.state.ui.element.hidden).toBe(true);
        current.app.sessionChanged(null, false);
        expect(current.app.state.getControls()).toEqual({ edit: true, trajectory: true });
        expect(current.app.state.ui.element.hidden).toBe(false);
        expect(current.viewer.setClickable).toHaveBeenLastCalledWith({}, false);
    });

    test('direct session URLs restore normal controls after joining, receiving snapshots, and leaving', async () => {
        current = await fixture({ url: 'https://example.test/viewer.html?session=Classroom' });
        expect(current.win.joinSession).toHaveBeenCalledWith('Classroom');
        expect(current.app.state.getControls()).toEqual({ edit: false, trajectory: false });
        const update = current.app.acceptSessionState(sessionState(1));
        await flushMicrotasks(); current.frames[0].resolve(); await update;
        current.app.sessionChanged('Classroom', true);
        expect(current.app.state.getControls()).toEqual({ edit: false, trajectory: false });
        current.app.sessionChanged(null, false);
        expect(current.app.state.getControls()).toEqual({ edit: true, trajectory: true });
        expect(current.app.state.ui.element.hidden).toBe(false);
    });

    test('reapplies cached appearance when only compact model coordinates change', async () => {
        current = await fixture();
        const style = [{ records: [{ style: { sphere: { color: 'red' } } }], patterns: [[0]], atoms: 0, frames: [] }];
        const first = current.app.acceptSessionState({ ...sessionState(1), modelStyles: style });
        await flushMicrotasks(); current.frames[0].resolve(); await first;
        const update = current.app.acceptSessionState({ ...sessionState(2), modelStyles: style }, ['models']);
        await flushMicrotasks(); current.frames[1].resolve(); await update;
        expect(current.viewer.setInternalState).toHaveBeenLastCalledWith(
            { models: sessionState(2).models, modelStyles: style }, { partial: true });
    });

    test('appearance-only patches preserve models and the disabled controls', async () => {
        current = await fixture();
        const first = current.app.acceptSessionState(sessionState(1));
        await flushMicrotasks(); current.frames[0].resolve(); await first;
        const controls = current.app.state;
        const next = { ...sessionState(1), viewSettings: { background: '#123456', backgroundAlpha: 0.5, projection: 'orthographic', slab: { near: -5, far: 5 } } };
        await current.app.acceptSessionState(next, ['viewSettings']);
        expect(current.applied).toEqual([1]);
        expect(current.frames).toHaveLength(1);
        expect(current.app.state).toBe(controls);
        expect(controls.getControls()).toEqual({ edit: false, trajectory: false });
        expect(controls.getViewSettings()).toEqual(next.viewSettings);
        expect(current.viewer.setSlab).toHaveBeenCalledWith(-5, 5);
    });

    test('coalesces different queued partial sections while a frame is loading', async () => {
        current = await fixture();
        const initial = current.app.acceptSessionState(sessionState(1));
        await flushMicrotasks(); current.frames[0].resolve(); await initial;
        const next = { ...sessionState(1), frame: 1, viewSettings: { background: '#123456' } };
        const first = current.app.acceptSessionState(next, ['frame']);
        await flushMicrotasks();
        const second = current.app.acceptSessionState(next, ['viewSettings']);
        const third = current.app.acceptSessionState({ ...next, modelTitles: [[0, 'Instructor model']] }, ['modelTitles']);
        current.frames[1].resolve(); await Promise.all([first, second, third]);
        expect(current.app.state.getViewSettings()).toEqual(next.viewSettings);
        expect(current.app.state.modelTitles.get(0)).toBe('Instructor model');
        expect(current.applied).toEqual([1]);
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
