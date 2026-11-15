/** @jest-environment node */
jest.mock('canvas', () => ({}));
const { JSDOM } = require('jsdom');
const { LearningSession } = require('../../build/3Dmol.ui.js');

function fixture() {
    const dom = new JSDOM('<div id="viewer"></div>', { url: 'https://classroom.test/lesson' });
    const handlers = new Map();
    const socket = {
        connected: true, on: jest.fn((event, fn) => handlers.set(event, fn)),
        off: jest.fn((event, fn) => { if (handlers.get(event) === fn) handlers.delete(event); }),
        emit: jest.fn(), disconnect: jest.fn()
    };
    let data = { models: [{ atoms: [{ elem: 'C', x: 0, y: 0, z: 0 }], frames: [] }],
        modelStyles: [], shapes: [], surfaces: [], labels: [], frame: 0,
        viewSettings: { background: '#ffffff' } };
    let view = [0, 0, 0, -30, 0, 0, 0, 1];
    let stateChange, sceneChange, adapterChange;
    const unsubscribeState = jest.fn(), unsubscribeScene = jest.fn(), unsubscribeAdapter = jest.fn();
    const viewer = {
        container: dom.window.document.querySelector('#viewer'),
        getView: () => view, getFrame: () => data.frame, getSlab: () => ({ near: -50, far: 50 }),
        subscribeStateChanges: fn => { sceneChange = fn; return unsubscribeScene; },
        setView: jest.fn(), setClickable: jest.fn(), setHoverable: jest.fn(), selectedAtoms: () => [],
        removeLabels: jest.fn(), render: jest.fn()
    };
    const state = {
        viewer, subscribe: fn => { stateChange = fn; return unsubscribeState; },
        getViewSettings: () => data.viewSettings, setDistancePicking: jest.fn(),
        getControls: () => ({ edit: true, trajectory: true }), setControls: jest.fn(),
        ui: { element: dom.window.document.createElement('div'), getTrajectorySettings: () => ({ playing: false }) }
    };
    const adapter = {
        getState: jest.fn(() => data), applyState: jest.fn(async () => {}), roleChanged: jest.fn(), showError: jest.fn(),
        subscribe: fn => { adapterChange = fn; return unsubscribeAdapter; }
    };
    const session = new LearningSession(state, { session: 'Classroom', socketFactory: () => socket, adapter });
    const receive = (event, message) => handlers.get(event)?.(message);
    return {
        session, socket, adapter, state, viewer, receive,
        unsubscribeState, unsubscribeScene, unsubscribeAdapter,
        get data() { return data; }, set data(value) { data = value; },
        set view(value) { view = value; },
        change() { stateChange(); sceneChange(); adapterChange(); },
        updates: () => socket.emit.mock.calls.filter(([event]) => event === 'viewer state change event').map(([, packet]) => packet),
        async create() { await session.start(); receive('create session response', 1234); },
        async join() { await session.join(); receive('join session response', 'student-1'); },
        close() { session.destroy(); dom.window.close(); }
    };
}
const settle = async () => { for (let i = 0; i < 8; i++) await Promise.resolve(); };

describe('UI classroom state transport', () => {
    let current;
    beforeEach(() => { jest.useFakeTimers(); current = fixture(); });
    afterEach(() => { current.close(); jest.useRealTimers(); });

    test('coalesces edits, sends only changed sections, and suppresses unchanged renders', async () => {
        await current.create(); jest.advanceTimersByTime(200);
        expect(current.updates()).toHaveLength(0);
        current.adapter.getState.mockClear();
        current.data.viewSettings = { background: '#000000' };
        for (let i = 0; i < 20; i++) current.change();
        expect(current.adapter.getState).not.toHaveBeenCalled();
        jest.advanceTimersByTime(100);
        expect(current.adapter.getState).toHaveBeenCalledTimes(1);
        expect(current.updates()).toEqual([{ name: 'Classroom', baseRevision: 0, revision: 1,
            patch: { viewSettings: { background: '#000000' } } }]);
        current.receive('viewer state revision', { revision: 1 });
        for (let i = 0; i < 10; i++) current.change();
        jest.advanceTimersByTime(200);
        expect(current.updates()).toHaveLength(1);
    });

    test('serializes revisions and flushes the latest edit after acknowledgement', async () => {
        await current.create(); jest.advanceTimersByTime(200);
        for (const frame of [1, 2, 3]) {
            current.data.frame = frame; current.change(); jest.advanceTimersByTime(100);
        }
        expect(current.updates()).toEqual([{ name: 'Classroom', baseRevision: 0, revision: 1, patch: { frame: 1 } }]);
        current.receive('viewer state revision', { revision: 1 }); jest.advanceTimersByTime(100);
        expect(current.updates()[1]).toEqual({ name: 'Classroom', baseRevision: 1, revision: 2, patch: { frame: 3 } });
    });

    test('camera changes do not serialize or send molecular state', async () => {
        await current.create(); jest.advanceTimersByTime(200);
        current.adapter.getState.mockClear(); current.socket.emit.mockClear();
        const view = [1, 2, 3, -40, 0, 0, 0, 1];
        current.view = view; jest.advanceTimersByTime(500);
        expect(current.socket.emit.mock.calls.filter(([event]) => event === 'viewer view change event')).toEqual([
            ['viewer view change event', { name: 'Classroom', view }]
        ]);
        expect(current.adapter.getState).not.toHaveBeenCalled();
        expect(current.updates()).toHaveLength(0);
    });

    test('students merge patches and request a full snapshot after a revision gap', async () => {
        await current.join();
        current.receive('viewer state change response', { state: current.data, revision: 4 }); await settle();
        current.receive('viewer state change response', { patch: { frame: 2 }, baseRevision: 4, revision: 5 }); await settle();
        expect(current.adapter.applyState).toHaveBeenLastCalledWith({ ...current.data, frame: 2 }, ['frame']);
        current.receive('viewer state change response', { patch: { frame: 3 }, baseRevision: 6, revision: 7 }); await settle();
        expect(current.adapter.applyState).toHaveBeenCalledTimes(2);
        expect(current.socket.emit).toHaveBeenCalledWith('request state', { name: 'Classroom' });
        current.receive('viewer state change response', { state: { ...current.data, frame: 3 }, revision: 7 }); await settle();
        expect(current.adapter.applyState).toHaveBeenLastCalledWith({ ...current.data, frame: 3 }, undefined);
    });

    test('reconnection sends a fresh full state using the server revision', async () => {
        await current.create(); jest.advanceTimersByTime(200); current.socket.emit.mockClear();
        current.socket.connected = false; current.receive('disconnect', 'transport close');
        current.data.frame = 3; current.change(); jest.advanceTimersByTime(200);
        expect(current.updates()).toHaveLength(0);
        current.socket.connected = true; current.receive('connect');
        expect(current.socket.emit.mock.calls.filter(([event]) => event === 'reconnection')).toEqual([
            ['reconnection', { name: 'Classroom', secret: 1234, sid: null }]
        ]);
        current.receive('reconnection response', { controller: true, sid: 'new-instructor', revision: 8 });
        jest.advanceTimersByTime(100);
        expect(current.updates()).toEqual([{ name: 'Classroom', baseRevision: 8, revision: 9, state: current.data }]);
    });

    test('leave cancels queued sends and detaches scene and UI subscriptions', async () => {
        await current.create(); current.data.frame = 3; current.change(); current.session.leave();
        jest.advanceTimersByTime(500);
        expect(current.unsubscribeState).toHaveBeenCalledTimes(1);
        expect(current.unsubscribeScene).toHaveBeenCalledTimes(1);
        expect(current.unsubscribeAdapter).toHaveBeenCalledTimes(1);
        expect(current.updates()).toHaveLength(0);
    });

    test('a new join applies its snapshot after an old asynchronous update finishes', async () => {
        await current.join();
        let complete;
        current.adapter.applyState.mockImplementationOnce(() => new Promise(resolve => { complete = resolve; }));
        current.receive('viewer state change response', { state: current.data, revision: 0 }); await settle();
        expect(complete).toBeDefined();
        current.session.leave(); await current.join();
        const next = { ...current.data, frame: 4 };
        current.receive('viewer state change response', { state: next, revision: 0 }); await settle();
        complete(); await settle();
        expect(current.adapter.applyState).toHaveBeenCalledTimes(2);
        expect(current.adapter.applyState).toHaveBeenLastCalledWith(next, undefined);
    });
});
