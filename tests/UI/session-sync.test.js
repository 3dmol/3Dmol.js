/** @jest-environment node */
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { JSDOM } = require('jsdom');
const jquery = require('jquery');
const source = fs.readFileSync(path.join(__dirname, '../../learning-environment/sessions.js'), 'utf8');

function fixture(url = 'https://classroom.test/viewer.html', server) {
    const dom = new JSDOM('<input id="session_name_input" value="Classroom"><button id="sessionconnectbutton">Create</button><button id="sessionleave"></button><button id="sessiondestroy"></button>', { url });
    const window = dom.window, $ = jquery(window), handlers = new Map();
    if (server !== undefined) window.viewerSessionServer = server;
    const socket = { connected: true, on: jest.fn((name, callback) => handlers.set(name, callback)), emit: jest.fn(), send: jest.fn() };
    const connect = jest.fn(() => socket);
    let data = { models: [{ atoms: Array.from({ length: 50 }, (_, index) => ({ index, x: index, y: 0, z: 0 })) }], modelStyles: [], shapes: [], surfaces: [], labels: [], frame: 0, viewSettings: { background: '#ffffff' } };
    let uiChange, sceneChange, viewChange;
    const unsubscribeUI = jest.fn(), unsubscribeScene = jest.fn();
    const viewer = {
        getView: () => [0, 0, 0, -30, 0, 0, 0, 1], getFrame: () => data.frame,
        getInternalState: () => data, setView: jest.fn(), setInternalState: jest.fn(),
        setViewChangeCallback: callback => { viewChange = callback; }, setStateChangeCallback: jest.fn(),
        subscribeStateChanges: callback => { sceneChange = callback; return unsubscribeScene; },
        setHoverable: jest.fn(), setClickable: jest.fn(), removeLabels: jest.fn()
    };
    const app = window.viewerApp = {
        getSessionState: jest.fn(() => data),
        subscribeSessionChanges: callback => { uiChange = callback; return unsubscribeUI; },
        acceptSessionState: jest.fn(() => Promise.resolve()), showError: jest.fn(), sessionChanged: jest.fn(),
        state: { ui: { element: window.document.createElement('div') } }
    };
    const context = vm.createContext({ window, io: { connect }, $, glviewer: viewer,
        console: { log() {} }, setTimeout, clearTimeout, alert: jest.fn() });
    vm.runInContext(source, context);
    context.initSessions();
    const receive = (event, payload) => handlers.get(event)(payload);
    return {
        window, $, socket, connect, app, viewer, receive, unsubscribeUI, unsubscribeScene,
        initialize: () => context.initSessions(), joinSession: name => context.joinSession(name),
        get data() { return data; }, set data(value) { data = value; },
        create() { $('#sessionconnectbutton').trigger('click'); receive('create session response', 1234); },
        join() { $('#sessionconnectbutton').text('Join').trigger('click'); receive('join session response', 'student-1'); },
        uiChange: () => uiChange(), sceneChange: () => sceneChange(), viewChange: view => viewChange(view),
        updates: () => socket.emit.mock.calls.filter(([event]) => event === 'viewer state change event').map(([, packet]) => packet),
        close: () => { window.leaveViewerSession(); dom.window.close(); }
    };
}

describe('classroom state transport', () => {
    let current;
    beforeEach(() => { jest.useFakeTimers(); current = fixture(); });
    afterEach(() => { current.close(); jest.useRealTimers(); });

    test.each([
        { url: 'http://localhost/viewer.html', expected: 'http://localhost' },
        { url: 'http://127.0.0.1/viewer.html', expected: 'http://127.0.0.1' },
        { url: 'http://localhost:5068/static/viewer.html', expected: 'http://localhost:5068' },
        { url: 'http://localhost/viewer.html', server: 'http://127.0.0.1:5050', expected: 'http://127.0.0.1:5050' }
    ])('$url connects to $expected', ({ url, server, expected }) => {
        current.close(); current = fixture(url, server);
        expect(current.connect).toHaveBeenCalledWith(expected);
    });

    test('repeated initialization retains one socket, one set of handlers, and the active role', () => {
        const handlerCount = current.socket.on.mock.calls.length;
        expect(current.initialize()).toBe(current.socket);
        current.create();
        expect(current.initialize()).toBe(current.socket);
        expect(current.connect).toHaveBeenCalledTimes(1);
        expect(current.socket.on).toHaveBeenCalledTimes(handlerCount);
        expect(current.socket.emit.mock.calls.filter(([event]) => event === 'create session event')).toHaveLength(1);
        current.data.frame = 1; current.uiChange(); jest.advanceTimersByTime(100);
        expect(current.updates()).toEqual([{ name: 'Classroom', baseRevision: 0, revision: 1, patch: { frame: 1 } }]);
        current.window.leaveViewerSession();
        expect(current.unsubscribeUI).toHaveBeenCalledTimes(1);
        expect(current.unsubscribeScene).toHaveBeenCalledTimes(1);
    });

    test('windows on the same origin retain independent roles through repeated connection actions', () => {
        const student = fixture();
        try {
            current.create(); student.join();
            expect(student.socket).not.toBe(current.socket);
            for (const client of [current, student]) {
                client.socket.emit.mockClear();
                client.joinSession('Other classroom');
                client.receive('join session response', 'unsolicited-student');
                client.$('#session_name_input').val('Other classroom');
                client.$('#sessionconnectbutton').text('Join').trigger('click');
                client.receive('join session response', 'unsolicited-student');
                client.$('#sessionconnectbutton').text('Create');
                client.$('#session_name_input').trigger(client.$.Event('keyup', { keyCode: 13 }));
                client.receive('create session response', 5678);
                expect(client.socket.emit).not.toHaveBeenCalled();
            }
            expect(current.app.sessionChanged).not.toHaveBeenCalledWith(expect.anything(), true);
            expect(current.app.state.ui.element.hidden).toBe(false);
            expect(student.app.sessionChanged).toHaveBeenCalledTimes(1);
            expect(student.app.sessionChanged).toHaveBeenCalledWith('Classroom', true);
            expect(student.app.state.ui.element.hidden).toBe(true);
            current.data.frame = 2; current.uiChange(); jest.advanceTimersByTime(100);
            const update = current.updates()[0];
            expect(update).toEqual({ name: 'Classroom', baseRevision: 0, revision: 1, patch: { frame: 2 } });
            student.receive('viewer state change response', { state: student.data, revision: 0 });
            student.receive('viewer state change response', update);
            expect(student.app.acceptSessionState).toHaveBeenLastCalledWith({ ...student.data, frame: 2 }, ['frame']);
            expect(student.updates()).toHaveLength(0);
        } finally { student.close(); }
    });

    test('coalesces edits, sends only changed sections, and suppresses unchanged renders', () => {
        current.create(); jest.advanceTimersByTime(100);
        expect(current.updates()).toHaveLength(0);
        current.app.getSessionState.mockClear();
        current.data.viewSettings = { background: '#000000' };
        for (let i = 0; i < 20; i++) { current.uiChange(); current.sceneChange(); }
        expect(current.app.getSessionState).not.toHaveBeenCalled();
        jest.advanceTimersByTime(100);
        expect(current.app.getSessionState).toHaveBeenCalledTimes(1);
        expect(current.updates()).toEqual([{ name: 'Classroom', baseRevision: 0, revision: 1, patch: { viewSettings: { background: '#000000' } } }]);
        current.receive('viewer state revision', { revision: 1 });
        for (let i = 0; i < 10; i++) current.sceneChange();
        jest.advanceTimersByTime(100);
        expect(current.updates()).toHaveLength(1);
    });

    test('serializes revisions and flushes the latest edit after acknowledgement', () => {
        current.create(); jest.advanceTimersByTime(100);
        current.data.frame = 1; current.uiChange(); jest.advanceTimersByTime(100);
        current.data.frame = 2; current.uiChange(); jest.advanceTimersByTime(100);
        current.data.frame = 3; current.uiChange(); jest.advanceTimersByTime(100);
        expect(current.updates()).toHaveLength(1);
        current.receive('viewer state revision', { revision: 1 }); jest.advanceTimersByTime(100);
        expect(current.updates()[1]).toEqual({ name: 'Classroom', baseRevision: 1, revision: 2, patch: { frame: 3 } });
    });

    test('camera changes do not serialize or send molecular state', () => {
        current.create(); jest.advanceTimersByTime(100); current.app.getSessionState.mockClear();
        const view = [1, 2, 3, -40, 0, 0, 0, 1];
        for (let i = 0; i < 20; i++) current.viewChange(view);
        jest.advanceTimersByTime(33);
        expect(current.socket.emit.mock.calls.filter(([event]) => event === 'viewer view change event')).toEqual([
            ['viewer view change event', { name: 'Classroom', view }]
        ]);
        current.viewChange(view); jest.advanceTimersByTime(100);
        expect(current.app.getSessionState).not.toHaveBeenCalled();
        expect(current.updates()).toHaveLength(0);
    });

    test('students merge patches and request a full snapshot after a revision gap', () => {
        current.join();
        expect(current.app.sessionChanged).toHaveBeenCalledWith('Classroom', true);
        expect(current.app.state.ui.element.hidden).toBe(true);
        current.receive('viewer state change response', { state: current.data, revision: 4 });
        current.receive('viewer state change response', { patch: { frame: 2 }, baseRevision: 4, revision: 5 });
        expect(current.app.acceptSessionState.mock.calls[1]).toEqual([{ ...current.data, frame: 2 }, ['frame']]);
        current.receive('viewer state change response', { patch: { frame: 3 }, baseRevision: 6, revision: 7 });
        expect(current.app.acceptSessionState).toHaveBeenCalledTimes(2);
        expect(current.socket.emit).toHaveBeenCalledWith('request state', { name: 'Classroom' });
        current.receive('leave session response');
        expect(current.app.sessionChanged).toHaveBeenLastCalledWith(null, false);
        expect(current.app.state.ui.element.hidden).toBe(false);
    });

    test('reconnect sends one fresh full state at the server revision', () => {
        current.create(); jest.advanceTimersByTime(100);
        current.socket.connected = false; current.receive('disconnect', 'transport close');
        current.data.frame = 3; current.sceneChange(); jest.advanceTimersByTime(100);
        expect(current.updates()).toHaveLength(0);
        current.socket.connected = true; current.receive('connect'); current.receive('reconnect');
        expect(current.socket.emit.mock.calls.filter(([event]) => event === 'reconnection')).toHaveLength(1);
        current.receive('reconnection response', { controller: true, sid: 'new-instructor', revision: 8 });
        jest.advanceTimersByTime(100);
        expect(current.updates()).toEqual([{ name: 'Classroom', baseRevision: 8, revision: 9, state: current.data }]);
    });

    test('leaving detaches state subscriptions and cancels queued sends', () => {
        current.create(); current.uiChange(); current.window.leaveViewerSession(); jest.advanceTimersByTime(100);
        expect(current.unsubscribeUI).toHaveBeenCalledTimes(1);
        expect(current.unsubscribeScene).toHaveBeenCalledTimes(1);
        expect(current.updates()).toHaveLength(0);
    });
});
