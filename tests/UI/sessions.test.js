/** @jest-environment node */
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(path.join(__dirname, '../../learning-environment/sessions.js'), 'utf8');

function fixture() {
    const session = { join: jest.fn(), leave: jest.fn() };
    const app = { enableSessions: jest.fn(() => Promise.resolve(session)) };
    const window = { viewerApp: app };
    const context = vm.createContext({ window, Promise, Error });
    vm.runInContext(source, context, { filename: 'sessions.js' });
    return { session, app, window, context };
}

describe('legacy classroom entry points', () => {
    test('initialization reuses the UI-owned client and join/leave delegate to it', async () => {
        const { session, app, window, context } = fixture();
        const first = context.initSessions();
        expect(context.initSessions()).toBe(first);
        await expect(first).resolves.toBe(session);
        await context.joinSession('Classroom');
        expect(app.enableSessions).toHaveBeenCalledTimes(1);
        expect(session.join).toHaveBeenCalledWith('Classroom');
        await window.leaveViewerSession();
        expect(session.leave).toHaveBeenCalledTimes(1);
    });

    test('leaving before initialization does not open a connection', () => {
        const { app, window } = fixture();
        window.leaveViewerSession();
        expect(app.enableSessions).not.toHaveBeenCalled();
    });

    test('leaving while initialization is pending cancels a queued join', async () => {
        const { session, app, window, context } = fixture();
        let resolve;
        app.enableSessions.mockImplementation(() => new Promise(done => { resolve = done; }));
        const joining = context.joinSession('Classroom');
        await Promise.resolve();
        const leaving = window.leaveViewerSession();
        resolve(session);
        await Promise.all([joining, leaving]);
        expect(session.join).not.toHaveBeenCalled();
        expect(session.leave).toHaveBeenCalledTimes(1);
    });

    test('initialization failures can be retried and a replaced viewer gets its own client', async () => {
        const { session, app, window, context } = fixture();
        app.enableSessions.mockRejectedValueOnce(new Error('Offline'));
        await expect(context.initSessions()).rejects.toThrow('Offline');
        await expect(context.initSessions()).resolves.toBe(session);
        const replacement = { join: jest.fn(), leave: jest.fn() };
        window.viewerApp = { enableSessions: jest.fn(() => replacement) };
        await expect(context.initSessions()).resolves.toBe(replacement);
        await context.joinSession('Other classroom');
        expect(replacement.join).toHaveBeenCalledWith('Other classroom');
        expect(session.join).not.toHaveBeenCalled();
    });

    test('calling before the viewer exists reports the missing initialization', async () => {
        const { window, context } = fixture();
        delete window.viewerApp;
        await expect(context.initSessions()).rejects.toThrow('Initialize the standalone viewer');
    });
});
