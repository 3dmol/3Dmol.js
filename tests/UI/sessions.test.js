/** @jest-environment node */
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const source = fs.readFileSync(path.join(__dirname, '../../learning-environment/sessions.js'), 'utf8');

function initialize(pageURL, server) {
    const socket = { on: jest.fn(), send: jest.fn(), emit: jest.fn() };
    const io = { connect: jest.fn(() => socket) };
    const element = { on: jest.fn(), html: jest.fn() };
    const window = { location: new URL(pageURL) };
    const context = vm.createContext({ window, io, $: () => element });
    vm.runInContext(source, context, { filename: 'sessions.js' });
    // Configuration may be supplied after loading the script, before initialization.
    if (server !== undefined) window.viewerSessionServer = server;
    context.initSessions();
    return io.connect;
}

describe('classroom session server routing', () => {
    test.each([
        ['http://localhost/viewer.html', 'http://localhost:5000'],
        ['http://localhost:8080/viewer.html', 'http://localhost:5000'],
        ['https://localhost:8443/viewer.html', 'https://localhost:5000'],
        ['https://classroom.example/viewer.html', 'https://classroom.example'],
        ['https://classroom.example:8443/static/viewer.html', 'https://classroom.example:8443'],
        ['http://classroom.example:8080/viewer.html', 'http://classroom.example:8080'],
        ['http://127.0.0.1:8080/viewer.html', 'http://127.0.0.1:8080']
    ])('connects from %s to the compatible default %s', (pageURL, expected) => {
        const connect = initialize(pageURL);
        expect(connect).toHaveBeenCalledTimes(1);
        expect(connect).toHaveBeenCalledWith(expected);
    });

    test.each([
        ['http://localhost:8080/viewer.html', 'http://localhost:6000'],
        ['http://localhost:8080/viewer.html', 'http://localhost:8080'],
        ['https://classroom.example/viewer.html', 'https://sessions.example:8443']
    ])('allows an explicit session server on %s', (pageURL, server) => {
        const connect = initialize(pageURL, server);
        expect(connect).toHaveBeenCalledTimes(1);
        expect(connect).toHaveBeenCalledWith(server);
    });
});
