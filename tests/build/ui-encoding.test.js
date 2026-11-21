/** @jest-environment node */
jest.mock('canvas', () => ({}));
const fs = require('node:fs');
const path = require('node:path');
const { TextDecoder } = require('node:util');
const { JSDOM } = require('jsdom');

// A classic external script can inherit Windows-1252 from a local HTML file
// without a charset declaration. Exercise the shipped bytes with that decoding.
test.each(['utf-8', 'windows-1252'])('classroom symbols survive %s script decoding', async encoding => {
    const dom = new JSDOM('<!doctype html><div id="viewer"></div>', {
        url: 'file:///classroom.html', runScripts: 'outside-only'
    });
    let session;
    try {
        const bytes = fs.readFileSync(path.join(__dirname, '../../build/3Dmol.ui-min.js'));
        dom.window.eval(new TextDecoder(encoding).decode(bytes));
        const { LearningSession } = dom.window['3Dmol.ui'];
        const socket = { on() {}, off() {}, emit() {}, disconnect() {} };
        session = new LearningSession({ viewer: { container: dom.window.document.getElementById('viewer') } }, {
            mode: 'hidden', session: 'encoding-test', socketFactory: () => socket
        });
        await session.ready;

        expect(session.element.querySelector('[aria-label="Close or leave classroom"]').textContent).toBe('\u00d7');
        expect(session.element.querySelector('[aria-label="Refresh responses"]').textContent).toBe('\u21bb');
        const rules = session.element.querySelector('style').sheet.cssRules;
        const handle = Array.from(rules).find(rule => rule.selectorText === '.molui-learning-drag::before');
        expect(handle.style.getPropertyValue('content')).toMatch(/^['"]\u283f['"]$/);
    } finally {
        session?.destroy();
        dom.window.close();
    }
});
