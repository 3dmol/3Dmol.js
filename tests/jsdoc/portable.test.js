/** @jest-environment node */
'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');
const os = require('node:os');
const vm = require('node:vm');
const { JSDOM } = require('jsdom');
const { makePortable } = require('../../jsdoc/portable');

let destination;
async function write(file, text = '') {
    await fs.mkdir(path.dirname(path.join(destination, file)), { recursive: true });
    await fs.writeFile(path.join(destination, file), text);
}
async function document(file) {
    return new JSDOM(await fs.readFile(path.join(destination, file), 'utf8')).window.document;
}
beforeEach(async () => { destination = await fs.mkdtemp(path.join(os.tmpdir(), '3dmol-portable-')); });
afterEach(async () => { await fs.rm(destination, { recursive: true, force: true }); });

test('makes styles, nested navigation and sibling assets relative with explicit HTML filenames', async () => {
    const content = '<link rel="stylesheet" href="/doc/_assets/style.css">'
        + '<a href="/doc">Home</a><a href="/doc/glviewer#addmodel">Viewer</a>'
        + '<a href="https://example.com/doc/">External</a><a href="#local">Local</a>'
        + '<script src="/build/3Dmol.js"></script><iframe src="/tests/example.html"></iframe>'
        + '<button data-pdb-url="/tutorials/doc-data/1ycr.pdb">Load</button>'
        + '<pre>&lt;a href="/doc/glviewer"&gt;Example&lt;/a&gt;</pre>';
    await write('index.html', content);
    await write('glviewer/index.html', content);
    await write('_assets/style.css');
    await makePortable(destination, { scriptPath: '_assets/runtime.js' });
    const home = await document('index.html');
    const viewer = await document('glviewer/index.html');
    expect(home.querySelector('link').getAttribute('href')).toBe('_assets/style.css');
    expect([...home.querySelectorAll('a')].map(node => node.getAttribute('href')))
        .toEqual(['index.html', 'glviewer/index.html#addmodel', 'https://example.com/doc/', '#local']);
    expect(viewer.querySelector('link').getAttribute('href')).toBe('../_assets/style.css');
    expect(viewer.querySelector('a').getAttribute('href')).toBe('../index.html');
    expect(viewer.querySelector('script').getAttribute('src')).toBe('../../build/3Dmol.js');
    expect(viewer.querySelector('iframe').getAttribute('src')).toBe('../../tests/example.html');
    expect(viewer.querySelector('button').getAttribute('data-pdb-url')).toBe('../../tutorials/doc-data/1ycr.pdb');
    expect(viewer.querySelector('pre').textContent).toBe('<a href="/doc/glviewer">Example</a>');
});

test('uses a classic runtime and gives hydrated navigation the same portable URLs', async () => {
    const props = { i0: { basePath: '/doc', currentSlug: 'glviewer', nav: [{ slug: 'glviewer' }, { slug: '' }, { href: '/doc/tutorials/home/' }], mdUrl: '/doc/glviewer/index.md' } };
    await write('index.html');
    await write('tutorials/home/index.html');
    await write('glviewer/index.html', `<script type="application/json" data-island-props>${JSON.stringify(props)}</script>`
        + `<script type="module">document.querySelectorAll('[data-island]');</script>`);
    await makePortable(destination, { scriptPath: '_assets/runtime.js' });
    const doc = await document('glviewer/index.html');
    const actual = JSON.parse(doc.querySelector('[data-island-props]').textContent).i0;
    expect(actual).toEqual({ basePath: '..', currentSlug: 'glviewer/index.html', nav: [{ slug: 'glviewer/index.html' }, { slug: 'index.html' }, { href: '../tutorials/home/index.html' }], mdUrl: 'index.md' });
    expect(doc.querySelector('script[src]').getAttribute('src')).toBe('../_assets/runtime.js');
    expect(doc.querySelector('script[src]').hasAttribute('defer')).toBe(true);
    expect(doc.querySelector('script[type="module"]')).toBeNull();
});

test('preserves legacy member redirects when opened as files', async () => {
    await write('glviewer/index.html');
    await write('GLViewer.html', '<title>Documentation moved</title>'
        + '<meta http-equiv="refresh" content="0;url=/doc/glviewer/">'
        + '<script>location.replace("/doc/glviewer/"+location.search+location.hash)</script>');
    await makePortable(destination, { scriptPath: '_assets/runtime.js' });
    const doc = await document('GLViewer.html');
    expect(doc.querySelector('meta').content).toBe('0;url=glviewer/index.html');
    const replace = jest.fn();
    vm.runInNewContext(doc.querySelector('script').textContent, { location: { search: '?test=1', hash: '#addmodel', replace } });
    expect(replace).toHaveBeenCalledWith('glviewer/index.html?test=1#addmodel');
});
