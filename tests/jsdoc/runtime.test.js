/** @jest-environment node */
'use strict';

const fs = require('node:fs/promises');
const os = require('node:os');
const path = require('node:path');
const vm = require('node:vm');
const { JSDOM } = require('jsdom');
const { prepareRuntime } = require('../../jsdoc/runtime');

let destination;
beforeEach(async () => {
    destination = await fs.mkdtemp(path.join(os.tmpdir(), '3dmol-doc-runtime-'));
    await fs.mkdir(path.join(destination, '_islands'));
    await fs.mkdir(path.join(destination, '_assets'));
    await fs.writeFile(path.join(destination, '_assets/search-index-test.json'), JSON.stringify([
        { slug: '', title: 'Home' },
        { slug: '#overview', title: 'Overview' },
        { slug: 'glviewer#addmodel', title: 'addModel' }
    ]));
    await fs.writeFile(path.join(destination, 'index.html'),
        '<script type="application/json" data-island-props>{"search":{"basePath":".","searchIndexUrl":"/doc/_assets/search-index-test.json"}}</script>');
});
afterEach(async () => {
    await fs.rm(destination, { recursive: true, force: true });
});

test('bundles island modules into a classic script with shared dependencies evaluated once', async () => {
    await fs.writeFile(path.join(destination, '_islands/chunk-shared.js'),
        'globalThis.sharedLoads++; export const value = 42;');
    for (const name of ['search', 'sidebar']) {
        await fs.writeFile(path.join(destination, `_islands/${name}-hash.js`),
            `import { value } from './chunk-shared.js'; globalThis.islands.push(['${name}', value]);`);
    }
    const { scriptPath } = await prepareRuntime(destination);
    expect(scriptPath).toMatch(/^_assets\/runtime-[a-f0-9]+\.js$/);
    const script = await fs.readFile(path.join(destination, scriptPath), 'utf8');
    const dom = new JSDOM(await fs.readFile(path.join(destination, 'index.html'), 'utf8'));
    const context = { document: dom.window.document, sharedLoads: 0, islands: [] };
    vm.runInNewContext(script, context);
    expect(context.sharedLoads).toBe(1);
    expect(context.islands).toEqual([['search', 42], ['sidebar', 42]]);
    const props = JSON.parse(context.document.querySelector('script').textContent);
    const index = JSON.parse(decodeURIComponent(props.search.searchIndexUrl.split(',', 2)[1]));
    expect(index.map(entry => entry.slug)).toEqual(['index.html', 'index.html#overview', 'glviewer/index.html#addmodel']);
    await expect(fs.access(path.join(destination, '_assets/.runtime-entry.cjs'))).rejects.toThrow();
    dom.window.close();
});

test('reports missing island entries instead of silently dropping interactive controls', async () => {
    await expect(prepareRuntime(destination)).rejects.toThrow('no island entry scripts');
});

test('normalizes dynamically generated Home links without changing external or member links', async () => {
    await fs.writeFile(path.join(destination, '_islands/mobile-nav-hash.js'), '');
    const { scriptPath } = await prepareRuntime(destination);
    const dom = new JSDOM('<body data-doc-root=".."><a href="../index.html#overview">Overview</a></body>', {
        url: 'file:///project/doc/glviewer/index.html', runScripts: 'outside-only'
    });
    dom.window.eval(await fs.readFile(path.join(destination, scriptPath), 'utf8'));
    dom.window.document.body.insertAdjacentHTML('beforeend',
        '<a id="home" href="..">Home</a><a id="external" href="https://example.com/">External</a>');
    await new Promise(resolve => setTimeout(resolve, 0));
    expect(dom.window.document.querySelector('#home').href).toBe('file:///project/doc/index.html');
    expect(dom.window.document.querySelector('#external').href).toBe('https://example.com/');
    expect(dom.window.document.querySelector('a').href).toBe('file:///project/doc/index.html#overview');
    dom.window.close();
});
