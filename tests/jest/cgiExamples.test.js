/** @jest-environment node */
'use strict';

const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const vm = require('node:vm');
const { spawnSync } = require('node:child_process');
const { JSDOM } = require('jsdom');

const root = path.resolve(__dirname, '../..');
let directory, python, manifestPath, manifest;

beforeAll(() => {
    const interpreter = spawnSync('python3', ['-c', 'import sys; print(sys.executable)'], { encoding: 'utf8' });
    if (interpreter.status !== 0) throw new Error(interpreter.stderr);
    python = interpreter.stdout.trim();
    directory = fs.mkdtempSync(path.join(os.tmpdir(), '3dmol-cgi-examples-'));
    for (const dir of ['src/nested', 'tests/auto/tests', 'unrelated-working-directory']) {
        fs.mkdirSync(path.join(directory, dir), { recursive: true });
    }
    for (const file of ['generate_tests.py', 'extract_examples.cjs', 'generate_test.cgi']) {
        fs.copyFileSync(path.join(root, 'tests/auto', file), path.join(directory, 'tests/auto', file));
    }
    fs.symlinkSync(path.join(root, 'node_modules'), path.join(directory, 'node_modules'), 'dir');
    fs.writeFileSync(path.join(directory, 'src/nested/demo.ts'), `/**
 * @example <caption>First example</caption>
 * viewer.value = 3 * 2;
 * viewer.render();
 * @example <caption>Second example</caption>
 * (async function () {
 *   await viewer.setFrame(1);
 *   viewer.value = 2 ** 3;
 *   viewer.render();
 * })();
 */
export function demo() {}
`);
    fs.writeFileSync(path.join(directory, 'tests/auto/tests/test49.js'), 'viewer.value = 49; viewer.render();');
    fs.writeFileSync(path.join(directory, 'tests/auto/tests/inline_embed.js'), `/* @data
<textarea id="inline-data">C 0 0 0</textarea>
*/
/* @div
<div class="viewer_3Dmoljs" data-element="inline-data" data-callback="window.checkFixture"></div>
*/
/* @script
window.checkFixture = function(viewer) { viewer.value = 42; viewer.render(); };
*/
window.fixtureScriptRan = true;
`);
    const generated = spawnSync(python, ['tests/auto/generate_tests.py'], { cwd: directory, encoding: 'utf8' });
    if (generated.status !== 0) throw new Error(generated.stderr || generated.stdout);
    manifestPath = path.join(directory, 'tests/auto/examples.json');
    manifest = fs.readFileSync(manifestPath, 'utf8');
}, 15000);

beforeEach(() => { fs.writeFileSync(manifestPath, manifest); });
afterAll(() => { if (directory) fs.rmSync(directory, { recursive: true, force: true }); });

function runCGI(query = '', withNode = false) {
    const result = spawnSync(python, [path.join(directory, 'tests/auto/generate_test.cgi')], {
        cwd: path.join(directory, 'unrelated-working-directory'), encoding: 'utf8',
        env: { ...process.env, PATH: withNode ? process.env.PATH : '', QUERY_STRING: query }
    });
    if (result.error) throw result.error;
    const boundary = result.stdout.search(/\r?\n\r?\n/);
    expect(boundary).toBeGreaterThan(0);
    const headers = result.stdout.slice(0, boundary);
    expect(headers).toMatch(/Content-Type: text\/html/i);
    const body = result.stdout.slice(boundary).trimStart();
    return { ...result, headers, body, document: new JSDOM(body).window.document };
}

function executeInline(document) {
    const order = [];
    const viewer = {
        render(callback) { order.push('render'); callback?.(this); },
        setFrame() { return Promise.resolve().then(() => { order.push('frame'); }); }
    };
    const context = { console, $: () => ({}), $3Dmol: { createViewer: () => viewer } };
    context.window = context;
    vm.createContext(context);
    for (const script of document.querySelectorAll('script:not([src])')) vm.runInContext(script.textContent, context);
    return { viewer, context, order };
}

test('generation writes all prepared examples and embed directives to the CGI manifest', () => {
    const data = JSON.parse(manifest);
    expect(data.version).toBe(1);
    expect(Object.keys(data.examples).sort()).toEqual([
        'inline_embed', 'nested_demo_ts_demo', 'nested_demo_ts_demo_2', 'test49'
    ]);
    expect(data.examples.nested_demo_ts_demo.script).toContain('3 * 2');
    expect(data.examples.nested_demo_ts_demo_2.script).toContain('viewer.render(callback)');
    expect(data.examples.nested_demo_ts_demo_2.script).not.toContain('<caption>');
    expect(data.examples.inline_embed.datas[0]).toContain('<textarea id="inline-data">');
    expect(data.examples.inline_embed.divs[0]).toContain('data-callback="window.checkFixture"');
    expect(data.examples.inline_embed.prescripts[0]).toContain('viewer.render(callback)');
});

test('default CGI request works without Node or source files and from another working directory', () => {
    const sources = path.join(directory, 'src');
    fs.renameSync(sources, sources + '.hidden');
    try {
        const result = runCGI();
        expect(result.status).toBe(0);
        expect(result.headers).not.toMatch(/Status: [45]/);
        expect(result.document.title).toBe('test49');
        const harness = executeInline(result.document);
        expect(harness.viewer.value).toBe(49);
        expect(harness.order).toEqual(['render']);
    } finally { fs.renameSync(sources + '.hidden', sources); }
});

test('CGI serves subsequent examples and preserves their asynchronous completion order', async () => {
    const result = runCGI('test=nested_demo_ts_demo_2');
    expect(result.status).toBe(0);
    expect(result.document.title).toBe('nested_demo_ts_demo_2');
    const harness = executeInline(result.document);
    expect(harness.order).toEqual([]);
    await new Promise(resolve => setImmediate(resolve));
    expect(harness.viewer.value).toBe(8);
    expect(harness.order).toEqual(['frame', 'render']);
});

test('CGI preserves inline data, embedded containers, prescripts, and the main script', () => {
    const result = runCGI('test=inline_embed');
    expect(result.status).toBe(0);
    expect(result.document.querySelector('#inline-data').textContent).toBe('C 0 0 0');
    expect(result.document.querySelector('.viewer_3Dmoljs').getAttribute('data-callback')).toBe('window.checkFixture');
    const harness = executeInline(result.document);
    expect(harness.context.fixtureScriptRan).toBe(true);
    expect(typeof harness.context.checkFixture).toBe('function');
    harness.context.checkFixture(harness.viewer);
    expect(harness.viewer.value).toBe(42);
    expect(harness.order).toEqual(['render']);
});

test('CGI escapes unknown names and returns a visible 404 response', () => {
    const name = '<script>alert("unknown")</script>';
    const result = runCGI('test=' + encodeURIComponent(name));
    expect(result.headers).toMatch(/Status: 404\b/);
    expect(result.document.body.textContent).toContain(name);
    expect(result.document.querySelectorAll('script')).toHaveLength(0);
    expect(result.body).toContain('&lt;script&gt;');
});

test('missing manifest and unavailable Node return actionable HTML instead of a blank response', () => {
    fs.unlinkSync(manifestPath);
    const result = runCGI();
    expect(result.headers).toMatch(/Status: 500\b/);
    expect(result.document.body.textContent).toContain('npm run generate:tests');
    expect(result.document.body.textContent.trim().length).toBeGreaterThan(30);
});

test('CGI can still extract live source when the manifest is absent and Node is available', () => {
    fs.unlinkSync(manifestPath);
    const result = runCGI('test=nested_demo_ts_demo', true);
    expect(result.status).toBe(0);
    expect(result.document.title).toBe('nested_demo_ts_demo');
    const harness = executeInline(result.document);
    expect(harness.viewer.value).toBe(6);
    expect(harness.order).toEqual(['render']);
});
