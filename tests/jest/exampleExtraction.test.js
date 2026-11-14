/** @jest-environment node */
'use strict';

const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const vm = require('node:vm');
const { spawnSync } = require('node:child_process');
const { extractExamples, prepareExample } = require('../auto/extract_examples.cjs');

const root = path.resolve(__dirname, '../..');
const declaration = comment => `${comment}\nexport function demo() {}`;
const example = code => declaration(`/**\n * @example\n${code}\n */`);

function runPython(directory, script) {
    const result = spawnSync('python3', [script], { cwd: directory, encoding: 'utf8' });
    if (result.status !== 0) throw new Error(result.stderr || result.stdout);
    return result;
}

function viewerContext() {
    const order = [];
    let completed = 0;
    const viewer = {
        render(callback) { order.push('render'); callback?.(this); },
        setFrame() { return Promise.resolve().then(() => { order.push('frame'); }); },
        surfacesFinished: () => true,
        isAnimated: () => false
    };
    const context = {
        viewer, callback: () => { completed++; }, console, setTimeout,
        $3Dmol: { createViewer: () => viewer, setSyncSurface() {}, viewers: { 0: viewer } },
        $: () => ({}), document: { body: {} }
    };
    context.window = context;
    context.global = context;
    return { context: vm.createContext(context), viewer, order, completions: () => completed };
}

const fixture = `/**
 * @example <caption>First view</caption>
 * viewer.value = 3 * 2;
 * viewer.literal = "viewer.render() @scope/pkg ** text";
 * viewer.render();
 * @param options - Trailing documentation is not executable.
 * @example <caption>Frame change</caption>
 * (async function () {
 *   await viewer.setFrame(1);
 *   viewer.value = 2 ** 3;
 *   viewer.render();
 * })();
 */
export function demo() {}
`;

let directory;
beforeEach(() => { directory = fs.mkdtempSync(path.join(os.tmpdir(), '3dmol-example-extraction-')); });
afterEach(() => { fs.rmSync(directory, { recursive: true, force: true }); });

test('extracts every example and excludes captions and intervening tags', () => {
    const examples = extractExamples(fixture, 'src/demo.ts');
    expect(examples.map(item => [item.name, item.ordinal])).toEqual([['demo', 1], ['demo', 2]]);
    expect(examples[0].code).toContain('3 * 2');
    expect(examples[1].code).toContain('2 ** 3');
    for (const item of examples) expect(item.code).not.toMatch(/<caption>|@param/);
});

test('handles legacy doubled tag stars, same-line code, typedefs, classes and nested declarations', () => {
    const source = `/** @typedef {Object} Shape
 * * @example viewer.render();
 */
export interface Shape {}
/** @example viewer.render(); */
export class Demo {
 /** @example viewer.render();
  * @example viewer.render();
  */
 public draw() {}
}
/** @example viewer.render(); */
export const colors = {};
`;
    expect(extractExamples(source, 'src/nested/demo.ts').map(item => [item.name, item.ordinal]))
        .toEqual([['Shape', 1], ['Demo', 1], ['draw', 1], ['draw', 2], ['colors', 1]]);
});

test('ignores fake JSDoc declarations in strings and template literals', () => {
    const source = 'const text = "/** @example viewer.render(); */ export function fake() {}";\n'
        + 'const other = `/** @example viewer.render(); */ export class AlsoFake {}`;\n' + fixture;
    expect(extractExamples(source, 'demo.ts')).toHaveLength(2);
});

test('keeps @ tags and operators inside strings, regexes and multiline templates', () => {
    const source = example(' * const email = "a@example.org";\n'
        + ' * const pattern = /a*b/;\n'
        + ' * const text = `first\n * @example still text ${2 * 3}\n * @returns also text\n * last`;\n'
        + ' * viewer.render();');
    const examples = extractExamples(source, 'demo.ts');
    expect(examples).toHaveLength(1);
    const code = examples[0].code;
    expect(code).toContain('a@example.org');
    expect(code).toContain('/a*b/');
    expect(code).toContain('@example still text ${2 * 3}');
    expect(code).toContain('@returns also text');
});

test.each(['*', '**'])('preserves unstarred continuation operator %s', operator => {
    const code = extractExamples(example(` const value = 3\n ${operator} 2;\n viewer.value = value;\n viewer.render();`), 'demo.ts')[0].code;
    const harness = viewerContext();
    vm.runInContext(prepareExample(code).script, harness.context);
    expect(harness.viewer.value).toBe(operator === '*' ? 6 : 9);
    expect(harness.completions()).toBe(1);
});

test('does not split multiline captions or swallow caption-looking code strings', () => {
    const caption = declaration('/**\n * @example <caption>First\n * @example still caption\n * </caption>\n * viewer.render();\n */');
    expect(extractExamples(caption, 'demo.ts')[0].code.trim()).toBe('viewer.render();');
    const strings = declaration('/**\n * @example const text="<caption>"; viewer.render();\n * @example const text="</caption>"; viewer.render();\n */');
    expect(extractExamples(strings, 'demo.ts')).toHaveLength(2);
});

test('keeps literal stars after template interpolation in unstarred examples', () => {
    const source = example('const text = `start ${1}\n *literal\nend`;\nviewer.value = text;\nviewer.render();');
    const code = extractExamples(source, 'demo.ts')[0].code;
    const harness = viewerContext();
    vm.runInContext(prepareExample(code).script, harness.context);
    expect(harness.viewer.value).toBe('start 1\n *literal\nend');
});

test('ignores fake caption openers inside templates before a later real example', () => {
    const source = declaration('/**\n * @example\n * const text = `start\n'
        + ' * @example <caption>fake text\n * end`;\n * viewer.render();\n'
        + ' * @example <caption>Second real caption</caption>\n * viewer.render();\n */');
    const examples = extractExamples(source, 'demo.ts');
    expect(examples).toHaveLength(2);
    expect(examples[0].code).toContain('@example <caption>fake text');
    expect(examples[1].code.trim()).toBe('viewer.render();');
});

test.each([
    ['', /empty @example/],
    ['const broken = ;', /invalid JavaScript/],
    ['const value: number = 2; viewer.render();', /invalid browser JavaScript/],
    ['await viewer.setFrame(1); viewer.render();', /wrap await\/return/]
])('rejects unusable examples rather than emitting empty passing tests: %s', (code, message) => {
    expect(() => extractExamples(example(code), 'demo.ts')).toThrow(message);
});

test('recognizes actual harness directives and preserves comment-looking data and token boundaries', () => {
    const code = '/* @data <textarea id="data">hello</textarea> */\n'
        + '/* @div <div class="viewer_3Dmoljs"></div> */\n'
        + '/* @script window.draw = function(viewer) { const text="viewer.render()"; viewer.render(); }; */\n'
        + 'const/* separator */value = 6;\nconst text="/* @div NOT HTML */";\nviewer.render(/*no callback*/);';
    const result = prepareExample(code);
    expect(result.datas).toHaveLength(1);
    expect(result.divs).toEqual(['<div class="viewer_3Dmoljs"></div> ']);
    expect(result.prescripts[0]).toContain('const text="viewer.render()"; viewer.render(callback);');
    expect(result.script).toContain('const/* separator */value = 6');
    expect(result.script).toContain('"/* @div NOT HTML */"');
    expect(result.script).toContain('viewer.render(/*no callback*/)');
});

test('only instruments real final renders and preserves intentional intermediate draws', async () => {
    const code = 'viewer.literal = "viewer.render()"; viewer.render( );\n'
        + '(async function () { await viewer.setFrame(1); viewer.render(); })();';
    const prepared = prepareExample(code);
    const harness = viewerContext();
    vm.runInContext(prepared.script, harness.context);
    expect(harness.order).toEqual(['render']);
    expect(harness.completions()).toBe(0);
    await new Promise(resolve => setImmediate(resolve));
    expect(harness.order).toEqual(['render', 'frame', 'render']);
    expect(harness.completions()).toBe(1);
    expect(harness.viewer.literal).toBe('viewer.render()');
});

test('a rejected awaited operation cannot signal a successful render', async () => {
    const code = '(async function () { await viewer.setFrame(1); viewer.render(); })();';
    const harness = viewerContext();
    harness.viewer.setFrame = () => Promise.reject(new Error('frame failed'));
    await expect(vm.runInContext(prepareExample(code).script, harness.context)).rejects.toThrow('frame failed');
    expect(harness.order).toEqual([]);
    expect(harness.completions()).toBe(0);
});

function prepareProject() {
    for (const dir of ['src/nested', 'tests/auto/tests', 'tests/jest', 'tests/glcheck', 'build']) {
        fs.mkdirSync(path.join(directory, dir), { recursive: true });
    }
    for (const file of ['tests/auto/generate_tests.py', 'tests/auto/extract_examples.cjs',
        'tests/jest/generate_jest_render_tests.py', 'tests/glcheck/generate_glcheck_render_tests.py']) {
        fs.copyFileSync(path.join(root, file), path.join(directory, file));
    }
    fs.symlinkSync(path.join(root, 'node_modules'), path.join(directory, 'node_modules'), 'dir');
    fs.writeFileSync(path.join(directory, 'src/nested/demo.ts'), fixture);
}

test('all three generators preserve legacy names, emit every example, and await async completion', async () => {
    prepareProject();
    for (const script of ['tests/auto/generate_tests.py', 'tests/jest/generate_jest_render_tests.py',
        'tests/glcheck/generate_glcheck_render_tests.py']) runPython(directory, script);
    const names = ['nested_demo_ts_demo', 'nested_demo_ts_demo_2'];
    const browser = fs.readFileSync(path.join(directory, 'tests/auto/build.js'), 'utf8');
    const browserHarness = viewerContext();
    vm.runInContext(browser, browserHarness.context);
    expect(Object.keys(browserHarness.context.system)).toEqual(names);
    browserHarness.context.system[names[1]](browserHarness.viewer, browserHarness.context.callback);
    expect(browserHarness.completions()).toBe(0);
    await new Promise(resolve => setImmediate(resolve));
    expect(browserHarness.order).toEqual(['frame', 'render']);
    expect(browserHarness.completions()).toBe(1);

    const jestSource = fs.readFileSync(path.join(directory, 'tests/jest/render.test.js'), 'utf8');
    const jestHarness = viewerContext();
    const cases = new Map();
    Object.assign(jestHarness.context, {
        beforeAll() {}, afterAll() {}, test: (name, fn) => cases.set(name, fn),
        jest: { setTimeout() {} }, WebGLRenderingContext: function () {},
        require: name => name === 'util' ? require('node:util')
            : name === 'jquery' ? jestHarness.context.$ : jestHarness.context.$3Dmol
    });
    vm.runInContext(jestSource, jestHarness.context);
    expect([...cases.keys()]).toEqual(names);
    let finished = 0;
    cases.get(names[1])(() => { finished++; });
    expect(finished).toBe(0);
    await new Promise(resolve => setImmediate(resolve));
    expect(jestHarness.order).toEqual(['frame', 'render']);
    expect(finished).toBe(1);

    for (const [index, name] of names.entries()) {
        const html = fs.readFileSync(path.join(directory, `tests/glcheck/render-tests/${name}.html`), 'utf8');
        expect(html).toContain('src="./jquery.min.js"');
        expect(html).not.toContain('https://code.jquery.com');
        const harness = viewerContext();
        for (const script of html.matchAll(/<script>([\s\S]*?)<\/script>/g)) vm.runInContext(script[1], harness.context);
        if (index) expect(harness.context.glcheck_renderDone).toBeUndefined();
        await new Promise(resolve => setImmediate(resolve));
        expect(harness.context.glcheck_renderDone).toBe(true);
        expect(harness.viewer.value).toBe(index ? 8 : 6);
        if (!index) expect(harness.viewer.literal).toBe('viewer.render() @scope/pkg ** text');
    }
    expect(fs.existsSync(path.join(directory, 'tests/glcheck/render-tests/jquery.min.js'))).toBe(true);
}, 15000);

test('duplicate names fail before generators can overwrite a fixture', () => {
    prepareProject();
    fs.writeFileSync(path.join(directory, 'tests/auto/tests/nested_demo_ts_demo_2.js'), 'viewer.render();');
    const result = spawnSync('python3', ['tests/auto/generate_tests.py'], { cwd: directory, encoding: 'utf8' });
    expect(result.status).not.toBe(0);
    expect(result.stderr).toContain("Duplicate render example name 'nested_demo_ts_demo_2'");
});

test('the interactive harness reports a rejected example and advances only once', () => {
    const source = fs.readFileSync(path.join(root, 'tests/auto/test.js'), 'utf8');
    const scheduled = [], ran = [];
    const chain = new Proxy({}, { get: (_target, name) => name === 'ready' ? callback => callback() : () => chain });
    const context = {
        document: {}, console: { log() {} }, $: () => chain,
        system: { first() { ran.push('first'); }, second() { ran.push('second'); } },
        $3Dmol: { createViewer: () => ({}) },
        resemble: { outputSettings() {} }, JSZip: function () {},
        setTimeout: callback => scheduled.push(callback)
    };
    context.window = context;
    vm.runInNewContext(source, context);
    expect(ran).toEqual(['first']);
    const event = { reason: new Error('async failure'), preventDefault: jest.fn() };
    context.onunhandledrejection(event);
    context.onerror('duplicate report');
    expect(event.preventDefault).toHaveBeenCalledTimes(1);
    expect(scheduled).toHaveLength(1);
    scheduled[0]();
    expect(ran).toEqual(['first', 'second']);
});
