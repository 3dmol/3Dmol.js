/** @jest-environment node */
'use strict';

const fs = require('node:fs/promises');
const os = require('node:os');
const path = require('node:path');
const vm = require('node:vm');
const { finish } = require('../../jsdoc/postprocess');

let destination;
const legacyLinks = new Map([['GLViewer', 'GLViewer.html'], ['Gradient.RWB', 'Gradient.RWB.html']]);

async function write(file, text) {
    await fs.mkdir(path.dirname(path.join(destination, file)), { recursive: true });
    await fs.writeFile(path.join(destination, file), text);
}

const read = file => fs.readFile(path.join(destination, file), 'utf8');

beforeEach(async () => {
    destination = await fs.mkdtemp(path.join(os.tmpdir(), '3dmol-jsdoc-postprocess-'));
    await write('glviewer/index.md', '---\nkind: class\nlongname: GLViewer\n---\n# GLViewer\n');
    await write('glviewer/index.html', '<h1>GLViewer</h1>');
    await write('global/index.md', '---\nkind: global\nlongname: Globals\n---\n');
    await write('global/index.html', '<h1>Globals</h1>');
});

afterEach(async () => {
    await fs.rm(destination, { recursive: true, force: true });
});

test('keeps source filename links and source pages, without inaccurate line links or captions', async () => {
    const source = '<pre><span id="L99">original source line</span></pre>';
    await write('source/glviewer-ts/index.html', source);
    await write('glviewer/index.html', [
        '<a href="/doc/source/glviewer-ts/#L99"><code>GLViewer.ts:99</code></a>',
        '<a class="source" href="/doc/source/glviewer-ts/#L42">GLViewer.ts:42</a>',
        '<a href="https://example.com/doc/source/file/#L9">external.ts:9</a>',
        '<a href="/doc/glviewer/#L9">unrelated:9</a>'
    ].join('\n'));
    await finish(destination, { legacyLinks });
    const html = await read('glviewer/index.html');
    expect(html).toContain('<a href="/doc/source/glviewer-ts/"><code>GLViewer.ts</code></a>');
    expect(html).toContain('href="/doc/source/glviewer-ts/">GLViewer.ts</a>');
    expect(html).toContain('href="https://example.com/doc/source/file/#L9">external.ts:9');
    expect(html).toContain('href="/doc/glviewer/#L9">unrelated:9');
    expect(await read('source/glviewer-ts/index.html')).toBe(source);
});

test('cleans source links in all emitted Markdown source components', async () => {
    await write('glviewer/index.md', (await read('glviewer/index.md')) + [
        '<SourceLink href="/source/glviewer-ts/#L99" label="GLViewer.ts:99" />',
        '<MemberMeta sourceHref="/source/glviewer-ts/#L42" sourceLabel="GLViewer.ts:42" />',
        '<MemberCard\n sourceHref="/source/glviewer-ts/#L7"\n sourceLabel="GLViewer.ts:7" />',
        '<SourceLink href="https://example.com/source/#L99" label="external.ts:99" />'
    ].join('\n'));
    await finish(destination, { legacyLinks });
    const markdown = await read('glviewer/index.md');
    expect(markdown).toContain('href="/source/glviewer-ts/" label="GLViewer.ts"');
    expect(markdown).toContain('sourceHref="/source/glviewer-ts/" sourceLabel="GLViewer.ts"');
    expect(markdown).toContain('sourceHref="/source/glviewer-ts/"\n sourceLabel="GLViewer.ts"');
    expect(markdown).toContain('href="https://example.com/source/#L99" label="external.ts:99"');
});

test('repairs restored HTML assets and known legacy links on nested pages', async () => {
    await write('tutorials/code/index.md', '---\nkind: guide\n---\n[Viewer](GLViewer.html#addModel)');
    await write('tutorials/code/index.html', [
        '<script src="../build/3Dmol.js">const text="../tests/unchanged";</script>',
        '<iframe src="../tests/embed.html?x=1#demo"></iframe>',
        '<button data-pdb-url="../tutorials/doc-data/1ycr.pdb">Load</button>',
        '<a href="GLViewer.html#addModel">API</a>',
        '<a href="tutorial-code.html">Tutorial</a>',
        '<a href="./tutorial-code.html#Example">Relative tutorial</a>',
        '<a href="unknown.html">Unknown</a>',
        '<a href="https://example.com/GLViewer.html">External</a>',
        '<a href="/doc/glviewer/">New URL</a>',
        '<img src="//example.com/image.png">'
    ].join('\n'));
    await finish(destination, { legacyLinks });
    const html = await read('tutorials/code/index.html');
    expect(html).toContain('src="/build/3Dmol.js">const text="../tests/unchanged"');
    expect(html).toContain('src="/tests/embed.html?x=1#demo"');
    expect(html).toContain('data-pdb-url="/tutorials/doc-data/1ycr.pdb"');
    expect(html).toContain('href="/doc/GLViewer.html#addModel"');
    expect(html).toContain('href="/doc/tutorial-code.html"');
    expect(html).toContain('href="/doc/tutorial-code.html#Example"');
    expect(html).toContain('href="unknown.html"');
    expect(html).toContain('href="https://example.com/GLViewer.html"');
    expect(html).toContain('href="/doc/glviewer/"');
    expect(html).toContain('src="//example.com/image.png"');
    expect(await read('tutorials/code/index.md')).toContain('[Viewer](/doc/glviewer/index.md#addmodel)');
});

test('preserves fenced Markdown code while repairing links and source components outside fences', async () => {
    await write('tutorials/code/index.md', '---\nkind: guide\n---\n');
    const example = [
        '<SourceLink href="/source/viewer-ts/#L42" label="viewer.ts:42" />',
        '<MemberMeta sourceHref="/source/viewer-ts/#L99" sourceLabel="viewer.ts:99" />',
        '<a href="../tests/example.html">Example</a>',
        '[Tutorial](./tutorial-code.html)',
        '<a href="/doc/source/viewer-ts/#L7">viewer.ts:7</a>'
    ].join('\n');
    const fences = ['```mdx', '~~~~html', '  ````javascript'].map(opening => {
        const closing = opening.match(/^\s*([`~]+)/)[1];
        return `${opening}\n${example}\n${closing}\n`;
    }).join('\n');
    await write('glviewer/index.md', (await read('glviewer/index.md'))
        + fences + '\n' + example);
    await finish(destination, { legacyLinks });
    const markdown = await read('glviewer/index.md');
    expect(markdown).toContain(fences);
    expect(markdown).toContain('[Tutorial](/doc/tutorials/code/index.md)');
    expect(markdown).toContain('href="/source/viewer-ts/" label="viewer.ts"');
    expect(markdown).toContain('sourceHref="/source/viewer-ts/" sourceLabel="viewer.ts"');
    expect(markdown).toContain('href="/tests/example.html"');
    expect(markdown).toContain('href="/doc/source/viewer-ts/">viewer.ts</a>');
});

test('emits API and tutorial redirects, preserving queries and translating legacy anchors', async () => {
    await write('gradient-rwb/index.md', '---\nlongname: "Gradient.RWB"\n---\n');
    await write('gradient-rwb/index.html', '<h1>RWB</h1>');
    await write('tutorials/learning-environment/index.md', '---\nkind: guide\n---\n');
    await write('tutorials/learning-environment/index.html', '<p>Learning</p>');
    await write('index.html', '<a href="tutorial-learning_environment.html">Learning</a>');
    await finish(destination, { legacyLinks });
    for (const [file, target] of [
        ['GLViewer.html', '/doc/glviewer/'],
        ['Gradient.RWB.html', '/doc/gradient-rwb/'],
        ['global.html', '/doc/global/'],
        ['tutorial-learning-environment.html', '/doc/tutorials/learning-environment/'],
        ['tutorial-learning_environment.html', '/doc/tutorials/learning-environment/']
    ]) {
        const html = await read(file);
        const replace = jest.fn();
        vm.runInNewContext(html.match(/<script>([\s\S]*?)<\/script>/)[1], {
            location: { hash: '#.addModel', search: '?example=1', replace }
        });
        expect(replace).toHaveBeenCalledWith(`${target}?example=1#addmodel`);
        vm.runInNewContext(html.match(/<script>([\s\S]*?)<\/script>/)[1], {
            location: { hash: '#hover_callback', search: '', replace }
        });
        expect(replace).toHaveBeenLastCalledWith(`${target}#hovercallback`);
        expect(html).toContain(`<link rel="canonical" href="${target}">`);
    }
    expect(await read('index.html')).toContain('href="/doc/tutorial-learning_environment.html"');
});

test('rejects external and parent-directory alias names and leaves the homepage intact', async () => {
    await write('index.html', '<h1>Original homepage</h1>');
    await finish(destination, { legacyLinks: new Map([['GLViewer', '../escape.html']]) });
    expect(await read('index.html')).toBe('<h1>Original homepage</h1>');
    await expect(fs.access(path.join(destination, 'GLViewer.html'))).rejects.toThrow();
    await finish(destination, { legacyLinks: new Map([['GLViewer', 'https://example.com/bad.html']]) });
    expect(await read('glviewer/index.html')).toBe('<h1>GLViewer</h1>');
});

test('redirects globals that moved to standalone pages without replacing the global page', async () => {
    await write('colorspec/index.md', '---\nkind: typedef\nlongname: ColorSpec\n---\n');
    await write('colorspec/index.html', '<h1>ColorSpec</h1>');
    await write('xyz/index.md', '---\nkind: function\nlongname: XYZ\n---\n');
    await write('xyz/index.html', '<h1>XYZ</h1>');
    await finish(destination, {
        legacyLinks: new Map([
            ...legacyLinks,
            ['ColorSpec', 'global.html#ColorSpec'],
            ['XYZ', 'global.html#%58YZ']
        ])
    });
    const html = await read('global.html');
    expect(html).toContain('<link rel="canonical" href="/doc/global/">');
    const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
    for (const [hash, target] of [
        ['', '/doc/global/?example=1'],
        ['#ColorSpec', '/doc/colorspec/?example=1'],
        ['#Color%53pec', '/doc/colorspec/?example=1'],
        ['#XYZ', '/doc/xyz/?example=1'],
        ['#%58%59%5A', '/doc/xyz/?example=1'],
        ['#createViewer', '/doc/global/?example=1#createviewer']
    ]) {
        const replace = jest.fn();
        vm.runInNewContext(script, { location: { hash, search: '?example=1', replace } });
        expect(replace).toHaveBeenCalledWith(target);
    }
});


test('links Markdown to emitted companion files under a custom base path', async () => {
    await write('source/glviewer-ts/index.html', '<pre>source</pre>');
    await write('tutorials/code/index.md', '---\nkind: guide\n---\n# Code\n');
    await write('tutorials/code/index.html', '<h1>Code</h1>');
    const original = await read('glviewer/index.md');
    await write('glviewer/index.md', original + [
        '[Root route](/global#createviewer)',
        '[Prefixed route](/reference/api/global/?mode=plain#createviewer)',
        '[Sibling route](../global#createviewer)',
        '[Sibling file](../global/index.html#createviewer)',
        '[Current file](./index.md#render)',
        '[Legacy](GLViewer.html?mode=plain#addModel)',
        '[Absolute legacy](/reference/api/GLViewer.html#addModel)',
        '[Current anchor](#render)',
        '<SourceLink href="/source/glviewer-ts/#L42" label="GLViewer.ts:42" />',
        '<MemberMeta sourceHref="../source/glviewer-ts/#L42" sourceLabel="GLViewer.ts:42" />',
        '<a href="/source/glviewer-ts/index.html#L42">GLViewer.ts:42</a>',
        '[External](https://example.com/global#createviewer)',
        '[Protocol relative](//example.com/glviewer)',
        '[Email](mailto:author@example.com)',
        '[Unknown](/unrelated-page#keepCase)',
        '<a href="https://example.com/source/glviewer-ts/#L42">external.ts:42</a>'
    ].join('\n'));
    await write('tutorials/code/index.md', (await read('tutorials/code/index.md'))
        + '[Nested relative](../../glviewer/#render)\n');
    await finish(destination, { legacyLinks, basePath: '/reference/api/' });
    const markdown = await read('glviewer/index.md');
    for (const expected of [
        '[Root route](/reference/api/global/index.md#createviewer)',
        '[Prefixed route](/reference/api/global/index.md?mode=plain#createviewer)',
        '[Sibling route](/reference/api/global/index.md#createviewer)',
        '[Sibling file](/reference/api/global/index.md#createviewer)',
        '[Current file](/reference/api/glviewer/index.md#render)',
        '[Legacy](/reference/api/glviewer/index.md?mode=plain#addmodel)',
        '[Absolute legacy](/reference/api/glviewer/index.md#addmodel)',
        '[Current anchor](#render)',
        'href="/reference/api/source/glviewer-ts/index.html" label="GLViewer.ts"',
        'sourceHref="/reference/api/source/glviewer-ts/index.html" sourceLabel="GLViewer.ts"',
        'href="/reference/api/source/glviewer-ts/index.html">GLViewer.ts</a>',
        '[External](https://example.com/global#createviewer)',
        '[Protocol relative](//example.com/glviewer)',
        '[Email](mailto:author@example.com)',
        '[Unknown](/unrelated-page#keepCase)',
        'href="https://example.com/source/glviewer-ts/#L42">external.ts:42</a>'
    ]) expect(markdown).toContain(expected);
    expect(await read('tutorials/code/index.md')).toContain(
        '[Nested relative](/reference/api/glviewer/index.md#render)');
    // Every rewritten local file exists, and repeating postprocessing is stable.
    for (const match of markdown.matchAll(/\/reference\/api\/([^#?"\s)]+\.(?:md|html))/g)) {
        await expect(fs.access(path.join(destination, match[1]))).resolves.toBeUndefined();
    }
    await finish(destination, { legacyLinks, basePath: '/reference/api/' });
    expect(await read('glviewer/index.md')).toBe(markdown);
});

test('resolves legacy Markdown links to symbols with their own emitted pages', async () => {
    await write('colorspec/index.md', '---\nkind: typedef\nlongname: ColorSpec\n---\n');
    await write('colorspec/index.html', '<h1>ColorSpec</h1>');
    await write('glviewer/index.md', (await read('glviewer/index.md'))
        + '[Color](global.html?mode=plain#ColorSpec)\n');
    await finish(destination, { legacyLinks: new Map([
        ...legacyLinks, ['ColorSpec', 'global.html#ColorSpec']
    ]) });
    expect(await read('glviewer/index.md')).toContain('[Color](/doc/colorspec/index.md?mode=plain)');
});
