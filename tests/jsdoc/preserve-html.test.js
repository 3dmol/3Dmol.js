/** @jest-environment node */

const { mkdtemp, mkdir, readFile, writeFile, rm } = require('node:fs/promises');
const { tmpdir } = require('node:os');
const { join } = require('node:path');
const { JSDOM } = require('jsdom');
const { prepare } = require('../../jsdoc/preserve-html');

const tokensIn = html => html.match(/THREEDMOLHTML[a-f0-9]+END/g) || [];

describe('preserving authored HTML through the documentation theme', () => {
    let destination;

    beforeEach(async () => {
        destination = await mkdtemp(join(tmpdir(), '3dmol-doc-html-'));
    });

    afterEach(async () => {
        await rm(destination, { recursive: true, force: true });
    });

    test('restores interactive fragments and source code without changing the inputs', async () => {
        const opts = {
            readme: '<h2>Overview</h2><style>.viewer {height:250px}</style>' +
                '<script src="../build/3Dmol-min.js"></script>' +
                '<div class="viewer_3Dmoljs" data-pdb="1UBQ"><button>Try it</button></div>' +
                '<script>window.rotateViewer = () => 1;</script>' +
                '<iframe src="../tests/example.html" title="Interactive example"></iframe>' +
                '<p>Example caption</p><pre><code>&lt;script&gt;sample()&lt;/script&gt;</code></pre>'
        };
        const originalReadme = opts.readme;
        const tutorial = {
            name: 'code', title: 'Using code', type: 2, content: 'raw markdown', children: [],
            parse() {
                return '<p>Tutorial caption</p><div id="container-01"><button id="btn-01">Try it</button></div>' +
                    '<script>document.querySelector("#btn-01").onclick = example;</script>' +
                    '<pre><code>&lt;div id="example"&gt;&lt;/div&gt;</code></pre>';
            }
        };
        const tutorials = { children: [tutorial] };
        const prepared = prepare(opts, tutorials);
        const readmeTokens = tokensIn(prepared.opts.readme);
        const tutorialTokens = tokensIn(prepared.tutorials.children[0].content);

        expect(opts.readme).toBe(originalReadme);
        expect(tutorial.content).toBe('raw markdown');
        expect(tutorial.type).toBe(2);
        expect(prepared.tutorials).not.toBe(tutorials);
        expect(prepared.tutorials.children[0]).not.toBe(tutorial);
        expect(prepared.tutorials.children[0].type).toBe(1);
        expect(readmeTokens).toHaveLength(5);
        expect(tutorialTokens).toHaveLength(2);
        expect(prepared.opts.readme).toContain('&lt;script&gt;sample()&lt;/script&gt;');
        expect(prepared.opts.readme).toContain('Example caption');

        await mkdir(join(destination, 'tutorials', 'code'), { recursive: true });
        await mkdir(join(destination, '_assets'));
        await writeFile(join(destination, 'index.html'),
            `<html><head><meta name="description" content="${readmeTokens.join(' ')}"></head>` +
            `<body><main>${prepared.opts.readme}</main></body></html>`);
        await writeFile(join(destination, 'index.md'), readmeTokens.join('\n\n'));
        await writeFile(join(destination, 'tutorials', 'code', 'index.html'), prepared.tutorials.children[0].content);
        await writeFile(join(destination, '_assets', 'search-index.json'), JSON.stringify({
            entries: [{ text: readmeTokens.join(' ') + tutorialTokens.join(' ') }]
        }));

        await prepared.restore(destination);

        const output = await readFile(join(destination, 'index.html'), 'utf8');
        const document = new JSDOM(output).window.document;
        expect(document.querySelector('style').textContent).toBe('.viewer {height:250px}');
        expect(document.querySelector('script[src]').getAttribute('src')).toBe('../build/3Dmol-min.js');
        expect(document.querySelector('.viewer_3Dmoljs').getAttribute('data-pdb')).toBe('1UBQ');
        expect(document.querySelector('.viewer_3Dmoljs button').textContent).toBe('Try it');
        expect(document.querySelector('iframe').getAttribute('src')).toBe('../tests/example.html');
        expect(document.querySelector('pre code').textContent).toBe('<script>sample()</script>');
        expect(document.querySelector('meta').getAttribute('content')).toContain('Try it');
        expect(document.querySelector('meta').getAttribute('content')).not.toContain('<script');
        expect(output).not.toContain('THREEDMOLHTML');
        expect(await readFile(join(destination, 'index.md'), 'utf8')).toContain('<iframe');
        const tutorialOutput = await readFile(join(destination, 'tutorials', 'code', 'index.html'), 'utf8');
        expect(tutorialOutput).toContain('id="btn-01"');
        expect(tutorialOutput).toContain('document.querySelector("#btn-01").onclick = example;');
        expect(tutorialOutput).toContain('Tutorial caption');
        const search = JSON.parse(await readFile(join(destination, '_assets', 'search-index.json'), 'utf8'));
        expect(search.entries[0].text).toContain('Try it');
        expect(search.entries[0].text).not.toMatch(/THREEDMOLHTML|rotateViewer|<script>/);
    });

    test('leaves prose and fenced HTML alone when no live elements exist', async () => {
        const html = '<p>A caption</p><pre><code>&lt;div&gt;Example&lt;/div&gt;</code></pre>';
        const tutorial = { name: 'plain', content: 'Original Markdown', type: 2, parse: () => html };
        const prepared = prepare({ readme: html }, { children: [tutorial] });
        expect(prepared.opts.readme).toBe(html);
        expect(prepared.tutorials.children[0].content).toBe('Original Markdown');
        expect(prepared.tutorials.children[0].type).toBe(2);
        await expect(prepared.restore(destination)).resolves.toBeUndefined();
    });

    test('fails the build if the theme silently drops a protected fragment', async () => {
        const prepared = prepare({ readme: '<div class="viewer_3Dmoljs"></div>' });
        await writeFile(join(destination, 'index.html'), '<p>Viewer is missing</p>');
        await expect(prepared.restore(destination)).rejects.toThrow('did not render 1 protected HTML fragment');
    });
});
