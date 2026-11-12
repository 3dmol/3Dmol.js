// Run after building the bundles and docs:
// node tests/jsdoc/browser.cjs [documentation-output-directory]
// Homepage assets/data are intercepted; tutorial demos use their real bundled
// library and protein. Analytics never leaves the browser.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const puppeteer = require('puppeteer');
const { checkTutorial } = require('./tutorial-browser.cjs');

const root = path.resolve(__dirname, '../..');
const destination = path.resolve(process.argv[2] || path.join(root, 'doc'));
const pdb = fs.readFileSync(path.join(root, 'tests/test_structs/1YCR.pdb'));
const mime = {
    '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css',
    '.json': 'application/json', '.svg': 'image/svg+xml', '.png': 'image/png',
    '.woff2': 'font/woff2', '.wasm': 'application/wasm', '.md': 'text/plain'
};
const bundles = new Map(['3Dmol.js', '3Dmol-min.js', '3Dmol.ui.js', '3Dmol.ui-min.js'].map(name => [
    '/build/' + name,
    Buffer.concat([
        fs.readFileSync(path.join(root, 'build', name)),
        Buffer.from(/^3Dmol(?:-min)?\.js$/.test(name) ? '\n;$3Dmol.setSyncSurface(true);\n' : '')
    ])
]));
assert.ok(fs.existsSync(path.join(destination, 'index.html')), `Build documentation first: ${destination}`);

const server = http.createServer((request, response) => {
    const url = new URL(request.url, 'http://localhost');
    let pathname;
    try { pathname = decodeURIComponent(url.pathname); }
    catch { response.writeHead(400).end(); return; }
    if (pathname === '/favicon.ico') { response.writeHead(204).end(); return; }
    if (pathname === '/doc') { response.writeHead(301, { Location: '/doc/' }).end(); return; }
    let base, relative;
    if (pathname.startsWith('/doc/')) {
        base = destination;
        relative = pathname.slice('/doc/'.length);
    } else if (pathname.startsWith('/build/') || pathname.startsWith('/tests/') || pathname.startsWith('/tutorials/doc-data/')) {
        base = root;
        relative = pathname.slice(1);
    } else {
        response.writeHead(404).end('Unknown test route');
        return;
    }
    let file = path.resolve(base, relative);
    if (file !== base && !file.startsWith(base + path.sep)) {
        response.writeHead(403).end();
        return;
    }
    try {
        if (fs.statSync(file).isDirectory()) file = path.join(file, 'index.html');
        const body = bundles.get(pathname) || fs.readFileSync(file);
        response.writeHead(200, { 'Content-Type': mime[path.extname(file)] || 'application/octet-stream' });
        response.end(request.method === 'HEAD' ? undefined : body);
    } catch {
        response.writeHead(404).end(`Missing ${pathname}`);
    }
});

(async () => {
    await new Promise((resolve, reject) => {
        server.once('error', reject);
        server.listen(0, '127.0.0.1', resolve);
    });
    const origin = `http://127.0.0.1:${server.address().port}`;
    let browser;
    try {
        browser = await puppeteer.launch({
            headless: true,
            args: ['--no-sandbox', '--enable-unsafe-swiftshader']
        });
        const page = await browser.newPage();
        await page.setViewport({ width: 1200, height: 900 });
        page.setDefaultTimeout(30000);
        page.setDefaultNavigationTimeout(60000);
        const errors = [], failedRequests = [], unexpectedRequests = [], missing = [];
        const assets = new Set(), links = new Set(), skippedFrames = new Set();
        let analyticsIntercepted = 0;
        page.on('pageerror', error => errors.push(`${page.url()}: ${error.message}`));
        page.on('requestfailed', request => {
            if (request.url().startsWith(origin + '/doc/')) {
                failedRequests.push(`${request.url()}: ${request.failure()?.errorText}`);
            }
        });
        page.on('response', response => {
            if (!response.url().startsWith(origin)) return;
            if (response.status() >= 400) missing.push(`${response.status()} ${response.url()}`);
            if (/\/doc\/(?:_assets|_islands)\//.test(response.url())) assets.add(response.url());
        });
        await page.setRequestInterception(true);
        page.on('request', async request => {
            try {
                const url = new URL(request.url());
                if (request.isNavigationRequest() && request.frame()?.parentFrame()) {
                    // The preserved example iframe pins jQuery 3.6.3 by SRI;
                    // validate its real URL separately without replacing that
                    // dependency with a different installed jQuery version.
                    skippedFrames.add(url.href);
                    await request.respond({ status: 200, contentType: 'text/html', body: '<!doctype html><title>Offline iframe fixture</title>' });
                } else if (url.origin === origin || ['blob:', 'data:'].includes(url.protocol)) {
                    await request.continue();
                } else if (url.hostname.toLowerCase() === '3dmol.org' && bundles.has(url.pathname)) {
                    await request.respond({ status: 200, contentType: 'text/javascript', body: bundles.get(url.pathname) });
                } else if (url.hostname === 'files.rcsb.org' && /\.pdb$/i.test(url.pathname)) {
                    await request.respond({ status: 200, contentType: 'text/plain', headers: { 'Access-Control-Allow-Origin': '*' }, body: pdb });
                } else if (/(^|\.)(google-analytics\.com|googletagmanager\.com)$/.test(url.hostname)) {
                    analyticsIntercepted++;
                    await request.respond({ status: 200, contentType: 'text/javascript', body: '/* Analytics disabled in documentation tests. */' });
                } else if (url.hostname === 'fonts.googleapis.com') {
                    await request.respond({ status: 200, contentType: 'text/css', body: '' });
                } else {
                    unexpectedRequests.push(url.href);
                    await request.abort();
                }
            } catch (error) {
                errors.push(`Request interception: ${error.message}`);
                if (!request.isInterceptResolutionHandled()) await request.abort();
            }
        });

        async function collectLinks() {
            for (const href of await page.$$eval('a[href]', nodes => nodes.map(node => node.href))) {
                const url = new URL(href);
                if (url.origin === origin && url.pathname.startsWith('/doc/')) {
                    url.hash = '';
                    links.add(url.href);
                }
            }
        }

        async function waitForViewers(count) {
            await page.waitForFunction(expected => {
                const viewers = Object.values(window.$3Dmol?.viewers || {});
                const hosts = [...document.querySelectorAll('.viewer_3Dmoljs')];
                return hosts.length === expected && viewers.length === expected &&
                    viewers.every(viewer => viewer.selectedAtoms({}).length > 0) &&
                    hosts.every(host => [...host.querySelectorAll('canvas')].some(canvas => canvas.width && canvas.height));
            }, {}, count);
        }

        await page.goto(origin + '/doc/', { waitUntil: 'networkidle0' });
        await waitForViewers(1);
        assert.equal(await page.$eval('.viewer_3Dmoljs', node => node.dataset.pdb), '1UBQ');
        const iframe = await page.$eval('iframe', node => node.src);
        assert.equal(iframe, origin + '/tests/example.html', 'Homepage iframe was not rebased');
        const iframeResponse = await fetch(iframe);
        assert.equal(iframeResponse.status, 200, 'Homepage iframe no longer resolves');
        assert.match(await iframeResponse.text(), /id="gldiv"/);
        assert.ok(skippedFrames.has(iframe));
        assert.ok(analyticsIntercepted > 0, 'Expected homepage analytics to be intercepted');
        await collectLinks();

        // Opening search exercises the theme's real hydrated controls and its
        // embedded search index instead of just its static HTML.
        await page.click('[data-island="cmdk"] button[aria-label="Search"]');
        await page.waitForSelector('input[type="search"]', { visible: true });
        await page.type('input[type="search"]', 'GLViewer');
        await page.waitForFunction(() => [...document.querySelectorAll('[role="dialog"] a')]
            .some(link => /glviewer/i.test(link.href)));
        await page.keyboard.press('Escape');

        await page.goto(origin + '/doc/tutorials/embeddable/', { waitUntil: 'networkidle0' });
        await waitForViewers(2);
        await page.waitForFunction(() => {
            const viewer = Object.values(window.$3Dmol.viewers).find(item => item.ui);
            return viewer?.ui?.ui?.element?.isConnected &&
                document.querySelector('.viewer_3Dmoljs[data-ui="true"] .molui');
        });
        assert.equal(await page.$$eval('.viewer_3Dmoljs[data-ui="true"] .molui', nodes => nodes.length), 1);
        assert.ok(await page.$$eval('pre code', nodes => nodes.some(node => node.textContent.includes('data-ui='))),
            'The embedding example source was lost');
        await collectLinks();

        await checkTutorial(page, origin + '/doc/', 'HTTP /doc/');
        await collectLinks();

        await page.goto(origin + '/doc/GLViewer.html#addModel', { waitUntil: 'networkidle0' });
        await page.waitForFunction(() => location.pathname === '/doc/glviewer/index.html' && location.hash === '#addmodel');
        assert.ok(await page.$('#addmodel'), 'Legacy member anchor did not resolve');
        const sources = await page.$$eval('a[href*="/source/"]', nodes => nodes.map(node => ({
            href: node.href, text: node.textContent.trim()
        })));
        assert.ok(sources.length > 0, 'The API page has no source links');
        for (const source of sources) {
            assert.doesNotMatch(source.href, /#L\d+$/, 'Source link still contains a line anchor');
            assert.doesNotMatch(source.text, /:\d+$/, 'Source label still contains a line number');
            links.add(source.href);
        }
        assert.ok(sources.some(source => /GLViewer\.ts/.test(source.text)), 'Source filenames are missing');
        await collectLinks();

        for (const [legacy, pathname, hash] of [
            ['AtomSpec.html#hover_callback', '/doc/atomspec/index.html', '#hovercallback'],
            ['global.html', '/doc/global/index.html', ''],
            ['global.html#createViewer', '/doc/global/index.html', '#createviewer'],
            ['global.html#ColorSpec', '/doc/colorspec/index.html', '']
        ]) {
            await page.goto(origin + '/doc/' + legacy, { waitUntil: 'networkidle0' });
            await page.waitForFunction((pathname, hash) => location.pathname === pathname && location.hash === hash,
                {}, pathname, hash);
            if (hash) assert.ok(await page.$(hash), `Legacy anchor did not resolve: ${legacy}`);
        }

        // Follow every internal link exposed on the visited pages, including
        // sidebar targets and source files. This checks nested routes without
        // loading unrelated source-editor CDN dependencies.
        await Promise.all([...links].map(async href => {
            const response = await fetch(href);
            assert.equal(response.status, 200, `Broken documentation link: ${href}`);
            await response.arrayBuffer();
        }));
        assert.ok([...assets].some(url => /\/_assets\/runtime-[^/]+\.js$/.test(url)), 'Theme runtime was never loaded');
        assert.deepEqual(unexpectedRequests, [], 'Unexpected external network request');
        assert.deepEqual(missing, [], 'Documentation requested missing local resources');
        assert.deepEqual(failedRequests, [], 'Documentation assets failed to load');
        assert.deepEqual(errors, [], 'Browser JavaScript errors');
        console.log(`Documentation browser smoke passed: homepage viewer, embedded viewer UI, search, source links, legacy member redirects, all three code tutorial buttons; ${links.size} links and ${assets.size} assets checked. Pinned-jQuery iframe execution is excluded.`);
    } finally {
        if (browser) await browser.close();
        await new Promise(resolve => server.close(resolve));
    }
})().catch(error => { console.error(error); process.exitCode = 1; });
