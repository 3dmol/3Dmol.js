// Run after npm run doc. No web-security or file-access browser overrides.
// node tests/jsdoc/local-browser.cjs [documentation-output-directory]
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { fileURLToPath, pathToFileURL } = require('node:url');
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

function readServedFile(file) {
    if (fs.statSync(file).isDirectory()) file = path.join(file, 'index.html');
    return { body: fs.readFileSync(file), type: mime[path.extname(file)] || 'application/octet-stream' };
}

async function run(browser, baseUrl, label) {
    const page = await browser.newPage();
    await page.setViewport({ width: 1200, height: 900 });
    page.setDefaultTimeout(15000);
    const failures = [], errors = [], unexpected = [];
    const local = url => url.protocol === 'file:' || url.origin === new URL(baseUrl).origin;
    page.on('pageerror', error => errors.push(error.message));
    page.on('requestfailed', request => {
        if (local(new URL(request.url()))) failures.push(`${request.url()}: ${request.failure()?.errorText}`);
    });
    page.on('response', response => {
        if (local(new URL(response.url())) && response.status() >= 400) failures.push(`${response.status()} ${response.url()}`);
    });
    await page.setRequestInterception(true);
    page.on('request', async request => {
        try {
            const url = new URL(request.url());
            if (request.isNavigationRequest() && request.frame()?.parentFrame()) {
                // This separate legacy example pins another jQuery version.
                await request.respond({ status: 200, contentType: 'text/html', body: '<!doctype html><title>Iframe fixture</title>' });
            } else if (/(^|\.)(google-analytics\.com|googletagmanager\.com)$/.test(url.hostname)) {
                // A protocol-relative legacy analytics URL becomes file://host
                // when the homepage is opened from disk. Intercept it too.
                await request.respond({ status: 200, contentType: 'text/javascript', body: '/* Analytics disabled. */' });
            } else if (local(url) || ['data:', 'blob:'].includes(url.protocol)) {
                await request.continue();
            } else if (url.hostname.toLowerCase() === '3dmol.org' && bundles.has(url.pathname)) {
                await request.respond({ status: 200, contentType: 'text/javascript', body: bundles.get(url.pathname) });
            } else if (url.hostname === 'files.rcsb.org' && /\.pdb$/i.test(url.pathname)) {
                await request.respond({ status: 200, contentType: 'text/plain', headers: { 'Access-Control-Allow-Origin': '*' }, body: pdb });
            } else if (url.hostname === 'fonts.googleapis.com') {
                await request.respond({ status: 200, contentType: 'text/css', body: '' });
            } else {
                unexpected.push(url.href);
                await request.abort();
            }
        } catch (error) {
            errors.push(`Request interception: ${error.message}`);
            if (!request.isInterceptResolutionHandled()) await request.abort();
        }
    });

    async function assertStyled() {
        const styling = await page.evaluate(() => ({
            margin: getComputedStyle(document.body).margin,
            sheets: [...document.querySelectorAll('link[rel="stylesheet"]')].map(link => ({ href: link.href, loaded: !!link.sheet }))
        }));
        assert.ok(styling.sheets.some(sheet => /_assets\/styles\./.test(sheet.href) && sheet.loaded),
            `${label}: theme stylesheet failed to load: ${JSON.stringify(styling.sheets)}; ${failures.join('; ')}`);
        assert.equal(styling.margin, '0px', `${label}: theme CSS did not style the page`);
    }

    async function assertTarget(url) {
        const target = new URL(url);
        assert.ok(target.href.startsWith(baseUrl), `${label}: documentation link escaped its deployment root: ${target.href}`);
        if (target.protocol === 'file:') {
            assert.equal(fs.statSync(fileURLToPath(target)).isFile(), true,
                `${label}: a file link points to a directory instead of index.html: ${target.href}`);
        } else {
            const response = await fetch(target);
            assert.equal(response.status, 200, `${label}: broken documentation link: ${target.href}`);
        }
    }

    async function assertMenu(selector) {
        const links = await page.$$eval(`${selector} a`, nodes => nodes.map(node => ({
            title: node.textContent.trim(), target: node.target, visible: !!node.getClientRects().length
        })));
        assert.deepEqual(links.map(link => link.title), ['Documentation', 'Tutorials', 'GitHub'], `${label}: menu links`);
        assert.ok(links.every(link => link.visible), `${label}: menu links must be visible`);
        assert.deepEqual(links.map(link => link.target), ['_self', '_self', '_blank'], `${label}: menu link targets`);
    }

    async function clickMenuLink(selector, index) {
        const links = await page.$$(`${selector} a`);
        const href = await links[index].evaluate(link => link.href);
        await assertTarget(href);
        await Promise.all([page.waitForNavigation({ waitUntil: 'networkidle0' }), links[index].click()]);
        assert.equal(page.url(), href, `${label}: menu link must navigate in the same tab`);
    }

    try {
        await page.goto(new URL('index.html', baseUrl).href, { waitUntil: 'networkidle0' });
        await assertStyled();
        const desktopMenu = 'header .docs-top-menu';
        await assertMenu(desktopMenu);
        assert.equal(await page.$(`${desktopMenu} svg`), null, `${label}: top menu should not show link-out icons`);
        assert.equal(await page.$eval('[data-island="sidebar"] nav > ul', list => getComputedStyle(list).display),
            'none', `${label}: desktop sidebar repeats the top menu`);
        const viewerLink = await page.$$eval('a[href]', nodes => nodes.find(node => /^GLViewer$/.test(node.textContent.trim()))?.href);
        assert.ok(viewerLink, `${label}: GLViewer navigation link is missing`);
        await assertTarget(viewerLink);
        await page.goto(viewerLink, { waitUntil: 'networkidle0' });
        await assertStyled();
        assert.ok(await page.$('#addmodel'), `${label}: GLViewer methods did not render`);
        const memberLink = await page.$$eval('a[href]', nodes => nodes.find(node => /#addmodel$/.test(node.href))?.href);
        assert.ok(memberLink, `${label}: member link is missing`);
        await assertTarget(memberLink);
        await page.goto(memberLink, { waitUntil: 'networkidle0' });
        assert.equal(new URL(page.url()).hash, '#addmodel');
        const source = await page.$$eval('a[href]', nodes => nodes.find(node => /\/source\//.test(node.href))?.href);
        assert.ok(source, `${label}: source link is missing`);
        await assertTarget(source);

        await page.click('[data-island="cmdk"] button[aria-label="Search"]');
        await page.waitForSelector('input[type="search"]', { visible: true });
        await page.type('input[type="search"]', 'GLModel');
        await page.waitForFunction(() => [...document.querySelectorAll('[role="dialog"] a')].some(link => /glmodel/i.test(link.href)));
        const result = await page.$$eval('[role="dialog"] a', nodes => nodes.find(node => /glmodel/i.test(node.href))?.href);
        await assertTarget(result);
        await page.keyboard.press('Escape');
        await page.waitForSelector('[role="dialog"]', { hidden: true });

        await assertMenu(desktopMenu);
        await clickMenuLink(desktopMenu, 1);
        await assertStyled();
        assert.ok(await page.$('main'), `${label}: tutorial did not render`);
        await clickMenuLink(desktopMenu, 0);
        assert.equal(page.url(), new URL('index.html', baseUrl).href, `${label}: Documentation should open the local homepage`);

        await page.setViewport({ width: 390, height: 844 });
        assert.equal(await page.$eval(desktopMenu, menu => getComputedStyle(menu).display), 'none',
            `${label}: desktop menu should be hidden on mobile`);
        const mobileMenu = '[role="dialog"] nav[aria-label="Documentation navigation"] > ul:first-child';
        for (const index of [1, 0]) {
            await page.click('[data-island="mobile-nav"] button[aria-haspopup="dialog"]');
            await page.waitForSelector(mobileMenu, { visible: true });
            await assertMenu(mobileMenu);
            assert.ok(await page.$eval('header', header => header.scrollWidth <= header.clientWidth),
                `${label}: mobile header overflows`);
            assert.ok(await page.$eval('[role="dialog"]', dialog => dialog.scrollWidth <= dialog.clientWidth),
                `${label}: mobile menu overflows`);
            await clickMenuLink(mobileMenu, index);
        }
        assert.equal(page.url(), new URL('index.html', baseUrl).href, `${label}: mobile Documentation should open the local homepage`);
        await page.setViewport({ width: 1200, height: 900 });

        await checkTutorial(page, baseUrl, label);

        await page.goto(new URL('GLViewer.html#addModel', baseUrl).href, { waitUntil: 'networkidle0' });
        await page.waitForFunction(() => location.hash === '#addmodel' && !!document.getElementById('addmodel'));
        await assertTarget(page.url());
        assert.deepEqual(failures, [], `${label}: local resources failed`);
        assert.deepEqual(errors, [], `${label}: browser JavaScript errors`);
        assert.deepEqual(unexpected, [], `${label}: unexpected external network request`);
        console.log(`${label}: stylesheet, desktop/mobile menus, member/source/tutorial links, search and legacy redirect passed.`);
    } finally {
        await page.close();
    }
}

(async () => {
    let browser;
    let server;
    try {
        browser = await puppeteer.launch({ headless: true, args: ['--no-sandbox', '--enable-unsafe-swiftshader'] });
        await run(browser, pathToFileURL(destination + path.sep).href, 'file://');
        for (const mount of ['/preview/docs/', '/']) {
            server = http.createServer((request, response) => {
                const url = new URL(request.url, 'http://localhost');
                if (url.pathname === '/favicon.ico') { response.writeHead(204).end(); return; }
                const relative = decodeURIComponent(url.pathname).slice(mount.length);
                const file = path.resolve(destination, relative);
                if (!url.pathname.startsWith(mount) || (file !== destination && !file.startsWith(destination + path.sep))) {
                    response.writeHead(404).end();
                    return;
                }
                try {
                    const { body, type } = readServedFile(file);
                    response.writeHead(200, { 'Content-Type': type }).end(body);
                } catch { response.writeHead(404).end(); }
            });
            await new Promise((resolve, reject) => {
                server.once('error', reject);
                server.listen(0, '127.0.0.1', resolve);
            });
            await run(browser, `http://127.0.0.1:${server.address().port}${mount}`, `HTTP ${mount}`);
            await new Promise(resolve => server.close(resolve));
            server = null;
        }
    } finally {
        if (browser) await browser.close();
        if (server) await new Promise(resolve => server.close(resolve));
    }
})().catch(error => { console.error(error); process.exitCode = 1; });
