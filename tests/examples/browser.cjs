// npm run test:examples:browser (after building the core and UI bundles).
// The Python example is executed directly from its README; no network is used.
'use strict';

const assert = require('node:assert/strict');
const { execFileSync } = require('node:child_process');
const fs = require('node:fs');
const path = require('node:path');
const puppeteer = require('puppeteer');

const root = path.resolve(__dirname, '../..');
const core = path.join(root, 'build/3Dmol.js');
const ui = path.join(root, 'build/3Dmol.ui.js');

async function pythonGrid(browser) {
    const exported = JSON.parse(execFileSync('python3', [path.join(__dirname, 'export_python_grid.py')], {
        cwd: root, encoding: 'utf8'
    }));
    const libraryURL = exported.html.match(/loadScriptAsync\('([^']+)'\)/)?.[1];
    assert.ok(libraryURL, 'Python export should load its configured JavaScript bundle');
    const page = await browser.newPage();
    try {
        const errors = [], requests = [];
        page.on('pageerror', error => errors.push(error.message));
        await page.setViewport({ width: 640, height: 320, deviceScaleFactor: 1 });
        await page.setRequestInterception(true);
        page.on('request', async request => {
            requests.push(request.url());
            if (request.url() === libraryURL) {
                await request.respond({ status: 200, contentType: 'text/javascript', body: fs.readFileSync(core) });
            } else await request.abort();
        });
        await page.setContent(exported.html);
        await page.waitForFunction(id => window['viewer_' + id]?.selectedAtoms({}).length === 2, {}, exported.id);
        const result = await page.evaluate(id => {
            const first = window['viewer_' + id];
            const [left, right] = first.getConfig().viewers[0];
            const canvas = left.getCanvas();
            const image = document.createElement('canvas');
            image.width = canvas.width; image.height = canvas.height;
            const context = image.getContext('2d');
            context.drawImage(canvas, 0, 0);
            const pixels = context.getImageData(0, 0, image.width, image.height).data;
            const colors = [{ red: 0, blue: 0 }, { red: 0, blue: 0 }];
            for (let y = 0; y < image.height; y++) for (let x = 0; x < image.width; x++) {
                const offset = 4 * (y * image.width + x);
                const [r, g, b, a] = pixels.subarray(offset, offset + 4);
                const half = colors[x < image.width / 2 ? 0 : 1];
                if (a > 200 && r > g + 50 && r > b + 50) half.red++;
                if (a > 200 && b > g + 50 && b > r + 50) half.blue++;
            }
            return {
                atoms: [left, right].map(viewer => viewer.selectedAtoms({}).length),
                styles: [left, right].map(viewer => viewer.selectedAtoms({}).map(atom => atom.style)),
                linked: [left, right].map(viewer => viewer.getConfig().control_all),
                canvases: document.querySelectorAll('canvas').length,
                sharedCanvas: left.getCanvas() === right.getCanvas(),
                sharedCapture: left.pngURI() === right.pngURI(), colors
            };
        }, exported.id);
        assert.deepEqual(result.atoms, [2, 2], 'Unqualified addModel should load both grid cells');
        assert.deepEqual(result.styles, [
            [{ sphere: { radius: 0.45, color: 'red' } }, { sphere: { radius: 0.45, color: 'red' } }],
            [{ stick: { radius: 0.15, color: 'blue' } }, { stick: { radius: 0.15, color: 'blue' } }]
        ], 'getModel(viewer=...) should style only its own cell');
        assert.deepEqual(result.linked, [false, false]);
        assert.equal(result.canvases, 1);
        assert.equal(result.sharedCanvas, true);
        assert.equal(result.sharedCapture, true);
        assert.ok(result.colors[0].red > 50 && result.colors[0].blue === 0,
            `Left cell should visibly render red spheres: ${JSON.stringify(result.colors)}`);
        assert.ok(result.colors[1].blue > 50 && result.colors[1].red === 0,
            `Right cell should visibly render blue sticks: ${JSON.stringify(result.colors)}`);
        assert.deepEqual(requests, [libraryURL], 'The export must not fetch molecule data or other external assets');
        assert.deepEqual(errors, []);
        console.log('Python README grid: exported HTML, cell-specific models/styles, shared canvas, and rendered colors passed.');
    } finally { await page.close(); }
}

async function viewerLifecycle(browser) {
    const page = await browser.newPage();
    try {
        const errors = [];
        page.on('pageerror', error => errors.push(error.message));
        await page.setViewport({ width: 800, height: 600, deviceScaleFactor: 1 });
        await page.setContent('<div id="viewer" style="width:400px;height:300px;position:static"></div>');
        await page.addScriptTag({ path: core });
        await page.addScriptTag({ path: ui });
        await page.evaluate(() => {
            // Track real polling timers so teardown cannot pass just by hiding controls.
            const set = window.setInterval.bind(window), clear = window.clearInterval.bind(window);
            window.intervals = new Set();
            window.setInterval = (callback, delay, ...args) => {
                const id = set(callback, delay, ...args); intervals.add(id); return id;
            };
            window.clearInterval = id => { intervals.delete(id); clear(id); };
            window.host = document.getElementById('viewer');
            window.viewer = $3Dmol.createViewer(host, { antialias: false, upscale: false });
            window.canvas = viewer.getCanvas();
            window.xyz = '1\nFrame 0\nC 0 0 0\n1\nFrame 1\nC 2 0 0\n';
            viewer.addModelsAsFrames(xyz, 'xyz');
            viewer.setStyle({}, { sphere: { color: 'blue' } });
            viewer.zoomTo(); viewer.render();
            window.controls = new $3Dmol.StateManager(viewer, { open: true });
            window.panel = controls.ui.element;
            window.initialIntervals = intervals.size;
            host.style.width = '520px'; host.style.height = '240px';
            viewer.resize();
        });
        await page.waitForFunction(() => canvas.width === 520 && canvas.height === 240);
        const reused = await page.evaluate(() => {
            viewer.clear();
            const empty = viewer.selectedAtoms({}).length === 0;
            viewer.addModelsAsFrames(xyz, 'xyz');
            viewer.setStyle({}, { sphere: { color: 'green' } });
            viewer.zoomTo(); viewer.render(); controls.refresh();
            return {
                empty, alive: !viewer.isDestroyed(), sameCanvas: canvas === viewer.getCanvas(),
                sameControls: panel === viewer.ui.ui.element, panels: host.querySelectorAll('.molui').length,
                frames: viewer.getNumFrames(), polling: intervals.size === initialIntervals && initialIntervals > 0
            };
        });
        assert.deepEqual(reused, {
            empty: true, alive: true, sameCanvas: true, sameControls: true, panels: 1, frames: 2, polling: true
        }, 'clear should allow reuse without replacing the canvas or duplicating UI resources');
        const detached = await page.evaluate(async () => {
            const setFrame = viewer.setFrame.bind(viewer);
            let release;
            viewer.setFrame = frame => new Promise(resolve => { release = () => setFrame(frame).then(resolve); });
            const pending = controls.ui.setTrajectorySettings({ frame: 1, playing: true });
            let draws = 0;
            viewer.setViewChangeCallback(() => { draws++; });
            controls.destroy();
            const retired = { disconnected: !panel.isConnected, timers: intervals.size, detached: viewer.ui === undefined };
            viewer.setFrame = setFrame;
            release();
            await pending;
            retired.lateDraws = draws;
            retired.playing = controls.ui.getTrajectorySettings().playing;
            viewer.setViewChangeCallback(undefined);
            window.replacementControls = new $3Dmol.StateManager(viewer, { open: true });
            retired.replacementPanels = host.querySelectorAll('.molui').length;
            retired.replacementTimers = intervals.size === initialIntervals;
            viewer.destroy();
            retired.destroyed = viewer.isDestroyed();
            retired.canvasRemoved = !canvas.isConnected;
            retired.controlsRemoved = !replacementControls.ui.element.isConnected;
            retired.allTimersReleased = intervals.size === 0;
            retired.restoredPosition = host.style.position;
            return retired;
        });
        assert.deepEqual(detached, {
            disconnected: true, timers: 0, detached: true, lateDraws: 0, playing: false,
            replacementPanels: 1, replacementTimers: true, destroyed: true,
            canvasRemoved: true, controlsRemoved: true, allTimersReleased: true, restoredPosition: 'static'
        }, 'Disposed controls must not resume pending playback or retain polling, and viewer destruction must release the replacement');
        await page.evaluate(() => {
            window.next = $3Dmol.createViewer(host, { antialias: false, upscale: false });
            next.addModel('1\nReplacement\nO 0 0 0\n', 'xyz');
            next.setStyle({}, { sphere: { color: 'red' } });
            next.zoomTo(); next.render();
            window.nextControls = new $3Dmol.StateManager(next);
            if (host.querySelectorAll('canvas').length !== 1 || host.querySelectorAll('.molui').length !== 1) {
                throw new Error('Recreating the viewer duplicated DOM resources');
            }
            host.remove();
        });
        await page.waitForFunction(() => next.isDestroyed() && intervals.size === 0);
        assert.equal(await page.evaluate(() => nextControls.ui.element.isConnected), false);
        assert.deepEqual(errors, []);
        console.log('Viewer lifecycle: CSS resize, clear/reuse, pending UI teardown, destroy/recreate, and automatic DOM cleanup passed.');
    } finally { await page.close(); }
}

async function reusedViewerBackground(browser) {
    // Run the same generated function and existing-viewer path as tests.html.
    // glcheck owns reference-image comparisons; this checks the background bug.
    execFileSync('python3', ['tests/auto/generate_tests.py'], { cwd: root });
    const page = await browser.newPage();
    try {
        const errors = [];
        page.on('pageerror', error => errors.push(error.message));
        await page.setViewport({ width: 400, height: 400, deviceScaleFactor: 1 });
        await page.setContent('<div id="viewer" style="width:400px;height:400px;position:relative"></div>');
        await page.addScriptTag({ path: path.join(root, 'node_modules/jquery/dist/jquery.min.js') });
        await page.addScriptTag({ path: core });
        await page.addScriptTag({ path: path.join(root, 'tests/auto/build.js') });
        const result = await page.evaluate(() => {
            const viewer = $3Dmol.createViewer(document.getElementById('viewer'));
            let captured;
            try {
                system.embed_style_order(viewer, rendered => {
                    const canvas = document.createElement('canvas');
                    canvas.width = rendered.getCanvas().width;
                    canvas.height = rendered.getCanvas().height;
                    const context = canvas.getContext('2d');
                    context.drawImage(rendered.getCanvas(), 0, 0);
                    const pixels = context.getImageData(0, 0, canvas.width, canvas.height).data;
                    let red = 0, green = 0;
                    for (let i = 0; i < pixels.length; i += 4) {
                        if (pixels[i] > pixels[i + 1] + 30 && pixels[i] > pixels[i + 2] + 30) red++;
                        if (pixels[i + 1] > pixels[i] + 30 && pixels[i + 1] > pixels[i + 2] + 30) green++;
                    }
                    captured = { reused: rendered === viewer, background: Array.from(pixels.slice(0, 4)), red, green };
                });
                return captured;
            } finally { viewer.destroy(); }
        });
        assert.ok(result?.reused, 'The inline embed must complete using the supplied viewer');
        assert.deepEqual(result.background, [255, 255, 255, 255],
            'data-config backgroundColor must reach the reused viewer renderer');
        assert.ok(result.red > 50 && result.green > 50, 'The captured scene must contain the embedded atoms');
        assert.deepEqual(errors, []);
        console.log('Embedded background: reused viewer renders white background and expected atom colors.');
    } finally { await page.close(); }
}

(async () => {
    const browser = await puppeteer.launch({ headless: true, args: ['--no-sandbox', '--enable-unsafe-swiftshader'] });
    try {
        await pythonGrid(browser);
        await viewerLifecycle(browser);
        await reusedViewerBackground(browser);
    } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exitCode = 1; });
