// Run after building: node tests/webgl/opacity-browser.cjs [path/to/3Dmol.js]
// Uses Puppeteer installed with glcheck; fixtures and bundles are entirely local.
const assert = require('node:assert/strict');
const path = require('node:path');
const puppeteer = require('puppeteer');
const bundle = path.resolve(process.argv[2] || path.join(__dirname, '../../build/3Dmol.js'));

(async () => {
    const browser = await puppeteer.launch({ headless: true, args: ['--no-sandbox', '--enable-unsafe-swiftshader'] });
    try {
        for (const preceding of ['none', 'stick', 'mixed-opacity spheres']) {
            const page = await browser.newPage();
            const errors = [];
            page.on('pageerror', error => errors.push(error.message));
            page.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
            await page.setViewport({ width: 640, height: 320 });
            await page.setContent('<style>body{margin:0;display:flex}.viewer{position:relative;width:320px;height:320px}</style>'
                + '<div id="first" class="viewer"></div><div id="second" class="viewer"></div>');
            await page.addScriptTag({ path: bundle });
            const result = await page.evaluate(preceding => {
                const first = $3Dmol.createViewer('first', { backgroundColor: 'white' });
                const gl = first.getRenderer().getContext();
                const snapshot = viewer => {
                    const canvas = document.createElement('canvas');
                    canvas.width = canvas.height = 320;
                    const context = canvas.getContext('2d');
                    context.drawImage(viewer.getCanvas(), 0, 0, 320, 320);
                    const pixels = context.getImageData(0, 0, 320, 320).data;
                    let opaqueBlue = 0, translucentBlue = 0;
                    for (let i = 0; i < pixels.length; i += 4) {
                        if (pixels[i + 2] <= pixels[i] + 20 || pixels[i + 2] <= pixels[i + 1] + 20) continue;
                        if (pixels[i + 3] > 240) opaqueBlue++;
                        if (pixels[i + 3] > 20 && pixels[i + 3] < 180) translucentBlue++;
                    }
                    return { center: Array.from(context.getImageData(160, 160, 1, 1).data), opaqueBlue, translucentBlue };
                };
                let mixedBefore;
                if (preceding === 'stick') {
                    first.addModel('2\nstick\nC -0.8 0 0\nO 0.8 0 0\n', 'xyz');
                    first.setStyle({}, { stick: { color: 'blue', radius: 0.25 } });
                    first.zoomTo(); first.render();
                } else if (preceding === 'mixed-opacity spheres') {
                    // Transparent black avoids color from the clear buffer contributing
                    // to the blended sphere when its bitmap is copied into a 2D canvas.
                    first.setBackgroundColor('black', 0);
                    first.addModel('2\nspheres\nC -3 0 0\nO 3 0 0\n', 'xyz');
                    first.setStyle({}, { sphere: { color: 'blue', radius: 1 } });
                    first.setStyle({ index: 0 }, { sphere: { color: 'blue', radius: 1, opacity: 0.5 } });
                    first.zoomTo(); first.render();
                    mixedBefore = snapshot(first);
                }
                // Separate Renderer instances reuse one OffscreenCanvas/WebGL context,
                // just as consecutive examples do in tests/auto/tests.html.
                const second = $3Dmol.createViewer('second', { backgroundColor: 'white' });
                second.addModel('1\nopaque sphere\nC 0 0 0\n', 'xyz');
                second.setStyle({}, { sphere: { color: 'red' } });
                second.zoomTo(); second.render();
                const initial = snapshot(second);
                // The auto harness renders again immediately before capturing its image.
                second.render();
                const repeated = snapshot(second);
                let mixedAfter;
                if (mixedBefore) { first.render(); mixedAfter = snapshot(first); }
                return {
                    shared: gl === second.getRenderer().getContext(),
                    imposters: second.getRenderer().supportedExtensions().supportsImposters,
                    initial, repeated, mixedBefore, mixedAfter, glError: gl.getError()
                };
            }, preceding);
            assert.equal(result.shared, true, 'Regression requires viewers sharing the WebGL context');
            assert.equal(result.imposters, true, 'Regression requires sphere imposters');
            for (const [render, snapshot] of [['first', result.initial], ['repeated', result.repeated]]) {
                const [red, green, blue, alpha] = snapshot.center;
                assert.ok(alpha > 250 && red > 150 && green < 40 && blue < 40,
                    `${render} opaque sphere render after ${preceding} lost its color/opacity: ${snapshot.center}`);
            }
            if (result.mixedBefore) {
                for (const [stage, snapshot] of [['before', result.mixedBefore], ['after', result.mixedAfter]]) {
                    assert.ok(snapshot.opaqueBlue > 50, `Mixed-opacity scene lost opaque spheres ${stage} switching viewers: ${JSON.stringify(snapshot)}`);
                    assert.ok(snapshot.translucentBlue > 50, `Per-atom sphere translucency was lost ${stage} switching viewers: ${JSON.stringify(snapshot)}`);
                }
            }
            assert.equal(result.glError, 0, `WebGL error after ${preceding}`);
            assert.deepEqual(errors, [], `Browser errors after ${preceding}`);
            await page.close();
        }
        console.log(`Sphere opacity browser regression passed (${path.basename(bundle)}): shared-context sticks, mixed opacity, repeated rendering.`);
    } finally {
        await browser.close();
    }
})().catch(error => { console.error(error); process.exitCode = 1; });
