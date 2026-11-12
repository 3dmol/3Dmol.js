const assert = require('node:assert/strict');

// Exercise the shipped tutorial assets without replacing its library or protein
// data. The caller may intercept unrelated homepage examples and analytics.
async function checkTutorial(page, baseUrl, label = 'Tutorial') {
    const requests = [];
    const record = request => requests.push(request.url());
    page.on('request', record);
    try {
        await page.goto(new URL('tutorials/code/index.html', baseUrl).href, { waitUntil: 'networkidle0' });
        assert.equal(await page.$('#btn-01-alt'), null, `${label}: background-color button remains`);
        assert.equal(await page.$$eval('.tutorial-try-button', nodes => nodes.length), 3,
            `${label}: expected the sphere and two protein buttons`);
        assert.equal(await page.evaluate(() => typeof window.jQuery), 'undefined',
            `${label}: live examples should not need a CDN jQuery dependency`);

        const appearance = await page.$$eval('.tutorial-try-button', nodes => nodes.map(node => {
            const style = getComputedStyle(node);
            return {
                id: node.id, background: style.backgroundColor, border: parseFloat(style.borderTopWidth),
                borderStyle: style.borderTopStyle, paddingX: parseFloat(style.paddingLeft),
                paddingY: parseFloat(style.paddingTop), cursor: style.cursor
            };
        }));
        for (const button of appearance) {
            assert.ok(button.background !== 'transparent' && !/rgba\([^)]*,\s*0\)$/.test(button.background),
                `${label}: ${button.id} has no visible button fill`);
            assert.ok(button.border >= 1 && button.borderStyle !== 'none', `${label}: ${button.id} has no button border`);
            assert.ok(button.paddingX >= 8 && button.paddingY >= 4, `${label}: ${button.id} has no button padding`);
            assert.equal(button.cursor, 'pointer');
        }

        await page.keyboard.press('Tab');
        await page.focus('#btn-01');
        const focus = await page.$eval('#btn-01', node => {
            const style = getComputedStyle(node);
            return { visible: node.matches(':focus-visible'), outline: parseFloat(style.outlineWidth), shadow: style.boxShadow };
        });
        assert.ok(focus.visible && (focus.outline >= 2 || focus.shadow !== 'none'), `${label}: keyboard focus is not visible`);

        // Capture real viewers when a demo renders; do not replace any drawing
        // or loading code. All clicks below still run the shipped handlers.
        await page.evaluate(() => {
            window.__documentationViewers = Object.create(null);
            window.__documentationRenders = Object.create(null);
            const original = $3Dmol.GLViewer.prototype.render;
            $3Dmol.GLViewer.prototype.render = function (...args) {
                const id = this.container?.id;
                if (id) {
                    window.__documentationViewers[id] = this;
                    window.__documentationRenders[id] = (window.__documentationRenders[id] || 0) + 1;
                }
                return original.apply(this, args);
            };
        });

        async function pixels(id) {
            return page.evaluate(async id => {
                const viewer = window.__documentationViewers[id];
                viewer.render();
                const image = new Image();
                image.src = viewer.pngURI();
                await image.decode();
                const canvas = document.createElement('canvas');
                canvas.width = image.width;
                canvas.height = image.height;
                const context = canvas.getContext('2d');
                context.drawImage(image, 0, 0);
                const rgba = context.getImageData(0, 0, canvas.width, canvas.height).data;
                const background = Array.from(rgba.slice(0, 3));
                let foreground = 0;
                for (let i = 0; i < rgba.length; i += 4) {
                    if (rgba[i + 3] > 0 && background.some((channel, index) => Math.abs(rgba[i + index] - channel) > 8)) foreground++;
                }
                return { background, foreground };
            }, id);
        }

        for (const [button, container] of [['btn-01', 'container-01'], ['btn-02', 'container-02'], ['btn-03', 'container-03']]) {
            if (button === 'btn-01') await page.keyboard.press('Enter');
            else await page.click('#' + button);
            await page.waitForFunction(id => !!window.__documentationViewers[id], {}, container);
            const rendered = await pixels(container);
            assert.deepEqual(rendered.background, button === 'btn-01' ? [255, 165, 0] : [255, 255, 255],
                `${label}: ${button} rendered the wrong background`);
            assert.ok(rendered.foreground > 200, `${label}: ${button} did not draw a visible molecule or sphere`);
            const before = await page.evaluate(id => ({
                atoms: window.__documentationViewers[id].selectedAtoms({}).length,
                renders: window.__documentationRenders[id]
            }), container);
            if (button !== 'btn-01') {
                assert.equal(before.atoms, 818, `${label}: ${button} did not load the bundled 1YCR protein`);
                assert.ok(await page.evaluate(id => window.__documentationViewers[id].selectedAtoms({})
                    .every(atom => !!atom.style.cartoon), container), `${label}: ${button} lost its cartoon style`);
            }
            if (button === 'btn-03') {
                assert.ok(await page.evaluate(() => window.__documentationViewers['container-03'].selectedAtoms({})
                    .every(atom => atom.style.cartoon?.colorfunc?.(atom) === (atom.resi % 2 === 0 ? 'white' : 'green'))),
                `${label}: dynamic coloring did not run`);
            }
            await page.click('#' + button);
            await page.waitForFunction((id, count) => window.__documentationRenders[id] > count, {}, container, before.renders);
            assert.equal(await page.evaluate(id => window.__documentationViewers[id].selectedAtoms({}).length, container), before.atoms,
                `${label}: repeated clicks duplicated the model`);
            assert.equal(await page.$$eval('#' + container + ' canvas', nodes => nodes.length), 1,
                `${label}: repeated clicks created extra viewers`);
        }

        const viewport = page.viewport();
        await page.setViewport({ width: 390, height: 844 });
        // Viewer canvases resize on the browser's resize event, after CSS
        // media queries have already changed the surrounding layout.
        await page.waitForFunction(() => document.documentElement.scrollWidth <= innerWidth + 1);
        await page.setViewport(viewport);

        for (const request of requests) {
            const url = new URL(request);
            assert.ok(request.startsWith(baseUrl) || ['data:', 'blob:'].includes(url.protocol)
                || (url.pathname === '/favicon.ico' && url.origin === new URL(baseUrl).origin)
                || ['fonts.googleapis.com', 'fonts.gstatic.com'].includes(url.hostname),
            `${label}: tutorial depends on a resource outside the generated docs: ${request}`);
        }
        console.log(`${label}: three styled/focusable buttons, visible sphere, both 818-atom proteins and repeat clicks passed.`);
    } finally {
        page.off('request', record);
    }
}

module.exports = { checkTutorial };
