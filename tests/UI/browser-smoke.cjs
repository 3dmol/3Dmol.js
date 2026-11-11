// Run after npm run build:dev. Uses the Puppeteer installed with glcheck.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const puppeteer = require('puppeteer');

(async () => {
    const browser = await puppeteer.launch({ headless: true, args: ['--no-sandbox', '--enable-unsafe-swiftshader'] });
    try {
        const page = await browser.newPage();
        const errors = [];
        page.on('pageerror', error => errors.push(error.message));
        await page.setViewport({ width: 1100, height: 760 });
        await page.setContent('<style>body{margin:0}#viewer{position:relative;width:100vw;height:100vh}</style><div id="viewer"></div>');
        await page.addScriptTag({ path: path.resolve('build/3Dmol.js') });
        await page.addScriptTag({ path: path.resolve('build/3Dmol.ui.js') });
        const pdb = fs.readFileSync('tests/test_structs/1YCR.pdb', 'utf8');
        await page.evaluate(data => {
            window.viewer = $3Dmol.createViewer('viewer', { backgroundColor: '#edf2f4' });
            const model = viewer.addModel(data, 'pdb');
            model.setStyle({}, { cartoon: { color: 'spectrum' } });
            viewer.zoomTo(); viewer.render();
            window.controls = new $3Dmol.StateManager(viewer, { open: true });
            controls.setModelTitle('1YCR · p53 / MDM2');
        }, pdb);
        assert.equal(await page.locator('.molui-count').map(el => el.textContent).wait(), '818 matching atoms');
        await page.evaluate(() => {
            const model = viewer.getModel();
            const atoms = model.selectedAtoms({});
            const coordinates = [];
            for (let frame = 0; frame < 6; frame++) for (const atom of atoms) {
                coordinates.push(atom.x + frame * 0.2, atom.y, atom.z);
            }
            model.setCoordinates(coordinates, 'array');
            controls.refresh();
        });
        // Exercise the selection card and material colors against actual geometry.
        await page.evaluate(() => {
            const root = controls.ui.element;
            const button = (scope, text) => [...scope.querySelectorAll('button')].find(node => node.textContent === text);
            button(root, 'Add selection').click();
            const add = root.querySelector('.molui-add-representation select');
            add.value = 'stick';
            root.querySelector('.molui-add-representation button').click();
            if ([...controls.selections.values()][0].styles.size !== 2) throw new Error('Additional representation missing');
            for (const style of root.querySelectorAll('.molui-style')) {
                const settings = style.querySelector('details.molui-geometry');
                if (!settings || settings.querySelector('summary').textContent !== 'Settings'
                    || !settings.querySelector('[data-focus$=":opacity"]')) {
                    throw new Error('Every style should contain Settings with opacity');
                }
            }
            if (root.textContent.includes('Show representation') || root.textContent.includes('Enable appearance')) {
                throw new Error('Redundant visibility controls remain');
            }
            if (root.textContent.includes('Clear labels') || root.textContent.includes('Export VRML')) {
                throw new Error('Removed actions remain');
            }
        });
        await page.evaluate(() => {
            const root = controls.ui.element;
            const selection = [...controls.selections.values()][0];
            const originalStyles = [...selection.styles.values()].map(style => style.id);
            const creator = root.querySelector('.molui-selection-create');
            const saved = root.querySelector('.molui-selection-list');
            if (!(saved.compareDocumentPosition(creator) & Node.DOCUMENT_POSITION_FOLLOWING)) {
                throw new Error('New selection form should follow existing selections');
            }
            const edit = () => root.querySelector('.molui-selection-summary-row button').click();
            const action = text => [...root.querySelectorAll('.molui-selection-edit button')].find(node => node.textContent === text).click();
            const field = key => root.querySelector(`.molui-selection-edit [data-field="${key}"]`);
            const chooseTwo = key => {
                const values = [...field(key).querySelectorAll('input[type=checkbox]')].slice(0, 2).map(input => input.dataset.value);
                if (values.length !== 2) throw new Error(`Insufficient fixture values for ${key}`);
                for (const value of values) {
                    const panel = field(key); panel.open = true;
                    [...panel.querySelectorAll('input[type=checkbox]')].find(input => input.dataset.value === value).click();
                }
                return values;
            };
            edit();
            chooseTwo('chain');
            action('Cancel');
            if (Object.keys(selection.spec).length) throw new Error('Cancel changed the saved selection');
            edit();
            for (const key of ['chain', 'resn']) chooseTwo(key);
            const residues = field('resi'); residues.open = true;
            if (residues.querySelector('input[type=checkbox]')) throw new Error('Residues should use range controls');
            const from = residues.querySelector('[data-focus$=":resi:from"]');
            const to = residues.querySelector('[data-focus$=":resi:to"]');
            const start = Number(from.value), end = Math.min(Number(to.value), start + 10);
            for (const [input, value] of [[from, start], [to, end]]) {
                input.value = value;
                input.dispatchEvent(new Event('input', { bubbles: true }));
            }
            field('resi').querySelector('[data-focus$=":resi:add"]').click();
            chooseTwo('elem');
            action('Save selection');
            for (const key of ['chain', 'resn', 'elem']) {
                if (!Array.isArray(selection.spec[key]) || selection.spec[key].length !== 2) {
                    throw new Error(`Multiple ${key} values were not saved`);
                }
            }
            if (selection.spec.resi !== (start === end ? start : `${start}-${end}`)) {
                throw new Error('Residue interval was not saved');
            }
            if (!viewer.selectedAtoms(selection.spec).length || [...selection.styles.keys()].join() !== originalStyles.join()) {
                throw new Error('Editing criteria lost matching atoms or representations');
            }
            const summary = root.querySelector('.molui-selection-summary').textContent;
            if (!summary.includes('Chains') || !summary.includes('Residues') || !summary.includes('Elements')) {
                throw new Error('Saved selection summary is incomplete');
            }
            edit();
            for (const key of ['chain', 'resn', 'resi', 'elem']) field(key).querySelector('[data-focus$=":all"]').click();
            action('Save selection');
            if (Object.keys(selection.spec).length) throw new Error('Clearing the filters did not restore all atoms');
        });
        await page.click('.molui-selection-name', { clickCount: 3 });
        await page.evaluate(() => {
            const input = document.querySelector('.molui-selection-name');
            const r = input.getBoundingClientRect();
            const target = document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2);
            if (target !== input || document.activeElement !== input) throw new Error('Selection input is obscured after scrolling');
        });
        await page.type('.molui-selection-name', 'Protein display');
        // Committing a name on blur must not swallow the following checkbox click.
        await page.click('input[data-focus$=":labels"]');
        if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: '/tmp/3dmol-ui-selection-controls.png' });
        await page.evaluate(() => {
            const entry = [...controls.selections.values()][0];
            if (entry.name !== 'Protein display' || !entry.labelsEnabled || !viewer.labels.length) {
                throw new Error(`Name editing or residue-label toggle failed: ${JSON.stringify({ name: entry.name, enabled: entry.labelsEnabled, labels: viewer.labels.length })}`);
            }
        });
        await page.click('input[data-focus$=":labels"]');
        // Exercise the size editors against actual rendered picking geometry.
        await page.evaluate(() => {
            const root = controls.ui.element;
            const selection = [...controls.selections.values()][0];
            const styles = () => [...selection.styles.values()];
            const styleFor = kind => styles().find(style => style.spec[kind]);
            const inputFor = (kind, name) => root.querySelector(`[data-focus="${selection.id}:${styleFor(kind).id}:${name}"]`);
            const edit = (kind, name, value) => {
                const input = inputFor(kind, name);
                input.closest('details.molui-geometry').open = true;
                input.focus();
                if (input.type === 'checkbox') input.checked = value;
                else input.value = value;
                input.dispatchEvent(new Event('input', { bubbles: true }));
                input.dispatchEvent(new Event('change', { bubbles: true }));
                const updated = inputFor(kind, name);
                if (!updated.closest('details.molui-geometry').open || document.activeElement !== updated) {
                    throw new Error('Shape editor lost expansion or keyboard focus');
                }
            };
            for (const kind of ['sphere', 'cross']) {
                const row = root.querySelector('.molui-add-representation');
                row.querySelector('select').value = kind;
                row.querySelector('button').click();
            }
            viewer.setClickable({ index: 0 }, true, () => {});
            edit('sphere', 'sphere-radius-mode', 'fixed');
            edit('sphere', 'sphere-radius', '1.25');
            edit('cross', 'cross-radius-mode', 'fixed');
            edit('cross', 'cross-radius', '0.75');
            const atom = viewer.selectedAtoms({ index: 0 })[0];
            if (atom.intersectionShape.sphere[0].radius !== 1.25) throw new Error('Sphere radius did not reach renderer');
            // Cross picking lines are intentionally 10% of their visible radius.
            if (Math.abs(atom.intersectionShape.line[0].x - atom.x - 0.075) > 1e-6) {
                throw new Error('Cross radius did not reach renderer');
            }
            const secondaryStructure = viewer.selectedAtoms({}).map(atom => atom.ss);
            for (const shape of ['oval', 'parabola', 'edged', 'default']) {
                edit('cartoon', 'cartoon-style', shape);
                if (atom.style.cartoon.style !== shape) throw new Error('Cartoon shape did not reach renderer');
            }
            edit('cartoon', 'cartoon-arrows', true);
            edit('cartoon', 'cartoon-tubes', true);
            edit('cartoon', 'cartoon-style', 'trace');
            edit('cartoon', 'cartoon-thickness', '0.65');
            if (atom.style.cartoon.thickness !== 0.65 || atom.style.cartoon.arrows || atom.style.cartoon.tubes) {
                throw new Error('Trace size or context-specific flags are wrong');
            }
            if (viewer.selectedAtoms({}).some((atom, index) => atom.ss !== secondaryStructure[index])) {
                throw new Error('Changing cartoon shape corrupted secondary structure');
            }
            edit('cartoon', 'cartoon-style', 'default');
            controls.removeStyle(selection.id, styleFor('sphere').id);
            controls.removeStyle(selection.id, styleFor('cross').id);
        });
        await page.evaluate(async () => {
            if (viewer.labels.length) throw new Error('Residue labels were not removed');
            const root = controls.ui.element;
            const field = (scope, text) => [...scope.querySelectorAll('label.molui-field')]
                .find(label => label.firstElementChild.textContent === text).querySelector('select, input');
            const change = (scope, text, value) => {
                const input = field(scope, text); input.value = value;
                input.dispatchEvent(new Event('change', { bubbles: true }));
            };
            change(root.querySelectorAll('.molui-style')[1], 'Color by', 'scheme');
            change(root.querySelectorAll('.molui-style')[1], 'Color scheme', 'greenCarbon');
            const style = [...[...controls.selections.values()][0].styles.values()][1];
            if (style.spec.stick.colorscheme !== 'greenCarbon') throw new Error('Representation scheme missing');
            const previousStyles = JSON.stringify(viewer.selectedAtoms({}).map(atom => atom.style));
            const creator = root.querySelector('.molui-selection-create');
            change(creator, 'Representation', 'none');
            [...creator.querySelectorAll('button')].find(node => node.textContent === 'Add selection').click();
            const reference = [...controls.selections.values()].pop();
            const card = () => root.querySelector(`[data-selection-id="${reference.id}"].molui-selection-card`);
            if (reference.styles.size || field(card(), 'Representation').value !== 'none') {
                throw new Error('None should create a selection without atom styles');
            }
            card().querySelector('[data-focus$=":labels"]').click();
            if (!viewer.labels.length) throw new Error('None selection could not show residue labels');
            card().querySelector('[data-focus$=":labels"]').click();
            if (viewer.labels.length || JSON.stringify(viewer.selectedAtoms({}).map(atom => atom.style)) !== previousStyles) {
                throw new Error('None selection changed other appearances or left labels behind');
            }
            const scope = controls.addSelection({ index: [0, 1, 2, 3, 4, 5, 6, 7] }, reference.id, 'Surface sample');
            $3Dmol.setSyncSurface(true);
            change(root.querySelector('.molui-surface-create'), 'Show surface for', scope);
            change(root.querySelector('.molui-surface-create'), 'Calculate using', '');
            change(root.querySelector('.molui-surface-create'), 'Color by', 'scheme');
            change(root.querySelector('.molui-surface-create'), 'Color scheme', 'greenCarbon');
            [...root.querySelectorAll('button')].find(node => node.textContent === 'Add surface').click();
        });
        await page.waitForFunction(() => controls.surfaces.size === 1);
        await page.evaluate(() => {
            const entry = [...controls.surfaces.values()][0];
            if (Object.keys(entry.calculationSpec).length || !entry.spec.index) {
                throw new Error('Surface display and calculation selections were not kept independent');
            }
            const actualSurface = viewer.getSurface(entry.surfaceId);
            if (Object.keys(actualSurface.allsel).length || actualSurface.atomsel.index.length !== 8) {
                throw new Error('Calculation atoms did not reach the surface renderer');
            }
            const geometry = viewer.getSurface(entry.surfaceId)[0].geo;
            const colors = geometry.geometryGroups.map(group => Array.from(group.colorArray));
            const card = controls.ui.element.querySelector('.molui-surface-card');
            const scheme = [...card.querySelectorAll('label.molui-field')]
                .find(label => label.firstElementChild.textContent === 'Color scheme').querySelector('select');
            scheme.value = 'orangeCarbon'; scheme.dispatchEvent(new Event('change', { bubbles: true }));
            if (viewer.getSurface(entry.surfaceId)[0].geo !== geometry || entry.style.opacity !== 0.7) {
                throw new Error('Recoloring rebuilt geometry or lost opacity');
            }
            if (!geometry.geometryGroups.some((group, i) => group.colorArray.some((value, j) => value !== colors[i][j]))) {
                throw new Error('Surface vertex colors did not change');
            }
            controls.removeSurface(entry.id);
            const sample = [...controls.selections.values()].find(selection => selection.name === 'Surface sample');
            controls.removeSelection(sample.id);
            const accent = getComputedStyle(controls.ui.element).getPropertyValue('--molui-accent').trim();
            if (accent !== '#000066') throw new Error('Navy palette missing');
        });
        await page.click('.molui-selection-collapse');
        await page.evaluate(() => {
            controls.refresh();
            const card = document.querySelector('.molui-selection-card');
            if (!card.querySelector('.molui-selection-content').hidden
                || card.querySelector('.molui-selection-collapse').getAttribute('aria-expanded') !== 'false'
                || document.activeElement !== card.querySelector('.molui-selection-collapse')) {
                throw new Error('Selection collapse or focus was lost after refresh');
            }
        });
        await page.keyboard.press('Enter');
        await page.evaluate(() => {
            if (document.querySelector('.molui-selection-content').hidden) throw new Error('Keyboard did not expand selection');
        });
        await page.click('button[aria-label="Play trajectory"]');
        await page.waitForFunction(() => viewer.getFrame() > 0);
        await page.click('button[aria-label="Pause trajectory"]');
        await page.click('button[aria-label="Last frame"]');
        await page.waitForFunction(() => viewer.getFrame() === 5);
        await page.click('button[aria-label="First frame"]');
        await page.waitForFunction(() => viewer.getFrame() === 0);
        const initialTransport = await page.evaluate(() => {
            const bar = document.querySelector('.molui-trajectory');
            const play = bar.querySelector('[aria-label="Play trajectory"]');
            if (bar.querySelector('.molui-trajectory-drag')) throw new Error('Playback should not need a grab button');
            const alpha = Number(getComputedStyle(bar).backgroundColor.match(/[\d.]+/g)[3]);
            if (!play.querySelector('svg') || play.textContent.trim()) throw new Error('Playback must use icons');
            if (!(alpha > 0 && alpha < 1)) throw new Error('Playback background should be translucent');
            if (bar.getBoundingClientRect().height > 100) throw new Error('Desktop playback should be compact');
            return { view: viewer.getView(), frame: viewer.getFrame(), left: bar.getBoundingClientRect().left, top: bar.getBoundingClientRect().top };
        });
        const dragStart = await page.locator('.molui-frame-count').map(node => {
            const box = node.getBoundingClientRect();
            return { x: box.x + box.width / 2, y: box.y + box.height / 2 };
        }).wait();
        await page.mouse.move(dragStart.x, dragStart.y);
        await page.mouse.down();
        await page.mouse.move(dragStart.x - 200, dragStart.y - 300, { steps: 8 });
        await page.mouse.up();
        const movedTransport = await page.evaluate(() => {
            const bar = document.querySelector('.molui-trajectory');
            if (!bar.classList.contains('molui-trajectory-moved') || !controls.ui.element.classList.contains('molui-trajectory-floating')) {
                throw new Error('Dragging did not float the controls');
            }
            controls.refresh(); controls.updateUI();
            const box = bar.getBoundingClientRect();
            return { view: viewer.getView(), frame: viewer.getFrame(), left: box.left, top: box.top };
        });
        assert.deepEqual(movedTransport.view, initialTransport.view, 'Moving playback rotated the molecule');
        assert.equal(movedTransport.frame, initialTransport.frame, 'Moving playback changed the frame');
        assert.ok(movedTransport.top < initialTransport.top - 200, 'Playback did not move or lost its position after refresh');
        await page.focus('.molui-trajectory');
        await page.keyboard.press('ArrowRight');
        assert.equal(await page.locator('.molui-trajectory').map(node => node.getBoundingClientRect().left).wait(), movedTransport.left + 10);
        const beforeScrub = await page.locator('.molui-trajectory').map(node => node.getBoundingClientRect().left).wait();
        await page.focus('.molui-frame-slider');
        await page.keyboard.press('ArrowRight');
        await page.waitForFunction(() => viewer.getFrame() === 1);
        assert.equal(await page.locator('.molui-trajectory').map(node => node.getBoundingClientRect().left).wait(), beforeScrub);
        await page.keyboard.press('Home');
        await page.waitForFunction(() => viewer.getFrame() === 0);
        assert.equal(await page.locator('.molui-trajectory').map(node => node.classList.contains('molui-trajectory-moved')).wait(), true);
        await page.focus('.molui-playback-speed');
        await page.keyboard.press('End');
        assert.equal(await page.locator('.molui-playback-speed').map(node => node.value).wait(), '60');
        assert.equal(await page.locator('.molui-trajectory').map(node => node.getBoundingClientRect().left).wait(), beforeScrub);
        await page.select('.molui-playback-speed', '10');
        // Padding also moves the widget; there is no dedicated drag target.
        const background = await page.locator('.molui-trajectory').map(node => {
            const box = node.getBoundingClientRect();
            return { x: box.left + 3, y: box.top + box.height / 2, top: box.top };
        }).wait();
        await page.mouse.move(background.x, background.y);
        await page.mouse.down();
        await page.mouse.move(background.x + 20, background.y - 30, { steps: 4 });
        await page.mouse.up();
        assert.equal(await page.locator('.molui-trajectory').map(node => node.getBoundingClientRect().top).wait(), background.top - 30);
        if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: '/tmp/3dmol-ui-floating-playback.png' });
        await page.setViewport({ width: 360, height: 280 });
        await page.evaluate(() => {
            viewer.resize(); controls.updateUI();
            const bar = document.querySelector('.molui-trajectory');
            const box = bar.getBoundingClientRect();
            const bounds = controls.ui.element.getBoundingClientRect();
            if (getComputedStyle(bar).visibility === 'hidden' || box.left < bounds.left || box.right > bounds.right
                || box.top < bounds.top || box.bottom > bounds.bottom) {
                throw new Error('Moved playback did not remain visible within the resized viewer');
            }
        });
        await page.setViewport({ width: 1100, height: 760 });
        await page.evaluate(() => { viewer.resize(); controls.updateUI(); });
        await page.focus('.molui-trajectory');
        await page.keyboard.press('Home');
        assert.equal(await page.locator('.molui-trajectory').map(node => node.classList.contains('molui-trajectory-moved')).wait(), false);
        if (process.env.UI_SCREENSHOTS) {
            await page.evaluate(() => document.querySelector('.molui-geometry').scrollIntoView({ block: 'center' }));
            await page.screenshot({ path: '/tmp/3dmol-ui-shape-controls.png' });
        }
        // Collapsing playback recovers sidebar space and leaves an accessible launcher.
        await page.setViewport({ width: 360, height: 640 });
        const expandedLayout = await page.evaluate(() => {
            viewer.resize(); controls.updateUI();
            const hamburger = document.querySelector('.molui-toggle').getBoundingClientRect();
            return {
                left: hamburger.left, top: hamburger.top, right: hamburger.right,
                panelBottom: document.querySelector('.molui-panel').getBoundingClientRect().bottom,
                view: viewer.getView(), frame: viewer.getFrame()
            };
        });
        await page.click('button[aria-label="Minimize trajectory controls"]');
        await page.evaluate(expected => {
            controls.refresh(); controls.updateUI();
            const launcher = document.querySelector('.molui-trajectory-toggle');
            const box = launcher.getBoundingClientRect();
            if (!document.querySelector('.molui-trajectory').hidden || launcher.hidden
                || !launcher.querySelector('svg') || launcher.getAttribute('aria-expanded') !== 'false'
                || document.activeElement !== launcher || !controls.ui.isTrajectoryMinimized) {
                throw new Error('Minimized trajectory launcher is missing, unfocused, or inaccessible');
            }
            if (box.top !== expected.top || box.left <= expected.right) throw new Error('Trajectory launcher overlaps the hamburger');
            if (document.querySelector('.molui-panel').getBoundingClientRect().bottom <= expected.panelBottom + 20) {
                throw new Error('Minimizing playback did not release reserved sidebar space');
            }
            if (viewer.getFrame() !== expected.frame || JSON.stringify(viewer.getView()) !== JSON.stringify(expected.view)) {
                throw new Error('Minimizing playback changed the molecule');
            }
        }, expandedLayout);
        if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: '/tmp/3dmol-ui-minimized-playback.png' });
        await page.keyboard.press('Enter');
        await page.evaluate(() => {
            if (document.querySelector('.molui-trajectory').hidden || !document.querySelector('.molui-trajectory-toggle').hidden
                || controls.ui.isTrajectoryMinimized) throw new Error('Keyboard did not restore trajectory controls');
            controls.setControls({ edit: false });
        });
        await page.click('button[aria-label="Minimize trajectory controls"]');
        await page.evaluate(expected => {
            const launcher = document.querySelector('.molui-trajectory-toggle');
            const box = launcher.getBoundingClientRect();
            if (!document.querySelector('.molui-toggle').hidden || document.querySelector('.molui-panel').getBoundingClientRect().height) {
                throw new Error('Trajectory-only mode still exposes the edit controls');
            }
            if (box.left !== expected.left || box.top !== expected.top) throw new Error('Trajectory-only launcher must occupy the hamburger position');
            controls.setControls({ trajectory: false });
            controls.refresh();
            if (!launcher.hidden || !document.querySelector('.molui-trajectory').hidden) throw new Error('Disabled trajectory controls reappeared');
            controls.setControls({ edit: true });
            if (document.querySelector('.molui-toggle').hidden || !launcher.hidden) throw new Error('Edit-only mode exposes trajectory controls');
            controls.setControls({ trajectory: true });
            if (launcher.hidden) throw new Error('Enabling trajectory controls lost the minimized state');
            controls.ui.setTrajectoryMinimized(false);
            controls.ui.setOpen(true);
        }, expandedLayout);
        await page.click('button[aria-label="Minimize trajectory controls"]');
        await page.setViewport({ width: 360, height: 280 });
        await page.evaluate(() => { viewer.resize(); controls.updateUI(); });
        await page.click('.molui-trajectory-toggle');
        await page.evaluate(() => {
            const bar = document.querySelector('.molui-trajectory');
            if (controls.ui.isOpen || bar.hidden || getComputedStyle(bar).visibility === 'hidden'
                || document.activeElement !== bar.querySelector('[aria-label="Minimize trajectory controls"]')) {
                throw new Error('Restoring playback in a tiny viewer must close the sidebar and focus visible controls');
            }
            controls.ui.setOpen(true);
        });
        for (const [width, height] of [[1100, 760], [640, 400], [360, 640], [320, 480], [360, 280]]) {
            await page.setViewport({ width, height });
            await page.evaluate(() => { viewer.resize(); controls.updateUI(); });
            await new Promise(resolve => setTimeout(resolve, 100));
            const layout = await page.evaluate(() => {
                const root = document.querySelector('.molui');
                const panel = document.querySelector('.molui-panel');
                const body = panel.querySelector('.molui-body');
                const header = panel.querySelector('.molui-header');
                const close = header.querySelector('button');
                const video = document.querySelector('.molui-trajectory');
                const rect = element => element.getBoundingClientRect();
                const closeIsReachable = () => {
                    const box = rect(close);
                    return close.contains(document.elementFromPoint(box.x + box.width / 2, box.y + box.height / 2));
                };
                body.scrollTop = body.scrollHeight;
                const closeAtBottom = closeIsReachable();
                body.scrollTop = 0;
                const p = rect(panel), v = rect(video);
                const speed = video.querySelector('.molui-playback-speed');
                const speedStyle = getComputedStyle(speed);
                const text = document.createElement('canvas').getContext('2d');
                text.font = speedStyle.font;
                const speedTextWidth = Math.max(...[...speed.options].map(option => text.measureText(option.text).width));
                const speedBox = rect(speed), loopBox = rect(video.querySelector('.molui-check'));
                return {
                    overflow: root.scrollWidth > root.clientWidth || body.scrollWidth > body.clientWidth,
                    panelWidth: p.width,
                    headerHeight: rect(header).height,
                    overlay: getComputedStyle(header).position === 'absolute' && rect(body).top - p.top <= 2,
                    closeReachable: closeAtBottom && closeIsReachable(),
                    speedFits: speed.clientWidth - parseFloat(speedStyle.paddingLeft) - parseFloat(speedStyle.paddingRight) >= speedTextWidth
                        && speedBox.right <= loopBox.left && speedBox.right <= v.right,
                    videoOverflow: video.scrollWidth > video.clientWidth,
                    videoVisible: getComputedStyle(video).display !== 'none' && getComputedStyle(video).visibility !== 'hidden',
                    overlap: root.classList.contains('molui-compact') && getComputedStyle(video).display !== 'none' && getComputedStyle(video).visibility !== 'hidden' && p.bottom > v.top
                };
            });
            assert.equal(layout.overflow, false, `${width}×${height} overflow`);
            assert.equal(layout.overlap, false, `${width}×${height} playback overlap`);
            assert.equal(layout.speedFits, true, `${width}×${height} FPS text needs room for the dropdown arrow`);
            assert.equal(layout.videoOverflow, false, `${width}×${height} playback overflow`);
            assert.ok(layout.panelWidth > 200);
            assert.ok(layout.headerHeight <= 44, `${width}×${height} header should remain compact`);
            assert.ok(layout.overlay, `${width}×${height} close button should not reserve a header row`);
            assert.ok(layout.closeReachable, `${width}×${height} close button should remain reachable while scrolling`);
            if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: `/tmp/3dmol-ui-${width}x${height}.png` });
            await page.evaluate(() => {
                document.querySelector('.molui-selection-summary-row button').click();
                const list = document.querySelector('.molui-selection-edit [data-field="resi"]');
                list.open = true;
                list.scrollIntoView({ block: 'center' });
                const body = document.querySelector('.molui-body');
                if (body.scrollWidth > body.clientWidth) throw new Error('Expanded selection editor overflows');
            });
            if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: `/tmp/3dmol-ui-selection-${width}x${height}.png` });
            await page.evaluate(() => [...document.querySelectorAll('.molui-selection-edit button')].find(node => node.textContent === 'Cancel').click());
            await page.evaluate(() => {
                const card = document.querySelector('.molui-selection-card');
                card.querySelector('.molui-selection-collapse').click();
                card.scrollIntoView({ block: 'center' });
                for (const selector of ['.molui-selection-name', '.molui-selection-summary']) {
                    if (!card.querySelector(selector).getBoundingClientRect().height) throw new Error('Collapsed selection lost its name or summary');
                }
                if (card.querySelector('.molui-selection-content').getBoundingClientRect().height
                    || card.querySelector('.molui-selection-summary-row button').getBoundingClientRect().height
                    || card.scrollWidth > card.clientWidth) throw new Error('Collapsed selection is not compact');
            });
            if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: `/tmp/3dmol-ui-collapsed-${width}x${height}.png` });
            await page.evaluate(() => document.querySelector('.molui-selection-collapse').click());
        }
        await page.click('button[aria-label="Close controls"]');
        await page.evaluate(() => {
            if (document.activeElement !== document.querySelector('.molui-toggle')) throw new Error('Close did not return focus');
            controls.destroy();
            if (document.querySelector('.molui')) throw new Error('Controls were not removed');
        });
        // Exercise the existing declarative entry point and numbered style layers.
        await page.evaluate(() => {
            const data = document.createElement('script');
            data.id = 'inline-molecule'; data.type = 'chemical/x-xyz';
            data.textContent = '2\nMolecule\nC 0 0 0\nO 1.2 0 0\n';
            document.body.append(data);
            const host = document.getElementById('viewer');
            host.className = 'viewer_3Dmoljs';
            Object.assign(host.dataset, {
                ui: 'true', element: 'inline-molecule', type: 'xyz', style: 'stick:color=blue',
                select1: 'elem:O', style1: 'sphere:color=red'
            });
            $3Dmol.autoload(viewer);
            const atoms = viewer.getModel().selectedAtoms({});
            if (!atoms[0].style.stick || !atoms[1].style.sphere || atoms[1].style.stick) {
                throw new Error('Numbered autoload styles were not applied independently');
            }
            if (viewer.ui.selections.size !== 2) throw new Error('Autoload layers missing');
            viewer.ui.destroy();
        });
        await page.evaluate(() => {
            const host = document.getElementById('viewer');
            const data = document.getElementById('inline-molecule');
            data.textContent = '1\nCarbon\nC 0 0 0\n1\nCarbon\nC 1 0 0\n';
            host.dataset.options = 'multimodel:true;frames:true';
            for (const [flags, edit, trajectory] of [
                [{ ui: 'true' }, true, true],
                [{ uiTrajectory: 'true' }, false, true],
                [{ ui: 'false', uiTrajectory: 'true' }, false, true],
                [{ ui: 'true', uiEdit: 'false' }, false, true],
                [{ ui: 'true', uiTrajectory: 'false' }, true, false],
                [{ uiEdit: 'true' }, true, false],
                [{ ui: 'false' }, false, false]
            ]) {
                viewer.removeAllModels();
                for (const key of ['ui', 'uiEdit', 'uiTrajectory']) delete host.dataset[key];
                Object.assign(host.dataset, flags);
                $3Dmol.autoload(viewer);
                if (viewer.getNumFrames() !== 2) throw new Error('Embedded trajectory fixture did not load frames');
                if (!edit && !trajectory) {
                    if (viewer.ui) throw new Error('Disabled embedded UI was created');
                    continue;
                }
                const root = viewer.ui.ui.element;
                if (root.querySelector('.molui-toggle').hidden === edit
                    || root.querySelector('.molui-trajectory').hidden === trajectory) {
                    throw new Error(`Embedded control flags were not respected: ${JSON.stringify(flags)}`);
                }
                viewer.ui.destroy();
            }
            viewer.removeAllModels();
            delete host.dataset.ui;
            host.dataset.uiTrajectory = 'true';
            data.textContent = '1\nCarbon\nC 0 0 0\n';
            $3Dmol.autoload(viewer);
            const controls = viewer.ui;
            const root = controls.ui.element;
            controls.ui.setTrajectoryMinimized(true);
            if (!root.querySelector('.molui-trajectory').hidden || !root.querySelector('.molui-trajectory-toggle').hidden) {
                throw new Error('A single-frame embedding exposes trajectory controls');
            }
            viewer.addModelsAsFrames('1\nCarbon\nC 0 0 0\n1\nCarbon\nC 1 0 0\n', 'xyz');
            controls.refresh();
            if (root.querySelector('.molui-trajectory-toggle').hidden) throw new Error('Adding frames did not reveal the trajectory launcher');
            viewer.removeAllModels();
            controls.refresh();
            if (!root.querySelector('.molui-trajectory-toggle').hidden) throw new Error('Removing frames left the trajectory launcher visible');
            controls.destroy();
        });
        // Compact desktop spacing must not shrink the interactive touch targets.
        const touch = await browser.newPage();
        touch.on('pageerror', error => errors.push(error.message));
        await touch.setViewport({ width: 360, height: 640, isMobile: true, hasTouch: true });
        await touch.setContent('<meta name="viewport" content="width=device-width,initial-scale=1"><style>body{margin:0}#viewer{position:relative;width:100vw;height:100vh}</style><div id="viewer"></div>');
        await touch.addScriptTag({ path: path.resolve('build/3Dmol.js') });
        await touch.addScriptTag({ path: path.resolve('build/3Dmol.ui.js') });
        await touch.evaluate(() => {
            const viewer = window.viewer = $3Dmol.createViewer('viewer');
            viewer.addModelsAsFrames('1\nCarbon\nC 0 0 0\n1\nCarbon\nC 1 0 0\n', 'xyz');
            const controls = window.controls = new $3Dmol.StateManager(viewer, { open: true });
            controls.createSelectionAndStyle({}, { sphere: {} });
            if (!matchMedia('(pointer: coarse)').matches) throw new Error('Touch fixture is not using coarse pointers');
            for (const selector of ['.molui-header button', '.molui-toggle', '.molui-style-header button', '.molui-selection-name', '.molui-selection-visibility', '.molui-selection-collapse', '.molui-trajectory-transport button']) {
                if (document.querySelector(selector).getBoundingClientRect().height < 40) {
                    throw new Error(`Touch target is too small: ${selector}`);
                }
            }
        });
        const touchInitial = await touch.evaluate(() => {
            const counter = document.querySelector('.molui-frame-count').getBoundingClientRect();
            const bar = document.querySelector('.molui-trajectory').getBoundingClientRect();
            if (bar.height > 140) throw new Error('Touch playback should stay compact');
            return { x: counter.x + counter.width / 2, y: counter.y + counter.height / 2, top: bar.top, view: viewer.getView(), frame: viewer.getFrame() };
        });
        const touchSession = await touch.createCDPSession();
        await touchSession.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x: touchInitial.x, y: touchInitial.y }] });
        await touchSession.send('Input.dispatchTouchEvent', { type: 'touchMove', touchPoints: [{ x: touchInitial.x, y: touchInitial.y - 120 }] });
        await touchSession.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
        const touchMoved = await touch.evaluate(() => {
            const bar = document.querySelector('.molui-trajectory');
            if (!bar.classList.contains('molui-trajectory-moved')) throw new Error('Touch drag did not move playback');
            return { top: bar.getBoundingClientRect().top, view: viewer.getView(), frame: viewer.getFrame() };
        });
        assert.ok(touchMoved.top < touchInitial.top - 100, 'Touch playback did not move');
        assert.deepEqual(touchMoved.view, touchInitial.view, 'Touch drag rotated the molecule');
        assert.equal(touchMoved.frame, touchInitial.frame, 'Touch drag changed frame');
        const touchSlider = await touch.locator('.molui-frame-slider').map(node => {
            const box = node.getBoundingClientRect();
            return { left: box.left + 8, right: box.right - 8, y: box.top + box.height / 2 };
        }).wait();
        await touchSession.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x: touchSlider.left, y: touchSlider.y }] });
        await touchSession.send('Input.dispatchTouchEvent', { type: 'touchMove', touchPoints: [{ x: touchSlider.right, y: touchSlider.y }] });
        await touchSession.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
        await touch.waitForFunction(() => viewer.getFrame() === 1);
        assert.equal(await touch.locator('.molui-trajectory').map(node => node.getBoundingClientRect().top).wait(), touchMoved.top, 'Touch scrubbing moved the widget');
        if (process.env.UI_SCREENSHOTS) await touch.screenshot({ path: '/tmp/3dmol-ui-touch.png' });
        await touch.tap('button[aria-label="Minimize trajectory controls"]');
        await touch.evaluate(() => {
            const launcher = document.querySelector('.molui-trajectory-toggle');
            const box = launcher.getBoundingClientRect();
            const hamburger = document.querySelector('.molui-toggle').getBoundingClientRect();
            if (box.width < 40 || box.height < 40 || box.left <= hamburger.right || box.top !== hamburger.top) {
                throw new Error('Touch trajectory launcher is too small or overlaps the hamburger');
            }
            controls.setControls({ edit: false });
            if (launcher.getBoundingClientRect().left !== hamburger.left) throw new Error('Touch trajectory-only launcher is misplaced');
        });
        await touch.tap('.molui-trajectory-toggle');
        assert.equal(await touch.locator('.molui-trajectory').map(node => node.getBoundingClientRect().top).wait(), touchMoved.top, 'Restoring playback lost its moved position');
        await touch.close();
        assert.deepEqual(errors, []);
        console.log('WebGL UI smoke passed: selection editing, labels, style settings, surfaces, icon playback, minimizing, independent controls, mouse/touch dragging, keyboard placement, resizing, five viewport sizes, cleanup, autoload.');
    } finally {
        await browser.close();
    }
})().catch(error => { console.error(error); process.exitCode = 1; });
