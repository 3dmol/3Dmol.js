// Run after npm run build:dev and npm run build:prod. All molecular data and assets are served locally;
// external structure requests are intercepted so this checks the URL API offline.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const puppeteer = require('puppeteer');

const root = path.resolve(__dirname, '../..');
const pdb = fs.readFileSync(path.join(root, 'tests/test_structs/1YCR.pdb'), 'utf8');
const sdf = fs.readFileSync(path.join(root, 'tests/test_structs/benzene.sdf'), 'utf8');
const cif = fs.readFileSync(path.join(root, 'tests/test_structs/1ubq.cif'), 'utf8');
const xyz = '2\nfirst\nC 0 0 0\nO 1 0 0\n2\nsecond\nC 10 0 0\nO 11 0 0\n';
// Negative coordinates avoid the core parser's existing ambiguous SDF heuristic.
const extensionlessXYZ = '2\nfirst\nC -1 -1 -1\nO -2 -1 -1\n2\nsecond\nC -10 -1 -1\nO -11 -1 -1\n';
const tutorial = fs.readFileSync(path.join(root, 'tutorials/url.md'), 'utf8');
const examples = [...new Set([...tutorial.matchAll(/(?:href=["']|\]\()([^"'\s)]*viewer\.html\?[^"'\s)]*)/g)].map(match => match[1].replace(/&amp;/g, '&')))];
const sessionSocket = `window.__sessionEvents = []; window.__sessionHandlers = {};
window.io = { connect: function () { return {
    on: function (name, callback) { window.__sessionHandlers[name] = callback; return this; },
    emit: function (name, data) { window.__sessionEvents.push({name: name, data: data}); return this; },
    send: function () {}, disconnect: function () {}
}; } };`;

// Observe the real renderer before viewer.js starts, including renders issued
// inside UI appearance updates, surface generation, and trajectory restoration.
const observeDrawing = `
window.__viewerDraws = [];
for (const method of ['render', 'zoomTo']) {
    const original = $3Dmol.GLViewer.prototype[method];
    $3Dmol.GLViewer.prototype[method] = function (...args) {
        const record = () => {
            if (this.selectedAtoms({}).length) window.__viewerDraws.push({
                method, camera: this.getView(), slab: this.getSlab(), frame: this.getFrame()
            });
        };
        if (method === 'render') record();
        const result = original.apply(this, args);
        if (method === 'zoomTo') record();
        return result;
    };
}
`;

const mime = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json', '.png': 'image/png', '.svg': 'image/svg+xml' };
const proxyRequests = [];
const server = http.createServer((request, response) => {
    const address = new URL(request.url, 'http://localhost');
    const pathname = decodeURIComponent(address.pathname);
    if (pathname === '/echo.cgi') {
        const remote = address.searchParams.get('url');
        proxyRequests.push(remote);
        if (/\/9ZZZ\.(?:pdb|cif)$/i.test(remote || '') || /\/9ZZY\.cif$/i.test(remote || '')) {
            response.writeHead(200, { 'Content-Type': 'text/plain' }).end(/\.cif$/i.test(remote) ? cif : pdb);
        } else if (/\/9ZZY\.pdb$/i.test(remote || '')) {
            // echo.cgi relays upstream error bodies without their HTTP status.
            response.writeHead(200, { 'Content-Type': 'text/html' }).end('<html><body>404: PDB file not found</body></html>');
        } else response.writeHead(404).end('Unexpected proxy request');
        return;
    }
    if (pathname === '/structure-without-extension' || pathname === '/fixture.mcif') {
        response.writeHead(200, { 'Content-Type': 'text/plain' }).end(pathname.endsWith('.mcif') ? cif : pdb);
        return;
    }
    if (pathname === '/literal+structure.xyz' || pathname === '/trajectory-without-extension') {
        response.writeHead(200, { 'Content-Type': 'text/plain' }).end(pathname.endsWith('without-extension') ? extensionlessXYZ : xyz);
        return;
    }
    if (pathname === '/slow.xyz') {
        setTimeout(() => { if (!response.destroyed) response.writeHead(200, { 'Content-Type': 'text/plain' }).end('1\nslow\nC 0 0 0\n'); }, 250);
        return;
    }
    if (pathname === '/fast.xyz') {
        response.writeHead(200, { 'Content-Type': 'text/plain' }).end('2\nfast\nC 0 0 0\nO 1 0 0\n');
        return;
    }
    if (pathname === '/frames.xyz') {
        response.writeHead(200, { 'Content-Type': 'text/plain' }).end(xyz);
        return;
    }
    const file = path.resolve(root, '.' + pathname);
    if (!file.startsWith(root + path.sep)) { response.writeHead(403).end(); return; }
    try {
        let body = fs.readFileSync(file);
        // Deterministic surface timing; the geometry code is the production code.
        if (/\/build\/3Dmol(?:-min)?\.js$/.test(pathname)) body = Buffer.concat([body, Buffer.from('\n;$3Dmol.setSyncSurface(true);\n' + observeDrawing)]);
        response.writeHead(200, { 'Content-Type': mime[path.extname(file)] || 'application/octet-stream' }).end(body);
    } catch { response.writeHead(404).end('Not found'); }
});

(async () => {
    await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
    const origin = `http://127.0.0.1:${server.address().port}`;
    let browser;
    try {
        browser = await puppeteer.launch({ headless: true, args: ['--no-sandbox', '--enable-unsafe-swiftshader'] });
        const page = await browser.newPage();
        page.setDefaultNavigationTimeout(60000);
        const errors = [], unexpectedRequests = [], fallbackRequests = [];
        page.on('pageerror', error => errors.push(error.message));
        await page.setViewport({ width: 1100, height: 760 });
        await page.setRequestInterception(true);
        page.on('request', request => {
            const url = new URL(request.url());
            if (url.origin === origin || ['data:', 'blob:'].includes(url.protocol)) return void request.continue();
            if (/jquery(?:-[\d.]+)?(?:\.min)?\.js$/.test(url.pathname)) {
                return void request.respond({ status: 200, contentType: 'text/javascript', body: fs.readFileSync(path.join(root, 'node_modules/jquery/dist/jquery.min.js')) });
            }
            if (/socket\.io(?:\.min)?\.js$/.test(url.pathname)) {
                return void request.respond({ status: 200, contentType: 'text/javascript', body: sessionSocket });
            }
            if (url.hostname === 'files.rcsb.org') {
                if (/\/9ZZ[ZY]\.(?:pdb|cif)$/i.test(url.pathname)) {
                    fallbackRequests.push(request.url());
                    return void request.respond({ status: 404, contentType: 'text/plain', headers: { 'Access-Control-Allow-Origin': '*' }, body: 'Fixture requires the proxy' });
                }
                return void request.respond({ status: 200, contentType: 'text/plain', headers: { 'Access-Control-Allow-Origin': '*' }, body: /\.cif$/.test(url.pathname) ? cif : pdb });
            }
            if (url.hostname === 'pubchem.ncbi.nlm.nih.gov') {
                return void request.respond({ status: 200, contentType: 'chemical/x-mdl-sdfile', headers: { 'Access-Control-Allow-Origin': '*' }, body: sdf });
            }
            if (url.hostname.toLowerCase() === '3dmol.org' && url.pathname === '/tests/test_structs/benzene.sdf') {
                return void request.respond({ status: 200, contentType: 'chemical/x-mdl-sdfile', headers: { 'Access-Control-Allow-Origin': '*' }, body: sdf });
            }
            unexpectedRequests.push(request.url());
            void request.abort();
        });

        const ready = async () => {
            await page.waitForFunction(() => !!window.viewerApp?.ready);
            await page.evaluate(() => viewerApp.ready);
            await page.waitForFunction(() => !!window.viewerApp?.state?.ui?.element?.isConnected);
            await page.evaluate(() => {
                const form = viewerApp.state.ui.element.querySelector('form');
                const sources = [...form.querySelector('select').options].map(option => option.value);
                if (sources.join(',') !== 'pdb,cid,url' || form.querySelector('input[type=file]')) {
                    throw new Error('Standalone viewer should offer only PDB, PubChem, and URL loading');
                }
            });
        };
        const navigate = async suffix => {
            try { await page.goto(origin + '/viewer.html' + suffix, { waitUntil: 'domcontentloaded' }); }
            catch (error) { error.message = `${suffix}: ${error.message}`; throw error; }
            await ready();
        };
        const waitURL = async previous => {
            await page.waitForFunction(previous => location.href !== previous && location.hash.startsWith('#state='), {}, previous);
            return page.url();
        };
        const assertFittedFirstRender = async (frame = 0) => {
            const draws = await page.evaluate(() => window.__viewerDraws);
            const fit = draws.find(entry => entry.method === 'zoomTo');
            const first = draws.find(entry => entry.method === 'render');
            assert.ok(fit && first, 'Load did not fit and render a molecule');
            assert.ok(draws.indexOf(fit) < draws.indexOf(first), 'Molecule rendered before its camera was fitted');
            assert.equal(fit.frame, frame, 'Camera was fitted to a different trajectory frame');
            assert.equal(first.frame, frame, 'First molecule render used a different trajectory frame');
            assert.deepEqual(first.camera, fit.camera, 'First molecule render did not use the fitted camera');
            assert.deepEqual(first.slab, fit.slab, 'First molecule render did not use the fitted clipping planes');
            return first;
        };

        // Every actual URL linked by the tutorial is loaded through viewer.html.
        for (const example of examples) {
            const query = new URL(example, origin + '/viewer.html').search;
            await navigate(query);
            await assertFittedFirstRender();
            await page.evaluate(({ data, format, query }) => {
                const { viewer, state } = viewerApp;
                if (document.querySelectorAll('.molui').length !== 1 || document.querySelector('#selection_list')) {
                    throw new Error('Standalone viewer did not replace its old editor with the UI');
                }
                const expected = new $3Dmol.GLModel(0);
                expected.addMolData(data, format); expected.setStyle({}, { line: {} });
                let selection = {};
                const surfaces = [], labels = [];
                for (const [command, value] of new URLSearchParams(query)) {
                    const spec = $3Dmol.specStringToObject(value);
                    if (command === 'select') selection = spec;
                    else if (command === 'style') expected.setStyle(selection, spec);
                    else if (command === 'addstyle') expected.setStyle(selection, spec, true);
                    else if (command === 'surface') surfaces.push({ spec: selection, style: spec });
                    else if (command === 'labelres') labels.push({ spec: selection, style: spec });
                }
                const atoms = viewer.selectedAtoms({});
                if (!atoms.length || atoms.length !== expected.selectedAtoms({}).length) throw new Error('Documented source did not load');
                const actualStyles = atoms.map(atom => atom.style);
                const expectedStyles = expected.selectedAtoms({}).map(atom => atom.style);
                if (JSON.stringify(actualStyles) !== JSON.stringify(expectedStyles)) throw new Error('Documented style order changed');
                if (state.surfaces.size !== surfaces.length) throw new Error('Documented surfaces missing');
                for (const [index, entry] of [...state.surfaces.values()].entries()) {
                    if (JSON.stringify(entry.spec) !== JSON.stringify(surfaces[index].spec)
                        || JSON.stringify(entry.calculationSpec) !== JSON.stringify(surfaces[index].spec)
                        || entry.type !== 1) throw new Error('Documented surface scope/type changed');
                }
                for (const item of labels) {
                    const entry = [...state.selections.values()].find(entry => entry.labelsEnabled && JSON.stringify(entry.spec) === JSON.stringify(item.spec));
                    if (!entry || JSON.stringify(entry.labelStyle) !== JSON.stringify(item.style)) throw new Error('Residue label settings were lost');
                    if (!viewer.labels.length) throw new Error('Documented labels were not rendered');
                    for (const [key, value] of Object.entries(item.style)) {
                        if (!viewer.labels.some(label => label.stylespec[key] === value)) throw new Error(`Label setting ${key} did not reach renderer`);
                    }
                }
            }, { data: query.includes('url=') ? sdf : pdb, format: query.includes('url=') ? 'sdf' : 'pdb', query });
            if (process.env.UI_SCREENSHOTS && example === examples[0]) await page.screenshot({ path: '/tmp/3dmol-viewer-url-example.png' });
        }

        await navigate('?url=tests/test_structs/benzene.sdf.gz&style=stick');
        assert.equal(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length), 12, 'Compressed relative URL did not load');
        await navigate('?cid=241&style=stick');
        assert.equal(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length), 12, 'CID URL did not load');
        await navigate('?cif=1UBQ&style=cartoon');
        assert.ok(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length > 100), 'CIF URL did not load');
        await navigate('?mmtf=1YCR&style=cartoon');
        assert.ok(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length > 100), 'Legacy MMTF identifier alias did not load');
        await navigate('?url=tests/test_structs/benzene-homo.cube&style=stick');
        assert.equal(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length), 12, 'Supported CUBE structure URL did not load');
        await page.evaluate(() => viewerApp.loadURL(viewerApp.syncURL()));
        assert.equal(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length), 12, 'CUBE source did not survive share/restore');
        await navigate('?url=tests/auto/data/mg.vasp&style=sphere');
        assert.deepEqual(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).map(atom => atom.elem).sort()), ['Mg', 'O'], 'Supported VASP structure URL did not load');
        await navigate('?url=fixture.mcif&style=cartoon');
        assert.ok(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length > 100), 'Supported MCIF alias did not load');
        await navigate('?url=structure-without-extension&style=cartoon');
        const detectedProtein = await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length);
        assert.ok(detectedProtein > 100, 'Extensionless PDB response was not detected from its contents');
        await page.evaluate(() => viewerApp.loadURL(viewerApp.syncURL()));
        assert.equal(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length), detectedProtein, 'Detected extensionless source did not survive share/restore');
        await navigate('?url=trajectory-without-extension&style=stick');
        assert.equal(await page.evaluate(() => viewerApp.viewer.getNumFrames()), 2, 'Extensionless XYZ trajectory was not detected from its contents');
        await navigate('?url=literal+structure.xyz&style=stick');
        assert.equal(await page.evaluate(() => viewerApp.getSnapshot().models[0].value), 'literal+structure.xyz', 'A literal plus in a legacy source URL was changed to a space');
        assert.equal(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length), 2, 'Literal-plus source path did not load');
        await navigate('?pdb=9ZZZ&style=stick');
        assert.ok(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length > 100), 'PDB proxy fallback did not load');
        assert.ok(fallbackRequests.some(url => url.endsWith('/9ZZZ.pdb')), 'PDB fallback fixture did not try the remote PDB source');
        assert.ok(proxyRequests.some(url => /\/9ZZZ\.(?:pdb|cif)$/i.test(url || '')), 'PDB source failure did not use the same-origin proxy');
        await navigate('?pdb=9ZZY&style=cartoon');
        assert.ok(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length > 100), 'CIF proxy fallback did not load after the PDB proxy failed');
        assert.deepEqual(proxyRequests.filter(url => /\/9ZZY\./.test(url || '')).map(url => new URL(url).pathname),
            ['/view/9ZZY.pdb', '/view/9ZZY.cif'], 'PDB fallback did not try both PDB and CIF through the proxy');
        await navigate('?pdb=1YCR&select=chain:B&style=%7B%7D');
        await page.evaluate(() => {
            if (viewerApp.viewer.selectedAtoms({ chain: 'B' }).some(atom => Object.keys(atom.style).length)) throw new Error('Explicit empty style did not clear atoms');
            if (viewerApp.viewer.selectedAtoms({ chain: 'A' }).some(atom => !atom.style.line)) throw new Error('Explicit empty style changed other atoms');
        });
        await navigate('?pdb=1YCR&select=chain:A&colorbyelement=colorscheme:greenCarbon&select=chain:B&colorbyproperty=prop:b;scheme:RWB');
        await page.evaluate(() => {
            const viewer = viewerApp.viewer;
            const carbon = viewer.selectedAtoms({ chain: 'A', elem: 'C' });
            if (!carbon.length || carbon.some(atom => atom.color !== $3Dmol.elementColors.greenCarbon.C)) throw new Error('Element color command did not reach atoms');
            const atoms = viewer.selectedAtoms({ chain: 'B' });
            const range = [Math.min(...atoms.map(atom => atom.b)), Math.max(...atoms.map(atom => atom.b))];
            const gradient = new $3Dmol.Gradient.RWB();
            if (atoms.some(atom => atom.color !== gradient.valueToHex(atom.b, range))) throw new Error('Property color command did not reach atoms');
        });
        const coloredAtoms = await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).map(atom => ({ color: atom.color, style: atom.style })));
        await page.evaluate(() => viewerApp.loadURL(viewerApp.syncURL()));
        assert.deepEqual(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).map(atom => ({ color: atom.color, style: atom.style }))), coloredAtoms,
            'Legacy element/property colors changed when sharing and restoring');

        // Default surface coloring is a snapshot of atom colors at this command,
        // so a later color command must not retroactively recolor that surface.
        const surfaceQuery = '?pdb=1YCR&select=index:0,1,2,3&surface=opacity:0.7';
        const surfaceColors = () => page.evaluate(() => {
            const surface = [...viewerApp.state.surfaces.values()][0];
            return viewerApp.viewer.getSurface(surface.surfaceId)[0].geo.geometryGroups.flatMap(group => Array.from(group.colorArray));
        });
        await navigate(surfaceQuery);
        const originalSurfaceColors = await surfaceColors();
        await navigate(surfaceQuery + '&colorbyelement=colorscheme:greenCarbon');
        assert.deepEqual(await surfaceColors(), originalSurfaceColors, 'Later atom coloring changed an earlier surface');
        await page.evaluate(() => viewerApp.loadURL(viewerApp.syncURL()));
        assert.deepEqual(await surfaceColors(), originalSurfaceColors, 'Surface color snapshot changed after URL restoration');
        await navigate('?pdb=1YCR&select=index:0,1,2,3&colorbyelement=colorscheme:greenCarbon&surface=opacity:0.7');
        assert.notDeepEqual(await surfaceColors(), originalSurfaceColors, 'Surface fixture did not exercise atom color order');

        // Reject malformed saved settings before changing a working scene, and
        // retain usable controls if a downloaded/embedded model contains no atoms.
        await page.evaluate(async () => {
            const count = viewerApp.viewer.selectedAtoms({}).length;
            try { await viewerApp.loadURL('#state=%7Bbroken'); throw new Error('Malformed state was accepted'); }
            catch (error) { if (error.message === 'Malformed state was accepted') throw error; }
            if (viewerApp.viewer.selectedAtoms({}).length !== count || !viewerApp.state.ui.element.isConnected) throw new Error('Malformed settings destroyed the current scene');
            const bad = { models: [{ kind: 'file', value: 'bad.xyz', format: 'xyz', frames: false, data: { encoding: 'text', value: 'not a molecule' } }], selections: [], surfaces: [], view: {} };
            try { await viewerApp.loadURL(ViewerURL.serialize(bad, location.href)); throw new Error('Empty model was accepted'); }
            catch (error) { if (error.message === 'Empty model was accepted') throw error; }
            if (!viewerApp.state.ui.element.isConnected || !viewerApp.state.ui.element.querySelector('form')) throw new Error('Model load error left unusable controls');
        });

        // Rapid URL navigation must not allow a late old request to replace the
        // newer scene. The real popstate path performs both asynchronous loads.
        const slowRequest = page.waitForRequest(request => new URL(request.url()).pathname === '/slow.xyz');
        await page.evaluate(() => {
            history.pushState(null, '', '?url=slow.xyz&style=stick');
            dispatchEvent(new PopStateEvent('popstate'));
            window.slowRestore = viewerApp.ready;
        });
        await slowRequest;
        await page.evaluate(async () => {
            history.pushState(null, '', '?url=fast.xyz&style=sphere');
            dispatchEvent(new PopStateEvent('popstate'));
            await viewerApp.ready;
            await window.slowRestore;
        });
        assert.equal(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).length), 2, 'Stale URL load replaced the latest model');
        assert.equal(await page.evaluate(() => viewerApp.getSnapshot().models[0].value), 'fast.xyz');
        assert.equal(await page.evaluate(() => document.querySelectorAll('.molui').length), 1);

        const colorsByFrame = () => page.evaluate(async () => {
            const viewer = viewerApp.viewer, originalFrame = viewer.getFrame(), colors = [];
            for (let frame = 0; frame < viewer.getNumFrames(); frame++) {
                await viewer.setFrame(frame);
                colors.push(viewer.selectedAtoms({}).map(atom => atom.color));
            }
            await viewer.setFrame(originalFrame); viewer.render();
            return colors;
        });
        await navigate('?url=frames.xyz&style=sphere&frame=-3');
        const firstTrajectoryRender = await assertFittedFirstRender(1);
        assert.ok(Math.abs(firstTrajectoryRender.camera[0] + 10.5) < 1e-6,
            'The initial camera was centered on frame zero instead of the requested frame');
        await page.evaluate(async () => {
            const { viewer, state } = viewerApp;
            if (viewer.getFrame() !== 1 || viewer.selectedAtoms({})[0].x !== 10) throw new Error('Negative URL frame did not wrap to the final trajectory frame');
            await viewer.setFrame(-2); viewer.render();
            if (viewer.getFrame() !== 0 || viewer.selectedAtoms({})[0].x !== 0) throw new Error('Direct negative frame did not wrap to the second-to-last frame');
            await state.ui.setTrajectorySettings({ frame: 3 });
            if (viewer.getFrame() !== 1 || viewer.selectedAtoms({})[0].x !== 10) throw new Error('Restored positive frame did not wrap modulo the trajectory length');
            const url = viewerApp.syncURL();
            await viewerApp.loadURL(url);
            if (viewerApp.viewer.getFrame() !== 1) throw new Error('Wrapped frame did not round-trip through the URL');
        });
        await navigate('?url=frames.xyz&style=sphere&colorbyelement=colorscheme:greenCarbon');
        const legacyFrameColors = await colorsByFrame();
        assert.notEqual(legacyFrameColors[0][0], legacyFrameColors[1][0], 'Legacy color fixture did not have distinct frame colors');
        await page.evaluate(() => viewerApp.loadURL(viewerApp.syncURL()));
        assert.deepEqual(await colorsByFrame(), legacyFrameColors, 'Sharing a legacy color command changed colors on another frame');
        await page.evaluate(async data => {
            const model = new $3Dmol.GLModel(0);
            model.addMolData(data, 'xyz', { multimodel: true, frames: true });
            model.setStyle({}, { sphere: {} });
            const internal = model.getInternalState();
            internal.frames[0][0].color = 0;
            internal.frames[0][1].color = 0x112233;
            internal.frames[1][0].color = 0xabcdef;
            internal.frames[1][1].color = 0;
            const scene = { models: [{ kind: 'snapshot', value: 'Colored snapshot', data: { encoding: 'text', value: JSON.stringify(internal) } }], selections: [], surfaces: [], view: { frame: 0 } };
            await viewerApp.loadURL(ViewerURL.serialize(scene, location.href));
        }, xyz);
        const expectedSnapshotColors = [[0, 0x112233], [0xabcdef, 0]];
        assert.deepEqual(await colorsByFrame(), expectedSnapshotColors, 'Snapshot import replaced black or distinct frame colors');
        await page.evaluate(() => viewerApp.loadURL(viewerApp.syncURL()));
        assert.deepEqual(await colorsByFrame(), expectedSnapshotColors, 'Snapshot share/restore lost per-frame colors');

        // Preserve imported residue-label options through real editing and hiding.
        await navigate('?pdb=1YCR&style=cartoon&select=chain:B&style=stick&labelres=fontSize:14;backgroundOpacity:0.8');
        const firstURL = page.url();
        const originalName = await page.evaluate(() => {
            const { state } = viewerApp;
            const entry = [...state.selections.values()].find(entry => entry.styles.size && entry.spec.chain === 'B');
            window.testSelectionId = entry.id;
            const card = state.ui.element.querySelector(`[data-selection-id="${entry.id}"].molui-selection-card`);
            const opacity = card.querySelector('[data-focus$=":opacity"]');
            opacity.closest('details').open = true;
            opacity.value = '0.35'; opacity.dispatchEvent(new Event('change', { bubbles: true }));
            return entry.name;
        });
        const opacityURL = await waitURL(firstURL);
        await page.evaluate(() => {
            const scene = ViewerURL.parse(location.href);
            if (!scene.selections.some(entry => entry.styles.some(style => style.spec.stick?.opacity === 0.35))) throw new Error('Style editing did not update URL');
            const name = viewerApp.state.ui.element.querySelector(`[data-selection-id="${testSelectionId}"].molui-selection-name`)
                || viewerApp.state.ui.element.querySelector(`[data-selection-id="${testSelectionId}"] .molui-selection-name`);
            name.value = 'Named & shareable + pocket'; name.dispatchEvent(new Event('change', { bubbles: true }));
        });
        const renamedURL = await waitURL(opacityURL);
        await page.evaluate(() => history.back());
        await page.waitForFunction(name => [...viewerApp.state.selections.values()].some(entry => entry.spec.chain === 'B' && entry.styles.size && entry.name === name), {}, originalName);
        assert.equal(page.url(), opacityURL, 'Back did not restore the prior URL');
        await page.evaluate(() => history.forward());
        await page.waitForFunction(() => [...viewerApp.state.selections.values()].some(entry => entry.name === 'Named & shareable + pocket'));
        assert.equal(page.url(), renamedURL, 'Forward did not restore renamed selection');
        await page.evaluate(() => {
            const { state, viewer } = viewerApp;
            const labelEntry = [...state.selections.values()].find(entry => entry.labelsEnabled);
            state.setSelectionVisible(labelEntry.id, false);
            state.setSelectionVisible(labelEntry.id, true);
            if (labelEntry.labelStyle.fontSize !== 14 || labelEntry.labelStyle.backgroundOpacity !== 0.8
                || !viewer.labels.some(label => label.stylespec.fontSize === 14)) throw new Error('Label toggling lost imported options');
        });

        // Load a trajectory URL through the UI; it must remain a reproducible source.
        const trajectoryURL = origin + '/frames.xyz';
        await page.evaluate(url => {
            const form = viewerApp.state.ui.element.querySelector('form');
            const source = form.querySelector('select');
            source.value = 'url'; source.dispatchEvent(new Event('change', { bubbles: true }));
            form.querySelector('input[type=text]').value = url;
            form.requestSubmit();
        }, trajectoryURL);
        await page.waitForFunction(() => viewerApp.viewer.getNumFrames() === 2);
        const trajectorySnapshot = await page.evaluate(() => viewerApp.getSnapshot());
        const trajectory = trajectorySnapshot.models.find(source => source.kind === 'url' && source.value === trajectoryURL);
        assert.ok(trajectory, 'UI-loaded trajectory lost its source URL');
        assert.equal(trajectory.format, 'xyz');
        assert.equal(trajectory.frames, true);
        assert.equal(trajectory.data, undefined, 'URL-loaded structures should not embed file data');
        assert.ok(await page.evaluate(() => !!viewerApp.viewer.getModel().selectedAtoms({})[0].style.stick),
            'An earlier global URL style replaced the newly loaded model default');

        // Capture one surface at frame zero, then share while viewing frame one.
        const surfaceAtZero = await page.evaluate(async () => {
            const { viewer, state } = viewerApp;
            const model = viewer.getModel();
            await viewer.setFrame(0);
            const id = await state.addSurface(1, { color: '#000066', opacity: 0.6 }, { model: model.getID(), index: 0 }, 'Motion surface', { model: model.getID() });
            const entry = state.surfaces.get(id);
            const surface = viewer.getSurface(entry.surfaceId);
            return surface[0].geo.geometryGroups.flatMap(group => Array.from(group.vertexArray));
        });
        await page.click('button[aria-label="Last frame"]');
        await page.waitForFunction(() => viewerApp.viewer.getFrame() === 1);
        await page.select('.molui-playback-speed', '24');
        await page.evaluate(() => {
            const root = viewerApp.state.ui.element;
            const loop = root.querySelector('.molui-trajectory-options input[type=checkbox]');
            loop.checked = false; loop.dispatchEvent(new Event('change', { bubbles: true }));
            const field = name => [...root.querySelectorAll('label.molui-field')].find(label => label.firstElementChild.textContent === name).querySelector('input,select');
            const background = field('Background'); background.value = '#000066'; background.dispatchEvent(new Event('input', { bubbles: true }));
            const projection = field('Projection'); projection.value = 'orthographic'; projection.dispatchEvent(new Event('change', { bubbles: true }));
            viewerApp.viewer.rotate(20, 'y'); viewerApp.viewer.render();
        });
        const sharedURL = await page.evaluate(() => viewerApp.syncURL());
        const shared = await page.evaluate(() => viewerApp.getSnapshot());
        const sharedAtoms = await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).map(atom => ({ color: atom.color, style: atom.style })));
        assert.equal(shared.view.frame, 1);
        assert.equal(shared.view.fps, 24);
        assert.equal(shared.view.loop, false);
        assert.equal(shared.view.background, '#000066');
        assert.equal(shared.view.projection, 'orthographic');
        assert.equal(shared.surfaces.find(surface => surface.name === 'Motion surface').frame, 0);
        assert.equal(shared.surfaces.find(surface => surface.name === 'Motion surface').spec.index, 0);
        assert.equal(shared.surfaces.find(surface => surface.name === 'Motion surface').calculationSpec.index, undefined);
        await page.evaluate(url => { window.__viewerDraws = []; return viewerApp.loadURL(url); }, sharedURL || page.url());
        await ready();
        const savedViewDraws = await page.evaluate(() => window.__viewerDraws);
        const firstSavedRender = savedViewDraws.find(entry => entry.method === 'render');
        assert.ok(firstSavedRender, 'Restoring the saved scene did not render a molecule');
        assert.equal(firstSavedRender.frame, 1, 'Saved scene first rendered at the wrong trajectory frame');
        for (const draw of savedViewDraws.filter(entry => entry.method === 'render')) {
            assert.deepEqual(draw.camera, shared.view.camera, 'Saved camera was applied after a molecule render');
            assert.deepEqual(draw.slab, shared.view.slab, 'Saved clipping planes were applied after a molecule render');
        }
        const restored = await page.evaluate(() => {
            const { viewer, state } = viewerApp;
            const entry = [...state.surfaces.values()].find(surface => surface.name === 'Motion surface');
            return {
                snapshot: viewerApp.getSnapshot(),
                vertices: viewer.getSurface(entry.surfaceId)[0].geo.geometryGroups.flatMap(group => Array.from(group.vertexArray)),
                frame: viewer.getFrame(),
                atoms: viewer.selectedAtoms({}).map(atom => ({ color: atom.color, style: atom.style })),
                panels: document.querySelectorAll('.molui').length
            };
        });
        assert.equal(restored.panels, 1, 'Restoring URL left duplicate UI instances');
        assert.equal(restored.frame, 1);
        assert.deepEqual(restored.vertices, surfaceAtZero, 'Surface geometry was rebuilt at the wrong trajectory frame');
        assert.deepEqual(restored.snapshot.view, shared.view, 'View/playback settings did not round-trip');
        assert.deepEqual(restored.snapshot.models, shared.models, 'Model source URLs did not round-trip');
        assert.deepEqual(restored.atoms, sharedAtoms, 'Adding a model changed its appearance only after restoring the URL');
        if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: '/tmp/3dmol-viewer-url.png' });

        // A removed first model leaves a hole in GLViewer IDs. Saved selections
        // and surfaces must target the surviving model's new index on restore.
        const afterRemoval = await page.evaluate(() => {
            viewerApp.state.removeModel(viewerApp.viewer.getModel(0));
            return {
                url: viewerApp.syncURL(),
                snapshot: viewerApp.getSnapshot(),
                atoms: viewerApp.viewer.selectedAtoms({}).map(atom => ({ color: atom.color, style: atom.style }))
            };
        });
        assert.equal(afterRemoval.snapshot.models.length, 1);
        assert.equal(afterRemoval.snapshot.surfaces.find(surface => surface.name === 'Motion surface').spec.model, 0);
        assert.ok(afterRemoval.snapshot.selections.some(selection => selection.spec.model === 0));
        await page.evaluate(url => viewerApp.loadURL(url), afterRemoval.url);
        assert.deepEqual(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).map(atom => ({ color: atom.color, style: atom.style }))), afterRemoval.atoms,
            'Removing a model changed remaining appearance after sharing');
        assert.equal(await page.evaluate(() => viewerApp.getSnapshot().surfaces[0].spec.model), 0);

        await page.setViewport({ width: 320, height: 480 });
        await page.evaluate(() => { viewerApp.viewer.resize(); viewerApp.state.updateUI(); });
        assert.equal(await page.evaluate(() => {
            const root = viewerApp.state.ui.element, body = root.querySelector('.molui-body');
            return root.scrollWidth > root.clientWidth || body.scrollWidth > body.clientWidth;
        }), false, 'Standalone viewer overflowed at mobile width');
        if (process.env.UI_SCREENSHOTS) await page.screenshot({ path: '/tmp/3dmol-viewer-url-mobile.png' });
        const sessionPayload = await page.evaluate(() => JSON.parse(JSON.stringify(viewerApp.viewer.getInternalState())));
        const sessionFrame = await page.evaluate(() => viewerApp.viewer.getFrame());

        // The learning-environment tutorial's session URLs remain usable. The
        // transport is mocked, while the page loads its real session adapter.
        await navigate('?session=URL+classroom%20room');
        await page.waitForFunction(() => window.__sessionEvents?.some(event => event.name === 'join session event'));
        await page.evaluate(async ({ payload, frame }) => {
            const join = __sessionEvents.find(event => event.name === 'join session event');
            if (join.data.name !== 'URL+classroom room') throw new Error('Session URL did not preserve a literal plus and decode an encoded space');
            __sessionHandlers['join session response']('test-session-client');
            if (getComputedStyle(document.getElementById('sessionmonitor')).display === 'none') throw new Error('Session response did not update controls');
            await __sessionHandlers['viewer state change response'](payload);
            if (viewerApp.viewer.getNumFrames() !== 2 || !viewerApp.viewer.selectedAtoms({}).length) throw new Error('Session model state was not applied');
            if (viewerApp.viewer.getFrame() !== frame) throw new Error('Older session state did not infer the displayed trajectory frame');
            await __sessionHandlers['viewer state change response']({ ...payload, frame });
            if (viewerApp.viewer.getFrame() !== frame) throw new Error('Session state did not retain its explicit trajectory frame');
            const coordinates = atoms => atoms.map(atom => [atom.x, atom.y, atom.z]);
            if (JSON.stringify(coordinates(viewerApp.viewer.selectedAtoms({}))) !== JSON.stringify(coordinates(payload.models.filter(Boolean).flatMap(model => model.atoms)))) {
                throw new Error('Session state changed the displayed trajectory coordinates');
            }
            // Socket updates can arrive in one event-loop turn. The final state
            // must be applied even while an earlier update is awaiting setFrame.
            const earlier = __sessionHandlers['viewer state change response']({ ...payload, frame: frame === 0 ? 1 : 0 });
            const latest = __sessionHandlers['viewer state change response']({ ...payload, frame });
            await Promise.all([earlier, latest]);
            if (viewerApp.viewer.getFrame() !== frame) throw new Error('A classroom update burst dropped the latest state');
            if (JSON.stringify(coordinates(viewerApp.viewer.selectedAtoms({}))) !== JSON.stringify(coordinates(payload.models.filter(Boolean).flatMap(model => model.atoms)))) {
                throw new Error('A classroom update burst left stale atom coordinates');
            }
            document.getElementById('sessionleave').click();
            if (viewerApp.state.ui.element.hidden) throw new Error('Leaving a session did not restore controls');
        }, { payload: sessionPayload, frame: sessionFrame });
        const sessionScene = await page.evaluate(() => ({
            snapshot: viewerApp.getSnapshot(), url: viewerApp.syncURL(),
            atoms: viewerApp.viewer.selectedAtoms({}).map(atom => ({ x: atom.x, y: atom.y, z: atom.z, color: atom.color, style: atom.style }))
        }));
        assert.equal(sessionScene.snapshot.models[0].kind, 'snapshot', 'Session models need a self-contained source after leaving');
        assert.equal(sessionScene.snapshot.view.frame, sessionFrame, 'Leaving a session lost its trajectory frame');
        await page.evaluate(url => viewerApp.loadURL(url), sessionScene.url);
        assert.deepEqual(await page.evaluate(() => viewerApp.viewer.selectedAtoms({}).map(atom => ({ x: atom.x, y: atom.y, z: atom.z, color: atom.color, style: atom.style }))), sessionScene.atoms,
            'A model received through a session did not survive share/restore after leaving');
        await page.evaluate(async payload => {
            await viewerApp.loadURL('?cid=241&style=stick');
            const before = viewerApp.viewer.selectedAtoms({}).length;
            await __sessionHandlers['viewer state change response'](payload);
            if (viewerApp.viewer.selectedAtoms({}).length !== before) throw new Error('Late session state replaced standalone model');
        }, sessionPayload);

        // A room update already in flight must also be invalidated when the user
        // opens a standalone URL. Hold its final frame promise deterministically.
        await navigate('?session=Pending+room');
        await page.waitForFunction(() => window.__sessionEvents?.some(event => event.name === 'join session event'));
        await page.evaluate(async payload => {
            __sessionHandlers['join session response']('pending-session-client');
            const viewer = viewerApp.viewer, originalSetFrame = viewer.setFrame;
            let release, entered, held = false;
            const enteredFrame = new Promise(resolve => { entered = resolve; });
            viewer.setFrame = function (...args) {
                const applied = originalSetFrame.apply(this, args);
                if (held) return applied;
                held = true;
                return applied.then(() => new Promise(resolve => { release = resolve; entered(); }));
            };
            const pending = __sessionHandlers['viewer state change response']({ ...payload, frame: 1 });
            try {
                await enteredFrame;
                await viewerApp.loadURL('?cid=241&style=stick');
                release();
                await pending;
                const snapshot = viewerApp.getSnapshot();
                if (snapshot.models.length !== 1 || snapshot.models[0].kind !== 'cid' || snapshot.models[0].value !== '241') {
                    throw new Error('An in-flight session update replaced the standalone URL source');
                }
                if (viewer.selectedAtoms({}).length !== 12 || viewer.getFrame() !== 0 || !viewerApp.state.ui.element.isConnected
                    || viewerApp.state.ui.element.hidden || document.querySelectorAll('.molui').length !== 1) {
                    throw new Error('An in-flight session update left stale coordinates or controls after URL navigation: ' + JSON.stringify({
                        atoms: viewer.selectedAtoms({}).length, frame: viewer.getFrame(), connected: viewerApp.state.ui.element.isConnected,
                        hidden: viewerApp.state.ui.element.hidden, panels: document.querySelectorAll('.molui').length
                    }));
                }
                await __sessionHandlers['viewer state change response'](payload);
                if (viewer.selectedAtoms({}).length !== 12) throw new Error('A detached room still accepted new state');
            } finally {
                release?.();
                viewer.setFrame = originalSetFrame;
            }
        }, sessionPayload);

        assert.deepEqual(unexpectedRequests, [], 'Viewer requested an unmocked external resource');
        assert.deepEqual(errors, [], 'Viewer produced browser errors');
        console.log(`Standalone viewer URL browser checks passed: ${examples.length} tutorial examples, source formats and detection, literal-plus URLs, proxy fallbacks, UI→URL, Back/Forward, cancellation, labels, trajectory URLs, view settings, model remapping, surface snapshots, per-frame colors, mobile layout, and classroom session races.`);
    } finally {
        if (browser) await browser.close();
        await new Promise(resolve => server.close(resolve));
    }
})().catch(error => { console.error(error); process.exitCode = 1; });
