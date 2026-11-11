/* Standalone URL adapter. URL syntax and browser history stay outside the UI library. */
(function (root, factory) {
    const api = factory(root);
    if (typeof module === 'object' && module.exports) module.exports = api;
    else {
        root.ViewerURL = api;
        const start = () => {
            if (!root.document.getElementById('gldiv')) return;
            try {
                root.viewerApp = api.createApp();
                root.viewerApp.ready.catch(error => root.viewerApp.showError(error));
            } catch (error) {
                const message = root.document.getElementById('viewer-message');
                message.hidden = false; message.textContent = error.message;
            }
        };
        if (root.document.readyState === 'loading') root.document.addEventListener('DOMContentLoaded', start, { once: true });
        else start();
    }
})(typeof window !== 'undefined' ? window : globalThis, function (root) {
    'use strict';
    const baseURL = 'https://3dmol.org/viewer.html';
    const sourceKinds = new Set(['pdb', 'cif', 'cid', 'mmtf', 'url', 'file', 'snapshot']);
    const forbiddenKeys = new Set(['__proto__', 'prototype', 'constructor']);
    const copy = value => JSON.parse(JSON.stringify(value));
    function assertObject(value, name) {
        if (!value || typeof value !== 'object' || Array.isArray(value)) throw new Error(`${name} must be an object.`);
    }
    function safeJSON(value) {
        if (value && typeof value === 'object') {
            for (const key of Object.keys(value)) {
                if (forbiddenKeys.has(key)) throw new Error(`Unsupported setting: ${key}.`);
                safeJSON(value[key]);
            }
        } else if (typeof value === 'number' && !Number.isFinite(value)) throw new Error('Settings must contain finite numbers.');
    }
    function asURL(input, base = root.location?.href || baseURL) {
        let value = typeof input === 'string' ? input : input?.href || '';
        if (/^(?:pdb|cif|cid|mmtf|url|session|state|select|style|background)=/i.test(value)) value = '?' + value;
        return new URL(value || base, base);
    }
    function inferFormat(value) {
        const pathname = new URL(value, baseURL).pathname.toLowerCase();
        const compressed = pathname.endsWith('.gz');
        const match = (compressed ? pathname.slice(0, -3) : pathname).match(/\.([^./]+)$/);
        // The core owns its parser registry and can detect a format from content
        // when the extension is missing or unknown. Keep compression for either case.
        return (match?.[1] || '') + (compressed ? '.gz' : '');
    }
    function queryEntries(query) {
        return query.split('&').filter(Boolean).map(part => {
            const separator = part.indexOf('=');
            // Legacy links treat + literally. Decode each value only after splitting
            // so encoded URL delimiters remain part of the source or specification.
            return [decodeURIComponent(separator < 0 ? part : part.slice(0, separator)),
                decodeURIComponent(separator < 0 ? '' : part.slice(separator + 1))];
        });
    }
    function validateModel(model) {
        assertObject(model, 'Model source');
        if (!sourceKinds.has(model.kind) || typeof model.value !== 'string' || !model.value.trim()) throw new Error('Choose a valid model source.');
        if (['pdb', 'cif', 'mmtf'].includes(model.kind) && !/^[a-z0-9]{4}$/i.test(model.value)) throw new Error('A PDB identifier must contain four letters or digits.');
        if (model.kind === 'cid' && !/^\d+$/.test(model.value)) throw new Error('A PubChem CID must be numeric.');
        if (model.kind === 'url' && !['https:', 'http:'].includes(new URL(model.value, baseURL).protocol)) throw new Error('Structure URLs must use http or https.');
        if (model.format !== undefined && typeof model.format !== 'string') throw new Error('A structure format must be a string.');
        if (['file', 'snapshot'].includes(model.kind) && (!model.data || !['text', 'base64'].includes(model.data.encoding) || typeof model.data.value !== 'string')) throw new Error('A local model URL must include its structure data.');
        if (model.options !== undefined) assertObject(model.options, 'Parser options');
    }
    function validateScene(scene) {
        assertObject(scene, 'Viewer state'); safeJSON(scene);
        for (const key of ['models', 'selections', 'surfaces']) if (!Array.isArray(scene[key])) throw new Error(`Viewer state needs a ${key} list.`);
        scene.models.forEach(validateModel);
        for (const selection of scene.selections) {
            assertObject(selection.spec, 'Atom selection');
            if (!Array.isArray(selection.styles)) throw new Error('A selection needs a styles list.');
            selection.styles.forEach(style => assertObject(style.spec, 'Representation'));
            if (selection.labelStyle !== undefined) assertObject(selection.labelStyle, 'Residue label style');
        }
        for (const surface of scene.surfaces) {
            assertObject(surface.spec, 'Surface selection'); assertObject(surface.style, 'Surface style');
            if (![1, 2, 3, 4].includes(surface.type)) throw new Error('Unsupported surface type.');
            if (surface.calculationSpec !== undefined) assertObject(surface.calculationSpec, 'Surface calculation selection');
        }
        if (scene.view !== undefined) {
            assertObject(scene.view, 'View settings');
            if (scene.view.projection !== undefined && !['perspective', 'orthographic'].includes(scene.view.projection)) throw new Error('Unsupported projection.');
            if (scene.view.camera && (!Array.isArray(scene.view.camera) || scene.view.camera.length !== 8 || scene.view.camera.some(value => typeof value !== 'number' || !Number.isFinite(value)))) throw new Error('The camera needs eight finite numbers.');
            if (scene.view.frame !== undefined && !Number.isInteger(scene.view.frame)) throw new Error('The trajectory frame must be an integer.');
            if (scene.view.backgroundAlpha !== undefined && (typeof scene.view.backgroundAlpha !== 'number' || scene.view.backgroundAlpha < 0 || scene.view.backgroundAlpha > 1)) throw new Error('Background opacity must be between zero and one.');
            if (scene.view.fps !== undefined && (typeof scene.view.fps !== 'number' || scene.view.fps <= 0 || scene.view.fps > 240)) throw new Error('Playback speed must be between zero and 240 fps.');
            if (scene.view.slab !== undefined) {
                assertObject(scene.view.slab, 'Clipping settings');
                if (!Number.isFinite(scene.view.slab.near) || !Number.isFinite(scene.view.slab.far)) throw new Error('Clipping distances must be finite numbers.');
            }
        }
        return scene;
    }
    /** Read established ordered query commands or a complete shared scene. */
    function parse(input, parseSpec = root.$3Dmol?.specStringToObject) {
        const url = asURL(input);
        const fragment = new URLSearchParams(url.hash.slice(1).replace(/\+/g, '%2B')).get('state');
        const commands = fragment === null ? queryEntries(url.search.slice(1)) : [];
        const shared = fragment ?? commands.find(([key]) => key === 'state')?.[1];
        if (shared !== undefined) {
            let scene;
            try { scene = JSON.parse(shared); } catch { throw new Error('The shared viewer state is not valid JSON.'); }
            return validateScene(scene);
        }
        const scene = { models: [], selections: [], surfaces: [], view: { background: '#ffffff', backgroundAlpha: 1, projection: 'perspective' } };
        let selection = {}, current;
        const readSpec = value => {
            if (!value || value === 'all') return {};
            let parsed;
            if (/^[\[{]/.test(value.trim())) {
                try { parsed = JSON.parse(value); } catch { throw new Error('Invalid JSON selection or style.'); }
            } else {
                if (!parseSpec) throw new Error('The 3Dmol specification parser is unavailable.');
                parsed = parseSpec(value);
            }
            assertObject(parsed, 'Selection or style'); safeJSON(parsed); return parsed;
        };
        const newSelection = () => {
            const entry = { name: `Selection ${scene.selections.length + 1}`, spec: copy(selection), styles: [] };
            scene.selections.push(entry); current = entry; return entry;
        };
        for (const [rawKey, value] of commands) {
            const key = rawKey.toLowerCase();
            if (sourceKinds.has(key) && !['file', 'snapshot'].includes(key)) scene.models.push({ kind: key, value, frames: true });
            else if (key === 'type' || key === 'frames') {
                const model = scene.models[scene.models.length - 1];
                if (!model) throw new Error(`${key} must follow a structure source.`);
                if (key === 'type') model.format = value.toLowerCase();
                else model.frames = !['0', 'false'].includes(value.toLowerCase());
            } else if (key === 'select') { selection = readSpec(value); current = undefined; }
            else if (key === 'style' || key === 'addstyle') {
                const entry = newSelection(); entry.additive = key === 'addstyle';
                entry.styles = Object.entries(readSpec(value)).map(([kind, settings]) => ({ spec: { [kind]: settings } }));
                if (!entry.styles.length) entry.styles.push({ spec: {} });
            } else if (key === 'labelres') {
                const entry = current && !current.labelsEnabled ? current : newSelection();
                entry.labelsEnabled = true; entry.labelStyle = readSpec(value);
            } else if (key === 'surface') scene.surfaces.push({ name: 'Molecular surface', type: 1, spec: copy(selection), calculationSpec: copy(selection), style: readSpec(value), colorOperationCount: scene.atomColors?.length || 0 });
            else if (key === 'colorbyelement' || key === 'colorbyproperty') {
                const settings = readSpec(value);
                (scene.atomColors ||= []).push(key === 'colorbyelement'
                    ? { type: 'element', spec: copy(selection), colorscheme: settings.colorscheme || value }
                    : { type: 'property', spec: copy(selection), prop: settings.prop, scheme: settings.scheme, ...(settings.range ? { range: settings.range } : {}) });
            } else if (key === 'session') {
                if (!value.trim()) throw new Error('A session name is required.'); scene.session = value;
            } else if (['background', 'backgroundalpha', 'projection', 'frame', 'fps', 'loop', 'playing', 'view', 'open'].includes(key)) {
                if (key === 'view') scene.view.camera = value.trim().startsWith('[') ? JSON.parse(value) : value.split(',').map(Number);
                else if (key === 'open') scene.open = !['0', 'false'].includes(value.toLowerCase());
                else if (['frame', 'fps', 'backgroundalpha'].includes(key)) scene.view[key === 'backgroundalpha' ? 'backgroundAlpha' : key] = Number(value);
                else if (key === 'loop' || key === 'playing') scene.view[key] = !['0', 'false'].includes(value.toLowerCase());
                else scene.view[key] = value;
            } else throw new Error(`Unknown viewer URL parameter: ${rawKey}.`);
        }
        for (const model of scene.models) {
            model.format ||= model.kind === 'url' ? inferFormat(model.value) : model.kind === 'cid' ? 'sdf' : model.kind === 'mmtf' ? 'cif' : model.kind;
            if (model.kind === 'url' && new URL(model.value, baseURL).pathname.endsWith('.gz') && !model.format.endsWith('.gz')) model.format += '.gz';
        }
        if (scene.models.length && !scene.selections.some(entry => entry.styles.length)) scene.selections.unshift({ name: 'All atoms', spec: {}, styles: [{ spec: { line: {} } }] });
        if (!scene.atomColors?.length) scene.surfaces.forEach(surface => delete surface.colorOperationCount);
        return validateScene(scene);
    }
    /** Fragments do not send local file data or long settings to the web server. */
    function serialize(scene, base = root.location?.href || baseURL) {
        validateScene(scene);
        const url = new URL(base, baseURL);
        url.search = ''; url.hash = 'state=' + encodeURIComponent(JSON.stringify(scene));
        return url.href;
    }
    function encodeData(data) {
        if (typeof data === 'string') return { encoding: 'text', value: data };
        const bytes = new Uint8Array(data); let binary = '';
        for (let i = 0; i < bytes.length; i += 8192) binary += String.fromCharCode(...bytes.subarray(i, i + 8192));
        return { encoding: 'base64', value: root.btoa(binary) };
    }
    function decodeData(data) {
        return data.encoding === 'text' ? data.value : Uint8Array.from(root.atob(data.value), character => character.charCodeAt(0)).buffer;
    }
    function modelList(viewer) {
        const models = [];
        for (let id = 0; id < viewer.getNextModelId(); id++) {
            const model = viewer.getModel(id); if (model && model.getID() === id) models.push(model);
        }
        return models;
    }
    function baselineStyles(model) {
        const groups = new Map();
        for (const atom of model.selectedAtoms({})) {
            const key = JSON.stringify(atom.style || {});
            if (!groups.has(key)) groups.set(key, { indices: [], style: copy(atom.style || {}) });
            groups.get(key).indices.push(atom.index);
        }
        return [...groups.values()];
    }
    function normalizeSelection(spec, indices, modelCount) {
        if (Array.isArray(spec)) return spec.map(value => normalizeSelection(value, indices, modelCount));
        if (!spec || typeof spec !== 'object') return spec;
        const result = {};
        for (const [key, value] of Object.entries(spec)) {
            if (key === 'model') {
                const index = item => {
                    const id = typeof item === 'number' ? item < 0 ? modelCount + item : item : item.getID();
                    if (!indices.has(id)) throw new Error('A selection refers to a model that is no longer loaded.');
                    return indices.get(id);
                };
                result.model = Array.isArray(value) ? value.map(index) : index(value);
            } else result[key] = normalizeSelection(value, indices, modelCount);
        }
        return result;
    }
    function createApp(options = {}) {
        const win = options.window || root, document = win.document, lib = options.library || win.$3Dmol;
        if (!lib?.StateManager) throw new Error('Load the 3Dmol core and UI bundles before viewer.js.');
        const viewer = options.viewer || lib.createViewer('gldiv', { defaultcolors: lib.rasmolElementColors });
        win.glviewer = viewer;
        const resources = new Map(), cache = new Map(), surfaceColors = new Map(), message = document.getElementById('viewer-message');
        let state, unsubscribe, controller, generation = 0, suspended = true, disposed = false;
        let restoringURL = false, sessionUpdate;
        let syncTimer, pendingPush = false, liveKey = '', lastURL = win.location.href, loadingURL;
        let session = null, colorOperations = [], sessionScripts;
        let clearPicking = () => {};
        const app = { viewer, get state() { return state; }, ready: null, getSnapshot, loadURL, syncURL, showError, acceptSessionState, sessionChanged, destroy };
        function showError(error) {
            if (!message || disposed) return;
            message.hidden = !error; message.dataset.error = String(!!error);
            message.textContent = error ? (error.message || String(error)) : '';
        }
        function makeState(scene = {}) {
            state = new lib.StateManager(viewer, {
                open: scene.open ?? (!modelList(viewer).length && !scene.session ? true : undefined),
                allowLocalFiles: false,
                view: scene.view, onOpenChange: () => scheduleSync(true),
                onModelLoaded(model, source) {
                    const descriptor = { kind: source.kind, value: source.value || source.name, format: source.format, frames: source.frames, title: source.name, baseStyles: baselineStyles(model) };
                    if (source.kind === 'file') descriptor.data = encodeData(source.data);
                    resources.set(model, descriptor);
                    cache.set(resourceKey(descriptor), { data: source.data, format: source.format });
                    // Make the loader's default appearance an explicit last layer.
                    // Earlier all-atom selections must not override it only on reload.
                    state.batch(() => {
                        descriptor.baseStyles.forEach((group, index) => {
                            const selection = { model: model.getID() };
                            if (descriptor.baseStyles.length > 1) selection.index = group.indices;
                            const id = state.createSelectionAndStyle(selection, group.style);
                            state.renameSelection(id, descriptor.baseStyles.length > 1 ? `${source.name} · ${index + 1}` : source.name);
                        });
                    });
                    installPicking();
                }
            });
            unsubscribe = state.subscribe(() => { rememberSurfaceColors(); scheduleSync(true); });
        }
        function captureColors() { return modelList(viewer).map(model => ({ model, colors: model.selectedAtoms({}).map(atom => atom.color) })); }
        function applyColors(records) {
            for (const { model, colors } of records) if (model) {
                for (const atom of model.selectedAtoms({})) if (colors[atom.index] !== undefined) atom.color = colors[atom.index];
            }
        }
        function rememberSurfaceColors() {
            for (const id of surfaceColors.keys()) if (!state.surfaces.has(id)) surfaceColors.delete(id);
            for (const id of state.surfaces.keys()) if (!surfaceColors.has(id)) surfaceColors.set(id, captureColors());
        }
        function getSnapshot() {
            if (session) return { models: [], selections: [], surfaces: [], session };
            const models = modelList(viewer), indices = new Map(models.map((model, index) => [model.getID(), index]));
            const spec = value => normalizeSelection(value, indices, viewer.getNextModelId());
            return {
                models: models.map(model => {
                    const source = resources.get(model);
                    if (!source) throw new Error('This model has no shareable source. Load it through the viewer controls.');
                    return { ...copy(source), title: state.modelTitles.get(model.getID()) || source.title || source.value, hidden: model.isHidden() };
                }),
                selections: [...state.selections.values()].map(entry => ({
                    name: entry.name, spec: spec(entry.spec), styles: [...entry.styles.values()].map(style => ({ spec: copy(style.spec), hidden: style.hidden })),
                    hidden: entry.hidden, additive: !!entry.additive, labelsEnabled: entry.labelsEnabled,
                    ...(entry.labelStyle !== undefined ? { labelStyle: copy(entry.labelStyle) } : {}),
                    ...(entry.labelFrame !== undefined ? { labelFrame: entry.labelFrame } : {})
                })),
                surfaces: [...state.surfaces.values()].map(surface => ({
                    name: surface.name, spec: spec(surface.spec), calculationSpec: spec(surface.calculationSpec), style: copy(surface.style), type: surface.type, frame: surface.frame ?? 0,
                    atomColors: (surfaceColors.get(surface.id) || captureColors()).filter(record => indices.has(record.model.getID()))
                        .map(record => ({ model: indices.get(record.model.getID()), colors: copy(record.colors) }))
                })),
                view: { ...state.getViewSettings(), camera: viewer.getView(), slab: viewer.getSlab(), ...state.ui.getTrajectorySettings() }, open: state.ui.isOpen
            };
        }
        function scheduleSync(push) {
            if (suspended || disposed) return;
            pendingPush ||= !!push;
            if (!syncTimer) syncTimer = win.setTimeout(() => {
                syncTimer = undefined;
                try { syncURL(pendingPush); } catch (error) { showError(error); }
                pendingPush = false;
            }, 180);
        }
        function syncURL(push = pendingPush) {
            if (suspended || disposed) return win.location.href;
            const url = serialize(getSnapshot(), win.location.href);
            if (url !== lastURL) { win.history[push ? 'pushState' : 'replaceState'](null, '', url); lastURL = url; }
            return url;
        }
        function resourceKey(source) { return JSON.stringify([source.kind, source.value, source.format, source.data]); }
        async function fetchSource(source, signal) {
            const key = resourceKey(source);
            if (cache.has(key)) return cache.get(key);
            if (['file', 'snapshot'].includes(source.kind)) return { data: decodeData(source.data), format: source.kind === 'snapshot' ? 'snapshot' : source.format || inferFormat(source.value) };
            const id = source.value.toUpperCase();
            let format = source.format || (source.kind === 'url' ? inferFormat(source.value) : source.kind === 'cid' ? 'sdf' : source.kind === 'mmtf' ? 'cif' : source.kind);
            let url = source.kind === 'url' ? new URL(source.value, win.location.href).href
                : source.kind === 'cid' ? `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${source.value}/SDF?record_type=3d`
                    : `https://files.rcsb.org/view/${id}.${source.kind === 'mmtf' || source.kind === 'cif' || format === 'cif' ? 'cif' : 'pdb'}`;
            const request = options.fetch || win.fetch.bind(win);
            const read = async address => {
                const response = await request(address, { signal });
                if (!response.ok) throw new Error(`Unable to load ${source.value} (HTTP ${response.status}).`);
                const data = await (/(?:\.gz|bcif|mmtf)$/.test(format) ? response.arrayBuffer() : response.text());
                // echo.cgi returns HTTP 200 even when RCSB supplies an error page.
                // Treat missing atom records as a failed PDB candidate so CIF can load.
                if (source.kind === 'pdb' && format === 'pdb' && !/^(?:ATOM  |HETATM)/m.test(data)) {
                    throw new Error(`No PDB atom records found for ${source.value}.`);
                }
                return data;
            };
            const readWithProxy = async address => {
                try { return await read(address); }
                catch (error) {
                    if (signal.aborted) throw error;
                    const proxy = new URL('echo.cgi', win.location.href); proxy.searchParams.set('url', address);
                    return read(proxy.href);
                }
            };
            let data;
            try { data = await readWithProxy(url); }
            catch (error) {
                if (signal.aborted) throw error;
                if (source.kind !== 'pdb' || format !== 'pdb') throw error;
                format = 'cif'; url = `https://files.rcsb.org/view/${id}.cif`; data = await readWithProxy(url);
            }
            const result = { data, format }; cache.set(key, result); return result;
        }
        function loadURL(input = win.location.href) {
            const url = asURL(input, win.location.href).href;
            // History navigation can emit both popstate and hashchange.
            if (restoringURL && loadingURL === url && app.ready) return app.ready;
            loadingURL = url;
            const promise = restore(url); app.ready = promise; return promise;
        }
        async function restore(input) {
            const token = ++generation;
            restoringURL = true; sessionUpdate = undefined;
            let scene = {};
            controller?.abort(); controller = new win.AbortController();
            const signal = controller.signal;
            suspended = true; pendingPush = false; win.clearTimeout(syncTimer); syncTimer = undefined; showError(null);
            try {
                scene = parse(input, lib.specStringToObject);
                const sources = await Promise.all(scene.models.map(source => fetchSource(source, signal)));
                if (token !== generation || disposed) return;
                win.leaveViewerSession?.();
                unsubscribe?.(); state?.destroy(); state = undefined; viewer.stopAnimate(); viewer.clear(); resources.clear(); surfaceColors.clear();
                session = scene.session || null; colorOperations = copy(scene.atomColors || []);
                for (let i = 0; i < scene.models.length; i++) {
                    const source = scene.models[i], loaded = sources[i];
                    let model;
                    if (source.kind === 'snapshot') {
                        const internal = JSON.parse(loaded.data); safeJSON(internal);
                        if (!Array.isArray(internal.atoms) || !Array.isArray(internal.frames)) throw new Error('Invalid model snapshot.');
                        model = viewer.addModel(); model.setInternalState(internal);
                        for (const atoms of [internal.atoms, ...internal.frames]) if (Array.isArray(atoms)) {
                            const colors = atoms.map(atom => atom.color);
                            model.setAtomDefaults(atoms);
                            atoms.forEach((atom, index) => { if (colors[index] !== undefined) atom.color = colors[index]; });
                        }
                    } else model = source.frames === false ? viewer.addModel(loaded.data, loaded.format, copy(source.options || {}))
                        : viewer.addModelsAsFrames(loaded.data, loaded.format, copy(source.options || {}));
                    if (!model.selectedAtoms({}).length) { viewer.removeModel(model); throw new Error(`No atoms found in ${source.value}. Check its file format.`); }
                    if (source.kind !== 'snapshot') model.setStyle({}, { line: {} });
                    for (const group of source.baseStyles || []) model.setStyle({ index: group.indices }, group.style);
                    if (source.colors) {
                        const internal = model.getInternalState();
                        for (const atoms of [internal.atoms, ...internal.frames]) if (Array.isArray(atoms)) {
                            for (const atom of atoms) if (source.colors[atom.index] !== undefined) atom.color = source.colors[atom.index];
                        }
                    }
                    if (source.frameColors) model.getInternalState().frames.forEach((atoms, index) => {
                        const colors = source.frameColors[index];
                        if (colors && Array.isArray(atoms)) for (const atom of atoms) {
                            if (colors[atom.index] !== undefined) atom.color = colors[atom.index];
                        }
                    });
                    if (source.hidden) model.hide(); resources.set(model, { ...copy(source), format: loaded.format });
                }
                const initialColors = captureColors();
                applyAtomColors();
                // Layers, surfaces and trajectory controls render as they are
                // restored. Fit the requested frame before any of those renders,
                // then keep its camera fixed while rebuilding frame-specific data.
                if (scene.models.length) {
                    await viewer.setFrame(scene.view?.frame ?? 0);
                    if (token !== generation || disposed) return;
                }
                viewer.zoomTo();
                if (scene.view?.slab) viewer.setSlab(scene.view.slab.near, scene.view.slab.far);
                if (scene.view?.camera) viewer.setView(scene.view.camera);
                makeState(scene);
                modelList(viewer).forEach((model, i) => state.modelTitles.set(model.getID(), scene.models[i].title || scene.models[i].value));
                const labels = [];
                state.batch(() => {
                    for (const entry of scene.selections) {
                        const id = state.addSelection(entry.spec, undefined, entry.name);
                        if (entry.additive) state.setSelectionAdditive(id, true);
                        for (const style of entry.styles) { const styleId = state.addStyle(style.spec, id); if (style.hidden) state.toggleHideStyle(id, styleId); }
                        if (entry.hidden) state.setSelectionVisible(id, false);
                        if (entry.labelStyle !== undefined || entry.labelsEnabled) labels.push({ id, entry });
                    }
                });
                for (const surface of scene.surfaces) {
                    if (token !== generation || disposed) return;
                    await viewer.setFrame(surface.frame ?? scene.view?.frame ?? 0);
                    if (token !== generation || disposed) return;
                    const colors = captureColors();
                    if (surface.atomColors) applyColors(surface.atomColors.map(record => ({ model: viewer.getModel(record.model), colors: record.colors })));
                    else if (surface.colorOperationCount !== undefined) {
                        applyColors(initialColors);
                        applyAtomColors(colorOperations.slice(0, surface.colorOperationCount), false);
                    }
                    try {
                        if (viewer.selectedAtoms(surface.spec).length) await state.addSurface(surface.type, surface.style, surface.spec, surface.name, surface.calculationSpec ?? surface.spec);
                    } finally { applyColors(colors); }
                }
                if (token !== generation || disposed) return;
                for (const { id, entry } of labels) {
                    await viewer.setFrame(entry.labelFrame ?? scene.view?.frame ?? 0);
                    if (token !== generation || disposed) return;
                    state.setResidueLabels(id, !!entry.labelsEnabled, entry.labelStyle);
                    if (entry.labelFrame !== undefined) state.selections.get(id).labelFrame = entry.labelFrame;
                }
                await state.ui.setTrajectorySettings({ ...scene.view, frame: scene.view?.frame ?? 0 });
                if (token !== generation || disposed) return;
                installPicking(); viewer.render(); state.refresh();
                if (scene.session) { await enableSessions(); if (token === generation) win.joinSession(scene.session); }
                if (token !== generation || disposed) return;
                lastURL = win.location.href; liveKey = readLiveKey();
            } catch (error) {
                if (token === generation && !signal.aborted) { if (!state) makeState({ open: true }); showError(error); throw error; }
            } finally { if (token === generation) { suspended = false; restoringURL = false; } }
        }
        function applyAtomColors(operations = colorOperations, persist = true) {
            for (const operation of operations) {
                if (operation.type === 'element') {
                    const colors = lib.elementColors[operation.colorscheme] || lib.builtinColorSchemes[operation.colorscheme];
                    if (!colors) throw new Error(`Unknown element color scheme: ${operation.colorscheme}.`);
                    const map = colors.map || colors;
                    for (const atom of viewer.selectedAtoms(operation.spec)) if (map[atom.elem] !== undefined) atom.color = map[atom.elem];
                } else {
                    const Constructor = lib.Gradient[operation.scheme] || lib.Gradient.builtinGradients[operation.scheme];
                    if (!Constructor || !operation.prop) throw new Error('Choose a property and supported gradient for colorbyproperty.');
                    const atoms = viewer.selectedAtoms(operation.spec), gradient = new Constructor();
                    const values = atoms.map(atom => Number(lib.getAtomProperty(atom, operation.prop))).filter(Number.isFinite);
                    const range = operation.range || [Math.min(...values), Math.max(...values)];
                    for (const atom of atoms) { const value = lib.getAtomProperty(atom, operation.prop); if (value !== undefined && Number.isFinite(Number(value))) atom.color = gradient.valueToHex(Number(value), range); }
                }
            }
            // Persist the result per source, rather than replaying old selection
            // commands on models loaded later or on shifted model identifiers.
            if (persist && operations.length) for (const model of modelList(viewer)) {
                const source = resources.get(model), internal = model.getInternalState();
                if (internal.frames.length) source.frameColors = internal.frames.map(atoms => Array.isArray(atoms) ? atoms.map(atom => atom.color) : null);
                else source.colors = internal.atoms.map(atom => atom.color);
            }
        }
        function installPicking() {
            clearPicking();
            let hover, clicked;
            clearPicking = () => {
                const labels = [hover, clicked?.label].filter(Boolean);
                if (labels.length) viewer.removeLabels(labels);
                hover = undefined; clicked = undefined;
            };
            const label = atom => atom.resn && atom.resi !== undefined ? `${atom.resn}${atom.resi}:${atom.atom || atom.elem}` : atom.atom || atom.elem;
            viewer.setHoverable({}, true, atom => {
                if (hover) viewer.removeLabel(hover);
                if (clicked?.atom !== atom) hover = viewer.addLabel(label(atom), { position: atom, fontSize: 12, backgroundColor: 'black', backgroundOpacity: 0.5, alignment: 'bottomCenter' });
            }, () => { if (hover) viewer.removeLabel(hover); hover = undefined; });
            viewer.setClickable({}, true, atom => {
                if (hover) viewer.removeLabel(hover); hover = undefined;
                if (!clicked) clicked = { atom, label: viewer.addLabel(label(atom), { position: atom, fontSize: 12, backgroundColor: 'blue', backgroundOpacity: 0.75 }) };
                else {
                    const previous = clicked.atom; viewer.removeLabel(clicked.label); clicked = undefined;
                    if (previous !== atom) {
                        const start = new lib.Vector3(previous.x, previous.y, previous.z), end = new lib.Vector3(atom.x, atom.y, atom.z);
                        const distance = lib.GLShape.distance_from(start, end); let distanceLabel;
                        viewer.addCylinder({ dashed: true, radius: 0.05, dashLength: 0.25, gapLength: 0.15, start, end, fromCap: 2, toCap: 2, color: 'blue', clickable: true,
                            callback(shape) { viewer.removeShape(shape); viewer.removeLabel(distanceLabel); viewer.render(); } });
                        distanceLabel = viewer.addLabel(distance.toFixed(3), { position: start.clone().add(end).multiplyScalar(0.5), fontSize: 12, backgroundColor: 'blue', backgroundOpacity: 0.75 });
                    }
                }
                viewer.render();
            });
        }
        function readLiveKey() { return state ? JSON.stringify([viewer.getView(), viewer.getSlab(), state.ui.getTrajectorySettings()]) : ''; }
        const poll = win.setInterval(() => {
            if (suspended || disposed || !state || document.hidden || session) return;
            const next = readLiveKey(); if (next !== liveKey) { liveKey = next; scheduleSync(false); }
        }, 200);
        const navigate = () => { loadURL(win.location.href).catch(showError); };
        win.addEventListener('popstate', navigate); win.addEventListener('hashchange', navigate);
        function script(src) {
            return new Promise((resolve, reject) => {
                const element = document.createElement('script'); element.src = src;
                element.onload = resolve; element.onerror = () => reject(new Error('Unable to load classroom session support.'));
                document.head.append(element);
            });
        }
        async function enableSessions() {
            if (!sessionScripts) sessionScripts = (async () => {
                if (!win.jQuery) await script('https://code.jquery.com/jquery-3.7.1.min.js');
                if (!win.io) await script('https://cdnjs.cloudflare.com/ajax/libs/socket.io/4.6.2/socket.io.js');
                win.initSessions();
            })().catch(error => { sessionScripts = undefined; throw error; });
            return sessionScripts;
        }
        function acceptSessionState(internal) {
            if (disposed || restoringURL) return Promise.resolve();
            if (sessionUpdate) {
                sessionUpdate.pending = internal;
                return sessionUpdate.promise;
            }
            const update = { token: ++generation, pending: internal, promise: null };
            sessionUpdate = update; suspended = true;
            // Coalesce a burst of packets and keep the newest state while a frame
            // is loading. All callers wait until the newest queued state is ready.
            update.promise = Promise.resolve().then(() => applySessionUpdates(update));
            return update.promise;
        }
        async function applySessionUpdates(update) {
            try {
                while (update.pending && update.token === generation && !disposed) {
                    const internal = update.pending, camera = viewer.getView(); update.pending = undefined;
                    unsubscribe?.(); state?.destroy(); viewer.stopAnimate();
                    viewer.setInternalState(internal); viewer.setView(camera); resources.clear(); surfaceColors.clear();
                    modelList(viewer).forEach((model, index) => resources.set(model, {
                        kind: 'snapshot', value: `Session model ${index + 1}`, title: `Session model ${index + 1}`,
                        data: encodeData(JSON.stringify(model.getInternalState()))
                    }));
                    // Older classroom senders omit the frame index. Their current
                    // atoms still identify the displayed locally stored frame.
                    let frame = internal.frame;
                    if (!Number.isInteger(frame) || frame < 0) {
                        const longest = internal.models.filter(Boolean).sort((a, b) => b.frames.length - a.frames.length)[0];
                        frame = longest?.frames.findIndex(atoms => Array.isArray(atoms) && atoms.length === longest.atoms.length
                            && atoms.every((atom, index) => ['x', 'y', 'z', 'elem', 'color'].every(key => atom[key] === longest.atoms[index][key]))) ?? 0;
                    }
                    try { await viewer.setFrame(Math.max(0, frame)); }
                    catch (error) {
                        if (update.token !== generation || disposed) return;
                        if (update.pending) continue;
                        throw error;
                    }
                    if (update.token !== generation || disposed) return;
                    if (update.pending) continue;
                    makeState({ open: false }); if (session) state.ui.element.hidden = true;
                    else installPicking();
                    viewer.render();
                }
            } finally {
                if (sessionUpdate === update) sessionUpdate = undefined;
                if (update.token === generation) { suspended = false; if (!session) scheduleSync(true); }
            }
        }
        function sessionChanged(name, joined) {
            session = joined ? name : null;
            if (!session && sessionUpdate) sessionUpdate.pending = undefined;
            if (!session && !disposed) installPicking();
            if (!suspended) scheduleSync(true);
        }
        document.getElementById('copy-link')?.addEventListener('click', async () => {
            try { const url = syncURL(); if (win.navigator.clipboard?.writeText) await win.navigator.clipboard.writeText(url); else win.prompt('Copy this viewer URL:', url); }
            catch (error) { showError(error); }
        });
        document.getElementById('sessionbutton')?.addEventListener('click', async event => {
            if (sessionScripts) return; event.preventDefault();
            try { await enableSessions(); win.jQuery('#sessionbutton').trigger('click'); } catch (error) { showError(error); }
        });
        function destroy() {
            disposed = true; generation++; sessionUpdate = undefined; controller?.abort(); win.leaveViewerSession?.(); clearPicking();
            win.clearInterval(poll); win.clearTimeout(syncTimer); unsubscribe?.(); state?.destroy();
            win.removeEventListener('popstate', navigate); win.removeEventListener('hashchange', navigate);
        }
        makeState({ open: true });
        app.ready = loadURL(options.url || win.location.href); return app;
    }
    return { parse, serialize, createApp, inferFormat };
});
