/** @jest-environment node */
// DOM interactions use real model selection without a WebGL/native canvas renderer.
jest.mock('canvas', () => ({}));
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!doctype html><html><body></body></html>', { pretendToBeVisual: true });
const globals = ['window', 'document', 'Option', 'FileReader', 'getComputedStyle'];
const original = new Map(globals.map(key => [key, global[key]]));
Object.assign(global, {
    window: dom.window, document: dom.window.document, Option: dom.window.Option,
    FileReader: dom.window.FileReader, getComputedStyle: dom.window.getComputedStyle.bind(dom.window)
});
global.URL.createObjectURL = jest.fn();
dom.window.URL.createObjectURL = jest.fn();
Object.defineProperty(global, '$3Dmol', {
    configurable: true, get: () => dom.window.$3Dmol, set: value => { dom.window.$3Dmol = value; }
});
const $3Dmol = require('../../build/3Dmol.js');
const { StateManager } = require('../../build/3Dmol.ui.js');
const { Event, File } = dom.window;

afterAll(() => {
    dom.window.close();
    delete global.$3Dmol;
    for (const key of globals) {
        if (original.get(key) === undefined) delete global[key];
        else global[key] = original.get(key);
    }
});

function fixture() {
    const container = document.createElement('div');
    container.style.position = 'static';
    Object.defineProperties(container, { clientWidth: { value: 900, configurable: true }, clientHeight: { value: 600, configurable: true } });
    const canvas = document.createElement('canvas');
    container.append(canvas);
    document.body.append(container);
    let sequence = 0;
    const models = [];
    let currentFrame = 0;
    let nextSurface = 0;
    let nextLabel = 0;
    const surfaces = new Map();
    const unrelatedLabel = { id: 'host-label' };
    const labels = new Set([unrelatedLabel]);
    const addAtoms = atoms => {
        const model = new $3Dmol.GLModel(sequence++);
        model.addAtoms(atoms.map((atom, index) => ({ x: index, y: 0, z: 0, bonds: [], ...atom })));
        models.push(model);
        return model;
    };
    addAtoms([
        { elem: 'C', atom: 'CA', chain: '', resn: 'GLY', resi: 0, hetflag: false },
        { elem: 'N', atom: 'N', chain: '', resn: 'GLY', resi: 0, hetflag: false },
        { elem: 'O', atom: 'O', chain: 'A', resn: 'HOH', resi: 4, hetflag: true },
        { elem: 'S', atom: 'SG', chain: 'A', resn: 'CYS', resi: 2, hetflag: false }
    ]);
    addAtoms([
        { elem: 'C', atom: 'C1', chain: 'B', resn: 'LIG', resi: 10, hetflag: true },
        { elem: 'Cl', atom: 'CL', chain: 'B', resn: 'LIG', resi: 10, hetflag: true }
    ]);
    const modelMatches = (model, spec) => spec.model === undefined || (Array.isArray(spec.model) ? spec.model : [spec.model])
        .some(item => item === model || item === model.getID());
    const parse = (data, format, options) => {
        const model = new $3Dmol.GLModel(sequence++);
        model.addMolData(data, format, options);
        models.push(model);
        return model;
    };
    const viewer = {
        container,
        getCanvas: () => canvas,
        getNextModelId: () => sequence,
        getModel: id => id === undefined ? models[models.length - 1] : models.find(model => model.getID() === id),
        getInternalState: () => ({ models: models.map(model => model.getInternalState()) }),
        selectedAtoms: spec => models.filter(model => modelMatches(model, spec)).flatMap(model => model.selectedAtoms(spec)),
        setStyle: (spec, style) => models.filter(model => modelMatches(model, spec)).forEach(model => model.setStyle(spec, style)),
        addStyle: (spec, style) => models.filter(model => modelMatches(model, spec)).forEach(model => model.setStyle(spec, style, true)),
        getNumFrames: () => Math.max(0, ...models.map(model => model.getNumFrames())),
        getFrame: () => currentFrame,
        setFrame: jest.fn(frame => Promise.all(models.map(model => model.setFrame(frame))).then(() => { currentFrame = frame; })),
        isAnimated: () => false,
        stopAnimate: jest.fn(),
        getSurface: id => surfaces.get(id),
        addSurface: jest.fn((type, style, spec, calculationSpec) => {
            const surfid = nextSurface++;
            surfaces.set(surfid, { type, style, spec, calculationSpec });
            return Object.assign(Promise.resolve(surfid), { surfid });
        }),
        removeSurface: jest.fn(id => surfaces.delete(id)),
        setSurfaceMaterialStyle: jest.fn((id, style) => { surfaces.get(id).style = style; }),
        addModelsAsFrames: jest.fn((data, format) => parse(data, format, { multimodel: true, frames: true })),
        addModel: jest.fn((data, format) => parse(data, format, {})),
        removeModel: jest.fn(model => { models.splice(models.indexOf(model), 1); }),
        render: jest.fn(), zoomTo: jest.fn(), removeLabels: jest.fn(list => list.forEach(label => labels.delete(label))),
        addResLabels: jest.fn(() => {
            const label = { id: `residue-${nextLabel++}` };
            labels.add(label);
            return [label];
        }), setBackgroundColor: jest.fn(), setProjection: jest.fn()
    };
    return { viewer, models, container, surfaces, labels, unrelatedLabel };
}

const selectionFields = { Chain: 'chain', 'Residue name': 'resn', 'Residue number': 'resi', Element: 'elem' };
function fieldScope(root, text) {
    if (root.classList.contains('molui') && [...Object.keys(selectionFields), 'Model', 'Atoms', 'Representation'].includes(text)) {
        return root.querySelector('.molui-selection-create') || root;
    }
    return root;
}
function multiselect(root, text) {
    const key = selectionFields[text];
    return key && fieldScope(root, text).querySelector(`.molui-multiselect[data-field="${key}"]`);
}
function field(root, text) {
    const choices = multiselect(root, text);
    if (choices) return choices.querySelector('summary');
    const label = [...fieldScope(root, text).querySelectorAll('label.molui-field')].find(item => item.firstElementChild.textContent === text);
    if (!label) throw new Error(`Missing field: ${text}`);
    return label.querySelector('select, input');
}
function chooseMany(root, text, values) {
    if (text === 'Residue number') {
        button(multiselect(root, text), 'All').click();
        for (const value of values) {
            const match = String(value).match(/^(-?\d+)(?:-(-?\d+))?$/);
            if (!match) throw new Error(`Invalid residue interval: ${value}`);
            addResidueRange(root, Number(match[1]), Number(match[2] ?? match[1]));
        }
        return;
    }
    const valuesAsStrings = values.map(String);
    const checked = [...multiselect(root, text).querySelectorAll('input[type=checkbox][data-value]:checked')].map(input => input.dataset.value);
    for (const value of checked) {
        if (!valuesAsStrings.includes(value)) {
            const input = [...multiselect(root, text).querySelectorAll('input[data-value]')].find(input => input.dataset.value === value);
            toggle(input, false);
        }
    }
    for (const value of valuesAsStrings) {
        const input = [...multiselect(root, text).querySelectorAll('input[data-value]')].find(input => input.dataset.value === value);
        if (!input) throw new Error(`Missing ${text} choice: ${value}`);
        if (!input.checked) toggle(input, true);
    }
}
function choose(root, text, value) {
    if (multiselect(root, text)) return chooseMany(root, text, value === '*' ? [] : [value]);
    const input = field(root, text);
    input.value = value;
    input.dispatchEvent(new Event('change', { bubbles: true }));
}
function button(root, text) {
    return [...root.querySelectorAll('button')].find(item => item.textContent === text);
}
function options(root, label) {
    const choices = multiselect(root, label);
    // Residue numbers expose the available interval rather than discrete choices.
    if (label === 'Residue number') return [...new Set([choices.dataset.min, choices.dataset.max].filter(value => value !== undefined))];
    return choices ? [...choices.querySelectorAll('input[type=checkbox][data-value]')].map(input => input.dataset.value)
        : [...field(root, label).options].map(option => option.value);
}
function checkedOptions(root, label) {
    if (label === 'Residue number') return [...multiselect(root, label).querySelectorAll('.molui-residue-chip')]
        .map(chip => chip.dataset.from === chip.dataset.to ? chip.dataset.from : `${chip.dataset.from}-${chip.dataset.to}`);
    return [...multiselect(root, label).querySelectorAll('input[type=checkbox][data-value]:checked')].map(input => input.dataset.value);
}
function addResidueRange(root, from, to) {
    const controls = () => multiselect(root, 'Residue number');
    setInput(controls().querySelector('[data-focus$=":resi:from"]'), String(from));
    setInput(controls().querySelector('[data-focus$=":resi:to"]'), String(to));
    button(controls(), 'Add range').click();
}
function createCount(root) { return root.querySelector('.molui-selection-create .molui-count'); }
function editStyle(root, name, value) {
    const input = root.querySelector(`.molui-style [data-focus$=":${name}"]`);
    if (!input) throw new Error(`Missing style control: ${name}`);
    input.value = String(value);
    input.dispatchEvent(new Event('input', { bubbles: true }));
    input.dispatchEvent(new Event('change', { bubbles: true }));
}
function selectionCard(root) { return root.querySelector('.molui-selection-header').closest('.molui-card'); }
function cardFor(root, id) { return root.querySelector(`.molui-selection-card[data-selection-id="${id}"]`); }
function checkbox(root, text) {
    const label = [...root.querySelectorAll('label')].find(item => item.textContent.trim() === text);
    if (!label) throw new Error(`Missing checkbox: ${text}`);
    return label.querySelector('input[type=checkbox]');
}
function toggle(input, checked) {
    input.checked = checked;
    input.dispatchEvent(new Event('change', { bubbles: true }));
}
function setInput(input, value) {
    input.value = value;
    input.dispatchEvent(new Event('input', { bubbles: true }));
    input.dispatchEvent(new Event('change', { bubbles: true }));
}
function addRepresentation(root, representation) {
    const row = root.querySelector('.molui-add-representation');
    row.querySelector('select').value = representation;
    button(row, 'Add').click();
}
async function settle() { for (let index = 0; index < 6; index++) await Promise.resolve(); }

describe('guided molecular viewer UI', () => {
    let states;
    let fetchOriginal;
    const mount = (config = {}) => {
        const data = fixture();
        const state = new StateManager(data.viewer, { open: true, ...config });
        states.push(state);
        return { ...data, state, root: state.ui.element };
    };

    beforeEach(() => {
        jest.useFakeTimers();
        states = [];
        fetchOriginal = global.fetch;
    });
    afterEach(() => {
        states.forEach(state => state.destroy());
        document.body.replaceChildren();
        global.fetch = fetchOriginal;
        jest.useRealTimers();
    });

    test('multiple viewers have independent controls, selections, and panel IDs', () => {
        const first = mount();
        const second = mount();
        const panels = [...document.querySelectorAll('.molui-panel')];
        expect(panels).toHaveLength(2);
        expect(new Set(panels.map(panel => panel.id)).size).toBe(2);
        for (const { root } of [first, second]) {
            const toggle = root.querySelector('.molui-toggle');
            expect(root.querySelector(`#${toggle.getAttribute('aria-controls')}`)).not.toBeNull();
            expect(toggle.textContent).toBe('☰');
            expect(toggle.getAttribute('aria-label')).toBe('Toggle controls');
            expect(toggle.getAttribute('aria-expanded')).toBe('true');
            const header = root.querySelector('.molui-header');
            expect(header.children).toHaveLength(1);
            expect(header.textContent).toBe('×');
            expect(header.querySelector('button').getAttribute('aria-label')).toBe('Close controls');
        }
        choose(first.root, 'Chain', 'B');
        button(first.root, 'Add selection').click();
        expect(first.state.selections.size).toBe(1);
        expect(second.state.selections.size).toBe(0);
        expect(createCount(second.root).textContent).toBe('6 matching atoms');
        first.root.querySelector('button[aria-label="Close controls"]').click();
        expect(first.root.classList.contains('molui-closed')).toBe(true);
        expect(first.root.querySelector('.molui-toggle').getAttribute('aria-expanded')).toBe('false');
        expect(second.root.classList.contains('molui-closed')).toBe(false);
    });

    test('embeddings can restore view and panel settings without manipulating control DOM', () => {
        const onViewChange = jest.fn(), onOpenChange = jest.fn();
        const { root, state, viewer } = mount({
            view: { background: '#000066', projection: 'orthographic' }, onViewChange, onOpenChange
        });
        expect(field(root, 'Background').value).toBe('#000066');
        expect(field(root, 'Projection').value).toBe('orthographic');
        expect(onOpenChange).not.toHaveBeenCalled();
        state.setViewSettings({ background: '#112233', projection: 'perspective' });
        expect(field(root, 'Background').value).toBe('#112233');
        choose(root, 'Projection', 'orthographic');
        expect(viewer.setProjection).toHaveBeenLastCalledWith('orthographic');
        expect(onViewChange).toHaveBeenLastCalledWith({ background: '#112233', backgroundAlpha: 1, projection: 'orthographic' });
        state.ui.setOpen(false);
        expect(state.ui.isOpen).toBe(false);
        expect(onOpenChange).toHaveBeenLastCalledWith(false);
        root.querySelector('.molui-toggle').click();
        expect(state.ui.isOpen).toBe(true);
        expect(onOpenChange).toHaveBeenLastCalledWith(true);
    });

    test.each([
        [{}, true, true],
        [{ edit: false }, false, true],
        [{ trajectory: false }, true, false],
        [{ edit: false, trajectory: false }, false, false]
    ])('control availability is independent for %j', (config, edit, trajectory) => {
        const { root, state, viewer } = mount(config);
        viewer.addModelsAsFrames('1\nfirst\nC 0 0 0\n1\nsecond\nC 1 0 0\n', 'xyz');
        state.refresh();
        expect(state.getControls()).toEqual({ edit, trajectory });
        expect(root.querySelector('.molui-toggle').hidden).toBe(!edit);
        expect(root.querySelector('.molui-panel').hidden).toBe(!edit);
        expect(root.querySelector('.molui-trajectory').hidden).toBe(!trajectory);
        expect(root.querySelector('.molui-trajectory-toggle').hidden).toBe(true);
        expect(root.classList.contains('molui-no-edit')).toBe(!edit);
        expect(root.classList.contains('molui-has-trajectory')).toBe(trajectory);
        state.ui.setOpen(true);
        expect(state.ui.isOpen).toBe(edit);
        state.ui.setTrajectoryMinimized(true);
        expect(root.querySelector('.molui-trajectory').hidden).toBe(true);
        expect(root.querySelector('.molui-trajectory-toggle').hidden).toBe(!trajectory);
        expect(root.classList.contains('molui-has-trajectory')).toBe(false);
    });

    test('availability can change without recreating controls or losing appearance state', () => {
        const onOpenChange = jest.fn();
        const { root, state, viewer } = mount({ onOpenChange });
        viewer.addModelsAsFrames('1\nfirst\nC 0 0 0\n1\nsecond\nC 1 0 0\n', 'xyz');
        state.refresh();
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][0];
        state.ui.setTrajectoryMinimized(true);
        state.setControls({ edit: false });
        expect(state.getControls()).toEqual({ edit: false, trajectory: true });
        expect(state.ui.isOpen).toBe(false);
        expect(onOpenChange).toHaveBeenLastCalledWith(false);
        expect(root.querySelector('.molui-toggle').hidden).toBe(true);
        expect(root.querySelector('.molui-trajectory-toggle').hidden).toBe(false);
        state.setControls({ trajectory: false });
        jest.advanceTimersByTime(1000);
        state.refresh();
        expect(root.querySelector('.molui-trajectory').hidden).toBe(true);
        expect(root.querySelector('.molui-trajectory-toggle').hidden).toBe(true);
        state.setControls({ edit: true, trajectory: true });
        expect(state.ui.isTrajectoryMinimized).toBe(true);
        root.querySelector('.molui-toggle').click();
        root.querySelector('.molui-trajectory-toggle').click();
        expect(state.ui.isOpen).toBe(true);
        expect(state.ui.isTrajectoryMinimized).toBe(false);
        expect(root.querySelector('.molui-trajectory').hidden).toBe(false);
        expect(root.classList.contains('molui-has-trajectory')).toBe(true);
        expect(state.ui.element).toBe(root);
        expect([...state.selections.values()]).toEqual([selection]);
    });

    test('trajectory-only controls follow data added and removed by the host while minimized', () => {
        const { root, state, viewer } = mount({ edit: false });
        const panel = root.querySelector('.molui-trajectory');
        const launcher = root.querySelector('.molui-trajectory-toggle');
        state.ui.setTrajectoryMinimized(true);
        expect(panel.hidden).toBe(true);
        expect(launcher.hidden).toBe(true);
        const model = viewer.addModelsAsFrames('1\nfirst\nC 0 0 0\n1\nsecond\nC 1 0 0\n', 'xyz');
        jest.advanceTimersByTime(500);
        expect(panel.hidden).toBe(true);
        expect(launcher.hidden).toBe(false);
        viewer.removeModel(model);
        jest.advanceTimersByTime(500);
        expect(panel.hidden).toBe(true);
        expect(launcher.hidden).toBe(true);
        expect(root.classList.contains('molui-has-trajectory')).toBe(false);
        expect(root.querySelector('.molui-toggle').hidden).toBe(true);
    });

    test('restoring playback in a tiny viewer closes the sidebar to reveal its controls', () => {
        const { root, state, viewer, container } = mount();
        Object.defineProperties(container, { clientWidth: { value: 360 }, clientHeight: { value: 280 } });
        viewer.addModelsAsFrames('1\nfirst\nC 0 0 0\n1\nsecond\nC 1 0 0\n', 'xyz');
        state.refresh();
        state.ui.resize();
        state.ui.setTrajectoryMinimized(true);
        expect(state.ui.isOpen).toBe(true);
        root.querySelector('.molui-trajectory-toggle').click();
        expect(state.ui.isOpen).toBe(false);
        expect(root.querySelector('.molui-trajectory').hidden).toBe(false);
        expect(document.activeElement.getAttribute('aria-label')).toBe('Minimize trajectory controls');
        state.ui.setTrajectoryMinimized(true);
        state.ui.setOpen(true);
        state.ui.setTrajectoryMinimized(false);
        expect(state.ui.isOpen).toBe(false);
    });

    test('dependent choices preserve blank chains and residue zero without broadening a draft', () => {
        const { root, state } = mount();
        expect(options(root, 'Chain')).toEqual(['', 'A', 'B']);
        choose(root, 'Chain', '');
        expect(options(root, 'Residue name')).toEqual(['GLY']);
        expect(options(root, 'Residue number')).toEqual(['0']);
        expect(options(root, 'Element')).toEqual(['C', 'N']);
        choose(root, 'Residue number', '0');
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][0];
        expect(selection.spec).toEqual({ chain: '', resi: 0 });
        expect(createCount(root).textContent).toBe('2 matching atoms');
        choose(root, 'Chain', 'B');
        expect(checkedOptions(root, 'Residue number')).toEqual(['0']);
        expect(options(root, 'Residue number')).toEqual(['10']);
        expect(createCount(root).textContent).toBe('0 matching atoms');
        choose(root, 'Residue number', '*');
        expect(createCount(root).textContent).toBe('2 matching atoms');
    });

    test('model filtering restricts matching atoms and representation choices', () => {
        const { root, state } = mount();
        expect(options(root, 'Representation')).toContain('cartoon');
        choose(root, 'Model', '1');
        expect(options(root, 'Chain')).toEqual(['B']);
        expect(createCount(root).textContent).toBe('2 matching atoms');
        expect(options(root, 'Representation')).not.toContain('cartoon');
        button(root, 'Add selection').click();
        expect([...state.selections.values()][0].spec).toEqual({ model: 1 });
    });

    test('empty matching groups disable appearance actions and dependent selectors', () => {
        const { root, state } = mount();
        choose(root, 'Model', '1');
        choose(root, 'Atoms', 'water');
        expect(createCount(root).textContent).toBe('0 matching atoms');
        expect(button(root, 'Add selection').disabled).toBe(true);
        expect(field(root, 'Representation').disabled).toBe(true);
        expect(options(root, 'Chain')).toEqual([]);
        button(root, 'Add selection').click();
        expect(state.selections.size).toBe(0);
        choose(root, 'Atoms', 'ligand');
        expect(button(root, 'Add selection').disabled).toBe(false);
    });

    test('multiple values use OR within each field and AND across selection fields', () => {
        const { root, state, viewer } = mount();
        chooseMany(root, 'Chain', ['', 'A']);
        expect(checkedOptions(root, 'Chain')).toEqual(['', 'A']);
        expect(createCount(root).textContent).toBe('4 matching atoms');
        chooseMany(root, 'Residue name', ['GLY', 'CYS']);
        chooseMany(root, 'Residue number', [0, 2]);
        chooseMany(root, 'Element', ['C', 'S']);
        expect(createCount(root).textContent).toBe('2 matching atoms');
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][0];
        expect(selection.spec.chain).toEqual(expect.arrayContaining(['', 'A']));
        expect(selection.spec.resn).toEqual(expect.arrayContaining(['GLY', 'CYS']));
        expect(selection.spec.resi).toEqual(expect.arrayContaining([0, 2]));
        expect(selection.spec.elem).toEqual(expect.arrayContaining(['C', 'S']));
        expect(viewer.selectedAtoms(selection.spec).map(atom => `${atom.chain}:${atom.resi}:${atom.elem}`)).toEqual([':0:C', 'A:2:S']);
        const summary = selectionCard(root).querySelector('.molui-selection-summary').textContent;
        expect(summary).toMatch(/chain/i);
        expect(summary).toContain('GLY');
        expect(summary).toContain('CYS');
        expect(summary).toContain('2 atoms');
    });

    test('residue intervals stay drafts until added and combine with disjoint ranges', () => {
        const { root, state, viewer } = mount();
        const controls = () => multiselect(root, 'Residue number');
        controls().open = true;
        setInput(controls().querySelector('[data-focus$=":resi:from"]'), '0');
        setInput(controls().querySelector('[data-focus$=":resi:to"]'), '2');
        expect(createCount(root).textContent).toBe('6 matching atoms');
        button(controls(), 'Add range').click();
        expect(createCount(root).textContent).toBe('3 matching atoms');
        expect(checkedOptions(root, 'Residue number')).toEqual(['0-2']);
        addResidueRange(root, 10, 10);
        expect(createCount(root).textContent).toBe('5 matching atoms');
        expect(checkedOptions(root, 'Residue number')).toEqual(['0-2', '10']);
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][0];
        expect(selection.spec.resi).toEqual(['0-2', 10]);
        expect(viewer.selectedAtoms(selection.spec).map(atom => atom.resi)).toEqual([0, 0, 2, 10, 10]);
        expect(selectionCard(root).querySelector('.molui-selection-summary').textContent).toContain('5 atoms');
        button(controls(), 'All').click();
        expect(createCount(root).textContent).toBe('6 matching atoms');
        expect(selection.spec.resi).toEqual(['0-2', 10]);
    });

    test('existing selections precede the new-selection form and expose criteria summaries and Edit', () => {
        const { root } = mount();
        choose(root, 'Chain', 'B');
        choose(root, 'Representation', 'stick');
        button(root, 'Add selection').click();
        const list = root.querySelector('.molui-selection-list');
        const create = root.querySelector('.molui-selection-create');
        expect(list.compareDocumentPosition(create) & dom.window.Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
        const card = selectionCard(root);
        expect(card.querySelector('.molui-selection-summary').textContent).toMatch(/chain.*B.*2 atoms/i);
        expect(button(card.querySelector('.molui-selection-summary-row'), 'Edit')).toBeDefined();
        expect(button(card.querySelector('.molui-add-representation'), 'Add')).toBeDefined();
        expect(button(card, 'Add representation')).toBeUndefined();
    });

    test('individual selections collapse without changing appearance and retain their state across refreshes', () => {
        const { root, state, viewer } = mount();
        const first = state.createSelectionAndStyle({ chain: 'B' }, { stick: {} });
        const second = state.createSelectionAndStyle({ chain: 'A' }, { sphere: {} });
        const card = () => cardFor(root, first);
        const collapse = () => card().querySelector('.molui-selection-collapse');
        const content = () => card().querySelector('.molui-selection-content');
        expect(collapse().getAttribute('aria-expanded')).toBe('true');
        expect(document.getElementById(collapse().getAttribute('aria-controls'))).toBe(content());
        expect(content().hidden).toBe(false);
        card().querySelector('.molui-geometry').open = true;
        const appearances = viewer.selectedAtoms({}).map(atom => JSON.stringify(atom.style));
        const renders = viewer.render.mock.calls.length;
        collapse().focus();
        collapse().click();
        expect(content().hidden).toBe(true);
        expect(collapse().getAttribute('aria-expanded')).toBe('false');
        expect(card().querySelector('.molui-selection-name').closest('[hidden]')).toBeNull();
        expect(card().querySelector('.molui-selection-summary').closest('[hidden]')).toBeNull();
        expect(cardFor(root, second).querySelector('.molui-selection-content').hidden).toBe(false);
        expect(viewer.render).toHaveBeenCalledTimes(renders);
        expect(viewer.selectedAtoms({}).map(atom => JSON.stringify(atom.style))).toEqual(appearances);

        setInput(card().querySelector('.molui-selection-name'), 'Binding pocket');
        expect(collapse().getAttribute('aria-label')).toBe('Expand selection Binding pocket');
        toggle(card().querySelector('.molui-selection-visibility input'), false);
        expect(state.selections.get(first).hidden).toBe(true);
        expect(content().hidden).toBe(true);
        collapse().focus();
        state.refresh();
        expect(content().hidden).toBe(true);
        expect(document.activeElement).toBe(collapse());
        expect(card().querySelector('.molui-selection-summary').textContent).toMatch(/chain.*B.*2 atoms/i);
        collapse().click();
        expect(content().hidden).toBe(false);
        expect(collapse().getAttribute('aria-label')).toBe('Collapse selection Binding pocket');
        expect(card().querySelector('.molui-geometry').open).toBe(true);
        expect(state.selections.get(first).hidden).toBe(true);
    });

    test('collapsed selections hide Edit and retain unsaved drafts when reopened', () => {
        const { root, state } = mount();
        const id = state.createSelectionAndStyle({ chain: 'B' }, { stick: {} });
        const card = () => cardFor(root, id);
        const collapse = () => card().querySelector('.molui-selection-collapse');
        collapse().click();
        expect(button(card(), 'Edit').hidden).toBe(true);
        collapse().click();
        expect(button(card(), 'Edit').hidden).toBe(false);
        button(card(), 'Edit').click();
        expect(card().querySelector('.molui-selection-content').hidden).toBe(false);
        expect(collapse().getAttribute('aria-expanded')).toBe('true');
        chooseMany(card().querySelector('.molui-selection-edit'), 'Chain', ['', 'A']);
        collapse().focus();
        collapse().click();
        state.refresh();
        expect(card().querySelector('.molui-selection-content').hidden).toBe(true);
        expect(card().querySelector('.molui-selection-summary').textContent).toMatch(/chain.*B.*2 atoms/i);
        expect(state.selections.get(id).spec).toEqual({ chain: 'B' });
        expect(button(card(), 'Edit').hidden).toBe(true);
        collapse().click();
        const edit = card().querySelector('.molui-selection-edit');
        expect(card().querySelector('.molui-selection-content').hidden).toBe(false);
        expect(checkedOptions(edit, 'Chain')).toEqual(['', 'A']);
        expect(document.activeElement).toBe(collapse());
        button(edit, 'Save selection').click();
        expect(state.selections.get(id).spec.chain).toEqual(['', 'A']);
        expect(card().querySelector('.molui-selection-edit')).toBeNull();
        expect(button(card(), 'Edit').hidden).toBe(false);
        expect(card().querySelector('.molui-selection-summary').textContent).toContain('4 atoms');
    });

    test('selection edits are drafts until saved and retain names, representations, and label settings', () => {
        const { root, state, viewer, labels, unrelatedLabel } = mount();
        choose(root, 'Chain', 'B');
        choose(root, 'Representation', 'stick');
        button(root, 'Add selection').click();
        addRepresentation(root, 'sphere');
        setInput(selectionCard(root).querySelector('.molui-selection-name'), 'Ligand');
        toggle(checkbox(selectionCard(root), 'Residue labels'), true);
        const selection = [...state.selections.values()][0];
        const originalSpec = { ...selection.spec };
        const originalStyles = [...selection.styles].map(([id, style]) => [id, JSON.parse(JSON.stringify(style))]);
        button(selectionCard(root), 'Edit').click();
        let edit = root.querySelector('.molui-selection-edit');
        chooseMany(edit, 'Chain', ['', 'A']);
        expect(selection.spec).toEqual(originalSpec);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        expect(labels.size).toBe(2);
        button(edit, 'Cancel').click();
        expect(root.querySelector('.molui-selection-edit')).toBeNull();
        expect(selection.spec).toEqual(originalSpec);
        expect([...selection.styles]).toEqual(originalStyles);
        expect(selection.name).toBe('Ligand');
        button(selectionCard(root), 'Edit').click();
        edit = root.querySelector('.molui-selection-edit');
        chooseMany(edit, 'Chain', ['', 'A']);
        chooseMany(edit, 'Element', ['C', 'N']);
        state.refresh();
        edit = root.querySelector('.molui-selection-edit');
        expect(checkedOptions(edit, 'Chain')).toEqual(['', 'A']);
        expect(checkedOptions(edit, 'Element')).toEqual(['C', 'N']);
        expect(edit.querySelector('.molui-count').textContent).toContain('2');
        button(edit, 'Save selection').click();
        expect(root.querySelector('.molui-selection-edit')).toBeNull();
        expect(selection.spec.chain).toEqual(expect.arrayContaining(['', 'A']));
        expect(selection.spec.elem).toEqual(expect.arrayContaining(['C', 'N']));
        expect([...selection.styles]).toEqual(originalStyles);
        expect(selection.name).toBe('Ligand');
        expect(selection.labelsEnabled).toBe(true);
        expect(checkbox(selectionCard(root), 'Residue labels').checked).toBe(true);
        expect(viewer.addResLabels).toHaveBeenLastCalledWith(selection.spec, expect.any(Object));
        expect(labels.size).toBe(2);
        expect(labels.has(unrelatedLabel)).toBe(true);
        expect(selectionCard(root).querySelector('.molui-selection-summary').textContent).toContain('2 atoms');
    });

    test('editing guided fields preserves imported model lists and custom predicates', () => {
        const { root, state, viewer } = mount();
        const predicate = atom => atom.elem !== 'S';
        const spec = { model: [0, 1], chain: ['', 'A'], not: { resn: 'HOH' }, predicate };
        const id = state.createSelectionAndStyle(spec, { sphere: { color: 'red', scale: 0.4 } });
        const selection = state.selections.get(id);
        button(selectionCard(root), 'Edit').click();
        const edit = root.querySelector('.molui-selection-edit');
        choose(edit, 'Chain', '');
        choose(edit, 'Element', 'C');
        button(edit, 'Save selection').click();
        expect(selection.spec).toMatchObject({ model: [0, 1], chain: '', elem: 'C', not: { resn: 'HOH' } });
        expect(selection.spec.predicate).toBe(predicate);
        expect(viewer.selectedAtoms(selection.spec)).toHaveLength(1);
        expect([...selection.styles.values()][0].spec).toEqual({ sphere: { color: 'red', scale: 0.4 } });
    });

    test('failed remote downloads retain the existing scene and allow retry', async () => {
        const { root, models, viewer } = mount();
        const originalModels = [...models];
        global.fetch = jest.fn().mockResolvedValue({ ok: false, status: 404 });
        choose(root, 'Source', 'pdb');
        field(root, 'PDB ID').value = '1CRN';
        root.querySelector('.molui-body').scrollTop = 200;
        root.querySelector('form').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
        await settle();
        expect(models).toEqual(originalModels);
        expect(viewer.removeModel).not.toHaveBeenCalled();
        expect(viewer.addModelsAsFrames).not.toHaveBeenCalled();
        expect(root.querySelector('[role="status"]').textContent).toContain('HTTP 404');
        expect(root.querySelector('.molui-body').scrollTop).toBe(0);
        expect(button(root, 'Load structure').disabled).toBe(false);
        expect(root.querySelector('form').hasAttribute('aria-busy')).toBe(false);
    });

    test('invalid local data removes only the unsuccessful new model', async () => {
        const { root, models, viewer } = mount();
        const originalModels = [...models];
        // FileReader completion is isolated from the fake timer clock.
        const read = jest.spyOn(FileReader.prototype, 'readAsText').mockImplementation(function () {
            Object.defineProperty(this, 'result', { configurable: true, value: 'This is not a structure.' });
            this.onload(new Event('load'));
        });
        try {
            Object.defineProperty(field(root, 'Structure file'), 'files', {
                configurable: true, value: [new File(['invalid'], 'invalid.xyz')]
            });
            root.querySelector('form').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
            await settle();
            expect(models).toEqual(originalModels);
            expect(viewer.removeModel).toHaveBeenCalledTimes(1);
            expect(viewer.removeModel.mock.calls[0][0]).not.toBe(originalModels[0]);
            expect(root.querySelector('[role="status"]').textContent).toContain('No atoms were found');
            expect(button(root, 'Load structure').disabled).toBe(false);
        } finally { read.mockRestore(); }
    });

    test('loading a local trajectory adds video controls without autoplay or replacing models', async () => {
        const events = [];
        const onModelLoaded = jest.fn(() => events.push('loaded'));
        const { root, models, viewer, state } = mount({ onModelLoaded });
        state.subscribe(() => events.push('changed'));
        const originalModels = [...models];
        const xyz = '2\nfirst\nC 0 0 0\nO 1 0 0\n2\nsecond\nC 0 1 0\nO 1 1 0\n';
        const read = jest.spyOn(FileReader.prototype, 'readAsText').mockImplementation(function () {
            Object.defineProperty(this, 'result', { configurable: true, value: xyz });
            this.onload(new Event('load'));
        });
        try {
            Object.defineProperty(field(root, 'Structure file'), 'files', {
                configurable: true, value: [new File([xyz], 'motion.xyz')]
            });
            root.querySelector('form').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
            await settle();
            expect(models).toHaveLength(3);
            expect(models.slice(0, 2)).toEqual(originalModels);
            expect(models[2].getNumFrames()).toBe(2);
            expect(onModelLoaded).toHaveBeenCalledWith(models[2], {
                kind: 'file', value: 'motion.xyz', name: 'motion.xyz', format: 'xyz', data: xyz, frames: true
            });
            expect(events[0]).toBe('loaded');
            expect(root.querySelector('.molui-trajectory').hidden).toBe(false);
            expect(root.querySelector('.molui-frame-count').textContent).toBe('Frame 1 / 2');
            expect(root.classList.contains('molui-has-trajectory')).toBe(true);
            jest.advanceTimersByTime(1000);
            await settle();
            expect(viewer.setFrame).not.toHaveBeenCalled();
        } finally { read.mockRestore(); }
    });

    test('keyboard focus survives dependent filters and appearance edits', () => {
        const { root } = mount();
        const chainChoice = [...multiselect(root, 'Chain').querySelectorAll('input[data-value]')].find(input => input.dataset.value === 'B');
        chainChoice.focus();
        toggle(chainChoice, true);
        expect(multiselect(root, 'Chain').contains(document.activeElement)).toBe(true);
        expect(document.activeElement.dataset.value).toBe('B');
        button(root, 'Add selection').click();
        let editor = root.querySelector('.molui-style');
        field(editor, 'Color by').focus();
        choose(editor, 'Color by', 'custom');
        editor = root.querySelector('.molui-style');
        expect(document.activeElement === field(editor, 'Color by')).toBe(true);
        expect(field(editor, 'Solid color').closest('label').hidden).toBe(false);
        field(editor, 'Representation').focus();
        choose(editor, 'Representation', 'sphere');
        editor = root.querySelector('.molui-style');
        expect(document.activeElement === field(editor, 'Representation')).toBe(true);
    });

    function importedAppearance() {
        const data = mount();
        const spec = {
            cartoon: { color: 'red', arrows: true, thickness: 0.4, tubes: true },
            stick: { colorscheme: { C: 0x112233, O: 0x445566 }, radius: 0.4, singleBonds: true },
            sphere: { color: 'navy', scale: 0.6 }
        };
        const id = data.state.createSelectionAndStyle({}, spec);
        return {
            ...data, spec,
            current: () => [...data.state.selections.get(id).styles.values()][0].spec
        };
    }

    test.each(['stick', 'ballstick', 'sphere', 'line', 'cross', 'cartoon'])('%s keeps opacity and geometry inside Settings with stable focus', kind => {
        const { root, state } = mount();
        expect(options(root, 'Representation')).toContain(kind);
        choose(root, 'Representation', kind);
        button(root, 'Add selection').click();
        const editor = () => root.querySelector('.molui-style');
        const settings = () => editor().querySelector('details.molui-geometry');
        const opacity = () => editor().querySelector('[data-focus$=":opacity"]');
        const selection = [...state.selections.values()][0];
        const style = [...selection.styles.values()][0];
        const before = JSON.parse(JSON.stringify(style.spec));
        const geometryKey = settings().dataset.geometry;
        expect(settings().firstElementChild.textContent).toBe('Settings');
        expect(settings().open).toBe(false);
        expect(settings().contains(opacity())).toBe(true);
        expect(editor().querySelector('.molui-style-grid').contains(opacity())).toBe(false);
        const radius = editor().querySelector('[data-focus$=":radius"]');
        expect(settings().contains(radius)).toBe(true);
        expect(radius.closest('label').hidden).toBe(!['stick', 'ballstick'].includes(kind));
        expect(settings().contains(field(editor(), 'Representation'))).toBe(false);
        expect(settings().contains(field(editor(), 'Color by'))).toBe(false);
        settings().open = true;
        opacity().focus();
        editStyle(root, 'opacity', 0.45);
        expect(document.activeElement).toBe(opacity());
        expect(settings().open).toBe(true);
        expect(settings().dataset.geometry).toBe(geometryKey);
        const expected = Object.fromEntries(Object.entries(before).map(([name, properties]) => [name, { ...properties, opacity: 0.45 }]));
        expect([...selection.styles.values()][0].spec).toEqual(expected);
        settings().open = false;
        state.refresh();
        expect(settings().open).toBe(false);
        expect(settings().dataset.geometry).toBe(geometryKey);
    });

    test('Settings stays open across representation changes including lines and sticks', () => {
        const { root } = mount();
        choose(root, 'Representation', 'sphere');
        button(root, 'Add selection').click();
        const settings = () => root.querySelector('.molui-style .molui-geometry');
        const geometryKey = settings().dataset.geometry;
        settings().open = true;
        for (const kind of ['line', 'stick', 'cartoon']) {
            const editor = root.querySelector('.molui-style');
            field(editor, 'Representation').focus();
            choose(editor, 'Representation', kind);
            expect(settings().open).toBe(true);
            expect(settings().dataset.geometry).toBe(geometryKey);
            expect(settings().firstElementChild.textContent).toBe('Settings');
            expect(document.activeElement).toBe(field(root.querySelector('.molui-style'), 'Representation'));
            expect(settings().querySelector('[data-focus$=":opacity"]')).not.toBeNull();
        }
    });

    test('opacity edits preserve imported combined representations and their custom options', () => {
        const { root, spec, current } = importedAppearance();
        const settings = () => root.querySelector('.molui-geometry');
        expect(settings().firstElementChild.textContent).toBe('Settings');
        settings().open = true;
        root.querySelector('[data-focus$=":opacity"]').focus();
        editStyle(root, 'opacity', 0.45);
        expect(current()).toEqual({
            cartoon: { ...spec.cartoon, opacity: 0.45 },
            stick: { ...spec.stick, opacity: 0.45 },
            sphere: { ...spec.sphere, opacity: 0.45 }
        });
        expect(document.activeElement).toBe(settings().querySelector('[data-focus$=":opacity"]'));
        expect(settings().open).toBe(true);
    });

    test('bond radius edits only change the stick radius in an imported combined appearance', () => {
        const { root, spec, current } = importedAppearance();
        const radius = root.querySelector('.molui-style [data-focus$=":radius"]');
        expect(radius.closest('label').hidden).toBe(false);
        expect(Number(radius.value)).toBe(0.4);
        editStyle(root, 'radius', 0.55);
        expect(current()).toEqual({ ...spec, stick: { ...spec.stick, radius: 0.55 } });
    });

    test('spheres start at atomic radii and can switch between scaled atoms and a fixed radius', () => {
        const { root, state } = mount();
        expect([...field(root, 'Representation').options].find(option => option.value === 'sphere').text).toBe('Spheres');
        choose(root, 'Representation', 'sphere');
        button(root, 'Add selection').click();
        const current = () => [...state.selections.values()][0].styles.values().next().value.spec.sphere;
        expect(current()).toEqual({ colorscheme: 'default', opacity: 1 });
        expect(root.querySelector('.molui-geometry').open).toBe(false);
        root.querySelector('.molui-geometry').open = true;
        editStyle(root, 'sphere-scale', 0.6);
        expect(current().scale).toBe(0.6);
        expect(current().radius).toBeUndefined();
        editStyle(root, 'sphere-radius-mode', 'fixed');
        editStyle(root, 'sphere-radius', 2.25);
        expect(current()).toEqual({ colorscheme: 'default', opacity: 1, radius: 2.25 });
        expect(root.querySelector('.molui-geometry').open).toBe(true);
        expect(root.querySelector('[data-focus$=":sphere-scale"]').closest('label').hidden).toBe(true);
        editStyle(root, 'sphere-radius-mode', 'atomic');
        expect(current()).toEqual({ colorscheme: 'default', opacity: 1 });
        expect(root.querySelector('[data-focus$=":sphere-radius"]').closest('label').hidden).toBe(true);
    });

    test('sphere and cross radii edit independently and use effective imported sizes', () => {
        const { root, state } = mount();
        const spec = {
            sphere: { color: 'red', radius: 2, scale: 0.5, opacity: 0.6 },
            cross: { colorscheme: 'greenCarbon', radius: 3, scale: 0.2, opacity: 0.7 }
        };
        const id = state.createSelectionAndStyle({}, spec);
        const current = () => [...state.selections.get(id).styles.values()][0].spec;
        expect(Number(root.querySelector('[data-focus$=":sphere-radius"]').value)).toBe(1);
        expect(Number(root.querySelector('[data-focus$=":cross-radius"]').value)).toBeCloseTo(0.6);
        root.querySelector('.molui-geometry').open = true;
        root.querySelector('[data-focus$=":sphere-radius"]').focus();
        editStyle(root, 'sphere-radius', 0.8);
        expect(current()).toEqual({ ...spec, sphere: { color: 'red', radius: 0.8, opacity: 0.6 } });
        expect(document.activeElement).toBe(root.querySelector('[data-focus$=":sphere-radius"]'));
        expect(root.querySelector('.molui-geometry').open).toBe(true);
        editStyle(root, 'cross-radius', 0.75);
        expect(current().cross).toEqual({ colorscheme: 'greenCarbon', radius: 0.75, opacity: 0.7 });
        expect(current().sphere.radius).toBe(0.8);
        editStyle(root, 'cross-radius-mode', 'atomic');
        editStyle(root, 'cross-scale', 1.25);
        expect(current().cross).toEqual({ colorscheme: 'greenCarbon', scale: 1.25, opacity: 0.7 });
    });

    test('ball and stick exposes sphere scale independently from bond radius', () => {
        const { root, state } = mount();
        choose(root, 'Representation', 'ballstick');
        button(root, 'Add selection').click();
        const current = () => [...state.selections.values()][0].styles.values().next().value.spec;
        expect(Number(root.querySelector('[data-focus$=":sphere-scale"]').value)).toBe(0.3);
        editStyle(root, 'sphere-scale', 0.5);
        expect(current().sphere.scale).toBe(0.5);
        expect(current().stick.radius).toBe(0.25);
        editStyle(root, 'radius', 0.55);
        expect(current().stick.radius).toBe(0.55);
        expect(current().sphere.scale).toBe(0.5);
    });

    test.each(['default', 'oval', 'parabola', 'edged'])('cartoon shape %s preserves other imported representations and options', shape => {
        const { root, spec, current } = importedAppearance();
        editStyle(root, 'cartoon-style', shape);
        expect(current()).toEqual({ ...spec, cartoon: { ...spec.cartoon, style: shape } });
        editStyle(root, 'cartoon-thickness', 0.65);
        expect(current().cartoon).toEqual({ ...spec.cartoon, style: shape, thickness: 0.65 });
    });

    test('cartoon controls follow trace, ribbon and width semantics without losing unrelated settings', () => {
        const { root, state } = mount();
        const cartoon = { color: 'blue', arrows: true, tubes: true, width: 1.8, gapcutoff: 4, opacity: 0.8 };
        const id = state.createSelectionAndStyle({}, { cartoon });
        const current = () => [...state.selections.get(id).styles.values()][0].spec.cartoon;
        const control = name => root.querySelector(`[data-focus$=":cartoon-${name}"]`);
        expect(control('arrows').closest('label').hidden).toBe(true);
        editStyle(root, 'cartoon-width-mode', 'auto');
        expect(current().width).toBeUndefined();
        expect(control('arrows').closest('label').hidden).toBe(false);
        toggle(control('ribbon'), true);
        expect(control('width-mode').closest('label').hidden).toBe(true);
        expect(control('arrows').closest('label').hidden).toBe(true);
        expect(current().arrows).toBe(true);
        toggle(control('ribbon'), false);
        editStyle(root, 'cartoon-width-mode', 'custom');
        editStyle(root, 'cartoon-width', 1.6);
        expect(current().width).toBe(1.6);
        editStyle(root, 'cartoon-style', 'trace');
        expect(current()).toEqual({ ...cartoon, width: 1.6, ribbon: false, style: 'trace', arrows: false, tubes: false });
        for (const name of ['ribbon', 'arrows', 'tubes', 'width-mode', 'width']) {
            expect(control(name).closest('label').hidden).toBe(true);
        }
        expect(control('thickness').closest('label').textContent).toContain('Trace radius');
        editStyle(root, 'cartoon-thickness', 0.3);
        expect(current().thickness).toBe(0.3);
        expect(current().gapcutoff).toBe(4);
    });

    test('solid-color edits preserve all imported geometry options and display named colors correctly', () => {
        const { root, spec, current } = importedAppearance();
        const color = root.querySelector('.molui-style [data-focus$=":color"]');
        expect(color.value).toBe('#ff0000');
        editStyle(root, 'coloring', 'custom');
        expect(root.querySelector('.molui-style [data-focus$=":color"]').closest('label').hidden).toBe(false);
        editStyle(root, 'color', '#123456');
        expect(current()).toEqual({
            cartoon: { ...spec.cartoon, color: '#123456' },
            stick: { radius: 0.4, singleBonds: true, color: '#123456' },
            sphere: { ...spec.sphere, color: '#123456' }
        });
    });

    test('a pending surface stays disabled across refreshes and cannot be submitted twice', async () => {
        const { root, state, viewer } = mount();
        let finish;
        viewer.addSurface = jest.fn(() => Object.assign(new Promise(resolve => { finish = resolve; }), { surfid: 7 }));
        const pendingSurface = {};
        viewer.getSurface = id => id === 7 ? pendingSurface : undefined;
        viewer.removeSurface = jest.fn();
        button(root, 'Add surface').click();
        expect(button(root, 'Add surface').disabled).toBe(true);
        state.refresh();
        expect(button(root, 'Add surface').disabled).toBe(true);
        button(root, 'Add surface').click();
        expect(viewer.addSurface).toHaveBeenCalledTimes(1);
        finish(7);
        await settle();
        expect(state.surfaces.size).toBe(1);
        expect(button(root, 'Add surface').disabled).toBe(false);
    });

    test('surface calculation defaults to the shown atoms and can explicitly include all atoms', async () => {
        const { root, state, viewer } = mount();
        const id = state.addSelection({ chain: 'B' }, undefined, 'Ligand');
        const form = root.querySelector('.molui-surface-create');
        const calculation = () => field(form, 'Calculate using');
        expect(calculation().value).toBe('shown');
        expect([...calculation().options].find(option => option.value === 'shown').text).toBe('Same as shown');
        expect([...calculation().options].find(option => option.value === '').text).toBe('All atoms');
        choose(form, 'Show surface for', id);
        button(form, 'Add surface').click();
        await settle();
        expect(viewer.addSurface).toHaveBeenNthCalledWith(1, 1, expect.any(Object), { chain: 'B' }, { chain: 'B' });
        const first = [...state.surfaces.values()][0];
        expect(first.spec).toEqual({ chain: 'B' });
        expect(first.calculationSpec).toEqual({ chain: 'B' });
        choose(form, 'Calculate using', '');
        state.refresh();
        expect(calculation().value).toBe('');
        expect(field(form, 'Show surface for').value).toBe(id);
        expect(button(form, 'Add surface').disabled).toBe(false);
        button(form, 'Add surface').click();
        await settle();
        expect(viewer.addSurface).toHaveBeenNthCalledWith(2, 1, expect.any(Object), { chain: 'B' }, {});
        const second = [...state.surfaces.values()][1];
        expect(second.spec).toEqual({ chain: 'B' });
        expect(second.calculationSpec).toEqual({});
    });

    test('surface calculation can use a distinct None selection and preserves choices through renaming and refresh', async () => {
        const { root, state, viewer } = mount();
        choose(root, 'Chain', 'B');
        choose(root, 'Representation', 'none');
        button(root, 'Add selection').click();
        const shown = [...state.selections.values()][0];
        chooseMany(root, 'Chain', ['A', 'B']);
        button(root, 'Add selection').click();
        const context = [...state.selections.values()][1];
        expect(shown.styles.size).toBe(0);
        expect(context.styles.size).toBe(0);
        expect(context.spec.chain).toEqual(expect.arrayContaining(['A', 'B']));
        expect(context.spec.chain).toHaveLength(2);
        const form = root.querySelector('.molui-surface-create');
        choose(form, 'Show surface for', shown.id);
        choose(form, 'Calculate using', context.id);
        setInput(cardFor(root, shown.id).querySelector('.molui-selection-name'), 'Shown ligand');
        setInput(cardFor(root, context.id).querySelector('.molui-selection-name'), 'Surrounding assembly');
        state.refresh();
        expect(field(form, 'Show surface for').value).toBe(shown.id);
        expect(field(form, 'Calculate using').value).toBe(context.id);
        for (const label of ['Show surface for', 'Calculate using']) {
            const choices = [...field(form, label).options];
            expect(choices.find(option => option.value === shown.id).text).toBe('Shown ligand');
            expect(choices.find(option => option.value === context.id).text).toBe('Surrounding assembly');
        }
        expect(button(form, 'Add surface').disabled).toBe(false);
        button(form, 'Add surface').click();
        await settle();
        expect(viewer.addSurface).toHaveBeenCalledWith(1, expect.any(Object), { chain: 'B' }, context.spec);
        const surface = [...state.surfaces.values()][0];
        expect(surface.spec).toEqual({ chain: 'B' });
        expect(surface.calculationSpec).toEqual(context.spec);
        button(cardFor(root, context.id), 'Remove selection').click();
        expect(field(form, 'Calculate using').value).toBe('shown');
        expect(field(form, 'Show surface for').value).toBe(shown.id);
        expect(options(form, 'Show surface for')).not.toContain(context.id);
        expect(options(form, 'Calculate using')).not.toContain(context.id);
        button(cardFor(root, shown.id), 'Remove selection').click();
        expect(field(form, 'Show surface for').value).toBe('');
        expect(field(form, 'Calculate using').value).toBe('shown');
        expect(button(form, 'Add surface').disabled).toBe(false);
        expect(surface.spec).toEqual({ chain: 'B' });
        expect(surface.calculationSpec).toEqual(context.spec);
    });

    test('invalid calculation selections show inline guidance and prevent surface creation', () => {
        const { root, state, viewer } = mount();
        const shown = state.addSelection({ chain: 'B' }, undefined, 'Ligand');
        const partial = state.addSelection({ chain: 'B', elem: 'C' }, undefined, 'One ligand atom');
        const empty = state.addSelection({ chain: 'missing' }, undefined, 'No matching atoms');
        const form = root.querySelector('.molui-surface-create');
        const feedback = () => form.querySelector('.molui-surface-validation');
        choose(form, 'Show surface for', shown);
        choose(form, 'Calculate using', partial);
        expect(button(form, 'Add surface').disabled).toBe(true);
        expect(feedback().hidden).toBe(false);
        expect(feedback().dataset.error).toBe('true');
        expect(feedback().textContent).toBe('Calculation atoms must include all shown atoms.');
        button(form, 'Add surface').click();
        expect(viewer.addSurface).not.toHaveBeenCalled();
        state.refresh();
        expect(field(form, 'Calculate using').value).toBe(partial);
        expect(button(form, 'Add surface').disabled).toBe(true);
        choose(form, 'Calculate using', empty);
        expect(feedback().textContent).toBe('Choose atoms to calculate the surface.');
        expect(feedback().dataset.error).toBe('true');
        expect(button(form, 'Add surface').disabled).toBe(true);
        button(form, 'Add surface').click();
        expect(viewer.addSurface).not.toHaveBeenCalled();
        choose(form, 'Show surface for', empty);
        choose(form, 'Calculate using', '');
        expect(feedback().textContent).toBe('Choose atoms to show before adding a surface.');
        expect(button(form, 'Add surface').disabled).toBe(true);
        button(form, 'Add surface').click();
        expect(viewer.addSurface).not.toHaveBeenCalled();
        choose(form, 'Show surface for', shown);
        expect(button(form, 'Add surface').disabled).toBe(false);
        expect(feedback().dataset.error).not.toBe('true');
        expect(feedback().textContent).toContain('Calculation atoms shape the surface');
    });

    test('surface completion keeps Add disabled if calculation choices became invalid while pending', async () => {
        const { root, state, viewer } = mount();
        const partial = state.addSelection({ chain: 'B' }, undefined, 'Only the ligand');
        const form = root.querySelector('.molui-surface-create');
        let finish;
        viewer.addSurface = jest.fn(() => Object.assign(new Promise(resolve => { finish = resolve; }), { surfid: 9 }));
        const pendingSurface = {};
        viewer.getSurface = id => id === 9 ? pendingSurface : undefined;
        button(form, 'Add surface').click();
        expect(button(form, 'Add surface').disabled).toBe(true);
        choose(form, 'Calculate using', partial);
        state.refresh();
        expect(button(form, 'Add surface').disabled).toBe(true);
        finish(9);
        await settle();
        expect(state.surfaces.size).toBe(1);
        expect(field(form, 'Calculate using').value).toBe(partial);
        expect(button(form, 'Add surface').disabled).toBe(true);
        expect(form.querySelector('.molui-surface-validation').dataset.error).toBe('true');
        button(form, 'Add surface').click();
        expect(viewer.addSurface).toHaveBeenCalledTimes(1);
        choose(form, 'Calculate using', 'shown');
        expect(button(form, 'Add surface').disabled).toBe(false);
    });

    test('selection names can be edited without losing focus or dependent surface choices', () => {
        const { root, state } = mount();
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][0];
        const name = selectionCard(root).querySelector('.molui-selection-name');
        expect(name.getAttribute('aria-label')).toBe('Selection name');
        name.focus();
        setInput(name, '  Binding pocket  ');
        expect(selection.name).toBe('Binding pocket');
        expect(document.activeElement).toBe(selectionCard(root).querySelector('.molui-selection-name'));
        for (const label of ['Show surface for', 'Calculate using']) {
            const choices = [...field(root, label).options];
            expect(choices.find(option => option.value === selection.id).text).toBe('Binding pocket');
        }
        setInput(selectionCard(root).querySelector('.molui-selection-name'), '   ');
        expect(selection.name).toBe('Binding pocket');
        expect(selectionCard(root).querySelector('.molui-selection-name').value).toBe('Binding pocket');
    });

    test('one checkbox beside the name controls every representation in a selection', () => {
        const { root, state, models } = mount();
        choose(root, 'Representation', 'stick');
        button(root, 'Add selection').click();
        addRepresentation(root, 'sphere');
        const selection = [...state.selections.values()][0];
        const original = [...selection.styles.values()].map(style => style.spec);
        expect(root.textContent).not.toContain('Enable appearance');
        expect(root.textContent).not.toContain('Show representation');
        expect(root.querySelectorAll('.molui-style input[type=checkbox]')).toHaveLength(0);
        expect(selectionCard(root).querySelectorAll('.molui-selection-header input[type=checkbox]')).toHaveLength(1);
        toggle(selectionCard(root).querySelector('.molui-selection-header input[type=checkbox]'), false);
        expect(selection.hidden).toBe(true);
        expect(models[0].selectedAtoms({})[0].style).toEqual({ line: {} });
        toggle(selectionCard(root).querySelector('.molui-selection-header input[type=checkbox]'), true);
        expect(selection.hidden).toBe(false);
        expect([...selection.styles.values()].map(style => style.spec)).toEqual(original);
        expect(models[0].selectedAtoms({})[0].style).toHaveProperty('stick');
        expect(models[0].selectedAtoms({})[0].style).toHaveProperty('sphere');
    });

    test('residue labels are a selection toggle with no duplicates or unrelated-label cleanup', () => {
        const { root, state, viewer, labels, unrelatedLabel } = mount();
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][0];
        expect(button(root, 'Clear labels')).toBeUndefined();
        expect(button(root, 'Export VRML')).toBeUndefined();
        toggle(checkbox(selectionCard(root), 'Residue labels'), true);
        expect(selection.labelsEnabled).toBe(true);
        expect(labels.size).toBe(2);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        state.refresh();
        expect(checkbox(selectionCard(root), 'Residue labels').checked).toBe(true);
        expect(labels.size).toBe(2);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        toggle(selectionCard(root).querySelector('.molui-selection-header input[type=checkbox]'), false);
        expect(labels).toEqual(new Set([unrelatedLabel]));
        toggle(selectionCard(root).querySelector('.molui-selection-header input[type=checkbox]'), true);
        expect(labels.size).toBe(2);
        toggle(checkbox(selectionCard(root), 'Residue labels'), false);
        expect(selection.labelsEnabled).toBe(false);
        expect(labels).toEqual(new Set([unrelatedLabel]));
        toggle(checkbox(selectionCard(root), 'Residue labels'), true);
        button(selectionCard(root), 'Remove selection').click();
        expect(state.selections.size).toBe(0);
        expect(labels).toEqual(new Set([unrelatedLabel]));
    });

    test('None creates a selection without replacing baseline or overlapping selection appearances', () => {
        const { root, state, viewer } = mount();
        viewer.setStyle({}, { line: { color: 'cyan', linewidth: 2 } });
        const earlierId = state.createSelectionAndStyle({ chain: 'A' }, { sphere: { color: 'red', radius: 1 } });
        const before = JSON.parse(JSON.stringify(viewer.selectedAtoms({}).map(atom => atom.style)));
        chooseMany(root, 'Chain', ['', 'A']);
        expect([...field(root, 'Representation').options].find(option => option.value === 'none').text).toBe('None');
        choose(root, 'Representation', 'none');
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()].find(entry => entry.id !== earlierId);
        expect(selection.spec.chain).toEqual(['', 'A']);
        expect(selection.styles.size).toBe(0);
        expect(viewer.selectedAtoms({}).map(atom => atom.style)).toEqual(before);
        const card = () => cardFor(root, selection.id);
        expect(card().querySelector('.molui-style')).toBeNull();
        expect(card().querySelector('.molui-add-representation')).toBeNull();
        expect(field(card(), 'Representation').value).toBe('none');
        expect(card().textContent).toMatch(/surfaces/i);
        expect(card().textContent).toMatch(/residue labels/i);
        toggle(card().querySelector('.molui-selection-header input[type=checkbox]'), false);
        expect(viewer.selectedAtoms({}).map(atom => atom.style)).toEqual(before);
        toggle(card().querySelector('.molui-selection-header input[type=checkbox]'), true);
        expect(viewer.selectedAtoms({}).map(atom => atom.style)).toEqual(before);
        expect(state.selections.get(earlierId).styles.size).toBe(1);
    });

    test('a selection with representation None supports residue labels and surfaces', async () => {
        const { root, state, viewer, labels, unrelatedLabel } = mount();
        const before = JSON.parse(JSON.stringify(viewer.selectedAtoms({}).map(atom => atom.style)));
        choose(root, 'Chain', 'B');
        choose(root, 'Representation', 'none');
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][0];
        setInput(selectionCard(root).querySelector('.molui-selection-name'), 'Ligand context');
        toggle(checkbox(selectionCard(root), 'Residue labels'), true);
        state.refresh();
        expect(selection.styles.size).toBe(0);
        expect(selection.labelsEnabled).toBe(true);
        expect(labels.size).toBe(2);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        expect(viewer.addResLabels).toHaveBeenLastCalledWith({ chain: 'B' }, expect.any(Object));
        const surfaceForm = root.querySelector('.molui-surface-create');
        expect([...field(surfaceForm, 'Show surface for').options].find(option => option.value === selection.id).text).toBe('Ligand context');
        choose(surfaceForm, 'Show surface for', selection.id);
        button(surfaceForm, 'Add surface').click();
        await settle();
        expect(viewer.addSurface).toHaveBeenCalledTimes(1);
        expect(viewer.addSurface).toHaveBeenCalledWith(1, expect.any(Object), { chain: 'B' }, { chain: 'B' });
        expect(state.surfaces.size).toBe(1);
        expect(selection.styles.size).toBe(0);
        expect(viewer.selectedAtoms({}).map(atom => atom.style)).toEqual(before);
        toggle(checkbox(selectionCard(root), 'Residue labels'), false);
        expect(labels).toEqual(new Set([unrelatedLabel]));
        expect(selection.styles.size).toBe(0);
        expect(state.surfaces.size).toBe(1);
    });

    test('switching the last representation to None keeps labels and restores the underlying appearance', () => {
        const { root, state, viewer, labels } = mount();
        state.createSelectionAndStyle({}, { line: { color: 'blue', linewidth: 3 } });
        const underlying = JSON.parse(JSON.stringify(viewer.selectedAtoms({}).map(atom => atom.style)));
        choose(root, 'Chain', 'B');
        choose(root, 'Representation', 'stick');
        button(root, 'Add selection').click();
        const selection = [...state.selections.values()][1];
        const card = () => cardFor(root, selection.id);
        setInput(card().querySelector('.molui-selection-name'), 'Ligand labels');
        toggle(checkbox(card(), 'Residue labels'), true);
        expect(viewer.selectedAtoms({ chain: 'B' })[0].style).toHaveProperty('stick');
        choose(card().querySelector('.molui-style'), 'Representation', 'none');
        expect(state.selections.size).toBe(2);
        expect(selection.styles.size).toBe(0);
        expect(selection.spec).toEqual({ chain: 'B' });
        expect(selection.name).toBe('Ligand labels');
        expect(selection.labelsEnabled).toBe(true);
        expect(labels.size).toBe(2);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        expect(viewer.selectedAtoms({}).map(atom => atom.style)).toEqual(underlying);
        expect(card().querySelector('.molui-style')).toBeNull();
        expect(card().querySelector('.molui-add-representation')).toBeNull();
        expect(field(card(), 'Representation').value).toBe('none');
        expect(field(card(), 'Representation').dataset.focus).toBe(`${selection.id}:representation`);
        choose(card(), 'Representation', 'sphere');
        expect(selection.styles.size).toBe(1);
        expect([...selection.styles.values()][0].spec).toEqual({ sphere: { colorscheme: 'default', opacity: 1 } });
        expect(viewer.selectedAtoms({ chain: 'B' }).every(atom => atom.style.sphere && !atom.style.line)).toBe(true);
        expect(viewer.selectedAtoms({ chain: 'A' }).every(atom => atom.style.line.color === 'blue')).toBe(true);
        expect(selection.name).toBe('Ligand labels');
        expect(selection.labelsEnabled).toBe(true);
        expect(labels.size).toBe(2);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        expect(card().querySelectorAll('.molui-style')).toHaveLength(1);
        expect(options(card().querySelector('.molui-style'), 'Representation')).toContain('none');
        expect([...card().querySelector('.molui-add-representation select').options].map(option => option.value)).not.toContain('none');
    });

    test('choosing None removes only its representation from a selection with multiple styles', () => {
        const { root, state, viewer, labels } = mount();
        choose(root, 'Representation', 'stick');
        button(root, 'Add selection').click();
        addRepresentation(root, 'sphere');
        const selection = [...state.selections.values()][0];
        const [stickId, sphereId] = [...selection.styles.keys()];
        const sphere = JSON.parse(JSON.stringify(selection.styles.get(sphereId)));
        setInput(selectionCard(root).querySelector('.molui-selection-name'), 'Two styles');
        toggle(checkbox(selectionCard(root), 'Residue labels'), true);
        choose(selectionCard(root).querySelector('.molui-style'), 'Representation', 'none');
        expect(selection.styles.size).toBe(1);
        expect(selection.styles.has(stickId)).toBe(false);
        expect(selection.styles.get(sphereId)).toEqual(sphere);
        expect(selection.name).toBe('Two styles');
        expect(selection.labelsEnabled).toBe(true);
        expect(labels.size).toBe(2);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        expect(viewer.selectedAtoms({}).every(atom => atom.style.sphere && !atom.style.stick)).toBe(true);
        expect(options(selectionCard(root).querySelector('.molui-style'), 'Representation')).toContain('none');
        expect([...selectionCard(root).querySelector('.molui-add-representation select').options].map(option => option.value)).not.toContain('none');
    });

    test('multiple representations can be edited and removed independently', () => {
        const { root, state, models } = mount();
        choose(root, 'Representation', 'stick');
        button(root, 'Add selection').click();
        addRepresentation(root, 'sphere');
        const selection = [...state.selections.values()][0];
        expect(selection.styles.size).toBe(2);
        expect(root.querySelectorAll('.molui-style')).toHaveLength(2);
        const offered = [...root.querySelector('.molui-add-representation select').options].map(option => option.value);
        expect(offered).not.toContain('stick');
        expect(offered).not.toContain('sphere');
        expect(offered).not.toContain('ballstick');
        const [stickId, sphereId] = [...selection.styles.keys()];
        const stickBefore = JSON.parse(JSON.stringify(selection.styles.get(stickId).spec));
        let sphere = root.querySelectorAll('.molui-style')[1];
        sphere.querySelector('.molui-geometry').open = true;
        sphere.querySelector('[data-focus$=":opacity"]').focus();
        setInput(sphere.querySelector('[data-focus$=":opacity"]'), '0.5');
        sphere = root.querySelectorAll('.molui-style')[1];
        expect(selection.styles.get(sphereId).spec.sphere.opacity).toBe(0.5);
        expect(selection.styles.get(stickId).spec).toEqual(stickBefore);
        expect(root.querySelectorAll('.molui-style')[0].querySelector('.molui-geometry').open).toBe(false);
        expect(sphere.querySelector('.molui-geometry').open).toBe(true);
        expect(document.activeElement).toBe(sphere.querySelector('[data-focus$=":opacity"]'));
        choose(sphere, 'Color by', 'custom');
        sphere = root.querySelectorAll('.molui-style')[1];
        setInput(field(sphere, 'Solid color'), '#123456');
        expect(selection.styles.get(stickId).spec).toEqual(stickBefore);
        expect(selection.styles.get(sphereId).spec.sphere.color).toBe('#123456');
        button(root.querySelectorAll('.molui-style')[0], 'Remove representation').click();
        expect(selection.styles.size).toBe(1);
        expect(selection.styles.has(stickId)).toBe(false);
        expect(models[0].selectedAtoms({})[0].style).not.toHaveProperty('stick');
        expect(models[0].selectedAtoms({})[0].style.sphere.color).toBe('#123456');
    });

    test('named color schemes apply to representations and survive unrelated edits', () => {
        const { root, state } = mount();
        choose(root, 'Representation', 'stick');
        button(root, 'Add selection').click();
        choose(root.querySelector('.molui-style'), 'Color by', 'scheme');
        let editor = root.querySelector('.molui-style');
        expect(field(editor, 'Color scheme').closest('label').hidden).toBe(false);
        expect(options(editor, 'Color scheme')).toEqual(expect.arrayContaining(['greenCarbon', 'orangeCarbon', 'chainHetatm']));
        choose(editor, 'Color scheme', 'greenCarbon');
        const selection = [...state.selections.values()][0];
        const current = () => [...selection.styles.values()][0].spec.stick;
        expect(current().colorscheme).toBe('greenCarbon');
        expect(current().color).toBeUndefined();
        editStyle(root, 'opacity', 0.5);
        expect(current().colorscheme).toBe('greenCarbon');
        expect(current().opacity).toBe(0.5);
        editor = root.querySelector('.molui-style');
        expect(field(editor, 'Color by').value).toBe('scheme');
        expect(field(editor, 'Color scheme').value).toBe('greenCarbon');
        // These names also have shorthand Color by modes, but selecting them
        // through Color scheme must not close that dropdown after a refresh.
        for (const scheme of ['default', 'chain', 'ssJmol']) {
            choose(root.querySelector('.molui-style'), 'Color scheme', scheme);
            state.refresh();
            editor = root.querySelector('.molui-style');
            expect(current().colorscheme).toBe(scheme);
            expect(field(editor, 'Color by').value).toBe('scheme');
            expect(field(editor, 'Color scheme').closest('label').hidden).toBe(false);
            expect(field(editor, 'Color scheme').value).toBe(scheme);
        }
    });

    test('surface color schemes can be selected when creating and recoloring a surface', async () => {
        const { root, state, viewer } = mount();
        const surfaces = root.querySelector('.molui-surface-create');
        choose(surfaces, 'Color by', 'scheme');
        choose(surfaces, 'Color scheme', 'greenCarbon');
        button(root, 'Add surface').click();
        await settle();
        expect(viewer.addSurface).toHaveBeenCalledWith(1, expect.objectContaining({ opacity: 0.7, colorscheme: 'greenCarbon' }), {}, {});
        const surface = [...state.surfaces.values()][0];
        expect(surface.style.colorscheme).toBe('greenCarbon');
        let card = root.querySelector('.molui-surface-card');
        expect(field(card, 'Color scheme').value).toBe('greenCarbon');
        choose(card, 'Color scheme', 'orangeCarbon');
        await settle();
        expect(viewer.setSurfaceMaterialStyle).toHaveBeenLastCalledWith(surface.surfaceId,
            expect.objectContaining({ opacity: 0.7, colorscheme: 'orangeCarbon' }));
        expect(viewer.addSurface).toHaveBeenCalledTimes(1);
        expect(surface.style.colorscheme).toBe('orangeCarbon');
        card = root.querySelector('.molui-surface-card');
        choose(card, 'Color by', 'custom');
        card = root.querySelector('.molui-surface-card');
        setInput(field(card, 'Solid color'), '#000066');
        await settle();
        expect(surface.style.color).toBe('#000066');
        expect(surface.style.colorscheme).toBeUndefined();
        expect(surface.style.opacity).toBe(0.7);
    });

    test('Escape closes the panel and returns focus to its toggle', () => {
        const { root } = mount();
        field(root, 'Chain').focus();
        field(root, 'Chain').dispatchEvent(new dom.window.KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
        expect(root.classList.contains('molui-closed')).toBe(true);
        expect(document.activeElement).toBe(root.querySelector('.molui-toggle'));
        expect(document.activeElement.getAttribute('aria-expanded')).toBe('false');
    });

    test('initialization and destruction are idempotent and release all polling timers', () => {
        const { state, root, viewer, container } = mount();
        // JSDOM queues native details-toggle events independently of UI polling.
        jest.advanceTimersByTime(0);
        const initialTimers = jest.getTimerCount();
        expect(state.initiateUI()).toBe(state.ui);
        expect(container.querySelectorAll('.molui')).toHaveLength(1);
        expect(jest.getTimerCount()).toBe(initialTimers);
        const replacement = new StateManager(viewer, { open: true });
        states.push(replacement);
        jest.advanceTimersByTime(0);
        expect(root.isConnected).toBe(false);
        expect(container.querySelectorAll('.molui')).toHaveLength(1);
        expect(jest.getTimerCount()).toBe(initialTimers);
        state.destroy();
        replacement.destroy();
        state.destroy();
        jest.advanceTimersByTime(0);
        expect(jest.getTimerCount()).toBe(0);
        expect(root.isConnected).toBe(false);
        expect(viewer.ui).toBeUndefined();
        expect(container.style.position).toBe('static');
    });
});
