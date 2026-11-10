/** @jest-environment node */
// Exercise the guided editor with real atom matching and a DOM-only viewer.
jest.mock('canvas', () => ({}));
const fs = require('node:fs');
const path = require('node:path');
const ts = require('typescript');
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!doctype html><html><body></body></html>');
global.window = dom.window;
global.document = dom.window.document;
global.URL.createObjectURL = jest.fn();
dom.window.URL.createObjectURL = jest.fn();
Object.defineProperty(global, '$3Dmol', {
    configurable: true, get: () => dom.window.$3Dmol, set: value => { dom.window.$3Dmol = value; }
});
const { GLModel } = require('../../build/3Dmol.js');
const source = fs.readFileSync(path.join(__dirname, '../../src/ui/selection.ts'), 'utf8');
const compiled = ts.transpileModule(source, {
    compilerOptions: { module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2019 }
}).outputText;
const editorModule = {};
new Function('exports', 'document', 'Option', 'HTMLInputElement', compiled)(
    editorModule, dom.window.document, dom.window.Option, dom.window.HTMLInputElement
);
const { SelectionEditor, selectionSummary } = editorModule;
const { Event, KeyboardEvent } = dom.window;

afterAll(() => {
    dom.window.close();
    delete global.window;
    delete global.document;
    delete global.$3Dmol;
});
afterEach(() => document.body.replaceChildren());

function fixture(spec = {}) {
    const model = new GLModel(0);
    model.addAtoms(Array.from({ length: 20 }, (_, index) => ({
        x: index, y: 0, z: 0, bonds: [], chain: 'A', resi: index,
        elem: index % 2 ? 'N' : 'C', atom: index % 2 ? 'N' : 'CA', resn: 'ALA', hetflag: false
    })));
    const onChange = jest.fn();
    const editor = new SelectionEditor({
        viewer: { selectedAtoms: selection => model.selectedAtoms(selection) },
        models: () => [{ id: 0, title: 'Protein' }], key: 'test-selection', spec, onChange
    });
    document.body.append(editor.element);
    const list = key => editor.element.querySelector(`.molui-multiselect[data-field="${key}"]`);
    const choice = (key, value) => [...list(key).querySelectorAll('input[data-value]')]
        .find(input => input.dataset.value === String(value));
    const toggle = (key, value, checked) => {
        const input = choice(key, value);
        input.checked = checked;
        input.dispatchEvent(new Event('change', { bubbles: true }));
    };
    return { editor, model, onChange, list, choice, toggle };
}

describe('selection editor', () => {
    test('search filters categorical checklists and preserves disclosure, caret, and choices across refresh', () => {
        const { editor, model, list, choice, toggle, onChange } = fixture();
        model.selectedAtoms({}).forEach((atom, index) => atom.chain = `A${index}`);
        editor.refresh();
        list('chain').open = true;
        let search = list('chain').querySelector('input[type=search]');
        expect(search.hidden).toBe(false);
        search.focus();
        search.value = '19'; search.setSelectionRange(1, 1);
        search.dispatchEvent(new Event('input', { bubbles: true }));
        expect(editor.getSpec()).toEqual({});
        expect(onChange).not.toHaveBeenCalled();
        expect(choice('chain', 'A19').closest('label').hidden).toBe(false);
        expect(choice('chain', 'A1').closest('label').hidden).toBe(true);
        editor.refresh();
        search = list('chain').querySelector('input[type=search]');
        expect(list('chain').open).toBe(true);
        expect(search.value).toBe('19');
        expect(document.activeElement).toBe(search);
        expect(search.selectionStart).toBe(1);
        choice('chain', 'A19').focus();
        toggle('chain', 'A19', true);
        expect(document.activeElement).toBe(choice('chain', 'A19'));
        expect(editor.getSpec()).toEqual({ chain: 'A19' });
        search = list('chain').querySelector('input[type=search]');
        search.value = ''; search.dispatchEvent(new Event('input', { bubbles: true }));
        expect(choice('chain', 'A1').closest('label').hidden).toBe(false);
        expect(choice('chain', 'A19').checked).toBe(true);
        toggle('chain', 'A1', true);
        expect(editor.getSpec().chain).toEqual(expect.arrayContaining(['A1', 'A19']));
        expect(editor.getAtoms()).toHaveLength(2);
    });

    const input = (list, name) => list('resi').querySelector(`[data-focus$=":resi:${name}"]`);
    const typeRange = (list, from, to) => {
        for (const [name, value] of [['from', from], ['to', to]]) {
            const control = input(list, name); control.value = String(value);
            control.dispatchEvent(new Event('input', { bubbles: true }));
        }
    };
    const addRange = (list, from, to) => { typeRange(list, from, to); input(list, 'add').click(); };

    test('Escape closes the active residue picker and All clears only its criterion', () => {
        const { editor, list } = fixture({ chain: 'A', resi: 0 });
        list('resi').open = true;
        input(list, 'from').focus();
        input(list, 'from').dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
        expect(list('resi').open).toBe(false);
        expect(document.activeElement).toBe(list('resi').querySelector('summary'));
        list('resi').open = true;
        addRange(list, 1, 3);
        input(list, 'all').click();
        expect(editor.getSpec()).toEqual({ chain: 'A' });
        expect(editor.getAtoms()).toHaveLength(20);
        expect(list('resi').querySelectorAll('.molui-residue-chip')).toHaveLength(0);
        expect(list('resi').querySelectorAll('input[type=checkbox]')).toHaveLength(0);
    });

    test('dense intervals merge overlaps and adjacency while separate spans stay compact', () => {
        const { editor, list } = fixture();
        expect(list('resi').dataset.min).toBe('0');
        expect(list('resi').dataset.max).toBe('19');
        addRange(list, 3, 12);
        expect(editor.getSpec()).toEqual({ resi: '3-12' });
        expect(editor.getAtoms()).toHaveLength(10);
        addRange(list, 10, 15);
        addRange(list, 16, 16);
        addRange(list, 19, 19);
        expect(editor.getSpec().resi).toEqual(['3-16', 19]);
        expect(editor.getAtoms()).toHaveLength(15);
        expect(selectionSummary(editor.getSpec())).toBe('Residues 3–16, 19');
        expect(list('resi').querySelectorAll('.molui-residue-chip')).toHaveLength(2);
    });

    test('range chips can be updated with Enter and removed without silently selecting everything', () => {
        const { editor, list } = fixture({ chain: 'A', resi: ['0-3', 10] });
        list('resi').querySelector('.molui-residue-chip button').click();
        expect(input(list, 'add').textContent).toBe('Update range');
        typeRange(list, 1, 4);
        input(list, 'to').dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
        expect(editor.getSpec()).toEqual({ chain: 'A', resi: ['1-4', 10] });
        list('resi').querySelector('.molui-residue-chip button[aria-label^="Remove"]').click();
        expect(editor.getSpec().resi).toBe(10);
        list('resi').querySelector('.molui-residue-chip button[aria-label^="Remove"]').click();
        expect(editor.getSpec().resi).toEqual([]);
        expect(editor.getAtoms()).toHaveLength(0);
        expect(list('resi').querySelector('summary').textContent).toContain('None');
        input(list, 'all').click();
        expect(editor.getSpec()).toEqual({ chain: 'A' });
    });

    test('negative and huge intervals match atoms without expanding a dense range into an array', () => {
        const { editor, model, list } = fixture();
        model.selectedAtoms({}).forEach(atom => atom.resi -= 5);
        editor.refresh();
        addRange(list, -4, 2);
        expect(editor.getSpec().resi).toBe('-4-2');
        expect(editor.getAtoms()).toHaveLength(7);
        addRange(list, -1000000000, 1000000000);
        expect(editor.getSpec().resi).toBe('-1000000000-1000000000');
        expect(editor.getAtoms()).toHaveLength(20);
        expect(list('resi').querySelectorAll('.molui-residue-chip')).toHaveLength(1);
        expect(list('resi').querySelectorAll('input[type=checkbox]')).toHaveLength(0);
    });

    test('invalid and unfinished bounds never commit, while drafts survive unrelated refreshes', () => {
        const { editor, list, toggle, onChange } = fixture({ chain: 'A' });
        list('resi').open = true;
        for (const [from, to] of [['', 5], [9, 2], [1.5, 5], [0, Number.MAX_SAFE_INTEGER + 1]]) {
            addRange(list, from, to);
            expect(editor.getSpec()).toEqual({ chain: 'A' });
            expect(list('resi').querySelector('.molui-residue-error').hidden).toBe(false);
        }
        expect(onChange).not.toHaveBeenCalled();
        typeRange(list, 3, 9);
        input(list, 'from').focus();
        editor.refresh();
        expect(input(list, 'from').value).toBe('3');
        expect(input(list, 'to').value).toBe('9');
        expect(document.activeElement).toBe(input(list, 'from'));
        expect(list('resi').open).toBe(true);
        toggle('elem', 'C', true);
        expect(input(list, 'from').value).toBe('3');
        expect(input(list, 'to').value).toBe('9');
        expect(editor.getSpec()).toEqual({ chain: 'A', elem: 'C' });
        input(list, 'add').click();
        expect(editor.getSpec()).toEqual({ chain: 'A', elem: 'C', resi: '3-9' });
        expect(editor.getAtoms()).toHaveLength(3);
    });

    test('unavailable chosen intervals stay visible and numeric drafts extend the sliders without clamping', () => {
        const { editor, list } = fixture({ resi: '100-110' });
        expect(editor.getAtoms()).toHaveLength(0);
        expect(list('resi').querySelector('.molui-residue-chip').dataset.from).toBe('100');
        expect(list('resi').dataset.max).toBe('19');
        list('resi').querySelector('.molui-residue-chip button').click();
        expect(input(list, 'from').value).toBe('100');
        expect(input(list, 'to-slider').value).toBe('110');
        typeRange(list, -10, 200);
        expect(input(list, 'from-slider').min).toBe('-10');
        expect(input(list, 'to-slider').max).toBe('200');
        expect(input(list, 'to-slider').value).toBe('200');
        expect(editor.getSpec().resi).toBe('100-110');
        input(list, 'cancel').click();
        expect(editor.getSpec().resi).toBe('100-110');
    });

    test('unexposed predicates, nested logic, ranges and explicit empty lists survive unrelated edits', () => {
        const predicate = atom => atom.x < 12;
        const spec = {
            model: [0], chain: ['A'], resi: '0-9', byres: true,
            and: [{ not: { elem: 'O' } }], predicate
        };
        const { editor, toggle } = fixture(spec);
        expect(editor.getSpec()).toEqual(spec);
        expect(editor.getAtoms()).toHaveLength(10);
        expect(selectionSummary(spec)).toContain('Additional constraints');
        expect(selectionSummary(spec)).toContain('0–9');
        toggle('elem', 'C', true);
        const next = editor.getSpec();
        expect(next).toEqual({ ...spec, elem: 'C' });
        expect(next.predicate).toBe(predicate);
        expect(editor.getAtoms()).toHaveLength(5);
        next.and[0].not.elem = 'C';
        expect(editor.getSpec().and[0].not.elem).toBe('O');
        expect(spec).not.toHaveProperty('elem');
        editor.setSpec({ chain: [], predicate });
        editor.refresh();
        expect(editor.getSpec()).toEqual({ chain: [], predicate });
        expect(editor.getAtoms()).toHaveLength(0);
        expect(selectionSummary(editor.getSpec())).toContain('(none)');
    });
});
