import type { AtomSelectionSpec, AtomSpec } from '../specs';
import type { AtomStyleSpec, GLModel } from '../GLModel';
import type { StateManager } from './state';
import { builtinColorSchemes, CC } from '../colors';
import type { ColorSpec, ColorschemeSpec } from '../colors';
import { TrajectoryControls } from './trajectory';
import { uiStyles } from './styles';
import { SelectionEditor, selectionSummary } from './selection';

export interface UIConfig {
    /** Start with the sidebar open. Small embedded viewers default to collapsed. */
    open?: boolean;
}

type Choice = [string, string];
type UIColors = { color?: ColorSpec; colorscheme?: ColorschemeSpec };
const schemeChoices: Choice[] = Object.keys(builtinColorSchemes).map(name => [name, name]);
const formats = ['pdb', 'cif', 'bcif', 'sdf', 'mol2', 'xyz', 'pqr', 'pdbqt', 'gro', 'lammpstrj', 'prmtop', 'mmtf'];
const water = ['HOH', 'WAT', 'H2O'];
let nextId = 0;

function el<K extends keyof HTMLElementTagNameMap>(tag: K, className = '', text?: string): HTMLElementTagNameMap[K] {
    const node = document.createElement(tag);
    node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
}
function button(text: string, action: () => void, className = ''): HTMLButtonElement {
    const node = el('button', className, text);
    node.type = 'button';
    node.addEventListener('click', action);
    return node;
}
function select(choices: Choice[], value?: string): HTMLSelectElement {
    const node = el('select');
    choices.forEach(([key, text]) => node.add(new Option(text, key)));
    if (value !== undefined) node.value = value;
    return node;
}
function field(text: string, input: HTMLElement): HTMLLabelElement {
    const label = el('label', 'molui-field');
    label.append(el('span', '', text), input);
    return label;
}
function check(text: string, checked = false): { label: HTMLLabelElement; input: HTMLInputElement } {
    const input = el('input'); input.type = 'checkbox'; input.checked = checked;
    const label = el('label', 'molui-check'); label.append(input, document.createTextNode(text));
    return { label, input };
}
function section(title: string, open = true): HTMLDetailsElement {
    const node = el('details', 'molui-section'); node.open = open;
    node.append(el('summary', '', title));
    return node;
}
function range(min: number, max: number, step: number, value: number): HTMLInputElement {
    const input = el('input'); input.type = 'range';
    Object.assign(input, { min: String(min), max: String(max), step: String(step), value: String(value) });
    return input;
}
function isPolymer(atoms: AtomSpec[]): boolean {
    return atoms.some(atom => !atom.hetflag && ['CA', 'P'].includes(atom.atom));
}
function errorText(error: unknown): string { return error instanceof Error ? error.message : String(error); }

/** Guided, instance-scoped controls for an existing viewer. Use StateManager to attach. */
export class UI {
    public readonly element = el('div', 'molui');
    private readonly panel = el('aside', 'molui-panel');
    private readonly toggle: HTMLButtonElement;
    private readonly status = el('div', 'molui-status');
    private readonly modelList = el('div', 'molui-list');
    private readonly selectionList = el('div', 'molui-list molui-selection-list');
    private readonly surfaceList = el('div', 'molui-list');
    private readonly selectionEditor: SelectionEditor;
    private readonly selectionEditors = new Map<string, SelectionEditor>();
    private readonly collapsedSelections = new Set<string>();
    private readonly count = el('p', 'molui-count');
    private readonly trajectory: TrajectoryControls;
    private readonly appearance = select([]);
    private readonly addSelectionButton: HTMLButtonElement;
    private readonly surfaceSelection = select([]);
    private readonly surfaceCalculation = select([['shown', 'Same as shown'], ['', 'All atoms']]);
    private readonly surfaceValidation = el('p', 'molui-muted molui-surface-validation');
    private readonly surfaceButton: HTMLButtonElement;
    private selection: AtomSelectionSpec = {};
    private selectionName = 'All atoms';
    private readonly colorModes = new Map<string, string>();
    private observer?: ResizeObserver;
    private unsubscribe: () => void;
    private poll: ReturnType<typeof setInterval>;
    private signature = '';
    private destroyed = false;
    private loading = false;
    private surfacePending = false;
    private renamingSelection = false;
    private abort?: AbortController;
    private originalPosition: string;
    private setPosition = false;
    private readonly onResize = () => this.resize();

    constructor(public readonly state: StateManager, config: UIConfig = {}) {
        const viewer = state.viewer;
        const container = viewer.container;
        this.originalPosition = container.style.position;
        if (getComputedStyle(container).position === 'static') {
            container.style.position = 'relative'; this.setPosition = true;
        }
        const style = el('style'); style.textContent = uiStyles;
        this.element.append(style);
        this.element.setAttribute('aria-label', 'Molecular viewer controls');
        this.panel.id = `molui-panel-${++nextId}`;
        this.panel.setAttribute('aria-label', 'Molecule settings');
        this.toggle = button('☰', () => this.setOpen(this.element.classList.contains('molui-closed')),
            'molui-toggle molui-icon');
        this.toggle.setAttribute('aria-label', 'Toggle controls');
        this.toggle.title = 'Toggle controls';
        this.toggle.setAttribute('aria-controls', this.panel.id);
        const header = el('div', 'molui-header');
        const close = button('×', () => this.setOpen(false), 'molui-icon');
        close.setAttribute('aria-label', 'Close controls');
        close.title = 'Close controls';
        header.append(close);
        const body = el('div', 'molui-body');
        this.status.setAttribute('role', 'status'); this.status.setAttribute('aria-live', 'polite');
        this.status.hidden = true;
        const models = section('Models');
        models.append(this.modelList, this.loadForm());
        const selections = section('Selections & appearance');
        this.addSelectionButton = button('Add selection', () => {
            if (!this.state.checkAtoms(this.selection)) return;
            const id = this.state.addSelection(this.selection, undefined, this.selectionName);
            if (this.appearance.value !== 'none') this.state.addStyle(this.styleSpec(this.appearance.value), id);
            this.showError();
            this.focusControl(`${id}:name`);
        }, 'molui-primary');
        this.selectionEditor = this.makeSelectionEditor('new-selection', {}, () => this.updateSelectionDraft());
        const create = el('div', 'molui-selection-create');
        create.append(el('h3', '', 'New selection'), this.selectionEditor.element, this.count,
            field('Representation', this.appearance), this.addSelectionButton);
        selections.append(this.selectionList, create);
        const surfaces = section('Surfaces', false);
        const surfaceForm = el('div', 'molui-surface-create');
        const surfaceColors = this.colorControls({}, 'surface-new');
        const type = select([['1', 'Van der Waals'], ['2', 'Solvent excluded'], ['3', 'Solvent accessible'], ['4', 'Molecular']]);
        const opacity = range(0.1, 1, 0.05, 0.7);
        const opacityField = field('Opacity · 70%', opacity);
        opacity.addEventListener('input', () => opacityField.firstElementChild.textContent = `Opacity · ${Math.round(Number(opacity.value) * 100)}%`);
        this.surfaceButton = button('Add surface', async () => {
            const entry = this.state.selections.get(this.surfaceSelection.value);
            const spec = entry?.spec || {};
            this.updateSurfaceSelection();
            if (this.surfaceButton.disabled) return;
            const calculation = this.surfaceCalculation.value === 'shown' ? spec
                : this.state.selections.get(this.surfaceCalculation.value)?.spec || {};
            this.surfacePending = true;
            this.surfaceButton.disabled = true;
            this.showError();
            try {
                await this.state.addSurface(Number(type.value), { ...surfaceColors.value(), opacity: Number(opacity.value) },
                    spec, entry?.name || 'All atoms', calculation);
            } catch (error) { this.showError(errorText(error)); }
            finally { this.surfacePending = false; if (!this.destroyed) this.updateSurfaceSelection(); }
        }, 'molui-primary');
        this.surfaceValidation.setAttribute('aria-live', 'polite');
        this.surfaceSelection.addEventListener('change', () => this.updateSurfaceSelection());
        this.surfaceCalculation.addEventListener('change', () => this.updateSurfaceSelection());
        surfaceForm.append(field('Show surface for', this.surfaceSelection), field('Calculate using', this.surfaceCalculation),
            this.surfaceValidation, field('Surface type', type), ...surfaceColors.fields, opacityField, this.surfaceButton);
        surfaces.append(surfaceForm, this.surfaceList);
        body.append(this.status, models, selections, surfaces, this.viewControls());
        this.panel.append(header, body);
        this.trajectory = new TrajectoryControls(viewer, message => this.showError(message));
        this.element.append(this.toggle, this.panel, this.trajectory.element);
        // Prevent viewer gesture handlers from treating controls as canvas interaction.
        for (const event of ['mousedown', 'pointerdown', 'touchstart', 'wheel', 'dblclick', 'contextmenu']) {
            this.element.addEventListener(event, e => e.stopPropagation());
        }
        this.element.addEventListener('keydown', event => {
            if (event.key === 'Escape' && !this.element.classList.contains('molui-closed')) {
                this.setOpen(false); this.toggle.focus(); event.stopPropagation();
            }
        });
        container.append(this.element);
        this.setOpen(config.open ?? (container.clientWidth >= 640 && container.clientHeight >= 340));
        this.unsubscribe = state.subscribe(() => this.refresh());
        window.addEventListener('resize', this.onResize);
        if (typeof ResizeObserver !== 'undefined') {
            this.observer = new ResizeObserver(this.onResize); this.observer.observe(container);
        }
        // Model/frame APIs can also be called by host applications. Only rebuild on topology changes.
        this.poll = setInterval(() => {
            if (this.destroyed) return;
            if (!this.element.isConnected) { this.state.destroy(); return; }
            if (document.hidden) return;
            const signature = this.modelSignature();
            if (signature !== this.signature) this.refresh();
            this.element.classList.toggle('molui-has-trajectory', viewer.getNumFrames() > 1);
            this.resize();
        }, 500);
        this.refresh(); this.resize();
    }

    public resize(): void {
        this.element.classList.toggle('molui-compact', this.state.viewer.container.clientWidth < 640);
        this.element.classList.toggle('molui-small', this.state.viewer.container.clientHeight < 440);
        this.element.classList.toggle('molui-tiny', this.state.viewer.container.clientHeight < 320);
        this.trajectory.resize();
        this.element.style.setProperty('--molui-trajectory-space', `${this.trajectory.element.offsetHeight + 24}px`);
    }

    public refresh(): void {
        if (this.destroyed) return;
        if (this.renamingSelection) {
            // A name's blur event precedes a click on the next control. Keep
            // those DOM nodes mounted so committing a name does not eat the click.
            for (const input of [this.surfaceSelection, this.surfaceCalculation]) {
                for (const option of Array.from(input.options)) {
                    const selection = this.state.selections.get(option.value);
                    if (selection) option.text = selection.name;
                }
            }
            return;
        }
        const active = document.activeElement as HTMLElement;
        const focusKey = this.element.contains(active) ? active?.dataset.focus : undefined;
        this.signature = this.modelSignature();
        this.renderModels(); this.renderFilters(); this.renderSelections(); this.renderSurfaces();
        this.trajectory.refresh();
        if (focusKey) {
            const control = Array.from(this.element.querySelectorAll<HTMLElement>('[data-focus]'))
                .find(node => node.dataset.focus === focusKey);
            control?.focus({ preventScroll: true });
        }
        this.element.classList.toggle('molui-has-trajectory', this.state.viewer.getNumFrames() > 1);
    }

    public destroy(): void {
        if (this.destroyed) return;
        this.destroyed = true;
        this.abort?.abort(); clearInterval(this.poll); this.unsubscribe?.(); this.observer?.disconnect();
        window.removeEventListener('resize', this.onResize);
        this.trajectory.destroy();
        this.element.remove();
        if (this.setPosition && this.state.viewer.container.style.position === 'relative') {
            this.state.viewer.container.style.position = this.originalPosition;
        }
    }

    private setOpen(open: boolean): void {
        if (!open && this.panel.contains(document.activeElement)) this.toggle.focus();
        this.element.classList.toggle('molui-closed', !open);
        this.toggle.setAttribute('aria-expanded', String(open));
        this.resize();
    }

    private showError(text = ''): void {
        if (this.destroyed) return;
        this.status.hidden = !text;
        this.status.textContent = text;
        this.status.dataset.error = String(!!text);
        if (text && this.status.parentElement) this.status.parentElement.scrollTop = 0;
    }

    private models(): GLModel[] {
        const models: GLModel[] = [];
        for (let id = 0; id < this.state.viewer.getNextModelId(); id++) {
            const model = this.state.viewer.getModel(id);
            if (model && model.getID() === id) models.push(model);
        }
        return models;
    }

    private modelSignature(): string {
        return this.models().map(model => `${model.getID()}:${model.selectedAtoms({}).length}:${model.getNumFrames()}:${model.isHidden()}`).join('|');
    }

    private renderModels(): void {
        this.modelList.replaceChildren();
        const models = this.models();
        if (!models.length) this.modelList.append(el('p', 'molui-empty', 'Load a structure to explore its atoms and appearance.'));
        for (const model of models) {
            const id = model.getID();
            const card = el('div', 'molui-card');
            const header = el('div', 'molui-card-header');
            const title = this.state.modelTitles.get(id) || `Model ${id + 1}`;
            header.append(el('strong', '', title));
            const actions = el('div', 'molui-actions');
            const visible = check('Show model', !model.isHidden());
            visible.input.addEventListener('change', () => {
                if (visible.input.checked) model.show(); else model.hide();
                this.state.viewer.render();
            });
            actions.append(button('Focus', () => { this.state.viewer.zoomTo({ model: id }); this.state.viewer.render(); }),
                button('Remove model', () => {
                    this.trajectory.pause();
                    this.state.removeModel(model);
                }, 'molui-danger'));
            card.append(visible.label);
            card.append(header, el('p', 'molui-muted', `${model.selectedAtoms({}).length.toLocaleString()} atoms · ${model.getNumFrames()} frame${model.getNumFrames() === 1 ? '' : 's'}`), actions);
            this.modelList.append(card);
        }
    }

    private loadForm(): HTMLDetailsElement {
        const details = section('Load structure', !this.models().length);
        const form = el('form');
        const source = select([['file', 'Local file'], ['pdb', 'PDB ID'], ['cid', 'PubChem CID'], ['url', 'URL']]);
        const file = el('input'); file.type = 'file';
        file.accept = formats.map(format => `.${format},.${format}.gz`).join(',');
        const identifier = el('input'); identifier.type = 'text'; identifier.autocomplete = 'off';
        const format = select([['auto', 'Detect from filename'], ...formats.map(f => [f, f.toUpperCase()] as Choice)]);
        const trajectory = check('Load multiple models as trajectory frames', true);
        const identifierField = field('PDB ID', identifier);
        const fileField = field('Structure file', file);
        const formatField = field('File format', format);
        const help = el('p', 'molui-muted', 'Files stay in your browser. New structures are added to the scene.');
        const submit = el('button', 'molui-primary', 'Load structure'); submit.type = 'submit';
        const update = () => {
            fileField.hidden = source.value !== 'file'; identifierField.hidden = source.value === 'file';
            formatField.hidden = !['file', 'url'].includes(source.value);
            identifierField.firstElementChild.textContent = source.value === 'pdb' ? 'PDB ID' : source.value === 'cid' ? 'PubChem CID' : 'Structure URL';
            identifier.placeholder = source.value === 'pdb' ? 'e.g. 1CRN' : source.value === 'cid' ? 'e.g. 2244' : 'https://…/structure.pdb';
            help.textContent = source.value === 'url' ? 'The server must allow cross-origin requests. Choose a format if the URL has no extension.'
                : source.value === 'file' ? 'Files stay in your browser. New structures are added to the scene.' : 'New structures are added to the scene.';
        };
        source.addEventListener('change', update); update();
        form.append(field('Source', source), fileField, identifierField, formatField, trajectory.label, help, submit);
        form.addEventListener('submit', async event => {
            event.preventDefault();
            if (this.loading) return;
            this.loading = true; submit.disabled = true; submit.textContent = 'Loading…';
            form.setAttribute('aria-busy', 'true');
            this.abort = new AbortController();
            try {
                let name: string; let data: string | ArrayBuffer; let fileFormat = format.value;
                if (source.value === 'file') {
                    const selected = file.files?.[0];
                    if (!selected) throw new Error('Choose a structure file first.');
                    name = selected.name;
                    if (fileFormat === 'auto') fileFormat = this.detectFormat(name);
                    else if (/\.gz$/i.test(name)) fileFormat += '.gz';
                    data = await this.readFile(selected, /(?:bcif|mmtf|gz)$/.test(fileFormat));
                } else {
                    const value = identifier.value.trim();
                    let url: URL;
                    if (source.value === 'pdb') {
                        if (!/^[a-z0-9]{4}$/i.test(value)) throw new Error('Enter a four-character PDB ID, such as 1CRN.');
                        name = value.toUpperCase(); fileFormat = 'pdb';
                        url = new URL(`https://files.rcsb.org/view/${name}.pdb`);
                    } else if (source.value === 'cid') {
                        if (!/^\d+$/.test(value)) throw new Error('Enter a numeric PubChem compound ID.');
                        name = `CID ${value}`; fileFormat = 'sdf';
                        url = new URL(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${value}/SDF?record_type=3d`);
                    } else {
                        try { url = new URL(value); } catch { throw new Error('Enter a complete https:// or http:// structure URL.'); }
                        if (!['https:', 'http:'].includes(url.protocol)) throw new Error('Use an https:// or http:// URL.');
                        name = decodeURIComponent(url.pathname.split('/').pop()) || url.hostname;
                        if (fileFormat === 'auto') fileFormat = this.detectFormat(url.pathname);
                        else if (/\.gz$/i.test(url.pathname)) fileFormat += '.gz';
                    }
                    const response = await fetch(url.href, { signal: this.abort.signal });
                    if (!response.ok) throw new Error(`The structure could not be downloaded (HTTP ${response.status}).`);
                    data = /(?:bcif|mmtf|gz)$/.test(fileFormat) ? await response.arrayBuffer() : await response.text();
                }
                if (this.destroyed) return;
                const model = trajectory.input.checked ? this.state.viewer.addModelsAsFrames(data, fileFormat)
                    : this.state.viewer.addModel(data, fileFormat);
                if (!model.selectedAtoms({}).length) {
                    this.state.viewer.removeModel(model);
                    throw new Error('No atoms were found. Check the file format and try again.');
                }
                this.state.modelTitles.set(model.getID(), name);
                const atoms = model.selectedAtoms({});
                model.setStyle({}, isPolymer(atoms) ? { cartoon: { color: 'spectrum' } } : { stick: {} });
                if (isPolymer(atoms)) model.setStyle({ hetflag: true, not: { or: water.map(resn => ({ resn })) } }, { stick: {} });
                this.state.viewer.zoomTo({ model: model.getID() }); this.state.viewer.render();
                this.state.refresh();
                details.open = false;
                this.showError();
            } catch (error) {
                if (!this.destroyed) this.showError(`Unable to load structure. ${errorText(error)}`);
            } finally {
                this.loading = false; submit.disabled = false; submit.textContent = 'Load structure'; form.removeAttribute('aria-busy');
            }
        });
        details.append(form);
        return details;
    }

    private detectFormat(name: string): string {
        const match = name.toLowerCase().match(/\.([a-z0-9]+)(\.gz)?$/);
        if (!match || !formats.includes(match[1])) throw new Error('Choose the file format; it could not be detected from the filename.');
        return match[1] + (match[2] || '');
    }

    private readFile(file: File, binary: boolean): Promise<string | ArrayBuffer> {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = () => resolve(reader.result as string | ArrayBuffer);
            reader.onerror = () => reject(new Error('The file could not be read.'));
            if (binary) reader.readAsArrayBuffer(file); else reader.readAsText(file);
        });
    }

    private makeSelectionEditor(key: string, spec: AtomSelectionSpec, onChange: () => void): SelectionEditor {
        return new SelectionEditor({
            viewer: this.state.viewer,
            models: () => this.models().map(model => ({
                id: model.getID(), title: this.state.modelTitles.get(model.getID()) || `Model ${model.getID() + 1}`
            })),
            key, spec, onChange
        });
    }

    private renderFilters(): void {
        this.selectionEditor.refresh();
        this.updateSelectionDraft();
    }

    private updateSelectionDraft(): void {
        const atoms = this.selectionEditor.getAtoms();
        this.selection = this.selectionEditor.getSpec();
        this.selectionName = selectionSummary(this.selection, this.state.modelTitles);
        this.count.textContent = `${atoms.length.toLocaleString()} matching atom${atoms.length === 1 ? '' : 's'}`;
        this.addSelectionButton.disabled = atoms.length === 0;
        const previous = this.appearance.value;
        this.appearance.replaceChildren();
        this.representations(atoms).forEach(([key, label]) => this.appearance.add(new Option(label, key)));
        this.appearance.value = Array.from(this.appearance.options).some(o => o.value === previous) ? previous : isPolymer(atoms) ? 'cartoon' : 'stick';
        this.appearance.disabled = !atoms.length;
    }

    private focusControl(key: string): void {
        Array.from(this.element.querySelectorAll<HTMLElement>('[data-focus]'))
            .find(node => node.dataset.focus === key)?.focus();
    }

    private editSelection(id: string): void {
        const entry = this.state.selections.get(id);
        if (!entry) return;
        this.collapsedSelections.delete(id);
        if (!this.selectionEditors.has(id)) {
            this.selectionEditors.set(id, this.makeSelectionEditor(`selection:${id}:edit`, entry.spec,
                () => this.updateSelectionEdit(id)));
        }
        this.renderSelections();
        this.selectionEditors.get(id)?.element.querySelector<HTMLElement>('select')?.focus();
    }

    private updateSelectionEdit(id: string): void {
        const editor = this.selectionEditors.get(id);
        const panel = Array.from(this.selectionList.querySelectorAll<HTMLElement>('.molui-selection-edit'))
            .find(node => node.dataset.selectionId === id);
        if (!editor || !panel) return;
        const count = editor.getAtoms().length;
        panel.querySelector('.molui-count').textContent = `${count.toLocaleString()} matching atom${count === 1 ? '' : 's'}`;
        panel.querySelector<HTMLButtonElement>('[data-action="save-selection"]').disabled = count === 0;
    }

    private finishSelectionEdit(id: string, save: boolean): void {
        const editor = this.selectionEditors.get(id);
        if (!editor || save && !editor.getAtoms().length) return;
        this.selectionEditors.delete(id);
        if (save) this.state.addSelection(editor.getSpec(), id);
        else this.renderSelections();
        this.focusControl(`${id}:edit`);
    }

    private representations(atoms: AtomSpec[]): Choice[] {
        const choices: Choice[] = [['stick', 'Sticks'], ['ballstick', 'Ball & stick'], ['sphere', 'Spheres'], ['line', 'Lines'], ['cross', 'Crosses']];
        if (isPolymer(atoms)) choices.unshift(['cartoon', 'Cartoon']);
        return [['none', 'None'], ...choices];
    }

    /** Shared color controls keep representation and surface choices consistent. */
    private colorControls(properties: UIColors, key: string, onChange?: (colors: UIColors) => void,
        spectrum = false, mixed = false) {
        const schemeName = typeof properties.colorscheme === 'string' ? properties.colorscheme : '';
        const knownScheme = schemeChoices.some(([name]) => name === schemeName);
        const mode = this.colorModes.get(key) === 'scheme' && knownScheme ? 'scheme' : mixed ? 'existing' : properties.color === 'spectrum' ? (spectrum ? 'spectrum' : 'existing')
            : properties.color !== undefined ? 'custom' : schemeName === 'chain' ? 'chain'
                : schemeName === 'ssJmol' ? 'secondary' : !properties.colorscheme || schemeName === 'default' ? 'element'
                    : knownScheme ? 'scheme' : 'existing';
        const coloring = select([['existing', 'Existing colors'], ['element', 'Element'], ['chain', 'Chain'],
            ['secondary', 'Secondary structure'], ['scheme', 'Color scheme'], ['custom', 'Solid color'],
            ...(spectrum ? [['spectrum', 'Spectrum']] as Choice[] : [])], mode);
        const scheme = select(schemeChoices, mode === 'scheme' && knownScheme ? schemeName : 'greenCarbon');
        const color = el('input'); color.type = 'color';
        color.value = properties.color !== undefined && properties.color !== 'spectrum'
            ? `#${CC.color(properties.color).getHex().toString(16).padStart(6, '0')}` : '#000066';
        const schemeField = field('Color scheme', scheme);
        const colorField = field('Solid color', color);
        const value = (): UIColors => {
            if (coloring.value === 'existing') {
                const result: UIColors = {};
                if (properties.color !== undefined) result.color = properties.color;
                if (properties.colorscheme !== undefined) result.colorscheme = properties.colorscheme;
                return result;
            }
            if (coloring.value === 'custom') return { color: color.value };
            if (coloring.value === 'spectrum') return { color: 'spectrum' };
            return { colorscheme: coloring.value === 'scheme' ? scheme.value
                : coloring.value === 'chain' ? 'chain' : coloring.value === 'secondary' ? 'ssJmol' : 'default' };
        };
        const update = () => {
            schemeField.hidden = coloring.value !== 'scheme';
            colorField.hidden = coloring.value !== 'custom';
        };
        for (const [name, input] of Object.entries({ coloring, scheme, color })) {
            input.dataset.focus = `${key}:${name}`;
            input.addEventListener('change', () => {
                this.colorModes.set(key, coloring.value);
                update(); onChange?.(value());
            });
        }
        update();
        return { fields: [field('Color by', coloring), schemeField, colorField], value, coloring };
    }

    private styleSpec(kind: string, colors: UIColors = { colorscheme: 'default' }, opacity = 1, radius = 0.25): AtomStyleSpec {
        // Spectrum is a cartoon color mode, not an atom color understood by other styles.
        if (kind !== 'cartoon' && colors.color === 'spectrum') colors = { colorscheme: 'default' };
        const base = { ...colors, opacity };
        if (kind === 'ballstick') return { stick: { ...base, radius }, sphere: { ...base, scale: 0.3 } };
        if (kind === 'stick') return { stick: { ...base, radius } };
        if (kind === 'sphere') return { sphere: base };
        if (kind === 'cartoon') return { cartoon: base };
        if (kind === 'cross') return { cross: base };
        return { line: base };
    }

    /** Settings belong to one representation, including imported combined styles. */
    private styleSettings(spec: AtomStyleSpec, selectionID: string, styleID: string, open: boolean,
        common: HTMLElement[]): HTMLDetailsElement {
        const key = `${selectionID}:${styleID}`;
        const details = section('Settings', open);
        details.classList.add('molui-geometry'); details.dataset.geometry = key;
        (details.firstElementChild as HTMLElement).dataset.focus = `${key}:geometry`;
        details.append(...common);
        const update = (changes: AtomStyleSpec) => this.state.addStyle({ ...spec, ...changes }, selectionID, styleID);
        const bind = (input: HTMLInputElement | HTMLSelectElement, name: string, change: () => void) => {
            input.dataset.focus = `${key}:${name}`;
            input.addEventListener('change', change);
        };
        const slider = (title: string, min: number, max: number, step: number, value: number, unit: string) => {
            const input = range(Math.min(min, value), Math.max(max, value), step, value);
            const text = () => `${title} · ${input.value}${unit}`;
            const label = field(text(), input);
            input.addEventListener('input', () => label.firstElementChild.textContent = text());
            return { input, label };
        };
        for (const kind of ['sphere', 'cross'] as const) {
            const properties = spec[kind];
            if (!properties) continue;
            const title = kind === 'sphere' ? 'Sphere' : 'Cross';
            const mode = select([['atomic', 'Atomic (van der Waals)'], ['fixed', 'Fixed radius']], properties.radius === undefined ? 'atomic' : 'fixed');
            const atomic = slider(`${title} scale`, 0.1, 2, 0.05, properties.scale ?? 1, '×');
            // The renderer multiplies explicit radius by scale. Display the effective
            // radius and remove scale when committing an absolute radius in Å.
            const fixed = slider(`${title} radius`, 0.05, 5, 0.05, (properties.radius ?? 1.5) * (properties.scale ?? 1), ' Å');
            atomic.label.hidden = mode.value !== 'atomic';
            fixed.label.hidden = mode.value !== 'fixed';
            bind(mode, `${kind}-radius-mode`, () => {
                const next = { ...properties };
                delete next.scale;
                if (mode.value === 'atomic') delete next.radius;
                else next.radius = Number(fixed.input.value);
                update({ [kind]: next });
            });
            bind(atomic.input, `${kind}-scale`, () => {
                const next = { ...properties, scale: Number(atomic.input.value) };
                delete next.radius;
                update({ [kind]: next });
            });
            bind(fixed.input, `${kind}-radius`, () => {
                const next = { ...properties, radius: Number(fixed.input.value) };
                delete next.scale;
                update({ [kind]: next });
            });
            details.append(field(`${title} size`, mode), atomic.label, fixed.label);
        }
        if (spec.cartoon) {
            const properties = spec.cartoon;
            const choices: Choice[] = [['default', 'Rectangle (default)'], ['oval', 'Oval'], ['parabola', 'Parabola'],
                ['edged', 'Edged'], ['trace', 'Backbone trace']];
            const current = properties.style === 'rectangle' ? 'default' : properties.style || 'default';
            if (!choices.some(([value]) => value === current)) choices.push([current, current]);
            const shape = select(choices, current);
            const trace = current === 'trace';
            bind(shape, 'cartoon-style', () => {
                const next = { ...properties, style: shape.value };
                // These flags modify secondary structure before the renderer enters
                // its trace branch, which cannot use them or undo that modification.
                if (shape.value === 'trace') { next.arrows = false; next.tubes = false; }
                update({ cartoon: next });
            });
            const thickness = slider(trace ? 'Trace radius' : 'Thickness',
                0.05, 2, 0.05, properties.thickness ?? 0.4, ' Å');
            bind(thickness.input, 'cartoon-thickness', () => update({ cartoon: { ...properties, thickness: Number(thickness.input.value) } }));
            const widthMode = select([['auto', 'By secondary structure'], ['custom', 'Custom width']], properties.width ? 'custom' : 'auto');
            const widthModeField = field('Strand width', widthMode);
            widthModeField.hidden = trace || !!properties.ribbon;
            const width = slider('Strand half-width', 0.1, 3, 0.05, properties.width || 1.3, ' Å');
            width.label.hidden = widthModeField.hidden || widthMode.value !== 'custom';
            bind(widthMode, 'cartoon-width-mode', () => {
                const next = { ...properties };
                if (widthMode.value === 'auto') delete next.width;
                else next.width = Number(width.input.value);
                update({ cartoon: next });
            });
            bind(width.input, 'cartoon-width', () => update({ cartoon: { ...properties, width: Number(width.input.value) } }));
            const options = el('div', 'molui-style-options');
            for (const [property, title] of [['ribbon', 'Uniform ribbon'], ['arrows', 'Sheet arrows'], ['tubes', 'Helix tubes']] as const) {
                const control = check(title, !!properties[property]);
                control.label.hidden = trace || property === 'arrows' && (!!properties.ribbon || !!properties.width);
                bind(control.input, `cartoon-${property}`, () => update({ cartoon: { ...properties, [property]: control.input.checked } }));
                options.append(control.label);
            }
            details.append(field('Cartoon style', shape), thickness.label, widthModeField, width.label, options);
        }
        return details;
    }

    private renderSelections(): void {
        const expanded = new Set(Array.from(this.selectionList.querySelectorAll<HTMLDetailsElement>('details[data-geometry]'))
            .filter(details => details.open).map(details => details.dataset.geometry));
        for (const id of this.selectionEditors.keys()) {
            if (!this.state.selections.has(id)) this.selectionEditors.delete(id);
        }
        for (const id of this.collapsedSelections) {
            if (!this.state.selections.has(id)) this.collapsedSelections.delete(id);
        }
        this.selectionList.replaceChildren();
        for (const entry of this.state.selections.values()) {
            const atoms = this.state.viewer.selectedAtoms(entry.spec);
            const card = el('div', 'molui-card molui-selection-card');
            card.dataset.selectionId = entry.id;
            const header = el('div', 'molui-selection-header');
            const content = el('div', 'molui-selection-content');
            content.id = `${this.panel.id}-selection-${encodeURIComponent(entry.id)}`;
            content.hidden = this.collapsedSelections.has(entry.id);
            const collapse = button('', () => {
                content.hidden = !content.hidden;
                if (content.hidden) this.collapsedSelections.add(entry.id);
                else this.collapsedSelections.delete(entry.id);
                updateCollapse();
                if (content.hidden && content.contains(document.activeElement)) collapse.focus();
            }, 'molui-selection-collapse');
            collapse.dataset.focus = `${entry.id}:collapse`;
            collapse.setAttribute('aria-controls', content.id);
            const updateCollapse = () => {
                collapse.setAttribute('aria-expanded', String(!content.hidden));
                collapse.setAttribute('aria-label', `${content.hidden ? 'Expand' : 'Collapse'} selection ${entry.name}`);
                collapse.title = `${content.hidden ? 'Expand' : 'Collapse'} selection`;
                edit.hidden = content.hidden || this.selectionEditors.has(entry.id);
            };
            const visible = check('', !entry.hidden);
            visible.label.classList.add('molui-selection-visibility');
            visible.input.setAttribute('aria-label', `Show selection ${entry.name}`);
            visible.input.title = 'Show all representations and labels for this selection';
            visible.input.dataset.focus = `${entry.id}:enabled`;
            visible.input.addEventListener('change', () => this.state.setSelectionVisible(entry.id, visible.input.checked));
            const name = el('input', 'molui-selection-name');
            name.type = 'text'; name.value = entry.name; name.maxLength = 160;
            name.setAttribute('aria-label', 'Selection name'); name.title = 'Edit selection name';
            name.dataset.focus = `${entry.id}:name`;
            const rename = () => {
                this.renamingSelection = true;
                try {
                    name.value = this.state.renameSelection(entry.id, name.value);
                    visible.input.setAttribute('aria-label', `Show selection ${name.value}`);
                    edit.setAttribute('aria-label', `Edit selection ${name.value}`);
                    updateCollapse();
                } finally { this.renamingSelection = false; }
            };
            name.addEventListener('change', rename);
            name.addEventListener('keydown', event => {
                if (event.key === 'Enter') { event.preventDefault(); rename(); }
                if (event.key === 'Escape') { event.stopPropagation(); name.value = entry.name; name.blur(); }
            });
            header.append(visible.label, name, collapse);
            const summaryRow = el('div', 'molui-selection-summary-row');
            const summary = el('p', 'molui-selection-summary',
                `${selectionSummary(entry.spec, this.state.modelTitles)} · ${atoms.length.toLocaleString()} atom${atoms.length === 1 ? '' : 's'}`);
            const edit = button('Edit', () => this.editSelection(entry.id));
            edit.setAttribute('aria-label', `Edit selection ${entry.name}`);
            edit.dataset.focus = `${entry.id}:edit`;
            updateCollapse();
            summaryRow.append(summary, edit);
            card.append(header, summaryRow, content);
            const draft = this.selectionEditors.get(entry.id);
            if (draft) {
                draft.refresh();
                const panel = el('div', 'molui-selection-edit');
                panel.dataset.selectionId = entry.id;
                panel.append(el('h4', '', 'Edit selection'), draft.element, el('p', 'molui-count'));
                const save = button('Save selection', () => this.finishSelectionEdit(entry.id, true), 'molui-primary');
                save.dataset.focus = `${entry.id}:save-selection`; save.dataset.action = 'save-selection';
                const cancel = button('Cancel', () => this.finishSelectionEdit(entry.id, false));
                cancel.dataset.focus = `${entry.id}:cancel-selection`;
                const draftActions = el('div', 'molui-actions'); draftActions.append(save, cancel);
                panel.append(draftActions); content.append(panel);
            }
            const labels = check('Residue labels', entry.labelsEnabled);
            labels.input.dataset.focus = `${entry.id}:labels`;
            labels.input.addEventListener('change', () => this.state.setResidueLabels(entry.id, labels.input.checked));
            const actions = el('div', 'molui-actions');
            actions.append(button('Focus', () => { this.state.viewer.zoomTo(entry.spec); this.state.viewer.render(); }),
                button('Remove selection', () => this.state.removeSelection(entry.id), 'molui-danger'));
            content.append(labels.label);
            const available = (exclude?: string): Choice[] => {
                const used = new Set(Array.from(entry.styles.values()).filter(style => style.id !== exclude)
                    .flatMap(style => Object.keys(style.spec)));
                return this.representations(atoms).filter(([kind]) => kind !== 'none' && (kind === 'ballstick'
                    ? !used.has('stick') && !used.has('sphere') : !used.has(kind)));
            };
            if (!entry.styles.size) {
                const representation = select(this.representations(atoms), 'none');
                representation.dataset.focus = `${entry.id}:representation`;
                representation.disabled = !atoms.length;
                representation.addEventListener('change', () => {
                    if (representation.value === 'none') return;
                    const id = this.state.addStyle(this.styleSpec(representation.value), entry.id);
                    this.focusControl(`${entry.id}:${id}:representation`);
                });
                const empty = el('div', 'molui-no-representation');
                empty.append(field('Representation', representation),
                    el('p', 'molui-muted', 'Use this selection for surfaces or residue labels.'));
                content.append(empty);
            }
            for (const style of entry.styles.values()) {
                const editor = el('div', 'molui-style');
                const styleHeader = el('div', 'molui-style-header');
                styleHeader.append(el('strong', '', 'Style'), button('Remove representation', () => {
                    this.state.removeStyle(entry.id, style.id);
                    this.focusControl(`${entry.id}:${entry.styles.size ? 'add-representation' : 'representation'}`);
                }, 'molui-danger'));
                const keys = Object.keys(style.spec);
                const kind = keys.length === 2 && keys.includes('stick') && keys.includes('sphere') ? 'ballstick'
                    : keys.length > 1 ? 'combined' : keys[0] || 'line';
                const choices: Choice[] = [['none', 'None'], ...available(style.id)];
                if (!choices.some(([key]) => key === kind)) choices.unshift([kind, kind === 'combined' ? `Combined (${keys.join(' + ')})` : kind]);
                const representation = select(choices, kind);
                const properties = (style.spec[keys[0]] || {}) as UIColors & { opacity?: number; colorfunc?: unknown };
                const mixedColors = !!properties.colorfunc || keys.some(key => style.spec[key]?.color !== properties.color
                    || style.spec[key]?.colorscheme !== properties.colorscheme || style.spec[key]?.colorfunc);
                const colors = this.colorControls(properties, `${entry.id}:${style.id}`, value => update('color', value), kind === 'cartoon', mixedColors);
                const opacity = range(0.05, 1, 0.05, properties.opacity ?? 1);
                const opacityField = field(`Opacity · ${Math.round(Number(opacity.value) * 100)}%`, opacity);
                opacity.addEventListener('input', () => opacityField.firstElementChild.textContent = `Opacity · ${Math.round(Number(opacity.value) * 100)}%`);
                const radius = range(0.1, 0.6, 0.05, style.spec.stick?.radius ?? 0.25);
                const radiusField = field(`Bond radius · ${radius.value} Å`, radius); radiusField.hidden = !keys.includes('stick');
                radius.addEventListener('input', () => radiusField.firstElementChild.textContent = `Bond radius · ${radius.value} Å`);
                for (const [name, input] of Object.entries({ representation, opacity, radius })) input.dataset.focus = `${entry.id}:${style.id}:${name}`;
                // Preserve imported representations and parameters not exposed by this editor.
                const update = (change: 'color' | 'opacity' | 'radius', colorStyle?: UIColors) => {
                    const spec: AtomStyleSpec = { ...style.spec };
                    for (const key of keys) {
                        if (!['cartoon', 'stick', 'sphere', 'line', 'cross'].includes(key)) continue;
                        const properties = { ...style.spec[key] };
                        if (change === 'opacity') properties.opacity = Number(opacity.value);
                        if (change === 'radius' && key === 'stick') properties.radius = Number(radius.value);
                        if (change === 'color' && colors.coloring.value !== 'existing') {
                            delete properties.color; delete properties.colorscheme; delete properties.colorfunc;
                            Object.assign(properties, colorStyle);
                        }
                        spec[key] = properties;
                    }
                    this.state.addStyle(spec, entry.id, style.id);
                };
                representation.addEventListener('change', () => {
                    if (representation.value === 'none') {
                        this.state.removeStyle(entry.id, style.id);
                        this.focusControl(`${entry.id}:${entry.styles.size ? 'add-representation' : 'representation'}`);
                        return;
                    }
                    this.state.addStyle(this.styleSpec(representation.value, colors.value(), Number(opacity.value), Number(radius.value)), entry.id, style.id);
                });
                opacity.addEventListener('change', () => update('opacity')); radius.addEventListener('change', () => update('radius'));
                const fields = el('div', 'molui-style-grid');
                fields.append(field('Representation', representation), ...colors.fields);
                const settings = this.styleSettings(style.spec, entry.id, style.id,
                    expanded.has(`${entry.id}:${style.id}`), [opacityField, radiusField]);
                editor.append(styleHeader, fields, settings);
                content.append(editor);
            }
            const addRow = el('div', 'molui-add-representation');
            const newRepresentation = select(available());
            newRepresentation.setAttribute('aria-label', 'New representation');
            newRepresentation.dataset.focus = `${entry.id}:add-representation`;
            const add = button('Add', () => {
                if (!newRepresentation.value) return;
                const id = this.state.addStyle(this.styleSpec(newRepresentation.value), entry.id);
                this.selectionList.querySelector<HTMLElement>(`[data-focus="${entry.id}:${id}:representation"]`)?.focus();
            });
            add.setAttribute('aria-label', 'Add representation');
            add.title = 'Add representation';
            add.disabled = newRepresentation.disabled = !newRepresentation.options.length || !atoms.length;
            addRow.append(field('New representation', newRepresentation), add);
            if (entry.styles.size) content.append(addRow);
            content.append(actions); this.selectionList.append(card);
            if (draft) this.updateSelectionEdit(entry.id);
        }
    }

    private updateSurfaceSelection(): void {
        const shown = this.state.selections.get(this.surfaceSelection.value)?.spec || {};
        const calculation = this.surfaceCalculation.value === 'shown' ? shown
            : this.state.selections.get(this.surfaceCalculation.value)?.spec || {};
        const shownAtoms = this.state.viewer.selectedAtoms(shown);
        const calculationAtoms = new Set(this.state.viewer.selectedAtoms(calculation));
        const error = !shownAtoms.length ? 'Choose atoms to show before adding a surface.'
            : !calculationAtoms.size ? 'Choose atoms to calculate the surface.'
                : shownAtoms.some(atom => !calculationAtoms.has(atom)) ? 'Calculation atoms must include all shown atoms.' : '';
        this.surfaceButton.disabled = this.surfacePending || !!error;
        const text = error || 'Calculation atoms shape the surface; only the shown atoms’ surface is displayed.';
        if (this.surfaceValidation.textContent !== text) this.surfaceValidation.textContent = text;
        this.surfaceValidation.dataset.error = String(!!error);
    }

    private renderSurfaces(): void {
        const selected = this.surfaceSelection.value;
        const calculation = this.surfaceCalculation.value;
        this.surfaceSelection.replaceChildren(new Option('All atoms', ''));
        this.surfaceCalculation.replaceChildren(new Option('Same as shown', 'shown'), new Option('All atoms', ''));
        this.state.selections.forEach(entry => {
            this.surfaceSelection.add(new Option(entry.name, entry.id));
            this.surfaceCalculation.add(new Option(entry.name, entry.id));
        });
        if (this.state.selections.has(selected)) this.surfaceSelection.value = selected;
        if (calculation === '' || calculation === 'shown' || this.state.selections.has(calculation)) this.surfaceCalculation.value = calculation;
        this.updateSurfaceSelection();
        this.surfaceList.replaceChildren();
        this.state.surfaces.forEach(entry => {
            const card = el('div', 'molui-card molui-surface-card');
            const colors = this.colorControls(entry.style, `${entry.id}:surface`, value => this.state.updateSurface(entry.id, value));
            card.append(el('strong', '', entry.name),
                el('p', 'molui-muted', `Calculated using: ${selectionSummary(entry.calculationSpec, this.state.modelTitles)}`),
                ...colors.fields, button('Remove surface', () => this.state.removeSurface(entry.id), 'molui-danger'));
            this.surfaceList.append(card);
        });
    }

    private viewControls(): HTMLDetailsElement {
        const details = section('View & export', false);
        const actions = el('div', 'molui-actions');
        actions.append(button('Fit all', () => { this.state.viewer.zoomTo(); this.state.viewer.render(); }),
            button('Save PNG', () => this.download(this.state.viewer.pngURI(), 'molecule.png')));
        const background = el('input'); background.type = 'color'; background.value = '#ffffff';
        background.addEventListener('input', () => { this.state.viewer.setBackgroundColor(background.value, 1); this.state.viewer.render(); });
        const projection = select([['perspective', 'Perspective'], ['orthographic', 'Orthographic']]);
        projection.addEventListener('change', () => { this.state.viewer.setProjection(projection.value); this.state.viewer.render(); });
        details.append(actions, field('Background', background), field('Projection', projection),
            el('p', 'molui-muted', 'Drag to rotate. Scroll or pinch to zoom. Shift-drag or use two fingers to move the molecule.'));
        return details;
    }

    private download(url: string, name: string): void {
        const link = el('a'); link.href = url; link.download = name; link.click();
    }
}
