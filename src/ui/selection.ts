import type { AtomSelectionSpec, AtomSpec } from '../specs';
import type { GLViewer } from '../GLViewer';

type FilterKey = 'chain' | 'resn' | 'resi' | 'elem';
type Group = '' | 'polymer' | 'ligand' | 'water' | 'hetero';
type FilterValue = string | number;
const filters: [FilterKey, string][] = [['chain', 'Chain'], ['resn', 'Residue name'], ['resi', 'Residue number'], ['elem', 'Element']];
const water = ['HOH', 'WAT', 'H2O'];

export interface SelectionEditorConfig {
    viewer: Pick<GLViewer, 'selectedAtoms'>;
    models: () => { id: number; title: string }[];
    key: string;
    spec?: AtomSelectionSpec;
    onChange?: (spec: AtomSelectionSpec, atoms: AtomSpec[]) => void;
}

/** Copy selection data without cloning model objects or retaining selection caches. */
function copy<T>(value: T): T {
    if (Array.isArray(value)) return value.map(copy) as T;
    if (value && typeof value === 'object' && Object.getPrototypeOf(value) === Object.prototype) {
        const result = {};
        Object.keys(value).filter(key => !key.startsWith('__cache')).forEach(key => result[key] = copy(value[key]));
        return result as T;
    }
    return value;
}

function isWater(spec: AtomSelectionSpec): boolean {
    return !!spec && Object.keys(spec).filter(key => !key.startsWith('__cache')).length === 1 && Array.isArray(spec.or)
        && spec.or.length === water.length && water.every(name => spec.or.some(part =>
            Object.keys(part).filter(key => !key.startsWith('__cache')).length === 1 && part.resn === name));
}

/** Remove only a recognized atom-group clause, leaving imported predicates intact. */
function splitGroup(spec: AtomSelectionSpec): { spec: AtomSelectionSpec; group: Group } {
    const base = copy(spec);
    let group: Group = '';
    if (isWater({ or: base.or })) {
        group = 'water'; delete base.or;
    } else if (base.hetflag === true && isWater(base.not)) {
        group = 'ligand'; delete base.hetflag; delete base.not;
    } else if (base.hetflag === false) {
        group = 'polymer'; delete base.hetflag;
    } else if (base.hetflag === true) {
        group = 'hetero'; delete base.hetflag;
    } else if (Array.isArray(base.and)) {
        for (let i = 0; i < base.and.length; i++) {
            const part = splitGroup(base.and[i]);
            if (part.group && Object.keys(part.spec).length === 0) {
                group = part.group;
                base.and.splice(i, 1);
                if (!base.and.length) delete base.and;
                break;
            }
        }
    }
    return { spec: base, group };
}

function groupSpec(group: Group): AtomSelectionSpec {
    const waterSpec = { or: water.map(resn => ({ resn })) };
    if (group === 'water') return waterSpec;
    if (group === 'ligand') return { hetflag: true, not: waterSpec };
    if (group === 'polymer') return { hetflag: false };
    if (group === 'hetero') return { hetflag: true };
    return {};
}

function values(value: unknown): FilterValue[] {
    return value === undefined ? [] : (Array.isArray(value) ? value : [value]) as FilterValue[];
}

type ResidueInterval = { from: number; to: number };
type ResidueDraft = { from: string; to: string; editing?: ResidueInterval };

function residueIntervals(items: FilterValue[]): { ranges: ResidueInterval[]; other: FilterValue[] } {
    const ranges: ResidueInterval[] = [];
    const other: FilterValue[] = [];
    for (const item of items) {
        const match = typeof item === 'string' ? item.match(/^\s*(-?\d+)\s*-\s*(-?\d+)\s*$/) : undefined;
        const from = match ? Number(match[1]) : typeof item === 'number' || /^-?\d+$/.test(String(item)) ? Number(item) : NaN;
        const to = match ? Number(match[2]) : from;
        if (Number.isSafeInteger(from) && Number.isSafeInteger(to) && from <= to) ranges.push({ from, to });
        else other.push(item);
    }
    ranges.sort((a, b) => a.from - b.from || a.to - b.to);
    const merged: ResidueInterval[] = [];
    for (const interval of ranges) {
        const last = merged[merged.length - 1];
        if (last && interval.from <= last.to + 1) last.to = Math.max(last.to, interval.to);
        else merged.push({ ...interval });
    }
    return { ranges: merged, other };
}

function intervalText(interval: ResidueInterval): string {
    return interval.from === interval.to ? String(interval.from) : `${interval.from}–${interval.to}`;
}

function numberList(items: FilterValue[]): string {
    const { ranges, other } = residueIntervals(items);
    return [...ranges.map(intervalText), ...other.map(String)].join(', ');
}

/** A readable description of the criteria, independent of the selection's name. */
export function selectionSummary(spec: AtomSelectionSpec, modelTitles?: ReadonlyMap<number, string>): string {
    const parts: string[] = [];
    const { spec: remaining, group } = splitGroup(spec);
    if (spec.model !== undefined) {
        const models = Array.isArray(spec.model) ? spec.model : [spec.model];
        parts.push(models.map(model => {
            const id = typeof model === 'number' ? model : model.getID();
            return modelTitles?.get(id) || `Model ${id + 1}`;
        }).join(', '));
    }
    if (group) parts.push({ polymer: 'Protein / nucleic acid', ligand: 'Hetero atoms (no water)', water: 'Water', hetero: 'Hetero atoms' }[group]);
    if (spec.chain !== undefined) parts.push(`${values(spec.chain).length > 1 ? 'Chains' : 'Chain'} ${values(spec.chain).map(value => value === '' ? '(blank)' : value).join(', ') || '(none)'}`);
    if (spec.resn !== undefined) parts.push(`Residue names ${values(spec.resn).join(', ') || '(none)'}`);
    if (spec.resi !== undefined) parts.push(`Residues ${numberList(values(spec.resi)) || '(none)'}`);
    if (spec.elem !== undefined) parts.push(`Elements ${values(spec.elem).join(', ') || '(none)'}`);
    for (const key of ['model', ...filters.map(([key]) => key)]) delete remaining[key];
    if (Object.keys(remaining).length) parts.push('Additional constraints');
    return parts.join(' · ') || 'All atoms';
}

function element<K extends keyof HTMLElementTagNameMap>(tag: K, className = '', text?: string): HTMLElementTagNameMap[K] {
    const node = document.createElement(tag);
    node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
}

function field(title: string, control: HTMLElement): HTMLLabelElement {
    const label = element('label', 'molui-field');
    label.append(element('span', '', title), control);
    return label;
}

/** Searchable checklists keep multi-value selection usable with touch and keyboards. */
export class SelectionEditor {
    public readonly element = element('div', 'molui-selection-editor');
    private spec: AtomSelectionSpec;
    private readonly search = new Map<FilterKey, string>();
    private residueDraft?: ResidueDraft;
    private residueError = '';

    constructor(private readonly config: SelectionEditorConfig) {
        this.spec = copy(config.spec || {});
        this.refresh();
    }

    public getSpec(): AtomSelectionSpec { return copy(this.spec); }
    public getAtoms(): AtomSpec[] { return this.config.viewer.selectedAtoms(this.getSpec()); }

    public setSpec(spec: AtomSelectionSpec): void {
        this.spec = copy(spec);
        this.residueDraft = undefined; this.residueError = '';
        this.refresh();
    }

    public refresh(): void {
        const expanded = new Set(Array.from(this.element.querySelectorAll<HTMLDetailsElement>('details[open]')).map(details => details.dataset.field));
        const active = document.activeElement as HTMLInputElement;
        const focus = this.element.contains(active) ? active.dataset.focus : undefined;
        const caret = active?.type === 'search' ? [active.selectionStart, active.selectionEnd] : undefined;
        this.element.replaceChildren();
        const grid = element('div', 'molui-filters');
        this.element.append(grid);
        const model = element('select');
        model.dataset.focus = `${this.config.key}:model`;
        model.add(new Option('All models', ''));
        this.config.models().forEach(item => model.add(new Option(item.title, String(item.id))));
        if (this.spec.model !== undefined) {
            if (typeof this.spec.model === 'number' && Array.from(model.options).some(option => option.value === String(this.spec.model))) model.value = String(this.spec.model);
            else { model.add(new Option('Existing models', 'existing')); model.value = 'existing'; }
        }
        model.addEventListener('change', () => {
            if (model.value === '') delete this.spec.model;
            else if (model.value !== 'existing') this.spec.model = Number(model.value);
            this.changed();
        });
        const group = element('select');
        group.dataset.focus = `${this.config.key}:group`;
        [['', 'All atoms'], ['polymer', 'Protein / nucleic acid'], ['ligand', 'Hetero atoms (no water)'], ['water', 'Water'], ['hetero', 'Hetero atoms']]
            .forEach(([key, title]) => group.add(new Option(title, key)));
        group.value = splitGroup(this.spec).group;
        group.addEventListener('change', () => {
            const base = splitGroup(this.spec).spec;
            const chosen = groupSpec(group.value as Group);
            if (Object.keys(chosen).some(key => key in base)) base.and = [...(base.and || []), chosen];
            else Object.assign(base, chosen);
            this.spec = base;
            this.changed();
        });
        grid.append(field('Model', model), field('Atoms', group));

        const prefix = this.getSpec();
        filters.forEach(([key]) => delete prefix[key]);
        for (const [key, title] of filters) {
            const available = Array.from(new Set(this.config.viewer.selectedAtoms(copy(prefix)).map(atom => atom[key])
                .filter(value => value !== undefined && value !== null).map(String))).sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
            if (key === 'resi') {
                grid.append(this.residueControls(available.map(Number).filter(Number.isSafeInteger), expanded.has(key)));
                if (this.spec.resi !== undefined) prefix.resi = copy(this.spec.resi);
                continue;
            }
            const chosen = values(this.spec[key]).map(String);
            const options = Array.from(new Set([...available, ...chosen]));
            const details = element('details', 'molui-multiselect');
            details.dataset.field = key;
            details.open = expanded.has(key);
            const summary = element('summary');
            summary.dataset.focus = `${this.config.key}:${key}:toggle`;
            const description = !chosen.length ? (this.spec[key] === undefined ? 'All' : 'None')
                : chosen.length > 3 ? `${chosen.length} selected` : chosen.map(value => value || '(blank)').join(', ');
            summary.append(element('span', 'molui-multiselect-label', title), element('span', 'molui-multiselect-value', description));
            summary.title = `${title}: ${chosen.length ? chosen.map(value => value || '(blank)').join(', ') : description}`;
            details.addEventListener('keydown', event => {
                if (event.key === 'Escape' && details.open) {
                    event.stopPropagation(); details.open = false; summary.focus();
                }
            });
            const body = element('div', 'molui-multiselect-body');
            const search = element('input', 'molui-multiselect-search');
            search.type = 'search'; search.placeholder = `Find ${title.toLowerCase()}`;
            search.setAttribute('aria-label', `Find ${title.toLowerCase()}`);
            search.dataset.focus = `${this.config.key}:${key}:search`;
            search.value = this.search.get(key) || '';
            search.hidden = options.length <= 8 && !search.value;
            const actions = element('div', 'molui-multiselect-actions');
            const all = element('button', '', 'All'); all.type = 'button';
            all.setAttribute('aria-label', `All ${title.toLowerCase()} values`);
            all.dataset.focus = `${this.config.key}:${key}:all`;
            all.addEventListener('click', () => { delete this.spec[key]; this.changed(); });
            actions.append(all);
            const list = element('div', 'molui-multiselect-options');
            list.setAttribute('role', 'group');
            list.setAttribute('aria-label', `${title} values`);
            for (const value of options) {
                const label = element('label', `molui-check${available.includes(value) ? '' : ' molui-unavailable'}`);
                const input = element('input'); input.type = 'checkbox'; input.value = value; input.checked = chosen.includes(value);
                input.dataset.value = value;
                input.dataset.focus = `${this.config.key}:${key}:value:${value}`;
                const text = value || '(blank)';
                label.append(input, element('span', '', `${text}${available.includes(value) ? '' : ' (no matches)'}`));
                if (!available.includes(value)) label.title = 'No matching atoms with the preceding filters';
                input.addEventListener('change', () => {
                    const next = chosen.filter(item => item !== value);
                    if (input.checked) next.push(value);
                    if (!next.length) delete this.spec[key];
                    else Object.assign(this.spec, { [key]: next.length === 1 ? next[0] : next });
                    this.changed();
                });
                list.append(label);
            }
            const empty = element('p', 'molui-muted', 'No matching values.');
            const filter = () => {
                this.search.set(key, search.value);
                const query = search.value.trim().toLowerCase();
                let shown = 0;
                for (const label of Array.from(list.children) as HTMLLabelElement[]) {
                    label.hidden = !label.textContent.toLowerCase().includes(query);
                    if (!label.hidden) shown++;
                }
                empty.hidden = shown > 0;
            };
            search.addEventListener('input', filter);
            filter();
            body.append(search, actions, list, empty);
            details.append(summary, body);
            grid.append(details);
            if (this.spec[key] !== undefined) Object.assign(prefix, { [key]: copy(this.spec[key]) });
        }
        if (focus) {
            const target = Array.from(this.element.querySelectorAll<HTMLElement>('[data-focus]')).find(node => node.dataset.focus === focus);
            target?.focus({ preventScroll: true });
            if (caret && target?.tagName === 'INPUT') (target as HTMLInputElement).setSelectionRange(caret[0], caret[1]);
        }
    }

    /** Commit whole intervals instead of allocating one checkbox per residue. */
    private residueControls(available: number[], open: boolean): HTMLDetailsElement {
        const key = `${this.config.key}:resi`;
        const details = element('details', 'molui-multiselect molui-residue-ranges');
        details.dataset.field = 'resi'; details.open = open;
        let min: number | undefined; let max: number | undefined;
        for (const value of available) {
            min = min === undefined ? value : Math.min(min, value);
            max = max === undefined ? value : Math.max(max, value);
        }
        if (min !== undefined) { details.dataset.min = String(min); details.dataset.max = String(max); }
        const { ranges, other } = residueIntervals(values(this.spec.resi));
        const description = this.spec.resi === undefined ? 'All' : numberList(values(this.spec.resi)) || 'None';
        const summary = element('summary'); summary.dataset.focus = `${key}:toggle`;
        summary.append(element('span', 'molui-multiselect-label', 'Residue number'),
            element('span', 'molui-multiselect-value', description));
        summary.title = `Residue number: ${description}`;
        details.addEventListener('keydown', event => {
            if (event.key === 'Escape' && details.open) {
                event.stopPropagation(); details.open = false; summary.focus();
            }
        });
        const body = element('div', 'molui-multiselect-body');
        body.append(element('p', 'molui-muted molui-residue-available', available.length
            ? `Available: ${intervalText({ from: min, to: max })} · ${available.length} residue number${available.length === 1 ? '' : 's'}`
            : 'No residue numbers match the preceding filters.'));
        const chips = element('div', 'molui-residue-chips');
        const save = (next: ResidueInterval[], extras = other) => {
            const normalized = residueIntervals([...next.map(interval => interval.from === interval.to ? interval.from : `${interval.from}-${interval.to}`), ...extras]);
            const result = [...normalized.ranges.map(interval => interval.from === interval.to ? interval.from : `${interval.from}-${interval.to}`), ...normalized.other];
            this.spec.resi = (result.length === 1 ? result[0] : result) as AtomSelectionSpec['resi'];
            this.residueDraft = undefined; this.residueError = '';
            this.changed();
        };
        const focus = (suffix: string) => Array.from(this.element.querySelectorAll<HTMLElement>('[data-focus]'))
            .find(node => node.dataset.focus === `${key}:${suffix}`)?.focus();
        const chip = (text: string, suffix: string, remove: () => void, edit?: () => void) => {
            const node = element('span', 'molui-residue-chip');
            if (edit) {
                const label = element('button', '', text); label.type = 'button';
                label.setAttribute('aria-label', `Edit residues ${text}`); label.dataset.focus = `${key}:edit:${suffix}`;
                label.addEventListener('click', edit); node.append(label);
            } else node.append(element('span', '', text));
            const close = element('button', '', '×'); close.type = 'button';
            close.setAttribute('aria-label', `Remove residues ${text}`); close.dataset.focus = `${key}:remove:${suffix}`;
            close.addEventListener('click', () => { remove(); focus('from'); }); node.append(close);
            chips.append(node); return node;
        };
        ranges.forEach((interval, index) => {
            const node = chip(intervalText(interval), String(index), () => save(ranges.filter((_, i) => i !== index)), () => {
                this.residueDraft = { from: String(interval.from), to: String(interval.to), editing: interval };
                this.residueError = ''; this.refresh(); focus('from');
            });
            node.dataset.from = String(interval.from); node.dataset.to = String(interval.to);
        });
        other.forEach((value, index) => {
            const node = chip(String(value), `other:${index}`, () => save(ranges, other.filter((_, i) => i !== index)));
            node.dataset.value = String(value);
        });
        body.append(chips);
        const draft: ResidueDraft = this.residueDraft || { from: min === undefined ? '' : String(min), to: max === undefined ? '' : String(max) };
        const from = element('input'); from.type = 'number'; from.step = '1'; from.value = draft.from;
        const to = element('input'); to.type = 'number'; to.step = '1'; to.value = draft.to;
        from.dataset.focus = `${key}:from`; to.dataset.focus = `${key}:to`;
        const message = element('p', 'molui-residue-error', this.residueError);
        message.setAttribute('role', 'status'); message.hidden = !this.residueError;
        const valid = () => from.value.trim() !== '' && to.value.trim() !== ''
            && Number.isSafeInteger(Number(from.value)) && Number.isSafeInteger(Number(to.value)) && Number(from.value) <= Number(to.value);
        const readDraft = () => {
            this.residueDraft = { from: from.value, to: to.value, editing: draft.editing };
            if (valid()) { this.residueError = ''; message.hidden = true; }
        };
        const endpoints = element('div', 'molui-residue-endpoints');
        const sliders: HTMLInputElement[] = [];
        const syncSliders = () => {
            const bounds = [min, max, ...[from, to].filter(input => input.value.trim() !== '').map(input => Number(input.value))]
                .filter(value => Number.isSafeInteger(value));
            const lower = bounds.length ? Math.min(...bounds) : 0;
            const upper = bounds.length ? Math.max(...bounds) : 0;
            sliders.forEach((slider, index) => {
                slider.min = String(lower); slider.max = String(upper); slider.disabled = lower === upper;
                slider.value = [from, to][index].value || String(index === 0 ? lower : upper);
            });
        };
        for (const [name, title, input] of [['from', 'From', from], ['to', 'To', to]] as const) {
            const slider = element('input'); slider.type = 'range'; slider.step = '1';
            slider.setAttribute('aria-label', `${title} residue number`); slider.dataset.focus = `${key}:${name}-slider`;
            slider.addEventListener('input', () => { input.value = slider.value; readDraft(); syncSliders(); });
            input.addEventListener('input', () => { readDraft(); syncSliders(); });
            const control = element('div', 'molui-residue-endpoint'); control.append(field(title, input), slider);
            endpoints.append(control); sliders.push(slider);
        }
        syncSliders();
        const commit = () => {
            readDraft();
            if (!valid()) {
                this.residueError = 'Enter whole numbers with From less than or equal to To.';
                message.textContent = this.residueError; message.hidden = false;
                return;
            }
            const previous = ranges.filter(interval => !draft.editing || interval.from !== draft.editing.from || interval.to !== draft.editing.to);
            save([...previous, { from: Number(from.value), to: Number(to.value) }]);
            focus('from');
        };
        for (const input of [from, to]) input.addEventListener('keydown', event => {
            if (event.key === 'Enter') { event.preventDefault(); event.stopPropagation(); commit(); }
        });
        const actions = element('div', 'molui-residue-actions');
        const add = element('button', 'molui-primary', draft.editing ? 'Update range' : 'Add range'); add.type = 'button'; add.dataset.focus = `${key}:add`;
        add.addEventListener('click', commit);
        const all = element('button', '', 'All'); all.type = 'button'; all.dataset.focus = `${key}:all`;
        all.setAttribute('aria-label', 'All residue number values');
        all.addEventListener('click', () => {
            delete this.spec.resi; this.residueDraft = undefined; this.residueError = ''; this.changed(); focus('all');
        });
        actions.append(add, all);
        if (draft.editing) {
            const cancel = element('button', '', 'Cancel edit'); cancel.type = 'button'; cancel.dataset.focus = `${key}:cancel`;
            cancel.addEventListener('click', () => { this.residueDraft = undefined; this.residueError = ''; this.refresh(); focus('from'); });
            actions.append(cancel);
        }
        body.append(endpoints, element('p', 'molui-muted', draft.editing ? 'Update range to apply these endpoints.' : 'Add range to apply these endpoints. Add more ranges for separate intervals.'), message, actions);
        details.append(summary, body);
        return details;
    }

    private changed(): void {
        this.refresh();
        this.config.onChange?.(this.getSpec(), this.getAtoms());
    }
}
