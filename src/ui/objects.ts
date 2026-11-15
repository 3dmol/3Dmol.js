import type { StateManager } from './state';
import type { UIConfig } from './ui';
import { CC } from '../colors';

function node<K extends keyof HTMLElementTagNameMap>(tag: K, text = ''): HTMLElementTagNameMap[K] {
    const element = document.createElement(tag); element.textContent = text; return element;
}
function field(title: string, input: HTMLElement): HTMLLabelElement {
    const label = node('label'); label.className = 'molui-field';
    label.append(node('span', title), input); return label;
}
function choice(values: [string, string][]): HTMLSelectElement {
    const select = node('select'); values.forEach(([key, title]) => select.add(new Option(title, key))); return select;
}
function number(value: number, min?: number): HTMLInputElement {
    const input = node('input'); input.type = 'number'; input.step = 'any'; input.value = String(value);
    if (min !== undefined) input.min = String(min); return input;
}
function section(title: string): HTMLDetailsElement {
    const result = node('details'); result.className = 'molui-section'; result.append(node('summary', title)); return result;
}

/** Map and annotation controls share the manager's change notifications. */
export class SceneObjectControls {
    readonly volumes = section('Volumetric data');
    readonly measurements = section('Measurements');
    private readonly volumeList = node('div');
    private readonly distanceList = node('div');
    private readonly measureToggle = node('input');
    private readonly clearDistances = this.action('Clear distances', () => this.state.clearDistances(), true);
    private readonly distanceHint = node('p');
    private enabled = true;
    private generation = 0;
    private abort?: AbortController;

    constructor(private readonly state: StateManager, config: UIConfig,
        private readonly error: (message?: string) => void) {
        this.volumeList.className = 'molui-list molui-volumes';
        this.volumes.append(this.volumeList, this.volumeForm(config));
        this.distanceControls();
    }

    setEnabled(enabled: boolean): void {
        this.enabled = enabled;
        if (!enabled) { this.generation++; this.abort?.abort(); this.state.setDistancePicking(false); }
    }

    destroy(): void { this.setEnabled(false); }

    private action(title: string, callback: () => void, danger = false): HTMLButtonElement {
        const button = node('button', title); button.type = 'button';
        if (danger) button.className = 'molui-danger';
        button.addEventListener('click', () => {
            if (!this.enabled) return;
            try { callback(); this.error(); } catch (error) { this.error(String(error.message || error)); }
        });
        return button;
    }

    private distanceControls(): void {
        this.distanceList.className = 'molui-list molui-distances';
        this.distanceHint.className = 'molui-muted'; this.distanceHint.setAttribute('aria-live', 'polite');
        this.measureToggle.type = 'checkbox'; this.measureToggle.dataset.focus = 'distance-picking';
        const enabled = node('label'); enabled.className = 'molui-check';
        enabled.append(this.measureToggle, node('span', 'Enable measuring'));
        this.measureToggle.addEventListener('change', () => {
            if (this.enabled) this.state.setDistancePicking(this.measureToggle.checked);
        });
        const actions = node('div'); actions.className = 'molui-row molui-distance-actions'; actions.append(this.clearDistances);
        this.measurements.append(enabled, this.distanceHint, this.distanceList, actions);
    }

    private volumeForm(config: UIConfig): HTMLElement {
        const form = node('form'); form.className = 'molui-create molui-volume-create';
        const source = choice([...(config.allowLocalFiles === false ? [] : [['file', 'Local file']] as [string, string][]), ['url', 'URL']]);
        const file = node('input'); file.type = 'file'; file.accept = '.dx,.cube,.cub,.vasp,.ccp4,.map,.mrc,.gz';
        const url = node('input'); url.type = 'url'; url.placeholder = 'https://…/density.dx';
        const format = choice([['auto', 'Detect from filename'], ['dx', 'OpenDX'], ['cube', 'Gaussian cube'], ['ccp4', 'CCP4 / MRC'], ['vasp', 'VASP']]);
        const mode = choice([['isosurface', 'Isosurface'], ['volume', 'Volume rendering (WebGL 2)']]);
        const color = node('input'); color.type = 'color'; color.value = '#3366cc';
        const fileField = field('Map file', file), urlField = field('Map URL', url);
        const submit = node('button', 'Load map'); submit.type = 'submit'; submit.className = 'molui-primary';
        const update = () => {
            fileField.hidden = source.value !== 'file'; urlField.hidden = source.value === 'file';
            url.disabled = source.value !== 'url';
        };
        source.addEventListener('change', update); update();
        file.addEventListener('change', () => {
            if (source.value === 'file' && file.files?.length) form.requestSubmit();
        });
        form.append(node('h3', 'New map'), field('Map source', source), ...(config.allowLocalFiles === false ? [] : [fileField]), urlField,
            field('Map format', format), field('Display as', mode), field('Map color', color), submit);
        form.addEventListener('submit', async event => {
            event.preventDefault(); if (!this.enabled || submit.disabled) return;
            const token = ++this.generation;
            this.abort?.abort(); this.abort = new AbortController(); submit.disabled = true;
            try {
                const selected = file.files?.[0];
                if (source.value === 'file' && !selected) throw new Error('Choose a map file first.');
                const address = source.value === 'url' ? new URL(url.value) : undefined;
                if (address && !['http:', 'https:'].includes(address.protocol)) throw new Error('Map URLs must use http or https.');
                const name = selected && source.value === 'file' ? selected.name : address.pathname.split('/').pop() || 'Volume';
                const extension = name.toLowerCase().match(/\.([a-z0-9]+)(\.gz)?$/);
                let kind = format.value === 'auto' ? extension?.[1] : format.value;
                kind = ({ map: 'ccp4', mrc: 'ccp4', cub: 'cube' } as Record<string, string>)[kind] || kind;
                if (!['dx', 'cube', 'ccp4', 'vasp'].includes(kind)) throw new Error('Choose the map format; the filename is ambiguous.');
                if (/\.gz$/i.test(name)) kind += '.gz';
                const binary = kind === 'ccp4' || kind.endsWith('.gz');
                let data: string | ArrayBuffer;
                if (source.value === 'file') {
                    data = await new Promise<string | ArrayBuffer>((resolve, reject) => {
                        const reader = new FileReader();
                        reader.onload = () => resolve(reader.result as string | ArrayBuffer);
                        reader.onerror = () => reject(new Error('The map file could not be read.'));
                        if (binary) reader.readAsArrayBuffer(selected); else reader.readAsText(selected);
                    });
                } else {
                    const response = await fetch(address.href, { signal: this.abort.signal });
                    if (!response.ok) throw new Error(`Map download failed (HTTP ${response.status}).`);
                    data = binary ? await response.arrayBuffer() : await response.text();
                }
                if (!this.enabled || token !== this.generation) return;
                this.state.addVolume(data, kind, { color: color.value, opacity: 0.8 }, name, mode.value as 'isosurface' | 'volume');
                this.state.viewer.zoomTo(); this.state.viewer.render(); this.error();
            } catch (error) {
                if (this.enabled && token === this.generation) this.error(`Unable to load map. ${error.message || error}`);
            } finally { submit.disabled = false; }
        });
        return form;
    }

    private styleFields(id: string, spec, update: (patch) => void): HTMLElement[] {
        const color = node('input'); color.type = 'color';
        const value = Array.isArray(spec.color) ? spec.color[0] : spec.color;
        const hex = typeof value?.getHex === 'function' ? value.getHex() : CC.color(value ?? '#3366cc').getHex();
        color.value = '#' + hex.toString(16).padStart(6, '0');
        const opacity = number(spec.opacity ?? 1, 0); opacity.max = '1'; opacity.step = '0.05';
        const visible = node('input'); visible.type = 'checkbox'; visible.checked = !spec.hidden;
        for (const [key, input] of Object.entries({ color, opacity, hidden: visible })) {
            input.dataset.focus = `${id}:${key}`;
            input.addEventListener('change', () => {
                if (!this.enabled) return;
                try {
                    if (!input.checkValidity()) throw new Error('Opacity must be between zero and one.');
                    update({ [key]: key === 'color' ? color.value : key === 'opacity' ? Number(opacity.value) : !visible.checked });
                    this.error();
                } catch (error) { this.error(String(error.message || error)); }
            });
        }
        return [field('Color', color), field('Opacity', opacity), field('Visible', visible)];
    }

    refresh(): void {
        this.volumeList.replaceChildren();
        const expanded = new Set(Array.from(this.distanceList.querySelectorAll<HTMLDetailsElement>('details[open]'))
            .map(card => card.dataset.distanceId));
        this.distanceList.replaceChildren();
        for (const entry of this.state.distances.values()) {
            const card = node('details'); card.className = 'molui-card molui-measurement-card'; card.dataset.distanceId = entry.id;
            card.open = expanded.has(entry.id);
            const summary = node('summary'); summary.dataset.focus = `${entry.id}:collapse`;
            const heading = node('span'); heading.className = 'molui-measurement-summary';
            const distance = `${entry.distance.toFixed(3)} Å`, name = entry.labelText.trim() || distance;
            const caption = node('span'); caption.append(node('strong', name));
            if (name !== distance) {
                const value = node('span', `(${distance})`); value.className = 'molui-measurement-distance';
                caption.append(' ', value);
            }
            const title = node('span', entry.name); title.className = 'molui-muted';
            heading.append(caption, title); summary.append(heading);
            const coordinates = (point: typeof entry.start) => [point.x, point.y, point.z].map(value => value.toFixed(3)).join(', ');
            const start = node('p', `Start (Å): ${coordinates(entry.start)}`), end = node('p', `End (Å): ${coordinates(entry.end)}`);
            start.className = end.className = 'molui-muted';
            const color = node('input'); color.type = 'color'; color.value = entry.color; color.dataset.focus = `${entry.id}:color`;
            color.addEventListener('change', () => {
                if (!this.enabled) return;
                try { this.state.updateDistance(entry.id, { color: color.value }); this.error(); }
                catch (error) { this.error(String(error.message || error)); }
            });
            const text = node('input'); text.type = 'text'; text.value = entry.labelText; text.dataset.focus = `${entry.id}:label-text`;
            text.addEventListener('change', () => {
                if (!this.enabled) return;
                this.state.updateDistance(entry.id, { labelText: text.value });
            });
            const visible = node('input'); visible.type = 'checkbox'; visible.checked = entry.showLabel; visible.dataset.focus = `${entry.id}:show-label`;
            visible.addEventListener('change', () => {
                if (!this.enabled) return;
                this.state.updateDistance(entry.id, { showLabel: visible.checked });
            });
            const label = node('label'); label.className = 'molui-check'; label.append(visible, node('span', 'Show label'));
            const appearance = node('div'); appearance.className = 'molui-measurement-appearance';
            appearance.append(field('Color', color), field('Label text', text));
            const actions = node('div'); actions.className = 'molui-row molui-measurement-actions';
            actions.append(label, this.action('Remove distance', () => this.state.removeDistance(entry.id), true));
            card.append(summary, start, end, appearance, actions);
            this.distanceList.append(card);
        }
        const picking = this.state.getDistancePicking();
        this.measureToggle.checked = picking !== 'off';
        this.measureToggle.disabled = !this.enabled;
        this.distanceHint.textContent = picking === 'off' ? 'Enable measuring, then select two atoms. Measurements stay at their original coordinates.'
            : picking === 'first' ? 'Select the first atom.' : 'Select the second atom. Select the first atom again to cancel.';
        this.clearDistances.hidden = this.state.distances.size === 0;
        for (const entry of this.state.volumes.values()) {
            const card = node('div'); card.className = 'molui-card';
            card.append(node('strong', entry.name), node('p', `Values: ${entry.min.toPrecision(4)} to ${entry.max.toPrecision(4)}`));
            if (entry.mode === 'isosurface') {
                const level = number(entry.spec.isoval); level.dataset.focus = `${entry.id}:isoval`;
                level.addEventListener('change', () => {
                    if (!this.enabled) return;
                    try { this.state.updateVolume(entry.id, { isoval: Number(level.value) }); this.error(); }
                    catch (error) { this.error(String(error.message || error)); }
                });
                card.append(field('Contour level', level));
            }
            card.append(...this.styleFields(entry.id, entry.spec, patch => this.state.updateVolume(entry.id, patch)),
                this.action('Remove map', () => this.state.removeVolume(entry.id), true));
            this.volumeList.append(card);
        }
    }
}
