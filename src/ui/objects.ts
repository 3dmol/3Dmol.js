import type { StateManager, UIShape } from './state';
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
    readonly shapes = section('Shapes');
    private readonly volumeList = node('div');
    private readonly shapeList = node('div');
    private enabled = true;
    private generation = 0;
    private abort?: AbortController;

    constructor(private readonly state: StateManager, config: UIConfig,
        private readonly error: (message?: string) => void) {
        this.volumeList.className = 'molui-list molui-volumes';
        this.shapeList.className = 'molui-list molui-shapes';
        this.volumes.append(this.volumeList, this.volumeForm(config));
        this.shapes.append(this.shapeList, this.shapeForm());
    }

    setEnabled(enabled: boolean): void {
        this.enabled = enabled;
        if (!enabled) { this.generation++; this.abort?.abort(); }
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
        const update = () => { fileField.hidden = source.value !== 'file'; urlField.hidden = source.value === 'file'; };
        source.addEventListener('change', update); update();
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

    private shapeForm(): HTMLElement {
        const form = node('form'); form.className = 'molui-create molui-shape-create';
        const kind = choice([['sphere', 'Sphere'], ['arrow', 'Arrow'], ['cylinder', 'Cylinder'], ['line', 'Line'], ['box', 'Box']]);
        const name = node('input'); name.type = 'text'; name.placeholder = 'Annotation';
        const start = [number(0), number(0), number(0)], end = [number(3), number(0), number(0)];
        const radius = number(0.5, 0.001), width = number(2, 0.001), height = number(2, 0.001), depth = number(2, 0.001);
        const color = node('input'); color.type = 'color'; color.value = '#ff8800';
        const endpoints = node('div'), dimensions = node('div');
        endpoints.append(...end.map((input, index) => field(`End ${'XYZ'[index]}`, input)));
        dimensions.append(field('Width', width), field('Height', height), field('Depth', depth));
        const radiusField = field('Radius', radius);
        const update = () => {
            endpoints.hidden = !['arrow', 'cylinder', 'line'].includes(kind.value);
            dimensions.hidden = kind.value !== 'box'; radiusField.hidden = ['box', 'line'].includes(kind.value);
        };
        kind.addEventListener('change', update); update();
        const submit = node('button', 'Add shape'); submit.type = 'submit'; submit.className = 'molui-primary';
        form.append(node('h3', 'New shape'), field('Shape', kind), field('Shape name', name), ...start.map((input, index) => field(`${'XYZ'[index]} (Å)`, input)),
            endpoints, radiusField, dimensions, field('Shape color', color), submit);
        form.addEventListener('submit', event => {
            event.preventDefault(); if (!this.enabled) return;
            try {
                const xyz = (inputs: HTMLInputElement[]) => ({ x: Number(inputs[0].value), y: Number(inputs[1].value), z: Number(inputs[2].value) });
                if ([...start, ...end, radius, width, height, depth].some(input => !Number.isFinite(Number(input.value)))) throw new Error('Shape coordinates and dimensions must be finite.');
                const spec: any = { color: color.value, opacity: 1 };
                if (kind.value === 'sphere') Object.assign(spec, { center: xyz(start), radius: Number(radius.value) });
                else if (kind.value === 'box') Object.assign(spec, { center: xyz(start), dimensions: { w: Number(width.value), h: Number(height.value), d: Number(depth.value) } });
                else Object.assign(spec, { start: xyz(start), end: xyz(end), radius: Number(radius.value) });
                this.state.addShape(kind.value as UIShape['kind'], spec, name.value.trim() || kind.options[kind.selectedIndex].text);
                if (!this.state.viewer.selectedAtoms({}).length) { this.state.viewer.zoomTo(); this.state.viewer.render(); }
                this.error();
            } catch (error) { this.error(String(error.message || error)); }
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
        this.shapeList.replaceChildren(); this.volumeList.replaceChildren();
        for (const entry of this.state.shapes.values()) {
            const card = node('div'); card.className = 'molui-card';
            card.append(node('strong', entry.name), ...this.styleFields(entry.id, entry.spec, patch => this.state.updateShape(entry.id, patch)),
                this.action('Remove shape', () => this.state.removeShape(entry.id), true));
            this.shapeList.append(card);
        }
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
