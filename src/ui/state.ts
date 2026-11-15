import type { AtomStyleSpec, GLModel } from "../GLModel";
import type { GLViewer, SurfaceStyleSpec } from "../GLViewer";
import type { Label, LabelSpec } from "../Label";
import type { AtomSelectionSpec, AtomSpec } from "../specs";
import type { GLShape, IsoSurfaceSpec, ShapeSpec } from "../GLShape";
import type { VolumeData } from "../VolumeData";
import type { GLVolumetricRender } from "../VolumetricRender";
import { CC } from "../colors";
import { UI, UIConfig, UIControlOptions } from "./ui";
import { LearningSession } from './learning';
import type { LearningConfig } from './learning';

export interface UIStyle {
    id: string;
    spec: AtomStyleSpec;
    hidden: boolean;
}

export interface UISelection {
    id: string;
    name: string;
    spec: AtomSelectionSpec;
    styles: Map<string, UIStyle>;
    hidden: boolean;
    labelsEnabled: boolean;
    labelStyle?: LabelSpec;
    /** Coordinate frame used for this selection's currently owned static labels. */
    labelFrame?: number;
    /** Preserve representations from earlier layers for matching atoms. */
    additive?: boolean;
}

export interface UIViewSettings {
    background: string;
    backgroundAlpha: number;
    projection: "perspective" | "orthographic";
}

export interface UISurface {
    id: string;
    name: string;
    surfaceId: number;
    spec: AtomSelectionSpec;
    calculationSpec: AtomSelectionSpec;
    style: SurfaceStyleSpec;
    type: number;
    /** Coordinate frame at which this static surface was calculated. */
    frame?: number;
}

export interface UIShape {
    id: string;
    name: string;
    kind: 'sphere' | 'arrow' | 'cylinder' | 'line' | 'box';
    shape: GLShape;
    spec: ShapeSpec;
}

export interface UIVolume {
    id: string;
    name: string;
    data: VolumeData;
    shape: GLShape | GLVolumetricRender;
    spec: IsoSurfaceSpec;
    mode: 'isosurface' | 'volume';
    min: number;
    max: number;
}

export interface UIDistance {
    id: string;
    name: string;
    start: { x: number; y: number; z: number };
    end: { x: number; y: number; z: number };
    distance: number;
    color: string;
    showLabel: boolean;
    labelText: string;
    shape: GLShape;
    label: Label;
}

// Match the public SurfaceType enum without importing the surface engine into the UI bundle.
const surfaceTypes: Record<string, number> = { VDW: 1, SES: 4, SAS: 3, MS: 2 };

/** Copy configuration objects without losing callbacks, colors or gradient instances. */
function copy<T>(value: T): T {
    if (Array.isArray(value)) return value.map(item => copy(item)) as T;
    if (value && typeof value === "object" && Object.getPrototypeOf(value) === Object.prototype) {
        const result: Record<string, unknown> = {};
        for (const key of Object.keys(value)) result[key] = copy(value[key]);
        return result as T;
    }
    return value;
}

/** Portable UI object settings; core and UI bundles may have distinct Color classes. */
function objectSettings(value: any): any {
    return JSON.parse(JSON.stringify(value, (key, item) => {
        if (['callback', 'hover_callback', 'unhover_callback', '__proto__', 'constructor', 'prototype'].includes(key)
            || typeof item === 'function') return undefined;
        if (item && typeof item.getHex === 'function' && ['r', 'g', 'b'].every(axis => Number.isFinite(item[axis]))) {
            return item.getHex();
        }
        return item;
    }));
}

/** Compare the values copied above, retaining identity for callbacks and class instances. */
function sameStyleValue(left: unknown, right: unknown): boolean {
    if (Object.is(left, right)) return true;
    if (!left || !right || typeof left !== "object" || typeof right !== "object") return false;
    if (Array.isArray(left)) {
        if (!Array.isArray(right) || left.length !== right.length) return false;
    } else if (Array.isArray(right) || Object.getPrototypeOf(left) !== Object.prototype ||
        Object.getPrototypeOf(right) !== Object.prototype) {
        return false;
    }
    const keys = Object.keys(left);
    return keys.length === Object.keys(right).length && keys.every(key =>
        Object.prototype.hasOwnProperty.call(right, key) && sameStyleValue(left[key], right[key]));
}

/**
 * Owns UI appearance layers without taking ownership of the rest of the scene.
 * Each enabled selection replaces the styles on its atoms, in list order; its
 * own styles are additive. Disabling/removing a layer restores earlier layers
 * or their underlying atom styles. Programmatic style changes on unaffected atoms
 * become their latest baseline; UI layers retain priority on atoms they restyle.
 * Coordinates, callbacks and labels are preserved.
 *
 * Available from the optional `3Dmol.ui` bundle, loaded after the core browser
 * bundle. Construction attaches controls and sets `viewer.ui` to this manager,
 * replacing any previous manager. Use {@link StateManager#refresh} after direct
 * model changes and {@link StateManager#destroy} to detach the controls.
 * @see UIConfig
 */
export class StateManager {
    public readonly selections = new Map<string, UISelection>();
    public readonly surfaces = new Map<string, UISurface>();
    public readonly shapes = new Map<string, UIShape>();
    public readonly volumes = new Map<string, UIVolume>();
    public readonly distances = new Map<string, UIDistance>();
    public readonly modelTitles = new Map<number, string>();
    public ui: UI;
    /** Optional classroom connection and controls, owned by this manager. */
    public learning?: LearningSession;
    private readonly listeners = new Set<() => void>();
    private readonly baseline = new WeakMap<AtomSpec, AtomStyleSpec>();
    private readonly styledByUI = new WeakSet<AtomSpec>();
    private readonly inheritedStyles = new WeakMap<AtomStyleSpec, AtomStyleSpec>();
    private readonly selectionLabels = new Map<string, Label[]>();
    // Weak keys avoid retaining removed geometry while an old UI entry still exists.
    private readonly surfaceOwners = new WeakMap<ReturnType<GLViewer['getSurface']>, WeakSet<UISurface>>();
    private readonly pendingSurfaces = new Map<number, {
        spec: AtomSelectionSpec; calculationSpec: AtomSelectionSpec; removed: boolean;
        surface: ReturnType<GLViewer['getSurface']>
    }>();
    private sequence = 0;
    private disposed = false;
    private viewSettings: UIViewSettings;
    private batchDepth = 0;
    private pendingLayerRender = false;
    private pendingRefresh = false;
    private distancePicking = false;
    private readonly distanceAtoms = new Map<AtomSpec, Pick<AtomSpec, 'clickable' | 'callback'>>();
    private distanceStart?: { atom: AtomSpec; label: Label };

    constructor(public readonly viewer: GLViewer, private readonly config: UIConfig = {}) {
        // Reinitializing an embedded viewer must not leave duplicate panels/timers.
        if (viewer.ui && viewer.ui !== this) viewer.ui.destroy();
        viewer.ui = this;
        const viewerConfig = viewer.getConfig?.();
        this.viewSettings = {
            background: this.backgroundHex(config.view?.background ?? viewerConfig?.backgroundColor ?? '#ffffff'),
            backgroundAlpha: config.view?.backgroundAlpha ?? viewerConfig?.backgroundAlpha ?? 1,
            projection: config.view?.projection ?? (viewerConfig?.orthographic ? 'orthographic' : 'perspective')
        };
        if (config.view) {
            viewer.setBackgroundColor(this.viewSettings.background, this.viewSettings.backgroundAlpha);
            viewer.setProjection(this.viewSettings.projection);
        }
        this.initiateUI();
        if (config.learning) this.enableLearning(config.learning);
    }

    /** Enable classroom sharing, replacing any existing connection. Pass false to disable it. */
    public enableLearning(config: boolean | LearningConfig = true): LearningSession | undefined {
        this.learning?.destroy();
        this.learning = undefined;
        if (config && !this.disposed) this.learning = new LearningSession(this, config === true ? {} : config);
        return this.learning;
    }

    /** Subscribe to state refreshes; returns a function that removes the subscription. */
    public subscribe(listener: () => void): () => void {
        if (this.disposed) return () => undefined;
        this.listeners.add(listener);
        return () => { this.listeners.delete(listener); };
    }

    /** Apply related synchronous layer edits before rendering and notifying controls. */
    public batch(callback: () => void): void {
        this.batchDepth++;
        try {
            callback();
        } finally {
            this.batchDepth--;
            if (!this.batchDepth) {
                const render = this.pendingLayerRender;
                const refresh = this.pendingRefresh;
                this.pendingLayerRender = this.pendingRefresh = false;
                if (render) this.renderLayers();
                else if (refresh) this.refresh();
            }
        }
    }

    public addSelection(spec: AtomSelectionSpec = {}, sid?: string, name?: string): string {
        const id = sid || this.nextID("selection", this.selections);
        const existing = this.selections.get(id);
        if (existing) {
            this.clearSelectionLabels(id);
            existing.spec = copy(spec);
            if (name !== undefined) existing.name = name;
        } else {
            this.selections.set(id, {
                id, name: name || `Selection ${this.selections.size + 1}`,
                spec: copy(spec), styles: new Map(), hidden: false, labelsEnabled: false
            });
        }
        this.renderLayers();
        return id;
    }

    public removeSelection(id: string): void {
        this.clearSelectionLabels(id);
        if (this.selections.delete(id)) this.renderLayers();
    }

    /** Rename a selection without changing its atom or representation settings. */
    public renameSelection(id: string, name: string): string {
        const selection = this.selection(id);
        selection.name = name.trim() || selection.name;
        this.refresh();
        return selection.name;
    }

    /** One visibility control governs every representation and label in a selection. */
    public setSelectionVisible(id: string, visible: boolean): void {
        const selection = this.selection(id);
        selection.hidden = !visible;
        if (visible) {
            for (const style of selection.styles.values()) style.hidden = false;
        }
        this.renderLayers();
    }

    public toggleHide(id: string): void {
        this.setSelectionVisible(id, this.selection(id).hidden);
    }

    /** Choose whether this layer replaces or adds to earlier representations. */
    public setSelectionAdditive(id: string, additive: boolean): void {
        this.selection(id).additive = additive;
        this.renderLayers();
    }

    /** Toggle residue labels owned by this selection, leaving other labels alone. */
    public setResidueLabels(id: string, enabled: boolean, style?: LabelSpec): void {
        const selection = this.selection(id);
        if (style !== undefined) {
            this.clearSelectionLabels(id);
            selection.labelStyle = copy(style);
        }
        selection.labelsEnabled = enabled;
        this.syncSelectionLabels(selection);
        this.refresh();
    }

    /** Application-readable appearance settings maintained by the view controls. */
    public getViewSettings(): UIViewSettings {
        return { ...this.viewSettings };
    }

    /** Apply view settings and keep the corresponding UI controls in sync. */
    public setViewSettings(patch: Partial<UIViewSettings>): void {
        const next = { ...this.viewSettings, ...patch };
        if (!Number.isFinite(next.backgroundAlpha) || next.backgroundAlpha < 0 || next.backgroundAlpha > 1) {
            throw new Error('Background opacity must be between zero and one.');
        }
        if (next.projection !== 'perspective' && next.projection !== 'orthographic') {
            throw new Error('Choose perspective or orthographic projection.');
        }
        next.background = this.backgroundHex(next.background);
        if (Object.keys(next).every(key => next[key] === this.viewSettings[key])) return;
        this.viewSettings = next;
        this.viewer.setBackgroundColor(next.background, next.backgroundAlpha);
        this.viewer.setProjection(next.projection);
        this.viewer.render();
        this.config.onViewChange?.(this.getViewSettings());
        this.refresh();
    }

    public addStyle(spec: AtomStyleSpec, sid: string, stid?: string): string {
        const selection = this.selection(sid);
        const id = stid || this.nextID("style", selection.styles);
        const previous = selection.styles.get(id);
        selection.styles.set(id, { id, spec: copy(spec), hidden: previous?.hidden || false });
        this.renderLayers();
        return id;
    }

    public removeStyle(sid: string, stid: string): void {
        if (this.selection(sid).styles.delete(stid)) this.renderLayers();
    }

    public toggleHideStyle(sid: string, stid: string): void {
        const style = this.selection(sid).styles.get(stid);
        if (!style) throw new Error("This appearance no longer exists.");
        style.hidden = !style.hidden;
        this.renderLayers();
    }

    public checkAtoms(spec: AtomSelectionSpec): boolean {
        return this.viewer.selectedAtoms(spec).length > 0;
    }

    public getSelectionList(): string[] {
        return Array.from(this.selections.keys());
    }

    /** Import an appearance from an embedding's data-select/data-style options. */
    public createSelectionAndStyle(spec: AtomSelectionSpec, style: AtomStyleSpec): string {
        // Separate layers preserve order, including overlapping numbered selections.
        const id = this.nextID("selection", this.selections);
        const styles = new Map<string, UIStyle>();
        if (Object.keys(style).length) {
            const styleId = this.nextID("style", styles);
            styles.set(styleId, { id: styleId, spec: copy(style), hidden: false });
        }
        this.selections.set(id, {
            id, name: `Selection ${this.selections.size + 1}`,
            spec: copy(spec), styles, hidden: false, labelsEnabled: false
        });
        this.renderLayers();
        return id;
    }

    public setModelTitle(title: string): void {
        const model = this.viewer.getModel();
        if (model && title) this.modelTitles.set(model.getID(), title);
        this.refresh();
    }

    /** Remove model-dependent UI entries before the viewer can reuse its model ID. */
    public removeModel(model: GLModel): void {
        const id = model.getID();
        if (this.viewer.getModel(id) !== model) return;
        const modelCount = this.viewer.getNextModelId();
        // Unscoped selections can include residues from this model too. Rebuild
        // their owned labels against the surviving scene after removal.
        this.clearSelectionLabels(...this.selectionLabels.keys());
        for (const [selectionId, selection] of this.selections) {
            if (this.referencesModel(selection.spec, model, modelCount)) this.selections.delete(selectionId);
        }
        for (const [surfaceId, surface] of this.surfaces) {
            if (this.surfaceIncludesModel(surface.spec, model, modelCount) ||
                this.surfaceIncludesModel(surface.calculationSpec, model, modelCount)) {
                if (this.ownsSurface(surface)) this.viewer.removeSurface(surface.surfaceId);
                this.surfaces.delete(surfaceId);
            }
        }
        for (const [surfaceId, pending] of this.pendingSurfaces) {
            if (this.surfaceIncludesModel(pending.spec, model, modelCount) ||
                this.surfaceIncludesModel(pending.calculationSpec, model, modelCount)) {
                pending.removed = true;
                if (pending.surface && this.viewer.getSurface(surfaceId) === pending.surface) this.viewer.removeSurface(surfaceId);
            }
        }
        this.modelTitles.delete(id);
        this.viewer.removeModel(model);
        // Removed compound selections may also have styled surviving models.
        this.renderLayers();
    }

    /** Calculate with all calculationSpec atoms, showing only the surface for spec. */
    public async addSurface(type: number, style: SurfaceStyleSpec, spec: AtomSelectionSpec,
        name = "Molecular surface", calculationSpec: AtomSelectionSpec = spec): Promise<string> {
        const surfaceType = this.surfaceType(type);
        const surfaceStyle = copy(style);
        const selection = copy(spec);
        const calculation = copy(calculationSpec);
        const frame = this.viewer.getFrame?.() || 0;
        const shownAtoms = this.viewer.selectedAtoms(selection);
        if (!shownAtoms.length) throw new Error("Select atoms before creating a surface.");
        const calculationAtoms = new Set(this.viewer.selectedAtoms(calculation));
        if (!calculationAtoms.size) throw new Error("Select atoms to calculate the surface.");
        if (shownAtoms.some(atom => !calculationAtoms.has(atom))) {
            throw new Error("Calculation atoms must include all shown atoms.");
        }
        // Capture the ID before completion so removal and ID reuse can be tracked
        // while surface generation is pending.
        const pending = this.viewer.addSurface(surfaceType, surfaceStyle, selection, calculation) as
            Promise<unknown> & { surfid: number };
        const entry = {
            spec: selection, calculationSpec: calculation, removed: false,
            surface: this.viewer.getSurface(pending.surfid)
        };
        const removeOwnedSurface = () => {
            if (entry.surface && this.viewer.getSurface(pending.surfid) === entry.surface) {
                this.viewer.removeSurface(pending.surfid);
            }
        };
        this.pendingSurfaces.set(pending.surfid, entry);
        try {
            await pending;
        } catch (error) {
            removeOwnedSurface();
            throw error;
        } finally {
            if (this.pendingSurfaces.get(pending.surfid) === entry) this.pendingSurfaces.delete(pending.surfid);
        }
        if (this.disposed || entry.removed || !entry.surface || this.viewer.getSurface(pending.surfid) !== entry.surface) {
            removeOwnedSurface();
            throw new Error(this.disposed ? "The viewer controls were closed before the surface finished."
                : entry.removed ? "A model used by this surface was removed before it finished."
                    : "This surface was removed before it finished.");
        }
        const id = this.nextID("surface", this.surfaces);
        this.registerSurface({
            id, name, surfaceId: pending.surfid, spec: selection, calculationSpec: calculation,
            style: surfaceStyle, type: surfaceType, frame
        }, entry.surface);
        this.viewer.render();
        this.refresh();
        return id;
    }

    /** Change colors or material settings without recomputing surface geometry. */
    public updateSurface(id: string, patch: SurfaceStyleSpec): void {
        const surface = this.surfaces.get(id);
        if (!surface) return;
        if (!this.ownsSurface(surface)) {
            this.surfaces.delete(id);
            this.refresh();
            return;
        }
        const style = { ...copy(surface.style), ...copy(patch) };
        if (patch.colorscheme !== undefined || patch.color !== undefined) {
            delete style.voldata; delete style.volscheme; delete style.volformat; delete style.map;
            if (patch.color !== undefined) delete style.colorscheme;
            else delete style.color;
        }
        this.viewer.setSurfaceMaterialStyle(surface.surfaceId, copy(style));
        surface.style = style;
        this.viewer.render();
        this.refresh();
    }

    public removeSurface(id: string): void {
        const surface = this.surfaces.get(id);
        if (!surface) return;
        const owned = this.ownsSurface(surface);
        if (owned) this.viewer.removeSurface(surface.surfaceId);
        this.surfaces.delete(id);
        if (owned) this.viewer.render();
        this.refresh();
    }

    /**
     * Add an editable annotation shape. Coordinates are in molecular angstroms.
     * @param {string} kind Sphere, arrow, cylinder, line, or box, in lowercase.
     * @param {Object} spec Shape settings for the chosen kind.
     * @param {string} [name="Shape"] Annotation name.
     */
    public addShape(kind: UIShape['kind'], spec: any, name = 'Shape'): string {
        const constructors = {
            sphere: value => this.viewer.addSphere(value), arrow: value => this.viewer.addArrow(value),
            cylinder: value => this.viewer.addCylinder(value), line: value => this.viewer.addLine(value),
            box: value => this.viewer.addBox(value)
        };
        if (!constructors[kind]) throw new Error('Choose a supported shape.');
        const shape = constructors[kind](copy(spec));
        const id = this.nextID('shape', this.shapes);
        this.shapes.set(id, { id, name, kind, shape, spec: copy(spec) });
        this.viewer.render(); this.refresh();
        return id;
    }

    public updateShape(id: string, patch: ShapeSpec): void {
        const entry = this.shapes.get(id);
        if (!entry) return;
        entry.spec = { ...entry.spec, ...copy(patch) };
        entry.shape.updateStyle(copy(patch));
        this.viewer.render(); this.refresh();
    }

    public removeShape(id: string): void {
        const entry = this.shapes.get(id);
        if (!entry) return;
        this.viewer.removeShape(entry.shape); this.shapes.delete(id);
        this.viewer.render(); this.refresh();
    }

    /**
     * Add a removable, fixed-coordinate distance annotation in angstroms.
     * @param {XYZ} start First endpoint.
     * @param {XYZ} end Second endpoint.
     * @param {string} [name="Distance"] Atom information displayed in the summary.
     * @param {Object} [options] Optional color, showLabel, and labelText settings.
     */
    public addDistance(start: UIDistance['start'], end: UIDistance['end'], name = 'Distance',
        options: Partial<Pick<UIDistance, 'color' | 'showLabel' | 'labelText'>> = {}): string {
        if (typeof name !== 'string' || ![start, end].every(point => point && ['x', 'y', 'z'].every(axis => Number.isFinite(point[axis])))) {
            throw new Error('Distance endpoints must contain finite coordinates.');
        }
        start = { x: start.x, y: start.y, z: start.z }; end = { x: end.x, y: end.y, z: end.z };
        const distance = Math.hypot(end.x - start.x, end.y - start.y, end.z - start.z);
        if (!Number.isFinite(distance)) throw new Error('Distance is too large.');
        const { color, showLabel, labelText } = this.distanceAppearance(options, distance);
        const shape = this.viewer.addCylinder({ dashed: true, radius: 0.05, dashLength: 0.25, gapLength: 0.15,
            start, end, fromCap: 2, toCap: 2, color });
        const label = this.viewer.addLabel(labelText, {
            position: { x: (start.x + end.x) / 2, y: (start.y + end.y) / 2, z: (start.z + end.z) / 2 },
            fontSize: 12, backgroundColor: color, backgroundOpacity: 0.75
        });
        if (!showLabel) label.hide();
        const id = this.nextID('distance', this.distances);
        this.distances.set(id, { id, name, start, end, distance, color, showLabel, labelText, shape, label });
        this.viewer.render(); this.refresh();
        return id;
    }

    /** Update measurement appearance and label text without rebuilding its line. */
    public updateDistance(id: string, patch: Partial<Pick<UIDistance, 'color' | 'showLabel' | 'labelText'>>): void {
        const entry = this.distances.get(id);
        if (!entry) return;
        if (this.viewer.getShape && this.viewer.getShape(entry.shape.shapePosition) !== entry.shape) {
            this.viewer.removeLabels([entry.label]); this.distances.delete(id); this.refresh();
            return;
        }
        const { color, showLabel, labelText } = this.distanceAppearance({
            color: entry.color, showLabel: entry.showLabel, labelText: entry.labelText, ...patch
        }, entry.distance);
        if (color === entry.color && showLabel === entry.showLabel && labelText === entry.labelText) return;
        if (color !== entry.color) {
            entry.shape.updateStyle({ color });
            this.viewer.setLabelStyle(entry.label, { ...entry.label.stylespec, backgroundColor: color });
        }
        if (labelText !== entry.labelText) this.viewer.setLabelText(entry.label, labelText);
        entry.color = color; entry.showLabel = showLabel; entry.labelText = labelText;
        if (showLabel) entry.label.show(); else entry.label.hide();
        this.viewer.render(); this.refresh();
    }

    private distanceAppearance(options: Partial<Pick<UIDistance, 'color' | 'showLabel' | 'labelText'>>,
        distance: number): Pick<UIDistance, 'color' | 'showLabel' | 'labelText'> {
        if (options.color !== undefined && typeof options.color !== 'string') throw new Error('Measurement color must be a string.');
        if (options.showLabel !== undefined && typeof options.showLabel !== 'boolean') throw new Error('Measurement label visibility must be true or false.');
        if (options.labelText !== undefined && typeof options.labelText !== 'string') throw new Error('Measurement label text must be a string.');
        return { color: this.backgroundHex(options.color ?? 'blue'), showLabel: options.showLabel ?? true,
            labelText: options.labelText ?? `${distance.toFixed(3)} Å` };
    }

    public removeDistance(id: string): void {
        const entry = this.distances.get(id);
        if (!entry) return;
        if (!this.viewer.getShape || this.viewer.getShape(entry.shape.shapePosition) === entry.shape) this.viewer.removeShape(entry.shape);
        this.viewer.removeLabels([entry.label]); this.distances.delete(id);
        this.viewer.render(); this.refresh();
    }

    public clearDistances(): void {
        this.batch(() => { for (const id of this.distances.keys()) this.removeDistance(id); });
    }

    /** Atom distance picking is opt-in and never restored from a shared URL. */
    public setDistancePicking(enabled: boolean): void {
        if (enabled && (this.disposed || this.ui?.getControls().edit === false)) return;
        if (this.distancePicking === enabled) return;
        this.distancePicking = enabled;
        if (!enabled) {
            this.clearDistanceStart();
            for (const [atom, previous] of this.distanceAtoms) this.restoreDistanceAtom(atom, previous);
            this.distanceAtoms.clear();
        }
        this.refresh(); this.viewer.render();
    }

    public getDistancePicking(): 'off' | 'first' | 'second' {
        return !this.distancePicking ? 'off' : this.distanceStart ? 'second' : 'first';
    }

    private clearDistanceStart(): void {
        if (this.distanceStart) this.viewer.removeLabels([this.distanceStart.label]);
        this.distanceStart = undefined;
    }

    private restoreDistanceAtom(atom: AtomSpec, previous: Pick<AtomSpec, 'clickable' | 'callback'>): void {
        if (atom.callback !== this.pickDistanceAtom) return;
        atom.callback = previous.callback;
        // Preserve an explicit disable made by the host or classroom controls.
        if (atom.clickable) atom.clickable = previous.clickable;
    }

    private readonly pickDistanceAtom = (atom: AtomSpec): void => {
        if (!this.distancePicking || this.disposed || this.ui?.getControls().edit === false) return;
        const title = (value: AtomSpec) => [value.chain, `${value.resn || ''}${value.resi ?? ''}`, value.atom || value.elem].filter(Boolean).join(':');
        const first = this.distanceStart?.atom;
        this.clearDistanceStart();
        if (first && first !== atom) {
            this.addDistance({ x: first.x, y: first.y, z: first.z }, { x: atom.x, y: atom.y, z: atom.z }, `${title(first)} – ${title(atom)}`);
        } else if (!first) {
            this.distanceStart = { atom, label: this.viewer.addLabel(title(atom), {
                _sessionTransient: true, position: { x: atom.x, y: atom.y, z: atom.z },
                fontSize: 12, backgroundColor: 'blue', backgroundOpacity: 0.75
            } as LabelSpec) };
        }
        this.refresh(); this.viewer.render();
    };

    private refreshDistancePicking(): void {
        if (!this.distancePicking) return;
        const atoms = new Set(this.viewer.selectedAtoms({}));
        if (this.distanceStart && !atoms.has(this.distanceStart.atom)) this.clearDistanceStart();
        for (const [atom, previous] of this.distanceAtoms) if (!atoms.has(atom)) {
            this.restoreDistanceAtom(atom, previous); this.distanceAtoms.delete(atom);
        }
        const added = new Set<AtomSpec>();
        for (const atom of atoms) if (!this.distanceAtoms.has(atom)) {
            this.distanceAtoms.set(atom, { clickable: atom.clickable, callback: atom.callback });
            added.add(atom);
        }
        if (added.size) this.viewer.setClickable({ predicate: atom => added.has(atom) }, true, this.pickDistanceAtom);
    }

    /**
     * Parse a scalar grid once, retaining it for subsequent contour edits.
     * @param {string|ArrayBuffer} data Volumetric file contents.
     * @param {string} format Volumetric file format.
     * @param {IsoSurfaceSpec} [spec] Display settings.
     * @param {string} [name="Volume"] Map name.
     * @param {string} [mode="isosurface"] Either isosurface or volume.
     */
    public addVolume(data: string | ArrayBuffer, format: string, spec: IsoSurfaceSpec = {},
        name = 'Volume', mode: UIVolume['mode'] = 'isosurface'): string {
        const Constructor = (window.$3Dmol as unknown as { VolumeData: typeof VolumeData })?.VolumeData;
        if (!Constructor) throw new Error('The core 3Dmol bundle must be loaded before volume controls.');
        const volume = new Constructor(data, format);
        return this.addVolumeData(volume, spec, name, mode);
    }

    private addVolumeData(volume: VolumeData, spec: IsoSurfaceSpec, name: string, mode: UIVolume['mode']): string {
        if (!['isosurface', 'volume'].includes(mode)) throw new Error('Choose an isosurface or volume display.');
        const expected = volume.size.x * volume.size.y * volume.size.z;
        if (!expected || volume.data.length !== expected) throw new Error('No complete scalar grid was found. Check the map format.');
        let min = Infinity, max = -Infinity;
        for (const value of volume.data) {
            if (!Number.isFinite(value)) throw new Error('The map contains non-finite values.');
            min = Math.min(min, value); max = Math.max(max, value);
        }
        if (min === max) throw new Error('The map must contain a nonconstant scalar field.');
        const id = this.nextID('volume', this.volumes);
        const style = { color: '#3366cc', opacity: 0.8, ...copy(spec), isoval: spec.isoval ?? (min + max) / 2 };
        const entry: UIVolume = { id, name, data: volume, spec: style, mode, min, max, shape: undefined };
        entry.shape = this.volumeShape(entry);
        this.volumes.set(id, entry);
        this.viewer.render(); this.refresh();
        return id;
    }

    /** Portable source settings for editable objects in a viewer URL. */
    public exportObjects() {
        return {
            shapes: Array.from(this.shapes.values(), entry => ({ kind: entry.kind, name: entry.name, spec: objectSettings(entry.spec) })),
            ...(this.distances.size ? { distances: Array.from(this.distances.values(), entry => ({ name: entry.name,
                start: { ...entry.start }, end: { ...entry.end }, color: entry.color, showLabel: entry.showLabel,
                labelText: entry.labelText })) } : {}),
            volumes: Array.from(this.volumes.values(), entry => ({ name: entry.name, mode: entry.mode, spec: objectSettings(entry.spec),
                data: { values: Array.from(entry.data.data), size: { ...entry.data.size }, origin: { ...entry.data.origin },
                    unit: { ...entry.data.unit }, matrix: entry.data.matrix ? Array.from(entry.data.matrix.elements) : null } }))
        };
    }

    /**
     * Restore data-only object settings without evaluating callback source.
     * @param {Object} records A snapshot returned by {@link StateManager#exportObjects}.
     */
    public restoreObjects(records: ReturnType<StateManager['exportObjects']>): void {
        if (!records || !Array.isArray(records.shapes) || !Array.isArray(records.volumes)) throw new Error('Invalid scene objects.');
        if (records.distances !== undefined && !Array.isArray(records.distances)) throw new Error('Invalid distance annotations.');
        for (const record of records.shapes) this.addShape(record.kind, objectSettings(record.spec), record.name);
        for (const record of records.distances || []) this.addDistance(record.start, record.end, record.name, record);
        const core = window.$3Dmol as unknown as { VolumeData: typeof VolumeData; Matrix4: typeof import('../WebGL/math').Matrix4 };
        for (const record of records.volumes) {
            const data = record.data;
            if (!data || !Array.isArray(data.values) || !data.values.every(Number.isFinite) || !['x', 'y', 'z'].every(axis =>
                Number.isInteger(data.size?.[axis]) && data.size[axis] > 0 && Number.isFinite(data.origin?.[axis]) && Number.isFinite(data.unit?.[axis]))) {
                throw new Error('Invalid map grid dimensions.');
            }
            const volume: VolumeData = Object.create(core.VolumeData.prototype);
            Object.assign(volume, { data: new Float32Array(data.values), size: { ...data.size }, origin: { ...data.origin },
                unit: { ...data.unit }, matrix: null, inversematrix: null });
            if (data.matrix) {
                if (data.matrix.length !== 16 || !data.matrix.every(Number.isFinite)) throw new Error('Invalid map coordinate transform.');
                volume.matrix = new core.Matrix4(); volume.matrix.elements.set(data.matrix);
            }
            this.addVolumeData(volume, objectSettings(record.spec), record.name, record.mode);
        }
    }

    private volumeShape(entry: UIVolume): GLShape | GLVolumetricRender {
        if (entry.mode === 'volume') {
            if (!this.viewer.hasVolumetricRender()) throw new Error('Volume rendering requires WebGL 2. Use an isosurface on this device.');
            const color = Array.isArray(entry.spec.color) ? entry.spec.color[0] : entry.spec.color;
            return this.viewer.addVolumetricRender(entry.data, { hidden: entry.spec.hidden,
                transferfn: [{ value: entry.min, color, opacity: 0 },
                    { value: entry.max, color, opacity: entry.spec.opacity }] });
        }
        return this.viewer.addIsosurface(entry.data, copy(entry.spec));
    }

    public updateVolume(id: string, patch: IsoSurfaceSpec): void {
        const entry = this.volumes.get(id);
        if (!entry) return;
        if (patch.isoval !== undefined && !Number.isFinite(patch.isoval)) throw new Error('Contour level must be finite.');
        const spec = { ...entry.spec, ...copy(patch) };
        if (entry.mode === 'isosurface' && (patch.isoval === undefined || patch.isoval === entry.spec.isoval)) {
            (entry.shape as GLShape).updateStyle(copy(patch));
        } else if (entry.mode === 'volume' && Object.keys(patch).every(key => key === 'hidden')) {
            entry.shape.hidden = !!spec.hidden;
        } else {
            const next = this.volumeShape({ ...entry, spec });
            this.viewer.removeShape(entry.shape);
            entry.shape = next;
        }
        entry.spec = spec;
        this.viewer.render(); this.refresh();
    }

    public removeVolume(id: string): void {
        const entry = this.volumes.get(id);
        if (!entry) return;
        this.viewer.removeShape(entry.shape); this.volumes.delete(id);
        this.viewer.render(); this.refresh();
    }

    /** Register a surface already created by autoload or an embedding. */
    public loadSurface(type: string | number, spec: AtomSelectionSpec,
        style: SurfaceStyleSpec, surfaceId: number, calculationSpec?: AtomSelectionSpec): string {
        const id = this.nextID("surface", this.surfaces);
        const surface = this.viewer.getSurface(surfaceId);
        this.registerSurface({
            id, name: "Molecular surface", surfaceId, spec: copy(spec), style: copy(style),
            calculationSpec: copy(calculationSpec ?? surface?.allsel ?? spec),
            type: this.surfaceType(type), frame: this.viewer.getFrame?.() || 0
        }, surface);
        this.refresh();
        return id;
    }

    public createSurface(type: string | number, spec: AtomSelectionSpec,
        style: SurfaceStyleSpec, surfaceId: number, calculationSpec?: AtomSelectionSpec): string {
        return this.loadSurface(type, spec, style, surfaceId, calculationSpec);
    }

    /** Notify controls after models or coordinates are changed by application code. */
    public refresh(): void {
        if (this.disposed) return;
        if (this.batchDepth) {
            this.pendingRefresh = true;
            return;
        }
        for (const [id, surface] of this.surfaces) {
            if (!this.ownsSurface(surface)) this.surfaces.delete(id);
        }
        for (const entries of [this.shapes, this.volumes]) {
            for (const [id, entry] of entries) {
                if (this.viewer.getShape && this.viewer.getShape(entry.shape.shapePosition) !== entry.shape) entries.delete(id);
            }
        }
        for (const [id, entry] of this.distances) {
            if (this.viewer.getShape && this.viewer.getShape(entry.shape.shapePosition) !== entry.shape) {
                this.viewer.removeLabels([entry.label]); this.distances.delete(id);
            }
        }
        this.refreshDistancePicking();
        for (const listener of this.listeners) listener();
    }

    public initiateUI(): UI {
        if (this.disposed) throw new Error("These viewer controls have been destroyed.");
        if (!this.ui) this.ui = new UI(this, this.config);
        else this.refresh();
        return this.ui;
    }

    public showUI(): UI {
        return this.initiateUI();
    }

    /** Enable the edit dialog and trajectory controls independently. */
    public setControls(options: UIControlOptions): void {
        if (this.disposed) return;
        if (options.edit === false) {
            for (const [id, pending] of this.pendingSurfaces) {
                pending.removed = true;
                if (pending.surface && this.viewer.getSurface(id) === pending.surface) this.viewer.removeSurface(id);
            }
        }
        this.ui?.setControls(options);
    }

    public getControls(): Required<UIControlOptions> {
        return this.ui.getControls();
    }

    public updateUI(): void {
        if (!this.disposed) this.ui?.resize();
    }

    /**
     * Detach controls, remove labels owned by UI selections, and release subscribers.
     * Leaves the viewer, models, styles, and completed surfaces in place. Pending
     * UI-created surfaces are discarded when their calculation finishes.
     * Safe to call repeatedly. Use GLViewer.destroy() to dispose the entire viewer.
     */
    public destroy(): void {
        if (this.disposed) return;
        this.learning?.destroy(); this.learning = undefined;
        this.setDistancePicking(false);
        this.disposed = true;
        this.clearSelectionLabels(...this.selectionLabels.keys());
        this.ui?.destroy();
        this.listeners.clear();
        if (this.viewer.ui === this) this.viewer.ui = undefined;
    }

    private clearSelectionLabels(...ids: string[]): void {
        const labels: Label[] = [];
        for (const id of ids) {
            for (const label of this.selectionLabels.get(id) || []) labels.push(label);
            this.selectionLabels.delete(id);
        }
        if (labels.length) this.viewer.removeLabels(labels);
    }

    private syncSelectionLabels(selection: UISelection): void {
        if (!selection.labelsEnabled || selection.hidden) {
            this.clearSelectionLabels(selection.id);
        } else if (!this.selectionLabels.has(selection.id)) {
            const frame = this.viewer.getFrame?.() ?? 0;
            this.selectionLabels.set(selection.id, this.viewer.addResLabels(selection.spec, copy(selection.labelStyle ?? {
                fontSize: 12, showBackground: false, fontColor: '#000066'
            })));
            selection.labelFrame = frame;
        }
    }

    private backgroundHex(color: string): string {
        return `#${CC.color(color).getHex().toString(16).padStart(6, '0')}`;
    }

    private nextID(kind: string, entries: ReadonlyMap<string, unknown>): string {
        let id: string;
        do {
            id = `${kind}-${++this.sequence}`;
        } while (entries.has(id));
        return id;
    }

    private selection(id: string): UISelection {
        const selection = this.selections.get(id);
        if (!selection) throw new Error("This selection no longer exists.");
        return selection;
    }

    private referencesModel(spec: AtomSelectionSpec, model: GLModel, modelCount: number,
        visited = new Set<AtomSelectionSpec>()): boolean {
        if (!spec || visited.has(spec)) return false;
        visited.add(spec);
        const references = Array.isArray(spec.model) ? spec.model : [spec.model];
        if (references.some(reference => reference === model ||
            (typeof reference === "number" && (reference < 0 ? modelCount + reference : reference) === model.getID()))) {
            return true;
        }
        const nested = [...(spec.and || []), ...(spec.or || []), spec.not, spec.within?.sel];
        return nested.some(selection => this.referencesModel(selection, model, modelCount, visited));
    }

    private surfaceIncludesModel(spec: AtomSelectionSpec, model: GLModel, modelCount: number): boolean {
        // Surfaces without an explicit model scope can contain atoms from every
        // model. Conservatively remove them instead of leaving stale geometry.
        return spec.model === undefined || this.referencesModel(spec, model, modelCount);
    }

    /** Bind an entry to its current surface object, since numeric IDs can be reused.
     * @param {UISurface} entry - UI entry that owns the surface.
     * @param {Surface} surface - Surface object returned by GLViewer.getSurface.
     */
    private registerSurface(entry: UISurface, surface: ReturnType<GLViewer['getSurface']>): void {
        if (surface) {
            let owners = this.surfaceOwners.get(surface);
            if (!owners) this.surfaceOwners.set(surface, owners = new WeakSet<UISurface>());
            owners.add(entry);
        }
        this.surfaces.set(entry.id, entry);
    }

    private ownsSurface(entry: UISurface): boolean {
        const surface = this.viewer.getSurface(entry.surfaceId);
        return !!surface && !!this.surfaceOwners.get(surface)?.has(entry);
    }

    private surfaceType(type: string | number): number {
        const result = typeof type === "string" ? surfaceTypes[type.toUpperCase()] : type;
        if (![1, 2, 3, 4].includes(result)) throw new Error("Choose a supported surface type.");
        return result;
    }

    private atoms(): Set<AtomSpec> {
        const atoms = new Set<AtomSpec>();
        // Appearance layers need live atoms; scene snapshots are detached data.
        const models = (this.viewer.getModels ? this.viewer.getModels().map(model => model.getInternalState())
            : this.viewer.getInternalState().models) as
            { atoms: AtomSpec[]; frames: AtomSpec[][] }[];
        for (const model of models) {
            if (!model) continue;
            for (const atom of model.atoms) if (atom) atoms.add(atom);
            // All locally stored frames must be restored, not only the current
            // frame. Remote trajectories reuse the model's current atoms.
            for (const frame of model.frames) {
                if (Array.isArray(frame)) for (const atom of frame) if (atom) atoms.add(atom);
            }
        }
        return atoms;
    }

    private captureStyles(): Set<AtomSpec> {
        const atoms = this.atoms();
        for (const atom of atoms) {
            if (!this.baseline.has(atom)) {
                // setCoordinates() creates new frame atoms by shallow-copying the
                // current atoms. Their shared style objects already contain UI
                // appearances; inherit the original baseline in that case.
                const inherited = atom.style && this.inheritedStyles.get(atom.style);
                this.baseline.set(atom, copy(inherited || atom.style || {}));
            } else if (!this.styledByUI.has(atom)) {
                // An atom outside the UI's current appearance layers can be
                // restyled by its embedding between any two UI edits.
                this.baseline.set(atom, copy(atom.style || {}));
            }
        }
        return atoms;
    }

    private renderLayers(): void {
        if (this.disposed) return;
        if (this.batchDepth) {
            this.pendingLayerRender = true;
            return;
        }
        const atoms = this.captureStyles();
        for (const atom of atoms) atom.style = copy(this.baseline.get(atom));
        // Invalidate model geometry after restoring style objects, including when
        // the final layer was removed. An empty additive style changes no styles.
        this.viewer.addStyle({}, {});
        for (const selection of this.selections.values()) {
            if (selection.hidden) continue;
            const styles = Array.from(selection.styles.values()).filter(style => !style.hidden);
            if (!styles.length) continue;
            if (!selection.additive) this.viewer.setStyle(selection.spec, {});
            for (const style of styles) this.viewer.addStyle(selection.spec, copy(style.spec));
        }
        for (const atom of atoms) {
            const baseline = this.baseline.get(atom);
            if (!sameStyleValue(atom.style, baseline)) {
                this.styledByUI.add(atom);
                if (atom.style) this.inheritedStyles.set(atom.style, baseline);
            } else {
                this.styledByUI.delete(atom);
                if (atom.style) this.inheritedStyles.delete(atom.style);
            }
        }
        for (const selection of this.selections.values()) this.syncSelectionLabels(selection);
        this.viewer.render();
        this.refresh();
    }
}
