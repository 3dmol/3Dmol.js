import type { AtomStyleSpec, GLModel } from "../GLModel";
import type { GLViewer, SurfaceStyleSpec } from "../GLViewer";
import type { Label, LabelSpec } from "../Label";
import type { AtomSelectionSpec, AtomSpec } from "../specs";
import { CC } from "../colors";
import { UI, UIConfig, UIControlOptions } from "./ui";

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
 */
export class StateManager {
    public readonly selections = new Map<string, UISelection>();
    public readonly surfaces = new Map<string, UISurface>();
    public readonly modelTitles = new Map<number, string>();
    public ui: UI;
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
    }

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
        if (!this.disposed) this.ui?.setControls(options);
    }

    public getControls(): Required<UIControlOptions> {
        return this.ui.getControls();
    }

    public updateUI(): void {
        if (!this.disposed) this.ui?.resize();
    }

    public destroy(): void {
        if (this.disposed) return;
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

    /** Bind an entry to its current surface object, since numeric IDs can be reused. */
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
        const models = this.viewer.getInternalState().models as
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
