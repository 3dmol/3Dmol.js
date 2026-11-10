import type { AtomStyleSpec, GLModel } from "../GLModel";
import type { GLViewer, SurfaceStyleSpec } from "../GLViewer";
import type { Label } from "../Label";
import type { AtomSelectionSpec, AtomSpec } from "../specs";
import { UI, UIConfig } from "./ui";

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
}

export interface UISurface {
    id: string;
    name: string;
    surfaceId: number;
    spec: AtomSelectionSpec;
    calculationSpec: AtomSelectionSpec;
    style: SurfaceStyleSpec;
    type: number;
}

const surfaceTypes: Record<string, number> = { VDW: 1, SES: 2, SAS: 3, MS: 4 };

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

/**
 * Owns UI appearance layers without taking ownership of the rest of the scene.
 * Each enabled selection replaces the styles on its atoms, in list order; its
 * own styles are additive. Disabling/removing a layer restores earlier layers
 * or the original atom styles. Coordinates, callbacks and labels are preserved.
 */
export class StateManager {
    public readonly selections = new Map<string, UISelection>();
    public readonly surfaces = new Map<string, UISurface>();
    public readonly modelTitles = new Map<number, string>();
    public ui: UI;
    private readonly listeners = new Set<() => void>();
    private readonly baseline = new WeakMap<AtomSpec, AtomStyleSpec>();
    private readonly inheritedStyles = new WeakMap<AtomStyleSpec, AtomStyleSpec>();
    private readonly selectionLabels = new Map<string, Label[]>();
    private readonly pendingSurfaces = new Map<number, {
        spec: AtomSelectionSpec; calculationSpec: AtomSelectionSpec; removed: boolean
    }>();
    private sequence = 0;
    private disposed = false;

    constructor(public readonly viewer: GLViewer, private readonly config: UIConfig = {}) {
        // Reinitializing an embedded viewer must not leave duplicate panels/timers.
        if (viewer.ui && viewer.ui !== this) viewer.ui.destroy();
        viewer.ui = this;
        this.initiateUI();
    }

    public subscribe(listener: () => void): () => void {
        if (this.disposed) return () => undefined;
        this.listeners.add(listener);
        return () => { this.listeners.delete(listener); };
    }

    public addSelection(spec: AtomSelectionSpec = {}, sid?: string, name?: string): string {
        const id = sid || this.nextID("selection");
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

    /** Toggle residue labels owned by this selection, leaving other labels alone. */
    public setResidueLabels(id: string, enabled: boolean): void {
        const selection = this.selection(id);
        selection.labelsEnabled = enabled;
        this.syncSelectionLabels(selection);
        this.refresh();
    }

    public addStyle(spec: AtomStyleSpec, sid: string, stid?: string): string {
        const selection = this.selection(sid);
        const id = stid || this.nextID("style");
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
        const id = this.nextID("selection");
        const styles = new Map<string, UIStyle>();
        if (Object.keys(style).length) {
            const styleId = this.nextID("style");
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
                if (this.viewer.getSurface(surface.surfaceId)) this.viewer.removeSurface(surface.surfaceId);
                this.surfaces.delete(surfaceId);
            }
        }
        for (const [surfaceId, pending] of this.pendingSurfaces) {
            if (this.surfaceIncludesModel(pending.spec, model, modelCount) ||
                this.surfaceIncludesModel(pending.calculationSpec, model, modelCount)) {
                pending.removed = true;
                if (this.viewer.getSurface(surfaceId)) this.viewer.removeSurface(surfaceId);
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
        const shownAtoms = this.viewer.selectedAtoms(selection);
        if (!shownAtoms.length) throw new Error("Select atoms before creating a surface.");
        const calculationAtoms = new Set(this.viewer.selectedAtoms(calculation));
        if (!calculationAtoms.size) throw new Error("Select atoms to calculate the surface.");
        if (shownAtoms.some(atom => !calculationAtoms.has(atom))) {
            throw new Error("Calculation atoms must include all shown atoms.");
        }
        // The completion result differs for symmetry-expanded surfaces. surfid on
        // the returned promise is the stable ID in both cases.
        const pending = this.viewer.addSurface(surfaceType, surfaceStyle, selection, calculation) as
            Promise<unknown> & { surfid: number };
        const entry = { spec: selection, calculationSpec: calculation, removed: false };
        this.pendingSurfaces.set(pending.surfid, entry);
        try {
            await pending;
        } catch (error) {
            if (this.viewer.getSurface(pending.surfid)) this.viewer.removeSurface(pending.surfid);
            throw error;
        } finally {
            this.pendingSurfaces.delete(pending.surfid);
        }
        if (this.disposed || entry.removed) {
            if (this.viewer.getSurface(pending.surfid)) this.viewer.removeSurface(pending.surfid);
            throw new Error(this.disposed ? "The viewer controls were closed before the surface finished."
                : "A model used by this surface was removed before it finished.");
        }
        const id = this.nextID("surface");
        this.surfaces.set(id, {
            id, name, surfaceId: pending.surfid, spec: selection, calculationSpec: calculation,
            style: surfaceStyle, type: surfaceType
        });
        this.viewer.render();
        this.refresh();
        return id;
    }

    /** Change colors or material settings without recomputing surface geometry. */
    public updateSurface(id: string, patch: SurfaceStyleSpec): void {
        const surface = this.surfaces.get(id);
        if (!surface || !this.viewer.getSurface(surface.surfaceId)) return;
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
        if (this.viewer.getSurface(surface.surfaceId)) this.viewer.removeSurface(surface.surfaceId);
        this.surfaces.delete(id);
        this.viewer.render();
        this.refresh();
    }

    /** Register a surface already created by autoload or an embedding. */
    public loadSurface(type: string | number, spec: AtomSelectionSpec,
        style: SurfaceStyleSpec, surfaceId: number, calculationSpec?: AtomSelectionSpec): string {
        const id = this.nextID("surface");
        this.surfaces.set(id, {
            id, name: "Molecular surface", surfaceId, spec: copy(spec), style: copy(style),
            calculationSpec: copy(calculationSpec ?? this.viewer.getSurface(surfaceId)?.allsel ?? spec),
            type: this.surfaceType(type)
        });
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
        for (const [id, surface] of this.surfaces) {
            if (!this.viewer.getSurface(surface.surfaceId)) this.surfaces.delete(id);
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
            this.selectionLabels.set(selection.id, this.viewer.addResLabels(selection.spec, {
                fontSize: 12, showBackground: false, fontColor: '#000066'
            }));
        }
    }

    private nextID(kind: string): string {
        return `${kind}-${++this.sequence}`;
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
            }
        }
        return atoms;
    }

    private renderLayers(): void {
        if (this.disposed) return;
        const atoms = this.captureStyles();
        for (const atom of atoms) atom.style = copy(this.baseline.get(atom));
        // Invalidate model geometry after restoring style objects, including when
        // the final layer was removed. An empty additive style changes no styles.
        this.viewer.addStyle({}, {});
        for (const selection of this.selections.values()) {
            if (selection.hidden) continue;
            const styles = Array.from(selection.styles.values()).filter(style => !style.hidden);
            if (!styles.length) continue;
            this.viewer.setStyle(selection.spec, {});
            for (const style of styles) this.viewer.addStyle(selection.spec, copy(style.spec));
        }
        for (const atom of atoms) {
            if (atom.style) this.inheritedStyles.set(atom.style, this.baseline.get(atom));
        }
        for (const selection of this.selections.values()) this.syncSelectionLabels(selection);
        this.viewer.render();
        this.refresh();
    }
}
