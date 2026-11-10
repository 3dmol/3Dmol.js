/** @jest-environment node */
// Model and state tests need DOM globals but never use native canvas rendering.
jest.mock("canvas", () => ({}));
const { JSDOM } = require("jsdom");
const dom = new JSDOM("<!doctype html><html><body></body></html>");
global.window = dom.window;
global.document = dom.window.document;
global.URL.createObjectURL = jest.fn();
dom.window.URL.createObjectURL = jest.fn();
Object.defineProperty(global, "$3Dmol", {
    configurable: true,
    get: () => dom.window.$3Dmol,
    set: value => { dom.window.$3Dmol = value; }
});
const $3Dmol = require("../../build/3Dmol.js");
const { StateManager } = require("../../build/3Dmol.ui.js");
afterAll(() => {
    dom.window.close();
    delete global.window;
    delete global.document;
    delete global.$3Dmol;
});

// Exercise real model selection and per-frame styling without a WebGL renderer.
function fixture() {
    const model = new $3Dmol.GLModel(0);
    model.addMolData("2\nfirst\nC 0 0 0\nO 1 0 0\n2\nsecond\nC 0 1 0\nO 1 1 0\n", "xyz", {
        multimodel: true, frames: true
    });
    const surfaceObjects = new Map();
    const externalLabel = { owner: 'host' };
    const labels = new Set([externalLabel]);
    const viewer = {
        getInternalState: () => ({ models: [model.getInternalState()] }),
        getModel: () => model,
        selectedAtoms: spec => model.selectedAtoms(spec),
        setStyle: (spec, style) => model.setStyle(spec, style),
        addStyle: (spec, style) => model.setStyle(spec, style, true),
        render: jest.fn(),
        addResLabels: jest.fn(spec => model.selectedAtoms(spec).map(atom => {
            const label = { atom };
            labels.add(label);
            return label;
        })),
        removeLabels: jest.fn(batch => batch.forEach(label => labels.delete(label))),
        setSurfaceMaterialStyle: jest.fn(),
        getSurface: id => surfaceObjects.get(id),
        removeSurface: jest.fn(id => surfaceObjects.delete(id)),
        addSurface: jest.fn(() => {
            surfaceObjects.set(7, {});
            // Symmetry surfaces resolve an array; the promise ID is authoritative.
            return Object.assign(Promise.resolve([19, 23]), { surfid: 7 });
        })
    };
    return { model, viewer, surfaceObjects, labels, externalLabel };
}

function modelRemovalFixture() {
    const makeModel = id => {
        const model = new $3Dmol.GLModel(id);
        model.addAtoms([{ elem: "C", x: id, y: 0, z: 0, bonds: [] }]);
        return model;
    };
    const other = makeModel(0);
    const target = makeModel(1);
    const viewer = Object.assign(Object.create($3Dmol.GLViewer.prototype), {
        models: [other, target], modelGroup: { remove: jest.fn() }, surfaces: {},
        render: jest.fn(), removeSurface: jest.fn(function (id) { delete this.surfaces[id]; })
    });
    return { viewer, target, other, makeModel };
}

describe("UI appearance state", () => {
    let mount;
    beforeEach(() => {
        mount = jest.spyOn(StateManager.prototype, "initiateUI").mockImplementation(function () {
            this.ui = { destroy: jest.fn(), resize: jest.fn() };
            return this.ui;
        });
    });
    afterEach(() => mount.mockRestore());

    test("overlapping appearances restore previous layers and original styles on every frame", () => {
        const { model, viewer } = fixture();
        model.setStyle({}, { line: { color: "blue" } });
        const state = new StateManager(viewer);
        const base = state.addSelection({});
        state.addStyle({ stick: { radius: 0.2 } }, base);
        const carbon = state.addSelection({ elem: "C" });
        state.addStyle({ sphere: { color: "red" } }, carbon);

        const frames = model.getInternalState().frames;
        expect(frames).toHaveLength(2);
        for (const frame of frames) {
            expect(frame[0].style).toEqual({ sphere: { color: "red" } });
            expect(frame[1].style).toEqual({ stick: { radius: 0.2 } });
        }
        state.toggleHide(carbon);
        for (const frame of frames) expect(frame[0].style).toEqual({ stick: { radius: 0.2 } });
        state.removeSelection(base);
        for (const frame of frames) {
            expect(frame[0].style).toEqual({ line: { color: "blue" } });
            expect(frame[1].style).toEqual({ line: { color: "blue" } });
        }
        state.destroy();
    });

    test("editing a selection restores atoms that no longer match", () => {
        const { model, viewer } = fixture();
        const state = new StateManager(viewer);
        const selection = state.addSelection({ elem: "C" });
        state.addStyle({ sphere: { color: "red" } }, selection);
        state.addSelection({ elem: "O" }, selection);
        for (const frame of model.getInternalState().frames) {
            expect(frame[0].style).toEqual({ line: {} });
            expect(frame[1].style).toEqual({ sphere: { color: "red" } });
        }
        state.destroy();
    });

    test("renaming trims text and keeps atom and representation settings intact", () => {
        const { model, viewer } = fixture();
        const state = new StateManager(viewer);
        const id = state.addSelection({ elem: "C" });
        state.addStyle({ stick: { radius: 0.3 } }, id);
        expect(state.renameSelection(id, "  Binding site  ")).toBe("Binding site");
        expect(state.renameSelection(id, "  ")).toBe("Binding site");
        expect(state.selections.get(id).spec).toEqual({ elem: "C" });
        expect(model.selectedAtoms({ elem: "C" })[0].style).toEqual({ stick: { radius: 0.3 } });
        state.destroy();
    });

    test("one selection visibility setting enables all its representations on every frame", () => {
        const { model, viewer } = fixture();
        const state = new StateManager(viewer);
        const id = state.addSelection({});
        const stick = state.addStyle({ stick: { radius: 0.2 } }, id);
        state.addStyle({ sphere: { scale: 0.3 } }, id);
        state.toggleHideStyle(id, stick);
        state.setSelectionVisible(id, true);
        for (const frame of model.getInternalState().frames) {
            expect(frame[0].style).toEqual({ stick: { radius: 0.2 }, sphere: { scale: 0.3 } });
        }
        state.setSelectionVisible(id, false);
        for (const frame of model.getInternalState().frames) expect(frame[0].style).toEqual({ line: {} });
        state.setSelectionVisible(id, true);
        for (const frame of model.getInternalState().frames) {
            expect(frame[0].style).toEqual({ stick: { radius: 0.2 }, sphere: { scale: 0.3 } });
        }
        state.destroy();
    });

    test("residue-label toggles avoid duplicates and keep labels scoped to their selection", () => {
        const { viewer, labels, externalLabel } = fixture();
        const state = new StateManager(viewer);
        const carbon = state.addSelection({ elem: "C" });
        const oxygen = state.addSelection({ elem: "O" });
        state.setResidueLabels(carbon, true);
        const carbonLabel = viewer.addResLabels.mock.results[0].value[0];
        state.setResidueLabels(carbon, true);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(1);
        state.setResidueLabels(oxygen, true);
        const oxygenLabel = viewer.addResLabels.mock.results[1].value[0];
        expect(labels).toEqual(new Set([externalLabel, carbonLabel, oxygenLabel]));

        state.setSelectionVisible(carbon, false);
        expect(state.selections.get(carbon).labelsEnabled).toBe(true);
        expect(labels).toEqual(new Set([externalLabel, oxygenLabel]));
        state.setSelectionVisible(carbon, true);
        expect(viewer.addResLabels).toHaveBeenCalledTimes(3);
        state.setResidueLabels(carbon, false);
        expect(labels).toEqual(new Set([externalLabel, oxygenLabel]));
        state.removeSelection(oxygen);
        expect(labels).toEqual(new Set([externalLabel]));
        state.setResidueLabels(carbon, true);
        state.destroy();
        expect(labels).toEqual(new Set([externalLabel]));
        expect(viewer.removeLabels.mock.calls.flatMap(([batch]) => batch)).not.toContain(externalLabel);
    });

    test("clearing many residue labels draws once and preserves other selections and host labels", () => {
        const { viewer } = fixture();
        const makeLabel = () => ({ sprite: {}, dispose: jest.fn() });
        const hostLabel = makeLabel();
        const owned = Array.from({ length: 500 }, makeLabel);
        const otherSelectionLabels = [makeLabel(), makeLabel()];
        const viewerLabels = [hostLabel];
        Object.assign(viewer, {
            labels: viewerLabels,
            modelGroup: { remove: jest.fn() },
            show: jest.fn(),
            removeLabels: $3Dmol.GLViewer.prototype.removeLabels,
            render: jest.fn(function () { this.show(); })
        });
        viewer.addResLabels.mockImplementationOnce(() => {
            viewerLabels.push(...owned);
            return owned;
        }).mockImplementationOnce(() => {
            viewerLabels.push(...otherSelectionLabels);
            return otherSelectionLabels;
        });
        const state = new StateManager(viewer);
        const first = state.addSelection({ elem: "C" });
        const second = state.addSelection({ elem: "O" });
        state.setResidueLabels(first, true);
        state.setResidueLabels(second, true);
        viewer.show.mockClear();
        viewer.render.mockClear();

        state.setResidueLabels(first, false);

        expect(viewer.show).toHaveBeenCalledTimes(1);
        expect(viewer.render).not.toHaveBeenCalled();
        expect(viewer.labels).toBe(viewerLabels);
        expect(viewerLabels).toEqual([hostLabel, ...otherSelectionLabels]);
        for (const label of owned) {
            expect(label.dispose).toHaveBeenCalledTimes(1);
            expect(viewer.modelGroup.remove).toHaveBeenCalledWith(label.sprite);
        }
        expect(hostLabel.dispose).not.toHaveBeenCalled();
        otherSelectionLabels.forEach(label => expect(label.dispose).not.toHaveBeenCalled());
        state.destroy();
        expect(viewerLabels).toEqual([hostLabel]);
    });

    test("enabled residue labels follow selection edits and wait until a hidden selection is shown", () => {
        const { viewer, labels, externalLabel } = fixture();
        const state = new StateManager(viewer);
        const id = state.addSelection({ elem: "C" });
        state.setSelectionVisible(id, false);
        state.setResidueLabels(id, true);
        expect(viewer.addResLabels).not.toHaveBeenCalled();
        state.setSelectionVisible(id, true);
        const original = viewer.addResLabels.mock.results[0].value[0];
        state.addSelection({ elem: "O" }, id);
        expect(viewer.removeLabels).toHaveBeenCalledWith([original]);
        expect(viewer.addResLabels).toHaveBeenLastCalledWith({ elem: "O" }, expect.any(Object));
        expect(labels.has(original)).toBe(false);
        state.destroy();
        expect(labels).toEqual(new Set([externalLabel]));
    });

    test("surface color-scheme edits preserve material settings and clear conflicting color modes", async () => {
        const { viewer } = fixture();
        const state = new StateManager(viewer);
        const id = await state.addSurface(1, { opacity: 0.55, wireframe: true, color: "red" }, {});
        state.updateSurface(id, { colorscheme: "greenCarbon" });
        expect(viewer.setSurfaceMaterialStyle).toHaveBeenLastCalledWith(7, {
            opacity: 0.55, wireframe: true, colorscheme: "greenCarbon"
        });
        expect(state.surfaces.get(id).style).toEqual({ opacity: 0.55, wireframe: true, colorscheme: "greenCarbon" });
        state.updateSurface(id, { color: "#000066" });
        expect(viewer.setSurfaceMaterialStyle).toHaveBeenLastCalledWith(7, {
            opacity: 0.55, wireframe: true, color: "#000066"
        });
        expect(viewer.addSurface).toHaveBeenCalledTimes(1);
        state.destroy();
    });

    test("styles applied by an embedding before the first UI edit become the baseline", () => {
        const { model, viewer } = fixture();
        const state = new StateManager(viewer);
        model.setStyle({}, { sphere: { color: "orange" } });
        const selection = state.addSelection({ elem: "C" });
        state.addStyle({ stick: {} }, selection);
        state.removeSelection(selection);
        for (const frame of model.getInternalState().frames) {
            expect(frame[0].style).toEqual({ sphere: { color: "orange" } });
            expect(frame[1].style).toEqual({ sphere: { color: "orange" } });
        }
        state.destroy();
    });

    test("trajectory frames added after a UI appearance restore the original baseline", () => {
        const { model, viewer } = fixture();
        const state = new StateManager(viewer);
        const selection = state.addSelection({});
        state.addStyle({ sphere: { color: "red" } }, selection);
        model.setCoordinates([[[0, 0, 0], [1, 0, 0]], [[0, 2, 0], [1, 2, 0]]], "array");
        state.removeSelection(selection);
        for (const frame of model.getInternalState().frames) {
            expect(frame[0].style).toEqual({ line: {} });
            expect(frame[1].style).toEqual({ line: {} });
        }
        state.destroy();
    });

    test("tracks the actual surface ID for symmetry-expanded surface results", async () => {
        const { viewer } = fixture();
        const state = new StateManager(viewer);
        const id = await state.addSurface(1, { opacity: 0.5 }, {});
        expect(state.surfaces.get(id).surfaceId).toBe(7);
        state.removeSurface(id);
        expect(viewer.removeSurface).toHaveBeenCalledWith(7);
        expect(state.surfaces.size).toBe(0);
        state.destroy();
    });

    test("surface display and calculation selections are independent snapshots", async () => {
        const { viewer } = fixture();
        const state = new StateManager(viewer);
        const shown = { elem: ["C"] };
        const calculation = { elem: ["C", "O"] };
        const pending = state.addSurface(1, { opacity: 0.5 }, shown, "Carbon surface", calculation);
        shown.elem.push("O");
        calculation.elem.splice(0, 1);
        const id = await pending;

        expect(viewer.addSurface).toHaveBeenCalledWith(1, { opacity: 0.5 },
            { elem: ["C"] }, { elem: ["C", "O"] });
        expect(state.surfaces.get(id)).toMatchObject({
            name: "Carbon surface", spec: { elem: ["C"] }, calculationSpec: { elem: ["C", "O"] }
        });
        state.destroy();
    });

    test("surface calculation defaults to shown atoms and preserves explicit All atoms", async () => {
        const { viewer } = fixture();
        const state = new StateManager(viewer);
        const shown = { elem: "C" };
        const defaultId = await state.addSurface(1, {}, shown);
        expect(viewer.addSurface).toHaveBeenLastCalledWith(1, {}, shown, shown);
        expect(state.surfaces.get(defaultId).calculationSpec).toEqual(shown);
        expect(state.surfaces.get(defaultId).calculationSpec).not.toBe(state.surfaces.get(defaultId).spec);

        const allId = await state.addSurface(1, {}, shown, "Carbon in molecule", {});
        expect(viewer.addSurface).toHaveBeenLastCalledWith(1, {}, shown, {});
        expect(state.surfaces.get(allId).calculationSpec).toEqual({});
        state.destroy();
    });

    test.each([
        [{ elem: "N" }, {}, "Select atoms before creating a surface"],
        [{ elem: "C" }, { elem: "N" }, "Select atoms to calculate the surface"],
        [{}, { elem: "C" }, "Calculation atoms must include all shown atoms"],
        [{ elem: "C" }, { elem: "O" }, "Calculation atoms must include all shown atoms"]
    ])("rejects invalid surface selections before starting a calculation: %j, %j", async (shown, calculation, message) => {
        const { viewer } = fixture();
        const state = new StateManager(viewer);
        await expect(state.addSurface(1, {}, shown, "Invalid surface", calculation)).rejects.toThrow(message);
        expect(viewer.addSurface).not.toHaveBeenCalled();
        expect(state.surfaces.size).toBe(0);
        state.destroy();
    });

    test("calculation containment distinguishes atoms with the same index in different models", async () => {
        const { viewer, target, other } = modelRemovalFixture();
        viewer.addSurface = jest.fn();
        const state = new StateManager(viewer);
        expect(target.selectedAtoms({})[0].index).toBe(other.selectedAtoms({})[0].index);
        await expect(state.addSurface(1, {}, { model: other }, "Invalid surface", { model: target }))
            .rejects.toThrow("Calculation atoms must include all shown atoms");
        expect(viewer.addSurface).not.toHaveBeenCalled();
        state.destroy();
    });

    test.each(["loadSurface", "createSurface"])("%s preserves imported calculation selections and explicit overrides", method => {
        const { viewer, surfaceObjects } = fixture();
        const state = new StateManager(viewer);
        const calculation = { elem: ["C", "O"] };
        surfaceObjects.set(7, { allsel: calculation });
        const inferred = state[method](1, { elem: "C" }, {}, 7);
        const overridden = state[method](1, { elem: "C" }, {}, 7, {});
        surfaceObjects.set(8, {});
        const fallback = state[method](1, { elem: "O" }, {}, 8);
        calculation.elem.splice(0, 1);

        expect(state.surfaces.get(inferred).calculationSpec).toEqual({ elem: ["C", "O"] });
        expect(state.surfaces.get(overridden).calculationSpec).toEqual({});
        expect(state.surfaces.get(fallback).calculationSpec).toEqual({ elem: "O" });
        state.destroy();
    });

    test("surfaces finishing after destruction are removed and failures propagate", async () => {
        const { viewer, surfaceObjects } = fixture();
        let resolve;
        viewer.addSurface.mockImplementation(() => {
            surfaceObjects.set(9, {});
            return Object.assign(new Promise(done => { resolve = done; }), { surfid: 9 });
        });
        const state = new StateManager(viewer);
        const pending = state.addSurface(1, {}, {});
        state.destroy();
        resolve(9);
        await expect(pending).rejects.toThrow("closed");
        expect(viewer.removeSurface).toHaveBeenCalledWith(9);
    });

    test("surface computation failures clean up the partial surface and remain visible to the caller", async () => {
        const { viewer, surfaceObjects } = fixture();
        viewer.addSurface.mockImplementation(() => {
            surfaceObjects.set(8, {});
            return Object.assign(Promise.reject(new Error("Surface calculation failed")), { surfid: 8 });
        });
        const state = new StateManager(viewer);
        await expect(state.addSurface(1, {}, {})).rejects.toThrow("Surface calculation failed");
        expect(viewer.removeSurface).toHaveBeenCalledWith(8);
        expect(state.surfaces.size).toBe(0);
        state.destroy();
    });

    test("destroy releases subscriptions and repeated installation disposes the old UI", () => {
        const { viewer } = fixture();
        const first = new StateManager(viewer);
        const listener = jest.fn();
        const unsubscribe = first.subscribe(listener);
        first.refresh();
        expect(listener).toHaveBeenCalledTimes(1);
        unsubscribe();
        first.refresh();
        expect(listener).toHaveBeenCalledTimes(1);
        const second = new StateManager(viewer);
        expect(first.ui.destroy).toHaveBeenCalledTimes(1);
        expect(viewer.ui).toBe(second);
        second.destroy();
        second.destroy();
        expect(second.ui.destroy).toHaveBeenCalledTimes(1);
        expect(viewer.ui).toBeUndefined();
    });

    test("removing a model prunes direct and nested references before its ID can be reused", () => {
        const { viewer, target, other, makeModel } = modelRemovalFixture();
        const state = new StateManager(viewer);
        const retained = [state.addSelection({}), state.addSelection({ model: other })];
        const affected = [
            { model: 1 }, { model: target }, { model: [0, 1] }, { model: -1 },
            { and: [{ model: 1 }] }, { or: [{ model: target }] },
            { not: { model: 1 } }, { within: { distance: 5, sel: { model: target } } },
            { and: [{ or: [{ within: { distance: 5, sel: { not: { model: 1 } } } }] }] }
        ].map(spec => state.addSelection(spec));
        state.addStyle({ sphere: { color: "red" } }, affected[2]);
        state.modelTitles.set(1, "Deleted molecule");
        expect(other.selectedAtoms({})[0].style).toEqual({ sphere: { color: "red" } });

        state.removeModel(target);
        expect(state.getSelectionList()).toEqual(retained);
        expect(state.modelTitles.has(1)).toBe(false);
        expect(viewer.getNextModelId()).toBe(1);
        expect(other.selectedAtoms({})[0].style).toEqual({ line: {} });

        const replacement = makeModel(1);
        viewer.models.push(replacement);
        state.toggleHide(retained[0]);
        expect(replacement.selectedAtoms({})[0].style).toEqual({ line: {} });
        state.destroy();
    });

    test("model removal discards shared/dependent surfaces and preserves other-model surfaces", () => {
        const { viewer, target } = modelRemovalFixture();
        const state = new StateManager(viewer);
        viewer.surfaces = { 41: [], 42: [], 43: [] };
        state.loadSurface(1, {}, {}, 41);
        state.loadSurface(1, { model: 1 }, {}, 42);
        const retained = state.loadSurface(1, { model: 0 }, {}, 43);
        state.removeModel(target);
        expect(viewer.removeSurface).toHaveBeenCalledWith(41);
        expect(viewer.removeSurface).toHaveBeenCalledWith(42);
        expect(viewer.removeSurface).not.toHaveBeenCalledWith(43);
        expect(Array.from(state.surfaces.keys())).toEqual([retained]);
        state.destroy();
    });

    test("model removal discards surfaces that only use the model for calculation", () => {
        const { viewer, target } = modelRemovalFixture();
        const state = new StateManager(viewer);
        viewer.surfaces = { 41: [], 42: [], 43: [] };
        state.loadSurface(1, { model: 0 }, {}, 41, { model: [0, 1] });
        state.loadSurface(1, { model: 0 }, {}, 42, {});
        const retained = state.loadSurface(1, { model: 0 }, {}, 43, { model: 0 });
        state.removeModel(target);
        expect(viewer.removeSurface).toHaveBeenCalledWith(41);
        expect(viewer.removeSurface).toHaveBeenCalledWith(42);
        expect(viewer.removeSurface).not.toHaveBeenCalledWith(43);
        expect(Array.from(state.surfaces.keys())).toEqual([retained]);
        state.destroy();
    });

    test("a pending surface is discarded when a calculation-only model is removed", async () => {
        const { viewer, target, other } = modelRemovalFixture();
        let finish;
        viewer.addSurface = () => {
            viewer.surfaces[44] = [];
            return Object.assign(new Promise(resolve => { finish = resolve; }), { surfid: 44 });
        };
        const state = new StateManager(viewer);
        const pending = state.addSurface(1, {}, { model: other }, "Context surface", { model: [other, target] });
        state.removeModel(target);
        expect(viewer.getModel(0)).toBe(other);
        expect(viewer.getSurface(44)).toBeUndefined();
        finish(44);
        await expect(pending).rejects.toThrow("model used by this surface was removed");
        expect(state.surfaces.size).toBe(0);
        state.destroy();
    });

    test("a surface finishing after its model is removed cannot reappear", async () => {
        const { viewer, target } = modelRemovalFixture();
        let finish;
        viewer.addSurface = () => {
            viewer.surfaces[44] = [];
            return Object.assign(new Promise(resolve => { finish = resolve; }), { surfid: 44 });
        };
        const state = new StateManager(viewer);
        const pending = state.addSurface(1, {}, { model: target });
        state.removeModel(target);
        expect(viewer.getSurface(44)).toBeUndefined();
        finish(44);
        await expect(pending).rejects.toThrow("model used by this surface was removed");
        expect(state.surfaces.size).toBe(0);
        state.destroy();
    });
});

describe("viewer trajectory state", () => {
    function viewer(models = [], shapes = [], labels = []) {
        return Object.assign(Object.create($3Dmol.GLViewer.prototype), {
            models, shapes, labels, viewer_frame: 0
        });
    }

    test("setFrame commits only after model updates finish and propagates failures", async () => {
        let finish;
        const model = { getNumFrames: () => 3, setFrame: jest.fn(() => new Promise(resolve => { finish = resolve; })) };
        const instance = viewer([model]);
        const pending = instance.setFrame(2);
        expect(instance.getFrame()).toBe(0);
        finish();
        await pending;
        expect(instance.getFrame()).toBe(2);
        model.setFrame.mockRejectedValue(new Error("Frame unavailable"));
        await expect(instance.setFrame(1)).rejects.toThrow("Frame unavailable");
        expect(instance.getFrame()).toBe(2);
    });

    test("clamps frame requests and skips removed models", async () => {
        const model = { getNumFrames: () => 3, setFrame: jest.fn(() => Promise.resolve()) };
        const instance = viewer([undefined, model]);
        await instance.setFrame(50);
        expect(instance.getFrame()).toBe(2);
        expect(model.setFrame).toHaveBeenLastCalledWith(2);
        await instance.setFrame(-2);
        expect(instance.getFrame()).toBe(2);
        expect(model.setFrame).toHaveBeenLastCalledWith(2);
    });

    test("counts shape and label frames in sparse scenes, including frame zero", () => {
        expect(viewer([], [undefined, { frame: 0 }]).getNumFrames()).toBe(1);
        expect(viewer([undefined], [undefined, { frame: 2 }], [undefined, { frame: 4 }]).getNumFrames()).toBe(5);
    });

    test("selections referring to removed or missing models are empty", () => {
        const instance = viewer([undefined]);
        expect(instance.selectedAtoms({ model: 99 })).toEqual([]);
        expect(instance.selectedAtoms({ model: [0, 99] })).toEqual([]);
    });

    test("setting an empty model frame preserves its empty atom list", async () => {
        const model = new $3Dmol.GLModel(0);
        await model.setFrame(0);
        expect(model.getInternalState().atoms).toEqual([]);
    });
});


describe("batch label removal", () => {
    const makeLabel = () => ({ sprite: {}, dispose: jest.fn() });
    const makeViewer = labels => Object.assign(Object.create($3Dmol.GLViewer.prototype), {
        labels, modelGroup: { remove: jest.fn() }, show: jest.fn()
    });

    test("removes only requested labels once, preserving order, shared array and caller list", () => {
        const labels = Array.from({ length: 1000 }, makeLabel);
        const original = [...labels];
        const viewer = makeViewer(labels);
        const foreign = makeLabel();
        const requested = [original[990], original[2], original[600], original[2], foreign];
        const snapshot = [...requested];
        // A batch must not call the single-label operation repeatedly.
        viewer.removeLabel = jest.fn(() => { throw new Error("Unexpected single removal"); });
        const splice = jest.spyOn(labels, "splice");
        expect(viewer.removeLabels(requested)).toBe(viewer);
        expect(viewer.labels).toBe(labels);
        expect([...labels]).toEqual(original.filter(label => !snapshot.includes(label)));
        expect(requested).toEqual(snapshot);
        expect(splice).not.toHaveBeenCalled();
        expect(viewer.removeLabel).not.toHaveBeenCalled();
        expect(viewer.show).toHaveBeenCalledTimes(1);
        expect(viewer.modelGroup.remove).toHaveBeenCalledTimes(3);
        for (const label of original) expect(label.dispose).toHaveBeenCalledTimes(snapshot.includes(label) ? 1 : 0);
        expect(foreign.dispose).not.toHaveBeenCalled();
        splice.mockRestore();
    });

    test("empty, unknown and already removed labels do not redraw or dispose twice", () => {
        const label = makeLabel();
        const foreign = makeLabel();
        const viewer = makeViewer([label]);
        viewer.removeLabels([]);
        viewer.removeLabels([foreign]);
        expect(viewer.show).not.toHaveBeenCalled();
        viewer.removeLabels([label, label]);
        viewer.removeLabels([label]);
        expect(viewer.show).toHaveBeenCalledTimes(1);
        expect(label.dispose).toHaveBeenCalledTimes(1);
        expect(foreign.dispose).not.toHaveBeenCalled();
    });

    test("supports the viewer label array itself as input", () => {
        const labels = Array.from({ length: 3 }, makeLabel);
        const original = [...labels];
        const viewer = makeViewer(labels);
        viewer.removeLabels(labels);
        expect(viewer.labels).toBe(labels);
        expect(labels).toEqual([]);
        original.forEach(label => expect(label.dispose).toHaveBeenCalledTimes(1));
        expect(viewer.show).toHaveBeenCalledTimes(1);
    });

    test("single-label removal delegates and continues to display immediately", () => {
        const label = makeLabel();
        const viewer = makeViewer([label]);
        const batch = jest.spyOn(viewer, "removeLabels");
        expect(viewer.removeLabel(label)).toBe(viewer);
        expect(batch).toHaveBeenCalledWith([label]);
        expect(viewer.labels).toEqual([]);
        expect(viewer.modelGroup.remove).toHaveBeenCalledWith(label.sprite);
        expect(label.dispose).toHaveBeenCalledTimes(1);
        expect(viewer.show).toHaveBeenCalledTimes(1);
    });
});
