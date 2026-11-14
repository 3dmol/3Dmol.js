global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

let viewers, nodes;
const settle = async () => {
    // Deliver real MutationObserver callbacks before allowing the deferred cleanup to run.
    await Promise.resolve();
    await new Promise(resolve => setTimeout(resolve, 10));
};

function element(parent = document.body) {
    const node = document.createElement("div");
    if (parent) parent.appendChild(node);
    nodes.push(node);
    return node;
}

function viewer(container = element(), config = {}) {
    const v = new $3Dmol.GLViewer(container, config);
    viewers.push(v);
    return v;
}

beforeEach(() => {
    jest.useRealTimers();
    viewers = []; nodes = [];
    jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        const canvas = this.config.canvas || document.createElement("canvas");
        this.renderer = {
            domElement: canvas, getCanvas: () => canvas, getAspect: () => 1,
            getXRatio: () => 1, getYRatio: () => 1, isLost: () => false,
            setSize: jest.fn(), setClearColorHex: jest.fn(), setViewport: jest.fn(),
            supportedExtensions: () => ({}), render: jest.fn(), dispose: jest.fn()
        };
    });
});

afterEach(async () => {
    for (const v of viewers) v.destroy();
    for (const node of nodes) node.remove();
    await Promise.resolve();
    jest.useRealTimers();
    jest.restoreAllMocks();
});

describe("automatic viewer cleanup after DOM removal", () => {
    test("removing the container destroys its viewer and graphics resources once", async () => {
        const v = viewer(), renderer = v.renderer, container = v.container, canvas = v.getCanvas();
        container.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
        expect(renderer.dispose).toHaveBeenCalledTimes(1);
        expect(canvas.parentNode).toBeNull();
        document.body.appendChild(container);
        await settle();
        expect(v.isDestroyed()).toBe(true);
        expect(renderer.dispose).toHaveBeenCalledTimes(1);
    });

    test("removing an ancestor also destroys the contained viewer", async () => {
        const ancestor = element(), container = element(element(ancestor)), v = viewer(container);
        ancestor.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("same-task reparenting survives and monitoring follows the new ancestors", async () => {
        const original = element(), destination = element(), container = element(original), v = viewer(container);
        destination.appendChild(container);
        await settle();
        expect(v.isDestroyed()).toBe(false);
        original.remove();
        await settle();
        expect(v.isDestroyed()).toBe(false);
        destination.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("reinsertion in a queued microtask cancels deferred destruction", async () => {
        const v = viewer(), container = v.container;
        container.remove();
        await Promise.resolve();
        document.body.appendChild(container);
        await settle();
        expect(v.isDestroyed()).toBe(false);
        container.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("initially detached containers remain usable until attached and subsequently removed", async () => {
        const container = element(null), v = viewer(container);
        await settle();
        expect(v.isDestroyed()).toBe(false);
        document.body.appendChild(container);
        await settle();
        expect(v.isDestroyed()).toBe(false);
        container.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("attachment and removal in one mutation batch still triggers cleanup", async () => {
        const container = element(null), v = viewer(container);
        document.body.appendChild(container);
        container.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("remembers attachment when the child and then its former parent are removed in one task", async () => {
        const parent = element(), container = element(null), v = viewer(container);
        parent.appendChild(container);
        container.remove();
        parent.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("staging inside an already removed subtree does not count as attachment", async () => {
        const parent = element(), container = element(null), v = viewer(container);
        parent.remove();
        parent.appendChild(container);
        container.remove();
        await settle();
        expect(v.isDestroyed()).toBe(false);
    });

    test("resize starts monitoring an initially detached viewer mounted into an existing shadow root", async () => {
        const root = element().attachShadow({ mode: "open" });
        const container = element(null), v = viewer(container);
        await settle();
        root.appendChild(container);
        v.resize();
        container.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("setContainer ignores the old element and monitors the replacement", async () => {
        const v = viewer(), old = v.container, replacement = element();
        v.setContainer(replacement);
        old.remove();
        await settle();
        expect(v.isDestroyed()).toBe(false);
        replacement.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("setContainer can deliberately stage a viewer in a detached element", async () => {
        const v = viewer(), old = v.container, replacement = element(null);
        v.setContainer(replacement);
        old.remove();
        await settle();
        expect(v.isDestroyed()).toBe(false);
        document.body.appendChild(replacement);
        await settle();
        replacement.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("setContainer cancels cleanup already queued for the previous container", async () => {
        const v = viewer(), replacement = element();
        v.container.remove();
        await Promise.resolve();
        v.setContainer(replacement);
        await settle();
        expect(v.isDestroyed()).toBe(false);
        replacement.remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("hiding an attached container does not destroy the viewer", async () => {
        const v = viewer();
        v.container.style.display = "none";
        v.container.hidden = true;
        await settle();
        expect(v.isDestroyed()).toBe(false);
    });

    test("autoDestroy false retains a detached viewer for explicit cleanup", async () => {
        const v = viewer(element(), { autoDestroy: false }), renderer = v.renderer;
        v.container.remove();
        await settle();
        expect(v.isDestroyed()).toBe(false);
        expect(renderer.dispose).not.toHaveBeenCalled();
        v.destroy();
        expect(renderer.dispose).toHaveBeenCalledTimes(1);
    });

    test.each([
        [false, "container"], [false, "host"],
        [true, "container"], [true, "host"]
    ])("detects %s nested shadow DOM removal of the %s", async (nested, removal) => {
        const host = element();
        let root = host.attachShadow({ mode: "open" });
        if (nested) root = element(root).attachShadow({ mode: "open" });
        const container = element(root), v = viewer(container);
        (removal === "host" ? host : container).remove();
        await settle();
        expect(v.isDestroyed()).toBe(true);
    });

    test("removing a grid container destroys every viewer and releases the shared canvas", async () => {
        const container = element(), grid = $3Dmol.createViewerGrid(container, { rows: 1, cols: 2 });
        const [a, b] = grid[0]; viewers.push(a, b);
        const canvas = a.getCanvas(), renderers = [a.renderer, b.renderer];
        container.remove();
        await settle();
        expect(a.isDestroyed()).toBe(true);
        expect(b.isDestroyed()).toBe(true);
        expect(grid[0]).toEqual([null, null]);
        expect(canvas.parentNode).toBeNull();
        for (const renderer of renderers) expect(renderer.dispose).toHaveBeenCalledTimes(1);
    });

    test("manual destruction cancels queued automatic cleanup", async () => {
        jest.useFakeTimers();
        const v = viewer(), renderer = v.renderer;
        v.container.remove();
        await Promise.resolve();
        await Promise.resolve();
        expect(jest.getTimerCount()).toBeGreaterThan(0);
        v.destroy();
        v.destroy();
        expect(jest.getTimerCount()).toBe(0);
        await jest.advanceTimersByTimeAsync(10);
        expect(renderer.dispose).toHaveBeenCalledTimes(1);
    });
});
