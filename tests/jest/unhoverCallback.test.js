global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function viewer(current, hits = []) {
    // Exercise the real hover handlers without allocating rendering resources.
    const result = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(result, {
        container: document.createElement("div"), current_hover: current,
        hoverables: [current], targetedObjects: () => hits
    });
    result._viewer = result;
    return result;
}

describe("optional unhover callbacks", () => {
    test.each([undefined, null])("hover exit tolerates callback %s", unhover_callback => {
        const current = { unhover_callback };
        const instance = viewer(current);
        instance.handleHoverContinue(0, 0);
        expect(instance.current_hover).toBeNull();
    });

    test.each([undefined, null])("the next target can be hovered after callback %s", unhover_callback => {
        const current = { unhover_callback };
        const next = { hover_callback: jest.fn() };
        const hits = [{ clickable: next, distance: 1 }];
        const instance = viewer(current, hits);
        const event = new MouseEvent("mousemove");
        instance.handleHoverSelection(0, 0, event);
        expect(instance.current_hover).toBe(next);
        expect(next.hover_callback).toHaveBeenCalledTimes(1);
        expect(next.hover_callback).toHaveBeenCalledWith(next, instance, event, instance.container, hits);
    });

    test.each([
        ["mouse", () => new MouseEvent("mousemove")],
        ["touch", () => new TouchEvent("touchmove")]
    ])("%s callbacks preserve arguments, receiver, and exit-before-enter order", (_name, makeEvent) => {
        const calls = [];
        const current = { unhover_callback(...args) { calls.push(["exit", this, args]); } };
        const next = { hover_callback(...args) { calls.push(["enter", this, args]); } };
        const hits = [{ clickable: next, distance: 1 }];
        const instance = viewer(current, hits);
        const event = makeEvent();
        instance.handleHoverSelection(0, 0, event);
        expect(calls).toEqual([
            ["exit", current, [current, instance, event, instance.container, hits]],
            ["enter", next, [next, instance, event, instance.container, hits]]
        ]);
        expect(instance.current_hover).toBe(next);
    });

    test("hover exit still passes undefined event and intersections", () => {
        const current = { unhover_callback: jest.fn() };
        const instance = viewer(current);
        instance.handleHoverContinue(0, 0);
        expect(current.unhover_callback).toHaveBeenCalledTimes(1);
        expect(current.unhover_callback).toHaveBeenCalledWith(current, instance, undefined, instance.container, undefined);
        expect(instance.current_hover).toBeNull();
    });

    test("string callbacks are converted and invoked", () => {
        const current = {
            unhover_callback: "function(atom, viewer, event, container, intersections) { this.received = [atom, viewer, event, container, intersections]; }"
        };
        const instance = viewer(current);
        instance.handleHoverContinue(0, 0);
        expect(typeof current.unhover_callback).toBe("function");
        expect(current.received).toEqual([current, instance, undefined, instance.container, undefined]);
        expect(instance.current_hover).toBeNull();
    });

    test("remaining on the same target invokes neither callback", () => {
        const current = { hover_callback: jest.fn(), unhover_callback: jest.fn() };
        const instance = viewer(current, [{ clickable: current, distance: 1 }]);
        instance.handleHoverSelection(0, 0, new MouseEvent("mousemove"));
        expect(current.hover_callback).not.toHaveBeenCalled();
        expect(current.unhover_callback).not.toHaveBeenCalled();
        expect(instance.current_hover).toBe(current);
    });
});
