global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const lists = [
    ["clickable", "clickables"],
    ["hoverable", "hoverables"],
    ["contextMenuEnabled", "contextMenuEnabledObjects"]
];
function model(atoms, frames = [atoms]) {
    const m = new $3Dmol.GLModel(0);
    m.setInternalState({ atoms, frames });
    return m;
}
function viewer(models, shapes = []) {
    return Object.assign(Object.create($3Dmol.GLViewer.prototype), {
        models, shapes, clickables: [], hoverables: [], contextMenuEnabledObjects: []
    });
}
const enabled = name => ({ name, clickable: true, hoverable: true, contextMenuEnabled: true });

afterEach(() => jest.restoreAllMocks());

describe("viewer interaction lists", () => {
    test.each(lists)("%s keeps atom-selection equality without invoking selectedAtoms", (flag, list) => {
        const atoms = [true, 1, "1", false, 0, "0", 2, "true", null, undefined].map(value => ({ [flag]: value }));
        const m = model(atoms), v = viewer([m]);
        const expected = m.selectedAtoms({ [flag]: true });
        expect(expected).toEqual(atoms.slice(0, 3));
        const selection = jest.spyOn(m, "selectedAtoms");
        v.updateClickables();
        expect(v[list]).toEqual(expected);
        expected.forEach((atom, i) => expect(v[list][i]).toBe(atom));
        expect(selection).not.toHaveBeenCalled();
    });

    test.each(lists)("%s retains truthiness for shapes", (flag, list) => {
        const shapes = [true, 1, "1", 2, "true", "0", false, 0, "", null, undefined].map(value => ({ [flag]: value }));
        const v = viewer([], shapes);
        v.updateClickables();
        expect(v[list]).toEqual(shapes.slice(0, 6));
    });

    test("preserves model/atom/shape order, duplicates and hidden entries while skipping holes", () => {
        const a = enabled("a"), b = { name: "b", clickable: true, contextMenuEnabled: true };
        const c = { name: "c", hoverable: true };
        const first = model([a, , null, b, a]), second = model([c]);
        first.hide();
        const shape = { ...enabled("shape"), hidden: true }, menu = { contextMenuEnabled: true };
        const v = viewer([null, first, , second, first], [undefined, shape, , menu, shape]);
        v.updateClickables();
        expect(v.clickables).toEqual([a, b, a, a, b, a, shape, shape]);
        expect(v.hoverables).toEqual([a, a, c, a, a, shape, shape]);
        expect(v.contextMenuEnabledObjects).toEqual([a, b, a, a, b, a, shape, menu, shape]);
    });

    test("observes direct flag edits, clears stale entries and retains the list objects", () => {
        const old = enabled("old"), current = { name: "current" };
        const m = model([old, current]), shape = enabled("shape"), v = viewer([m], [shape]);
        const references = lists.map(([, list]) => v[list]);
        v.updateClickables();
        for (const [flag] of lists) { old[flag] = false; current[flag] = true; shape[flag] = false; }
        v.updateClickables();
        lists.forEach(([, list], i) => {
            expect(v[list]).toBe(references[i]);
            expect(v[list]).toEqual([current]);
        });
        m.getInternalState().atoms.length = 0;
        v.shapes.length = 0;
        v.updateClickables();
        lists.forEach(([, list], i) => {
            expect(v[list]).toBe(references[i]);
            expect(v[list]).toEqual([]);
        });
    });

    test("frame changes use the current frame's atom objects", async () => {
        const first = enabled("first"), second = enabled("second"), frames = [[first], [second]];
        const m = model(frames[0], frames), v = viewer([m]);
        for (const frame of [0, 1, 0]) {
            await m.setFrame(frame);
            v.updateClickables();
            for (const [, list] of lists) {
                expect(v[list]).toHaveLength(1);
                expect(v[list][0]).toBe(frames[frame][0]);
            }
        }
    });
});
