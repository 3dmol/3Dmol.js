global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

let raycaster;

beforeAll(() => {
    // Exercise the viewer's real raycaster without allocating rendering resources.
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"),
            getAspect: (width, height) => width / height,
            setSize() {},
            setClearColorHex() {}
        };
    });
    const container = document.createElement("div");
    Object.defineProperties(container, {
        offsetWidth: { value: 500 }, offsetHeight: { value: 500 }
    });
    try {
        raycaster = new $3Dmol.GLViewer(container, {}).raycaster;
    } finally {
        setup.mockRestore();
    }
});

function line(points) {
    return {
        intersectionShape: {
            line: points.map(point => new $3Dmol.Vector3(...point)),
            sphere: [], cylinder: [], triangle: []
        }
    };
}

function pick(objects) {
    raycaster.set(new $3Dmol.Vector3(0, 0, 0), new $3Dmol.Vector3(0, 0, 1));
    return raycaster.intersectObjects({ matrixWorld: new $3Dmol.Matrix4() }, objects);
}

describe("line picking", () => {
    test.each([
        ["before the first endpoint", [[1, 0, 5], [3, 0, 5]]],
        ["after the second endpoint", [[-3, 0, 5], [-1, 0, 5]]],
        ["behind the ray origin", [[-1, 0, -5], [1, 0, -5]]]
    ])("rejects hits %s in either endpoint order", (_name, points) => {
        expect(pick([line(points)])).toEqual([]);
        expect(pick([line([...points].reverse())])).toEqual([]);
    });

    test.each([
        ["interior", [[-1, 0, 5], [1, 0, 5]]],
        ["first endpoint", [[0, 0, 5], [2, 0, 5]]],
        ["second endpoint", [[-2, 0, 5], [0, 0, 5]]]
    ])("accepts the %s in either endpoint order", (_name, points) => {
        for (const ordered of [points, [...points].reverse()]) {
            const object = line(ordered);
            const hits = pick([object]);
            expect(hits).toHaveLength(1);
            expect(hits[0].clickable).toBe(object);
            expect(hits[0].distance).toBeCloseTo(5);
        }
    });

    test("preserves line precision and excludes zero-length segments", () => {
        const nearby = line([[-1, 0.1, 5], [1, 0.1, 5]]);
        const distant = line([[-1, 0.3, 5], [1, 0.3, 5]]);
        const point = line([[0, 0, 5], [0, 0, 5]]);
        expect(pick([nearby, distant, point]).map(hit => hit.clickable)).toEqual([nearby]);
    });

    test("sorts forward hits without selecting a line behind the camera", () => {
        const near = line([[-1, 0, 3], [1, 0, 3]]);
        const far = line([[-1, 0, 8], [1, 0, 8]]);
        const behind = line([[-1, 0, -2], [1, 0, -2]]);
        const hits = pick([far, behind, near]);
        expect(hits.map(hit => hit.clickable)).toEqual([near, far]);
        expect(hits.map(hit => hit.distance)).toEqual([3, 8]);
    });
});
