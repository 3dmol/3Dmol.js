const fs = require("fs");
const path = require("path");
const { Worker: NodeWorker } = require("worker_threads");
global.$ = require("jquery");
global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const coordinates = [[[0, 0, 0], [1.3, 0, 0]], [[2, 3, 4], [3.3, 3, 4]]];
const positions = model => model.selectedAtoms({}).map(({ x, y, z }) => [x, y, z]);
function modelWithTrajectory() {
    const model = new $3Dmol.GLModel(0);
    model.addMolData("2\ninitial\nC 0 0 0\nO 1.3 0 0\n", "xyz");
    model.setCoordinates(coordinates, "array");
    return model;
}

describe("trajectory coordinate validation", () => {
    let warn;
    beforeEach(() => { warn = jest.spyOn(console, "warn").mockImplementation(() => {}); });
    afterEach(() => { jest.restoreAllMocks(); });

    test("an empty model returns without hanging or changing state", async () => {
        // A worker timeout keeps the original infinite-loop regression from hanging Jest.
        const script = `
            const assert = require('assert/strict');
            const { parentPort, workerData } = require('worker_threads');
            global.window = globalThis;
            global.URL.createObjectURL = () => '';
            global.document = { readyState: 'loading' };
            const { GLModel } = require(workerData);
            const model = new GLModel(0);
            const before = model.getInternalState();
            const warnings = [];
            console.warn = (...args) => warnings.push(args);
            assert.deepEqual(model.setCoordinates([[[1, 2, 3]]], 'array'), []);
            assert.equal(warnings.length, 1);
            assert.strictEqual(model.getInternalState().atoms, before.atoms);
            assert.strictEqual(model.getInternalState().frames, before.frames);
            parentPort.postMessage('empty-model validation passed');
        `;
        const worker = new NodeWorker(script, {
            eval: true, workerData: require.resolve("../../build/3Dmol.js")
        });
        let timeout;
        try {
            const result = await new Promise((resolve, reject) => {
                timeout = setTimeout(() => reject(new Error("Empty-model coordinate loading did not return.")), 5000);
                worker.once("message", resolve);
                worker.once("error", reject);
            });
            expect(result).toBe("empty-model validation passed");
        } finally {
            clearTimeout(timeout);
            await worker.terminate();
        }
    }, 10000);

    test.each([
        ["empty trajectory", [], "array"],
        ["missing atom", [[[1, 2, 3]]], "array"],
        ["ragged triples with a valid flattened length", [[[1, 2], [3, 4, 5, 6]]], "array"],
        ["ragged frames after a complete frame", [coordinates[0], [[1, 2, 3]]], "array"],
        ["sparse trajectory", new Array(1), "array"],
        ["sparse atom list", [new Array(2)], "array"],
        ["sparse coordinates before a complete frame", [[new Array(3), new Array(3)], coordinates[0]], "array"],
        ["NaN coordinate", [[[NaN, 0, 0], [1, 2, 3]]], "array"],
        ["infinite coordinate", [[[Infinity, 0, 0], [1, 2, 3]]], "array"],
        ["numeric string coordinate", [[["1", 0, 0], [1, 2, 3]]], "array"],
        ["non-array input", "invalid", "array"],
        ["ragged array accepted by the MDCRD parser", [[[1, 2], [3, 4, 5, 6]]], "mdcrd"],
        ["incomplete text frame", "title\n1 2 3 4 5", "mdcrd"],
        ["nonfinite text coordinate", "title\n1 2 3 4 5 Infinity", "mdcrd"],
        ["text parser exception", "title\n", "mdcrd"],
        ["no parsed coordinates", "no atom records", "pdb"],
        ["binary parser exception", new ArrayBuffer(0), "netcdf"],
        ["decompression failure", new ArrayBuffer(0), "mdcrd.gz"]
    ])("%s warns and preserves the existing trajectory", async (_name, input, format) => {
        const model = modelWithTrajectory();
        await model.setFrame(1);
        const before = model.getInternalState();
        const snapshot = JSON.parse(JSON.stringify(before));

        expect(model.setCoordinates(input, format)).toEqual([]);
        expect(warn).toHaveBeenCalledTimes(1);
        expect(model.getInternalState().atoms).toBe(before.atoms);
        expect(model.getInternalState().frames).toBe(before.frames);
        expect(model.getInternalState()).toEqual(snapshot);
        expect(positions(model)).toEqual(coordinates[1]);
        await model.setFrame(0);
        expect(positions(model)).toEqual(coordinates[0]);
    });

    test.each([[null, "array"], ["", "mdcrd"], ["data", "unsupported"]])(
        "missing or unsupported input (%p, %s) remains a silent no-op", (input, format) => {
            const model = modelWithTrajectory();
            const before = model.getInternalState();
            expect(model.setCoordinates(input, format)).toEqual([]);
            expect(model.getInternalState().atoms).toBe(before.atoms);
            expect(model.getInternalState().frames).toBe(before.frames);
            expect(warn).not.toHaveBeenCalled();
        }
    );

    test.each(["array", "mdcrd"])("valid arrays load all frames with format %s", async format => {
        const model = modelWithTrajectory();
        const frames = model.setCoordinates(coordinates, format);
        expect(frames).toHaveLength(2);
        expect(positions(model)).toEqual(coordinates[0]);
        await model.setFrame(1);
        expect(positions(model)).toEqual(coordinates[1]);
        expect(model.selectedAtoms({}).map(atom => atom.elem)).toEqual(["C", "O"]);
        expect(model.selectedAtoms({})[0].bonds).toEqual([1]);
        expect(warn).not.toHaveBeenCalled();
    });

    test.each([
        ["TRPcage.prmtop", "prmtop", "TRPcage.inpcrd", "inpcrd", 1],
        ["TC5b.prmtop", "prmtop", "heat1.mdcrd", "mdcrd", 200],
        ["TC5b.prmtop", "prmtop", "heat1.mdcrd.gz", "mdcrd.gz", 200],
        ["4csv.pdb", "pdb", "4csv.pdb", "pdb", 1],
        ["model1.prmtop", "prmtop", "model1_md2.nc", "netcdf", 10]
    ])("loads %s (%s) with %s (%s)", async (topology, topFormat, file, format, frameCount) => {
        const fixture = name => path.join(__dirname, "../auto/data", name);
        const model = new $3Dmol.GLModel(0);
        model.addMolData(fs.readFileSync(fixture(topology), "utf8"), topFormat);
        const atomCount = model.selectedAtoms({}).length;
        const binary = format === "netcdf" || format.endsWith(".gz");
        let data = fs.readFileSync(fixture(file), binary ? undefined : "utf8");
        if (binary) data = data.buffer.slice(data.byteOffset, data.byteOffset + data.byteLength);

        const frames = model.setCoordinates(data, format);
        expect(frames).toHaveLength(frameCount);
        expect(frames.every(frame => frame.length === atomCount)).toBe(true);
        expect(warn).not.toHaveBeenCalled();
        await model.setFrame(frameCount - 1);
        expect(positions(model).every(point => point.every(Number.isFinite))).toBe(true);
    });
});
