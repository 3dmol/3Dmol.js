/** @jest-environment node */
const path = require("path");
const ts = require("typescript");

test("published coordinate types accept nested arrays and existing trajectory inputs", () => {
    const filename = path.join(__dirname, "coordinate-types-consumer.ts");
    const source = `
        import type { GLModel } from '../../build/types/GLModel';
        declare const model: GLModel;

        model.setCoordinates([[[0, 1, 2]], [[3, 4, 5]]], 'array');
        const inferred = [[[0, 1, 2]]];
        const explicit: number[][][] = [[[0, 1, 2]]];
        const tuples: [number, number, number][][] = [[[0, 1, 2]]];
        model.setCoordinates(inferred, 'array');
        model.setCoordinates(explicit, 'array');
        model.setCoordinates(tuples, 'array');
        model.setCoordinates('trajectory text', 'mdcrd');
        model.setCoordinates(new ArrayBuffer(24), 'netcdf');
        declare const input: string | ArrayBuffer | number[][][];
        model.setCoordinates(input, 'array');

        // @ts-expect-error coordinates need frames x atoms x coordinates
        model.setCoordinates([[0, 1, 2]], 'array');
        // @ts-expect-error coordinates have three array dimensions
        model.setCoordinates([[[[0, 1, 2]]]], 'array');
        // @ts-expect-error coordinate values must be numbers
        model.setCoordinates([[['0', '1', '2']]], 'array');
        // @ts-expect-error all coordinate values must be numbers
        model.setCoordinates([[[0, 1, '2']]], 'array');
        // @ts-expect-error coordinate points are arrays, not objects
        model.setCoordinates([[{ x: 0, y: 1, z: 2 }]], 'array');
        // @ts-expect-error a scalar is not trajectory data
        model.setCoordinates(42, 'array');
        // @ts-expect-error flat typed arrays are not nested coordinate arrays
        model.setCoordinates(new Float32Array([0, 1, 2]), 'array');
        // @ts-expect-error a format is still required
        model.setCoordinates(inferred);
    `;
    const options = {
        noEmit: true,
        strict: true,
        skipLibCheck: true,
        target: ts.ScriptTarget.ES2015,
        module: ts.ModuleKind.ESNext,
        moduleResolution: ts.ModuleResolutionKind.Node10,
        baseUrl: path.resolve(__dirname, "../../build/types"),
        types: []
    };
    const host = ts.createCompilerHost(options);
    const readFile = host.readFile;
    host.readFile = file => file === filename ? source : readFile(file);
    const program = ts.createProgram([filename], options, host);
    const diagnostics = ts.getPreEmitDiagnostics(program);
    expect(ts.formatDiagnostics(diagnostics, host)).toBe("");
});
