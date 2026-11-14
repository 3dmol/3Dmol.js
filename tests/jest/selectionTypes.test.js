/** @jest-environment node */
const path = require("path");
const ts = require("typescript");

test("published selection types accept attribute lists while preserving atom and control types", () => {
    const filename = path.join(__dirname, "selection-types-consumer.ts");
    const source = `
        import type { AtomSelectionSpec, AtomSpec } from '../../build/types/specs';
        import type { GLModel } from '../../build/types/GLModel';
        import type { GLViewer } from '../../build/types/GLViewer';
        declare const model: GLModel;
        declare const viewer: GLViewer;

        const scalar: AtomSelectionSpec = { chain: 'A', serial: 1, hetflag: false };
        const lists: AtomSelectionSpec = {
            chain: ['A', 'B'], elem: ['C', 'O'], resn: ['GLY', 'ALA'], atom: ['CA'],
            serial: [1, 2], x: [0, 1], b: [20, 30], hetflag: [true, false],
            clickable: [true], contextMenuEnabled: [true, false],
            color: ['red', 0xff0000], surfaceColor: [0xff0000, 'blue']
        };
        const controls: AtomSelectionSpec = {
            model: [model], frame: 0, index: [1, 2], bonds: 0, resi: [1, '5-10'],
            predicate: atom => atom.elem === 'C', invert: false, byres: true, expand: 3,
            within: { distance: 4, sel: { chain: ['B', 'C'] } },
            and: [{ elem: ['C', 'O'] }], or: [{ hetflag: [false] }], not: { resn: ['HOH'] }
        };
        const structured: AtomSelectionSpec = {
            bondOrder: [1, 2], bondStyles: [{ color1: 'red' }], properties: { score: 1 },
            style: { stick: {} }, symmetries: [], color: { r: 1, g: 0, b: 0 }
        };
        viewer.selectedAtoms(scalar);
        model.selectedAtoms(lists);
        viewer.setStyle({ chain: ['A', 'B'] }, { stick: {} });
        model.setClickable({ elem: ['C', 'O'] }, true, () => undefined);
        viewer.addSurface('VDW', {}, controls);
        model.selectedAtoms(structured);

        // @ts-expect-error chain values are strings
        const badChain: AtomSelectionSpec = { chain: [1] };
        // @ts-expect-error serial values are numbers
        const badSerial: AtomSelectionSpec = { serial: ['one'] };
        // @ts-expect-error hetflag values are booleans
        const badHetflag: AtomSelectionSpec = { hetflag: ['true'] };
        // @ts-expect-error invert is one boolean
        const badInvert: AtomSelectionSpec = { invert: [true] };
        // @ts-expect-error frame is one number
        const badFrame: AtomSelectionSpec = { frame: [0, 1] };
        // @ts-expect-error byres is one boolean
        const badByres: AtomSelectionSpec = { byres: [true] };
        // @ts-expect-error expand is one distance
        const badExpand: AtomSelectionSpec = { expand: [1, 2] };
        // @ts-expect-error bonds is one bond count
        const badBonds: AtomSelectionSpec = { bonds: [0, 1] };
        // @ts-expect-error array fields do not gain nested-list semantics
        const badBondOrder: AtomSelectionSpec = { bondOrder: [[1, 2]] };
        // @ts-expect-error object fields do not gain list semantics
        const badStyle: AtomSelectionSpec = { style: [{ stick: {} }] };
        // @ts-expect-error callback lists are not primitive attribute lists
        const badCallback: AtomSelectionSpec = { callback: [() => undefined] };
        // @ts-expect-error only primitive color values gain list semantics
        const badColor: AtomSelectionSpec = { color: [{ r: 1, g: 0, b: 0 }] };
        // @ts-expect-error atom data retains scalar strings
        const badAtomChain: AtomSpec = { chain: ['A', 'B'] };
        // @ts-expect-error atom data retains scalar numbers
        const badAtomSerial: AtomSpec = { serial: [1, 2] };
        // @ts-expect-error atom data retains scalar booleans
        const badAtomHetflag: AtomSpec = { hetflag: [true] };
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
