/** @jest-environment node */
const path = require("path");
const ts = require("typescript");

test("published label coordinate types accept plain points and vectors", () => {
    const filename = path.join(__dirname, "label-coordinate-types-consumer.ts");
    const source = `
        import type { LabelSpec, Label } from '../../build/types/Label';
        import type { GLViewer } from '../../build/types/GLViewer';
        import { Vector2 } from '../../build/types/WebGL/math/Vector2';
        declare const viewer: GLViewer;
        declare const label: Label;

        const plain: LabelSpec = { screenOffset: { x: 10, y: -15 }, alignment: { x: 0.5, y: -1 } };
        const vectors: LabelSpec = { screenOffset: new Vector2(10, -15), alignment: new Vector2(0.5, -1) };
        const preset: LabelSpec = { alignment: 'topLeft' };
        viewer.addLabel('plain', { screenOffset: { x: 0, y: 0 }, alignment: { x: 1, y: -1 } });
        viewer.addLabel('vectors', vectors);
        viewer.setLabelStyle(label, { screenOffset: { x: 10, y: 15 } });
        viewer.addResLabels({}, { screenOffset: { x: 0, y: 0 } });
        viewer.addPropertyLabels('b', {}, { alignment: { x: 0, y: 0 } });

        const stored: LabelSpec = Math.random() > 0.5 ? plain : vectors;
        if (stored.screenOffset instanceof Vector2) stored.screenOffset.clone();
        if (stored.alignment instanceof Vector2) stored.alignment.clone();
        if (stored.screenOffset) {
            stored.screenOffset.x.toFixed();
            // @ts-expect-error a plain coordinate object need not have Vector2 methods
            stored.screenOffset.clone();
        }
        if (stored.alignment && typeof stored.alignment !== 'string') stored.alignment.y.toFixed();

        // @ts-expect-error both offset coordinates are required
        const missingOffset: LabelSpec = { screenOffset: { x: 1 } };
        // @ts-expect-error offset coordinates must be numbers
        const badOffset: LabelSpec = { screenOffset: { x: '1', y: 2 } };
        // @ts-expect-error screenOffset is two-dimensional
        const extraOffset: LabelSpec = { screenOffset: { x: 1, y: 2, z: 3 } };
        // @ts-expect-error both alignment coordinates are required
        const missingAlignment: LabelSpec = { alignment: { y: 1 } };
        // @ts-expect-error alignment coordinates must be numbers
        const badAlignment: LabelSpec = { alignment: { x: 1, y: '2' } };
        // @ts-expect-error scalar values are not alignment vectors
        const scalarAlignment: LabelSpec = { alignment: 1 };
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
