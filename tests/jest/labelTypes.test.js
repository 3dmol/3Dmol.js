/** @jest-environment node */
const path = require("path");
const ts = require("typescript");

test("published label types accept supported style options", () => {
    const filename = path.join(__dirname, "label-types-consumer.ts");
    const source = `
        import type { LabelSpec, Label } from '../../build/types/Label';
        import type { GLViewer } from '../../build/types/GLViewer';
        import { RWB, GradientType } from '../../build/types/Gradient';
        declare const viewer: GLViewer;
        declare const label: Label;
        declare const canvas: HTMLCanvasElement;

        const bold: LabelSpec = { bold: true };
        const padding: LabelSpec = { padding: 8 };
        const width: LabelSpec = { backgroundWidth: 160 };
        const height: LabelSpec = { backgroundHeight: 40 };
        const gradientSpec: LabelSpec = { backgroundGradient: { gradient: 'RWB', min: -1, max: 1 } };
        const gradientObject: LabelSpec = { backgroundGradient: new RWB(-1, 1) };
        const combined: LabelSpec = {
            bold: false, padding: 0, backgroundWidth: 0, backgroundHeight: 0,
            backgroundImage: canvas,
            backgroundGradient: { gradient: 'linear', colors: ['red', 'blue'] },
            fontSize: 18, fontColor: 'white', position: { x: 1, y: 2, z: 3 }
        };
        class CustomGradient extends GradientType {
            range() { return [0, 1]; }
            valueToHex(value: number) { return value > 0.5 ? 0xffffff : 0x000000; }
        }
        viewer.addLabel('combined', combined);
        viewer.addLabel('inline', { bold: true, padding: 6, backgroundWidth: 100, backgroundHeight: 40 });
        viewer.addLabel('custom gradient', { backgroundGradient: new CustomGradient() });
        viewer.setLabelStyle(label, { backgroundGradient: { gradient: 'linear_red_blue' } });
        viewer.addResLabels({}, { bold: true });
        viewer.addPropertyLabels('b', {}, { padding: 5 });
        const previous: LabelSpec = { font: 'serif', useScreen: true, frame: 0, alignment: 'center' };

        // @ts-expect-error bold is boolean
        const badBold: LabelSpec = { bold: 'true' };
        // @ts-expect-error padding is numeric
        const badPadding: LabelSpec = { padding: '4' };
        // @ts-expect-error width is numeric
        const badWidth: LabelSpec = { backgroundWidth: '160' };
        // @ts-expect-error height is numeric
        const badHeight: LabelSpec = { backgroundHeight: false };
        // @ts-expect-error gradients require a specification or gradient object
        const badGradient: LabelSpec = { backgroundGradient: 'RWB' };
        // @ts-expect-error gradient range endpoints are numeric
        const badRange: LabelSpec = { backgroundGradient: { gradient: 'RWB', min: 'zero' } };
        // @ts-expect-error unsupported properties remain errors
        const unknown: LabelSpec = { madeUp: true };
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
