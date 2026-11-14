/** @jest-environment node */
const path = require("path");
const ts = require("typescript");

test("published volume transfer types describe supported color and numeric inputs", () => {
    const filename = path.join(__dirname, "volume-transfer-types-consumer.ts");
    const source = `
        import type { VolumetricRendererSpec } from '../../build/types/VolumetricRender';
        import type { GLViewer } from '../../build/types/GLViewer';
        import type { VolumeData } from '../../build/types/VolumeData';
        import { Color, ColorSpec } from '../../build/types/colors';
        declare const viewer: GLViewer;
        declare const data: VolumeData;

        const numeric: VolumetricRendererSpec = { transferfn: [
            { color: 'red', opacity: 0, value: -1 },
            { color: 0x00ff00, opacity: 1, value: 1 }
        ] };
        const strings: VolumetricRendererSpec = { transferfn: [
            { color: '#000000', opacity: '0', value: '-1.5' },
            { color: 'rgb(255,255,255)', opacity: '1', value: '1.5' }
        ] };
        const objects: VolumetricRendererSpec = { transferfn: [
            { color: new Color(0xffffff), opacity: 0.5, value: 0 },
            { color: { r: 1, g: 0, b: 0 }, opacity: 1, value: 2 },
            { color: { r: 0, g: 1, b: 0, a: 1 }, opacity: '0.75', value: '3' }
        ] };
        const omitted: VolumetricRendererSpec = {};
        const empty: VolumetricRendererSpec = { transferfn: [] };
        viewer.addVolumetricRender(data, numeric);
        viewer.addVolumetricRender(data, strings);
        viewer.addVolumetricRender(data, objects);
        viewer.addVolumetricData('volume file', 'cube', numeric);

        const stop = numeric.transferfn![0];
        const color: ColorSpec = stop.color;
        const opacity: number | string = stop.opacity;
        const value: number | string = stop.value;
        if (typeof stop.value === 'number') stop.value.toFixed();
        else stop.value.trim();

        // @ts-expect-error color must be a supported color representation
        const badColor: VolumetricRendererSpec = { transferfn: [{ color: true, opacity: 1, value: 0 }] };
        // @ts-expect-error one point requires one color rather than a list
        const colorList: VolumetricRendererSpec = { transferfn: [{ color: ['red'], opacity: 1, value: 0 }] };
        // @ts-expect-error RGB colors require all three channels
        const missingChannel: VolumetricRendererSpec = { transferfn: [{ color: { r: 1, g: 0 }, opacity: 1, value: 0 }] };
        // @ts-expect-error opacity must be numeric or string
        const badOpacity: VolumetricRendererSpec = { transferfn: [{ color: 'red', opacity: true, value: 0 }] };
        // @ts-expect-error value must be numeric or string
        const badValue: VolumetricRendererSpec = { transferfn: [{ color: 'red', opacity: 1, value: {} }] };
        // @ts-expect-error numeric arrays are not a stop value
        const arrayValue: VolumetricRendererSpec = { transferfn: [{ color: 'red', opacity: 1, value: [0] }] };
        // @ts-expect-error opacity remains required for each point
        const missingOpacity: VolumetricRendererSpec = { transferfn: [{ color: 'red', value: 0 }] };
        declare const unchecked: unknown;
        // @ts-expect-error unchecked input must be validated before use
        const uncheckedValue: VolumetricRendererSpec = { transferfn: [{ color: 'red', opacity: 1, value: unchecked }] };
        if (typeof unchecked === 'number' || typeof unchecked === 'string') {
            const checkedValue: VolumetricRendererSpec = { transferfn: [{ color: 'red', opacity: 1, value: unchecked }] };
        }
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
