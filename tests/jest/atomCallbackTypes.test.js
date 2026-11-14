/** @jest-environment node */
const path = require("path");
const ts = require("typescript");

test("published atom callback types describe the complete runtime arguments", () => {
    const filename = path.join(__dirname, "atom-callback-types-consumer.ts");
    const source = `
        import type { AtomSpec, PickingIntersection } from '../../build/types/specs';
        import type { GLViewer } from '../../build/types/GLViewer';
        import type { GLShape } from '../../build/types/GLShape';
        declare const viewer: GLViewer;
        declare const shape: GLShape;
        declare const mouse: MouseEvent;
        declare const touch: TouchEvent;
        declare const container: HTMLElement;

        const atom: AtomSpec = {
            callback(atom, viewer, event, container, intersects) {
                const x: number | undefined = atom.x;
                viewer.render();
                if ('changedTouches' in event) event.changedTouches.item(0);
                else event.clientX.toFixed();
                container?.getBoundingClientRect();
                const target: AtomSpec | GLShape = intersects[0].clickable;
                intersects[0].distance.toFixed();
                // @ts-expect-error touch events do not have clientX
                event.clientX;
                // @ts-expect-error the container can be null
                container.getBoundingClientRect();
                // @ts-expect-error intersections are records, not raw atoms
                intersects[0].serial;
            },
            hover_callback(atom, viewer, event, container, intersects) {
                event.preventDefault();
                const hit: PickingIntersection = intersects[0];
            },
            unhover_callback(atom, viewer, event, container, intersects) {
                event?.preventDefault();
                intersects?.forEach(hit => hit.distance.toFixed());
                // @ts-expect-error unhover can lack an event
                event.preventDefault();
                // @ts-expect-error unhover can lack intersections
                intersects.length;
            }
        };
        const noArguments: AtomSpec = {
            callback: () => {}, hover_callback: () => {}, unhover_callback: () => {}
        };
        const twoArguments: AtomSpec = {
            callback: (atom, viewer) => viewer.render(),
            hover_callback: (atom, viewer) => viewer.render(),
            unhover_callback: (atom, viewer) => viewer.render()
        };
        const hits: PickingIntersection[] = [
            { clickable: atom, distance: 1 }, { clickable: shape, distance: 2 }
        ];
        atom.callback!(atom, viewer, mouse, container, hits);
        atom.callback!(atom, viewer, touch, null, hits);
        atom.hover_callback!(atom, viewer, touch, container, hits);
        atom.unhover_callback!(atom, viewer, undefined, container, undefined);
        // @ts-expect-error click callbacks always receive an event
        atom.callback!(atom, viewer, undefined, container, hits);
        // @ts-expect-error intersections cannot refer to a string
        const badHit: PickingIntersection = { clickable: 'atom', distance: 1 };
        // @ts-expect-error incorrect event type
        const badCallback: AtomSpec = { callback: (atom, viewer, event: string) => {} };
        // @ts-expect-error callbacks cannot require a sixth argument
        const tooManyArguments: AtomSpec = { callback: (atom, viewer, event, container, hits, extra: number) => {} };
        // @ts-expect-error unhover handlers must accept absent events
        const narrowUnhover: AtomSpec = { unhover_callback: (atom, viewer, event: MouseEvent | TouchEvent) => {} };
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
