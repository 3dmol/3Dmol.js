global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const outlineKeys = ["_outlineMaterial", "_outlineSphereImposterMaterial", "_outlineStickImposterMaterial"];
const outlines = renderer => outlineKeys.map(key => renderer[key]);
const references = renderer => renderer._programs.reduce((sum, entry) => sum + entry.usedTimes, 0);
const flush = async () => { await Promise.resolve(); await Promise.resolve(); };

function fixture() {
    // Obtain the real renderer through the viewer without allocating a WebGL context.
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: {}, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer;
    const programs = new Set();
    let constant = 100;
    // Exercise real shader caching and material disposal with tracked GPU handles.
    const gl = new Proxy({
        createProgram: jest.fn(() => { const program = {}; programs.add(program); return program; }),
        deleteProgram: jest.fn(program => { expect(programs.delete(program)).toBe(true); }),
        createShader: () => ({}),
        getShaderParameter: () => true,
        getProgramParameter: () => true,
        getUniformLocation: () => null,
        getAttribLocation: () => 0
    }, {
        get(target, key) {
            if (key in target) return target[key];
            if (/^[A-Z0-9_]+$/.test(key)) return target[key] = constant++;
            return target[key] = () => {};
        }
    });
    renderer._gl = gl;
    renderer._webglversion = 2;
    return { renderer, gl, programs };
}

function initialize(renderer, materials = outlines(renderer)) {
    for (const material of materials) {
        renderer.initMaterial(material, [], null, {});
        material.needsUpdate = false;
    }
}

function disposeCurrent(renderer) {
    for (const material of outlines(renderer)) material.dispose();
}

describe("outline material ownership", () => {
    test("unused replaced materials retire once at the next microtask", async () => {
        const { renderer, gl } = fixture();
        const old = outlines(renderer);
        const disposed = old.map(material => {
            const listener = jest.fn();
            material.addEventListener("dispose", listener);
            return listener;
        });
        renderer.enableOutline({ width: 0.3 });
        disposed.forEach(listener => expect(listener).not.toHaveBeenCalled());
        outlines(renderer).forEach((material, i) => expect(material).not.toBe(old[i]));
        await flush();
        disposed.forEach(listener => expect(listener).toHaveBeenCalledTimes(1));
        expect(gl.createProgram).not.toHaveBeenCalled();
        disposeCurrent(renderer);
    });

    test("repeated immediate redraws reuse all three shaders without retaining old references", async () => {
        const { renderer, gl, programs } = fixture();
        initialize(renderer);
        const handles = outlines(renderer).map(material => material.program);
        for (let i = 0; i < 20; i++) {
            const old = outlines(renderer);
            renderer.enableOutline({ width: 0.1 + i / 100, color: 0xabcdef });
            initialize(renderer);
            expect(outlines(renderer).map(material => material.program)).toEqual(handles);
            expect(gl.deleteProgram).not.toHaveBeenCalled();
            expect(references(renderer)).toBe(6);
            await flush();
            old.forEach(material => expect(material.program).toBeUndefined());
            expect(references(renderer)).toBe(3);
            expect(gl.createProgram).toHaveBeenCalledTimes(3);
        }
        expect(programs.size).toBe(3);
        disposeCurrent(renderer);
        expect(programs.size).toBe(0);
        expect(references(renderer)).toBe(0);
    });

    test("several replacements before redraw dispose every abandoned material", async () => {
        const { renderer, gl, programs } = fixture();
        initialize(renderer);
        const retired = [];
        for (let i = 0; i < 4; i++) {
            outlines(renderer).forEach(material => {
                const listener = jest.fn();
                material.addEventListener("dispose", listener);
                retired.push({ material, listener });
            });
            renderer.enableOutline({});
        }
        initialize(renderer);
        await flush();
        retired.forEach(({ material, listener }) => {
            expect(listener).toHaveBeenCalledTimes(1);
            expect(material.program).toBeUndefined();
        });
        expect(references(renderer)).toBe(3);
        expect(gl.createProgram).toHaveBeenCalledTimes(3);
        disposeCurrent(renderer);
        expect(programs.size).toBe(0);
    });

    test("disabling retains current materials while delayed replacements release unused shaders", async () => {
        const { renderer, gl, programs } = fixture();
        initialize(renderer);
        const old = outlines(renderer);
        renderer.disableOutline();
        await flush();
        expect(renderer._outlineEnabled).toBe(false);
        expect(outlines(renderer)).toEqual(old);
        expect(references(renderer)).toBe(3);
        expect(gl.deleteProgram).not.toHaveBeenCalled();
        renderer.enableOutline({});
        await flush();
        expect(programs.size).toBe(0);
        expect(references(renderer)).toBe(0);
        initialize(renderer);
        expect(gl.createProgram).toHaveBeenCalledTimes(6);
        expect(programs.size).toBe(3);
        disposeCurrent(renderer);
        expect(programs.size).toBe(0);
    });

    test("releasing outlines preserves a shader still used by another material", async () => {
        const { renderer, programs } = fixture();
        initialize(renderer);
        const other = renderer._outlineMaterial.clone();
        initialize(renderer, [other]);
        const shared = other.program;
        renderer.enableOutline({});
        await flush();
        expect(references(renderer)).toBe(1);
        expect(programs.has(shared)).toBe(true);
        expect(programs.size).toBe(1);
        initialize(renderer);
        expect(renderer._outlineMaterial.program).toBe(shared);
        disposeCurrent(renderer);
        expect(programs.size).toBe(1);
        other.dispose();
        expect(programs.size).toBe(0);
    });
});
