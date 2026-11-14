global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function fixture() {
    // Obtain the real renderer without allocating a browser WebGL context.
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: {}, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer;
    const programs = new Set(), shaders = new Set();
    const state = { failure: "", failProgram: 0, framebuffer: null };
    let constant = 100, programCount = 0;
    const validProgram = program => {
        expect(programs.has(program)).toBe(true);
        expect(program.linked).toBe(true);
    };
    const gl = new Proxy({
        drawingBufferWidth: 640, drawingBufferHeight: 480,
        isContextLost: () => false,
        createProgram: jest.fn(() => {
            if (state.failure === "program allocation") return null;
            const program = { number: ++programCount, attached: new Set() };
            programs.add(program);
            return program;
        }),
        deleteProgram: jest.fn(program => expect(programs.delete(program)).toBe(true)),
        createShader: jest.fn(type => {
            const stage = type === gl.FRAGMENT_SHADER ? "fragment" : "vertex";
            if (state.failure === stage + " allocation") return null;
            const shader = { stage };
            shaders.add(shader);
            return shader;
        }),
        deleteShader: jest.fn(shader => expect(shaders.delete(shader)).toBe(true)),
        compileShader: shader => { shader.compiled = state.failure !== shader.stage + " compilation"; },
        getShaderParameter: shader => shader.compiled,
        getShaderInfoLog: () => "controlled compilation failure",
        attachShader: (program, shader) => program.attached.add(shader),
        detachShader: (program, shader) => program.attached.delete(shader),
        linkProgram: program => {
            program.linked = state.failure !== "linking" && program.number !== state.failProgram &&
                program.attached.size === 2 && [...program.attached].every(shader => shader.compiled);
        },
        getProgramParameter: program => program.linked,
        getProgramInfoLog: () => "controlled link failure",
        getUniformLocation: jest.fn(program => { validProgram(program); return null; }),
        getAttribLocation: jest.fn(program => { validProgram(program); return 0; }),
        useProgram: jest.fn(validProgram),
        createBuffer: () => ({}), createTexture: () => ({}), createFramebuffer: () => ({}),
        bindFramebuffer: jest.fn((target, framebuffer) => { state.framebuffer = framebuffer; })
    }, {
        get(target, key) {
            if (key in target) return target[key];
            if (/^[A-Z0-9_]+$/.test(key)) return target[key] = constant++;
            return target[key] = jest.fn();
        }
    });
    renderer._gl = gl;
    renderer._webglversion = 2;
    return { renderer, gl, programs, shaders, state };
}

function build(renderer) {
    return renderer.buildProgram(
        "void main() { gl_FragColor = vec4(1.0); }",
        "void main() { gl_Position = vec4(0.0); }", {}, {}
    );
}

function mesh() {
    const shape = new $3Dmol.GLShape({});
    shape.addCustom({
        vertexArr: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }],
        faceArr: [0, 1, 2]
    });
    let root;
    shape.globj({ add(object) { root = object; } });
    const object = root.children[0];
    object._modelViewMatrix = new $3Dmol.Matrix4();
    return object;
}

describe("shader failure recovery", () => {
    let errors;
    beforeEach(() => { errors = jest.spyOn(console, "error").mockImplementation(() => {}); });
    afterEach(() => { errors.mockRestore(); });

    test.each([
        "fragment compilation", "vertex compilation", "linking",
        "program allocation", "fragment allocation", "vertex allocation"
    ])("%s releases resources without querying or caching an invalid program", failure => {
        const { renderer, gl, programs, shaders, state } = fixture();
        state.failure = failure;
        expect(build(renderer)).toBeNull();
        expect(errors).toHaveBeenCalled();
        expect(programs.size).toBe(0);
        expect(shaders.size).toBe(0);
        expect(renderer._programs).toHaveLength(0);
        expect(renderer.info.memory.programs).toBe(0);
        expect(gl.getUniformLocation).not.toHaveBeenCalled();
        expect(gl.getAttribLocation).not.toHaveBeenCalled();
    });

    test("successful programs still share handles until the last reference is released", () => {
        const { renderer, gl, programs, shaders } = fixture();
        const program = build(renderer);
        expect(build(renderer)).toBe(program);
        expect(gl.createShader).toHaveBeenCalledTimes(2);
        expect(shaders.size).toBe(0);
        expect(renderer._programs[0].usedTimes).toBe(2);
        renderer.deallocateMaterial({ program, id: 1 });
        expect(programs.size).toBe(1);
        renderer.deallocateMaterial({ program, id: 2 });
        expect(programs.size).toBe(0);
        expect(renderer.info.memory.programs).toBe(0);
    });

    test("a failed material skips drawing until needsUpdate requests a retry", () => {
        const { renderer, gl, programs, state } = fixture();
        const object = mesh(), material = object.material, camera = {};
        const fog = { color: { r: 1, g: 1, b: 1 }, near: 1, far: 100 };
        state.failure = "fragment compilation";
        renderer.renderBuffer(camera, [], fog, material, {}, object);
        expect(material.program).toBeNull();
        expect(material.needsUpdate).toBe(false);
        for (let i = 0; i < 3; i++) renderer.renderBuffer(camera, [], fog, material, {}, object);
        expect(gl.createProgram).toHaveBeenCalledTimes(1);
        expect(gl.useProgram).not.toHaveBeenCalled();
        expect(gl.drawElements).not.toHaveBeenCalled();
        expect(gl.drawArrays).not.toHaveBeenCalled();

        state.failure = "";
        material.needsUpdate = true;
        expect(renderer.setProgram(camera, [], fog, material, object, renderer)).toBeTruthy();
        expect(gl.createProgram).toHaveBeenCalledTimes(2);
        expect(programs.size).toBe(1);
        material.dispose();
        expect(programs.size).toBe(0);
    });

    test.each([
        [1, "_screenshader"], [2, "_AOshader"], [3, "_blurshader"]
    ])("postprocessing failure %i leaves other programs usable and skips AO", (number, key) => {
        const { renderer, gl, programs, shaders, state } = fixture();
        state.failProgram = number;
        renderer.setSize(100, 80);
        expect(renderer[key]).toBeNull();
        expect(programs.size).toBe(2);
        expect(shaders.size).toBe(0);
        expect(renderer._programs).toHaveLength(2);
        renderer.setShading({}, {}, "opaque");
        expect(gl.framebufferTexture2D).not.toHaveBeenCalled();
        expect(gl.useProgram).not.toHaveBeenCalled();

        const object = mesh();
        object.material.hasAO = true;
        const record = { object, buffer: {}, render: true };
        renderer.unrollBufferMaterial(record);
        renderer.renderBuffer = jest.fn();
        renderer.renderObjects([record], true, "opaque", {}, [], {}, false);
        expect(renderer.renderBuffer).toHaveBeenCalledTimes(1);
        expect(renderer.renderBuffer.mock.calls[0][3]).toBe(object.material);
    });

    test.each([false, true])("screen failure preserves direct drawing (grid: %s)", grid => {
        const { renderer, gl, state } = fixture();
        if (grid) Object.assign(renderer, { rows: 2, cols: 2, row: 1, col: 1 });
        state.failProgram = 1;
        renderer.setSize(100, 80);
        renderer.setFrameBuffer();
        expect(state.framebuffer).toBeNull();
        const viewport = grid ? [320, 240, 320, 240] : [0, 0, 640, 480];
        expect(gl.viewport).toHaveBeenLastCalledWith(...viewport);
        expect(gl.scissor).toHaveBeenLastCalledWith(...viewport);
        expect(gl.framebufferTexture2D).not.toHaveBeenCalled();
        renderer.renderFrameBuffertoScreen();
        expect(gl.clear).not.toHaveBeenCalled();
        expect(gl.useProgram).not.toHaveBeenCalled();
        expect(gl.drawArrays).not.toHaveBeenCalled();
    });
});
