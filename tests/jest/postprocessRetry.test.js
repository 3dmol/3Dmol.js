global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const programKeys = ["_screenshader", "_AOshader", "_blurshader"];
const resourceKeys = ["_fb", "_screenQuadVBO", "_targetTexture", "_depthTexture", "_shadingTexture", "_scratchTexture"];

function fixture(version = 2) {
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: {}, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer;
    const programs = new Set(), shaders = new Set();
    const state = { failProgram: 0, failLink: false, framebuffer: null };
    let constant = 100, programCount = 0;
    const validProgram = program => {
        expect(programs.has(program)).toBe(true);
        expect(program.linked).toBe(true);
    };
    const gl = new Proxy({
        drawingBufferWidth: 640, drawingBufferHeight: 480,
        isContextLost: () => false,
        createProgram: jest.fn(() => {
            const program = { number: ++programCount, attached: new Set() };
            programs.add(program);
            return program;
        }),
        deleteProgram: jest.fn(program => expect(programs.delete(program)).toBe(true)),
        createShader: () => { const shader = {}; shaders.add(shader); return shader; },
        deleteShader: shader => expect(shaders.delete(shader)).toBe(true),
        compileShader: shader => { shader.compiled = true; },
        getShaderParameter: shader => shader.compiled,
        attachShader: (program, shader) => program.attached.add(shader),
        detachShader: (program, shader) => program.attached.delete(shader),
        linkProgram: program => {
            program.linked = !state.failLink && program.number !== state.failProgram &&
                program.attached.size === 2 && [...program.attached].every(shader => shader.compiled);
        },
        getProgramParameter: program => program.linked,
        getProgramInfoLog: () => "controlled link failure",
        getUniformLocation: program => { validProgram(program); return null; },
        getAttribLocation: jest.fn(program => { validProgram(program); return 0; }),
        useProgram: jest.fn(validProgram),
        createBuffer: jest.fn(() => ({})),
        createTexture: jest.fn(() => ({})),
        createFramebuffer: jest.fn(() => ({})),
        bindFramebuffer: (target, framebuffer) => { state.framebuffer = framebuffer; }
    }, {
        get(target, key) {
            if (key in target) return target[key];
            if (/^[A-Z0-9_]+$/.test(key)) return target[key] = constant++;
            return target[key] = jest.fn();
        }
    });
    renderer._gl = gl;
    renderer._webglversion = version;
    return { renderer, gl, programs, shaders, state };
}

function camera() {
    // Use the public viewer to obtain a real Camera for Renderer.render's type check.
    const setup = jest.spyOn($3Dmol.GLViewer.prototype, "setupRenderer").mockImplementation(function () {
        this.renderer = {
            domElement: document.createElement("canvas"), getAspect: () => 1,
            setSize() {}, setClearColorHex() {}, setViewport() {}, supportedExtensions: () => ({})
        };
    });
    try {
        return new $3Dmol.GLViewer(document.createElement("div"), { nomouse: true }).camera;
    } finally {
        setup.mockRestore();
    }
}

function expectResources(env, handles) {
    resourceKeys.forEach((key, index) => expect(env.renderer[key]).toBe(handles[index]));
    expect(env.gl.createTexture).toHaveBeenCalledTimes(4);
    expect(env.gl.createFramebuffer).toHaveBeenCalledTimes(1);
    expect(env.gl.createBuffer).toHaveBeenCalledTimes(1);
}

function expectPrograms(env, count) {
    expect(env.programs.size).toBe(count);
    expect(env.renderer._programs).toHaveLength(count);
    expect(env.renderer.info.memory.programs).toBe(count);
    for (const program of env.renderer._programs) expect(program.usedTimes).toBe(1);
    expect(env.shaders.size).toBe(0);
}

describe("postprocessing shader retries", () => {
    let errors;
    beforeEach(() => { errors = jest.spyOn(console, "error").mockImplementation(() => {}); });
    afterEach(() => { errors.mockRestore(); });

    test.each(programKeys)("resize recovers %s while retaining successful programs and resources", key => {
        const env = fixture(), { renderer, gl, state } = env;
        state.failProgram = programKeys.indexOf(key) + 1;
        renderer.setSize(100, 80);
        const handles = resourceKeys.map(name => renderer[name]);
        const successful = programKeys.filter(name => name !== key).map(name => [name, renderer[name]]);
        expect(renderer[key]).toBeNull();
        expectPrograms(env, 2);

        state.failProgram = 0;
        renderer.setSize(200, 120);
        expect(renderer[key].linked).toBe(true);
        for (const [name, program] of successful) expect(renderer[name]).toBe(program);
        expect(gl.createProgram).toHaveBeenCalledTimes(4);
        expectResources(env, handles);
        expectPrograms(env, 3);
        expect(renderer.canRenderAmbientOcclusion()).toBe(true);
        renderer.setFrameBuffer();
        renderer.renderFrameBuffertoScreen();
        expect(gl.useProgram).toHaveBeenCalledWith(renderer._screenshader);

        const queries = gl.getAttribLocation.mock.calls.length;
        renderer.initFrameBuffer();
        expect(gl.getAttribLocation).toHaveBeenCalledTimes(queries);
        expect(gl.createProgram).toHaveBeenCalledTimes(4);
        expectPrograms(env, 3);
    });

    test.each(programKeys)("persistent %s failure retains fallbacks without retrying during render", key => {
        const env = fixture(), { renderer, gl, state } = env;
        state.failProgram = programKeys.indexOf(key) + 1;
        renderer.setSize(100, 80);
        const handles = resourceKeys.map(name => renderer[name]);
        const successful = programKeys.filter(name => name !== key).map(name => [name, renderer[name]]);
        state.failLink = true;
        for (let i = 0; i < 5; i++) {
            renderer.setSize(200 + i, 120 + i);
            expect(renderer[key]).toBeNull();
            for (const [name, program] of successful) expect(renderer[name]).toBe(program);
            expectResources(env, handles);
            expectPrograms(env, 2);
        }
        expect(gl.createProgram).toHaveBeenCalledTimes(8);
        expect(renderer.canRenderAmbientOcclusion()).toBe(false);
        renderer.setFrameBuffer();
        expect(state.framebuffer).toBe(key === "_screenshader" ? null : renderer._fb);

        renderer.autoUpdateScene = false;
        renderer.autoUpdateObjects = false;
        renderer.renderObjects = jest.fn();
        renderer.renderSprites = jest.fn();
        renderer.render({ __lights: [], __webglObjects: [], fog: {} }, camera());
        expect(renderer.renderObjects).toHaveBeenCalled();
        expect(gl.createProgram).toHaveBeenCalledTimes(8);
        expectPrograms(env, 2);
    });

    test("successful resizes preserve every handle and program reference count", () => {
        const env = fixture(), { renderer, gl } = env;
        renderer.setSize(100, 80);
        const handles = resourceKeys.map(name => renderer[name]);
        const programs = programKeys.map(name => renderer[name]);
        for (let i = 0; i < 10; i++) renderer.setSize(200 + i, 120 + i);
        expectResources(env, handles);
        programKeys.forEach((key, index) => expect(renderer[key]).toBe(programs[index]));
        expect(gl.createProgram).toHaveBeenCalledTimes(3);
        expectPrograms(env, 3);
    });

    test("all failed programs recover together on explicit initialization", () => {
        const env = fixture(), { renderer, gl, state } = env;
        state.failLink = true;
        renderer.setSize(100, 80);
        const handles = resourceKeys.map(name => renderer[name]);
        expectPrograms(env, 0);
        state.failLink = false;
        renderer.initFrameBuffer();
        expectPrograms(env, 3);
        expectResources(env, handles);
        expect(gl.createProgram).toHaveBeenCalledTimes(6);
        expect(renderer.canRenderAmbientOcclusion()).toBe(true);
    });

    test("WebGL1 skips postprocessing resources and retries", () => {
        const env = fixture(1), { renderer, gl } = env;
        renderer.setSize(100, 80);
        renderer.initFrameBuffer();
        renderer.setSize(200, 120);
        for (const key of resourceKeys) expect(renderer[key]).toBeNull();
        expect(gl.createProgram).not.toHaveBeenCalled();
        expect(gl.createTexture).not.toHaveBeenCalled();
        expect(gl.createFramebuffer).not.toHaveBeenCalled();
        expect(gl.createBuffer).not.toHaveBeenCalled();
        expectPrograms(env, 0);
    });
});
