global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const textureKeys = ["_targetTexture", "_depthTexture", "_shadingTexture", "_scratchTexture"];

function environment() {
    const textures = new Set(), storage = new Map(), attachments = new Map(), clears = [];
    const state = { texture: null, framebuffer: null, clearDepth: 0, depthWrite: true };
    let constant = 100;
    const gl = new Proxy({
        drawingBufferWidth: 201, drawingBufferHeight: 161,
        isContextLost: () => false,
        createTexture: jest.fn(() => { const texture = {}; textures.add(texture); return texture; }),
        deleteTexture: jest.fn(texture => { expect(textures.delete(texture)).toBe(true); storage.delete(texture); }),
        bindTexture: jest.fn((target, texture) => { state.texture = texture; }),
        texImage2D: jest.fn((target, level, internalFormat, width, height, border, format, type) => {
            expect(textures.has(state.texture)).toBe(true);
            storage.set(state.texture, { width: Math.trunc(width), height: Math.trunc(height), internalFormat, format, type });
        }),
        createFramebuffer: jest.fn(() => ({})),
        bindFramebuffer: jest.fn((target, framebuffer) => { state.framebuffer = framebuffer; }),
        framebufferTexture2D: jest.fn((target, attachment, textarget, texture) => {
            if (!attachments.has(state.framebuffer)) attachments.set(state.framebuffer, new Map());
            attachments.get(state.framebuffer).set(attachment, texture);
        }),
        clearDepth: jest.fn(value => { state.clearDepth = value; }),
        depthMask: jest.fn(value => { state.depthWrite = value; }),
        clear: jest.fn(bits => clears.push({
            bits, value: state.clearDepth, depthWrite: state.depthWrite,
            depth: attachments.get(state.framebuffer)?.get(gl.DEPTH_ATTACHMENT)
        })),
        createBuffer: jest.fn(() => ({})),
        getAttribLocation: () => 0
    }, {
        get(target, key) {
            if (key in target) return target[key];
            if (/^[A-Z0-9_]+$/.test(key)) return target[key] = constant++;
            return target[key] = jest.fn();
        }
    });
    return { gl, textures, storage, attachments, clears, state };
}

function fixture(env = environment(), config = {}) {
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: { antialias: false, upscale: false, ...config }, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer;
    renderer._gl = renderer.context = env.gl;
    renderer._webglversion = 2;
    renderer.buildProgram = jest.fn(() => ({ uniforms: {} }));
    return { ...env, viewer, renderer };
}

function expectStorage(env, width, height) {
    const { renderer, gl, storage } = env;
    textureKeys.forEach((key, i) => expect(storage.get(renderer[key])).toEqual({
        width, height,
        internalFormat: i === 0 ? gl.RGBA : gl.DEPTH_COMPONENT32F,
        format: i === 0 ? gl.RGBA : gl.DEPTH_COMPONENT,
        type: i === 0 ? gl.UNSIGNED_BYTE : gl.FLOAT
    }));
}

function expectRouting(env) {
    const { renderer, gl, attachments, clears, state } = env;
    expect(state.framebuffer).toBe(renderer._fb);
    expect(attachments.get(renderer._fb).get(gl.COLOR_ATTACHMENT0)).toBe(renderer._targetTexture);
    expect(attachments.get(renderer._fb).get(gl.DEPTH_ATTACHMENT)).toBe(renderer._depthTexture);
    expect(clears[clears.length - 1]).toEqual({
        bits: gl.DEPTH_BUFFER_BIT, value: 1, depthWrite: true, depth: renderer._shadingTexture
    });
}

describe("framebuffer storage reuse", () => {
    test.each([false, true])("steady-size redraws retain storage and clear AO before drawing (AO: %s)", ao => {
        const env = fixture(), { renderer, gl, clears } = env;
        renderer._AOEnabled = ao;
        renderer.setSize(100, 80);
        renderer.setFrameBuffer();
        expectStorage(env, 100, 80);
        const textures = textureKeys.map(key => renderer[key]);
        expect(gl.texImage2D).toHaveBeenCalledTimes(4);
        const parameters = gl.texParameteri.mock.calls.length;
        for (let i = 0; i < 10; i++) {
            // Foreground/AO/volume passes can change depth writes and attachments.
            renderer.setDepthWrite(false);
            gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.TEXTURE_2D, renderer._scratchTexture, 0);
            renderer.setFrameBuffer();
            expectRouting(env);
        }
        expect(textureKeys.map(key => renderer[key])).toEqual(textures);
        expect(gl.createTexture).toHaveBeenCalledTimes(4);
        expect(gl.deleteTexture).not.toHaveBeenCalled();
        expect(gl.texImage2D).toHaveBeenCalledTimes(4);
        expect(gl.texParameteri).toHaveBeenCalledTimes(parameters);
        expect(clears).toHaveLength(11);
    });

    test("resizing reallocates storage once while preserving texture handles", () => {
        const env = fixture(), { renderer, gl } = env;
        renderer.setSize(100, 80);
        renderer.setFrameBuffer();
        const textures = textureKeys.map(key => renderer[key]);
        renderer.setSize(200, 120);
        renderer.setFrameBuffer();
        renderer.setFrameBuffer();
        expectStorage(env, 200, 120);
        expectRouting(env);
        expect(gl.texImage2D).toHaveBeenCalledTimes(8);
        expect(textureKeys.map(key => renderer[key])).toEqual(textures);
        expect(gl.createTexture).toHaveBeenCalledTimes(4);
        expect(gl.deleteTexture).not.toHaveBeenCalled();
    });

    test("a zero-sized viewport can become visible at a new size", () => {
        const env = fixture(), { renderer, gl } = env;
        renderer.setSize(0, 0);
        renderer.setFrameBuffer();
        const before = gl.texImage2D.mock.calls.length;
        renderer.setSize(160, 90);
        renderer.setFrameBuffer();
        expectStorage(env, 160, 90);
        expect(gl.texImage2D).toHaveBeenCalledTimes(before + 4);
        expectRouting(env);
    });

    test("grid viewport dimensions are cached independently of the full canvas", () => {
        const env = fixture(undefined, { rows: 2, cols: 2, row: 1, col: 1 });
        const { renderer, gl } = env;
        renderer.setSize(201, 161);
        renderer.setViewport();
        renderer.setFrameBuffer();
        renderer.setViewport();
        renderer.setFrameBuffer();
        expectStorage(env, 100, 80); // WebGL truncates half-pixel cell dimensions.
        expect(gl.texImage2D).toHaveBeenCalledTimes(4);
        gl.drawingBufferWidth = 240;
        gl.drawingBufferHeight = 180;
        renderer.setViewport();
        renderer.setFrameBuffer();
        expectStorage(env, 120, 90);
        expect(gl.texImage2D).toHaveBeenCalledTimes(8);
    });

    test("shader fallback does not mark storage allocated before recovery", () => {
        const env = fixture(), { renderer, gl, state } = env;
        renderer.setSize(100, 80);
        renderer._screenshader = null;
        renderer.setFrameBuffer();
        expect(state.framebuffer).toBeNull();
        expect(gl.texImage2D).not.toHaveBeenCalled();
        renderer.initFrameBuffer();
        renderer.setFrameBuffer();
        renderer.setFrameBuffer();
        expectStorage(env, 100, 80);
        expectRouting(env);
        expect(gl.texImage2D).toHaveBeenCalledTimes(4);
    });

    test("a replacement renderer allocates fresh storage at the same dimensions", () => {
        const env = environment();
        const first = fixture(env);
        first.renderer.setSize(100, 80);
        first.renderer.setFrameBuffer();
        const old = textureKeys.map(key => first.renderer[key]);
        first.renderer.dispose();
        expect(env.textures.size).toBe(0);
        const second = fixture(env);
        second.renderer.setSize(100, 80);
        second.renderer.setFrameBuffer();
        expectStorage(second, 100, 80);
        expect(env.gl.texImage2D).toHaveBeenCalledTimes(8);
        textureKeys.forEach(key => expect(old).not.toContain(second.renderer[key]));
        second.renderer.dispose();
        expect(env.textures.size).toBe(0);
    });

    test("shared offscreen sizing follows the actual canvas when viewers alternate", () => {
        const env = environment(), first = fixture(env), second = fixture(env);
        const size = { width: 0, height: 0 }, shared = {};
        const writes = { width: jest.fn(), height: jest.fn() };
        for (const key of Object.keys(size)) Object.defineProperty(shared, key, {
            get: () => size[key], set: value => { writes[key](value); size[key] = value; }
        });
        first.renderer._offscreen = second.renderer._offscreen = shared;
        first.renderer.setSize(100, 80);
        second.renderer.setSize(200, 120);
        first.renderer.setViewport();
        first.renderer.setFrameBuffer();
        for (let i = 0; i < 3; i++) first.renderer.setViewport();
        expect(writes.width.mock.calls).toEqual([[100]]);
        expect(writes.height.mock.calls).toEqual([[80]]);
        second.renderer.setViewport();
        second.renderer.setFrameBuffer();
        first.renderer.setViewport();
        first.renderer.setFrameBuffer();
        expect(writes.width.mock.calls).toEqual([[100], [200], [100]]);
        expect(writes.height.mock.calls).toEqual([[80], [120], [80]]);
        expectStorage(first, 100, 80);
        expectStorage(second, 200, 120);
        expect(env.gl.texImage2D).toHaveBeenCalledTimes(8);
    });
});
