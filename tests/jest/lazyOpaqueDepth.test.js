global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const flush = async () => { await Promise.resolve(); await Promise.resolve(); };
let environments;

function camera() {
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

function record(properties = {}, { volume = false, alpha = false, visible = true } = {}) {
    let shape;
    if (volume) {
        const data = new $3Dmol.VolumeData("", "cube");
        Object.assign(data, {
            size: { x: 2, y: 2, z: 2 }, unit: { x: 1, y: 1, z: 1 },
            origin: { x: 0, y: 0, z: 0 }, data: new Float32Array(8).fill(1)
        });
        shape = new $3Dmol.GLVolumetricRender(data, { transferfn: [
            { value: 0, color: "black", opacity: 0 }, { value: 1, color: "white", opacity: 1 }
        ] });
    } else {
        shape = new $3Dmol.GLShape({});
        shape.addCustom({
            vertexArr: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }],
            faceArr: [0, 1, 2]
        });
    }
    let group;
    shape.globj({ add(object) { group = object; }, remove() {} });
    const object = group.children[0];
    Object.assign(object.material, properties);
    object.visible = visible;
    object.geometry.alpha = alpha;
    object._modelViewMatrix = new $3Dmol.Matrix4();
    object._normalMatrix = new $3Dmol.Matrix3();
    return { object, buffer: {}, render: false, clone: jest.spyOn(object.material, "clone") };
}

function environment(records, { ao = false, failProgram = 0 } = {}) {
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: {}, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer;
    const state = { failProgram, framebuffer: null, depth: null };
    let constant = 100, programCount = 0;
    const gl = new Proxy({
        drawingBufferWidth: 640, drawingBufferHeight: 480, isContextLost: () => false,
        createProgram: () => ({ number: ++programCount }),
        createShader: () => ({}),
        getShaderParameter: () => true,
        linkProgram: program => { program.linked = program.number !== state.failProgram; },
        getProgramParameter: program => program.linked,
        getProgramInfoLog: () => "controlled link failure",
        getUniformLocation: () => null, getAttribLocation: () => 0,
        createBuffer: () => ({}), createTexture: () => ({}), createFramebuffer: () => ({}),
        useProgram: program => expect(program.linked).toBe(true),
        bindFramebuffer: (target, framebuffer) => { state.framebuffer = framebuffer; },
        framebufferTexture2D: (target, attachment, textureTarget, texture) => {
            if (attachment === gl.DEPTH_ATTACHMENT) state.depth = texture;
        }
    }, {
        get(target, key) {
            if (key in target) return target[key];
            if (/^[A-Z0-9_]+$/.test(key)) return target[key] = constant++;
            return target[key] = jest.fn();
        }
    });
    renderer._gl = gl;
    renderer._webglversion = 2;
    renderer.setSize(100, 80);
    renderer._AOEnabled = ao;
    renderer.autoUpdateScene = false;
    renderer.autoUpdateObjects = false;
    renderer.renderSprites = jest.fn();
    const scene = {
        __lights: [], __webglObjects: records, __webglSprites: [],
        fog: { color: { r: 1, g: 1, b: 1 }, near: 1, far: 100 }
    };
    const viewCamera = camera(), draws = [];
    // Keep actual render ordering, AO setup, and framebuffer routing; record the final draws.
    renderer.renderBuffer = (camera, lights, fog, material, buffer, object) => {
        draws.push({ material, object, depth: state.depth, framebuffer: state.framebuffer });
    };
    const env = {
        renderer, state, scene, records: [...records], draws, camera: viewCamera,
        draw() { draws.length = 0; renderer.render(scene, viewCamera); }
    };
    environments.push(env);
    return env;
}

const drawsFor = (env, record) => env.draws.filter(draw => draw.object === record.object);

function expectMainAndVolume(env, opaque, volume, shaded = false) {
    const main = drawsFor(env, opaque).filter(draw => draw.material.opacity !== 0);
    expect(main).toHaveLength(1);
    expect(main[0].material).toBe(shaded ? opaque.opaqueShaded : opaque.object.material);
    expect(main[0].depth).toBe(env.renderer._depthTexture);
    const volumeDraws = drawsFor(env, volume);
    expect(volumeDraws).toHaveLength(1);
    expect(volumeDraws[0].material).toBe(volume.object.material);
    expect(volumeDraws[0].depth).toBeNull();
    expect(volumeDraws[0].framebuffer).toBe(env.renderer._fb);
}

describe("opaque depth material allocation", () => {
    beforeEach(() => {
        environments = [];
        jest.spyOn(console, "error").mockImplementation(() => {});
    });
    afterEach(async () => {
        try {
            for (const env of environments) {
                for (const entry of env.records) {
                    env.renderer.removeObject(entry.object, env.scene);
                    entry.object.material.dispose();
                }
            }
            await flush();
        } finally {
            jest.restoreAllMocks();
        }
    });

    test("ordinary opaque and volume frames allocate no unused depth materials", async () => {
        const opaque = record(), volume = record({}, { volume: true }), env = environment([opaque, volume]);
        for (let i = 0; i < 20; i++) {
            env.draw();
            expectMainAndVolume(env, opaque, volume);
            await flush();
        }
        expect(opaque.clone).not.toHaveBeenCalled();
        expect(opaque.opaqueDepth).toBeNull();
        expect(volume.clone).not.toHaveBeenCalled();
    });

    test.each([true, false])("AO includes every opaque occluder (global AO: %s)", globalAO => {
        const selected = record({ hasAO: true }), other = record(), volume = record({}, { volume: true });
        const env = environment([selected, other, volume], { ao: globalAO });
        env.draw();
        for (const entry of [selected, other]) {
            const depth = drawsFor(env, entry).filter(draw => draw.material.opacity === 0);
            expect(depth).toHaveLength(1);
            expect(depth[0].material).toBe(entry.opaqueDepth);
            expect(depth[0].depth).toBe(env.renderer._shadingTexture);
            expect(entry.opaqueDepth).not.toBe(entry.object.material);
        }
        expect(selected.clone).toHaveBeenCalledTimes(2);
        expect(other.clone).toHaveBeenCalledTimes(globalAO ? 2 : 1);
        expectMainAndVolume(env, selected, volume, true);
        expectMainAndVolume(env, other, volume, globalAO);
        // A second AO pass in the same frame uses the existing depth materials.
        env.renderer.setShading(env.scene, env.camera, "opaque");
        expect(selected.clone).toHaveBeenCalledTimes(2);
        expect(other.clone).toHaveBeenCalledTimes(globalAO ? 2 : 1);
    });

    test("failed AO leaves depth unallocated until shader recovery", async () => {
        const opaque = record(), volume = record({}, { volume: true });
        const env = environment([opaque, volume], { ao: true, failProgram: 2 });
        env.draw();
        expect(opaque.opaqueDepth).toBeNull();
        expect(opaque.clone).toHaveBeenCalledTimes(1); // Existing shaded pass only.
        expectMainAndVolume(env, opaque, volume);
        await flush();
        env.state.failProgram = 0;
        env.renderer.setSize(100, 80);
        env.draw();
        expect(opaque.opaqueDepth.opacity).toBe(0);
        expectMainAndVolume(env, opaque, volume, true);
    });

    test("transparent priming and per-vertex alpha retain their depth routes", () => {
        const opaque = record(), transparent = record({ transparent: true, opacity: 0.4 });
        const alpha = record({ transparent: true, opacity: 0.4 }, { alpha: true });
        const volume = record({}, { volume: true }), env = environment([opaque, transparent, alpha, volume]);
        env.draw();
        expect(opaque.clone).not.toHaveBeenCalled();
        expect(transparent.clone).toHaveBeenCalledTimes(1);
        expect(alpha.clone).not.toHaveBeenCalled();
        const transparentDraws = drawsFor(env, transparent);
        expect(transparentDraws.map(draw => draw.material)).toEqual([transparent.transparentDepth, transparent.object.material]);
        expect(transparent.transparentDepth.opacity).toBe(0);
        for (const draw of transparentDraws) expect(draw.depth).toBe(env.renderer._depthTexture);
        expect(drawsFor(env, alpha).map(draw => draw.material)).toEqual([alpha.object.material]);
        expectMainAndVolume(env, opaque, volume);
    });

    test("AO toggles retire passes and subsequent clones copy current settings", async () => {
        const opaque = record(), volume = record({}, { volume: true });
        const env = environment([opaque, volume], { ao: true });
        env.draw();
        const firstDepth = opaque.opaqueDepth, firstShaded = opaque.opaqueShaded;
        const depthDisposed = jest.fn(), shadedDisposed = jest.fn();
        firstDepth.addEventListener("dispose", depthDisposed);
        firstShaded.addEventListener("dispose", shadedDisposed);
        env.renderer._AOEnabled = false;
        env.draw();
        expect(opaque.opaqueDepth).toBeNull();
        expect(opaque.opaqueShaded).toBeNull();
        await flush();
        expect(depthDisposed).toHaveBeenCalledTimes(1);
        expect(shadedDisposed).toHaveBeenCalledTimes(1);

        opaque.object.material.color.setHex(0x123456);
        opaque.object.material.depthWrite = false;
        env.renderer._AOEnabled = true;
        env.draw();
        expect(opaque.opaqueDepth).not.toBe(firstDepth);
        expect(opaque.opaqueDepth.color.getHex()).toBe(0x123456);
        expect(opaque.opaqueDepth.depthWrite).toBe(false);
        expect(opaque.opaqueDepth.opacity).toBe(0);
        expectMainAndVolume(env, opaque, volume, true);

        const removed = jest.fn();
        opaque.opaqueDepth.addEventListener("dispose", removed);
        env.renderer.removeObject(opaque.object, env.scene);
        await flush();
        expect(removed).toHaveBeenCalledTimes(1);
        env.draw();
        expect(drawsFor(env, opaque)).toHaveLength(0);
        expect(drawsFor(env, volume)).toHaveLength(1);
        expect(depthDisposed).toHaveBeenCalledTimes(1);
    });

    test("hidden and wireframe objects do not acquire AO depth materials", () => {
        const hidden = record({}, { visible: false }), wireframe = record({ wireframe: true });
        const solid = record(), env = environment([hidden, wireframe, solid], { ao: true });
        env.draw();
        expect(hidden.clone).not.toHaveBeenCalled();
        expect(hidden.opaqueDepth).toBeUndefined();
        expect(wireframe.opaqueDepth).toBeNull();
        expect(wireframe.clone).toHaveBeenCalledTimes(1); // Existing shaded pass only.
        expect(solid.opaqueDepth.opacity).toBe(0);
    });
});
