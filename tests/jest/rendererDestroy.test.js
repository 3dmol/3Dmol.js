global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const outlineKeys = ["_outlineMaterial", "_outlineSphereImposterMaterial", "_outlineStickImposterMaterial"];
const flush = async () => { await Promise.resolve(); await Promise.resolve(); };

function fakeGL() {
    const live = { Buffer: new Set(), Texture: new Set(), Framebuffer: new Set(), Program: new Set(), Shader: new Set() };
    const state = { lost: false, compile: true };
    let constant = 100;
    const methods = {
        drawingBufferWidth: 320, drawingBufferHeight: 240,
        isContextLost: () => state.lost,
        getShaderParameter: () => state.compile,
        getShaderInfoLog: () => "controlled compile failure",
        getProgramParameter: () => true,
        getUniformLocation: () => null,
        getAttribLocation: () => 0,
        getExtension: jest.fn(() => ({ loseContext: jest.fn() }))
    };
    for (const kind of Object.keys(live)) {
        methods["create" + kind] = jest.fn(() => {
            const handle = { kind };
            live[kind].add(handle);
            return handle;
        });
        methods["delete" + kind] = jest.fn(handle => expect(live[kind].delete(handle)).toBe(true));
    }
    const gl = new Proxy(methods, {
        get(target, key) {
            if (key in target) return target[key];
            if (/^[A-Z0-9_]+$/.test(key)) return target[key] = constant++;
            return target[key] = jest.fn();
        }
    });
    return { gl, live, state };
}

function renderer(env, canvas = document.createElement("canvas")) {
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: { canvas }, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const renderer = viewer.renderer;
    if (env) {
        renderer._gl = renderer.context = env.gl;
        renderer._webglversion = 2;
    }
    return renderer;
}

function populate(renderer) {
    renderer.setSize(100, 80);
    renderer.sprites.init(renderer);
    for (const key of outlineKeys) renderer.initMaterial(renderer[key], [], null, {});
    const parent = { add: jest.fn(), remove: jest.fn() };
    const shape = new $3Dmol.GLShape({ opacity: 0.5 });
    shape.addSphere({ center: { x: 1, y: 2, z: 3 }, radius: 1 });
    shape.globj(parent);
    const mesh = parent.add.mock.calls[0][0].children[0];
    const scene = { __webglObjects: [], __webglSprites: [] };
    renderer.addObject(mesh, scene);
    renderer.initMaterial(mesh.material, [], null, mesh);
    renderer.unrollBufferMaterial(scene.__webglObjects[0]);
    const depth = scene.__webglObjects[0].transparentDepth;
    if (depth) renderer.initMaterial(depth, [], null, mesh);
    const label = new $3Dmol.Label("dispose", {});
    label.setContext();
    renderer.setTexture(label.sprite.material.map, 0, false);
    return { shape, mesh, parent, scene, label };
}

function expectEmpty(env) {
    for (const handles of Object.values(env.live)) expect(handles.size).toBe(0);
}

describe("renderer destruction", () => {
    test("releases geometry, labels, sprite, outline and postprocessing allocations exactly once", async () => {
        const env = fakeGL(), r = renderer(env), objects = populate(r);
        expect(env.live.Buffer.size).toBeGreaterThan(2);
        expect(env.live.Texture.size).toBe(5);
        expect(env.live.Program.size).toBeGreaterThan(4);
        r.dispose();
        expectEmpty(env);
        expect(r.info.memory).toEqual({ programs: 0, geometries: 0, textures: 0 });
        expect(r.getContext()).toBeNull();
        expect(r.isLost()).toBe(true);
        for (const resource of [objects.mesh.geometry, objects.mesh.material, objects.label.sprite.material.map]) {
            expect(resource.listeners.dispose || []).toHaveLength(0);
        }
        r.dispose();
        objects.shape.removegl(objects.parent);
        objects.label.dispose();
        await flush();
        expectEmpty(env);
    });

    test("already released owners and deferred material disposal cannot double-delete resources", async () => {
        const env = fakeGL(), r = renderer(env), objects = populate(r);
        const callback = jest.fn();
        objects.mesh.material.addEventListener("dispose", callback);
        objects.shape.removegl(objects.parent);
        objects.label.dispose();
        r.removeObject(objects.mesh, objects.scene);
        r.enableOutline({ width: 0.3 });
        r.dispose();
        expectEmpty(env);
        await flush();
        expectEmpty(env);
        expect(callback).toHaveBeenCalledTimes(1);
    });

    test("disposing one renderer preserves peer allocations and its shared canvas/context", async () => {
        const env = fakeGL(), canvas = document.createElement("canvas");
        document.body.appendChild(canvas);
        const first = renderer(env, canvas), second = renderer(env, canvas);
        populate(first);
        const peer = populate(second);
        const peerPrograms = second._programs.map(entry => entry.program);
        const peerTexture = peer.label.sprite.material.map.__webglTexture;
        const width = canvas.width, height = canvas.height;
        first.dispose();
        await flush();
        expect(canvas.parentNode).toBe(document.body);
        expect([canvas.width, canvas.height]).toEqual([width, height]);
        expect(second.getContext()).toBe(env.gl);
        expect(second.isLost()).toBe(false);
        peerPrograms.forEach(program => expect(env.live.Program.has(program)).toBe(true));
        expect(env.live.Texture.has(peerTexture)).toBe(true);
        expect(env.gl.getExtension).not.toHaveBeenCalledWith("WEBGL_lose_context");
        second.setSize(200, 120);
        second.setTexture(peer.label.sprite.material.map, 0, false);
        second.dispose();
        expectEmpty(env);
        canvas.remove();
    });

    test("late rendering and allocation entry points do not recreate disposed resources", () => {
        const env = fakeGL(), r = renderer(env), objects = populate(r);
        const canvas = r.getCanvas(), width = canvas.width, height = canvas.height;
        r.dispose();
        const creates = Object.keys(env.live).map(kind => env.gl["create" + kind].mock.calls.length);
        r.render(null, null);
        r.setSize(400, 300);
        r.setViewport();
        r.initFrameBuffer();
        r.setFrameBuffer();
        r.renderFrameBuffertoScreen();
        r.initMaterial(objects.mesh.material, [], null, objects.mesh);
        r.setTexture(objects.label.sprite.material.map, 0, false);
        r.addObject(objects.mesh, objects.scene);
        r.sprites.init(r);
        r.sprites.render(null, null, 100, 80);
        expect(Object.keys(env.live).map(kind => env.gl["create" + kind].mock.calls.length)).toEqual(creates);
        expect([canvas.width, canvas.height]).toEqual([width, height]);
        expectEmpty(env);
    });

    test("zero-size initialization and lost contexts can be disposed repeatedly", () => {
        const uninitialized = renderer();
        expect(() => { uninitialized.dispose(); uninitialized.dispose(); }).not.toThrow();
        const env = fakeGL(), r = renderer(env);
        populate(r);
        env.state.lost = true;
        r.dispose();
        r.dispose();
        expectEmpty(env);
    });

    test("sprite initialization releases shader objects after linking", () => {
        const env = fakeGL(), r = renderer(env);
        r.sprites.init(r);
        expect(env.live.Shader.size).toBe(0);
        expect(env.live.Program.size).toBe(1);
        expect(env.live.Buffer.size).toBe(2);
        r.dispose();
        expectEmpty(env);
    });

    test("failed sprite compilation leaves no program, shader or buffer allocations", () => {
        const env = fakeGL(), r = renderer(env);
        env.state.compile = false;
        expect(() => r.sprites.init(r)).toThrow("Error compiling shader");
        expectEmpty(env);
        r.dispose();
        expectEmpty(env);
    });
});
