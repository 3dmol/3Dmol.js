global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function fixture() {
    const viewer = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(viewer, { config: {}, WIDTH: 0, HEIGHT: 0 });
    viewer.setupRenderer();
    const plugin = viewer.renderer.sprites, textures = [];
    let constant = 100;
    const gl = new Proxy({}, { get(target, key) {
        if (!(key in target)) target[key] = /^[A-Z0-9_]+$/.test(key) ? constant++ : jest.fn();
        return target[key];
    } });
    const renderer = {
        devicePixelRatio: 2, setDepthTest: jest.fn(), setDepthWrite: jest.fn(),
        setTexture: jest.fn(texture => textures.push(texture.name))
    };
    plugin.gl = gl; plugin.renderer = renderer;
    plugin.sprite.program = {};
    plugin.sprite.attributes = { position: 0, uv: 1 };
    plugin.sprite.uniforms = new Proxy({}, { get: (target, key) => key });
    const camera = { projectionMatrix: new $3Dmol.Matrix4(), matrixWorldInverse: new $3Dmol.Matrix4() };
    return { plugin, gl, renderer, camera, textures,
        draw(sprites, inFront = false, width = 400, height = 200) {
            plugin.render({ __webglSprites: sprites }, camera, width, height, inFront);
        }
    };
}
function sprite(id, z = 0, material = {}) {
    const world = new $3Dmol.Matrix4(); world.elements[14] = z;
    return {
        id, visible: true, position: { x: 20, y: 30, z }, scale: { x: 1, y: 1 }, rotation: 0,
        matrixWorld: world, _modelViewMatrix: new $3Dmol.Matrix4(),
        material: { depthTest: true, depthWrite: true, opacity: 1, useScreenCoordinates: true,
            color: { r: 1, g: 1, b: 1 }, map: { name: id, image: { width: 64, height: 32 } }, ...material }
    };
}
afterEach(() => jest.restoreAllMocks());

test("both passes reuse their list without retaining labels or sorting the scene array", () => {
    const env = fixture(), a = sprite(1, -1), b = sprite(2, -5), front = sprite(3, 0, { depthTest: false });
    const source = [a, front, b], lists = [];
    const original = env.plugin.drawSprites;
    jest.spyOn(env.plugin, "drawSprites").mockImplementation(function (list, ...args) {
        lists.push(list);
        return original.call(this, list, ...args);
    });
    const scratch = env.plugin.renderList;
    env.draw(source); env.draw(source, true); env.draw([]);
    expect(env.textures).toEqual([2, 1, 3]);
    expect(source).toEqual([a, front, b]);
    expect(lists).toHaveLength(2);
    for (const list of lists) expect(list).toBe(scratch);
    expect(scratch).toHaveLength(0);
    expect(env.plugin.rendering).toBe(false);
});

test("sparse lists and legacy depth-test matching select the same foreground and embedded labels", () => {
    const env = fixture(), source = [];
    source[1] = sprite(1, 0, { depthTest: true });
    source[3] = sprite(2, 0, { depthTest: 0 });
    source[5] = sprite(3, 0, { depthTest: undefined });
    env.draw(source); env.draw(source, true);
    expect(env.textures).toEqual([1, 2]);
    expect(source.length).toBe(6);
    expect(0 in source).toBe(false);
    expect(env.plugin.renderList).toHaveLength(0);
});

test("hidden labels retain their sorting participation and opacity-zero labels are not drawn", () => {
    const env = fixture(), near = sprite(1, -1), hidden = sprite(2), far = sprite(3, -9);
    hidden.visible = false;
    const invisible = sprite(4, 0, { opacity: 0 });
    env.draw([near, hidden, far, invisible]);
    // An undefined hidden-label depth separates these entries in the existing stable sort.
    expect(env.textures).toEqual([1, 3]);
    expect(hidden.z).toBeUndefined();
    expect(invisible.z).toBeUndefined();
    env.textures.length = 0;
    env.draw([near, far]);
    expect(env.textures).toEqual([3, 1]);
});

test.each([false, true])("scale uniforms preserve viewport, device-pixel ratio and sprite scaling (screen: %s)", screen => {
    const env = fixture(), label = sprite(7, -3, { useScreenCoordinates: screen, scaleByViewport: true });
    label.scale = { x: 3, y: 4 };
    env.draw([label]);
    expect(env.gl.uniform2f).toHaveBeenCalledWith("scale", (64 * 2 / 400) * (3 / 200), (32 * 2 / 200) * (4 / 200));
    expect(env.gl.uniform2fv).not.toHaveBeenCalled();
    if (screen) expect(env.gl.uniform3f).toHaveBeenCalledWith("screenPosition", -0.8, 0.4, 0);
    else expect(env.gl.uniformMatrix4fv).toHaveBeenCalledWith("modelViewMatrix", false, label._modelViewMatrix.elements);
    env.gl.uniform2f.mockClear();
    env.renderer.devicePixelRatio = 1;
    label.material.scaleByViewport = false;
    env.draw([label], false, 200, 100);
    expect(env.gl.uniform2f).toHaveBeenCalledWith("scale", 64 / 200 * 3, 32 / 100 * 4);
});

test("zero-size viewport keeps the previous nonfinite scale values", () => {
    const env = fixture(), label = sprite(1);
    label.scale.x = 0;
    env.draw([label], false, 0, 0);
    const scale = env.gl.uniform2f.mock.calls.find(call => call[0] === "scale");
    expect(Number.isNaN(scale[1])).toBe(true);
    expect(scale[2]).toBe(Infinity);
    expect(env.plugin.renderList).toHaveLength(0);
});

test.each([false, true])("higher-resolution textures preserve sprite size and screen offsets (screen: %s)", screen => {
    const env = fixture(), label = sprite(7, -3, {
        useScreenCoordinates: screen, alignment: { x: -1, y: 1 }, screenOffset: { x: 12, y: -8 },
        map: { name: 7, image: { width: 192, height: 96 }, logicalWidth: 64, logicalHeight: 32 }
    });
    label.scale = { x: 3, y: 4 };
    env.draw([label]);
    expect(env.gl.uniform2f).toHaveBeenCalledWith("scale", (64 * 2 / 400) * 3, (32 * 2 / 200) * 4);
    expect(env.gl.uniform2f).toHaveBeenCalledWith("alignment", -0.625, 0.5);
    if (screen) expect(env.gl.uniform3f).toHaveBeenCalledWith("screenPosition", -0.8, 0.4, 0);
    else expect(env.gl.uniformMatrix4fv).toHaveBeenCalledWith("modelViewMatrix", false, label._modelViewMatrix.elements);
});

test("empty and nonmatching passes do no GL work", () => {
    const env = fixture();
    env.draw([]); env.draw([sprite(1, 0, { depthTest: false })]);
    env.plugin.render(null, env.camera, 400, 200);
    expect(env.gl.useProgram).not.toHaveBeenCalled();
    expect(env.gl.drawElements).not.toHaveBeenCalled();
    expect(env.renderer.setTexture).not.toHaveBeenCalled();
    expect(env.plugin.renderList).toHaveLength(0);
});

test("errors release label references and leave the scratch list reusable", () => {
    const env = fixture(), label = sprite(1);
    env.renderer.setTexture.mockImplementationOnce(() => { throw new Error("texture upload failed"); });
    expect(() => env.draw([label])).toThrow("texture upload failed");
    expect(env.plugin.renderList).toHaveLength(0);
    expect(env.plugin.rendering).toBe(false);
    env.draw([sprite(2)]);
    expect(env.textures).toEqual([2]);
    expect(env.plugin.renderList).toHaveLength(0);
});

test("nested texture callbacks preserve the outer pass list", () => {
    const env = fixture(), outer = [sprite(1, -4), sprite(2, -1)], nested = sprite(3);
    let entered = false;
    env.renderer.setTexture.mockImplementation(texture => {
        env.textures.push(texture.name);
        if (!entered) { entered = true; env.draw([nested]); }
    });
    env.draw(outer);
    expect(env.textures).toEqual([1, 3, 2]);
    expect(env.gl.drawElements).toHaveBeenCalledTimes(3);
    expect(env.plugin.renderList).toHaveLength(0);
    expect(env.plugin.rendering).toBe(false);
});

test("dispose clears scratch references and prevents subsequent passes", () => {
    const env = fixture();
    env.plugin.renderList.push(sprite(1));
    env.plugin.dispose();
    expect(env.plugin.renderList).toHaveLength(0);
    expect(() => env.draw([sprite(2)])).not.toThrow();
    expect(env.gl.drawElements).not.toHaveBeenCalled();
});
