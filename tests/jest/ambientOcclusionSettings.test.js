global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function viewer(ambientOcclusion) {
    const v = Object.create($3Dmol.GLViewer.prototype);
    Object.assign(v, { config: { ambientOcclusion }, WIDTH: 0, HEIGHT: 0 });
    // Use the real renderer constructor without creating a WebGL context.
    v.setupRenderer();
    return v;
}

const routes = {
    setConfig: (v, settings) => v.setConfig({ ambientOcclusion: settings }),
    setViewStyle: (v, settings) => v.setViewStyle({ style: "ambientOcclusion", ...settings })
};
const state = v => [v.renderer._AOEnabled, v.renderer._AOstrength, v.renderer._AOradius];

describe("ambient occlusion settings", () => {
    test.each(Object.keys(routes))("%s matches constructor settings", route => {
        const settings = { strength: 0.6, radius: 12 };
        const v = viewer();
        routes[route](v, settings);
        expect(state(v)).toEqual(state(viewer(settings)));
        expect(state(v)).toEqual([true, 0.6, 12]);
    });

    test.each(Object.keys(routes))("%s accepts explicit zeros and can re-enable shading", route => {
        const v = viewer({ strength: 1.5, radius: 9 });
        routes[route](v, { strength: 0, radius: 0 });
        expect(state(v)).toEqual(state(viewer({ strength: 0, radius: 0 })));
        expect(state(v)).toEqual([false, 0, 0]);
        routes[route](v, { radius: 11 });
        expect(state(v)).toEqual([false, 0, 11]);
        routes[route](v, { strength: 0.8 });
        expect(state(v)).toEqual([true, 0.8, 11]);
    });

    test.each(Object.keys(routes))("%s retains omitted settings", route => {
        const v = viewer({ strength: 0.7, radius: 8 });
        routes[route](v, {});
        expect(state(v)).toEqual([true, 0.7, 8]);
        routes[route](v, { strength: undefined, radius: undefined });
        expect(state(v)).toEqual([true, 0.7, 8]);
    });

    test("direct renderer calls retain the legacy scale alias", () => {
        const v = viewer();
        v.renderer.enableAmbientOcclusion({ scale: 17 });
        expect(state(v)).toEqual([true, 1, 17]);
        v.renderer.disableAmbientOcclusion();
        v.renderer.enableAmbientOcclusion();
        expect(state(v)).toEqual([true, 1, 17]);
        v.setConfig({ ambientOcclusion: { scale: 19 } });
        expect(state(v)).toEqual([true, 1, 19]);
    });

    test("radius takes priority over the legacy scale alias, including zero", () => {
        const v = viewer();
        v.renderer.enableAmbientOcclusion({ radius: 13, scale: 17 });
        expect(v.renderer._AOradius).toBe(13);
        v.renderer.enableAmbientOcclusion({ radius: 0, scale: 19 });
        expect(v.renderer._AOradius).toBe(0);
        v.renderer.enableAmbientOcclusion({ radius: undefined, scale: 21 });
        expect(v.renderer._AOradius).toBe(21);
    });

    test.each(Object.keys(routes))("%s supplies the requested values to the shading uniforms", route => {
        const v = viewer(), renderer = v.renderer;
        routes[route](v, { strength: 0.4, radius: 14 });
        const gl = new Proxy({}, {
            get(target, key) {
                if (!(key in target)) target[key] = /^[A-Z0-9_]+$/.test(key) ? 1 : jest.fn();
                return target[key];
            }
        });
        Object.assign(renderer, {
            _gl: gl, _fb: {}, _screenshader: {}, _blurshader: {},
            _AOshader: { uniforms: { total_strength: "strength", radius: "radius" } },
            renderObjects: jest.fn()
        });
        const scene = { __lights: [], fog: null, __webglObjects: [{
            render: true, opaque: {}, opaqueDepth: {},
            object: { _modelViewMatrix: new $3Dmol.Matrix4() }
        }] };
        renderer.setShading(scene, { projectionMatrix: new $3Dmol.Matrix4() }, "opaque");
        expect(renderer.renderObjects).toHaveBeenCalledTimes(1);
        expect(gl.uniform1f).toHaveBeenCalledWith("strength", 0.4);
        expect(gl.uniform1f).toHaveBeenCalledWith("radius", 14);
    });
});
