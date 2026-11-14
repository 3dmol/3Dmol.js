global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const Geometry = new $3Dmol.GLShape({}).geo.constructor;
const methods = ["drawAtomSphere", "drawAtomClickSphere", "drawAtomInstanced", "drawAtomImposter"];

function atom(x, y, z) {
    return {
        x, y, z, elem: "C", color: 0xff0000, clickable: true,
        style: { sphere: { radius: 1.25 }, clicksphere: { radius: 1.25 } },
        intersectionShape: { sphere: [], cylinder: [], line: [], triangle: [] }
    };
}
const xyz = v => [v.x, v.y, v.z];

test.each(methods)("%s creates independent picking centers with the original coordinate defaults", method => {
    const model = new $3Dmol.GLModel(0);
    for (const coords of [[1, -2, 3], [NaN, -0, Infinity], [undefined, null, 0]]) {
        const a = atom(...coords), geo = new Geometry(true, true, true);
        const expected = new $3Dmol.Sphere(new $3Dmol.Vector3(...coords), 1.25);
        model[method](a, geo);
        expect(a.intersectionShape.sphere).toHaveLength(1);
        const sphere = a.intersectionShape.sphere[0];
        expect(sphere.center).not.toBe(a);
        expect(sphere.center).toBeInstanceOf($3Dmol.Vector3);
        expect(xyz(sphere.center)).toEqual(xyz(expected.center));
        expect(sphere.radius).toBe(expected.radius);
        // Geometry owns a snapshot, including across later atom updates and rebuilds.
        a.x = 8; a.y = 9; a.z = 10;
        model[method](a, geo);
        const next = a.intersectionShape.sphere[1];
        expect(next.center).not.toBe(sphere.center);
        expect(xyz(next.center)).toEqual([8, 9, 10]);
        expect(xyz(sphere.center)).toEqual(xyz(expected.center));
        sphere.center.set(11, 12, 13);
        expect(xyz(a)).toEqual([8, 9, 10]);
        expect(xyz(next.center)).toEqual([8, 9, 10]);
    }
});

test.each(methods)("%s preserves hover-only picking and excludes noninteractive or hidden atoms", method => {
    const model = new $3Dmol.GLModel(0);
    const a = atom(2, 3, 4), geo = new Geometry(true, true, true);
    a.clickable = false;
    model[method](a, geo);
    expect(a.intersectionShape.sphere).toEqual([]);
    a.hoverable = true;
    model[method](a, geo);
    expect(a.intersectionShape.sphere).toHaveLength(1);
    a.style.sphere.hidden = a.style.clicksphere.hidden = true;
    model[method](a, geo);
    expect(a.intersectionShape.sphere).toHaveLength(1);
});
