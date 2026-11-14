global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");
const raw = values => new $3Dmol.Vector3().set(...values);
const values = v => [v.x, v.y, v.z];
const exact = (actual, expected) => values(actual).forEach((value, i) => {
    expect(Object.is(value, values(expected)[i])).toBe(true);
});

// Previous implementations serve as arithmetic references, including clone defaults.
function oldApplyQuaternion(v, q) {
    const { x, y, z } = v, qx = q.x, qy = q.y, qz = q.z, qw = q.w;
    const t = { x: 2 * (y * qz - z * qy), y: 2 * (z * qx - x * qz), z: 2 * (x * qy - y * qx) };
    const t2 = { x: t.y * qz - t.z * qy, y: t.z * qx - t.x * qz, z: t.x * qy - t.y * qx };
    v.x = x + qw * t.x + t2.x;
    v.y = y + qw * t.y + t2.y;
    v.z = z + qw * t.z + t2.z;
    return v;
}
function oldRotateAboutVector(v, axis, angle) {
    axis.normalize();
    const c = Math.cos(angle), s = Math.sin(angle);
    const t1 = v.clone().multiplyScalar(c);
    const t2 = axis.clone().cross(v).multiplyScalar(s);
    const t3 = axis.clone().multiplyScalar(axis.clone().dot(v)).multiplyScalar(1 - c);
    const result = t1.add(t2).add(t3);
    return v.set(result.x, result.y, result.z);
}
function cases() {
    let seed = 0x3d012345;
    const random = () => { seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0; return seed / 4294967296 * 200 - 100; };
    const result = Array.from({ length: 250 }, () => ({
        vector: [random(), random(), random()], axis: [random(), random(), random()], angle: random(), w: random()
    }));
    for (const vector of [[0, -0, 0], [1, 2, 3], [NaN, 1, -0], [Infinity, -Infinity, 1]]) {
        for (const axis of [[0, 0, 0], [1, 0, -0], [NaN, 2, 1], [Infinity, 1, 0]]) {
            for (const angle of [0, -0, Math.PI / 3, NaN, Infinity]) result.push({ vector, axis, angle, w: angle });
        }
    }
    return result;
}
afterEach(() => jest.restoreAllMocks());

test("quaternion application matches previous arithmetic exactly for finite and exceptional values", () => {
    for (const input of cases()) {
        const actual = raw(input.vector), expected = raw(input.vector);
        const q = { x: input.axis[0], y: input.axis[1], z: input.axis[2], w: input.w };
        const before = { ...q };
        oldApplyQuaternion(expected, q);
        expect(actual.applyQuaternion(q)).toBe(actual);
        exact(actual, expected);
        expect(q).toEqual(before);
    }
});

test("quaternion application captures every component before writing an aliased input", () => {
    const actual = raw([0.5, -0.25, 0.75]), expected = raw(values(actual));
    actual.w = expected.w = 0.625;
    oldApplyQuaternion(expected, expected);
    expect(actual.applyQuaternion(actual)).toBe(actual);
    exact(actual, expected);
    expect(actual.w).toBe(expected.w);
});

test.each([false, true])("axis rotation preserves exact results and normalization side effects (alias: %s)", alias => {
    for (const input of cases()) {
        const actual = raw(input.vector), expected = raw(input.vector);
        const actualAxis = alias ? actual : raw(input.axis);
        const expectedAxis = alias ? expected : raw(input.axis);
        oldRotateAboutVector(expected, expectedAxis, input.angle);
        expect(actual.rotateAboutVector(actualAxis, input.angle)).toBe(actual);
        exact(actual, expected);
        exact(actualAxis, expectedAxis);
    }
});

test("repeated rotations normalize the supplied axis while allocating no cloned vectors", () => {
    const vector = raw([2, 3, 4]), axis = raw([0, 0, 8]);
    const clone = jest.spyOn($3Dmol.Vector3.prototype, "clone");
    for (let i = 0; i < 1000; i++) vector.rotateAboutVector(axis, 0.001);
    expect(clone).not.toHaveBeenCalled();
    expect(values(axis)).toEqual([0, 0, 1]);
    expect(vector.length()).toBeCloseTo(Math.sqrt(29), 10);
});
