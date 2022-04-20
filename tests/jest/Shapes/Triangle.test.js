/**
*@jest-environment jsdom
*/
const $3Dmol = require("../../../build/3Dmol.js");

let range = (start, end) => [...Array.from(new Array(end + 1).keys()).slice(start)];

describe('Triangle Tests', function triangleTests() {
  test("Triangle constructor: a is all 1's, b is all 2's, and c is all 3's", () => {
    let a = new $3Dmol.Vector3(1, 1, 1);
    let b = new $3Dmol.Vector3(2, 2, 2);
    let c = new $3Dmol.Vector3(3, 3, 3);
    let tri = new $3Dmol.Triangle(a, b, c);
    expect(tri.a).toEqual({x: 1, y: 1, z: 1});
    expect(tri.b).toEqual({x: 2, y: 2, z: 2});
    expect(tri.c).toEqual({x: 3, y: 3, z: 3});
  });
  test('Triangle copy: same properties as previous test', () => {
    let a = new $3Dmol.Vector3(1, 1, 1);
    let b = new $3Dmol.Vector3(2, 2, 2);
    let c = new $3Dmol.Vector3(3, 3, 3);
    let triA = new $3Dmol.Triangle(a, b, c);
    let triB = triA.copy(triA);
    expect(triA.a).toEqual(triB.a);
    expect(triA.b).toEqual(triB.b);
    expect(triA.c).toEqual(triB.c);
  });
  test(`Triangle applyMatrix4: same properties as previous test, mat: [1-16]
		  Return: a: <10, 26, 42>, b: <16, 44, 72>, c: <22, 62, 102>`, () => {
    let a = new $3Dmol.Vector3(1, 1, 1);
    let b = new $3Dmol.Vector3(2, 2, 2);
    let c = new $3Dmol.Vector3(3, 3, 3);
    let tri = new $3Dmol.Triangle(a, b, c);
    let mat = new $3Dmol.Matrix4(...range(1, 16));
    tri = tri.applyMatrix4(mat);
    expect(tri.a).toEqual({x: 10, y: 26, z: 42});
    expect(tri.b).toEqual({x: 16, y: 44, z: 72});
    expect(tri.c).toEqual({x: 22, y: 62, z: 102});
  });
});
