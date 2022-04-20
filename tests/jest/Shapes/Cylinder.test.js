/**
*@jest-environment jsdom
*/
const $3Dmol = require("../../../build/3Dmol.js");

let range = (start, end) => [...Array.from(Array(end + 1).keys()).slice(start)];

describe('Cylinder Tests', () => {
  test(`Cylinder constructor, c1: <1, 2, 3>, c2: <4, 5, 6>, radius: 9, 
		return with direction of: <3, 3, 3> all divided by sqrt(27)`, () => {
    let c1 = new $3Dmol.Vector3(1, 2, 3);
    let c2 = new $3Dmol.Vector3(4, 5, 6);
    let radius = 9;
    let cyl = new $3Dmol.Cylinder(c1, c2, radius);
    expect(cyl.c1).toEqual({x: 1, y: 2, z: 3});
    expect(cyl.c2).toEqual({x: 4, y: 5, z: 6});
    expect(cyl.direction).toEqual({
      x: 3 / Math.sqrt(27),
      y: 3 / Math.sqrt(27),
      z: 3 / Math.sqrt(27),
    });
    expect(cyl.radius).toBe(9);
  });
  test('Cylinder copy, same properties as previous test', () => {
    let c1 = new $3Dmol.Vector3(1, 2, 3);
    let c2 = new $3Dmol.Vector3(4, 5, 6);
    let radius = 9;
    let cyl1 = new $3Dmol.Cylinder(c1, c2, radius);
    let cyl2 = new $3Dmol.Cylinder(c1, c2, radius);
    expect(cyl1.c1).toEqual(cyl2.c1);
    expect(cyl1.c2).toEqual(cyl2.c2);
    expect(cyl1.direction).toEqual(cyl2.direction);
    expect(cyl1.radius).toBe(cyl2.radius);
  });
  test(`Cylinder applyMatrix4, same properties as previous test, Matrix: [1-16], 
		Return: c1: <18, 46, 74>, c2: <36, 100, 164>, radius: 9 * sqrt(179)
		and direction: <18, 54, 90> all divided by sqrt(11340) `, () => {
    let c1 = new $3Dmol.Vector3(1, 2, 3);
    let c2 = new $3Dmol.Vector3(4, 5, 6);
    let radius = 9;
    let cyl = new $3Dmol.Cylinder(c1, c2, radius);
    let mat = new $3Dmol.Matrix4(...range(1, 16));
    cyl = cyl.applyMatrix4(mat);
    expect(cyl.c1).toEqual({x: 18, y: 46, z: 74});
    expect(cyl.c2).toEqual({x: 36, y: 100, z: 164});
    expect(cyl.direction.x).toBeCloseTo(18 / Math.sqrt(11340));
    expect(cyl.direction.y).toBeCloseTo(54 / Math.sqrt(11340));
    expect(cyl.direction.z).toBeCloseTo(90 / Math.sqrt(11340));
    expect(cyl.radius).toBe(9 * Math.sqrt(179));
  });
});
