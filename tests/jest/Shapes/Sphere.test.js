/**
 *@jest-environment jsdom
 */
const $3Dmol = require('../../../build/3Dmol.js');
let range = (start, end) => [...Array.from(new Array(end + 1).keys()).slice(start)];
describe('Sphere class', () => {
  test('Sphere constructor, center: <1, 2, 3>, radius: 9', () => {
    let center = new $3Dmol.Vector3(1, 2, 3);
    let radius = 9;
    let sphere = new $3Dmol.Sphere(center, radius);
    expect(sphere.center).toEqual({x: 1, y: 2, z: 3});
    expect(sphere.radius).toBe(9);
  });
  test('Sphere set, center: <1, 2, 3>, radius: 9', () => {
    let sphere = new $3Dmol.Sphere();
    let center = new $3Dmol.Vector3(1, 2, 3);
    let radius = 9;
    sphere = sphere.set(center, radius);
    expect(sphere.center).toEqual({x: 1, y: 2, z: 3});
    expect(sphere.radius).toBe(9);
  });
  test('Sphere copy, center: <4, 5, 6>, radius: 4', () => {
    let sphereA = new $3Dmol.Sphere();
    let center = new $3Dmol.Vector3(4, 5, 6);
    let radius = 4;
    let sphereB = new $3Dmol.Sphere(center, radius);
    sphereA = sphereA.copy(sphereB);
    expect(sphereA.center).toEqual({x: 4, y: 5, z: 6});
    expect(sphereA.radius).toBe(4);
  });
  test(`Sphere apply matrix 4, center: <1, 2, 3>, radius: 9, 
		  matrix: [1-16], Return center: <18, 46, 74>, 
		  radius: 9 * sqrt(179)`, () => {
    let center = new $3Dmol.Vector3(1, 2, 3);
    let radius = 9;
    let sphere = new $3Dmol.Sphere(center, radius);
    let mat = new $3Dmol.Matrix4(...range(1, 16));
    sphere = sphere.applyMatrix4(mat);
    expect(sphere.center).toEqual({x: 18, y: 46, z: 74});
    expect(sphere.radius).toBe(9 * Math.sqrt(179));
  });
  test(`Sphere translate, center: <1, 2, 3>, radius: 9, 
		  offset: <1, 1, 1>, new center will be: <2, 3, 4>`, () => {
    let center = new $3Dmol.Vector3(1, 2, 3);
    let radius = 9;
    let sphere = new $3Dmol.Sphere(center, radius);
    let offset = new $3Dmol.Vector3(1, 1, 1);
    sphere = sphere.translate(offset);
    expect(sphere.center).toEqual({x: 2, y: 3, z: 4});
    expect(sphere.radius).toBe(9);
  });
  test(`Sphere equals, center: <1, 2, 3>, radius: 9, two of. 
		  Center's and radii will be equal`, () => {
    let center = new $3Dmol.Vector3(1, 2, 3);
    let radius = 9;
    let sphereA = new $3Dmol.Sphere(center, radius);
    let sphereB = new $3Dmol.Sphere(center, radius);
    expect(sphereA.center).toEqual(sphereB.center);
    expect(sphereA.radius).toEqual(sphereB.radius);
  });
  test(`Sphere clone, center: <1, 2, 3>, radius: 9, two of.
	      Center's and radii will be equal, different refs`, () => {
    let center = new $3Dmol.Vector3(1, 2, 3);
    let radius = 9;
    let sphereA = new $3Dmol.Sphere(center, radius);
    let sphereB = new $3Dmol.Sphere(center, radius);
    expect(sphereA.center).toEqual(sphereB.center);
    expect(sphereA.radius).toEqual(sphereB.radius);
    expect(sphereA).not.toBe(sphereB);
  });
});
