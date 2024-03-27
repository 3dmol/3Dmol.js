/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function () { };
let $3Dmol = require("../../build/3Dmol.js");

//A poorman's JS version of Python's range
let range = (start, end) => [...Array.from(Array(end + 1).keys()).slice(start)];
let nDigits = (digit, n) => Array.from(`${digit}`.repeat(n), Number)


describe("Cylinder Tests", () => {
	test(`Cylinder constructor, c1: <1, 2, 3>, c2: <4, 5, 6>, radius: 9, 
		return with direction of: <3, 3, 3> all divided by sqrt(27)`, () => {
			let c1 = new $3Dmol.Vector3(1, 2, 3);
			let c2 = new $3Dmol.Vector3(4, 5, 6);
			let radius = 9;
			let cyl = new $3Dmol.Cylinder(c1, c2, radius);
			expect(cyl.c1).toEqual({ x: 1, y: 2, z: 3 });
			expect(cyl.c2).toEqual({ x: 4, y: 5, z: 6 });
			expect(cyl.direction).toEqual({ x: 3 / Math.sqrt(27), y: 3 / Math.sqrt(27), z: 3 / Math.sqrt(27) });
			expect(cyl.radius).toBe(9);
		})
	test("Cylinder copy, same properties as previous test", () => {
		let c1 = new $3Dmol.Vector3(1, 2, 3);
		let c2 = new $3Dmol.Vector3(4, 5, 6);
		let radius = 9;
		let cyl1 = new $3Dmol.Cylinder(c1, c2, radius);
		let cyl2 = new $3Dmol.Cylinder(c1, c2, radius);
		expect(cyl1.c1).toEqual(cyl2.c1);
		expect(cyl1.c2).toEqual(cyl2.c2);
		expect(cyl1.direction).toEqual(cyl2.direction);
		expect(cyl1.radius).toBe(cyl2.radius);
	})
	test(`Cylinder applyMatrix4, same properties as previous test, Matrix: [1-16], 
		Return: c1: <18, 46, 74>, c2: <36, 100, 164>, radius: 9 * sqrt(179)
		and direction: <18, 54, 90> all divided by sqrt(11340) `, () => {
			let c1 = new $3Dmol.Vector3(1, 2, 3);
			let c2 = new $3Dmol.Vector3(4, 5, 6);
			let radius = 9;
			let cyl = new $3Dmol.Cylinder(c1, c2, radius);
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			cyl = cyl.applyMatrix4(mat);
			expect(cyl.c1).toEqual({ x: 18, y: 46, z: 74 });
			expect(cyl.c2).toEqual({ x: 36, y: 100, z: 164 });
			expect(cyl.direction.x).toBeCloseTo(18 / Math.sqrt(11340));
			expect(cyl.direction.y).toBeCloseTo(54 / Math.sqrt(11340));
			expect(cyl.direction.z).toBeCloseTo(90 / Math.sqrt(11340));
			expect(cyl.radius).toBe(9 * Math.sqrt(179));
	})
})

describe("Sphere Tests", () => {
	test("Sphere constructor, center: <1, 2, 3>, radius: 9", () => {
		let center = new $3Dmol.Vector3(1, 2, 3);
		let radius = 9;
		let sphere = new $3Dmol.Sphere(center, radius);
		expect(sphere.center).toEqual({ x: 1, y: 2, z: 3 });
		expect(sphere.radius).toBe(9);
	})
	test("Sphere set, center: <1, 2, 3>, radius: 9", () => {
		let sphere = new $3Dmol.Sphere();
		let center = new $3Dmol.Vector3(1, 2, 3);
		let radius = 9;
		sphere = sphere.set(center, radius);
		expect(sphere.center).toEqual({ x: 1, y: 2, z: 3 });
		expect(sphere.radius).toBe(9);
	})
	test("Sphere copy, center: <4, 5, 6>, radius: 4", () => {
		let sphereA = new $3Dmol.Sphere();
		let center = new $3Dmol.Vector3(4, 5, 6);
		let radius = 4;
		let sphereB = new $3Dmol.Sphere(center, radius);
		sphereA = sphereA.copy(sphereB);
		expect(sphereA.center).toEqual({ x: 4, y: 5, z: 6 });
		expect(sphereA.radius).toBe(4);
	})
	test(`Sphere apply matrix 4, center: <1, 2, 3>, radius: 9, 
		  matrix: [1-16], Return center: <18, 46, 74>, 
		  radius: 9 * sqrt(179)`, () => {
		let center = new $3Dmol.Vector3(1, 2, 3);
		let radius = 9;
		let sphere = new $3Dmol.Sphere(center, radius);
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		sphere = sphere.applyMatrix4(mat);
		expect(sphere.center).toEqual({ x: 18, y: 46, z: 74 });
		expect(sphere.radius).toBe(9 * Math.sqrt(179));
	})
	test(`Sphere translate, center: <1, 2, 3>, radius: 9, 
		  offset: <1, 1, 1>, new center will be: <2, 3, 4>`, () => {
		let center = new $3Dmol.Vector3(1, 2, 3);
		let radius = 9;
		let sphere = new $3Dmol.Sphere(center, radius);
		let offset = new $3Dmol.Vector3(1, 1, 1);
		sphere = sphere.translate(offset);
		expect(sphere.center).toEqual({ x: 2, y: 3, z: 4 });
		expect(sphere.radius).toBe(9);
	})
	test(`Sphere equals, center: <1, 2, 3>, radius: 9, two of. 
		  Center's and radii will be equal`, () => {
		let center = new $3Dmol.Vector3(1, 2, 3);
		let radius = 9;
		let sphereA = new $3Dmol.Sphere(center, radius);
		let sphereB = new $3Dmol.Sphere(center, radius);
		expect(sphereA.center).toEqual(sphereB.center);
		expect(sphereA.radius).toEqual(sphereB.radius);
	})
	test(`Sphere clone, center: <1, 2, 3>, radius: 9, two of.
	      Center's and radii will be equal, different refs`, () => {
		let center = new $3Dmol.Vector3(1, 2, 3);
		let radius = 9;
		let sphereA = new $3Dmol.Sphere(center, radius);
		let sphereB = new $3Dmol.Sphere(center, radius);
		expect(sphereA.center).toEqual(sphereB.center);
		expect(sphereA.radius).toEqual(sphereB.radius);
		expect(sphereA).not.toBe(sphereB);
	})
})

describe('Triangle Tests', () => {
	const a = new $3Dmol.Vector3(1, 1, 1);
	const b = new $3Dmol.Vector3(2, 2, 2);
	const c = new $3Dmol.Vector3(3, 3, 3);

	test("Triangle constructor: a is all 1's, b is all 2's, and c is all 3's", () => {
		const tri = new $3Dmol.Triangle(a, b, c);
		expect(tri.a).toEqual(a);
		expect(tri.b).toEqual(b);
		expect(tri.c).toEqual(c);
	})

	test("Triangle copy: same properties as previous test", () => {
		const triA = new $3Dmol.Triangle(a, b, c);
		const triB = triA.copy(triA);
		expect(triA).toEqual(triB);
	})

	test(`Triangle applyMatrix4: same properties as previous test, mat: [1-16]
		  Return: a: <10, 26, 42>, b: <16, 44, 72>, c: <22, 62, 102>`, () => {
		const tri = new $3Dmol.Triangle(a, b, c);
		const mat = new $3Dmol.Matrix4(...range(1, 16));
		tri.applyMatrix4(mat);
		expect(tri.a).toEqual(new $3Dmol.Vector3(10, 26, 42));
		expect(tri.b).toEqual(new $3Dmol.Vector3(16, 44, 72));
		expect(tri.c).toEqual(new $3Dmol.Vector3(22, 62, 102));
	})
})