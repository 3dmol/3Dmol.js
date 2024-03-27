/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function () { };
var $3Dmol = require("3dmol");

//A poorman's JS version of Python's range
var range = (start, end) => [...Array.from(Array(end + 1).keys()).slice(start)];
var nDigits = (digit, n) => Array.from(`${digit}`.repeat(n), Number);
const piValue = Math.PI

describe('Math Tests', () => {
	describe('Clamp Tests', () => {
		describe.each([
			[12, 6, 25, 12],
			[1, 2, 3, 2],
			[100, 25, 75, 75],
		])('clamp(%i, %i, %i)', (x, min, max, expected) => {
			test(`returns ${expected}`, () => {
			expect($3Dmol.clamp(x, min, max)).toBe(expected);
			});
		});
	});

	describe('DegtoRad Tests', () => {
		const testCases = [
			[180, piValue],
			[45, piValue / 4],
		];

		test.each(testCases)('degToRad(%i) returns %f', (input, expected) => {
			expect($3Dmol.degToRad(input)).toBeCloseTo(expected);
		});
	});
})

describe("Quaternion Tests", () => {
	describe("Quat Constructor Tests", () => {
		test("Constructor all args given", () => {
			const quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			expect(quat).toEqual({ x: 1, y: 2, z: 3, w: 4 });
		})
		test("Constructor no args given", () => {
			const quat = new $3Dmol.Quaternion();
			expect(quat).toEqual({ x: 0, y: 0, z: 0, w: 1 });
		})
	});
	describe("Quat Methods Tests", () => {
		test("Set with 1, 2, 3, 4", () => {
			let quat = new $3Dmol.Quaternion();
			quat = quat.set(1, 2, 3, 4);
			expect(quat).toEqual({ x: 1, y: 2, z: 3, w: 4 });
		})
		test("Copy, make sure values are equal but references are not", () => {
			let q = new $3Dmol.Quaternion();
			q = q.copy(new $3Dmol.Quaternion(1, 2, 3, 4));
			expect(q).toEqual({ x: 1, y: 2, z: 3, w: 4 });
		})
		test("Conjugate with 1, -2, -3, 4: Signs should flip but 4 shouldn't change", () => {
			let quat = new $3Dmol.Quaternion(1, -2, -3, 4);
			quat = quat.conjugate();
			expect(quat).toEqual({ x: -1, y: 2, z: 3, w: 4 });
		})
		test("Inverse with 1, 2, 3, 4: Flip x-z's signs and multiply by 1/sqrt(30)", () => {
			let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			quat = quat.inverse();
			//nf is normalizing factor
			const nf = 1 / Math.sqrt(30);
			expect(quat).toEqual({ x: -1 * nf, y: -2 * nf, z: -3 * nf, w: 4 * nf });
		})
		test("Length with 1, 2, 3, 4: Should be sqrt(30)", () => {
			const quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			const ans = quat.length();
			expect(ans).toBe(Math.sqrt(30));
		})
		test("Lengthxyz with 1, 2, 3,4: Should be sqrt(14)", () => {
			const quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			const ans = quat.lengthxyz();
			expect(ans).toBe(Math.sqrt(14));
		})
		test("Normalize, length=0, x,y,z = 0, w =1", () => {
			let quat = new $3Dmol.Quaternion(0, 0, 0, 0);
			quat = quat.normalize();
			expect(quat).toEqual({ x: 0, y: 0, z: 0, w: 1 });
		})
		test("Normalize, non-zero length: All's values multiplied by length", () => {
			let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			quat = quat.normalize();

			//nf is normalizing factor
			const nf = 1 / Math.sqrt(30);
			expect(quat).toEqual({ x: nf * 1, y: nf * 2, z: nf * 3, w: nf * 4 });
		})
		test("Multiply Quaternions: Test the function", () => {
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let p = new $3Dmol.Quaternion(5, 6, 7, 8);
			q = q.multiply(p);
			expect(q).toEqual({ x: 24, y: 48, z: 48, w: -6 });
		})
		test("Multiply Scalar: Should multiply all values by 2", () => {
			let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			const scalar = 2;
			quat = quat.multiplyScalar(2);
			expect(quat).toEqual({ x: 2, y: 4, z: 6, w: 8 });
		})
		test("Multiply Quaternions: Should do proper operation", () => {
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let p = new $3Dmol.Quaternion(5, 6, 7, 8);
			let resultQuat = q.multiplyQuaternions(q, p);
			expect(resultQuat).toEqual({ x: 24, y: 48, z: 48, w: -6 })
		})
		test("Subtract Quaternions: Should do prop - argProp for each property", () => {
			let p = new $3Dmol.Quaternion(5, 6, 7, 8);
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			p = p.sub(q);
			expect(p).toEqual({ x: 4, y: 4, z: 4, w: 4 });
		})
		test("Clone a Quaternion, Equals props but different references", () => {
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let p = q.clone();
			expect(q).toEqual(p);
			expect(q).not.toBe(p);
		})
		test("Set Quaternion from Euler, Converts an euler x,y,z to a Quaternion", () => {
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let e = { x: 2, y: 4, z: 6 };
			q = q.setFromEuler(e);

			const c1 = Math.cos(1);
			const c2 = Math.cos(2);
			const c3 = Math.cos(3);
			const s1 = Math.sin(1);
			const s2 = Math.sin(2);
			const s3 = Math.sin(3);

			const newX = s1 * c2 * c3 + c1 * s2 * s3;
			const newY = c1 * s2 * c3 - s1 * c2 * s3;
			const newZ = c1 * c2 * s3 + s1 * s2 * c3;
			const newW = c1 * c2 * c3 - s1 * s2 * s3;

			expect(q).toEqual({ x: newX, y: newY, z: newZ, w: newW });
		})
	})
})

describe("Vectors Tests", () => {
	describe("2D Vec Tests", () => {
		test("Vector constructor: 1,2", () => {
			const vec = new $3Dmol.Vector2(1, 2);
			expect(vec).toEqual({ x: 1, y: 2 });
		});

		test("Vector constructor no args", () => {
			const vec = new $3Dmol.Vector2();
			expect(vec).toEqual({ x: 0, y: 0 });
		});

		test("Vector set: 1,2", () => {
			const vec = new $3Dmol.Vector2();
			vec.set(1, 2);
			expect(vec).toEqual({ x: 1, y: 2 });
		});

		test("Vector subtraction: <3, 4> - <1, 2>", () => {
			const a = new $3Dmol.Vector2(3, 4);
			const b = new $3Dmol.Vector2(1, 2);
			let c = new $3Dmol.Vector2();
			c = c.subVectors(a, b);
			expect(c).toEqual({ x: 2, y: 2 })
		});

		test("Vector copy: <1, 2> Same values", () => {
			const a = new $3Dmol.Vector2();
			const b = new $3Dmol.Vector2(1, 2);
			a.copy(b);
			expect(a).toEqual({ x: 1, y: 2 });
		});

		test("Vector clone: <1, 2> Same values, different reference", () => {
			const a = new $3Dmol.Vector2(1, 2);
			const aRef = a;
			const aClone = a.clone();
			expect(aClone).toEqual({ x: 1, y: 2 });
			expect(aClone).not.toBe(aRef);
		});
	})

	describe("3D Vec Tests", () => {
		test("Vector constructor: 1,2", () => {
			const vec = new $3Dmol.Vector3(1, 2, 3);
			expect(vec).toEqual({ x: 1, y: 2, z: 3 });
		})
		test("Vector constructor no args", () => {
			const vec = new $3Dmol.Vector3();
			expect(vec).toEqual({ x: 0.0, y: 0.0, z: 0.0 });
		})
		test("Vector set: 1,2,3", () => {
			const vec = new $3Dmol.Vector3();
			vec.set(1, 2, 3);
			expect(vec).toEqual({ x: 1, y: 2, z: 3 });
		})
		test("Vector copy: <1, 2, 3> Same values", () => {
			let a = new $3Dmol.Vector3();
			const b = new $3Dmol.Vector3(1, 2, 3);
			a = a.copy(b);
			expect(a).toEqual({ x: 1, y: 2, z: 3 });
		})
		test("Vector add: <1, 2, 3> + <1, 1, 1>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(1, 1, 1);
			a = a.add(b);
			expect(a).toEqual({ x: 2, y: 3, z: 4 });
		})
		test("Add vectors: <1, 2, 3> + <1, 1, 1>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(1, 1, 1);
			let c = new $3Dmol.Vector3();
			c = c.addVectors(a, b);
			expect(c).toEqual({ x: 2, y: 3, z: 4 });
		})
		test("Multiply vectors elementwise: <1, 2, 3> * <4, 5, 6>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(4, 5, 6);
			let c = new $3Dmol.Vector3();
			c = c.multiplyVectors(a, b);
			expect(c).toEqual({ x: 4, y: 10, z: 18 });
		})
		test("Subtract vectors: <4, 5, 6> - <1, 2, 3> ", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(4, 5, 6);
			let c = new $3Dmol.Vector3();
			c = c.multiplyVectors(a, b);
			expect(c).toEqual({ x: 4, y: 10, z: 18 });
		})
		test("Subtract vectors: <4, 5, 6> - <1, 2, 3> ", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(4, 5, 6);
			let c = new $3Dmol.Vector3();
			c = c.subVectors(b, a);
			expect(c).toEqual({ x: 3, y: 3, z: 3 });
		})
		test("Scale vector: <1, 2, 3> * 2", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			a = a.multiplyScalar(2);
			expect(a).toEqual({ x: 2, y: 4, z: 6 });
		})
		test("Shrink vector: <1, 2, 3> / 0", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			a = a.divideScalar(0);
			expect(a).toEqual({ x: 0, y: 0, z: 0 });
		})
		test("Shrink vector: <2, 4, 6> / 2", () => {
			let a = new $3Dmol.Vector3(2, 4, 6);
			a = a.divideScalar(2);
			expect(a).toEqual({ x: 1, y: 2, z: 3 });
		})
		test("Vector max: <1, 2, 3>, <4, 5, 6>: Return <4, 5, 6>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			a = a.max(b);
			expect(a).toEqual({ x: 4, y: 5, z: 6 });
		})
		test("Vector min: <1, 2, 3>, <4, 5, 6>: Return <1, 2, 3>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			a = a.min(b);
			expect(a).toEqual({ x: 1, y: 2, z: 3 });
		})
		test("Vector distance to from: <1, 2, 3> to <4, 5, 6>, return sqrt(27)", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			const dist = a.distanceTo(b);
			expect(dist).toBe(Math.sqrt(27));
		})
		test("Vector distance to squared from: <1, 2, 3> to <4, 5, 6>, Return 27", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			const distSquared = a.distanceToSquared(b);
			expect(distSquared).toBe(27);
		})
		test("Apply matrix3 [1 2 3 4 5 6 7 8 9] to vector <1, 2, 3> Return <14, 32, 50>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix3(...range(1, 9));
			a = a.applyMatrix3(m);
			expect(a).toEqual({ x: 14, y: 32, z: 50 });
		})
		test("Apply matrix4 [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] to vector <1, 2, 3> Return <7, 7, 7>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(...nDigits(1, 16));
			a = a.applyMatrix4(m);
			expect(a).toEqual({ x: 7, y: 7, z: 7 });
		})
		test("Apply projection [1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0] to vector <1, 2, 3> Return <INF, INF, INF>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0);
			a = a.applyProjection(m);
			expect(a).toEqual({ x: Infinity, y: Infinity, z: Infinity });
		})
		test("Apply quaternion 1, 1, 1, 1 to vector <1, 2, 3> Return <5, 6, -5> ", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let q = new $3Dmol.Quaternion(1, 1, 1, 1);
			a = a.applyQuaternion(q);
			expect(a).toEqual({ x: 5, y: 6, z: -5 });
		})
		test("Vector negation: <1, 2, 3>, Return <-1, -2, -3>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			a = a.negate();
			expect(a).toEqual({ x: -1, y: -2, z: -3 });
		})
		test("Vector dot product: <1, 2, 3> . <4, 5, 6>, Return 32", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			const c = a.dot(b);
			expect(c).toBe(32);
		})
		test("Vector length: <3, 4, 12>, Return 13", () => {
			let a = new $3Dmol.Vector3(3, 4, 12);
			const c = a.length();
			expect(c).toBe(13);
		})
		test("Vector length squared: <1, 2, 3>, Return 14", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const c = a.lengthSq();
			expect(c).toBe(14);
		})
		test("Vector normalize: <3, 4, 12>, Return <3/13, 4/13, 12/13>", () => {
			let a = new $3Dmol.Vector3(3, 4, 12);
			a = a.normalize();
			expect(a).toEqual({ x: 3 / 13, y: 4 / 13, z: 12 / 13 })
		})
		test("Vector cross product: <1, 2, 3> x <4, 5, 6>. Return <-3, 6, -3>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			a = a.cross(b);
			expect(a).toEqual({ x: -3, y: 6, z: -3 });
		})
		test("Vector cross product new obj: <1, 2, 3> x <4, 5, 6>. Return <-3, 6, -3>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			let c = new $3Dmol.Vector3();
			c = c.crossVectors(a, b);
			expect(c).toEqual({ x: -3, y: 6, z: -3 })
		})
		test("Vector get position from matrix: [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15], Return <3, 7, 11>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let matArr = range(0, 16)
			let m = new $3Dmol.Matrix4(...matArr);
			a = a.getPositionFromMatrix(m);
			expect(a).toEqual({ x: 3, y: 7, z: 11 });
		})
		test("Vector setEulerFromRotationMatrix bogus order Call console.error", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let matArr = range(0, 16)
			let m = new $3Dmol.Matrix4(...matArr);
			a = a.setFromMatrixPosition(m);
			expect(a).toEqual({ x: 3, y: 7, z: 11 });
		})
		test("Vector setEulerFromRotationMatrix element 13 < 0.9999999, Return <0, 0, PI/2>", () => {
			//Necessary wrapper to prevent a console.log wall
			let prevConsoleError = console.error;
			console.error = (() => null);

			let a = new $3Dmol.Vector3(1, 2, 3);
			let matArr = range(0, 16)
			let m = new $3Dmol.Matrix4(...matArr);
			const consoleSpy = jest.spyOn(console, 'error');
			a.setEulerFromRotationMatrix(m, "ZZZ");
			expect(consoleSpy).toHaveBeenCalledWith("Error with vector's setEulerFromRotationMatrix: Unknown order: ZZZ");

			//Necessary wrapper to prevent a console.log wall
			console.error = prevConsoleError
		})
		test("Vector setEulerFromRotationMatrix element 13 > 0.9999999, Return <PI/2, 1, 0>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(0, -1, 0, 24, 1, 0, 0, 24, 0, 0, 1, 24, 24, 24, 24, 24)
			a = a.setEulerFromRotationMatrix(m, "XYZ");
			expect(a).toEqual({ x: -0, y: 0, z: piValue / 2 });
		})
		test("Vector rotateAboutVector: <1, 2, 3>, ang = PI/2 rads, Return same vector with all divided by sqrt(14)", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(0, 0, 1, 24, 1, 0, 0, 24, 0, 1, 0, 24, 24, 24, 24, 24)
			a = a.setEulerFromRotationMatrix(m, "XYZ");
			expect(a).toEqual({ x: piValue / 2, y: piValue / 2, z: 0 });
		})
		test("Vector set position from matrix: [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15], Return <3, 7, 11>", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const ang = piValue / 2;
			a = a.rotateAboutVector(a, ang);
			const sqrt14 = Math.sqrt(14);
			expect(a.x).toBeCloseTo(1 / sqrt14);
			expect(a.y).toBeCloseTo(2 / sqrt14);
			expect(a.z).toBeCloseTo(3 / sqrt14);
		})
		test("Vector transfromDirection from 0-15 4x4 matrix, Return <1, 2, 7> all divided by sqrt(66)", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			let matArr = range(0, 16);
			let m = new $3Dmol.Matrix4(...matArr);
			a = a.transformDirection(m);
			const sqrt66 = Math.sqrt(66);
			expect(a.x).toBeCloseTo(1 / sqrt66);
			expect(a.y).toBeCloseTo(4 / sqrt66);
			expect(a.z).toBeCloseTo(7 / sqrt66);
		})
		test("Vector clone: <1, 2, 3> Same values, different reference", () => {
			let a = new $3Dmol.Vector3(1, 2, 3);
			const aRef = a;
			a = a.clone(a);
			expect(a).toEqual({ x: 1, y: 2, z: 3 });
			expect(a).not.toBe(aRef);
		})
	})
})

describe("Matrix Tests", () => {
	describe("3x3 Matrix Tests", () => {
			test("Matrix 3 Constructor, [1 2 3 4 5 6 7 8 9]", () => {
			let mat = new $3Dmol.Matrix3(...range(1, 9));
			expect([...mat.elements]).toEqual([1, 4, 7, 2, 5, 8, 3, 6, 9]);
		})
		test("Matrix 3 Set, [1 1 1 1 1 1 1 1 1]", () => {
			let mat = new $3Dmol.Matrix3();
			mat = mat.set(...nDigits(1, 9))
			expect([...mat.elements]).toEqual(nDigits(1, 9));
		})
		test("Matrix 3 Identity, [1 0 0 0 1 0 0 1]", () => {
			let mat = new $3Dmol.Matrix3();
			mat = mat.identity();
			expect([...mat.elements]).toEqual([1, 0, 0, 0, 1, 0, 0, 0, 1]);
		})
		test("Matrix 3 Copy, [0 0 0 0 0 0 0 0 0]", () => {
			let mat = new $3Dmol.Matrix3();
			let elems = nDigits(0, 9);
			let from = new $3Dmol.Matrix3(...elems);
			mat.copy(from);
			expect([...mat.elements]).toEqual(elems);
		})
		test("Matrix 3 Multiply Scalar, [1 1 1 1 1 1 1 1 1] by 2", () => {
			let mat = new $3Dmol.Matrix3(...nDigits(1, 9));
			mat = mat.multiplyScalar(2);	
			expect([...mat.elements]).toEqual(nDigits(2, 9));
		})
		test("Matrix 3 getInverse3, [1, 2, 3, 4, 5, 6, 7, 8, 24]", () => {
			let mat = new $3Dmol.Matrix3();
			let given = new $3Dmol.Matrix3(1, 2, 3, 4, 5, 6, 7, 8, 24);
			mat = mat.getInverse3(given);
			expect(mat.elements[0]).toBeCloseTo(-8 / 5);
			expect(mat.elements[1]).toBeCloseTo(1.2);
			expect(mat.elements[2]).toBeCloseTo(1 / 15);
			expect(mat.elements[3]).toBeCloseTo(8 / 15);
			expect(mat.elements[4]).toBeCloseTo(-1 / 15);
			expect(mat.elements[5]).toBeCloseTo(-2 / 15);
			expect(mat.elements[6]).toBeCloseTo(1 / 15);
			expect(mat.elements[7]).toBeCloseTo(-2 / 15);
			expect(mat.elements[8]).toBeCloseTo(1 / 15);
		})
		test("Matrix 3 getInverse Inputs: Random stuff that works for happy path Throw on invertible = false", () => {
			let mat3 = new $3Dmol.Matrix3(9, 5, 7, 3, 2, 9, 5, 7, 2);
			let mat4 = new $3Dmol.Matrix4(7, 10, 6, 8, 9, 2, 4, 6, 2, 7, 4, 1, 3, 4, 10, 8);
			let toI = false;
			mat3 = mat3.getInverse(mat4, toI);
			let expectedElems = [-20, -28, 59, 2, 16, -29, 28, 27, -76]
			expectedElems = expectedElems.map(x => x / -66);
			for (let i = 0; i < expectedElems.length; i++) {
				expect(mat3.elements[i]).toBeCloseTo(expectedElems[i], 1);
			}
		})
		test("Matrix 3 getInverse Inputs: 3x3 of all 1's and 4x4 of all 1's, throw on invertible = false", () => {
			//Necessary wrapper to prevent a console.log wall
			let prevConsoleWarn = console.warn;
			console.warn = (() => null);

			let mat3 = new $3Dmol.Matrix3(...nDigits(1, 9));
			let mat4 = new $3Dmol.Matrix4(...nDigits(1, 16));
			let toI = false;
			const consoleSpy = jest.spyOn(console, 'warn');
			mat3 = mat3.getInverse(mat4, toI);
			expect(consoleSpy).toHaveBeenCalledWith("Matrix3.getInverse(): can't invert matrix, determinant is 0");

			//Necessary wrapper to prevent a console.log wall
			console.warn = prevConsoleWarn
		})
		test("Matrix 3 getInverse Inputs: 3x3 of all 1's and 4x4 of all 1's, throw on invertible = true", () => {
			let mat3 = new $3Dmol.Matrix3(...nDigits(1, 9));
			let mat4 = new $3Dmol.Matrix4(...nDigits(1, 16));
			let toI = true;
			expect(() => mat3 = mat3.getInverse(mat4, toI)).toThrow("Matrix3.getInverse(): can't invert matrix, determinant is 0");
		})
		test("Matrix 3 Determinant, [1, 2, 3, 4, 5, 6, 7, 8, 24], Return -45", () => {
			let mat = new $3Dmol.Matrix3(1, 2, 3, 4, 5, 6, 7, 8, 24);
			const det = mat.getDeterminant();
			expect(det).toBe(-45);
		})
		test("Matrix 3 Get Matrix 4, [1 1 1 1 1 1 1 1 1], Throw in 0s to pad out", () => {
			let mat = new $3Dmol.Matrix3(...nDigits(1, 9));
			let mat4 = mat.getMatrix4();
			expect([...mat4.elements]).toEqual([1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1])
		})
		test("Matrix 3 Transpose, [1 2 3 4 5 6 7 8 9], Return [1 4 7 2 5 8 3 6 9]", () => {
			/*These look not different, but you write to the object in row-major order
			* and read from the object in column-major order
			*/
			let mat = new $3Dmol.Matrix3(...range(1, 9));
			mat = mat.transpose();
			expect([...mat.elements]).toEqual(range(1, 9));
		})
		test("Matrix 3 Clone, [1 2 3 4 5 6 7 8 9], Return new Matrix with exact same values", () => {
			let mat = new $3Dmol.Matrix3(...range(1, 9));
			let matNew = mat.clone();
			expect(mat).not.toBe(matNew);
			expect([...mat.elements]).toEqual([...matNew.elements]);
		})
		test("Matrix 3 conversionMatrix3, some wacky stuff happens, Return [1 0 0 sqrt(2) sqrt(2) 0 1.5 (1.5 * (sqrt(6) -1)) really long arg]", () => {
			let mat3 = $3Dmol.conversionMatrix3(1, 2, 3, 30, 60, 45);
			let lastElem = (3 * Math.sqrt(2) * Math.sqrt(Math.sqrt(6) - 2)) / 2;
			let expected = [1, 0, 0, Math.sqrt(2), Math.sqrt(2), 0, 1.5, (1.5 * (Math.sqrt(6) - 1)), lastElem];
			for (let i = 0; i < mat3.elements.length; i++) {
				expect(mat3.elements[i]).toBeCloseTo(expected[i], 1);
			}
		})
	})

	describe("4x4 Matrix Tests", () => {
		test("Matrix 4 Constructor, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]", () => {
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
		})
		test("Matrix 4 Set, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]", () => {
			let mat = new $3Dmol.Matrix4();
			mat = mat.set(...range(1, 16));
			expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
		})
		test("Matrix 4 Identity, [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1]", () => {
			let mat = new $3Dmol.Matrix4();
			mat = mat.identity();
			expect([...mat.elements]).toEqual([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
		})
		test("Matrix 4 Copy, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]", () => {
			let mat = new $3Dmol.Matrix4();
			mat = mat.copy(new $3Dmol.Matrix4(...range(1, 16)));
			expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
		})
		test("Matrix 4 Matrix3 from Top left, [1 - 16] -> [1 5 9 2 6 10 3 7 11] ", () => {
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			let mat3 = mat.matrix3FromTopLeft();
			expect([...mat3.elements]).toEqual([1, 5, 9, 2, 6, 10, 3, 7, 11]);
		})
		test("Matrix 4 setRotationFromEuler bogus order, call console.error", () => {
			//Necessary wrapper to prevent a console.log wall
			let prevConsoleError = console.error;
			console.error = (() => null);

			let mat = new $3Dmol.Matrix4(...range(1, 16));
			let vec = new $3Dmol.Vector3(1, 2, 3);
			const consoleSpy = jest.spyOn(console, 'error');
			mat.setRotationFromEuler(vec, "ZZZ");
			expect(consoleSpy).toHaveBeenCalledWith("Error with matrix4 setRotationFromEuler. Order: ZZZ");
			//Necessary wrapper to prevent a console.log wall
			console.error = prevConsoleError
		})
		test("Matrix 4 setRotationFromEuler with [1 - 16] and vec: <1, 2, 3>", () => {
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			let vec = new $3Dmol.Vector3(piValue / 2, piValue / 2, piValue / 2);
			mat = mat.setRotationFromEuler(vec, "XYZ");

			//Pretreat some math weirdness
			mat.elements[0] = Math.round(mat.elements[0]);
			mat.elements[1] = Math.round(mat.elements[1]);
			mat.elements[4] = Math.abs(Math.round(mat.elements[4]));
			mat.elements[6] = Math.round(mat.elements[6]);
			mat.elements[9] = Math.abs(Math.round(mat.elements[9]));
			mat.elements[10] = Math.round(mat.elements[10]);

			expect([...mat.elements]).toEqual([0, 0, 1, 13, 0, -1, 0, 14, 1, 0, 0, 15, 4, 8, 12, 16]);
		})
		test("Matrix 4 setRotationFromQuaternion with [1 - 16] and quat: (1, 2, 3, 4)", () => {
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			mat = mat.setRotationFromQuaternion(quat);
			expect([...mat.elements]).toEqual([-25, 28, -10, 13, -20, -19, 20, 14, 22, 4, -9, 15, 4, 8, 12, 16]);
		})
		//test("Matrix 4 lookAt, with [1-16], Return ", mat4lookAt)
		test("Matrix 4 multplication wtih [1 - 16] and [1 - 16', Return [90 100 110 120 202 228 254 280 314 356 398 440 426 484 542 600]", () => {
			let mat = new $3Dmol.Matrix4();
			let a = new $3Dmol.Matrix4(...range(1, 16));
			let b = new $3Dmol.Matrix4(...range(1, 16));
			mat = mat.multiplyMatrices(a, b);
			expect([...mat.elements]).toEqual([90, 202, 314, 426, 100, 228, 356, 484, 110, 254, 398, 542, 120, 280, 440, 600]);
		})
		test("Matrix 4 multiply scalar all 1's with 2, Return all 2's", () => {
			let mat = new $3Dmol.Matrix4(nDigits(1, 16));
			const scalar = 2;
			mat = mat.multiplyScalar(scalar);
			expect([...mat.elements]).toEqual(nDigits(2, 16));
		})
		test("Matrix 4 makeTranslation [1-16]", () => {
			let mat = new $3Dmol.Matrix4();
			mat = mat.makeTranslation(1, 2, 3);
			expect([...mat.elements]).toEqual([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 2, 3, 1]);
		})
		test("Matrix 4 snap [1-16] with all elements + 0.00001, Return [1-16]", () => {
			let elems = range(1, 16);
			elems = elems.map(x => x + 0.00001);
			let mat = new $3Dmol.Matrix4(elems);
			mat = mat.snap(4);
			expect([...mat.elements]).toEqual(range(1, 16));
		})
		test("Matrix 4 transpose [1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16], Return [1-16]", () => {
			let mat = new $3Dmol.Matrix4([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
			mat = mat.transpose();
			expect([...mat.elements]).toEqual(range(1, 16));
		})
		//test("Matrix 4 getPostion, not sure how to test closures", mat4getPosition)
		test("Matrix 4 setPosition, [1 - 16] with <1, 2, 3>, Return [1 2 3 1 5 6 7 2 9 10 11 3 13 14 15 16]", () => {
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			let vec = new $3Dmol.Vector3(1, 2, 3);
			mat = mat.setPosition(vec);
			expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 1, 2, 3, 16]);
		})
		test("Matrix 4 translate, [1 - 16] with <1, 2, 3>, Return [1 2 3 5 5 6 7 10 9 10 11 15 13 14 15 16]", () => {
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			let vec = new $3Dmol.Vector3(1, 2, 3);
			mat = mat.translate(vec);
			expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 5, 10, 15, 16]);
		})
		test("Matrix 4 getInverse, random inputs that make the happy path pass, throw on invertible = false", () => {
			let matA = new $3Dmol.Matrix4(7, 10, 6, 8, 9, 2, 4, 6, 2, 7, 4, 1, 3, 4, 10, 8);
			let matB = new $3Dmol.Matrix4(7, 10, 6, 8, 9, 2, 4, 6, 2, 7, 4, 1, 3, 4, 10, 8);
			matA = matA.getInverse(matB, false);
			let elems = [160, -194, 364, -418, -340, 140, -130, 220, -200, -104, -356, 572, 120, 102, -222, -66];
			elems = elems.map(x => x / -1980)
			for (let i = 0; i < elems.length; i++) {
				expect(matA.elements[i]).toBeCloseTo(elems[i], 1);
			}
		})
		test("Matrix 4 getInverse, two [1-16] matrixes with throw on invertible = false", () => {
			//Necessary wrapper to prevent a console.log wall
			let prevConsoleWarn = console.warn;
			console.warn = (() => null);

			let matA = new $3Dmol.Matrix4(...range(1, 16));
			let matB = new $3Dmol.Matrix4(...range(1, 16));
			const consoleSpy = jest.spyOn(console, 'warn');
			matA = matA.getInverse(matB, false);

			expect(consoleSpy).toHaveBeenCalledWith("Matrix4.getInverse(): can't invert matrix, determinant is 0");
			//Necessary wrapper to prevent a console.log wall
			console.warn = prevConsoleWarn
		})
		test("Matrix 4 getInverse, two [1-16] matrixes with throw on invertible = true", () => {
			let matA = new $3Dmol.Matrix4(...range(1, 16));
			let matB = new $3Dmol.Matrix4(...range(1, 16));
			expect(() => matA = matA.getInverse(matB, true)).toThrow("Matrix4.getInverse(): can't invert matrix, determinant is 0");
		})
		test("Matrix 4 isReflected, upperLeft corner  of 4x4 matrix is [1 2 3 4 5 6 7 8 10], Return true", () => {
			let mat = new $3Dmol.Matrix4(1, 2, 3, 0, 4, 5, 6, 0, 7, 8, 10, 0, 0, 0, 0, 0);
			const reflected = mat.isReflected();
			expect(reflected).toBe(true);
		})
		test("Matrix 4 isReflected, upperLeft corner  of 4x4 matrix is [1 2 3 4 5 6 7 8 10], Return false", () => {
			let mat = new $3Dmol.Matrix4(1, 2, 3, 0, 4, 5, 6, 0, 7, 8, -9, 0, 0, 0, 0, 0);
			const reflected = mat.isReflected();
			expect(reflected).toBe(false);
		})
		test("Matrix 4 scale, scale a matrix of all 1's with vector <2, 2, 2>, Return all 2's expect last 4 are 1's", () => {
			let mat4 = new $3Dmol.Matrix4(...nDigits(1, 16));
			let v = new $3Dmol.Vector3(2, 2, 2);
			mat4 = mat4.scale(v);
			let expected = nDigits(2, 16);
			for (let i = 12; i < 16; i++) {
				expected[i] = 1;
			}
			expect([...mat4.elements]).toEqual(expected);
		})
		test("Matrix 4 getMaxScaleOnAxis, a single [1-16] matrix, Return sqrt(179)", () => {
			let mat4 = new $3Dmol.Matrix4(...range(1, 16));
			const ans = mat4.getMaxScaleOnAxis();
			expect(ans).toEqual(Math.sqrt(179));
		})
		test("Matrix 4 makeFrustum, a single [1-16] matrix, args in order: 1, 2, 3, 4, 5, 6, Return [10 0 0 0 0 10 0 0 3 7 -11 -1 0 0 -60 0]", () => {
			let mat4 = new $3Dmol.Matrix4(...range(1, 16));
			mat4 = mat4.makeFrustum(1, 2, 3, 4, 5, 6);
			expect([...mat4.elements]).toEqual([10, 0, 0, 0, 0, 10, 0, 0, 3, 7, -11, -1, 0, 0, -60, 0]);
		})
		test("Matrix 4 makePerspective, a single [1-16] matrix, args in order: 90, 1, 2, 3, Return [NaN 0 0 0 0 NaN 0 0 -1.0224719047546387 5 NaN -1 0 0 NaN 0]", () => {
			let mat4 = new $3Dmol.Matrix4(...range(1, 16));
			mat4 = mat4.makePerspective(90, 1, 2, 3);
			expect([...mat4.elements]).toEqual([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -5, -1, 0, 0, -12, 0]);
		})
		test("Matrix 4 makeOrthographic, a single [1-16] matrix, args in order: 1, 2, 3, 4, 5, 6, Return: [2 0 0 0 0 -2 0 0 0 0 -2 0 -3 7 -11 1]", () => {
			let mat4 = new $3Dmol.Matrix4(...range(1, 16));
			mat4 = mat4.makeOrthographic(1, 2, 3, 4, 5, 6);
			expect([...mat4.elements]).toEqual([2, 0, 0, 0, 0, -2, 0, 0, 0, 0, -2, 0, -3, 7, -11, 1]);
		})
		test("Matrix 4 isEqual, two all 1's 4x4 matricies, Return true", () => {
				let a = new $3Dmol.Matrix4(...nDigits(1, 16));
			let b = new $3Dmol.Matrix4(...nDigits(1, 16));
			let equal = a.isEqual(b);
			expect(equal).toBe(true);
		})
		test("Matrix 4 isEqual, one all 1's 4x4 matricies and other one is all 2's 4x4 matricies, Return false", () => {
			let a = new $3Dmol.Matrix4(...nDigits(1, 16));
			let b = new $3Dmol.Matrix4(...nDigits(2, 16));
			let equal = a.isEqual(b);
			expect(equal).toBe(false);
		})
		test("Matrix 4 clone, all 1's 4x4 matricies, Return identical matrix with different reference", () => {
			let a = new $3Dmol.Matrix4(...nDigits(1, 16));
			let b = a.clone();
			expect([...a.elements]).toEqual([...b.elements]);
			expect(a).not.toBe(b);
		})
		test("Matrix 4 isIdentity, all 0's expect 1's down the diagonal, Return true", () => {
			let mat = new $3Dmol.Matrix4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
			let ident = mat.isIdentity();
			expect(ident).toEqual(true);
		})
		test("Matrix 4 isIdentity, [1-16], Return false", () => {
			let mat = new $3Dmol.Matrix4(...range(1, 16));
			let ident = mat.isIdentity();
			expect(ident).toEqual(false);
		})
		test("Matrix 4 isNearlyIdentity, all 0's expect 1's down the diagonal, add 0.00001 to all elements, Return true", () => {
			let elems = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
			elems = elems.map(x => x + 0.00001);
			let mat = new $3Dmol.Matrix4(...elems);
			let ident = mat.isNearlyIdentity(4);
			expect(ident).toEqual(true);
		})
	})
})


describe("Ray Tests", () => {
	const origin = new $3Dmol.Vector3(1, 2, 3);
	const direction = new $3Dmol.Vector3(4, 5, 6);

	test("Ray constructor, origin: <1, 2, 3>, direction: <4, 5, 6>", () => {
		const ray = new $3Dmol.Ray(origin, direction);
		expect(ray.origin).toEqual(origin);
		expect(ray.direction).toEqual(direction);
	});

	test("Ray set, origin: <1, 2, 3>, direction: <4, 5, 6>", () => {
		const ray = new $3Dmol.Ray().set(origin, direction);
		expect(ray.origin).toEqual(origin);
		expect(ray.direction).toEqual(direction);
	});

	test("Ray copy, origin: <1, 2, 3> direction: <4, 5, 6>", () => {
		const ray = new $3Dmol.Ray(origin, direction);
		const newRay = ray.copy(ray);
		expect(newRay.origin).toEqual(origin);
		expect(newRay.direction).toEqual(direction);
	});

	test("Ray at, scalar: 24, origin: <1, 2, 3>, direction: <4, 5, 6>, No optional target, Return <97, 122, 147>", () => {
		const ray = new $3Dmol.Ray(origin, direction);
		const resultVec = ray.at(24);
		expect(resultVec).toEqual(new $3Dmol.Vector3(97, 122, 147));
	});

	test("Ray at, scalar: 24, origin: <1, 2, 3>, direction: <4, 5, 6>, optional target <1, 1, 1>, Return <97, 122, 147>", () => {
		const ray = new $3Dmol.Ray(origin, direction);
		const opt = new $3Dmol.Vector3(1, 2, 3);
		const resultVec = ray.at(24, opt);
		expect(resultVec).toEqual(new $3Dmol.Vector3(97, 122, 147));
	});

	test("Ray closestPointToPoint, origin: <1, 2, 3>, direction: <4, 5, 6>, point: <1, 1, 1>, Return <-67, -83, -99>", () => {
		const ray = new $3Dmol.Ray(origin, direction);
		const point = new $3Dmol.Vector3(1, 1, 1);
		const result = ray.closestPointToPoint(point);
		expect(result).toEqual(new $3Dmol.Vector3(-67, -83, -99));
	});

	test("Ray isIntersectionCylinder, Empty function, Return undefined", () => {
		const ray = new $3Dmol.Ray(origin, direction);
		const result = ray.isIntersectionCylinder();
		expect(result).toBeUndefined();
	});

	test("Ray isIntersectionSphere, origin: <1, 2, 3>, direction: <4, 5, 6>, center: <1, 2, 3>, radius: 9, Return true", () => {
		const ray = new $3Dmol.Ray(origin, direction);
		const sphCent = new $3Dmol.Vector3(1, 2, 3);
		const sph = new $3Dmol.Sphere(sphCent, 9);
		const result = ray.isIntersectionSphere(sph);
		expect(result).toBe(true);
	});
	/*
	* isIntersectionPlane
	* distanceToPlane
	* intersectPlane
	* The plane class no longer exists?
	*/
	test("Ray applyMatrix4, origin: <1, 2, 3>, direction: <4, 5, 6>, Matrix4: [1-16], newOrig: <18, 46, 74> newDir: <32, 92, 152>", () => {
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray = new $3Dmol.Ray(origin, direction);
		let mat4 = new $3Dmol.Matrix4(...range(1, 16));
		ray = ray.applyMatrix4(mat4);
		expect({ ...ray.origin }).toEqual({ x: 18, y: 46, z: 74 });
		expect({ ...ray.direction }).toEqual({ x: 32, y: 92, z: 152 });
	})
	test("Ray equals, origin: <1, 2, 3>, direction: <4, 5, 6> for both, Throw an error because Vector3 doesn't have an equals method", () => {
		/*
		* As it stands, Vector3 doesn't have an equals function, so this can't be tested
		*/
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray1 = new $3Dmol.Ray(origin, direction);
		let ray2 = new $3Dmol.Ray(origin, direction);
		// the error message can change we just test that this function doesn't exist
		// follow up question, why does this function not exist?
		expect(() => ray1.equals(ray2)).toThrow();
	})
	test("Ray clone, origin: <1, 2, 3>, direction: <4, 5, 6>, Identical values but different references", () => {
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray1 = new $3Dmol.Ray(origin, direction);
		let ray2 = new $3Dmol.Ray(origin, direction);
		expect({ ...ray1.origin }).toEqual({ ...ray2.origin });
		expect({ ...ray1.direction }).toEqual({ ...ray2.direction });
		expect(ray1).not.toBe(ray2);
	})
})
