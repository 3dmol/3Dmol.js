/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");

//A poorman's JS version of Python's range
let range = (start, end) => [...Array.from(Array(end + 1).keys()).slice(start)];
let nDigits = (digit, n) => Array.from(`${digit}`.repeat(n), Number)

describe("Math Tests", mathTests)
describe("Quaternion Tests", quatTests)
describe("Vectors Tests", vecTests)
describe("Matrix Tests", matTests)
describe("Ray Tests", rayTests)


function mathTests()
{
	describe("Clamp Tests", clampTests)
	describe("DegtoRad Tests", degToRadTests)
	describe("Test square", squareTest)
}

	function clampTests()
	{
		test("X in the range of min and max" , clampInRange);
		test("X less than min" , clampLessThanMin);
		test("X greater than max" , clampGreaterThanMax);
	}

		function clampInRange()
		{
			expect($3Dmol.clamp(12, 6, 25)).toBe(12);
		}

		function clampLessThanMin()
		{
			expect($3Dmol.clamp(1, 2, 3)).toBe(2);
		}

		function clampGreaterThanMax()
		{
			expect($3Dmol.clamp(100, 25, 75)).toBe(75);
		}

	function degToRadTests()
	{
		test("180 degrees", degToRadPiRadian)
		test("45 degrees", degToRadQuarterPiRadian)
	}

		function degToRadPiRadian()
		{
			expect($3Dmol.degToRad(180)).toBe(Math.PI);
		}
		function degToRadQuarterPiRadian()
		{
			expect($3Dmol.degToRad(45)).toBe(Math.PI / 4);
		}

	function squareTest()
	{
		test("2 squared", twoSquared);
	}
		function twoSquared()
		{
			expect($3Dmol.square(2)).toBe(4);
		}

function quatTests()
{
	describe("Quat Constructor Tests", quatConTests)
	describe("Quat Methods Tests", quatMethodTests)
}

	function quatConTests()
	{
		test("Constructor all args given", quatConstructorBasic)
		test("Constructor no args given", quatConstructorNoArgs)
	}

		function quatConstructorBasic()
		{
			const quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			expect(quat).toEqual({x: 1, y: 2, z: 3, w:4});
		}

		function quatConstructorNoArgs()
		{
			const quat = new $3Dmol.Quaternion();
			expect(quat).toEqual({x: 0, y: 0, z: 0, w:1});
		}

	function quatMethodTests()
	{
		test("Set with 1, 2, 3, 4", quatSet)
		test("Copy, make sure values are equal but references are not", quatCopy)
		test("Conjugate with 1, -2, -3, 4: Signs should flip but 4 shouldn't change", quatConjugate)
		test("Inverse with 1, 2, 3, 4: Flip x-z's signs and multiply by 1/sqrt(30)", quatInverse)
		test("Length with 1, 2, 3, 4: Should be sqrt(30)", quatLength)
		test("Lengthxyz with 1, 2, 3,4: Should be sqrt(14)", quatLengthXYZ)
		test("Normalize, length=0, x,y,z = 0, w =1", quatNormalizeZero)
		test("Normalize, non-zero length: All's values multiplied by length", quatNormalizeNonZero)
		test("Multiply Quaternions: Test the function", quatMultiply)
		test("Multiply Scalar: Should multiply all values by 2", quatMultiplyScalar)
		test("Multiply Quaternions: Should do proper operation", quatMultiplyQuats)
		test("Subtract Quaternions: Should do prop - argProp for each property", quatSubtract)
		test("Clone a Quaternion, Equals props but different references", quatClone)
		test("Set Quaternion from Euler, Converts an euler x,y,z to a Quaternion", quatEulerConversion)
	}
		function quatSet()
		{
			let quat = new $3Dmol.Quaternion();
			quat = quat.set(1, 2, 3, 4);
			expect(quat).toEqual({x: 1, y: 2, z: 3, w:4});
		}
		function quatCopy()
		{
			let q = new $3Dmol.Quaternion();
			q = q.copy(new $3Dmol.Quaternion(1, 2, 3, 4));
			expect(q).toEqual({x: 1, y: 2, z: 3, w:4});		}
		function quatConjugate()
		{
			let quat = new $3Dmol.Quaternion(1, -2, -3, 4);
			quat = quat.conjugate();
			expect(quat).toEqual({x: -1, y: 2, z: 3, w:4});
		}
		function quatInverse()
		{
			let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			quat = quat.inverse();

			//nf is normalizing factor
			const nf = 1 / Math.sqrt(30);
			expect(quat).toEqual({x: -1 * nf, y: -2 * nf, z: -3 * nf, w:4 * nf});

		}
		function quatLength()
		{
			const quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			const ans = quat.length();
			expect(ans).toBe(Math.sqrt(30));
		}
		function quatLengthXYZ()
		{
			const quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			const ans = quat.lengthxyz();
			expect(ans).toBe(Math.sqrt(14));
		}
		function quatNormalizeZero()
		{
			let quat = new $3Dmol.Quaternion(0, 0, 0, 0);
			quat = quat.normalize();
			expect(quat).toEqual({x: 0, y: 0, z:0, w:1});
		}
		function quatNormalizeNonZero()
		{
			let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			quat = quat.normalize();

			//nf is normalizing factor
			const nf = 1 / Math.sqrt(30);
			expect(quat).toEqual({x: nf * 1, y: nf * 2, z:nf * 3, w:nf * 4});
		}
		function quatMultiply()
		{
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let p = new $3Dmol.Quaternion(5, 6, 7, 8);
			q = q.multiply(p);
			expect(q).toEqual({x: 24, y: 48, z: 48, w:-6});

		}
		function quatMultiplyScalar()
		{
			let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
			const scalar = 2;
			quat = quat.multiplyScalar(2);
			expect(quat).toEqual({x: 2, y: 4, z: 6, w:8});
		}
		function quatMultiplyQuats()
		{
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let p = new $3Dmol.Quaternion(5, 6, 7, 8);
			let resultQuat = q.multiplyQuaternions(q, p);
			expect(resultQuat).toEqual({x: 24, y: 48, z: 48, w:-6})
		}
		function quatSubtract()
		{
			let p = new $3Dmol.Quaternion(5, 6, 7, 8);
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			p = p.sub(q);
			expect(p).toEqual({x: 4, y: 4, z: 4, w:4});
		}
		function quatClone()
		{
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let p = q.clone();
			expect(q).toEqual(p);
			expect(q).not.toBe(p); 
		}
		function quatEulerConversion()
		{
			let q = new $3Dmol.Quaternion(1, 2, 3, 4);
			let e = {x: 2, y: 4, z: 6};
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

			expect(q).toEqual({x: newX, y: newY, z: newZ, w: newW});
		}

function vecTests()
{
	describe("2D Vec Tests", vec2Tests)
	describe("3D Vec Tests", vec3Tests)
}

	function vec2Tests()
	{
		test("Vector constructor: 1,2" , vecConstructorBasic)
		test("Vector constructor no args", vecConstructorNoArgs)
		test("Vector set: 1,2", vecSet)
		test("Vector subtraction: <3, 4> - <1, 2>", vecSub)
		test("Vector copy: <1, 2> Same values", vecCopy)
		test("Vector clone: <1, 2> Same values, different reference", vecClone)
	}
		function vecConstructorBasic()
		{
			const vec = new $3Dmol.Vector2(1, 2);
			expect(vec).toEqual({x: 1, y:2});
		}
		function vecConstructorNoArgs()
		{
			const vec = new $3Dmol.Vector2();
			expect(vec).toEqual({x: 0.0, y:0.0})
		}
		function vecSet()
		{
			const vec = new $3Dmol.Vector2();
			vec.set(1, 2);
			expect(vec).toEqual({x: 1, y:2});
		}
		function vecSub()
		{
			const a = new $3Dmol.Vector2(3, 4);
			const b = new $3Dmol.Vector2(1, 2);
			let c = new $3Dmol.Vector2();
			c = c.subVectors(a, b);
			expect(c).toEqual({x: 2, y: 2})
		}
		function vecCopy()
		{
			let a = new $3Dmol.Vector2();
			const b = new $3Dmol.Vector2(1, 2);
			a = a.copy(b);
			expect(a).toEqual({x: 1, y:2});
		}
		function vecClone()
		{
			let a = new $3Dmol.Vector2(1, 2);
			const aRef = a;
			a = a.clone(a);
			expect(a).toEqual({x: 1, y:2});
			expect(a).not.toBe(aRef);
		}

	function vec3Tests()
	{
		test("Vector constructor: 1,2" , vec3ConstructorBasic)
		test("Vector constructor no args", vec3ConstructorNoArgs)
		test("Vector set: 1,2,3", vec3Set)
		test("Vector copy: <1, 2, 3> Same values", vec3Copy)
		test("Vector add: <1, 2, 3> + <1, 1, 1>", vecAdd)
		test("Add vectors: <1, 2, 3> + <1, 1, 1>", addVects)
		test("Multiply vectors elementwise: <1, 2, 3> * <4, 5, 6>", vecMul)
		test("Subtract vectors: <4, 5, 6> - <1, 2, 3> ", subVects)
		test("Subtract vectors: <4, 5, 6> - <1, 2, 3> ", subVects2)
		test("Scale vector: <1, 2, 3> * 2", vecMulScalar)
		test("Shrink vector: <1, 2, 3> / 0", vecDivScalarZero)
		test("Shrink vector: <2, 4, 6> / 2", vecDivScalar)
		test("Vector max: <1, 2, 3>, <4, 5, 6>: Return <4, 5, 6>", vecMax)
		test("Vector min: <1, 2, 3>, <4, 5, 6>: Return <1, 2, 3>", vecMin)
		test("Vector distance to from: <1, 2, 3> to <4, 5, 6>, return sqrt(27)", vecDist)
		test("Vector distance to squared from: <1, 2, 3> to <4, 5, 6>, Return 27", vecDistSq)
		test("Apply matrix3 [1 2 3 4 5 6 7 8 9] to vector <1, 2, 3> Return <14, 32, 50>", appMat3)
		test("Apply matrix4 [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] to vector <1, 2, 3> Return <7, 7, 7>", appMat4)
		test("Apply projection [1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0] to vector <1, 2, 3> Return <INF, INF, INF>", appProj)
		test("Apply quaternion 1, 1, 1, 1 to vector <1, 2, 3> Return <5, 6, -5> ", appQuat)
		test("Vector negation: <1, 2, 3>, Return <-1, -2, -3>", vecNeg)
		test("Vector dot product: <1, 2, 3> . <4, 5, 6>, Return 32", vecDot)
		test("Vector length: <3, 4, 12>, Return 13", vecLength)
		test("Vector length squared: <1, 2, 3>, Return 14", vecLengthSq)
		test("Vector normalize: <3, 4, 12>, Return <3/13, 4/13, 12/13>", vecNormalize)
		test("Vector cross product: <1, 2, 3> x <4, 5, 6>. Return <-3, 6, -3>", vecCrossProd)
		test("Vector cross product new obj: <1, 2, 3> x <4, 5, 6>. Return <-3, 6, -3>", vecCrossProdNew)
		test("Vector get position from matrix: [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15], Return <3, 7, 11>", vecGetMatPos)
		test("Vector setEulerFromRotationMatrix bogus order Call console.error", vecSetEulerBogusOrder)
		test("Vector setEulerFromRotationMatrix element 13 < 0.9999999, Return <0, 0, PI/2>", vec13LessThanOne)
		test("Vector setEulerFromRotationMatrix element 13 > 0.9999999, Return <PI/2, 1, 0>", vec13GThanEqOne)
		test("Vector rotateAboutVector: <1, 2, 3>, ang = PI/2 rads, Return same vector with all divided by sqrt(14)", vecRotate)
		test("Vector set position from matrix: [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15], Return <3, 7, 11>", vecSetMatPos)
		test("Vector transfromDirection from 0-15 4x4 matrix, Return <1, 2, 7> all divided by sqrt(66)", vecTransformDirection)
		test("Vector clone: <1, 2, 3> Same values, different reference", vec3Clone)

	}
		function vec3ConstructorBasic()
		{
			const vec = new $3Dmol.Vector3(1, 2, 3);
			expect(vec).toEqual({x: 1, y:2, z: 3});
		}
		function vec3ConstructorNoArgs()
		{
			const vec = new $3Dmol.Vector3();
			expect(vec).toEqual({x: 0.0, y:0.0, z:0.0})
		}
		function vec3Set()
		{
			const vec = new $3Dmol.Vector3();
			vec.set(1, 2, 3);
			expect(vec).toEqual({x: 1, y:2, z:3});
		}
		function vec3Copy()
		{
			let a = new $3Dmol.Vector3();
			const b = new $3Dmol.Vector3(1, 2, 3);
			a = a.copy(b);
			expect(a).toEqual({x: 1, y:2, z: 3});
		}
		function vecAdd()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(1, 1, 1);
			a = a.add(b);
			expect(a).toEqual({x: 2, y:3, z: 4});
		}
		function addVects()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(1, 1, 1);
			let c = new $3Dmol.Vector3();
			c = c.addVectors(a, b);
			expect(c).toEqual({x: 2, y:3, z: 4});
		}
		function vecMul()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(4, 5, 6);
			let c = new $3Dmol.Vector3();
			c = c.multiplyVectors(a, b);
			expect(c).toEqual({x: 4, y:10, z: 18});
		}
		function subVects()
		{
			const a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			b = b.sub(a);
			expect(b).toEqual({x: 3, y:3, z: 3});
		}
		function subVects2()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			const b = new $3Dmol.Vector3(4, 5, 6);
			let c = new $3Dmol.Vector3();
			c = c.subVectors(b, a);
			expect(c).toEqual({x: 3, y:3, z:3});
		}
		function vecMulScalar()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			a = a.multiplyScalar(2);
			expect(a).toEqual({x: 2, y:4, z:6});
		}
		function vecDivScalarZero()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			a = a.divideScalar(0);
			expect(a).toEqual({x: 0, y:0, z:0});
		}
		function vecDivScalar()
		{
			let a = new $3Dmol.Vector3(2, 4, 6);
			a = a.divideScalar(2);
			expect(a).toEqual({x: 1, y:2, z:3});
		}
		function vecMax()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			a = a.max(b);
			expect(a).toEqual({x: 4, y:5, z:6});
		}
		function vecMin()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			a = a.min(b);
			expect(a).toEqual({x: 1, y:2, z:3});
		}
		function vecDist()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			const dist = a.distanceTo(b);
			expect(dist).toBe(Math.sqrt(27));
		}
		function vecDistSq()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			const distSquared = a.distanceToSquared(b);
			expect(distSquared).toBe(27);
		}
		function appMat3()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix3(...range(1, 9));
			a = a.applyMatrix3(m);
			expect(a).toEqual({x: 14, y:32, z:50});
		}
		function appMat4()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(...nDigits(1, 16));
			a = a.applyMatrix4(m);
			expect(a).toEqual({x: 7, y:7, z:7});
		}
		function appProj()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0);
			a = a.applyProjection(m);
			expect(a).toEqual({x: Infinity, y: Infinity, z: Infinity});
		}
		function appQuat()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let q = new $3Dmol.Quaternion(1, 1, 1, 1);
			a = a.applyQuaternion(q);
			expect(a).toEqual({x: 5, y: 6, z: -5});
		}
		function vecNeg()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			a = a.negate();
			expect(a).toEqual({x: -1, y: -2, z: -3});
		}
		function vecDot()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			const c = a.dot(b);
			expect(c).toBe(32);
		}
		function vecLength()
		{
			let a = new $3Dmol.Vector3(3, 4, 12);
			const c = a.length();
			expect(c).toBe(13);
		}
		function vecLengthSq()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			const c = a.lengthSq();
			expect(c).toBe(14);
		}
		function vecNormalize()
		{
			let a = new $3Dmol.Vector3(3, 4, 12);
			a = a.normalize();
			expect(a).toEqual({x: 3/13, y: 4/13, z: 12/13})
		}
		function vecCrossProd()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			a = a.cross(b);
			expect(a).toEqual({x: -3, y: 6, z: -3});
		}
		function vecCrossProdNew()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let b = new $3Dmol.Vector3(4, 5, 6);
			let c = new $3Dmol.Vector3();
			c  = c.crossVectors(a, b);
			expect(c).toEqual({x: -3, y: 6, z: -3})
		}
		function vecGetMatPos()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let matArr = range(0, 16)
			let m = new $3Dmol.Matrix4(...matArr);
			a = a.getPositionFromMatrix(m);
			expect(a).toEqual({x: 3, y: 7, z: 11});
		}
		function vecSetMatPos()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let matArr = range(0, 16)
			let m = new $3Dmol.Matrix4(...matArr);
			a = a.setFromMatrixPosition(m);
			expect(a).toEqual({x: 3, y: 7, z: 11});
		}
		function vecSetEulerBogusOrder()
		{
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
		}
		function vec13LessThanOne()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(0, -1, 0, 24, 1, 0, 0, 24, 0, 0, 1, 24, 24, 24, 24, 24)
			a = a.setEulerFromRotationMatrix(m, "XYZ");
			expect(a).toEqual({x: -0, y: 0, z: Math.PI / 2});
		}
		function vec13GThanEqOne()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let m = new $3Dmol.Matrix4(0, 0, 1, 24, 1, 0, 0, 24, 0, 1, 0, 24, 24, 24, 24, 24)
			a = a.setEulerFromRotationMatrix(m, "XYZ");
			expect(a).toEqual({x: Math.PI / 2, y: Math.PI / 2, z: 0});
		}
		function vecRotate()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			const ang = Math.PI / 2;
			a = a.rotateAboutVector(a, ang);
			const sqrt14 = Math.sqrt(14);
			expect(a.x).toBeCloseTo(1 / sqrt14);
			expect(a.y).toBeCloseTo(2 / sqrt14);
			expect(a.z).toBeCloseTo(3 / sqrt14);
		}
		function vecTransformDirection()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			let matArr = range(0, 16);
			let m = new $3Dmol.Matrix4(...matArr);
			a = a.transformDirection(m);
			const sqrt66 = Math.sqrt(66);
			expect(a.x).toBeCloseTo(1 / sqrt66);
			expect(a.y).toBeCloseTo(4 / sqrt66);
			expect(a.z).toBeCloseTo(7 / sqrt66);
		}
		function vec3Clone()
		{
			let a = new $3Dmol.Vector3(1, 2, 3);
			const aRef = a;
			a = a.clone(a);
			expect(a).toEqual({x: 1, y:2, z: 3});
			expect(a).not.toBe(aRef);
		}

function matTests()
{
	describe("3x3 Matrix Tests", mat3Tests)
	describe("4x4 Matrix Tests", mat4Tests)
}
	function mat3Tests()
	{
		test("Matrix 3 Constructor, [1 2 3 4 5 6 7 8 9]", mat3Constructor)
		test("Matrix 3 Set, [1 1 1 1 1 1 1 1 1]", mat3Set)
		test("Matrix 3 Identity, [1 0 0 0 1 0 0 1]", mat3Id)
		test("Matrix 3 Copy, [0 0 0 0 0 0 0 0 0]", mat3Cop)
		test("Matrix 3 Multiply Scalar, [1 1 1 1 1 1 1 1 1] by 2", mat3MultScal)
		test("Matrix 3 getInverse3, [1, 2, 3, 4, 5, 6, 7, 8, 24]", mat3Inv3)
		test("Matrix 3 getInverse Inputs: Random stuff that works for happy path Throw on invertible = false", mat3Inv)
		test("Matrix 3 getInverse Inputs: 3x3 of all 1's and 4x4 of all 1's, throw on invertible = false", mat3InvWarn)
		test("Matrix 3 getInverse Inputs: 3x3 of all 1's and 4x4 of all 1's, throw on invertible = true", mat3InvThrow)
		test("Matrix 3 Determinant, [1, 2, 3, 4, 5, 6, 7, 8, 24], Return -45", mat3Det)
		test("Matrix 3 Get Matrix 4, [1 1 1 1 1 1 1 1 1], Throw in 0s to pad out", mat3GetMat4)
		test("Matrix 3 Transpose, [1 2 3 4 5 6 7 8 9], Return [1 4 7 2 5 8 3 6 9]", mat3Transpose)
		test("Matrix 3 Clone, [1 2 3 4 5 6 7 8 9], Return new Matrix with exact same values", mat3Clone)
		test("Matrix 3 conversionMatrix3, some wacky stuff happens, Return [1 0 0 sqrt(2) sqrt(2) 0 1.5 (1.5 * (sqrt(6) -1)) really long arg]", convertMatrix3)
	}
		function mat3Constructor()
		{
			let mat = new $3Dmol.Matrix3(...range(1, 9));
			expect([...mat.elements]).toEqual([1, 4, 7, 2, 5, 8, 3, 6, 9]);
		}
		function mat3Set()
		{
			let mat = new $3Dmol.Matrix3();
			mat = mat.set(...nDigits(1, 9))
			expect([...mat.elements]).toEqual(nDigits(1, 9));
		}
		function mat3Id()
		{
			let mat = new $3Dmol.Matrix3();
			mat = mat.identity();
			expect([...mat.elements]).toEqual([1, 0, 0, 0, 1, 0, 0, 0, 1]);
		}
		function mat3Cop()
		{
			let mat = new $3Dmol.Matrix3();
			let elems = nDigits(0, 9);
			let from = new $3Dmol.Matrix3(...elems);
			mat.copy(from);
			expect([...mat.elements]).toEqual(elems);
		}
		function mat3MultScal()
		{
			let mat = new $3Dmol.Matrix3(...nDigits(1, 9));
			mat = mat.multiplyScalar(2);
			expect([...mat.elements]).toEqual(nDigits(2, 9));
		}
		function mat3Inv3()
		{
			let mat = new $3Dmol.Matrix3();
			let given = new $3Dmol.Matrix3(1, 2, 3, 4, 5, 6, 7, 8, 24);
			mat = mat.getInverse3(given);
			expect(mat.elements[0]).toBeCloseTo(-8/5);
			expect(mat.elements[1]).toBeCloseTo(1.2);
			expect(mat.elements[2]).toBeCloseTo(1/15);
			expect(mat.elements[3]).toBeCloseTo(8/15);
			expect(mat.elements[4]).toBeCloseTo(-1/15);
			expect(mat.elements[5]).toBeCloseTo(-2/15);
			expect(mat.elements[6]).toBeCloseTo(1/15);
			expect(mat.elements[7]).toBeCloseTo(-2/15);
			expect(mat.elements[8]).toBeCloseTo(1/15);
		}
		function mat3Inv()
		{
			let mat3 = new $3Dmol.Matrix3(9, 5, 7, 3, 2, 9, 5, 7, 2);
			let mat4 = new $3Dmol.Matrix4(7, 10, 6, 8, 9, 2, 4, 6, 2, 7, 4, 1, 3, 4, 10, 8);
			let toI = false;
			mat3 = mat3.getInverse(mat4, toI);
			let expectedElems = [-20, -28, 59, 2, 16, -29, 28, 27, -76]
			expectedElems = expectedElems.map(x => x / -66);
			for(let i = 0; i < expectedElems.length; i++)
			{
				expect(mat3.elements[i]).toBeCloseTo(expectedElems[i], 1);
			}
		}
		function mat3InvWarn()
		{
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


		}
		function mat3InvThrow()
		{
			let mat3 = new $3Dmol.Matrix3(...nDigits(1, 9));
			let mat4 = new $3Dmol.Matrix4(...nDigits(1, 16));
			let toI = true;
			expect(() => mat3 = mat3.getInverse(mat4, toI)).toThrow("Matrix3.getInverse(): can't invert matrix, determinant is 0");
		}
		function mat3Det()
		{
			let mat = new $3Dmol.Matrix3(1, 2, 3, 4, 5, 6, 7, 8, 24);
			const det = mat.getDeterminant();
			expect(det).toBe(-45);
		}
		function mat3GetMat4()
		{
			let mat = new $3Dmol.Matrix3(...nDigits(1, 9));
			let mat4 = mat.getMatrix4();
			expect([...mat4.elements]).toEqual([1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1])
		}
		function mat3Transpose()
		{
			/*These look not different, but you write to the object in row-major order
			 * and read from the object in column-major order
			 */ 
			let mat = new $3Dmol.Matrix3(...range(1, 9));
			mat = mat.transpose();
			expect([...mat.elements]).toEqual(range(1, 9));
		}
		function mat3Clone()
		{
			let mat = new $3Dmol.Matrix3(...range(1, 9));
			let matNew = mat.clone();
			expect(mat).not.toBe(matNew);
			expect([...mat.elements]).toEqual([...matNew.elements]);
		}
		function convertMatrix3()
		{
			let mat3 = $3Dmol.conversionMatrix3(1, 2, 3, 30, 60, 45);
			let lastElem = (3 * Math.sqrt(2) * Math.sqrt(Math.sqrt(6) - 2)) / 2;
			let expected = [1, 0, 0, Math.sqrt(2), Math.sqrt(2), 0, 1.5, (1.5 * (Math.sqrt(6) - 1)), lastElem];
			for(let i = 0; i < mat3.elements.length; i++)
			{
				expect(mat3.elements[i]).toBeCloseTo(expected[i], 1);
			}
		}

function mat4Tests()
{
	test("Matrix 4 Constructor, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]", mat4Constructor)
	test("Matrix 4 Set, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]", mat4Set)
	test("Matrix 4 Identity, [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1]", mat4Identity)
	test("Matrix 4 Copy, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]", mat4Copy)
	test("Matrix 4 Matrix3 from Top left, [1 - 16] -> [1 5 9 2 6 10 3 7 11] ", mat4Mat3)
	test("Matrix 4 setRotationFromEuler bogus order, call console.error", mat4SetRotFromEulerBogus)
	test("Matrix 4 setRotationFromEuler with [1 - 16] and vec: <1, 2, 3>", mat4SetRotFromEuler)
	test("Matrix 4 setRotationFromQuaternion with [1 - 16] and quat: (1, 2, 3, 4)", mat4SetRotFromQuat)
	//test("Matrix 4 lookAt, with [1-16], Return ", mat4lookAt)
	test("Matrix 4 multplication wtih [1 - 16] and [1 - 16', Return [90 100 110 120 202 228 254 280 314 356 398 440 426 484 542 600]", mat4Mult)
	test("Matrix 4 multiply scalar all 1's with 2, Return all 2's", mat4MultiplyScalar)
	test("Matrix 4 makeTranslation [1-16]", mat4MakeTranslation)
	test("Matrix 4 snap [1-16] with all elements + 0.00001, Return [1-16]", mat4Snap)
	test("Matrix 4 transpose [1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16], Return [1-16]", mat4Transpose)
	//test("Matrix 4 getPostion, not sure how to test closures", mat4getPosition)
	test("Matrix 4 setPosition, [1 - 16] with <1, 2, 3>, Return [1 2 3 1 5 6 7 2 9 10 11 3 13 14 15 16]", mat4SetPostion)
	test("Matrix 4 translate, [1 - 16] with <1, 2, 3>, Return [1 2 3 5 5 6 7 10 9 10 11 15 13 14 15 16]", mat4Translate)
	test("Matrix 4 getInverse, random inputs that make the happy path pass, throw on invertible = false", mat4Inv)
	test("Matrix 4 getInverse, two [1-16] matrixes with throw on invertible = false", mat4InvWarn)
	test("Matrix 4 getInverse, two [1-16] matrixes with throw on invertible = true", mat4InvThrow)
	test("Matrix 4 isReflected, upperLeft corner  of 4x4 matrix is [1 2 3 4 5 6 7 8 10], Return true", mat4IsReflected)
	test("Matrix 4 isReflected, upperLeft corner  of 4x4 matrix is [1 2 3 4 5 6 7 8 10], Return false", mat4IsNotReflected)
	test("Matrix 4 scale, scale a matrix of all 1's with vector <2, 2, 2>, Return all 2's expect last 4 are 1's", mat4Scale)
	test("Matrix 4 getMaxScaleOnAxis, a single [1-16] matrix, Return sqrt(179)", mat4gmscoa)
	test("Matrix 4 makeFrustum, a single [1-16] matrix, args in order: 1, 2, 3, 4, 5, 6, Return [10 0 0 0 0 10 0 0 3 7 -11 -1 0 0 -60 0]", mat4MakeFrustum)
	test("Matrix 4 makePerspective, a single [1-16] matrix, args in order: 90, 1, 2, 3, Return [NaN 0 0 0 0 NaN 0 0 -1.0224719047546387 5 NaN -1 0 0 NaN 0]", mat4MakePerspective)
	test("Matrix 4 makeOrthographic, a single [1-16] matrix, args in order: 1, 2, 3, 4, 5, 6, Return: [2 0 0 0 0 -2 0 0 0 0 -2 0 -3 7 -11 1]", mat4MakeOrthographic)
	test("Matrix 4 isEqual, two all 1's 4x4 matricies, Return true", mat4AreEqual)
	test("Matrix 4 isEqual, one all 1's 4x4 matricies and other one is all 2's 4x4 matricies, Return false", mat4AreNotEqual)
	test("Matrix 4 clone, all 1's 4x4 matricies, Return identical matrix with different reference", mat4Clone)
	test("Matrix 4 isIdentity, all 0's expect 1's down the diagonal, Return true", mat4IsIdentity)
	test("Matrix 4 isIdentity, [1-16], Return false", mat4IsNotIdentity)
	test("Matrix 4 isNearlyIdentity, all 0's expect 1's down the diagonal, add 0.00001 to all elements, Return true", mat4IsNearlyIdentity)
}
	function mat4Constructor()
	{
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
	}
	function mat4Set()
	{
		let mat = new $3Dmol.Matrix4();
		mat = mat.set(...range(1, 16));
		expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
	}
	function mat4Identity()
	{
		let mat = new $3Dmol.Matrix4();
		mat = mat.identity();
		expect([...mat.elements]).toEqual([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
	}
	function mat4Copy()
	{
		let mat = new $3Dmol.Matrix4();
		mat = mat.copy(new $3Dmol.Matrix4(...range(1, 16)));
		expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
	}
	function mat4Mat3()
	{
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		let mat3 = mat.matrix3FromTopLeft();
		expect([...mat3.elements]).toEqual([1, 5, 9, 2, 6, 10, 3, 7, 11]);
	}
	function mat4SetRotFromEulerBogus()
	{
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
	}
	function mat4SetRotFromEuler()
	{
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		let vec = new $3Dmol.Vector3(Math.PI/2, Math.PI /2, Math.PI /2);
		mat = mat.setRotationFromEuler(vec, "XYZ");

		//Pretreat some math weirdness
		mat.elements[0] = Math.round(mat.elements[0]);
		mat.elements[1] = Math.round(mat.elements[1]);
		mat.elements[4] = Math.abs(Math.round(mat.elements[4]));
		mat.elements[6] = Math.round(mat.elements[6]);
		mat.elements[9] = Math.abs(Math.round(mat.elements[9]));
		mat.elements[10] = Math.round(mat.elements[10]);

		expect([...mat.elements]).toEqual([0, 0, 1, 13, 0, -1, 0, 14, 1, 0, 0, 15, 4, 8, 12, 16]);
	}
	function mat4SetRotFromQuat()
	{
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		let quat = new $3Dmol.Quaternion(1, 2, 3, 4);
		mat = mat.setRotationFromQuaternion(quat);
		expect([...mat.elements]).toEqual([-25, 28, -10, 13, -20, -19, 20, 14, 22, 4, -9, 15, 4, 8, 12, 16]);
	}
	
	 /*function mat4lookAt()
	{
		
	}
	*/
	
	function mat4Mult()
	{
		let mat = new $3Dmol.Matrix4();
		let a = new $3Dmol.Matrix4(...range(1, 16));
		let b = new $3Dmol.Matrix4(...range(1, 16));
		mat = mat.multiplyMatrices(a, b);
		expect([...mat.elements]).toEqual([90, 202, 314, 426, 100, 228, 356, 484, 110, 254, 398, 542, 120, 280, 440, 600]);
	}
	function mat4MultiplyScalar()
	{
		let mat = new $3Dmol.Matrix4(nDigits(1, 16));
		const scalar = 2;
		mat = mat.multiplyScalar(scalar);
		expect([...mat.elements]).toEqual(nDigits(2, 16));
	}
	function mat4MakeTranslation()
	{
		let mat = new $3Dmol.Matrix4();
		mat = mat.makeTranslation(1, 2, 3);
		expect([...mat.elements]).toEqual([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 2, 3, 1]);
	}
	function mat4Snap()
	{
		let elems = range(1, 16);
		elems = elems.map(x => x + 0.00001);
		let mat = new $3Dmol.Matrix4(elems);
		mat = mat.snap(4);
		expect([...mat.elements]).toEqual(range(1, 16));
	}
	function mat4Transpose()
	{
		let mat = new $3Dmol.Matrix4([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]);
		mat = mat.transpose();
		expect([...mat.elements]).toEqual(range(1, 16));
	}
	/*function mat4getPosition()
	{

	}
	*/
	function mat4SetPostion()
	{
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		let vec = new $3Dmol.Vector3(1, 2, 3);
		mat = mat.setPosition(vec);
		expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 1, 2, 3, 16]);
	}
	function mat4Translate()
	{
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		let vec = new $3Dmol.Vector3(1, 2, 3);
		mat = mat.translate(vec);
		expect([...mat.elements]).toEqual([1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 5, 10, 15, 16]);
	}
	function mat4Inv()
	{
		let matA = new $3Dmol.Matrix4(7, 10, 6, 8, 9, 2, 4, 6, 2, 7, 4, 1, 3, 4, 10, 8);
		let matB = new $3Dmol.Matrix4(7, 10, 6, 8, 9, 2, 4, 6, 2, 7, 4, 1, 3, 4, 10, 8);
		matA = matA.getInverse(matB, false);
		let elems = [160, -194, 364, -418, -340, 140, -130, 220, -200, -104, -356, 572, 120, 102, -222, -66];
		elems = elems.map(x => x / -1980)
		for(let i = 0; i < elems.length; i++)
		{
			expect(matA.elements[i]).toBeCloseTo(elems[i], 1);
		}
	}
	function mat4InvWarn()
	{
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
	}
	function mat4InvThrow()
	{
		let matA = new $3Dmol.Matrix4(...range(1, 16));
		let matB = new $3Dmol.Matrix4(...range(1, 16));
		expect(() => matA = matA.getInverse(matB, true)).toThrow("Matrix4.getInverse(): can't invert matrix, determinant is 0");
	}
	function mat4IsReflected()
	{
		let mat = new $3Dmol.Matrix4(1, 2, 3, 0, 4, 5, 6, 0, 7, 8, 10, 0, 0, 0, 0, 0);
		const reflected = mat.isReflected();
		expect(reflected).toBe(true);
	}
	function mat4IsNotReflected()
	{
		let mat = new $3Dmol.Matrix4(1, 2, 3, 0, 4, 5, 6, 0, 7, 8, -9, 0, 0, 0, 0, 0);
		const reflected = mat.isReflected();
		expect(reflected).toBe(false);
	}
	function mat4Scale()
	{
		let mat4 = new $3Dmol.Matrix4(...nDigits(1, 16));
		let v = new $3Dmol.Vector3(2, 2, 2);
		mat4 = mat4.scale(v);
		let expected = nDigits(2, 16);
		for(let i = 12; i < 16; i++)
		{
			expected[i] = 1;
		}
		expect([...mat4.elements]).toEqual(expected);
	}
	function mat4gmscoa()
	{
		let mat4 = new $3Dmol.Matrix4(...range(1, 16));
		const ans = mat4.getMaxScaleOnAxis();
		expect(ans).toEqual(Math.sqrt(179));
	}
	function mat4MakeFrustum()
	{
		let mat4 = new $3Dmol.Matrix4(...range(1, 16));
		mat4 = mat4.makeFrustum(1, 2, 3, 4, 5, 6);
		expect([...mat4.elements]).toEqual([10, 0, 0, 0, 0, 10, 0, 0, 3, 7, -11, -1, 0, 0, -60, 0]);
	}
	function mat4MakePerspective()
	{
		let mat4 = new $3Dmol.Matrix4(...range(1, 16));
		mat4 = mat4.makePerspective(90, 1, 2, 3);
		expect([...mat4.elements]).toEqual([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -5, -1, 0, 0, -12, 0]);
	}
	function mat4MakeOrthographic()
	{
		let mat4 = new $3Dmol.Matrix4(...range(1, 16));
		mat4 = mat4.makeOrthographic(1, 2, 3, 4, 5, 6);
		expect([...mat4.elements]).toEqual([2, 0, 0, 0, 0, -2, 0, 0, 0, 0, -2, 0, -3, 7, -11, 1]);
	}
	function mat4AreEqual()
	{
		let a = new $3Dmol.Matrix4(...nDigits(1, 16));
		let b = new $3Dmol.Matrix4(...nDigits(1, 16));

		let equal = a.isEqual(b);
		expect(equal).toBe(true);
	}
	function mat4AreNotEqual()
	{
		let a = new $3Dmol.Matrix4(...nDigits(1, 16));
		let b = new $3Dmol.Matrix4(...nDigits(2, 16));

		let equal = a.isEqual(b);
		expect(equal).toBe(false);
	}
	function mat4Clone()
	{
		let a = new $3Dmol.Matrix4(...nDigits(1, 16));
		let b = a.clone();
		expect([...a.elements]).toEqual([...b.elements]);
		expect(a).not.toBe(b);
	}
	function mat4IsIdentity()
	{
		let mat = new $3Dmol.Matrix4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
		let ident = mat.isIdentity();
		expect(ident).toEqual(true);
	}
	function mat4IsNotIdentity()
	{
		let mat = new $3Dmol.Matrix4(...range(1, 16));
		let ident = mat.isIdentity();
		expect(ident).toEqual(false);
	}
	function mat4IsNearlyIdentity()
	{
		let elems = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
		elems = elems.map(x => x + 0.00001);
		let mat = new $3Dmol.Matrix4(...elems);
		let ident = mat.isNearlyIdentity(4);
		expect(ident).toEqual(true);
	}

function rayTests()
{
	test("Ray constructor, origin: <1, 2, 3>, direction: <4, 5, 6>", rayCon)
	test("Ray set, origin: <1, 2, 3>, direction: <4, 5, 6>", raySet)
	test("Ray copy, origin: <1, 2, 3> direction: <4, 5, 6>", rayCopy)
	test("Ray at, scalar: 24, origin: <1, 2, 3>, direction: <4, 5, 6>, No optional target, Return <97, 122, 147>", rayAt)
	test("Ray at, scalar: 24, origin: <1, 2, 3>, direction: <4, 5, 6>, optional target <1, 1, 1>, Return <97, 122, 147>", rayAtOpt)
	//Recast: Closure, not sure how to test
	test("Ray closestPointToPoint, origin: <1, 2, 3>, direction: <4, 5, 6>, point: <1, 1, 1>, Return <-67, -83, -99>", rayCPTP)
	//distanceToPoint: Closure, not sure how to test
	test("Ray isIntersectionCylinder, Empty function, Return undefined", rayiic)
	test("Ray isIntersectionSphere, origin: <1, 2, 3>, direction: <4, 5, 6>, center: <1, 2, 3>, radius: 9, Return true", rayIsIntersectionSphere)
	/*
	* isIntersectionPlane
	* distanceToPlane
	* intersectPlane
	* The plane class no longer exists?
	*/
	test("Ray applyMatrix4, origin: <1, 2, 3>, direction: <4, 5, 6>, Matrix4: [1-16], newOrig: <18, 46, 74> newDir: <32, 92, 152>", rayApplyMatrix4)
	test("Ray equals, origin: <1, 2, 3>, direction: <4, 5, 6> for both, Throw an error because Vector3 doesn't have an equals method", rayEquals)
	test("Ray clone, origin: <1, 2, 3>, direction: <4, 5, 6>, Identical values but different references", rayClone)
}
	function rayCon()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		expect({...ray.origin}).toEqual({x: 1, y: 2, z: 3});
		expect({...ray.direction}).toEqual({x: 4, y: 5, z: 6});
	}
	function raySet()
	{
		let ray = new $3Dmol.Ray();
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		ray = ray.set(origin, direction);
		expect({...ray.origin}).toEqual({x: 1, y: 2, z: 3});
		expect({...ray.direction}).toEqual({x: 4, y: 5, z: 6});
	}
	function rayCopy()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		let newRay = ray.copy(ray);
		expect({...newRay.origin}).toEqual({x: 1, y: 2, z: 3});
		expect({...newRay.direction}).toEqual({x: 4, y: 5, z: 6});
	}
	function rayAt()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		let resultVec = ray.at(24);
		expect(resultVec).toEqual({x: 97, y: 122, z: 147});
	}
	function rayAtOpt()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		let opt = new $3Dmol.Vector3(1, 2, 3);
		let resultVec = ray.at(24, opt);
		expect(resultVec).toEqual({x: 97, y: 122, z: 147});
	}
	function rayCPTP()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		let point = new $3Dmol.Vector3(1, 1, 1);
		let result = ray.closestPointToPoint(point);
		expect({...result}).toEqual({x: -67, y: -83, z: -99})
	}
	function rayiic()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		expect(ray.isIntersectionCylinder()).toEqual(undefined);
	}
	function rayIsIntersectionSphere()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		let sphCent = new $3Dmol.Vector3(1, 2, 3);
		let sph = new $3Dmol.Sphere(sphCent, 9);
		const res = ray.isIntersectionSphere(sph);
		expect(res).toEqual(true);
	}
	function rayApplyMatrix4()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray =  new $3Dmol.Ray(origin, direction);
		let mat4 = new $3Dmol.Matrix4(...range(1, 16));
		ray = ray.applyMatrix4(mat4);
		expect({...ray.origin}).toEqual({x: 18, y: 46, z: 74});
		expect({...ray.direction}).toEqual({x: 32, y: 92, z: 152});
	}
	function rayEquals()
	{
		/*
		* As it stands, Vector3 doesn't have an equals function, so this can't be tested
		*/
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray1 =  new $3Dmol.Ray(origin, direction);
		let ray2 = new $3Dmol.Ray(origin, direction);
		// the error message can change we just test that this function doesn't exist
		// follow up question, why does this function not exist?
		expect(() => ray1.equals(ray2)).toThrow();
	}
	function rayClone()
	{
		let origin = new $3Dmol.Vector3(1, 2, 3);
		let direction = new $3Dmol.Vector3(4, 5, 6);
		let ray1 =  new $3Dmol.Ray(origin, direction);
		let ray2 = new $3Dmol.Ray(origin, direction);
		expect({...ray1.origin}).toEqual({...ray2.origin});
		expect({...ray1.direction}).toEqual({...ray2.direction});
		expect(ray1).not.toBe(ray2);
	}
