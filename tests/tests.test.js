/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../build/3Dmol.js");

describe('Tests Clamp, degToRad and Quaternion', () => {
	// Clamp Tests
	test("Test clamp with x being in the range of min and max", () => {
		expect($3Dmol.clamp(12, 6, 25)).toBe(12);
	});
	test("Test clamp with x being less than min", () => {
		expect($3Dmol.clamp(1, 2, 3)).toBe(2);
	});
	test("Test clamp with x being greater than max", () => {
		expect($3Dmol.clamp(100, 25, 75)).toBe(75);
	});

	// degToRad Tests
	test("Test degToRad on 180 degrees", () => {
		expect($3Dmol.degToRad(180)).toBe(Math.PI);
	});
	test("Test degToRad on 45 degrees", () => {
		expect($3Dmol.degToRad(45)).toBe(Math.PI / 4);
	});

	// Quaternion Tests
	test("Test Quaternion constructor with all args given", () => {
		const quat = new $3Dmol.Quaternion(1, 2, 3, 4);
		expect(quat).toEqual({ x: 1, y: 2, z: 3, w: 4 });
	});
	test("Test Quaternion constructor with no args given", () => {
		const quat = new $3Dmol.Quaternion();
		expect(quat).toEqual({ x: 0, y: 0, z: 0, w: 1 });
	});
});
