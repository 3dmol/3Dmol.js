/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../build/3Dmol.js");

test("Test clamp with x being in the range of min and max" , clampInRange);
test("Test clamp with x being less than min" , clampLessThanMin);
test("Test clamp with x being greater than max" , clampGreaterThanMax);
test("Test degToRad on 180 degrees", degToRadPiRadian)
test("Test degToRad on 45 degrees", degToRadQuarterPiRadian)
test("Test Quaternion constructor with all args given", quatConstructorBasic)
test("Test Quaternion constructor with no args given", quatConstructorNoArgs)

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
function degToRadPiRadian()
{
	expect($3Dmol.degToRad(180)).toBe(Math.PI);
}
function degToRadQuarterPiRadian()
{
	expect($3Dmol.degToRad(45)).toBe(Math.PI / 4);
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
