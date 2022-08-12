/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");

//A poorman's JS version of Python's range
let range = (start, end) => [...Array.from(Array(end + 1).keys()).slice(start)];
let nDigits = (digit, n) => Array.from(`${digit}`.repeat(n), Number)

describe("Common Tests", commonTests)
describe("RWB Tests", rwbTests)
describe("ROYGB Tests", roygbTests)
describe("Sinebow Tests", sinebowTests)

function commonTests()
{
	test("Gradient normalize value, 0, 100, 50", gradNorm)
	test("Gradient normalize value, 0, 100, -1 Return: 0, 100, 0", gradNormLT)
	test("Gradient normalize value, 0, 100, 101 Return: 0, 100, 100", gradNormGT)
	test("Gradient normalize value, 100, 0, 50", gradNormFlip)
	test("Gradient normalize value, 100, 0, -1 Return: 0, 100, 100", gradNormFlipLT)
	test("Gradient normalize value, 100, 0, 101, Return: 0, 100, 0", gradNormFlipGT)
	test("Gradient get happy path", gradGet)
	test("Gradient get bogus", gradGetBogus)
	test("Gradient get RWB", gradGetPredef)
	test("Gradient get RWB predefine mid, max, and min", gradGetPredefSet)
}
	function gradNorm()
	{
		let normFunc = $3Dmol.Gradient.normalizeValue
		let obj = normFunc(0, 100, 50)
		expect(obj).toEqual({lo: 0, hi: 100, val: 50});
	}
	function gradNormLT()
	{
		let normFunc = $3Dmol.Gradient.normalizeValue
		let obj = normFunc(0, 100, -1)
		expect(obj).toEqual({lo: 0, hi: 100, val: 0});
	}
	function gradNormGT()
	{
		let normFunc = $3Dmol.Gradient.normalizeValue
		let obj = normFunc(0, 100, 101)
		expect(obj).toEqual({lo: 0, hi: 100, val: 100});
	}
	function gradNormFlip()
	{
		let normFunc = $3Dmol.Gradient.normalizeValue
		let obj = normFunc(100, 0, 50)
		expect(obj).toEqual({lo: 0, hi: 100, val: 50});
	}
	function gradNormFlipLT()
	{
		let normFunc = $3Dmol.Gradient.normalizeValue
		let obj = normFunc(100, 0, -1)
		expect(obj).toEqual({lo: 0, hi: 100, val: 100});
	}
	function gradNormFlipGT()
	{
		let normFunc = $3Dmol.Gradient.normalizeValue
		let obj = normFunc(100, 0, 101)
		expect(obj).toEqual({lo: 0, hi: 100, val: 0});
	}
	function gradGet()
	{
		let grad = new $3Dmol.Gradient();
		let obj = $3Dmol.Gradient.getGradient(grad);
		expect(obj).toEqual({});
	}
	function gradGetBogus()
	{
		let grad = 2;
		expect($3Dmol.Gradient.getGradient(grad)).toBe(2);
	}
	function gradGetPredef()
	{
		let grad = {};
		grad.gradient = "rwb";
		grad = $3Dmol.Gradient.getGradient(grad);
		expect(grad).toBeInstanceOf($3Dmol.Gradient.RWB);
	}
	function gradGetPredefSet()
	{
		let grad = {};
		grad.gradient = "rwb";
		grad.min = 0;
		grad.max = 100;
		grad.mid = 50;
		grad = $3Dmol.Gradient.getGradient(grad);
		expect(grad).toBeInstanceOf($3Dmol.Gradient.RWB);
	}
function rwbTests()
{
	test("RWB no args", rwbNargs)
	test("RWB Constructor with range", rwbRange)
	test("RWB vth normal", rwbVTHNorm)
	test("RWB vth min max mid args", rwbVTHThreeArgs)
	test("RWB vth min max mid args no range", rwbVTHThreeArgsNoRange)
	test("RWB vth min max mid args, val > mid", rwbVTHGT)
	test("RWB vth min max mid args, val == mid", rwbVTHEQ)
}
	function rwbNargs()
	{
		let grad = new $3Dmol.Gradient.RWB();
		expect(grad.range()).toEqual(null);
	}
	function rwbRange()
	{
		let grad = new $3Dmol.Gradient.RWB([0, 100]);
		expect(grad.range()).toEqual([0, 100]);
	}
	function rwbVTHNorm()
	{
		let grad = new $3Dmol.Gradient.RWB([0, 100]);
		let hexVal = grad.valueToHex(10, grad.range());
		expect(hexVal).toEqual(0xff7272);
	}
	function rwbVTHThreeArgs()
	{
		let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		let hexVal = grad.valueToHex(10, [0, 50, 100]);
		expect(hexVal).toEqual(0xff5050);
	}
	function rwbVTHThreeArgsNoRange()
	{
		let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		let hexVal = grad.valueToHex(10);
		expect(hexVal).toEqual(0xff7272);
	}
	function rwbVTHGT()
	{
		let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		let hexVal = grad.valueToHex(60, [0, 100, 50]);
		expect(hexVal).toEqual(0xe4e4ff);
	}
	function rwbVTHEQ()
	{
		let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		let hexVal = grad.valueToHex(50, [0, 100, 50]);
		expect(hexVal).toEqual(0xffffff);
	}

function roygbTests()
{
	test("ROYGB no args", roygbNargs)
	test("ROYGB Constructor with range", roygbRange)
	test("ROYGB vth normal", roygbVTHNorm)
	test("ROYGB vth min max mid args", roygbVTHThreeArgs)
	test("ROYGB vth min max mid args no range", roygbVTHThreeArgsNoRange)
	test("ROYGB vth min max mid args, val < q3", roygbVTHLT)
	test("ROYGB vth min max mid args, val < mid", roygbVTHLTMin)
	test("ROYGB vth min max mid args, val > max", roygbVTHGTMax)
}
	function roygbNargs()
	{
		let grad = new $3Dmol.Gradient.ROYGB();
		expect(grad.range()).toEqual(null);
	}
	function roygbRange()
	{
		let grad = new $3Dmol.Gradient.ROYGB([0, 100]);
		expect(grad.range()).toEqual([0, 100]);
	}
	function roygbVTHNorm()
	{
		let grad = new $3Dmol.Gradient.ROYGB([0, 100]);
		let hexVal = grad.valueToHex(10, grad.range());
		expect(hexVal).toEqual(0xffa100);
	}
	function roygbVTHThreeArgs()
	{
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(10, [0, 50, 100]);
		expect(hexVal).toEqual(0xffe400);
	}
	function roygbVTHThreeArgsNoRange()
	{
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(10);
		expect(hexVal).toEqual(0xffa100);
	}
	function roygbVTHLT()
	{
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(60, [0, 100, 50]);
		expect(hexVal).toEqual(0xffa1);
	}
	function roygbVTHLTMin()
	{
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(49, [0, 100, 50]);
		expect(hexVal).toEqual(0x33ff00);
	}
	function roygbVTHGTMax()
	{
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(101, [0, 100, 50]);
		expect(hexVal).toEqual(0xff);
	}

function sinebowTests()
{
	test("Sinebow no args", sinebowNargs)
	test("Sinebow Constructor with range", sinebowRange)
	test("Sinebow vth normal", sinebowVTHNorm)
	test("Sinebow vth min max mid args no range", sinebowVTHThreeArgsNoRange)
}
	function sinebowNargs()
	{
		let grad = new $3Dmol.Gradient.Sinebow();
		expect(grad.range()).toEqual(null);
	}
	function sinebowRange()
	{
		let grad = new $3Dmol.Gradient.Sinebow([0, 100]);
		expect(grad.range()).toEqual([0, 100]);
	}
	function sinebowVTHNorm()
	{
		let grad = new $3Dmol.Gradient.Sinebow([100, 0]);
		let hexVal = grad.valueToHex(10, grad.range());
		expect(hexVal).toEqual(0x7f11ed);
	}
	function sinebowVTHThreeArgsNoRange()
	{
		let grad = new $3Dmol.Gradient.Sinebow(0, 100, 50);
		expect(grad).toBeDefined();
		let hexVal = grad.valueToHex(10);
		expect(hexVal).toEqual(0xed7f11);
	}