/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function () { };
let $3Dmol = require("3dmol");

// Common Tests
describe("Common Tests", () => {
	test("Gradient normalize value, 0, 100, 50", () => {
		let obj = $3Dmol.Gradient.normalizeValue(0, 100, 50);
		expect(obj).toEqual({ lo: 0, hi: 100, val: 50 });
	});

	test("Gradient normalize value, 0, 100, -1 Return: 0, 100, 0", () => {
		let obj = $3Dmol.Gradient.normalizeValue(0, 100, -1);
		expect(obj).toEqual({ lo: 0, hi: 100, val: 0 });
	});

	test("Gradient normalize value, 0, 100, 101 Return: 0, 100, 100", () => {
		let obj = $3Dmol.Gradient.normalizeValue(0, 100, 101);
		expect(obj).toEqual({ lo: 0, hi: 100, val: 100 });
	});

	test("Gradient normalize value, 100, 0, 50", () => {
		let obj = $3Dmol.Gradient.normalizeValue(100, 0, 50);
		expect(obj).toEqual({ lo: 0, hi: 100, val: 50 });
	});

	test("Gradient normalize value, 100, 0, -1 Return: 0, 100, 100", () => {
		let obj = $3Dmol.Gradient.normalizeValue(100, 0, -1);
		expect(obj).toEqual({ lo: 0, hi: 100, val: 100 });
	});

	test("Gradient normalize value, 100, 0, 101, Return: 0, 100, 0", () => {
		let obj = $3Dmol.Gradient.normalizeValue(100, 0, 101);
		expect(obj).toEqual({ lo: 0, hi: 100, val: 0 });
	});

	test("Gradient get happy path", () => {
		let grad = new $3Dmol.Gradient();
		let obj = $3Dmol.Gradient.getGradient(grad);
		expect(obj).toEqual({});
	});

	test("Gradient get bogus", () => {
		let grad = 2;
		expect($3Dmol.Gradient.getGradient(grad)).toBe(2);
	});

	test("Gradient get RWB", () => {
		let grad = {};
		grad.gradient = "rwb";
		grad = $3Dmol.Gradient.getGradient(grad);
		expect(grad).toBeInstanceOf($3Dmol.Gradient.RWB);
	});

	test("Gradient get RWB predefine mid, max, and min", () => {
		let grad = {};
		grad.gradient = "rwb";
		grad.min = 0;
		grad.max = 100;
		grad.mid = 50;
		grad = $3Dmol.Gradient.getGradient(grad);
		expect(grad).toBeInstanceOf($3Dmol.Gradient.RWB);
	});
})

// RWB Tests
describe("RWB Tests", () => {
	test("RWB no args", () => {
		const grad = new $3Dmol.Gradient.RWB();
		expect(grad.range()).toBeNull();
	});

	test("RWB Constructor with range", () => {
		const grad = new $3Dmol.Gradient.RWB([0, 100]);
		expect(grad.range()).toEqual([0, 100]);
	});

	test("RWB vth normal", () => {
		let grad = new $3Dmol.Gradient.RWB([0, 100]);
		let hexVal = grad.valueToHex(10, grad.range());
		expect(hexVal).toEqual(0xff7272);
	});

	test("RWB vth min max mid args", () => {
		const grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		const hexVal = grad.valueToHex(10, [0, 50, 100]);
		expect(hexVal).toEqual(0xff5050);
	});

	test("RWB vth min max mid args no range", () => {
		const grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		const hexVal = grad.valueToHex(10);
		expect(hexVal).toEqual(0xff7272);
	});

	test("RWB vth min max mid args, val > mid", () => {
		const grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		const hexVal = grad.valueToHex(60, [0, 100, 50]);
		expect(hexVal).toEqual(0xe4e4ff);
	});

	test("RWB vth min max mid args, val == mid", () => {
		const grad = new $3Dmol.Gradient.RWB(0, 100, 50);
		const hexVal = grad.valueToHex(50, [0, 100, 50]);
		expect(hexVal).toEqual(0xffffff);
	});
});

// ROYGB Tests
describe("ROYGB Tests", () => {
	test("ROYGB no args", () => {
		let grad = new $3Dmol.Gradient.ROYGB();
		expect(grad.range()).toEqual(null);
	})
	test("ROYGB Constructor with range", () => {
		let grad = new $3Dmol.Gradient.ROYGB([0, 100]);
		expect(grad.range()).toEqual([0, 100]);
	})
	test("ROYGB vth normal", () => {
		let grad = new $3Dmol.Gradient.ROYGB([0, 100]);
		let hexVal = grad.valueToHex(10, grad.range());
		expect(hexVal).toEqual(0xffa100);
	})
	test("ROYGB vth min max mid args", () => {
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(10, [0, 50, 100]);
		expect(hexVal).toEqual(0xffe400);
	})
	test("ROYGB vth min max mid args no range", () => {
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(10);
		expect(hexVal).toEqual(0xffa100);
	})
	test("ROYGB vth min max mid args, val < q3", () => {
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(60, [0, 100, 50]);
		expect(hexVal).toEqual(0xffa1);
	})
	test("ROYGB vth min max mid args, val < mid", () => {
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(49, [0, 100, 50]);
		expect(hexVal).toEqual(0x33ff00);
	})
	test("ROYGB vth min max mid args, val > max", () => {
		let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
		let hexVal = grad.valueToHex(101, [0, 100, 50]);
		expect(hexVal).toEqual(0xff);
	})
})

// Sinebow Tests
describe("Sinebow Tests", () => {
	test("Sinebow no args", () => {
		let grad = new $3Dmol.Gradient.Sinebow();
		expect(grad.range()).toEqual(null);
	})
	test("Sinebow Constructor with range", () => {
		let grad = new $3Dmol.Gradient.Sinebow([0, 100]);
		expect(grad.range()).toEqual([0, 100]);
	})
	test("Sinebow vth normal", () => {
		let grad = new $3Dmol.Gradient.Sinebow([100, 0]);
		let hexVal = grad.valueToHex(10, grad.range());
		expect(hexVal).toEqual(0x7f11ed);
	})
	test("Sinebow vth min max mid args no range", () => {
		let grad = new $3Dmol.Gradient.Sinebow(0, 100, 50);
		expect(grad).toBeDefined();
		let hexVal = grad.valueToHex(10);
		expect(hexVal).toEqual(0xed7f11);
	})
})
