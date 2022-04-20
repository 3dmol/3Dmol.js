/**
*@jest-environment jsdom
*/
const $3Dmol = require("../../../build/3Dmol.js");

describe('Gradient.ROYGB',() => {
  test('ROYGB no args', () => {
    let grad = new $3Dmol.Gradient.ROYGB();
    expect(grad.range()).toEqual(null);
  });
  test('ROYGB Constructor with range', () => {
    let grad = new $3Dmol.Gradient.ROYGB([0, 100]);
    let hexVal = grad.valueToHex(10, grad.range());
    expect(hexVal).toEqual(0xffa100);
  });
  test('ROYGB vth normal', () => {
    let grad = new $3Dmol.Gradient.ROYGB([0, 100]);
    expect(grad.range()).toEqual([0, 100]);
  });
  test('ROYGB vth min max mid args', () => {
    let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
    let hexVal = grad.valueToHex(10, [0, 50, 100]);
    expect(hexVal).toEqual(0xffe400);
  });
  test('ROYGB vth min max mid args no range', () => {
    let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
    let hexVal = grad.valueToHex(10);
    expect(hexVal).toEqual(0xffa100);
  });
  test('ROYGB vth min max mid args, val < q3', () => {
    let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
    let hexVal = grad.valueToHex(60, [0, 100, 50]);
    expect(hexVal).toEqual(0xffa1);
  });
  test('ROYGB vth min max mid args, val < mid', () => {
    let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
    let hexVal = grad.valueToHex(49, [0, 100, 50]);
    expect(hexVal).toEqual(0x33ff00);
  });
  test('ROYGB vth min max mid args, val > max', () => {
    let grad = new $3Dmol.Gradient.ROYGB(0, 100, 50);
    let hexVal = grad.valueToHex(101, [0, 100, 50]);
    expect(hexVal).toEqual(0xff);
  });
});
