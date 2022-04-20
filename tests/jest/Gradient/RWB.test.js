/**
*@jest-environment jsdom
*/
const $3Dmol = require("../../../build/3Dmol.js");

describe('Gradient.RWB',() => {
  test('RWB no args', () => {
    let grad = new $3Dmol.Gradient.RWB();
    expect(grad.range()).toEqual(null);
  });
  test('RWB Constructor with range', () => {
    let grad = new $3Dmol.Gradient.RWB([0, 100]);
    expect(grad.range()).toEqual([0, 100]);
  });
  test('RWB vth normal', () => {
    let grad = new $3Dmol.Gradient.RWB([0, 100]);
    let hexVal = grad.valueToHex(10, grad.range());
    expect(hexVal).toEqual(0xff7272);
  });
  test('RWB vth min max mid args', () => {
    let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
    let hexVal = grad.valueToHex(10, [0, 50, 100]);
    expect(hexVal).toEqual(0xff5050);
  });
  test('RWB vth min max mid args no range', () => {
    let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
    let hexVal = grad.valueToHex(10);
    expect(hexVal).toEqual(0xff7272);
  });
  test('RWB vth min max mid args, val > mid', () => {
    let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
    let hexVal = grad.valueToHex(60, [0, 100, 50]);
    expect(hexVal).toEqual(0xe4e4ff);
  });
  test('RWB vth min max mid args, val == mid', () =>{
    let grad = new $3Dmol.Gradient.RWB(0, 100, 50);
    let hexVal = grad.valueToHex(50, [0, 100, 50]);
    expect(hexVal).toEqual(0xffffff);
  });
})
