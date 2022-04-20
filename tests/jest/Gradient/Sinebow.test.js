/**
*@jest-environment jsdom
*/
const $3Dmol = require("../../../build/3Dmol.js");

describe('Gradient.Sinebow', () => {
  test('Sinebow no args', () => {
    let grad = new $3Dmol.Gradient.Sinebow();
    expect(grad.range()).toEqual(null);
  });
  test('Sinebow Constructor with range', () => {
    let grad = new $3Dmol.Gradient.Sinebow([0, 100]);
    expect(grad.range()).toEqual([0, 100]);
  });
  test('Sinebow vth normal', () => {
    let grad = new $3Dmol.Gradient.Sinebow([100, 0]);
    let hexVal = grad.valueToHex(10, grad.range());
    expect(hexVal).toEqual(0x7f11ed);
  });
  test('Sinebow vth min max mid args no range', () => {
    let grad = new $3Dmol.Gradient.Sinebow(0, 100, 50);
    let hexVal = grad.valueToHex(10);
    expect(hexVal).toEqual(0xed7f11);
  });
});
