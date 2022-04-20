/**
*@jest-environment jsdom
*/
const $3Dmol = require("../../../build/3Dmol.js");

describe('Gradient', () => {
  test('Gradient normalize value, 0, 100, 50', () => {
    let normFunc = $3Dmol.Gradient.normalizeValue;
    let obj = normFunc(0, 100, 50);
    expect(obj).toEqual({lo: 0, hi: 100, val: 50});
  });
  test('Gradient normalize value, 0, 100, -1 Return: 0, 100, 0', () => {
    let normFunc = $3Dmol.Gradient.normalizeValue;
    let obj = normFunc(0, 100, -1);
    expect(obj).toEqual({lo: 0, hi: 100, val: 0});
  });
  test('Gradient normalize value, 0, 100, 101 Return: 0, 100, 100', () => {
    let normFunc = $3Dmol.Gradient.normalizeValue;
    let obj = normFunc(0, 100, 101);
    expect(obj).toEqual({lo: 0, hi: 100, val: 100});
  });
  test('Gradient normalize value, 100, 0, 50', () => {
    let normFunc = $3Dmol.Gradient.normalizeValue;
    let obj = normFunc(100, 0, 50);
    expect(obj).toEqual({lo: 0, hi: 100, val: 50});
  });
  test('Gradient normalize value, 100, 0, -1 Return: 0, 100, 100', () => {
    let normFunc = $3Dmol.Gradient.normalizeValue;
    let obj = normFunc(100, 0, -1);
    expect(obj).toEqual({lo: 0, hi: 100, val: 100});
  });
  test('Gradient normalize value, 100, 0, 101, Return: 0, 100, 0', () => {
    let normFunc = $3Dmol.Gradient.normalizeValue;
    let obj = normFunc(100, 0, 101);
    expect(obj).toEqual({lo: 0, hi: 100, val: 0});
  });
  test('Gradient get happy path', () => {
    let grad = new $3Dmol.Gradient();
    let obj = $3Dmol.Gradient.getGradient(grad);
    expect(obj).toEqual({});
  });
  test('Gradient get bogus', () => {
    let grad = 2;
    expect($3Dmol.Gradient.getGradient(grad)).toBe(2);
  });
  test('Gradient get RWB', () => {
    let grad = new $3Dmol.Gradient.RWB();
    grad.gradient = 'rwb';
    grad = $3Dmol.Gradient.getGradient(grad);
    expect(grad).toBeInstanceOf($3Dmol.Gradient.RWB);
  });
  test('Gradient get RWB predefine mid, max, and min', () => {
    let grad = new $3Dmol.Gradient.RWB();
    grad.gradient = 'rwb';
    grad.min = 0;
    grad.max = 100;
    grad.mid = 50;
    grad = $3Dmol.Gradient.getGradient(grad);
    expect(grad).toHaveProperty('min', 0);
    expect(grad).toHaveProperty('max', 100);
    expect(grad).toHaveProperty('mid', 50);
    expect(grad).toBeInstanceOf($3Dmol.Gradient.RWB);
  });
});
