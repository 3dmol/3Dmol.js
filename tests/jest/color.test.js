/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function () { };
let $3Dmol = require("../../build/3Dmol.js");


describe("Color", () => {
    test("Has an identity constructor", () => {
        expect(new $3Dmol.Color()).toBeInstanceOf($3Dmol.Color);
    })

    test("Has a 0 to 1 based number constructor", () => {
        expect(new $3Dmol.Color(1, 0, 0)).toBeInstanceOf($3Dmol.Color);
    })

    test("Has a copy constructor", () => {
        let color = new $3Dmol.Color(1, 0, 0);
        let copy = new $3Dmol.Color(color)
        expect(copy).toBeInstanceOf($3Dmol.Color);
        expect(copy).toMatchObject(color);
    })

    test("Has a hex constructor", () => {
        let color = new $3Dmol.Color(0xff0000);
        expect(color).toBeInstanceOf($3Dmol.Color);
        expect(color).toMatchObject({ r: 1, g: 0, b: 0 });
    })

    test("Has Colored object constructor", () => {
        let color = new $3Dmol.Color({ r: 1, g: 0, b: 0 });
        expect(color).toBeInstanceOf($3Dmol.Color);
        expect(color).toMatchObject({ r: 1, g: 0, b: 0 });
    })

    test("Has a hex getter", () => {
        let color = new $3Dmol.Color(0xff0000);
        expect(color.getHex()).toBe(0xff0000);
    })

    test("Has a hex setter", () => {
        let color = new $3Dmol.Color();
        expect(color.getHex()).toBe(0x000000);
        color.setHex(0x00ff00);
        expect(color.getHex()).toBe(0x00ff00);
    })

    test("Can copy a colored obj", () => {
        let color = new $3Dmol.Color();
        let other = new $3Dmol.Color({ r: 1, g: 0, b: 0 });
        color.copy(other);
        expect(color).toBeInstanceOf($3Dmol.Color);
        expect(color).toMatchObject({ r: 1, g: 0, b: 0 });
    })

    test("Has a clone fn", () => {
        let color = new $3Dmol.Color(0xff0000);
        let other = color.clone();
        expect(other).toBeInstanceOf($3Dmol.Color);
        expect(other).toMatchObject({ r: 1, g: 0, b: 0 });
        expect(other).not.toBe(color);
    })

    test("Can produce a 255 scaled color", () => {
        let color = new $3Dmol.Color(0xff0000);
        let scaled = color.scaled();
        expect(scaled).toEqual(expect.objectContaining({ r: 255, g: 0, b: 0, a: 1 }));
    })

    test("Has r, g, and b properties that are zero initialized", () => {
        let color = new $3Dmol.Color();
        expect(color).toMatchObject({ r: 0, g: 0, b: 0 });
    })
})
