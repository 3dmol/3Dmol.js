/**
*@jest-environment jsdom
*/

global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");


describe("Color", () => {
    it("has an identity constructor", () => {
        expect(new $3Dmol.Color()).toBeInstanceOf($3Dmol.Color);
    })

    it("has a 0 to 1 based number constructor", () => {
        expect(new $3Dmol.Color(1,0,0)).toBeInstanceOf($3Dmol.Color);
    })

    it("has a copy constructor", () => {
        let color = new $3Dmol.Color(1,0,0);
        let copy = new $3Dmol.Color(color)
        expect(copy).toBeInstanceOf($3Dmol.Color);
        expect(copy).toHaveProperty("r", 1);
        expect(copy).toHaveProperty("g", 0);
        expect(copy).toHaveProperty("b", 0);
    })

    it("has a hex constructor", () => {
        let color = new $3Dmol.Color(0xff0000);
        expect(color).toBeInstanceOf($3Dmol.Color);
        expect(color).toHaveProperty("r", 1);
        expect(color).toHaveProperty("g", 0);
        expect(color).toHaveProperty("b", 0);
    })

    it("has Colored object constructor", () => {
        let color = new $3Dmol.Color({r:1, g:0, b:0});
        expect(color).toBeInstanceOf($3Dmol.Color);
        expect(color).toHaveProperty("r", 1);
        expect(color).toHaveProperty("g", 0);
        expect(color).toHaveProperty("b", 0);
    })

    it("has a hex getter", () => {
        let color = new $3Dmol.Color(0xff0000);
        expect(color.getHex()).toBe(0xff0000);
    })

    it("has a hex setter", () => {
        let color = new $3Dmol.Color();
        expect(color.getHex()).toBe(0x000000);
        color.setHex(0x00ff00);
        expect(color.getHex()).toBe(0x00ff00);
    })

    it("can copy a colored obj", () => {
        let color = new $3Dmol.Color();
        let other = new $3Dmol.Color({r:1, g:0, b:0});
        color.copy(other);
        expect(color).toBeInstanceOf($3Dmol.Color);
        expect(color).toHaveProperty("r", 1);
        expect(color).toHaveProperty("g", 0);
        expect(color).toHaveProperty("b", 0);
    })

    it("has a clone fn", () => {
        let color = new $3Dmol.Color(0xff0000);
        let other = color.clone();
        expect(other).toBeInstanceOf($3Dmol.Color);
        expect(other).toHaveProperty("r", 1);
        expect(other).toHaveProperty("g", 0);
        expect(other).toHaveProperty("b", 0);
        expect(other).not.toBe(color);
    })

    it("can produce a 255 scaled color", () => {
        let color = new $3Dmol.Color(0xff0000);
        let scaled = color.scaled();
        console.log(scaled);
        expect(JSON.stringify(scaled)).toEqual(JSON.stringify({r:255, g:0, b:0, a: 1}));
    })

    it("has r, g, and b properties that are zero initialized", () => {
        let color = new $3Dmol.Color();
        expect(color).toHaveProperty("r", 0);
        expect(color).toHaveProperty("g", 0);
        expect(color).toHaveProperty("b", 0);
    })
})