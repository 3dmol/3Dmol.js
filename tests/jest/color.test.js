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

    test("should set the RGB values from the object when val is an object", () => {
        const color = new $3Dmol.Color();
        const obj = { r: 0.5, g: 0.8, b: 1.0 };

        color.set(obj);

        expect(color.r).toBe(0.5);
        expect(color.g).toBe(0.8);
        expect(color.b).toBe(1.0);
    });

    test("should return an array of hex values when input is an array", () => {
        const color = new $3Dmol.Color();
        const hexArray = ["#ff0000", "#00ff00", "#0000ff"];
        const expectedResult = [0xff0000, 0x00ff00, 0x0000ff];

        const result =  color.getHex(hexArray);

        expect(result).toEqual(expectedResult);
    });

    test("should expand short hex strings to full hex number", () => {
        const color = new $3Dmol.Color();
        const hexString = "#abc";
        const expectedResult = "#aabbcc";

        const result = color.getHex(hexString);

        expect(result).toBe(expectedResult);
    });

    test("should return the color value from htmlColors object when input is a string", () => {
        const color = new $3Dmol.Color();
        // Mocking the htmlColors object
        const htmlColors = {
            red: 0xff0000,
            green: 0x00ff00,
            blue: 0x0000ff,
        };

        // Mocking the window object
        const mockWindow = {
            $3Dmol: {
                htmlColors,
            },
        };

        // Providing the mock window object to the global scope
        global.window = mockWindow;

        const colorString = "green";
        const expectedResult = 0x00ff00;

        const result = color.getHex(colorString);

        expect(result).toBe(expectedResult);
    });
})
