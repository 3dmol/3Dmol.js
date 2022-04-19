/**
 * @jest-environment jsdom
 */
const $3Dmol = require("../build/3Dmol.js");

describe("Sanity check", function() {
    it("should exist", () => {
        expect($3Dmol).toBeDefined();
    })

    it("should have properties", () => {
        expect(Object.keys($3Dmol).length).toBeGreaterThan(0);
    })
})