
global.$ = require("jquery");
global.URL.createObjectURL = function() {};
let $3Dmol = require("../../build/3Dmol.js");
const fs = require('fs');

describe('Function json |', ()=>{
    const data = '{"m":[{"a":[{"x":85,"y":144},{"x":102.32050807568878,"y":134},{"x":67.67949192431124,"y":134},{"x":85,"y":164}],"b":[{"b":0,"e":1},{"b":0,"e":2},{"b":0,"e":3}]},{"a":[{"x":213,"y":131},{"x":230.32050807568876,"y":121},{"x":247.64101615137753,"y":131},{"x":264.9615242270663,"y":121}],"b":[{"b":0,"e":1},{"b":1,"e":2},{"b":2,"e":3}]}]}';
    let atoms = $3Dmol.Parsers.json(data, {});

    test("Atoms should match the snapshot", ()=>{
        expect(atoms).toMatchSnapshot();
    });

});
