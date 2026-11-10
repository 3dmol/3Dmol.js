global.$ = require("jquery");
global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

function oldShapes() {
    return { sphere: ["stale"], cylinder: ["stale"], line: ["stale"], triangle: ["stale"] };
}

describe("picking-shape rebuilds", () => {
    test("rebuilds bounded, current picking shapes across revisited frames", async () => {
        const model = new $3Dmol.GLModel(0);
        model.addMolData("3\nfirst\nC 0 0 0\nO 1.4 0 0\nN 8 0 0\n3\nsecond\nC 0 2 0\nO 1.4 2 0\nN 8 2 0\n", "xyz", {
            multimodel: true, frames: true
        });
        model.setStyle({}, { sphere: { radius: 0.5 }, stick: { radius: 0.2 }, line: {} });

        const frames = [];
        const noninteractiveShapes = [];
        for (let frame = 0; frame < model.getNumFrames(); frame++) {
            await model.setFrame(frame);
            model.setClickable({ index: 0 }, true, () => {});
            model.setHoverable({ index: 1 }, true, () => {}, () => {});
            const atoms = model.selectedAtoms({});
            expect(atoms[0].bonds).toEqual([1]);
            expect(atoms[1].bonds).toEqual([0]);
            atoms.forEach(atom => { atom.intersectionShape = oldShapes(); });
            frames.push(atoms);
            noninteractiveShapes.push(atoms[2].intersectionShape);
        }

        // Real geometry generation, without a WebGL renderer or retained scene meshes.
        const group = { add() {}, remove() {} };

        for (let rebuild = 0; rebuild < 200; rebuild++) {
            const frame = rebuild % frames.length;
            await model.setFrame(frame);
            const atoms = model.selectedAtoms({});
            // Changing coordinates on each revisit also catches retained stale hit regions.
            atoms[0].y = atoms[1].y = rebuild;
            model.globj(group, { supportsImposters: true, supportsAIA: false });

            for (const atom of atoms.slice(0, 2)) {
                const shapes = atom.intersectionShape;
                expect(shapes.sphere).toHaveLength(2); // atom sphere and bond endpoint
                expect(shapes.cylinder).toHaveLength(1);
                expect(shapes.line).toHaveLength(2);
                expect(shapes.triangle).toHaveLength(0);
                for (const sphere of shapes.sphere) {
                    expect(sphere.center).toMatchObject({ x: atom.x, y: rebuild, z: 0 });
                }
                expect(shapes.cylinder[0].c1).toMatchObject({ x: atom.x, y: rebuild, z: 0 });
                expect(shapes.cylinder[0].c2).toMatchObject({ x: 0.7, y: rebuild, z: 0 });
                expect(shapes.line.map(point => point.y)).toEqual([rebuild, rebuild]);
            }
            expect(atoms[2].intersectionShape).toBe(noninteractiveShapes[frame]);
            expect(atoms[2].intersectionShape).toEqual(oldShapes());
        }
    });
});
