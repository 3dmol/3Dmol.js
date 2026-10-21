const { resetIntersectionShapes } = require("../../src/intersectionShapes");

function oldShapes() {
    return { sphere: ["stale"], cylinder: ["stale"], line: ["stale"], triangle: ["stale"] };
}

describe("picking-shape rebuilds", () => {
    test("keeps clickable and hoverable frame shapes bounded on revisits", () => {
        const frames = [
            [
                { clickable: true, intersectionShape: oldShapes() },
                { hoverable: true, intersectionShape: oldShapes() },
                { intersectionShape: { sphere: ["unchanged"] } },
            ],
            [{ clickable: true, intersectionShape: oldShapes() }],
        ];
        const noninteractiveShapes = frames[0][2].intersectionShape;

        for (let rebuild = 0; rebuild < 200; rebuild++) {
            const frame = frames[rebuild % frames.length];
            resetIntersectionShapes(frame);
            for (const atom of frame) {
                if (atom.clickable || atom.hoverable) {
                    atom.intersectionShape.sphere.push(rebuild);
                    atom.intersectionShape.cylinder.push(rebuild);
                }
            }
        }

        for (const frame of frames) {
            for (const atom of frame) {
                if (atom.clickable || atom.hoverable) {
                    expect(atom.intersectionShape.sphere).toHaveLength(1);
                    expect(atom.intersectionShape.cylinder).toHaveLength(1);
                    expect(atom.intersectionShape.line).toHaveLength(0);
                    expect(atom.intersectionShape.triangle).toHaveLength(0);
                }
            }
        }
        expect(frames[0][2].intersectionShape).toBe(noninteractiveShapes);
    });
});
