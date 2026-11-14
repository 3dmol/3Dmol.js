global.URL.createObjectURL = function () {};
const $3Dmol = require("../../build/3Dmol.js");

const labels = [];
function label(style = {}) {
    // An image fixes the logical dimensions independently of platform font metrics.
    const image = document.createElement("canvas");
    image.width = 101;
    image.height = 41;
    const result = new $3Dmol.Label("Hi", { backgroundImage: image, ...style });
    labels.push(result);
    result.setContext();
    return result;
}
afterEach(() => { labels.splice(0).forEach(item => item.dispose()); });

test.each([
    [undefined, 202, 82], [1, 101, 41], [1.5, 152, 62], [3, 303, 123]
])("resolution %p changes raster density without changing logical size", (resolution, width, height) => {
    const result = label({ resolution });
    expect([result.canvas.width, result.canvas.height]).toEqual([width, height]);
    const texture = result.sprite.material.map;
    expect([texture.logicalWidth, texture.logicalHeight]).toEqual([101, 41]);
    expect(texture.image).toBe(result.canvas);
    const transform = result.context.getTransform();
    expect([transform.a, transform.d]).toEqual([width / 101, height / 41]);
    expect([result.sprite.scale.x, result.sprite.scale.y]).toEqual([1, 1]);
});

test.each([0, -1, 0.5, NaN, Infinity, "3"])("invalid resolution %p uses the default density", resolution => {
    const result = label({ resolution });
    expect([result.canvas.width, result.canvas.height]).toEqual([202, 82]);
});

test("redrawing and changing resolution replace the canvas transform instead of accumulating it", () => {
    const result = label();
    for (const resolution of [2, 3, 1, 2]) {
        result.dispose();
        result.stylespec.resolution = resolution;
        result.text = "New";
        result.setContext();
        const transform = result.context.getTransform();
        expect([transform.a, transform.b, transform.c, transform.d, transform.e, transform.f])
            .toEqual([resolution, 0, 0, resolution, 0, 0]);
    }
});

test("cloning a label texture retains its logical dimensions and image", () => {
    const result = label({ resolution: 3 });
    const texture = result.sprite.material.map;
    const copy = texture.clone();
    expect(copy).not.toBe(texture);
    expect(copy.image).toBe(result.canvas);
    expect([copy.logicalWidth, copy.logicalHeight]).toEqual([101, 41]);
    copy.dispose();
});
