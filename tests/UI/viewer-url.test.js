/** @jest-environment node */
// The URL adapter stays in viewer.js; verify it against actual model selection
// and style semantics, not a second implementation of the new UI state.
jest.mock('canvas', () => ({}));
const fs = require('node:fs');
const path = require('node:path');
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!doctype html><html><body></body></html>', { url: 'https://example.test/viewer.html' });
const previousWindow = global.window;
const previousDocument = global.document;
const previous3Dmol = Object.getOwnPropertyDescriptor(global, '$3Dmol');
const previousCreateObjectURL = global.URL.createObjectURL;
global.window = dom.window;
global.document = dom.window.document;
global.URL.createObjectURL = jest.fn();
dom.window.URL.createObjectURL = jest.fn();
Object.defineProperty(global, '$3Dmol', { configurable: true, get: () => dom.window.$3Dmol, set: value => { dom.window.$3Dmol = value; } });
const $3Dmol = require('../../build/3Dmol.js');
const ViewerURL = require('../../viewer.js');
const parse = url => ViewerURL.parse(url, $3Dmol.specStringToObject);
const baseURL = 'https://example.test/viewer.html';
const pdb = fs.readFileSync(path.join(__dirname, '../test_structs/1YCR.pdb'), 'utf8');
const tutorial = fs.readFileSync(path.join(__dirname, '../../tutorials/url.md'), 'utf8');
const examples = [...new Set([...tutorial.matchAll(/(?:href=["']|\]\()([^"'\s)]*viewer\.html\?[^"'\s)]*)/g)].map(match => match[1].replace(/&amp;/g, '&')))];

afterAll(() => {
    dom.window.close();
    if (previousWindow === undefined) delete global.window;
    else global.window = previousWindow;
    if (previousDocument === undefined) delete global.document;
    else global.document = previousDocument;
    if (previous3Dmol) Object.defineProperty(global, '$3Dmol', previous3Dmol);
    else delete global.$3Dmol;
    global.URL.createObjectURL = previousCreateObjectURL;
});

function model() {
    const result = new $3Dmol.GLModel(0);
    result.addMolData(pdb, 'pdb');
    result.setStyle({}, { line: {} });
    return result;
}

function applySelections(scene) {
    const result = model();
    for (const selection of scene.selections) {
        if (selection.hidden) continue;
        const styles = selection.styles.filter(style => !style.hidden);
        if (!styles.length) continue;
        if (!selection.additive) result.setStyle(selection.spec, {});
        for (const style of styles) result.setStyle(selection.spec, style.spec, true);
    }
    return result;
}

function legacyReference(url) {
    const result = model();
    let selection = {};
    const surfaces = [], labels = [];
    for (const [command, text] of new URL(url, baseURL).searchParams) {
        const spec = $3Dmol.specStringToObject(text);
        if (command === 'select') selection = spec;
        else if (command === 'style') result.setStyle(selection, spec);
        else if (command === 'addstyle') result.setStyle(selection, spec, true);
        else if (command === 'surface') surfaces.push({ spec: selection, style: spec });
        else if (command === 'labelres') labels.push({ spec: selection, style: spec });
    }
    return { result, surfaces, labels };
}

const appearances = model => model.selectedAtoms({}).map(atom => atom.style);

describe('standalone viewer URL API', () => {
    test('extracts the actual linked tutorial examples', () => {
        expect(examples.length).toBeGreaterThanOrEqual(9);
    });

    test.each(examples)('preserves documented selection/style/label/surface semantics: %s', url => {
        const scene = parse(url);
        const expected = legacyReference(url);
        expect(scene.models).toHaveLength(1);
        expect(appearances(applySelections(scene))).toEqual(appearances(expected.result));
        expect(scene.surfaces).toHaveLength(expected.surfaces.length);
        for (const [index, surface] of scene.surfaces.entries()) {
            expect(surface.type).toBe(1);
            expect(surface.spec).toEqual(expected.surfaces[index].spec);
            expect(surface.calculationSpec).toEqual(expected.surfaces[index].spec);
            expect(surface.style).toEqual(expected.surfaces[index].style);
        }
        const labels = scene.selections.filter(selection => selection.labelsEnabled);
        expect(labels).toHaveLength(expected.labels.length);
        for (const [index, selection] of labels.entries()) {
            expect(selection.spec).toEqual(expected.labels[index].spec);
            expect(selection.labelStyle).toEqual(expected.labels[index].style);
            expect(expected.result.selectedAtoms(selection.spec).length).toBeGreaterThan(0);
        }
        const restored = parse(ViewerURL.serialize(scene, baseURL));
        expect(appearances(applySelections(restored))).toEqual(appearances(expected.result));
        expect(restored.surfaces).toEqual(scene.surfaces);
        expect(restored.selections).toEqual(scene.selections);
    });

    test('style replaces while addstyle layers onto earlier representations, in query order', () => {
        const url = '?pdb=1YCR&style=cartoon&select=chain:B&addstyle=stick:radius~0.2&style=sphere:scale~0.5&addstyle=line&select=chain:A&addstyle=cross';
        const actual = applySelections(parse(url));
        expect(appearances(actual)).toEqual(appearances(legacyReference(url).result));
        expect(actual.selectedAtoms({ chain: 'B' })[0].style).toEqual({ sphere: { scale: 0.5 }, line: {} });
        expect(actual.selectedAtoms({ chain: 'A' })[0].style).toEqual({ cartoon: {}, cross: {} });
    });

    test.each(['?pdb=1YCR&style=%7B%7D', '?pdb=1YCR&style=cartoon&select=chain:B&style=%7B%7D', '?pdb=1YCR&style=cartoon&select=chain:B&style=%7B%7D&addstyle=line'])('an explicitly empty style clears the selected atoms: %s', url => {
        expect(appearances(applySelections(parse(url)))).toEqual(appearances(legacyReference(url).result));
    });

    test('retains repeated surfaces and residue-label layers on the same selection', () => {
        const scene = parse('?pdb=1YCR&select=chain:B&surface=color:red;opacity:0.2&surface=color:blue;opacity:0.4&labelres=fontSize:12&labelres=fontSize:18');
        expect(scene.surfaces.map(surface => surface.style)).toEqual([{ color: 'red', opacity: 0.2 }, { color: 'blue', opacity: 0.4 }]);
        expect(scene.selections.filter(selection => selection.labelsEnabled).map(selection => selection.labelStyle.fontSize)).toEqual([12, 18]);
    });

    test('decodes source URLs once without splitting encoded query values or losing plus signs', () => {
        const source = 'https://files.example.test/molecule.sdf?token=a=b&name=A+B%2FC#part';
        const scene = parse('?url=' + encodeURIComponent(source) + '&style=stick:color~%23000066');
        expect(scene.models[0].value).toBe(source);
        expect(scene.models[0].format).toBe('sdf');
        expect(scene.selections.flatMap(selection => selection.styles).some(style => style.spec.stick?.color === '#000066')).toBe(true);
        expect(parse(ViewerURL.serialize(scene, baseURL)).models[0].value).toBe(source);
    });

    test('retains literal plus signs from legacy source, selection, and session URLs', () => {
        const scene = parse('?url=A+B%20C.sdf&select=chain:A+B&style=stick');
        expect(scene.models[0].value).toBe('A+B C.sdf');
        expect(scene.selections[0].spec.chain).toBe('A+B');
        expect(parse('?session=Class+One').session).toBe('Class+One');
        expect(parse('?session=Class%2BOne').session).toBe('Class+One');
        expect(parse('?session=Class%20One').session).toBe('Class One');
        expect(parse(ViewerURL.serialize(scene, baseURL))).toEqual(scene);
    });

    test.each([
        ['cube', '../test_structs/h-bn.cube'],
        ['vasp', '../auto/data/mg.vasp'],
        ['mcif', '../test_structs/1ubq.cif']
    ])('loads the core-supported %s format with inferred or explicit type', (format, filename) => {
        const data = fs.readFileSync(path.join(__dirname, filename), 'utf8');
        for (const query of [`?url=molecule.${format}`, `?url=download&type=${format}`]) {
            const scene = parse(query);
            const model = new $3Dmol.GLModel(0);
            model.addMolData(data, scene.models[0].format);
            expect(model.selectedAtoms({}).length).toBeGreaterThan(0);
            expect(parse(ViewerURL.serialize(scene, baseURL))).toEqual(scene);
        }
    });

    test('allows parsers registered by the application without a viewer format whitelist', () => {
        $3Dmol.Parsers.viewer_test_format = $3Dmol.Parsers.xyz;
        try {
            for (const query of ['?url=download&type=viewer_test_format', '?url=molecule.viewer_test_format']) {
                const scene = parse(query);
                expect(scene.models[0].format).toBe('viewer_test_format');
                const model = new $3Dmol.GLModel(0);
                model.addMolData('1\ncustom parser\nC 0 0 0\n', scene.models[0].format);
                expect(model.selectedAtoms({})).toHaveLength(1);
            }
        } finally { delete $3Dmol.Parsers.viewer_test_format; }
    });

    test.each(['download?id=1YCR', 'molecule.unknown', 'download.gz', 'molecule.unknown.gz'])('retains core content detection for %s', source => {
        const scene = parse('?url=' + encodeURIComponent(source));
        const compressed = source.endsWith('.gz');
        const data = compressed ? require('node:zlib').gzipSync(pdb) : pdb;
        const model = new $3Dmol.GLModel(0);
        const log = jest.spyOn(console, 'log').mockImplementation(() => {});
        try { model.addMolData(data, scene.models[0].format); }
        finally { log.mockRestore(); }
        expect(model.selectedAtoms({})).toHaveLength(818);
        expect(scene.models[0].format.endsWith('.gz')).toBe(compressed);
        expect(parse(ViewerURL.serialize(scene, baseURL))).toEqual(scene);
    });

    test('infers compressed formats from paths and accepts explicit format away from the source', () => {
        expect(parse('?url=tests/test_structs/benzene.sdf.gz?download=1').models[0].format).toBe('sdf.gz');
        const explicit = parse('?url=download?id=12&style=stick&type=sdf');
        expect(explicit.models[0].value).toBe('download?id=12');
        expect(explicit.models[0].format).toBe('sdf');
        expect(parse('?url=molecule.sdf&style=stick').models[0].value).toBe('molecule.sdf');
    });

    test.each([['pdb', '1YCR'], ['cif', '1UBQ'], ['cid', '241'], ['mmtf', '1YCR']])('recognizes %s structure identifiers', (kind, value) => {
        expect(parse(`?${kind}=${value}`).models[0]).toMatchObject({ kind, value });
    });

    test('reads view and trajectory options without losing boolean false or frame zero', () => {
        const scene = parse('?url=motion.xyz&frames=false&background=%23000066&backgroundalpha=0.5&projection=orthographic&frame=0&fps=24&loop=false&playing=false&open=false&view=1,2,3,-75,0,0,0,1');
        expect(scene.models[0].frames).toBe(false);
        expect(scene.open).toBe(false);
        expect(scene.view).toEqual({ background: '#000066', backgroundAlpha: 0.5, projection: 'orthographic', frame: 0, fps: 24, loop: false, playing: false, camera: [1, 2, 3, -75, 0, 0, 0, 1] });
    });

    test.each([-1, -2, -11, 7, 100])('preserves integer frame request %i for trajectory wrapping', frame => {
        const scene = parse(`?url=motion.xyz&frame=${frame}`);
        expect(scene.view.frame).toBe(frame);
        expect(parse(ViewerURL.serialize(scene, baseURL))).toEqual(scene);
    });

    test('supports encoded nested JSON atom selections and complete style objects', () => {
        const spec = { and: [{ chain: ['A', 'B'] }, { not: { resn: 'HOH' } }], within: { distance: 4, sel: { resn: 'TRP' } } };
        const style = { stick: { radius: 0.2, colorscheme: { prop: 'b', gradient: 'rwb', min: 0, max: 40 } }, sphere: { scale: 0.25 } };
        const scene = parse(`?pdb=1YCR&select=${encodeURIComponent(JSON.stringify(spec))}&style=${encodeURIComponent(JSON.stringify(style))}`);
        const layer = scene.selections.find(selection => JSON.stringify(selection.spec) === JSON.stringify(spec));
        expect(layer).toBeDefined();
        expect(Object.assign({}, ...layer.styles.map(item => item.spec))).toEqual(style);
        expect(appearances(applySelections(scene))).toEqual(appearances(legacyReference(`?pdb=1YCR&select=${encodeURIComponent(JSON.stringify(spec))}&style=${encodeURIComponent(JSON.stringify(style))}`).result));
    });

    test('retains legacy element and property color commands with their atom scopes', () => {
        const scene = parse('?pdb=1YCR&select=chain:A&colorbyelement=colorscheme:greenCarbon&select=chain:B&colorbyproperty=prop:b;scheme:RWB');
        expect(scene.atomColors).toEqual([
            expect.objectContaining({ type: 'element', spec: { chain: 'A' }, colorscheme: 'greenCarbon' }),
            expect.objectContaining({ type: 'property', spec: { chain: 'B' }, prop: 'b', scheme: 'RWB' })
        ]);
        expect(parse(ViewerURL.serialize(scene, baseURL)).atomColors).toEqual(scene.atomColors);
    });

    test('canonical scene links retain every UI setting, source, and surface snapshot frame', () => {
        const scene = {
            models: [
                { kind: 'pdb', value: '1YCR', format: 'pdb', frames: false, title: 'Protein & ligand', hidden: true },
                { kind: 'file', value: 'two frames.xyz', format: 'xyz', frames: true, title: 'Local motion', data: { encoding: 'text', value: '1\nfirst\nC 0 0 0\n1\nsecond\nC 1 0 0\n' }, baseStyles: [{ indices: [0], style: { sphere: { scale: 0.5 } } }] },
                { kind: 'file', value: 'benzene.sdf.gz', format: 'sdf.gz', frames: false, data: { encoding: 'base64', value: Buffer.from([0, 255, 128, 43]).toString('base64') } }
            ],
            selections: [{ name: 'A+B & "pocket"', spec: { model: [0, 1], resi: ['10-20', 25], not: { elem: 'H' } }, styles: [{ spec: { sphere: { scale: 0.7, opacity: 0.3 } }, hidden: false }, { spec: { cross: { radius: 0.2 } }, hidden: true }], hidden: true, additive: true, labelsEnabled: true, labelStyle: { fontSize: 14, backgroundOpacity: 0.8, fontColor: '#000066' } }],
            surfaces: [{ name: 'Pocket surface', spec: { model: 0, chain: 'B' }, calculationSpec: { model: 0 }, style: { colorscheme: 'greenCarbon', opacity: 0.65 }, type: 2, frame: 0 }],
            view: { background: '#010266', backgroundAlpha: 0.5, projection: 'orthographic', camera: [1, 2, 3, -75, 0, 0, 0, 1], frame: 1, fps: 24, loop: false, playing: false },
            open: false
        };
        const url = ViewerURL.serialize(scene, baseURL + '?pdb=obsolete#old');
        expect(new URL(url).hash.startsWith('#state=')).toBe(true);
        expect(new URL(url).search).toBe('');
        expect(parse(url)).toEqual(scene);
    });

    test('a saved scene takes precedence over stale query settings and does not pollute prototypes', () => {
        const scene = parse('?pdb=1YCR&style=stick');
        const serialized = new URL(ViewerURL.serialize(scene, baseURL));
        serialized.search = '?pdb=4KW4&style=cartoon';
        expect(parse(serialized.href)).toEqual(scene);
        serialized.search = '?url=malformed%escape';
        expect(parse(serialized.href)).toEqual(scene);
        expect(() => parse('#state=' + encodeURIComponent('{"__proto__":{"polluted":true}}'))).toThrow();
        expect({}.polluted).toBeUndefined();
    });

    test.each(['#state=%7Bbroken', '#state=' + encodeURIComponent('{"models":[]}'), '?url=javascript:alert(1)', '?pdb=not-a-pdb', '?cid=one', '?pdb=1YCR&notACommand=anything'])('reports invalid URL state or structure identifiers: %s', url => {
        expect(() => parse(url)).toThrow();
    });

    test.each([{ fps: 0 }, { frame: 1.5 }, { backgroundAlpha: 2 }, { projection: 'fisheye' }, { slab: { near: 'bad', far: 10 } }])('rejects invalid saved view settings before applying a scene: %j', view => {
        const scene = { models: [], selections: [], surfaces: [], view };
        expect(() => parse('#state=' + encodeURIComponent(JSON.stringify(scene)))).toThrow();
    });
});
