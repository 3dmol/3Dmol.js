/** @jest-environment node */
// Build cartoon geometry without a WebGL context or native canvas renderer.
jest.mock('canvas', () => ({}));
const fs = require('fs');
const path = require('path');
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!doctype html><html><body></body></html>');
global.window = dom.window;
global.document = dom.window.document;
global.URL.createObjectURL = jest.fn();
dom.window.URL.createObjectURL = jest.fn();
Object.defineProperty(global, '$3Dmol', {
    configurable: true,
    get: () => dom.window.$3Dmol,
    set: value => { dom.window.$3Dmol = value; }
});
const $3Dmol = require('../../build/3Dmol.js');
const protein = fs.readFileSync(path.join(__dirname, '../test_structs/1YCR.pdb'), 'utf8');
afterAll(() => {
    dom.window.close();
    delete global.window;
    delete global.document;
    delete global.$3Dmol;
});

function fixture() {
    const model = new $3Dmol.GLModel(0);
    model.addMolData(protein, 'pdb');
    const atoms = model.selectedAtoms({});
    const group = { add: jest.fn(), remove: jest.fn() };
    return {
        model, atoms, group,
        draw: cartoon => {
            model.setStyle({}, { cartoon });
            model.globj(group, { supportsImposters: false, supportsAIA: false });
        }
    };
}

describe('cartoon secondary-structure preservation', () => {
    test('arrows, helix tubes and trace can be redrawn without changing the molecule', () => {
        const { atoms, group, draw } = fixture();
        const secondaryStructure = atoms.map(atom => atom.ss);
        for (const style of [
            { style: 'oval', arrows: true },
            { style: 'oval', arrows: true, tubes: true },
            { style: 'trace', arrows: true, tubes: true },
            { style: 'rectangle' }
        ]) {
            draw(style);
            expect(atoms.map(atom => atom.ss)).toEqual(secondaryStructure);
        }
        expect(group.add).toHaveBeenCalledTimes(4);
    });

    test('a failed geometry build also restores temporary secondary-structure markers', () => {
        const { atoms, draw } = fixture();
        const secondaryStructure = atoms.map(atom => atom.ss);
        // Fail after arrow/tube preprocessing, while reading the first residue.
        atoms.find(atom => atom.atom === 'CA').resn = undefined;
        expect(() => draw({ arrows: true, tubes: true })).toThrow();
        expect(atoms.map(atom => atom.ss)).toEqual(secondaryStructure);
    });
});
