/** @jest-environment node */
jest.mock('canvas', () => ({}));
const fs = require('fs');
const path = require('path');
const ts = require('typescript');
const { JSDOM } = require('jsdom');

const source = fs.readFileSync(path.join(__dirname, '../../src/autoload.ts'), 'utf8');
const compiled = ts.transpileModule(source, {
    compilerOptions: { module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2019 }
}).outputText;

describe('embedded control availability', () => {
    let dom, host, viewer, StateManager, autoload;

    beforeEach(() => {
        dom = new JSDOM('<!doctype html><div id="viewer" class="viewer_3Dmoljs"></div>');
        host = dom.window.document.getElementById('viewer');
        viewer = {
            setBackgroundColor: jest.fn(), setConfig: jest.fn(), setStyle: jest.fn(),
            addStyle: jest.fn(), render: jest.fn(), zoomTo: jest.fn(), addModel: jest.fn()
        };
        StateManager = jest.fn(function (attachedViewer, config) {
            this.config = config;
            this.setControls = jest.fn();
            this.createSelectionAndStyle = jest.fn();
            attachedViewer.ui = this;
        });
        dom.window.$3Dmol = { StateManager };
        const modules = {
            './GLViewer': { createViewer: jest.fn(() => viewer) },
            './ProteinSurface4': { SurfaceType: { VDW: 1 } },
            './colors': { CC: { color: () => 'white' } },
            './utilities': {
                specStringToObject: value => value === 'multimodel:true;frames:true'
                    ? { multimodel: true, frames: true } : {}
            }
        };
        const exports = {};
        new Function('exports', 'require', 'window', 'document', compiled)(
            exports, name => modules[name], dom.window, dom.window.document
        );
        // Exercise explicit initialization, independent of JSDOM's ready-state events.
        dom.window.document.onreadystatechange = null;
        autoload = exports.autoload;
    });

    afterEach(() => dom.window.close());

    test.each([
        [{ ui: 'true' }, { edit: true, trajectory: true }],
        [{ ui: 'true', uiEdit: 'false' }, { edit: false, trajectory: true }],
        [{ ui: 'true', uiTrajectory: 'false' }, { edit: true, trajectory: false }],
        [{ uiTrajectory: 'true' }, { edit: false, trajectory: true }],
        [{ uiEdit: 'true' }, { edit: true, trajectory: false }],
        [{ ui: 'false', uiTrajectory: 'true' }, { edit: false, trajectory: true }],
        [{ ui: 'false', uiEdit: 'true' }, { edit: true, trajectory: false }],
        [{ uiEdit: 'true', uiTrajectory: 'true' }, { edit: true, trajectory: true }]
    ])('initializes requested controls for %j', (attributes, expected) => {
        Object.assign(host.dataset, attributes);
        autoload();
        expect(StateManager).toHaveBeenCalledTimes(1);
        expect(StateManager).toHaveBeenCalledWith(viewer, expected);
    });

    test.each([
        {}, { ui: 'false' }, { ui: '' },
        { ui: 'true', uiEdit: 'false', uiTrajectory: 'false' },
        { uiEdit: 'false', uiTrajectory: 'false' }
    ])('does not create controls when all are disabled: %j', attributes => {
        Object.assign(host.dataset, attributes);
        autoload();
        expect(StateManager).not.toHaveBeenCalled();
    });

    test('updates an existing UI and can disable both components', () => {
        host.dataset.ui = 'true';
        autoload(viewer);
        const controls = viewer.ui;
        host.dataset.uiEdit = 'false';
        autoload(viewer);
        expect(controls.setControls).toHaveBeenLastCalledWith({ edit: false, trajectory: true });
        host.dataset.uiTrajectory = 'false';
        autoload(viewer);
        expect(controls.setControls).toHaveBeenLastCalledWith({ edit: false, trajectory: false });
        expect(StateManager).toHaveBeenCalledTimes(1);
        expect(viewer.ui).toBe(controls);
    });

    test('leaves a programmatically installed UI alone without data-ui options', () => {
        const controls = viewer.ui = { setControls: jest.fn() };
        autoload(viewer);
        expect(controls.setControls).not.toHaveBeenCalled();
        expect(StateManager).not.toHaveBeenCalled();
    });

    test('trajectory-only embeds still pass frame parser options when loading inline data', () => {
        const data = dom.window.document.createElement('script');
        data.id = 'frames';
        data.type = 'chemical/x-xyz';
        data.textContent = '1\nFirst\nC 0 0 0\n1\nSecond\nC 1 0 0\n';
        dom.window.document.body.append(data);
        Object.assign(host.dataset, {
            uiTrajectory: 'true', element: 'frames', type: 'xyz', options: 'multimodel:true;frames:true'
        });
        autoload();
        expect(viewer.addModel).toHaveBeenCalledWith(data.textContent, 'xyz', { multimodel: true, frames: true });
        expect(StateManager).toHaveBeenCalledWith(viewer, { edit: false, trajectory: true });
    });

    test('embeds work without the optional UI bundle', () => {
        delete dom.window.$3Dmol.StateManager;
        host.dataset.uiTrajectory = 'true';
        autoload();
        expect(StateManager).not.toHaveBeenCalled();
        expect(viewer.render).toHaveBeenCalled();
        expect(host.textContent).toBe('');
    });
});
