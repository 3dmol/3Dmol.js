/** @jest-environment node */
// The controls only need a DOM; they do not render through native canvas.
jest.mock('canvas', () => ({}));
const fs = require('fs');
const path = require('path');
const ts = require('typescript');
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!doctype html><html><body></body></html>');
const { document, Event } = dom.window;
afterAll(() => dom.window.close());

// Exercise the source module without depending on a generated viewer bundle.
const source = fs.readFileSync(path.join(__dirname, '../../src/ui/trajectory.ts'), 'utf8');
const compiled = ts.transpileModule(source, {
    compilerOptions: { module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2019 }
}).outputText;
const uiModule = {};
new Function('exports', compiled)(uiModule);
const { TrajectoryControls } = uiModule;

describe('trajectory controls', () => {
    let viewer;
    let controls;
    let frames;
    let current;
    let animated;
    const settle = async () => {
        await Promise.resolve();
        await Promise.resolve();
    };
    const button = label => controls.element.querySelector(`button[aria-label="${label}"]`);
    const scrub = frame => {
        const slider = controls.element.querySelector('input[type="range"]');
        slider.value = String(frame);
        slider.dispatchEvent(new Event('input', { bubbles: true }));
    };
    const pointer = (target, type, clientX, clientY, extra = {}) => {
        const event = new Event(type, { bubbles: true, cancelable: true });
        Object.assign(event, { pointerId: 1, button: 0, isPrimary: true, clientX, clientY, ...extra });
        target.dispatchEvent(event);
        return event;
    };
    const moveKey = (key, shiftKey = false) => {
        controls.element.dispatchEvent(new dom.window.KeyboardEvent('keydown', {
            key, shiftKey, bubbles: true, cancelable: true
        }));
    };
    const geometry = (scale = 1) => {
        const root = document.createElement('div');
        root.className = 'molui';
        document.body.append(root);
        root.append(controls.element);
        const bounds = { left: 100, top: 50, width: 600, height: 400 };
        const rect = ({ left, top, width, height }) => ({ left, top, width, height,
            x: left, y: top, right: left + width, bottom: top + height });
        Object.defineProperties(root, {
            clientWidth: { get: () => bounds.width },
            clientHeight: { get: () => bounds.height }
        });
        root.getBoundingClientRect = () => rect({ ...bounds, width: bounds.width * scale, height: bounds.height * scale });
        controls.element.getBoundingClientRect = () => rect({
            left: bounds.left + (parseFloat(controls.element.style.left) || 240) * scale,
            top: bounds.top + (parseFloat(controls.element.style.top) || 290) * scale,
            width: Math.min(root.classList.contains('molui-trajectory-floating') ? 320 : 280, bounds.width - 16) * scale,
            height: Math.min(84, bounds.height - 16) * scale
        });
        const widget = controls.element;
        const counter = widget.querySelector('output');
        counter.getBoundingClientRect = () => {
            const bar = controls.element.getBoundingClientRect();
            return rect({ left: bar.right - 105 * scale, top: bar.top + 6 * scale, width: 90 * scale, height: 30 * scale });
        };
        widget.setPointerCapture = jest.fn();
        widget.releasePointerCapture = jest.fn();
        const start = () => {
            const box = counter.getBoundingClientRect();
            pointer(counter, 'pointerdown', box.left + 6 * scale, box.top + 6 * scale);
        };
        const moveTo = (left, top) => pointer(document, 'pointermove',
            bounds.left + (left + Math.min(320, bounds.width - 16) - 105 + 6) * scale,
            bounds.top + (top + 12) * scale);
        return { root, bounds, widget, counter, start, moveTo };
    };

    beforeEach(() => {
        jest.useFakeTimers();
        Object.defineProperty(document, 'hidden', { configurable: true, value: false });
        frames = 5;
        current = 0;
        animated = false;
        const canvas = document.createElement('canvas');
        document.body.append(canvas);
        viewer = {
            getCanvas: () => canvas,
            getNumFrames: () => frames,
            getFrame: () => current,
            isAnimated: () => animated,
            stopAnimate: jest.fn(() => { animated = false; }),
            setFrame: jest.fn(frame => { current = frame; return Promise.resolve(); }),
            render: jest.fn()
        };
        controls = new TrajectoryControls(viewer);
        document.body.append(controls.element);
    });

    afterEach(() => {
        controls.destroy();
        document.body.replaceChildren();
        delete document.hidden;
        jest.useRealTimers();
    });

    test('appears only for multiple frames and never starts playback on load', () => {
        expect(controls.element.hidden).toBe(false);
        expect(button('First frame').disabled).toBe(true);
        expect(viewer.setFrame).not.toHaveBeenCalled();
        frames = 1;
        controls.refresh();
        expect(controls.element.hidden).toBe(true);
        frames = 7;
        controls.refresh();
        expect(controls.element.hidden).toBe(false);
        expect(controls.element.querySelector('output').textContent).toBe('Frame 1 / 7');
    });

    test('steps and seeks with one-based labels and correct boundary controls', async () => {
        button('Next frame').click();
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(1);
        expect(viewer.render).toHaveBeenCalledTimes(1);
        expect(controls.element.querySelector('output').textContent).toBe('Frame 2 / 5');
        button('Last frame').click();
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(4);
        expect(button('Next frame').disabled).toBe(true);
        expect(button('First frame').disabled).toBe(false);
    });

    test('uses play/pause icons while keeping accessible names and pressed state', () => {
        const play = button('Play trajectory');
        expect(play.textContent).toBe('');
        expect(play.querySelector('svg').getAttribute('aria-hidden')).toBe('true');
        const playPath = play.querySelector('path').getAttribute('d');
        play.click();
        const pause = button('Pause trajectory');
        expect(pause).toBe(play);
        expect(pause.textContent).toBe('');
        expect(pause.querySelector('path').getAttribute('d')).not.toBe(playPath);
        expect(pause.title).toBe('Pause trajectory');
        expect(pause.getAttribute('aria-pressed')).toBe('true');
        pause.click();
        expect(play.getAttribute('aria-pressed')).toBe('false');
        expect(play.querySelector('path').getAttribute('d')).toBe(playPath);
    });

    test('clicking the widget does not undock or move it and needs no grab button', () => {
        const { root, widget, counter, start } = geometry();
        const box = counter.getBoundingClientRect();
        expect(widget.querySelector('.molui-trajectory-drag')).toBeNull();
        expect(widget.tabIndex).toBe(0);
        start();
        pointer(document, 'pointermove', box.left + 8, box.top + 7);
        pointer(document, 'pointerup', box.left + 8, box.top + 7);
        expect(root.classList.contains('molui-trajectory-floating')).toBe(false);
        expect(controls.element.style.left).toBe('');
        expect(widget.setPointerCapture).toHaveBeenCalledWith(1);
        expect(widget.releasePointerCapture).toHaveBeenCalledWith(1);
    });

    test('drags the frame counter, preserves its grab point as the bar widens, and blocks viewer gestures', () => {
        const { root, start, moveTo } = geometry();
        const viewerGesture = jest.fn();
        root.addEventListener('pointerdown', viewerGesture);
        start();
        expect(viewerGesture).not.toHaveBeenCalled();
        pointer(document, 'pointermove', 0, 0, { pointerId: 2 });
        expect(controls.element.style.left).toBe('');
        moveTo(80, 70);
        expect(controls.element.style.left).toBe('80px');
        expect(controls.element.style.top).toBe('70px');
        expect(root.classList.contains('molui-trajectory-floating')).toBe(true);
        expect(controls.element.classList.contains('molui-trajectory-dragging')).toBe(true);
        pointer(document, 'pointerup', 0, 0);
        expect(controls.element.classList.contains('molui-trajectory-dragging')).toBe(false);
        controls.refresh();
        expect(controls.element.style.left).toBe('80px');
        expect(viewer.setFrame).not.toHaveBeenCalled();
    });

    test.each(['', '.molui-trajectory-header', '.molui-trajectory-transport', '.molui-trajectory-options'])(
        'drags from the widget background or a control gap: %s', selector => {
            const { widget } = geometry();
            const target = selector ? widget.querySelector(selector) : widget;
            if (target !== widget) target.getBoundingClientRect = widget.getBoundingClientRect;
            const box = target.getBoundingClientRect();
            pointer(target, 'pointerdown', box.left + 6, box.top + 6);
            pointer(document, 'pointermove', box.left - 44, box.top - 34);
            pointer(document, 'pointerup', box.left - 44, box.top - 34);
            expect(widget.style.left).toBe('190px');
            expect(widget.style.top).toBe('250px');
            expect(document.activeElement).toBe(widget);
            expect(viewer.setFrame).not.toHaveBeenCalled();
        }
    );

    test('clamps the floating bar to viewer bounds and reclamps after resize or reappearance', () => {
        const { bounds, start, moveTo } = geometry();
        start();
        moveTo(10000, 10000);
        expect(controls.element.style.left).toBe('272px');
        expect(controls.element.style.top).toBe('308px');
        moveTo(-10000, -10000);
        expect(controls.element.style.left).toBe('8px');
        expect(controls.element.style.top).toBe('8px');
        moveTo(200, 200);
        pointer(document, 'pointerup', 0, 0);
        bounds.width = 350;
        bounds.height = 150;
        controls.resize();
        expect(controls.element.style.left).toBe('22px');
        expect(controls.element.style.top).toBe('58px');
        frames = 1;
        controls.refresh();
        bounds.width = 200;
        bounds.height = 60;
        controls.resize();
        expect(controls.element.style.left).toBe('22px');
        frames = 5;
        controls.refresh();
        expect(controls.element.style.left).toBe('8px');
        expect(controls.element.style.top).toBe('8px');
    });

    test('uses viewer coordinates for dragging, bounds, and keyboard movement in a scaled embed', () => {
        const { start, moveTo } = geometry(0.5);
        moveKey('ArrowLeft');
        expect(controls.element.style.left).toBe('230px');
        moveKey('Home');
        start();
        moveTo(80, 70);
        expect(controls.element.style.left).toBe('80px');
        expect(controls.element.style.top).toBe('70px');
        moveTo(10000, 10000);
        expect(controls.element.style.left).toBe('272px');
        expect(controls.element.style.top).toBe('308px');
        pointer(document, 'pointerup', 0, 0);
        moveKey('ArrowLeft');
        expect(controls.element.style.left).toBe('262px');
    });

    test('buttons, sliders, speed, and checkbox labels retain their native actions without dragging', async () => {
        const { root } = geometry();
        const slider = controls.element.querySelector('input[type="range"]');
        const loop = controls.element.querySelector('input[type="checkbox"]');
        const speed = controls.element.querySelector('select');
        for (const target of [slider, button('Next frame'), button('Play trajectory').querySelector('path'), speed, speed.parentElement, loop, loop.parentElement]) {
            expect(pointer(target, 'pointerdown', 300, 300).defaultPrevented).toBe(false);
            pointer(document, 'pointermove', 40, 40);
            pointer(document, 'pointerup', 40, 40);
        }
        expect(root.classList.contains('molui-trajectory-floating')).toBe(false);
        scrub(3);
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(2);
        button('Next frame').click();
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(3);
        loop.parentElement.click();
        expect(loop.checked).toBe(false);
    });

    test('touch adjustment to a nearby button still drags the physically pressed background', () => {
        const { widget, counter, moveTo } = geometry();
        const originalHitTest = document.elementFromPoint;
        document.elementFromPoint = jest.fn(() => counter);
        try {
            const box = counter.getBoundingClientRect();
            const event = pointer(button('Next frame'), 'pointerdown', box.left + 6, box.top + 6, { pointerType: 'touch' });
            expect(document.elementFromPoint).toHaveBeenCalledWith(box.left + 6, box.top + 6);
            expect(event.defaultPrevented).toBe(true);
            moveTo(80, 70);
            expect(widget.style.left).toBe('80px');
            expect(widget.style.top).toBe('70px');
            pointer(document, 'pointerup', 0, 0);
            expect(viewer.setFrame).not.toHaveBeenCalled();
        } finally { document.elementFromPoint = originalHitTest; }
    });

    test('physical hits on controls keep their native behavior even with a different event target', () => {
        const { widget } = geometry();
        const originalHitTest = document.elementFromPoint;
        document.elementFromPoint = jest.fn(() => button('Next frame'));
        try {
            const event = pointer(widget, 'pointerdown', 300, 300, { pointerType: 'touch' });
            expect(event.defaultPrevented).toBe(false);
            pointer(document, 'pointermove', 100, 100);
            pointer(document, 'pointerup', 100, 100);
            expect(widget.setPointerCapture).not.toHaveBeenCalled();
            expect(widget.classList.contains('molui-trajectory-moved')).toBe(false);
        } finally { document.elementFromPoint = originalHitTest; }
    });

    test('supports keyboard movement and restores docking with Home or double-click', () => {
        const { root, counter } = geometry();
        moveKey('ArrowLeft');
        expect(controls.element.style.left).toBe('230px');
        expect(controls.element.style.top).toBe('290px');
        moveKey('ArrowUp', true);
        expect(controls.element.style.top).toBe('250px');
        moveKey('Home');
        expect(controls.element.style.left).toBe('');
        expect(controls.element.style.top).toBe('');
        expect(root.classList.contains('molui-trajectory-floating')).toBe(false);
        moveKey('ArrowUp');
        counter.dispatchEvent(new Event('dblclick', { bubbles: true, cancelable: true }));
        expect(controls.element.classList.contains('molui-trajectory-moved')).toBe(false);
        expect(root.classList.contains('molui-trajectory-floating')).toBe(false);
    });

    test('keyboard and double-click events from playback inputs do not move or reset the widget', () => {
        const { widget } = geometry();
        moveKey('ArrowLeft');
        for (const target of widget.querySelectorAll('input, select, button, label')) {
            for (const key of ['ArrowLeft', 'ArrowRight', 'Home']) {
                const event = new dom.window.KeyboardEvent('keydown', { key, bubbles: true, cancelable: true });
                target.dispatchEvent(event);
                expect(event.defaultPrevented).toBe(false);
            }
            const event = new Event('dblclick', { bubbles: true, cancelable: true });
            target.dispatchEvent(event);
            expect(event.defaultPrevented).toBe(false);
        }
        expect(widget.style.left).toBe('230px');
        expect(widget.classList.contains('molui-trajectory-moved')).toBe(true);
    });

    test('pointer cancellation and lost capture end a drag without changing playback', () => {
        const { widget, start, moveTo } = geometry();
        button('Play trajectory').click();
        start();
        moveTo(80, 70);
        pointer(document, 'pointercancel', 0, 0);
        moveTo(180, 170);
        expect(controls.element.style.left).toBe('80px');
        expect(controls.isPlaying).toBe(true);
        start();
        pointer(widget, 'lostpointercapture', 0, 0);
        moveTo(180, 170);
        expect(controls.element.style.left).toBe('80px');
        expect(widget.releasePointerCapture).toHaveBeenCalledTimes(2);
    });

    test('losing a child implicit capture does not end the widget drag', () => {
        const { widget, counter, start, moveTo } = geometry();
        start();
        pointer(counter, 'lostpointercapture', 0, 0);
        moveTo(80, 70);
        expect(widget.style.left).toBe('80px');
        expect(widget.releasePointerCapture).not.toHaveBeenCalled();
        pointer(widget, 'lostpointercapture', 0, 0);
        moveTo(180, 170);
        expect(widget.style.left).toBe('80px');
        expect(widget.releasePointerCapture).toHaveBeenCalledTimes(1);
    });

    test('destroy during a drag removes document listeners, capture, and floating layout state', () => {
        const { root, widget, start, moveTo } = geometry();
        const added = jest.spyOn(document, 'addEventListener');
        const removed = jest.spyOn(document, 'removeEventListener');
        start();
        moveTo(80, 70);
        const listeners = added.mock.calls.filter(([type]) => type.startsWith('pointer'));
        controls.destroy();
        for (const listener of listeners) expect(removed).toHaveBeenCalledWith(...listener);
        expect(listeners).toHaveLength(3);
        expect(widget.releasePointerCapture).toHaveBeenCalledWith(1);
        expect(root.classList.contains('molui-trajectory-floating')).toBe(false);
        moveTo(180, 170);
        expect(controls.element.style.left).toBe('');
        added.mockRestore();
        removed.mockRestore();
    });

    test('serializes slow loads and coalesces rapid scrub requests', async () => {
        frames = 10;
        controls.refresh();
        const complete = [];
        viewer.setFrame.mockImplementation(frame => new Promise(resolve => {
            complete.push(() => { current = frame; resolve(); });
        }));
        scrub(3);
        scrub(6);
        scrub(9);
        expect(viewer.setFrame.mock.calls).toEqual([[2]]);
        complete.shift()();
        await settle();
        expect(viewer.setFrame.mock.calls).toEqual([[2], [8]]);
        complete.shift()();
        await settle();
        expect(controls.element.querySelector('output').textContent).toBe('Frame 9 / 10');
        expect(controls.element.getAttribute('aria-busy')).toBe('false');
    });

    test('plays at the selected speed and stops at the end when loop is off', async () => {
        frames = 3;
        controls.refresh();
        controls.element.querySelector('input[type="checkbox"]').checked = false;
        const speed = controls.element.querySelector('select');
        speed.value = '5';
        speed.dispatchEvent(new Event('change'));
        button('Play trajectory').click();
        jest.advanceTimersByTime(199);
        await settle();
        expect(viewer.setFrame).not.toHaveBeenCalled();
        jest.advanceTimersByTime(1);
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(1);
        jest.advanceTimersByTime(200);
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(2);
        jest.advanceTimersByTime(200);
        await settle();
        expect(button('Play trajectory')).not.toBeNull();
        expect(viewer.setFrame).toHaveBeenCalledTimes(2);
    });

    test('loops and pauses before manual stepping', async () => {
        current = 4;
        controls.refresh();
        button('Play trajectory').click();
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(0);
        button('Next frame').click();
        await settle();
        expect(viewer.setFrame).toHaveBeenLastCalledWith(1);
        jest.advanceTimersByTime(1000);
        await settle();
        expect(viewer.setFrame).toHaveBeenCalledTimes(2);
    });

    test('resumes from a pending seek rather than the previously displayed last frame', async () => {
        current = 4;
        controls.refresh();
        let complete;
        viewer.setFrame.mockImplementationOnce(frame => new Promise(resolve => {
            complete = () => { current = frame; resolve(); };
        }));
        scrub(2);
        button('Play trajectory').click();
        complete();
        await settle();
        jest.advanceTimersByTime(100);
        await settle();
        expect(viewer.setFrame.mock.calls).toEqual([[1], [2]]);
    });

    test('reflects external animation and can pause it', () => {
        animated = true;
        current = 2;
        controls.refresh();
        expect(controls.element.querySelector('output').textContent).toBe('Frame 3 / 5');
        button('Pause trajectory').click();
        expect(viewer.stopAnimate).toHaveBeenCalledTimes(1);
        expect(button('Play trajectory')).not.toBeNull();
    });

    test('pauses when the page is hidden', async () => {
        button('Play trajectory').click();
        Object.defineProperty(document, 'hidden', { configurable: true, value: true });
        document.dispatchEvent(new Event('visibilitychange'));
        jest.advanceTimersByTime(500);
        await settle();
        expect(viewer.setFrame).not.toHaveBeenCalled();
        expect(controls.isPlaying).toBe(false);
    });

    test('surfaces a failed load and releases the controls for retry', async () => {
        viewer.setFrame.mockRejectedValueOnce(new Error('Connection lost'));
        scrub(3);
        await settle();
        expect(controls.element.querySelector('[role="alert"]').textContent).toContain('Connection lost');
        expect(controls.element.getAttribute('aria-busy')).toBe('false');
        button('Next frame').click();
        await settle();
        expect(viewer.render).toHaveBeenCalledTimes(1);
    });

    test('destroy cancels timers and prevents a late frame from rendering', async () => {
        let complete;
        viewer.setFrame.mockImplementation(() => new Promise(resolve => { complete = resolve; }));
        scrub(3);
        controls.destroy();
        complete();
        await settle();
        jest.advanceTimersByTime(1000);
        expect(viewer.render).not.toHaveBeenCalled();
        expect(jest.getTimerCount()).toBe(0);
        expect(controls.element.isConnected).toBe(false);
    });

    test('releases its polling timer when mounted controls are removed', () => {
        controls.refresh();
        controls.element.remove();
        jest.advanceTimersByTime(250);
        expect(jest.getTimerCount()).toBe(0);
    });
});
