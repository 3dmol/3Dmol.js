import type { GLViewer } from "../GLViewer";

export interface TrajectorySettings {
    fps: number;
    loop: boolean;
    playing: boolean;
    /** Zero-based index; negative and overflowing inputs wrap like GLViewer.setFrame(). */
    frame: number;
}

type TransportIcon = "first" | "previous" | "play" | "pause" | "next" | "last";

const iconPaths: Record<TransportIcon, string> = {
    first: "M6 5v14M19 5l-9 7 9 7V5Z",
    previous: "M15 5l-9 7 9 7",
    play: "M8 5v14l11-7Z",
    pause: "M6 5h4v14H6ZM14 5h4v14h-4Z",
    next: "M9 5l9 7-9 7",
    last: "M18 5v14M5 5l9 7-9 7V5Z"
};

/** Accessible transport controls for the frames of all models in a viewer. */
export class TrajectoryControls {
    readonly element: HTMLElement;

    private readonly document: Document;
    private readonly playButton: HTMLButtonElement;
    private readonly firstButton: HTMLButtonElement;
    private readonly previousButton: HTMLButtonElement;
    private readonly nextButton: HTMLButtonElement;
    private readonly lastButton: HTMLButtonElement;
    private readonly slider: HTMLInputElement;
    private readonly frameCount: HTMLOutputElement;
    private readonly speed: HTMLSelectElement;
    private readonly loop: HTMLInputElement;
    private readonly error: HTMLParagraphElement;
    private readonly refreshTimer: ReturnType<typeof setInterval>;
    private playbackTimer?: ReturnType<typeof setTimeout>;
    private playing = false;
    private destroyed = false;
    private wasConnected = false;
    private pending = false;
    private queuedFrame?: number;
    private frameLoad?: Promise<void>;
    private settingsRevision = 0;
    private totalFrames = 0;
    private currentFrame = 0;
    private position?: { left: number; top: number };
    private floatingRoot?: HTMLElement;
    private drag?: { pointerId: number; startX: number; startY: number; anchor: Element; grabX: number; grabY: number; moving: boolean };

    constructor(private readonly viewer: GLViewer, private readonly onError?: (message: string) => void) {
        this.document = viewer.getCanvas().ownerDocument;
        const document = this.document;
        this.element = document.createElement("section");
        this.element.className = "molui-trajectory";
        this.element.setAttribute("aria-label", "Trajectory controls");
        this.element.tabIndex = 0;
        const moveInstructions = "Drag to move. When the widget is focused, arrow keys move; Shift moves farther. Home or double-click the background resets position.";
        this.element.setAttribute("aria-description", moveInstructions);
        this.element.setAttribute("aria-keyshortcuts", "ArrowUp ArrowDown ArrowLeft ArrowRight Home");
        this.element.addEventListener("pointerdown", this.onPointerDown);
        this.element.addEventListener("lostpointercapture", this.onPointerEnd);
        this.element.addEventListener("keydown", this.onMoveKey);
        this.element.addEventListener("dblclick", event => {
            if (this.isControlTarget(event.target)) return;
            event.preventDefault();
            event.stopPropagation();
            this.resetPosition();
        });

        const heading = document.createElement("div");
        heading.className = "molui-trajectory-header";
        heading.title = moveInstructions;
        this.frameCount = document.createElement("output");
        this.frameCount.className = "molui-frame-count";
        // Frame updates must not interrupt a screen reader on every animation tick.
        this.frameCount.setAttribute("aria-live", "off");
        heading.append(this.frameCount);

        const transport = document.createElement("div");
        transport.className = "molui-trajectory-transport";
        this.firstButton = this.button("first", "First frame", () => this.seek(0));
        this.previousButton = this.button("previous", "Previous frame", () => this.seek(this.requestedFrame - 1));
        this.playButton = this.button("play", "Play trajectory", () => this.togglePlayback());
        this.playButton.classList.add("molui-primary");
        this.nextButton = this.button("next", "Next frame", () => this.seek(this.requestedFrame + 1));
        this.lastButton = this.button("last", "Last frame", () => this.seek(this.totalFrames - 1));
        transport.append(this.firstButton, this.previousButton, this.playButton, this.nextButton, this.lastButton);

        this.slider = document.createElement("input");
        this.slider.className = "molui-frame-slider";
        this.slider.type = "range";
        this.slider.min = "1";
        this.slider.step = "1";
        this.slider.setAttribute("aria-label", "Trajectory frame");
        this.slider.addEventListener("input", () => this.seek(Number(this.slider.value) - 1));

        const options = document.createElement("div");
        options.className = "molui-trajectory-options";
        const speedLabel = document.createElement("label");
        speedLabel.textContent = "Speed ";
        this.speed = document.createElement("select");
        this.speed.className = "molui-playback-speed";
        this.speed.setAttribute("aria-label", "Playback speed");
        for (const fps of [1, 5, 10, 15, 24, 30, 60]) {
            const option = document.createElement("option");
            option.value = String(fps);
            option.textContent = `${fps} fps`;
            this.speed.append(option);
        }
        this.speed.value = "10";
        this.speed.addEventListener("change", () => {
            if (this.playing && !this.pending) this.scheduleFrame();
        });
        speedLabel.append(this.speed);
        const loopLabel = document.createElement("label");
        loopLabel.className = "molui-check";
        this.loop = document.createElement("input");
        this.loop.type = "checkbox";
        this.loop.checked = true;
        loopLabel.append(this.loop, document.createTextNode(" Loop"));
        options.append(speedLabel, loopLabel);

        this.error = document.createElement("p");
        this.error.className = "molui-trajectory-error";
        this.error.setAttribute("role", "alert");
        this.error.hidden = true;
        this.element.append(heading, transport, this.slider, options, this.error);
        this.document.addEventListener("visibilitychange", this.onVisibilityChange);
        // Also reflects frames added and animation started by the host application.
        this.refreshTimer = setInterval(() => this.refresh(), 250);
        this.refresh();
    }

    get isPlaying(): boolean {
        return this.playing || this.viewer.isAnimated();
    }

    /** Read transport settings without exposing or depending on its DOM. */
    getSettings(): TrajectorySettings {
        return {
            fps: Number(this.speed.value), loop: this.loop.checked,
            playing: this.isPlaying, frame: this.clamp(this.viewer.getFrame())
        };
    }

    /** Restore transport settings after models have been loaded. */
    async setSettings(settings: Partial<TrajectorySettings>): Promise<void> {
        if (this.destroyed) return;
        if (settings.fps !== undefined && (!Number.isFinite(settings.fps) || settings.fps <= 0 || settings.fps > 240)) {
            throw new Error('Playback speed must be between zero and 240 frames per second.');
        }
        if (settings.frame !== undefined && !Number.isInteger(settings.frame)) {
            throw new Error('Choose an integer trajectory frame.');
        }
        const revision = ++this.settingsRevision;
        const playing = settings.playing ?? this.isPlaying;
        this.pause();
        this.refresh();
        if (settings.fps !== undefined) {
            const value = String(settings.fps);
            if (!Array.from(this.speed.options).some(option => option.value === value)) {
                const option = this.document.createElement('option');
                option.value = value; option.textContent = `${value} fps`;
                this.speed.append(option);
            }
            this.speed.value = value;
        }
        if (settings.loop !== undefined) this.loop.checked = settings.loop;
        if (settings.frame !== undefined) await this.queueFrame(settings.frame, true);
        else if (this.pending) await this.frameLoad;
        if (this.destroyed || revision !== this.settingsRevision) return;
        if (!this.error.hidden) throw new Error(this.error.textContent || 'Unable to load trajectory frame.');
        if (playing) this.togglePlayback();
        else this.updateControls();
    }

    private get requestedFrame(): number {
        return this.queuedFrame ?? (this.pending ? Number(this.slider.value) - 1 : this.currentFrame);
    }

    private button(icon: TransportIcon, label: string, action?: () => void): HTMLButtonElement {
        const button = this.document.createElement("button");
        button.type = "button";
        button.className = "molui-button";
        this.setIcon(button, icon);
        button.title = label;
        button.setAttribute("aria-label", label);
        if (action) button.addEventListener("click", action);
        return button;
    }

    private setIcon(button: HTMLButtonElement, icon: TransportIcon): void {
        if (button.dataset.icon === icon) return;
        const svg = this.document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.setAttribute("viewBox", "0 0 24 24");
        svg.setAttribute("width", "16");
        svg.setAttribute("height", "16");
        svg.setAttribute("aria-hidden", "true");
        svg.setAttribute("focusable", "false");
        svg.setAttribute("fill", icon === "play" || icon === "pause" ? "currentColor" : "none");
        svg.setAttribute("stroke", "currentColor");
        svg.setAttribute("stroke-width", "1.8");
        svg.setAttribute("stroke-linecap", "round");
        svg.setAttribute("stroke-linejoin", "round");
        const path = this.document.createElementNS("http://www.w3.org/2000/svg", "path");
        path.setAttribute("d", iconPaths[icon]);
        svg.append(path);
        button.replaceChildren(svg);
        button.dataset.icon = icon;
    }

    /** Keep a moved transport reachable when its viewer changes size. */
    resize(): void {
        if (!this.destroyed && this.position && !this.element.hidden) {
            this.place(this.position.left, this.position.top);
        }
    }

    /** Restore the default docked position. */
    resetPosition(): void {
        this.endDrag();
        this.position = undefined;
        this.floatingRoot?.classList.remove("molui-trajectory-floating");
        this.floatingRoot = undefined;
        this.element.classList.remove("molui-trajectory-moved");
        this.element.style.removeProperty("left");
        this.element.style.removeProperty("top");
    }

    private float(): void {
        this.floatingRoot = this.element.parentElement ?? undefined;
        this.floatingRoot?.classList.add("molui-trajectory-floating");
        this.element.classList.add("molui-trajectory-moved");
    }

    private layoutBounds() {
        const parent = this.element.parentElement;
        if (!parent) return undefined;
        const bounds = parent.getBoundingClientRect();
        // Pointer coordinates are viewport pixels; positions use the viewer's CSS
        // pixels, which differ when a host page scales its embedded viewer.
        const scaleX = bounds.width / (parent.offsetWidth || parent.clientWidth || bounds.width) || 1;
        const scaleY = bounds.height / (parent.offsetHeight || parent.clientHeight || bounds.height) || 1;
        return {
            left: bounds.left + parent.clientLeft * scaleX,
            top: bounds.top + parent.clientTop * scaleY,
            width: parent.clientWidth || bounds.width / scaleX,
            height: parent.clientHeight || bounds.height / scaleY,
            scaleX, scaleY
        };
    }

    private place(left: number, top: number): void {
        const bounds = this.layoutBounds();
        if (!bounds) return;
        const bar = this.element.getBoundingClientRect();
        // Hidden viewers have no layout; retain their chosen position until shown.
        if (!bounds.width || !bounds.height || !bar.width || !bar.height) return;
        const maxLeft = Math.max(0, bounds.width - bar.width / bounds.scaleX);
        const maxTop = Math.max(0, bounds.height - bar.height / bounds.scaleY);
        const insetX = Math.min(8, maxLeft / 2);
        const insetY = Math.min(8, maxTop / 2);
        this.position = {
            left: Math.max(insetX, Math.min(maxLeft - insetX, left)),
            top: Math.max(insetY, Math.min(maxTop - insetY, top))
        };
        this.element.style.left = `${this.position.left}px`;
        this.element.style.top = `${this.position.top}px`;
    }

    private isControlTarget(target: EventTarget | null): boolean {
        return !!(target as Element | null)?.closest?.("button, input, select, textarea, label, a, [contenteditable]");
    }

    private readonly onPointerDown = (event: PointerEvent): void => {
        if (this.destroyed || this.drag || event.button !== 0 || event.isPrimary === false) return;
        // Touch browsers may adjust a press on nearby background to target a
        // button. Use the physical hit so that background drags remain possible.
        const hit = this.document.elementFromPoint?.(event.clientX, event.clientY);
        const anchor = hit && this.element.contains(hit) ? hit : event.target as Element;
        if (this.isControlTarget(anchor)) return;
        event.preventDefault();
        event.stopPropagation();
        this.element.focus({ preventScroll: true });
        const box = anchor.getBoundingClientRect();
        this.drag = {
            pointerId: event.pointerId, startX: event.clientX, startY: event.clientY, anchor,
            grabX: event.clientX - box.left, grabY: event.clientY - box.top, moving: false
        };
        // Document listeners also cover browsers without pointer capture support.
        this.document.addEventListener("pointermove", this.onPointerMove, true);
        this.document.addEventListener("pointerup", this.onPointerEnd, true);
        this.document.addEventListener("pointercancel", this.onPointerEnd, true);
        this.document.defaultView?.addEventListener("blur", this.onDragBlur);
        try { this.element.setPointerCapture(event.pointerId); } catch { /* Capture can be unavailable. */ }
    };

    private readonly onPointerMove = (event: PointerEvent): void => {
        const drag = this.drag;
        if (!drag || event.pointerId !== drag.pointerId) return;
        event.preventDefault();
        event.stopPropagation();
        if (!drag.moving) {
            if (Math.hypot(event.clientX - drag.startX, event.clientY - drag.startY) < 4) return;
            drag.moving = true;
            this.float();
            this.element.classList.add("molui-trajectory-dragging");
        }
        const bounds = this.layoutBounds();
        if (bounds) {
            // Floating controls can be wider than the docked bar. Keep the clicked
            // text or background point under the pointer when that layout changes.
            const bar = this.element.getBoundingClientRect();
            const anchor = drag.anchor.getBoundingClientRect();
            this.place((event.clientX - bounds.left - (anchor.left - bar.left) - drag.grabX) / bounds.scaleX,
                (event.clientY - bounds.top - (anchor.top - bar.top) - drag.grabY) / bounds.scaleY);
        }
    };

    private readonly onPointerEnd = (event: PointerEvent): void => {
        if (event.pointerId !== this.drag?.pointerId) return;
        // A child can lose implicit touch capture as it transfers to the widget.
        if (event.type === "lostpointercapture" && event.target !== this.element) return;
        event.stopPropagation();
        this.endDrag();
    };

    private readonly onDragBlur = (): void => this.endDrag();

    private endDrag(): void {
        const drag = this.drag;
        this.drag = undefined;
        this.document.removeEventListener("pointermove", this.onPointerMove, true);
        this.document.removeEventListener("pointerup", this.onPointerEnd, true);
        this.document.removeEventListener("pointercancel", this.onPointerEnd, true);
        this.document.defaultView?.removeEventListener("blur", this.onDragBlur);
        this.element.classList.remove("molui-trajectory-dragging");
        if (drag) {
            try { this.element.releasePointerCapture(drag.pointerId); } catch { /* Capture may already be released. */ }
        }
    }

    private readonly onMoveKey = (event: KeyboardEvent): void => {
        if (this.destroyed || event.target !== this.element) return;
        if (event.key === "Home") {
            event.preventDefault();
            event.stopPropagation();
            this.resetPosition();
            return;
        }
        const direction: Record<string, [number, number]> = {
            ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1]
        };
        const delta = direction[event.key];
        if (!delta || event.altKey || event.ctrlKey || event.metaKey) return;
        event.preventDefault();
        event.stopPropagation();
        const bounds = this.layoutBounds();
        if (!bounds) return;
        const bar = this.element.getBoundingClientRect();
        const left = this.position?.left ?? (bar.left - bounds.left) / bounds.scaleX;
        const top = this.position?.top ?? (bar.top - bounds.top) / bounds.scaleY;
        this.endDrag();
        this.float();
        const step = event.shiftKey ? 40 : 10;
        this.place(left + delta[0] * step, top + delta[1] * step);
    };

    /** Reconcile the controls with models and externally selected frames. */
    refresh(): void {
        if (this.destroyed) return;
        if (this.element.isConnected) this.wasConnected = true;
        else if (this.wasConnected) {
            // Removing a mounted viewer must not retain it through the poll timer.
            this.destroy();
            return;
        }
        const count = Number(this.viewer.getNumFrames());
        this.totalFrames = Number.isFinite(count) ? Math.max(0, Math.floor(count)) : 0;
        this.element.hidden = this.totalFrames <= 1;
        if (this.playing && (this.totalFrames <= 1 || this.document.hidden || !this.element.isConnected)) {
            this.stopPlayback();
        }
        // An externally started viewer animation owns frame changes until paused.
        if (this.playing && this.viewer.isAnimated()) this.stopPlayback();
        if (!this.pending) this.currentFrame = this.clamp(this.viewer.getFrame());
        this.updateControls();
        this.resize();
    }

    /** Pause playback, including an animation started through viewer.animate(). */
    pause(): void {
        this.stopPlayback();
        if (this.viewer.isAnimated()) this.viewer.stopAnimate();
        this.updateControls();
    }

    private stopPlayback(): void {
        this.playing = false;
        if (this.playbackTimer !== undefined) clearTimeout(this.playbackTimer);
        this.playbackTimer = undefined;
    }

    private clamp(frame: number): number {
        return Math.max(0, Math.min(this.totalFrames - 1, Math.floor(Number(frame) || 0)));
    }

    private updateControls(): void {
        const available = this.totalFrames > 1;
        const requested = this.clamp(this.requestedFrame);
        const active = available && this.isPlaying;
        this.slider.max = String(Math.max(1, this.totalFrames));
        this.slider.value = String(requested + 1);
        this.slider.disabled = !available;
        this.slider.setAttribute("aria-valuetext", `Frame ${requested + 1} of ${this.totalFrames}`);
        this.frameCount.textContent = this.totalFrames ? `Frame ${this.currentFrame + 1} / ${this.totalFrames}` : "No frames";
        this.element.setAttribute("aria-busy", String(this.pending));
        this.setIcon(this.playButton, active ? "pause" : "play");
        this.playButton.title = active ? "Pause trajectory" : "Play trajectory";
        this.playButton.setAttribute("aria-label", this.playButton.title);
        this.playButton.setAttribute("aria-pressed", String(active));
        this.playButton.disabled = !available;
        this.firstButton.disabled = this.previousButton.disabled = !available || requested <= 0;
        this.nextButton.disabled = this.lastButton.disabled = !available || requested >= this.totalFrames - 1;
        this.speed.disabled = this.loop.disabled = !available;
    }

    private togglePlayback(): void {
        if (this.isPlaying) {
            this.pause();
            return;
        }
        this.refresh();
        if (this.totalFrames <= 1 || this.document.hidden) return;
        this.playing = true;
        this.error.hidden = true;
        if (this.requestedFrame >= this.totalFrames - 1) {
            this.queueFrame(0);
        } else if (!this.pending) {
            this.scheduleFrame();
        }
        this.updateControls();
    }

    private seek(frame: number): void {
        this.pause();
        if (this.destroyed || this.totalFrames <= 1) return;
        this.queueFrame(frame);
    }

    private queueFrame(frame: number, wrap = false): Promise<void> {
        // Explicit settings use the viewer's modulo semantics; pointer/keyboard
        // scrubbing stays within the transport's displayed bounds.
        this.queuedFrame = wrap ? frame : this.clamp(frame);
        const shownFrame = wrap && this.totalFrames ? ((frame % this.totalFrames) + this.totalFrames) % this.totalFrames : this.clamp(frame);
        this.slider.value = String(shownFrame + 1);
        if (!this.pending) this.frameLoad = this.loadQueuedFrames();
        return this.frameLoad ?? Promise.resolve();
    }

    private async loadQueuedFrames(): Promise<void> {
        this.pending = true;
        this.error.hidden = true;
        try {
            // Keep only the latest scrub request while a remote frame is loading.
            // Rendering one frame completes before another request is issued.
            while (this.queuedFrame !== undefined && !this.destroyed) {
                const frame = this.queuedFrame;
                this.queuedFrame = undefined;
                this.updateControls();
                await this.viewer.setFrame(frame);
                if (this.destroyed) return;
                this.currentFrame = this.clamp(this.viewer.getFrame());
                this.viewer.render();
            }
        } catch (error) {
            this.stopPlayback();
            this.queuedFrame = undefined;
            if (!this.destroyed) {
                const detail = error instanceof Error ? error.message : String(error);
                const message = `Unable to load trajectory frame: ${detail}`;
                this.error.textContent = message;
                this.error.hidden = false;
                this.onError?.(message);
            }
        } finally {
            this.pending = false;
            if (!this.destroyed) {
                this.refresh();
                if (this.playing) this.scheduleFrame();
            }
        }
    }

    private scheduleFrame(): void {
        if (this.playbackTimer !== undefined) clearTimeout(this.playbackTimer);
        this.playbackTimer = setTimeout(() => {
            this.playbackTimer = undefined;
            this.refresh();
            if (!this.playing || this.pending || this.destroyed) return;
            let next = this.currentFrame + 1;
            if (next >= this.totalFrames) {
                if (!this.loop.checked) {
                    this.stopPlayback();
                    this.updateControls();
                    return;
                }
                next = 0;
            }
            this.queueFrame(next);
        }, 1000 / Number(this.speed.value));
    }

    private readonly onVisibilityChange = (): void => {
        if (this.document.hidden) this.pause();
    };

    /** Release timers and listeners; an in-flight load will no longer render. */
    destroy(): void {
        if (this.destroyed) return;
        this.pause();
        this.destroyed = true;
        this.resetPosition();
        this.queuedFrame = undefined;
        clearInterval(this.refreshTimer);
        this.document.removeEventListener("visibilitychange", this.onVisibilityChange);
        this.element.remove();
    }
}
