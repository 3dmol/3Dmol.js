import type { StateManager } from './state';
import type { UIControlOptions } from './ui';
import { learningStyles } from './learning-styles';

/** The viewer's current classroom role. */
export type LearningRole = 'idle' | 'instructor' | 'student';

/** Optional host integration for URL state and application-owned scene loading. */
export interface LearningAdapter {
    getManager?: () => StateManager;
    getState?: () => any;
    applyState?: (snapshot: any, sections?: string[]) => Promise<void> | void;
    roleChanged?: (name: string | null, following: boolean) => void;
    subscribe?: (listener: () => void) => (() => void);
    showError?: (error: Error) => void;
}

/** Options for classroom controls in the optional UI bundle. */
export interface LearningConfig {
    /** Socket.IO server origin. Defaults to https://3dmol.csb.pitt.edu. */
    server?: string;
    /** Show a classroom button, or wait invisibly for an instructor. Defaults to button. */
    mode?: 'button' | 'hidden';
    /** Classroom name. Required in hidden mode; prefills the connection dialog in button mode. */
    session?: string;
    /** Hidden-mode instructor shortcut while the viewer is focused. Defaults to Ctrl+Alt+I. */
    hotkey?: string;
    /** Socket.IO browser-client script URL; use a local URL when required by CSP. */
    clientURL?: string;
    /** Optional Socket.IO-compatible factory, useful when the host already manages dependencies. */
    socketFactory?: (server: string, options: any) => any;
    /** Host application hooks; ordinary embedded viewers do not need an adapter. */
    adapter?: LearningAdapter;
}

const clientLoads = new WeakMap<Document, Map<string, Promise<any>>>();
function loadClient(document: Document, url: string): Promise<any> {
    const win = document.defaultView as any;
    if (win.io) return Promise.resolve(win.io);
    let loads = clientLoads.get(document);
    if (!loads) clientLoads.set(document, loads = new Map());
    if (!loads.has(url)) {
        const promise = new Promise((resolve, reject) => {
            const script = document.createElement('script'); script.src = url; script.async = true;
            const timer = setTimeout(() => finish(new Error('Loading the classroom connection timed out.')), 15000);
            const finish = (error?: Error) => {
                clearTimeout(timer); script.onload = script.onerror = null;
                if (error) { script.remove(); reject(error); }
                else if (win.io) resolve(win.io);
                else reject(new Error('The classroom client did not provide Socket.IO.'));
            };
            script.onload = () => finish(); script.onerror = () => finish(new Error('Unable to load the classroom connection.'));
            document.head.append(script);
        });
        loads.set(url, promise);
        promise.catch(() => loads.delete(url));
    }
    return loads.get(url);
}

/**
 * Classroom synchronization and atom queries for one viewer.
 * Instances own their socket, dialogs, shortcuts, and subscriptions. Dispose with
 * destroy() when removing the viewer. Hidden mode's shortcut is a convenience,
 * not instructor authentication: use a private session name/server when needed.
 */
export class LearningSession {
    /** Resolves after the Socket.IO client and event handlers are ready. */
    readonly ready: Promise<void>;
    /** The scoped classroom controls, mounted inside the viewer container. */
    readonly element: HTMLDivElement;
    private readonly container: HTMLElement;
    private readonly document: Document;
    private readonly window: Window;
    private readonly adapter: LearningAdapter;
    private readonly button: HTMLButtonElement;
    private readonly dialog: HTMLDivElement;
    private readonly heading: HTMLSpanElement;
    private readonly connectPanel: HTMLDivElement;
    private readonly input: HTMLInputElement;
    private readonly connectButton: HTMLButtonElement;
    private readonly status: HTMLDivElement;
    private readonly controls: HTMLDivElement;
    private readonly ask: HTMLButtonElement;
    private readonly users: HTMLSpanElement;
    private readonly responses: HTMLSpanElement;
    private readonly refreshButton: HTMLButtonElement;
    private readonly monitor: HTMLDivElement;
    private readonly errorText: HTMLParagraphElement;
    private socket: any;
    private role: LearningRole = 'idle';
    private name: string | null = null;
    private secret: number | null = null;
    private member: string | null = null;
    private pending: 'create' | 'join' | null = null;
    private pendingSent = false;
    private activation = 0;
    private instructorRequested = false;
    private disposed = false;
    private userActivated = false;
    private watching = true;
    private revision = 0;
    private lastSent: any;
    private inFlight: any;
    private received: any;
    private dirty = false;
    private full = false;
    private timer?: ReturnType<typeof setTimeout>;
    private poll?: ReturnType<typeof setInterval>;
    private lastView = '';
    private lastSettings = '';
    private query = false;
    private showingResults = false;
    private results: any[] = [];
    private choices = new Map<any, any>();
    private interactions = new Map<any, any>();
    private savedControls?: UIControlOptions;
    private savedHidden = false;
    private applyGeneration = 0;
    private pendingState?: { snapshot: any; sections?: string[] };
    private applying?: Promise<void>;
    private cleanups: Array<() => void> = [];
    private subscriptions: Array<() => void> = [];
    private oldTabIndex: string | null;
    private drag?: { x: number; y: number; pointer: number };

    constructor(private readonly initialState: StateManager, public readonly config: LearningConfig = {}) {
        if (config.mode === 'hidden' && !config.session?.trim()) throw new Error('Hidden classroom mode requires a session identifier.');
        const state = initialState;
        this.adapter = config.adapter || {};
        this.container = state.viewer.container;
        this.document = this.container.ownerDocument; this.window = this.document.defaultView;
        const node = <K extends keyof HTMLElementTagNameMap>(tag: K, text = '') => {
            const result = this.document.createElement(tag); result.textContent = text; return result;
        };
        const button = (text: string, action: () => void) => {
            const result = node('button', text); result.type = 'button';
            result.addEventListener('click', action); return result;
        };
        this.element = node('div'); this.element.className = 'molui-learning';
        const style = node('style'); style.textContent = learningStyles; this.element.append(style);
        this.button = button('', () => this.open()); this.button.className = 'molui-learning-launch';
        this.button.title = 'Classroom session'; this.button.setAttribute('aria-label', 'Classroom session');
        this.button.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="7" r="3"/><path d="M5 21v-3a7 7 0 0 1 14 0v3M4 5a3 3 0 0 0 0 6M20 5a3 3 0 0 1 0 6"/></svg>';
        this.dialog = node('div'); this.dialog.className = 'molui-learning-dialog'; this.dialog.hidden = true;
        const header = node('div'); header.className = 'molui-learning-header';
        this.heading = node('span', 'Classroom'); this.heading.className = 'molui-learning-drag'; this.heading.tabIndex = 0;
        this.heading.setAttribute('aria-label', 'Move classroom dialog');
        const close = button('×', () => { if (this.role !== 'idle' || this.pending) this.leave(); else this.hide(); });
        close.setAttribute('aria-label', 'Close or leave classroom'); header.append(this.heading, close);
        const body = node('div'); body.className = 'molui-learning-body';
        this.connectPanel = node('div'); const label = node('label', 'Session name');
        this.input = node('input'); this.input.type = 'text'; this.input.autocomplete = 'off'; this.input.value = config.session || '';
        this.input.setAttribute('aria-label', 'Session name'); label.append(this.input);
        this.connectButton = button('Create', () => this.connectByName()); this.connectButton.disabled = true;
        const connectRow = node('div'); connectRow.className = 'molui-learning-row molui-learning-connect'; connectRow.append(label, this.connectButton);
        this.connectPanel.append(connectRow);
        this.status = node('div'); this.status.setAttribute('role', 'status');
        this.input.addEventListener('input', () => this.checkName());
        this.input.addEventListener('keydown', event => { if (event.key === 'Enter') this.connectByName(); });
        this.controls = node('div'); this.controls.className = 'molui-learning-row';
        this.ask = button('Query Atoms', () => this.queryAction());
        this.users = node('span', 'Users: 0'); this.responses = node('span');
        const stats = node('div'); stats.className = 'molui-learning-stats'; stats.setAttribute('aria-live', 'polite'); stats.append(this.users, this.responses);
        this.refreshButton = button('↻', () => this.socket?.emit('query fetch', { name: this.name }));
        this.refreshButton.title = 'Refresh responses'; this.refreshButton.setAttribute('aria-label', 'Refresh responses'); this.refreshButton.hidden = true;
        this.controls.append(this.ask, stats, this.refreshButton);
        this.monitor = node('div', 'Following instructor'); this.monitor.setAttribute('role', 'status');
        this.errorText = node('p'); this.errorText.className = 'molui-learning-error'; this.errorText.hidden = true; this.errorText.setAttribute('role', 'alert');
        body.append(this.connectPanel, this.status, this.controls, this.monitor, this.errorText); this.dialog.append(header, body); this.element.append(this.button, this.dialog);
        this.container.append(this.element); this.renderRole();
        this.oldTabIndex = this.container.getAttribute('tabindex');
        if (this.oldTabIndex === null) this.container.tabIndex = 0;
        const listen = (target: EventTarget, type: string, handler: EventListener) => {
            target.addEventListener(type, handler); this.cleanups.push(() => target.removeEventListener(type, handler));
        };
        ['pointerdown', 'mousedown', 'touchstart', 'wheel', 'dblclick'].forEach(type => listen(this.element, type, event => event.stopPropagation()));
        listen(this.container, 'pointerdown', event => {
            if (!this.element.contains(event.target as Node) && !(event.target as HTMLElement).closest('input,select,button,textarea,[contenteditable]')) this.container.focus({ preventScroll: true });
        });
        listen(this.container, 'keydown', event => this.keydown(event as KeyboardEvent));
        listen(this.window, 'beforeunload', () => this.leave());
        listen(this.heading, 'pointerdown', event => {
            const e = event as PointerEvent; if (e.button !== 0) return;
            const bounds = this.element.getBoundingClientRect(); this.drag = { x: e.clientX - bounds.left, y: e.clientY - bounds.top, pointer: e.pointerId };
            e.preventDefault(); this.heading.focus();
        });
        listen(this.document, 'pointermove', event => {
            const e = event as PointerEvent; if (!this.drag || e.pointerId !== this.drag.pointer) return;
            const bounds = this.container.getBoundingClientRect();
            this.place(e.clientX - bounds.left - this.drag.x, e.clientY - bounds.top - this.drag.y); e.preventDefault();
        });
        ['pointerup', 'pointercancel'].forEach(type => listen(this.document, type, () => { this.drag = undefined; }));
        listen(this.window, 'blur', () => { this.drag = undefined; });
        listen(this.window, 'resize', () => {
            if (this.element.style.left) this.place(this.element.offsetLeft, this.element.offsetTop);
        });
        listen(this.heading, 'keydown', event => {
            const e = event as KeyboardEvent, delta = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1] }[e.key];
            if (!delta || e.ctrlKey || e.altKey || e.metaKey) return;
            e.preventDefault(); e.stopPropagation();
            this.place(this.element.offsetLeft + delta[0] * (e.shiftKey ? 1 : 10), this.element.offsetTop + delta[1] * (e.shiftKey ? 1 : 10));
        });
        this.ready = this.connect(); this.ready.catch(error => this.showError(error));
    }

    private get state(): StateManager { return this.adapter.getManager?.() || this.initialState; }
    private get viewer() { return this.initialState.viewer; }
    private room(): string {
        if (this.config.session) return this.config.session;
        return this.input.value.trim();
    }
    private async connect(): Promise<void> {
        const io = this.config.socketFactory || await loadClient(this.document, this.config.clientURL || 'https://cdn.socket.io/4.8.1/socket.io.min.js');
        if (this.disposed) return;
        this.socket = (typeof io === 'function' ? io : io.connect)(this.config.server || 'https://3dmol.csb.pitt.edu', { forceNew: true, multiplex: false });
        const on = (event: string, callback: (message?: any) => void) => {
            const guarded = (message?: any) => { if (!this.disposed) return callback(message); };
            this.socket.on(event, guarded); this.cleanups.push(() => this.socket?.off?.(event, guarded));
        };
        on('connect', () => {
            this.errorText.hidden = true;
            if (this.config.mode === 'hidden' && this.watching) this.socket.emit('watch session event', { name: this.room() });
            if (this.role !== 'idle' && this.name) this.socket.emit('reconnection', { name: this.name, secret: this.secret, sid: this.member });
            this.sendPending();
            if (!this.dialog.hidden && this.role === 'idle') this.checkName();
        });
        on('disconnect', () => {
            if (this.pending) {
                const creating = this.pending === 'create';
                this.activation++; this.pending = null; this.pendingSent = false; this.name = null; this.renderRole();
                if (creating) this.showError(new Error('The classroom connection was interrupted. Try starting the classroom again.'));
            }
        });
        on('connect_error', error => {
            if (this.pending) {
                this.activation++; this.pending = null; this.pendingSent = false; this.name = null; this.renderRole();
            }
            this.showError(new Error(`Classroom connection failed: ${error?.message || error}`));
        });
        on('session available', message => {
            if (this.config.mode === 'hidden' && this.watching && !this.instructorRequested && message?.name === this.room() && this.role === 'idle' && !this.pending) this.join(message.name).catch(error => this.showError(error));
        });
        on('check session name response', message => {
            if (this.role !== 'idle' || this.pending) return;
            this.connectButton.disabled = !this.input.value.trim(); this.connectButton.textContent = message === 'exists' ? 'Join' : 'Create';
        });
        on('create session response', secret => {
            if (this.pending !== 'create' || !this.name) return;
            this.pending = null; this.pendingSent = false;
            if (!secret) {
                this.name = null; this.renderRole();
                this.showError(new Error('The classroom name is already in use. End the existing session or choose a different session name.'));
                return;
            }
            this.errorText.hidden = true;
            this.savedControls = this.state.getControls(); this.savedHidden = this.state.ui.element.hidden;
            this.state.setControls({ edit: true, trajectory: true }); this.state.ui.element.hidden = false;
            this.secret = Number(secret); this.role = 'instructor'; this.revision = 0; this.inFlight = undefined;
            this.subscriptions.push(this.state.subscribe(() => this.schedule()));
            if (this.viewer.subscribeStateChanges) this.subscriptions.push(this.viewer.subscribeStateChanges(() => this.schedule()));
            if (this.adapter.subscribe) this.subscriptions.push(this.adapter.subscribe(() => this.schedule()));
            this.poll = setInterval(() => this.pollScene(), 100); this.lastView = ''; this.lastSettings = ''; this.pollScene(); this.schedule(); this.renderRole();
        });
        on('join session response', member => {
            if (this.pending !== 'join' || !this.name) return;
            this.pending = null; this.pendingSent = false;
            if (!member || member === '0') { this.name = null; this.renderRole(); this.showError(new Error('The classroom is no longer available.')); return; }
            this.errorText.hidden = true;
            this.member = member; this.role = 'student'; this.revision = 0; this.received = undefined;
            this.state.setDistancePicking(false); this.savedControls = this.state.getControls(); this.savedHidden = this.state.ui.element.hidden;
            this.adapter.roleChanged?.(this.name, true);
            this.state.setControls({ edit: false, trajectory: false }); this.state.ui.element.hidden = true;
            this.saveInteractions(); this.disableInteractions(); this.renderRole();
        });
        on('viewer state revision', message => this.acknowledge(message));
        on('viewer state resync', message => { if (this.role === 'instructor') { this.revision = message.revision; this.inFlight = undefined; this.full = true; this.schedule(); } });
        on('reconnection response', message => {
            this.revision = message.revision;
            if (this.role === 'instructor' && message.controller) { this.inFlight = undefined; this.full = true; this.lastView = ''; this.schedule(); }
            else if (this.role === 'student') this.member = message.sid;
        });
        on('viewer view change response', view => { if (this.role === 'student') this.viewer.setView(view); });
        on('viewer state change response', message => this.receive(message));
        on('leave session response', () => this.reset());
        on('delete session response', () => this.reset());
        on('session count', count => { if (this.role === 'instructor') this.users.textContent = `Users: ${count}`; });
        on('query start response', () => {
            if (this.role === 'idle') return;
            this.query = true; this.showingResults = false; this.ask.textContent = 'Show Results'; this.responses.textContent = 'Responses: 0';
            if (this.role === 'student') { this.clearChoices(); this.restoreQuery(); }
        });
        on('query end response', () => { this.query = false; this.showingResults = false; this.ask.textContent = 'Query Atoms'; this.responses.textContent = ''; this.refreshButton.hidden = true; this.clearResults(); this.clearChoices(); this.restoreQuery(); });
        on('query update response', count => { if (this.role === 'instructor') this.responses.textContent = `Responses: ${count}`; });
        on('query fetch response', results => this.showResults(results));
        if (this.socket.connected && this.config.mode === 'hidden' && this.watching) this.socket.emit('watch session event', { name: this.room() });
    }

    /** Open the classroom connection dialog. */
    open(): void {
        if (this.disposed) return;
        this.userActivated = true; this.dialog.hidden = false; this.button.hidden = true; this.checkName();
    }
    /** Queue an instructor creation request; omitted names use the configured session. */
    async start(name = this.config.session || this.room()): Promise<void> {
        if (this.disposed || this.role !== 'idle' || this.pending) return;
        if (!name.trim()) throw new Error('A classroom name is required.');
        const activation = ++this.activation;
        this.userActivated = true; this.instructorRequested = true;
        this.pending = 'create'; this.pendingSent = false; this.name = name; this.lastSent = undefined; this.renderRole();
        try {
            await this.ready;
            if (this.disposed || activation !== this.activation) return;
            this.lastSent = this.snapshot();
            if (this.config.mode === 'hidden') {
                this.watching = true; this.socket.emit('watch session event', { name: this.room() });
            }
            this.sendPending();
        } catch (error) {
            if (this.disposed || activation !== this.activation) return;
            this.pending = null; this.pendingSent = false; this.name = null; this.renderRole(); throw error;
        }
    }
    /** Queue a request to join an existing instructor session. */
    async join(name = this.config.session || this.room()): Promise<void> {
        if (this.disposed || this.role !== 'idle' || this.pending) return;
        if (!name.trim()) throw new Error('A classroom name is required.');
        const activation = ++this.activation;
        this.instructorRequested = false; this.pending = 'join'; this.pendingSent = false; this.name = name; this.renderRole();
        try {
            await this.ready;
            if (this.disposed || activation !== this.activation) return;
            this.sendPending();
        } catch (error) {
            if (this.disposed || activation !== this.activation) return;
            this.pending = null; this.pendingSent = false; this.name = null; this.renderRole(); throw error;
        }
    }
    private sendPending(): void {
        if (!this.socket || this.socket.connected === false || !this.pending || this.pendingSent) return;
        if (this.pending === 'create' && !this.lastSent) return;
        this.pendingSent = true;
        if (this.pending === 'create') this.socket.emit('create session event', { name: this.name, state: this.lastSent, view: this.viewer.getView() });
        else this.socket.emit('join session event', { name: this.name });
    }
    /** Leave a classroom; instructors end it for all connected students. */
    leave(): void {
        this.watching = false;
        if (this.config.mode === 'hidden') this.socket?.emit('unwatch session event', { name: this.room() });
        if (this.name && (this.role !== 'idle' || this.pendingSent)) this.socket?.emit(this.role === 'instructor' || this.pending === 'create' ? 'delete session event' : 'leave session event', { name: this.name });
        this.reset();
    }
    /** Return idle, instructor, or student. */
    getRole(): LearningRole { return this.role; }
    /** Permanently release sockets, subscriptions, DOM, shortcuts, and query markers. */
    destroy(): void {
        if (this.disposed) return;
        this.leave(); this.disposed = true; this.cleanups.splice(0).forEach(cleanup => cleanup()); this.socket?.disconnect(); this.element.remove();
        if (this.oldTabIndex === null) this.container.removeAttribute('tabindex'); else this.container.setAttribute('tabindex', this.oldTabIndex);
    }
    private hide(): void { this.userActivated = false; this.dialog.hidden = true; this.button.hidden = this.config.mode === 'hidden'; }
    private checkName(): void {
        if (!this.socket || this.role !== 'idle' || this.pending) return;
        this.connectButton.disabled = true;
        if (this.input.value.trim()) this.socket.emit('check session name event', { name: this.input.value.trim() });
    }
    private connectByName(): void {
        if (this.connectButton.disabled) return;
        const action = this.connectButton.textContent === 'Join' ? this.join(this.input.value.trim()) : this.start(this.input.value.trim());
        action.catch(error => this.showError(error));
    }
    private renderRole(): void {
        this.connectPanel.hidden = this.role !== 'idle' || !!this.pending; this.controls.hidden = this.role !== 'instructor'; this.monitor.hidden = this.role !== 'student';
        this.status.hidden = !this.pending;
        this.status.textContent = this.pending === 'create' ? 'Starting classroom...' : 'Joining classroom...';
        this.heading.textContent = this.name || 'Classroom';
        if (this.role === 'idle' && !this.userActivated) this.hide(); else { this.dialog.hidden = false; this.button.hidden = true; }
    }
    private showError(error: Error): void {
        if (this.disposed) return;
        this.errorText.textContent = error.message || String(error); this.errorText.hidden = false;
        if (this.config.mode !== 'hidden' || this.role !== 'idle' || this.userActivated) {
            this.dialog.hidden = false; this.button.hidden = true;
        }
        this.adapter.showError?.(error);
    }
    private place(left: number, top: number): void {
        this.element.style.left = `${Math.max(0, Math.min(left, this.container.clientWidth - this.element.offsetWidth))}px`;
        this.element.style.top = `${Math.max(0, Math.min(top, this.container.clientHeight - this.element.offsetHeight))}px`; this.element.style.right = 'auto';
    }
    private keydown(event: KeyboardEvent): void {
        if (this.config.mode !== 'hidden' || event.repeat || (event.target as HTMLElement).closest('input,textarea,select,[contenteditable]')) return;
        const keys = (this.config.hotkey || 'Ctrl+Alt+I').toLowerCase().split('+').map(key => key.trim());
        const key = keys.pop();
        const matches = event.key.toLowerCase() === key || (/^[a-z]$/.test(key) && event.code === `Key${key.toUpperCase()}`);
        if (!matches || event.ctrlKey !== keys.includes('ctrl') || event.altKey !== keys.includes('alt') || event.shiftKey !== keys.includes('shift') || event.metaKey !== keys.includes('meta')) return;
        event.preventDefault(); event.stopPropagation();
        if (this.role === 'idle') this.start().catch(error => this.showError(error)); else this.open();
    }
    private snapshot(): any {
        const value = this.adapter.getState?.() || {
            ...this.viewer.getInternalState({ compact: true }), frame: this.viewer.getFrame(),
            viewSettings: { ...this.state.getViewSettings(), slab: this.viewer.getSlab() },
            trajectory: this.state.ui.getTrajectorySettings(), modelTitles: Array.from(this.state.modelTitles)
        };
        if (value.labels) value.labels = value.labels.filter(label => !label.style?._sessionTransient);
        return JSON.parse(JSON.stringify(value));
    }
    private schedule(): void {
        if (this.role !== 'instructor' || this.disposed) return;
        this.dirty = true; if (!this.timer) this.timer = setTimeout(() => this.flush(), 100);
    }
    private flush(): void {
        this.timer = undefined;
        if (this.role !== 'instructor' || !this.dirty || this.inFlight || this.socket.connected === false) return;
        let next; try { next = this.snapshot(); } catch (error) { this.showError(error); return; }
        const patch = {}; Object.keys(next).forEach(key => { if (!this.lastSent || JSON.stringify(next[key]) !== JSON.stringify(this.lastSent[key])) patch[key] = next[key]; });
        this.dirty = false;
        if (!this.full && this.lastSent && !Object.keys(patch).length) return;
        const message: any = { name: this.name, baseRevision: this.revision, revision: this.revision + 1 };
        if (this.full || !this.lastSent) message.state = next; else message.patch = patch;
        this.full = false; this.inFlight = next; this.socket.emit('viewer state change event', message);
    }
    private acknowledge(message): void {
        if (this.role !== 'instructor' || !Number.isInteger(message.revision) || message.revision < this.revision) return;
        this.revision = message.revision;
        if (this.inFlight) { this.lastSent = this.inFlight; this.inFlight = undefined; }
        if (this.dirty) this.schedule();
    }
    private pollScene(): void {
        if (this.role !== 'instructor' || this.socket.connected === false) return;
        const view = this.viewer.getView(), key = JSON.stringify(view);
        if (key !== this.lastView) { this.lastView = key; this.socket.emit('viewer view change event', { name: this.name, view }); }
        const settings = JSON.stringify([this.viewer.getFrame(), this.state.getViewSettings(), this.state.ui.getTrajectorySettings(), this.viewer.getSlab()]);
        if (settings !== this.lastSettings) { this.lastSettings = settings; this.schedule(); }
    }
    private receive(message) {
        if (this.role === 'instructor') { this.acknowledge(message); return; }
        if (this.role !== 'student') return;
        let sections: string[];
        if (message.patch) {
            if (message.revision <= this.revision) return;
            if (!this.received || message.baseRevision !== this.revision) { this.socket.emit('request state', { name: this.name }); return; }
            sections = Object.keys(message.patch); this.received = { ...this.received, ...message.patch };
        } else {
            if (Number.isInteger(message.revision) && message.revision < this.revision) return;
            this.received = message.state || message;
        }
        this.revision = message.revision || 0;
        return this.queueState(this.received, sections).catch(error => this.showError(error));
    }
    private queueState(snapshot: any, sections?: string[]): Promise<void> {
        this.pendingState = { snapshot, sections: this.pendingState ? (!sections || !this.pendingState.sections ? undefined : Array.from(new Set([...sections, ...this.pendingState.sections]))) : sections };
        if (!this.applying) {
            const generation = this.applyGeneration;
            this.applying = Promise.resolve().then(async () => {
                while (this.pendingState && generation === this.applyGeneration && !this.disposed) {
                    const current = this.pendingState; this.pendingState = undefined;
                    if (this.adapter.applyState) await this.adapter.applyState(current.snapshot, current.sections);
                    else await this.applyState(current.snapshot, current.sections, generation);
                    if (generation === this.applyGeneration && this.role === 'student') this.restoreQuery();
                }
            }).finally(() => {
                this.applying = undefined;
                if (this.pendingState && this.role === 'student' && !this.disposed) {
                    const pending = this.pendingState; this.pendingState = undefined;
                    this.queueState(pending.snapshot, pending.sections).catch(error => this.showError(error));
                }
            });
        }
        return this.applying;
    }
    private async applyState(snapshot: any, sections: string[], generation: number): Promise<void> {
        const full = !sections, changed = (key: string) => full || sections.includes(key), camera = this.viewer.getView();
        const scene = sections ? Object.fromEntries(sections.filter(key => !['frame', 'viewSettings', 'trajectory', 'modelTitles'].includes(key)).map(key => [key, snapshot[key]])) : { ...snapshot };
        if (changed('models')) {
            this.viewer.stopAnimate();
            if (snapshot.modelStyles) scene.modelStyles = snapshot.modelStyles;
            this.state.selections.clear(); this.state.surfaces.clear(); this.state.shapes.clear(); this.state.volumes.clear(); this.state.distances.clear();
        }
        this.viewer.setInternalState(scene, { partial: !full }); this.viewer.setView(camera);
        if (changed('models') || changed('frame')) await this.viewer.setFrame(Math.max(0, Number.isInteger(snapshot.frame) ? snapshot.frame : 0));
        if (generation !== this.applyGeneration || this.disposed || this.role !== 'student') return;
        if (changed('viewSettings') && snapshot.viewSettings) {
            this.state.setViewSettings(snapshot.viewSettings);
            if (snapshot.viewSettings.slab) this.viewer.setSlab(snapshot.viewSettings.slab.near, snapshot.viewSettings.slab.far);
        }
        if (changed('modelTitles') && snapshot.modelTitles) { this.state.modelTitles.clear(); snapshot.modelTitles.forEach(([id, title]) => this.state.modelTitles.set(id, title)); }
        if (changed('trajectory') && snapshot.trajectory) {
            const trajectory = { ...snapshot.trajectory, playing: false }; delete trajectory.frame;
            await this.state.ui.setTrajectorySettings(trajectory);
        }
        if (generation !== this.applyGeneration || this.disposed || this.role !== 'student') return;
        this.state.refresh(); this.viewer.render();
    }
    private reset(): void {
        const following = this.role === 'student', active = this.role !== 'idle';
        this.activation++; this.instructorRequested = false;
        this.role = 'idle'; this.name = null; this.secret = null; this.member = null; this.pending = null; this.pendingSent = false; this.userActivated = false;
        this.applyGeneration++; this.pendingState = undefined;
        clearTimeout(this.timer); clearInterval(this.poll); this.timer = this.poll = undefined;
        this.dirty = this.full = false; this.inFlight = this.lastSent = this.received = undefined; this.revision = 0;
        this.subscriptions.splice(0).forEach(remove => remove()); this.query = false; this.showingResults = false; this.clearResults(); this.clearChoices();
        if (following) this.restoreInteractions();
        if (active) {
            this.state.setControls(this.savedControls || { edit: true, trajectory: true }); this.state.ui.element.hidden = this.savedHidden;
        }
        if (following) this.adapter.roleChanged?.(null, false);
        this.savedControls = undefined; this.ask.textContent = 'Query Atoms'; this.responses.textContent = ''; this.refreshButton.hidden = true; this.renderRole();
    }
    private saveInteractions(): void {
        for (const atom of this.viewer.selectedAtoms({})) if (!this.interactions.has(atom)) this.interactions.set(atom, {
            clickable: atom.clickable, callback: atom.callback, hoverable: atom.hoverable, hover_callback: atom.hover_callback, unhover_callback: atom.unhover_callback
        });
    }
    private disableInteractions(): void { this.viewer.setClickable({}, false, undefined); this.viewer.setHoverable({}, false, undefined, undefined); this.viewer.render(); }
    private restoreInteractions(): void {
        const atoms = new Set(this.viewer.selectedAtoms({}));
        this.viewer.setClickable({}, false, undefined); this.viewer.setHoverable({}, false, undefined, undefined);
        for (const [atom, flags] of this.interactions) if (atoms.has(atom)) Object.assign(atom, flags);
        this.interactions.clear(); this.viewer.render();
    }
    private clearChoices(): void { if (this.choices.size) this.viewer.removeLabels(Array.from(this.choices.values())); this.choices.clear(); }
    private clearResults(): void { if (this.results.length) this.viewer.removeLabels(this.results); this.results = []; }
    private restoreQuery(): void {
        if (this.role !== 'student') return;
        this.monitor.textContent = this.query ? 'Select atom(s)' : 'Following instructor';
        if (!this.query) { this.disableInteractions(); return; }
        this.viewer.setHoverable({}, false, undefined, undefined);
        this.viewer.setClickable({}, true, atom => {
            if (this.role !== 'student' || !this.query) return;
            if (this.choices.has(atom)) { this.viewer.removeLabels([this.choices.get(atom)]); this.choices.delete(atom); }
            else {
                const text = [atom.chain, `${atom.resn || ''}${atom.resi ?? ''}`, atom.atom || atom.elem].filter(Boolean).join(':');
                this.choices.set(atom, this.viewer.addLabel(text, { _sessionTransient: true, position: atom, backgroundColor: '#800080', backgroundOpacity: 1 } as any));
            }
            this.socket.emit('query update', { name: this.name, selected: Array.from(this.choices.keys(), atom => [atom.x, atom.y, atom.z]) });
            this.viewer.render();
        }); this.viewer.render();
    }
    private queryAction(): void {
        if (this.role !== 'instructor') return;
        if (!this.query) this.socket.emit('query start', { name: this.name, kind: 'atoms' });
        else if (!this.showingResults) { this.showingResults = true; this.ask.textContent = 'Clear Labels'; this.refreshButton.hidden = false; this.socket.emit('query fetch', { name: this.name }); }
        else this.socket.emit('query end', { name: this.name });
    }
    private showResults(results): void {
        if (this.role !== 'instructor' || !Array.isArray(results)) return;
        this.clearResults(); const max = Math.max(1, ...results.map(item => Number(item[1])));
        this.results = results.map(([position, count]) => this.viewer.addLabel(String(count), {
            _sessionTransient: true, position: { x: position[0], y: position[1], z: position[2] }, backgroundOpacity: 0.5 + Number(count) / max * 0.5
        } as any)); this.viewer.render();
    }
}
