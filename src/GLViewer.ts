//a molecular viewer based on GLMol

import { Geometry, Renderer, Camera, Raycaster, Projector, Light, Fog, Scene, Coloring, FrontSide, Material } from "./WebGL";
import { Vector3, Matrix4, Matrix3, Quaternion, XYZ } from "./WebGL/math";
import { MeshLambertMaterial, Object3D, Mesh, LineBasicMaterial, Line } from "./WebGL";
import { elementColors, CC, ColorSpec, ColorschemeSpec } from "./colors";
import { extend, getExtent, makeFunction, getPropertyRange, isEmptyObject, adjustVolumeStyle, mergeGeos, PausableTimer, getColorFromStyle, getElement } from "./utilities";
import { getGradient, Gradient } from "./Gradient";
import { AtomStyleSpec, GLModel, LineStyleSpec } from "./GLModel";
import { Label, LabelSpec } from "./Label";
import { ArrowSpec, BoxSpec, CurveSpec, CustomShapeSpec, CylinderSpec, GLShape, IsoSurfaceSpec, LineSpec, ShapeSpec, SphereSpec, splitMesh } from "./GLShape";
import { VolumeData } from "./VolumeData";
import { ProteinSurface, SurfaceType, syncSurface } from "./ProteinSurface4";
import { GLVolumetricRender, VolumetricRendererSpec } from "./VolumetricRender";
import { AtomSelectionSpec, AtomSpec } from "./specs";
import { decode, toRGBA8, encode } from 'upng-js'


/**
 * WebGL-based 3Dmol.js viewer
 * Note: The preferred method of instantiating a GLViewer is through {@link createViewer}
 *
 * @class
*/
export class GLViewer {
    // private class variables
    private static numWorkers = 4; // number of threads for surface generation
    private static maxVolume = 64000; // how much to break up surface calculations

    private callback: any;
    private defaultcolors: any;
    private config: any;
    private nomouse = false;
    private bgColor: any;
    private camerax: number;
    private _viewer: GLViewer;
    private glDOM: HTMLCanvasElement | null = null;

    private models: GLModel[] = []; // atomistic molecular models
    private surfaces: any = {};
    private shapes = []; // Generic shapes
    private labels: Label[] = [];
    private clickables = []; //things you can click on
    private hoverables = []; //things you can hover over
    private contextMenuEnabledAtoms = []; // atoms with context menu
    private current_hover: any = null;
    private hoverDuration = 500;
    private viewer_frame = 0;
    private WIDTH: number;
    private HEIGHT: number;
    private viewChangeCallback: any = null;
    private stateChangeCallback: any = null;

    private NEAR = 1;
    private FAR = 800;
    private CAMERA_Z = 150;
    private fov = 20;

    private linkedViewers = [];
    private renderer: any = null;

    private row: number;
    private col: number;
    private cols: number;
    private rows: number;
    private viewers: any;
    private control_all = false;

    private ASPECT: any;
    private camera: Camera;
    private lookingAt: Vector3;

    private raycaster: Raycaster;
    private projector: Projector;

    private scene: any = null;
    private rotationGroup: any = null; // which contains modelGroup
    private modelGroup: any = null;

    private fogStart = 0.4;
    private slabNear = -50; // relative to the center of rotationGroup
    private slabFar = 50;

    public container: HTMLElement | null;

    static readonly surfaceTypeMap = {
        "VDW": SurfaceType.VDW,
        "MS": SurfaceType.MS,
        "SAS": SurfaceType.SAS,
        "SES": SurfaceType.SES
    };

    private cq = new Quaternion(0, 0, 0, 1);
    private dq = new Quaternion(0, 0, 0, 1);
    private animated = 0;
    private animationTimers = new Set<PausableTimer>();
    private isDragging = false;
    private mouseStartX = 0;
    private mouseStartY = 0;
    private touchDistanceStart = 0;
    private touchHold = false;
    private currentModelPos = 0;
    private cz = 0;
    private cslabNear = 0;
    private cslabFar = 0;

    private mouseButton: any;
    private hoverTimeout: any;

    private divwatcher: any;
    private spinInterval: any;


    //reimplement jquery getwidth/height
    private getRect() {
        let div = this.container;
        let rect = div.getBoundingClientRect();
        if (rect.width == 0 && rect.height == 0 && div.style.display === 'none') {
            let oldpos = div.style.position;
            let oldvis = div.style.visibility;
            div.style.display = 'block';
            div.style.visibility = 'hidden';
            div.style.position = 'absolute';
            rect = div.getBoundingClientRect();
            div.style.display = 'none';
            div.style.visibility = oldvis;
            div.style.position = oldpos;
        }
        return rect;
    };

    private getWidth() {
        return this.getRect().width;
    };

    private getHeight() {
        return this.getRect().height;
    };

    private setupRenderer() {

        this.renderer = new Renderer({
            antialias: this.config.antialias,
            preserveDrawingBuffer: true, //so we can export images
            premultipliedAlpha: false,/* more traditional compositing with background */
            id: this.config.id,
            row: this.config.row,
            col: this.config.col,
            rows: this.config.rows,
            cols: this.config.cols,
            canvas: this.config.canvas,
            //cannot initialize with zero size
            containerWidth: this.WIDTH || 1,
            containerHeight: this.HEIGHT || 1,
        });
        this.renderer.domElement.style.width = "100%";
        this.renderer.domElement.style.height = "100%";
        this.renderer.domElement.style.padding = "0";
        this.renderer.domElement.style.position = "absolute"; //TODO: get rid of this
        this.renderer.domElement.style.top = "0px";
        this.renderer.domElement.style.left = "0px";
        this.renderer.domElement.style.zIndex = "0";
    }

    private initializeScene() {

        this.scene = new Scene();
        this.scene.fog = new Fog(this.bgColor, 100, 200);

        this.modelGroup = new Object3D();
        this.rotationGroup = new Object3D();
        this.rotationGroup.useQuaternion = true;
        this.rotationGroup.quaternion = new Quaternion(0, 0, 0, 1);
        this.rotationGroup.add(this.modelGroup);

        this.scene.add(this.rotationGroup);

        // setup lights
        var directionalLight = new Light(0xFFFFFF);
        directionalLight.position = new Vector3(0.2, 0.2, 1)
            .normalize();
        directionalLight.intensity = 1.0;
        this.scene.add(directionalLight);
    };

    private initContainer(element) {
        this.container = element;
        this.WIDTH = this.getWidth();
        this.HEIGHT = this.getHeight();
        this.ASPECT = this.renderer.getAspect(this.WIDTH, this.HEIGHT);
        this.renderer.setSize(this.WIDTH, this.HEIGHT);
        this.container.append(this.renderer.domElement);
        this.glDOM = this.renderer.domElement;

        if (!this.nomouse) {
            // user can request that the mouse handlers not be installed
            this.glDOM.addEventListener('mousedown', this._handleMouseDown.bind(this), { passive: false });
            this.glDOM.addEventListener('touchstart', this._handleMouseDown.bind(this), { passive: false });
            this.glDOM.addEventListener('wheel', this._handleMouseScroll.bind(this), { passive: false });
            this.glDOM.addEventListener('mousemove', this._handleMouseMove.bind(this), { passive: false });
            this.glDOM.addEventListener('touchmove', this._handleMouseMove.bind(this), { passive: false });
            this.glDOM.addEventListener("contextmenu", this._handleContextMenu.bind(this), { passive: false });
        }

    };

    private decAnim() {
        //decrement the number of animations currently
        this.animated--;
        if (this.animated < 0) this.animated = 0;
    };

    private incAnim() {
        this.animated++;
    };

    private nextSurfID() {
        //compute the next highest surface id directly from surfaces
        //this is necessary to support linking of model data
        var max = 0;
        for (let i in this.surfaces) { // this is an object with possible holes
            if (!this.surfaces.hasOwnProperty(i)) continue;
            var val = parseInt(i);
            if (!isNaN(val)) {
                if (val > max)
                    max = val;
            }
        }
        return max + 1;
    };

    private setSlabAndFog() {

        let center = this.camera.position.z - this.rotationGroup.position.z;
        if (center < 1)
            center = 1;
        this.camera.near = center + this.slabNear;
        if (this.camera.near < 1)
            this.camera.near = 1;
        this.camera.far = center + this.slabFar;
        if (this.camera.near + 1 > this.camera.far)
            this.camera.far = this.camera.near + 1;

        this.camera.fov = this.fov;
        this.camera.right = center * Math.tan(Math.PI / 180 * this.fov);
        this.camera.left = -this.camera.right;
        this.camera.top = this.camera.right / this.ASPECT;
        this.camera.bottom = -this.camera.top;

        this.camera.updateProjectionMatrix();

        this.scene.fog.near = this.camera.near + this.fogStart * (this.camera.far - this.camera.near);
        // if (scene.fog.near > center) scene.fog.near = center;
        this.scene.fog.far = this.camera.far;

        if (this.config.disableFog) {
            this.scene.fog.near = this.scene.fog.far;
        }
    };

    // display scene
    //if nolink is set/true, don't propagate changes to linked viewers
    private show(nolink?) {
        this.renderer.setViewport();
        if (!this.scene)
            return;
        // var time = new Date();
        this.setSlabAndFog();
        this.renderer.render(this.scene, this.camera);
        // console.log("rendered in " + (+new Date() - time) + "ms");

        //have any scene change trigger a callback
        if (this.viewChangeCallback) this.viewChangeCallback(this._viewer.getView());

        if (!nolink && this.linkedViewers.length > 0) {
            var view = this._viewer.getView();
            for (var i = 0; i < this.linkedViewers.length; i++) {
                var other = this.linkedViewers[i];
                other.setView(view, true);
            }
        }
    };


    //regenerate the list of clickables
    //also updates hoverables
    private updateClickables() {
        this.clickables.splice(0, this.clickables.length);
        this.hoverables.splice(0, this.hoverables.length);
        this.contextMenuEnabledAtoms.splice(0, this.contextMenuEnabledAtoms.length);

        for (let i = 0, il = this.models.length; i < il; i++) {
            let model = this.models[i];
            if (model) {
                let atoms = model.selectedAtoms({
                    clickable: true
                });

                let hoverable_atoms = model.selectedAtoms({
                    hoverable: true
                });

                let contextMenuEnabled_atom = model.selectedAtoms({ contextMenuEnabled: true });
                // Array.prototype.push.apply(hoverables,hoverable_atoms);
                for (let n = 0; n < hoverable_atoms.length; n++) {
                    this.hoverables.push(hoverable_atoms[n]);
                }

                // Array.prototype.push.apply(clickables, atoms); //add atoms into clickables
                for (let m = 0; m < atoms.length; m++) {
                    this.clickables.push(atoms[m]);
                }

                // add atoms into contextMenuEnabledAtoms
                for (let m = 0; m < contextMenuEnabled_atom.length; m++) {
                    this.contextMenuEnabledAtoms.push(contextMenuEnabled_atom[m]);
                }

            }
        }
        for (let i = 0, il = this.shapes.length; i < il; i++) {

            let shape = this.shapes[i];
            if (shape && shape.clickable) {
                this.clickables.push(shape);
            }
            if (shape && shape.hoverable) {
                this.hoverables.push(shape);
            }
        }
    };

    // Checks for selection intersects on mousedown
    private handleClickSelection(mouseX: number, mouseY: number, event) {
        let intersects = this.targetedObjects(mouseX, mouseY, this.clickables);
        // console.log('handleClickSelection', mouseX, mouseY, intersects);
        if (intersects.length) {
            var selected = intersects[0].clickable;
            if (selected.callback !== undefined) {
                if (typeof (selected.callback) != "function") {
                    selected.callback = makeFunction(selected.callback);
                }
                if (typeof (selected.callback) === "function") {
                    selected.callback(selected, this._viewer, event, this.container, intersects);
                }
            }
        }
    };


    //return offset of container
    private canvasOffset() {
        let canvas = this.glDOM;
        let rect = canvas.getBoundingClientRect();
        let doc = canvas.ownerDocument;
        let docElem = doc.documentElement;
        let win = doc.defaultView;
        return {
            top: rect.top + win.pageYOffset - docElem.clientTop,
            left: rect.left + win.pageXOffset - docElem.clientLeft
        };
    };

    //set current_hover to sel (which can be null), calling appropraite callbacks
    private setHover(selected, event?, intersects?) {
        if (this.current_hover == selected) return;
        if (this.current_hover) {
            if (typeof (this.current_hover.unhover_callback) != "function") {
                this.current_hover.unhover_callback = makeFunction(this.current_hover.unhover_callback);
            }
            this.current_hover.unhover_callback(this.current_hover, this._viewer, event, this.container, intersects);
        }
        this.current_hover = selected;

        if (selected && selected.hover_callback !== undefined) {
            if (typeof (selected.hover_callback) != "function") {
                selected.hover_callback = makeFunction(selected.hover_callback);
            }
            if (typeof (selected.hover_callback) === "function") {
                selected.hover_callback(selected, this._viewer, event, this.container, intersects);
            }
        }

    };

    //checks for selection intersects on hover
    private handleHoverSelection(mouseX, mouseY, event) {
        if (this.hoverables.length == 0) return;
        let intersects = this.targetedObjects(mouseX, mouseY, this.hoverables);
        if (intersects.length) {
            var selected = intersects[0].clickable;
            this.setHover(selected, event, intersects);
            this.current_hover = selected;
        }
        else {
            this.setHover(null);
        }
    };

    //sees if the mouse is still on the object that invoked a hover event and if not then the unhover callback is called
    private handleHoverContinue(mouseX: number, mouseY: number) {
        let intersects = this.targetedObjects(mouseX, mouseY, this.hoverables);
        if (intersects.length == 0 || intersects[0] === undefined) {
            this.setHover(null);
        }
        if (intersects[0] !== undefined && intersects[0].clickable !== this.current_hover) {
            this.setHover(null);
        }
    };

    /**
     * Determine if a positioned event is "close enough" to mouseStart to be considered a click.
     * With a mouse, the position should be exact, but allow a slight delta for a touch interface.
     * @param {Event} event
     * @param {{ allowTolerance, tolerance: number }} options
     */
    private closeEnoughForClick(event, { allowTolerance=event.targetTouches, tolerance=5}={}) {
        const x = this.getX(event);
        const y = this.getY(event);
        if (allowTolerance) {
            const deltaX = Math.abs(x - this.mouseStartX);
            const deltaY = Math.abs(y - this.mouseStartY);
            return deltaX <= tolerance && deltaY <= tolerance;
        } else {
            return x === this.mouseStartX && y === this.mouseStartY;
        }
    }

    private calcTouchDistance(ev) { // distance between first two
        // fingers
        var xdiff = ev.targetTouches[0].pageX -
            ev.targetTouches[1].pageX;
        var ydiff = ev.targetTouches[0].pageY -
            ev.targetTouches[1].pageY;
        return Math.sqrt(xdiff * xdiff + ydiff * ydiff);
    };

    //check targetTouches as well
    private getX(ev) {
        var x = ev.pageX;
        if (x == undefined) x = ev.pageX; //firefox
        if (ev.targetTouches &&
            ev.targetTouches[0]) {
            x = ev.targetTouches[0].pageX;
        }
        else if (ev.changedTouches &&
            ev.changedTouches[0]) {
            x = ev.changedTouches[0].pageX;
        }
        return x;
    };

    private getY(ev) {
        var y = ev.pageY;
        if (y == undefined) y = ev.pageY;
        if (ev.targetTouches &&
            ev.targetTouches[0]) {
            y = ev.targetTouches[0].pageY;
        }
        else if (ev.changedTouches &&
            ev.changedTouches[0]) {
            y = ev.changedTouches[0].pageY;
        }
        return y;
    };

    //for grid viewers, return true if point is in this viewer
    private isInViewer(x: number, y: number) {
        if (this.viewers != undefined && !this.control_all) {
            var width = this.WIDTH / this.cols;
            var height = this.HEIGHT / this.rows;
            var offset = this.canvasOffset();
            var relx = (x - offset.left);
            var rely = (y - offset.top);

            var r = this.rows - Math.floor(rely / height) - 1;
            var c = Math.floor(relx / width);

            if (r != this.row || c != this.col)
                return false;
        }
        return true;
    };

    //if the user has specify zoom limits, readjust to fit within them
    //also, make sure we don't go past CAMERA_Z
    private adjustZoomToLimits(z: number) {
        //a lower limit of 0 is at CAMERA_Z
        if (this.config.lowerZoomLimit && this.config.lowerZoomLimit > 0) {
            let lower = this.CAMERA_Z - this.config.lowerZoomLimit;
            if (z > lower) z = lower;
        }

        if (this.config.upperZoomLimit && this.config.upperZoomLimit > 0) {
            let upper = this.CAMERA_Z - this.config.upperZoomLimit;
            if (z < upper) z = upper;
        }

        if (z > this.CAMERA_Z) {
            z = this.CAMERA_Z * 0.999; //avoid getting stuck
        }
        return z;
    };
    //interpolate between two normalized quaternions (t between 0 and 1)
    //https://en.wikipedia.org/wiki/Slerp
    private static slerp(v0: Quaternion, v1: Quaternion, t: number) {
        // Compute the cosine of the angle between the two vectors.
        //dot product
        if (t == 1) return v1.clone();
        else if (t == 0) return v0.clone();
        let dot = v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w * v1.w;
        if (dot > 0.9995) {
            // If the inputs are too close for comfort, linearly interpolate
            // and normalize the result.
            let result = new Quaternion(
                v0.x + t * (v1.x - v0.x),
                v0.y + t * (v1.y - v0.y),
                v0.z + t * (v1.z - v0.z),
                v0.w + t * (v1.w - v0.w));

            result.normalize();
            return result;
        }

        // If the dot product is negative, the quaternions
        // have opposite handed-ness and slerp won't take
        // the shorted path. Fix by reversing one quaternion.
        if (dot < 0.0) {
            v1 = v1.clone().multiplyScalar(-1);
            dot = -dot;
        }

        if (dot > 1) dot = 1.0;
        else if (dot < -1) dot = -1.0;

        var theta_0 = Math.acos(dot);  // theta_0 = angle between input vectors
        var theta = theta_0 * t;    // theta = angle between v0 and result

        var v2 = v1.clone();
        v2.sub(v0.clone().multiplyScalar(dot));
        v2.normalize();              // { v0, v2 } is now an orthonormal basis

        var c = Math.cos(theta);
        var s = Math.sin(theta);
        var ret = new Quaternion(
            v0.x * c + v2.x * s,
            v0.y * c + v2.y * s,
            v0.z * c + v2.z * s,
            v0.w * c + v2.w * s
        );
        ret.normalize();
        return ret;
    };


    /* @param {Object} element HTML element within which to create viewer
     * @param {ViewerSpec} config Object containing optional configuration for the viewer
     */
    constructor(element, c: ViewerSpec = {}) {
        // set variables
        this.config = c;
        this.callback = this.config.callback;
        this.defaultcolors = this.config.defaultcolors;
        if (!this.defaultcolors)
            this.defaultcolors = elementColors.defaultColors;
        this.nomouse = this.config.nomouse;
        this.bgColor = 0;
        this.config.backgroundColor = this.config.backgroundColor || "#ffffff";
        if (typeof (this.config.backgroundColor) != 'undefined') {
            this.bgColor = CC.color(this.config.backgroundColor).getHex();
        }
        this.config.backgroundAlpha = this.config.backgroundAlpha == undefined ? 1.0 : this.config.backgroundAlpha;

        this.camerax = 0;
        if (typeof (this.config.camerax) != 'undefined') {
            this.camerax = parseFloat(this.config.camerax);
        }
        this._viewer = this;
        this.container = element; //we expect container to be HTMLElement

        if (this.config.hoverDuration != undefined) {
            this.hoverDuration = this.config.hoverDuration;
        }
        if (this.config.antialias === undefined) this.config.antialias = true;
        if (this.config.cartoonQuality === undefined) this.config.cartoonQuality = 10;

        this.WIDTH = this.getWidth();
        this.HEIGHT = this.getHeight();

        this.setupRenderer();

        this.row = this.config.row;
        this.col = this.config.col;
        this.cols = this.config.cols;
        this.rows = this.config.rows;
        this.viewers = this.config.viewers;
        this.control_all = this.config.control_all;

        this.ASPECT = this.renderer.getAspect(this.WIDTH, this.HEIGHT);

        this.camera = new Camera(this.fov, this.ASPECT, this.NEAR, this.FAR, this.config.orthographic);
        this.camera.position = new Vector3(this.camerax, 0, this.CAMERA_Z);
        this.lookingAt = new Vector3();
        this.camera.lookAt(this.lookingAt);

        this.raycaster = new Raycaster(new Vector3(0, 0, 0), new Vector3(0, 0, 0));
        this.projector = new Projector();

        this.initializeScene();
        this.renderer.setClearColorHex(this.bgColor, this.config.backgroundAlpha);
        this.scene.fog.color = CC.color(this.bgColor);

        // this event is bound to the body element, not the container,
        // so no need to put it inside initContainer()
        document.body.addEventListener('mouseup', this._handleMouseUp.bind(this));
        document.body.addEventListener('touchend', this._handleMouseUp.bind(this));

        this.initContainer(this.container);
        if (this.config.style) { //enable setting style in constructor
            this.setViewStyle(this.config);
        }

        window.addEventListener("resize", this.resize.bind(this));

        if (typeof (window.ResizeObserver) !== "undefined") {
            this.divwatcher = new window.ResizeObserver(this.resize.bind(this));
            this.divwatcher.observe(this.container);
        }

        try {
            if (typeof (this.callback) === "function")
                this.callback(this);
        } catch (e) {
            // errors in callback shouldn't invalidate the viewer
            console.log("error with glviewer callback: " + e);
        }

    };

    /**
    * Return a list of objects that intersect that at the specified viewer position.
    *
    * @param x - x position in screen coordinates
    * @param y - y position in screen coordinates
    * @param {Object[]} - list of objects or selection object specifying what object to check for targeting
    */
    public targetedObjects(x: number, y: number, objects) {
        var mouse = {
            x: x,
            y: y,
            z: -1.0
        };
        if (!Array.isArray(objects)) { //assume selection object
            objects = this.selectedAtoms(objects);
        }
        if (objects.length == 0) return [];
        this.raycaster.setFromCamera(mouse, this.camera);
        return this.raycaster.intersectObjects(this.modelGroup, objects);
    };

    /** Convert model coordinates to screen coordinates.
     * @param {object | list} - an object or list of objects with x,y,z attributes (e.g. an atom)
     * @return {object | list} - and object or list of {x: screenX, y: screenY}
     */
    public modelToScreen(coords) {
        let returnsingle = false;
        if (!Array.isArray(coords)) {
            coords = [coords];
            returnsingle = true;
        }

        let results = [];
        let offset = this.canvasOffset();
        coords.forEach(coord => {
            let t = new Vector3(coord.x, coord.y, coord.z);
            t.applyMatrix4(this.modelGroup.matrixWorld);
            this.projector.projectVector(t, this.camera);
            let screenX = this.WIDTH * (t.x + 1) / 2.0 + offset.left;
            let screenY = -this.HEIGHT * (t.y - 1) / 2.0 + offset.top;
            results.push({ x: screenX, y: screenY });
        });
        if (returnsingle) results = results[0];
        return results;
    };

    /**
     * For a given screen (x,y) displacement return model displacement
     * @param{x} x displacement in screen coordinates
     * @param{y} y displacement in screen corodinates
     * @param{modelz} z coordinate in model coordinates to compute offset for, default is model axis
    */
    public screenOffsetToModel(x: number, y: number, modelz?) {
        var dx = x / this.WIDTH;
        var dy = y / this.HEIGHT;
        var zpos = (modelz === undefined ? this.rotationGroup.position.z : modelz);
        var q = this.rotationGroup.quaternion;
        var t = new Vector3(0, 0, zpos);
        this.projector.projectVector(t, this.camera);
        t.x += dx * 2;
        t.y -= dy * 2;
        this.projector.unprojectVector(t, this.camera);
        t.z = 0;
        t.applyQuaternion(q);
        return t;
    };

    /**
     * Distance from screen coordinate to model coordinate assuming screen point
     * is projected to the same depth as model coordinate
     * @param{screen} xy screen coordinate
     * @param{model} xyz model coordinate
    */
    public screenToModelDistance(screen: XYZ, model) {
        let offset = this.canvasOffset();

        //convert model to screen to get screen z
        let mvec = new Vector3(model.x, model.y, model.z);
        mvec.applyMatrix4(this.modelGroup.matrixWorld);
        let m = mvec.clone();
        this.projector.projectVector(mvec, this.camera);

        let t = new Vector3((screen.x - offset.left) * 2 / this.WIDTH - 1, (screen.y - offset.top) * 2 / -this.HEIGHT + 1, mvec.z);
        this.projector.unprojectVector(t, this.camera);

        return t.distanceTo(m);
    };

    /**
     * Set a callback to call when the view has potentially changed.
     *
    */
    public setViewChangeCallback(callback) {
        if (typeof (callback) === 'function' || callback == null)
            this.viewChangeCallback = callback;
    };

    /**
     * Set a callback to call when the view has potentially changed.
     *
    */
    public setStateChangeCallback(callback) {
        if (typeof (callback) === 'function' || callback == null)
            this.stateChangeCallback = callback;
    };

    /**
     * Return configuration of viewer
     */
    public getConfig() {
        return this.config;
    };

    /**
     * Set the configuration object.  Note that some setting may only
     * have an effect at viewer creation time.
     */
    public setConfig(c) {
        this.config = c;
    };

    /**
     * Return object representing internal state of
     * the viewer appropriate for passing to setInternalState
     *
    */
    public getInternalState() {
        var ret = { 'models': [], 'surfaces': [], 'shapes': [], 'labels': [] };
        for (let i = 0; i < this.models.length; i++) {
            if (this.models[i]) {
                ret.models[i] = this.models[i].getInternalState();
            }
        }

        //todo: labels, shapes, surfaces

        return ret;
    };

    /**
     * Overwrite internal state of the viewer with passed  object
     * which should come from getInternalState.
     *
    */
    public setInternalState(state) {

        //clear out current viewer
        this.clear();

        //set model state
        var newm = state.models;
        for (let i = 0; i < newm.length; i++) {
            if (newm[i]) {
                this.models[i] = new GLModel(i);
                this.models[i].setInternalState(newm[i]);
            }
        }

        //todo: labels, shapes, surfaces
        this.render();
    };

    /**
     * Set lower and upper limit stops for zoom.
     *
     * @param {lower} - limit on zoom in (positive number).  Default 0.
     * @param {upper} - limit on zoom out (positive number).  Default infinite.
     * @example
      $3Dmol.get("data/set1_122_complex.mol2", function(moldata) {
            var m = viewer.addModel(moldata);
            viewer.setStyle({stick:{colorscheme:"Jmol"}});
            viewer.setZoomLimits(100,200);
            viewer.zoomTo();
            viewer.zoom(10); //will not zoom all the way
            viewer.render();
        });
    */
    public setZoomLimits(lower, upper) {
        if (typeof (lower) !== 'undefined') this.config.lowerZoomLimit = lower;
        if (upper) this.config.upperZoomLimit = upper;
        this.rotationGroup.position.z = this.adjustZoomToLimits(this.rotationGroup.position.z);
        this.show();
    };

    /**
     * Set camera parameters (distance to the origin and field of view)
     *
     * @param {parameters} - new camera parameters, with possible fields
     *                       being fov for the field of view, z for the
     *                       distance to the origin, and orthographic (boolean)
     *                       for kind of projection (default false).
     * @example
      $3Dmol.get("data/set1_122_complex.mol2", function(data) {
            var m = viewer.addModel(data);
            viewer.setStyle({stick:{}});
            viewer.zoomTo();
            viewer.setCameraParameters({ fov: 10 , z: 300 });
            viewer.render();
        });
    */
    public setCameraParameters(parameters) {
        if (parameters.fov !== undefined) {
            this.fov = parameters.fov;
            this.camera.fov = this.fov;
        }

        if (parameters.z !== undefined) {
            this.CAMERA_Z = parameters.z;
            this.camera.z = this.CAMERA_Z;
        }
        if (parameters.orthographic !== undefined) {
            this.camera.ortho = parameters.orthographic;
        }
    };

    public _handleMouseDown(ev) {
        ev.preventDefault();
        if (!this.scene)
            return;
        var x = this.getX(ev);
        var y = this.getY(ev);
        if (x === undefined)
            return;
        this.isDragging = true;
        this.mouseButton = ev.which;
        this.mouseStartX = x;
        this.mouseStartY = y;
        this.touchHold = true;
        this.touchDistanceStart = 0;
        if (ev.targetTouches &&
            ev.targetTouches.length == 2) {
            this.touchDistanceStart = this.calcTouchDistance(ev);
        }
        this.cq = this.rotationGroup.quaternion.clone();
        this.cz = this.rotationGroup.position.z;
        this.currentModelPos = this.modelGroup.position.clone();
        this.cslabNear = this.slabNear;
        this.cslabFar = this.slabFar;

        let self = this;
        setTimeout(function () {
            if (ev.targetTouches) {
                if (self.touchHold == true) {
                    // console.log('Touch hold', x,y);
                    self.glDOM = self.renderer.domElement;
                    self.glDOM.dispatchEvent(new Event('contextmenu'));
                }
                else {
                    // console.log('Touch hold ended earlier');

                }
            }
        }, 1000);

    };

    public _handleMouseUp(ev) {
        // handle touch
        this.touchHold = false;

        // handle selection
        if (this.isDragging && this.scene) { //saw mousedown, haven't moved
            var x = this.getX(ev);
            var y = this.getY(ev);
            if (this.closeEnoughForClick(ev)) {
                var offset = this.canvasOffset();
                var mouseX = ((x - offset.left) / this.WIDTH) * 2 - 1;
                var mouseY = -((y - offset.top) / this.HEIGHT) * 2 + 1;
                this.handleClickSelection(mouseX, mouseY, ev);
            }
        }

        this.isDragging = false;
    }

    public _handleMouseScroll(ev) { // Zoom
        ev.preventDefault();
        if (!this.scene)
            return;

        var x = this.getX(ev);
        var y = this.getY(ev);
        if (x === undefined)
            return;
        if (!this.isInViewer(x, y)) {
            return;
        }

        var scaleFactor = (this.CAMERA_Z - this.rotationGroup.position.z) * 0.85;
        var mult = 1.0;
        if (ev.ctrlKey) {
            mult = -1.0; //this is a pinch event turned into a wheel event (or they're just holding down the ctrl)
        }
        if (ev.detail) {
            this.rotationGroup.position.z += mult * scaleFactor * ev.detail / 10;
        } else if (ev.wheelDelta) {
            this.rotationGroup.position.z -= mult * scaleFactor * ev.wheelDelta / 400;
        }
        this.rotationGroup.position.z = this.adjustZoomToLimits(this.rotationGroup.position.z);
        this.show();
    };

    /**
     * Return image URI of viewer contents (base64 encoded).     *
     */
    public pngURI() {
        return this.getCanvas().toDataURL('image/png');
    };

    /**
     * Return a promise that resolves to an animated PNG image URI of
     viewer contents (base64 encoded) for nframes of viewer changes.
     * @return {Promise}
     */
    public apngURI(nframes: number) {
        let viewer = this;
        nframes = nframes ? nframes : 1;
        return new Promise(function (resolve) {
            let framecnt = 0;
            let oldcb = viewer.viewChangeCallback;
            let bufpromise = [];
            let delays = [];
            let lasttime = Date.now();
            viewer.viewChangeCallback = function () {
                delays.push(Date.now() - lasttime);
                lasttime = Date.now();
                bufpromise.push(new Promise(resolve => {
                    viewer.getCanvas().toBlob(function (blob) {
                        blob.arrayBuffer().then(resolve);
                    }, "image/png");
                }));
                framecnt += 1;
                if (framecnt == nframes) {
                    viewer.viewChangeCallback = oldcb;

                    Promise.all(bufpromise).then((buffers) => {
                        //convert to apng
                        let rgbas = [];
                        //have to convert png to rgba, before creating the apng
                        for (let i = 0; i < buffers.length; i++) {
                            let img = decode(buffers[i]);
                            rgbas.push(toRGBA8(img)[0]);
                        }
                        let width = viewer.getCanvas().width;
                        let height = viewer.getCanvas().height;
                        let apng = encode(rgbas, width, height, 0, delays);
                        let blob = new Blob([apng], { type: 'image/png' });
                        let fr = new FileReader();
                        fr.onload = function (e) {
                            resolve(e.target.result);
                        };
                        fr.readAsDataURL(blob);
                    });
                }
            };
        });

    };


    /**
     * Return underlying canvas element.
     */
    public getCanvas(): HTMLCanvasElement {
        return this.glDOM;
    };

    /**
     * Return renderer element.
     */
    public getRenderer() {
        return this.renderer;
    };

    /**
         * Set the duration of the hover delay
         *
         * @param {number}
         *            [hoverDuration] - an optional parameter that denotes
         *            the duration of the hover delay (in milliseconds) before the hover action is called
         *
     */
    public setHoverDuration(duration?: number) {
        this.hoverDuration = duration;
    };

    public _handleMouseMove(ev) { // touchmove

        clearTimeout(this.hoverTimeout);
        var offset = this.canvasOffset();
        var mouseX = ((this.getX(ev) - offset.left) / this.WIDTH) * 2 - 1;
        var mouseY = -((this.getY(ev) - offset.top) / this.HEIGHT) * 2 + 1;
        let self = this;
        // hover timeout
        if (this.current_hover !== null) {
            this.handleHoverContinue(mouseX, mouseY);
        }

        if (this.hoverables.length > 0) {
            this.hoverTimeout = setTimeout(
                function () {
                    self.handleHoverSelection(mouseX, mouseY, ev);
                },
                this.hoverDuration);
        }

        ev.preventDefault();
        if (!this.scene)
            return;
        if (!this.isDragging)
            return;
        var mode = 0;

        var x = this.getX(ev);
        var y = this.getY(ev);
        if (x === undefined)
            return;

        if (!this.isInViewer(x, y)) {
            return;
        }


        var dx = (x - this.mouseStartX) / this.WIDTH;
        var dy = (y - this.mouseStartY) / this.HEIGHT;
        // check for pinch
        if (this.touchDistanceStart != 0 &&
            ev.targetTouches &&
            ev.targetTouches.length == 2) {
            var newdist = this.calcTouchDistance(ev);
            // change to zoom
            mode = 2;
            dy = (newdist - this.touchDistanceStart) * 2 / (this.WIDTH + this.HEIGHT);
        } else if (ev.targetTouches &&
            ev.targetTouches.length == 3) {
            // translate
            mode = 1;
        }
        var ratioX = this.renderer.getXRatio();
        var ratioY = this.renderer.getYRatio();
        dx *= ratioX;
        dy *= ratioY;
        var r = Math.sqrt(dx * dx + dy * dy);
        var scaleFactor;
        if (mode == 3 || (this.mouseButton == 3 && ev.ctrlKey)) { // Slab
            this.slabNear = this.cslabNear + dx * 100;
            this.slabFar = this.cslabFar - dy * 100;
        } else if (mode == 2 || this.mouseButton == 3 || ev.shiftKey) { // Zoom
            scaleFactor = (this.CAMERA_Z - this.rotationGroup.position.z) * 0.85;
            if (scaleFactor < 80)
                scaleFactor = 80;
            this.rotationGroup.position.z = this.cz + dy * scaleFactor;
            this.rotationGroup.position.z = this.adjustZoomToLimits(this.rotationGroup.position.z);
        } else if (mode == 1 || this.mouseButton == 2 || ev.ctrlKey) { // Translate
            var t = this.screenOffsetToModel(ratioX * (x - this.mouseStartX), ratioY * (y - this.mouseStartY));
            this.modelGroup.position.addVectors(this.currentModelPos, t);

        } else if ((mode === 0 || this.mouseButton == 1) && r !== 0) { // Rotate
            var rs = Math.sin(r * Math.PI) / r;
            this.dq.x = Math.cos(r * Math.PI);
            this.dq.y = 0;
            this.dq.z = rs * dx;
            this.dq.w = -rs * dy;
            this.rotationGroup.quaternion.set(1, 0, 0, 0);
            this.rotationGroup.quaternion.multiply(this.dq);
            this.rotationGroup.quaternion.multiply(this.cq);
        }
        this.show();
    };

    /** User specified function for handling a context menu event.
     * Handler is passed the selected object, x and y in canvas coordinates,
     * and original event.
     */
    public userContextMenuHandler: Function | null = null;

    public _handleContextMenu(ev) {
        ev.preventDefault();
        var newX = this.getX(ev);
        var newY = this.getY(ev);

        if (newX != this.mouseStartX || newY != this.mouseStartY) {
            return;
        } else {
            var x = this.mouseStartX;
            var y = this.mouseStartY;
            var offset = this.canvasOffset();
            var mouseX = ((x - offset.left) / this.WIDTH) * 2 - 1;
            var mouseY = -((y - offset.top) / this.HEIGHT) * 2 + 1;
            let intersects = this.targetedObjects(mouseX, mouseY, this.contextMenuEnabledAtoms);
            var selected = null;
            if (intersects.length) {
                selected = intersects[0].clickable;
            }

            var offset = this.canvasOffset();
            var x = this.mouseStartX - offset.left;
            var y = this.mouseStartY - offset.top;
            if (this.userContextMenuHandler) {
                this.userContextMenuHandler(selected, x, y,intersects);
            }
        }
    };


    /**
     * Change the viewer's container element
     * Also useful if the original container element was removed from the DOM.
     *
     * @param {Object | string} element
     *            Either HTML element or string identifier. Defaults to the element used to initialize the viewer.

     */
    public setContainer(element) {
        let elem = getElement(element) || this.container;
        this.initContainer(elem);
        return this;
    };

    /**
     * Set the background color (default white)
     *
     * @param {number}
     *            hex Hexcode specified background color, or standard color spec
     * @param {number}
     *            a Alpha level (default 1.0)
     *
     * @example
     *
     * viewer.setBackgroundColor(0x000000);


     *
     */
    public setBackgroundColor(hex: ColorSpec, a: number) {
        if (typeof (a) == "undefined") {
            a = 1.0;
        }
        else if (a < 0 || a > 1.0) {
            a = 1.0;
        }
        var c = CC.color(hex);
        this.scene.fog.color = c;
        this.bgColor = c.getHex();
        this.renderer.setClearColorHex(c.getHex(), a);
        this.show();

        return this;
    };

    /**
     * Set view projection scheme.  Either orthographic or perspective.
     * Default is perspective.  Orthographic can also be enabled on viewer creation
     * by setting orthographic to true in the config object.
     *
     *
     * @example
     viewer.setViewStyle({style:"outline"});
          $3Dmol.get('data/1fas.pqr', function(data){
              viewer.addModel(data, "pqr");
              $3Dmol.get("data/1fas.cube",function(volumedata){
                  viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
              });
              viewer.zoomTo();

              viewer.setProjection("orthographic");
              viewer.render(callback);
          });
     *
     */
    public setProjection(proj) {
        this.camera.ortho = (proj === "orthographic");
        this.setSlabAndFog();
    };

    /**
     * Set global view styles.
     *
     * @example
     *   viewer.setViewStyle({style:"outline"});
          $3Dmol.get('data/1fas.pqr', function(data){
              viewer.addModel(data, "pqr");
              $3Dmol.get("data/1fas.cube",function(volumedata){
                  viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
              });
              viewer.zoomTo();
              viewer.render(callback);
          });
     *
     */
    public setViewStyle(parameters) {
        if (parameters.style === "outline") {
            var params: any = {};
            if (parameters.color) params.color = CC.color(parameters.color);
            if (parameters.width) params.width = parameters.width;
            this.renderer.enableOutline(params);
        } else {
            this.renderer.disableOutline();
        }
        return this;
    };


    private updateSize() {
        this.renderer.setSize(this.WIDTH, this.HEIGHT);
        this.ASPECT = this.renderer.getAspect(this.WIDTH, this.HEIGHT);
        this.renderer.setSize(this.WIDTH, this.HEIGHT);
        this.camera.aspect = this.ASPECT;
        this.camera.updateProjectionMatrix();
    }
    /**
     * Set viewer width independently of the HTML container.  This is probably not what you want.
     *
     * @param {number} w Width in pixels
     */
    public setWidth(w: number) {
        this.WIDTH = w || this.WIDTH;
        this.updateSize();
        return this;
    };

    /**
     * Set viewer height independently of the HTML container.  This is probably not what you want.
     *
     * @param {number} h Height in pixels
     */
    public setHeight(h: number) {
        this.HEIGHT = h || this.HEIGHT;
        this.updateSize();
        return this;
    };

    /**
     * Resize viewer according to containing HTML element's dimensions
     *
     */
    public resize() {
        this.WIDTH = this.getWidth();
        this.HEIGHT = this.getHeight();
        let regen = false;
        if (this.renderer.isLost() && this.WIDTH > 0 && this.HEIGHT > 0) {
            //create new context
            this.container.querySelector('canvas').remove(); //remove existing
            this.setupRenderer();
            this.initContainer(this.container);
            regen = true;
        }
        if (this.WIDTH == 0 || this.HEIGHT == 0) {
            if (this.animated) this._viewer.pauseAnimate();
        } else if (this.animated) {
            this._viewer.resumeAnimate();
        }

        this.updateSize();

        if (regen) { //restored rendere, need to regenerate scene
            let options = this.renderer.supportedExtensions();
            options.regen = true;
            this._viewer.render(null, options);
        } else {
            this.show();
        }

        return this;
    };


    /**
     * Return specified model
     *
     * @param {number}
     *            [id=last model id] - Retrieve model with specified id
     * @default Returns last model added to viewer or null if there are no models
     * @return {GLModel}
     *
     * @example // Retrieve reference to first GLModel added var m =
     *    $3Dmol.download("pdb:1UBQ",viewer,{},function(m1){
              $3Dmol.download("pdb:1UBI", viewer,{}, function(m2) {
                viewer.zoomTo();
                m1.setStyle({cartoon: {color:'green'}});
                //could use m2 here as well
                viewer.getModel().setStyle({cartoon: {color:'blue'}});
                viewer.render();
            })
          });
     */
    public getModel(id?: number | GLModel) {
        if (id === undefined) {
            return this.models.length == 0 ? null : this.models[this.models.length - 1];
        }
        if (id instanceof GLModel) {
            return id;
        }
        if (!(id in this.models)) {
            if (this.models.length == 0)
                return null;
            else
                return this.models[this.models.length - 1]; //get last model if no (or invalid) id specified
        }
        return this.models[id];
    };


    /**
     * Continuously rotate a scene around the specified axis.
     *
     * Call `spin(false)` to stop spinning.
     *
     * @param  {string|boolean|Array} axis
     *            [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around.
     *            Default "y".  View relative (rather than model relative) axes are prefixed with v.
     * @param  {number} speed
     *            [speed] - Speed multiplier for spinning the viewer. 1 is default and a negative
     *             value reverses the direction of the spin.
     *
     */
    public spin(axis, speed: number = 1) {
        clearInterval(this.spinInterval);
        if (typeof axis == 'undefined')
            axis = 'y';
        if (typeof axis == "boolean") {
            if (!axis)
                return;
            else
                axis = 'y';
        }

        if (Array.isArray(axis)) {
            axis = { x: axis[0], y: axis[1], z: axis[2] };
        }
        //out of bounds check

        var viewer = this;

        this.spinInterval = setInterval(
            function () {
                if (!viewer.getCanvas().isConnected && viewer.renderer.isLost()) {
                    clearInterval(viewer.spinInterval);
                }
                viewer.rotate(1 * speed, axis);
            }, 25);

    };

    //animate motion between current position and passed position
    // can set some parameters to null
    //if fixed is true will enforce the request animation, otherwise
    //does relative updates
    //positions objects have modelggroup position, rotation group position.z,
    //and rotationgroup quaternion
    //return array includes final position, but not current
    //the returned array includes an animate method
    private animateMotion(duration: number, fixed: boolean, mpos: Vector3, rz: number, rot: Quaternion, cam: Vector3) {
        var interval = 20;
        var nsteps: number = Math.ceil(duration / interval);
        if (nsteps < 1) nsteps = 1;
        this.incAnim();

        var curr = {
            mpos: this.modelGroup.position.clone(),
            rz: this.rotationGroup.position.z,
            rot: this.rotationGroup.quaternion.clone(),
            cam: this.lookingAt.clone()
        };

        if (fixed) { //precompute path and stick to it
            let steps = new Array(nsteps);
            for (let i = 0; i < nsteps; i++) {
                let frac = (i + 1) / nsteps;
                let next: any = { mpos: curr.mpos, rz: curr.rz, rot: curr.rot };
                next.mpos = mpos.clone().sub(curr.mpos).multiplyScalar(frac).add(curr.mpos);
                next.rz = curr.rz + frac * (rz - curr.rz);
                next.rot = GLViewer.slerp(curr.rot, rot, frac);
                next.cam = cam.clone().sub(curr.cam).multiplyScalar(frac).add(curr.cam);
                steps[i] = next;
            }

            let step = 0;
            let self = this;
            let callback = function () {
                var p = steps[step];
                step += 1;
                self.modelGroup.position = p.mpos;
                self.rotationGroup.position.z = p.rz;
                self.rotationGroup.quaternion = p.rot;
                self.camera.lookAt(p.cam);

                if (step < steps.length) {
                    setTimeout(callback, interval);
                } else {
                    self.decAnim();
                }
                self.show();
            };
            setTimeout(callback, interval);

        } else { //relative update
            var delta: any = {};
            let frac = 1.0 / nsteps;
            if (mpos) {
                delta.mpos = mpos.clone().sub(curr.mpos).multiplyScalar(frac);
            }
            if (typeof (rz) != 'undefined' && rz != null) {
                delta.rz = frac * (rz - curr.rz);
            }
            if (rot) {
                var next = GLViewer.slerp(curr.rot, rot, frac);
                //comptute step delta rotation
                delta.rot = curr.rot.clone().inverse().multiply(next);
            }
            if (cam) {
                delta.cam = cam.clone().sub(curr.cam).multiplyScalar(frac);
            }
            let step = 0.0;
            let self = this;
            let callback = function () {
                step += 1;
                if (delta.mpos) {
                    self.modelGroup.position.add(delta.mpos);
                }
                if (delta.rz) {
                    self.rotationGroup.position.z += delta.rz;
                }
                if (delta.rot) {
                    self.rotationGroup.quaternion.multiply(delta.rot);
                }
                if (delta.cam) {
                    self.lookingAt.add(delta.cam);
                    self.camera.lookAt(self.lookingAt);
                }

                if (step < nsteps) {
                    setTimeout(callback, interval);
                } else {
                    self.decAnim();
                }
                self.show();
            };
            setTimeout(callback, interval);
        }
    };

    /**
     * Rotate scene by angle degrees around axis
     *
     * @param {number}
     *            [angle] - Angle, in degrees, to rotate by.
     * @param {string}
     *            [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around.
     *            Default "y".  View relative (rather than model relative) axes are prefixed with v.
     *            Axis can also be specified as a vector.
     * @param {number}
     *            [animationDuration] - an optional parameter that denotes
     *            the duration of the rotation animation. Default 0 (no animation)
     * @param {boolean} [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     * @example     $3Dmol.download('cid:4000', viewer, {}, function() {
    viewer.setStyle({stick:{}});
    viewer.zoomTo();
    viewer.rotate(90,'y',1);
    viewer.render(callback);
    });

     *
     */
    public rotate(angle: number, axis: any = "y", animationDuration: number = 0, fixedPath: boolean = false) {

        if (axis == "x") {
            axis = { x: 1, y: 0, z: 0 };
        } else if (axis == "y") {
            axis = { x: 0, y: 1, z: 0 };
        } else if (axis == "z") {
            axis = { x: 0, y: 0, z: 1 };
        }

        //support rotating with respect to view axis, not model
        if (axis == "vx") {
            axis = { vx: 1, vy: 0, vz: 0 };
        } else if (axis == "vy") {
            axis = { vx: 0, vy: 1, vz: 0 };
        } else if (axis == "vz") {
            axis = { vx: 0, vy: 0, vz: 1 };
        }

        if (typeof (axis.vx) !== 'undefined') {
            var vaxis = new Vector3(axis.vx, axis.vy, axis.vz);
            vaxis.applyQuaternion(this.rotationGroup.quaternion);
            axis = { x: vaxis.x, y: vaxis.y, z: vaxis.z };
        }

        var qFromAngle = function (rangle) {
            var s = Math.sin(rangle / 2.0);
            var c = Math.cos(rangle / 2.0);
            var i = 0, j = 0, k = 0;

            i = axis.x * s;
            j = axis.y * s;
            k = axis.z * s;

            return new Quaternion(i, j, k, c).normalize();
        };

        var rangle = Math.PI * angle / 180.0;
        var q = qFromAngle(rangle);

        if (animationDuration) {
            var final = new Quaternion().copy(this.rotationGroup.quaternion).multiply(q);//final
            this.animateMotion(animationDuration, fixedPath,
                this.modelGroup.position,
                this.rotationGroup.position.z,
                final,
                this.lookingAt);
        } else { //not animated
            this.rotationGroup.quaternion.multiply(q);
            this.show();
        }
        return this;

    };

    public surfacesFinished() {
        for (var key in this.surfaces) {
            if (!this.surfaces[key][0].done) {
                return false;
            }
        }
        return true;
    };

    /** Returns an array representing the current viewpoint.
     * Translation, zoom, and rotation quaternion.
     * @returns {Array.<number>} [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y, q.z, q.w ]
     *  */
    public getView() {
        if (!this.modelGroup)
            return [0, 0, 0, 0, 0, 0, 0, 1];
        var pos = this.modelGroup.position;
        var q = this.rotationGroup.quaternion;
        return [pos.x, pos.y, pos.z, this.rotationGroup.position.z, q.x, q.y,
        q.z, q.w];
    };

    /** Sets the view to the specified translation, zoom, and rotation.
     *
     * @param {Array.<number>} arg Array formatted identically to the return value of getView */
    public setView(arg, nolink?) {

        if (arg === undefined ||
            !(arg instanceof Array || arg.length !== 8))
            return this;

        if (!this.modelGroup || !this.rotationGroup)
            return this;
        this.modelGroup.position.x = arg[0];
        this.modelGroup.position.y = arg[1];
        this.modelGroup.position.z = arg[2];
        this.rotationGroup.position.z = arg[3];
        this.rotationGroup.quaternion.x = arg[4];
        this.rotationGroup.quaternion.y = arg[5];
        this.rotationGroup.quaternion.z = arg[6];
        this.rotationGroup.quaternion.w = arg[7];
        if (typeof (arg[8]) != "undefined") {
            this.rotationGroup.position.x = arg[8];
            this.rotationGroup.position.y = arg[9];
        }

        this.show(nolink);
        return this;

    };

    // apply styles, models, etc in viewer
    /**
     * Render current state of viewer, after
     * adding/removing models, applying styles, etc.
     *
     */
    public render(callback?, exts?) {
        this.renderer.setViewport();
        this.updateClickables(); //must render for clickable styles to take effect
        var view = this.getView();

        if (this.stateChangeCallback) {
            //todo: have ability to only send delta updates
            this.stateChangeCallback(this.getInternalState());
        }

        var i, n;
        if (!exts) exts = this.renderer.supportedExtensions();
        for (i = 0; i < this.models.length; i++) {
            if (this.models[i]) {
                this.models[i].globj(this.modelGroup, exts);
            }
        }

        for (i = 0; i < this.shapes.length; i++) {
            if (this.shapes[i]) { //exists
                if ((typeof (this.shapes[i].frame) === 'undefined' || this.viewer_frame < 0 ||
                    this.shapes[i].frame < 0 || this.shapes[i].frame == this.viewer_frame)) {
                    this.shapes[i].globj(this.modelGroup, exts);
                } else { //should not be displayed in current frame
                    this.shapes[i].removegl(this.modelGroup);
                }
            }
        }

        for (i = 0; i < this.labels.length; i++) {
            if (this.labels[i] && typeof (this.labels[i].frame) != 'undefined' && this.labels[i].frame >= 0) { //exists and has frame specifier
                this.modelGroup.remove(this.labels[i].sprite);
                if (this.viewer_frame < 0 || this.labels[i].frame == this.viewer_frame) {
                    this.modelGroup.add(this.labels[i].sprite);
                }
            }
        }

        for (i in this.surfaces) { // this is an object with possible holes
            if (!this.surfaces.hasOwnProperty(i)) continue;
            var surfArr = this.surfaces[i];
            for (n = 0; n < surfArr.length; n++) {
                if (surfArr.hasOwnProperty(n)) {
                    var geo = surfArr[n].geo;
                    // async surface generation can cause
                    // the geometry to be webgl initialized before it is fully
                    // formed; force various recalculations until full surface
                    // is
                    // available
                    if (!surfArr[n].finished || exts.regen) {
                        geo.verticesNeedUpdate = true;
                        geo.elementsNeedUpdate = true;
                        geo.normalsNeedUpdate = true;
                        geo.colorsNeedUpdate = true;
                        geo.buffersNeedUpdate = true;
                        geo.boundingSphere = null;

                        if (surfArr[n].done)
                            surfArr[n].finished = true;

                        // remove partially rendered surface
                        if (surfArr[n].lastGL)
                            this.modelGroup.remove(surfArr[n].lastGL);

                        // create new surface
                        var smesh = null;

                        if (surfArr[n].mat instanceof LineBasicMaterial) {
                            //special case line meshes
                            smesh = new Line(geo, surfArr[n].mat);
                        }
                        else {
                            smesh = new Mesh(geo, surfArr[n].mat);
                        }
                        if (surfArr[n].mat.transparent && surfArr[n].mat.opacity == 0) {
                            //don't bother with hidden surfaces
                            smesh.visible = false;
                        } else {
                            smesh.visible = true;
                        }
                        if (surfArr[n].symmetries.length > 1 ||
                            (surfArr[n].symmetries.length == 1 &&
                                !(surfArr[n].symmetries[n].isIdentity()))) {
                            var j;
                            var tmeshes = new Object3D(); //transformed meshes
                            for (j = 0; j < surfArr[n].symmetries.length; j++) {
                                var tmesh = smesh.clone();
                                tmesh.matrix = surfArr[n].symmetries[j];
                                tmesh.matrixAutoUpdate = false;
                                tmeshes.add(tmesh);
                            }
                            surfArr[n].lastGL = tmeshes;
                            this.modelGroup.add(tmeshes);
                        }
                        else {
                            surfArr[n].lastGL = smesh;
                            this.modelGroup.add(smesh);
                        }
                    } // else final surface already there
                }
            }
        }

        this.setView(view); // Calls show() => renderer render
        if (typeof callback === 'function') {
            callback(this);
        }
        return this;
    };

    /* @param {AtomSelectionSpec|any} sel
     * @return list of models specified by sel
     */
    private getModelList(sel: any): GLModel[] {
        let ms: GLModel[] = [];
        if (typeof sel === 'undefined' || typeof sel.model === "undefined") {
            for (let i = 0; i < this.models.length; i++) {
                if (this.models[i])
                    ms.push(this.models[i]);
            }
        } else { // specific to some models
            let selm: any = sel.model;
            if (!Array.isArray(selm))
                selm = [selm];

            for (let i = 0; i < selm.length; i++) {
                //allow referencing models by order of creation
                if (typeof selm[i] === 'number') {
                    var index = selm[i];
                    //support python backward indexing
                    if (index < 0) index += this.models.length;
                    ms.push(this.models[index]);
                } else {
                    ms.push(selm[i]);
                }
            }
        }

        return ms;
    }
    /**
     *
     * @param {AtomSelectionSpec}
     *            sel
     * @return {AtomSpec[]}
     */
    private getAtomsFromSel(sel: AtomSelectionSpec): AtomSpec[] {
        var atoms = [];
        if (typeof (sel) === "undefined")
            sel = {};

        var ms = this.getModelList(sel);

        for (let i = 0; i < ms.length; i++) {
            atoms = atoms.concat(ms[i].selectedAtoms(sel));
        }

        return atoms;
    }

    /**
     *
     * @param {AtomSpec}
     *            atom
     * @param {AtomSelectionSpec}
     *            sel
     * @return {boolean}
     */
    private atomIsSelected(atom: AtomSpec, sel: AtomSelectionSpec) {
        if (typeof (sel) === "undefined")
            sel = {};

        var ms = this.getModelList(sel);

        for (var i = 0; i < ms.length; i++) {
            if (ms[i].atomIsSelected(atom, sel))
                return true;
        }

        return false;
    }


    /** return list of atoms selected by sel
     *
     * @param {AtomSelectionSpec} sel
     * @return {AtomSpec[]}
     */
    public selectedAtoms(sel: AtomSelectionSpec): AtomSpec[] {
        return this.getAtomsFromSel(sel);
    };

    /**
    * Returns valid values for the specified attribute in the given selection
    * @param {string} attribute
    * @param {AtomSelectionSpec} sel
    * @return {Array.<Object>}
    *
    */
    public getUniqueValues(attribute: string, sel?: AtomSelectionSpec) {
        if (typeof (sel) === "undefined")
            sel = {};
        var atoms = this.getAtomsFromSel(sel);
        var values = {};

        for (var atom in atoms) {
            if (atoms[atom].hasOwnProperty(attribute)) {
                var value = atoms[atom][attribute];
                values[value] = true;
            }
        }

        return Object.keys(values);
    };

    /**
     * Return pdb output of selected atoms (if atoms from pdb input)
     *
     * @param {AtomSelectionSpec} sel - Selection specification specifying model and atom properties to select.  Default: all atoms in viewer
     * @return {string} PDB string of selected atoms
     */
    public pdbData(sel: AtomSelectionSpec) {
        var atoms = this.getAtomsFromSel(sel);
        var ret = "";
        for (var i = 0, n = atoms.length; i < n; ++i) {
            ret += atoms[i].pdbline + "\n";
        }
        return ret;
    };


    /**
     * Zoom current view by a constant factor
     *
     * @param {number}
     *            [factor] - Magnification factor. Values greater than 1
     *            will zoom in, less than one will zoom out. Default 2.
     * @param {number}
     *            [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param {Boolean} [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation
     * @example
    $3Dmol.get('data/4csv.pdb', function(data) {
    viewer.addModel(data,'pdb');
    viewer.setStyle({cartoon:{},stick:{}});
    viewer.zoomTo()
    viewer.zoom(2,1000);
    viewer.render();
    });

         */
    public zoom(factor: number = 2, animationDuration: number = 0, fixedPath: boolean = false) {
        var scale = (this.CAMERA_Z - this.rotationGroup.position.z) / factor;
        var final_z = this.CAMERA_Z - scale;

        if (animationDuration > 0) {
            this.animateMotion(animationDuration, fixedPath,
                this.modelGroup.position,
                this.adjustZoomToLimits(final_z),
                this.rotationGroup.quaternion,
                this.lookingAt);
        } else { //no animation
            this.rotationGroup.position.z = this.adjustZoomToLimits(final_z);
            this.show();
        }
        return this;
    };

    /**
     * Translate current view by x,y screen coordinates
     * This pans the camera rather than translating the model.
     *
     * @param {number} x Relative change in view coordinates of camera
     * @param {number} y Relative change in view coordinates of camera
     * @param {number}
     *            [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param {Boolean} [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     * @example     $3Dmol.get('data/4csv.pdb', function(data) {
    viewer.addModel(data,'pdb');
    viewer.setStyle({cartoon:{},stick:{}});
    viewer.zoomTo();
    viewer.translate(200,50);
    viewer.rotate(90,'z');
    viewer.render(callback);
    });
     */
    public translate(x: number, y: number, animationDuration: number = 0, fixedPath: boolean = false) {
        var dx = x / this.WIDTH;
        var dy = y / this.HEIGHT;
        var v = new Vector3(0, 0, -this.CAMERA_Z);

        this.projector.projectVector(v, this.camera);
        v.x -= dx;
        v.y -= dy;
        this.projector.unprojectVector(v, this.camera);
        v.z = 0;

        var final_position = this.lookingAt.clone().add(v);
        if (animationDuration > 0) {
            this.animateMotion(animationDuration, fixedPath,
                this.modelGroup.position,
                this.rotationGroup.position.z,
                this.rotationGroup.quaternion,
                final_position);
        } else { //no animation
            this.lookingAt = final_position;
            this.camera.lookAt(this.lookingAt);
            this.show();
        }
        return this;
    };

    /**
     * Translate current models by x,y screen coordinates
     * This translates the models relative to the current view. It does
     * not change the center of rotation.
     *
     * @param {number} x Relative change in x screen coordinate
     * @param {number} y Relative change in y screen coordinate
     * @param {number}
     *            [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param {Boolean} [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     * @example     $3Dmol.get('data/4csv.pdb', function(data) {
    viewer.addModel(data,'pdb');
    viewer.setStyle({cartoon:{},stick:{}});
    viewer.zoomTo();
    viewer.translateScene(200,50);
    viewer.rotate(90,'z'); // will no longer be around model center
    viewer.render(callback);
    });
     */
    public translateScene(x: number, y: number, animationDuration: number = 0, fixedPath = false) {

        var t = this.screenOffsetToModel(x, y);
        var final_position = this.modelGroup.position.clone().add(t);

        if (animationDuration > 0) {
            this.animateMotion(animationDuration, fixedPath,
                this.modelGroup.position,
                this.rotationGroup.position.z,
                this.rotationGroup.quaternion,
                this.lookingAt);
        } else { //no animation
            this.modelGroup.position = final_position;
            this.show();
        }
        return this;
    };

    /**
     * Adjust slab to fully enclose selection (default everything).
     *
     * @param {AtomSelectionSpec} sel
     *            Selection specification specifying model and atom
     *            properties to select. Default: all atoms in viewer
     */
    public fitSlab(sel: AtomSelectionSpec) {
        sel = sel || {};
        var atoms = this.getAtomsFromSel(sel);
        var tmp = getExtent(atoms);

        // fit to bounding box
        var x = tmp[1][0] - tmp[0][0],
            y = tmp[1][1] - tmp[0][1],
            z = tmp[1][2] - tmp[0][2];

        var maxD = Math.sqrt(x * x + y * y + z * z);
        if (maxD < 5)
            maxD = 5;

        // use full bounding box for slab/fog
        this.slabNear = -maxD / 1.9;
        this.slabFar = maxD / 2;

        return this;
    };

    /**
     * Re-center the viewer around the provided selection (unlike zoomTo, does not zoom).
     *
     * @param {AtomSelectionSpec}
     *            [sel] - Selection specification specifying model and atom
     *            properties to select. Default: all atoms in viewer
     * @param {number}
     *            [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param {Boolean} [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
     * @example // if the user were to pass the animationDuration value to
     *           // the function like so viewer.zoomTo({resn:'STI'},1000);
     *         //   the program would center on resn 'STI' over the course
     *         //   of 1 second(1000 milleseconds).
     *  // Reposition to centroid of all atoms of all models in this
     * //viewer glviewer.center();
    $3Dmol.get('data/4csv.pdb', function(data) {
    viewer.addModel(data,'pdb');
    viewer.setStyle({cartoon:{},stick:{}});
    viewer.center();
    viewer.render(callback);
    });
     */
    public center(sel: AtomSelectionSpec = {}, animationDuration: number = 0, fixedPath: boolean = false) {
        var allatoms, alltmp;
        var atoms = this.getAtomsFromSel(sel);
        var tmp = getExtent(atoms);

        if (isEmptyObject(sel)) {
            //include shapes when zooming to full scene
            //TODO: figure out a good way to specify shapes as part of a selection
            this.shapes.forEach((shape) => {
                if (shape && shape.boundingSphere && shape.boundingSphere.center) {
                    var c = shape.boundingSphere.center;
                    var r = shape.boundingSphere.radius;
                    if (r > 0) {
                        //make sure full shape is visible
                        atoms.push(new Vector3(c.x + r, c.y, c.z));
                        atoms.push(new Vector3(c.x - r, c.y, c.z));
                        atoms.push(new Vector3(c.x, c.y + r, c.z));
                        atoms.push(new Vector3(c.x, c.y - r, c.z));
                        atoms.push(new Vector3(c.x, c.y, c.z + r));
                        atoms.push(new Vector3(c.x, c.y, c.z - r));
                    } else {
                        atoms.push(c);
                    }
                }
            });
            tmp = getExtent(atoms);
            allatoms = atoms;
            alltmp = tmp;

        }
        else {
            allatoms = this.getAtomsFromSel({});
            alltmp = getExtent(allatoms);
        }

        // use selection for center
        var center = new Vector3(tmp[2][0], tmp[2][1], tmp[2][2]);

        // but all for bounding box
        var x = alltmp[1][0] - alltmp[0][0], y = alltmp[1][1] -
            alltmp[0][1], z = alltmp[1][2] - alltmp[0][2];

        var maxD = Math.sqrt(x * x + y * y + z * z);
        if (maxD < 5)
            maxD = 5;

        // use full bounding box for slab/fog
        this.slabNear = -maxD / 1.9;
        this.slabFar = maxD / 2;

        // for zoom, use selection box
        x = tmp[1][0] - tmp[0][0];
        y = tmp[1][1] - tmp[0][1];
        z = tmp[1][2] - tmp[0][2];
        maxD = Math.sqrt(x * x + y * y + z * z);
        if (maxD < 5)
            maxD = 5;

        //find the farthest atom from center to get max distance needed for view
        var maxDsq = 25;
        for (var i = 0; i < atoms.length; i++) {
            if (atoms[i]) {
                var dsq = center.distanceToSquared(atoms[i] as XYZ);
                if (dsq > maxDsq)
                    maxDsq = dsq;
            }
        }

        maxD = Math.sqrt(maxDsq) * 2;
        var finalpos = center.clone().multiplyScalar(-1);
        if (animationDuration > 0) {
            this.animateMotion(animationDuration, fixedPath,
                finalpos,
                this.rotationGroup.position.z,
                this.rotationGroup.quaternion,
                this.lookingAt);
        } else { //no animation
            this.modelGroup.position = finalpos;
            this.show();
        }
        return this;
    };

    /**
     * Zoom to center of atom selection.  The slab will be set appropriately for
     * the selection, unless an empty selection is provided, in which case there will be no slab.
     *
     * @param {Object}
     *            [sel] - Selection specification specifying model and atom
     *            properties to select. Default: all atoms in viewer
     * @param {number}
     *            [animationDuration] - an optional parameter that denotes
     *            the duration of a zoom animation
     * @param {Boolean} [fixedPath] - if true animation is constrained to
     *      requested motion, overriding updates that happen during the animation         *
      * @example


          $3Dmol.get('data/1fas.pqr', function(data){
              viewer.addModel(data, "pqr");
              $3Dmol.get("data/1fas.cube",function(volumedata){
                  viewer.addSurface($3Dmol.SurfaceType.VDW, {
                      opacity:0.85,
                      voldata: new $3Dmol.VolumeData(volumedata, "cube"),
                      volscheme: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'charge'))
                  },{});

              viewer.render();
              });
              viewer.zoomTo();
            });
     */
    public zoomTo(sel: AtomSelectionSpec = {}, animationDuration: number = 0, fixedPath: boolean = false) {
        let atoms = this.getAtomsFromSel(sel);
        let atombox = getExtent(atoms);
        let allbox = atombox;

        if (isEmptyObject(sel)) {
            //include shapes when zooming to full scene
            //TODO: figure out a good way to specify shapes as part of a selection
            let natoms = atoms && atoms.length;
            this.shapes.forEach((shape) => {
                if (shape && shape.boundingSphere) {
                    if (shape.boundingSphere.box) {
                        let box = shape.boundingSphere.box;
                        atoms.push(new Vector3(box.min.x, box.min.y, box.min.z));
                        atoms.push(new Vector3(box.max.x, box.max.y, box.max.z));
                    } else if (shape.boundingSphere.center) {
                        var c = shape.boundingSphere.center;
                        var r = shape.boundingSphere.radius;
                        if (r > 0) {
                            //make sure full shape is visible
                            atoms.push(new Vector3(c.x + r, c.y, c.z));
                            atoms.push(new Vector3(c.x - r, c.y, c.z));
                            atoms.push(new Vector3(c.x, c.y + r, c.z));
                            atoms.push(new Vector3(c.x, c.y - r, c.z));
                            atoms.push(new Vector3(c.x, c.y, c.z + r));
                            atoms.push(new Vector3(c.x, c.y, c.z - r));
                        } else {
                            atoms.push(c);
                        }
                    }
                }
            });
            allbox = getExtent(atoms);
            if (!natoms) { //if no atoms, use shapes for center
                for (let i = 0; i < 3; i++) { //center of bounding box
                    atombox[2][i] = (allbox[0][i] + allbox[1][i]) / 2;
                }
            }
        } else { //include all atoms in slab calculation
            let allatoms = this.getAtomsFromSel({});
            allbox = getExtent(allatoms);
        }

        // use selection for center
        var center = new Vector3(atombox[2][0], atombox[2][1], atombox[2][2]);

        // but all for bounding box
        var x = allbox[1][0] - allbox[0][0], y = allbox[1][1]
            - allbox[0][1], z = allbox[1][2] - allbox[0][2];

        var maxD = Math.sqrt(x * x + y * y + z * z);
        if (maxD < 5)
            maxD = 5;

        // use full bounding box for slab/fog
        this.slabNear = -maxD / 1.9;
        this.slabFar = maxD / 2;

        //if we are selecting everything, have ver permissive slab
        //can't do "infinity" size since this will break orthographic
        if (Object.keys(sel).length === 0) {
            this.slabNear = Math.min(-maxD * 2, -50);
            this.slabFar = Math.max(maxD * 2, 50);
        }

        // keep at least this much space in view
        var MAXD = this.config.minimumZoomToDistance || 5;
        // for zoom, use selection box
        x = atombox[1][0] - atombox[0][0];
        y = atombox[1][1] - atombox[0][1];
        z = atombox[1][2] - atombox[0][2];
        maxD = Math.sqrt(x * x + y * y + z * z);
        if (maxD < MAXD)
            maxD = MAXD;

        //find the farthest atom from center to get max distance needed for view
        var maxDsq = MAXD * MAXD;
        for (var i = 0; i < atoms.length; i++) {
            if (atoms[i]) {
                var dsq = center.distanceToSquared(atoms[i] as XYZ);
                if (dsq > maxDsq)
                    maxDsq = dsq;
            }
        }

        maxD = Math.sqrt(maxDsq) * 2;
        var finalpos = center.clone().multiplyScalar(-1);
        var finalz = -(maxD * 0.5
            / Math.tan(Math.PI / 180.0 * this.camera.fov / 2) - this.CAMERA_Z);

        finalz = this.adjustZoomToLimits(finalz);
        if (animationDuration > 0) {
            this.animateMotion(animationDuration, fixedPath,
                finalpos,
                finalz,
                this.rotationGroup.quaternion,
                this.lookingAt);
        } else {
            this.modelGroup.position = finalpos;
            this.rotationGroup.position.z = finalz;
            this.show();
        }
        return this;

    };

    /**
     * Set slab of view (contents outside of slab are clipped).
     * Must call render to update.
     *
     * @param {number} near near clipping plane distance
     * @param {number} far far clipping plane distance
     */
    public setSlab(near: number, far: number) {
        this.slabNear = near;
        this.slabFar = far;
    };

    /**
     * Get slab of view (contents outside of slab are clipped).
     *
     * @return {Object}
     *      @property {number} near - near clipping plane distance
     *      @property {number} far - far clipping plane distance
     */
    public getSlab() {
        return { near: this.slabNear, far: this.slabFar };
    };

    /**
     * Add label to viewer
     *
     * @param {string}
     *            text - Label text
     * @param {LabelSpec}
     *            options - Label style specification
      @param {AtomSelection}
     *            sel - Set position of label to center of this selection
     * @param {boolean} noshow - if true, do not immediately display label - when adding multiple labels this is more efficient
     * @return {Label}
     *
     * @example
     *  $3Dmol.download("pdb:2EJ0",viewer,{},function(){

              viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
              viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                       borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                       position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
              viewer.setStyle({chain:'A'},{cross:{hidden:true}});
              viewer.setStyle({chain:'B'},{cross:{hidden:false,
                                                  linewidth:1.0,
                                                  colorscheme:'greenCarbon'}});
              viewer.setStyle({chain:'C'},{cross:{hidden:false,
                                                  linewidth:1.0,
                                                  radius:0.5}});
              viewer.setStyle({chain:'D'},{cross:{hidden:false,
                                                  linewidth:10.0}});
              viewer.setStyle({chain:'E'},{cross:{hidden:false,
                                                  linewidth:1.0,
                                                  color:'black'}});

              viewer.render();


            });

     */
    public addLabel(text: string, options: LabelSpec = {}, sel?: AtomSelectionSpec, noshow: boolean = false) {
        if (sel) {
            var extent = getExtent(this.getAtomsFromSel(sel));
            options.position = { x: extent[2][0], y: extent[2][1], z: extent[2][2] };
        }
        var label = new Label(text, options);
        label.setContext();
        this.modelGroup.add(label.sprite);
        this.labels.push(label);

        if (!noshow) this.show();
        return label;
    };



    /** Add residue labels.  This will generate one label per a
     * residue within the selected atoms.  The label will be at the
     * centroid of the atoms and styled according to the passed style.
     * The label text will be [resn][resi]
     *
     * @param {AtomSelectionSpec} sel
     * @param {AtomStyleSpec} style
     * @param {boolean} byframe - if true, create labels for every individual frame, not just current
     *
     * @example
         $3Dmol.download("mmtf:2ll5",viewer,{},function(){
              viewer.setStyle({stick:{radius:0.15},cartoon:{}});
              viewer.addResLabels({hetflag:false}, {font: 'Arial', fontColor:'black',showBackground:false, screenOffset: {x:0,y:0}});
              viewer.zoomTo();
              viewer.render();
            });
     */
    public addResLabels(sel: AtomSelectionSpec, style: AtomStyleSpec, byframe: boolean = false) {
        let start = this.labels.length;
        this.applyToModels("addResLabels", sel, this, style, byframe);
        this.show();
        return this.labels.slice(start);
    };

    /** Add property labels.  This will generate one label per a selected
     * atom at the atom's coordinates with the property value as the label text.
     *
     * @param {string} prop - property name
     * @param {AtomSelectionSpec} sel
     * @param {AtomStyleSpec} style
     *
     * * @example
         $3Dmol.download("cid:5291",viewer,{},function(){
              viewer.setStyle({stick: {radius:.2}});
              viewer.addPropertyLabels("index",{not:{elem:'H'}}, {fontColor:'black',font: 'sans-serif', fontSize: 28, showBackground:false,alignment:'center'});
              viewer.zoomTo();
              viewer.render();
            });
     */
    public addPropertyLabels(prop: string, sel: AtomSelectionSpec, style: AtomStyleSpec) {
        this.applyToModels("addPropertyLabels", prop, sel, this, style);
        this.show();
        return this;
    };

    /**
     * Remove label from viewer
     *
     * @param {Label} label - $3Dmol label
     *
     * @example // Remove labels created in
     $3Dmol.download("pdb:2EJ0",viewer,{},function(){
              var toremove = viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
              viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                       borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                       position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
              viewer.removeLabel(toremove);
              viewer.render();


            });

     */
    public removeLabel(label: Label) {
        //todo: don't do the linear search
        for (var i = 0; i < this.labels.length; i++) {
            if (this.labels[i] == label) {
                this.labels.splice(i, 1);
                label.dispose();
                this.modelGroup.remove(label.sprite);
                break;
            }
        }
        this.show();
        return this;
    };

    /**
     * Remove all labels from viewer
     *
     *         @example
    $3Dmol.download("pdb:1ubq",viewer,{},function(){

           viewer.addResLabels();
           viewer.setStyle({},{stick:{}});
           viewer.render( ); //show labels

           viewer.removeAllLabels();
           viewer.render(); //hide labels
    });
     */
    public removeAllLabels() {
        for (var i = 0; i < this.labels.length; i++) {
            if (this.labels[i] && this.labels[i].sprite) {
                this.modelGroup.remove(this.labels[i].sprite);
            }
        }
        this.labels.splice(0, this.labels.length); //don't overwrite in case linked
        this.show();
        return this;
    };

    // Modify label style
    /**
     * Modify existing label's style
     *
     * @param {Label} label - $3Dmol label
     * @param {LabelSpec}
     *            stylespec - Label style specification
     * @return {Label}
     */
    public setLabelStyle(label: Label, stylespec: LabelSpec) {
        this.modelGroup.remove(label.sprite);
        label.dispose();
        label.stylespec = stylespec;
        label.setContext();
        this.modelGroup.add(label.sprite);
        this.show();
        return label;

    };

    // Change label text
    /**
     * Modify existing label's text
     *
     * @param {Label}  label - $3Dmol label
     * @param {String}
     *            text - Label text
     * @return {Label}
     */
    public setLabelText(label: Label, text: string) {
        this.modelGroup.remove(label.sprite);
        label.dispose();
        label.text = text;
        label.setContext();
        this.modelGroup.add(label.sprite);
        this.show();
        return label;

    };

    /**
     * Add shape object to viewer
     * @see {GLShape}
     *
     * @param {ShapeSpec} shapeSpec - style specification for label
     * @return {GLShape}
     */
    public addShape(shapeSpec: ShapeSpec) {
        shapeSpec = shapeSpec || {};
        var shape = new GLShape(shapeSpec);
        shape.shapePosition = this.shapes.length;
        this.shapes.push(shape);

        return shape;
    };

    /**
     * Remove shape object from viewer
     *
     * @param {GLShape} shape - Reference to shape object to remove
     */
    public removeShape(shape: GLShape) {
        if (!shape)
            return this;
        shape.removegl(this.modelGroup);
        delete this.shapes[shape.shapePosition];
        // clear off back of model array
        while (this.shapes.length > 0
            && typeof (this.shapes[this.shapes.length - 1]) === "undefined")
            this.shapes.pop();
        return this;
    };

    /**
     * Remove all shape objects from viewer
     */
    public removeAllShapes() {
        for (var i = 0; i < this.shapes.length; i++) {
            var shape = this.shapes[i];
            if (shape) shape.removegl(this.modelGroup);
        }
        this.shapes.splice(0, this.shapes.length);
        return this;
    };

    //gets the center of the selection
    private getSelectionCenter(spec: AtomSelectionSpec): XYZ {
        if (spec.hasOwnProperty("x") && spec.hasOwnProperty("y") && spec.hasOwnProperty("z"))
            return spec as XYZ;
        var atoms = this.getAtomsFromSel(spec);
        if (atoms.length == 0)
            return { x: 0, y: 0, z: 0 };

        var extent = getExtent(atoms);
        return { x: extent[0][0] + (extent[1][0] - extent[0][0]) / 2, y: extent[0][1] + (extent[1][1] - extent[0][1]) / 2, z: extent[0][2] + (extent[1][2] - extent[0][2]) / 2 };
    };

    /**
     * Create and add sphere shape. This method provides a shorthand
     * way to create a spherical shape object
     *
     * @param {SphereShapeSpec} spec - Sphere shape style specification
     * @return {GLShape}
     @example

     viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});

     viewer.render();
     */
    public addSphere(spec: SphereSpec) {
        spec = spec || {};

        spec.center = this.getSelectionCenter(spec.center);

        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        s.addSphere(spec);
        this.shapes.push(s);
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended
        return s;
    };

    /**
     * Create and add box shape. This method provides a shorthand
     * way to create a box shape object
     *
     * @param {BoxSpec} spec - Box shape style specification
     * @return {GLShape}
     @example

     viewer.addLine({color:'red',start:{x:0,y:0,z:0},end:{x:5,y:0,z:0}});
     viewer.addLine({color:'blue',start:{x:0,y:0,z:0},end:{x:0,y:5,z:0}});
     viewer.addLine({color:'green',start:{x:0,y:0,z:0},end:{x:0,y:0,z:5}});

     viewer.addBox({center:{x:0,y:0,z:0},dimensions: {w:3,h:4,d:2},color:'magenta'});
     viewer.zoomTo();
     viewer.rotate(45, {x:1,y:1,z:1});
     viewer.render();
     */
    public addBox(spec: BoxSpec = {}) {

        if (spec.corner != undefined) {
            spec.corner = this.getSelectionCenter(spec.corner);
        }
        if (spec.center != undefined) {
            spec.center = this.getSelectionCenter(spec.center);
        }

        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        s.addBox(spec);
        this.shapes.push(s);
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

        return s;
    };

    /**
     * Create and add arrow shape
     *
     * @param {ArrowSpec} spec - Style specification
     * @return {GLShape}
     @example
    $3Dmol.download("pdb:4DM7",viewer,{},function(){

              viewer.setBackgroundColor(0xffffffff);
              viewer.addArrow({
                  start: {x:-10.0, y:0.0, z:0.0},
                  end: {x:0.0, y:-10.0, z:0.0},
                  radius: 1.0,
                  radiusRadio:1.0,
                  mid:1.0,
                  clickable:true,
                  callback:function(){
                      this.color.setHex(0xFF0000FF);
                      viewer.render( );
                  }
              });
              viewer.render();
            });
     */
    public addArrow(spec: ArrowSpec = {}) {

        spec.start = this.getSelectionCenter(spec.start);
        spec.end = this.getSelectionCenter(spec.end);

        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        s.addArrow(spec);
        this.shapes.push(s);
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

        return s;
    };

    /**
     * Create and add cylinder shape
     *
     * @param {CylinderSpec} spec - Style specification
     * @return {GLShape}

      @example
     viewer.setBackgroundColor(0xffffffff);
          viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
                              end:{x:10.0,y:0.0,z:0.0},
                              radius:1.0,
                              fromCap:1,
                              toCap:2,
                              color:'red',
                              hoverable:true,
                              clickable:true,
                              callback:function(){ this.color.setHex(0x00FFFF00);viewer.render( );},
                              hover_callback: function(){ viewer.render( );},
                              unhover_callback: function(){ this.color.setHex(0xFF000000);viewer.render( );}
                             });
          viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
                              end:{x:0.0,y:10.0,z:0.0},
                              radius:0.5,
                              fromCap:false,
                              toCap:true,
                              color:'teal'});
          viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
                              end:{x:20.0,y:0.0,z:0.0},
                              radius:1.0,
                              color:'black',
                              fromCap:false,
                              toCap:false});
          viewer.render();
     */
    public addCylinder(spec: CylinderSpec = {}) {

        spec.start = this.getSelectionCenter(spec.start);
        spec.end = this.getSelectionCenter(spec.end);

        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        if (spec.dashed)
            s.addDashedCylinder(spec);
        else
            s.addCylinder(spec);
        this.shapes.push(s);
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

        return s;
    };

    /**
     * Create and add Curve shape
     *
     * @param {CurveSpec} spec - Style specification
     * @return {GLShape}

     @example
          viewer.addCurve({points: [{x:0.0,y:0.0,z:0.0}, {x:5.0,y:3.0,z:0.0}, {x:5.0,y:7.0,z:0.0}, {x:0.0,y:10.0,z:0.0}],
                              radius:0.5,
                              smooth: 10,
                              fromArrow:false,
                              toArrow: true,
                              color:'orange',
                              });
          viewer.addCurve({points: [{x:-1,y:0.0,z:0.0}, {x:-5.0,y:5.0,z:0.0}, {x:-2,y:10.0,z:0.0}],
                              radius:1,
                              fromArrow:true,
                              toArrow: false,
                              color:'purple',
                              });
          viewer.zoomTo();
          viewer.render();
     */
    public addCurve(spec: CurveSpec = {}) {
        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        s.addCurve(spec);
        this.shapes.push(s);
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

        return s;
    };


    /**
     * Create and add line shape
     *
     * @param {LineSpec} spec - Style specification, can specify dashed, dashLength, and gapLength
     * @return {GLShape}
     @example
     $3Dmol.download("pdb:2ABJ",viewer,{},function(){

              viewer.setViewStyle({style:"outline"});
              viewer.setStyle({chain:'A'},{sphere:{hidden:true}});
              viewer.setStyle({chain:'D'},{sphere:{radius:3.0}});
              viewer.setStyle({chain:'G'},{sphere:{colorscheme:'greenCarbon'}});
              viewer.setStyle({chain:'J'},{sphere:{color:'blue'}});
              viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
              viewer.render();
          });

     */
    public addLine(spec: LineSpec = {}) {

        spec.start = this.getSelectionCenter(spec.start);
        spec.end = this.getSelectionCenter(spec.end);

        spec.wireframe = true;
        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        if (spec.dashed)
            s = this.addLineDashed(spec, s);
        else
            s.addLine(spec);
        this.shapes.push(s);
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

        return s;
    };


    /**
     * Create and add unit cell visualization.
     *
     * @param {GLModel|number} model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
     * @param {UnitCellStyleSpec} spec - visualization style
       @example

            $3Dmol.get('data/1jpy.cif', function(data) {
              let m = viewer.addModel(data);
              viewer.addUnitCell(m, {box:{color:'purple'},alabel:'X',blabel:'Y',clabel:'Z',alabelstyle: {fontColor: 'black',backgroundColor:'white',inFront:true,fontSize:40},astyle:{color:'darkred', radius:5,midpos: -10}});
              viewer.zoomTo();
              viewer.render();
    });
     */
    public addUnitCell(model?: GLModel | number, spec?: UnitCellStyleSpec) {
        model = this.getModel(model);
        spec = spec || { alabel: 'a', blabel: 'b', clabel: 'c' };

        spec.box = spec.box || {};
        spec.astyle = spec.astyle || { color: 'red', radius: 0.1, midpos: -1 };
        spec.bstyle = spec.bstyle || { color: 'green', radius: 0.1, midpos: -1 };
        spec.cstyle = spec.cstyle || { color: 'blue', radius: 0.1, midpos: -1 };
        spec.alabelstyle = spec.alabelstyle || { fontColor: 'red', showBackground: false, alignment: 'center', inFront: false };
        spec.blabelstyle = spec.blabelstyle || { fontColor: 'green', showBackground: false, alignment: 'center', inFront: false };
        spec.clabelstyle = spec.clabelstyle || { fontColor: 'blue', showBackground: false, alignment: 'center', inFront: false };

        //clear any previous box
        if (model.unitCellObjects) {
            this.removeUnitCell(model);
        }
        model.unitCellObjects = { shapes: [], labels: [] };
        //calculate points
        var data = model.getCrystData();
        var matrix = null;
        if (data) {

            if (data.matrix) {
                matrix = data.matrix;
            } else {
                var a = data.a, b = data.b, c = data.c, alpha = data.alpha, beta = data.beta, gamma = data.gamma;
                alpha = alpha * Math.PI / 180.0;
                beta = beta * Math.PI / 180.0;
                gamma = gamma * Math.PI / 180.0;

                var u, v, w;

                u = Math.cos(beta);
                v = (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma)) / Math.sin(gamma);
                w = Math.sqrt(Math.max(0, 1 - u * u - v * v));

                matrix = new Matrix3(a, b * Math.cos(gamma), c * u,
                    0, b * Math.sin(gamma), c * v,
                    0, 0, c * w);
            }

            var points = [new Vector3(0, 0, 0),
            new Vector3(1, 0, 0),
            new Vector3(0, 1, 0),
            new Vector3(0, 0, 1),
            new Vector3(1, 1, 0),
            new Vector3(0, 1, 1),
            new Vector3(1, 0, 1),
            new Vector3(1, 1, 1)];

            // console.log('Matrix4', data.matrix4, data.matrix);
            if (data.matrix4) {
                for (let i = 0; i < points.length; i++) {
                    if (data.size) points[i].multiplyVectors(points[i], data.size); //matrix is for unit vectors, not whole box
                    points[i] = points[i].applyMatrix4(data.matrix4);
                }
            } else {
                for (let i = 0; i < points.length; i++) {
                    points[i] = points[i].applyMatrix3(matrix);
                }
            }

            //draw box
            if (spec.box && !spec.box.hidden) {
                spec.box.wireframe = true;
                var s = new GLShape(spec.box);
                s.shapePosition = this.shapes.length;

                s.addLine({ start: points[0], end: points[1] });
                s.addLine({ start: points[0], end: points[2] });
                s.addLine({ start: points[1], end: points[4] });
                s.addLine({ start: points[2], end: points[4] });

                s.addLine({ start: points[0], end: points[3] });
                s.addLine({ start: points[3], end: points[5] });
                s.addLine({ start: points[2], end: points[5] });

                s.addLine({ start: points[1], end: points[6] });
                s.addLine({ start: points[4], end: points[7] });
                s.addLine({ start: points[6], end: points[7] });

                s.addLine({ start: points[3], end: points[6] });
                s.addLine({ start: points[5], end: points[7] });

                this.shapes.push(s);
                model.unitCellObjects.shapes.push(s);
                s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended
            }

            //draw arrows
            if (!spec.astyle.hidden) {
                spec.astyle.start = points[0];
                spec.astyle.end = points[1];
                let arrow = this.addArrow(spec.astyle);
                model.unitCellObjects.shapes.push(arrow);
            }

            if (!spec.bstyle.hidden) {
                spec.bstyle.start = points[0];
                spec.bstyle.end = points[2];
                let arrow = this.addArrow(spec.bstyle);
                model.unitCellObjects.shapes.push(arrow);
            }

            if (!spec.cstyle.hidden) {
                spec.cstyle.start = points[0];
                spec.cstyle.end = points[3];
                let arrow = this.addArrow(spec.cstyle);
                model.unitCellObjects.shapes.push(arrow);
            }

            if (spec.alabel) {
                spec.alabelstyle.position = points[1];
                let label = this.addLabel(spec.alabel, spec.alabelstyle);
                model.unitCellObjects.labels.push(label);

            }
            if (spec.blabel) {
                spec.blabelstyle.position = points[2];
                let label = this.addLabel(spec.blabel, spec.blabelstyle);
                model.unitCellObjects.labels.push(label);
            }
            if (spec.clabel) {
                spec.clabelstyle.position = points[3];
                let label = this.addLabel(spec.clabel, spec.clabelstyle);
                model.unitCellObjects.labels.push(label);
            }

        }

    };

    /**
    * Remove unit cell visualization from model.
    *
    * @param {GLModel|number} model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
      @example
           $3Dmol.get('data/icsd_200866.cif', function(data) {
             let m = viewer.addModel(data);
             viewer.setStyle({sphere:{}})
             viewer.addUnitCell();
             viewer.zoomTo();
             viewer.removeUnitCell();
             viewer.render();
       });
    */
    public removeUnitCell(model?: GLModel | number) {
        model = this.getModel(model);
        if (model.unitCellObjects) {
            let viewer = this;
            model.unitCellObjects.shapes.forEach(function (s) { viewer.removeShape(s); });
            model.unitCellObjects.labels.forEach(function (l) { viewer.removeLabel(l); });
        }
        delete model.unitCellObjects;
    };

    /**
    * Replicate atoms in model to form a super cell of the specified dimensions.
    * Original cell will be centered as much as possible.
    *
    * @param {integer} A - number of times to replicate cell in X dimension.
    * @param {integer} B - number of times to replicate cell in Y dimension.  If absent, X value is used.
    * @param {integer} C - number of times to replicate cell in Z dimension.  If absent, Y value is used.
    * @param {GLModel} model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
    * @param {boolean} addBonds - Create bonds between unit cells based on distances.
      @example
           $3Dmol.get('data/icsd_200866.cif', function(data) {
             let m = viewer.addModel(data);
             viewer.setStyle({sphere:{scale:.25}})
             viewer.addUnitCell();
             viewer.zoomTo();
             viewer.replicateUnitCell(3,2,1,m);
             viewer.render();
       });
    */
    public replicateUnitCell(A: number = 3, B: number = A, C: number = B, model?: GLModel | number, addBonds?: boolean) {
        model = this.getModel(model);
        let cryst = model.getCrystData();
        if (cryst) {
            const atoms = model.selectedAtoms({});
            const matrix = cryst.matrix;
            let makeoff = function (I) {
                //alternate around zero: 1,-1,2,-2...
                if (I % 2 == 0) return -I / 2;
                else return Math.ceil(I / 2);
            };

            for (let i = 0; i < A; i++) {
                for (let j = 0; j < B; j++) {
                    for (let k = 0; k < C; k++) {
                        if (i == 0 && j == 0 && k == 0) continue; //actual unit cell
                        let offset = new Vector3(makeoff(i), makeoff(j), makeoff(k));
                        offset.applyMatrix3(matrix);

                        let newatoms = [];
                        for (let a = 0; a < atoms.length; a++) {
                            let newAtom: any = {};
                            for (let p in atoms[a]) {
                                newAtom[p] = atoms[a][p];
                            }
                            newAtom.x += offset.x;
                            newAtom.y += offset.y;
                            newAtom.z += offset.z;
                            newatoms.push(newAtom);
                        }
                        model.addAtoms(newatoms);
                    }
                }
            }

            if(addBonds) {
                model.assignBonds();
            }
        }
    };

    /** Add dashed line to shape */
    public addLineDashed(spec: CylinderSpec, s: GLShape) {
        spec.dashLength = spec.dashLength || 0.5;
        spec.gapLength = spec.gapLength || 0.5;

        var p1: Vector3;
        if (!spec.start) {
            p1 = new Vector3(0, 0, 0);
        } else {
            p1 = new Vector3(spec.start.x || 0,
                spec.start.y || 0, spec.start.z || 0);
        }

        var p2: Vector3;
        if(!spec.end) p2 = new Vector3(0,0,0);
        else p2 = new Vector3(spec.end.x, spec.end.y || 0, spec.end.z || 0);

        var dir = new Vector3();
        var dash = new Vector3();
        var gap = new Vector3();
        var length, dashAmt, gapAmt;
        var temp = p1.clone();
        var drawn = 0;

        dir.subVectors(p2, p1);
        length = dir.length();
        dir.normalize();
        dash = dir.clone();
        gap = dir.clone();
        dash.multiplyScalar(spec.dashLength);
        gap.multiplyScalar(spec.gapLength);
        dashAmt = dash.length();
        gapAmt = gap.length();

        while (drawn < length) {
            if ((drawn + dashAmt) > length) {
                spec.start = p1;
                spec.end = p2;
                s.addLine(spec);
                break;
            }
            temp.addVectors(p1, dash);
            spec.start = p1;
            spec.end = temp;
            s.addLine(spec);
            p1 = temp.clone();
            drawn += dashAmt;

            temp.addVectors(p1, gap);
            p1 = temp.clone();
            drawn += gapAmt;
        }
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

        return s;
    }

    /**
     * Add custom shape component from user supplied function
     *
     * @param {CustomSpec} spec - Style specification
     * @return {GLShape}
     @example
     function triangle(viewer) {
        var vertices = [];
        var normals = [];
        var colors = [];
        var r = 20;
        //triangle
        vertices.push(new $3Dmol.Vector3(0,0,0));
        vertices.push(new $3Dmol.Vector3(r,0,0));
        vertices.push(new $3Dmol.Vector3(0,r,0));

        normals.push(new $3Dmol.Vector3(0,0,1));
        normals.push(new $3Dmol.Vector3(0,0,1));
        normals.push(new $3Dmol.Vector3(0,0,1));

        colors.push({r:1,g:0,b:0});
        colors.push({r:0,g:1,b:0});
        colors.push({r:0,g:0,b:1});

        var faces = [ 0,1,2 ];

        var spec = {vertexArr:vertices, normalArr: normals, faceArr:faces,color:colors};
        viewer.addCustom(spec);
    }
    triangle(viewer);
    viewer.render();
     */
    public addCustom(spec:CustomShapeSpec) {
        spec = spec || {};
        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        s.addCustom(spec);
        this.shapes.push(s);
        s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

        return s;
    };

    /**
     * Construct isosurface from volumetric data in gaussian cube format
     * @param {String} data - Input file contents
     * @param {String} format - Input file format
     * @param {VolumetricRendererSpec|IsoSurfaceSpec} spec - Shape style specification
     * @return {GLShape}
     *
     * @example


    $3Dmol.get('data/bohr.cube', function(data) {

    viewer.addVolumetricData(data, "cube", {isoval: -0.01, color: "red", opacity: 0.95});
    viewer.setStyle({cartoon:{},stick:{}});
    viewer.zoomTo();
    viewer.render();
    });


     */
    public addVolumetricData(data, format:string, spec:VolumetricRendererSpec|IsoSurfaceSpec={}) {

        var voldata = new VolumeData(data, format);
        if (spec.hasOwnProperty('transferfn')) { //volumetric rendering
            return this.addVolumetricRender(voldata, spec as VolumetricRendererSpec);
        } else {
            return this.addIsosurface(voldata, spec as IsoSurfaceSpec);
        }
    };

    /**
     * Construct isosurface from volumetric data.  This is more flexible
    * than addVolumetricData, but can not be used with py3Dmol.
     * @param {VolumeData} data - volumetric data
     * @param {IsoSurfaceSpec} spec - Shape style specification
     * @return {GLShape}
     *
     @example
     $3Dmol.get('../test_structs/benzene-homo.cube', function(data){
              var voldata = new $3Dmol.VolumeData(data, "cube");
              viewer.addIsosurface(voldata, {isoval: 0.01,
                                             color: "blue"});
              viewer.addIsosurface(voldata, {isoval: -0.01,
                                             color: "red"});
              viewer.zoomTo();
              viewer.render();
            });
     */
    public addIsosurface(data, spec:IsoSurfaceSpec={}, callback?) {
        var s = new GLShape(spec);
        s.shapePosition = this.shapes.length;
        s.addIsosurface(data, spec, callback);
        this.shapes.push(s);
        return s;
    };

    /**
     * Create volumetric renderer for volumetricData
     * @param {VolumeData} data - volumetric data
     * @param {VolumetricRenderSpec} spec - specification of volumetric render
     *
     * @return {GLShape}
     *
     */
    public addVolumetricRender(data, spec:VolumetricRendererSpec) {
        spec = spec || {};
        var s = new GLVolumetricRender(data, spec);
        s.shapePosition = this.shapes.length;
        this.shapes.push(s);
        return s;
    };

    /**
     * Return true if volumetric rendering is supported (WebGL 2.0 required)
     *
     * @return {boolean}
     */
    public hasVolumetricRender() {
        return this.renderer.supportsVolumetric();
    };

    /**
     * Enable/disable fog for content far from the camera
     *
     * @param {boolean} fog whether to enable or disable the fog
     */
    public enableFog(fog:boolean) {
        if (fog) {
            this.scene.fog = new Fog(this.bgColor, 100, 200);
        } else {
            this.config.disableFog = true;
            this.show();
        }
    };

    /**
     * Sets the atomlists of all models in the viewer to specified frame.
     * Shapes and labels can also be displayed by frame.
     * Sets to last frame if framenum out of range
     *
     * @param {number} framenum - fame index to use, starts at zero
     * @return {Promise}
     */
    public setFrame(framenum:number) {
        this.viewer_frame = framenum;
        let viewer = this;
        return new Promise<void>(function (resolve) {
            var modelMap = viewer.models.map(function (model) {
                return model.setFrame(framenum, viewer);
            });
            Promise.all(modelMap)
                .then(function () { resolve(); });
        });
    };

    /**
     * Gets the current viewer frame.
     *
     */
    public getFrame() {
        return this.viewer_frame;
    };

    /**
     * Returns the number of frames that the model with the most frames in the viewer has
     *
     * @return {number}
     */
    public getNumFrames() {
        var mostFrames = 0;
        for (let i = 0; i < this.models.length; i++) {
            if (this.models[i].getNumFrames() > mostFrames) {
                mostFrames = this.models[i].getNumFrames();
            }
        }
        for (let i = 0; i < this.shapes.length; i++) {
            if (this.shapes[i].frame && this.shapes[i].frame >= mostFrames) {
                mostFrames = this.shapes[i].frame + 1;
            }
        }
        for (let i = 0; i < this.labels.length; i++) {
            if (this.labels[i].frame && this.labels[i].frame >= mostFrames) {
                mostFrames = this.labels[i].frame + 1;
            }
        }
        return mostFrames;
    };


    /**
     * Animate all models in viewer from their respective frames
     * @param {Object} options - can specify interval (speed of animation), loop (direction
     * of looping, 'backward', 'forward' or 'backAndForth'), step interval between frames ('step'), startFrame, and reps (numer of repetitions, 0 indicates infinite loop)
     *
     */

    public animate(options) {
        this.incAnim();
        var interval = 100;
        var loop = "forward";
        var reps = 0;
        options = options || {};
        if (options.interval) {
            interval = options.interval;
        }
        if (options.loop) {
            loop = options.loop;
        }
        if (options.reps) {
            reps = options.reps;
        }
        var mostFrames = this.getNumFrames();
        var self = this;
        var currFrame = 0;
        if (options.startFrame) {
            currFrame = options.startFrame % mostFrames;
        }
        var inc = 1;
        if (options.step) {
            inc = options.step;
            reps /= inc;
        }
        var displayCount = 0;
        var displayMax = mostFrames * reps;
        var time = new Date();
        var resolve, timer;
        var display = function (direction) {
            time = new Date();
            if (direction == "forward") {
                self.setFrame(currFrame)
                    .then(function () {
                        currFrame = (currFrame + inc) % mostFrames;
                        resolve();
                    });
            }
            else if (direction == "backward") {
                self.setFrame((mostFrames - 1) - currFrame)
                    .then(function () {
                        currFrame = (currFrame + inc) % mostFrames;
                        resolve();
                    });
            }
            else { //back and forth
                self.setFrame(currFrame)
                    .then(function () {
                        currFrame += inc;
                        inc *= (((currFrame % (mostFrames - 1)) == 0) ? -1 : 1);
                        resolve();
                    });
            }
        };

        resolve = function () {
            self.render();
            if (!self.getCanvas().isConnected) {
                //we no longer exist as part of the DOM
                self.stopAnimate();
            }
            else if (++displayCount == displayMax || !self.isAnimated()) {
                timer.cancel();
                self.animationTimers.delete(timer);
                self.decAnim();
            }
            else {
                var newInterval = interval - (new Date().getTime() - time.getTime());
                newInterval = (newInterval > 0) ? newInterval : 0;
                self.animationTimers.delete(timer);
                timer = new PausableTimer(display, newInterval, loop);
                self.animationTimers.add(timer);
            }
        };

        timer = new PausableTimer(display, 0, loop);
        this.animationTimers.add(timer);
        return this;
    };

    /**
     * Stop animation of all models in viewer
     */
    public stopAnimate() {
        this.animated = 0;
        this.animationTimers.forEach(function (timer: PausableTimer) { timer.cancel(); });
        this.animationTimers = new Set();
        return this;
    };

    /**
     * Pause animation of all models in viewer
     */
    public pauseAnimate() {
        this.animationTimers.forEach(function (timer) { timer.pause(); });
        return this;
    };

    /**
     * Resume animation of all models in viewer
     */
    public resumeAnimate() {
        this.animationTimers.forEach(function (timer) { timer.resume(); });
        return this;
    };


    /**
     * Return true if viewer is currently being animated, false otherwise
     * @return {boolean}
     */
    public isAnimated() {
        return this.animated > 0;
    };


    /**
     * Create and add model to viewer, given molecular data and its format
     *
     * @param {string} data - Input data
     * @param {string} format - Input format ('pdb', 'sdf', 'xyz', 'pqr', or 'mol2')
     * @param {ParserOptionsSpec} options - format dependent options. Attributes depend on the input file format.
     * @example


          viewer.setViewStyle({style:"outline"});
          $3Dmol.get('data/1fas.pqr', function(data){
              viewer.addModel(data, "pqr");
              $3Dmol.get("data/1fas.cube",function(volumedata){
                  viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});

              viewer.render();
              });
              viewer.zoomTo();
          });
     *
     * @return {GLModel}
     */
    public addModel(data?, format="", options?) {
        if (options && !options.defaultcolors) {
            options.defaultcolors = this.defaultcolors;
            options.cartoonQuality = options.cartoonQuality || this.config.cartoonQuality;
        } else if (typeof (options) === 'undefined') {
            options = { defaultcolors: this.defaultcolors, cartoonQuality: this.config.cartoonQuality };
        }
        var m = new GLModel(this.models.length, options);
        m.addMolData(data, format, options);
        this.models.push(m);

        return m;
    };

    /**
     * Given multimodel file and its format, add atom data to the viewer as separate models
     * and return list of these models
     *
     * @param {string} data - Input data
     * @param {string} format - Input format (see {@link FileFormats})
     * @return {Array<GLModel>}
     */
    public addModels(data, format:string, options?) {
        options = options || {};
        options.multimodel = true;
        options.frames = true;

        var modelatoms = GLModel.parseMolData(data, format, options);

        for (var i = 0; i < modelatoms.length; i++) {
            var newModel = new GLModel(this.models.length, this.defaultcolors);
            newModel.setAtomDefaults(modelatoms[i]);
            newModel.addFrame(modelatoms[i]);
            newModel.setFrame(0);
            if (modelatoms.modelData)
                newModel.setModelData(modelatoms.modelData[i]);
            newModel.setDontDuplicateAtoms(!options.duplicateAssemblyAtoms);
            this.models.push(newModel);
        }

        return this.models;
    };

    /**
     * Create and add model to viewer. Given multimodel file and its format,
     * different atomlists are stored in model's frame
     * property and model's atoms are set to the 0th frame
     *
     * @param {string} data - Input data
     * @param {string} format - Input format (see {@link FileFormats})
     * @return {GLModel}
     *
     * @example
            $3Dmol.get('../test_structs/multiple2.xyz', function(data){
              viewer.addModelsAsFrames(data, "xyz");
              viewer.animate({loop: "forward",reps: 1});
              viewer.setStyle({stick:{colorscheme:'magentaCarbon'}});
              viewer.zoomTo();
              viewer.render();
          });
     */
    public addModelsAsFrames(data, format:string, options?) {
        options = options || {};
        options.multimodel = true;
        options.frames = true;
        var m = new GLModel(this.models.length, this.defaultcolors);
        m.addMolData(data, format, options);
        this.models.push(m);

        return m;
    };

    /**
     * Create and add model to viewer. Given multimodel file and its format,
     * all atoms are added to one model
     *
     * @param {string} data - Input data
     * @param {string} format - Input format (see {@link FileFormats})
     * @return {GLModel}
     @example


          $3Dmol.get('../test_structs/multiple.sdf', function(data){
              viewer.addAsOneMolecule(data, "sdf");
              viewer.zoomTo();
              viewer.render();
          });
     */
    public addAsOneMolecule(data, format:string, options?) {
        options = options || {};
        options.multimodel = true;
        options.onemol = true;
        var m = new GLModel(this.models.length, this.defaultcolors);
        m.addMolData(data, format, options);
        this.models.push(m);

        return m;
    };


    /**
     * Delete specified model from viewer
     *
     * @param {GLModel|number} model
     */
    public removeModel(model?:GLModel|number) {
        model = this.getModel(model);
        if (!model)
            return;
        model.removegl(this.modelGroup);
        delete this.models[model.getID()];
        // clear off back of model array
        while (this.models.length > 0
            && typeof (this.models[this.models.length - 1]) === "undefined")
            this.models.pop();
        return this;
    };

    /**
     * Delete all existing models
     */
    public removeAllModels() {
        for (var i = 0; i < this.models.length; i++) {
            var model = this.models[i];
            if (model) model.removegl(this.modelGroup);

        }
        this.models.splice(0, this.models.length); //don't simply overwrite array in case linked
        return this;
    };

    /**
     * Export one or all of the loaded models into ChemDoodle compatible JSON.
     * @param {boolean} includeStyles - Whether or not to include style information.
     * @param {number} modelID - Optional parameter for which model to export. If left out, export all of them.
     * @return {string}
     */
    public exportJSON(includeStyles:boolean, modelID:number) {
        var object: any = {};
        if (modelID === undefined) {
            object.m = this.models.map(function (model) {
                return model.toCDObject(includeStyles);
            });
        } else {
            object.m = [this.models[modelID].toCDObject()];
        }
        return JSON.stringify(object);
    };

    /** return a VRML string representation of the scene.  Include VRML header information
     * @return VRML
     */
    public exportVRML() {
        var savedmodelGroup = this.modelGroup;
        this.applyToModels("removegl", this.modelGroup); //cleanup
        this.modelGroup = new Object3D();
        //rendering with plain mesh
        this.render(null, { supportsImposters: false, supportsAIA: false, regen: true });
        var ret = '#VRML V2.0 utf8\n' + this.modelGroup.vrml() + '\n';
        this.applyToModels("removegl", this.modelGroup); //cleanup
        this.modelGroup = savedmodelGroup;
        return ret;
    };

    /**
     * Create a new model from atoms specified by sel.
     * If extract, removes selected atoms from existing models
     *
     * @param {AtomSelectionSpec} sel - Atom selection specification
     * @param {boolean=} extract - If true, remove selected atoms from existing models
     * @return {GLModel}
     */
    public createModelFrom(sel:AtomSelectionSpec, extract:boolean=false) {
        var m = new GLModel(this.models.length, this.defaultcolors);
        for (var i = 0; i < this.models.length; i++) {
            if (this.models[i]) {
                var atoms = this.models[i].selectedAtoms(sel);
                m.addAtoms(atoms);
                if (extract)
                    this.models[i].removeAtoms(atoms);
            }
        }
        this.models.push(m);
        return m;
    };

    private applyToModels(func:string, sel:any, value1?, value2?, value3?, value4?, value5?) {

        //apply func to all models that are selected by sel with value1 and 2
        //sel might not be a selection, in which case getModelList returns everything
        var ms = this.getModelList(sel);
        for (var i = 0; i < ms.length; i++) {
            ms[i][func](sel, value1, value2, value3, value4, value5);
        }
    }

    /**
     * Set style properties to all selected atoms
     *
     * @param {AtomSelectionSpec} sel - Atom selection specification.  Can be omitted to select all.
     * @param {AtomStyleSpec} style - Style spec to apply to specified atoms
     *
     * @example
     viewer.setBackgroundColor(0xffffffff);
    $3Dmol.download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
    m.setStyle({chain:'A'},{'cartoon':{color:'spectrum'}});
    m.setStyle({chain:'C'},{'cartoon':{style:'trace',color:'blue'}});
    m.setStyle({chain:'E'},{'cartoon':{tubes:true,arrows:true,color:'green',opacity:0.75}});
    m.setStyle({chain:'B'},{'cartoon':{color:'red',opacity:0.5}});
    m.setStyle({chain:'D'},{'cartoon':{style:'trace',color:'grey',opacity:0.75}});
    m.setStyle({chain:'F'},{'cartoon':{arrows:true,color:'white'}});
    // viewer.addStyle({chain:'B'},{line:{}});
    viewer.zoomTo();
    viewer.render();
    });
     */
    public setStyle(sel:AtomSelectionSpec, style:AtomStyleSpec);
    public setStyle(sel:AtomStyleSpec);
    public setStyle(sel:unknown, style?:unknown) {
        if (typeof (style) === 'undefined') {
            //if a single argument is provided, assume it is a style and select all
            style = sel as AtomStyleSpec;
            sel = {};
        }

        this.applyToModels("setStyle", sel, style, false);
        return this;
    };

    /**
     * Add style properties to all selected atoms
     *
     * @param {AtomSelectionSpec} sel - Atom selection specification.  Can be omitted to select all
     * @param {AtomStyleSpec} style - style spec to add to specified atoms
     @example

    $3Dmol.download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
    viewer.setStyle({cartoon:{}});
    //keep cartoon style, but show thick sticks for chain A
    viewer.addStyle({chain:'A'},{stick:{radius:.5,colorscheme:"magentaCarbon"}});
    viewer.zoomTo();
    viewer.render();
    });
     */
    public addStyle(sel:AtomSelectionSpec, style:AtomStyleSpec);
    public addStyle(sel:AtomStyleSpec);
    public addStyle(sel:unknown, style?:unknown) {
        if (typeof (style) === 'undefined') {
            //if a single argument is provided, assume it is a style and select all
            style = sel;
            sel = {};
        }
        this.applyToModels("setStyle", sel, style, true);
        return this;
    };


    /**
     * Set click-handling properties to all selected atoms. *Important*: render must be called for this to take effect.
     *
     * @param {AtomSelectionSpec} sel - atom selection to apply clickable settings to
     * @param {boolean} clickable - whether click-handling is enabled for the selection
     * @param {function} callback - function called when an atom in the selection is clicked. The function is passed
     * the selected (foremost) object, the viewer, the triggering event, the associated container, and a list
     * of all intersecting objects with their distances from the viewer.
     *
     * @example
        $3Dmol.download("cid:307900",viewer,{},function(){

               viewer.setStyle({},{sphere:{}});
               viewer.setClickable({},true,function(atom,viewer,event,container) {
                   viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8});
               });
               viewer.render();
    });
     */
    public setClickable(sel:AtomSelectionSpec, clickable:boolean, callback) {
        this.applyToModels("setClickable", sel, clickable, callback);
        return this;
    };
    /** Set hoverable and callback of selected atoms
     *
     * @param {AtomSelectionSpec} sel - atom selection to apply hoverable settings to
     * @param {boolean} hoverable - whether hover-handling is enabled for the selection
     * @param {function} hover_callback - function called when an atom in the selection is hovered over.  The function has the same signature as a click handler.
     * @param {function} unhover_callback - function called when the mouse moves out of the hover area
    @example
    $3Dmol.download("pdb:1ubq",viewer,{},function(){

           viewer.setHoverable({},true,function(atom,viewer,event,container) {
               if(!atom.label) {
                atom.label = viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
               }
           },
           function(atom) {
               if(atom.label) {
                viewer.removeLabel(atom.label);
                delete atom.label;
               }
            }
           );
           viewer.setStyle({},{stick:{}});
           viewer.render();
    });

     */
    public setHoverable(sel:AtomSelectionSpec, hoverable:boolean, hover_callback, unhover_callback) {
        this.applyToModels("setHoverable", sel, hoverable, hover_callback, unhover_callback);
        return this;
    };

    /** enable context menu and callback of selected atoms
     *
     * @param {AtomSelectionSpec} sel - atom selection to apply hoverable settings to
     * @param {boolean} contextMenuEnabled - whether contextMenu-handling is enabled for the selection

     */
    public enableContextMenu(sel:AtomSelectionSpec, contextMenuEnabled:boolean) {
        this.applyToModels("enableContextMenu", sel, contextMenuEnabled);
        return this;
    };

    /**
     * If  atoms have dx, dy, dz properties (in some xyz files), vibrate populates each model's frame property based on parameters.
     * Models can then be animated
     *
     * @param {number} numFrames - number of frames to be created, default to 10
     * @param {number} amplitude - amplitude of distortion, default to 1 (full)
     * @param {boolean} bothWays - if true, extend both in positive and negative directions by numFrames
     * @param {ArrowSpec} arrowSpec - specification for drawing animated arrows. If color isn't specified, atom color (sphere, stick, line preference) is used.
     */
    public vibrate(numFrames:number, amplitude:number, bothways:boolean, arrowSpec:ArrowSpec) {
        this.applyToModels("vibrate", numFrames, amplitude, bothways, this, arrowSpec);
        return this;
    };

    /**
     * @param {AtomSelectionSpec} sel
     * @param {string} prop
     * @param {Gradient|string} scheme
     * @param {object} range
     */
    public setColorByProperty(sel:AtomSelectionSpec, prop:string, scheme:Gradient|string, range) {
        this.applyToModels("setColorByProperty", sel, prop, scheme, range);
        return this;
    };

    /**
     * @param {AtomSelectionSpec} sel
     * @param {object} colors
     */
    public setColorByElement(sel:AtomSelectionSpec, colors) {
        this.applyToModels("setColorByElement", sel, colors);
        return this;
    };

    /**
     *
     * @param {AtomSpec[]} atomlist
     * @param {Array}
     *            extent
     * @return {Array}
     */
    private static getAtomsWithin(atomlist:AtomSpec[], extent) {
        var ret = [];

        for (let i = 0; i < atomlist.length; i++) {
            var atom = atomlist[i];
            if (typeof (atom) == "undefined")
                continue;

            if (atom.x < extent[0][0] || atom.x > extent[1][0])
                continue;
            if (atom.y < extent[0][1] || atom.y > extent[1][1])
                continue;
            if (atom.z < extent[0][2] || atom.z > extent[1][2])
                continue;
            ret.push(atom);
        }
        return ret;
    };

    // return volume of extent
    private static volume(extent) {
        var w = extent[1][0] - extent[0][0];
        var h = extent[1][1] - extent[0][1];
        var d = extent[1][2] - extent[0][2];
        return w * h * d;
    }; // volume
    /*
     * Break up bounding box/atoms into smaller pieces so we can parallelize
     * with webworkers and also limit the size of the working memory Returns
     * a list of bounding boxes with the corresponding atoms. These extents
     * are expanded by 4 angstroms on each side.
     */
    /**
     *
     * @param {Array}
     *            extent
     * @param {AtomSpec[]} atomlist
     * @param {AtomSpec[]} atomstoshow
     * @return {Array}
     */
    private carveUpExtent(extent, atomlist: AtomSpec[], atomstoshow: AtomSpec[]) {
        let ret = [];

        let index2atomlist = {}; //map from atom.index to position in atomlist
        for (let i = 0, n = atomlist.length; i < n; i++) {
            index2atomlist[atomlist[i].index] = i;
        }

        let atomsToListIndex = function (atoms) {
            //return a list of indices into atomlist
            let ret = [];
            for (let i = 0, n = atoms.length; i < n; i++) {
                if (atoms[i].index in index2atomlist)
                    ret.push(index2atomlist[atoms[i].index]);
            }
            return ret;
        };
        let copyExtent = function (extent) {
            // copy just the dimensions
            let ret = [];
            ret[0] = [extent[0][0], extent[0][1], extent[0][2]];
            ret[1] = [extent[1][0], extent[1][1], extent[1][2]];
            return ret;
        }; // copyExtent
        let splitExtentR = function (extent) {
            // recursively split until volume is below maxVol
            if (GLViewer.volume(extent) < GLViewer.maxVolume) {
                return [extent];
            } else {
                // find longest edge
                var w = extent[1][0] - extent[0][0];
                var h = extent[1][1] - extent[0][1];
                var d = extent[1][2] - extent[0][2];

                var index;

                if (w > h && w > d) {
                    index = 0;
                } else if (h > w && h > d) {
                    index = 1;
                } else {
                    index = 2;
                }

                // create two halves, splitting at index
                var a = copyExtent(extent);
                var b = copyExtent(extent);
                var mid = (extent[1][index] - extent[0][index]) / 2
                    + extent[0][index];
                a[1][index] = mid;
                b[0][index] = mid;

                var alist = splitExtentR(a);
                var blist = splitExtentR(b);
                return alist.concat(blist);
            }
        }; // splitExtentR

        // divide up extent
        let splits = splitExtentR(extent);
        // now compute atoms within expanded (this could be more efficient)
        let off = 6; // enough for water and 2*r, also depends on scale
        // factor
        for (let i = 0, n = splits.length; i < n; i++) {
            let e = copyExtent(splits[i]);
            e[0][0] -= off;
            e[0][1] -= off;
            e[0][2] -= off;
            e[1][0] += off;
            e[1][1] += off;
            e[1][2] += off;

            let atoms = GLViewer.getAtomsWithin(atomlist, e);
            let toshow = GLViewer.getAtomsWithin(atomstoshow, splits[i]);

            // ultimately, divide up by atom for best meshing
            ret.push({
                extent: splits[i],
                atoms: atomsToListIndex(atoms),
                toshow: atomsToListIndex(toshow)
            });
        }

        return ret;
    };

    // create a mesh defined from the passed vertices and faces and material
    // Just create a single geometry chunk - broken up whether sync or not
    /**
     *
     * @param {AtomSpec[]} atoms
     * @param {{vertices:number,faces:number}}
     *            VandF
     * @param {MeshLambertMaterial}
     *            mat
     * @return {Mesh}
     */
    private static generateSurfaceMesh(atoms: AtomSpec[], VandF, mat: MeshLambertMaterial) {
        var geo = new Geometry(true);
        // Only one group per call to generate surface mesh (addSurface
        // should split up mesh render)
        var geoGroup = geo.updateGeoGroup(0);

        // set colors for vertices
        var colors = [];
        for (let i = 0, il = atoms.length; i < il; i++) {
            var atom = atoms[i];
            if (atom) {
                if (typeof (atom.surfaceColor) != "undefined") {
                    colors[i] = atom.surfaceColor;
                } else if (atom.color) // map from atom
                    colors[i] = CC.color(atom.color);
            }
        }

        var vertexArray = geoGroup.vertexArray;
        // reconstruct vertices and faces
        var v = VandF.vertices;
        for (let i = 0, il = v.length; i < il; i++) {
            let offset = geoGroup.vertices * 3;
            vertexArray[offset] = v[i].x;
            vertexArray[offset + 1] = v[i].y;
            vertexArray[offset + 2] = v[i].z;
            geoGroup.vertices++;
        }

        //set colorArray of there are per-atom colors
        var colorArray = geoGroup.colorArray;

        if (mat.voldata && mat.volscheme) {
            //convert volumetric data into colors
            var scheme = mat.volscheme;
            var voldata = mat.voldata;
            var range = scheme.range() || [-1, 1];
            for (let i = 0, il = v.length; i < il; i++) {
                let val = voldata.getVal(v[i].x, v[i].y, v[i].z);
                let col = CC.color(scheme.valueToHex(val, range));
                let offset = i * 3;
                colorArray[offset] = col.r;
                colorArray[offset + 1] = col.g;
                colorArray[offset + 2] = col.b;
            }
        }
        else if (colors.length > 0) { //have atom colors
            for (let i = 0, il = v.length; i < il; i++) {
                let A = v[i].atomid;
                let offsetA = i * 3;

                colorArray[offsetA] = colors[A].r;
                colorArray[offsetA + 1] = colors[A].g;
                colorArray[offsetA + 2] = colors[A].b;
            }
        }

        var faces = VandF.faces;
        geoGroup.faceidx = faces.length;// *3;
        geo.initTypedArrays();

        var verts = geoGroup.vertexArray;
        var normalArray = geoGroup.normalArray;
        var vA, vB, vC, norm;

        // Setup colors, faces, and normals
        for (let i = 0, il = faces.length; i < il; i += 3) {

            // var a = faces[i].a, b = faces[i].b, c = faces[i].c;
            var a = faces[i], b = faces[i + 1], c = faces[i + 2];
            var offsetA = a * 3, offsetB = b * 3, offsetC = c * 3;

            // setup Normals
            // todo - calculate normals in parallel code
            vA = new Vector3(verts[offsetA], verts[offsetA + 1],
                verts[offsetA + 2]);
            vB = new Vector3(verts[offsetB], verts[offsetB + 1],
                verts[offsetB + 2]);
            vC = new Vector3(verts[offsetC], verts[offsetC + 1],
                verts[offsetC + 2]);

            vC.subVectors(vC, vB);
            vA.subVectors(vA, vB);
            vC.cross(vA);

            // face normal
            norm = vC;
            norm.normalize();

            normalArray[offsetA] += norm.x;
            normalArray[offsetB] += norm.x;
            normalArray[offsetC] += norm.x;
            normalArray[offsetA + 1] += norm.y;
            normalArray[offsetB + 1] += norm.y;
            normalArray[offsetC + 1] += norm.y;
            normalArray[offsetA + 2] += norm.z;
            normalArray[offsetB + 2] += norm.z;
            normalArray[offsetC + 2] += norm.z;

        }
        geoGroup.faceArray = new Uint16Array(faces);
        var mesh = new Mesh(geo, mat as Material);
        return mesh;
    };

    // do same thing as worker in main thread
    /**
     *
     * @param {SurfaceType}
     *            type
     * @param {Array}
     *            expandedExtent
     * @param {AtomSpec[]}
     *            extendedAtoms
     * @param {AtomSpec[]}
     *            atomsToShow
     * @param {AtomSpec[]} atoms
     * @param {number}
     *            vol
     * @return {Object}
     */
    private static generateMeshSyncHelper(type: SurfaceType, expandedExtent,
        extendedAtoms: AtomSpec[], atomsToShow: AtomSpec[], atoms: AtomSpec[], vol: number) {
        //            var time = new Date();
        var ps = new ProteinSurface();
        ps.initparm(expandedExtent, (type === 1) ? false : true, vol);

        //            var time2 = new Date();
        //console.log("initialize " + (time2 - time) + "ms");

        ps.fillvoxels(atoms, extendedAtoms);

        //            var time3 = new Date();
        //console.log("fillvoxels " + (time3 - time2) + "  " + (time3 - time) + "ms");

        ps.buildboundary();

        if (type == SurfaceType.SES || type == SurfaceType.MS) {
            ps.fastdistancemap();
            ps.boundingatom(false);
            ps.fillvoxelswaals(atoms, extendedAtoms);
        }

        //            var time4 = new Date();
        //console.log("buildboundaryetc " + (time4 - time3) + "  " + (time4 - time) + "ms");

        ps.marchingcube(type);

        //            var time5 = new Date();
        //console.log("marching cube " + (time5 - time4) + "  "+ (time5 - time) + "ms");

        return ps.getFacesAndVertices(atomsToShow);
    };

    /*
     *
     * @param {SurfaceStyleSpec}
     *            style
     * @return {MeshLambertMaterial}
     */
    private static getMatWithStyle(style: SurfaceStyleSpec) {
        var mat = new MeshLambertMaterial();
        mat.vertexColors = Coloring.VertexColors;

        for (var prop in style) {
            if (prop === "color" || prop === "map") {
                // ignore
            } else if (style.hasOwnProperty(prop))
                mat[prop] = style[prop];
        }
        if (style.opacity !== undefined) {
            if (style.opacity === 1)
                mat.transparent = false;
            else
                mat.transparent = true;
        }

        return mat;
    }


    /**
     * Adds an explicit mesh as a surface object.
     * @param {Mesh}
     *            mesh
     * @param {Object}
     *            style
     * @returns {number} surfid
     */
    public addMesh(mesh: Mesh) {
        var surfobj = {
            geo: mesh.geometry,
            mat: mesh.material,
            done: true,
            finished: false //the rendered finishes surfaces when they are done
        };
        var surfid = this.nextSurfID();
        this.surfaces[surfid] = surfobj;
        return surfid;
    };

    //return a shallow copy of list l, e.g., for atoms so we can
    //ignore superficial changes (ie surfacecolor, position) that happen
    //while we're surface building
    private static shallowCopy(l) {
        var ret = [];
        let length = l.length;
        for (let i = 0; i < length; i++) {
            ret[i] = extend({}, l[i]);
        }
        return ret;
    };



    /**
     * Add surface representation to atoms
     * @param {SurfaceType|string} type - Surface type (VDW, MS, SAS, or SES)
     * @param {SurfaceStyleSpec} style - optional style specification for surface material (e.g. for different coloring scheme, etc)
     * @param {AtomSelectionSpec} atomsel - Show surface for atoms in this selection
     * @param {AtomSelectionSpec} allsel - Use atoms in this selection to calculate surface; may be larger group than 'atomsel'
     * @param {AtomSelectionSpec} focus - Optionally begin rendering surface specified atoms
     * @param {function} surfacecallback - function to be called after setting the surface
     * @return {Promise} promise - Returns a promise that ultimately resovles to the surfid.  Returns surfid immediately if surfacecallback is specified.  Returned promise has a [surfid, GLViewer, style, atomsel, allsel, focus] fields for immediate access.
     */
    public addSurface(stype:SurfaceType|string, style:SurfaceStyleSpec={}, atomsel:AtomSelectionSpec={},
        allsel?:AtomSelectionSpec, focus?: AtomSelectionSpec, surfacecallback?) {
        // type 1: VDW 3: SAS 4: MS 2: SES
        // if sync is true, does all work in main thread, otherwise uses
        // workers
        // with workers, must ensure group is the actual modelgroup since
        // surface
        // will get added asynchronously
        // all atoms in atomlist are used to compute surfaces, but only the
        // surfaces
        // of atomsToShow are displayed (e.g., for showing cavities)
        // if focusSele is specified, will start rending surface around the

        //surfacecallback gets called when done
        let surfid = this.nextSurfID();
        let mat = null;
        let self = this;
        let type: SurfaceType | 0 = SurfaceType.VDW;

        if (typeof stype == "string") {
            if (GLViewer.surfaceTypeMap[stype.toUpperCase()] !== undefined)
                type = GLViewer.surfaceTypeMap[stype];
            else {
                console.log("Surface type : " + stype + " is not recognized");
            }
        } else if(typeof stype == "number") {
            type = stype;
        }

        // atoms specified by this selection
        var atomlist = null, focusSele = null;
        //TODO: currently generating a shallow copy to avoid problems when atoms are chagned
        //during surface generation - come up with a better solution
        var atomsToShow = GLViewer.shallowCopy(this.getAtomsFromSel(atomsel));
        if (!allsel) {
            atomlist = atomsToShow;
        }
        else {
            atomlist = GLViewer.shallowCopy(this.getAtomsFromSel(allsel));
        }

        adjustVolumeStyle(style);
        var symmetries = false;
        var n;
        for (n = 0; n < this.models.length; n++) {
            if (this.models[n]) {
                var symMatrices = this.models[n].getSymmetries();
                if (symMatrices.length > 1 || (symMatrices.length == 1 && !(symMatrices[0].isIdentity()))) {
                    symmetries = true;
                    break;
                }
            }
        }

        var addSurfaceHelper = function addSurfaceHelper(surfobj, atomlist: AtomSpec[], atomsToShow: AtomSpec[]) {
            //function returns promise with surfid resolved
            if (!focus) {
                focusSele = atomsToShow;
            } else {
                focusSele = GLViewer.shallowCopy(self.getAtomsFromSel(focus));
            }

            var atom;
            //                var time = new Date();
            var extent = getExtent(atomsToShow, true);
            if (style.map && style.map.prop) {
                // map color space using already set atom properties
                var prop = style.map.prop;
                let scheme = getGradient(style.map.scheme || style.map.gradient || new Gradient.RWB());
                let range = scheme.range();
                if (!range) {
                    range = getPropertyRange(atomsToShow, prop);
                }
                style.colorscheme = { prop: prop as string, gradient: scheme };

            }

            //cache surface color on each atom
            for (let i = 0, il = atomlist.length; i < il; i++) {
                atom = atomlist[i];
                atom.surfaceColor = getColorFromStyle(atom, style);
            }

            var totalVol = GLViewer.volume(extent); // used to scale resolution
            var extents = self.carveUpExtent(extent, atomlist, atomsToShow);

            if (focusSele && focusSele.length && focusSele.length > 0) {
                var seleExtent = getExtent(focusSele, true);
                // sort by how close to center of seleExtent
                var sortFunc = function (a, b) {
                    var distSq = function (ex, sele) {
                        // distance from e (which has no center of mass) and
                        // sele which does
                        var e = ex.extent;
                        var x = e[1][0] - e[0][0];
                        var y = e[1][1] - e[0][1];
                        var z = e[1][2] - e[0][2];
                        var dx = (x - sele[2][0]);
                        dx *= dx;
                        var dy = (y - sele[2][1]);
                        dy *= dy;
                        var dz = (z - sele[2][2]);
                        dz *= dz;

                        return dx + dy + dz;
                    };
                    var d1 = distSq(a, seleExtent);
                    var d2 = distSq(b, seleExtent);
                    return d1 - d2;
                };
                extents.sort(sortFunc);
            }

            var reducedAtoms = [];
            // to reduce amount data transfered, just pass x,y,z,serial and elem
            for (let i = 0, il = atomlist.length; i < il; i++) {
                atom = atomlist[i];
                reducedAtoms[i] = {
                    x: atom.x,
                    y: atom.y,
                    z: atom.z,
                    serial: i,
                    elem: atom.elem
                };
            }

            var sync = !!(syncSurface);
            if (sync) { // don't use worker, still break up for memory purposes

                // to keep the browser from locking up, call through setTimeout
                var callSyncHelper = function callSyncHelper(i) {
                    return new Promise<void>(function (resolve) {
                        var VandF = GLViewer.generateMeshSyncHelper(type as SurfaceType, extents[i].extent,
                            extents[i].atoms, extents[i].toshow, reducedAtoms,
                            totalVol);
                        //complicated surfaces sometimes have > 2^16 vertices
                        var VandFs = splitMesh({ vertexArr: VandF.vertices, faceArr: VandF.faces });
                        for (var vi = 0, vl = VandFs.length; vi < vl; vi++) {
                            VandF = {
                                vertices: VandFs[vi].vertexArr,
                                faces: VandFs[vi].faceArr
                            };
                            var mesh = GLViewer.generateSurfaceMesh(atomlist, VandF, mat);
                            mergeGeos(surfobj.geo, mesh);
                        }
                        self.render();
                        resolve();
                    });
                };
                var promises = [];
                for (let i = 0; i < extents.length; i++) {
                    promises.push(callSyncHelper(i));
                }
                return Promise.all(promises)
                    .then(function () {
                        surfobj.done = true;
                        return Promise.resolve(surfid);
                    });

                // TODO: Asynchronously generate geometryGroups (not separate
                // meshes) and merge them into a single geometry
            } else { // use worker

                var workers = [];
                if (type < 0)
                    type = 0; // negative reserved for atom data
                for (let i = 0, il = GLViewer.numWorkers; i < il; i++) {
                    var w = new Worker($3Dmol.SurfaceWorker);
                    workers.push(w);
                    w.postMessage({
                        'type': -1,
                        'atoms': reducedAtoms,
                        'volume': totalVol
                    });
                }

                return new Promise(function (resolve, reject) {
                    var cnt = 0;

                    var releaseMemory = function () {
                        if (!workers || !workers.length) return;
                        workers.forEach(function (worker) {
                            if (worker && worker.terminate) {
                                worker.terminate();
                            }
                        });
                    };

                    var rfunction = function (event) {
                        var VandFs = splitMesh({
                            vertexArr: event.data.vertices,
                            faceArr: event.data.faces
                        });
                        for (var i = 0, vl = VandFs.length; i < vl; i++) {
                            var VandF = {
                                vertices: VandFs[i].vertexArr,
                                faces: VandFs[i].faceArr
                            };
                            var mesh = GLViewer.generateSurfaceMesh(atomlist, VandF, mat);
                            mergeGeos(surfobj.geo, mesh);
                        }
                        self.render();

                        //    console.log("async mesh generation " + (+new Date() - time) + "ms");
                        cnt++;
                        if (cnt == extents.length) {
                            surfobj.done = true;
                            releaseMemory();
                            resolve(surfid); //caller of helper will resolve callback if present
                        }
                    };

                    var efunction = function (event) {
                        releaseMemory();
                        console.log(event.message + " (" + event.filename + ":" + event.lineno + ")");
                        reject(event);
                    };

                    for (let i = 0; i < extents.length; i++) {
                        var worker = workers[i % workers.length];
                        worker.onmessage = rfunction;
                        worker.onerror = efunction;

                        worker.postMessage({
                            'type': type,
                            'expandedExtent': extents[i].extent,
                            'extendedAtoms': extents[i].atoms,
                            'atomsToShow': extents[i].toshow
                        });
                    }
                });
            }
        };

        style = style || {};
        mat = GLViewer.getMatWithStyle(style);
        var surfobj: any = [];
        //save configuration of surface
        surfobj.style = style;
        surfobj.atomsel = atomsel;
        surfobj.allsel = allsel;
        surfobj.focus = focus;
        var promise = null;
        if (symmetries) { //do preprocessing
            var modelsAtomList = {};
            var modelsAtomsToShow = {};
            for (n = 0; n < this.models.length; n++) {
                modelsAtomList[n] = [];
                modelsAtomsToShow[n] = [];
            }
            for (n = 0; n < atomlist.length; n++) {
                modelsAtomList[atomlist[n].model].push(atomlist[n]);
            }
            for (n = 0; n < atomsToShow.length; n++) {
                modelsAtomsToShow[atomsToShow[n].model].push(atomsToShow[n]);
            }
            var promises = [];
            for (n = 0; n < this.models.length; n++) {
                if (modelsAtomsToShow[n].length > 0) {
                    surfobj.push({
                        geo: new Geometry(true),
                        mat: mat,
                        done: false,
                        finished: false,
                        symmetries: this.models[n].getSymmetries()
                        // also webgl initialized
                    });
                    promises.push(addSurfaceHelper(surfobj[surfobj.length - 1], modelsAtomList[n], modelsAtomsToShow[n]));
                }
            }
            promise = Promise.all(promises);
        }
        else {
            surfobj.push({
                geo: new Geometry(true),
                mat: mat,
                done: false,
                finished: false,
                symmetries: [new Matrix4()]
            });
            promise = addSurfaceHelper(surfobj[surfobj.length - 1], atomlist, atomsToShow);
        }
        this.surfaces[surfid] = surfobj;
        promise.surfid = surfid;

        if (surfacecallback && typeof (surfacecallback) == "function") {
            promise.then(function (surfid) {
                surfacecallback(surfid);
            });
            return surfid;
        }
        else {
            return promise;
        }
    };

    /**
     * Set the surface material to something else, must render change
     * @param {number} surf - Surface ID to apply changes to
     * @param {SurfaceStyleSpec} style - new material style specification
     @example
     $3Dmol.get("data/9002806.cif",function(data){
        viewer.addModel(data);
        viewer.setStyle({stick:{}});
        let surf = viewer.addSurface("SAS");
        surf.then(function() {
            viewer.setSurfaceMaterialStyle(surf.surfid, {color:'blue',opacity:0.5});
            viewer.render();
            });
       });
     */
    public setSurfaceMaterialStyle(surf:number, style:SurfaceStyleSpec) {
        adjustVolumeStyle(style);
        if (this.surfaces[surf]) {
            var surfArr = this.surfaces[surf];
            surfArr.style = style;
            for (var i = 0; i < surfArr.length; i++) {
                var mat = surfArr[i].mat = GLViewer.getMatWithStyle(style);
                surfArr[i].mat.side = FrontSide;
                if (style.color) {
                    surfArr[i].mat.color = style.color;
                    surfArr[i].geo.colorsNeedUpdate = true;
                    const c = CC.color(style.color);
                    surfArr[i].geo.setColors(function () { return c; });
                }
                else if (mat.voldata && mat.volscheme) {
                    //convert volumetric data into colors
                    const scheme = mat.volscheme;
                    const voldata = mat.voldata;
                    const cc = CC;
                    const range = scheme.range() || [-1, 1];
                    surfArr[i].geo.setColors(function (x, y, z) {
                        let val = voldata.getVal(x, y, z);
                        let col = cc.color(scheme.valueToHex(val, range));
                        return col;
                    });
                }
                surfArr[i].finished = false; // trigger redraw
            }
        }
        return this;
    };

    /**
     * Return surface object
     * @param {number} surf - surface id
     */
    public getSurface(surf:number) {
        return this.surfaces[surf];
    };

    /**
     * Remove surface with given ID
     * @param {number} surf - surface id
     */
    public removeSurface(surf:number) {
        var surfArr = this.surfaces[surf];
        for (var i = 0; i < surfArr.length; i++) {
            if (surfArr[i] && surfArr[i].lastGL) {
                if (surfArr[i].geo !== undefined)
                    surfArr[i].geo.dispose();
                if (surfArr[i].mat !== undefined)
                    surfArr[i].mat.dispose();
                this.modelGroup.remove(surfArr[i].lastGL); // remove from scene
            }
        }
        delete this.surfaces[surf];
        this.show();
        return this;
    };

    /** Remove all surfaces.
     **/
    public removeAllSurfaces() {
        for (var n in this.surfaces) {
            if (!this.surfaces.hasOwnProperty(n)) continue;
            var surfArr = this.surfaces[n];
            for (var i = 0; i < surfArr.length; i++) {
                if (surfArr[i] && surfArr[i].lastGL) {
                    if (surfArr[i].geo !== undefined)
                        surfArr[i].geo.dispose();
                    if (surfArr[i].mat !== undefined)
                        surfArr[i].mat.dispose();
                    this.modelGroup.remove(surfArr[i].lastGL); // remove from scene
                }
            }
            delete this.surfaces[n];
        }
        this.show();
        return this;
    };

    /** return Jmol moveto command to position this scene */
    public jmolMoveTo() {
        var pos = this.modelGroup.position;
        // center on same position
        var ret = "center { " + (-pos.x) + " " + (-pos.y) + " " + (-pos.z)
            + " }; ";
        // apply rotation
        var q = this.rotationGroup.quaternion;
        ret += "moveto .5 quaternion { " + q.x + " " + q.y + " " + q.z
            + " " + q.w + " };";
        // zoom is tricky.. maybe i would be best to let callee zoom on
        // selection?
        // can either do a bunch of math, or maybe zoom to the center with a
        // fixed
        // but reasonable percentage

        return ret;
    };

    /** Clear scene of all objects
     * */
    public clear() {
        this.removeAllSurfaces();
        this.removeAllModels();
        this.removeAllLabels();
        this.removeAllShapes();
        this.show();
        return this;
    };

    // props is a list of objects that select certain atoms and enumerate
    // properties for those atoms
    /**
     * Add specified properties to all atoms matching input argument
     * @param {Object} props, either array of atom selectors with associated props, or function that takes atom and sets its properties
     * @param {AtomSelectionSpec} sel  - subset of atoms to work on - model selection must be specified here
         @example
         $3Dmol.get('../test_structs/b.sdf', function(data){
                  viewer.addModel(data,'sdf');
                  let props = [];
                  //make the atom index a property x
                  for(let i = 0; i < 8; i++) {
                    props.push({index:i,props:{'x':i}});
                  }
                  viewer.mapAtomProperties(props);
                  viewer.setStyle({sphere:{colorscheme:{gradient:'roygb',prop:'x',min:0,max:8}}});
                  viewer.zoomTo();
                  viewer.render();
                });
     */
    public mapAtomProperties(props, sel:AtomSelectionSpec) {
        sel = sel || {};
        var atoms = this.getAtomsFromSel(sel);

        if (typeof (props) == "function") {
            for (let a = 0, numa = atoms.length; a < numa; a++) {
                let atom = atoms[a];
                props(atom);
            }
        }
        else {
            for (let a = 0, numa = atoms.length; a < numa; a++) {
                var atom = atoms[a];
                for (let i = 0, n = props.length; i < n; i++) {
                    let prop = props[i];
                    if (prop.props) {
                        for (var p in prop.props) {
                            if (prop.props.hasOwnProperty(p)) {
                                // check the atom
                                if (this.atomIsSelected(atom, prop)) {
                                    if (!atom.properties)
                                        atom.properties = {};
                                    atom.properties[p] = prop.props[p];
                                }
                            }
                        }
                    }
                }
            }
        }
        return this;
    };


    /**
     * Synchronize this view matrix of this viewer to the passed viewer.
     * When the viewpoint of this viewer changes, the other viewer will
     * be set to this viewer's view.
     * @param {$3Dmol.GLViewer} otherview
     */
    public linkViewer(otherviewer:GLViewer) {
        this.linkedViewers.push(otherviewer);
        return this;
    };

    /**
     * Return the z distance between the model and the camera
     * @return {number} distance
     */
    public getPerceivedDistance() {
        return this.CAMERA_Z - this.rotationGroup.position.z;
    };

    /**
     * Set the distance between the model and the camera
     * Essentially zooming. Useful while stereo rendering.
     */
    public setPerceivedDistance(dist:number) {
        this.rotationGroup.position.z = this.CAMERA_Z - dist;
    };

    /**
     * Used for setting an approx value of eyeSeparation. Created for calling by StereoViewer object
     * @return {number} camera x position
     */
    public setAutoEyeSeparation(isright:boolean, x:number) {
        var dist = this.getPerceivedDistance();
        if (!x) x = 5.0;
        if (isright || this.camera.position.x > 0) //setting a value of dist*tan(x)
            this.camera.position.x = dist * Math.tan(Math.PI / 180.0 * x);
        else
            this.camera.position.x = -dist * Math.tan(Math.PI / 180.0 * x);
        this.camera.lookAt(new Vector3(0, 0, this.rotationGroup.position.z));
        return this.camera.position.x;
    };

    /**
     * Set the default cartoon quality for newly created models.  Default is 5.
     * Current models are not affected.
     * @number quality, higher results in higher resolution renders
     */
    public setDefaultCartoonQuality(val:number) {
        this.config.cartoonQuality = val;
    };

}


/**
 * Create and initialize an appropriate viewer at supplied HTML element using specification in config
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {ViewerSpec} [config] Viewer configuration
 * @return {GLViewer} GLViewer, null if unable to instantiate WebGL
 * @example
   var viewer = $3Dmol.createViewer(
     'gldiv', //id of div to create canvas in
     {
       defaultcolors: $3Dmol.elementColors.rasmol,
       backgroundColor: 'black'
     }
   );
 *
 */
export function createViewer(element, config?:ViewerSpec) {
    element = getElement(element);
    if (!element) return;

    config = config || {};
    //try to create the  viewer
    try {
        var viewer = new GLViewer(element, config);
        return viewer;
    }
    catch (e) {
        throw "error creating viewer: " + e;
    }
};

/**
 * Create and initialize an appropriate a grid of viewers that share a WebGL canvas
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {GridSpec} [config] - grid configuration
 * @param {ViewerGridSpec} [viewer_config] - Viewer specification to apply to all subviewers
 * @return [[GLViewer]] 2D array of GLViewers
 * @example
   var viewers = $3Dmol.createViewerGrid(
     'gldiv', //id of div to create canvas in
     {
       rows: 2,
       cols: 2,
       control_all: true  //mouse controls all viewers
     },
     { backgroundColor: 'lightgrey' }
   );
   $3Dmol.get('data/1jpy.cif', function(data) {
     var viewer = viewers[0][0];
     viewer.addModel(data,'cif');
     viewer.setStyle({sphere:{}});
     viewer.zoomTo();
     viewer.render( );

     viewer = viewers[0][1];
     viewer.addModel(data,'cif');
     viewer.setStyle({stick:{}});
     viewer.zoomTo();
     viewer.render( );

     viewer = viewers[1][0];
     viewer.addModel(data,'cif');
     viewer.setStyle({cartoon:{color:'spectrum'}});
     viewer.zoomTo();
     viewer.render( );

     viewer = viewers[1][1];
     viewer.addModel(data,'cif');
     viewer.setStyle({cartoon:{colorscheme:'chain'}});
     viewer.zoomTo();
     viewer.render();


   });

 */
export function createViewerGrid(element, config:ViewerGridSpec={}, viewer_config:ViewerSpec={}) {
    element = getElement(element);
    if (!element) return;

    var viewers = [];
    //create canvas
    var canvas = document.createElement('canvas');

    viewer_config.rows = config.rows;
    viewer_config.cols = config.cols;
    viewer_config.control_all = config.control_all != undefined ? config.control_all : false;
    element.appendChild(canvas);

    //try to create the  viewer
    try {
        for (var r = 0; r < config.rows; r++) {
            var row = [];
            for (var c = 0; c < config.cols; c++) {
                viewer_config.row = r;
                viewer_config.col = c;
                viewer_config.canvas = canvas;
                viewer_config.viewers = viewers;
                viewer_config.control_all = config.control_all;
                var viewer = createViewer(element, viewer_config);
                row.push(viewer);
            }
            viewers.unshift(row); //compensate for weird ordering in renderer
        }
    } catch (e) {
        throw "error creating viewer grid: " + e;
    }

    return viewers;
};


/* StereoViewer for stereoscopic viewing
* @param {Object | string} element - Either HTML element or string identifier
*
*/
export function createStereoViewer(element) {
    var that = this;
    element = getElement(element);
    if (!element) return;

    var viewers = createViewerGrid(element, { rows: 1, cols: 2, control_all: true });

    this.glviewer1 = viewers[0][0];
    this.glviewer2 = viewers[0][1];

    this.glviewer1.setAutoEyeSeparation(false);
    this.glviewer2.setAutoEyeSeparation(true);

    this.glviewer1.linkViewer(this.glviewer2);
    this.glviewer2.linkViewer(this.glviewer1);

    var methods = Object.getOwnPropertyNames(this.glviewer1.__proto__) //get all methods of glviewer object
        .filter(function (property) {
            return typeof that.glviewer1[property] == 'function';
        });

    for (var i = 0; i < methods.length; i++) { //create methods of the same name
        this[methods[i]] = (function (method) {
            return function () {
                var m1 = this.glviewer1[method].apply(this.glviewer1, arguments);
                var m2 = this.glviewer2[method].apply(this.glviewer2, arguments);
                return [m1, m2];
            };
        })(methods[i]);
    }

    //special cased methods
    this.setCoordinates = function (models, data, format) { //for setting the coordinates of the models
        for (var i = 0; i < models.length; i++) {
            models[i].setCoordinates(data, format);
        }
    };

    this.surfacesFinished = function () {
        return this.glviewer1.surfacesFinished() && this.glviewer2.surfacesFinished();
    };

    this.isAnimated = function () {
        return this.glviewer1.isAnimated() || this.glviewer2.isAnimated();
    };

    this.render = function (callback) {
        this.glviewer1.render();
        this.glviewer2.render();
        if (callback) {
            callback(this); //call only once
        }
    };

    this.getCanvas = function () {
        return this.glviewer1.getCanvas(); //same for both
    };

};

/**
 * GLViewer input specification
 */
export interface ViewerSpec {
    /** Callback function to be executed with this viewer after setup is complete */
    callback?: (viewer: ViewerSpec) => void;
    /** Object defining default atom colors as atom => color property value pairs for all models within this viewer */
    defaultcolors?: Record<string, ColorSpec>;
    /**
     * Whether to disable disable handling of mouse events.
     * If you want to use your own mouse handlers, set this then bind your handlers to the canvas object.
                  The default 3Dmol.js handlers are available for use:
                  'mousedown touchstart': viewer._handleMouseDown,
                  'DOMMouseScroll mousewheel': viewer._handleMouseScroll
                  'mousemove touchmove': viewer._handleMouseMove
     */
    nomouse?: boolean | string;
    /** Color of the canvas background */
    backgroundColor?: string;
    /** Alpha transparency of canvas background */
    backgroundAlpha?: number;
    /** */
    camerax?: number;
    /** */
    hoverDuration?: number;
    /** id of the canvas */
    id?: string;
    /** default 5 */
    cartoonQuality?: number;
    /** */
    row?: number;
    /** */
    col?: number;
    /** */
    rows?: number;
    /** */
    cols?: number;
    /** */
    canvas?: HTMLCanvasElement;
    viewers?: GLViewer[];
    /** */
    minimumZoomToDistance?: number;
    /** */
    lowerZoomLimit?: number;
    /** */
    upperZoomLimit?: number;
    /** */
    antialias?: boolean;
    /** */
    control_all?: boolean;
    /** */
    orthographic?: boolean;
    /** Disable fog, default to false */
    disableFog?: boolean;

};

/**
 * Grid GLViewer input specification
 */
export interface ViewerGridSpec {
    /** number of rows in grid */
    rows?: number;
    /** number of columns in grid */
    cols?: number;
    /** if true, mouse events are linked */
    control_all?: boolean;
};


/**
 * @example
 * var setStyles = function(volumedata){
 *  var data = new $3Dmol.VolumeData(volumedata, "cube");
 *  viewer.addSurface("VDW", {opacity:0.85, voldata: data, volscheme: new $3Dmol.Gradient.RWB(-10,10)},{chain:'A'});
 *  viewer.mapAtomProperties($3Dmol.applyPartialCharges);
 *  viewer.addSurface($3Dmol.SurfaceType.SAS, {map:{prop:'partialCharge',scheme:new $3Dmol.Gradient.RWB(-.05,.05)}, opacity:1.0},{chain:'B'});
 *  viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: data, color:'red'},{chain:'C'});
 *  viewer.addSurface($3Dmol.SurfaceType.SAS, {opacity:0.85,voldata: data, colorscheme:'greenCarbon'},{chain:'D'});
 *  viewer.render();
 * };
 * $3Dmol.download("pdb:4DLN",viewer,{},function(){
 *   $.get("data/1fas.cube",setStyles);
 * });
 */
export interface SurfaceStyleSpec {
    /** sets the transparency: 0 to hide, 1 for fully opaque */
    opacity?: number;
    /** element based coloring */
    colorscheme?: ColorschemeSpec;
    /** fixed coloring, overrides colorscheme */
    color?: ColorSpec;
    /** volumetric data for vertex coloring, can be VolumeData object or raw data if volformat is specified */
    voldata?: VolumeData;
    /** coloring scheme for mapping volumetric data to vertex color, if not a Gradient object, show describe a builtin gradient one by providing an object with gradient, min, max, and (optionally) mid fields. */
    volscheme?: Gradient;
    /** format of voldata if not a {VolumeData} object */
    volformat?: string;
    /* specifies a numeric atom property (prop) and color mapping (scheme) such as {@link $3Dmol.Gradient.RWB}.  Deprecated, use colorscheme instead. */
    map?: Record<string, unknown>
};


/** Style specification ofr unit cell shape.  */
export interface UnitCellStyleSpec {
    /** line style used to draw box */
    box?: LineStyleSpec;
    /** arrow specification of the "a" axis */
    astyle?: ArrowSpec;
    /** arrow specification of the "b" axis */
    bstyle?: ArrowSpec;
    /** arrow specification of the "c" axis */
    cstyle?: ArrowSpec;
    /** label for "a" axis */
    alabel?: string;
    /** label style for a axis */
    alabelstyle?: LabelSpec;
    /** label for "b" axis */
    blabel?: string;
    /** label style for b axis */
    blabelstyle?: LabelSpec;
    /** label for "c" axis */
    clabel?: string;
    /** label style for c axis */
    clabelstyle?: LabelSpec;
}
