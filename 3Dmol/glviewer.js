/* eslint-disable no-underscore-dangle */
// @ts-check
/* eslint-disable no-multi-assign */
/* eslint-disable eqeqeq */
import UPNG from 'upng-js';
import $ from 'jquery';
import Renderer from './WebGL/Renderer';
import Camera from './WebGL/Camera';
import {Object3D, Light, Raycaster, Scene, Projection, Geometry} from './WebGL/core';
import Fog from './WebGL/Fog';
import {Vector3, Quaternion, Matrix3, Matrix4} from './WebGL/math';
import {
  LineBasicMaterial,
  FrontSide,
  MeshLambertMaterial,
  VertexColors,
} from './WebGL/materials';
import {Line, Mesh} from './WebGL/objects';
import {elementColors, CC, getColorFromStyle} from './colors';
import Gradient from './Gradient';
import GLShape, {splitMesh} from './glshape';
import makeFunction from './util/makeFunction';
import getExtent from './util/getExtent';
import isEmptyObject from './util/isEmptyObject';
import SurfaceType, { normalizeSurfaceType } from './enum/SurfaceType';
import mergeGeos from './util/mergeGeos';
import SyncSurface from './util/SyncSurface';
import extend from './util/extend';
import GLModel from './GLModel';
import getPropertyRange from './util/getPropertyRange';
import adjustVolumeStyle from './util/adjustVolumeStyle';
import Label from './Label';
import ProteinSurface from './ProteinSurface4';
import {VolumeData, GLVolumetricRender} from './volume';
import StateManager from './ui/StateManager';
import Viewers from './singletons/Viewers';
import { SurfaceWorker } from './SurfaceWorker';

// a molecular viewer based on GLMol

// private class variables
const numWorkers = 4; // number of threads for surface generation
const maxVolume = 64000; // how much to break up surface calculations

// touch event math helpers
function calcTouchDistance(ev) {
  // distance between first two
  // fingers
  const xdiff = ev.originalEvent.targetTouches[0].pageX - ev.originalEvent.targetTouches[1].pageX;
  const ydiff = ev.originalEvent.targetTouches[0].pageY - ev.originalEvent.targetTouches[1].pageY;
  return Math.sqrt(xdiff * xdiff + ydiff * ydiff);
}

// check targetTouches as well
function getX(ev) {
  let x = ev.pageX;
  if (x == undefined) x = ev.originalEvent.pageX; // firefox
  if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
    x = ev.originalEvent.targetTouches[0].pageX;
  } else if (ev.originalEvent.changedTouches && ev.originalEvent.changedTouches[0]) {
    x = ev.originalEvent.changedTouches[0].pageX;
  }
  return x;
}

function getY(ev) {
  let y = ev.pageY;
  if (y == undefined) y = ev.originalEvent.pageY;
  if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
    y = ev.originalEvent.targetTouches[0].pageY;
  } else if (ev.originalEvent.changedTouches && ev.originalEvent.changedTouches[0]) {
    y = ev.originalEvent.changedTouches[0].pageY;
  }
  return y;
}

// interpolate between two normalized quaternions (t between 0 and 1)
// https://en.wikipedia.org/wiki/Slerp
function slerp(v0, v1, t) {
  // Compute the cosine of the angle between the two vectors.
  // dot product
  if (t == 1) return v1;
  if (t == 0) return v0;
  let dot = v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w * v1.w;
  if (dot > 0.9995) {
    // If the inputs are too close for comfort, linearly interpolate
    // and normalize the result.
    const result = new Quaternion(
      v0.x + t * (v1.x - v0.x),
      v0.y + t * (v1.y - v0.y),
      v0.z + t * (v1.z - v0.z),
      v0.w + t * (v1.w - v0.w)
    );

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

  const theta0 = Math.acos(dot); // theta_0 = angle between input vectors
  const theta = theta0 * t; // theta = angle between v0 and result

  const v2 = v1.clone();
  v2.sub(v0.clone().multiplyScalar(dot));
  v2.normalize(); // { v0, v2 } is now an orthonormal basis

  const c = Math.cos(theta);
  const s = Math.sin(theta);
  const ret = new Quaternion(
    v0.x * c + v2.x * s,
    v0.y * c + v2.y * s,
    v0.z * c + v2.z * s,
    v0.w * c + v2.w * s
  );
  ret.normalize();
  return ret;
}

/**
 * WebGL-based 3Dmol.js viewer
 * Note: The preferred method of instantiating a GLViewer is through {@link createViewer}
 */
export default class GLViewer {
  // Scene data
  /** @type {Array<GLModel>} */
  models = []; // atomistic molecular models
  surfaces = {};
  shapes = []; // Generic shapes
  labels = [];
  fixed_labels = [];
  clickables = []; // things you can click on
  hoverables = []; // things you can hover over
  contextMenuEnabledAtoms = []; // atoms with context menu
  scene = new Scene();
  rotationGroup = new Object3D(); // which contains modelGroup
  modelGroup = new Object3D();

  config;

  // State data
  clickedAtom = null;
  glDOM;
  _stateManager;
  callback;
  /** @type {any} */
  current_hover = null;
  /** @type {number | undefined} */
  hoverDuration = 500;
  viewer_frame = 0;
  bgColor = 0;
  linkedViewers = [];

  // Child Compoenents
  raycaster = new Raycaster(new Vector3(0, 0, 0), new Vector3(0, 0, 0));
  projector = new Projection();

  // camera const
  NEAR = 1;
  FAR = 800;
  CAMERA_Z = 150;
  // camera dynamics
  /** @type {number|undefined} */
  camerax = 0;
  fov = 20;
  lookingAt = new Vector3();
  fogStart = 0.4;
  slabNear = -50; // relative to the center of rotationGroup
  slabFar = 50;

  // UI variables
  cq = new Quaternion(0, 0, 0, 1);
  dq = new Quaternion(0, 0, 0, 1);
  animated = 0;
  animationTimers = new Set();
  isDragging = false;
  mouseStartX = 0;
  mouseStartY = 0;
  touchDistanceStart = 0;
  touchHold = false;
  /** @type {number|Vector3} */
  currentModelPos = 0;
  cz = 0;
  cslabNear = 0;
  cslabFar = 0;

  // uninitialized members
  control_all;
  cols;
  rows;
  row;
  col;
  nomouse;

  /**
   * @param {Object} element - HTML element within which to create viewer
   * @param {import('./specs').ViewerSpec} cfg - Object containing optional this.configuration for the viewer
   * @param {any} [svr]
   */
  constructor(element, cfg, svr) {
    // config defaults
    this.config = cfg || {};
    this.config.defaultcolors = this.config.defaultColors || elementColors.defaultColors;
    this.config.backgroundColor = this.config.backgroundColor || '#ffffff';
    // this.config.disableFog= this.config.disableFog || false;
    this.bgColor = CC.color(this.config.backgroundColor || 0).getHex();
    this.config.backgroundAlpha = this.config.backgroundAlpha || 1.0;
    this.camerax = this.config.camerax
      ? parseFloat(/** @type {string} */ (this.config.camerax))
      : undefined;
    this.container = $(element); // we expect container to be jquery
    this.hoverDuration = this.config.hoverDuration || undefined;
    this.config.antialias = this.config.antialias || true;
    this.config.cartoonQuality = this.config.cartoonQuality || 5;

    this.renderer = new Renderer({
      antialias: this.config.antialias,
      preserveDrawingBuffer: true,
      premultipliedAlpha: false /* more traditional compositing with background */,
      id: this.config.id,
      row: this.config.row,
      col: this.config.col,
      rows: this.config.rows,
      cols: this.config.cols,
      canvas: this.config.canvas,
      // cannot initialize with zero size
      containerWidth: this.WIDTH || 1,
      containerHeight: this.HEIGHT || 1,
    });
    this.renderer.domElement.style.width = '100%';
    this.renderer.domElement.style.height = '100%';
    this.renderer.domElement.style.padding = '0';
    this.renderer.domElement.style.position = 'absolute'; // TODO: get rid of this
    this.renderer.domElement.style.top = '0px';
    this.renderer.domElement.style.left = '0px';
    this.renderer.domElement.style.zIndex = '0';

    this.ASPECT = this.renderer.getAspect(this.WIDTH, this.HEIGHT);

    this.camera = new Camera(this.fov, this.ASPECT, this.NEAR, this.FAR, this.config.orthographic);
    this.camera.position = new Vector3(this.camerax, 0, this.CAMERA_Z);
    this.camera.lookAt(this.lookingAt);

    // Initialize the scene
    this.scene.fog = new Fog(this.bgColor, 100, 200);
    this.rotationGroup.useQuaternion = true;
    this.rotationGroup.quaternion = new Quaternion(0, 0, 0, 1);
    this.rotationGroup.add(this.modelGroup);
    this.scene.add(this.rotationGroup);

    // setup lights
    const directionalLight = new Light(0xffffff);
    directionalLight.position = new Vector3(0.2, 0.2, 1).normalize();
    directionalLight.intensity = 1.0;
    this.scene.add(directionalLight);

    this.renderer.setClearColorHex(this.bgColor, this.config.backgroundAlpha);
    this.scene.fog.color = CC.color(this.bgColor);

    // this event is bound to the body element, not the container
    $('body').on('mouseup touchend', ev => {
      // handle touch
      this.touchHold = false;

      // handle selection
      if (this.isDragging && this.scene) {
        // saw mousedown, haven't moved
        const x = getX(ev);
        const y = getY(ev);
        if (x == this.mouseStartX && y == this.mouseStartY) {
          const offset = this.canvasOffset();
          const mouseX = ((x - offset.left) / this.WIDTH) * 2 - 1;
          const mouseY = -((y - offset.top) / this.HEIGHT) * 2 + 1;
          this.handleClickSelection(mouseX, mouseY, ev);
        }
      }

      this.isDragging = false;
    });

    // eslint-disable-next-line no-underscore-dangle
    this._stateManager = new StateManager(this, this.config); // Creates the UI state management tool

    this.initContainer(this.container);

    if (this.config.style) {
      // enable setting style in constructor
      this.setViewStyle(this.config);
    }

    $(window).resize(this.resize);

    if (typeof window.ResizeObserver !== 'undefined') {
      const divwatcher = new window.ResizeObserver(this.resize);
      divwatcher.observe(this.container[0]);
    }

    try {
      if (typeof this.callback === 'function') this.callback(this);
    } catch (e) {
      // errors in callback shouldn't invalidate the viewer
      console.log(`error with glviewer callback: ${e}`);
    }

    /**
     * TODO: extract this into an internal class
     */
    this.ui = {};

    // State Management function
    /**
     * Calls StateManager to add selection and style on the ui
     * @function GLViewer#ui.loadSelectionStyle
     * @param {import('./specs').AtomSelectionSpec} sel Atom Selection spec
     * @param {import('./specs').AtomStyleSpec} style Atom Style spec
     */
    this.ui.loadSelectionStyle = (sel, style) => {
      this._stateManager.createSelectionAndStyle(sel, style);
    };

    /**
     * Calls StateManager to add surface and selection on the ui
     *
     * @function GLViewer#ui.loadSurface
     * @param {String} surfaceType Name of the surface type
     * @param {import('./specs').AtomSelectionSpec} sel Atom Selection Spec
     * @param {import('./specs').AtomStyleSpec} style Atom Style Spec
     * @param {string} sid Id of the surface that is already added
     */
    this.ui.loadSurface = (surfaceType, sel, style, sid) => {
      this._stateManager.createSurface(surfaceType, sel, style, sid);
    };

    /**
     * Sets the name of the file as title in the ui
     *
     * @function GlViewer#ui.setModelTitle
     * @param {String} title Name of file loaded
     */
    this.ui.setModelTitle = title => {
      this._stateManager.setModelTitle(title);
    };

    /**
     * Calls StateManager to start the
     * @function GLViewer#ui.initiateUI
     */
    this.ui.initiateUI = () => {
      this._stateManager.initiateUI();
    };
  }

  /**
   * Set a callback to call when the view has potentially changed.
   *
   * @function GLViewer#setViewChangeCallback
   */
  setViewChangeCallback(callback) {
    if (typeof callback === 'function' || callback == null) this.viewChangeCallback = callback;
  }

  /**
   * Set a callback to call when the view has potentially changed.
   *
   * @function GLViewer#setStateChangeCallback
   */
  setStateChangeCallback(callback) {
    if (typeof callback === 'function' || callback == null) this.stateChangeCallback = callback;
  }

  /**
   * Return this.configuration of viewer
   */
  getConfig() {
    return this.config;
  }

  /**
   * Set the this.configuration object.  Note that some setting may only
   * have an effect at viewer creation time.
   */
  setConfig(c) {
    this.config = c;
  }

  /**
   * Return object representing internal state of
   * the viewer appropriate for passing to setInternalState
   *
   * @function GLViewer#getInternalState
   */
  getInternalState() {
    /** @type {{models: any[], surfaces: any[], shapes: any[], labels: any[]}} */
    const ret = {models: [], surfaces: [], shapes: [], labels: []};
    for (let i = 0; i < this.models.length; i++) {
      if (this.models[i]) {
        ret.models[i] = this.models[i].getInternalState();
      }
    }

    // todo: labels, shapes, surfaces
    return ret;
  }

  /**
   * Overwrite internal state of the viewer with passed  object
   * which should come from getInternalState.
   *
   * @function GLViewer#setInternalState
   */
  setInternalState(state) {
    // clear out current viewer
    this.clear();

    // set model state
    const newm = state.models;
    for (let i = 0; i < newm.length; i++) {
      if (newm[i]) {
        this.models[i] = new GLModel(i);
        this.models[i].setInternalState(newm[i]);
      }
    }

    // todo: labels, shapes, surfaces
    this.render();
  }

  /**
           * Set lower and upper limit stops for zoom.
           *
           * @function GLViewer#setZoomLimits
           * @param {number} lower - limit on zoom in (positive number).  Default 0.
           * @param {number} upper - limit on zoom out (positive number).  Default infinite.
           * @example
            $.get("data/set1_122_complex.mol2", function(moldata) {
                  var m = viewer.addModel(moldata);
                  viewer.setStyle({stick:{colorscheme:"Jmol"}});
                  viewer.setZoomLimits(100,200);
                  viewer.zoomTo();
                  viewer.zoom(10); //will not zoom all the way
                  viewer.render();
              });
          */
  setZoomLimits(lower, upper) {
    if (typeof lower !== 'undefined') this.config.lowerZoomLimit = lower;
    if (upper) this.config.upperZoomLimit = upper;
    this.rotationGroup.position.z = this.adjustZoomToLimits(this.rotationGroup.position.z);
    this.show();
  }

  /**
           * Set camera parameters (distance to the origin and field of view)
           *
           * @function GLViewer#setCameraParameters
           * @param {{fov: number, z: number, orthographic?: boolean}} parameters - new camera parameters, with possible fields
           *                       being fov for the field of view, z for the
           *                       distance to the origin, and orthographic (boolean)
           *                       for kind of projection (default false).
           * @example
            $.get("data/set1_122_complex.mol2", function(data) {
                  var m = viewer.addModel(data);
                  viewer.setStyle({stick:{}});
                  viewer.zoomTo();
                  viewer.setCameraParameters({ fov: 10 , z: 300 });
                  viewer.render();
              });
          */
  setCameraParameters(parameters) {
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
  }

  // reimplement jquery getwidth/height
  getRect() {
    const div = this.container[0];
    let rect = div.getBoundingClientRect();
    if (rect.width == 0 && rect.height == 0 && div.style.display === 'none') {
      const oldpos = div.style.position;
      const oldvis = div.style.visibility;
      div.style.display = 'block';
      div.style.visibility = 'hidden';
      div.style.position = 'absolute';
      rect = div.getBoundingClientRect();
      div.style.display = 'none';
      div.style.visibility = oldvis;
      div.style.position = oldpos;
    }
    return rect;
  }

  getWidth() {
    return this.getRect().width;
  }

  get WIDTH() {
    return this.getRect().width;
  }

  getHeight() {
    return this.getRect().height;
  }

  get HEIGHT() {
    return this.getRect().height;
  }

  decAnim() {
    // decrement the number of animations currently
    this.animated--;
    if (this.animated < 0) this.animated = 0;
  }

  incAnim() {
    this.animated++;
  }

  nextSurfID() {
    // compute the next highest surface id directly from surfaces
    // this is necessary to support linking of model data
    let max = 0;
    for (let i in this.surfaces) {
      // this is an object with possible holes
      if (!this.surfaces[i]) continue;
      const val = Number.parseInt(i, 10);
      // @ts-ignore
      if (!Number.isNaN(val)) i = val;
      // @ts-ignore
      if (i > max) max = i;
    }
    return max + 1;
  }

  setSlabAndFog() {
    let center = this.camera.position.z - this.rotationGroup.position.z;
    if (center < 1) center = 1;
    this.camera.near = center + this.slabNear;
    if (this.camera.near < 1) this.camera.near = 1;
    this.camera.far = center + this.slabFar;
    if (this.camera.near + 1 > this.camera.far) this.camera.far = this.camera.near + 1;

    this.camera.fov = this.fov;
    this.camera.right = center * Math.tan((Math.PI / 180) * this.fov);
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
  }

  // display scene
  // if nolink is set/true, don't propagate changes to linked viewers
  show(nolink) {
    this.renderer.setViewport();
    if (!this.scene) return;
    // var time = new Date();
    this.setSlabAndFog();
    this.renderer.render(this.scene, this.camera);
    // console.log("rendered in " + (+new Date() - time) + "ms");
    // have any scene change trigger a callback
    if (this.viewChangeCallback) this.viewChangeCallback(this.getView());

    if (!nolink && this.linkedViewers.length > 0) {
      const view = this.getView();
      for (let i = 0; i < this.linkedViewers.length; i++) {
        const other = this.linkedViewers[i];
        other.setView(view, true);
      }
    }
  }

  // enable mouse support
  // regenerate the list of clickables
  // also updates hoverables
  updateClickables() {
    this.clickables.splice(0, this.clickables.length);
    this.hoverables.splice(0, this.hoverables.length);
    this.contextMenuEnabledAtoms.splice(0, this.contextMenuEnabledAtoms.length);

    for (let i = 0, il = this.models.length; i < il; i++) {
      const model = this.models[i];
      if (model) {
        const atoms = model.selectedAtoms({
          clickable: true,
        });

        const hoverableAtoms = model.selectedAtoms({
          hoverable: true,
        });

        const contextMenuEnabledAtom = model.selectedAtoms({contextMenuEnabled: true});
        // Array.prototype.push.apply(hoverables,hoverable_atoms);
        for (let n = 0; n < hoverableAtoms.length; n++) {
          this.hoverables.push(hoverableAtoms[n]);
        }

        // Array.prototype.push.apply(clickables, atoms); //add atoms into clickables
        for (let m = 0; m < atoms.length; m++) {
          this.clickables.push(atoms[m]);
        }

        // add atoms into contextMenuEnabledAtoms
        for (let m = 0; m < contextMenuEnabledAtom.length; m++) {
          this.contextMenuEnabledAtoms.push(contextMenuEnabledAtom[m]);
        }
      }
    }
    for (let i = 0, il = this.shapes.length; i < il; i++) {
      const shape = this.shapes[i];
      if (shape && shape.clickable) {
        this.clickables.push(shape);
      }
      if (shape && shape.hoverable) {
        this.hoverables.push(shape);
      }
    }
  }

  /**
   * Return a list of objects that intersect that at the specified viewer position.
   *
   * @function GLViewer#targetedObjects
   * @param {number} x - x position in screen coordinates
   * @param {number} y - y position in screen coordinates
   * @param {any[]} objects - list of objects or selection object specifying what object to check for targeting
   */
  targetedObjects(x, y, objects) {
    const mouse = {
      x,
      y,
      z: -1.0,
    };
    if (!Array.isArray(objects)) {
      // assume selection object
      objects = this.selectedAtoms(objects);
    }
    if (objects.length == 0) return [];
    this.raycaster.setFromCamera(mouse, this.camera);
    return this.raycaster.intersectObjects(this.modelGroup, objects);
  }

  // return offset of container
  canvasOffset() {
    const canvas = this.glDOM.get(0);
    const rect = canvas.getBoundingClientRect();
    const doc = canvas.ownerDocument;
    const docElem = doc.documentElement;
    const win = doc.defaultView;
    return {
      top: rect.top + win.pageYOffset - docElem.clientTop,
      left: rect.left + win.pageXOffset - docElem.clientLeft,
    };
  }

  /** Convert model coordinates to screen coordinates.
   * @function GLViewer#modelToScreen
   * @param {import('./specs').Vector3Like | Array<import('./specs').Vector3Like>} coords - an object or list of objects with x,y,z attributes (e.g. an atom)
   * @return {import("./WebGL/math").Vector2 | Array<import("./WebGL/math").Vector2>} - and object or list of {x: screenX, y: screenY}
   */
  modelToScreen(coords) {
    let returnsingle = false;
    if (!Array.isArray(coords)) {
      coords = [coords];
      returnsingle = true;
    }

    let results = [];
    const offset = this.canvasOffset();
    coords.forEach(coord => {
      const t = new Vector3(coord.x, coord.y, coord.z);
      t.applyMatrix4(this.modelGroup.matrixWorld);
      this.projector.projectVector(t, this.camera);
      const screenX = (this.WIDTH * (t.x + 1)) / 2.0 + offset.left;
      const screenY = (-this.HEIGHT * (t.y - 1)) / 2.0 + offset.top;
      results.push({x: screenX, y: screenY});
    });
    if (returnsingle) results = results[0];
    return results;
  }

  // Checks for selection intersects on mousedown
  handleClickSelection(mouseX, mouseY, event) {
    const intersects = this.targetedObjects(mouseX, mouseY, this.clickables);
    // console.log('handleClickSelection', mouseX, mouseY, intersects);
    if (intersects.length) {
      const selected = intersects[0].clickable;
      if (selected.callback !== undefined) {
        if (typeof selected.callback != 'function') {
          selected.callback = makeFunction(selected.callback);
        }
        if (typeof selected.callback === 'function') {
          selected.callback(selected, this, event, this.container);
        }
      }
    }
  }

  // set current_hover to sel (which can be null), calling appropraite callbacks
  setHover(selected, event) {
    if (this.current_hover == selected) return;
    if (this.current_hover) {
      if (typeof this.current_hover.unhover_callback != 'function') {
        this.current_hover.unhover_callback = makeFunction(this.current_hover.unhover_callback);
      }
      this.current_hover.unhover_callback(this.current_hover, this, event, this.container);
    }
    this.current_hover = selected;

    if (selected && selected.hover_callback !== undefined) {
      if (typeof selected.hover_callback != 'function') {
        selected.hover_callback = makeFunction(selected.hover_callback);
      }
      if (typeof selected.hover_callback === 'function') {
        selected.hover_callback(selected, this, event, this.container);
      }
    }
  }

  // checks for selection intersects on hover
  handleHoverSelection(mouseX, mouseY) {
    if (this.hoverables.length == 0) return;
    const intersects = this.targetedObjects(mouseX, mouseY, this.hoverables);
    if (intersects.length) {
      const selected = intersects[0].clickable;
      this.setHover(selected);
      this.current_hover = selected;
    } else {
      this.setHover(null);
    }
  }

  // sees if the mouse is still on the object that invoked a hover event and if not then the unhover callback is called
  handleHoverContinue(mouseX, mouseY) {
    const intersects = this.targetedObjects(mouseX, mouseY, this.hoverables);
    if (intersects.length == 0 || intersects[0] === undefined) {
      this.setHover(null);
    }
    if (intersects[0] !== undefined && intersects[0].clickable !== this.current_hover) {
      this.setHover(null);
    }
  }

  /**
   * For a given screen (x,y) displacement return model displacement
   * @param {number} x - displacement in screen coordinates
   * @param {number} y - displacement in screen corodinates
   * @param {number|undefined} [modelz] - coordinate in model coordinates to compute offset for, default is model axis
   * @returns {Vector3} - displacement in model coordinates
   */
  screenOffsetToModel(x, y, modelz) {
    const dx = x / this.WIDTH;
    const dy = y / this.HEIGHT;
    const zpos = modelz === undefined ? this.rotationGroup.position.z : modelz;
    const q = this.rotationGroup.quaternion;
    const t = new Vector3(0, 0, zpos);
    this.projector.projectVector(t, this.camera);
    t.x += dx * 2;
    t.y -= dy * 2;
    this.projector.unprojectVector(t, this.camera);
    t.z = 0;
    t.applyQuaternion(q);
    return t;
  }

  /**
   * Distance from screen coordinate to model coordinate assuming screen point
   * is projected to the same depth as model coordinate
   * @param {import("./WebGL/math").Vector2} screen - xy screen coordinate
   * @param {Vector3} model - xyz model coordinate
   *
   */
  screenToModelDistance(screen, model) {
    const offset = this.canvasOffset();

    // convert model to screen to get screen z
    const mvec = new Vector3(model.x, model.y, model.z);
    mvec.applyMatrix4(this.modelGroup.matrixWorld);
    const m = mvec.clone();
    this.projector.projectVector(mvec, this.camera);

    const t = new Vector3(
      ((screen.x - offset.left) * 2) / this.WIDTH - 1,
      ((screen.y - offset.top) * 2) / -this.HEIGHT + 1,
      mvec.z
    );
    this.projector.unprojectVector(t, this.camera);

    return t.distanceTo(m);
  }

  // for grid viewers, return true if point is in this viewer
  isInViewer(x, y) {
    if (Viewers != undefined && !this.control_all) {
      const width = this.WIDTH / this.cols;
      const height = this.HEIGHT / this.rows;
      const offset = this.canvasOffset();
      const relx = x - offset.left;
      const rely = y - offset.top;

      const r = this.rows - Math.floor(rely / height) - 1;
      const c = Math.floor(relx / width);

      if (r != this.row || c != this.col) return false;
    }
    return true;
  }

  // if the user has specify zoom limits, readjust to fit within them
  // also, make sure we don't go past CAMERA_Z
  adjustZoomToLimits(z) {
    // a lower limit of 0 is at CAMERA_Z
    if (this.config.lowerZoomLimit && this.config.lowerZoomLimit > 0) {
      const lower = this.CAMERA_Z - this.config.lowerZoomLimit;
      if (z > lower) z = lower;
    }

    if (this.config.upperZoomLimit && this.config.upperZoomLimit > 0) {
      const upper = this.CAMERA_Z - this.config.upperZoomLimit;
      if (z < upper) z = upper;
    }

    if (z > this.CAMERA_Z) {
      z = this.CAMERA_Z * 0.999; // avoid getting stuck
    }
    return z;
  }

  mouseButton;
  _handleMouseDown(ev) {
    ev.preventDefault();
    if (!this.scene) return;
    const x = getX(ev);
    const y = getY(ev);
    if (x === undefined) return;
    this.isDragging = true;
    this.clickedAtom = null;
    this.mouseButton = ev.which;
    this.mouseStartX = x;
    this.mouseStartY = y;
    this.touchHold = true;
    this.touchDistanceStart = 0;
    if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches.length == 2) {
      this.touchDistanceStart = calcTouchDistance(ev);
    }
    this.cq = this.rotationGroup.quaternion.clone();
    this.cz = this.rotationGroup.position.z;
    this.currentModelPos = this.modelGroup.position.clone();
    this.cslabNear = this.slabNear;
    this.cslabFar = this.slabFar;

    setTimeout(() => {
      if (ev.originalEvent.targetTouches) {
        if (this.touchHold == true) {
          // console.log('Touch hold', x,y);
          this.glDOM = $(this.renderer.domElement);
          this.glDOM.trigger('contextmenu');
        } else {
          // console.log('Touch hold ended earlier');
        }
      }
    }, 1000);

    // Exiting context menu on next mouse event
    this._stateManager.exitContextMenu();
  }

  _handleMouseScroll(ev) {
    // Zoom
    ev.preventDefault();
    if (!this.scene) return;

    const x = getX(ev);
    const y = getY(ev);
    if (x === undefined) return;
    if (!this.isInViewer(x, y)) {
      return;
    }

    const scaleFactor = (this.CAMERA_Z - this.rotationGroup.position.z) * 0.85;
    let mult = 1.0;
    if (ev.originalEvent.ctrlKey) {
      mult = -1.0; // this is a pinch event turned into a wheel event (or they're just holding down the ctrl)
    }
    if (ev.originalEvent.detail) {
      this.rotationGroup.position.z += (mult * scaleFactor * ev.originalEvent.detail) / 10;
    } else if (ev.originalEvent.wheelDelta) {
      this.rotationGroup.position.z -= (mult * scaleFactor * ev.originalEvent.wheelDelta) / 400;
    }
    this.rotationGroup.position.z = this.adjustZoomToLimits(this.rotationGroup.position.z);
    this.show();

    // Exiting context menu on next mouse event
    this._stateManager.exitContextMenu();
  }

  /**
   * Return image URI of viewer contents (base64 encoded).
   * @returns {string}
   *
   */
  pngURI() {
    return this.getCanvas().toDataURL('image/png');
  }

  /**
   * Return a promise that resolves to an animated PNG image URI of
   viewer contents (base64 encoded) for nframes of viewer changes.
   * @function GLViewer#apngURI
   * @return {Promise<any>}
   */
  apngURI(nframes) {
    const viewer = this;
    nframes = nframes || 1;
    return new Promise(resolve => {
      let framecnt = 0;
      const oldcb = this.viewChangeCallback;
      const bufpromise = [];
      const delays = [];
      let lasttime = Date.now();
      this.viewChangeCallback = function () {
        delays.push(Date.now() - lasttime);
        lasttime = Date.now();
        bufpromise.push(
          new Promise(resolve => {
            viewer.getCanvas().toBlob(blob => {
              if (blob) blob.arrayBuffer().then(resolve);
            }, 'image/png');
          })
        );
        framecnt += 1;
        if (framecnt == nframes) {
          this.viewChangeCallback = oldcb;

          Promise.all(bufpromise).then(buffers => {
            // convert to apng
            const rgbas = [];
            // have to convert png to rgba, before creating the apng
            for (let i = 0; i < buffers.length; i++) {
              const img = UPNG.decode(buffers[i]);
              rgbas.push(UPNG.toRGBA8(img)[0]);
            }
            const {width} = viewer.getCanvas();
            const {height} = viewer.getCanvas();
            const apng = UPNG.encode(rgbas, width, height, 0, delays);
            const blob = new Blob([apng], {type: 'image/png'});
            const fr = new FileReader();
            fr.onload = e => {
              resolve((e.target || {}).result);
            };
            fr.readAsDataURL(blob);
          });
        }
      };
    });
  }

  /**
   * @returns {HTMLCanvasElement}
   */
  getCanvas() {
    return this.glDOM.get(0);
  }

  /**
   * @returns {Renderer}
   */
  getRenderer() {
    return this.renderer;
  }

  /**
   * Set the duration of the hover delay
   *
   * @function GLViewer#setHoverDuration
   * @param {number} [duration] - the duration of the hover delay (in milliseconds) before the hover action is called
   *
   */
  setHoverDuration(duration) {
    this.hoverDuration = duration;
  }

  hoverTimeout;
  _handleMouseMove(ev) {
    // touchmove
    clearTimeout(this.hoverTimeout);
    const offset = this.canvasOffset();
    const mouseX = ((getX(ev) - offset.left) / this.WIDTH) * 2 - 1;
    const mouseY = -((getY(ev) - offset.top) / this.HEIGHT) * 2 + 1;

    // hover timeout
    if (this.current_hover !== null) {
      this.handleHoverContinue(mouseX, mouseY);
    }

    if (this.hoverables.length > 0) {
      this.hoverTimeout = setTimeout(() => {
        this.handleHoverSelection(mouseX, mouseY);
      }, this.hoverDuration);
    }

    ev.preventDefault();
    if (!this.scene) return;
    if (!this.isDragging) return;
    let mode = 0;

    const x = getX(ev);
    const y = getY(ev);
    if (x === undefined) return;

    if (!this.isInViewer(x, y)) {
      return;
    }

    let dx = (x - this.mouseStartX) / this.WIDTH;
    let dy = (y - this.mouseStartY) / this.HEIGHT;
    // check for pinch
    if (
      this.touchDistanceStart != 0 &&
      ev.originalEvent.targetTouches &&
      ev.originalEvent.targetTouches.length == 2
    ) {
      const newdist = calcTouchDistance(ev);
      // change to zoom
      mode = 2;
      dy = ((newdist - this.touchDistanceStart) * 2) / (this.WIDTH + this.HEIGHT);
    } else if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches.length == 3) {
      // translate
      mode = 1;
    }
    const ratioX = this.renderer.getXRatio();
    const ratioY = this.renderer.getYRatio();
    dx *= ratioX;
    dy *= ratioY;
    const r = Math.sqrt(dx * dx + dy * dy);
    let scaleFactor;
    if (mode == 3 || (this.mouseButton == 3 && ev.ctrlKey)) {
      // Slab
      this.slabNear = this.cslabNear + dx * 100;
      this.slabFar = this.cslabFar - dy * 100;
    } else if (mode == 2 || this.mouseButton == 3 || ev.shiftKey) {
      // Zoom
      scaleFactor = (this.CAMERA_Z - this.rotationGroup.position.z) * 0.85;
      if (scaleFactor < 80) scaleFactor = 80;
      this.rotationGroup.position.z = this.cz + dy * scaleFactor;
      this.rotationGroup.position.z = this.adjustZoomToLimits(this.rotationGroup.position.z);
    } else if (mode == 1 || this.mouseButton == 2 || ev.ctrlKey) {
      // Translate
      const t = this.screenOffsetToModel(
        ratioX * (x - this.mouseStartX),
        ratioY * (y - this.mouseStartY)
      );
      this.modelGroup.position.addVectors(this.currentModelPos, t);
    } else if ((mode === 0 || this.mouseButton == 1) && r !== 0) {
      // Rotate
      const rs = Math.sin(r * Math.PI) / r;
      this.dq.x = Math.cos(r * Math.PI);
      this.dq.y = 0;
      this.dq.z = rs * dx;
      this.dq.w = -rs * dy;
      this.rotationGroup.quaternion.set(1, 0, 0, 0);
      this.rotationGroup.quaternion.multiply(this.dq);
      this.rotationGroup.quaternion.multiply(this.cq);
    }
    this.show();
  }

  handleContextMenuSelection(mouseX, mouseY) {
    const intersects = this.targetedObjects(mouseX, mouseY, this.contextMenuEnabledAtoms);
    // console.log('Intersected Objects',mouseX, mouseY, contextMenuEnabledAtoms,  intersects[0]);
    let selected = null;
    if (intersects.length) {
      selected = intersects[0].clickable;
      // console.log('intersects and selected', selected);
    }

    const offset = this.canvasOffset();
    const x = this.mouseStartX - offset.left;
    const y = this.mouseStartY - offset.top;
    this._stateManager.openContextMenu(selected, x, y);
  }

  // eslint-disable-next-line no-underscore-dangle
  _handleContextMenu(ev) {
    ev.preventDefault();
    const newX = getX(ev);
    const newY = getY(ev);

    if (!(newX != this.mouseStartX || newY != this.mouseStartY)) {
      // console.log('Context Menu Called', ev);
      const x = this.mouseStartX;
      const y = this.mouseStartY;
      const offset = this.canvasOffset();
      const mouseX = ((x - offset.left) / this.WIDTH) * 2 - 1;
      const mouseY = -((y - offset.top) / this.HEIGHT) * 2 + 1;
      this.handleContextMenuSelection(mouseX, mouseY);
    }
  }

  initContainer(element) {
    this.container = element;
    this.ASPECT = this.renderer.getAspect(this.WIDTH, this.HEIGHT);
    this.renderer.setSize(this.WIDTH, this.HEIGHT);
    this.container.append(this.renderer.domElement);
    this.glDOM = $(this.renderer.domElement);

    if (!this.nomouse) {
      // user can request that the mouse handlers not be installed
      this.glDOM.on('mousedown touchstart', this._handleMouseDown);
      this.glDOM.on('wheel', this._handleMouseScroll);
      this.glDOM.on('mousemove touchmove', this._handleMouseMove);

      this.glDOM.on('contextmenu', this._handleContextMenu);
    }
  }

  /**
   * Change the viewer's container element
   * Also useful if the original container element was removed from the DOM.
   *
   * @function GLViewer#setContainer
   *
   * @param {HTMLElement | string} element - Either HTML element or string identifier. Defaults to the element used to initialize the viewer.
   */
  setContainer(element) {
    if (typeof element === 'string') element = $(`#${element}`);
    if (!element) {
      element = this.container;
    }
    this.initContainer(element);
    return this;
  }

  /**
   * Set the background color (default white)
   *
   * @function GLViewer#setBackgroundColor
   * @param {number} hex - Hexcode specified background color, or standard color spec
   * @param {number} a - Alpha level (default 1.0)
   *
   * @example
      viewer.setBackgroundColor(0x000000);
   *
   */
  setBackgroundColor(hex, a) {
    if (typeof a == 'undefined') {
      a = 1.0;
    } else if (a < 0 || a > 1.0) {
      a = 1.0;
    }
    const c = CC.color(hex);
    this.scene.fog.color = c;
    this.bgColor = c.getHex();
    this.renderer.setClearColorHex(c.getHex(), a);
    this.show();

    return this;
  }

  /**
  * Set view projection scheme.  Either orthographic or perspective.
  * Default is perspective.  Orthographic can also be enabled on viewer creation
  * by setting orthographic to true in the this.config object.
  *
  * @param {string} proj
  *
  * @example
    viewer.setViewStyle({style:"outline"});
         $.get('data/1fas.pqr', function(data){
             viewer.addModel(data, "pqr");
             $.get("data/1fas.cube",function(volumedata){
                 viewer.addSurface(SurfaceType.VDW, {opacity:0.85,voldata: new VolumeData(volumedata, "cube"), volscheme: new Gradient.RWB(-10,10)},{});
             });
             viewer.zoomTo();
             viewer.setProjection("orthographic");
             viewer.render(callback);
         });
  *
  */
  setProjection(proj) {
    this.camera.ortho = proj === 'orthographic';
    this.setSlabAndFog();
  }

  /**
   * Set global view styles.
   * @function GLViewer#setViewStyle
   *
   * @example
   *   viewer.setViewStyle({style:"outline"});
        $.get('data/1fas.pqr', function(data){
            viewer.addModel(data, "pqr");
            $.get("data/1fas.cube",function(volumedata){
                viewer.addSurface(SurfaceType.VDW, {opacity:0.85,voldata: new VolumeData(volumedata, "cube"), volscheme: new Gradient.RWB(-10,10)},{});
            });
            viewer.zoomTo();
            viewer.render(callback);
        });
   *
   */
  setViewStyle(parameters) {
    if (parameters.style === 'outline') {
      const params = {};
      if (parameters.color) params.color = CC.color(parameters.color);
      if (parameters.width) params.width = parameters.width;
      this.renderer.enableOutline(params);
    } else {
      this.renderer.disableOutline();
    }
    return this;
  }

  /**
   * Set viewer width
   *
   * @function GLViewer#setWidth
   * @param {number}
   *            w Width in pixels
   */
  setWidth(w) {
    this.renderer.setSize(w, this.HEIGHT);
    return this;
  }

  /**
   * Set viewer height
   *
   * @function GLViewer#setHeight
   * @param {number}
   *            h Height in pixels
   */
  setHeight(h) {
    this.renderer.setSize(this.WIDTH, h);
    return this;
  }

  /**
   * Resize viewer according to containing HTML element's dimensions
   *
   * @function GLViewer#resize
   */
  resize() {
    let regen = false;
    if (this.renderer.isLost() && this.WIDTH > 0 && this.HEIGHT > 0) {
      // create new context
      this.container.children('canvas').remove(); // remove existing
      this.setupRenderer();
      this.initContainer(this.container);
      regen = true;
    }
    this.ASPECT = this.renderer.getAspect(this.WIDTH, this.HEIGHT);
    this.renderer.setSize(this.WIDTH, this.HEIGHT);
    this.camera.aspect = this.ASPECT;
    this.camera.updateProjectionMatrix();

    if (regen) {
      // restored rendere, need to regenerate scene
      const options = this.renderer.supportedExtensions();
      options.regen = true;
      this.render(null, options);
    } else {
      this.show();
    }

    // UI Edition
    this._stateManager.updateUI();
    return this;
  }

  // eslint-disable-next-line class-methods-use-this
  setupRenderer() {
    throw new Error('Method not implemented.');
  }

  /**
   * Return specified model
   *
   * @function GLViewer#getModel
   * @param {number|GLModel} [id] - Retrieve model with specified id
   * @default Returns last model added to viewer or null if there are no models
   * @return {GLModel|null}
   *
   * @example // Retrieve reference to first GLModel added var m =
   *    download("pdb:1UBQ",viewer,{},function(m1){
            download("pdb:1UBI", viewer,{}, function(m2) {
              viewer.zoomTo();
              m1.setStyle({cartoon: {color:'green'}});
              //could use m2 here as well
              viewer.getModel().setStyle({cartoon: {color:'blue'}});
              viewer.render();
          })
        });
   */
  getModel(id) {
    if (id === undefined) {
      return this.models.length == 0 ? null : this.models[this.models.length - 1];
    }
    if (id instanceof GLModel) {
      return id;
    }
    if (!(id in this.models)) {
      if (this.models.length == 0) return null;
      return this.models[this.models.length - 1]; // get last model if no (or invalid) id specified
    }
    return this.models[id];
  }

  spinInterval;
  /**
   * Continuously rotate a scene around the specified axis.
   *
   * Call `GLViewer.spin(false)` to stop spinning.
   *
   * @function GLViewer#spin
   * @param {string| import('./specs').Vector3Like} [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around. Default "y".  View relative (rather than model relative) axes are prefixed with v.
   * @param {number} [speed] - Speed multiplier for spinning the viewer. 1 is default and a negative value reverses the direction of the spin.
   *
   */
  spin(axis, speed) {
    clearInterval(this.spinInterval);
    if (typeof axis == 'undefined') axis = 'y';
    if (typeof axis == 'boolean') {
      if (!axis) return;
      axis = 'y';
    }

    if (Array.isArray(axis)) {
      axis = {x: axis[0], y: axis[1], z: axis[2]};
    }
    // out of bounds check
    const viewer = this;

    this.spinInterval = setInterval(() => {
      if (!viewer.getCanvas().isConnected && this.renderer.isLost()) {
        clearInterval(this.spinInterval);
      }
      viewer.rotate(1 * (speed || 1), axis);
    }, 25);
  }

  // animate motion between current position and passed position
  // can set some parameters to null
  // if fixed is true will enforce the request animation, otherwise
  // does relative updates
  // positions objects have modelggroup position, rotation group position.z,
  // and rotationgroup quaternion
  // return array includes final position, but not current
  // the returned array includes an animate method
  animateMotion(duration, fixed, mpos, rz, rot, cam) {
    const interval = 20;
    let stepLen = Math.ceil(duration / interval);
    if (stepLen < 1) stepLen = 1;
    this.incAnim();

    const curr = {
      mpos: this.modelGroup.position.clone(),
      rz: this.rotationGroup.position.z,
      rot: this.rotationGroup.quaternion.clone(),
      cam: this.lookingAt.clone(),
    };

    if (fixed) {
      // precompute path and stick to it
      let steps = new Array(stepLen);
      const n = steps.length;
      for (let i = 0; i < n; i++) {
        const frac = (i + 1) / n;
        const next = {mpos: curr.mpos, rz: curr.rz, rot: curr.rot};
        if (mpos) {
          next.mpos = mpos.clone().sub(curr.mpos).multiplyScalar(frac).add(curr.mpos);
        }
        if (typeof rz != 'undefined' && rz != null) {
          next.rz = curr.rz + frac * (rz - curr.rz);
        }
        if (rot) {
          next.rot = slerp(curr.rot, rot, frac);
        }
        if (cam) {
          next.cam = cam.clone().sub(curr.cam).multiplyScalar(frac).add(curr.cam);
        }

        steps[i] = next;
      }

      let step = 0;
      const callback = () => {
        const p = steps[step];
        step += 1;
        if (p.mpos) {
          this.modelGroup.position = p.mpos;
        }
        if (p.rz) {
          this.rotationGroup.position.z = p.rz;
        }
        if (p.rot) {
          this.rotationGroup.quaternion = p.rot;
        }
        if (p.cam) {
          this.camera.lookAt(p.cam);
        }

        if (step < steps.length) {
          setTimeout(callback, interval);
        } else {
          this.decAnim();
        }
        this.show();
      };
      setTimeout(callback, interval);
    } else {
      // relative update
      const delta = {};
      const frac = 1.0 / stepLen;
      if (mpos) {
        delta.mpos = mpos.clone().sub(curr.mpos).multiplyScalar(frac);
      }
      if (typeof rz != 'undefined' && rz != null) {
        delta.rz = frac * (rz - curr.rz);
      }
      if (rot) {
        const next = slerp(curr.rot, rot, frac);
        // comptute step delta rotation
        delta.rot = curr.rot.clone().inverse().multiply(next);
      }
      if (cam) {
        delta.cam = cam.clone().sub(curr.cam).multiplyScalar(frac);
      }
      let step = 0.0;
      const callback = () => {
        step += 1;
        if (delta.mpos) {
          this.modelGroup.position.add(delta.mpos);
        }
        if (delta.rz) {
          this.rotationGroup.position.z += delta.rz;
        }
        if (delta.rot) {
          this.rotationGroup.quaternion.multiply(delta.rot);
        }
        if (delta.cam) {
          this.lookingAt.add(delta.cam);
          this.camera.lookAt(this.lookingAt);
        }

        if (step < stepLen) {
          setTimeout(callback, interval);
        } else {
          this.decAnim();
        }
        this.show();
      };
      setTimeout(callback, interval);
    }
  }

  /**
   * Rotate scene by angle degrees around axis
   *
   * @function GLViewer#rotate
   * @param {number} angle - Angle, in degrees, to rotate by.
   * @param {string | import('./specs').Vector3Like} [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around. Default "y".  View relative (rather than model relative) axes are prefixed with v. Axis can also be specified as a vector.
   * @param {number} [animationDuration] - an optional parameter that denotes the duration of the rotation animation. Default 0 (no animation)
   * @param {boolean} [fixedPath] - if true animation is constrained to requested motion, overriding updates that happen during the animation
   * @example     
   *  download('cid:4000', viewer, {}, function() {
        viewer.setStyle({stick:{}});
        viewer.zoomTo();
        viewer.rotate(90,'y',1);
        viewer.render(callback);
      });
   */
  rotate(angle, axis, animationDuration, fixedPath) {
    animationDuration = animationDuration !== undefined ? animationDuration : 0;

    if (typeof axis === 'undefined') {
      axis = 'y';
    }

    let rAxis;

    if (axis == 'x') {
      rAxis = {x: 1, y: 0, z: 0};
    } else if (axis == 'y') {
      rAxis = {x: 0, y: 1, z: 0};
    } else if (axis == 'z') {
      rAxis = {x: 0, y: 0, z: 1};
    }

    // support rotating with respect to view axis, not model
    if (axis == 'vx') {
      rAxis = {vx: 1, vy: 0, vz: 0};
    } else if (axis == 'vy') {
      rAxis = {vx: 0, vy: 1, vz: 0};
    } else if (axis == 'vz') {
      rAxis = {vx: 0, vy: 0, vz: 1};
    }

    if (rAxis && typeof rAxis.vx !== 'undefined') {
      const vaxis = new Vector3(rAxis.vx, rAxis.vy, rAxis.vz);
      vaxis.applyQuaternion(this.rotationGroup.quaternion);
      rAxis = {x: vaxis.x, y: vaxis.y, z: vaxis.z};
    }

    const qFromAngle = function (rangle) {
      const s = Math.sin(rangle / 2.0);
      const c = Math.cos(rangle / 2.0);
      let i = 0;
      let j = 0;
      let k = 0;

      i = rAxis.x * s;
      j = rAxis.y * s;
      k = rAxis.z * s;

      return new Quaternion(i, j, k, c).normalize();
    };

    const rangle = (Math.PI * angle) / 180.0;
    const q = qFromAngle(rangle);

    if (animationDuration) {
      const final = new Quaternion().copy(this.rotationGroup.quaternion).multiply(q); // final
      this.animateMotion(
        animationDuration,
        fixedPath,
        this.modelGroup.position,
        this.rotationGroup.position.z,
        final,
        this.lookingAt
      );
    } else {
      // not animated
      this.rotationGroup.quaternion.multiply(q);
      this.show();
    }
    return this;
  }

  surfacesFinished() {
    for (const key in this.surfaces) {
      if (!this.surfaces[key][0].done) {
        return false;
      }
    }
    return true;
  }

  /** Returns an array representing the current viewpoint.
   * Translation, zoom, and rotation quaternion.
   * @function GLViewer#getView
   * @returns {Array.<number>} [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y, q.z, q.w ]
   *  */
  getView() {
    if (!this.modelGroup) return [0, 0, 0, 0, 0, 0, 0, 1];
    const pos = this.modelGroup.position;
    const q = this.rotationGroup.quaternion;
    return [pos.x, pos.y, pos.z, this.rotationGroup.position.z, q.x, q.y, q.z, q.w];
  }

  /** Sets the view to the specified translation, zoom, and rotation.
   *
   * @function GLViewer#setView
   * @param {Array<number>} [arg] - Array formatted identically to the return value of getView
   * @param {boolean} [nolink]
   * */
  setView(arg, nolink) {
    if (arg === undefined || (Array.isArray(arg) && arg.length >= 8)) return this;

    if (!this.modelGroup || !this.rotationGroup) return this;
    this.modelGroup.position.x = arg[0];
    this.modelGroup.position.y = arg[1];
    this.modelGroup.position.z = arg[2];
    this.rotationGroup.position.z = arg[3];
    this.rotationGroup.quaternion.x = arg[4];
    this.rotationGroup.quaternion.y = arg[5];
    this.rotationGroup.quaternion.z = arg[6];
    this.rotationGroup.quaternion.w = arg[7];
    if (typeof arg[8] != 'undefined') {
      this.rotationGroup.position.x = arg[8];
      this.rotationGroup.position.y = arg[9];
    }

    this.show(nolink);
    return this;
  }

  // apply styles, models, etc in viewer
  /**
   * Render current state of viewer, after
   * adding/removing models, applying styles, etc.
   *
   * @function GLViewer#render
   */
  render(callback, exts) {
    this.renderer.setViewport();
    this.updateClickables(); // must render for clickable styles to take effect
    const view = this.getView();

    if (this.stateChangeCallback) {
      // todo: have ability to only send delta updates
      this.stateChangeCallback(this.getInternalState());
    }

    let i;
    let n;
    if (!exts) exts = this.renderer.supportedExtensions();
    for (i = 0; i < this.models.length; i++) {
      if (this.models[i]) {
        this.models[i].globj(this.modelGroup, exts);
      }
    }

    for (i = 0; i < this.shapes.length; i++) {
      if (this.shapes[i]) {
        // exists
        if (
          typeof this.shapes[i].frame === 'undefined' ||
          this.viewer_frame < 0 ||
          this.shapes[i].frame < 0 ||
          this.shapes[i].frame == this.viewer_frame
        ) {
          this.shapes[i].globj(this.modelGroup, exts);
        } else {
          // should not be displayed in current frame
          this.shapes[i].removegl(this.modelGroup);
        }
      }
    }

    for (i = 0; i < this.labels.length; i++) {
      if (
        this.labels[i] &&
        typeof this.labels[i].frame != 'undefined' &&
        this.labels[i].frame >= 0
      ) {
        // exists and has frame specifier
        this.modelGroup.remove(this.labels[i].sprite);
        if (this.viewer_frame < 0 || this.labels[i].frame == this.viewer_frame) {
          this.modelGroup.add(this.labels[i].sprite);
        }
      }
    }

    for (i in this.surfaces) {
      // this is an object with possible holes
      if (!this.surfaces.hasOwnProperty(i)) continue;
      const surfArr = this.surfaces[i];
      for (n = 0; n < surfArr.length; n++) {
        if (surfArr.hasOwnProperty(n)) {
          const {geo} = surfArr[n];
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

            if (surfArr[n].done) surfArr[n].finished = true;

            // remove partially rendered surface
            if (surfArr[n].lastGL) this.modelGroup.remove(surfArr[n].lastGL);

            // create new surface
            let smesh = null;

            if (surfArr[n].mat instanceof LineBasicMaterial) {
              // special case line meshes
              smesh = new Line(geo, surfArr[n].mat);
            } else {
              smesh = new Mesh(geo, surfArr[n].mat);
            }
            if (surfArr[n].mat.transparent && surfArr[n].mat.opacity == 0) {
              // don't bother with hidden surfaces
              smesh.visible = false;
            } else {
              smesh.visible = true;
            }
            if (
              surfArr[n].symmetries.length > 1 ||
              (surfArr[n].symmetries.length == 1 && !surfArr[n].symmetries[n].isIdentity())
            ) {
              const tmeshes = new Object3D(); // transformed meshes
              for (let j = 0; j < surfArr[n].symmetries.length; j++) {
                const tmesh = smesh.clone();
                tmesh.matrix = surfArr[n].symmetries[j];
                tmesh.matrixAutoUpdate = false;
                tmeshes.add(tmesh);
              }
              surfArr[n].lastGL = tmeshes;
              this.modelGroup.add(tmeshes);
            } else {
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
  }

  /** @param {import('./specs').AtomSelectionSpec} sel
   * @return list of models specified by sel
   */
  getModelList(sel) {
    /** @type {GLModel|GLModel[]} */
    let ms = [];
    if (typeof sel === 'undefined' || typeof sel.model === 'undefined') {
      for (let i = 0; i < this.models.length; i++) {
        if (this.models[i]) ms.push(this.models[i]);
      }
    } else {
      // specific to some models
      ms = sel.model;
      if (!Array.isArray(ms)) ms = [ms];

      for (let i = 0; i < ms.length; i++) {
        // allow referencing models by order of creation
        if (typeof ms[i] === 'number') {
          let index = ms[i];
          // support python backward indexing
          if (index < 0) index += this.models.length;
          ms[i] = this.models[index];
        }
      }
    }

    return ms;
  }

  /**
   *
   * @param {import('./specs').AtomSelectionSpec}
   *            sel
   * @return {import('./specs').AtomSpec[]}
   */
  getAtomsFromSel(sel) {
    let atoms = [];
    if (typeof sel === 'undefined') sel = {};

    const ms = this.getModelList(sel);

    for (let i = 0; i < ms.length; i++) {
      atoms = atoms.concat(ms[i].selectedAtoms(sel));
    }

    return atoms;
  }

  /**
   *
   * @param {import('./specs').AtomSpec} atom
   * @param {import('./specs').AtomSelectionSpec} sel
   * @return {boolean}
   */
  atomIsSelected(atom, sel) {
    if (typeof sel === 'undefined') sel = {};

    const ms = this.getModelList(sel);

    for (let i = 0; i < ms.length; i++) {
      if (ms[i].atomIsSelected(atom, sel)) return true;
    }

    return false;
  }

  /** return list of atoms selected by sel
   *
   * @function GLViewer#selectedAtoms
   * @param {import('./specs').AtomSelectionSpec} sel
   * @return {Array.<Object>}
   */
  selectedAtoms(sel) {
    return this.getAtomsFromSel(sel);
  }

  /**
   * Returns valid values for the specified attribute in the given selection
   * @function GlViewer#getUniqueValues
   * @param {string} attribute
   * @param {import('./specs').AtomSelectionSpec} sel
   * @return {Array.<Object>}
   *
   */
  getUniqueValues(attribute, sel) {
    if (typeof sel === 'undefined') sel = {};
    const atoms = this.getAtomsFromSel(sel);
    const values = {};

    for (const atom in atoms) {
      if (atoms[atom].hasOwnProperty(attribute)) {
        const value = atoms[atom][attribute];
        values[value] = true;
      }
    }

    return Object.keys(values);
  }

  /**
   * Return pdb output of selected atoms (if atoms from pdb input)
   *
   * @function GLViewer#pdbData
   * @param {Object=} [sel] - Selection specification specifying model and atom properties to select.  Default: all atoms in viewer
   * @return {string} PDB string of selected atoms
   */
  pdbData(sel) {
    const atoms = this.getAtomsFromSel(sel);
    let ret = '';
    for (let i = 0, n = atoms.length; i < n; ++i) {
      ret += `${atoms[i].pdbline}\n`;
    }
    return ret;
  }

  /**
   * Zoom current view by a constant factor
   *
   * @function GLViewer#zoom
   * @param {number}
   *            [factor] - Magnification factor. Values greater than 1
   *            will zoom in, less than one will zoom out. Default 2.
   * @param {number}
   *            [animationDuration] - an optional parameter that denotes
   *            the duration of a zoom animation
   * @param {Boolean} [fixedPath] - if true animation is constrained to
   *      requested motion, overriding updates that happen during the animation
   * @example
    $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo()
      viewer.zoom(2,1000);
      viewer.render();
    });
              
   */
  zoom(factor, animationDuration, fixedPath) {
    factor = factor || 2;
    animationDuration = animationDuration !== undefined ? animationDuration : 0;
    const scale = (this.CAMERA_Z - this.rotationGroup.position.z) / factor;
    const finalZ = this.CAMERA_Z - scale;

    if (animationDuration > 0) {
      this.animateMotion(
        animationDuration,
        fixedPath,
        this.modelGroup.position,
        this.adjustZoomToLimits(finalZ),
        this.rotationGroup.quaternion,
        this.lookingAt
      );
    } else {
      // no animation
      this.rotationGroup.position.z = this.adjustZoomToLimits(finalZ);
      this.show();
    }
    return this;
  }

  /**
   * Translate current view by x,y screen coordinates
   * This pans the camera rather than translating the model.
   *
   * @function GLViewer#translate
   * @param {number} x Relative change in view coordinates of camera
   * @param {number} y Relative change in view coordinates of camera
   * @param {number}
   *            [animationDuration] - an optional parameter that denotes
   *            the duration of a zoom animation
   * @param {Boolean} [fixedPath] - if true animation is constrained to
   *      requested motion, overriding updates that happen during the animation         *
   * @example     
      $.get('data/4csv.pdb', function(data) {
        viewer.addModel(data,'pdb');
        viewer.setStyle({cartoon:{},stick:{}});
        viewer.zoomTo();
        viewer.translate(200,50);
        viewer.rotate(90,'z');
        viewer.render(callback);
      });
   */
  translate(x, y, animationDuration, fixedPath) {
    animationDuration = animationDuration !== undefined ? animationDuration : 0;
    const dx = x / this.WIDTH;
    const dy = y / this.HEIGHT;
    const v = new Vector3(0, 0, -this.CAMERA_Z);

    this.projector.projectVector(v, this.camera);
    v.x -= dx;
    v.y -= dy;
    this.projector.unprojectVector(v, this.camera);
    v.z = 0;

    const finalPosition = this.lookingAt.clone().add(v);
    if (animationDuration > 0) {
      this.animateMotion(
        animationDuration,
        fixedPath,
        this.modelGroup.position,
        this.rotationGroup.position.z,
        this.rotationGroup.quaternion,
        finalPosition
      );
    } else {
      // no animation
      this.lookingAt = finalPosition;
      this.camera.lookAt(this.lookingAt);
      this.show();
    }
    return this;
  }

  /**
   * Translate current models by x,y screen coordinates
   * This translates the models relative to the current view. It does
   * not change the center of rotation.
   *
   * @function GLViewer#translateScene
   * @param {number} x Relative change in x screen coordinate
   * @param {number} y Relative change in y screen coordinate
   * @param {number}
   *            [animationDuration] - an optional parameter that denotes
   *            the duration of a zoom animation
   * @param {Boolean} [fixedPath] - if true animation is constrained to
   *      requested motion, overriding updates that happen during the animation         *
   * @example     
      $.get('data/4csv.pdb', function(data) {
        viewer.addModel(data,'pdb');
        viewer.setStyle({cartoon:{},stick:{}});
        viewer.zoomTo();
        viewer.translateScene(200,50);
        viewer.rotate(90,'z'); // will no longer be around model center
        viewer.render(callback);
      });
   */
  translateScene(x, y, animationDuration, fixedPath) {
    animationDuration = animationDuration !== undefined ? animationDuration : 0;

    const t = this.screenOffsetToModel(x, y);
    const finalPosition = this.modelGroup.position.clone().add(t);

    if (animationDuration > 0) {
      this.animateMotion(
        animationDuration,
        fixedPath,
        this.modelGroup.position,
        this.rotationGroup.position.z,
        this.rotationGroup.quaternion,
        this.lookingAt
      );
    } else {
      // no animation
      this.modelGroup.position = finalPosition;
      this.show();
    }
    return this;
  }

  /**
   * Adjust slab to fully enclose selection (default everything).
   *
   * @function GLViewer#fitSlab
   * @param {Object}
   *            [sel] - Selection specification specifying model and atom
   *            properties to select. Default: all atoms in viewer
   */
  fitSlab(sel) {
    sel = sel || {};
    const atoms = this.getAtomsFromSel(sel);
    const tmp = getExtent(atoms);

    // fit to bounding box
    const x = tmp[1][0] - tmp[0][0];
    const y = tmp[1][1] - tmp[0][1];
    const z = tmp[1][2] - tmp[0][2];

    let maxD = Math.sqrt(x * x + y * y + z * z);
    if (maxD < 5) maxD = 5;

    // use full bounding box for slab/fog
    this.slabNear = -maxD / 1.9;
    this.slabFar = maxD / 2;

    return this;
  }

  /**
   * Re-center the viewer around the provided selection (unlike zoomTo, does not zoom).
   *
   * @function GLViewer#center
   * @param {Object}
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
      $.get('data/4csv.pdb', function(data) {
        viewer.addModel(data,'pdb');
        viewer.setStyle({cartoon:{},stick:{}});
        viewer.center();
        viewer.render(callback);
      });
   */
  center(sel, animationDuration, fixedPath) {
    animationDuration = animationDuration !== undefined ? animationDuration : 0;
    let allatoms;
    let alltmp;
    sel = sel || {};
    /** @type {Array<import('./specs').AtomSpec|Vector3>} */
    const atoms = this.getAtomsFromSel(sel);
    let tmp = getExtent(/** @type {import('./specs').AtomSpec[]} */ (atoms));

    if (isEmptyObject(sel)) {
      // include shapes when zooming to full scene
      // TODO: figure out a good way to specify shapes as part of a selection
      this.shapes.forEach(shape => {
        if (shape && shape.boundingSphere && shape.boundingSphere.center) {
          const c = shape.boundingSphere.center;
          const r = shape.boundingSphere.radius;
          if (r > 0) {
            // make sure full shape is visible
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
      tmp = getExtent(/** @type {import('./specs').AtomSpec[]} */ (atoms));
      allatoms = atoms;
      alltmp = tmp;
    } else {
      allatoms = this.getAtomsFromSel({});
      alltmp = getExtent(allatoms);
    }

    // use selection for center
    const center = new Vector3(tmp[2][0], tmp[2][1], tmp[2][2]);

    // but all for bounding box
    let x = alltmp[1][0] - alltmp[0][0];
    let y = alltmp[1][1] - alltmp[0][1];
    let z = alltmp[1][2] - alltmp[0][2];

    let maxD = Math.sqrt(x * x + y * y + z * z);
    if (maxD < 5) maxD = 5;

    // use full bounding box for slab/fog
    this.slabNear = -maxD / 1.9;
    this.slabFar = maxD / 2;

    // for zoom, use selection box
    x = tmp[1][0] - tmp[0][0];
    y = tmp[1][1] - tmp[0][1];
    z = tmp[1][2] - tmp[0][2];
    maxD = Math.sqrt(x * x + y * y + z * z);
    if (maxD < 5) maxD = 5;

    // find the farthest atom from center to get max distance needed for view
    let maxDsq = 25;
    for (let i = 0; i < atoms.length; i++) {
      if (atoms[i]) {
        const dsq = center.distanceToSquared(atoms[i]);
        if (dsq > maxDsq) maxDsq = dsq;
      }
    }

    maxD = Math.sqrt(maxDsq) * 2;
    const finalpos = center.clone().multiplyScalar(-1);
    if (animationDuration > 0) {
      this.animateMotion(
        animationDuration,
        fixedPath,
        finalpos,
        this.rotationGroup.position.z,
        this.rotationGroup.quaternion,
        this.lookingAt
      );
    } else {
      // no animation
      this.modelGroup.position = finalpos;
      this.show();
    }
    return this;
  }

  /**
   * Zoom to center of atom selection.  The slab will be set appropriately for
   * the selection, unless an empty selection is provided, in which case there will be no slab.
   *
   * @function GLViewer#zoomTo
   * @param {Object}
   *            [sel] - Selection specification specifying model and atom
   *            properties to select. Default: all atoms in viewer
   * @param {number}
   *            [animationDuration] - an optional parameter that denotes
   *            the duration of a zoom animation
   * @param {Boolean} [fixedPath] - if true animation is constrained to
   *      requested motion, overriding updates that happen during the animation         *
   * @example
        $.get('data/1fas.pqr', function(data){
            viewer.addModel(data, "pqr");
            $.get("data/1fas.cube",function(volumedata){
                viewer.addSurface(SurfaceType.VDW, {
                    opacity:0.85,
                    voldata: new VolumeData(volumedata, "cube"),
                    volscheme: new Gradient.Sinebow(getPropertyRange(viewer.selectedAtoms(),'charge'))
                },{});
            viewer.render();
            });
            viewer.zoomTo();
          });
   */
  zoomTo(sel, animationDuration, fixedPath) {
    animationDuration = animationDuration !== undefined ? animationDuration : 0;
    sel = sel || {};
    /** @type {Array<import('./specs').AtomSpec|Vector3>} */
    const atoms = this.getAtomsFromSel(sel);
    const atombox = getExtent(/** @type {import('./specs').AtomSpec[]} */ (atoms));
    let allbox = atombox;

    if (isEmptyObject(sel)) {
      // include shapes when zooming to full scene
      // TODO: figure out a good way to specify shapes as part of a selection
      const natoms = atoms && atoms.length;
      this.shapes.forEach(shape => {
        if (shape && shape.boundingSphere) {
          if (shape.boundingSphere.box) {
            const {box} = shape.boundingSphere;
            atoms.push(new Vector3(box.min.x, box.min.y, box.min.z));
            atoms.push(new Vector3(box.max.x, box.max.y, box.max.z));
          } else if (shape.boundingSphere.center) {
            const c = shape.boundingSphere.center;
            const r = shape.boundingSphere.radius;
            if (r > 0) {
              // make sure full shape is visible
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
      allbox = getExtent(/** @type {import('./specs').AtomSpec[]} */ (atoms));
      if (!natoms) {
        // if no atoms, use shapes for center
        for (let i = 0; i < 3; i++) {
          // center of bounding box
          atombox[2][i] = (allbox[0][i] + allbox[1][i]) / 2;
        }
      }
    } else {
      // include all atoms in slab calculation
      const allatoms = this.getAtomsFromSel({});
      allbox = getExtent(allatoms);
    }

    // use selection for center
    const center = new Vector3(atombox[2][0], atombox[2][1], atombox[2][2]);

    // but all for bounding box
    let x = allbox[1][0] - allbox[0][0];
    let y = allbox[1][1] - allbox[0][1];
    let z = allbox[1][2] - allbox[0][2];

    let maxD = Math.sqrt(x * x + y * y + z * z);
    if (maxD < 5) maxD = 5;

    // use full bounding box for slab/fog
    this.slabNear = -maxD / 1.9;
    this.slabFar = maxD / 2;

    // if we are selecting everything, have ver permissive slab
    // can't do "infinity" size since this will break orthographic
    if (Object.keys(sel).length === 0) {
      this.slabNear = Math.min(-maxD * 2, -50);
      this.slabFar = Math.max(maxD * 2, 50);
    }

    // keep at least this much space in view
    const MAXD = this.config.minimumZoomToDistance || 5;
    // for zoom, use selection box
    x = atombox[1][0] - atombox[0][0];
    y = atombox[1][1] - atombox[0][1];
    z = atombox[1][2] - atombox[0][2];
    maxD = Math.sqrt(x * x + y * y + z * z);
    if (maxD < MAXD) maxD = MAXD;

    // find the farthest atom from center to get max distance needed for view
    let maxDsq = MAXD * MAXD;
    for (let i = 0; i < atoms.length; i++) {
      if (atoms[i]) {
        const dsq = center.distanceToSquared(atoms[i]);
        if (dsq > maxDsq) maxDsq = dsq;
      }
    }

    maxD = Math.sqrt(maxDsq) * 2;
    const finalpos = center.clone().multiplyScalar(-1);
    let finalz = -(
      (maxD * 0.5) / Math.tan(((Math.PI / 180.0) * this.camera.fov) / 2) -
      this.CAMERA_Z
    );

    finalz = this.adjustZoomToLimits(finalz);
    if (animationDuration > 0) {
      this.animateMotion(
        animationDuration,
        fixedPath,
        finalpos,
        finalz,
        this.rotationGroup.quaternion,
        this.lookingAt
      );
    } else {
      this.modelGroup.position = finalpos;
      this.rotationGroup.position.z = finalz;
      this.show();
    }
    return this;
  }

  /**
   * Set slab of view (contents outside of slab are clipped).
   * Must call render to update.
   *
   * @function GLViewer#setSlab
   * @param {number} near near clipping plane distance
   * @param {number} far far clipping plane distance
   */
  setSlab(near, far) {
    this.slabNear = near;
    this.slabFar = far;
  }

  /**
   * Get slab of view (contents outside of slab are clipped).
   *
   * @function GLViewer#getSlab
   * @return {Object}
   *      @property {number} near - near clipping plane distance
   *      @property {number} far - far clipping plane distance
   */
  getSlab() {
    return {near: this.slabNear, far: this.slabFar};
  }

  /**
   * Add label to viewer
   *
   * @function GLViewer#addLabel
   * @param {string} text - Label text
   * @param {import('./specs').LabelSpec} options - Label style specification
   * @param {import('./specs').AtomSelectionSpec} [sel] - Set position of label to center of this selection
   * @param {boolean} [noshow] - if true, do not immediately display label - when adding multiple labels this is more efficient
   * @return {import("./Label").default}
   * @example
      download("pdb:2EJ0",viewer,{},function(){
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
  addLabel(text, options, sel, noshow) {
    options = options || {};
    if (sel) {
      const extent = getExtent(this.getAtomsFromSel(sel));
      options.position = {x: extent[2][0], y: extent[2][1], z: extent[2][2]};
    }
    const label = new Label(text, options);
    label.setContext();
    this.modelGroup.add(label.sprite);
    if (options.fixed) this.fixed_labels.push(this.labels.length);
    this.labels.push(label);

    if (!noshow) this.show();
    return label;
  }

  /** Add residue labels.  This will generate one label per a
   * residue within the selected atoms.  The label will be at the
   * centroid of the atoms and styled according to the passed style.
   * The label text will be [resn][resi]
   *
   * @function GLViewer#addResLabels
   * @param {Object} sel
   * @param {Object} style
   * @param {boolean} byframe - if true, create labels for every individual frame, not just current
   *
   * * @example
       download("mmtf:2ll5",viewer,{},function(){
            viewer.setStyle({stick:{radius:0.15},cartoon:{}});
            viewer.addResLabels({hetflag:false}, {font: 'Arial', fontColor:'black',showBackground:false, screenOffset: {x:0,y:0}});
            viewer.zoomTo();
            viewer.render();
          });
   */
  addResLabels(sel, style, byframe) {
    const start = this.labels.length;
    this.applyToModels('addResLabels', sel, this, style, byframe);
    this.show();
    return this.labels.slice(start);
  }

  /** Add property labels.  This will generate one label per a selected
   * atom at the atom's coordinates with the property value as the label text.
   *
   * @function GLViewer#addPropertyLabels
   * @param {string} prop - property name
   * @param {Object} sel
   * @param {Object} style
   *
   * * @example
       download("cid:5291",viewer,{},function(){
            viewer.setStyle({stick: {radius:.2}});
            viewer.addPropertyLabels("index",{not:{elem:'H'}}, {fontColor:'black',font: 'sans-serif', fontSize: 28, showBackground:false,alignment:'center'});
            viewer.zoomTo();
            viewer.render();
          });
   */
  addPropertyLabels(prop, sel, style) {
    this.applyToModels('addPropertyLabels', prop, sel, this, style);
    this.show();
    return this;
  }

  /**
   * Remove label from viewer
   *
   * @function GLViewer#removeLabel
   * @param {Label}
   *            label - $3Dmol label
   *
   * @example // Remove labels created in
   download("pdb:2EJ0",viewer,{},function(){
            var toremove = viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
            viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                     borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                     position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
            viewer.removeLabel(toremove);
            viewer.render();
          });
   */
  removeLabel(label) {
    // todo: don't do the linear search
    for (let i = 0; i < this.labels.length; i++) {
      if (this.labels[i] == label) {
        this.labels.splice(i, 1);
        label.dispose();
        this.modelGroup.remove(label.sprite);
        break;
      }
    }
    this.show();
    return this;
  }

  /**
   * Remove all labels from viewer
   *
   * @function GLViewer#removeAllLabels
   *         @example
  download("pdb:1ubq",viewer,{},function(){
         viewer.addResLabels();
         viewer.setStyle({},{stick:{}});
         viewer.render( ); //show labels
         viewer.removeAllLabels();
         viewer.render(); //hide labels
  });
   */
  removeAllLabels() {
    for (let i = 0; i < this.labels.length; i++) {
      if (this.labels[i] && this.labels[i].sprite) {
        this.modelGroup.remove(this.labels[i].sprite);
      }
    }
    this.labels.splice(0, this.labels.length); // don't overwrite in case linked
    this.show();
    return this;
  }

  // Modify label style
  /**
   * Modify existing label's style
   *
   * @function GLViewer#setLabelStyle
   * @param {Label}
   *            label - $3Dmol label
   * @param {Object}
   *            stylespec - Label style specification
   * @return {Label}
   */
  setLabelStyle(label, stylespec) {
    this.modelGroup.remove(label.sprite);
    label.dispose();
    label.stylespec = stylespec;
    label.setContext();
    this.modelGroup.add(label.sprite);
    this.show();
    return label;
  }

  // Change label text
  /**
   * Modify existing label's text
   *
   * @function GLViewer#setLabelText
   * @param {Label}
   *            label - $3Dmol label
   * @param {String}
   *            text - Label text
   * @return {Label}
   */
  setLabelText(label, text) {
    this.modelGroup.remove(label.sprite);
    label.dispose();
    label.text = text;
    label.setContext();
    this.modelGroup.add(label.sprite);
    this.show();
    return label;
  }

  /**
   * Add shape object to viewer
   * @see {@link GLShape}
   *
   * @function GLViewer#addShape
   * @param {import('./specs').ShapeSpec} shapeSpec - style specification for label
   * @return {GLShape}
   */
  addShape(shapeSpec) {
    shapeSpec = shapeSpec || {};
    const shape = new GLShape(shapeSpec);
    shape.shapePosition = this.shapes.length;
    this.shapes.push(shape);

    return shape;
  }

  /**
   * Remove shape object from viewer
   *
   * @function GLViewer#removeShape
   * @param {GLShape} shape - Reference to shape object to remove
   */
  removeShape(shape) {
    if (!shape) return this;
    shape.removegl(this.modelGroup);
    delete this.shapes[shape.shapePosition];
    // clear off back of model array
    while (this.shapes.length > 0 && typeof this.shapes[this.shapes.length - 1] === 'undefined')
      this.shapes.pop();
    return this;
  }

  /**
   * Remove all shape objects from viewer
   * @function GLViewer#removeAllShapes
   */
  removeAllShapes() {
    for (let i = 0; i < this.shapes.length; i++) {
      const shape = this.shapes[i];
      if (shape) shape.removegl(this.modelGroup);
    }
    this.shapes.splice(0, this.shapes.length);
    return this;
  }

  // gets the center of the selection
  getSelectionCenter(spec) {
    if (spec.hasOwnProperty('x') && spec.hasOwnProperty('y') && spec.hasOwnProperty('z'))
      return spec;
    const atoms = this.getAtomsFromSel(spec);
    if (atoms.length == 0) return {x: 0, y: 0, z: 0};

    const extent = getExtent(atoms);
    return {
      x: extent[0][0] + (extent[1][0] - extent[0][0]) / 2,
      y: extent[0][1] + (extent[1][1] - extent[0][1]) / 2,
      z: extent[0][2] + (extent[1][2] - extent[0][2]) / 2,
    };
  }

  /**
   * Create and add sphere shape. This method provides a shorthand
   * way to create a spherical shape object
   *
   * @function GLViewer#addSphere
   * @param {import('./specs').SphereShapeSpec} spec - Sphere shape style specification
   * @return {GLShape}
   @example
   viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
   viewer.render();
   */
  addSphere(spec) {
    spec = spec || {};

    spec.center = this.getSelectionCenter(spec.center);

    const s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    s.addSphere(spec);
    this.shapes.push(s);
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended
    return s;
  }

  /**
   * Create and add box shape. This method provides a shorthand
   * way to create a box shape object
   *
   * @function GLViewer#addBox
   * @param {import('./specs').BoxSpec} spec - Box shape style specification
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
  addBox(spec) {
    spec = spec || {};

    if (spec.corner != undefined) {
      spec.corner = this.getSelectionCenter(spec.corner);
    }
    if (spec.center != undefined) {
      spec.center = this.getSelectionCenter(spec.center);
    }

    const s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    s.addBox(spec);
    this.shapes.push(s);
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended

    return s;
  }

  /**
   * Create and add arrow shape
   *
   * @function GLViewer#addArrow
   * @param {import('./specs').ArrowSpec} spec - Style specification
   * @return {GLShape}
   @example
    download("pdb:4DM7",viewer,{},function(){
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
  addArrow(spec) {
    spec = spec || {};

    spec.start = this.getSelectionCenter(spec.start);
    spec.end = this.getSelectionCenter(spec.end);

    const s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    s.addArrow(spec);
    this.shapes.push(s);
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended

    return s;
  }

  /**
   * Create and add cylinder shape
   *
   * @function GLViewer#addCylinder
   * @param {import('./specs').CylinderSpec} spec - Style specification
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
  addCylinder(spec) {
    spec = spec || {};

    spec.start = this.getSelectionCenter(spec.start);
    spec.end = this.getSelectionCenter(spec.end);

    const s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    if (spec.dashed) s.addDashedCylinder(spec);
    else s.addCylinder(spec);
    this.shapes.push(s);
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended

    return s;
  }

  /**
   * Create and add Curve shape
   *
   * @function GLViewer#addCurve
   * @param {import('./specs').CurveSpec} spec - Style specification
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
  addCurve(spec) {
    spec = spec || {};
    const s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    s.addCurve(spec);
    this.shapes.push(s);
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended

    return s;
  }

  /**
   * Create and add line shape
   *
   * @function GLViewer#addLine
   * @param {import('./specs').LineSpec} spec - Style specification, can specify dashed, dashLength, and gapLength
   * @return {GLShape}
   @example
   download("pdb:2ABJ",viewer,{},function(){
            viewer.setViewStyle({style:"outline"});
            viewer.setStyle({chain:'A'},{sphere:{hidden:true}});
            viewer.setStyle({chain:'D'},{sphere:{radius:3.0}});
            viewer.setStyle({chain:'G'},{sphere:{colorscheme:'greenCarbon'}});
            viewer.setStyle({chain:'J'},{sphere:{color:'blue'}});
            viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
            viewer.render();
        });
   */
  addLine(spec) {
    spec = spec || {};

    spec.start = this.getSelectionCenter(spec.start);
    spec.end = this.getSelectionCenter(spec.end);

    spec.wireframe = true;
    let s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    if (spec.dashed) s = GLViewer.addLineDashed(spec, s);
    else s.addLine(spec);
    this.shapes.push(s);
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended

    return s;
  }

  /**
   * Create and add unit cell visualization.
   *
   * @function GLViewer#addUnitCell
   * @param {GLModel} modelIn - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
   * @param {import('./specs').UnitCellStyleSpec} [specIn] - visualization style
     @example
          $.get('data/1jpy.cif', function(data) {
            let m = viewer.addModel(data);
            viewer.addUnitCell(m, {box:{color:'purple'},alabel:'X',blabel:'Y',clabel:'Z',alabelstyle: {fontColor: 'black',backgroundColor:'white',inFront:true,fontSize:40},astyle:{color:'darkred', radius:5,midpos: -10}});
            viewer.zoomTo();
            viewer.render();
          });
   */
  addUnitCell(modelIn, specIn) {
    const model = this.getModel(modelIn);
    if (!model) return;
    const spec = specIn || {alabel: 'a', blabel: 'b', clabel: 'c'};

    spec.box = spec.box || {};
    spec.astyle = spec.astyle || {color: 'red', radius: 0.1, midpos: -1};
    spec.bstyle = spec.bstyle || {color: 'green', radius: 0.1, midpos: -1};
    spec.cstyle = spec.cstyle || {color: 'blue', radius: 0.1, midpos: -1};
    spec.alabelstyle = spec.alabelstyle || {
      fontColor: 'red',
      showBackground: false,
      alignment: 'center',
      inFront: false,
    };
    spec.blabelstyle = spec.blabelstyle || {
      fontColor: 'green',
      showBackground: false,
      alignment: 'center',
      inFront: false,
    };
    spec.clabelstyle = spec.clabelstyle || {
      fontColor: 'blue',
      showBackground: false,
      alignment: 'center',
      inFront: false,
    };

    // clear any previous box
    if (model.unitCellObjects) {
      this.removeUnitCell(model);
    }
    model.unitCellObjects = {shapes: [], labels: []};
    // calculate points
    const data = model.getCrystData();
    let matrix = null;
    if (data) {
      if (data.matrix) {
        matrix = data.matrix;
      } else {
        const {a} = data;
        const {b} = data;
        const {c} = data;
        let {alpha} = data;
        let {beta} = data;
        let {gamma} = data;
        alpha = (alpha * Math.PI) / 180.0;
        beta = (beta * Math.PI) / 180.0;
        gamma = (gamma * Math.PI) / 180.0;

        let u;
        let v;
        let w;

        u = Math.cos(beta);
        v = (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma)) / Math.sin(gamma);
        w = Math.sqrt(Math.max(0, 1 - u * u - v * v));

        matrix = new Matrix3(
          a,
          b * Math.cos(gamma),
          c * u,
          0,
          b * Math.sin(gamma),
          c * v,
          0,
          0,
          c * w
        );
      }

      const points = [
        new Vector3(0, 0, 0),
        new Vector3(1, 0, 0),
        new Vector3(0, 1, 0),
        new Vector3(0, 0, 1),
        new Vector3(1, 1, 0),
        new Vector3(0, 1, 1),
        new Vector3(1, 0, 1),
        new Vector3(1, 1, 1),
      ];

      // console.log('Matrix4', data.matrix4, data.matrix);
      if (data.matrix4) {
        for (let i = 0; i < points.length; i++) {
          if (data.size) points[i].multiplyVectors(points[i], data.size); // matrix is for unit vectors, not whole box
          points[i] = points[i].applyMatrix4(data.matrix4);
        }
      } else {
        for (let i = 0; i < points.length; i++) {
          points[i] = points[i].applyMatrix3(matrix);
        }
      }

      // draw box
      if (spec.box && !spec.box.hidden) {
        spec.box.wireframe = true;
        const s = new GLShape(spec.box);
        s.shapePosition = this.shapes.length;

        s.addLine({start: points[0], end: points[1]});
        s.addLine({start: points[0], end: points[2]});
        s.addLine({start: points[1], end: points[4]});
        s.addLine({start: points[2], end: points[4]});

        s.addLine({start: points[0], end: points[3]});
        s.addLine({start: points[3], end: points[5]});
        s.addLine({start: points[2], end: points[5]});

        s.addLine({start: points[1], end: points[6]});
        s.addLine({start: points[4], end: points[7]});
        s.addLine({start: points[6], end: points[7]});

        s.addLine({start: points[3], end: points[6]});
        s.addLine({start: points[5], end: points[7]});

        this.shapes.push(s);
        model.unitCellObjects.shapes.push(s);
        s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended
      }

      // draw arrows
      if (!spec.astyle.hidden) {
        spec.astyle.start = points[0];
        spec.astyle.end = points[1];
        const arrow = this.addArrow(spec.astyle);
        model.unitCellObjects.shapes.push(arrow);
      }

      if (!spec.bstyle.hidden) {
        spec.bstyle.start = points[0];
        spec.bstyle.end = points[2];
        const arrow = this.addArrow(spec.bstyle);
        model.unitCellObjects.shapes.push(arrow);
      }

      if (!spec.cstyle.hidden) {
        spec.cstyle.start = points[0];
        spec.cstyle.end = points[3];
        const arrow = this.addArrow(spec.cstyle);
        model.unitCellObjects.shapes.push(arrow);
      }

      if (spec.alabel) {
        spec.alabelstyle.position = points[1];
        const label = this.addLabel(spec.alabel, spec.alabelstyle);
        model.unitCellObjects.labels.push(label);
      }
      if (spec.blabel) {
        spec.blabelstyle.position = points[2];
        const label = this.addLabel(spec.blabel, spec.blabelstyle);
        model.unitCellObjects.labels.push(label);
      }
      if (spec.clabel) {
        spec.clabelstyle.position = points[3];
        const label = this.addLabel(spec.clabel, spec.clabelstyle);
        model.unitCellObjects.labels.push(label);
      }
    }
  }

  /**
   * Remove unit cell visualization from model.
   *
   * @function GLViewer#removeUnitCell
   * @param {GLModel} modelIn - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
   * @returns {boolean} - True if unit cell was removed. false if model was not found.
     @example
          $.get('data/icsd_200866.cif', function(data) {
            let m = viewer.addModel(data);
            viewer.setStyle({sphere:{}})
            viewer.addUnitCell();
            viewer.zoomTo();
            viewer.removeUnitCell();
            viewer.render();
      });
   */
  removeUnitCell(modelIn) {
    const model = this.getModel(modelIn);
    if (!model) return false;
    if (model.unitCellObjects && typeof model.unitCellObjects === 'object') {
      const viewer = this;
      model.unitCellObjects.shapes.forEach(s => {
        viewer.removeShape(s);
      });
      model.unitCellObjects.labels.forEach(l => {
        viewer.removeLabel(l);
      });
    }
    delete model.unitCellObjects;
    return true;
  }

  /**
   * Replicate atoms in model to form a super cell of the specified dimensions.
   * Original cell will be centered as much as possible.
   *
   * @function GLViewer#replicateUnitCell
   * @param {number} A - number of times to replicate cell in X dimension.
   * @param {number} B - number of times to replicate cell in Y dimension.  If absent, X value is used.
   * @param {number} C - number of times to replicate cell in Z dimension.  If absent, Y value is used.
   * @param {GLModel} modelIn - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
   * @returns {boolean} - True if unit cell was replicated. false if model was not found.
     @example
          $.get('data/icsd_200866.cif', function(data) {
            let m = viewer.addModel(data);
            viewer.setStyle({sphere:{scale:.25}})
            viewer.addUnitCell();
            viewer.zoomTo();
            viewer.replicateUnitCell(3,2,1,m);
            viewer.render();
      });
   */
  replicateUnitCell(A, B, C, modelIn) {
    const model = this.getModel(modelIn);
    if (!model) return false;
    A = A || 3;
    B = B || A;
    C = C || B;
    const cryst = model.getCrystData();
    if (cryst) {
      const atoms = model.selectedAtoms({});
      const {matrix} = cryst;
      const makeoff = function (I) {
        // alternate around zero: 1,-1,2,-2...
        if (I % 2 == 0) return -I / 2;
        return Math.ceil(I / 2);
      };

      for (let i = 0; i < A; i++) {
        for (let j = 0; j < B; j++) {
          for (let k = 0; k < C; k++) {
            if (i == 0 && j == 0 && k == 0) continue; // actual unit cell
            const offset = new Vector3(makeoff(i), makeoff(j), makeoff(k));
            offset.applyMatrix3(matrix);

            const newatoms = [];
            for (let a = 0; a < atoms.length; a++) {
              /** @type {Partial<import('./specs').AtomSpec>} */
              const newAtom = {};
              for (const p in atoms[a]) {
                newAtom[p] = atoms[a][p];
              }
              newAtom.x = (offset.x || 0) + (newAtom.x || 0);
              newAtom.y = (offset.y || 0) + (newAtom.y || 0);
              newAtom.z = (offset.z || 0) + (newAtom.z || 0);
              newatoms.push(newAtom);
            }
            model.addAtoms(/** @type {Array<import("./specs").AtomSpec>} */ (newatoms));
          }
        }
      }
    }
    return true;
  }

  static addLineDashed(spec, s) {
    spec.dashLength = spec.dashLength || 0.5;
    spec.gapLength = spec.gapLength || 0.5;
    spec.start = spec.start || {};
    spec.end = spec.end || {};

    let p1 = new Vector3(spec.start.x || 0, spec.start.y || 0, spec.start.z || 0);
    const p2 = new Vector3(spec.end.x, spec.end.y || 0, spec.end.z || 0);

    const dir = new Vector3();
    let dash = new Vector3();
    let gap = new Vector3();
    let length;
    let dashAmt;
    let gapAmt;
    const temp = p1.clone();
    let drawn = 0;

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
      if (drawn + dashAmt > length) {
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
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended

    return s;
  }

  /**
   * Add custom shape component from user supplied function
   *
   * @function GLViewer#addCustom
   * @param {import("./specs").CustomShapeSpec} spec - Style specification
   * @return {GLShape}
   @example
   function triangle(viewer) {
      var vertices = [];
      var normals = [];
      var colors = [];
      var r = 20;
      //triangle
      vertices.push(new Vector3(0,0,0));
      vertices.push(new Vector3(r,0,0));
      vertices.push(new Vector3(0,r,0));
  
      normals.push(new Vector3(0,0,1));
      normals.push(new Vector3(0,0,1));
      normals.push(new Vector3(0,0,1));
  
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
  addCustom(spec) {
    const s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    s.addCustom(spec);
    this.shapes.push(s);
    s.finalize(); // finalize shape for memory efficiency, assume shape won't be extended

    return s;
  }

  /**
   * Construct isosurface from volumetric data in gaussian cube format
   * @function GLViewer#addVolumetricData
   * @param {String} data - Input file contents
   * @param {String} format - Input file format
   * @param {import('./specs').VolumetricRendererSpec} spec - Shape style specification
   * @return {GLShape|import("./volume").GLVolumetricRender}
   *
   * @example
              
      $.get('data/bohr.cube', function(data) {
        viewer.addVolumetricData(data, "cube", {isoval: -0.01, color: "red", opacity: 0.95});
        viewer.setStyle({cartoon:{},stick:{}});
        viewer.zoomTo();
        viewer.render();
      });
     */
  addVolumetricData(data, format, spec) {
    spec = spec || {};

    const voldata = new VolumeData(data, format);
    if (spec.transferfn) {
      // volumetric rendering
      return this.addVolumetricRender(voldata, spec);
    }
    return this.addIsosurface(voldata, spec);
  }

  /**
   * Construct isosurface from volumetric data.  This is more flexible
   * than addVolumetricData, but can not be used with py3Dmol.
   * @function GLViewer#addIsosurface
   * @param {VolumeData} data - volumetric data
   * @param {import('./specs').IsoSurfaceSpec} spec - Shape style specification
   * @param {Function} [callback] - Callback function to be called when isosurface is constructed
   * @return {GLShape}
   *
   @example
   $.get('../test_structs/benzene-homo.cube', function(data){
            var voldata = new VolumeData(data, "cube");
            viewer.addIsosurface(voldata, {isoval: 0.01,
                                           color: "blue"});
            viewer.addIsosurface(voldata, {isoval: -0.01,
                                           color: "red"});
            viewer.zoomTo();
            viewer.render();
          });
   */
  addIsosurface(data, spec, callback) {
    const s = new GLShape(spec);
    s.shapePosition = this.shapes.length;
    s.addIsosurface(data, spec, callback);
    this.shapes.push(s);
    return s;
  }

  /**
   * Create volumetric renderer for volumetricData
   * @function GLViewer#addVolumetricRender
   * @param {VolumeData} data - volumetric data
   * @param {import("./specs").VolumetricRendererSpec} spec - specification of volumetric render
   * @return {GLVolumetricRender}
   *
   */
  addVolumetricRender(data, spec) {
    const s = new GLVolumetricRender(data, spec);
    s.shapePosition = this.shapes.length;
    this.shapes.push(s);
    return s;
  }

  /**
   * Return true if volumetric rendering is supported (WebGL 2.0 required)
   *
   * @function GLViewer#hasVolumetricRender
   * @return {boolean}
   */
  hasVolumetricRender() {
    return this.renderer.supportsVolumetric();
  }

  /**
   * Enable/disable fog for content far from the camera
   *
   * @function GLViewer#enableFog
   * @param {boolean} fog whether to enable or disable the fog
   */
  enableFog(fog) {
    if (fog) {
      this.scene.fog = new Fog(this.bgColor, 100, 200);
    } else {
      this.config.disableFog = true;
      this.show();
    }
  }

  /**
   * Sets the atomlists of all models in the viewer to specified frame.
   * Shapes and labels can also be displayed by frame.
   * Sets to last frame if framenum out of range
   *
   * @function GLViewer#setFrame
   * @param {number} framenum - fame index to use, starts at zero
   * @return {Promise<void[]>}
   */
  async setFrame(framenum) {
    this.viewer_frame = framenum;
    const modelMap = this.models.map(model => model.setFrame(framenum, this));
    return Promise.all(modelMap);
  }

  /**
   * Gets the current viewer frame.
   *
   * @function GLViewer#getFrame
   */
  getFrame() {
    return this.viewer_frame;
  }

  /**
   * Returns the number of frames that the model with the most frames in the viewer has
   *
   * @function GLViewer#getNumFrames
   * @return {number}
   */
  getNumFrames() {
    let mostFrames = 0;
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
  }

  /**
   * Animate all models in viewer from their respective frames
   * @function GLViewer#animate
   * @param {Object} options - can specify interval (speed of animation), loop (direction
   * of looping, 'backward', 'forward' or 'backAndForth'), step interval between frames ('step'), and reps (numer of repetitions, 0 indicates infinite loop)
   *
   */
  animate(options) {
    this.incAnim();
    let interval = 100;
    let loop = 'forward';
    let reps = 0;
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
    const mostFrames = this.getNumFrames();
    const that = this;
    let currFrame = 0;
    let inc = 1;
    if (options.step) {
      inc = options.step;
      reps /= inc;
    }
    let displayCount = 0;
    const displayMax = mostFrames * reps;
    let time = Date.now();
    let intervalID;

    const resolve = () => {
      this.render();
      if (!this.getCanvas().isConnected && this.renderer.isLost()) {
        // we no longer exist
        this.stopAnimate();
      } else if (++displayCount == displayMax || !that.isAnimated()) {
        clearTimeout(intervalID);
        this.animationTimers.delete(intervalID);
        this.decAnim();
      } else {
        let newInterval = interval - (Date.now() - time);
        newInterval = newInterval > 0 ? newInterval : 0;
        this.animationTimers.delete(intervalID);
        intervalID = setTimeout(display, newInterval, loop);
        this.animationTimers.add(intervalID);
      }
    };

    const display = direction => {
      time = Date.now();
      if (direction == 'forward') {
        this.setFrame(currFrame).then(() => {
          currFrame = (currFrame + inc) % mostFrames;
          resolve();
        });
      } else if (direction == 'backward') {
        this.setFrame(mostFrames - 1 - currFrame).then(() => {
          currFrame = (currFrame + inc) % mostFrames;
          resolve();
        });
      } else {
        // back and forth
        this.setFrame(currFrame).then(() => {
          currFrame += inc;
          inc *= currFrame % (mostFrames - 1) == 0 ? -1 : 1;
          resolve();
        });
      }
    };

    intervalID = setTimeout(display, 0, loop);
    this.animationTimers.add(intervalID);
    return this;
  }

  /**
   * Stop animation of all models in viewer
   * @function GLViewer#stopAnimate
   */
  stopAnimate() {
    this.animated = 0;
    this.animationTimers.forEach(timer => {
      clearTimeout(timer);
    });
    this.animationTimers = new Set();
    return this;
  }

  /**
   * Return true if viewer is currently being animated, false otherwise
   * @function GLViewer#isAnimated
   * @return {boolean}
   */
  isAnimated() {
    return this.animated > 0;
  }

  /**
   * Create and add model to viewer, given molecular data and its format
   *
   * @function GLViewer#addModel
   * @param {string} [data] - Input data
   * @param {string} [format] - Input format ('pdb', 'sdf', 'xyz', 'pqr', or 'mol2')
   * @param {import('./specs').ParserOptionsSpec} [options] - format dependent options. Attributes depend on the input file format.
   * @example
        viewer.setViewStyle({style:"outline"});
        $.get('data/1fas.pqr', function(data){
            viewer.addModel(data, "pqr");
            $.get("data/1fas.cube",function(volumedata){
                viewer.addSurface(SurfaceType.VDW, {opacity:0.85,voldata: new VolumeData(volumedata, "cube"), volscheme: new Gradient.RWB(-10,10)},{});
            viewer.render();
            });
            viewer.zoomTo();
        });
   *
   * @return {GLModel}
   */
  addModel(data, format, options) {
    if (options && !options.defaultcolors) {
      options.defaultcolors = this.config.defaultcolors;
      options.cartoonQuality = this.config.cartoonQuality;
    } else if (typeof options === 'undefined') {
      options = {
        defaultcolors: elementColors.defaultcolors,
        cartoonQuality: this.config.cartoonQuality,
      };
    }
    const m = new GLModel(this.models.length, options);
    m.addMolData(data, format, options);
    this.models.push(m);

    return m;
  }

  /**
   * Given multimodel file and its format, add atom data to the viewer as separate models
   * and return list of these models
   *
   * @function GLViewer#addModels
   * @param {string} data - Input data
   * @param {string} format - Input format (see {@link FileFormats})
   * @return {Array<GLModel>}
   */
  addModels(data, format, options) {
    options = options || {};
    options.multimodel = true;
    options.frames = true;

    const modelatoms = GLModel.parseMolData(data, format, options);

    for (let i = 0; i < modelatoms.length; i++) {
      const newModel = new GLModel(this.models.length, options.defaultcolors);
      newModel.setAtomDefaults(modelatoms[i]);
      newModel.addFrame(modelatoms[i]);
      newModel.setFrame(0);
      if (modelatoms.modelData) newModel.setModelData(modelatoms.modelData[i]);
      newModel.setDontDuplicateAtoms(!options.duplicateAssemblyAtoms);
      this.models.push(newModel);
    }

    return this.models;
  }

  /**
   * Create and add model to viewer. Given multimodel file and its format,
   * different atomlists are stored in model's frame
   * property and model's atoms are set to the 0th frame
   *
   * @function GLViewer#addModelsAsFrames
   * @param {string} data - Input data
   * @param {string} format - Input format (see {@link FileFormats})
   * @return {GLModel}
   *
   * @example
          $.get('../test_structs/multiple2.xyz', function(data){
            viewer.addModelsAsFrames(data, "xyz");
            viewer.animate({loop: "forward",reps: 1});
            viewer.setStyle({stick:{colorscheme:'magentaCarbon'}});
            viewer.zoomTo();
            viewer.render();
        });
   */
  addModelsAsFrames(data, format, options) {
    options = options || {};
    options.multimodel = true;
    options.frames = true;
    const m = new GLModel(this.models.length, options.defaultcolors);
    m.addMolData(data, format, options);
    this.models.push(m);

    return m;
  }

  /**
   * Create and add model to viewer. Given multimodel file and its format,
   * all atoms are added to one model
   *
   * @function GLViewer#addAsOneMolecule
   * @param {string} data - Input data
   * @param {string} format - Input format (see {@link FileFormats})
   * @return {GLModel}
   @example
        $.get('../test_structs/multiple.sdf', function(data){
            viewer.addAsOneMolecule(data, "sdf");
            viewer.zoomTo();
            viewer.render();
        });
   */
  addAsOneMolecule(data, format, options) {
    options = options || {};
    options.multimodel = true;
    options.onemol = true;
    const m = new GLModel(this.models.length, options.defaultcolors);
    m.addMolData(data, format, options);
    this.models.push(m);

    return m;
  }

  /**
   * Delete specified model from viewer
   *
   * @function GLViewer#removeModel
   * @param {GLModel} modelIn
   * @returns {boolean} true if model was deleted, false if model was not found
   */
  removeModel(modelIn) {
    const model = this.getModel(modelIn);
    if (!model) return false;
    model.removegl(this.modelGroup);
    delete this.models[model.getID()];
    // clear off back of model array
    while (this.models.length > 0 && typeof this.models[this.models.length - 1] === 'undefined')
      this.models.pop();
    return true;
  }

  /**
   * Delete all existing models
   * @function GLViewer#removeAllModels
   * @returns {boolean} true if models were deleted, false if no models were found
   */
  removeAllModels() {
    if (this.models.length === 0) return false;
    for (let i = 0; i < this.models.length; i++) {
      const model = this.models[i];
      if (model) model.removegl(this.modelGroup);
    }
    this.models.splice(0, this.models.length); // don't simply overwrite array in case linked
    return true;
  }

  /**
   * Export one or all of the loaded models into ChemDoodle compatible JSON.
   * @function GLViewer#exportJSON
   * @param {boolean} includeStyles - Whether or not to include style information.
   * @param {number} modelID - Optional parameter for which model to export. If left out, export all of them.
   * @return {string}
   */
  exportJSON(includeStyles, modelID) {
    const object = {};
    if (modelID === undefined) {
      object.m = this.models.map(model => model.toCDObject(includeStyles));
    } else {
      object.m = [this.models[modelID].toCDObject()];
    }
    return JSON.stringify(object);
  }

  /** return a VRML string representation of the scene.  Include VRML header information
   * @function GLViewer#exportVRML
   * @returns {string}
   */
  exportVRML() {
    const savedmodelGroup = this.modelGroup;
    this.applyToModels('removegl', this.modelGroup); // cleanup
    this.modelGroup = new Object3D();
    // rendering with plain mesh
    this.render(null, {supportsImposters: false, supportsAIA: false, regen: true});
    const ret = `#VRML V2.0 utf8\n${this.modelGroup.vrml()}\n`;
    this.applyToModels('removegl', this.modelGroup); // cleanup
    this.modelGroup = savedmodelGroup;
    return ret;
  }

  /**
   * Create a new model from atoms specified by sel.
   * If extract, removes selected atoms from existing models
   *
   * @function GLViewer#createModelFrom
   * @param {Object} sel - Atom selection specification
   * @param {boolean=} extract - If true, remove selected atoms from existing models
   * @return {GLModel}
   */
  createModelFrom(sel, extract) {
    const m = new GLModel(this.models.length, this.config.defaultColors);
    for (let i = 0; i < this.models.length; i++) {
      if (this.models[i]) {
        const atoms = this.models[i].selectedAtoms(sel);
        m.addAtoms(atoms);
        if (extract) this.models[i].removeAtoms(atoms);
      }
    }
    this.models.push(m);
    return m;
  }

  applyToModels(func, sel, value1, value2, value3, value4, value5) {
    // apply func to all models that are selected by sel with value1 and 2
    const ms = this.getModelList(sel);
    for (let i = 0; i < ms.length; i++) {
      ms[i][func](sel, value1, value2, value3, value4, value5);
    }
  }

  /**
   * Set style properties to all selected atoms
   *
   * @function GLViewer#setStyle
   * @param {import('./specs').AtomSelectionSpec|import("./specs").AtomStyleSpec} sel - Atom selection specification
   * @param {import('./specs').AtomStyleSpec} [style] - Style spec to apply to specified atoms
   *
   * @example
      viewer.setBackgroundColor(0xffffffff);
         download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
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
  setStyle(sel, style) {
    if (typeof style === 'undefined') {
      // if a single argument is provided, assume it is a style and select all
      style = /** @type {import('./specs').AtomStyleSpec} */ (sel);
      sel = {};
    }

    this.applyToModels('setStyle', sel, style, false);
    return this;
  }

  /**
   * Add style properties to all selected atoms
   * @function GLViewer#addStyle
   * @param {import('./specs').AtomSelectionSpec|import('./specs').AtomStyleSpec} sel - Atom selection specification
   * @param {import('./specs').AtomStyleSpec} style - style spec to add to specified atoms
   * @example
    download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
      viewer.setStyle({cartoon:{}});
      //keep cartoon style, but show thick sticks for chain A
      viewer.addStyle({chain:'A'},{stick:{radius:.5,colorscheme:"magentaCarbon"}});
      viewer.zoomTo();
      viewer.render();
    });
   */
  addStyle(sel, style) {
    if (typeof style === 'undefined') {
      // if a single argument is provided, assume it is a style and select all
      style = /** @type {import('./specs').AtomStyleSpec} */ (sel);
      sel = {};
    }
    this.applyToModels('setStyle', sel, style, true);
    return this;
  }

  /**
   * Set click-handling properties to all selected atomsthis.
   *
   * @function GLViewer#setClickable
   * @param {import('./specs').AtomSelectionSpec} sel - atom selection to apply clickable settings to
   * @param {boolean} clickable - whether click-handling is enabled for the selection
   * @param {function} callback - function called when an atom in the selection is clicked
   *
   * @example
      download("cid:307900",viewer,{},function(){
             viewer.setStyle({},{sphere:{}});
             viewer.setClickable({},true,function(atom,viewer,event,container) {
                 viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8});
             });
             viewer.render();
      });
   */
  setClickable(sel, clickable, callback) {
    this.applyToModels('setClickable', sel, clickable, callback);
    return this;
  }

  /** Set hoverable and callback of selected atoms
   *
   * @function GLViewer#setHoverable
   * @param {import('./specs').AtomSelectionSpec} sel - atom selection to apply hoverable settings to
   * @param {boolean} hoverable - whether hover-handling is enabled for the selection
   * @param {()=>any} hoverCallback - function called when an atom in the selection is hovered over
   * @param {()=>any} unhoverCallback - function called when the mouse moves out of the hover area
   * @returns {GLViewer}
   * @example
    download("pdb:1ubq",viewer,{},function(){
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
  setHoverable(sel, hoverable, hoverCallback, unhoverCallback) {
    this.applyToModels('setHoverable', sel, hoverable, hoverCallback, unhoverCallback);
    return this;
  }

  /** enable context menu and callback of selected atoms
   *
   * @function GLViewer#enableContextMenu
   * @param {import('./specs').AtomSelectionSpec} sel - atom selection to apply hoverable settings to
   * @param {boolean} contextMenuEnabled - whether contextMenu-handling is enabled for the selection
   * @returns {GLViewer}
   */
  enableContextMenu(sel, contextMenuEnabled) {
    this.applyToModels('enableContextMenu', sel, contextMenuEnabled);
    return this;
  }

  /**
   * If  atoms have dx, dy, dz properties (in some xyz files), vibrate populates each model's frame property based on parameters.
   * Models can then be animated
   *
   * @function GLViewer#vibrate
   * @param {number} numFrames - number of frames to be created, default to 10
   * @param {number} amplitude - amplitude of distortion, default to 1 (full)
   * @param {boolean} bothWays - if true, extend both in positive and negative directions by numFrames
   * @param {import('./specs').ArrowSpec} arrowSpec - specification for drawing animated arrows. If color isn't specified, atom color (sphere, stick, line preference) is used.
   */
  vibrate(numFrames, amplitude, bothWays, arrowSpec) {
    this.applyToModels('vibrate', numFrames, amplitude, bothWays, this, arrowSpec);
    return this;
  }

  /**
   * @function GLViewer#setColorByProperty
   * @param {import('./specs').AtomSelectionSpec} sel
   * @param {any} prop
   * @param {any} scheme
   */
  setColorByProperty(sel, prop, scheme, range) {
    this.applyToModels('setColorByProperty', sel, prop, scheme, range);
    return this;
  }

  /**
   * @function GLViewer#setColorByElement
   * @param {import('./specs').AtomSelectionSpec} sel
   * @param {any} colors
   */
  setColorByElement(sel, colors) {
    this.applyToModels('setColorByElement', sel, colors);
    return this;
  }

  /**
   *
   * @param {import('./specs').AtomSpec[]} atomlist
   * @param {number[][]} extent
   * @return {Array<import('./specs').AtomSpec>}
   */
  static getAtomsWithin(atomlist, extent) {
    const ret = [];

    for (let i = 0; i < atomlist.length; i++) {
      const atom = atomlist[i];
      if (typeof atom == 'undefined') continue;

      if (atom.x < extent[0][0] || atom.x > extent[1][0]) continue;
      if (atom.y < extent[0][1] || atom.y > extent[1][1]) continue;
      if (atom.z < extent[0][2] || atom.z > extent[1][2]) continue;
      ret.push(atom);
    }
    return ret;
  }

  /**
   * return volume of extent
   * @param {number[][]} extent
   * @return {number}
   */
  static volume(extent) {
    const w = extent[1][0] - extent[0][0];
    const h = extent[1][1] - extent[0][1];
    const d = extent[1][2] - extent[0][2];
    return w * h * d;
  } // volume

  /*
   * Break up bounding box/atoms into smaller pieces so we can parallelize
   * with webworkers and also limit the size of the working memory Returns
   * a list of bounding boxes with the corresponding atoms. These extents
   * are expanded by 4 angstroms on each side.
   *
   * @param {number[][]} extent
   * @param {import('./specs').AtomSpec[]} atomlist
   * @param {import('./specs').AtomSpec[]} atomstoshow
   * @return {Array}
   */
  static carveUpExtent(extent, atomlist, atomstoshow) {
    const ret = [];

    const index2atomlist = {}; // map from atom.index to position in atomlist
    for (let i = 0, n = atomlist.length; i < n; i++) {
      index2atomlist[atomlist[i].index] = i;
    }

    const atomsToListIndex = function (atoms) {
      // return a list of indices into atomlist
      const ret = [];
      for (let i = 0, n = atoms.length; i < n; i++) {
        if (atoms[i].index in index2atomlist) ret.push(index2atomlist[atoms[i].index]);
      }
      return ret;
    };
    const copyExtent = function (extent) {
      // copy just the dimensions
      const ret = [];
      ret[0] = [extent[0][0], extent[0][1], extent[0][2]];
      ret[1] = [extent[1][0], extent[1][1], extent[1][2]];
      return ret;
    }; // copyExtent
    const splitExtentR = function (extent) {
      // recursively split until volume is below maxVol
      if (GLViewer.volume(extent) < maxVolume) {
        return [extent];
      }
      // find longest edge
      const w = extent[1][0] - extent[0][0];
      const h = extent[1][1] - extent[0][1];
      const d = extent[1][2] - extent[0][2];

      let index;

      if (w > h && w > d) {
        index = 0;
      } else if (h > w && h > d) {
        index = 1;
      } else {
        index = 2;
      }

      // create two halves, splitting at index
      const a = copyExtent(extent);
      const b = copyExtent(extent);
      const mid = (extent[1][index] - extent[0][index]) / 2 + extent[0][index];
      a[1][index] = mid;
      b[0][index] = mid;

      const alist = splitExtentR(a);
      const blist = splitExtentR(b);
      return alist.concat(blist);
    }; // splitExtentR

    // divide up extent
    const splits = splitExtentR(extent);
    // now compute atoms within expanded (this could be more efficient)
    const off = 6; // enough for water and 2*r, also depends on scale

    // factor
    for (let i = 0, n = splits.length; i < n; i++) {
      const e = copyExtent(splits[i]);
      e[0][0] -= off;
      e[0][1] -= off;
      e[0][2] -= off;
      e[1][0] += off;
      e[1][1] += off;
      e[1][2] += off;

      const atoms = GLViewer.getAtomsWithin(atomlist, e);
      const toshow = GLViewer.getAtomsWithin(atomstoshow, splits[i]);

      // ultimately, divide up by atom for best meshing
      ret.push({
        extent: splits[i],
        atoms: atomsToListIndex(atoms),
        toshow: atomsToListIndex(toshow),
      });
    }

    return ret;
  }

  // create a mesh defined from the passed vertices and faces and material
  // Just create a single geometry chunk - broken up whether sync or not
  /**
   *
   * @param {import('./specs').AtomSpec[]} atoms
   * @param {{vertices:Vector3[]; faces:number[]}} VandF
   * @param {import("./WebGL/materials").Material} mat
   * @return {Mesh}
   */
  static generateSurfaceMesh(atoms, VandF, mat) {
    const geo = new Geometry(true);
    // Only one group per call to generate surface mesh (addSurface
    // should split up mesh render)
    const geoGroup = geo.updateGeoGroup(0);

    // set colors for vertices
    const colors = [];
    for (let i = 0, il = atoms.length; i < il; i++) {
      const atom = atoms[i];
      if (atom) {
        if (typeof atom.surfaceColor != 'undefined') {
          colors[i] = atom.surfaceColor;
        } else if (atom.color)
          // map from atom
          colors[i] = CC.color(atom.color);
      }
    }

    const {vertexArray} = geoGroup;

    // reconstruct vertices and faces
    const v = VandF.vertices;
    for (let i = 0, il = v.length; i < il; i++) {
      const offset = geoGroup.vertices * 3;
      vertexArray[offset] = v[i].x;
      vertexArray[offset + 1] = v[i].y;
      vertexArray[offset + 2] = v[i].z;
      geoGroup.vertices++;
    }

    // set colorArray of there are per-atom colors
    const {colorArray} = geoGroup;

    if (mat.voldata && mat.volscheme) {
      // convert volumetric data into colors
      const scheme = mat.volscheme;
      const {voldata} = mat;
      const range = scheme.range() || [-1, 1];
      for (let i = 0, il = v.length; i < il; i++) {
        const val = voldata.getVal(v[i].x, v[i].y, v[i].z);
        const col = CC.color(scheme.valueToHex(val, range));
        const offset = i * 3;
        colorArray[offset] = col.r;
        colorArray[offset + 1] = col.g;
        colorArray[offset + 2] = col.b;
      }
    } else if (colors.length > 0) {
      // have atom colors
      for (let i = 0, il = v.length; i < il; i++) {
        const A = v[i].atomid;
        const offsetA = i * 3;

        colorArray[offsetA] = colors[A].r;
        colorArray[offsetA + 1] = colors[A].g;
        colorArray[offsetA + 2] = colors[A].b;
      }
    }

    const {faces} = VandF;
    geoGroup.faceidx = faces.length; // *3;
    geo.initTypedArrays();

    const verts = geoGroup.vertexArray;
    const {normalArray} = geoGroup;
    let vA;
    let vB;
    let vC;
    let norm;

    // Setup colors, faces, and normals
    for (let i = 0, il = faces.length; i < il; i += 3) {
      // var a = faces[i].a, b = faces[i].b, c = faces[i].c;
      const a = faces[i];
      const b = faces[i + 1];
      const c = faces[i + 2];
      const offsetA = a * 3;
      const offsetB = b * 3;
      const offsetC = c * 3;

      // setup Normals
      // todo - calculate normals in parallel code
      vA = new Vector3(verts[offsetA], verts[offsetA + 1], verts[offsetA + 2]);
      vB = new Vector3(verts[offsetB], verts[offsetB + 1], verts[offsetB + 2]);
      vC = new Vector3(verts[offsetC], verts[offsetC + 1], verts[offsetC + 2]);

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
    const mesh = new Mesh(geo, mat);
    mesh.doubleSided = true;
    return mesh;
  }

  // do same thing as worker in main thread
  /**
   * @param {number} surfaceType
   * @param {Array} expandedExtent
   * @param {Array} extendedAtoms
   * @param {Array} atomsToShow
   * @param {import('./specs').AtomSpec[]} atoms
   * @param {number} vol
   * @return {Object}
   */
  static generateMeshSyncHelper(surfaceType, expandedExtent, extendedAtoms, atomsToShow, atoms, vol) {
    // var time = new Date();
    const ps = new ProteinSurface();
    ps.initparm(expandedExtent, surfaceType !== 1, vol);

    // var time2 = new Date();
    // console.log("initialize " + (time2 - time) + "ms");
    ps.fillvoxels(atoms, extendedAtoms);

    // var time3 = new Date();
    // console.log("fillvoxels " + (time3 - time2) + "  " + (time3 - time) + "ms");
    ps.buildboundary();

    if (surfaceType == SurfaceType.SES || surfaceType == SurfaceType.MS) {
      ps.fastdistancemap();
      ps.boundingatom(false);
      ps.fillvoxelswaals(atoms, extendedAtoms);
    }

    //            var time4 = new Date();
    // console.log("buildboundaryetc " + (time4 - time3) + "  " + (time4 - time) + "ms");
    ps.marchingcube(surfaceType);

    //            var time5 = new Date();
    // console.log("marching cube " + (time5 - time4) + "  "+ (time5 - time) + "ms");
    return ps.getFacesAndVertices(atomsToShow);
  }

  /**
   *
   * @param {any} style
   * @return {MeshLambertMaterial}
   */
  static getMatWithStyle(style) {
    const mat = new MeshLambertMaterial();
    mat.vertexColors = VertexColors;

    for (const prop in style) {
      if (!(prop === 'color' || prop === 'map') && style[prop]) {
        mat[prop] = style[prop];
      }
    }
    if (style.opacity !== undefined) {
      if (style.opacity === 1) {
        mat.transparent = false;
      } else {
        mat.transparent = true;
      }
    }

    return mat;
  }

  /**
   * Adds an explicit mesh as a surface object.
   * @function GLViewer#addMesh
   * @param {Mesh} mesh
   * @returns {number} surfid
   */
  addMesh(mesh) {
    const surfobj = {
      geo: mesh.geometry,
      mat: mesh.material,
      done: true,
      finished: false, // the rendered finishes surfaces when they are done
    };
    const surfid = this.nextSurfID();
    this.surfaces[surfid] = surfobj;
    return surfid;
  }

  // return a shallow copy of list l, e.g., for atoms so we can
  // ignore superficial changes (ie surfacecolor, position) that happen
  // while we're surface building
  static shallowCopy(l) {
    const ret = [];
    const {length} = l;
    for (let i = 0; i < length; i++) {
      ret[i] = extend({}, l[i]);
    }
    return ret;
  }


  /**
   * Add surface representation to atoms
   * @function GLViewer#addSurface
   * @param {number | string} surfaceTypeIn - Surface type (VDW, MS, SAS, or SES)
   * @param {import('./specs').SurfaceStyleSpec} style - optional style specification for surface material (e.g. for different coloring scheme, etc)
   * @param {import('./specs').AtomSelectionSpec} atomsel - Show surface for atoms in this selection
   * @param {import('./specs').AtomSelectionSpec} allsel - Use atoms in this selection to calculate surface; may be larger group than 'atomsel'
   * @param {import('./specs').AtomSelectionSpec} focus - Optionally begin rendering surface specified atoms
   * @param {() => any} surfacecallback - function to be called after setting the surface
   * @return {Promise} promise - Returns a promise that ultimately resovles to the surfid.  
   * Returns surfid immediately if surfacecallback is specified. 
   *  Returned promise has a [surfid, GLViewer, style, atomsel, allsel, focus] fields for immediate access.
   */
  async addSurface(surfaceTypeIn, style, atomsel, allsel, focus, surfacecallback) {
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
    // surfacecallback gets called when done
    const surfid = this.nextSurfID();
    let surfaceType = normalizeSurfaceType(surfaceTypeIn);
    if (!surfaceType) throw new Error(`Surface type not recognized: ${surfaceTypeIn}`);
    let mat = null;
    // atoms specified by this selection
    let atomlist = null;
    let focusSele = null;
    // TODO: currently generating a shallow copy to avoid problems when atoms are chagned
    // during surface generation - come up with a better solution
    const atomsToShow = GLViewer.shallowCopy(this.getAtomsFromSel(atomsel));
    if (!allsel) {
      atomlist = atomsToShow;
    } else {
      atomlist = GLViewer.shallowCopy(this.getAtomsFromSel(allsel));
    }

    adjustVolumeStyle(style);
    let symmetries = false;
    let n;
    for (n = 0; n < this.models.length; n++) {
      if (this.models[n]) {
        const symMatrices = this.models[n].getSymmetries();
        if (symMatrices.length > 1 || (symMatrices.length == 1 && !symMatrices[0].isIdentity())) {
          symmetries = true;
          break;
        }
      }
    }
    const self = this;

    const addSurfaceHelper = (surfobj, atomlist, atomsToShow) => {
      // function returns promise with surfid resolved
      if (!focus) {
        focusSele = atomsToShow;
      } else {
        focusSele = GLViewer.shallowCopy(this.getAtomsFromSel(focus));
      }

      let atom;
      //                var time = new Date();
      const extent = getExtent(atomsToShow, true);
      if (style.map && style.map.prop) {
        // map color space using already set atom properties
        const {prop} = style.map;
        /** @type {Gradient} */
        const scheme = style.map.scheme || style.map.gradient || new Gradient.RWB();
        let range = scheme.range();
        if (!range) {
          range = getPropertyRange(atomsToShow, prop);
        }
        style.colorscheme = {prop, gradient: scheme};
      }

      // cache surface color on each atom
      for (let i = 0, il = atomlist.length; i < il; i++) {
        atom = atomlist[i];
        atom.surfaceColor = getColorFromStyle(atom, style);
      }

      const totalVol = GLViewer.volume(extent); // used to scale resolution
      const extents = GLViewer.carveUpExtent(extent, atomlist, atomsToShow);

      if (focusSele && focusSele.length && focusSele.length > 0) {
        const seleExtent = getExtent(focusSele, true);
        // sort by how close to center of seleExtent
        const sortFunc = (a, b) => {
          const distSq = (ex, sele) => {
            // distance from e (which has no center of mass) and
            // sele which does
            const e = ex.extent;
            const x = e[1][0] - e[0][0];
            const y = e[1][1] - e[0][1];
            const z = e[1][2] - e[0][2];
            let dx = x - sele[2][0];
            dx *= dx;
            let dy = y - sele[2][1];
            dy *= dy;
            let dz = z - sele[2][2];
            dz *= dz;

            return dx + dy + dz;
          };
          const d1 = distSq(a, seleExtent);
          const d2 = distSq(b, seleExtent);
          return d1 - d2;
        };
        extents.sort(sortFunc);
      }

      // console.log("Extents " + extents.length + "  "+ (+new Date() - time) + "ms");
      const reducedAtoms = [];
      // to reduce amount data transfered, just pass x,y,z,serial and elem
      for (let i = 0, il = atomlist.length; i < il; i++) {
        atom = atomlist[i];
        reducedAtoms[i] = {
          x: atom.x,
          y: atom.y,
          z: atom.z,
          serial: i,
          elem: atom.elem,
        };
      }

      const sync = !!SyncSurface;
      if (sync) {
        // don't use worker, still break up for memory purposes
        // to keep the browser from locking up, call through setTimeout
        const callSyncHelper = function callSyncHelper(i) {
          return new Promise(resolve => {
            let VandF = GLViewer.generateMeshSyncHelper(
              surfaceType,
              extents[i].extent,
              extents[i].atoms,
              extents[i].toshow,
              reducedAtoms,
              totalVol
            );
            // complicated surfaces sometimes have > 2^16 vertices
            const VandFs = splitMesh({vertexArr: VandF.vertices, faceArr: VandF.faces});
            for (let vi = 0, vl = VandFs.length; vi < vl; vi++) {
              VandF = {vertices: VandFs[vi].vertexArr, faces: VandFs[vi].faceArr};
              const mesh = GLViewer.generateSurfaceMesh(atomlist, VandF, mat);
              mergeGeos(surfobj.geo, mesh);
            }
            self.render();
            resolve(undefined);
          });
        };
        const promises = [];
        for (let i = 0; i < extents.length; i++) {
          promises.push(callSyncHelper(i));
        }
        return Promise.all(promises).then(() => {
          surfobj.done = true;
          return Promise.resolve(surfid);
        });

        // TODO: Asynchronously generate geometryGroups (not separate
        // meshes) and merge them into a single geometry
      } // use worker

      const workers = [];
      if (surfaceType < 0) surfaceType = 0; // negative reserved for atom data
      for (let i = 0, il = numWorkers; i < il; i++) {
        const w = new Worker((/** @type {string} */(SurfaceWorker)));
        workers.push(w);
        w.postMessage({
          type: -1,
          atoms: reducedAtoms,
          volume: totalVol,
        });
      }

      return new Promise((resolve, reject) => {
        let cnt = 0;

        const releaseMemory = () => {
          if (!workers || !workers.length) return;
          workers.forEach(worker => {
            if (worker && worker.terminate) {
              worker.terminate();
            }
          });
        };

        const rfunction = (event) => {
          const VandFs = splitMesh({
            vertexArr: event.data.vertices,
            faceArr: event.data.faces,
          });
          for (let i = 0, vl = VandFs.length; i < vl; i++) {
            const VandF = {vertices: VandFs[i].vertexArr, faces: VandFs[i].faceArr};
            const mesh = GLViewer.generateSurfaceMesh(atomlist, VandF, mat);
            mergeGeos(surfobj.geo, mesh);
          }
          this.render();

          //    console.log("async mesh generation " + (+new Date() - time) + "ms");
          cnt++;
          if (cnt == extents.length) {
            surfobj.done = true;
            releaseMemory();
            resolve(surfid); // caller of helper will resolve callback if present
          }
        };

        const efunction = function (event) {
          releaseMemory();
          console.log(`${event.message} (${event.filename}:${event.lineno})`);
          reject(event);
        };

        for (let i = 0; i < extents.length; i++) {
          const worker = workers[i % workers.length];
          worker.onmessage = rfunction;

          worker.onerror = efunction;

          worker.postMessage({
            type: surfaceType,
            expandedExtent: extents[i].extent,
            extendedAtoms: extents[i].atoms,
            atomsToShow: extents[i].toshow,
          });
        }
      });
    };

    style = style || {};
    mat = GLViewer.getMatWithStyle(style);
    const surfobj = {};
    // save this.configuration of surface
    surfobj.style = style;
    surfobj.atomsel = atomsel;
    surfobj.allsel = allsel;
    surfobj.focus = focus;
    /** @type {Promise<any[]> & {surfid?:number} | null} */
    let promise = null;
    if (symmetries) {
      // do preprocessing
      const modelsAtomList = {};
      const modelsAtomsToShow = {};
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
      const promises = [];
      for (n = 0; n < this.models.length; n++) {
        if (modelsAtomsToShow[n].length > 0) {
          surfobj.push({
            geo: new Geometry(true),
            mat,
            done: false,
            finished: false,
            symmetries: this.models[n].getSymmetries(),
            // also webgl initialized
          });
          promises.push(
            addSurfaceHelper(surfobj[surfobj.length - 1], modelsAtomList[n], modelsAtomsToShow[n])
          );
        }
      }

      promise = Promise.all(promises);
    } else {
      surfobj.push({
        geo: new Geometry(true),
        mat,
        done: false,
        finished: false,
        symmetries: [new Matrix4()],
      });
      promise = addSurfaceHelper(surfobj[surfobj.length - 1], atomlist, atomsToShow);
    }
    this.surfaces[surfid] = surfobj;
    promise.surfid = surfid;

    if (surfacecallback && typeof surfacecallback == 'function') {
      promise.then(surfacecallback);
      return surfid;
    }
    return promise;
  }

  /**
   * Set the surface material to something else, must render change
   * @function GLViewer#setSurfaceMaterialStyle
   * @param {number} surf - Surface ID to apply changes to
   * @param {import('./specs').SurfaceStyleSpec} style - new material style specification
   @example
   $.get("data/9002806.cif",function(data){
      viewer.addModel(data);
      viewer.setStyle({stick:{}});
      let surf = viewer.addSurface("SAS");
      surf.then(function() {
          viewer.setSurfaceMaterialStyle(surf.surfid, {color:'blue',opacity:0.5});
          viewer.render();
          });
     });
   */
  setSurfaceMaterialStyle(surf, style) {
    this.adjustVolumeStyle(style);
    if (this.surfaces[surf]) {
      const surfArr = this.surfaces[surf];
      surfArr.style = style;
      for (let i = 0; i < surfArr.length; i++) {
        const mat = (surfArr[i].mat = GLViewer.getMatWithStyle(style));
        surfArr[i].mat.side = FrontSide;
        if (style.color) {
          surfArr[i].mat.color = style.color;
          surfArr[i].geo.colorsNeedUpdate = true;
          const c = CC.color(style.color);
          surfArr[i].geo.setColors(() => c);
        } else if (mat.voldata && mat.volscheme) {
          // convert volumetric data into colors
          const scheme = mat.volscheme;
          const {voldata} = mat;
          const cc = CC;
          const range = scheme.range() || [-1, 1];
          surfArr[i].geo.setColors((x, y, z) => {
            const val = voldata.getVal(x, y, z);
            const col = cc.color(scheme.valueToHex(val, range));
            return col;
          });
        }
        surfArr[i].finished = false; // trigger redraw
      }
    }
    return this;
  }

  // eslint-disable-next-line class-methods-use-this
  adjustVolumeStyle(style) {
    throw new Error('Method not implemented.');
  }

  /**
   * Return surface object
   * @function GLViewer#getSurface
   * @param {number} surf - surface id
   */
  getSurface(surf) {
    return this.surfaces[surf];
  }

  /**
   * Remove surface with given ID
   * @function GLViewer#removeSurface
   * @param {number} surf - surface id
   */
  removeSurface(surf) {
    const surfArr = this.surfaces[surf];
    for (let i = 0; i < surfArr.length; i++) {
      if (surfArr[i] && surfArr[i].lastGL) {
        if (surfArr[i].geo !== undefined) surfArr[i].geo.dispose();
        if (surfArr[i].mat !== undefined) surfArr[i].mat.dispose();
        this.modelGroup.remove(surfArr[i].lastGL); // remove from scene
      }
    }
    delete this.surfaces[surf];
    this.show();
    return this;
  }

  /** Remove all surfaces.
   * @function GLViewer#removeAllSurfaces */
  removeAllSurfaces() {
    for (const n in this.surfaces) {
      if (!this.surfaces[n]) continue;
      const surfArr = this.surfaces[n];
      for (let i = 0; i < surfArr.length; i++) {
        if (surfArr[i] && surfArr[i].lastGL) {
          if (surfArr[i].geo !== undefined) surfArr[i].geo.dispose();
          if (surfArr[i].mat !== undefined) surfArr[i].mat.dispose();
          this.modelGroup.remove(surfArr[i].lastGL); // remove from scene
        }
      }
      delete this.surfaces[n];
    }
    this.show();
    return this;
  }

  /** return Jmol moveto command to position this scene */
  jmolMoveTo() {
    const pos = this.modelGroup.position;
    // center on same position
    let ret = `center { ${-pos.x} ${-pos.y} ${-pos.z} }; `;
    // apply rotation
    const q = this.rotationGroup.quaternion;
    ret += `moveto .5 quaternion { ${q.x} ${q.y} ${q.z} ${q.w} };`;
    // zoom is tricky.. maybe i would be best to let callee zoom on
    // selection?
    // can either do a bunch of math, or maybe zoom to the center with a
    // fixed
    // but reasonable percentage
    return ret;
  }

  /** Clear scene of all objects
   * @function GLViewer#clear
   * */
  clear() {
    this.removeAllSurfaces();
    this.removeAllModels();
    this.removeAllLabels();
    this.removeAllShapes();
    this.show();
    return this;
  }

  // props is a list of objects that select certain atoms and enumerate
  // properties for those atoms
  /**
   * Add specified properties to all atoms matching input argument
   * @function GLViewer#mapAtomProperties
   * @param {Object} props, either array of atom selectors with associated props, or function that takes atom and sets its properties
   * @param {import('./specs').AtomSelectionSpec} sel  - subset of atoms to work on - model selection must be specified here
       @example
       $.get('../test_structs/b.sdf', function(data){
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
  mapAtomProperties(props, sel) {
    sel = sel || {};
    const atoms = this.getAtomsFromSel(sel);

    if (typeof props == 'function') {
      for (let a = 0, numa = atoms.length; a < numa; a++) {
        const atom = atoms[a];
        props(atom);
      }
    } else {
      for (let a = 0, numa = atoms.length; a < numa; a++) {
        const atom = atoms[a];
        for (let i = 0, n = props.length; i < n; i++) {
          const prop = props[i];
          if (prop.props) {
            for (const p in prop.props) {
              if (prop.props.hasOwnProperty(p)) {
                // check the atom
                if (this.atomIsSelected(atom, prop)) {
                  if (!atom.properties) atom.properties = {};
                  atom.properties[p] = prop.props[p];
                }
              }
            }
          }
        }
      }
    }
    return this;
  }

  /**
   * Synchronize this view matrix of this viewer to the passed viewer.
   * When the viewpoint of this viewer changes, the other viewer will
   * be set to this viewer's view.
   * @function GLViewer#linkViewer
   * @param {GLViewer} otherviewer
   */
  linkViewer(otherviewer) {
    this.linkedViewers.push(otherviewer);
    return this;
  }

  /**
   * Return the z distance between the model and the camera
   * @function GLViewer#getPerceivedDistance
   * @return {number} distance
   */
  getPerceivedDistance() {
    return this.CAMERA_Z - this.rotationGroup.position.z;
  }

  /**
   * Set the distance between the model and the camera
   * Essentially zooming. Useful while stereo rendering.
   * @function GLViewer#setPerceivedDistance
   */
  setPerceivedDistance(dist) {
    this.rotationGroup.position.z = this.CAMERA_Z - dist;
  }

  /**
   * Used for setting an approx value of eyeSeparation. Created for calling by StereoViewer object
   * @function GLViewer#setAutoEyeSeparation
   * @return {number} camera x position
   */
  setAutoEyeSeparation(isright, x) {
    const dist = this.getPerceivedDistance();
    if (!x) x = 5.0;
    if (isright || this.camera.position.x > 0)
      // setting a value of dist*tan(x)
      this.camera.position.x = dist * Math.tan((Math.PI / 180.0) * x);
    else this.camera.position.x = -dist * Math.tan((Math.PI / 180.0) * x);
    this.camera.lookAt(new Vector3(0, 0, this.rotationGroup.position.z));
    return this.camera.position.x;
  }

  /**
   * Set the default cartoon quality for newly created models.  Default is 5.
   * Current models are not affected.
   * @number quality, higher results in higher resolution renders
   * @function GLViewer#ui.setDefaultCartoonQuality
   */
  setDefaultCartoonQuality(val) {
    this.config.cartoonQuality = val;
  }
}

export const glmolViewer = GLViewer;
