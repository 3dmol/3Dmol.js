//@ts-check

import { Camera } from "./objects/Camera";
import { Color, Geometry, GeometryGroup, Scene } from "./core";
import { SpritePlugin } from "./SpritePlugin";
import {
  BackSide,
  DoubleSide,
  ImposterMaterial,
  Material,
  MeshOutlineMaterial,
  NearestFilter,
  RGBAFormat,
  SphereImposterMaterial,
  SphereImposterOutlineMaterial,
  StickImposterOutlineMaterial,
  UnsignedByteType,
} from "./materials";
import { Matrix3, Matrix4, Vector3 } from "./math";
import { Light, Line, Mesh, Sprite } from "./objects";
import { ShaderLib, ShaderUtils } from "./ShaderUtils";
import { Fog } from "./Fog";

export class Renderer {
  _alpha;
  _antialias;
  /**
   * @type {HTMLCanvasElement}
   */
  _canvas;
  _clearAlpha;
  _clearColor;
  /**
   * @type {Camera | null}
   */
  _currentCamera = null;
  /**
   * @type {number | undefined}
   */
  _currentGeometryGroupHash;
  _currentHeight = 0;
  _currentMaterialId = -1;
  /**
   * @type {WebGLProgram | null}
   */
  _currentProgram = null;
  _currentWidth = 0;
  /**
   * @type {WebGLTexture | null}
   */
  _depthTexture = null;
  _direction = new Vector3();
  /**
   * @type {Record<number, boolean>}
   */
  _enabledAttributes = {};
  _extColorBufferFloat;
  _extFloatLinear;
  _extFragDepth;
  _extInstanced;
  /**
   * @type {WebGLFramebuffer | null}
   */
  _fb = null;
  _geometryGroupCounter = 0;
  /**
   * @type {WebGLRenderingContext | WebGL2RenderingContext | null}
   */
  _gl;
  _lights = {
    ambient: [0, 0, 0],
    directional: {
      length: 0,
      colors: [],
      positions: [],
    },
    point: {
      length: 0,
      colors: [],
      positions: [],
      distances: [],
    },
    spot: {
      length: 0,
      colors: [],
      positions: [],
      distances: [],
      directions: [],
      anglesCos: [],
      exponents: [],
    },
    hemi: {
      length: 0,
      skyColors: [],
      groundColors: [],
      positions: [],
    },
  };
  _lightsNeedUpdate = true;
  /**
   * @type {number | boolean}
   */
  _oldBlending = -1;
  /**
   * @type {number | boolean}
   */
  _oldDepthTest = -1;
  /**
   * @type {number | boolean}
   */
  _oldDepthWrite = -1;
  /**
   * @type {number | boolean}
   */
  _oldDoubleSided = -1;
  /**
   * @type {number | boolean}
   */
  _oldFlipSided = -1;
  _oldLineWidth = null;
  _oldPolygonOffset = null;
  _outlineEnabled;
  _outlineMaterial;
  _outlineSphereImposterMaterial;
  _outlineStickImposterMaterial;
  /**
   * @type {string}
   */
  _precision;
  _premultipliedAlpha;
  _preserveDrawingBuffer;
  /**
   * @type {any[]}
   */
  _programs = [];
  _programs_counter = 0;
  _projInverse = new Matrix4();
  _projScreenMatrix = new Matrix4();
  /**
   * @type {WebGLBuffer | null}
   */
  _screenQuadVBO = null;
  /**
   * @type {WebGLProgram | null}
   */
  _screenshader = null;
  _stencil;
  /**
   * @type {WebGLTexture | null}
   */
  _targetTexture = null;
  _textureMatrix = new Matrix4();
  _vector3 = new Vector3();
  /**
   * @type {any}
   */
  _vertexattribpos;
  _viewportHeight = 0;
  _viewportWidth = 0;
  _webglversion = 1;
  /**
   * @type {Matrix4|undefined}
   */
  _worldInverse;
  autoClear = true;
  autoClearColor = true;
  autoClearDepth = true;
  autoClearStencil = true;
  autoUpdateObjects = true;
  autoUpdateScene = true;
  devicePixelRatio = 1.0; //set in setSize
  info = {
    memory: {
      programs: 0,
      geometries: 0,
      textures: 0,
    },
    render: {
      calls: 0,
      vertices: 0,
      faces: 0,
      points: 0,
    },
  };
  sortObjects = true;
  sprites;
  /**
   * @param {{ row?: any; col?: any; rows?: any; cols?: any; canvas?: any; precision?: any; alpha?: any; premultipliedAlpha?: any; antialias?: any; stencil?: any; preserveDrawingBuffer?: any; clearColor?: any; clearAlpha?: any; outline?: any; id?: any; }} parameters
   */
  constructor(parameters) {
    parameters = parameters || {};
    this.row = parameters.row;
    this.col = parameters.col;
    this.rows = parameters.rows;
    this.cols = parameters.cols;
    this._canvas = parameters.canvas || document.createElement("canvas");
    this._precision = parameters.precision || "highp";
    this._alpha = parameters.alpha || true;
    this._premultipliedAlpha = parameters.premultipliedAlpha || true;
    this._antialias = parameters.antialias || false;
    this._stencil = parameters.stencil || true;
    this._preserveDrawingBuffer = parameters.preserveDrawingBuffer || false;
    this._clearColor =
      parameters.clearColor !== undefined
        ? new Color(parameters.clearColor)
        : new Color(0x000000);
    this._clearAlpha = parameters.clearAlpha || 0;
    this._outlineMaterial = new MeshOutlineMaterial(parameters.outline);
    this._outlineSphereImposterMaterial = new SphereImposterOutlineMaterial(
      parameters.outline
    );
    this._outlineStickImposterMaterial = new StickImposterOutlineMaterial(
      parameters.outline
    );
    this._outlineEnabled = !!parameters.outline;
    this._canvas.id = parameters.id;
    // initialize
    //note setting antialis to true doesn't seem to do much and
    //causes problems on iOS Safari
    // note: duck typing RenderingContext -> WebGLRenderingContext
    try {
      // @ts-ignore
      this._gl = this._canvas.getContext("webgl2", {
        alpha: this._alpha,
        premultipliedAlpha: this._premultipliedAlpha,
        antialias: false,
        stencil: this._stencil,
        preserveDrawingBuffer: this._preserveDrawingBuffer,
      });
      if (!this._gl)
        // @ts-ignore
        this._gl = this._canvas.getContext("experimental-webgl", {
          alpha: this._alpha,
          premultipliedAlpha: this._premultipliedAlpha,
          antialias: false,
          stencil: this._stencil,
          preserveDrawingBuffer: this._preserveDrawingBuffer,
        });
      if (!this._gl)
        // @ts-ignore
        this._gl = this._canvas.getContext("webgl", {
          alpha: this._alpha,
          premultipliedAlpha: this._premultipliedAlpha,
          antialias: false,
          stencil: this._stencil,
          preserveDrawingBuffer: this._preserveDrawingBuffer,
        });
      if (!this._gl) throw "Error creating WebGL context.";

      let vers = this._gl.getParameter(this._gl.VERSION);
      this._webglversion = parseInt(vers[6]);
    } catch (error) {
      console.error(error);
      throw error;
    }

    // set default GL state
    this._gl.clearColor(0, 0, 0, 1);
    this._gl.clearDepth(1);
    this._gl.clearStencil(0);

    this._gl.enable(this._gl.DEPTH_TEST);
    this._gl.depthFunc(this._gl.LEQUAL);

    this._gl.frontFace(this._gl.CCW);
    this._gl.cullFace(this._gl.BACK);
    this._gl.enable(this._gl.CULL_FACE);

    this._gl.enable(this._gl.BLEND);
    this._gl.blendEquation(this._gl.FUNC_ADD);
    this._gl.blendFunc(this._gl.SRC_ALPHA, this._gl.ONE_MINUS_SRC_ALPHA);

    this._gl.clearColor(
      this._clearColor.r,
      this._clearColor.g,
      this._clearColor.b,
      this._clearAlpha
    );
    this._extInstanced = this._gl.getExtension("ANGLE_instanced_arrays");
    this._extFragDepth = this._gl.getExtension("EXT_frag_depth");
    this._extFloatLinear = this._gl.getExtension("OES_texture_float_linear");
    this._extColorBufferFloat = this._gl.getExtension("EXT_color_buffer_float");
    this.context = this._gl;
    this.sprites = new SpritePlugin(this);
  }

  get domElement() {
    return this._canvas;
  }

  supportsVolumetric() {
    return !this._isWebGL1();
  }

  supportedExtensions() {
    return {
      supportsAIA: Boolean(this._extInstanced),
      supportsImposters: Boolean(this._extFragDepth) || !this._isWebGL1(),
    };
  }

  getContext() {
    return this._gl;
  }

  isLost() {
    if (!this._gl) return true;
    return this._gl.isContextLost();
  }

  getPrecision() {
    return this._precision;
  }

  /**
   * @param {number} hex
   * @param {any} alpha
   */
  setClearColorHex(hex, alpha) {
    if (!this._gl) throw new Error("gl not initialized");
    this._clearColor.setHex(hex);
    this._clearAlpha = alpha;

    this._gl.clearColor(
      this._clearColor.r,
      this._clearColor.g,
      this._clearColor.b,
      this._clearAlpha
    );
  }

  /**
   * @param {{ color?: any; width?: any; pushback?: any; } | undefined} parameters
   */
  enableOutline(parameters) {
    this._outlineMaterial = new MeshOutlineMaterial(parameters);
    this._outlineSphereImposterMaterial = new SphereImposterOutlineMaterial(
      parameters
    );
    this._outlineStickImposterMaterial = new StickImposterOutlineMaterial(
      parameters
    );
    this._outlineEnabled = true;
  }

  disableOutline() {
    this._outlineEnabled = false;
  }

  setViewport() {
    if (
      this.rows != undefined &&
      this.cols != undefined &&
      this.row != undefined &&
      this.col != undefined &&
      this._gl
    ) {
      var wid = this._canvas.width / this.cols;
      var hei = this._canvas.height / this.rows;

      this._viewportWidth = wid;
      this._viewportHeight = hei;
      this._gl.enable(this._gl.SCISSOR_TEST);
      this._gl.scissor(wid * this.col, hei * this.row, wid, hei);
      this._gl.viewport(wid * this.col, hei * this.row, wid, hei);
    }
  }

  /**
   * @param {number} width
   * @param {number} height
   */
  setSize(width, height) {
    //zooming (in the browser) changes the pixel ratio and width/height
    this.devicePixelRatio =
      window.devicePixelRatio !== undefined ? window.devicePixelRatio : 1;
    //with antialiasing on (which doesn't seem to do much), render at double rsolution to eliminate jaggies
    //my iphone crashes if we do though, so as a hacky workaround, don't do it with retina displays
    if (this._antialias && this.devicePixelRatio < 2.0)
      this.devicePixelRatio *= 2.0;

    if (
      this.rows != undefined &&
      this.cols != undefined &&
      this.row != undefined &&
      this.col != undefined
    ) {
      var wid = width / this.cols;
      var hei = height / this.rows;
      this._canvas.width = width * this.devicePixelRatio;
      this._canvas.height = height * this.devicePixelRatio;

      this._viewportWidth = wid * this.devicePixelRatio;
      this._viewportHeight = hei * this.devicePixelRatio;

      this._canvas.style.width = width + "px";
      this._canvas.style.height = height + "px";

      this.setViewport();
    } else if (this._gl) {
      this._viewportWidth = this._canvas.width = width * this.devicePixelRatio;
      this._viewportHeight = this._canvas.height =
        height * this.devicePixelRatio;

      this._canvas.style.width = width + "px";
      this._canvas.style.height = height + "px";

      this._gl.viewport(
        0,
        0,
        this._gl.drawingBufferWidth,
        this._gl.drawingBufferHeight
      );
    }
    this.initFrameBuffer();
  }

  /**
   * @param {boolean | undefined} color
   * @param {boolean | undefined} depth
   * @param {boolean | undefined} stencil
   */
  clear(color, depth, stencil) {
    var bits = 0;
    if (!this._gl) throw new Error("WebGL not initialized");
    if (color === undefined || color) bits |= this._gl.COLOR_BUFFER_BIT;
    if (depth === undefined || depth) bits |= this._gl.DEPTH_BUFFER_BIT;
    if (stencil === undefined || stencil) bits |= this._gl.STENCIL_BUFFER_BIT;
    this._gl.clear(bits);
  }

  clearTarget = this.clear;

  /**
   * @param {{ side: number; imposter: any; }} material
   * @param {any} reflected
   */
  setMaterialFaces(material, reflected) {
    if (!this._gl) throw new Error("WebGL not initialized");
    var doubleSided = material.side === DoubleSide;
    var flipSided = material.side === BackSide;

    if (!material.imposter)
      //ignore reflection with imposters
      flipSided = reflected ? !flipSided : flipSided;

    // @ts-ignore
    if (this._oldDoubleSided !== doubleSided) {
      if (doubleSided) {
        this._gl.disable(this._gl.CULL_FACE);
      } else {
        this._gl.enable(this._gl.CULL_FACE);
      }

      // @ts-ignore
      this._oldDoubleSided = doubleSided;
    }

    // @ts-ignore
    if (this._oldFlipSided !== flipSided) {
      if (flipSided) {
        this._gl.frontFace(this._gl.CW);
      } else {
        this._gl.frontFace(this._gl.CCW);
      }

      // @ts-ignore
      this._oldFlipSided = flipSided;
    }

    this._gl.cullFace(this._gl.BACK);
  }

  /**
   * @param {boolean | number} depthTest
   */
  setDepthTest(depthTest) {
    if (!this._gl) throw new Error("WebGL not initialized");
    if (this._oldDepthTest !== depthTest) {
      if (depthTest) {
        this._gl.enable(this._gl.DEPTH_TEST);
      } else {
        this._gl.disable(this._gl.DEPTH_TEST);
      }

      this._oldDepthTest = depthTest;
    }
  }

  /**
   * @param {number | boolean} depthWrite
   */
  setDepthWrite(depthWrite) {
    if (this._oldDepthWrite !== depthWrite && this._gl) {
      this._gl.depthMask(!!depthWrite);
      this._oldDepthWrite = depthWrite;
    }
  }

  /**
   * @param {boolean} blending
   */
  setBlending(blending) {
    if (!this._gl) throw new Error("WebGL not initialized");
    if (!blending) {
      this._gl.disable(this._gl.BLEND);
    } else {
      this._gl.enable(this._gl.BLEND);
      this._gl.blendEquationSeparate(this._gl.FUNC_ADD, this._gl.FUNC_ADD);
      this._gl.blendFuncSeparate(
        this._gl.SRC_ALPHA,
        this._gl.ONE_MINUS_SRC_ALPHA,
        this._gl.ONE,
        this._gl.ONE_MINUS_SRC_ALPHA
      );
    }
    this._oldBlending = blending;
  }

  // TODO: need to set up shader attributes and uniforms as attributes on
  // material object after attaching prgm
  // We need to attach appropriate uniform variables to material after shaders
  // have been chosen
  /**
   * @param {any} material
   */
  initMaterial(material) {
    material.addEventListener("dispose", this._onMaterialDispose);

    var parameters, shaderID;

    shaderID = material.shaderID;

    if (shaderID) {
      var shader = ShaderLib[shaderID];
      material.vertexShader = shader.vertexShader;
      material.fragmentShader = shader.fragmentShader;
      material.uniforms = ShaderUtils.clone(shader.uniforms);
      // TODO: set material uniforms to shader uniform variables
    }

    parameters = {
      wireframe: material.wireframe,
      fragdepth: material.imposter,
      volumetric: material.volumetric,
    };

    material.program = this._buildProgram(
      material.fragmentShader,
      material.vertexShader,
      material.uniforms,
      parameters
    );
  }

  /**
   * @param {Scene} scene
   * @param {Camera} camera
   * @param {boolean} forceClear
   */
  render(scene, camera, forceClear) {
    if (camera instanceof Camera === false) {
      console.error("Renderer.render: camera is not an instance of Camera.");
      return;
    }

    var i,
      il,
      webglObject,
      object,
      renderList,
      lights = scene.__lights,
      fog = scene.fog;

    // reset caching for this frame

    this._currentMaterialId = -1;
    this._lightsNeedUpdate = true;

    // update scene graph

    if (this.autoUpdateScene) scene.updateMatrixWorld();

    // update camera matrices
    // Pretty sure camera's parent is always going to be undefined for our
    // purposes...
    if (camera.parent === undefined) camera.updateMatrixWorld();

    camera.matrixWorldInverse.getInverse(camera.matrixWorld);

    this._projScreenMatrix.multiplyMatrices(
      camera.projectionMatrix,
      camera.matrixWorldInverse
    );

    // update WebGL objects

    if (this.autoUpdateObjects) this.initWebGLObjects(scene);

    this.info.render.calls = 0;
    this.info.render.vertices = 0;
    this.info.render.faces = 0;
    this.info.render.points = 0;

    this._currentWidth = this._viewportWidth;
    this._currentHeight = this._viewportHeight;
    this.setViewport();
    this.setFrameBuffer();

    if ((this.autoClear || forceClear) && this._gl) {
      this._gl.clearColor(
        this._clearColor.r,
        this._clearColor.g,
        this._clearColor.b,
        this._clearAlpha
      );
      this.clear(
        this.autoClearColor,
        this.autoClearDepth,
        this.autoClearStencil
      );
    }

    // set matrices for regular objects (frustum culled)

    renderList = scene.__webglObjects;
    var hasvolumetric = false;
    for (i = 0, il = renderList.length; i < il; i++) {
      webglObject = renderList[i];
      object = webglObject.object;

      webglObject.render = false;

      if (object.visible) {
        this._setupMatrices(object, camera);
        this._unrollBufferMaterial(webglObject);
        webglObject.render = true;
        if (webglObject.volumetric) hasvolumetric = true;
      }
    }

    // set matrices for immediate objects

    var material = null;

    // opaque pass (front-to-back order)

    this.setBlending(false);

    this._renderObjects(
      scene.__webglObjects,
      true,
      "opaque",
      camera,
      lights,
      fog,
      false,
      // @ts-ignore
      material
    );

    // Render embedded labels (sprites)
    this._renderSprites(scene, camera, false);

    // prime depth buffer
    this._renderObjects(
      scene.__webglObjects,
      true,
      "blank",
      camera,
      lights,
      fog,
      true,
      // @ts-ignore
      material
    );

    // transparent pass (back-to-front order)

    this._renderObjects(
      scene.__webglObjects,
      false,
      "transparent",
      camera,
      lights,
      fog,
      true,
      // @ts-ignore
      material
    );

    //volumetric is separate
    if (hasvolumetric && this._fb) {
      //disconnect framebuffer to get depth texture
      this.reinitFrameBuffer();
      this._renderObjects(
        scene.__webglObjects,
        false,
        "volumetric",
        camera,
        lights,
        fog,
        true,
        // @ts-ignore
        material
      );
    }

    this.renderFrameBuffertoScreen();
    this.setDepthTest(true);
    this.setDepthWrite(true);

    // Render floating labels (sprites)
    this._renderSprites(scene, camera, true);
  }

  //reinitialize framebuffer without the depth texture attached so we can read to it
  //do not allocate new textures
  reinitFrameBuffer() {
    // only needed/works with webgl2
    if (this._isWebGL1()) return;
    if (!this._gl) return;

    // Create and bind the framebuffer
    this._fb = this._gl.createFramebuffer();
    this._gl.bindFramebuffer(this._gl.FRAMEBUFFER, this._fb);
    this._gl.framebufferTexture2D(
      this._gl.FRAMEBUFFER,
      this._gl.COLOR_ATTACHMENT0,
      this._gl.TEXTURE_2D,
      this._targetTexture,
      0
    );
  }

  //setup framebuffer for drawing into, assumes buffers already allocated
  setFrameBuffer() {
    if (this._isWebGL1() || !this._fb || !this._gl) return;
    let width = this._viewportWidth;
    let height = this._viewportHeight;

    //when using framebuffer, always draw from origin, will shift the viewport when we render
    this._gl.enable(this._gl.SCISSOR_TEST);
    this._gl.scissor(0, 0, width, height);
    this._gl.viewport(0, 0, width, height);

    //color texture
    this._gl.bindTexture(this._gl.TEXTURE_2D, this._targetTexture);
    this._gl.texImage2D(
      this._gl.TEXTURE_2D,
      0,
      this._gl.RGBA,
      width,
      height,
      0,
      this._gl.RGBA,
      this._gl.UNSIGNED_BYTE,
      null
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_MIN_FILTER,
      this._gl.LINEAR
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_MAG_FILTER,
      this._gl.LINEAR
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_WRAP_S,
      this._gl.CLAMP_TO_EDGE
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_WRAP_T,
      this._gl.CLAMP_TO_EDGE
    );

    //depth texture
    this._gl.bindTexture(this._gl.TEXTURE_2D, this._depthTexture);
    this._gl.texImage2D(
      this._gl.TEXTURE_2D,
      0,
      // @ts-ignore
      this._gl.DEPTH_COMPONENT32F,
      width,
      height,
      0,
      this._gl.DEPTH_COMPONENT,
      this._gl.FLOAT,
      null
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_MIN_FILTER,
      this._gl.NEAREST
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_MAG_FILTER,
      this._gl.NEAREST
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_WRAP_S,
      this._gl.CLAMP_TO_EDGE
    );
    this._gl.texParameteri(
      this._gl.TEXTURE_2D,
      this._gl.TEXTURE_WRAP_T,
      this._gl.CLAMP_TO_EDGE
    );

    //bind fb
    this._gl.bindFramebuffer(this._gl.FRAMEBUFFER, this._fb);
    this._gl.framebufferTexture2D(
      this._gl.FRAMEBUFFER,
      this._gl.COLOR_ATTACHMENT0,
      this._gl.TEXTURE_2D,
      this._targetTexture,
      0
    );
    this._gl.framebufferTexture2D(
      this._gl.FRAMEBUFFER,
      this._gl.DEPTH_ATTACHMENT,
      this._gl.TEXTURE_2D,
      this._depthTexture,
      0
    );
  }

  //allocate buffers for framebuffer, needs to be called with every resize
  initFrameBuffer() {
    // only needed/works with webgl2
    if (this._isWebGL1() || !this._gl) return;

    let width = this._viewportWidth;
    let height = this._viewportHeight;

    //when using framebuffer, always draw from origin, will shift the viewport when we render
    this._gl.enable(this._gl.SCISSOR_TEST);
    this._gl.scissor(0, 0, width, height);
    this._gl.viewport(0, 0, width, height);

    //create textures and frame buffer, will be initialized in setFrameBuffer
    this._targetTexture = this._gl.createTexture();
    this._depthTexture = this._gl.createTexture();
    this._fb = this._gl.createFramebuffer();

    // build screenshader
    var screenshader = ShaderLib.screen;

    this._screenshader = this._buildProgram(
      screenshader.fragmentShader,
      screenshader.vertexShader,
      screenshader.uniforms,
      {}
    );
    if (this._screenshader)
      this._vertexattribpos = this._gl.getAttribLocation(
        this._screenshader,
        "vertexPosition"
      );
    // create the vertex array and attrib array for the full screenquad
    var verts = [
      // First triangle:
      1.0, 1.0, -1.0, 1.0, -1.0, -1.0,
      // Second triangle:
      -1.0, -1.0, 1.0, -1.0, 1.0, 1.0,
    ];
    this._screenQuadVBO = this._gl.createBuffer();
    this._gl.bindBuffer(this._gl.ARRAY_BUFFER, this._screenQuadVBO);
    this._gl.bufferData(
      this._gl.ARRAY_BUFFER,
      new Float32Array(verts),
      this._gl.STATIC_DRAW
    );
  }

  renderFrameBuffertoScreen() {
    // only needed/works with webgl2
    if (this._isWebGL1() || this._fb === null || !this._gl) return;

    this.setViewport(); //draw texture in correct viewport

    // bind default framebuffer
    this._gl.bindFramebuffer(this._gl.FRAMEBUFFER, null);
    this._gl.clear(this._gl.COLOR_BUFFER_BIT | this._gl.DEPTH_BUFFER_BIT);
    this._gl.frontFace(this._gl.CCW);
    this._gl.cullFace(this._gl.BACK);

    // set screen shader and use it
    this._gl.useProgram(this._screenshader);
    this._currentProgram = this._screenshader;
    // disable depth test
    this.setDepthTest(-1);
    this.setDepthWrite(-1);

    // bind vertexarray buffer and texture
    this._gl.bindBuffer(this._gl.ARRAY_BUFFER, this._screenQuadVBO);
    this._gl.enableVertexAttribArray(this._vertexattribpos);
    this._gl.vertexAttribPointer(
      this._vertexattribpos,
      2,
      this._gl.FLOAT,
      false,
      0,
      0
    );

    this._gl.activeTexture(this._gl.TEXTURE0);
    this._gl.bindTexture(this._gl.TEXTURE_2D, this._targetTexture);

    // Draw 6 vertexes => 2 triangles:
    this._gl.drawArrays(this._gl.TRIANGLES, 0, 6);
  }
  getYRatio() {
    if (this.rows !== undefined && this.row !== undefined) return this.rows;
    return 1;
  }

  getXRatio() {
    if (this.cols !== undefined && this.col !== undefined) return this.cols;
    return 1;
  }

  /**
   * @param {number | undefined} width
   * @param {number | undefined} height
   */
  getAspect(width, height) {
    if (width == undefined || height == undefined) {
      width = this._canvas.width;
      height = this._canvas.height;
    }
    var aspect = width / height;
    if (
      this.rows != undefined &&
      this.cols != undefined &&
      this.row != undefined &&
      this.col != undefined
    ) {
      var wid = width / this.cols;
      var hei = height / this.rows;
      aspect = wid / hei;
    }
    return aspect;
  }

  /**
   * @param {any} texture
   * @param {any} slot
   * @param {boolean} [is3D]
   */
  setTexture(texture, slot, is3D) {
    if (!this._gl) return;
    if (texture.needsUpdate) {
      if (!texture.__webglInit) {
        texture.__webglInit = true;
        texture.addEventListener("dispose", this._onTextureDispose);
        texture.__webglTexture = this._gl.createTexture();
        this.info.memory.textures++;
      }

      this._gl.activeTexture(this._gl.TEXTURE0 + slot);
      // @ts-ignore
      var gltextureType = is3D ? this._gl.TEXTURE_3D : this._gl.TEXTURE_2D;
      this._gl.bindTexture(gltextureType, texture.__webglTexture);
      this._gl.pixelStorei(this._gl.UNPACK_FLIP_Y_WEBGL, texture.flipY);
      this._gl.pixelStorei(
        this._gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL,
        texture.premultiplyAlpha
      );
      this._gl.pixelStorei(this._gl.UNPACK_ALIGNMENT, texture.unpackAlignment);
      this._gl.pixelStorei(this._gl.PACK_ALIGNMENT, texture.unpackAlignment);

      var glFormat = this._paramToGL(texture.format),
        glType = this._paramToGL(texture.type);

      if (!is3D) {
        var image = texture.image;
        var width = image.width; //might not be defined
        var height = image.height;
        if (typeof width === "undefined") {
          //if no width,
          width = image.length;
          if (glFormat == this._gl.RGBA) {
            width /= 4; //each element takes up 4 bytes
          }
          height = 1;
        }
        this._setTextureParameters(this._gl.TEXTURE_2D, texture);
        if (
          !this._isWebGL1() &&
          glFormat !== undefined &&
          glType !== undefined
        ) {
          //webgl2
          this._gl.texImage2D(
            this._gl.TEXTURE_2D,
            0,
            glFormat,
            width,
            height,
            0,
            glFormat,
            glType,
            texture.image
          );
        } else if (glFormat !== undefined && glType !== undefined) {
          this._gl.texImage2D(
            this._gl.TEXTURE_2D,
            0,
            glFormat,
            glFormat,
            glType,
            texture.image
          );
        }
      } else {
        //3D
        this._setTextureParameters(
          /** @type {WebGL2RenderingContext} */ (this._gl).TEXTURE_3D,
          texture
        );
        /** @type {WebGL2RenderingContext} */ (this._gl).texImage3D(
          /** @type {WebGL2RenderingContext} */ (this._gl).TEXTURE_3D,
          0,
          /** @type {WebGL2RenderingContext} */ (this._gl).R32F,
          texture.image.size.z,
          texture.image.size.y,
          texture.image.size.x,
          0,
          /** @type {WebGL2RenderingContext} */ (this._gl).RED,
          this._gl.FLOAT,
          texture.image.data
        );
      }

      texture.needsUpdate = false;

      if (texture.onUpdate) texture.onUpdate();
    } else {
      this._gl.activeTexture(this._gl.TEXTURE0 + slot);
      if (is3D)
        this._gl.bindTexture(
          /** @type {WebGL2RenderingContext} */ (this._gl).TEXTURE_3D,
          texture.__webglTexture
        );
      else this._gl.bindTexture(this._gl.TEXTURE_2D, texture.__webglTexture);
    }
  }

  /**
   * @param {Scene} scene
   */
  initWebGLObjects(scene) {
    if (!scene.__webglObjects) {
      scene.__webglObjects = [];
      scene.__webglObjectsImmediate = [];
      scene.__webglSprites = [];
      scene.__webglFlares = [];
    }

    // Add objects; this sets up buffers for each geometryGroup
    if (scene.__objectsAdded.length) {
      while (scene.__objectsAdded.length) {
        this._addObject(scene.__objectsAdded[0], scene);
        scene.__objectsAdded.splice(0, 1);
      }

      // Force buffer update during render
      // Hackish fix for initial cartoon-render-then-transparent-surface
      // bug
      this._currentGeometryGroupHash = -1;
    }

    while (scene.__objectsRemoved.length) {
      this._removeObject(scene.__objectsRemoved[0], scene);
      scene.__objectsRemoved.splice(0, 1);
    }

    // update must be called after objects adding / removal
    // This sends typed arrays to GL buffers for each geometryGroup
    for (var o = 0, ol = scene.__webglObjects.length; o < ol; o++) {
      this._updateObject(scene.__webglObjects[o].object);
    }
  }

  /**
   * @param {Camera} camera
   * @param {Light[]} lights
   * @param {Fog} fog
   * @param {any} material
   * @param {GeometryGroup} geometryGroup
   * @param {any} object
   */
  renderBuffer(camera, lights, fog, material, geometryGroup, object) {
    if (!material.visible || !this._gl) return;

    var program, attributes;

    // Sets up proper vertex and fragment shaders and attaches them to webGL
    // program
    // Also sets appropriate uniform variables
    program = this._setProgram(camera, lights, fog, material, object, this);

    attributes = program.attributes;

    var updateBuffers = false,
      wireframeBit = material.wireframe ? 1 : 0,
      geometryGroupHash =
        geometryGroup.id * 0xffffff + program.id * 2 + wireframeBit;

    if (geometryGroupHash !== this._currentGeometryGroupHash) {
      this._currentGeometryGroupHash = geometryGroupHash;
      updateBuffers = true;
    }

    // rebind shader attributes to appropriate (and already initialized) gl
    // buffers
    if (updateBuffers) {
      this._disableAttributes();

      // Vertices
      if (attributes.position >= 0) {
        this._gl.bindBuffer(
          this._gl.ARRAY_BUFFER,
          geometryGroup.__webglVertexBuffer
        );
        this._enableAttribute(attributes.position);
        this._gl.vertexAttribPointer(
          attributes.position,
          3,
          this._gl.FLOAT,
          false,
          0,
          0
        );
      }

      // Colors
      if (attributes.color >= 0) {
        this._gl.bindBuffer(
          this._gl.ARRAY_BUFFER,
          geometryGroup.__webglColorBuffer
        );
        this._enableAttribute(attributes.color);
        this._gl.vertexAttribPointer(
          attributes.color,
          3,
          this._gl.FLOAT,
          false,
          0,
          0
        );
      }

      // Normals
      if (attributes.normal >= 0) {
        this._gl.bindBuffer(
          this._gl.ARRAY_BUFFER,
          geometryGroup.__webglNormalBuffer
        );
        this._enableAttribute(attributes.normal);
        this._gl.vertexAttribPointer(
          attributes.normal,
          3,
          this._gl.FLOAT,
          false,
          0,
          0
        );
      }

      // Offsets (Instanced only)
      if (attributes.offset >= 0) {
        this._gl.bindBuffer(
          this._gl.ARRAY_BUFFER,
          geometryGroup.__webglOffsetBuffer
        );
        this._enableAttribute(attributes.offset);
        this._gl.vertexAttribPointer(
          attributes.offset,
          3,
          this._gl.FLOAT,
          false,
          0,
          0
        );
      }

      // Radii (Instanced only)
      if (attributes.radius >= 0) {
        this._gl.bindBuffer(
          this._gl.ARRAY_BUFFER,
          geometryGroup.__webglRadiusBuffer
        );
        this._enableAttribute(attributes.radius);
        this._gl.vertexAttribPointer(
          attributes.radius,
          1,
          this._gl.FLOAT,
          false,
          0,
          0
        );
      }
    }

    // Render
    var faceCount, lineCount;
    // lambert shaders - draw triangles
    // TODO: make sure geometryGroup's face count is setup correctly
    if (object instanceof Mesh) {
      if (material.shaderID === "instanced" && this._extInstanced) {
        var sphereGeometryGroup = material.sphere.geometryGroups[0];
        if (updateBuffers && this._gl) {
          this._gl.bindBuffer(
            this._gl.ARRAY_BUFFER,
            geometryGroup.__webglVertexBuffer
          );
          this._gl.bufferData(
            this._gl.ARRAY_BUFFER,
            sphereGeometryGroup.vertexArray,
            this._gl.STATIC_DRAW
          );
          this._gl.bindBuffer(
            this._gl.ARRAY_BUFFER,
            geometryGroup.__webglNormalBuffer
          );
          this._gl.bufferData(
            this._gl.ARRAY_BUFFER,
            sphereGeometryGroup.normalArray,
            this._gl.STATIC_DRAW
          );
          this._gl.bindBuffer(
            this._gl.ELEMENT_ARRAY_BUFFER,
            geometryGroup.__webglFaceBuffer
          );
          this._gl.bufferData(
            this._gl.ELEMENT_ARRAY_BUFFER,
            sphereGeometryGroup.faceArray,
            this._gl.STATIC_DRAW
          );
        }

        faceCount = sphereGeometryGroup.faceidx;

        this._extInstanced.vertexAttribDivisorANGLE(attributes.offset, 1);
        this._extInstanced.vertexAttribDivisorANGLE(attributes.radius, 1);
        this._extInstanced.vertexAttribDivisorANGLE(attributes.color, 1);

        this._extInstanced.drawElementsInstancedANGLE(
          this._gl.TRIANGLES,
          faceCount,
          this._gl.UNSIGNED_SHORT,
          0,
          geometryGroup.radiusArray.length
        );

        this._extInstanced.vertexAttribDivisorANGLE(attributes.offset, 0);
        this._extInstanced.vertexAttribDivisorANGLE(attributes.radius, 0);
        this._extInstanced.vertexAttribDivisorANGLE(attributes.color, 0);
      } else if (material.wireframe) {
        lineCount = geometryGroup.lineidx;
        this._setLineWidth(material.wireframeLinewidth);

        if (updateBuffers)
          this._gl.bindBuffer(
            this._gl.ELEMENT_ARRAY_BUFFER,
            geometryGroup.__webglLineBuffer
          );

        this._gl.drawElements(
          this._gl.LINES,
          lineCount,
          this._gl.UNSIGNED_SHORT,
          0
        );
      } else {
        faceCount = geometryGroup.faceidx;

        if (updateBuffers)
          this._gl.bindBuffer(
            this._gl.ELEMENT_ARRAY_BUFFER,
            geometryGroup.__webglFaceBuffer
          );
        this._gl.drawElements(
          this._gl.TRIANGLES,
          faceCount,
          this._gl.UNSIGNED_SHORT,
          0
        );
      }

      this.info.render.calls++;
      this.info.render.vertices += faceCount;
      this.info.render.faces += faceCount / 3;
    }

    // basic shaders - draw lines
    else if (object instanceof Line) {
      lineCount = geometryGroup.vertices;

      this._setLineWidth(material.linewidth);
      this._gl.drawArrays(this._gl.LINES, 0, lineCount);

      this.info.render.calls++;
    }
  }

  /**
   * @param {Camera} camera
   * @param {Light[]} lights
   * @param {Fog} fog
   * @param {any} material
   * @param {{ _modelViewMatrix: Matrix4; _normalMatrix: { elements: any; }; material: { texmatrix: Matrix4; unit: number; maxdepth: number; transfermax: any; transfermin: any; subsamples: any; transferfn: any; map: any; }; }} object
   * @param {this} renderer
   */
  _setProgram(camera, lights, fog, material, object, renderer) {
    if (!this._gl) return;
    if (material.needsUpdate) {
      if (material.program) this._deallocateMaterial(material);

      // @ts-ignore
      this.initMaterial(material, lights, fog, object);
      material.needsUpdate = false;
    }

    var refreshMaterial = false;

    // p_uniforms: uniformVarName => uniformLocation
    // m_uniforms: uniformVarName => uniformJsVal
    var program = material.program,
      p_uniforms = program.uniforms,
      m_uniforms = material.uniforms;

    if (program != this._currentProgram) {
      this._gl.useProgram(program);
      this._currentProgram = program;

      refreshMaterial = true;
    }

    if (material.id != this._currentMaterialId) {
      this._currentMaterialId = material.id;
      refreshMaterial = true;
    }

    if (camera != this._currentCamera) {
      this._currentCamera = camera;
      refreshMaterial = true;
    }

    this._gl.uniformMatrix4fv(
      p_uniforms.projectionMatrix,
      false,
      camera.projectionMatrix.elements
    );
    this._gl.uniformMatrix4fv(
      p_uniforms.modelViewMatrix,
      false,
      object._modelViewMatrix.elements
    );
    this._gl.uniformMatrix3fv(
      p_uniforms.normalMatrix,
      false,
      object._normalMatrix.elements
    );

    // Send projection matrix to uniform variable in shader
    if (refreshMaterial) {
      // Load projection, model-view matrices for perspective

      // Set up correct fog uniform vals
      m_uniforms.fogColor.value = fog.color;
      m_uniforms.fogNear.value = fog.near;
      m_uniforms.fogFar.value = fog.far;

      // Set up lights for lambert shader
      if (
        material.shaderID.startsWith("lambert") ||
        material.shaderID === "instanced" ||
        material.shaderID.endsWith("imposter")
      ) {
        // load view and normal matrices for directional and object
        // lighting
        this._gl.uniformMatrix4fv(
          p_uniforms.viewMatrix,
          false,
          camera.matrixWorldInverse.elements
        );

        if (this._lightsNeedUpdate) {
          this._setupLights(program, lights);
          this._lightsNeedUpdate = false;
        }

        // Set up correct light uniform var vals
        m_uniforms.directionalLightColor.value =
          this._lights.directional.colors;
        m_uniforms.directionalLightDirection.value =
          this._lights.directional.positions;
      } else if (material.shaderID.endsWith("outline")) {
        m_uniforms.outlineColor.value = material.outlineColor;
        m_uniforms.outlineWidth.value = material.outlineWidth;
        m_uniforms.outlinePushback.value = material.outlinePushback;
      } else if (material.shaderID === "volumetric" && this._worldInverse) {
        //need a matrix that maps back from model coordinates to texture coordinates
        //  textureMat*modelInv*position
        object._modelViewMatrix.getScale(this._direction); //scale factor of conversion
        this._worldInverse.getInverse(object._modelViewMatrix);
        this._projInverse.getInverse(camera.projectionMatrix);
        this._textureMatrix.multiplyMatrices(
          object.material.texmatrix,
          this._worldInverse
        );
        this._gl.uniformMatrix4fv(
          p_uniforms.textmat,
          false,
          this._textureMatrix.elements
        );
        this._gl.uniformMatrix4fv(
          p_uniforms.projinv,
          false,
          this._projInverse.elements
        );

        //  need the resolution (step size of ray in viewer coordinates)
        let invscale = Math.min(
          Math.min(this._direction.x, this._direction.y),
          this._direction.z
        );
        m_uniforms.step.value = object.material.unit * invscale;
        m_uniforms.maxdepth.value = object.material.maxdepth * invscale;
        m_uniforms.transfermax.value = object.material.transfermax;
        m_uniforms.transfermin.value = object.material.transfermin;
        m_uniforms.subsamples.value = object.material.subsamples;

        renderer.setTexture(object.material.transferfn, 4, false);
        renderer.setTexture(object.material.map, 3, true);
        //depth texture from the renderbuffer, for volumetric integration with surfaces
        this._gl.activeTexture(this._gl.TEXTURE5);
        this._gl.bindTexture(this._gl.TEXTURE_2D, this._depthTexture);
      }

      // opacity, diffuse, emissive, etc
      m_uniforms.opacity.value = material.opacity;

      // Load any other material specific uniform variables to gl shaders
      this._loadMaterialUniforms(p_uniforms, m_uniforms);
    }

    return program;
  }

  /**
   * @param {{ [x: string]: any; }} p_uniforms
   * @param {{ [x: string]: { value: any; type: any; }; }} m_uniforms
   */
  _loadMaterialUniforms(p_uniforms, m_uniforms) {
    if (!this._gl) return;
    var uniformVar, type, uniformVal, uniformLoc;

    for (uniformVar in m_uniforms) {
      if (!p_uniforms[uniformVar]) continue;

      type = m_uniforms[uniformVar].type;
      uniformVal = m_uniforms[uniformVar].value;
      uniformLoc = p_uniforms[uniformVar];

      // single float
      if (type === "f") this._gl.uniform1f(uniformLoc, uniformVal);
      // single integer
      else if (type === "i") this._gl.uniform1i(uniformLoc, uniformVal);
      // array of floats
      else if (type === "fv") this._gl.uniform3fv(uniformLoc, uniformVal);
      // color - r,g,b floats
      else if (type === "c")
        this._gl.uniform3f(
          uniformLoc,
          uniformVal.r,
          uniformVal.g,
          uniformVal.b
        );
      else if (type === "f4")
        this._gl.uniform4f(
          uniformLoc,
          uniformVal[0],
          uniformVal[1],
          uniformVal[2],
          uniformVal[3]
        );
    }
  }

  _isWebGL1() {
    return this._webglversion == 1;
  }

  /**
   * @param {number} attribute
   */
  _enableAttribute(attribute) {
    if (!this._enabledAttributes[attribute] && this._gl) {
      this._gl.enableVertexAttribArray(attribute);
      this._enabledAttributes[attribute] = true;
    }
  }

  _disableAttributes() {
    if (!this._gl) return;
    for (let attribute in this._enabledAttributes) {
      if (this._enabledAttributes[attribute]) {
        this._gl.disableVertexAttribArray(Number.parseInt(attribute));
        this._enabledAttributes[attribute] = false;
      }
    }
  }

  /**
   * @param {null} polygonOffset
   */
  _setPolygonOffset(polygonOffset) {
    if (!this._gl) return;
    if (this._oldPolygonOffset !== polygonOffset) {
      if (polygonOffset) this._gl.enable(this._gl.POLYGON_OFFSET_FILL);
      else this._gl.disable(this._gl.POLYGON_OFFSET_FILL);
    }
  }

  /**
   * @param {null} width
   */
  _setLineWidth(width) {
    if (!this._gl) return;
    if (width !== this._oldLineWidth) {
      this._gl.lineWidth(width);
      this._oldLineWidth = width;
    }
  }

  /**
   * @param {Geometry} geometry
   */
  _deallocateGeometry(geometry) {
    if (!this._gl) return;
    geometry.__webglInit = undefined;

    if (geometry.__webglVertexBuffer !== undefined)
      this._gl.deleteBuffer(geometry.__webglVertexBuffer);

    if (geometry.__webglColorBuffer !== undefined)
      this._gl.deleteBuffer(geometry.__webglColorBuffer);

    if (geometry.geometryGroups !== undefined) {
      for (var g = 0, gl = geometry.groups; g < gl; g++) {
        var geometryGroup = geometry.geometryGroups[g];

        if (geometryGroup.__webglVertexBuffer !== undefined)
          this._gl.deleteBuffer(geometryGroup.__webglVertexBuffer);

        if (geometryGroup.__webglColorBuffer !== undefined)
          this._gl.deleteBuffer(geometryGroup.__webglColorBuffer);

        if (geometryGroup.__webglNormalBuffer !== undefined)
          this._gl.deleteBuffer(geometryGroup.__webglNormalBuffer);

        if (geometryGroup.__webglFaceBuffer !== undefined)
          this._gl.deleteBuffer(geometryGroup.__webglFaceBuffer);

        if (geometryGroup.__webglLineBuffer !== undefined)
          this._gl.deleteBuffer(geometryGroup.__webglLineBuffer);
      }
    }
  }

  /**
   * @param {any} material
   */
  _deallocateMaterial(material) {
    if (!this._gl) return;
    var program = material.program;

    if (program === undefined) return;

    material.program = undefined;

    // only deallocate GL program if this was the last use of shared program
    // assumed there is only single copy of any program in the _programs
    // list
    // (that's how it's constructed)

    var i, il, programInfo;
    var deleteProgram = false;

    for (i = 0, il = this._programs.length; i < il; i++) {
      programInfo = this._programs[i];

      if (programInfo.program === program) {
        programInfo.usedTimes--;

        if (programInfo.usedTimes === 0) {
          deleteProgram = true;
        }

        break;
      }
    }

    if (deleteProgram === true) {
      // avoid using array.splice, this is costlier than creating new
      // array from scratch

      var newPrograms = [];

      for (i = 0, il = this._programs.length; i < il; i++) {
        programInfo = this._programs[i];

        if (programInfo.program !== program) {
          newPrograms.push(programInfo);
        }
      }

      this._programs = newPrograms;
      this._gl.deleteProgram(program);
      this.info.memory.programs--;
    }
  }

  /**
   * @param {{ image: { __webglTextureCube: any; }; __webglInit: boolean; __webglTexture: any; }} texture
   */
  _deallocateTexture(texture) {
    if (!this._gl) return;
    if (texture.image && texture.image.__webglTextureCube) {
      // cube texture

      this._gl.deleteTexture(texture.image.__webglTextureCube);
    } else {
      // 2D texture

      if (!texture.__webglInit) return;

      texture.__webglInit = false;
      this._gl.deleteTexture(texture.__webglTexture);
    }
  }

  /**
   * @param {{ target: any; }} event
   */
  _onGeometryDispose(event) {
    var geometry = event.target;
    geometry.removeEventListener("dispose", this._onGeometryDispose);

    this._deallocateGeometry(geometry);

    this.info.memory.geometries--;
  }

  /**
   * @param {{ target: any; }} event
   */
  _onTextureDispose(event) {
    var texture = event.target;

    texture.removeEventListener("dispose", this._onTextureDispose);

    this._deallocateTexture(texture);

    this.info.memory.textures--;
  }

  /**
   * @param {{ target: any; }} event
   */
  _onMaterialDispose(event) {
    var material = event.target;
    material.removeEventListener("dispose", this._onMaterialDispose);

    this._deallocateMaterial(material);
  }

  // Compile and return shader
  /**
   * @param {string} type
   * @param {string} str
   */
  _getShader(type, str) {
    if (!this._gl) return;
    var shader;

    if (!this._isWebGL1() && !str.startsWith("_version")) {
      //convert webgl1 to webgl2, unless a version is already explicit
      str = str.replace(/gl_FragDepthEXT/g, "gl_FragDepth");
      if (type == "fragment") {
        str = str.replace(/varying/g, "in");
      } else {
        str = str.replace(/varying/g, "out");
      }
      str = str.replace(/attribute/g, "in");
      str = str.replace(/texture2D/g, "texture");
      str = str.replace(/\/\/DEFINEFRAGCOLOR/g, "out vec4 glFragColor;");
      str = str.replace(/gl_FragColor/g, "glFragColor");
      str = "_version 300 es\n" + str;
    }
    if (type === "fragment")
      shader = this._gl.createShader(this._gl.FRAGMENT_SHADER);
    else if (type === "vertex")
      shader = this._gl.createShader(this._gl.VERTEX_SHADER);

    this._gl.shaderSource(/** @type {WebGLShader} */ (shader), str);
    this._gl.compileShader(/** @type {WebGLShader} */ (shader));

    if (
      !this._gl.getShaderParameter(
        /** @type {WebGLShader} */ (shader),
        this._gl.COMPILE_STATUS
      )
    ) {
      console.error(
        this._gl.getShaderInfoLog(/** @type {WebGLShader} */ (shader))
      );
      console.error("could not initialize shader");
      return null;
    }

    return shader;
  }

  // Compile appropriate shaders (if necessary) from source code and attach to
  // gl program.
  /**
   * @param {string} fragmentShader
   * @param {string} vertexShader
   * @param {{}} uniforms
   * @param {{ [x: string]: any; wireframe?: any; fragdepth?: any; volumetric?: any; }} parameters
   */
  _buildProgram(fragmentShader, vertexShader, uniforms, parameters) {
    if (!this._gl) return;
    var p, pl, program, code;
    var chunks = [];

    chunks.push(fragmentShader);
    chunks.push(vertexShader);

    for (p in parameters) {
      chunks.push(p);
      chunks.push(parameters[p]);
    }

    code = chunks.join();

    // check if program has already been compiled

    for (p = 0, pl = this._programs.length; p < pl; p++) {
      var programInfo = this._programs[p];

      if (programInfo.code === code) {
        programInfo.usedTimes++;

        return programInfo.program;
      }
    }

    // check if program requires webgl2
    if (this._isWebGL1()) {
      if (parameters.volumetric)
        throw new Error(
          "Volumetric rendering requires webgl2 which is not supported by your hardware."
        );
    }

    // Set up new program and compile shaders

    program = this._gl.createProgram();
    if (!program) return null;
    // set up precision
    var precision = this._precision;
    var prefix = "precision " + precision + " float;";

    var prefix_vertex = [
      parameters.volumetric ? "_version 300 es" : "",
      prefix,
    ].join("\n");

    var prefix_fragment = [
      parameters.volumetric ? "_version 300 es" : "",
      parameters.fragdepth && this._isWebGL1()
        ? "_extension GL_EXT_frag_depth: enable"
        : "",
      parameters.wireframe ? "_define WIREFRAME 1" : "",
      prefix,
    ].join("\n");

    var glFragmentShader = this._getShader(
      "fragment",
      prefix_fragment + fragmentShader
    );
    var glVertexShader = this._getShader(
      "vertex",
      prefix_vertex + vertexShader
    );

    if (glVertexShader) this._gl.attachShader(program, glVertexShader);
    if (glFragmentShader) this._gl.attachShader(program, glFragmentShader);
    this._gl.linkProgram(program);

    if (!this._gl.getProgramParameter(program, this._gl.LINK_STATUS))
      console.error("Could not initialize shader");

    // gather and cache uniform variables and attributes

    // @ts-ignore
    program.uniforms = {};
    // @ts-ignore
    program.attributes = {};

    var identifiers, u, i;

    // uniform vars
    identifiers = [
      "viewMatrix",
      "modelViewMatrix",
      "projectionMatrix",
      "normalMatrix",
    ];

    // custom uniform vars
    for (u in uniforms) identifiers.push(u);

    for (i = 0; i < identifiers.length; i++) {
      var uniformVar = identifiers[i];
      // @ts-ignore
      program.uniforms[uniformVar] = this._gl.getUniformLocation(
        program,
        uniformVar
      );
    }

    // attributes
    identifiers = [
      "position",
      "normal",
      "color",
      "lineDistance",
      "offset",
      "radius",
    ];

    /*
     * for (a in attributes) identifiers.push(a);
     */

    for (i = 0; i < identifiers.length; i++) {
      var attributeVar = identifiers[i];
      // @ts-ignore
      program.attributes[attributeVar] = this._gl.getAttribLocation(
        program,
        attributeVar
      );
    }

    // @ts-ignore
    program.id = this._programs_counter++;
    this._programs.push({
      program: program,
      code: code,
      usedTimes: 1,
    });
    this.info.memory.programs = this._programs.length;

    return program;
  }

  // rendering
  /**
   * @param {string | any[]} renderList
   * @param {boolean} reverse
   * @param {string} materialType
   * @param {any} camera
   * @param {any} lights
   * @param {any} fog
   * @param {boolean} useBlending
   */
  _renderObjects(
    renderList,
    reverse,
    materialType,
    camera,
    lights,
    fog,
    useBlending
  ) {
    var webglObject, object, buffer, material, start, end, delta;

    // Forward or backward render

    if (reverse) {
      start = renderList.length - 1;
      end = -1;
      delta = -1;
    } else {
      start = 0;
      end = renderList.length;
      delta = 1;
    }

    for (var i = start; i !== end; i += delta) {
      webglObject = renderList[i];

      if (webglObject.render) {
        object = webglObject.object;
        buffer = webglObject.buffer;
        material = webglObject[materialType];

        if (!material) continue;

        if (useBlending) this.setBlending(true);

        this.setDepthTest(material.depthTest);
        this.setDepthWrite(material.depthWrite);
        this._setPolygonOffset(
          material.polygonOffset,
          // @ts-ignore
          material.polygonOffsetFactor,
          material.polygonOffsetUnits
        );

        var reflected = object._modelViewMatrix.isReflected();

        this.setMaterialFaces(material, reflected);
        this.renderBuffer(camera, lights, fog, material, buffer, object);
        if (this._outlineEnabled || material.outline) {
          if (material.shaderID == "sphereimposter") {
            this.renderBuffer(
              camera,
              lights,
              fog,
              this._outlineSphereImposterMaterial,
              buffer,
              object
            );
          } else if (material.shaderID == "stickimposter") {
            this.renderBuffer(
              camera,
              lights,
              fog,
              this._outlineStickImposterMaterial,
              buffer,
              object
            );
          } else if (
            !material.wireframe &&
            material.shaderID !== "basic" &&
            material.opacity !== 0.0
          ) {
            this.renderBuffer(
              camera,
              lights,
              fog,
              this._outlineMaterial,
              buffer,
              object
            );
          }
        }
      }
    }
  }

  /**
   * @param {Scene} scene
   * @param {Camera} camera
   * @param {boolean} inFront
   */
  _renderSprites(scene, camera, inFront) {
    // Reset state once regardless
    // This should also fix cartoon render bug (after transparent surface
    // render)

    this._currentGeometryGroupHash = -1;
    this._currentProgram = null;
    this._currentCamera = null;
    this._oldBlending = -1;
    this._oldDepthWrite = -1;
    this._oldDepthTest = -1;
    this._oldDoubleSided = -1;
    this._currentMaterialId = -1;
    this._oldFlipSided = -1;
    this._lightsNeedUpdate = true;

    this.sprites.render(
      scene,
      camera,
      this._currentWidth,
      this._currentHeight,
      inFront
    );

    // Reset state a
    this._currentGeometryGroupHash = -1;
    this._currentProgram = null;
    this._currentCamera = null;
    this._oldBlending = -1;
    this._oldDepthWrite = -1;
    this._oldDepthTest = -1;
    this._oldDoubleSided = -1;
    this._currentMaterialId = -1;
    this._oldFlipSided = -1;
  }

  /**
   * @param {{ __webglInit: boolean; _modelViewMatrix: Matrix4; _normalMatrix: Matrix3; geometry: { __webglInit: boolean | undefined; addEventListener: (arg0: string, arg1: (event: any) => void) => void; } | undefined; material: any; __webglActive: boolean; }} object
   * @param {{ __webglObjects: any; __webglSprites: Sprite[]; }} scene
   */
  _addObject(object, scene) {
    // @ts-ignore
    var g, gl, geometry, material, geometryGroup;

    if (!object.__webglInit) {
      object.__webglInit = true;

      object._modelViewMatrix = new Matrix4();
      object._normalMatrix = new Matrix3();

      if (
        object.geometry !== undefined &&
        object.geometry.__webglInit === undefined
      ) {
        object.geometry.__webglInit = true;
        object.geometry.addEventListener("dispose", this._onGeometryDispose);
      }

      if (object instanceof Mesh || object instanceof Line) {
        geometry = object.geometry;
        material = object.material;

        for (g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
          geometryGroup = geometry.geometryGroups[g];

          geometryGroup.id = this._geometryGroupCounter++;

          // initialise VBO on the first access

          if (!geometryGroup.__webglVertexBuffer) {
            if (object instanceof Mesh) {
              this._createMeshBuffers(geometryGroup);
              geometry.elementsNeedUpdate = true;
              geometry.normalsNeedUpdate = true;
            } else if (object instanceof Line)
              this._createLineBuffers(geometryGroup);

            geometry.verticesNeedUpdate = true;
            geometry.colorsNeedUpdate = true;
          }
        }
      }
    }

    if (!object.__webglActive) {
      if (object instanceof Mesh || object instanceof Line) {
        geometry = object.geometry;

        for (g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
          geometryGroup = geometry.geometryGroups[g];

          this._addBuffer(scene.__webglObjects, geometryGroup, object);
        }
      }

      // Sprite
      else if (object instanceof Sprite) scene.__webglSprites.push(object);

      object.__webglActive = true;
    }
  }

  /**
   * @param {{ geometry: any; }} object
   */
  _updateObject(object) {
    if (!this._gl) return;
    var geometry = object.geometry,
      geometryGroup;

    if (object instanceof Mesh || object instanceof Line) {
      for (var g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
        geometryGroup = geometry.geometryGroups[g];

        if (
          geometry.verticesNeedUpdate ||
          geometry.elementsNeedUpdate ||
          geometry.colorsNeedUpdate ||
          geometry.normalsNeedUpdate
        ) {
          this._setBuffers(geometryGroup, this._gl.STATIC_DRAW);
        }
      }

      geometry.verticesNeedUpdate = false;
      geometry.elementsNeedUpdate = false;
      geometry.normalsNeedUpdate = false;
      geometry.colorsNeedUpdate = false;

      geometry.buffersNeedUpdate = false;
    }
  }

  /**
   * @param {{ __webglActive: boolean; }} object
   * @param {{ __webglObjects: any; __webglSprites: any; }} scene
   */
  _removeObject(object, scene) {
    if (object instanceof Mesh || object instanceof Line)
      this._removeInstances(scene.__webglObjects, object);
    else if (object instanceof Sprite)
      this._removeInstancesDirect(scene.__webglSprites, object);

    object.__webglActive = false;
  }

  /**
   * @param {any[]} objList
   * @param {Mesh | Line} object
   */
  _removeInstances(objList, object) {
    for (var o = objList.length - 1; o >= 0; --o) {
      if (objList[o].object === object) objList.splice(o, 1);
    }
  }

  /**
   * @param {any[]} objList
   * @param {Sprite} object
   */
  _removeInstancesDirect(objList, object) {
    for (var o = objList.length - 1; o >= 0; --o) {
      if (objList[o] === object) objList.splice(o, 1);
    }
  }

  /**
   * @param {{ object: any; opaque: null; transparent: null; volumetric: null; blank: any; }} globject
   */
  _unrollBufferMaterial(globject) {
    var object = globject.object;
    var material = object.material;

    if (material.volumetric) {
      globject.opaque = null;
      globject.transparent = null;
      globject.volumetric = material;
    } else if (material.transparent) {
      globject.opaque = null;
      globject.volumetric = null;
      globject.transparent = material;
      if (!material.wireframe) {
        var blankMaterial = material.clone();
        blankMaterial.opacity = 0.0;
        globject.blank = blankMaterial;
      }
    } else {
      globject.opaque = material;
      globject.transparent = null;
      globject.volumetric = null;
    }
  }

  /**
   * @param {{ vertexArray: any; colorArray: any; __webglOffsetBuffer: undefined; __webglVertexBuffer: any; __webglColorBuffer: any; normalArray: any; __webglNormalBuffer: undefined; radiusArray: any; __webglRadiusBuffer: undefined; faceArray: any; __webglFaceBuffer: undefined; lineArray: any; __webglLineBuffer: undefined; }} geometryGroup
   * @param {any} hint
   */
  _setBuffers(geometryGroup, hint) {
    if (!this._gl) return;
    var vertexArray = geometryGroup.vertexArray;
    var colorArray = geometryGroup.colorArray;

    // offset buffers
    if (geometryGroup.__webglOffsetBuffer !== undefined) {
      this._gl.bindBuffer(
        this._gl.ARRAY_BUFFER,
        geometryGroup.__webglOffsetBuffer
      );
      this._gl.bufferData(this._gl.ARRAY_BUFFER, vertexArray, hint);
    } else {
      //normal, non-instanced case
      this._gl.bindBuffer(
        this._gl.ARRAY_BUFFER,
        geometryGroup.__webglVertexBuffer
      );
      this._gl.bufferData(this._gl.ARRAY_BUFFER, vertexArray, hint);
    }
    // color buffers
    this._gl.bindBuffer(
      this._gl.ARRAY_BUFFER,
      geometryGroup.__webglColorBuffer
    );
    this._gl.bufferData(this._gl.ARRAY_BUFFER, colorArray, hint);

    // normal buffers
    if (
      geometryGroup.normalArray &&
      geometryGroup.__webglNormalBuffer !== undefined
    ) {
      var normalArray = geometryGroup.normalArray;
      this._gl.bindBuffer(
        this._gl.ARRAY_BUFFER,
        geometryGroup.__webglNormalBuffer
      );
      this._gl.bufferData(this._gl.ARRAY_BUFFER, normalArray, hint);
    }

    // radius buffers
    if (
      geometryGroup.radiusArray &&
      geometryGroup.__webglRadiusBuffer !== undefined
    ) {
      this._gl.bindBuffer(
        this._gl.ARRAY_BUFFER,
        geometryGroup.__webglRadiusBuffer
      );
      this._gl.bufferData(
        this._gl.ARRAY_BUFFER,
        geometryGroup.radiusArray,
        hint
      );
    }

    // face (index) buffers
    if (
      geometryGroup.faceArray &&
      geometryGroup.__webglFaceBuffer !== undefined
    ) {
      var faceArray = geometryGroup.faceArray;
      this._gl.bindBuffer(
        this._gl.ELEMENT_ARRAY_BUFFER,
        geometryGroup.__webglFaceBuffer
      );
      this._gl.bufferData(this._gl.ELEMENT_ARRAY_BUFFER, faceArray, hint);
    }

    // line (index) buffers (for wireframe)
    if (
      geometryGroup.lineArray &&
      geometryGroup.__webglLineBuffer !== undefined
    ) {
      var lineArray = geometryGroup.lineArray;
      this._gl.bindBuffer(
        this._gl.ELEMENT_ARRAY_BUFFER,
        geometryGroup.__webglLineBuffer
      );
      this._gl.bufferData(this._gl.ELEMENT_ARRAY_BUFFER, lineArray, hint);
    }
  }

  // Creates appropriate gl buffers for geometry chunk
  // TODO: do we need line buffer for mesh objects?
  // Also, can we integrate this with _createLineBuffers?
  /**
   * @param {{ radiusArray: any; __webglRadiusBuffer: any; useOffset: any; __webglOffsetBuffer: any; __webglVertexBuffer: any; __webglNormalBuffer: any; __webglColorBuffer: any; __webglFaceBuffer: any; __webglLineBuffer: any; }} geometryGroup
   */
  _createMeshBuffers(geometryGroup) {
    if (!this._gl) return;
    if (geometryGroup.radiusArray) {
      geometryGroup.__webglRadiusBuffer = this._gl.createBuffer();
    }
    if (geometryGroup.useOffset) {
      geometryGroup.__webglOffsetBuffer = this._gl.createBuffer();
    }
    geometryGroup.__webglVertexBuffer = this._gl.createBuffer();
    geometryGroup.__webglNormalBuffer = this._gl.createBuffer();
    geometryGroup.__webglColorBuffer = this._gl.createBuffer();

    geometryGroup.__webglFaceBuffer = this._gl.createBuffer();
    geometryGroup.__webglLineBuffer = this._gl.createBuffer();

    this.info.memory.geometries++;
  }

  /**
   * @param {{ __webglVertexBuffer: any; __webglColorBuffer: any; }} geometry
   */
  _createLineBuffers(geometry) {
    if (!this._gl) return;
    geometry.__webglVertexBuffer = this._gl.createBuffer();
    geometry.__webglColorBuffer = this._gl.createBuffer();

    this.info.memory.geometries++;
  }

  /**
   * @param {{ buffer: any; object: any; opaque: null; transparent: null; }[]} objlist
   * @param {any} buffer
   * @param {Mesh | Line} object
   */
  _addBuffer(objlist, buffer, object) {
    objlist.push({
      buffer: buffer,
      object: object,
      opaque: null,
      transparent: null,
    });
  }

  /**
   * @param {{ _modelViewMatrix: { multiplyMatrices: (arg0: any, arg1: any) => void; }; matrixWorld: any; _normalMatrix: { getInverse: (arg0: any) => void; transpose: () => void; }; }} object
   * @param {{ matrixWorldInverse: any; }} camera
   */
  _setupMatrices(object, camera) {
    object._modelViewMatrix.multiplyMatrices(
      camera.matrixWorldInverse,
      object.matrixWorld
    );

    object._normalMatrix.getInverse(object._modelViewMatrix);
    object._normalMatrix.transpose();
  }

  // Fallback filters for non-power-of-2 textures
  _filterFallback() {
    if (!this._gl) return;
    return this._gl.LINEAR;
  }

  /**
   * @param {any} textureType
   * @param {{ magFilter: any; minFilter: any; }} texture
   */
  _setTextureParameters(textureType, texture) {
    if (!this._gl) return;
    if (textureType == this._gl.TEXTURE_2D) {
      this._gl.texParameteri(
        textureType,
        this._gl.TEXTURE_WRAP_S,
        this._gl.CLAMP_TO_EDGE
      );
      this._gl.texParameteri(
        textureType,
        this._gl.TEXTURE_WRAP_T,
        this._gl.CLAMP_TO_EDGE
      );
      this._gl.texParameteri(
        textureType,
        this._gl.TEXTURE_MAG_FILTER,
        // @ts-ignore
        this._filterFallback(texture.magFilter)
      );
      this._gl.texParameteri(
        textureType,
        this._gl.TEXTURE_MIN_FILTER,
        // @ts-ignore
        this._filterFallback(texture.minFilter)
      );
    } else {
      // 3Dtexture
      this._gl.texParameteri(
        textureType,
        this._gl.TEXTURE_WRAP_S,
        this._gl.CLAMP_TO_EDGE
      );
      this._gl.texParameteri(
        textureType,
        this._gl.TEXTURE_WRAP_T,
        this._gl.CLAMP_TO_EDGE
      );
      this._gl.texParameteri(
        textureType,
        /** @type {WebGL2RenderingContext} */ (this._gl).TEXTURE_WRAP_R,
        this._gl.CLAMP_TO_EDGE
      );

      // @ts-ignore
      if (this._extColorBufferFloat && this._extFloatLinear) {
        //linear interpolation isn't supported by default (despite being the default??)
        this._gl.texParameteri(
          textureType,
          this._gl.TEXTURE_MAG_FILTER,
          this._gl.LINEAR
        );
        this._gl.texParameteri(
          textureType,
          this._gl.TEXTURE_MIN_FILTER,
          this._gl.LINEAR
        );
      } else {
        this._gl.texParameteri(
          textureType,
          this._gl.TEXTURE_MAG_FILTER,
          this._gl.NEAREST
        );
        this._gl.texParameteri(
          textureType,
          this._gl.TEXTURE_MIN_FILTER,
          this._gl.NEAREST
        );
      }
    }
  }

  // Map constants to WebGL constants

  /**
   * @param {number} p
   */
  _paramToGL(p) {
    if (!this._gl) return;
    if (p === UnsignedByteType) return this._gl.UNSIGNED_BYTE;
    if (p === RGBAFormat) return this._gl.RGBA;
    if (p === NearestFilter) return this._gl.NEAREST;

    return 0;
  }

  // @ts-ignore
  _setupLights(program, lights) {
    /**
     * @type {number[]}
     */
    /**
     * @type {any[]}
     */
    var l,
      ll,
      light,
      r = 0,
      g = 0,
      b = 0,
      color,
      intensity,
      // @ts-ignore
      distance,
      zlights = this._lights,
      dirColors = zlights.directional.colors,
      dirPositions = zlights.directional.positions,
      // @ts-ignore
      dirCount = 0,
      dirLength = 0,
      dirOffset = 0;

    for (let l = 0, ll = lights.length; l < ll; l++) {
      light = lights[l];

      color = light.color;
      intensity = light.intensity;
      distance = light.distance;

      if (light instanceof Light) {
        dirCount++;

        this._direction.getPositionFromMatrix(light.matrixWorld);
        this._vector3.getPositionFromMatrix(light.target.matrixWorld);
        this._direction.sub(this._vector3);
        this._direction.normalize();

        if (
          this._direction.x === 0 &&
          this._direction.y === 0 &&
          this._direction.z === 0
        )
          continue;

        dirPositions[dirOffset] = this._direction.x;
        dirPositions[dirOffset + 1] = this._direction.y;
        dirPositions[dirOffset + 2] = this._direction.z;

        dirColors[dirOffset] = color.r * intensity;
        dirColors[dirOffset + 1] = color.g * intensity;
        dirColors[dirOffset + 2] = color.b * intensity;

        dirOffset += 3;

        dirLength++;
      }
    }

    zlights.ambient[0] = r;
    zlights.ambient[1] = g;
    zlights.ambient[2] = b;
    zlights.directional.length = dirLength;
  }
}
