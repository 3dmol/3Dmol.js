import { Fog } from "./Fog";
import type { Camera } from "./Camera";
import { Scene } from "./core";
//Render plugins go here
import { Sprite } from "./objects/Sprite";
import { ShaderLib } from "./shaders/index";
import { Shader } from './shaders/ShaderType';

export type Nullable<T> = T | null;
export type NullableWebGLUniformLocation = Nullable<WebGLUniformLocation>;

export type UniformsValueType = {
  uvOffset: NullableWebGLUniformLocation;
  uvScale: NullableWebGLUniformLocation;
  rotation: NullableWebGLUniformLocation;
  scale: NullableWebGLUniformLocation;
  alignment: NullableWebGLUniformLocation;
  color: NullableWebGLUniformLocation;
  map: NullableWebGLUniformLocation;
  opacity: NullableWebGLUniformLocation;
  useScreenCoordinates: NullableWebGLUniformLocation;
  screenPosition: NullableWebGLUniformLocation;
  modelViewMatrix: NullableWebGLUniformLocation;
  projectionMatrix: NullableWebGLUniformLocation;
  fogNear: NullableWebGLUniformLocation;
  fogFar: NullableWebGLUniformLocation;
  fogColor: NullableWebGLUniformLocation;
  alphaTest: NullableWebGLUniformLocation;
};

export type SpriteMeta = {
  vertices: Nullable<Float32Array>;
  faces: Nullable<Uint16Array>;
  vertexBuffer: Nullable<WebGLBuffer>;
  elementBuffer: Nullable<WebGLBuffer>;
  program: Nullable<WebGLProgram>;
  attributes: Record<string, number>;
  uniforms: Nullable<UniformsValueType>;
};

/*
 * Sprite render plugin
 */
export class SpritePlugin {
  private gl?: WebGLRenderingContext;
  private disposed = false;
  private readonly renderList: Sprite[] = [];
  private rendering = false;
  private renderer: any;
  private precision?: number;
  private sprite: SpriteMeta = {
    vertices: null,
    faces: null,
    vertexBuffer: null,
    elementBuffer: null,
    program: null,
    attributes: {},
    uniforms: null,
  };
  sprites?: Sprite[];

  init(renderer: any) {
    if (this.disposed || this.sprite.program) return;
    this.gl = renderer.context as WebGLRenderingContext;
    this.renderer = renderer;

    this.precision = renderer.getPrecision();

    this.sprite.vertices = new Float32Array(8 + 8);
    this.sprite.faces = new Uint16Array(6);

    var i = 0;

    this.sprite.vertices[i++] = -1;
    this.sprite.vertices[i++] = -1; // vertex 0
    this.sprite.vertices[i++] = 0;
    this.sprite.vertices[i++] = 0; // uv 0

    this.sprite.vertices[i++] = 1;
    this.sprite.vertices[i++] = -1; // vertex 1
    this.sprite.vertices[i++] = 1;
    this.sprite.vertices[i++] = 0; // uv 1

    this.sprite.vertices[i++] = 1;
    this.sprite.vertices[i++] = 1; // vertex 2
    this.sprite.vertices[i++] = 1;
    this.sprite.vertices[i++] = 1; // uv 2

    this.sprite.vertices[i++] = -1;
    this.sprite.vertices[i++] = 1; // vertex 3
    this.sprite.vertices[i++] = 0;
    this.sprite.vertices[i++] = 1; // uv 3

    i = 0;

    this.sprite.faces[i++] = 0;
    this.sprite.faces[i++] = 1;
    this.sprite.faces[i++] = 2;
    this.sprite.faces[i++] = 0;
    this.sprite.faces[i++] = 2;
    this.sprite.faces[i++] = 3;

    this.sprite.program = this.createProgram(
      ShaderLib.sprite,
      this.precision || 1 /** added default to single precision */
    );

    this.sprite.vertexBuffer = this.gl.createBuffer();
    this.sprite.elementBuffer = this.gl.createBuffer();

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.sprite.vertexBuffer);
    this.gl.bufferData(
      this.gl.ARRAY_BUFFER,
      this.sprite.vertices,
      this.gl.STATIC_DRAW
    );

    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.sprite.elementBuffer);
    this.gl.bufferData(
      this.gl.ELEMENT_ARRAY_BUFFER,
      this.sprite.faces,
      this.gl.STATIC_DRAW
    );

    this.sprite.attributes = {};

    this.sprite.attributes.position = this.gl.getAttribLocation(
      this.sprite.program,
      "position"
    );
    this.sprite.attributes.uv = this.gl.getAttribLocation(
      this.sprite.program,
      "uv"
    );

    this.sprite.uniforms = {
      uvOffset: this.gl.getUniformLocation(this.sprite.program, "uvOffset"),
      uvScale: this.gl.getUniformLocation(this.sprite.program, "uvScale"),
      rotation: this.gl.getUniformLocation(this.sprite.program, "rotation"),
      scale: this.gl.getUniformLocation(this.sprite.program, "scale"),
      alignment: this.gl.getUniformLocation(this.sprite.program, "alignment"),
      color: this.gl.getUniformLocation(this.sprite.program, "color"),
      map: this.gl.getUniformLocation(this.sprite.program, "map"),
      opacity: this.gl.getUniformLocation(this.sprite.program, "opacity"),
      useScreenCoordinates: this.gl.getUniformLocation(this.sprite.program, "useScreenCoordinates"),
      screenPosition: this.gl.getUniformLocation(this.sprite.program, "screenPosition"),
      modelViewMatrix: this.gl.getUniformLocation(this.sprite.program, "modelViewMatrix"),
      projectionMatrix: this.gl.getUniformLocation(this.sprite.program, "projectionMatrix"),
      fogNear: this.gl.getUniformLocation(this.sprite.program, "fogNear"),
      fogFar: this.gl.getUniformLocation(this.sprite.program, "fogFar"),
      fogColor: this.gl.getUniformLocation(this.sprite.program, "fogColor"),
      alphaTest: this.gl.getUniformLocation(this.sprite.program, "alphaTest"),
    };
  }

  render(
    scene: Scene,
    camera: Camera,
    viewportWidth: number,
    viewportHeight: number,
    inFront?: boolean
  ) {
    if (this.disposed) return;
    if (!this.gl) throw new Error("WebGLRenderer not initialized");
    // Texture update callbacks may render recursively; preserve the outer pass list.
    const nested = this.rendering;
    const sprites = nested ? [] : this.renderList;
    this.rendering = true;
    try {
      const sceneSprites = scene?.__webglSprites;
      for (let i = 0, n = sceneSprites?.length || 0; i < n; i++) {
        if (!(i in sceneSprites)) continue;
        const sprite = sceneSprites[i];
        // Keep pass membership and sorting order, including hidden sprites.
        if (inFront && sprite.material.depthTest == false) {
          sprites.push(sprite);
        } else if (!inFront && sprite.material.depthTest) {
          sprites.push(sprite);
        }
      }
      if (sprites.length) {
        this.drawSprites(sprites, scene, camera, viewportWidth, viewportHeight, inFront);
      }
    } finally {
      // Removed labels must not be retained between frames, including after errors.
      sprites.length = 0;
      this.rendering = nested;
    }
  }

  private drawSprites(
    sprites: Sprite[], scene: Scene, camera: Camera,
    viewportWidth: number, viewportHeight: number, inFront?: boolean
  ) {
    const nSprites = sprites.length;
    const attributes = this.sprite.attributes;
    const uniforms = this.sprite.uniforms;

    if (!uniforms) throw new Error("Uniforms not defined");

    var halfViewportWidth = viewportWidth * 0.5,
      halfViewportHeight = viewportHeight * 0.5;

    // setup gl

    this.gl.useProgram(this.sprite.program);

    this.gl.enableVertexAttribArray(attributes.position);
    this.gl.enableVertexAttribArray(attributes.uv);

    this.gl.disable(this.gl.CULL_FACE);
    this.gl.enable(this.gl.BLEND);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.sprite.vertexBuffer);
    this.gl.vertexAttribPointer(
      attributes.position,
      2,
      this.gl.FLOAT,
      false,
      2 * 8,
      0
    );
    this.gl.vertexAttribPointer(
      attributes.uv,
      2,
      this.gl.FLOAT,
      false,
      2 * 8,
      8
    );

    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.sprite.elementBuffer);

    this.gl.uniformMatrix4fv(
      uniforms.projectionMatrix,
      false,
      camera.projectionMatrix.elements
    );

    this.gl.activeTexture(this.gl.TEXTURE0);
    this.gl.uniform1i(uniforms.map, 0);

    var fog = scene.fog as Fog;

    if (fog) {
      this.gl.uniform3f(
        uniforms.fogColor,
        fog.color.r,
        fog.color.g,
        fog.color.b
      );

      this.gl.uniform1f(uniforms.fogNear, fog.near);
      this.gl.uniform1f(uniforms.fogFar, fog.far);
    } else {
      this.gl.uniform1f(uniforms.fogNear, 0);
      this.gl.uniform1f(uniforms.fogFar, 0);
    }

    // update positions and sort

    var i;
    let sprite: Sprite;
    let material;
    let size;

    for (i = 0; i < nSprites; i++) {
      sprite = sprites[i] as Sprite;
      material = sprite.material;
      if (!material) continue;
      if (material.depthTest == false && !inFront) continue;

      if (!sprite.visible || material.opacity === 0) continue;

      if (!material.useScreenCoordinates) {
        sprite._modelViewMatrix.multiplyMatrices(
          camera.matrixWorldInverse,
          sprite.matrixWorld
        );
        sprite.z = -sprite._modelViewMatrix.elements[14];
      } else {
        sprite.z = -sprite.position.z;
      }
    }

    sprites.sort(painterSortStable);

    // render all sprites
    for (i = 0; i < nSprites; i++) {
      sprite = sprites[i] as Sprite;
      material = sprite.material;
      if (!material) continue;
      if (!sprite.visible || material.opacity === 0) continue;

      if (material.map && material.map.image && material.map.image.width) {
        this.gl.uniform1f(uniforms?.alphaTest || null, material.alphaTest);
        var w = material.map.logicalWidth ?? material.map.image.width;
        var h = material.map.logicalHeight ?? material.map.image.height;

        let scaleX = (w * this.renderer.devicePixelRatio) / viewportWidth;
        let scaleY = (h * this.renderer.devicePixelRatio) / viewportHeight;

        if (material.useScreenCoordinates === true) {
          this.gl.uniform1i(uniforms.useScreenCoordinates, 1);
          this.gl.uniform3f(
            uniforms.screenPosition,
            (sprite.position.x * this.renderer.devicePixelRatio -
              halfViewportWidth) /
            halfViewportWidth,
            (halfViewportHeight -
              sprite.position.y * this.renderer.devicePixelRatio) /
            halfViewportHeight,
            Math.max(0, Math.min(1, sprite.position.z))
          );
        } else {
          this.gl.uniform1i(uniforms.useScreenCoordinates, 0);
          this.gl.uniformMatrix4fv(
            uniforms.modelViewMatrix,
            false,
            sprite._modelViewMatrix.elements
          );
        }

        size = 1 / (material.scaleByViewport ? viewportHeight : 1);

        scaleX *= size * sprite.scale.x;
        scaleY *= size * sprite.scale.y;

        let alignx = material?.alignment?.x,
          aligny = material?.alignment?.y;
        if (material.screenOffset) {
          //adjust alignment offset by screenOffset adjusted to sprite coords
          alignx = (alignx || 0) + (2.0 * material.screenOffset.x) / w;
          aligny = (aligny || 0) + (2.0 * material.screenOffset.y) / h;
        }

        this.gl.uniform2f(
          uniforms.uvScale,
          material?.uvScale?.x || 1,
          material?.uvScale?.y || 1
        );
        this.gl.uniform2f(
          uniforms.uvOffset,
          material?.uvOffset?.x || 0,
          material?.uvOffset?.y || 0
        );
        this.gl.uniform2f(uniforms.alignment, alignx || 0, aligny || 0);

        this.gl.uniform1f(uniforms.opacity, material.opacity);
        this.gl.uniform3f(
          uniforms.color,
          material?.color?.r || 0,
          material?.color?.g || 0,
          material?.color?.b || 0
        );

        this.gl.uniform1f(uniforms.rotation, sprite.rotation as number);
        this.gl.uniform2f(uniforms.scale, scaleX, scaleY);

        //this.renderer.setBlending( material.blending, material.blendEquation, material.blendSrc, material.blendDst );
        this.renderer.setDepthTest(material.depthTest);
        this.renderer.setDepthWrite(material.depthWrite);
        this.renderer.setTexture(material.map, 0);

        this.gl.drawElements(this.gl.TRIANGLES, 6, this.gl.UNSIGNED_SHORT, 0);
      }
    }

    // restore gl
    this.gl.enable(this.gl.CULL_FACE);
  }

  /** Release sprite resources while leaving the renderer's shared context usable. */
  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    if (this.gl) {
      if (this.sprite.vertexBuffer) this.gl.deleteBuffer(this.sprite.vertexBuffer);
      if (this.sprite.elementBuffer) this.gl.deleteBuffer(this.sprite.elementBuffer);
      if (this.sprite.program) this.gl.deleteProgram(this.sprite.program);
    }
    this.sprite.vertexBuffer = this.sprite.elementBuffer = this.sprite.program = null;
    this.sprite.vertices = this.sprite.faces = this.sprite.uniforms = null;
    this.sprite.attributes = {};
    this.sprites = undefined;
    this.renderList.length = 0;
    this.renderer = undefined;
    this.gl = undefined;
  }

  private createProgram(shader: Shader, precision: number): WebGLProgram {
    if (!this.gl) throw new Error("WebGL Rendering context not found");
    const program = this.gl.createProgram();
    if (!program) throw new Error("Error creating webgl program");
    let fragmentShader: WebGLShader = null;
    let vertexShader: WebGLShader = null;
    let attached = false;
    try {
      fragmentShader = this.gl.createShader(this.gl.FRAGMENT_SHADER);
      vertexShader = this.gl.createShader(this.gl.VERTEX_SHADER);
      if (!fragmentShader)
        throw new Error("Unable to create fragment shader SpritePlugin.createProgram");
      if (!vertexShader)
        throw new Error("Unable to create vertex shader SpritePlugin.createProgram");

      const prefix = "precision " + precision + " float;\n";
      this.gl.shaderSource(fragmentShader, prefix + shader.fragmentShader);
      this.gl.shaderSource(vertexShader, prefix + shader.vertexShader);
      this.gl.compileShader(fragmentShader);
      this.gl.compileShader(vertexShader);
      if (
        !this.gl.getShaderParameter(fragmentShader, this.gl.COMPILE_STATUS) ||
        !this.gl.getShaderParameter(vertexShader, this.gl.COMPILE_STATUS)
      ) {
        throw new Error(`Error compiling shader:
      ${this.gl.getShaderInfoLog(fragmentShader)}
      ${this.gl.getShaderInfoLog(vertexShader)}`);
      }
      this.gl.attachShader(program, fragmentShader);
      this.gl.attachShader(program, vertexShader);
      attached = true;
      this.gl.linkProgram(program);
      if (!this.gl.getProgramParameter(program, this.gl.LINK_STATUS))
        console.error("Could not initialize shader");
      return program;
    } catch (error) {
      this.gl.deleteProgram(program);
      attached = false;
      throw error;
    } finally {
      // Shader objects can be deleted after linking; the program retains its executable.
      if (attached) {
        this.gl.detachShader(program, fragmentShader);
        this.gl.detachShader(program, vertexShader);
      }
      if (fragmentShader) this.gl.deleteShader(fragmentShader);
      if (vertexShader) this.gl.deleteShader(vertexShader);
    }
  }

}
function painterSortStable(a: any, b: any): number {
  if (a.z !== b.z) {
    return b.z - a.z;
  } else {
    return b.id - a.id;
  }
}
