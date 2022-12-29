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
  fogType: NullableWebGLUniformLocation;
  fogDensity: NullableWebGLUniformLocation;
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

    this.sprite.program = this.createProgram(
      ShaderLib.sprite,
      this.precision || 1 /** added default to single precision */
    );

    this.sprite.attributes = {};
    const uniforms = {} as Partial<UniformsValueType>;

    this.sprite.attributes.position = this.gl.getAttribLocation(
      this.sprite.program,
      "position"
    );
    this.sprite.attributes.uv = this.gl.getAttribLocation(
      this.sprite.program,
      "uv"
    );

    uniforms.uvOffset = this.gl.getUniformLocation(
      this.sprite.program,
      "uvOffset"
    );
    uniforms.uvScale = this.gl.getUniformLocation(
      this.sprite.program,
      "uvScale"
    );
    uniforms.rotation = this.gl.getUniformLocation(
      this.sprite.program,
      "rotation"
    );
    uniforms.scale = this.gl.getUniformLocation(this.sprite.program, "scale");
    uniforms.alignment = this.gl.getUniformLocation(
      this.sprite.program,
      "alignment"
    );
    uniforms.color = this.gl.getUniformLocation(this.sprite.program, "color");
    uniforms.map = this.gl.getUniformLocation(this.sprite.program, "map");
    uniforms.opacity = this.gl.getUniformLocation(
      this.sprite.program,
      "opacity"
    );
    uniforms.useScreenCoordinates = this.gl.getUniformLocation(
      this.sprite.program,
      "useScreenCoordinates"
    );
    uniforms.screenPosition = this.gl.getUniformLocation(
      this.sprite.program,
      "screenPosition"
    );
    uniforms.modelViewMatrix = this.gl.getUniformLocation(
      this.sprite.program,
      "modelViewMatrix"
    );
    uniforms.projectionMatrix = this.gl.getUniformLocation(
      this.sprite.program,
      "projectionMatrix"
    );
    uniforms.fogType = this.gl.getUniformLocation(
      this.sprite.program,
      "fogType"
    );
    uniforms.fogDensity = this.gl.getUniformLocation(
      this.sprite.program,
      "fogDensity"
    );
    uniforms.fogNear = this.gl.getUniformLocation(
      this.sprite.program,
      "fogNear"
    );
    uniforms.fogFar = this.gl.getUniformLocation(this.sprite.program, "fogFar");
    uniforms.fogColor = this.gl.getUniformLocation(
      this.sprite.program,
      "fogColor"
    );
    uniforms.alphaTest = this.gl.getUniformLocation(
      this.sprite.program,
      "alphaTest"
    );

    this.sprite.uniforms = uniforms as UniformsValueType;
  }

  render(
    scene: Scene,
    camera: Camera,
    viewportWidth: number,
    viewportHeight: number,
    inFront?: boolean
  ) {
    if (!this.gl) throw new Error("WebGLRenderer not initialized");
    let sprites: unknown[] = [];
    scene?.__webglSprites?.forEach((sprite) => {
      //depthTest is false for inFront labels
      if (inFront && sprite.material.depthTest == false) {
        sprites.push(sprite);
      } else if (!inFront && sprite.material.depthTest) {
        sprites.push(sprite);
      }
    });

    let nSprites = sprites.length;

    if (!nSprites) return;

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

    var oldFogType = 0;
    var sceneFogType = 0;
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

      this.gl.uniform1i(uniforms.fogType, 1);
      oldFogType = 1;
      sceneFogType = 1;
    } else {
      this.gl.uniform1i(uniforms.fogType, 0);
      oldFogType = 0;
      sceneFogType = 0;
    }

    // update positions and sort

    var i;
    let sprite: Sprite;
    let material;
    let size;
    let fogType;
    let scale: number[] = [];

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
        var w = material.map.image.width;
        var h = material.map.image.height;

        scale[0] = (w * this.renderer.devicePixelRatio) / viewportWidth;
        scale[1] = (h * this.renderer.devicePixelRatio) / viewportHeight;

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

        if (scene.fog && material.fog) {
          fogType = sceneFogType;
        } else {
          fogType = 0;
        }

        if (oldFogType !== fogType) {
          this.gl.uniform1i(uniforms.fogType, fogType);
          oldFogType = fogType;
        }

        size = 1 / (material.scaleByViewport ? viewportHeight : 1);

        scale[0] *= size * sprite.scale.x;
        scale[1] *= size * sprite.scale.y;

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
        this.gl.uniform2fv(uniforms.scale, scale);

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

  private createProgram(shader: Shader, precision: number): WebGLProgram {
    if (!this.gl) throw new Error("WebGL Rendering context not found");
    var program = this.gl.createProgram();

    if (!program) throw new Error("Error creating webgl program");

    var fragmentShader = this.gl.createShader(this.gl.FRAGMENT_SHADER);
    var vertexShader = this.gl.createShader(this.gl.VERTEX_SHADER);

    if (!fragmentShader)
      throw new Error(
        "Unable to create fragment shader SpritePlugin.createProgram"
      );
    if (!vertexShader)
      throw new Error(
        "Unable to create vertex shader SpritePlugin.createProgram"
      );

    var prefix = "precision " + precision + " float;\n";

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

    this.gl.linkProgram(program);

    if (!this.gl.getProgramParameter(program, this.gl.LINK_STATUS))
      console.error("Could not initialize shader");

    return program;
  }
}
function painterSortStable(a: any, b: any): number {
  if (a.z !== b.z) {
    return b.z - a.z;
  } else {
    return b.id - a.id;
  }
}
