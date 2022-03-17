// @ts-check
//Render plugins go here

import { Scene } from "./core";
import { Material, SpriteMaterial } from "./materials";
import { Camera, Sprite } from "./objects";
import { Renderer } from "./Renderer";
import { ShaderLib } from "./ShaderUtils";

/**
 * @param {{ z: number; id: number; }} a
 * @param {{ z: number; id: number; }} b
 */
function painterSortStable(a, b) {
  if (a.z !== b.z) {
    return b.z - a.z;
  } else {
    return b.id - a.id;
  }
}

/**
 * Sprite render plugin
 * @this {SpritePlugin}
 */

export class SpritePlugin {
  /**
   * @type {WebGLRenderingContext}
   */
  _gl;
  /**
   * @type {Renderer}
   */
  _renderer;
  /**
   * @type {string}
   */
  _precision="";

  _sprite = {};

  /**
   * @param {Renderer} renderer
   */
   constructor(renderer) {
    if (!renderer._gl) throw new Error("gl context not found");
    this._gl = renderer._gl;
    this._renderer = renderer;

    this._precision = renderer.getPrecision();

    this._sprite.vertices = new Float32Array(8 + 8);
    this._sprite.faces = new Uint16Array(6);

    var i = 0;

    this._sprite.vertices[i++] = -1;
    this._sprite.vertices[i++] = -1; // vertex 0
    this._sprite.vertices[i++] = 0;
    this._sprite.vertices[i++] = 0; // uv 0

    this._sprite.vertices[i++] = 1;
    this._sprite.vertices[i++] = -1; // vertex 1
    this._sprite.vertices[i++] = 1;
    this._sprite.vertices[i++] = 0; // uv 1

    this._sprite.vertices[i++] = 1;
    this._sprite.vertices[i++] = 1; // vertex 2
    this._sprite.vertices[i++] = 1;
    this._sprite.vertices[i++] = 1; // uv 2

    this._sprite.vertices[i++] = -1;
    this._sprite.vertices[i++] = 1; // vertex 3
    this._sprite.vertices[i++] = 0;
    this._sprite.vertices[i++] = 1; // uv 3

    i = 0;

    this._sprite.faces[i++] = 0;
    this._sprite.faces[i++] = 1;
    this._sprite.faces[i++] = 2;
    this._sprite.faces[i++] = 0;
    this._sprite.faces[i++] = 2;
    this._sprite.faces[i++] = 3;

    this._sprite.vertexBuffer = this._gl.createBuffer();
    this._sprite.elementBuffer = this._gl.createBuffer();

    this._gl.bindBuffer(this._gl.ARRAY_BUFFER, this._sprite.vertexBuffer);
    this._gl.bufferData(
      this._gl.ARRAY_BUFFER,
      this._sprite.vertices,
      this._gl.STATIC_DRAW
    );

    this._gl.bindBuffer(
      this._gl.ELEMENT_ARRAY_BUFFER,
      this._sprite.elementBuffer
    );
    this._gl.bufferData(
      this._gl.ELEMENT_ARRAY_BUFFER,
      this._sprite.faces,
      this._gl.STATIC_DRAW
    );

    this._sprite.program = this.createProgram(
      ShaderLib.sprite,
      this._precision
    );
    if (!this._sprite.program) throw new Error("shader program created");

    this._sprite.attributes = {};
    this._sprite.uniforms = {};

    this._sprite.attributes.position = this._gl.getAttribLocation(
      this._sprite.program,
      "position"
    );
    this._sprite.attributes.uv = this._gl.getAttribLocation(
      this._sprite.program,
      "uv"
    );

    this._sprite.uniforms.uvOffset = this._gl.getUniformLocation(
      this._sprite.program,
      "uvOffset"
    );
    this._sprite.uniforms.uvScale = this._gl.getUniformLocation(
      this._sprite.program,
      "uvScale"
    );

    this._sprite.uniforms.rotation = this._gl.getUniformLocation(
      this._sprite.program,
      "rotation"
    );
    this._sprite.uniforms.scale = this._gl.getUniformLocation(
      this._sprite.program,
      "scale"
    );
    this._sprite.uniforms.alignment = this._gl.getUniformLocation(
      this._sprite.program,
      "alignment"
    );

    this._sprite.uniforms.color = this._gl.getUniformLocation(
      this._sprite.program,
      "color"
    );
    this._sprite.uniforms.map = this._gl.getUniformLocation(
      this._sprite.program,
      "map"
    );
    this._sprite.uniforms.opacity = this._gl.getUniformLocation(
      this._sprite.program,
      "opacity"
    );

    this._sprite.uniforms.useScreenCoordinates = this._gl.getUniformLocation(
      this._sprite.program,
      "useScreenCoordinates"
    );
    this._sprite.uniforms.screenPosition = this._gl.getUniformLocation(
      this._sprite.program,
      "screenPosition"
    );
    this._sprite.uniforms.modelViewMatrix = this._gl.getUniformLocation(
      this._sprite.program,
      "modelViewMatrix"
    );
    this._sprite.uniforms.projectionMatrix = this._gl.getUniformLocation(
      this._sprite.program,
      "projectionMatrix"
    );

    this._sprite.uniforms.fogType = this._gl.getUniformLocation(
      this._sprite.program,
      "fogType"
    );
    this._sprite.uniforms.fogDensity = this._gl.getUniformLocation(
      this._sprite.program,
      "fogDensity"
    );
    this._sprite.uniforms.fogNear = this._gl.getUniformLocation(
      this._sprite.program,
      "fogNear"
    );
    this._sprite.uniforms.fogFar = this._gl.getUniformLocation(
      this._sprite.program,
      "fogFar"
    );
    this._sprite.uniforms.fogColor = this._gl.getUniformLocation(
      this._sprite.program,
      "fogColor"
    );

    this._sprite.uniforms.alphaTest = this._gl.getUniformLocation(
      this._sprite.program,
      "alphaTest"
    );
  }

  /**
   * @param {Scene} scene
   * @param {Camera} camera
   * @param {number} viewportWidth
   * @param {number} viewportHeight
   * @param {boolean} inFront
   */
  render(scene, camera, viewportWidth, viewportHeight, inFront) {
    /**
     * @type {any[]}
     */
    let sprites = [];
    scene.__webglSprites.forEach((/** @type {{ material: { depthTest: boolean; }; }} */ sprite) => {
      //depthTest is false for inFront labels
      if (inFront && sprite.material.depthTest == false) {
        sprites.push(sprite);
      } else if (!inFront && sprite.material.depthTest) {
        sprites.push(sprite);
      }
    });

    let nSprites = sprites.length;

    if (!nSprites) return;

    var attributes = this._sprite.attributes,
      uniforms = this._sprite.uniforms;

    var halfViewportWidth = viewportWidth * 0.5,
      halfViewportHeight = viewportHeight * 0.5;

    // setup gl

    this._gl.useProgram(this._sprite.program);

    this._gl.enableVertexAttribArray(attributes.position);
    this._gl.enableVertexAttribArray(attributes.uv);

    this._gl.disable(this._gl.CULL_FACE);
    this._gl.enable(this._gl.BLEND);

    this._gl.bindBuffer(this._gl.ARRAY_BUFFER, this._sprite.vertexBuffer);
    this._gl.vertexAttribPointer(
      attributes.position,
      2,
      this._gl.FLOAT,
      false,
      2 * 8,
      0
    );
    this._gl.vertexAttribPointer(
      attributes.uv,
      2,
      this._gl.FLOAT,
      false,
      2 * 8,
      8
    );

    this._gl.bindBuffer(
      this._gl.ELEMENT_ARRAY_BUFFER,
      this._sprite.elementBuffer
    );

    this._gl.uniformMatrix4fv(
      uniforms.projectionMatrix,
      false,
      camera.projectionMatrix.elements
    );

    this._gl.activeTexture(this._gl.TEXTURE0);
    this._gl.uniform1i(uniforms.map, 0);

    var oldFogType = 0;
    var sceneFogType = 0;
    var fog = scene.fog;

    if (fog) {
      this._gl.uniform3f(
        uniforms.fogColor,
        fog.color.r,
        fog.color.g,
        fog.color.b
      );

      this._gl.uniform1f(uniforms.fogNear, fog.near);
      this._gl.uniform1f(uniforms.fogFar, fog.far);

      this._gl.uniform1i(uniforms.fogType, 1);
      oldFogType = 1;
      sceneFogType = 1;
    } else {
      this._gl.uniform1i(uniforms.fogType, 0);
      oldFogType = 0;
      sceneFogType = 0;
    }

    // update positions and sort

    var i,
      sprite,
      material,
      size,
      fogType,
      scale = [];

    for (i = 0; i < nSprites; i++) {
      sprite = sprites[i];
      material = sprite.material;
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
      sprite = sprites[i];
      material = sprite.material;

      if (!sprite.visible || material.opacity === 0) continue;

      if (material.map && material.map.image && material.map.image.width) {
        this._gl.uniform1f(uniforms.alphaTest, material.alphaTest);
        var w = material.map.image.width;
        var h = material.map.image.height;

        scale[0] = (w * this._renderer.devicePixelRatio) / viewportWidth;
        scale[1] = (h * this._renderer.devicePixelRatio) / viewportHeight;

        if (material.useScreenCoordinates === true) {
          this._gl.uniform1i(uniforms.useScreenCoordinates, 1);
          this._gl.uniform3f(
            uniforms.screenPosition,
            (sprite.position.x * this._renderer.devicePixelRatio -
              halfViewportWidth) /
              halfViewportWidth,
            (halfViewportHeight -
              sprite.position.y * this._renderer.devicePixelRatio) /
              halfViewportHeight,
            Math.max(0, Math.min(1, sprite.position.z))
          );
        } else {
          this._gl.uniform1i(uniforms.useScreenCoordinates, 0);
          this._gl.uniformMatrix4fv(
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
          this._gl.uniform1i(uniforms.fogType, fogType);
          oldFogType = fogType;
        }

        size = 1 / (material.scaleByViewport ? viewportHeight : 1);

        scale[0] *= size * sprite.scale.x;
        scale[1] *= size * sprite.scale.y;

        let alignx = material.alignment.x,
          aligny = material.alignment.y;
        if (material.screenOffset) {
          //adjust alignment offset by screenOffset adjusted to sprite coords
          alignx += (2.0 * material.screenOffset.x) / w;
          aligny += (2.0 * material.screenOffset.y) / h;
        }

        this._gl.uniform2f(
          uniforms.uvScale,
          material.uvScale.x,
          material.uvScale.y
        );
        this._gl.uniform2f(
          uniforms.uvOffset,
          material.uvOffset.x,
          material.uvOffset.y
        );
        this._gl.uniform2f(uniforms.alignment, alignx, aligny);

        this._gl.uniform1f(uniforms.opacity, material.opacity);
        this._gl.uniform3f(
          uniforms.color,
          material.color.r,
          material.color.g,
          material.color.b
        );

        this._gl.uniform1f(uniforms.rotation, sprite.rotation);
        this._gl.uniform2fv(uniforms.scale, scale);

        //this._renderer.setBlending( material.blending, material.blendEquation, material.blendSrc, material.blendDst );
        this._renderer.setDepthTest(material.depthTest);
        this._renderer.setDepthWrite(material.depthWrite);
        this._renderer.setTexture(material.map, 0);

        this._gl.drawElements(
          this._gl.TRIANGLES,
          6,
          this._gl.UNSIGNED_SHORT,
          0
        );
      }
    }

    // restore gl
    this._gl.enable(this._gl.CULL_FACE);
  }

  /**
   * @param {{ fragmentShader: any; vertexShader: any; uniforms?: {}; }} shader
   * @param {string} precision
   */
  createProgram(shader, precision) {
    var program = this._gl.createProgram();

    var fragmentShader = this._gl.createShader(this._gl.FRAGMENT_SHADER);
    var vertexShader = this._gl.createShader(this._gl.VERTEX_SHADER);

    var prefix = "precision " + precision + " float;\n";

    if (!vertexShader || !fragmentShader || !program) throw new Error("failed to initialize shader program")

    this._gl.shaderSource(fragmentShader, prefix + shader.fragmentShader);
    this._gl.shaderSource(vertexShader, prefix + shader.vertexShader);

    this._gl.compileShader(fragmentShader);
    this._gl.compileShader(vertexShader);

    if (
      !this._gl.getShaderParameter(fragmentShader, this._gl.COMPILE_STATUS) ||
      !this._gl.getShaderParameter(vertexShader, this._gl.COMPILE_STATUS)
    ) {
      console.error(this._gl.getShaderInfoLog(fragmentShader));
      console.error("could not initialize shader");
      return null;
    }

    this._gl.attachShader(program, fragmentShader);
    this._gl.attachShader(program, vertexShader);

    this._gl.linkProgram(program);

    if (!this._gl.getProgramParameter(program, this._gl.LINK_STATUS))
      console.error("Could not initialize shader");

    return program;
  }
}
