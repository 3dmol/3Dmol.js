/* eslint-disable no-param-reassign */
/* eslint-disable max-classes-per-file */
import Color from './core/Color';
import EventDispatcher from './core/EventDispatcher';
import {Vector2, Vector3} from './math';

// sides
export const FrontSide = 0;
export const BackSide = 1;
export const DoubleSide = 2;

// shading
export const NoShading = 0;
export const FlatShading = 1;
export const SmoothShading = 2;

// colors
export const NoColors = 0;
export const FaceColors = 1;
export const VertexColors = 2;

// Texture constants
// TODO: Which of these do I need (since I only use textures to display label sprites) ?
export const MultiplyOperation = 0;
export const MixOperation = 1;
export const AddOperation = 2;

// mapping modes

export function UVMapping() {}

// wrapping modes
export const ClampToEdgeWrapping = 1001;

// Filters
export const LinearFilter = 1006;
export const NearestFilter = 1007;
export const LinearMipMapLinearFilter = 1008;

// Data types
export const UnsignedByteType = 1009;
export const FloatType = 1010;

// Pixel formats
export const RGBAFormat = 1021;
export const RFormat = 1022;
export const R32Format = 1023;

// eslint-disable-next-line import/no-mutable-exports
export let MaterialIdCount = 0;

/**
 * Line and Mesh material types
 * @constructor
 */
/**
 * Line and Mesh material types
 * @constructor
 */
export class Material extends EventDispatcher {
  overdraw;
  voldata;
  volscheme;

  constructor() {
    super();

    this.id = MaterialIdCount += 1;

    this.name = '';

    // TODO: Which of these instance variables can I remove??
    this.side = FrontSide;

    this.opacity = 1;
    this.transparent = false;

    this.depthTest = true;
    this.depthWrite = true;

    this.stencilTest = true;

    this.polygonOffset = false;
    this.polygonOffsetFactor = 0;
    this.polygonOffsetUnits = 0;

    this.alphaTest = 0;

    this.visible = true;

    this.needsUpdate = true;
  }

  setValues(values) {
    if (values === undefined) return;

    for (const key in values) {
      const newValue = values[key];

      if (newValue === undefined) {
        console.warn(`Material: '${key}' parameter is undefined.`);
        continue;
      }

      if (key in this) {
        const currentValue = this[key];

        if (currentValue instanceof Color && newValue instanceof Color) {
          currentValue.copy(newValue);
        } else if (currentValue instanceof Color) {
          currentValue.set(newValue);
        } else if (currentValue instanceof Vector3 && newValue instanceof Vector3) {
          currentValue.copy(newValue);
        } else {
          this[key] = newValue;
        }
      }
    }
  }

  // TODO: might want to look into blending equations
  clone(material) {
    if (material === undefined) material = new Material();

    material.name = this.name;

    material.side = this.side;

    material.opacity = this.opacity;
    material.transparent = this.transparent;

    material.depthTest = this.depthTest;
    material.depthWrite = this.depthWrite;
    material.stencilTest = this.stencilTest;

    material.polygonOffset = this.polygonOffset;
    material.polygonOffsetFactor = this.polygonOffsetFactor;
    material.polygonOffsetUnits = this.polygonOffsetUnits;

    material.alphaTest = this.alphaTest;

    material.overdraw = this.overdraw;

    material.visible = this.visible;

    return material;
  }

  dispose() {
    this.dispatchEvent({type: 'dispose'});
  }
}

// Line basic material
/** @constructor */
// Line basic material
/** @constructor */
export class LineBasicMaterial extends Material {
  constructor(parameters) {
    super();

    this.color = new Color(0xffffff);

    this.linewidth = 1;
    this.linecap = 'round';
    this.linejoin = 'round';

    this.vertexColors = false;

    this.fog = true;
    this.shaderID = 'basic';
    this.setValues(parameters);
  }

  clone() {
    const material = new LineBasicMaterial();

    Material.prototype.clone.call(this, material);

    material.color.copy(this.color);
    return material;
  }
}

// Mesh Lambert material
/** @constructor */
// Mesh Lambert material
/** @constructor */
export class MeshLambertMaterial extends Material {
  combine;

  morphTargets;

  morphNormals;

  constructor(parameters) {
    super();

    this.color = new Color(0xffffff);
    this.ambient = new Color(0xfffff);
    this.emissive = new Color(0x000000);

    // TODO: Which of these instance variables do I really need?
    this.wrapAround = false;
    this.wrapRGB = new Vector3(1, 1, 1);

    this.map = null;

    this.lightMap = null;

    this.specularMap = null;

    this.envMap = null;
    this.reflectivity = 1;
    this.refractionRatio = 0.98;

    this.fog = true;

    this.wireframe = false;
    this.wireframeLinewidth = 1;
    this.wireframeLinecap = 'round';
    this.wireframeLinejoin = 'round';

    this.shading = SmoothShading;
    this.shaderID = 'lambert';
    this.vertexColors = NoColors;

    this.skinning = false;

    this.setValues(parameters);
  }

  clone(material) {
    if (typeof material === 'undefined') material = new MeshLambertMaterial();

    Material.prototype.clone.call(this, material);

    material.color.copy(this.color);
    material.ambient.copy(this.ambient);
    material.emissive.copy(this.emissive);

    material.wrapAround = this.wrapAround;
    material.wrapRGB.copy(this.wrapRGB);

    material.map = this.map;

    material.lightMap = this.lightMap;

    material.specularMap = this.specularMap;

    material.envMap = this.envMap;
    material.combine = this.combine;
    material.reflectivity = this.reflectivity;
    material.refractionRatio = this.refractionRatio;

    material.fog = this.fog;

    material.shading = this.shading;
    material.shaderID = this.shaderID;
    material.vertexColors = this.vertexColors;

    material.skinning = this.skinning;
    material.morphTargets = this.morphTargets;
    material.morphNormals = this.morphNormals;

    return material;
  }
}

// Double sided Mesh Lambert material
/** @constructor */
// Double sided Mesh Lambert material
/** @constructor */
export class MeshDoubleLambertMaterial extends MeshLambertMaterial {
  outline;

  constructor(parameters) {
    super(parameters);

    this.shaderID = 'lambertdouble';
    this.side = DoubleSide;
  }

  clone() {
    const material = new MeshDoubleLambertMaterial();

    MeshLambertMaterial.prototype.clone.call(this, material);

    return material;
  }
}

// Outlined Mesh Lamert material
export class MeshOutlineMaterial extends Material {
  constructor(parameters) {
    super();
    parameters = parameters || {};
    this.fog = true;
    this.shaderID = 'outline';
    this.wireframe = false;
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;
  }

  clone(material) {
    if (typeof material === 'undefined') material = new MeshOutlineMaterial();
    Material.prototype.clone.call(this, material);
    material.fog = this.fog;
    material.shaderID = this.shaderID;
    material.wireframe = this.wireframe;
    return material;
  }
}

// Imposter material
export class ImposterMaterial extends Material {
  combine;

  morphTargets;

  morphNormals;

  constructor(parameters) {
    super();

    this.color = new Color(0xffffff);
    this.ambient = new Color(0xfffff);
    this.emissive = new Color(0x000000);
    this.imposter = true;

    // TODO: Which of these instance variables do I really need?
    this.wrapAround = false;
    this.wrapRGB = new Vector3(1, 1, 1);

    this.map = null;

    this.lightMap = null;

    this.specularMap = null;

    this.envMap = null;
    this.reflectivity = 1;
    this.refractionRatio = 0.98;

    this.fog = true;

    this.wireframe = false;
    this.wireframeLinewidth = 1;
    this.wireframeLinecap = 'round';
    this.wireframeLinejoin = 'round';

    this.shading = SmoothShading;
    this.shaderID = null;
    this.vertexColors = NoColors;

    this.skinning = false;

    this.setValues(parameters);
  }

  clone(mat) {
    const material = mat || new ImposterMaterial();

    Material.prototype.clone.call(this, material);

    material.color.copy(this.color);
    material.ambient.copy(this.ambient);
    material.emissive.copy(this.emissive);

    material.wrapAround = this.wrapAround;
    material.wrapRGB.copy(this.wrapRGB);

    material.map = this.map;

    material.lightMap = this.lightMap;

    material.specularMap = this.specularMap;

    material.envMap = this.envMap;
    material.combine = this.combine;
    material.reflectivity = this.reflectivity;
    material.refractionRatio = this.refractionRatio;

    material.fog = this.fog;

    material.shading = this.shading;
    material.shaderID = this.shaderID;
    material.vertexColors = this.vertexColors;

    material.skinning = this.skinning;
    material.morphTargets = this.morphTargets;
    material.morphNormals = this.morphNormals;

    return material;
  }
}

export class SphereImposterMaterial extends ImposterMaterial {
  constructor(parameters) {
    super();

    this.shaderID = 'sphereimposter';
    this.setValues(parameters);
  }

  clone(mat) {
    const material = mat || new SphereImposterMaterial();
    ImposterMaterial.prototype.clone.call(this, material);
    return material;
  }
}

export class SphereImposterOutlineMaterial extends ImposterMaterial {
  constructor(parameters) {
    super();
    parameters = parameters || {};

    this.shaderID = 'sphereimposteroutline';
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;

    this.setValues(parameters);
  }

  clone() {
    const material = new SphereImposterOutlineMaterial();
    ImposterMaterial.prototype.clone.call(this, material);
    material.outlineColor = this.outlineColor;
    material.outlineWidth = this.outlineWidth;
    material.outlinePushback = this.outlinePushback;
    return material;
  }
}

export class StickImposterMaterial extends ImposterMaterial {
  constructor(parameters) {
    super();

    this.shaderID = 'stickimposter';
    this.setValues(parameters);
  }

  clone() {
    const material = new StickImposterMaterial();
    ImposterMaterial.prototype.clone.call(this, material);
    return material;
  }
}

export class StickImposterOutlineMaterial extends ImposterMaterial {
  constructor(parameters) {
    super();
    parameters = parameters || {};

    this.shaderID = 'stickimposteroutline';
    this.outlineColor = parameters.color || new Color(0.0, 0.0, 0.0);
    this.outlineWidth = parameters.width || 0.1;
    this.outlinePushback = parameters.pushback || 1.0;

    this.setValues(parameters);
  }

  clone() {
    const material = new StickImposterOutlineMaterial();
    ImposterMaterial.prototype.clone.call(this, material);
    material.outlineColor = this.outlineColor;
    material.outlineWidth = this.outlineWidth;
    material.outlinePushback = this.outlinePushback;
    return material;
  }
}

export class InstancedMaterial extends Material {
  combine;

  morphTargets;

  morphNormals;

  constructor(parameters) {
    super();

    this.color = new Color(0xffffff);
    this.ambient = new Color(0xfffff);
    this.emissive = new Color(0x000000);

    // TODO: Which of these instance variables do I really need?
    this.wrapAround = false;
    this.wrapRGB = new Vector3(1, 1, 1);

    this.map = null;

    this.lightMap = null;

    this.specularMap = null;

    this.envMap = null;
    this.reflectivity = 1;
    this.refractionRatio = 0.98;

    this.fog = true;

    this.wireframe = false;
    this.wireframeLinewidth = 1;
    this.wireframeLinecap = 'round';
    this.wireframeLinejoin = 'round';

    this.shading = SmoothShading;
    this.shaderID = 'instanced';
    this.vertexColors = NoColors;

    this.skinning = false;

    this.sphere = null;

    this.setValues(parameters);
  }

  clone() {
    const material = new InstancedMaterial();

    Material.prototype.clone.call(this, material);

    material.color.copy(this.color);
    material.ambient.copy(this.ambient);
    material.emissive.copy(this.emissive);

    material.wrapAround = this.wrapAround;
    material.wrapRGB.copy(this.wrapRGB);

    material.map = this.map;

    material.lightMap = this.lightMap;

    material.specularMap = this.specularMap;

    material.envMap = this.envMap;
    material.combine = this.combine;
    material.reflectivity = this.reflectivity;
    material.refractionRatio = this.refractionRatio;

    material.fog = this.fog;

    material.shading = this.shading;
    material.shaderID = this.shaderID;
    material.vertexColors = this.vertexColors;

    material.skinning = this.skinning;
    material.morphTargets = this.morphTargets;
    material.morphNormals = this.morphNormals;

    material.sphere = this.sphere;

    return material;
  }
}

// Volumetric material
export class VolumetricMaterial extends Material {
  constructor(parameters) {
    super();

    this.transparent = false;
    this.volumetric = true;

    this.color = new Color(0xffffff);
    this.transferfn = null;
    this.map = null;
    this.volumetric = true;
    this.extent = [];
    this.maxdepth = 100.0;
    this.unit = 0;
    this.texmatrix = null;
    this.transfermin = -1.0;
    this.transfermax = 1.0;
    this.subsamples = 5.0;

    // this.fog = true; // TODO: to integrate the new shader with the fog stuff
    this.shaderID = 'volumetric';
    this.side = FrontSide;

    this.setValues(parameters);
  }

  clone() {
    const material = Object.assign(new VolumetricMaterial(), this);

    Material.prototype.clone.call(this, material);
    return material;
  }
}

// Texture
// eslint-disable-next-line import/no-mutable-exports
export let TextureIdCount = 0;
// We really only create textures from 2d rendering contexts (to display text labels)
// edit: we can now create 3dtextures using volumetric data
export class Texture extends EventDispatcher {
  constructor(image, is3D) {
    super();

    this.id = TextureIdCount += 1;

    this.name = '';

    this.image = image;

    this.mapping = new UVMapping();

    this.wrapS = ClampToEdgeWrapping;
    this.wrapT = ClampToEdgeWrapping;

    this.anisotropy = 1;

    if (is3D) {
      this.format = RFormat;
      this.type = FloatType;

      this.premultiplyAlpha = false;
      this.flipY = false;

      this.unpackAlignment = 1;

      this.magFilter = NearestFilter;
      this.minFilter = NearestFilter;
    } else {
      this.format = RGBAFormat;
      this.type = UnsignedByteType;

      this.offset = new Vector2(0, 0);
      this.repeat = new Vector2(1, 1);

      this.premultiplyAlpha = false;
      this.flipY = true;
      this.unpackAlignment = 4;

      this.magFilter = LinearFilter;
      this.minFilter = LinearMipMapLinearFilter;
    }

    this.needsUpdate = false;
    this.onUpdate = null;
  }

  clone(tex) {
    const texture = tex || new Texture();

    texture.image = this.image;

    texture.mapping = this.mapping;

    texture.wrapS = this.wrapS;
    texture.wrapT = this.wrapT;

    texture.magFilter = this.magFilter;
    texture.minFilter = this.minFilter;

    texture.anisotropy = this.anisotropy;

    texture.format = this.format;
    texture.type = this.type;

    texture.offset.copy(this.offset);
    texture.repeat.copy(this.repeat);

    texture.premultiplyAlpha = this.premultiplyAlpha;
    texture.flipY = this.flipY;
    texture.unpackAlignment = this.unpackAlignment;

    return texture;
  }

  dispose() {
    this.dispatchEvent({type: 'dispose'});
  }
}

// Alignment for Sprites

export const SpriteAlignment = {};
SpriteAlignment.topLeft = new Vector2(1, -1);
SpriteAlignment.topCenter = new Vector2(0, -1);
SpriteAlignment.topRight = new Vector2(-1, -1);
SpriteAlignment.centerLeft = new Vector2(1, 0);
SpriteAlignment.center = new Vector2(0, 0);
SpriteAlignment.centerRight = new Vector2(-1, 0);
SpriteAlignment.bottomLeft = new Vector2(1, 1);
SpriteAlignment.bottomCenter = new Vector2(0, 1);
SpriteAlignment.bottomRight = new Vector2(-1, 1);

// Sprite material
export class SpriteMaterial extends Material {
  scaleByViewport;

  constructor(parameters) {
    super();

    this.color = new Color(0xffffff);
    this.map = new Texture();

    this.useScreenCoordinates = true;
    this.depthTest = !this.useScreenCoordinates;
    this.sizeAttenuation = !this.useScreenCoordinates;
    this.screenOffset = this.screenOffset;
    this.scaleByViewPort = !this.sizeAttenuation;
    this.alignment = SpriteAlignment.center.clone();

    this.fog = false; // use scene fog

    this.uvOffset = new Vector2(0, 0);
    this.uvScale = new Vector2(1, 1);

    this.setValues(parameters);

    parameters = parameters || {};

    if (parameters.depthTest === undefined) this.depthTest = !this.useScreenCoordinates;
    if (parameters.sizeAttenuation === undefined) this.sizeAttenuation = !this.useScreenCoordinates;
    if (parameters.scaleByViewPort === undefined) this.scaleByViewPort = !this.sizeAttenuation;
  }

  clone() {
    const material = new SpriteMaterial();

    Material.prototype.clone.call(this, material);

    material.color.copy(this.color);
    material.map = this.map;

    material.useScreenCoordinates = this.useScreenCoordinates;
    material.screenOffset = this.screenOffset;
    material.sizeAttenuation = this.sizeAttenuation;
    material.scaleByViewport = this.scaleByViewPort;
    material.alignment.copy(this.alignment);

    material.uvOffset.copy(this.uvOffset);

    return material;
  }
}


