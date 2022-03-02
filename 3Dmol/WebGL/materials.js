// @ts-check

import { Color, EventDispatcher } from "./core";
import { Vector2, Vector3 } from "./math";

/**
 * Line and Mesh material types
 * @constructor
 */
export class Material {
    constructor() {

        EventDispatcher.call( this );

        this.id = MaterialIdCount ++;

        this.name = '';

        //TODO: Which of these instance variables can I remove??
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

        if ( values === undefined ) return;

        for ( var key in values ) {

            var newValue = values[ key ];

            if ( newValue === undefined ) {

                console.warn( 'Material: \'' + key + '\' parameter is undefined.' );
                continue;

            }

            if ( key in this ) {

                var currentValue = this[ key ];

                if ( currentValue instanceof Color && newValue instanceof Color ) {

                    currentValue.copy( newValue );

                } else if ( currentValue instanceof Color ) {

                    currentValue.set( newValue );

                } else if ( currentValue instanceof Vector3 && newValue instanceof Vector3 ) {

                    currentValue.copy( newValue );

                } else {

                    this[ key ] = newValue;

                }

            }

        }

    }

    clone(material) {

        if ( material === undefined ) material = new Material();

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

        this.dispatchEvent( { type: 'dispose' } );

    }
}

export var MaterialIdCount = 0;

//Line basic material
/** @constructor */
export class LineBasicMaterial {
    constructor(parameters) {

        Material.call(this);

        this.color = new Color(0xffffff);

        this.linewidth = 1;
        this.linecap = 'round';
        this.linejoin = 'round';

        this.vertexColors = false;

        this.fog = true;
        this.shaderID = "basic";
        this.setValues(parameters);

    }

    clone() {

        var material = new LineBasicMaterial();

        Material.prototype.clone.call(this, material);

        material.color.copy(this.color);
        return material;
    }
}

LineBasicMaterial.prototype = Object.create(Material.prototype);

//Mesh Lambert material
/** @constructor */
export class MeshLambertMaterial {
    constructor(parameters) {

        Material.call(this);

        this.color = new Color(0xffffff);
        this.ambient = new Color(0xfffff);
        this.emissive = new Color(0x000000);

        //TODO: Which of these instance variables do I really need?
        this.wrapAround = false;
        this.wrapRGB = new Vector3(1,1,1);

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
        this.shaderID = "lambert";
        this.vertexColors = NoColors;

        this.skinning = false;

        this.setValues(parameters);

    }

    clone(material) {

        if ( typeof material === "undefined" ) material = new MeshLambertMaterial();

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

MeshLambertMaterial.prototype = Object.create(Material.prototype);

//Double sided Mesh Lambert material
/** @constructor */
export class MeshDoubleLambertMaterial {
    constructor(parameters) {

        MeshLambertMaterial.call(this, parameters);

        this.shaderID = "lambertdouble";
        this.side = DoubleSide;

    }

    clone() {

        var material = new MeshDoubleLambertMaterial();

        MeshLambertMaterial.prototype.clone.call(this, material);

        return material;

    }
}

MeshDoubleLambertMaterial.prototype = Object.create(MeshLambertMaterial.prototype);

//Outlined Mesh Lamert material
/** @constructor */
export class MeshOutlineMaterial {
    constructor(parameters) {
        Material.call(this);
        parameters = parameters || {};
        this.fog = true;
        this.shaderID = "outline";
        this.wireframe=false;
        this.outlineColor= parameters.color || new Color(0.0,0.0,0.0);
        this.outlineWidth= parameters.width || 0.1;
        this.outlinePushback= parameters.pushback || 1.0;

    }

    clone(material) {
        if ( typeof material === "undefined" ) material = new MeshOutlineMaterial();
        Material.prototype.clone.call(this, material);
        material.fog = this.fog;
        material.shaderID = this.shaderID;
        material.wireframe = this.wireframe;
        return material;
    }
}

MeshOutlineMaterial.prototype = Object.create(Material.prototype);


//Imposter material
/** @constructor */
export class ImposterMaterial {
    constructor(parameters) {

      Material.call(this);

      this.color = new Color(0xffffff);
      this.ambient = new Color(0xfffff);
      this.emissive = new Color(0x000000);
      this.imposter = true;

      //TODO: Which of these instance variables do I really need?
      this.wrapAround = false;
      this.wrapRGB = new Vector3(1,1,1);

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

    clone() {

      var material = new ImposterMaterial();

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

ImposterMaterial.prototype = Object.create(Material.prototype);


export class SphereImposterMaterial {
    constructor(parameters) {

        ImposterMaterial.call(this);

        this.shaderID = "sphereimposter";
        this.setValues(parameters);

    }

    clone() {

        var material = new SphereImposterMaterial();
        ImposterMaterial.prototype.clone.call(this, material);
        return material;
    }
}

SphereImposterMaterial.prototype = Object.create(ImposterMaterial.prototype);


export class SphereImposterOutlineMaterial {
    constructor(parameters) {

        ImposterMaterial.call(this);
        parameters = parameters || {};

        this.shaderID = "sphereimposteroutline";
        this.outlineColor= parameters.color || new Color(0.0,0.0,0.0);
        this.outlineWidth= parameters.width || 0.1;
        this.outlinePushback= parameters.pushback || 1.0;

        this.setValues(parameters);

    }

    clone() {

        var material = new SphereImposterOutlineMaterial();
        ImposterMaterial.prototype.clone.call(this, material);
        material.outlineColor = this.outlineColor;
        material.outlineWidth = this.outlineWidth;
        material.outlinePushback = this.outlinePushback;
        return material;
    }
}

SphereImposterOutlineMaterial.prototype = Object.create(ImposterMaterial.prototype);


export class StickImposterMaterial {
    constructor(parameters) {

        ImposterMaterial.call(this);

        this.shaderID = "stickimposter";
        this.setValues(parameters);

    }

    clone() {

        var material = new StickImposterMaterial();
        ImposterMaterial.prototype.clone.call(this, material);
        return material;
    }
}

StickImposterMaterial.prototype = Object.create(ImposterMaterial.prototype);


export class StickImposterOutlineMaterial {
    constructor(parameters) {

        ImposterMaterial.call(this);
        parameters = parameters || {};

        this.shaderID = "stickimposteroutline";
        this.outlineColor= parameters.color || new Color(0.0,0.0,0.0);
        this.outlineWidth= parameters.width || 0.1;
        this.outlinePushback= parameters.pushback || 1.0;

        this.setValues(parameters);

    }

    clone() {

        var material = new StickImposterOutlineMaterial();
        ImposterMaterial.prototype.clone.call(this, material);
        material.outlineColor = this.outlineColor;
        material.outlineWidth = this.outlineWidth;
        material.outlinePushback = this.outlinePushback;
        return material;
    }
}

StickImposterOutlineMaterial.prototype = Object.create(ImposterMaterial.prototype);


export class InstancedMaterial {
    constructor(parameters) {

        Material.call(this);

        this.color = new Color(0xffffff);
        this.ambient = new Color(0xfffff);
        this.emissive = new Color(0x000000);

        //TODO: Which of these instance variables do I really need?
        this.wrapAround = false;
        this.wrapRGB = new Vector3(1,1,1);

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
        this.shaderID = "instanced";
        this.vertexColors = NoColors;

        this.skinning = false;

        this.sphere = null;

        this.setValues(parameters);

    }

    clone() {

        var material = new InstancedMaterial();

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

InstancedMaterial.prototype = Object.create(Material.prototype);


//Volumetric material
/** @constructor */
export class VolumetricMaterial {
    constructor(parameters) {

        Material.call(this);

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

        this.shaderID = "volumetric";
        this.side = FrontSide;

        this.setValues(parameters);
    }

    clone() {

        var material = Object.assign(new VolumetricMaterial(),this);

        Material.prototype.clone.call(this, material);
        return material;

    }
}

VolumetricMaterial.prototype = Object.create(Material.prototype);

//Sprite material
/** @constructor */
export class SpriteMaterial {
    constructor(parameters) {

        Material.call(this);

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

        if (parameters.depthTest === undefined)
            this.depthTest = !this.useScreenCoordinates;
        if (parameters.sizeAttenuation === undefined)
            this.sizeAttenuation = !this.useScreenCoordinates;
        if (parameters.scaleByViewPort === undefined)
            this.scaleByViewPort = !this.sizeAttenuation;

    }

    clone() {

        var material = new SpriteMaterial();

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

SpriteMaterial.prototype = Object.create(Material.prototype);

//Alignment for Sprites
export var SpriteAlignment = {};
SpriteAlignment.topLeft = new Vector2(1, -1);
SpriteAlignment.topCenter = new Vector2(0, -1);
SpriteAlignment.topRight = new Vector2(-1, -1);
SpriteAlignment.centerLeft = new Vector2(1, 0);
SpriteAlignment.center = new Vector2(0, 0);
SpriteAlignment.centerRight = new Vector2(-1, 0);
SpriteAlignment.bottomLeft = new Vector2(1, 1);
SpriteAlignment.bottomCenter = new Vector2(0, 1);
SpriteAlignment.bottomRight = new Vector2(-1, 1);


//Texture
//We really only create textures from 2d rendering contexts (to display text labels)
//edit: we can now create 3dtextures using volumetric data
/** @constructor */
export class Texture {
    constructor(image, is3D) {

        EventDispatcher.call(this);

        this.id = TextureIdCount++;

        this.name = "";

        this.image = image;

        this.mapping = new UVMapping();

        this.wrapS = ClampToEdgeWrapping;
        this.wrapT = ClampToEdgeWrapping;

        this.anisotropy = 1;

        if (is3D){
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
}

Texture.prototype = {

    constructor : Texture,

    clone : function(texture) {

        if (texture === undefined)
            texture = new Texture();

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

    },

    dispose : function() {

        this.dispatchEvent( {type: 'dispose'});

    }

};

export var TextureIdCount = 0;


// sides
export var FrontSide = 0;
export var BackSide = 1;
export var DoubleSide = 2;

// shading
export var NoShading = 0;
export var FlatShading = 1;
export var SmoothShading = 2;

// colors
export var NoColors = 0;
export var FaceColors = 1;
export var VertexColors = 2;

//Texture constants
//TODO: Which of these do I need (since I only use textures to display label sprites) ?
export var MultiplyOperation = 0;
export var MixOperation = 1;
export var AddOperation = 2;

// mapping modes

export class UVMapping {
    constructor() {}
}

// wrapping modes
export var ClampToEdgeWrapping = 1001;

//Filters
export var LinearFilter = 1006;
export var NearestFilter = 1007;
export var LinearMipMapLinearFilter = 1008;

//Data types
export var UnsignedByteType = 1009;
export var FloatType = 1010;

//Pixel formats
export var RGBAFormat = 1021;
export var RFormat = 1022;
export var R32Format = 1023;
