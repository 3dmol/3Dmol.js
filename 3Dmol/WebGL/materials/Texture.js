// @ts-check

import { EventDispatcher } from "../core/EventDispatcher";
import { Vector2 } from "../math/Vector2";



//Texture constants
//TODO: Which of these do I need (since I only use textures to display label sprites) ?
export var MultiplyOperation = 0;
export var MixOperation = 1;
export var AddOperation = 2;

// mapping modes

export class UVMapping {}

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


//Texture
//We really only create textures from 2d rendering contexts (to display text labels)
//edit: we can now create 3dtextures using volumetric data
/** @constructor */
export class Texture extends EventDispatcher {
    constructor(image, is3D) {
  
      super();
  
      this.id = TextureIdCount++;
  
      this.name = "";
  
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
  
    clone(texture) {
  
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
    }
  
    dispose() {
  
      this.dispatchEvent({ type: 'dispose' });
  
    }
  }
  
  export var TextureIdCount = 0;