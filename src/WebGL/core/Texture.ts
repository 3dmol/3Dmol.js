//Texture
//We really only create textures from 2d rendering contexts (to display text labels)
//edit: we can now create 3dtextures using volumetric data

import {
  ClampToEdgeWrapping,
  RFormat,
  FloatType,
  NearestFilter,
  RGBAFormat,
  UnsignedByteType,
  LinearFilter,
  LinearMipMapLinearFilter,
} from "../constants/TextureConstants";
import { Vector2 } from "../math";
import { EventDispatcher } from "./EventDispatcher";
import { UVMapping } from "./UVMapping";

/* @constructor */
export class Texture extends EventDispatcher {
  id: number;
  name: string;
  image: any;
  mapping: any;
  wrapS: number;
  wrapT: number;
  anisotropy: number;
  format: number;
  type: number;
  premultiplyAlpha: boolean;
  flipY: boolean;
  unpackAlignment: number;
  magFilter: number;
  minFilter: number;
  offset: any;
  repeat: any;
  needsUpdate: boolean;
  onUpdate: null;
  constructor(image?: any, is3D?: boolean) {
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

  clone(texture = new Texture()): Texture {

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
    this.dispatchEvent({ type: "dispose" });
  }
}
export let TextureIdCount = 0;
