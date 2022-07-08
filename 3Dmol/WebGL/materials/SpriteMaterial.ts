import { SpriteAlignment } from '../constants/SpriteAlignment';
import { Texture } from '../core/Texture';
import { Color } from '../core/Color';
import { Vector2 } from '../math';
import { Material } from './Material';
export class SpriteMaterial extends Material {
  color: Color;
  map: any;
  useScreenCoordinates: boolean;
  sizeAttenuation: boolean;
  screenOffset: any;
  scaleByViewPort: boolean;
  alignment: any;
  fog: boolean;
  uvOffset: Vector2;
  uvScale: Vector2;
  scaleByViewport: any; 
  constructor(parameters?: any) {
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

    if (parameters.depthTest === undefined)
      this.depthTest = !this.useScreenCoordinates;
    if (parameters.sizeAttenuation === undefined)
      this.sizeAttenuation = !this.useScreenCoordinates;
    if (parameters.scaleByViewPort === undefined)
      this.scaleByViewPort = !this.sizeAttenuation;

  };

  clone() {

    var material = new SpriteMaterial();

    super.clone.call(this, material);

    material.color.copy(this.color);
    material.map = this.map;

    material.useScreenCoordinates = this.useScreenCoordinates;
    material.screenOffset = this.screenOffset;
    material.sizeAttenuation = this.sizeAttenuation;
    material.scaleByViewport = this.scaleByViewPort;
    material.alignment.copy(this.alignment);

    material.uvOffset.copy(this.uvOffset);

    return material;

  };
}