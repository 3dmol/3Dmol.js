import { SpriteAlignment } from "../constants/SpriteAlignment";
import { Texture } from "../core/Texture";
import { Color } from "../../colors";
import { Vector2 } from "../math";
import { Material } from "./Material";
export class SpriteMaterial extends Material {
  sizeAttenuation: boolean;
  screenOffset: any;
  scaleByViewPort: boolean;
  alignment: any;
  scaleByViewport: any;

  color = new Color(0xffffff);
  map = new Texture();
  useScreenCoordinates = true;
  fog = false; // use scene fog
  uvOffset = new Vector2(0, 0);
  uvScale = new Vector2(1, 1);
  
  constructor(parameters?: any) {
    super();
    this.depthTest = !this.useScreenCoordinates;
    this.sizeAttenuation = !this.useScreenCoordinates;
    this.screenOffset = this.screenOffset;
    this.scaleByViewPort = !this.sizeAttenuation;
    this.alignment = SpriteAlignment.center.clone();

    this.setValues(parameters);

    parameters = parameters || {};

    if (parameters.depthTest === undefined)
      this.depthTest = !this.useScreenCoordinates;
    if (parameters.sizeAttenuation === undefined)
      this.sizeAttenuation = !this.useScreenCoordinates;
    if (parameters.scaleByViewPort === undefined)
      this.scaleByViewPort = !this.sizeAttenuation;
  }

  clone<T extends this>(material = new SpriteMaterial() as T): T {
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
  }
}
