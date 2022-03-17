// @ts-check

import { Color } from "../core/Color";
import { Vector2 } from "../math/Vector2";
import { Material } from "./Material";
import { Texture } from "./Texture";


//Alignment for Sprites
export var SpriteAlignment = {
  topLeft: new Vector2(1, -1),
  topCenter: new Vector2(0, -1),
  topRight: new Vector2(-1, -1),
  centerLeft: new Vector2(1, 0),
  center: new Vector2(0, 0),
  centerRight: new Vector2(-1, 0),
  bottomLeft: new Vector2(1, 1),
  bottomCenter: new Vector2(0, 1),
  bottomRight: new Vector2(-1, 1),
}

//Sprite material
/** @constructor */
export class SpriteMaterial extends Material {
  scaleByViewport = true;
  /**
   * @param {Record<any, any> | undefined} [parameters]
   */
  constructor(parameters) {

    super();

    this.color = new Color(0xffffff);
    this.map = new Texture();

    this.useScreenCoordinates = true;
    this.depthTest = !this.useScreenCoordinates;
    this.sizeAttenuation = !this.useScreenCoordinates;
    /**
     * @type {undefined}
     */
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

  /**
   * 
   * @param {SpriteMaterial} [material] 
   * @returns SpriteMaterial
   */
  clone(material) {

    if (material === undefined) material = new SpriteMaterial();

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