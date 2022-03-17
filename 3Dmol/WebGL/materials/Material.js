// @ts-check

import { Color, EventDispatcher } from "../core/index";
import { Vector3 } from "../math/index";
import { FrontSide } from "./sides";

export var MaterialIdCount = 0;

/**
 * Line and Mesh material types
 * @constructor
 */
export class Material extends EventDispatcher {
  /**
   * @type {any}
   */
  combine;
  /**
   * @type {any}
   */
  morphTargets;
  /**
   * @type {any}
   */
  morphNormals;
  /**
   * @type {any}
   */
  overdraw;

  id = MaterialIdCount++;
  name = "";
  //TODO: Which of these instance variables can I remove??
  side = FrontSide;
  opacity = 1;
  transparent = false;
  depthTest = true;
  depthWrite = true;
  stencilTest = true;
  polygonOffset = false;
  polygonOffsetFactor = 0;
  polygonOffsetUnits = 0;
  alphaTest = 0;
  visible = true;
  needsUpdate = true;
  /**
   * @param {Record<any, any> | undefined} [values]
   */
  setValues(values) {
    if (values === undefined) return;

    for (var key in values) {
      var newValue = values[key];

      if (newValue === undefined) {
        console.warn("Material: '" + key + "' parameter is undefined.");
        continue;
      }

      if (key in this) {
        // @ts-ignore
        var currentValue = this[key];

        if (currentValue instanceof Color && newValue instanceof Color) {
          currentValue.copy(newValue);
        } else if (currentValue instanceof Color) {
          currentValue.set(newValue);
        } else if (
          currentValue instanceof Vector3 &&
          newValue instanceof Vector3
        ) {
          currentValue.copy(newValue);
        } else {
          // @ts-ignore
          this[key] = newValue;
        }
      }
    }
  }

  /**
   * @param {Material | undefined} material
   */
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
    this.dispatchEvent({ type: "dispose" });
  }
}
