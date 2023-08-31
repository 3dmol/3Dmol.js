import {FrontSide}  from "../constants/Sides";
import { EventDispatcher } from "../core";
import type { Texture } from "../core";
import { Vector2, Vector3 } from "../math";
import { Color } from "../../colors";
/*
 * Line and Mesh material types
 * @constructor
 */
export class Material extends EventDispatcher {
  id = MaterialIdCount++;
  name = "";
  overdraw: any;
  color?: Color;
  map?: Texture;
  useScreenCoordinates?: boolean;
  alignment?: Vector2;
  screenOffset?: Vector2;
  uvScale?: Vector2;
  uvOffset?: Vector2;
  scaleByViewport?: boolean;
  fog?: unknown;

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
  outline = false;

  setValues(
    values: Partial<Record<keyof Material, any>> = {} as any
  ) {
    if (values === undefined) return;

    for (var key in values) {
      var newValue: Material[keyof Material] = values[key as keyof Material];

      if (newValue === undefined) {
        console.warn("$3Dmol.Material: '" + key + "' parameter is undefined.");
        continue;
      }

      if (key in this) {
        var currentValue = this[key as keyof Material];

        if (currentValue instanceof Color && newValue instanceof Color) {
          currentValue.copy(newValue);
        } else if (currentValue instanceof Color) {
          currentValue.set(newValue as unknown as Color);
        } else if (
          currentValue instanceof Vector3 &&
          newValue instanceof Vector3
        ) {
          currentValue.copy(newValue);
        } else {
          this[key] = newValue;
        }
      }
    }
  }

  //TODO: might want to look into blending equations
  clone<T extends this>(material = new Material() as T): T {
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

export let MaterialIdCount = 0;
