import { Sides } from './../constants/Sides';
import { Color } from "../core/Color";
import { EventDispatcher } from "../core/EventDispatcher";
import { Vector3 } from "../math";
/**
 * Line and Mesh material types
 * @constructor
 */
export class Material extends EventDispatcher {
  copy(source: any) {
    throw new Error('Method not implemented.');
  }
  id: number;
  name: string;
  side: number;
  opacity: number;
  transparent: boolean;
  depthTest: boolean;
  depthWrite: boolean;
  stencilTest: boolean;
  polygonOffset: boolean;
  polygonOffsetFactor: number;
  polygonOffsetUnits: number;
  alphaTest: number;
  visible: boolean;
  needsUpdate: boolean;
  overdraw: any;
  constructor() {
    super();
    this.id = MaterialIdCount++;

    this.name = "";

    //TODO: Which of these instance variables can I remove??
    this.side = Sides.FrontSide;

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

    for (var key in values) {
      var newValue = values[key];

      if (newValue === undefined) {
        console.warn("$3Dmol.Material: '" + key + "' parameter is undefined.");
        continue;
      }

      if (key in this) {
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
          this[key] = newValue;
        }
      }
    }
  }

  //TODO: might want to look into blending equations
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

export let MaterialIdCount = 0;