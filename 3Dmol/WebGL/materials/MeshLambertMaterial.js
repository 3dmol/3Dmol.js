// @ts-check

import { Color } from "../core/Color";
import { Vector3 } from "../math/Vector3";
import { NoColors } from "./colors";
import { Material } from "./Material";
import { SmoothShading } from "./shading";

//Mesh Lambert material
export class MeshLambertMaterial extends Material {
  color = new Color(0xffffff);
  ambient = new Color(0xfffff);
  emissive = new Color(0x000000);
  //TODO: Which of these instance variables do I really need?
  wrapAround = false;
  wrapRGB = new Vector3(1, 1, 1);
  map = null;
  lightMap = null;
  specularMap = null;
  envMap = null;
  reflectivity = 1;
  refractionRatio = 0.98;
  fog = true;
  wireframe = false;
  wireframeLinewidth = 1;
  wireframeLinecap = "round";
  wireframeLinejoin = "round";
  shading = SmoothShading;
  shaderID = "lambert";
  vertexColors = NoColors;
  skinning = false;

  /**
   * @param {Record<any, any> | undefined} [parameters]
   */
  constructor(parameters) {
    super();
    this.setValues(parameters);
  }

  /**
   * @param {MeshLambertMaterial | undefined} material
   */
  clone(material) {
    if (typeof material === "undefined") material = new MeshLambertMaterial();

    super.clone.call(this, material);

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
