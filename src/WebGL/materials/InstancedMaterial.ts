import { Coloring } from "../constants/Coloring";
import { Shading } from "../constants/Shading";
import { Color } from "../../colors";
import { Vector3 } from "../math";
import { Material } from "./Material";

export class InstancedMaterial extends Material {
  combine: any;
  morphTargets: any;
  morphNormals: any;

  color = new Color(0xffffff);
  ambient = new Color(0xfffff);
  emissive = new Color(0x000000);

  //TODO: Which of these instance variables do I really need?
  wrapAround = false;
  wrapRGB = new Vector3(1, 1, 1);
  map = undefined;
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
  shading = Shading.SmoothShading;
  shaderID = "instanced";
  vertexColors = Coloring.NoColors;
  skinning = false;
  sphere = null;
  constructor(parameters?: any) {
    super();
    this.setValues(parameters);
  }

  clone<T extends this>(material: T = new InstancedMaterial() as T): T {
    super.clone.call(this, material as T);

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
