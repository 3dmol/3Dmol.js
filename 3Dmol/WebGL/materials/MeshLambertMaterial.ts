import { Colors } from './../constants/Colors';
import { Shading } from './../constants/Shading';
import { Material } from "./Material";
//Mesh Lambert material

import { Color } from "../core/Color";
import { Vector3 } from "../math";

/** @constructor */
export class MeshLambertMaterial extends Material {
  color: any;
  ambient: any;
  emissive: any;
  wrapAround: boolean;
  wrapRGB: any;
  map: null;
  lightMap: null;
  specularMap: null;
  envMap: null;
  reflectivity: number;
  refractionRatio: number;
  fog: boolean;
  wireframe: boolean;
  wireframeLinewidth: number;
  wireframeLinecap: string;
  wireframeLinejoin: string;
  shading: number;
  shaderID: string;
  vertexColors: number;
  skinning: boolean;
  combine: any;
  morphTargets: any;
  morphNormals: any;
  constructor(parameters?: any) {
    super();

    this.color = new Color(0xffffff);
    this.ambient = new Color(0xfffff);
    this.emissive = new Color(0x000000);

    //TODO: Which of these instance variables do I really need?
    this.wrapAround = false;
    this.wrapRGB = new Vector3(1, 1, 1);

    this.map = null;

    this.lightMap = null;

    this.specularMap = null;

    this.envMap = null;
    this.reflectivity = 1;
    this.refractionRatio = 0.98;

    this.fog = true;

    this.wireframe = false;
    this.wireframeLinewidth = 1;
    this.wireframeLinecap = "round";
    this.wireframeLinejoin = "round";

    this.shading = Shading.SmoothShading;
    this.shaderID = "lambert";
    this.vertexColors = Colors.NoColors;

    this.skinning = false;

    this.setValues(parameters);
  }

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
