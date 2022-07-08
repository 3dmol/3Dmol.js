import { Colors } from './../constants/Colors';
import { Shading } from './../constants/Shading';
import { Color } from "../core/Color";
import { Vector3 } from "../math";
import { Material } from "./Material";
//Imposter material
/** @constructor */
export class ImposterMaterial extends Material {
  color: any;
  ambient: any;
  emissive: any;
  imposter: boolean;
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
  shaderID: unknown;
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
    this.imposter = true;

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
    this.shaderID = null;
    this.vertexColors = Colors.NoColors;

    this.skinning = false;

    this.setValues(parameters);
  }
  clone() {
    var material = new ImposterMaterial();

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
