import { TextureOperations } from "./../constants/TextureOperations";
import { Color } from "../core/Color";
import { Material } from "./Material";

/**
 * This class is ripped right out of three.js.
 * It was missing in 3Dmol despite being referenced in Mesh.ts
 */
class MeshBasicMaterial extends Material {
  isMeshBasicMaterial = true;
  type = "MeshBasicMaterial";
  color = new Color(0xffffff); // emissive
  map = undefined;
  lightMap = null;
  lightMapIntensity = 1.0;
  aoMap = null;
  aoMapIntensity = 1.0;
  specularMap = null;
  alphaMap = null;
  envMap = null;
  combine = TextureOperations.MultiplyOperation;
  reflectivity = 1;
  refractionRatio = 0.98;
  wireframe = false;
  wireframeLinewidth = 1;
  wireframeLinecap = "round";
  wireframeLinejoin = "round";
  fog = true;

  constructor(values: Partial<Record<keyof MeshBasicMaterial, any>> = {} as any) {
    super();
    this.setValues(values);
  }

	clone<T extends this>(material = new MeshBasicMaterial() as T): T {
    super.clone.call(this, material);

    material.color.copy(this.color);
    material.map = this.map
    material.lightMap = this.lightMap
    material.lightMapIntensity = this.lightMapIntensity
    material.aoMap = this.aoMap
    material.aoMapIntensity = this.aoMapIntensity
    material.specularMap = this.specularMap
    material.alphaMap = this.alphaMap
    material.envMap = this.envMap
    material.combine = this.combine
    material.reflectivity = this.reflectivity
    material.refractionRatio = this.refractionRatio
    material.wireframe = this.wireframe
    material.wireframeLinewidth = this.wireframeLinewidth
    material.wireframeLinecap = this.wireframeLinecap
    material.wireframeLinejoin = this.wireframeLinejoin
    material.fog = this.fog

    return material;
  }
}

export { MeshBasicMaterial };
