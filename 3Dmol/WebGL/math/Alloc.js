import { Matrix4 } from "./Matrix4";
import { Vector3 } from "./Vector3";

// tmp allocation for lookAt
export const x = new Vector3();
export const y = new Vector3();
export const z = new Vector3();
// tmp allocation for getPosition
export const v1 = new Vector3();
// tmp allocation for decompose
// note uses the same private temp vectors as Matrix4.lookAt()
export const matrix = new Matrix4();
// tmp matrix allocation for compose
export const mRotation = new Matrix4();
export const mScale = new Matrix4();