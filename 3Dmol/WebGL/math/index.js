// @ts-check

import { Matrix3 } from "./Matrix3";

export * from "./Matrix3";
export * from "./Matrix4";
export * from "./Quaternion";
export * from "./Ray";
export * from "./Vector2";
export * from "./Vector3";

// patch the math object
/*
 * math-like functionality
 * quaternion, vector, matrix
 */

export function clamp(x, min, max) {
  return Math.min(Math.max(x, min), max);
}

let degreeToRadiansFactor = Math.PI / 180;
export function degToRad(deg) {
  return deg * degreeToRadiansFactor;
}

//return conversion matrix given crystal unit cell parameters
export function conversionMatrix3(a, b, c, alpha, beta, gamma) {
  //convert to radians
  alpha = (alpha * Math.PI) / 180;
  beta = (beta * Math.PI) / 180;
  gamma = (gamma * Math.PI) / 180;
  let sqr = (a) => a * a;
  let cos_alpha = Math.cos(alpha);
  let cos_beta = Math.cos(beta);
  let cos_gamma = Math.cos(gamma);
  let sin_gamma = Math.sin(gamma);
  let conversionMatrix = new Matrix3(
    a,
    b * cos_gamma,
    c * cos_beta,
    0,
    b * sin_gamma,
    (c * (cos_alpha - cos_beta * cos_gamma)) / sin_gamma,
    0,
    0,
    (c *
      Math.sqrt(
        1 -
          sqr(cos_alpha) -
          sqr(cos_beta) -
          sqr(cos_gamma) +
          2 * cos_alpha * cos_beta * cos_gamma
      )) /
      sin_gamma
  );
  return conversionMatrix;
}
