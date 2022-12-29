// return conversion matrix given crystal unit cell parameters

import { Matrix3 } from "../math";

/**
 *
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @param {number} alpha
 * @param {number} beta
 * @param {number} gamma
 * @returns {Matrix3}
 */
export function conversionMatrix3(
  a: any,
  b: number,
  c: number,
  alpha: number,
  beta: number,
  gamma: number
) {
  // convert to radians
  alpha = (alpha * Math.PI) / 180;
  beta = (beta * Math.PI) / 180;
  gamma = (gamma * Math.PI) / 180;
  const sqr = (x) => x*x;
  const cosAlpha = Math.cos(alpha);
  const cosBeta = Math.cos(beta);
  const cosGamma = Math.cos(gamma);
  const sinGamma = Math.sin(gamma);
  const conversionMatrix = new Matrix3(
    a,
    b * cosGamma,
    c * cosBeta,
    0,
    b * sinGamma,
    (c * (cosAlpha - cosBeta * cosGamma)) / sinGamma,
    0,
    0,
    (c *
      Math.sqrt(
        1 -
          sqr(cosAlpha) -
          sqr(cosBeta) -
          sqr(cosGamma) +
          2 * cosAlpha * cosBeta * cosGamma
      )) /
      sinGamma
  );
  return conversionMatrix;
}
