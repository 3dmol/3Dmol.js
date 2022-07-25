/*
 * math-like functionality
 * quaternion, vector, matrix
 */
export function clamp(x: number, min: number, max: number) {
  return Math.min(Math.max(x, min), max);
}