const degreeToRadiansFactor = Math.PI / 180;

export function degToRad(deg: number) {
  return deg * degreeToRadiansFactor;
}