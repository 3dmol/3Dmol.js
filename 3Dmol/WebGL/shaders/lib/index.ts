export * from "./basic"
export * from "./instanced"
export * from "./lambert"
export * from "./lambertdouble"
export * from "./outline"
export * from "./screen"
export * from "./screenaa"
export * from "./sphereimposter"
export * from "./sphereimposteroutline"
export * from "./sprite"
export * from "./stickimposter"
export * from "./stickimposteroutline"
export * from "./volumetric"

import { Shader } from '../ShaderType';
import { basic } from "./basic"
import { instanced } from "./instanced"
import { lambert } from "./lambert"
import { lambertdouble } from "./lambertdouble"
import { outline } from "./outline"
import { screen } from "./screen"
import { screenaa } from "./screenaa"
import { sphereimposter } from "./sphereimposter"
import { sphereimposteroutline } from "./sphereimposteroutline"
import { sprite } from "./sprite"
import { stickimposter } from "./stickimposter"
import { stickimposteroutline } from "./stickimposteroutline"
import { volumetric } from "./volumetric"


export const ShaderLib: Record<string, Shader> = {
    basic,
    instanced,
    lambert,
    lambertdouble,
    outline,
    screen,
    screenaa,
    sphereimposter,
    sphereimposteroutline,
    sprite,
    stickimposter,
    stickimposteroutline,
    volumetric,
}