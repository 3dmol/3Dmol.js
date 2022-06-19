import { ShaderUtils } from './shaders/utils/index';
import { ShaderLib } from './shaders/lib/index';
import { SpritePlugin } from './SpritePlugin';
import { Color } from './core/Color';
import { Sphere, Cylinder, Triangle } from "./shapes"
import {
  Matrix4,
  Matrix3,
  Quaternion,
  Ray,
  conversionMatrix3,
  Vector2,
  Vector3,
  clamp,
  degToRad,
} from "./math";
import { square } from "./math/utils/square";

// @ts-ignore
window.$3Dmol = {
  // @ts-ignore
  ...(window.$3Dmol || {}),
  // @ts-ignore
  Matrix3,
  Matrix4,
  Quaternion,
  Ray,
  conversionMatrix3,
  Vector2,
  Vector3,
  Math: {
    clamp,
    degToRad,
  },
  square,
  Cylinder,
  Sphere,
  Triangle,
  Color,
  SpritePlugin,
  ShaderLib,
  ShaderUtils
};
