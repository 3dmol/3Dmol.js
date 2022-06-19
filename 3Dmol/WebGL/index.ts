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
} from "./math/index";
import { square } from "./math/utils/square";
import { ShaderLib, ShaderUtils } from "./shaders/index";

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
