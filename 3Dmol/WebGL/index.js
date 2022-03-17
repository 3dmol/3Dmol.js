import {
  Color,
  EventDispatcher,
  Geometry,
  GeometryGroup,
  Object3D,
  Projector,
  Raycaster,
  Scene
} from "./core";
import { Fog } from "./Fog";
import {
  BackSide,
  DoubleSide, FaceColors, FlatShading, FrontSide, ImposterMaterial, InstancedMaterial, LineBasicMaterial, Material, MeshDoubleLambertMaterial, MeshLambertMaterial, MeshOutlineMaterial, NoColors, NoShading, SmoothShading, SphereImposterMaterial,
  SphereImposterOutlineMaterial, SpriteMaterial, StickImposterMaterial,
  StickImposterOutlineMaterial, VertexColors, VolumetricMaterial
} from "./materials";
import {
  clamp,
  degToRad,
  Matrix3,
  Matrix4,
  Quaternion,
  Ray,
  Vector2,
  Vector3
} from "./math";
import { Camera, Light, Line, Mesh, Sprite } from "./objects";
import { Renderer } from "./Renderer";
import { ShaderUtils } from "./ShaderUtils";
import { Cylinder, Sphere, Triangle } from "./shapes";
import { SpritePlugin } from "./SpritePlugin";

globalThis.$3Dmol = globalThis.$3Dmol || {};
globalThis.$3Dmol = {
  ...globalThis.$3Dmol,
  Object3D,
  Color,
  EventDispatcher,
  Geometry,
  GeometryGroup,
  Projector,
  Raycaster,
  Scene,
  NoColors,
  FaceColors,
  VertexColors,
  Material,
  LineBasicMaterial,
  MeshLambertMaterial,
  MeshDoubleLambertMaterial,
  MeshOutlineMaterial,
  ImposterMaterial,
  SphereImposterMaterial,
  SphereImposterOutlineMaterial,
  StickImposterMaterial,
  StickImposterOutlineMaterial,
  InstancedMaterial,
  VolumetricMaterial,
  SpriteMaterial,
  NoShading,
  FlatShading,
  SmoothShading,
  FrontSide,
  BackSide,
  DoubleSide,
  Light,
  Fog,
  SpritePlugin,
  Matrix3,
  Matrix4,
  Quaternion,
  Ray,
  Vector2,
  Vector3,
  Math: {
    clamp,
    degToRad,
  },
  Triangle,
  Sphere,
  Cylinder,
  Camera,
  Renderer,
  ShaderUtils,
  Line,
  Mesh,
  Sprite,
};
