import { Mesh } from './objects/Mesh';
import { Line, LineStyle } from './objects/Line';
import { Scene } from './Scene';
import { Light } from './Light';
import { Camera } from './Camera';
import { VolumetricMaterial } from './materials/VolumetricMaterial';
import { StickImposterOutlineMaterial } from './materials/StickImposterOutlineMaterial';
import { StickImposterMaterial } from './materials/StickImposterMaterial';
import { SpriteMaterial } from './materials/SpriteMaterial';
import { SphereImposterOutlineMaterial } from './materials/SphereImposterOutlineMaterial';
import { SphereImposterMaterial } from './materials/SphereImposterMaterial';
import { MeshOutlineMaterial } from './materials/MeshOutlineMaterial';
import { MeshLambertMaterial } from './materials/MeshLambertMaterial';
import { MeshDoubleLambertMaterial } from './materials/MeshDoubleLambertMaterial';
import { Material, MaterialIdCount } from './materials/Material';
import { LineBasicMaterial } from './materials/LineBasicMaterial';
import { InstancedMaterial } from './materials/InstancedMaterial';
import { ImposterMaterial } from './materials/ImposterMaterial';
import { UVMapping } from './core/UVMapping';
import { Texture, TextureIdCount } from './core/Texture';
import { Raycaster } from './core/Raycaster';
import { Projector } from './core/Projector';
import { Object3D, Object3DIDCount } from './core/Object3D';
import { Geometry, GeometryIDCount, GeometryGroup } from './core/Geometry';
import { EventDispatcher } from './core/EventDispatcher';
import { TextureOperations } from './constants/TextureOperations';
import { SpriteAlignment } from './constants/SpriteAlignment';
import { Sides } from './constants/Sides';
import { Shading } from './constants/Shading';
import { SpritePlugin } from './SpritePlugin';
import { Color, CC } from './core/Color';
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
import { Colors } from './constants/Colors';
import { ClampToEdgeWrapping, LinearFilter, NearestFilter, LinearMipMapLinearFilter, UnsignedByteType, FloatType, RGBAFormat, RFormat, R32Format } from './constants/TextureConstants';
import { intersectObject } from './core/intersectObject';
import { MeshBasicMaterial } from './materials/MeshBasicMaterial';
import { Sprite } from './objects/Sprite';


(window as Record<string, any>).$3Dmol = {
  ...((window as Record<string, any>).$3Dmol || {}),
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
  SpritePlugin,
  ShaderLib,
  ShaderUtils,
  NoColors: Colors.NoColors,
  FaceColors: Colors.FaceColors,
  VertexColors: Colors.VertexColors,
  NoShading: Shading.NoShading,
  FlatShading: Shading.FlatShading,
  SmoothShading: Shading.SmoothShading,
  FrontSide: Sides.FrontSide,
  BackSide: Sides.BackSide,
  DoubleSide: Sides.DoubleSide,
  SpriteAlignment,
  ClampToEdgeWrapping,
  LinearFilter,
  NearestFilter,
  LinearMipMapLinearFilter,
  UnsignedByteType,
  FloatType,
  RGBAFormat,
  RFormat,
  R32Format,
  MultiplyOperation: TextureOperations.MultiplyOperation,
  AddOperation: TextureOperations.AddOperation,
  MixOperation: TextureOperations.MixOperation,
  Color,
  CC,
  EventDispatcher,
  GeometryIDCount,
  Geometry,
  GeometryGroup,
  intersectObject,
  Object3D,
  Object3DIDCount,
  Projector,
  Raycaster,
  Texture,
  TextureIdCount,
  UVMapping,
  ImposterMaterial,
  InstancedMaterial,
  LineBasicMaterial,
  Material,
  MaterialIdCount,
  MeshDoubleLambertMaterial,
  MeshLambertMaterial,
  MeshOutlineMaterial,
  SphereImposterMaterial,
  SphereImposterOutlineMaterial,
  SpriteMaterial,
  StickImposterMaterial,
  StickImposterOutlineMaterial,
  VolumetricMaterial,
  MeshBasicMaterial,
  Camera,
  Light,
  Scene,
  Line,
  Mesh,
  Sprite,
  LineStrip: LineStyle.LineStrip,
  LinePieces: LineStyle.LinePieces,
};
