/* eslint-disable import/prefer-default-export */
import "./effects/binaryPlusAjaxExtension"
import "./effects/autoloadEffect";
import $ from "jquery"
import {CC, elementColors, htmlColors, ssColors} from './colors';
import drawCartoon from './drawCartoon';
import GLDraw from './GLDraw';
import GLModel from './GLModel';
import GLShape from './glshape';
import GLViewer from './GLViewer';
import Gradient from './Gradient';
import Label, {LabelCount} from './Label';
import Color from './WebGL/core/Color';
import Parsers from './Parsers';
import MarchingCube, {MarchingCubeInitializer} from './MarchingCube';
import {VolumeData, GLVolumetricRender} from './volume';
import Camera from './WebGL/Camera';
import Fog from './WebGL/Fog';
import ProteinSurface from './ProteinSurface4';
import {
  Geometry,
  Object3D,
  GeometryIDCount,
  Light,
  Object3DIDCount,
  Projection,
  Raycaster,
  Scene,
  geometryGroup,
} from './WebGL/core';
import {
  LineBasicMaterial,
  Material,
  MeshDoubleLambertMaterial,
  MeshLambertMaterial,
  SphereImposterMaterial,
  AddOperation,
  BackSide,
  ClampToEdgeWrapping,
  DoubleSide,
  FaceColors,
  FlatShading,
  FloatType,
  FrontSide,
  ImposterMaterial,
  InstancedMaterial,
  LinearFilter,
  LinearMipMapLinearFilter,
  MaterialIdCount,
  MeshOutlineMaterial,
  MixOperation,
  MultiplyOperation,
  NearestFilter,
  NoColors,
  NoShading,
  R32Format,
  RFormat,
  RGBAFormat,
  SmoothShading,
  SphereImposterOutlineMaterial,
  SpriteAlignment,
  SpriteMaterial,
  StickImposterMaterial,
  StickImposterOutlineMaterial,
  Texture,
  TextureIdCount,
  UVMapping,
  UnsignedByteType,
  VertexColors,
  VolumetricMaterial,
} from './WebGL/materials';

import {Cylinder, Sphere, Triangle} from './WebGL/shapes';
import {SurfaceWorker} from './SurfaceWorker';
import Renderer from './WebGL/Renderer';
import {ShaderUtils, ShaderLib} from './WebGL/shaders';
import SpritePlugin from './WebGL/SpritePlugin';
import {
  Matrix3,
  Matrix4,
  Quaternion,
  Ray,
  Vector2,
  Vector3,
  clamp,
  conversionMatrix3,
  degToRad,
  square,
} from './WebGL/math';
import EventDispatcher from './WebGL/core/EventDispatcher';
import createViewerGrid from "./util/createViewerGrid";
import createViewer from "./util/createViewer";
import Viewers from "./singletons/Viewers";
import download from "./util/download";
import MMTF from "./MMTF";
import getPropertyRange from "./util/getPropertyRange";
import SurfaceType from "./enum/SurfaceType";
import getbin from "./util/getbin";
import CAP from "./enum/CAP";
import createStereoViewer from "./util/createStereoViewer";
import extend from "./util/extend";

globalThis.$ = $;
globalThis.MMTF = MMTF;
export default {
  Color,
  GLViewer,
  GLModel,
  GLShape,
  Gradient,
  GLDraw,
  Label,
  CC,
  elementColors,
  drawCartoon,
  Parsers,
  MarchingCube,
  VolumeData,
  GLVolumetricRender,
  htmlColors,
  Camera,
  Fog,
  Object3D,
  Geometry,
  Material,
  LineBasicMaterial,
  MeshLambertMaterial,
  MeshDoubleLambertMaterial,
  SphereImposterMaterial,
  AddOperation,
  BackSide,
  ClampToEdgeWrapping,
  DoubleSide,
  FaceColors,
  FlatShading,
  FloatType,
  FrontSide,
  ImposterMaterial,
  InstancedMaterial,
  LinearFilter,
  LinearMipMapLinearFilter,
  MaterialIdCount,
  MeshOutlineMaterial,
  MixOperation,
  MultiplyOperation,
  NearestFilter,
  NoColors,
  NoShading,
  R32Format,
  RFormat,
  RGBAFormat,
  SmoothShading,
  SphereImposterOutlineMaterial,
  SpriteAlignment,
  SpriteMaterial,
  StickImposterMaterial,
  StickImposterOutlineMaterial,
  Texture,
  TextureIdCount,
  UVMapping,
  UnsignedByteType,
  VertexColors,
  VolumetricMaterial,
  GeometryIDCount,
  Light,
  Object3DIDCount,
  Projection,
  Raycaster,
  Scene,
  geometryGroup,
  MarchingCubeInitializer,
  LabelCount,
  ProteinSurface,
  Cylinder,
  Sphere,
  Triangle,
  SurfaceWorker,
  Renderer,
  ShaderUtils,
  ShaderLib,
  SpritePlugin,
  Matrix3,
  Matrix4,
  Quaternion,
  Ray,
  Vector2,
  Vector3,
  Math: {
    clamp,
    conversionMatrix3,
    degToRad,
    square,
  },
  square,
  conversionMatrix3,
  EventDispatcher,
  createViewer,
  createViewerGrid,
  viewers: Viewers,
  download,
  getPropertyRange,
  SurfaceType,
  getbin,
  ssColors,
  CAP,
  createStereoViewer,
  extend
};

