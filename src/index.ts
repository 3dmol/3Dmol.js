// Create the global $3Dmol "namespace"


export * from "./Gradient";
export * from "./colors";
export * from "./Label";
export * from "./partialCharges";
export * from "./parsers";
export * from "./WebGL/math";
export * from "./WebGL/shapes";
export * from "./utilities";
export * from "./ProteinSurface4"
export * from "./VolumeData"
export * from "./VolumetricRender"
export * from "./GLShape"
export * from "./GLDraw"
export * from "./glcartoon"
export * from "./GLModel"
export * from "./GLViewer"
export * from "./autoload"
export * from "./specs"

//import * as $ from 'jquery';
declare var __webpack_exports__: any;

if (window) {
    //this needs to be exported here so the webworker can see it
    window.$3Dmol = __webpack_exports__;
}
