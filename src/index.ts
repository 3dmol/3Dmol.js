// Create the global $3Dmol "namespace"


export * from "./3Dmol";

//import * as $ from 'jquery';
declare var __webpack_exports__: any;

if (window) {
    //this needs to be exported here so the webworker can see it
    window.$3Dmol = __webpack_exports__;
}
