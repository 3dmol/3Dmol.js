
var WebMol = {};

/**
 * 
 * @param {Object} element
 * @param {ViewerSpec} config
 */
WebMol.createViewer = function(element, config) {};
/**
 * 
 * @param {Object} query
 * @param {Object} viewer
 */
WebMol.download = function(query, viewer) {};

//Viewer
// The constructor
/**
 * WebGL WebMol viewer
 * Note: The preferred method of instantiating a GLViewer is through {@link WebMol.createViewer} 
 * 
 * @constructor 
 * @param {Object} element HTML element within which to create viewer
 * @param {function} callback - Callback function to be immediately executed on this viewer
 * @param {Object} defaultcolors - Object defining default atom colors as atom => color property value pairs for all models within this viewer
 */
WebMol.GLViewer = function(element, callback, defaultcolors) {};



WebMol.GLModel = {};

WebMol.JmolElementColors = {};
WebMol.rasmolElementColors = {};


//Specification arguments

//Viewer spec - for WebMol.createViewer

var ViewerSpec = {};
ViewerSpec.order = [];
ViewerSpec.defaultcolors = {};

/**
 * Input specification argument for WebMol.createViewer
 * @param {WebMol.GLViewer} viewer
 */
ViewerSpec.callback = function(viewer) {};
