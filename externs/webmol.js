/** @namespace */
var WebMol = {};

/** Object with x, y, and z properties
 * @typedef {{x:number, y:number, z:number}} */
var vectorlike;

/**
 * 
 * @param {string | Object} element
 * @param {ViewerSpec} config
 */
WebMol.createViewer = function(element, config) {};
/**
 * @param {Object} query
 * @param {WebMol.GLViewer} viewer
 */
WebMol.download = function(query, viewer) {};

/**
 * @constructor
 * @param {...number} color
 */
WebMol.Color = function(color) {};
WebMol.Color.r;
WebMol.Color.g;
WebMol.Color.b;

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
//TODO: flesh out the annotations

/**
 * 
 * @constructor
 * @struct
 */
var ViewerSpec = {};
ViewerSpec.order;
ViewerSpec.defaultcolors;
/** @type {function(WebMol.GLViewer)} */
ViewerSpec.callback;


