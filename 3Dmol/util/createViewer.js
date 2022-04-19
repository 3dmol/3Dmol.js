import $ from "jquery";
import GLViewer from "../GLViewer";


/**
 * Create and initialize an appropriate viewer at supplied HTML element using specification in config
 @function createViewer
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {import("../specs").ViewerSpec} config Viewer configuration
 * @param {Object} [sharedViewerResources] shared resources between viewers' renderers
 * @return {GLViewer|null} GLViewer, null if unable to instantiate WebGL
 * @example
   var viewer = createViewer(
     'gldiv', //id of div to create canvas in
     {
       defaultcolors: elementColors.rasmol,
       backgroundColor: 'black'
     }
   );
 *                        
 */

export default function createViewer(element, config, sharedViewerResources) {
  const el = typeof element == 'string' ? $(`#${element}`) : element;
  if (!el) return null;

  const cfg = config || {};
  const svr = sharedViewerResources || {};

  // try to create the  viewer
  try {
    const viewer = new GLViewer(el, cfg, svr);
    return viewer;
  } catch (e) {
    // throw `error creating viewer: ${e}`;
    throw new Error('error creating viewer');
  }

  // return null;
}
