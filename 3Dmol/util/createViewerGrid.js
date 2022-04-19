import $ from "jquery";
import createViewer from "./createViewer";
/**
 * Create and initialize an appropriate a grid of viewers that share a WebGL canvas
 @function createViewerGrid
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {import("../specs").GridSpec} config - grid configuration
 * @param {import("../specs").ViewerGridSpec} viewerConfig - Viewer specification to apply to all subviewers
 * @return [[GLViewer]] 2D array of GLViewers
 * @example                    
   var viewers = createViewerGrid(
     'gldiv', //id of div to create canvas in
     {
       rows: 2,
       cols: 2,
       controlAll: true  //mouse controls all viewers
     },
     { backgroundColor: 'lightgrey' }
   );
   $.get('data/1jpy.cif', function(data) {
     var viewer = viewers[0][0];
     viewer.addModel(data,'cif');
     viewer.setStyle({sphere:{}});
     viewer.zoomTo();
     viewer.render( );

     viewer = viewers[0][1];
     viewer.addModel(data,'cif');
     viewer.setStyle({stick:{}});
     viewer.zoomTo();     
     viewer.render( );
     
     viewer = viewers[1][0];
     viewer.addModel(data,'cif');
     viewer.setStyle({cartoon:{color:'spectrum'}});
     viewer.zoomTo();     
     viewer.render( );
     
     viewer = viewers[1][1];
     viewer.addModel(data,'cif');
     viewer.setStyle({cartoon:{colorscheme:'chain'}});
     viewer.zoomTo();     
     viewer.render();
     
     
   });
     
 */
export default function createViewerGrid(element, config, viewerConfig) {
  const el = typeof element == 'string' ? $(`#${element}`) : element;
  if (!el) return null;

  const cfg = config || {};
  const vCfg = viewerConfig || {};

  const viewers = [];
  // create canvas
  const canvas = document.createElement('canvas');

  vCfg.rows = cfg.rows;
  vCfg.cols = cfg.cols;
  vCfg.controlAll = cfg.controlAll !== undefined ? cfg.controlAll : false;
  $(el).append($(canvas));

  // try to create the  viewer
  try {
    for (let r = 0; r < (cfg.rows||0); r++) {
      const row = [];
      for (let c = 0; c < (cfg.cols||0); c++) {
        vCfg.row = r;
        vCfg.col = c;
        vCfg.canvas = canvas;
        vCfg.viewers = viewers;
        vCfg.controlAll = cfg.controlAll;
        // @ts-ignore
        const viewer = createViewer(el, vCfg);
        row.push(viewer);
      }
      viewers.unshift(row); // compensate for weird ordering in renderer
    }
  } catch (e) {
    // throw `error creating viewer grid: ${e}`;
    throw new Error('error creating viewer grid');
  }

  return viewers;
}
