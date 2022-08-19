/**
 * GLViewer input specification
 */
export type ViewerSpec = Partial<{
  /** Callback function to be immediately executed on this viewer */
  callback: (viewer: ViewerSpec) => void;
  /** Object defining default atom colors as atom => color property value pairs for all models within this viewer */
  defaultcolors: Record<string, string>;
  /** 
   * Whether to disable disable handling of mouse events.
   * If you want to use your own mouse handlers, set this then bind your handlers to the canvas object.  
                The default 3Dmol.js handlers are available for use: 
                'mousedown touchstart': viewer._handleMouseDown,
                'DOMMouseScroll mousewheel': viewer._handleMouseScroll
                'mousemove touchmove': viewer._handleMouseMove   
   */
  nomouse: boolean | string;
  /** Color of the canvas background */
  backgroundColor: string;
  /** Alpha transparency of canvas background */
  backgroundAlpha: number;
  /** */
  camerax: number;
  /** */
  hoverDuration: number;
  /** id of the canvas */
  id: string;
  /** default 5 */
  cartoonQuality: number;
  /** */
  row: number;
  /** */
  col: number;
  /** */
  rows: number;
  /** */
  cols: number;
  /** */
  canvas: HTMLCanvasElement;
  /** GLViewer has not been ported but this is a viewer object */
  viewers: unknown[];
  /** */
  minimumZoomToDistance: number;
  /** */
  lowerZoomLimit: number;
  /** */
  upperZoomLimit: number;
  /** */
  antialias: boolean;
  /** */
  control_all: boolean;
  /** */
  orthographic: boolean;
  /** Disable fog, default to false */
  disableFog: boolean;

}>

/**
 * Grid GLViewer input specification
 */
export type ViewerGridSpec = Partial<{
  /** number of rows in grid */
  rows: number;
  /** number of columns in grid */
  cols: number;
  /** if true, mouse events are linked */
  control_all: boolean;
}>;