// This defines the $3Dmol object which is used to create viewers
// and configure system-wide settings

/**
 * All of the functionality of $3Dmol.js is contained within the
 * $3Dmol global namespace
 * @namespace */
// eslint-disable-next-line no-var
var $3Dmol = (function threeDeeMollContructor(window) {
  const my = window.$3Dmol || {};

  return my;
})(window);

if (typeof module == 'object' && typeof module.exports == 'object') {
  // for node.js exporting
  module.exports = $3Dmol;
}

/* The following code "phones home" to register that an ip 
   address has loaded 3Dmol.js.  Being able track this usage
   is very helpful when reporting to funding agencies.  Please
   leave this code in if you would like to increase the 
   likelihood of 3Dmol.js remaining supported.
*/
if (!$3Dmol.notrack) {
  /* The https traffic is just too much for the current server
   * to handle, so disable for now */
  // $.get("https://3dmol.csb.pitt.edu/track/report.cgi");
}

/* shims for IE */
/*
IE Doesn't have a .startsWith 
*/
if (!String.prototype.startsWith) {
  String.prototype.startsWith = str => this.lastIndexOf(str, 0) === 0;
}

// or endsWith
if (!String.prototype.endsWith) {
  String.prototype.endsWith = suffix => this.indexOf(suffix, this.length - suffix.length) !== -1;
}

/**
 *
 * jquery.binarytransport.js
 *
 * @description. jQuery ajax transport for making binary data type requests.
 * @version 1.0
 * @author Henry Algus <henryalgus@gmail.com>
 *
 */

// use this transport for "binary" data type
$.ajaxTransport('+binary', (options, originalOptions, jqXHR) => {
  // check for conditions and support for blob / arraybuffer response type
  if (
    window.FormData &&
    ((options.dataType && options.dataType === 'binary') ||
      (options.data &&
        ((window.ArrayBuffer && options.data instanceof ArrayBuffer) ||
          (window.Blob && options.data instanceof Blob))))
  ) {
    return {
      // create new XMLHttpRequest
      send(headers, callback) {
        // setup all variables
        const xhr = new XMLHttpRequest();
        const {url} = options;
        const {type} = options;
        const async = options.async || true;
        // blob or arraybuffer. Default is blob
        const dataType = options.responseType || 'blob';
        let data = options.data || null;
        const username = options.username || null;
        const password = options.password || null;

        const xhrret = () => {
          data = {};
          data[options.dataType] = xhr.response;
          // make callback and send data
          callback(xhr.status, xhr.statusText, data, xhr.getAllResponseHeaders());
        };

        xhr.addEventListener('load', xhrret);
        xhr.addEventListener('error', xhrret);
        xhr.addEventListener('abort', xhrret);

        xhr.open(type, url, async, username, password);

        // setup custom headers
        for (const i in headers) {
          xhr.setRequestHeader(i, headers[i]);
        }

        xhr.responseType = dataType;
        xhr.send(data);
      },
      abort() {
        jqXHR.abort();
      },
    };
  }
  return null;
});

/**
 * Create and initialize an appropriate viewer at supplied HTML element using specification in config
 @function $3Dmol.createViewer
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {ViewerSpec} config Viewer configuration
 * @param {Object} sharedViewerResources shared resources between viewers' renderers
 * @return {$3Dmol.GLViewer} GLViewer, null if unable to instantiate WebGL
 * @example
   var viewer = $3Dmol.createViewer(
     'gldiv', //id of div to create canvas in
     {
       defaultcolors: $3Dmol.elementColors.rasmol,
       backgroundColor: 'black'
     }
   );
 *                        
 */

$3Dmol.createViewer = function createViewer(element, config, sharedViewerResources) {
  const el = typeof element == 'string' ? $(`#${element}`) : element;
  if (!el) return null;

  const cfg = config || {};
  const svr = sharedViewerResources || {};

  // try to create the  viewer
  try {
    const viewer = new $3Dmol.GLViewer(el, cfg, svr);
    return viewer;
  } catch (e) {
    // throw `error creating viewer: ${e}`;
    throw new Error('error creating viewer');
  }

  // return null;
};

/**
 * Create and initialize an appropriate a grid of viewers that share a WebGL canvas
 @function $3Dmol.createViewerGrid
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {GridSpec} config - grid configuration
 * @param {ViewerGridSpec} viewerConfig - Viewer specification to apply to all subviewers
 * @return [[$3Dmol.GLViewer]] 2D array of GLViewers
 * @example                    
   var viewers = $3Dmol.createViewerGrid(
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
$3Dmol.createViewerGrid = function createViewerGrid(element, config, viewerConfig) {
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
    for (let r = 0; r < cfg.rows; r++) {
      const row = [];
      for (let c = 0; c < cfg.cols; c++) {
        vCfg.row = r;
        vCfg.col = c;
        vCfg.canvas = canvas;
        vCfg.viewers = viewers;
        vCfg.controlAll = cfg.controlAll;
        const viewer = $3Dmol.createViewer(el, vCfg);
        row.push(viewer);
      }
      viewers.unshift(row); // compensate for weird ordering in renderer
    }
  } catch (e) {
    // throw `error creating viewer grid: ${e}`;
    throw new Error('error creating viewer grid');
  }

  return viewers;
};

/**
 * Contains a dictionary of embedded viewers created from HTML elements
 * with a the viewer_3Dmoljs css class indexed by their id (or numerically
 * if they do not have an id).
 */
$3Dmol.viewers = {};

/**
 * Download binary data (e.g. a gzipped file) into an array buffer and provide
 * arraybuffer to callback.
 * @function $3Dmol.getbin
 * @param {string} uri - location of data
 * @param {Function} callback - Function to call with arraybuffer as argument.
 * @param {string} request - type of request
 * @return {Promise}
 */
$3Dmol.getbin = function getbin(uri, callback, request, postdata) {
  const promise = new Promise((resolve, reject) => {
    $.ajax({
      url: uri,
      dataType: 'binary',
      method: request || 'GET',
      data: postdata,
      responseType: 'arraybuffer',
      processData: false,
    })
      .done(ret => {
        resolve(ret);
      })
      .fail((e, txt) => {
        console.log(txt);
        reject();
      });
  });
  if (callback) return promise.then(callback);
  return promise;
};

/**
 * Convert a base64 encoded string to a Uint8Array
 * @function $3Dmol.base64ToArray
 * @param {string} base64 encoded string
 */
$3Dmol.base64ToArray = function base64ToArray(base64) {
  const binaryString = window.atob(base64);
  const len = binaryString.length;
  const bytes = new Uint8Array(len);
  for (let i = 0; i < len; i++) {
    bytes[i] = binaryString.charCodeAt(i);
  }
  return bytes;
};

/**
 * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
 * @function $3Dmol.download
 * @param {string} query - String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
 * @param {$3Dmol.GLViewer} viewer - Add new model to existing viewer
 * @param {Object} options - Specify additional options
 *                           format: file format to download, if multiple are available, default format is pdb
 *                           pdbUri: URI to retrieve PDB files, default URI is http://www.rcsb.org/pdb/files/
 * @param {Function} callback - Function to call with model as argument after data is loaded.
  
 * @return {$3Dmol.GLModel} GLModel, Promise if callback is not provided
 * @example
 viewer.setBackgroundColor(0xffffffff);
       $3Dmol.download('pdb:2nbd',viewer,{onemol: true,multimodel: true},function(m) {
        m.setStyle({'cartoon':{colorscheme:{prop:'ss',map:$3Dmol.ssColors.Jmol}}});
       viewer.zoomTo();
       viewer.render(callback);
    });
 */
$3Dmol.download = function download(query, viewer, options, callback) {
  let type = '';
  let pdbUri = '';
  let mmtfUri = '';
  let uri = '';
  let promise = null;
  const m = viewer.addModel();

  if (query.indexOf(':') < 0) {
    // no type specifier, guess
    if (query.length === 4) {
      query = `pdb:${query}`;
    } else if (!Number.isNaN(query)) {
      query = `cid:${query}`;
    } else {
      query = `url:${query}`;
    }
  }
  if (query.substring(0, 5) === 'mmtf:') {
    pdbUri = options && options.pdbUri ? options.pdbUri : 'https://mmtf.rcsb.org/v1.0/full/';
    query = query.substring(5).toUpperCase();
    uri = pdbUri + query;
    if (options && typeof options.noComputeSecondaryStructure == 'undefined') {
      // when fetch directly from pdb, trust structure annotations
      options.noComputeSecondaryStructure = true;
    }
    promise = new Promise(resolve => {
      $3Dmol.getbin(uri).then(
        ret => {
          m.addMolData(ret, 'mmtf', options);
          viewer.zoomTo();
          viewer.render();
          resolve(m);
        },
        () => {
          console.log(`fetch of ${uri} failed.`);
        }
      );
    });
  } else {
    if (query.substr(0, 4) === 'pdb:') {
      type = 'mmtf';
      if (options && options.format) {
        type = options.format; // can override and require pdb
      }

      if (options && typeof options.noComputeSecondaryStructure == 'undefined') {
        // when fetch directly from pdb, trust structure annotations
        options.noComputeSecondaryStructure = true;
      }
      query = query.substr(4).toUpperCase();
      if (!query.match(/^[1-9][A-Za-z0-9]{3}$/)) {
        // eslint-disable-next-line no-alert
        alert('Wrong PDB ID');
        return null;
      }
      if (type === 'mmtf') {
        mmtfUri = options && options.mmtfUri ? options.mmtfUri : 'https://mmtf.rcsb.org/v1.0/full/';
        uri = mmtfUri + query.toUpperCase();
      } else {
        pdbUri = options && options.pdbUri ? options.pdbUri : 'https://files.rcsb.org/view/';
        uri = `${pdbUri + query}.${type}`;
      }
    } else if (query.substring(0, 4) === 'cid:') {
      type = 'sdf';
      query = query.substring(4);
      if (!query.match(/^[0-9]+$/)) {
        // eslint-disable-next-line no-alert
        alert('Wrong Compound ID');
        return null;
      }
      uri = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${query}/SDF?record_type=3d`;
    } else if (query.substring(0, 4) === 'url:') {
      uri = query.substring(4);
      type = uri;
    }

    const handler = ret => {
      m.addMolData(ret, type, options);
      viewer.zoomTo();
      viewer.render();
    };
    promise = new Promise(resolve => {
      if (type === 'mmtf') {
        // binary data
        $3Dmol
          .getbin(uri)
          .then(ret => {
            handler(ret);
            resolve(m);
          })
          .catch(() => {
            // if mmtf server is being annoying, fallback to text
            pdbUri = options && options.pdbUri ? options.pdbUri : 'https://files.rcsb.org/view/';
            uri = `${pdbUri + query}.pdb`;
            console.log('falling back to pdb format');
            $.get(uri, ret => {
              handler(ret);
              resolve(m);
            }).fail(e => {
              handler('');
              resolve(m);
              console.log(`fetch of ${uri} failed: ${e.statusText}`);
            });
          }); // an error msg has already been printed
      } else {
        $.get(uri, ret => {
          handler(ret);
          resolve(m);
        }).fail(e => {
          handler('');
          resolve(m);
          console.log(`fetch of ${uri} failed: ${e.statusText}`);
        });
      }
    });
  }
  if (callback) {
    promise.then(callback);
    return m;
  }
  return promise;
};

/**
 * $3Dmol surface types
 * @enum {number}
 */
$3Dmol.SurfaceType = {
  VDW: 1,
  MS: 2,
  SAS: 3,
  SES: 4,
};

// Miscellaneous functions and classes - to be incorporated into $3Dmol proper
/**
 *
 * @param {$3Dmol.Geometry} geometry
 * @param {$3Dmol.Mesh} mesh
 * @returns {undefined}
 */
$3Dmol.mergeGeos = function (geometry, mesh) {
  const meshGeo = mesh.geometry;

  if (meshGeo === undefined) return;

  geometry.geometryGroups.push(meshGeo.geometryGroups[0]);
};

$3Dmol.multiLineString = function (f) {
  return f
    .toString()
    .replace(/^[^\/]+\/\*!?/, '')
    .replace(/\*\/[^\/]+$/, '');
};

/**
 * Render surface synchronously if true
 * @param {boolean} [$3Dmol.SyncSurface=false]
 * @type {boolean} */
$3Dmol.syncSurface = false;

// Internet Explorer refuses to allow webworkers in data blobs.  I can find
// no way of checking for this feature directly, so must do a browser check
if (
  window.navigator.userAgent.indexOf('MSIE ') >= 0 ||
  window.navigator.userAgent.indexOf('Trident/') >= 0
) {
  $3Dmol.syncSurface = true; // can't use webworkers
}

/**
 * Parse a string that represents a style or atom selection and convert it
 * into an object.  The goal is to make it easier to write out these specifications
 * without resorting to json. Objects cannot be defined recursively.
 * ; - delineates fields of the object
 * : - if the field has a value other than an empty object, it comes after a colon
 * , - delineates key/value pairs of a value object
 *     If the value object consists of ONLY keys (no = present) the keys are
 *     converted to a list.  Otherwise a object of key/value pairs is created with
 *     any missing values set to null
 * = OR ~ - separates key/value pairs of a value object, if not provided value is null
 *     twiddle is supported since = has special meaning in URLs
 * @param (String) str
 * @returns {Object}
 */
$3Dmol.specStringToObject = function (str) {
  if (typeof str == 'object') {
    return str; // not string, assume was converted already
  }
  if (typeof str == 'undefined' || str == null) {
    return str;
  }

  str = str.replace(/%7E/, '~'); // copy/pasting urls sometimes does this
  // convert things that look like numbers into numbers
  const massage = val => {
    if ($3Dmol.isNumeric(val)) {
      // hexadecimal does not parse as float
      if (Math.floor(parseFloat(val, 16)) === parseInt(val, 16)) {
        return parseFloat(val);
      }
      if (val.indexOf('.') >= 0) {
        return parseFloat(val); // ".7" for example, does not parseInt
      }

      return parseInt(val, 10);
    }
    // boolean conversions
    if (val === 'true') {
      return true;
    }
    if (val === 'false') {
      return false;
    }
    return val;
  };

  const ret = {};
  if (str === 'all') return ret;
  const fields = str.split(';');
  for (let i = 0; i < fields.length; i++) {
    const fv = fields[i].split(':');
    const f = fv[0];
    let val = {};
    let vstr = fv[1];
    if (vstr) {
      vstr = vstr.replace(/~/g, '=');
      if (vstr.indexOf('=') !== -1) {
        // has key=value pairs, must be object
        const kvs = vstr.split(',');
        for (let j = 0; j < kvs.length; j++) {
          const kv = kvs[j].split('=', 2);
          val[kv[0]] = massage(kv[1]);
        }
      } else if (vstr.indexOf(',') !== -1) {
        // has multiple values, must list
        val = vstr.split(',');
      } else {
        val = massage(vstr); // value itself
      }
    }
    ret[f] = val;
  }

  return ret;
};

/**
 * computes the bounding box around the provided atoms
 * @param {AtomSpec[]} atomlist
 * @return {Array}
 */
$3Dmol.getExtent = function (atomlist, ignoreSymmetries) {
  let xmin;
  let ymin;
  let zmin;
  let xmax;
  let ymax;
  let zmax;
  let xsum;
  let ysum;
  let zsum;
  let cnt;
  const includeSym = !ignoreSymmetries;

  xmin = ymin = zmin = 9999;
  xmax = ymax = zmax = -9999;
  xsum = ysum = zsum = cnt = 0;

  if (atomlist.length === 0)
    return [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
  for (let i = 0; i < atomlist.length; i++) {
    const atom = atomlist[i];
    if (
      typeof atom == 'undefined' ||
      !Number.isFinite(atom.x) ||
      !Number.isFinite(atom.y) ||
      !Number.isFinite(atom.z)
    )
      continue;
    cnt += 1;
    xsum += atom.x;
    ysum += atom.y;
    zsum += atom.z;

    xmin = xmin < atom.x ? xmin : atom.x;
    ymin = ymin < atom.y ? ymin : atom.y;
    zmin = zmin < atom.z ? zmin : atom.z;
    xmax = xmax > atom.x ? xmax : atom.x;
    ymax = ymax > atom.y ? ymax : atom.y;
    zmax = zmax > atom.z ? zmax : atom.z;

    if (atom.symmetries && includeSym) {
      for (let n = 0; n < atom.symmetries.length; n++) {
        cnt += 1;
        xsum += atom.symmetries[n].x;
        ysum += atom.symmetries[n].y;
        zsum += atom.symmetries[n].z;
        xmin = xmin < atom.symmetries[n].x ? xmin : atom.symmetries[n].x;
        ymin = ymin < atom.symmetries[n].y ? ymin : atom.symmetries[n].y;
        zmin = zmin < atom.symmetries[n].z ? zmin : atom.symmetries[n].z;
        xmax = xmax > atom.symmetries[n].x ? xmax : atom.symmetries[n].x;
        ymax = ymax > atom.symmetries[n].y ? ymax : atom.symmetries[n].y;
        zmax = zmax > atom.symmetries[n].z ? zmax : atom.symmetries[n].z;
      }
    }
  }

  return [
    [xmin, ymin, zmin],
    [xmax, ymax, zmax],
    [xsum / cnt, ysum / cnt, zsum / cnt],
  ];
};

// return the value of an atom property prop, or null if non existent
// looks first in properties, then in the atom itself
$3Dmol.getAtomProperty = function (atom, prop) {
  let val = null;
  if (atom.properties && typeof atom.properties[prop] != 'undefined') {
    val = atom.properties[prop];
  } else if (typeof atom[prop] != 'undefined') {
    val = atom[prop];
  }
  return val;
};

/* get the min and max values of the specified property in the provided
 * @function $3Dmol.getPropertyRange
 * @param {AtomSpec[]} atomlist - list of atoms to evaluate
 * @param {string} prop - name of property
 * @return {Array} - [min, max] values
 */
$3Dmol.getPropertyRange = function (atomlist, prop) {
  let min = Number.POSITIVE_INFINITY;
  let max = Number.NEGATIVE_INFINITY;

  for (let i = 0, n = atomlist.length; i < n; i++) {
    const atom = atomlist[i];
    const val = $3Dmol.getAtomProperty(atom, prop);

    if (val != null) {
      if (val < min) min = val;
      if (val > max) max = val;
    }
  }

  if (!Number.isFinite(min) && !Number.isFinite(max)) min = max = 0;
  else if (!Number.isFinite(min)) min = max;
  else if (!Number.isFinite(max)) max = min;

  return [min, max];
};

// hackish way to work with requirejs - doesn't actually work yet
// since we don't use the require optimizer to combine modules
if (typeof _3dmolSavedDefine != 'undefined') {
  /** When pulling in external sources, disable amd to ensure they
   * populate the global namespace as expected.  Restore it so code
   * using amd still works. */
  /* global _3dmolSavedDefine, _3dmolSavedRequire, define:true, require:true */
  define = _3dmolSavedDefine;
  require = _3dmolSavedRequire;
}
if (typeof define == 'function' && define.amd) {
  define('$3Dmol', [], () => $3Dmol);
}

/* StereoViewer for stereoscopic viewing
  @function $3Dmol.createStereoViewer
* @param {Object | string} element - Either HTML element or string identifier
* 
*/

$3Dmol.createStereoViewer = function (element) {
  const that = this;
  element = typeof element == 'string' ? $(`#${element}`) : element;
  if (!element) return;

  const viewers = $3Dmol.createViewerGrid(element, {rows: 1, cols: 2, controlAll: true});

  this.glviewer1 = viewers[0][0];
  this.glviewer2 = viewers[0][1];

  this.glviewer1.setAutoEyeSeparation(false);
  this.glviewer2.setAutoEyeSeparation(true);

  this.glviewer1.linkViewer(this.glviewer2);
  this.glviewer2.linkViewer(this.glviewer1);

  const methods = Object.getOwnPropertyNames(this.glviewer1) // get all methods of glviewer object
    .filter(property => typeof that.glviewer1[property] == 'function');

  for (let i = 0; i < methods.length; i++) {
    // create methods of the same name
    this[methods[i]] = (function (method) {
      return function () {
        const m1 = this.glviewer1[method].apply(this.glviewer1, arguments);
        const m2 = this.glviewer2[method].apply(this.glviewer2, arguments);
        return [m1, m2];
      };
    })(methods[i]);
  }

  // special cased methods
  this.setCoordinates = (models, data, format) => {
    // for setting the coordinates of the models
    for (let i = 0; i < models.length; i++) {
      models[i].setCoordinates(data, format);
    }
  };

  this.surfacesFinished = () =>
    this.glviewer1.surfacesFinished() && this.glviewer2.surfacesFinished();

  this.isAnimated = () => this.glviewer1.isAnimated() || this.glviewer2.isAnimated();

  this.render = callback => {
    this.glviewer1.render();
    this.glviewer2.render();
    if (callback) {
      callback(this); // call only once
    }
  };

  this.getCanvas = () => this.glviewer1.getCanvas(); // same for both
};

// simplified version of $.extend
$3Dmol.extend = (obj1, src1) => {
  for (const key in src1) {
    if (src1.hasOwnProperty(key) && src1[key] !== undefined) {
      // if(Object.prototype.hasOwnProperty.call(src1,key) && src1[key] !== undefined){ // use Object.prototype
      obj1[key] = src1[key];
    }
  }
  return obj1;
};

// deep copy, cannot deal with circular refs; undefined input becomes an empty object
// https://medium.com/javascript-in-plain-english/how-to-deep-copy-objects-and-arrays-in-javascript-7c911359b089
$3Dmol.deepCopy = inObject => {
  let outObject;
  let value;
  let key;

  if (inObject === undefined) {
    return {};
  }
  if (typeof inObject != 'object' || inObject == null) {
    return inObject; // Return the value if inObject is not an object
  }

  // Create an array or object to hold the values
  // eslint-disable-next-line prefer-const
  outObject = Array.isArray(inObject) ? [] : {};

  for (key in inObject) {
    value = inObject[key];
    // Recursively (deep) copy for nested objects, including arrays
    outObject[key] = $3Dmol.deepCopy(value);
  }

  return outObject;
};

$3Dmol.isNumeric = obj => {
  const type = typeof obj;
  return (type === 'number' || type === 'string') && !isNaN(obj - parseFloat(obj));
};

$3Dmol.isEmptyObject = obj => {
  let name;
  for (name in obj) {
    return false;
  }
  return true;
};

$3Dmol.makeFunction = callback => {
  // for py3dmol let users provide callback as string
  if (callback && typeof callback == 'string') {
    // eslint-disable-next-line no-eval
    callback = eval(`(${callback})`);
  }
  // report to console if callback is not a valid function
  if (callback && typeof callback != 'function') {
    return null;
  }
  return callback;
};

// standardize voldata/volscheme in style
$3Dmol.adjustVolumeStyle = style => {
  if (style) {
    if (style.volformat && !(style.voldata instanceof $3Dmol.VolumeData)) {
      style.voldata = new $3Dmol.VolumeData(style.voldata, style.volformat);
    }
    if (style.volscheme) {
      style.volscheme = $3Dmol.Gradient.getGradient(style.volscheme);
    }
  }
};
