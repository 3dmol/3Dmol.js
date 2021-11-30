
//This defines the $3Dmol object which is used to create viewers
//and configure system-wide settings

/** 
 * All of the functionality of $3Dmol.js is contained within the
 * $3Dmol global namespace
 * @namespace */
$3Dmol = (function(window) {
    
    var my = window.$3Dmol || {};
    
    return my;

})(window);

if ( typeof module === "object" && typeof module.exports === "object" ) { 
	//for node.js exporting
	module.exports = $3Dmol; 
}

/* The following code "phones home" to register that an ip 
   address has loaded 3Dmol.js.  Being able track this usage
   is very helpful when reporting to funding agencies.  Please
   leave this code in if you would like to increase the 
   likelihood of 3Dmol.js remaining supported.
*/
if(!$3Dmol.notrack) {
/* The https traffic is just too much for the current server
 * to handle, so disable for now */
// $.get("https://3dmol.csb.pitt.edu/track/report.cgi");
}

/* shims for IE */
/*
IE Doesn't have a .startsWith 
*/
if (!String.prototype.startsWith) {
    String.prototype.startsWith = function (str){
        return this.lastIndexOf(str, 0) === 0;
    };
}

// or endsWith
if (!String.prototype.endsWith) {
    String.prototype.endsWith = function(suffix) {
        return this.indexOf(suffix, this.length - suffix.length) !== -1;
    };
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
$.ajaxTransport(
               "+binary",
               function(options, originalOptions, jqXHR) {
                   // check for conditions and support for blob / arraybuffer response type
                   if (window.FormData && ((options.dataType && (options.dataType == 'binary')) || 
                           (options.data && ((window.ArrayBuffer && options.data instanceof ArrayBuffer) || 
                                   (window.Blob && options.data instanceof Blob))))) {
                       return {
                           // create new XMLHttpRequest
                           send : function(headers, callback) {
                               // setup all variables
                               var xhr = new XMLHttpRequest(), url = options.url, type = options.type, async = options.async || true,
                               // blob or arraybuffer. Default is blob
                               dataType = options.responseType || "blob", 
                                           data = options.data || null, 
                                           username = options.username || null, 
                                           password = options.password || null;

                               var xhrret = function() {
                                   var data = {};
                                   data[options.dataType] = xhr.response;
                                   // make callback and send data
                                   callback(xhr.status, xhr.statusText,
                                           data,
                                           xhr.getAllResponseHeaders());
                               };
                               
                               xhr.addEventListener('load', xhrret);
                               xhr.addEventListener('error', xhrret);
                               xhr.addEventListener('abort', xhrret);
                               
                               xhr.open(type, url, async, username,
                                       password);

                               // setup custom headers
                               for ( var i in headers) {
                                   xhr.setRequestHeader(i, headers[i]);
                               }

                               xhr.responseType = dataType;
                               xhr.send(data);
                           },
                           abort : function() {
                               jqXHR.abort();
                           }
                       };
                   }
               });

    
/**
 * Create and initialize an appropriate viewer at supplied HTML element using specification in config
 @function $3Dmol.createViewer
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {ViewerSpec} config Viewer configuration
 * @param {Object} shared_viewer_resources shared resources between viewers' renderers
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

$3Dmol.createViewer = function(element, config, shared_viewer_resources)
{
    if(typeof(element) === "string")
    element = $("#"+element);
    if(!element) return;
    
    config = config || {}; 
    shared_viewer_resources = shared_viewer_resources || {};
    
    //try to create the  viewer
    try {
        var viewer = new $3Dmol.GLViewer(element, config, shared_viewer_resources);
        return viewer;
    }
    catch(e) {
        throw "error creating viewer: "+e;
    }
    
    return null;
};

/**
 * Create and initialize an appropriate a grid of viewers that share a WebGL canvas
 @function $3Dmol.createViewerGrid
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {GridSpec} config - grid configuration
 * @param {ViewerGridSpec} viewer_config - Viewer specification to apply to all subviewers
 * @return [[$3Dmol.GLViewer]] 2D array of GLViewers
 * @example                    
   var viewers = $3Dmol.createViewerGrid(
     'gldiv', //id of div to create canvas in
     {
       rows: 2,
       cols: 2,
       control_all: true  //mouse controls all viewers
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
$3Dmol.createViewerGrid  = function(element,config,viewer_config){
    if(typeof(element) === "string")
        element = $("#"+element);
    if(!element) return;

    config = config || {}; 
    viewer_config = viewer_config || {};
    
    var viewers = [];
    //create canvas
    var canvas = document.createElement('canvas');

    viewer_config.rows = config.rows;
    viewer_config.cols = config.cols;
    viewer_config.control_all = config.control_all != undefined ? config.control_all : false;
    $(element).append($(canvas));

      //try to create the  viewer
    try {  
      for(var r =0;r<config.rows;r++){
        var row = [];
        for(var c = 0;c<config.cols;c++){
          viewer_config.row = r;
          viewer_config.col = c;
          viewer_config.canvas = canvas;
          viewer_config.viewers = viewers;
          viewer_config.control_all = config.control_all;
          var viewer = $3Dmol.createViewer(element, viewer_config);
          row.push(viewer);
        }
        viewers.unshift(row); //compensate for weird ordering in renderer
      }
    }catch(e) {
        throw "error creating viewer grid: "+e;
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
$3Dmol.getbin = function(uri, callback, request,postdata) {
    var promise = new Promise(function(resolve, reject) {
        
        request = (request == undefined)?"GET":request;
        $.ajax({url:uri, 
            dataType: "binary",
            method: request,
            data: postdata,
            responseType: "arraybuffer",
            processData: false})
        .done(function(ret) {
            resolve(ret);
        })
        .fail(function(e,txt) { 
            console.log(txt);
            reject();
        });
    });
    if (callback) return promise.then(callback);
    else return promise;
};

/**
 * Convert a base64 encoded string to a Uint8Array
 * @function $3Dmol.base64ToArray
 * @param {string} base64 encoded string
 */
$3Dmol.base64ToArray = function(base64) {
    var binary_string =  window.atob(base64);
    var len = binary_string.length;
    var bytes = new Uint8Array( len );
    for (var i = 0; i < len; i++)        {
        bytes[i] = binary_string.charCodeAt(i);
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
$3Dmol.download = function(query, viewer, options, callback) {
    var type = "";
    var pdbUri = "";
    var mmtfUri = "";
    var uri = "";
    var promise = null;
    var m = viewer.addModel();
    if (query.substr(0, 5) === 'mmtf:') {
        pdbUri = options && options.pdbUri ? options.pdbUri : "https://mmtf.rcsb.org/v1.0/full/";
        query = query.substr(5).toUpperCase();
        uri = pdbUri + query;        
        if(options && typeof options.noComputeSecondaryStructure === 'undefined') {
                //when fetch directly from pdb, trust structure annotations
                options.noComputeSecondaryStructure = true;
        }
        promise = new Promise(function(resolve) {
            $3Dmol.getbin(uri)
            .then(function(ret) {
                m.addMolData(ret, 'mmtf',options);
                viewer.zoomTo();
                viewer.render();
                resolve(m);
            },function() {console.log("fetch of "+uri+" failed.");});
        });
    }
    else {
        if (query.substr(0, 4) === 'pdb:') {
            type = 'mmtf';
            if(options && options.format) {
                type = options.format; //can override and require pdb
            }
            
            if(options && typeof options.noComputeSecondaryStructure === 'undefined') {
                //when fetch directly from pdb, trust structure annotations
                options.noComputeSecondaryStructure = true;
            }
            query = query.substr(4).toUpperCase();
            if (!query.match(/^[1-9][A-Za-z0-9]{3}$/)) {
               alert("Wrong PDB ID"); return;
            }
            if(type == 'mmtf') {
                mmtfUri = options && options.mmtfUri ? options.mmtfUri : 'https://mmtf.rcsb.org/v1.0/full/';
                uri = mmtfUri + query.toUpperCase();
            }
            else  {
                pdbUri = options && options.pdbUri ? options.pdbUri : "https://files.rcsb.org/view/";
                uri = pdbUri + query + "." + type;
            }
    
        } else if (query.substr(0, 4) == 'cid:') {
            type = "sdf";
            query = query.substr(4);
            if (!query.match(/^[0-9]+$/)) {
               alert("Wrong Compound ID"); return;
            }
            uri = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + query + 
              "/SDF?record_type=3d";
        } else if (query.substr(0,4) == 'url:') {
            uri = query.substr(4);
            type = uri;
        }
    
        var handler = function(ret) {
            m.addMolData(ret, type, options);
            viewer.zoomTo();
            viewer.render();
        };
        promise = new Promise(function(resolve) {
            if(type == 'mmtf') { //binary data
                $3Dmol.getbin(uri)
                .then(function(ret) {
                    handler(ret);
                    resolve(m);
                }).catch(function(){
                    //if mmtf server is being annoying, fallback to text
                    pdbUri = options && options.pdbUri ? options.pdbUri : "https://files.rcsb.org/view/";
                    uri = pdbUri + query + ".pdb";
                    console.log("falling back to pdb format");
                    $.get(uri, function(ret) {
                        handler(ret);
                        resolve(m);
                    }).fail(function(e) {
                       handler("");
                       resolve(m);
                       console.log("fetch of "+uri+" failed: "+e.statusText);
                       });
                }); //an error msg has already been printed
            }
            else {        
               $.get(uri, function(ret) {
                   handler(ret);
                   resolve(m);
               }).fail(function(e) {
                   handler("");
                   resolve(m);
                   console.log("fetch of "+uri+" failed: "+e.statusText);
               });
            }
        });
    }
    if (callback) {
        promise.then(function(m){
            callback(m);
        });
        return m;
    }
    else return promise;
};
       

/**
 * $3Dmol surface types
 * @enum {number}
 */
$3Dmol.SurfaceType = {
    VDW : 1,
    MS : 2,
    SAS : 3,
    SES  : 4
};


//Miscellaneous functions and classes - to be incorporated into $3Dmol proper
/**
 * 
 * @param {$3Dmol.Geometry} geometry
 * @param {$3Dmol.Mesh} mesh
 * @returns {undefined}
 */
$3Dmol.mergeGeos = function(geometry, mesh) {
    
    var meshGeo = mesh.geometry;
    
    if (meshGeo === undefined) 
        return;
    
    geometry.geometryGroups.push( meshGeo.geometryGroups[0] );
    
};

$3Dmol.multiLineString = function(f) {
    return f.toString()
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
if(window.navigator.userAgent.indexOf('MSIE ') >= 0 ||
        window.navigator.userAgent.indexOf('Trident/') >= 0) {
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
$3Dmol.specStringToObject = function(str) {
    if(typeof(str) === "object") {
        return str; //not string, assume was converted already
    }
    else if(typeof(str) === "undefined" || str == null) {
        return str; 
    }
    
    str = str.replace(/%7E/,'~'); //copy/pasting urls sometimes does this
    //convert things that look like numbers into numbers
    var massage = function(val) {
        if($3Dmol.isNumeric(val)) {
           //hexadecimal does not parse as float
           if(Math.floor(parseFloat(val)) == parseInt(val)) {
              return parseFloat(val);
           }
           else if(val.indexOf('.') >= 0) {
               return parseFloat(val); // ".7" for example, does not parseInt
           }
           else{
               return parseInt(val);
           }
        }
        //boolean conversions
        else if(val === 'true') {
            return true;
        }
        else if(val === 'false') {
            return false;
        }
        return val;
    };
    
    var ret = {};
    if(str === 'all') return ret;
    var fields = str.split(';');
    for(var i = 0; i < fields.length; i++) {
        var fv = fields[i].split(':');
        var f = fv[0];
        var val = {};
        var vstr = fv[1];
        if(vstr) {
            vstr = vstr.replace(/~/g,"=");
            if(vstr.indexOf('=') !== -1) {
                //has key=value pairs, must be object
                var kvs = vstr.split(',');
                for(var j = 0; j < kvs.length; j++) {
                    var kv = kvs[j].split('=',2);
                    val[kv[0]] = massage(kv[1]);
                }
            }
            else if(vstr.indexOf(',') !== -1) {
                //has multiple values, must list
                val = vstr.split(',');
            }
            else {
                val = massage(vstr); //value itself
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
$3Dmol.getExtent = function(atomlist, ignoreSymmetries) {
    var xmin, ymin, zmin, xmax, ymax, zmax, xsum, ysum, zsum, cnt;
    var includeSym = !ignoreSymmetries;

    xmin = ymin = zmin = 9999;
    xmax = ymax = zmax = -9999;
    xsum = ysum = zsum = cnt = 0;
    
    if (atomlist.length === 0)
        return [ [ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ] ];
    for (var i = 0; i < atomlist.length; i++) {
        var atom = atomlist[i];
        if (typeof atom === 'undefined' || !isFinite(atom.x) ||
                !isFinite(atom.y) || !isFinite(atom.z))
            continue;
        cnt++;
        xsum += atom.x;
        ysum += atom.y;
        zsum += atom.z;
        
        xmin = (xmin < atom.x) ? xmin : atom.x;
        ymin = (ymin < atom.y) ? ymin : atom.y;
        zmin = (zmin < atom.z) ? zmin : atom.z;
        xmax = (xmax > atom.x) ? xmax : atom.x;
        ymax = (ymax > atom.y) ? ymax : atom.y;
        zmax = (zmax > atom.z) ? zmax : atom.z;
        
        if (atom.symmetries && includeSym) {
            for (var n = 0; n < atom.symmetries.length; n++) {
                cnt++;
                xsum += atom.symmetries[n].x;
                ysum += atom.symmetries[n].y;
                zsum += atom.symmetries[n].z;
                xmin = (xmin < atom.symmetries[n].x) ? xmin : atom.symmetries[n].x;
                ymin = (ymin < atom.symmetries[n].y) ? ymin : atom.symmetries[n].y;
                zmin = (zmin < atom.symmetries[n].z) ? zmin : atom.symmetries[n].z;
                xmax = (xmax > atom.symmetries[n].x) ? xmax : atom.symmetries[n].x;
                ymax = (ymax > atom.symmetries[n].y) ? ymax : atom.symmetries[n].y;
                zmax = (zmax > atom.symmetries[n].z) ? zmax : atom.symmetries[n].z; 
            }
        }  
    }

    return [ [ xmin, ymin, zmin ], [ xmax, ymax, zmax ],
            [ xsum / cnt, ysum / cnt, zsum / cnt ] ];
};


//return the value of an atom property prop, or null if non existent
// looks first in properties, then in the atom itself
$3Dmol.getAtomProperty = function(atom, prop) {
    var val = null;
    if (atom.properties &&
            typeof (atom.properties[prop]) != "undefined") {
        val = atom.properties[prop];
    } else if(typeof(atom[prop]) != 'undefined') {
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
    var min = Number.POSITIVE_INFINITY;
    var max = Number.NEGATIVE_INFINITY;

    for (var i = 0, n = atomlist.length; i < n; i++) {
        var atom = atomlist[i];
        var val = $3Dmol.getAtomProperty(atom, prop);
        
        if(val != null) {
            if (val < min)
                min = val;
            if (val > max)
                max = val;                
        }
    }

    if (!isFinite(min) && !isFinite(max))
        min = max = 0;
    else if (!isFinite(min))
        min = max;
    else if (!isFinite(max))
        max = min;

    return [ min, max ];
};

//hackish way to work with requirejs - doesn't actually work yet
//since we don't use the require optimizer to combine modules
if( typeof(define) === 'function' && define.amd) {
    define('$3Dmol',[], function() { return $3Dmol; });
}

/* StereoViewer for stereoscopic viewing
  @function $3Dmol.createStereoViewer
* @param {Object | string} element - Either HTML element or string identifier
* 
*/

$3Dmol.createStereoViewer = function(element) {
    var that = this;
    if(typeof(element) === "string")
        element = $("#"+element);
    if(!element) return;
    
    var viewers = $3Dmol.createViewerGrid(element, {rows: 1, cols: 2, control_all: true});
    
    this.glviewer1 = viewers[0][0];
    this.glviewer2 = viewers[0][1];
    
    this.glviewer1.setAutoEyeSeparation(false);
    this.glviewer2.setAutoEyeSeparation(true);    

    this.glviewer1.linkViewer(this.glviewer2);
    this.glviewer2.linkViewer(this.glviewer1);

    var methods = Object.getOwnPropertyNames(this.glviewer1) //get all methods of glviewer object
    .filter(function(property) {
        return typeof that.glviewer1[property] == 'function';
    });

    for (var i = 0; i < methods.length; i++) { //create methods of the same name
        this[methods[i]] = (function(method){
            return function(){
                var m1=this.glviewer1[method].apply(this.glviewer1,arguments);
                var m2=this.glviewer2[method].apply(this.glviewer2,arguments);
                return [m1,m2];
            };
        })(methods[i]);
    }
    
    //special cased methods
    this.setCoordinates = function (models, data, format) { //for setting the coordinates of the models
        for (var i = 0; i < models.length; i++) {
            models[i].setCoordinates(data, format);
        }
    };
    
    this.surfacesFinished = function() {
        return this.glviewer1.surfacesFinished() && this.glviewer2.surfacesFinished();
    };
    
    this.isAnimated = function() {
        return this.glviewer1.isAnimated() || this.glviewer2.isAnimated();
    };
    
    this.render = function(callback) {
        this.glviewer1.render();
        this.glviewer2.render();
        if(callback) {
            callback(this); //call only once
        }
    };
    
    this.getCanvas = function() {
        return this.glviewer1.getCanvas(); //same for both
    };

};

//simplified version of $.extend
$3Dmol.extend = function (obj1, src1) {
    for (var key in src1) {
        if(src1.hasOwnProperty(key) && src1[key] !== undefined) {
            obj1[key] = src1[key];            
        }
    }   
    return obj1;
}; 

//deep copy, cannot deal with circular refs; undefined input becomes an empty object
//https://medium.com/javascript-in-plain-english/how-to-deep-copy-objects-and-arrays-in-javascript-7c911359b089
$3Dmol.deepCopy = function(inObject)  {
  let outObject, value, key;

  if ( inObject == undefined) {
    return {};
  }
  if (typeof inObject !== "object" || inObject === null) {
    return inObject; // Return the value if inObject is not an object
  }

  // Create an array or object to hold the values
  outObject = Array.isArray(inObject) ? [] : {};

  for (key in inObject) {
    value = inObject[key];
    // Recursively (deep) copy for nested objects, including arrays
    outObject[key] = $3Dmol.deepCopy(value);
  }

  return outObject;
};

$3Dmol.isNumeric = function( obj ) {

    var type = typeof( obj );
    return ( type === "number" || type === "string" ) &&
        !isNaN( obj - parseFloat( obj ) );
};

$3Dmol.isEmptyObject = function( obj ) {
    var name;
    for ( name in obj ) {
        return false;
    }
    return true;
};

$3Dmol.makeFunction = function(callback) {
    //for py3dmol let users provide callback as string
    if (callback && typeof callback === "string") {
    /* jshint ignore:start */
        callback = eval("("+callback+")");
    /* jshint ignore:end */
    }
    // report to console if callback is not a valid function
    if (callback && typeof callback != "function") {
        return null;
    }    
    return callback;
};

//standardize voldata/volscheme in style
$3Dmol.adjustVolumeStyle = function(style) {
    if(style) {
        if(style.volformat && !(style.voldata instanceof $3Dmol.VolumeData)) {
            style.voldata = new $3Dmol.VolumeData(style.voldata, style.volformat);
        }
        if(style.volscheme) {
            style.volscheme = $3Dmol.Gradient.getGradient(style.volscheme);
        }
    }
};
