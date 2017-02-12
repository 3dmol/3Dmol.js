
//This defines the $3Dmol object which is used to create viewers
//and configure system-wide settings

/** 
 * All of the functionality of $3Dmol.js is contained within the
 * $3Dmol global namespace
 * @namespace */
$3Dmol = (function(window) {
    
    var my = window['$3Dmol'] || {};
    //var $ = window['jQuery'];
    
    return my;

})(window);

/* The following code "phones home" to register that an ip 
   address has loaded 3Dmol.js.  Being able track this usage
   is very helpful when reporting to funding agencies.  Please
   leave this code in if you would like to increase the 
   likelihood of 3Dmol.js remaining supported.
*/
$.get("http://3dmol.csb.pitt.edu/track/report.cgi");

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
                   if (window.FormData
                           && ((options.dataType && (options.dataType == 'binary')) || (options.data && ((window.ArrayBuffer && options.data instanceof ArrayBuffer) || (window.Blob && options.data instanceof Blob))))) {
                       return {
                           // create new XMLHttpRequest
                           send : function(headers, callback) {
                               // setup all variables
                               var xhr = new XMLHttpRequest(), url = options.url, type = options.type, async = options.async || true,
                               // blob or arraybuffer. Default is blob
                               dataType = options.responseType || "blob", data = options.data
                                       || null, username = options.username
                                       || null, password = options.password
                                       || null;

                               xhr.addEventListener('load', function() {
                                   var data = {};
                                   data[options.dataType] = xhr.response;
                                   // make callback and send data
                                   callback(xhr.status, xhr.statusText,
                                           data,
                                           xhr.getAllResponseHeaders());
                               });

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
 * @param {Object | string} element - Either HTML element or string identifier
 * @param {ViewerSpec} config Viewer specification
 * @return {$3Dmol.GLViewer} GLViewer, null if unable to instantiate WebGL
 * @example
   $3Dmol.download("pdb:4UAA",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  var colorAsSnake = function(atom) {
                    return atom.resi % 2 ? 'white': 'green'
                  };

                  viewer.setStyle( {}, { cartoon: {colorfunc: colorAsSnake }});

                  viewer.render();
              });                     
 *                        
 */
$3Dmol.createViewer = function(element, config)
{
    if($.type(element) === "string")
        element = $("#"+element);
    if(!element) return;

    config = config || {}; 

    //try to create the  viewer
    try {
        return new $3Dmol.GLViewer(element, config);
    }
    catch(e) {
        throw "error creating viewer: "+e;
    }
    
    return null;
};
   
/**
 * Contains a dictionary of embedded viewers created from HTML elements
 * with a the viewer_3Dmoljs css class indexed by their id (or numerically
 * if they do not have an id).
*/
$3Dmol.viewers = {};

/**
 * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
 * 
 * @param {string} query - String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
 * @param {$3Dmol.GLViewer} viewer - Add new model to existing viewer
 * @param {Object} options - Specify additional options
 *                           format: file format to download, if multiple are available, default format is pdb
 *                           pdbUri: URI to retrieve PDB files, default URI is http://www.rcsb.org/pdb/files/
 * @param {Function} callback - Function to call with model as argument after data is loaded.

 * @example
 viewer.setBackgroundColor(0xffffffff);
       $3Dmol.download('pdb:2nbd',viewer,{onemol: true,multimodel: true},function(m) {
        m.setStyle({'cartoon':{colorscheme:{prop:'ss',map:$3Dmol.ssColors.Jmol}}});
       viewer.zoomTo();
       viewer.render(callback);
    });
 * @return {$3Dmol.GLModel} GLModel
 */ 
$3Dmol.download = function(query, viewer, options, callback) {
    var baseURL = '';
    var type = "";
    var pdbUri = "";
    var mmtfUri = "";
    var m = viewer.addModel();
    
    if (query.substr(0, 5) === 'mmtf:') {
        pdbUri = options && options.pdbUri ? options.pdbUri : "https://mmtf.rcsb.org/v1.0/full/";
        query = query.substr(5).toUpperCase();
        var uri = pdbUri + query;        
        if(options && typeof options.noComputeSecondaryStructure === 'undefined') {
                //when fetch directly from pdb, trust structure annotations
                options.noComputeSecondaryStructure = true;
        }
            
        $.ajax({url:uri, 
            type: "GET",
            dataType: "binary",
            responseType: "arraybuffer",
            processData: false}).done(
                function(ret, txt, response) {
                    m.addMolData(ret, 'mmtf',options);
                    viewer.zoomTo();
                    viewer.render();
                    if(callback) callback(m);
                }).fail(function(e,txt) { 
                    console.log(txt);
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
          if(callback) callback(m);
        };
        
        if(type == 'mmtf') { //binary data
            $.ajax({url:uri, 
            type: "GET",
            dataType: "binary",
            responseType: "arraybuffer",
            processData: false}).done(handler).fail(function(e,txt) { 
                    console.log(txt);
            });
        }
        else {        
           $.get(uri, handler).fail(function(e) {
            console.log("fetch of "+uri+" failed: "+e.statusText);
           });
        }
   }
   
   return m;
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
// no way of checking for this feature directly, so must do a brower check
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
        if($.isNumeric(val)) {
           //hexadecimal does not parse as float
           if(Math.floor(parseFloat(val)) == parseInt(val)) {
              return parseFloat(val);
           }
           else if(val.indexOf('.') >= 0) {
               return parseFloat(val); // ".7" for example, does not parseInt
           }
           else {
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
    }
    
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
}


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
//looks first in properties, then in the atom itself
$3Dmol.getAtomProperty = function(atom, prop) {
    var val = null;
    if (atom.properties
            && typeof (atom.properties[prop]) != "undefined") {
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
}

//hackish way to work with requirejs - doesn't actually work yet
//since we doing use the require optimizer to combine modules
if( typeof(define) === 'function' && define.amd) {
    define('$3Dmol',$3Dmol);
}

