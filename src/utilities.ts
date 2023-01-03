//a collection of miscellaneous utility functions

import { builtinGradients, Gradient } from "./Gradient";
import { VolumeData } from "./VolumeData";
import { builtinColorSchemes, CC, elementColors, htmlColors, Color } from "./colors";
import { IsoSurfaceSpec } from "GLShape";

//simplified version of jquery extend
export function extend(obj1, src1) {
    for (var key in src1) {
        if (src1.hasOwnProperty(key) && src1[key] !== undefined) {
            obj1[key] = src1[key];
        }
    }
    return obj1;
};

//deep copy, cannot deal with circular refs; undefined input becomes an empty object
//https://medium.com/javascript-in-plain-english/how-to-deep-copy-objects-and-arrays-in-javascript-7c911359b089
export function deepCopy(inObject) {
    let outObject, value, key;

    if (inObject == undefined) {
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
        outObject[key] = deepCopy(value);
    }

    return outObject;
};

export function isNumeric(obj) {

    var type = typeof (obj);
    return (type === "number" || type === "string") &&
        !isNaN(obj - parseFloat(obj));
};

export function isEmptyObject(obj) {
    var name;
    for (name in obj) {
        return false;
    }
    return true;
};

export type Func = Function|string|undefined|null;

export function makeFunction(callback:Func) {
    //for py3dmol let users provide callback as string
    if (callback && typeof callback === "string") {
        /* jshint ignore:start */
        callback = eval("(" + callback + ")");
        /* jshint ignore:end */
    }
    // report to console if callback is not a valid function
    if (callback && typeof callback != "function") {
        return null;
    }
    return callback;
};

//standardize voldata/volscheme in style
export function adjustVolumeStyle(style: IsoSurfaceSpec) {
    if (style) {
        if (style.volformat && !(style.voldata instanceof VolumeData)) {
            style.voldata = new VolumeData(style.voldata, style.volformat);
        }
        if (style.volscheme) {
            style.volscheme = Gradient.getGradient(style.volscheme);
        }
    }
};


/*
 * computes the bounding box around the provided atoms
 * @param {AtomSpec[]} atomlist
 * @return {Array}
 */
export function getExtent(atomlist, ignoreSymmetries?) {
    var xmin, ymin, zmin, xmax, ymax, zmax, xsum, ysum, zsum, cnt;
    var includeSym = !ignoreSymmetries;

    xmin = ymin = zmin = 9999;
    xmax = ymax = zmax = -9999;
    xsum = ysum = zsum = cnt = 0;

    if (atomlist.length === 0)
        return [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
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

    return [[xmin, ymin, zmin], [xmax, ymax, zmax],
    [xsum / cnt, ysum / cnt, zsum / cnt]];
};


/* get the min and max values of the specified property in the provided
* @function $3Dmol.getPropertyRange
* @param {AtomSpec[]} atomlist - list of atoms to evaluate
* @param {string} prop - name of property 
* @return {Array} - [min, max] values
*/
export function getPropertyRange(atomlist, prop) {
    var min = Number.POSITIVE_INFINITY;
    var max = Number.NEGATIVE_INFINITY;

    for (var i = 0, n = atomlist.length; i < n; i++) {
        var atom = atomlist[i];
        var val = getAtomProperty(atom, prop);

        if (val != null) {
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

    return [min, max];
};


//adapted from https://stackoverflow.com/questions/3969475/javascript-pause-settimeout
export class PausableTimer {
    ident: any;
    total_time_run = 0;
    start_time: number;
    countdown: number;
    fn: any;
    arg: any;

    constructor(fn, countdown, arg?) {
        this.fn = fn;
        this.arg = arg;
        this.countdown = countdown;
        this.start_time = new Date().getTime();
        this.ident = setTimeout(fn, countdown, arg);
    }

    cancel() {
        clearTimeout(this.ident);
    }

    pause() {
        clearTimeout(this.ident);
        this.total_time_run = new Date().getTime() - this.start_time;
    }

    resume() {
        this.ident = setTimeout(this.fn, Math.max(0, this.countdown - this.total_time_run), this.arg);
    }

};

/*
 * Convert a base64 encoded string to a Uint8Array
 * @param {string} base64 encoded string
 */
export function base64ToArray(base64) {
    var binary_string = window.atob(base64);
    var len = binary_string.length;
    var bytes = new Uint8Array(len);
    for (var i = 0; i < len; i++) {
        bytes[i] = binary_string.charCodeAt(i);
    }
    return bytes;
};

//return the value of an atom property prop, or null if non existent
// looks first in properties, then in the atom itself
export function getAtomProperty(atom, prop) {
    var val = null;
    if (atom.properties &&
        typeof (atom.properties[prop]) != "undefined") {
        val = atom.properties[prop];
    } else if (typeof (atom[prop]) != 'undefined') {
        val = atom[prop];
    }
    return val;
};

//Miscellaneous functions and classes - to be incorporated into $3Dmol proper
/*
 * 
 * @param {$3Dmol.Geometry} geometry
 * @param {$3Dmol.Mesh} mesh
 * @returns {undefined}
 */
export function mergeGeos(geometry, mesh) {

    var meshGeo = mesh.geometry;

    if (meshGeo === undefined)
        return;

    geometry.geometryGroups.push(meshGeo.geometryGroups[0]);

};


/*
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
export function specStringToObject(str) {
    if (typeof (str) === "object") {
        return str; //not string, assume was converted already
    }
    else if (typeof (str) === "undefined" || str == null) {
        return str;
    }

    //if this is a json string, parse it directly
    try {
        let parsed = JSON.parse(str);
        return parsed;
    } catch (error) {

    }

    str = str.replace(/%7E/, '~'); //copy/pasting urls sometimes does this
    //convert things that look like numbers into numbers
    var massage = function (val) {
        if (isNumeric(val)) {
            //hexadecimal does not parse as float
            if (Math.floor(parseFloat(val)) == parseInt(val)) {
                return parseFloat(val);
            }
            else if (val.indexOf('.') >= 0) {
                return parseFloat(val); // ".7" for example, does not parseInt
            }
            else {
                return parseInt(val);
            }
        }
        //boolean conversions
        else if (val === 'true') {
            return true;
        }
        else if (val === 'false') {
            return false;
        }
        return val;
    };

    var ret = {};
    if (str === 'all') return ret;
    var fields = str.split(';');
    for (var i = 0; i < fields.length; i++) {
        var fv = fields[i].split(':');
        var f = fv[0];
        var val = {};
        var vstr = fv[1];
        if (vstr) {
            vstr = vstr.replace(/~/g, "=");
            if (vstr.indexOf('=') !== -1) {
                //has key=value pairs, must be object
                var kvs = vstr.split(',');
                for (var j = 0; j < kvs.length; j++) {
                    var kv = kvs[j].split('=', 2);
                    val[kv[0]] = massage(kv[1]);
                }
            }
            else if (vstr.indexOf(',') !== -1) {
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



function checkStatus(response) {
    if (!response.ok) {
        throw new Error(`HTTP ${response.status} - ${response.statusText}`);
    }
    return response;
}

/**
 * Fetch data from URL
 * 
 * @param uri URL
 * @param callback Function to call with data 
 */
export function get(uri, callback?) {
    var promise = fetch(uri).then(checkStatus).then((response) => response.text());
    if (callback)
        return promise.then(callback);
    else
        return promise;
}

/**
 * Download binary data (e.g. a gzipped file) into an array buffer and provide
 * arraybuffer to callback.
 * @param {string} uri - location of data
 * @param {Function} [callback] - Function to call with arraybuffer as argument.  
 * @param {string} [request] - type of request
 * @param {string} [postdata] - data for POST request
 * @return {Promise}
 */
export function getbin(uri, callback?, request?, postdata?) {
    var promise;
    if (request == "POST") {
        promise = fetch(uri, { method: 'POST', body: postdata })
            .then((response) => checkStatus(response))
            .then((response) => response.arrayBuffer());
    } else {
        promise = fetch(uri).then((response) => checkStatus(response))
            .then((response) => response.arrayBuffer());
    }

    if (callback) return promise.then(callback);
    else return promise;
};


/**
 * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
 * @param {string} query - String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
 * @param {GLViewer} viewer - Add new model to existing viewer
 * @param {Object} options - Specify additional options
 *                           format: file format to download, if multiple are available, default format is pdb
 *                           pdbUri: URI to retrieve PDB files, default URI is http://www.rcsb.org/pdb/files/
 * @param {Function} [callback] - Function to call with model as argument after data is loaded.
  
 * @return {GLModel} GLModel, Promise if callback is not provided
 * @example
 viewer.setBackgroundColor(0xffffffff);
       $3Dmol.download('pdb:2nbd',viewer,{onemol: true,multimodel: true},function(m) {
        m.setStyle({'cartoon':{colorscheme:{prop:'ss',map:$3Dmol.ssColors.Jmol}}});
       viewer.zoomTo();
       viewer.render(callback);
    });
 */
export function download(query, viewer, options, callback?) {
    var type = "";
    var pdbUri = "";
    var mmtfUri = "";
    var uri = "";
    var promise = null;
    var m = viewer.addModel();

    if (query.indexOf(':') < 0) {
        //no type specifier, guess
        if (query.length == 4) {
            query = 'pdb:' + query;
        } else if (!isNaN(query)) {
            query = 'cid:' + query;
        } else {
            query = 'url:' + query;
        }
    }
    if (query.substr(0, 5) === 'mmtf:') {
        pdbUri = options && options.pdbUri ? options.pdbUri : "https://mmtf.rcsb.org/v1.0/full/";
        query = query.substr(5).toUpperCase();
        uri = pdbUri + query;
        if (options && typeof options.noComputeSecondaryStructure === 'undefined') {
            //when fetch directly from pdb, trust structure annotations
            options.noComputeSecondaryStructure = true;
        }
        promise = new Promise(function (resolve) {
            getbin(uri)
                .then(function (ret) {
                    m.addMolData(ret, 'mmtf', options);
                    viewer.zoomTo();
                    viewer.render();
                    resolve(m);
                }, function () { console.log("fetch of " + uri + " failed."); });
        });
    }
    else {
        if (query.substr(0, 4) === 'pdb:') {
            type = 'mmtf';
            if (options && options.format) {
                type = options.format; //can override and require pdb
            }

            if (options && typeof options.noComputeSecondaryStructure === 'undefined') {
                //when fetch directly from pdb, trust structure annotations
                options.noComputeSecondaryStructure = true;
            }
            query = query.substr(4).toUpperCase();
            if (!query.match(/^[1-9][A-Za-z0-9]{3}$/)) {
                alert("Wrong PDB ID");
                return;
            }
            if (type == 'mmtf') {
                mmtfUri = options && options.mmtfUri ? options.mmtfUri : 'https://mmtf.rcsb.org/v1.0/full/';
                uri = mmtfUri + query.toUpperCase();
            }
            else {
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
        } else if (query.substr(0, 4) == 'url:') {
            uri = query.substr(4);
            type = uri;
        }

        var handler = function (ret) {
            m.addMolData(ret, type, options);
            viewer.zoomTo();
            viewer.render();
        };
        promise = new Promise(function (resolve) {
            if (type == 'mmtf') { //binary data
                getbin(uri)
                    .then(function (ret) {
                        handler(ret);
                        resolve(m);
                    }).catch(function () {
                        //if mmtf server is being annoying, fallback to text
                        pdbUri = options && options.pdbUri ? options.pdbUri : "https://files.rcsb.org/view/";
                        uri = pdbUri + query + ".pdb";
                        type = "pdb";
                        console.log("falling back to pdb format");
                        get(uri).then(function (data) {
                            handler(data);
                            resolve(m);
                        }).catch(function (e) {
                            handler("");
                            resolve(m);
                            console.log("fetch of " + uri + " failed: " + e.statusText);
                        });
                    }); //an error msg has already been printed
            }
            else {
                get(uri).then(function (data) {
                    handler(data);
                    resolve(m);
                }).catch(function (e) {
                    handler("");
                    resolve(m);
                    console.log("fetch of " + uri + " failed: " + e.statusText);
                });
            }
        });
    }
    if (callback) {
        promise.then(function (m) {
            callback(m);
        });
        return m;
    }
    else return promise;
};


/* Return proper color for atom given style
 * @param {AtomSpec} atom
 * @param {AtomStyle} style
 * @return {Color}
 */
export function getColorFromStyle(atom, style): Color {
    let scheme = style.colorscheme;
    if (typeof builtinColorSchemes[scheme] != "undefined") {
        scheme = builtinColorSchemes[scheme];
    } else if (typeof scheme == "string" && scheme.endsWith("Carbon")) {
        //any color you want of carbon
        let ccolor = scheme
            .substring(0, scheme.lastIndexOf("Carbon"))
            .toLowerCase();
        if (typeof htmlColors[ccolor] != "undefined") {
            let newscheme = { ...elementColors.defaultColors };
            newscheme.C = htmlColors[ccolor];
            builtinColorSchemes[scheme] = { prop: "elem", map: newscheme };
            scheme = builtinColorSchemes[scheme];
        }
    }

    let color = atom.color;
    if (typeof style.color != "undefined" && style.color != "spectrum")
        color = style.color;
    if (typeof scheme != "undefined") {
        let prop, val;
        if (typeof elementColors[scheme] != "undefined") {
            //name of builtin colorscheme
            scheme = elementColors[scheme];
            if (typeof scheme[atom[scheme.prop]] != "undefined") {
                color = scheme.map[atom[scheme.prop]];
            }
        } else if (typeof scheme[atom[scheme.prop]] != "undefined") {
            //actual color scheme provided
            color = scheme.map[atom[scheme.prop]];
        } else if (
            typeof scheme.prop != "undefined" &&
            typeof scheme.gradient != "undefined"
        ) {
            //apply a property mapping
            prop = scheme.prop;
            var grad = scheme.gradient; //redefining scheme
            if (typeof builtinGradients[grad] != "undefined") {
                grad = new builtinGradients[grad](
                    scheme.min,
                    scheme.max,
                    scheme.mid ? scheme.mid : scheme.colors
                );
            }

            let range = grad.range() || [-1, 1]; //sensible default
            val = getAtomProperty(atom, prop);
            if (val != null) {
                color = grad.valueToHex(val, range);
            }
        } else if (
            typeof scheme.prop != "undefined" &&
            typeof scheme.map != "undefined"
        ) {
            //apply a discrete property mapping
            prop = scheme.prop;
            val = getAtomProperty(atom, prop);
            if (typeof scheme.map[val] != "undefined") {
                color = scheme.map[val];
            }
        } else if (typeof style.colorscheme[atom.elem] != "undefined") {
            //actual color scheme provided
            color = style.colorscheme[atom.elem];
        } else {
            console.log("Could not interpret colorscheme " + scheme);
        }
    } else if (typeof style.colorfunc != "undefined") {
        //this is a user provided function for turning an atom into a color
        color = style.colorfunc(atom);
    }

    let C = CC.color(color);
    return C;
};

//given a string selector, element, or jquery object, return the HTMLElement
export function getElement(element): HTMLElement | null {
    let ret = element;
    if (typeof (element) === "string") {
        ret = document.querySelector("#" + element);
    } else if (typeof element === 'object' && element.get) { //jquery
        ret = element.get(0);
    }
    return ret;
}