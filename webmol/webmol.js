
//This defines the WebMol object which is used to create viewers
//and configure system-wide settings

//TODO: automated build system.  (doc gen system, too?)

//the presence of jquery is assumed
/** 
 * WebMol global namespace
 * 
 * @namespace
 * 
 * @property JmolElementColors - Jmol style atom colors (default color scheme)
 * @property rasmolElementColors - Rasmol style atom colors
 */
var WebMol = (function(window) {

    var my = {};
    var $ = window['jQuery'];
    //window['WebMol'] = my;
        
    return my;

})(window);
    
/**
 * Create and initialize an appropriate viewer at supplied HTML element using specification in config
 * @function WebMol.createViewer
 * @param {Object | string} element Either HTML element or string identifier
 * @param {Object} config Viewer specification
 * @returns {WebMol.GLViewer} GLViewer
 * 
 * @example
 * // Assume there exists an HTML div with id "gldiv"
 * var element = $("#gldiv");
 * 
 * // Viewer config - properties 'defaultcolors' and 'callback'
 * var config = {defaultcolors: WebMol.rasmolElementColors,
 *               callback : function(viewer) {
 *                            //'data' is a string containing molecule data in pdb format  
 *                            viewer.addModel(data, "pdb");
 *                            viewer.zoomTo();
 *                            viewer.render();
 *                          }  
 *                        
 *               };
 * 
 * // Create GLViewer within 'gldiv' and execute callback
 * var myviewer = WebMol.createViewer(element, config);
 *      
 */
WebMol.createViewer = function(element, config)
{
    if($.type(element) === "string")
        element = $("#"+element);
    if(!element) return;

    config = config || {};
    if(!config.order)
        config.order = ["glmol","jmol"];
    if(!config.defaultcolors)
        config.defaultcolors = WebMol.defaultElementColors;

    //try to create the appropriate viewer
    for(var i = 0; i < config.order.length; i++) {
        var kind = config.order[i];
        var fname =kind+"Viewer";

        if(typeof(this[fname]) === "function")
        {
            try {
                return new this[fname](element, config.callback, config.defaultcolors);
            }
            catch(e) {
                console.log("error with "+kind+":"+e);
            }
        }
    }
    alert("Unable to instantiate webmol viewer: "+config.order);
    return null;
};

/**
 * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
 * 
 * @function WebMol.download
 * @param {string} query String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
 * @param {Object} viewer Add new model to existing viewer
 * @example
 * var myviewer = WebMol.createViewer(gldiv);
 * 
 * // GLModel 'm' created and loaded into glviewer for PDB id 2POR
 * var m = WebMol.download('pdb: 2POR', myviewer);
 * 
 * @returns {WebMol.GLModel} GLModel
 */    
WebMol.download = function(query, viewer) {
    var baseURL = '';
    var type = "";
    var m = null;
    if (query.substr(0, 4) === 'pdb:') {
        type = "pdb";
        query = query.substr(4).toUpperCase();
        if (!query.match(/^[1-9][A-Za-z0-9]{3}$/)) {
           alert("Wrong PDB ID"); return;
        }
        uri = "http://www.pdb.org/pdb/files/" + query + ".pdb";
    } else if (query.substr(0, 4) == 'cid:') {
        type = "sdf";
        query = query.substr(4);
        if (!query.match(/^[1-9]+$/)) {
           alert("Wrong Compound ID"); return;
        }
        uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + query + 
          "/SDF?record_type=3d";
    }

   $.get(uri, function(ret) {
      viewer.addModel(ret, type);
      viewer.zoomTo();
      viewer.render();                            
   });
   
   return m;
};
       

WebMol.SurfaceType = {
    VDW : 1,
    MS : 2,
    SAS : 3,
    SES  : 4
};

// in an attempt to reduce memory overhead, cache all WebMol.Colors
//this makes things a little faster
WebMol.CC = {
    cache : {},
    color : function(hex) {
        if(typeof(this.cache[hex]) !== "undefined") {
            return this.cache[hex];
        }
        else {
            var c = new WebMol.Color(hex);
            this.cache[hex] = c;
            return c;
        }
    }
};

WebMol.Color = function( color ){
    
    if ( arguments.length > 1) {
            this.r = arguments[0] || 0.0;
            this.g = arguments[1] || 0.0;
            this.b = arguments[2] || 0.0;

            return this;
    }
    
    return this.set(color);
                
};

WebMol.Color.prototype = {
    
    constructor: WebMol.Color,
    
    r: 0.0, g: 0.0, b: 0.0,
    
    set : function(val) {
        
            if (val instanceof WebMol.Color) 
                return val.clone();

            else if (typeof val === 'number')
                this.setHex(val);
    },
    
    setHex: function(hex) {
        
            hex = Math.floor(hex);

            this.r = (hex >> 16 & 255) / 255;
            this.g = (hex >> 8 & 255) / 255;
            this.b = (hex & 255) / 255;                                                                                     
        
            return this;
    },
    
    clone : function() {
            return new WebMol.Color(this.r, this.g, this.b);
    },
        
    copy : function(color) {
        this.r = color.r;
        this.g = color.g;
        this.b = color.b;
        
        return this;
    }
    
};

//Miscellaneous functions and classes - to be incorporated into WebMol proper

WebMol.mergeGeos = function(geometry, mesh) {
    
    var meshGeo = mesh.geometry;
    
    if (meshGeo === undefined) 
        return;
    
    geometry.geometryGroups.push( meshGeo.geometryGroups[0] );
    
};

WebMol.multiLineString = function(f) {
    return f.toString()
            .replace(/^[^\/]+\/\*!?/, '')
            .replace(/\*\/[^\/]+$/, '');
            
};

//Synchronized (i.e. not threaded) surface gen? Used mainly for debugging
WebMol.syncSurface = false;
//WebMol constants (replaces needed THREE constants)

//material constants

// sides
WebMol.FrontSide = 0;
WebMol.BackSide = 1;
WebMol.DoubleSide = 2;

// blending modes
WebMol.NoBlending = 0;
WebMol.NormalBlending = 1;
WebMol.AdditiveBlending = 2;
WebMol.SubtractiveBlending = 3;
WebMol.MultiplyBlending = 4;
WebMol.CustomBlending = 5;

// shading
WebMol.NoShading = 0;
WebMol.FlatShading = 1;
WebMol.SmoothShading = 2;

// colors
WebMol.NoColors = 0;
WebMol.FaceColors = 1;
WebMol.VertexColors = 2;

//Texture constants
//TODO: Which of these do I need (since I only use textures to display label sprites) ?
WebMol.MultiplyOperation = 0;
WebMol.MixOperation = 1;
WebMol.AddOperation = 2;

// mapping modes

WebMol.UVMapping = function() {};

// wrapping modes
WebMol.ClampToEdgeWrapping = 1001;

//Filters
WebMol.LinearFilter = 1006;
WebMol.LinearMipMapLinearFilter = 1008;

//Data types
WebMol.UnsignedByteType = 1009;

//Pixel formats
WebMol.RGBAFormat = 1021;

