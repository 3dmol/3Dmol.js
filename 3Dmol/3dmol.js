
//This defines the $3Dmol object which is used to create viewers
//and configure system-wide settings

$3Dmol = (function(window) {
    
    var my = window['$3Dmol'] || {};
    //var $ = window['jQuery'];
    
    return my;

})(window);
    
$3Dmol.createViewer = function(element, config)
{
    if($.type(element) === "string")
        element = $("#"+element);
    if(!element) return;

    config = config || {};
 
    
    if(!config.defaultcolors)
        config.defaultcolors = $3Dmol.elementColors.defaultColors;

    //try to create the  viewer
    try {
    	return new $3Dmol.GLViewer(element, config.callback, config.defaultcolors, config.nomouse);
    }
    catch(e) {
    	console.log("error with "+kind+":"+e);
    }
    
    alert("Unable to instantiate 3Dmol viewer: "+config.order);
    return null;
};
   
$3Dmol.download = function(query, viewer) {
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
       

$3Dmol.SurfaceType = {
    VDW : 1,
    MS : 2,
    SAS : 3,
    SES  : 4
};

// in an attempt to reduce memory overhead, cache all $3Dmol.Colors
// this makes things a little faster
$3Dmol.CC = {
    cache : {},
    color : function(hex) {

        if(typeof(this.cache[hex]) !== "undefined") {
            return this.cache[hex];
        }
        else {
            hex = this.getHex(hex);
            var c = new $3Dmol.Color(hex);
            this.cache[hex] = c;
            return c;
        }
    },
    colorTab : {
        'white' : 0xFFFFFF,
        'silver' : 0xC0C0C0,
        'gray' : 0x808080,
        'grey' : 0x808080,
        'black' : 0x000000,
        'red' : 0xFF0000,
        'maroon' : 0x800000,
        'yellow' : 0xFFFF00,
        'orange' : 0xFF6600,
        'olive' : 0x808000,
        'lime' : 0x00FF00,
        'green' : 0x008000,
        'aqua' : 0x00FFFF,
        'cyan' : 0x00FFFF,
        'teal' : 0x008080,
        'blue' : 0x0000FF,
        'navy' : 0x000080,
        'fuchsia' : 0xFF00FF,
        'magenta' : 0xFF00FF,
        'purple' : 0x800080
    },    
    getHex : function(hex) {
        if (parseInt(hex))
            return parseInt(hex);
        
        else if (typeof(hex) === 'string') {
            
            return this.colorTab[hex.trim().toLowerCase()] || 0x000000;
        }
        
    }
    
};



$3Dmol['CC'] = $3Dmol.CC;
$3Dmol['CC']['color'] = $3Dmol.CC.color;

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

//Synchronized (i.e. not threaded) surface gen? Used mainly for debugging
$3Dmol.syncSurface = false;

// Internet Explorer refuses to allow webworkers in data blobs.  I can find
// no way of checking for this feature directly, so must do a brower check
if(window.navigator.userAgent.indexOf('MSIE ') >= 0 ||
		window.navigator.userAgent.indexOf('Trident/') >= 0) {
	$3Dmol.syncSurface = true; // can't use webworkers
}

/**
 * Parse a string that represents a style or atom selection and convert it
 * into an object.  
 */
$3Dmol.specStringToObject = function(str) {
	if(typeof(str) === "object") {
		return str; //not string, assume was converted already
	}
	else if(typeof(str) === "undefined" || str == null) {
		return str; 
	}
	var ret = {};
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
					val[kv[0]] = kv[1];
				}
			}
			else if(vstr.indexOf(',') !== -1) {
				//has multiple values, must list
				val = vstr.split(',');
			}
			else {
				val = vstr; //value itself
			}
		}
		ret[f] = val;
	}
	
	return ret;
}



