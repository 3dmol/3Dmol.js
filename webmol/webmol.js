
//This defines the WebMol object which is used to create viewers
//and configure system-wide settings

//the presence of jquery is assumed
var WebMol = (function() {
	var my = {};
	var $ = jQuery; //avoid any conflicts with the dollar
	
	
	
	//create the best viewer according to parameters in config within element
	//config.width - width of viewer, if unset use html elment width
	//config.height - height, if unset use html elment width
	//config.order - preference for types of viewer, glmol or jmol
	//config.callback - for intialization commands to immediately apply to viewer
	//element can either be the html element object or its identifier
	my.createViewer = function(element, config)
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
			
			if(typeof(my[fname]) === "function")
			{
				try {
					return new my[fname](element, config.callback, config.defaultcolors);
				}
				catch(e) {
					console.log("error with "+kind+":"+e);
				}
			}
		}
		alert("Unable to instantiate webmol viewer: "+config.order);
		return null;
	};
	
	//loads a pdb/pubchem structure into the provided viewer
	my.download = function(query, viewer) {
		   var baseURL = '';
		   var type = "";
		   if (query.substr(0, 4) == 'pdb:') {
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
		};
	return my;
})();

WebMol.SurfaceType = {
		VDW : 1,
		SAS : 3,
		SES : 2
	};

// in an attempt to reduce memory overhead, cache all THREE.Colors
//this makes things a little faster
WebMol.CC = {
	cache : {},
	color : function(hex) {
		if(typeof(this.cache[hex]) != "undefined") {
			return this.cache[hex];
		}
		else {
			var c = new THREE.Color(hex);
			this.cache[hex] = c;
			return c;
		}
	}
};
	