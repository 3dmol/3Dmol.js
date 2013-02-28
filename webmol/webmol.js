
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
		if(!config.width)
			config.width = $(element).width();
		if(!config.height)
			config.height = $(element).height();
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
					return new my[fname](element, config.width, config.height, config.callback, config.defaultcolors);
				}
				catch(e) {
					console.log("error with "+kind+":"+e);
				}
			}
		}
		alert("Unable to instantiate webmol viewer: "+config.order);
		return null;
	}
	
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
	