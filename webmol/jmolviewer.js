//a molecular viewer based on GLMol

var WebMol = WebMol || {};

// a webmol unified interace to gmol
WebMol.jmolViewer = (function() {
	// private class variables


	// private class helper functions
	
	//intialize jmol applet inside container which is a jquery element
	function initializeJMol(container)
	{
		// jmol initialization
		var usedSigned = (document.location.search.indexOf("SIGNEDJMOL") >= 0);
		jmolSetDocument(false);
		jmolInitialize("Jmol/appletweb", usedSigned);
		jmolSetAppletColor("#222222");
		var init = "set statusReporting off; frank off;";
		init += "set ambient 40; set specpower 40; wireframe only;";
		init += "color labels blue;";
		init += "set showHydrogens false;";
		init += "set antialiasDisplay on;";
		init += "set antialiasTranslucent on;";
		init += "set background [xf8f8f8];";
		init += "data \"model dummy\"|1|dummy|C 0 0 0end \"model dummy\"; hide model=1.1;";
		init += "javascript setAppletLoaded();"; //jmolSetCallback does not work

		var appletstr = jmolApplet([ 1, 1 ], init);
		container.html(appletstr);
	}

	// The constructor
	function jmolViewer(element, width, height, callback) {
		// check dependencies
		if (typeof (Jmol) === "undefined") {
			// three.js not loaded, take matters into our own hands
			throw "Missing JMol.js";
		}


		// set variables
		var container = element;
		var id = container.id;
		
		var models = []; //atomistic molecular models

		var Info = {
				  addSelectionOptions: false,
				  color: "#FFFFFF",
				  debug: false,
				  defaultModel: "",
				  height: height,
				  isSigned: false,
				  jarFile: "JmolApplet.jar",
				  jarPath: "Jmol/appletweb",
				  memoryLimit: 512,
				  readyFunction: null,
				  script: "frank off",
				  src: null,
				  use: "Java noWebGL noHTML5 noImage",
				  width: width
				};	 

		Jmol.setDocument(false);
		var japp = Jmol.getApplet("japp", Info) 

		container.html(Jmol.getAppletHtml(japp));

		// public methods
		this.setBackgroundColor = function(hex, a) {
	
		};

		this.setWidth = function(w) {

		};

		this.setHeight = function(h) {
	
		};


		// return specified model
		this.getModel = function(id) {
			return models[id];
		};

		this.getView = function() {

		}

		this.setView = function(arg) {

		}
		
		// apply styles, models, etc in viewer
		this.render = function() {

		};
		
		//zoom to atom selection
		this.zoomTo = function(sel) {

		};
		
		//given molecular data and its format (pdb, sdf or xyz)
		//create a model and add it, returning the model identifier
		this.addModel = function(data, format) {
			var m;
			return m;
		};
		
		this.removeModel = function(model) {

		}

		//apply sel to all models and apply style
		this.setStyle = function(style, sel) {

		};
		
		//add a surface
		this.addSurface = function(type, style, atomsel, allsel) {
			
		}

		this.removeSurface = function(surf) {
			
		}
	}

	return jmolViewer;
})();
