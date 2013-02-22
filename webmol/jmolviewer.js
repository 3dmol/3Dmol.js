//a molecular viewer based on GLMol

var WebMol = WebMol || {};

// a webmol unified interace to gmol
WebMol.jmolViewer = (function() {
	// private class variables


	// private class helper functions
	

	// The constructor
	function jmolViewer(element, width, height, callback) {
		// check dependencies
		if (typeof (jmolScript) === "undefined") {
			// three.js not loaded, take matters into our own hands
			throw "Missing JMol.js";
		}


		// set variables
		var container = element;
		var id = container.id;
		
		var models = []; //atomistic molecular models


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

	return GLViewer;
})();
