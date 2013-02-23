//a molecular viewer based on GLMol

var WebMol = WebMol || {};

// a webmol unified interace to gmol
WebMol.jmolViewer = (function() {
	// private class variables
	var instance = 0;
	// private class helper functions

	// The constructor
	function jmolViewer(element, width, height, callback) {
		var japp = null;
		var container = element;
		var models = []; // atomistic molecular models

		// check dependencies
		if (typeof (Jmol) === "undefined") {
			// three.js not loaded, take matters into our own hands
			throw "Missing JMol.js";
		}

		// public methods
		this.setBackgroundColor = function(hex, a) {
			var hex = hex.toString(16);
			// must have padding zeroes
			hex = "000000".substr(0, 6 - hex.length) + hex;
			Jmol.script(japp, "set backgroundColor \"[x" + hex + "]\"");
		};

		this.setWidth = function(w) {
			var h = container.height();
			Jmol.resizeApplet(japp, [ w, h ]);
		};

		this.setHeight = function(h) {
			var w = container.width();
			Jmol.resizeApplet(japp, [ w, h ]);
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
			for ( var i = 0; i < models.length; i++) {
				models[i].render();
			}

		};

		// zoom to atom selection
		this.zoomTo = function(sel) {
			var sel = sel || {};
			// apply to all models unless sell specifies a model by id
			var ms = [];
			if (typeof sel.model == "undefined") {
				for(var i = 0; i < models.length; i++) {
					if(models[i]) ms.push(i);
				}
			} else { // specific to some models
				var ms = sel.model;
				if (!$.isArray(ms))
					ms = [ ms ];
			}

			var ors = []; // combine all with or
			for ( var i = 0; i < ms.length; i++) {
				if (typeof models[ms[i]] != "undefined") {
					var m = models[ms[i]];
					ors.push("(" + m.jmolSelect(sel) + ")");
				}
			}

			var script = "zoomto " + ors.join(" or ");
			Jmol.script(japp, script);

		};

		// given molecular data and its format (pdb, sdf or xyz)
		// create a model and add it, returning the model identifier
		this.addModel = function(data, format) {
			var m = new WebMol.jmolModel(japp, models.length);
			m.addMolData(data, format);
			models.push(m);
			return m;
		};

		this.removeModel = function(model) {
			var script = "delete " + model.jmolSelect({});
			delete models[model.getID()];
			Jmol.script(japp, script);
		}

		// apply sel to all models and apply style
		this.setStyle = function(style, sel) {
			var sel = sel || {};
			var ms = [];
			if (typeof sel.model == "undefined") {
				for(var i = 0; i < models.length; i++) {
					if(models[i]) ms.push(i);
				}
			} else { // specific to some models
				var ms = sel.model;
				if (!$.isArray(ms))
					ms = [ ms ];
			}

			for ( var i = 0; i < ms.length; i++) {
				if (typeof models[ms[i]] != "undefined") {
					var m = models[ms[i]];
					m.setStyle(style, sel);
				}
			}
		};

		// add a surface
		this.addSurface = function(type, style, atomsel, allsel) {

		}

		this.removeSurface = function(surf) {

		}


		var Info = {
			addSelectionOptions : false,
			color : "#FFFFFF",
			debug : false,
			defaultModel : "",
			height : height,
			isSigned : false,
			jarFile : "JmolApplet.jar",
			jarPath : "Jmol/appletweb",
			memoryLimit : 512,
			readyFunction : null,
			script : "frank off; ",
			src : null,
			use : "Java noWebGL noHTML5 noImage",
			width : width
		};

		if (typeof (callback) === "function" && callback) {
			// have to let the java applet initialize before doing anything with
			// it
			// readyFunction doesn't work for multiple reasons
			var callbackname = "__jmolInitHack" + instance;
			var self = this;
			WebMol[callbackname] = function() {
				// delete WebMol[callbackname];
				callback(self);
			};
			Info.script += "javascript WebMol." + callbackname + "();";

		}
		Jmol.setDocument(false);
		japp = Jmol.getApplet("japp" + instance, Info)
		instance++;

		container.html(Jmol.getAppletHtml(japp));

		$(window).resize(function() { // only window can capture resize event
			var w = container.width();
			var h = container.height();
			Jmol.resizeApplet(japp, [ w, h ]);
		});

	}

	return jmolViewer;
})();
