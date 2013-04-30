//a molecular viewer based on GLMol

var WebMol = WebMol || {};

// a webmol unified interace to gmol
WebMol.jmolViewer = (function() {
	// private class variables
	var instance = 0;
	// private class helper functions

	function jmolColor(c) {
		var hex = c.toString(16);
		// must have padding zeroes
		hex = "000000".substr(0, 6 - hex.length) + hex;
		return "[x" + hex + "]";
	}

	// The constructor
	function jmolViewer(element, callback, defaultcolors) {
		var japp = null;
		var container = element;
		var models = []; // atomistic molecular models
		var surfaceCounter = 0;
		// check dependencies
		if (typeof (Jmol) === "undefined") {
			// three.js not loaded, take matters into our own hands
			throw "Missing JMol.js";
		}

		var defaultcolorsstr = "jmol";
		if (defaultcolors == WebMol.rasmolElementColors)
			defaultcolorsstr = "rasmol";

		// public methods
		this.setBackgroundColor = function(hex, a) {
			Jmol.script(japp, "set backgroundColor \"" + jmolColor(hex))+"\"";
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
			var orient = Jmol.getPropertyAsArray(japp,"orientationInfo");
			return orient;
		}

		this.setView = function(arg) {

		}

		// apply styles, models, etc in viewer
		this.render = function() {
			for ( var i = 0; i < models.length; i++) {
				models[i].render();
			}

		};

		// return atom expression representing sel in jmol
		function getJMolSel(sel) {
			var sel = sel || {};
			// apply to all models unless sell specifies a model by id
			var ms = [];
			if (typeof sel.model == "undefined") {
				for ( var i = 0; i < models.length; i++) {
					if (models[i])
						ms.push(i);
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
			return ors.join(" or ");
		}

		// zoom to atom selection
		this.zoomTo = function(sel) {
			var script = "zoomto 0 " + getJMolSel(sel);
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
		
		this.removeAllModels = function() {
			models = [];
			Jmol.script(japp, "delete");
		};
		
		var applyToModels = function(func, sel, value, value2) {
			var sel = sel || {};
			var ms = [];
			if (typeof sel.model == "undefined") {
				for ( var i = 0; i < models.length; i++) {
					if (models[i])
						ms.push(i);
				}
			} else { // specific to some models
				var ms = sel.model;
				if (!$.isArray(ms))
					ms = [ ms ];
			}

			for ( var i = 0; i < ms.length; i++) {
				if (typeof models[ms[i]] != "undefined") {
					var m = models[ms[i]];
					m[func](sel, value, value2);
				}
			}
		}
		
		this.setColorByProperty = function(sel, prop, scheme) {
			applyToModels("setColorByProperty", sel, prop, scheme);
		}
		
		this.setColorByElement = function(sel, colors) {
			applyToModels("setColorByElement", sel, colors);
		}

		// apply sel to all models and apply style
		this.setStyle = function(sel, style) {
			applyToModels("setStyle", sel, style);
		};

		// add a surface
		this.addSurface = function(type, style, atomsel, allsel) {
			var surfid = "id" + surfaceCounter++;
			var s = getJMolSel(atomsel);
			var a = getJMolSel(allsel);

			var ST = WebMol.SurfaceType;
			var t = "molecular";
			switch (type) {
			case ST.VDW:
				t = "vdw";
				break;
			case ST.SAS:
				t = "sasurface";
				break;
			case ST.SES:
				t = "solvent 1.4";
				break;
			}
			

			
			var script = "isosurface " + surfid + " select \"" + s
					+ "\" ignore \"not (" + a + ")\" " + t + " ";

			if(style.map && style.map.prop) {
				//map color space using already set atom properties
				var prop = style.map.prop;
				var scheme = style.map.scheme || new WebMol.RWB();
				script += " colorscheme \""+scheme.jmolID()+"\" ";
				var range = scheme.range();
				if(range) {
					script += " color range "+range[0]+" "+range[1] + " ";
				}
				script += "map property " + prop+" ";
			}
			
			if (style.color) {
				script += "color " + jmolColor(style.color);
				if(style.opacity)
					script += " translucent "+(1-style.opacity);
			}
			else if (style.opacity)
				script += "color translucent "
						+ (1 - style.opacity);

			console.log(script);
			Jmol.script(japp, script);
			return surfid;
		};

		this.setSurfaceMaterialStyle = function(surfid, style) {
			var script = "";
			if (style.color)
				script += "color $" + surfid + " " + jmolColor(style.color)
						+ ";";
			if (!isNaN(style.opacity))
			{
				if(style.opacity == 0) 
				{ //jmol won't do a fully translucent surface
					script += "isosurface "+surfid+" hide;";
				}
				else
				{
					script += "isosurface "+surfid+" display;";
					script += "color $" + surfid + " translucent "
						+ (1 - style.opacity) + ";";
				}
			}
			console.log(script);
			Jmol.script(japp, script);
		};

		this.removeSurface = function(surf) {
			Jmol.script(japp, "isosurface " + surf + " delete");
		};

		// take a list of property objects that define selections of atoms
		// and store the named properties in the atom
		this.mapAtomProperties = function(props) {
			var script = "";
			for ( var i = 0; i < props.length; i++) {
				var p = props[i];
				if (p.props) {
					var sel = "{" + getJMolSel(p) + "}.";
					for( var prop in p.props) {
						if(p.props.hasOwnProperty(prop)) {
							script += sel + prop + " = " + p.props[prop] + ";";
						}
					}
				}
			}
			Jmol.script(japp, script);
		};
		
		this.pdbData = function(sel) {
			Jmol.scriptWait(japp, "select "+getJMolSel(sel) + "; ");
			var data = Jmol.scriptEcho(japp,"write pdb");
			// jmol doesn't output the atom names correctly, the element symbol
			// is
			// suppose to
			// be right justified to column 14 (to be fair, this is far from
			// clear
			// in the documentation)
			var lines = data.split("\n");
			var newlines = [];
			for ( var i = 0, n = lines.length; i < n; i++) {
				var l = lines[i];
				if (l.charAt(12) != " " && !l.charAt(13).match(/[a-z]/)) {
					// need to shift
					l = l.substr(0, 12) + " " + l.substr(12, 3) + l.substr(16);
				}
				if (l.substr(0, 4) == "ATOM")
					newlines.push(l);
			}
			return newlines.join("\n");
		};

		var Info = {
			addSelectionOptions : false,
			color : "#FFFFFF",
			debug : false,
			defaultModel : "",
			height : "100%",
			isSigned : false,
			jarFile : "JmolApplet.jar",
			jarPath : "Jmol/appletweb",
			memoryLimit : 512,
			readyFunction : null,
			script : "frank off; set showHydrogens false;",
			src : null,
			use : "Java noWebGL noHTML5 noImage",
			width : "100%"
		};

		// have to let the java applet initialize before doing anything with
		// it
		// readyFunction doesn't work for multiple reasons
		var callbackname = "__jmolInitHack" + instance;
		var self = this;
		WebMol[callbackname] = function() {
			Jmol.script(japp, "set defaultcolorscheme " + defaultcolorsstr);
			delete WebMol[callbackname];
			if (typeof (callback) === "function" && callback)
				callback(self);
		};
		Info.script += "javascript WebMol." + callbackname + "();";

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
