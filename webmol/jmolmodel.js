// A model is a collection of related atoms.  Bonds are only allowed between
//atoms in the same model.  An atom is uniquely specified by its model id and
//its serial number.
//A jmolmodel corresponds to a model within jmol

var WebMol = WebMol || {};

var Jmol = Jmol || {};

//no idea why this didn't make the transition to jsmol
//put back ability to load inline data
Jmol.loadInline = function(jsapp, model) {
	  if (!model) return null;
	  var applet= jsapp._applet;
	  if (!applet) return null;
	  if (typeof(model) == "string")
		return applet.loadInlineString(model, "", false);
	  else
	    return applet.loadInlineArray(model, "", false);
};

WebMol.jmolModel = (function() {
	// class variables go here
	var defaultAtomStyle = {
		sphere : {},
		stick : null,
		line : null,
		cross : null,
		cartoon : null
	};

	defaultAtomStyle = {
		sphere : {}
	};


	function getCurJMolID(japp)
	{
		var modelInfo = Jmol.getPropertyAsArray(japp,"modelInfo");
		var cnt = modelInfo.modelCount;
		if(cnt == 0) { return "";}
		var models = modelInfo.models;
		var ret = models[modelInfo.modelCount-1].file_model;
		return ret;
	}

	//jmol color spec from hex value
	function jmolColor(c) {
		var hex = c.toString(16);
		//must have padding zeroes
		hex = "000000".substr(0, 6 - hex.length) + hex; 			  
		return "\"[x"+hex+"]\"";		
	}
	
	function jmolModel(japp, mid) {


		// private variables
		var atoms = [];
		var id = mid;
		var jmolid = null;

		this.getID = function() {
			return id;
		};

		// add atoms to this model from molecular data string
		this.addMolData = function(data, format) {
			Jmol.loadInline(japp, data);
			//figure out what model was just created
			jmolid = getCurJMolID();
		};

		//create an or statement if necessary from select= sel[i]
		//fn knows how to do the select
		function constructOrStatement(select, fn) {
			if($.isArray(select)) {
				var or = [];
				for(var j = 0; j < select.length; j++) {
					or.push(fn(select[j]));
				}
				return "("+or.join(" or ")+")";
			}
			else {
				return fn(select);
			}
		};
		
		// returns a jmol selection string (with the select keyword) for
		//the passed selection objection
		this.jmolSelect = function(sel) {
			var ret = ["model="+jmolid];
			for(var i in sel) {
				if(sel.hasOWnProperty(i)) {
					switch(i) {
					case "resn": //residue name
						if(typeof(sel[i]) != "undefined") {
							ret.push(constructOrStatement(sel[i], 
									function(x) {
										return "["+x+"]";
									}));
						}
						break;
					case "elem":
						if(typeof(sel[i]) != "undefined") {
							ret.push(constructOrStatement(sel[i], 
									function(x) {
										return "element=\""+x+"\"";
									}));
						}
						break;
					case "hetflag":
						if(typeof(sel[i]) != "undefined") {
							if(sel[i])
								ret.push("(hetero)");
							else
								ret.push("(not hetero)");
						}
						break;
					case "chain":
						if(typeof(sel[i]) != "undefined") {
							ret.push(constructOrStatement(sel[i], 
									function(x) {
										return ":"+x;
									}));
						}
						break;
					case "resi": //resid
						if(typeof(sel[i]) != "undefined") {
							ret.push(constructOrStatement(sel[i], 
									function(x) {
										return "resno="+x;
									}));
						}
						break;
					case "icode":
						if(typeof(sel[i]) != "undefined") {
							ret.push(constructOrStatement(sel[i], 
									function(x) {
										return "^"+x;
									}));
						}
						break;
					}
				}
			}
			return ret.join(" and ");
		}

		// style the select atoms with style
		this.setStyle = function(style, sel) {
			var select = "select " + jmolSelect(sel);
			var style = "";
			
			if(style.sphere) {
				style += "spacefill ";
				if(typeof(style.sphere.scale) != "undefined") {
					style += Math.round(style.sphere.scale*100)+"%";
				}
				else if(typeof(style.sphere.radius) != "undefined") {
					style += style.sphere.radius.toFixed(3);
				}
				style += ";";
				if(typeof(style.sphere.color) {
					style += "color "+jmolColor(style.sphere.color) + ";";
				}
			}
			else {
				style += "spacefill off;"; 
			}

			if(style.line || style.stick) {
				//ignore line styling of stick is set
				var c = null;
				if(style.stick) {
					var r = style.stick.radius || 0.25;
					style += "wireframe "+r + ";";
					if(style.stick.color) c = style.stick.color;
				}
				else {
					style += "wireframe;";
					if(style.line.color) c = style.line.color;
				}
				if(c != null) {
					style += "color wireframe "+jmolColor(c)+";";
				}
			}
			else {
				style += "wireframe off;"; 
			}
			
			if(style.cross) {
				style += "stars ";
				if(typeof(style.cross.scale) != "undefined") {
					style += Math.round(style.sphere.scale*100)+"%";
				}
				else if(typeof(style.cross.radius) != "undefined") {
					style += style.sphere.radius.toFixed(3);
				}
				style += ";";
				if(typeof(style.cross.color) {
					style += "color stars "+jmolColor(style.cross.color) + ";";
				}
			} else {
				style += "stars off";
			}
			
			if(style.cartoon) {
				style += "cartoon on;";
				if(style.cartoon.color) {
					style += "color cartoon "+jmolColor(style.cartoon.color) + ";";
				}
			} else {
				style += "cartoon off";
			}			
		};



	};

	return jmolModel;
})();