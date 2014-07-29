//  JSmolJME.js   Bob Hanson hansonr@stolaf.edu  6/14/2012 and 3/20/2013

// see http://peter-ertl.com/jsme/JSME_2013-10-13/api_javadoc/index.html


// BH 3/1/2014 4:31:18 PM fix for evaluate returning atom sets as arrays
// BH 1/27/2014 8:37:06 AM adding Info.viewSet  
// BH 12/4/2013 7:44:26 PM fix for JME independent search box

/*

	Only HTML5 version (JSME) is supported.

	JME 2D option -- use 

		Jmol.getJMEApplet(id, Info)
		Jmol.getJMEApplet(id, Info, linkedApplet)

	no option for getJMEAppletHtml(), but instead we indicate the
	target div using Info.divId.

	linkedApplet puts JME into INFO block for that applet; 
	use Jmol.showInfo(jme,true/false) to show/hide JME applet with id "jme"

	JME licensing: http://www.molinspiration.com/jme/doc/index.html
	note that required boilerplate: "JME Editor courtesy of Peter Ertl, Novartis"

	API includes:

	Jmol.jmeSmiles = function(jme, withStereoChemistry); 

		// returns SMILES string

	Jmol.jmeGetFile = function(jme, asJME)

		// retrieves structure as JME or MOL data

	Jmol.jmeReadMolecule = function(jme, jmeOrMolData); 

		// loads JME or MOL data into the app
		// JME data is recognized as a single line with no line ending


	Jmol.jmeReset = function(jme);

		// clears the app

	Jmol.jmeOptions = function(jme, options);

			 
	All other methods are private to JSmolJME.js




*/


;(function (Jmol, document) {

	Jmol._JMEApplet = function(id, Info, linkedApplet, checkOnly) {
		this._isJME = true;
		this._isJava = false;//(Info.use && (Info.use.toUpperCase() != "HTML5"))
		this._jmolType = "Jmol._JME" + (this._isJava ? "(JAVA)" : "(HTML5)");
		this._viewType = "JME";
		if (checkOnly)
			return this;
		window[id] = this;
		Jmol._setObject(this, id, Info);
		this._options = Info.options;
		(this._options.indexOf("autoez") < 0) && (this._options += ",autoez");
		this._setCheck(true);
		this._editEnabled = true;
		this._editOptions =  Info.editOptions || "editEnabled;";
		this._setEditOptions();
		this._editMol = "";
		this._highlightColor = Info.highlightColor || 1;
		if (this._viewSet != null) {
			this._options += ",star";
		}
		var jmol = this._linkedApplet = linkedApplet;
		this._hasOptions = Info.addSelectionOptions;
		this._readyFunction = Info.readyFunction;
		this._ready = false; 
 		this._jarPath = Info.jarPath;
		this._jarFile = Info.jarFile;
		if (jmol)
			this._console = jmol._console;
		Jmol._setConsoleDiv(this._console);
		this._isEmbedded = !Info.divId;
		this._divId = Info.divId || this._id + "_jmeappletdiv";
 		if (Jmol._document) {
			if (jmol) {
				jmol._2dapplet = this;
				var id = jmol._id + "_2dappletdiv";
				if (this._isEmbedded)
					this._divId = id;
				var d = Jmol._document;
				Jmol._document = null;
				Jmol.$html(this._divId, this.create());
				Jmol._document = d;
				this.__showContainer(false, false);
			} else {
				this.create();
			}
		}
		return this;
	}

	Jmol._JMEApplet._get = function(id, Info, linkedApplet, checkOnly) {

	// requires JmolJME.js and JME.jar
	// Jmol.getJMEApplet("jme", Info)
	// window.JME will be created as the return to this function

		Info || (Info = {});
		var DefaultInfo = {
			width: 300,
			height: 300,
			jarPath: "jme",
			jarFile: "JME.jar",
			use: "HTML5",
			structureChangedCallback: null, // could be myFunction(); first parameter will be reference to this object
			editOptions: "editEnabled",
			highlightColor: 1,  // 1-6
			options: "autoez"
			// see http://www2.chemie.uni-erlangen.de/services/fragment/editor/jme_functions.html
			// rbutton, norbutton - show / hide R button
			// hydrogens, nohydrogens - display / hide hydrogens
			// query, noquery - enable / disable query features
			// autoez, noautoez - automatic generation of SMILES with E,Z stereochemistry
			// nocanonize - SMILES canonicalization and detection of aromaticity supressed
			// nostereo - stereochemistry not considered when creating SMILES
			// reaction, noreaction - enable / disable reaction input
			// multipart - possibility to enter multipart structures
			// number - possibility to number (mark) atoms
			// depict - the applet will appear without editing butons,this is used for structure display only
		};		
		Jmol._addDefaultInfo(Info, DefaultInfo);
		var applet = new Jmol._JMEApplet(id, Info, linkedApplet, checkOnly);
		return (checkOnly ? applet : Jmol._registerApplet(id, applet));  
	}

	jsmeOnLoad = Jmol._JMEApplet.onload = function() {
		for (var i in Jmol._applets) {
			var app = Jmol._applets[i]
			if (app._isJME && !app._isJava && !app._ready) {
				app._applet = new JSApplet.JSME(app._divId, app.__Info);
				app._applet.options(app._options);
				var f = "";
				if (app._viewSet) {
					f = "(function(){" + app._id + "._myEditCallback()})";
				} else if (app.__Info.structureChangedCallback) {
					var m = app.__Info.structureChangedCallback;
					if (!m.endsWith(")"))
						m += "()";
					f = "(function(a,b,c,d){"+m.replace(/\(\)/, "(" + app._id +",a,b,c,d)") + "})";
				}
				if (f) 
					app._applet.setNotifyStructuralChangeJSfunction(f);
				app._ready = true;
				if (app._isEmbedded && app._linkedApplet._ready && app.__Info.visible)
					app._linkedApplet.show2d(true);
				Jmol._setReady(app);
				app._setSVG();
			}
		}
	}   

;(function(proto){

	proto.create = function() {
		var s = "";
		if (this._isJava) {
			var w = (this._linkedApplet ? "2px" : this._containerWidth);
			var h = (this._linkedApplet ? "2px" : this._containerHeight);
			s = '<applet code="JME.class" id="' + this._id + '_object" name="' + this._id 
				+ '_object" archive="' + this._jarFile + '" codebase="' + this._jarPath + '" width="'+w+'" height="'+h+'">'
				+ '<param name="options" value="' + this._options + '" />'	
				+ '</applet>';
		} else if (this._isEmbedded) {
			return this._code = "";
		}    
		if (this._hasOptions)
			s += Jmol._getGrabberOptions(this);      
		return this._code = Jmol._documentWrite(s);
	}

	proto._checkDeferred = function(script) {
		return false;
	}	

	proto._search = function(query){
		Jmol._search(this, query);
	}

	proto._searchDatabase = function(query, database, _jme_searchDatabase){
		this._showInfo(false);
		this._searchQuery = database + query;
		if (database == "$")
			query = "$" + query; // 2D variant;  will be $$caffeine
		var dm = database + query;
		this._loadFile(dm, {chemID: this._searchQuery});
	}

 	proto._loadFile = function(fileName, params, _jme_loadFile){
		var chemID = (params ? params.chemID : fileName);
		this._showInfo(false);
		this._thisJmolModel = "" + Math.random();
		var me = this;
		Jmol._loadFileData(this, fileName, function(data){me.__loadModel(data, chemID)}, function() {me.__loadModel(null, chemID)});
	}

	proto.__loadModel = function(jmeOrMolData, chemID, viewID, _jme__loadModel) {
		System.out.println("__loadModel")
		if (jmeOrMolData == null)
			return;
		Jmol.jmeReadMolecule(this, jmeOrMolData);
		if (this._viewSet != null)
			Jmol.View.updateView(this, {chemID:chemID, data:jmeOrMolData, viewID: viewID});      
	}

	proto._loadModelFromView = function(view, _jme_loadModelFromView) {
		// request from Jmol.View to update view with view.JME.data==null or needs changing
		var rec = view.JME;
		this._currentView = view;
		if (rec.data != null) {
			this.__loadModel(rec.data, view.info.chemID, view.info.viewID);
		} else if (rec.chemID != null) {
			Jmol._searchMol(this, view.info.chemID, null, false);
		} else {
			rec = view.Jmol;
			if (rec) {
				this._show2d(true, rec.applet);
			}
		}
		var me = this;
	}

	proto._updateView = function(_jme_updateView) {
		// called from model change without chemical identifier, possibly by user action and call to Jmol.updateView(applet)
		if (this._viewSet != null)
			this._search("$" + this._getSmiles())
		var me = this;
	}

	proto._setCheck = function(b, why) {
		this._checkEnabled = b;
	}
		 
	proto._setEditOptions = function(o) {
		o || (o = this._editOptions);
		o = o.toLowerCase();
		if (o.indexOf("edit") >= 0)
			this._editEnabled = (o.indexOf("editdisable") < 0);
	}
	
	proto._setSVG = function() {
		var me = this;
		Jmol._$(this._divId).find("div").each(function(){this._applet = me})
		this._starPressed = false;
	}
	
	proto._enableEdit = function(tf) {
	  this._editEnabled = tf;
		this._setSVG();
	}

	Jmol._jmeHook = function(a, e) {
		// called from within function C(a) in JSME code
		// a is the function to send this event to using apply
		var d;
		return (!(d = e.target) || !(d._applet || (d = Jmol.$closest(d, 'div')[0])) || !d || !d._applet ? a : d._applet.__allowEvent(a, e));
	}
	
	proto.__allowEvent = function(a, e){
		var t;
		if(!e || !(t = e.target) || "DIV rect line text polygon ellipse path".indexOf("" + t.tagName) < 0){
			// cancel only these specific tags
			return a;
		}
		
		var isKey = (e.type.indexOf("key") >= 0);
		if  (!this._editEnabled && isKey) {
			// cancel all keyboard events
			return 0;
		}
		if ("dblclick mousedown mouseup".indexOf(e.type) < 0) {
			// cancel only these specific events
			return a;
		}	
		var x = (t.textContent ? t.x.baseVal[0].value : t.points ? t.animatedPoints[0].x : isKey || !t.x && !t.x1 ? 200 : (t.x || t.x1).baseVal.value);
		var y = (t.textContent ? t.y.baseVal[0].value : t.points ? t.animatedPoints[0].y : isKey || !t.x && !t.x1 ? 200 : (t.y || t.y1).baseVal.value);
		// when editing is disabled, only the star key and main-panel clicking will be allowed
		var isStar = (x >= 100 && x < 124 && y < 24); // fifth icon from the left on top row
		var isMain = (x > 25 && y > 50 || x == 0 && y == 0 && t.width.baseVal.value > 100 && t.height.baseVal.value > 100);
		if (isStar || this._editEnabled && !isMain)
			this._starPressed = isStar;			
		var ok = (isStar || this._editEnabled || isMain && this._starPressed);
		return (ok ? a : 0);
	}
			
	proto._myEditCallback = function(_jme_myEditCallback) {
		// direct callback from JSME applet
		var data = this._applet.jmeFile().replace(/\:1/g,"");
		System.out.println("myEditCallback " + data)
		this._editMol = this._editMol.replace(/\:1/g,"");
		if (this._checkEnabled && !this._editEnabled) {
			// data is not null, and we don't allow editing
			// data is null, and we don't allow clearing
	 		this._editEnabled && data && (this._editMol = data); /// was .editEnabled
			if (data && data != this._editMol) {
				var me = this;
				var m = me._editMol;
				setTimeout(function(){
					me._setCheck(false, "sorry");
					(m ? me._applet.readMolecule(m) : me._applet.reset());
					me._setSVG();
					me._editMol = me._applet.jmeFile();
				},150);
				return;
			}
		}
		var me = this;
		setTimeout(function() {me._setCheck(true, "callback")}, 10);
		data && (this._editMol = data);
		if (this._viewSet == null)
			return;
		var a = this._applet.molFile().split("V2000")[1];
		var b = ("" + this._molData).split("V2000")[1];
		if (a != b) {
			this._molData = "<modified>";
			this._thisJmolModel = null;
			return this.__clearAtomSelection(false);
		}
		// not a structural change
		var data = this._applet.jmeFile();
		if (data.indexOf(":") < 0) 
			return this.__clearAtomSelection(true);
		data = data.split(" ");
		var n = parseInt(data[0]);
		var iAtom = 0;
		for (var i = 0; i < n; i++)
			if (!this.__atomSelection[i + 1] && data[i*3 + 2].indexOf(":") >= 0) {
				iAtom = i + 1;
				break;
			}
		if (iAtom <= 0)
			return;
		var A = [];
		var map = this._currentView.JME.atomMap;    
		A.push(map == null ? iAtom : map.toJmol[iAtom]);
		Jmol.View.updateAtomPick(this, A);
		this._updateAtomPick(A);
		if (this._atomPickCallback)
			setTimeout(this._atomPickCallback+"([" + iAtom + "])",10);    
	}

	proto.__clearAtomSelection = function(andUpdate) {
		System.out.println("clearAtomSelection");
		this.__atomSelection = [];
		this._applet.resetAtomColors(1);
		if (andUpdate)
			Jmol.View.updateAtomPick(this, []);
	}	

	proto._updateAtomPick = function(A, _jme_updateAtomPick) {
		this._applet.resetAtomColors(1);
		if (A.length == 0)
			return;
		var B = [];
		var C = [];
		var j;		
		var map = this._currentView.JME.atomMap;		
		for (var i = 0; i < A.length; i++) {
		 C[j = (map == null ? A[i] : map.fromJmol[A[i]])] = 1; 
		 B.push(j);
		 B.push(this._highlightColor);
		}
		System.out.println("JME setting atom colors " + B.join(","))
		this._applet.setAtomBackgroundColors(1, B.join(","));
		this.__atomSelection = C;
	}

	proto._showInfo = Jmol._Applet.prototype._showInfo;
	proto._show = Jmol._Applet.prototype._show;

	proto._show = function(tf) {
		var x = (!tf ? 2 : "100%");
		Jmol.$setSize(Jmol.$(this, "object"), x, x);
		if (!this._isJava) {
			Jmol.$setVisible(Jmol.$(this, "appletdiv"), tf);
			if (this._isEmbedded && !tf)
				Jmol.$setVisible(Jmol.$(this._linkedApplet, "2dappletdiv"), false);
		}
	}

	proto._show2d = function(toJME, jmol) {
		jmol || (jmol = this._linkedApplet);
		if (jmol) {
			var jme = this._applet;
			if (jme == null && this._isJava)
				jme = this._applet = Jmol._getElement(this, "object");
			var isOK = true;
			if (this._viewSet != null) {
			  isOK = false;
			} else if (jme != null) {
				var jmeSMILES = this._getSmiles();
				// testing here to see that we have the same structure as in the JMOL applet
				// feature change here --- evaluation of an atom set returns an array now, not an uninterpretable string
				var jmolAtoms = (jmeSMILES ? jmol._evaluate("{*}.find('SMILES', '" + jmeSMILES.replace(/\\/g,"\\\\")+ "')") : []);
				var isOK = (jmolAtoms.length > 0);
			}
			if (!isOK) {
				if (toJME) {
				  this._loadFromJmol(jmol);
				} else {
					// toJmol
					if (jmeSMILES)
						Jmol.script(jmol, "load \"$" + jmeSMILES + "\"");
				}
			}
		}
		if (this._linkedApplet) {
		 	this.__showContainer(toJME, true);
			this._showInfo(!toJME);
		}
	}
	
	proto._loadFromJmol = function(jmol) {
		this._molData = jmol._getMol2D();
		setTimeout(this._id + ".__readMolData()",10);
 }

	proto.__readMolData = function() {
		if (!this._applet)return;
		this._setCheck(false, "readmoldata");
		if (this._molData) {
			this._applet.readMolFile(this._molData);
			this._molData = this._applet.molFile();
			if (this._viewSet) {
			  var v = this._currentView;
			  v.JME.data = this._molData;
				v.JME.atomMap = (v.Jmol && v.Jmol.applet? v.Jmol.applet._getAtomCorrelation(this._molData) : null);
			}
		} else {
			this._applet.reset();
			this._molData = "<zapped>";
		}
		this._editMol = this._applet.jmeFile();
		this._setSVG();
	}
  
	proto.__showContainer = function(tf, andShow) {
		var jmol = this._linkedApplet;
		var mydiv = Jmol.$(jmol, "2dappletdiv");
		if (this._isJava) {
			var w = (tf ? "100%" : 2);
			var h = (tf ? "100%" : 2);
			Jmol.$setSize(mydiv, w, h);
			if (andShow)
				Jmol.$setSize(Jmol.$(this, "object"), w, h);
		} else {
			Jmol.$setVisible(mydiv, tf);
		}
	}

	proto._script = function(script) {
	 // a little scripting language for JME!
	 // print 
		var list = script.split(";");
		for (i = 0; i < list.length; i++) {
			switch (list[i].split(" ")[0].trim().toLowerCase()) {
			case "print":
				// only FF and Chrome right now. Maybe Safari
				var svg = Jmol.$("#" + this._divId).find("svg")[0];
				var img = new Image();
				var me = this;
				img.onload = function() {
					var canvas = document.createElement("canvas");
					var ctx = canvas.getContext('2d');
					ctx.clearRect( 0, 0, (canvas.width = svg.width.animVal.value - 5), (canvas.height = svg.height.animVal.value));
					ctx.drawImage(img, 0, 0);
					// throw out "data:image/png;base64," because we will reconstruct that if we need to, and we might not
					Jmol._saveFile(me._id + ".png", canvas.toDataURL("image/png").substring(22), "image/png", "base64");
				}
				img.src = "data:image/svg+xml;base64," + btoa(svg.outerHTML);
				break;
			}
		}
	}

  proto._getSmiles = function(withStereoChemistry) {
  	var s = (arguments.length == 0 || withStereoChemistry ? jme._applet.smiles() : jme._applet.nonisomericSmiles()).replace(/\:1/g,"");
		s = s.replace(/H/g,"").replace(/\[\]/g,"").replace(/@\]/g,"@H]").replace(/\(\)/g,"");
		if (s.indexOf("\\") == 0 || s.indexOf("/") == 0)
		  s= "[H]" + s;
		return s;
  }

  proto._getMol = function() {
		return this._applet.molFile();   
  }
  
})(Jmol._JMEApplet.prototype);

	//////  additional API for JME /////////

	// see also http://www2.chemie.uni-erlangen.de/services/fragment/editor/jme_functions.html

	// The final replacement here is to remove markings from star option.

	Jmol.jmeSmiles = function(jme, withStereoChemistry) {
		return jme._getSmiles();
	}

	Jmol.jmeReadMolecule = function(jme, jmeOrMolData) {
		// JME data is a single line with no line ending
		jme._setCheck(false, "readmolecule");
		if (jmeOrMolData.indexOf("\n") < 0 && jmeOrMolData.indexOf("\r") < 0)
			jme._applet.readMolecule(jmeOrMolData);
		else 
			jme._applet.readMolFile(jmeOrMolData);   
	 jme._molData = jme._applet.molFile();
	 jme._editMol = jme._applet.jmeFile();
	 jme._setSVG();
	}

	Jmol.jmeGetFile = function(jme, asJME) {
		jme._molData = jme._applet.molFile();
		return  (asJME ? jme._applet.jmeFile() : jme._molData);
	}

	Jmol.jmeReset = function(jme) {
		jme._applet.reset();
	}

	Jmol.jmeOptions = function(jme, options) {
		jme._applet.options(options);
	}

// doesn't work because of the way JSME is created using frames and SVG.
//	
//  Jmol.getJSVAppletHtml = function(applet, Info, linkedApplet) {
//    if (Info) {
//      var d = Jmol._document;
//      Jmol._document = null;
//      applet = Jmol.getJMEApplet(applet, Info, linkedApplet);
//      Jmol._document = d;
//    }  
//    return applet._code;
//	}


})(Jmol, document);

