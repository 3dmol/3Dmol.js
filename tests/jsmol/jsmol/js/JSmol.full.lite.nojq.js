// JSmoljQueryExt.js
// 9/2/2013 7:43:12 AM BH Opera/Safari fix for binary file reading
// 3/11/2014 6:31:01 AM BH fix for MSIE not working locally

;(function($) {

	function createXHR(isMSIE) {
		try {
			return (isMSIE ? new window.ActiveXObject( "Microsoft.XMLHTTP" ) : new window.XMLHttpRequest());
		} catch( e ) {}
	}

 $.ajaxSettings.xhr = (window.ActiveXObject === undefined ? createXHR :  
	function() {
		return (this.url == document.location || this.url.indexOf("http") == 0 || !this.isLocal) &&  // BH MSIE fix
			/^(get|post|head|put|delete|options)$/i.test( this.type ) &&
			createXHR() || createXHR(1);
	});
 
	// incorporates jquery.iecors MSIE asynchronous cross-domain request for MSIE < 10

	$.extend( $.support, { iecors: !!window.XDomainRequest });

	if ($.support.iecors) {
		// source: https://github.com/dkastner/jquery.iecors
		// author: Derek Kastner dkastner@gmail.com http://dkastner.github.com    
		$.ajaxTransport(function(s) {
			return {
				send: function( headers, complete ) {
					// Note that xdr is not synchronous.
					// This is only being used in JSmol for transport of java code packages.
					var xdr = new window.XDomainRequest();
					xdr.onload = function() {          
						var headers = { 'Content-Type': xdr.contentType };
						complete(200, 'OK', { text: xdr.responseText }, headers);
					};
					if ( s.xhrFields ) {
						xdr.onerror = s.xhrFields.error;
						xdr.ontimeout = s.xhrFields.timeout;
					}
					xdr.open( s.type, s.url );
					xdr.send( ( s.hasContent && s.data ) || null );
				},
				abort: function() {        
					xdr.abort();
				}
			};
		});

	} else {

	// adds support for synchronous binary file reading

		$.ajaxSetup({
			accepts: { binary: "text/plain; charset=x-user-defined" },
			responseFields: { binary: "response" }
		})

		$.ajaxTransport('binary', function(s) {
			var callback;
			return {
				// synchronous or asynchronous binary transfer only
				send: function( headers, complete ) {        
					var xhr = s.xhr();
					xhr.open( s.type, s.url, s.async );					
					var isOK = false;
					try {
						if (xhr.hasOwnProperty("responseType")) {
								xhr.responseType = "arraybuffer";
								isOK = true;
						} 
					} catch(e) {
					  //
					}
					try {
						if (!isOK && xhr.overrideMimeType) {
							xhr.overrideMimeType('text/plain; charset=x-user-defined');
						}
					} catch(e) {
							//
					}
					if ( !s.crossDomain && !headers["X-Requested-With"] ) {
						headers["X-Requested-With"] = "XMLHttpRequest";
					}
					try {
						for (var i in headers )
							xhr.setRequestHeader( i, headers[ i ] );
					} catch(_) {}

					xhr.send( ( s.hasContent && s.data ) || null );

 					// Listener
					callback = function( _, isAbort ) {

					var 
						status = xhr.status,
						statusText = "",
						responseHeaders = xhr.getAllResponseHeaders(),
						responses = {},
						xml;

					try {

						// Firefox throws exceptions when accessing properties
						// of an xhr when a network error occured
						// http://helpful.knobs-dials.com/index.php/Component_returned_failure_code:_0x80040111_(NS_ERROR_NOT_AVAILABLE)
						// Was never called and is aborted or complete
						if ( callback && ( xhr.readyState === 4 ) ) {

							// Only called once
							callback = undefined;

							// When requesting binary data, IE6-9 will throw an exception
							// on any attempt to access responseText (#11426)
							try {
								responses.text = (typeof xhr.responseText === "string" ? xhr.responseText : null);
							} catch( _ ) {
							}
							try {
								responses.binary = xhr.response;
							} catch( _ ) {
							}

							// Firefox throws an exception when accessing
							// statusText for faulty cross-domain requests
							try {
								statusText = xhr.statusText;
							} catch( _ ) {
								// We normalize with Webkit giving an empty statusText
								statusText = "";
							}
							// Filter status for non standard behaviors

							// If the request is local and we have data: assume a success
							// (success with no data won't get notified, that's the best we
							// can do given current implementations)
							if ( !status && s.isLocal && !s.crossDomain ) {
								status = (responses.text ? 200 : 404);
							// IE - #1450: sometimes returns 1223 when it should be 204
							} else if ( status === 1223 ) {
								status = 204;
							}
							complete( status, statusText, responses, responseHeaders );
						}
					} catch( e ) {
						alert(e)
						complete( -1, e );
					}
					};
					
					if ( !s.async ) {
						// if we're in sync mode we fire the callback
						callback();
					} else if ( xhr.readyState === 4 ) {
						// (IE6 & IE7) if it's in cache and has been
						// retrieved directly we need to fire the callback
						setTimeout( callback );
					} else {
						// Add to the list of active xhr callbacks
						xhr.onreadystatechange = callback;
					}
					
				},
				abort: function() {}
			};
		});
	}
})( jQuery );
	 
/*
 * jQuery outside events - v1.1 - 3/16/2010
 * http://benalman.com/projects/jquery-outside-events-plugin/
 * 
 * Copyright (c) 2010 "Cowboy" Ben Alman
 * Dual licensed under the MIT and GPL licenses.
 * http://benalman.com/about/license/
 * 
 * Modified by Bob Hanson for JSmol-specific events and to add parameter reference to actual jQuery event.
 * Used for closing the pop-up menu.
 *   
 */

;(function($,doc,eventList,id){  
	// was 'click dblclick mousemove mousedown mouseup mouseover mouseout change select submit keydown keypress keyup'
	$.map(
		eventList.split(' '),
		function( event_name ) { jq_addOutsideEvent( event_name ); }
	);
	jq_addOutsideEvent( 'focusin',  'focus' + id );
	jq_addOutsideEvent( 'focusout', 'blur' + id );
	function jq_addOutsideEvent( event_name, outside_event_name ) {
		outside_event_name = outside_event_name || event_name + id;
		var elems = $(),
			event_namespaced = event_name + '.' + outside_event_name + '-special-event';
		$.event.special[ outside_event_name ] = {    
			setup: function(){
				elems = elems.add( this );
				if ( elems.length === 1 ) {
					$(doc).bind( event_namespaced, handle_event );
				}
			},
			teardown: function(){
				self.Jmol && Jmol._setMouseOwner(null);
				elems = elems.not( this );
				if ( elems.length === 0 ) {
					$(doc).unbind( event_namespaced );
				}
			},
			add: function( handleObj ) {
				var old_handler = handleObj.handler;
				handleObj.handler = function( event, elem ) {
					event.target = elem;
					old_handler.apply( this, arguments );
				};
			}
		};
		function handle_event( event ) {
			$(elems).each(function(){
				self.Jmol && (outside_event_name.indexOf("mouseup") >= 0 || outside_event_name.indexOf("touchend") >= 0) && Jmol._setMouseOwner(null);
				var elem = $(this);
				if ( this !== event.target && !elem.has(event.target).length ) {
					//BH: adds event to pass that along to our handler as well.
					elem.triggerHandler( outside_event_name, [ event.target, event ] );
				}
			});
		};
	};
})(jQuery,document,"click mousemove mouseup touchmove touchend", "outjsmol");
// JSmolCore.js -- Jmol core capability  4/27/2014 6:32:42 PM

// see JSmolApi.js for public user-interface. All these are private functions

// BH 5/30/2014 7:20:07 AM better dragging for console and menu
// BH 4/27/2014 6:31:52 PM allows _USE=SIGNED HTML5 as well as _USE=JAVA HTML5
// BH 3/8/2014 5:50:51 PM adds support for dataURI download in FF and Chrome
// BH 3/8/2014 8:43:10 AM moves PubChem access to https
// BH 3/4/2014 8:40:15 PM adds Jmol.Cache for JSV/Jmol sharing files
// BH 2/10/2014 10:07:14 AM added Info.z and Info.zIndexBase
// BH 2/9/2014 9:56:06 PM updated JSmolCore.js with option to extend Viewer with code PRIOR to loading Viewer classes
// BH 2/6/2014 8:46:25 AM disabled Jmol._tracker for localhost and 127.x.x.x 
// BH 1/29/2014 8:02:23 AM Jmol.View and Info.viewSet
// BH 1/21/2014 12:06:59 PM adding Jmol.Info.cacheFiles (applet, true/false) and applet._cacheFiles and Jmol._fileCache
// BH 1/13/2014 2:12:38 PM adding "http://www.nmrdb.org/tools/jmol/predict.php":"%URL", to _DirectDatabaseCalls
// BH 12/21/2013 6:38:35 PM applet sync broken
// BH 12/6/2013 6:18:32 PM cover.htm and coverImage fix
// BH 12/4/2013 7:44:26 PM fix for JME independent search box
// BH 12/3/2013 6:30:08 AM fix for ready function returning Boolean instead of boolean in HTML5 version
// BH 11/30/2013 10:31:37 AM added type:"GET" for jQuery.ajax() requests instead of using defaults
// BH 11/30/2013 10:31:37 AM added cache:true for jQuery.ajax() requests; can override with cache:"NO", not cache:false
// BH 11/28/2013 11:09:27 AM added Jmol._alertNoBinary:true
// BH 11/26/2013 8:19:55 PM fix !xxxx search commmand entry and stop MSIE from duplicating command
// BH 11/25/2013 7:38:31 AM adds Jmol._tracker: option for GoogleAnalytics tracking
// BH 11/25/2013 7:39:03 AM adds URL options _J2S=  _JAR=  _USE=
// BH 11/23/2013 10:51:37 PM  adds JNLP support for local applet
// BH 11/2/2013 12:05:11 PM JSmolJSME fixes; https access fixed
// BH 10/31/2013 7:50:06 PM Jmol.Dialog as SwingController; Jmol._mouseOwner added
// BH 10/19/2013 7:05:04 AM adding Jmol._ajaxCall for Error Contacting Server; database POST method enabled
// BH 10/17/2013 1:40:51 PM  adding javajs/swing and Jmol.Dialog
// BH 9/30/2013 6:42:24 PM: pdb.gz switch  pdb should only be for www.rcsb.org
// BH 9/17/2013 10:17:51 AM: asynchronous file reading and saving
// BH 8/16/2013 12:02:20 PM: JSmoljQueryExt.js pulled out
// BH 8/16/2013 12:02:20 PM: Jmol._touching used properly

// BH 3/22/2013 5:53:02 PM: Adds noscript option, JSmol.min.core.js
// BH 1/17/2013 5:20:44 PM: Fixed problem with console not getting initial position if no first click
// 1/13/2013 BH: Fixed MSIE not-reading-local-files problem.
// 11/28/2012 BH: Fixed MacOS Safari binary ArrayBuffer problem
// 11/21/2012 BH: restructuring of files as JS... instead of J...
// 11/20/2012 BH: MSIE9 cannot do a synchronous file load cross-domain. See Jmol._getFileData
// 11/4/2012 BH: RCSB REST format change "<structureId>" to "<dimStructure.structureId>"
// 9/13/2012 BH: JmolCore.js changes for JSmol doAjax() method -- _3ata()
// 6/12/2012 BH: JmolApi.js: adds Jmol.setInfo(applet, info, isShown) -- third parameter optional 
// 6/12/2012 BH: JmolApi.js: adds Jmol.getInfo(applet) 
// 6/12/2012 BH: JmolApplet.js: Fixes for MSIE 8
// 6/5/2012  BH: fixes problem with Jmol "javascript" command not working and getPropertyAsArray not working
// 6/4/2012  BH: corrects problem with MSIE requiring mouse-hover to activate applet
// 5/31/2012 BH: added JSpecView interface and api -- see JmolJSV.js
//               also changed "jmolJarPath" to just "jarPath"
//               jmolJarFile->jarFile, jmolIsSigned->isSigned, jmolReadyFunction->readyFunction
//               also corrects a double-loading issue
// 5/14/2012 BH: added AJAX queue for ChemDoodle option with multiple canvases 
// 8/12/2012 BH: adds support for MSIE xdr cross-domain request (jQuery.iecors.js)

// allows Jmol applets to be created on a page with more flexibility and extendability
// provides an object-oriented interface for JSpecView and syncing of Jmol/JSpecView


// required/optional libraries (preferably in the following order):

//    jquery/jquery.js     -- at least jQuery.1.9
//    js/JSmoljQueryext.js -- required for binary file transfer; otherwise standard jQuery should be OK
//    js/JSmolCore.js      -- required
//    js/j2sjmol.js        -- required
//    js/JSmol.js          -- required
//    js/JSmolApplet.js    -- required; internal functions for _Applet and _Image; must be after JSmolCore
//    js/JSmolControls.js  -- optional; internal functions for buttons, links, menus, etc.; must be after JSmolCore
//    js/JSmolConsole.js   -- optional; for the pop-up console
//    js/JSmolApi.js       -- required; all user functions; must be after JSmolCore
//    js/JSmolTHREE.js     -- optional; WebGL library required for JSmolGLmol.js
//    js/JSmolGLmol.js     -- optional; WebGL version of JSmol.
//    js/JSmolJME.js       -- optional; JSME (2D editor)
//    jsme/jsme/jsme.nocache.js   --  required for JSME 
//    js/JSmolMenu.js      -- optional; required for menuing in JSV
//    js/JSmolJSV.js       -- optional; for creating and interacting with a JSpecView applet 

// most of these will be loaded automatically, and for most installations, all you need is JSmol.min.js


// Allows Jmol-like objects to be displayed on Java-challenged (iPad/iPhone)
// or applet-challenged (Android/iPhone) platforms, with automatic switching to 

// For your installation, you should consider putting JmolData.jar and jsmol.php 
// on your own server. Nothing more than these two files is needed on the server, and this 
// allows more options for MSIE and Chrome when working with cross-domain files (such as RCSB or pubChem) 

// The NCI and RCSB databases are accessed via direct AJAX if available (xhr2/xdr).


if(typeof(jQuery)=="undefined") alert ("Note -- JSmoljQuery is required for JSmol, but it's not defined.")

// An example of how to extend Jmol with code PRIOR to JSmolCore.js or JSmol.min.js follows:
//
// 
//	Jmol = {
//  	z:3000,
//		extend: function(what, obj) {if (what == "viewer") { obj._testing = true } }
//	}

self.Jmol || (Jmol = {});

if (!Jmol._version)
Jmol = (function(document) {
	var z=Jmol.z || 9000;
	var getZOrders = function(z) {
		return {
			header:z++,
			rear:z++,
			main:z++,
			image:z++,
			front:z++,
			fileOpener:z++,
			coverImage:z++,
			dialog:z++, // could be several of these, JSV only
			menu:z+90000, // way front
			console:z+91000, // even more front
			monitorZIndex:z+99999 // way way front
		}
	};
	var j = {
		_version: 'JSmol 14.1.14 Apr 27, 2014',
		_alertNoBinary: true,
		// this url is used to Google Analytics tracking of Jmol use. You may remove it or modify it if you wish. 
		_allowedJmolSize: [25, 2048, 300],   // min, max, default (pixels)
		/*  By setting the Jmol.allowedJmolSize[] variable in the webpage
				before calling Jmol.getApplet(), limits for applet size can be overriden.
				2048 standard for GeoWall (http://geowall.geo.lsa.umich.edu/home.html)
		*/
		_fileCache: null, // enabled by Jmol.setFileCaching(applet, true/false)
		_jarFile: null,  // can be set in URL using _JAR=
		_j2sPath: null,  // can be set in URL using _J2S=
		_use: null,      // can be set in URL using _USE=
		_j2sLoadMonitorOpacity: 90, // initial opacity for j2s load monitor message
		_applets: {},
		_asynchronous: true,
		_ajaxQueue: [],
		_getZOrders: getZOrders,
		_z:getZOrders(z),
		_debugCode: true,  // set false in process of minimization
		db: {
			_databasePrefixes: "$=:",
			_fileLoadScript: ";if (_loadScript = '' && defaultLoadScript == '' && _filetype == 'Pdb') { select protein or nucleic;cartoons Only;color structure; select * };",
			_nciLoadScript: ";n = ({molecule=1}.length < {molecule=2}.length ? 2 : 1); select molecule=n;display selected;center selected;",
			_pubChemLoadScript: "",
			_DirectDatabaseCalls:{
				"cactus.nci.nih.gov": "%URL",
				"www.rcsb.org": "%URL",
				"pubchem.ncbi.nlm.nih.gov":"%URL",
				"http://www.nmrdb.org/tools/jmol/predict.php":"%URL",
				"$": "http://cactus.nci.nih.gov/chemical/structure/%FILENCI/file?format=sdf&get3d=True",
				"$$": "http://cactus.nci.nih.gov/chemical/structure/%FILENCI/file?format=sdf",
				"=": "http://www.rcsb.org/pdb/files/%FILE.pdb",
				"==": "http://www.rcsb.org/pdb/files/ligand/%FILE.cif",
				":": "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/%FILE/SDF?record_type=3d"
			},
			_restQueryUrl: "http://www.rcsb.org/pdb/rest/search",
			_restQueryXml: "<orgPdbQuery><queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType><description>Text Search</description><keywords>QUERY</keywords></orgPdbQuery>",
			_restReportUrl: "http://www.pdb.org/pdb/rest/customReport?pdbids=IDLIST&customReportColumns=structureId,structureTitle"
		},
		_debugAlert: false,
		_document: document,
		_isXHTML: false,
		_lastAppletID: null,
		_mousePageX: null,
		_mouseOwner: null,
		_serverUrl: "http://your.server.here/jsmol.php",
		_syncId: ("" + Math.random()).substring(3),
		_touching: false,
		_XhtmlElement: null,
		_XhtmlAppendChild: false
	}
	
	var ref = document.location.href.toLowerCase();
	j._httpProto = (ref.indexOf("https") == 0 ? "https://" : "http://"); 
	j._isFile = (ref.indexOf("file:") == 0);
	j._ajaxTestSite = j._httpProto + "google.com";
	var isLocal = (j._isFile || ref.indexOf("http://localhost") == 0 || ref.indexOf("http://127.") == 0);
	j._tracker = (j._httpProto == "http://" && !isLocal && 'http://chemapps.stolaf.edu/jmol/JmolTracker.htm?id=UA-45940799-1');
	
	j._isChrome = (navigator.userAgent.toLowerCase().indexOf("chrome") >= 0);
	j._isSafari = (!j._isChrome && navigator.userAgent.toLowerCase().indexOf("safari") >= 0);
	j._isMsie = (window.ActiveXObject !== undefined);
	j._useDataURI = !j._isSafari && !j._isMsie; // safari may be OK here -- untested

	for(var i in Jmol) j[i] = Jmol[i]; // allows pre-definition
	return j;
})(document, Jmol);


(function (Jmol, $) {

// this library is organized into the following sections:

	// jQuery interface 
	// protected variables
	// feature detection
	// AJAX-related core functionality
	// applet start-up functionality
	// misc core functionality
	// mouse events


	////////////////////// jQuery interface ///////////////////////

	// hooks to jQuery -- if you have a different AJAX tool, feel free to adapt.
	// There should be no other references to jQuery in all the JSmol libraries.

	Jmol.$ = function(objectOrId, subdiv) {
		// if a subdivv, then return $("#objectOrId_subdiv") 
		// or if no subdiv, then just $(objectOrId)
		if (objectOrId == null)alert (subdiv + arguments.callee.caller.toString());
			return $(subdiv ? "#" + objectOrId._id + "_" + subdiv : objectOrId);
	} 

	Jmol._$ = function(id) {
		// either the object or $("#" + id)
		return (typeof id == "string" ? $("#" + id) : id);
	}

	/// special functions:

	Jmol.$ajax = function (info) {
		Jmol._ajaxCall = info.url;
		info.cache = (info.cache != "NO");
		if (info.url.indexOf("http://pubchem.ncbi.nlm.nih") == 0)
			info.url = "https://" + info.url.substring(7);
			// fluke at pubchem --- requires https now from all pages 3/8/2014
		// don't let jQuery add $_=... to URL unless we 
		// use cache:"NO"; other packages might use $.ajaxSetup() to set this to cache:false
		return $.ajax(info);
	}

	Jmol._getNCIInfo = function(identifier, what, fCallback) {
		if (what == "name")
			what = "names"
		url = "http://cactus.nci.nih.gov/chemical/structure/"+identifier +"/" + what; 
		return Jmol._getFileData(url);
	}
	


	Jmol.$appEvent = function(app, subdiv, evt, f) {
		var o = Jmol.$(app, subdiv); 
		o.off(evt) && f && o.on(evt, f);
	}   

	Jmol.$resize = function (f) {
		return $(window).resize(f);
	}

	//// full identifier expected (could be "body", for example):

	Jmol.$after = function (what, s) {
		return $(what).after(s);
	}

	Jmol.$bind = function(what, list, f) {
		return (f ? $(what).bind(list, f) : $(what).unbind(list));
	}

	Jmol.$closest = function(what, d) {
		return $(what).closest(d);
	}
	
	Jmol.$get = function(what, i) {
	return $(what).get(i);
	}
 
	// element id expected
			 
	Jmol.$documentOff = function(evt, id) {
		return $(document).off(evt, "#" + id);
	}

	Jmol.$documentOn = function(evt, id, f) {
		return $(document).on(evt, "#" + id, f);
	}

	Jmol.$getAncestorDiv = function(id, className) {
		return $("div." + className + ":has(#" + id + ")")[0];
	}

	Jmol.$supportsIECrossDomainScripting = function() {
		return $.support.iecors;
	}

	//// element id or jQuery object expected:

	Jmol.$attr = function (id, a, val) {
		return Jmol._$(id).attr(a, val);
	}

	Jmol.$css = function(id, style) {
		return Jmol._$(id).css(style);
	}
	 
	Jmol.$find = function(id, d) {
		return Jmol._$(id).find(d);
	}
	
	Jmol.$focus = function(id) {
		return Jmol._$(id).focus();
	}

	Jmol.$html = function(id, html) {
		return Jmol._$(id).html(html);
	}
	 
	Jmol.$offset = function(id) {
		return Jmol._$(id).offset();
	}

	Jmol.$windowOn = function(evt, f) {
		return $(window).on(evt, f);
	}

	Jmol.$prop = function(id, p, val) {
		var o = Jmol._$(id);
		return (arguments.length == 3 ? o.prop(p, val) : o.prop(p));
	}

	Jmol.$remove = function(id) {
		return Jmol._$(id).remove();
	}

	Jmol.$scrollTo = function (id, n) {
		var o = Jmol._$(id);
		return o.scrollTop(n < 0 ? o[0].scrollHeight : n);
	}

	Jmol.$setEnabled = function(id, b) {
		return Jmol._$(id).attr("disabled", b ? null : "disabled");  
	}

	Jmol.$setSize = function(id, w, h) {
		return Jmol._$(id).width(w).height(h);
	}

	Jmol.$setVisible = function(id, b) {
		var o = Jmol._$(id);
		return (b ? o.show() : o.hide());  
	}

	Jmol.$submit = function(id) {
		return Jmol._$(id).submit();
	}

	Jmol.$val = function (id, v) {
		var o = Jmol._$(id);
		return (arguments.length == 1 ? o.val() : o.val(v));
	}

	////////////// protected variables ///////////


	Jmol._clearVars = function() {
		// only on page closing -- appears to improve garbage collection

		delete jQuery;
		delete $;
		delete Jmol;
		delete SwingController;
		delete J;
		delete JM;
		delete JS;
		delete JSV;
		delete JU;
		delete JV;
		delete java;
		delete javajs;
		delete Clazz;
		delete c$; // used in p0p; could be gotten rid of
	}

	////////////// feature detection ///////////////

	Jmol.featureDetection = (function(document, window) {

		var features = {};
		features.ua = navigator.userAgent.toLowerCase()

		features.os = (function(){
			var osList = ["linux","unix","mac","win"]
			var i = osList.length;

			while (i--){
				if (features.ua.indexOf(osList[i])!=-1) return osList[i]
			}
			return "unknown";
		})();

		features.browser = function(){
			var ua = features.ua;
			var browserList = ["konqueror","webkit","omniweb","opera","webtv","icab","msie","mozilla"];
			for (var i = 0; i < browserList.length; i++)
			if (ua.indexOf(browserList[i])>=0) 
				return browserList[i];
			return "unknown";
		}
		features.browserName = features.browser();
		features.browserVersion= parseFloat(features.ua.substring(features.ua.indexOf(features.browserName)+features.browserName.length+1));
		features.supportsXhr2 = function() {return ($.support.cors || $.support.iecors)}
		features.allowDestroy = (features.browserName != "msie");
		features.allowHTML5 = (features.browserName != "msie" || navigator.appVersion.indexOf("MSIE 8") < 0);
		features.getDefaultLanguage = function() {
			return navigator.language || navigator.userLanguage || "en-US";
		};

		features._webGLtest = 0;

		features.supportsWebGL = function() {
		if (!Jmol.featureDetection._webGLtest) { 
			var canvas;
			Jmol.featureDetection._webGLtest = ( 
				window.WebGLRenderingContext 
					&& ((canvas = document.createElement("canvas")).getContext("webgl") 
				|| canvas.getContext("experimental-webgl")) ? 1 : -1);
		}
		return (Jmol.featureDetection._webGLtest > 0);
	};

	features.supportsLocalization = function() {
		//<meta charset="utf-8">                                     
		var metas = document.getElementsByTagName('meta'); 
		for (var i= metas.length; --i >= 0;) 
			if (metas[i].outerHTML.toLowerCase().indexOf("utf-8") >= 0) return true;
		return false;
		};

	features.supportsJava = function() {
		if (!Jmol.featureDetection._javaEnabled) {
			if (Jmol._isMsie) {
				if (!navigator.javaEnabled()) {
					Jmol.featureDetection._javaEnabled = -1;
				} else {
					//more likely -- would take huge testing
					Jmol.featureDetection._javaEnabled = 1;
				}
			} else {
				Jmol.featureDetection._javaEnabled = (navigator.javaEnabled() && (!navigator.mimeTypes || navigator.mimeTypes["application/x-java-applet"]) ? 1 : -1);
			}
		}
		return (Jmol.featureDetection._javaEnabled > 0);
	};

	features.compliantBrowser = function() {
		var a = !!document.getElementById;
		var os = features.os;
		// known exceptions (old browsers):
		if (features.browserName == "opera" && features.browserVersion <= 7.54 && os == "mac"
			|| features.browserName == "webkit" && features.browserVersion < 125.12
			|| features.browserName == "msie" && os == "mac"
			|| features.browserName == "konqueror" && features.browserVersion <= 3.3
		) a = false;
		return a;
	}

	features.isFullyCompliant = function() {
		return features.compliantBrowser() && features.supportsJava();
	}

	features.useIEObject = (features.os == "win" && features.browserName == "msie" && features.browserVersion >= 5.5);
	features.useHtml4Object = (features.browserName == "mozilla" && features.browserVersion >= 5) ||
		(features.browserName == "opera" && features.browserVersion >= 8) ||
		(features.browserName == "webkit"/* && features.browserVersion >= 412.2 && features.browserVersion < 500*/); // 500 is a guess; required for 536.3

	features.hasFileReader = (window.File && window.FileReader);

	return features;

})(document, window);


		////////////// AJAX-related core functionality //////////////

	Jmol._ajax = function(info) {
		if (!info.async) {
			return Jmol.$ajax(info).responseText;
		}
		Jmol._ajaxQueue.push(info)
		if (Jmol._ajaxQueue.length == 1)
			Jmol._ajaxDone()
	}
	Jmol._ajaxDone = function() {
		var info = Jmol._ajaxQueue.shift();
		info && Jmol.$ajax(info);
	}

	Jmol._grabberOptions = [
		["$", "NCI(small molecules)"],
		[":", "PubChem(small molecules)"],
		["=", "RCSB(macromolecules)"]
	];

	Jmol._getGrabberOptions = function(applet) {
		// feel free to adjust this look to anything you want
		if (Jmol._grabberOptions.length == 0)
			return ""


		var s = '<input type="text" id="ID_query" onfocus="jQuery(this).select()" onkeypress="if(13==event.which){Jmol._applets[\'ID\']._search();return false}" size="32" value="" />';
		var b = '<button id="ID_submit" onclick="Jmol._applets[\'ID\']._search()">Search</button></nobr>'
		if (Jmol._grabberOptions.length == 1) {
			s = '<nobr>' + s + '<span style="display:none">';
			b = '</span>' + b;
		} else {
			s += '<br /><nobr>'
		}
		s += '<select id="ID_select">'
		for (var i = 0; i < Jmol._grabberOptions.length; i++) {
			var opt = Jmol._grabberOptions[i];
			s += '<option value="' + opt[0] + '" ' + (i == 0 ? 'selected' : '') + '>' + opt[1] + '</option>';
		}
		s = (s + '</select>' + b).replace(/ID/g, applet._id);
		return '<br />' + s;
	}

	Jmol._getScriptForDatabase = function(database) {
		return (database == "$" ? Jmol.db._nciLoadScript : database == ":" ? Jmol.db._pubChemLoadScript : Jmol.db._fileLoadScript);
	}

	 //   <dataset><record><structureId>1BLU</structureId><structureTitle>STRUCTURE OF THE 2[4FE-4S] FERREDOXIN FROM CHROMATIUM VINOSUM</structureTitle></record><record><structureId>3EUN</structureId><structureTitle>Crystal structure of the 2[4Fe-4S] C57A ferredoxin variant from allochromatium vinosum</structureTitle></record></dataset>

	Jmol._setInfo = function(applet, database, data) {
		var info = [];
		var header = "";
		if (data.indexOf("ERROR") == 0)
			header = data;
		else
			switch (database) {
			case "=":
				var S = data.split("<dimStructure.structureId>");
				var info = ["<table>"];
				for (var i = 1; i < S.length; i++) {
					info.push("<tr><td valign=top><a href=\"javascript:Jmol.search(" + applet._id + ",'=" + S[i].substring(0, 4) + "')\">" + S[i].substring(0, 4) + "</a></td>");
					info.push("<td>" + S[i].split("Title>")[1].split("</")[0] + "</td></tr>");
				}
				info.push("</table>");
				header = (S.length - 1) + " matches";
				break;      
			case "$": // NCI
			case ":": // pubChem
			break;
			default:
				return;
		}
		applet._infoHeader = header;
		applet._info = info.join("");
		applet._showInfo(true);
	}

	Jmol._loadSuccess = function(a, fSuccess) {
		if (!fSuccess)
			return;
		Jmol._ajaxDone();
		fSuccess(a);
	}

	Jmol._loadError = function(fError){
		Jmol._ajaxDone();
		Jmol.say("Error connecting to server: " + Jmol._ajaxCall);  
		null!=fError&&fError()
	}

	Jmol._isDatabaseCall = function(query) {
		return (Jmol.db._databasePrefixes.indexOf(query.substring(0, 1)) >= 0);
	}

	
	Jmol._getDirectDatabaseCall = function(query, checkXhr2) {
		if (checkXhr2 && !Jmol.featureDetection.supportsXhr2())
			return query;
		var pt = 2;
		var db;
		var call = Jmol.db._DirectDatabaseCalls[query.substring(0,pt)];
		if (!call)
			call = Jmol.db._DirectDatabaseCalls[db = query.substring(0,--pt)];
		if (call) {
			if (db == ":") {
				var ql = query.toLowerCase();
				if (!isNaN(parseInt(query.substring(1)))) {
					query = "cid/" + query.substring(1);
				} else if (ql.indexOf(":smiles:") == 0) {
					call += "?POST?smiles=" + query.substring(8);
					query = "smiles";
				} else if (ql.indexOf(":cid:") == 0) {
					query = "cid/" + query.substring(5);
				} else {
					if (ql.indexOf(":name:") == 0)
						query = query.substring(5);
					else if (ql.indexOf(":cas:") == 0)
						query = query.substring(4);
					query = "name/" + encodeURIComponent(query.substring(pt));
				}
			} else {
				query = encodeURIComponent(query.substring(pt));		
			}
			if (call.indexOf("FILENCI") >= 0) {
				query = query.replace(/\%2F/g, "/");				
				query = call.replace(/\%FILENCI/, query);
			} else {
				query = call.replace(/\%FILE/, query);
			}
		}		
		return query;
	}

	Jmol._getRawDataFromServer = function(database,query,fSuccess,fError,asBase64,noScript){
	  // note that this method is now only enabled for "_"
	  // server-side processing of database queries was too slow and only useful for 
	  // the IMAGE option, which has been abandoned.
		var s = 
			"?call=getRawDataFromDatabase&database=" + database + (query.indexOf("?POST?") >= 0 ? "?POST?" : "")
				+ "&query=" + encodeURIComponent(query)
				+ (asBase64 ? "&encoding=base64" : "")
				+ (noScript ? "" : "&script=" + encodeURIComponent(Jmol._getScriptForDatabase(database)));
		return Jmol._contactServer(s, fSuccess, fError);
	}

	Jmol._checkFileName = function(applet, fileName, isRawRet) {
		if (Jmol._isDatabaseCall(fileName)) {
			if (isRawRet)
				Jmol._setQueryTerm(applet, fileName);
			fileName = Jmol._getDirectDatabaseCall(fileName, true);
			if (Jmol._isDatabaseCall(fileName)) {
				// xhr2 not supported (MSIE)
				fileName = Jmol._getDirectDatabaseCall(fileName, false);
				isRawRet && (isRawRet[0] = true);
			}
		}
		return fileName;
	}
	
	Jmol._checkCache = function(applet, fileName, fSuccess) {
		if (applet._cacheFiles && Jmol._fileCache && !fileName.endsWith(".js")) {
			var data = Jmol._fileCache[fileName];
			if (data) {
				System.out.println("using "  + data.length + " bytes of cached data for "  + fileName);
				fSuccess(data);
				return null;
			} else {
				fSuccess = function(fileName, data) { fSuccess(Jmol._fileCache[fileName] = data) };     
			}
		}
		return fSuccess;
	}
	
	Jmol._loadFileData = function(applet, fileName, fSuccess, fError){
		var isRaw = [];
		fileName = Jmol._checkFileName(applet, fileName, isRaw);
		fSuccess = Jmol._checkCache(applet, fileName, fSuccess);
		if (isRaw[0]) {
				Jmol._getRawDataFromServer("_",fileName,fSuccess,fError);   
				return;
		} 
		var info = {
			type: "GET",
			dataType: "text",
			url: fileName,
			async: Jmol._asynchronous,
			success: function(a) {Jmol._loadSuccess(a, fSuccess)},
			error: function() {Jmol._loadError(fError)}
		}
		Jmol._checkAjaxPost(info);
		Jmol._ajax(info);
	}

	Jmol._getInfoFromDatabase = function(applet, database, query){
		if (database == "====") {
			var data = Jmol.db._restQueryXml.replace(/QUERY/,query);

			var info = {
				dataType: "text",
				type: "POST",
				contentType:"application/x-www-form-urlencoded",
				url: Jmol.db._restQueryUrl,
				data: encodeURIComponent(data) + "&req=browser",
				success: function(data) {Jmol._ajaxDone();Jmol._extractInfoFromRCSB(applet, database, query, data)},
				error: function() {Jmol._loadError(null)},
				async: Jmol._asynchronous
			}
			return Jmol._ajax(info);
		}   
		query = "?call=getInfoFromDatabase&database=" + database
				+ "&query=" + encodeURIComponent(query);
		return Jmol._contactServer(query, function(data) {Jmol._setInfo(applet, database, data)});
	}

	Jmol._extractInfoFromRCSB = function(applet, database, query, output) {
		var n = output.length/5;
		if (n == 0)
			return; 
		if (query.length == 4 && n != 1) {
			var QQQQ = query.toUpperCase();
			var pt = output.indexOf(QQQQ);
			if (pt > 0 && "123456789".indexOf(QQQQ.substring(0, 1)) >= 0)
				output = QQQQ + "," + output.substring(0, pt) + output.substring(pt + 5);
			if (n > 50)
				output = output.substring(0, 250);
			output = output.replace(/\n/g,",");
			var url = Jmol._restReportUrl.replace(/IDLIST/,output);
			Jmol._loadFileData(applet, url, function(data) {Jmol._setInfo(applet, database, data) });   
		}
	}

	Jmol._checkAjaxPost = function(info) {
		var pt = info.url.indexOf("?POST?");
		if (pt > 0) {
			info.data = info.url.substring(pt + 6);
			info.url = info.url.substring(0, pt);
			info.type = "POST";
			info.contentType = "application/x-www-form-urlencoded";
		}
	}
	Jmol._contactServer = function(data,fSuccess,fError){
		var info = {
			dataType: "text",
			type: "GET",
			url: Jmol._serverUrl + data,
			success: function(a) {Jmol._loadSuccess(a, fSuccess)},
			error:function() { Jmol._loadError(fError) },
			async:fSuccess ? Jmol._asynchronous : false
		}
		Jmol._checkAjaxPost(info);
		return Jmol._ajax(info);
	}

	Jmol._setQueryTerm = function(applet, query) {
		if (!query || !applet._hasOptions || query.substring(0, 7) == "http://")
			return;
		if (Jmol._isDatabaseCall(query)) {
			var database = query.substring(0, 1);
			query = query.substring(1);
			if (query.substring(0,1) == database && "=$".indexOf(database) >= 0)
				query = query.substring(1);
			var d = Jmol._getElement(applet, "select");
			if (d && d.options)
				for (var i = 0; i < d.options.length; i++)
					if (d[i].value == database)
						d[i].selected = true;
		}
		Jmol.$val(Jmol.$(applet, "query"), query);
	}

	Jmol._search = function(applet, query, script) {
		arguments.length > 1 || (query = null);
		Jmol._setQueryTerm(applet, query);
		query || (query = Jmol.$val(Jmol.$(applet, "query")))
		if (query.indexOf("!") == 0) {
		// command prompt in this box as well
		// remove exclamation point "immediate" indicator
			applet._script(query.substring(1));
			return;
		} 
		query && (query = query.replace(/\"/g, ""));
		applet._showInfo(false);
		Jmol._searchMol(applet, query, script, true);
	}

	Jmol._searchMol = function(applet, query, script, checkView) {
		var database;
		if (Jmol._isDatabaseCall(query)) {
			database = query.substring(0, 1);
			query = query.substring(1);
		} else {
			database = (applet._hasOptions ? Jmol.$val(Jmol.$(applet, "select")) : "$");
		}
		if (database == "=" && query.length == 3)
			query = "=" + query; // this is a ligand      
		var dm = database + query;
		if (!query || dm.indexOf("?") < 0 && dm == applet._thisJmolModel) {
			return;    
		}
		applet._thisJmolModel = dm;
		var view;
		if (checkView && applet._viewSet != null && (view = Jmol.View.__findView(applet._viewSet, {chemID:dm})) != null) {
			Jmol.View.__setView(view, applet, false);
			return;
		}

		if (database == "$" || database == ":")
			applet._jmolFileType = "MOL";
		else if (database == "=")
			applet._jmolFileType = "PDB";
		applet._searchDatabase(query, database, script);
	}

	Jmol._searchDatabase = function(applet, query, database, script) {
		applet._showInfo(false);
		if (query.indexOf("?") >= 0) {
			Jmol._getInfoFromDatabase(applet, database, query.split("?")[0]);
			return true;
		}
		if (Jmol.db._DirectDatabaseCalls[database]) {
			applet._loadFile(database + query, script);
			return true;
		}
		return false;
	}

	Jmol._syncBinaryOK="?";

	Jmol._canSyncBinary = function(isSilent) {
		if (self.VBArray) return (Jmol._syncBinaryOK = false);
		if (Jmol._syncBinaryOK != "?") return Jmol._syncBinaryOK;
		Jmol._syncBinaryOK = true;
		try {
			var xhr = new window.XMLHttpRequest();
			xhr.open( "text", Jmol._ajaxTestSite, false );
			if (xhr.hasOwnProperty("responseType")) {
				xhr.responseType = "arraybuffer";
			} else if (xhr.overrideMimeType) {
				xhr.overrideMimeType('text/plain; charset=x-user-defined');
			}
		} catch( e ) {
			var s = "JSmolCore.js: synchronous binary file transfer is requested but not available";
			System.out.println(s);
			if (Jmol._alertNoBinary && !isSilent)
				alert (s)
			return Jmol._syncBinaryOK = false;
		}
		return true;  
	}

	Jmol._binaryTypes = [".gz",".jpg",".png",".zip",".jmol",".bin",".smol",".spartan",".mrc",".pse", ".map", ".omap"];

	Jmol._isBinaryUrl = function(url) {
		for (var i = Jmol._binaryTypes.length; --i >= 0;)
			if (url.indexOf(Jmol._binaryTypes[i]) >= 0) return true;
		return false;
	}

	Jmol._getFileData = function(fileName, fSuccess) {
		// use host-server PHP relay if not from this host
		var type = (Jmol._isBinaryUrl(fileName) ? "binary" : "text");
		var isPDB = (fileName.indexOf("pdb.gz") >= 0 && fileName.indexOf("http://www.rcsb.org/pdb/files/") == 0);
		var asBase64 = (type == "binary" && !Jmol._canSyncBinary(isPDB));
		if (asBase64 && isPDB) {
			// avoid unnecessary binary transfer
			fileName = fileName.replace(/pdb\.gz/,"pdb");
			asBase64 = false;
			type = "text";
		}
		var isPost = (fileName.indexOf("?POST?") >= 0);
		if (fileName.indexOf("file:/") == 0 && fileName.indexOf("file:///") != 0)
			fileName = "file://" + fileName.substring(5);      /// fixes IE problem
		var isMyHost = (fileName.indexOf("://") < 0 || fileName.indexOf(document.location.protocol) == 0 && fileName.indexOf(document.location.host) >= 0);
		var isDirectCall = Jmol._isDirectCall(fileName);
		//if (fileName.indexOf("http://pubchem.ncbi.nlm.nih.gov/") == 0)isDirectCall = false;

		var cantDoSynchronousLoad = (!isMyHost && Jmol.$supportsIECrossDomainScripting());
		if (!fSuccess || asBase64)
			if (cantDoSynchronousLoad || asBase64 || !isMyHost && !isDirectCall)
				return Jmol._getRawDataFromServer("_",fileName, fSuccess, fSuccess, asBase64, true);
		fileName = fileName.replace(/file:\/\/\/\//, "file://"); // opera
		var info = {dataType:type,async:!!fSuccess};
		if (isPost) {
			info.type = "POST";
			info.url = fileName.split("?POST?")[0]
			info.data = fileName.split("?POST?")[1]
		} else {
			info.type = "GET";
			info.url = fileName;
		}
		if (fSuccess) {
			info.success = function(data) { fSuccess(Jmol._xhrReturn(info.xhr))};
			info.error = function() { fSuccess(xhr.statusText)};
		}
		info.xhr = Jmol.$ajax(info);
		if (!fSuccess) 
			return Jmol._xhrReturn(info.xhr);			
	}
	
	Jmol._xhrReturn = function(xhr){
		if (!xhr.responseText || self.Clazz && Clazz.instanceOf(xhr.response, self.ArrayBuffer)) {
			// Safari or error 
			return xhr.response || xhr.statusText;
		} 
		return xhr.responseText;
	}

	Jmol._isDirectCall = function(url) {
		for (var key in Jmol.db._DirectDatabaseCalls) {
			if (key.indexOf(".") >= 0 && url.indexOf(key) >= 0)
				return true;
		}
		return false;
	}

	Jmol._cleanFileData = function(data) {
		if (data.indexOf("\r") >= 0 && data.indexOf("\n") >= 0) {
			return data.replace(/\r\n/g,"\n");
		}
		if (data.indexOf("\r") >= 0) {
			return data.replace(/\r/g,"\n");
		}
		return data;
	};

	Jmol._getFileType = function(name) {
		var database = name.substring(0, 1);
		if (database == "$" || database == ":")
			return "MOL";
		if (database == "=")
			return (name.substring(1,2) == "=" ? "LCIF" : "PDB");
		// just the extension, which must be PDB, XYZ..., CIF, or MOL
		name = name.split('.').pop().toUpperCase();
		return name.substring(0, Math.min(name.length, 3));
	};

	Jmol._getZ = function(applet, what) {
		return applet && applet._z && applet._z[what] || Jmol._z[what];
	}
	
	Jmol._incrZ = function(applet, what) {
		return applet && applet._z && ++applet._z[what] || ++Jmol._z[what];
	}
	
	Jmol._loadFileAsynchronously = function(fileLoadThread, applet, fileName, appData) {
		if (fileName.indexOf("?") != 0) {
			// LOAD ASYNC command
			fileName = Jmol._checkFileName(applet, fileName);
			var fSuccess = function(data) {Jmol._setData(fileLoadThread, fileName, data, appData)};
			fSuccess = Jmol._checkCache(applet, fileName, fSuccess);
			return (fSuccess == null ? null : Jmol._getFileData(fileName, fSuccess));		
		}
		// we actually cannot suggest a fileName, I believe.
		if (!Jmol.featureDetection.hasFileReader)
				return fileLoadThread.setData("Local file reading is not enabled in your browser", null, appData);
		if (!applet._localReader) {
			var div = '<div id="ID" style="z-index:'+Jmol._getZ(applet, "fileOpener") + ';position:absolute;background:#E0E0E0;left:10px;top:10px"><div style="margin:5px 5px 5px 5px;"><input type="file" id="ID_files" /><button id="ID_loadfile">load</button><button id="ID_cancel">cancel</button></div><div>'
			Jmol.$after("#" + applet._id + "_appletdiv", div.replace(/ID/g, applet._id + "_localReader"));
			applet._localReader = Jmol.$(applet, "localReader");
		}
		Jmol.$appEvent(applet, "localReader_loadfile", "click");
		Jmol.$appEvent(applet, "localReader_loadfile", "click", function(evt) {
			var file = Jmol.$(applet, "localReader_files")[0].files[0];   
			var reader = new FileReader();
			reader.onloadend = function(evt) {
				if (evt.target.readyState == FileReader.DONE) { // DONE == 2
					Jmol.$css(Jmol.$(applet, "localReader"), {display : "none"});
					Jmol._setData(fileLoadThread, file.name, evt.target.result, appData);
				}
			};
			reader.readAsArrayBuffer(file);
		});
		Jmol.$appEvent(applet, "localReader_cancel", "click");
		Jmol.$appEvent(applet, "localReader_cancel", "click", function(evt) {
			Jmol.$css(Jmol.$(applet, "localReader"), {display: "none"});
			fileLoadThread.setData(null, appData);
		});
		Jmol.$css(Jmol.$(applet, "localReader"), {display : "block"});
	}

  Jmol._setData = function(fileLoadThread, filename, data, appData) {
  	data = Jmol._strToBytes(data);
		if (filename.indexOf(".jdx") >= 0)
			Jmol.Cache.put("cache://" + filename, data);
		fileLoadThread.setData(filename, data, appData);
  }
  
	Jmol._toBytes = function(data) {
	if (typeof data == "string") 
		return data.getBytes();
	// ArrayBuffer assumed here
	data = new Uint8Array(data);
	var b = Clazz.newByteArray(data.length, 0);
		for (var i = data.length; --i >= 0;)
			b[i] = data[i];
	return b;
	}

	Jmol._doAjax = function(url, postOut, dataOut) {
		// called by org.jmol.awtjs2d.JmolURLConnection.doAjax()
		url = url.toString();

		if (dataOut != null) 
			return Jmol._saveFile(url, dataOut);
		if (postOut)
			url += "?POST?" + postOut;
		var data = Jmol._getFileData(url)
		return Jmol._processData(data, Jmol._isBinaryUrl(url));
	}

	// Jmol._localFileSaveFunction --  // do something local here; Maybe try the FileSave interface? return true if successful
	 
	Jmol._saveFile = function(filename, data, mimetype, encoding) {
		if (Jmol._localFileSaveFunction && Jmol._localFileSaveFunction(filename, data))
			return "OK";
		var filename = filename.substring(filename.lastIndexOf("/") + 1);
		mimetype || (mimetype = (filename.indexOf(".pdf") >= 0 ? "application/pdf" : filename.indexOf(".png") >= 0 ? "image/png" : filename.indexOf(".jpg") >= 0 ? "image/jpg" : ""));
		var isString = (typeof data == "string");
		if (!isString)
			data = (JU ? JU : J.util).Base64.getBase64(data).toString();
		encoding || (encoding = (isString ? "" : "base64"));
		var url = Jmol._serverUrl;
		url && url.indexOf("your.server") >= 0 && (url = "");
		if (Jmol._useDataURI || !url) {
			// Asynchronous output generated using an anchor tag
			encoding || (data = btoa(data));
			var a = document.createElement("a");
			a.href = "data:" + mimetype + ";base64," + data;
			a.type = mimetype || (mimetype = "text/plain");	
			a.download = filename;
			a.target = "_blank";
				$("body").append(a);
			a.click();
			a.remove();		
		} else {
		// Asynchronous outputto be reflected as a download
			if (!Jmol._formdiv) {
					var sform = '<div id="__jsmolformdiv__" style="display:none">\
						<form id="__jsmolform__" method="post" target="_blank" action="">\
						<input name="call" value="saveFile"/>\
						<input id="__jsmolmimetype__" name="mimetype" value=""/>\
						<input id="__jsmolencoding__" name="encoding" value=""/>\
						<input id="__jsmolfilename__" name="filename" value=""/>\
						<textarea id="__jsmoldata__" name="data"></textarea>\
						</form>\
						</div>'
				Jmol.$after("body", sform);
				Jmol._formdiv = "__jsmolform__";
			}
			Jmol.$attr(Jmol._formdiv, "action", url + "?" + (new Date()).getMilliseconds());
			Jmol.$val("__jsmoldata__", data);
			Jmol.$val("__jsmolfilename__", filename);
			Jmol.$val("__jsmolmimetype__", mimetype);
			Jmol.$val("__jsmolencoding__", encoding);
			Jmol.$submit("__jsmolform__");
			Jmol.$val("__jsmoldata__", "");
			Jmol.$val("__jsmolfilename__", "");
		}
		return "OK";
	}

	Jmol._processData = function(data, isBinary) {
		if (typeof data == "undefined") {
			data = "";
			isBinary = false;
		}
		if (isBinary)
			isBinary = Jmol._canSyncBinary();
		return (isBinary ? Jmol._strToBytes(data) : JU.SB.newS(data));
	};

	Jmol._strToBytes = function(s) {
		if (Clazz.instanceOf(s, self.ArrayBuffer))
			return Jmol._toBytes(s);
		var b = Clazz.newByteArray(s.length, 0);
		for (var i = s.length; --i >= 0;)
			b[i] = s.charCodeAt(i) & 0xFF;
		return b;
	}

	////////////// applet start-up functionality //////////////

	Jmol._setConsoleDiv = function (d) {
		if (!self.Clazz)return;
		Clazz.setConsoleDiv(d);
	}

	Jmol._setJmolParams = function(params, Info, isHashtable) {
		var availableValues = ";progressbar;progresscolor;boxbgcolor;boxfgcolor;allowjavascript;boxmessage;\
									;messagecallback;pickcallback;animframecallback;appletreadycallback;atommovedcallback;\
									;echocallback;evalcallback;hovercallback;language;loadstructcallback;measurecallback;\
									;minimizationcallback;resizecallback;scriptcallback;statusform;statustext;statustextarea;\
									;synccallback;usecommandthread;syncid;appletid;startupscript;menufile;";
		for (var i in Info)
			if(availableValues.indexOf(";" + i.toLowerCase() + ";") >= 0){
				if (i == "language" && !Jmol.featureDetection.supportsLocalization())
					continue;
				if (isHashtable)
					params.put(i, (Info[i] === true ? Boolean.TRUE: Info[i] === false ? Boolean.FALSE : Info[i]))
				else
					params[i] = Info[i];
			}
	}     
	 
	Jmol._registerApplet = function(id, applet) {
		return window[id] = Jmol._applets[id] = Jmol._applets[applet] = Jmol._applets[id + "__" + Jmol._syncId + "__"] = applet;
	} 

	Jmol._readyCallback = function (appId,fullId,isReady,jmolApp) {
		appId = appId.split("_object")[0];
		isReady = (isReady.booleanValue ? isReady.booleanValue() : isReady);
		// necessary for MSIE in strict mode -- apparently, we can't call 
		// jmol._readyCallback, but we can call Jmol._readyCallback. Go figure...

		Jmol._track(Jmol._applets[appId])._readyCallback(appId, fullId, isReady, jmolApp);
	}

	Jmol._getWrapper = function(applet, isHeader) {

			// id_appletinfotablediv
			//     id_appletdiv
			//     id_coverdiv
			//     id_infotablediv
			//       id_infoheaderdiv
			//          id_infoheaderspan
			//          id_infocheckboxspan
			//       id_infodiv


			// for whatever reason, without DOCTYPE, with MSIE, "height:auto" does not work, 
			// and the text scrolls off the page.
			// So I'm using height:95% here.
			// The table was a fix for MSIE with no DOCTYPE tag to fix the miscalculation
			// in height of the div when using 95% for height. 
			// But it turns out the table has problems with DOCTYPE tags, so that's out. 
			// The 95% is a compromise that we need until the no-DOCTYPE MSIE solution is found. 
			// (100% does not work with the JME linked applet)
		var s;
		// ... here are just for clarification in this code; they are removed immediately
		if (isHeader) {
			var img = ""; 
			if (applet._coverImage){
				var more = " onclick=\"Jmol.coverApplet(ID, false)\" title=\"" + (applet._coverTitle) + "\"";
				var play = "<image id=\"ID_coverclickgo\" src=\"" + applet._j2sPath + "/img/play_make_live.jpg\" style=\"width:25px;height:25px;position:absolute;bottom:10px;left:10px;"
					+ "z-index:" + Jmol._getZ(applet, "coverImage")+";opacity:0.5;\"" + more + " />"  
				img = "<div id=\"ID_coverdiv\" style=\"background-color:red;z-index:" + Jmol._getZ(applet, "coverImage")+";width:100%;height:100%;display:inline;position:absolute;top:0px;left:0px\"><image id=\"ID_coverimage\" src=\""
				 + applet._coverImage + "\" style=\"width:100%;height:100%\"" + more + "/>" + play + "</div>";
			}
			s = "\
...<div id=\"ID_appletinfotablediv\" style=\"width:Wpx;height:Hpx;position:relative;font-size:14px;text-align:left\">IMG\
......<div id=\"ID_appletdiv\" style=\"z-index:" + Jmol._getZ(applet, "header") + ";width:100%;height:100%;position:absolute;top:0px;left:0px;\">";
			var height = applet._height;
			var width = applet._width;
			if (typeof height !== "string" || height.indexOf("%") < 0) 
				height += "px";
			if (typeof width !== "string" || width.indexOf("%") < 0)
				width += "px";
			s = s.replace(/IMG/, img).replace(/Hpx/g, height).replace(/Wpx/g, width);
		} else {
			s = "\
......</div>\
......<div id=\"ID_2dappletdiv\" style=\"position:absolute;width:100%;height:100%;overflow:hidden;display:none\"></div>\
......<div id=\"ID_infotablediv\" style=\"width:100%;height:100%;position:absolute;top:0px;left:0px\">\
.........<div id=\"ID_infoheaderdiv\" style=\"height:20px;width:100%;background:yellow;display:none\"><span id=\"ID_infoheaderspan\"></span><span id=\"ID_infocheckboxspan\" style=\"position:absolute;text-align:right;right:1px;\"><a href=\"javascript:Jmol.showInfo(ID,false)\">[x]</a></span></div>\
.........<div id=\"ID_infodiv\" style=\"position:absolute;top:20px;bottom:0px;width:100%;height:100%;overflow:auto\"></div>\
......</div>\
...</div>";
		}
		return s.replace(/\.\.\./g,"").replace(/[\n\r]/g,"").replace(/ID/g, applet._id);
	}

	Jmol._documentWrite = function(text) {
		if (Jmol._document) {
			if (Jmol._isXHTML && !Jmol._XhtmlElement) {
				var s = document.getElementsByTagName("script");
				Jmol._XhtmlElement = s.item(s.length - 1);
				Jmol._XhtmlAppendChild = false;
			}
			if (Jmol._XhtmlElement)
				Jmol._domWrite(text);
			else
				Jmol._document.write(text);
		}
		return text;
	}

	Jmol._domWrite = function(data) {
		var pt = 0
		var Ptr = []
		Ptr[0] = 0
		while (Ptr[0] < data.length) {
			var child = Jmol._getDomElement(data, Ptr);
			if (!child)
				break;
			if (Jmol._XhtmlAppendChild)
				Jmol._XhtmlElement.appendChild(child);
			else
				Jmol._XhtmlElement.parentNode.insertBefore(child, _jmol.XhtmlElement);
		}
	}

	Jmol._getDomElement = function(data, Ptr, closetag, lvel) {

		// there is no "document.write" in XHTML

		var e = document.createElement("span");
		e.innerHTML = data;
		Ptr[0] = data.length;

/*
	// unnecessary ?  

		closetag || (closetag = "");
		lvel || (lvel = 0);
		var pt0 = Ptr[0];
		var pt = pt0;
		while (pt < data.length && data.charAt(pt) != "<") 
			pt++
		if (pt != pt0) {
			var text = data.substring(pt0, pt);
			Ptr[0] = pt;
			return document.createTextNode(text);
		}
		pt0 = ++pt;
		var ch;
		while (pt < data.length && "\n\r\t >".indexOf(ch = data.charAt(pt)) < 0) 
			pt++;
		var tagname = data.substring(pt0, pt);
		var e = (tagname == closetag  || tagname == "/" ? ""
			: document.createElementNS ? document.createElementNS('http://www.w3.org/1999/xhtml', tagname)
			: document.createElement(tagname));
		if (ch == ">") {
			Ptr[0] = ++pt;
			return e;
		}
		while (pt < data.length && (ch = data.charAt(pt)) != ">") {
			while (pt < data.length && "\n\r\t ".indexOf(ch = data.charAt(pt)) >= 0) 
				pt++;
			pt0 = pt;
			while (pt < data.length && "\n\r\t =/>".indexOf(ch = data.charAt(pt)) < 0) 
				pt++;
			var attrname = data.substring(pt0, pt).toLowerCase();
			if (attrname && ch != "=")
				e.setAttribute(attrname, "true");
			while (pt < data.length && "\n\r\t ".indexOf(ch = data.charAt(pt)) >= 0) 
				pt++;
			if (ch == "/") {
				Ptr[0] = pt + 2;
				return e;
			} else if (ch == "=") {
				var quote = data.charAt(++pt);
				pt0 = ++pt;
				while (pt < data.length && (ch = data.charAt(pt)) != quote) 
					pt++;
				var attrvalue = data.substring(pt0, pt);
				e.setAttribute(attrname, attrvalue);
				pt++;
			}
		}
		Ptr[0] = ++pt;
		while (Ptr[0] < data.length) {
			var child = Jmol._getDomElement(data, Ptr, "/" + tagname, lvel+1);
			if (!child)
				break;
			e.appendChild(child);
		}
*/
		return e;    
	}

	Jmol._setObject = function(obj, id, Info) {
		obj._id = id;
		obj.__Info = {};
		Info.z && Info.zIndexBase && (Jmol._z = Jmol._getZOrders(Info.zIndexBase));
		for (var i in Info)
			obj.__Info[i] = Info[i];
		(obj._z = Info.z) || Info.zIndexBase && (obj._z = obj.__Info.z = Jmol._getZOrders(Info.zIndexBase));
		obj._width = Info.width;
		obj._height = Info.height;
		obj._noscript = !obj._isJava && Info.noscript;
		obj._console = Info.console;
		obj._cacheFiles = !!Info.cacheFiles;
		obj._viewSet = (Info.viewSet == null || obj._isJava ? null : "Set" + Info.viewSet);
		if (obj._viewSet != null) {
			Jmol.View.__init(obj);
			obj._currentView = null;
		}
		!Jmol._fileCache && obj._cacheFiles && (Jmol._fileCache = {});
		if (!obj._console)
			obj._console = obj._id + "_infodiv";
		if (obj._console == "none")
			obj._console = null;

		obj._color = (Info.color ? Info.color.replace(/0x/,"#") : "#FFFFFF");
		obj._disableInitialConsole = Info.disableInitialConsole;
		obj._noMonitor = Info.disableJ2SLoadMonitor;
		Jmol._j2sPath && (Info.j2sPath = Jmol._j2sPath);
		obj._j2sPath = Info.j2sPath;
		obj._deferApplet = Info.deferApplet;
		obj._deferUncover = Info.deferUncover;
		obj._coverImage = !obj._isJava && Info.coverImage;
		obj._isCovered = !!obj._coverImage; 
		obj._coverScript = Info.coverScript;
		obj._coverTitle = Info.coverTitle;

		if (!obj._coverTitle)
			obj._coverTitle = (obj._deferApplet ? "activate 3D model" : "3D model is loading...")
		obj._containerWidth = obj._width + ((obj._width==parseFloat(obj._width))? "px":"");
		obj._containerHeight = obj._height + ((obj._height==parseFloat(obj._height))? "px":"");
		obj._info = "";
		obj._infoHeader = obj._jmolType + ' "' + obj._id + '"'
		obj._hasOptions = Info.addSelectionOptions;
		obj._defaultModel = Info.defaultModel;
		obj._readyScript = (Info.script ? Info.script : "");
		obj._readyFunction = Info.readyFunction;
		if (obj._coverImage && !obj._deferApplet)
			obj._readyScript += ";javascript " + id + "._displayCoverImage(false)";
		obj._src = Info.src;

	}

	Jmol._addDefaultInfo = function(Info, DefaultInfo) {
		for (var x in DefaultInfo)
			if (typeof Info[x] == "undefined")
				Info[x] = DefaultInfo[x];
		Jmol._use && (Info.use = Jmol._use);
		if (Info.use.indexOf("SIGNED") >= 0) {
			if (Info.jarFile.indexOf("Signed") < 0)
				Info.jarFile = Info.jarFile.replace(/Applet/,"AppletSigned");
			Info.use = Info.use.replace(/SIGNED/, "JAVA");
			Info.isSigned = true;
		}
	}

	Jmol._syncedApplets = [];
	Jmol._syncedCommands = [];
	Jmol._syncedReady = [];
	Jmol._syncReady = false;
	Jmol._isJmolJSVSync = false;

	Jmol._setReady = function(applet) {
		Jmol._syncedReady[applet] = 1;
		var n = 0;
		for (var i = 0; i < Jmol._syncedApplets.length; i++) {
			if (Jmol._syncedApplets[i] == applet._id) {
				Jmol._syncedApplets[i] = applet;
				Jmol._syncedReady[i] = 1;
			} else if (!Jmol._syncedReady[i]) {
				continue;
			}
			n++;
		}
		if (n != Jmol._syncedApplets.length)
			return;
		Jmol._setSyncReady();
	}

	Jmol._setDestroy = function(applet) {
		//MSIE bug responds to any link click even if it is just a JavaScript call

		if (Jmol.featureDetection.allowDestroy)
			Jmol.$windowOn('beforeunload', function () { Jmol._destroy(applet); } );
	}

	Jmol._destroy = function(applet) {
		try {
			if (applet._applet) applet._applet.destroy();
			applet._applet = null;
			Jmol._unsetMouse(applet._canvas)
			applet._canvas = null;
			var n = 0;
			for (var i = 0; i < Jmol._syncedApplets.length; i++) {
				if (Jmol._syncedApplets[i] == applet)
					Jmol._syncedApplets[i] = null;
				if (Jmol._syncedApplets[i])
					n++;
			}
			if (n > 0)
				return;
			Jmol._clearVars();
		} catch(e){}
	}

	////////////// misc core functionality //////////////

	Jmol._setSyncReady = function() {
		Jmol._syncReady = true;
		var s = ""
		for (var i = 0; i < Jmol._syncedApplets.length; i++)
			if (Jmol._syncedCommands[i])
				s += "Jmol.script(Jmol._syncedApplets[" + i + "], Jmol._syncedCommands[" + i + "]);"
		setTimeout(s, 50);  
	}

	Jmol._mySyncCallback = function(appFullName,msg) {
		app = Jmol._applets[appFullName];
		if (app._viewSet) {
			// when can we do this?
//			if (app._viewType == "JSV" && !app._currentView.JMOL)
				Jmol.View.updateFromSync(app, msg);
			return;
		}
		if(!Jmol._syncReady || !Jmol._isJmolJSVSync)
			return 1; // continue processing and ignore me
		for (var i = 0; i < Jmol._syncedApplets.length; i++) {
			if (msg.indexOf(Jmol._syncedApplets[i]._syncKeyword) >= 0)
				Jmol._syncedApplets[i]._syncScript(msg);
		}
		return 0 // prevents further Jmol sync processing 
	}              

	Jmol._getElement = function(applet, what) {
		var d = document.getElementById(applet._id + "_" + what);
		return (d || {});
	} 
	 
	Jmol._evalJSON = function(s,key){
		s = s + "";
		if(!s)
			return [];
		if(s.charAt(0) != "{") {
			if(s.indexOf(" | ") >= 0)
				s = s.replace(/\ \|\ /g, "\n");
			return s;
		}
		var A = (new Function( "return " + s ) )();
		return (!A ? null : key && A[key] != undefined ? A[key] : A);
	}

	Jmol._sortMessages = function(A){
		/*
		 * private function
		 */
		function _sortKey0(a,b){
			return (a[0]<b[0]?1:a[0]>b[0]?-1:0);
		}

		if(!A || typeof (A) != "object")
			return [];
		var B = [];
		for(var i = A.length - 1; i >= 0; i--)
			for(var j = 0, jj= A[i].length; j < jj; j++)
				B[B.length] = A[i][j];
		if(B.length == 0)
			return;
		B = B.sort(_sortKey0);
		return B;
	}

	//////////////////// mouse events //////////////////////

	Jmol._setMouseOwner = function(who, tf) {
		if (who == null || tf)
			Jmol._mouseOwner = who;
		else if (Jmol._mouseOwner == who)
			Jmol._mouseOwner = null;
	}

	Jmol._jsGetMouseModifiers = function(ev) {
		var modifiers = 0;
		switch (ev.button) {
		case 0:
			modifiers = 16;//J.api.Event.MOUSE_LEFT;
			break;
		case 1:
			modifiers = 8;//J.api.Event.MOUSE_MIDDLE;
			break;
		case 2:
			modifiers = 4;//J.api.Event.MOUSE_RIGHT;
			break;
		}
		if (ev.shiftKey)
			modifiers += 1;//J.api.Event.SHIFT_MASK;
		if (ev.altKey)
			modifiers += 8;//J.api.Event.ALT_MASK;
		if (ev.ctrlKey)
			modifiers += 2;//J.api.Event.CTRL_MASK;
		return modifiers;
	}

	Jmol._jsGetXY = function(canvas, ev) {
		if (!canvas.applet._ready || Jmol._touching && ev.type.indexOf("touch") < 0)
			return false;
		ev.preventDefault();
		var offsets = Jmol.$offset(canvas.id);
		var x, y;
		var oe = ev.originalEvent;
		// drag-drop jQuery event is missing pageX
		ev.pageX || (ev.pageX = oe.pageX);
		ev.pageY || (ev.pageY = oe.pageY);
		Jmol._mousePageX = ev.pageX;
		Jmol._mousePageY = ev.pageY;
		if (oe.targetTouches && oe.targetTouches[0]) {
			x = oe.targetTouches[0].pageX - offsets.left;
			y = oe.targetTouches[0].pageY - offsets.top;
		} else if (oe.changedTouches) {
			x = oe.changedTouches[0].pageX - offsets.left;
			y = oe.changedTouches[0].pageY - offsets.top;
		} else {
			x = ev.pageX - offsets.left;
			y = ev.pageY - offsets.top;
		}
		return (x == undefined ? null : [Math.round(x), Math.round(y), Jmol._jsGetMouseModifiers(ev)]);
	}

	Jmol._gestureUpdate = function(canvas, ev) {
		ev.stopPropagation();
		ev.preventDefault();
		var oe = ev.originalEvent;
		switch (ev.type) {
		case "touchstart":
			Jmol._touching = true;
			break;
		case "touchend":
			Jmol._touching = false;
			break;
		}
		if (!oe.touches || oe.touches.length != 2) return false;
		switch (ev.type) {
		case "touchstart":
			canvas._touches = [[],[]];
			break;
		case "touchmove":
			var offsets = Jmol.$offset(canvas.id);
			var t0 = canvas._touches[0];
			var t1 = canvas._touches[1];
			t0.push([oe.touches[0].pageX - offsets.left, oe.touches[0].pageY - offsets.top]);
			t1.push([oe.touches[1].pageX - offsets.left, oe.touches[1].pageY - offsets.top]);
			var n = t0.length;
			if (n > 3) {
				t0.shift();
				t1.shift();
			}
			if (n >= 2)
				canvas.applet._processGesture(canvas._touches);
			break;
		}
		return true;
	}

	Jmol._jsSetMouse = function(canvas) {
		Jmol.$bind(canvas, 'mousedown touchstart', function(ev) {
			Jmol._setMouseOwner(canvas, true);
			ev.stopPropagation();
			ev.preventDefault();
			canvas.isDragging = true;
			if ((ev.type == "touchstart") && Jmol._gestureUpdate(canvas, ev))
				return false;
			Jmol._setConsoleDiv(canvas.applet._console);
			var xym = Jmol._jsGetXY(canvas, ev);
			if(!xym)
				return false;
			if (ev.button != 2)
				Jmol.Swing.hideMenus(canvas.applet);

			canvas.applet._processEvent(501, xym); //J.api.Event.MOUSE_DOWN
			return false;
		});
		Jmol.$bind(canvas, 'mouseup touchend', function(ev) {
			Jmol._setMouseOwner(null);
			ev.stopPropagation();
			ev.preventDefault();
			canvas.isDragging = false;
			if (ev.type == "touchend" && Jmol._gestureUpdate(canvas, ev))
				return false;
			var xym = Jmol._jsGetXY(canvas, ev);
			if(!xym) return false;
			canvas.applet._processEvent(502, xym);//J.api.Event.MOUSE_UP
			return false;
		});
		Jmol.$bind(canvas, 'mousemove touchmove', function(ev) { // touchmove
		  // defer to console or menu when dragging within this canvas
			if (Jmol._mouseOwner && Jmol._mouseOwner != canvas && Jmol._mouseOwner.isDragging) {
				Jmol._mouseOwner.mouseMove(ev);
				return false;
			}
			return Jmol._drag(canvas, ev);
		});
		
		Jmol._drag = function(canvas, ev) {
			ev.stopPropagation();
			ev.preventDefault();
			var isTouch = (ev.type == "touchmove");
			if (isTouch && Jmol._gestureUpdate(canvas, ev))
				return false;
			var xym = Jmol._jsGetXY(canvas, ev);
			if(!xym) return false;
			if (!canvas.isDragging)
				xym[2] = 0;
			canvas.applet._processEvent((canvas.isDragging ? 506 : 503), xym); // J.api.Event.MOUSE_DRAG : J.api.Event.MOUSE_MOVE
			return false;
		}
		
		Jmol.$bind(canvas, 'DOMMouseScroll mousewheel', function(ev) { // Zoom
			ev.stopPropagation();
			ev.preventDefault();
			// Webkit or Firefox
			canvas.isDragging = false;
			var oe = ev.originalEvent;
			var scroll = (oe.detail ? oe.detail : (Jmol.featureDetection.os == "mac" ? 1 : -1) * oe.wheelDelta); // Mac and PC are reverse; but 
			var modifiers = Jmol._jsGetMouseModifiers(ev);
			canvas.applet._processEvent(-1,[scroll < 0 ? -1 : 1,0,modifiers]);
			return false;
		});

		// context menu is fired on mouse down, not up, and it's handled already anyway.

		Jmol.$bind(canvas, "contextmenu", function() {return false;});

		Jmol.$bind(canvas, 'mouseout', function(ev) {
			if (canvas.applet._applet)
				canvas.applet._applet.viewer.startHoverWatcher(false);
			//canvas.isDragging = false;
			var xym = Jmol._jsGetXY(canvas, ev);
			if (!xym)
				return false;
			//canvas.applet._processEvent(502, xym);//J.api.Event.MOUSE_UP
			//canvas.applet._processEvent(505, xym);//J.api.Event.MOUSE_EXITED
			return false;
		});

		Jmol.$bind(canvas, 'mouseenter', function(ev) {
			if (canvas.applet._applet)
				canvas.applet._applet.viewer.startHoverWatcher(true);
			if (ev.buttons === 0 || ev.which === 0) {
				canvas.isDragging = false;
				var xym = Jmol._jsGetXY(canvas, ev);
				if (!xym)
					return false;
				canvas.applet._processEvent(504, xym);//J.api.Event.MOUSE_ENTERED	
				canvas.applet._processEvent(502, xym);//J.api.Event.MOUSE_UP
				return false;
			}
		});

	Jmol.$bind(canvas, 'mousemoveoutjsmol', function(evspecial, target, ev) {
		if (canvas == Jmol._mouseOwner && canvas.isDragging) {
			return Jmol._drag(canvas, ev);
		}
	});

		if (canvas.applet._is2D)
			Jmol.$resize(function() {
				if (!canvas.applet)
					return;
				canvas.applet._resize();
			});
 
		Jmol.$bind('body', 'mouseup touchend', function(ev) {
			if (canvas.applet)
				canvas.isDragging = false;
			Jmol._setMouseOwner(null);
		});

	}

	Jmol._jsUnsetMouse = function(canvas) {
		canvas.applet = null;
		Jmol.$bind(canvas, 'mousedown touchstart mousemove touchmove mouseup touchend DOMMouseScroll mousewheel contextmenu mouseout mouseenter', null);
		Jmol._setMouseOwner(null);
	}


////// Jmol.Swing interface  for Javascript implementation of Swing dialogs and menus

Jmol.Swing = {
	// a static class
	count:0,
	menuInitialized:0,
	menuCounter:0,
	htDialogs:{}
};

(function(Swing) {

SwingController = Swing; // see javajs.api.SwingController

Swing.setDraggable = function(Obj) {
	
	var proto = Obj.prototype;
	if (proto.setContainer)
		return;
	
	// for menus, console, and 
	proto.setContainer = function(container) {
		this.container = container;
		container.obj = this;
		this.isDragging = false;
		this.ignoreMouse = false;
		var me = this;
		container.bind('mousedown touchstart', function(ev) {
			if (me.ignoreMouse) {
				me.ignoreMouse = false;
				return true;
			}
			Jmol._setMouseOwner(me, true);
			me.isDragging = true;
			me.pageX = ev.pageX;
			me.pageY = ev.pageY;
			return false;
		});
		container.bind('mousemove touchmove', function(ev) {
			if (me.isDragging && Jmol._mouseOwner == me) {
				me.mouseMove(ev);
				return false;
			}
		});
		container.bind('mouseup touchend', function(ev) {
			me.mouseUp(ev);
			Jmol._setMouseOwner(null);
		});
	};

	proto.mouseUp = function(ev) {
		if (this.isDragging && Jmol._mouseOwner == this) {
			this.pageX0 += (ev.pageX - this.pageX);
			this.pageY0 += (ev.pageY - this.pageY);
			this.isDragging = false;
			return false;
		}
		Jmol._setMouseOwner(null);
	}

	proto.setPosition = function() {
		if (Jmol._mousePageX === null) {
			var id = this.applet._id + "_" + (this.applet._is2D ? "canvas2d" : "canvas");
			var offsets = Jmol.$offset(id);
			Jmol._mousePageX = offsets.left;
			Jmol._mousePageY = offsets.top;
		}
		this.pageX0 = Jmol._mousePageX;
		this.pageY0 = Jmol._mousePageY;
		var pos = { top: Jmol._mousePageY + 'px', left: Jmol._mousePageX + 'px' };
		this.container.css(pos);
	};

	proto.mouseMove = function(ev) {
		if (this.isDragging && Jmol._mouseOwner == this) {
			this.timestamp = System.currentTimeMillis(); // used for menu closure
			var x = this.pageX0 + (ev.pageX - this.pageX);
			var y = this.pageY0 + (ev.pageY - this.pageY);
			this.container.css({ top: y + 'px', left: x + 'px' })
		}
	};

	proto.dragBind = function(isBind) {
		this.applet._ignoreMouse = !isBind;
		this.container.unbind('mousemoveoutjsmol');
		this.container.unbind('touchmoveoutjsmol');
		this.container.unbind('mouseupoutjsmol');
		this.container.unbind('touchendoutjsmol');
		Jmol._setMouseOwner(null);
		if (isBind) {
			var me = this;
			this.container.bind('mousemoveoutjsmol touchmoveoutjsmol', function(evspecial, target, ev) {
				me.mouseMove(ev);
			});
			this.container.bind('mouseupoutjsmol touchendoutjsmol', function(evspecial, target, ev) {
				me.mouseUp(ev);
			});
		}
	};
}

// Dialog //

Swing.JSDialog = function () {
}

Swing.setDraggable(Swing.JSDialog);

///// calls from javajs and other Java-derived packages /////

Swing.getScreenDimensions = function(d) {
	d.width = $(window).width();
	d.height = $(window).height();
}

Swing.dispose = function(dialog) {
	Jmol.$remove(dialog.id + "_mover");
	delete Swing.htDialogs[dialog.id]
	dialog.container.obj.dragBind(false);
//  var btns = $("#" + dialog.id + " *[id^='J']"); // add descendents with id starting with "J"
//  for (var i = btns.length; --i >= 0;)
//    delete Dialog.htDialogs[btns[i].id]
	//System.out.println("JSmolCore.js: dispose " + dialog.id)
}
 
Swing.register = function(dialog, type) {
	dialog.id = type + (++Swing.count);
	Swing.htDialogs[dialog.id] = dialog;
	//System.out.println("JSmolCore.js: register " + dialog.id)

}

Swing.setDialog = function(dialog) {
	Jmol._setMouseOwner(null);
	Jmol.$remove(dialog.id);
	//System.out.println("removed " + dialog.id)
	var id = dialog.id + "_mover";
	var container = Jmol._$(id);
	var jd;
	//System.out.println("JSmolCore.js: setDialog " + dialog.id);
	if (container[0]) {
		container.html(dialog.html);
		jd = container[0].jd;
	} else {
		Jmol.$after("body","<div id='" + id + "' style='position:absolute;left:0px;top:0px;'>" + dialog.html + "</div>");
		var jd = new Swing.JSDialog();
		container = Jmol._$(id);
		dialog.container = container;
		jd.applet = dialog.manager.vwr.applet;
		jd.setContainer(container);
		jd.dialog = dialog;
		jd.setPosition();  
		jd.dragBind(true);
		container[0].jd = jd; 
	}
	Jmol.$bind("#" + dialog.id + " .JButton", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JComboBox", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JCheckBox", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JTextField", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JTable", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JScrollPane", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JEditorPane", "mousedown touchstart", function(event) { jd.ignoreMouse=true });

}
 
Swing.setSelected = function(chk) {
 Jmol.$prop(chk.id, 'checked', !!chk.selected);
}

Swing.setSelectedIndex = function(cmb) {
 Jmol.$prop(cmb.id, 'selectedIndex', cmb.selectedIndex);
}

Swing.setText = function(btn) {
 Jmol.$prop(btn.id, 'value', btn.text);
}

Swing.setVisible = function(c) {
	Jmol.$setVisible(c.id, c.visible);
}

Swing.setEnabled = function(c) {
	Jmol.$setEnabled(c.id, c.enabled);
}

/// callbacks from the HTML elements ////
 
Swing.click = function(element, keyEvent) {
	var component = Swing.htDialogs[element.id];
	if (component) {
		//System.out.println("click " + element + " " + component)
		var info = component.toString();
		// table cells will have an id but are not registered
		if (info.indexOf("JCheck") >= 0) {
				component.selected = element.checked;
		} else if (info.indexOf("JCombo") >= 0) {
			component.selectedIndex = element.selectedIndex;
		} else if (component.text != null) {  // JButton, JTextField
			component.text = element.value;
			if (keyEvent && ((keyEvent.charCode || keyEvent.keyCode) != 13))
				return;
		}    
	}
	var dialog = Swing.htDialogs[Jmol.$getAncestorDiv(element.id, "JDialog").id];
	var key = (component ? component.name :  dialog.registryKey + "/" + element.id);
	//System.out.println("JSmolCore.js: click " + key); 
	dialog.manager.actionPerformed(key);
}

Swing.setFront = function(dialog) {
  var applet = dialog.manager.vwr.applet;
	if (dialog.zIndex != Jmol._getZ(applet, "dialog"))
	 dialog.zIndex = Jmol._incrZ(applet, "dialog");
	dialog.container && ((dialog.container[0] || dialog.container).style.zIndex = dialog.zIndex);
}

Swing.hideMenus = function(applet) {
	// called from LEFT-DOWN mouse event
	var menus = applet._menus;
	if (menus)
		for (var i in menus)
			if (menus[i].visible)
				Swing.hideMenu(menus[i]);
}

Swing.windowClosing = function(element) {
	var dialog = Swing.htDialogs[Jmol.$getAncestorDiv(element.id, "JDialog").id];
	if (dialog.registryKey) {
		//System.out.println("JSmolCore.js: windowClosing " + dialog.registryKey); 
		dialog.manager.processWindowClosing(dialog.registryKey);
	} else {
		//System.out.println("JSmolCore.js: windowClosing " + dialog.title); 
		dialog.dispose();
	}
}

})(Jmol.Swing);

Jmol._track = function(applet) {
	// this function inserts an iFrame that can be used to track your page's applet use. 
	// By default it tracks to a page at St. Olaf College, but you can change that. 
	// and you can use
	//
	// delete Jmol._tracker
	//
	// yourself to not have you page execute this 
	//
	if (Jmol._tracker){
		try {  
			var url = Jmol._tracker + "&applet=" + applet._jmolType + "&version=" + Jmol._version 
				+ "&appver=" + self.___JmolVersion + "&url=" + encodeURIComponent(document.location.href);
			var s = '<iframe style="display:none" width="0" height="0" frameborder="0" tabindex="-1" src="' + url + '"></iframe>'
			Jmol.$after("body", s);
		} catch (e) {
			// ignore
		}
		delete Jmol._tracker;
	}
	return applet;
}

Jmol.getProfile = function() {
	window["j2s.doProfile"] = true;
	if (self.Clazz && self.JSON) {
		Clazz._profile || (Clazz._profile = {});
		return Clazz.getProfile();
	}
}

Jmol._getInChIKey = function(applet, data) {
	if (data.indexOf("MOL=") >= 0)
		data = data.split("MOL=")[1].split("\"")[0];

}

Jmol._getAttr = function(s, a) {
	var pt = s.indexOf(a + "=");
	return (pt >= 0 && (pt = s.indexOf('"', pt)) >= 0 
		? s.substring(pt+1, s.indexOf('"', pt+1)) : null);
}

Jmol.User = {
	viewUpdatedCallback: null
}

Jmol.View = {

// The objective of Jmol.View is to coordinate
// asynchronous applet loading and atom/peak picking
// among any combination of Jmol, JME, and JSV.
// 
// basic element is a view object:
//   view = {
//     viewType1: viewRecord1,
//     viewType2: viewRecord2,
//     viewType3: viewRecord3
//   }
// where viewType is one of (Jmol, JME, JSV)
// and a viewRecord is an object
// with elements .chemID, .applet, .data, .script
//
// Jmol.View.sets is a list of cached views[0..n]
// for a given group of applet objects with common Info.viewSet
//
// Bob Hanson 1/22/2014 7:05:38 AM

	count: 0,
	applets: {},
	sets: {}
};

(function(View) {

// methods called from other modules have no "_" in their name

View.updateView = function(applet, Info, _View_updateView) {
	// Info.chemID, Info.data, possibly Info.viewID if no chemID
	// return from applet after asynchronous download of new data
	if (applet._viewSet == null)
		return;
	Info.chemID || (applet._searchQuery = null);
	Info.data || (Info.data = "N/A");
	Info.type = applet._viewType;
	if((applet._currentView = View.__findView(applet._viewSet, Info)) == null)
		applet._currentView = View.__createViewSet(applet._viewSet, Info.chemID, Info.viewID || Info.chemID);
	applet._currentView[Info.type].data = Info.data;
	applet._currentView[Info.type].smiles = applet._getSmiles();
	if (Jmol.User.viewUpdatedCallback)
		Jmol.User.viewUpdatedCallback(applet, "updateView");
	View.__setView(applet._currentView, applet, false);
}

View.updateFromSync = function(applet, msg) {
	applet._updateMsg = msg;
	var id = Jmol._getAttr(msg, "sourceID") || Jmol._getAttr(msg, "file");
	if (!id)
		return;
	var view = View.__findView(applet._viewSet, {viewID:id});
	if (view == null)
		return Jmol.updateView(applet, msg); // JSV has been updated internally
	if (view != applet._currentView)
		View.__setView(view, applet, true);
	var A = ((id = Jmol._getAttr(msg, "atoms")) && msg.indexOf("selectionhalos ON") >= 0  
		? eval("[" + id + "]") : []);
	setTimeout(function(){if (applet._currentView == view)View.updateAtomPick(applet, A)}, 10);
	if (Jmol.User.viewUpdatedCallback)
		Jmol.User.viewUpdatedCallback(applet, "updateFromSync");
}

View.updateAtomPick = function(applet, A, _View_updateAtomPick) {
	var view = applet._currentView;
	if (view == null)
		return;
	for (var viewType in view) {
		if (viewType != "info" && view[viewType].applet != applet)
			view[viewType].applet._updateAtomPick(A);
	}
	if (Jmol.User.viewUpdatedCallback)
		Jmol.User.viewUpdatedCallback(applet, "updateAtomPick");
}

View.dumpViews = function(setID) {
	var views = View.sets[setID];
	if (!views)
	  return;
	var s = "View set " + setID + ":\n";
	var applets = View.applets[setID];
	for (var i in applets)
		s += "\napplet " + applets[i]._id
			+ " currentView=" + (applets[i]._currentView ? applets[i]._currentView.info.viewID : null);
	for (var i = views.length; --i >= 0;) {
		var view = views[i];
		s += "\n\n<b>view=" + i 
			+ " viewID=" + view.info.viewID 
			+ " chemID=" + view.info.chemID + "</b>\n"
		var v;
		for (var viewType in view) 
			if (viewType != "info")
				s += "\nview=" + i + " type=" + viewType + " applet=" 
					+ ((v = view[viewType]).applet ? v.applet._id : null) 
					+ " SMILES=" + v.smiles + "\n"
					+ " atomMap=" + JSON.stringify(v.atomMap) + "\n"
					+ " data=\n" + v.data + "\n"
	}
	return s
}


// methods starting with "__" are "private" to JSmolCore.js

View.__init = function(applet) {
  var set = applet._viewSet;
	var a = View.applets;
	a[set] || (a[set] = {});
	a[set][applet._viewType] = applet;
}

View.__findView = function(set, Info) {
	var views = View.sets[set];
	if (views == null)
		views = View.sets[set] = [];
	for (var i = views.length; --i >= 0;) {
		var view = views[i];
		if (Info.viewID) {
			if (view.info.viewID == Info.viewID)
				return view;
		} else if (Info.chemID != null && Info.chemID == view.info.chemID) {
				return view;
		} else {
			for (var viewType in view) { 
				if (viewType != "info") {
					if (Info.data != null && view[viewType].data != null ? Info.data == view[viewType].data 
						: Info.type == viewType)
							return view;
				}
			}
		}
	}
	return null;  
}

View.__createViewSet = function(set, chemID, viewID, _createViewSet) {
	View.count++;
	var view = {info:{chemID: chemID, viewID: viewID || "model_" + View.count}};

	for (var id in Jmol._applets) {
		var a = Jmol._applets[id];
		if (a._viewSet == set)
			view[a._viewType] = {applet:a, data: null};
	}
	View.sets[set].push(view);
	return view;
}

View.__setView = function(view, applet, isSwitch, _setView) {
	// called from Jmol._searchMol and Jmol.View.setCurrentView 
	// to notify the applets in the set that there may be new data for them
	// skip the originating applet itself and cases where the data has not changed.
	// stop at first null data, because that will initiate some sort of asynchronous
	// call that will be back here afterward.

	for (var viewType in view) {
			if (viewType == "info") 
				continue;
		var rec = view[viewType];
		var a = rec.applet;
		var modified = (isSwitch || a != null && a._molData == "<modified>");

		if (a == null || a == applet && !modified)
			continue; // may be a mol3d required by JSV but not having a corresponding applet
		var wasNull = (rec.data == null);
		var haveView = (a._currentView != null);
		a._currentView = view; 
		if (haveView && view[viewType].data == rec.data && !wasNull & !modified)
			continue;
		a._loadModelFromView(view);
		if (wasNull)
			break;
	}

	// Either all are taken care of or one was null,
	// in which case we have started an asynchronous
	// process to get the data, and we can quit here.
	// In either case, we are done.
}

}) (Jmol.View);

Jmol.Cache = {fileCache: {}};

Jmol.Cache.get = function(filename) {
	return Jmol.Cache.fileCache[filename];
}

Jmol.Cache.put = function(filename, data) {
  Jmol.Cache.fileCache[filename] = data;
}

	Jmol.Cache.setDragDrop = function(me) {
		Jmol.$appEvent(me, "appletdiv", "dragover", function(e) {
			e = e.originalEvent;
			e.stopPropagation();
			e.preventDefault();
			e.dataTransfer.dropEffect = 'copy';
		});
		Jmol.$appEvent(me, "appletdiv", "drop", function(e) {
			var e = e.originalEvent;
			e.stopPropagation();
			e.preventDefault();
			var file = ev.dataTransfer.files[0];
			if (file == null) {
				// FF and Chrome will drop an image here
				// but it will be only a URL, not an actual file. 
				try {
				  file = "" + ev.dataTransfer.getData("text");
				  if (file.indexOf("file:/") == 0 || file.indexOf("http:/") == 0) {
				  	me._script("load \"" + file + "\"");
				  	return;
			  	}
				} catch(e) {
				  return;
				}
			  // some other format
			  return;
			}
			// MSIE will drop an image this way, though, and load it!
			var reader = new FileReader();
			reader.onloadend = function(evt) {
				if (evt.target.readyState == FileReader.DONE) {
					var cacheName = "cache://DROP_" + file.name;
					var bytes = Jmol._toBytes(evt.target.result);
					me._applet.viewer.cacheFileByName("cache://DROP_*",false);
					if (me._viewType == "JSV" || cacheName.endsWith(".jdx")) // shared by Jmol and JSV
						Jmol.Cache.put(cacheName, bytes);
					else
						me._applet.viewer.cachePut(cacheName, bytes);
					var xym = Jmol._jsGetXY(me._canvas, e);
					if(xym && (!me._applet.viewer.setStatusDragDropped || me._applet.viewer.setStatusDragDropped(0, xym[0], xym[1], cacheName))) {
						me._applet.viewer.openFileAsyncSpecial(cacheName, 1);
					}
				}
			};
			reader.readAsArrayBuffer(file);
		});
	}
  
})(Jmol, jQuery);
Jmol._debugCode = false;
// JSmolTM.js -- JSmol TwirlyMol implementation
//
// TM for "TwirlyMol" -- which this is derived from.
// Can't say there is much of TwirlyMol left in here, as
// we are not using dojo, there are no shadows, and I use a simpler
// rendering idea. Nonetheless, TwirlyMol is the inspiration.
//
// author: Bob Hanson, hansonr@stolaf.edu	3/23/2012

// This library requires
//
//	JSmoljQuery.js
//	JSmolCore.js
//
// and is combined with those in JSmol.lite.jq.js
// or packaged without jQuery -- because you are using it already -- in JSmol.lite.js
//


Jmol._grabberOptions = [
	["$", "NCI(small molecules)"],
	[":", "PubChem(small molecules)"]
];

Jmol.say = function(msg) {
	alert(msg);
}

Jmol._TMApplet = function(id, Info, checkOnly){
		this._uniqueId = ("" + Math.random()).substring(3);
		this._id = id;
		this._is2D = true;
		this._isJava = false;
		this._ready = true;
		this._mouseDown = false;
		this._jmolType = "Jmol._Canvas2D (TwirlyMol)";
		if (checkOnly)
			return this;
		this._createCanvas(id, Info);
		return this;
}

Jmol._TMApplet._getApplet = function(id, Info, checkOnly) {
	if (!Jmol.featureDetection.allowHTML5)
		return null;
				checkOnly || (checkOnly = false);
	Info || (Info = {});
	var DefaultInfo = {
		color: "#FFFFFF",
 		width: 300,
		height: 300,
		addSelectionOptions: false,
		serverURL: "http://your.server.here/jsmol.php",
		defaultModel: "",
		readyFunction: null,
		use: "HTML5",
		bondWidth: 5,
		shadeAtoms: false,
		zoomScaling: 1.5,
		pinchScaling: 2.0,
		mouseDragFactor: 0.5,
		touchDragFactor: 0.15,
		multipleBondSpacing: 4,
		spinRateX: 0.0,
		spinRateY: 0.5,
		spinFPS: 20,
		spin: false,
		noscript: true,
		debug: false
	};	 
	Jmol._addDefaultInfo(Info, DefaultInfo);
	var applet = new Jmol._TMApplet(id, Info, checkOnly);
	if (checkOnly)
		return applet;
	return Jmol._registerApplet(id, applet);
}

Jmol.getTMApplet = Jmol._TMApplet._getApplet;

;(function(tm){

tm._CPK = ["#FF1493", "#FFFFFF", "#D9FFFF", 
	"#CC80FF", "#C2FF00", "#FFB5B5", "#909090", "#3050F8", "#FF0D0D", "#90E050", "#B3E3F5", 
	"#AB5CF2", "#8AFF00", "#BFA6A6", "#F0C8A0", "#FF8000", "#FFFF30", "#1FF01F", "#80D1E3",
	"#8F40D4", "#3DFF00", "#E6E6E6", "#BFC2C7", "#A6A6AB", "#8A99C7", "#9C7AC7", "#E06633", 
		"#F090A0", "#50D050", "#C88033", "#7D80B0", "#C28F8F", "#668F8F", "#BD80E3", "#FFA100", "#A62929", "#5CB8D1", 
	"#702EB0", "#00FF00", "#94FFFF", "#94E0E0", "#73C2C9", "#54B5B5", "#3B9E9E", "#248F8F", 
		"#0A7D8C", "#006985", "#C0C0C0", "#FFD98F", "#A67573", "#668080", "#9E63B5", "#D47A00", "#940094", "#429EB0", 
	"#57178F", "#00C900"] // through barium
	//  "#70D4FF", "#FFFFC7", "#D9FFC7", "#C7FFC7", "#A3FFC7", "#8FFFC7", "#61FFC7", "#45FFC7", "#30FFC7", "#1FFFC7", "#00FF9C", "#00E675", "#00D452", "#00BF38", "#00AB24", "#4DC2FF", "#4DA6FF", "#2194D6", "#267DAB", "#266696", "#175487", "#D0D0E0", "#FFD123", "#B8B8D0", "#A6544D", "#575961", "#9E4FB5", "#AB5C00", "#754F45", "#428296", "#420066", "#007D00", "#70ABFA", "#00BAFF", "#00A1FF", "#008FFF", "#0080FF", "#006BFF", "#545CF2", "#785CE3", "#8A4FE3", "#A136D4", "#B31FD4", "#B31FBA", "#B30DA6", "#BD0D87", "#C70066", "#CC0059", "#D1004F", "#D90045", "#E00038", "#E6002E", "#EB0026"]; 
tm._elem = ['X', 'H', 'He', 
	'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 
	'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 
	'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
	'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 
	'Cs', 'Ba', 
		'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 
			 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb','Bi', 'Po', 'At', 'Rn', 
	'Fr', 'Ra', 
		'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es'];
tm._elemNo = {};

;(function(proto) {

// A user-available function

proto.spin = function(TF) {
	this.__Info.spin = TF;
	this._spin(TF);
}

proto._spin = function(TF) {
	if (this._spinThread)clearTimeout(this._spinThread);
	if (this.spinFPS == 0 || this.spinRateX == 0 && this.spinRateY == 0)
		TF = false;
	if (!TF) return;
	var me = this;
	var delay = 1000 / this.spinFPS;
	if (!this._mouseDown) {
		this._rotate(this.spinRateY, this.spinRateX);
		this._draw();
	}
	this._spinThread = setTimeout(function(){me._spin(true)}, delay);
}

proto._initParams = function() {
	this.zoom = this.__Info.defaultZoom || 100;
	this.doSpin = this.__Info.spin || false;
	this.center2D = [this._canvas.width/2, this._canvas.height/2, 0];
	this._getCenterAndRadius();
	this.rotation = new tm.M3();
	this.shadeAtoms = false;
	this._setParams();
}

proto._setParams = function() {
	this.bondWidth = this.__Info.bondWidth || 5;
	this.zoomScaling = this.__Info.zoomScaling || 1.5;
	this.pinchScaling = this.__Info.pinchScaling || 1;
	this.mouseDragFactor = this.__Info.mouseDragFactor || 0.5;
	this.touchDragFactor = this.__Info.touchDragFactor || 0.15;
	this.multipleBondSpacing = this.__Info.multipleBondSpacing || 4;
	this.spinRateX = this.__Info.spinRateX || 0;
	this.spinRateY = this.__Info.spinRateY || 0;
	this.spinFPS = this.__Info.spinFPS || 0;
	var p = this.shadeAtoms;
	this.shadeAtoms = this.__Info.shadeAtoms || false;
	if (this.shadeAtoms && !p)
		this._setAtomShades();
}

proto._setAtomShades = function() {
	if (!this.atoms)
		return;
	for (var i = this.atoms.length; --i >= 0;)
		this.atoms[i].color50 = this._getColor(this.atoms[i].color, 0.5);
}

proto._createCanvas = function(id, Info) {
	Jmol._setObject(this, id, Info);
	this._color = this._color.replace(/0x/,"#");
	var t = Jmol._getWrapper(this, true);
	if (Jmol._document) {
		Jmol._documentWrite(t);
		this._createCanvas2d(false);				
		t = "";
	} else {
		t += '<script type="text/javascript">' + id + '._createCanvas2d(false)</script>';
	}
	t += Jmol._getWrapper(this, false);
	if (Info.addSelectionOptions)
		t += Jmol._getGrabberOptions(this, "");
	if (Jmol._debugAlert && !Jmol._document)
		alert(t);    
	this._code = Jmol._documentWrite(t);
};

proto._readyCallback = function(id, fullid, isReady, applet) {
	if (!isReady)
		return; // ignore -- page is closing
	var me = this;
	Jmol._setDestroy(me);
	me._ready = true;
	var script = me._readyScript;
	me._applet = applet;
	if (me._defaultModel) 
		me._search(me._defaultModel, (script ? ";" + script : ""));
	else if (me._src)
		me._loadFile(me._src);
	me._showInfo(true);
	me._showInfo(false);
	me._readyFunction && me._readyFunction(me);
	Jmol._setReady(this);
};

proto._createCanvas2d = function(doReplace) {
	var canvas = document.createElement( 'canvas' );
	var container = Jmol.$(this, "appletdiv");
	if (doReplace) {
		container[0].removeChild(this._canvas);
		Jmol._jsUnsetMouse(this._canvas);
	}
	var w = Math.round(container.width());
	var h = Math.round(container.height());
	canvas.applet = this;
	this._canvas = canvas;
	canvas.style.width = "100%";
	canvas.style.height = "100%";
	canvas.width = w;
	canvas.height = h; // w and h used in setScreenDimension
	canvas.id = this._id + "_canvas2d";
	container.append(canvas);
	setTimeout(this._id+"._start()",10);
}

proto._start = function() {
	Jmol._jsSetMouse(this._canvas);
	if (this._defaultModel)
		Jmol._search(this, this._defaultModel);
	else if (this._src)
		this._loadFile(this._src);
	this._showInfo(true);
	this._showInfo(false);
}      

proto._search = function(query, script){
	Jmol._search(this, query, script);
}

proto._searchDatabase = function(query, database, script){
	this._showInfo(false);
	if (query.indexOf("?") >= 0) {
		Jmol._getInfoFromDatabase(this, database, query.split("?")[0]);
		return;
	}
	var dm = database + query; 
	this._loadFile(dm, script);
}

proto.__loadModel = function(mol) {
	this._spin(false);
	if (mol == "''")mol = this._mol;
	if (!mol) {
		alert("No model data.");
		return;
	}
	this._parse(mol);
	this._initParams();
	this._draw();
	this._showInfo(false);
	if (this.doSpin) 
		this._spin(true);
};

proto._showInfo = function(tf) {
	Jmol.$html(Jmol.$(this, "infoheaderspan"), this._infoHeader);
	if (this._info)
		Jmol.$html(Jmol.$(this, "infodiv"), this._info);
	if ((!this._isInfoVisible) == (!tf))
		return;
	this._isInfoVisible = tf;
	Jmol.$setVisible(Jmol.$(this, "infotablediv"), tf);
	Jmol.$setVisible(Jmol.$(this, "infoheaderdiv"), tf);
	this._show(!tf);
}

proto._show = function(tf) {
	Jmol.$setVisible(Jmol.$(this, "appletdiv"), tf);
	if (tf)
		this._draw();
}

proto._resize = function() {
	var s = "__resizeTimeout_" + this._id;
	// only at end
	if (Jmol[s])
		clearTimeout(Jmol[s]);
	var me = this;
	Jmol[s] = setTimeout(function() {me._draw();Jmol[s]=null}, 100);
}

proto._canScript = function(script) {return false};

proto._loadFile = function(fileName){
	this._showInfo(false);
	this._thisJmolModel = fileName;
	var app = this;
	Jmol._loadFileData(this, fileName, function(data){app.__loadModel(data)});
}

proto._processEvent = function(type, xym) {
	switch(type) {
	case 502: // mouse_up
	case 503: // mouse_move
		this._mouseDown = false;
	default:
		return;
	case 501: // mouse_down
		this._mouseX = xym[0];
		this._mouseY = xym[1];
		this._mouseDown = true;
		return;
	case 506: // mouse_drag
		var dx = xym[0] - this._mouseX;
		var dy = xym[1] - this._mouseY;
		dx = (dx < 0 ? -1 : dx > 0 ? 1 : 0)
		dy = (dy < 0 ? -1 : dy > 0 ? 1 : 0)

		switch (xym[2]) {
		case 17: // left-shift
			this._zoomBy(dy);
			break;
		case 24: // left-alt
			this.center2D[0] += dx;
			this.center2D[1] += dy;
			break;
		default:
			this._rotate(dx, dy);
			break;
		}
		this._mouseX = xym[0];
		this._mouseY = xym[1];
		break;
	case -1:  // mouse_scroll
		this._zoomBy(xym[1])
		break;
	}
	this._draw();
}

proto._processGesture = function (touches) {
	if (touches[0].length < 2) return;
	var t1 = touches[0];
	var t2 = touches[1];
	var t10 = t1[0];
	var t11 = t1[t2.length - 1];
	var x10 = t10[0];
	var x11 = t11[0];
	var dx1 = x11 - x10;
	var y10 = t10[1];
	var y11 = t11[1];
	var dy1 = y11 - y10;
	var v1 = [dx1, dy1, 0];
	var d1 = tm._length(v1);
	var t20 = t2[0];
	var t21 = t2[t2.length - 1];
	var x20 = t20[0];
	var x21 = t21[0];
	var dx2 = x21 - x20;
	var y20 = t20[1];
	var y21 = t21[1];
	var dy2 = y21 - y20;
	var v2 = [dx2, dy2, 0];
	var d2 = tm._length(v2);
	if (d1 < 3 || d2 < 3) return;
	tm._scale(v1, 1/d1);
	tm._scale(v2, 1/d2);
	var cos12 = tm._dot(v1,v2);
	if (cos12 > 0.8) {
		var deltaX = Math.round(x11 - t1[t1.length - 2][0]);
		var deltaY = Math.round(y11 - t1[t1.length - 2][1]);
		this.center2D[0] += deltaX;
		this.center2D[1] += deltaY;
	} else if (cos12 < -0.8) {
		v1 = [x20 - x10, y20 - y10, 0]
		v2 = [x21 - x11, y21 - y11, 0];
		var dx = tm._length(v2) - tm._length(v1);
		this.zoom += (dx < 0 ? -1 : 1) * this.pinchScaling;
	}
	this._draw();
}

proto._zoomBy = function(dz) {
	this.zoom += dz * this.zoomScaling;
	if (this.zoom < 5) this.zoom = 5;
	if (this.zoom > 500) this.zoom = 500;
}

proto._getRotationScaling = function() {
	return [1 / Math.max(this._canvas.width, 500) * this.mouseDragFactor * tm.deg_per_radian, 
		1 / Math.max(this._canvas.height, 500) * this.mouseDragFactor * tm.deg_per_radian];
}

proto._rotate = function(dx, dy) {
	var f = this._getRotationScaling();
	if (dy) {
		tm._m3.rotX(dy * f[1]);
		this.rotation.mul(tm._m3);
	}  
	if (dx) {
		tm._m3.rotY(dx * f[0]);
		this.rotation.mul(tm._m3);
	}
}

proto._getCenterAndRadius = function() {
	var p = this;
	var size = Math.min(p._canvas.width, p._canvas.height);
	var min = [999999, 999999, 999999];
	var max = [-999999, -999999, -999999];
	for (var i= p.atoms.length; --i >= 0;)
		for (var j = 3; --j >= 0;) {
			var v = p.atoms[i].xyz[j];
			if (v < min[j]) min[j]=v;
			if (v > max[j]) max[j]=v;
		}
	var c = [0, 0, 0];
	for (var j = 3; --j >= 0;)
		c[j] = (max[j] + min[j])/2;
	var r = 0;
	for (var i= p.atoms.length; --i >= 0;) {
		var a = p.atoms[i];
		var d = tm._distance(a.xyz, c) + (a.elemNo == 1 ? 1 : 1.5);//a.radius;
		if (d > r)
			r = d;
	}
	this.center = c;
	this.modelRadius = (r == 0 ? 10 : r);
}

proto._setScreenCoords = function() {
	var c = this.center;
	var mat = this.rotation;
	var sc = this.center2D;
	var f = Math.min(this._canvas.width, this._canvas.height) / 2 / this.modelRadius * this.zoom / 100;
	for (var i=this.atoms.length; --i >= 0;) {
		var a = this.atoms[i];
		this._transform(a.xyz, a.sxyz, c, mat, f, sc);
		a.srad = f * a.radius;
	}
	for (var i=this.bonds.length; --i >= 0;) {
		var b = this.bonds[i];
		this._transform(b.xyz, b.sxyz, c, mat, f, sc);
		this._transform(b.pts[0], b.spts[0], c, mat, f, sc);
		this._transform(b.pts[1], b.spts[1], c, mat, f, sc);
	}
}

proto._transform = function(pt, spt, c, mat, f, sc) {
	var a = tm._temp1;
	var b = tm._temp2;
	for (var i = 3; --i >= 0;)
		a[i] = (pt[i] - c[i]) * f;
	mat.transform(a, b);
	b[1] = -b[1];
	b[2] = -b[2];
	for (var i = 3; --i >= 0;)
		spt[i] = b[i] + sc[i];  
}


//  adaptable --- 

proto._parse = function(data) {
	// just MOL data for now
	this._parseSDF(data);
}

proto._parseSDF = function(sdf) {
	this._mol = sdf;
	if (!tm._elemNo["H"])
		for (var i = tm._elem.length; --i >= 0;)
			tm._elemNo[tm._elem[i]] = i;
	var lines = sdf.split("\n");
	var pt = 3;
	var line = lines[pt++];
	this.nAtoms = parseFloat(line.substring(0, 3));
	this.nBonds = parseFloat(line.substring(3, 6));
	this.atoms = Array(this.nAtoms);
	this.bonds = Array(this.nBonds);
	this.elements = Array(this.nAtoms + this.nBonds)
	var ept = 0;
	var x, y, z;
	for (var i = 0; i < this.nAtoms; i++) {
		line = lines[pt++];
		x = parseFloat(line.substring(0, 10));
		y = parseFloat(line.substring(10, 20));
		z = parseFloat(line.substring(20, 30));
		var e = line.substring(31, 34).replace(/\s+/g, '');
		var elemNo = tm._elemNo[e] || 0;
		var radius = (elemNo == 1 ? 0.23 : 0.35);
		var color = tm._CPK[elemNo] || tm._CPK[0];
		this.elements[ept++] = this.atoms[i] = {type:0, elemNo: elemNo, xyz:[x, y, z], sxyz:[0, 0, 0], radius: radius, color: color, color50: color};
	}
	for (var i = 0; i < this.nBonds; i++) {
		line = lines[pt++]
		var a1 = this.atoms[parseFloat(line.substring(0, 3)) - 1];
		var a2 = this.atoms[parseFloat(line.substring(3, 6)) - 1];
		var order = Math.abs(parseFloat(line.substring(6, 9)));
		switch (order) {
		case 4: // JmolAdapter.ORDER_AROMATIC;
		case 5: // JmolAdapter.ORDER_PARTIAL12;
		case 6: // JmolAdapter.ORDER_AROMATIC_SINGLE;
		case 8: // JmolAdapter.ORDER_PARTIAL01;
			order = 1;
			break;
		case 7: // JmolAdapter.ORDER_AROMATIC_DOUBLE;
			order = 2; 
			break;
		}
		var xyz = tm._getPointAlong(a1.xyz, a2.xyz, 0.5);
		var d = tm._distance(a1.xyz, xyz); 
		var p1 = (a1.radius < d ? tm._getPointAlong(a1.xyz, xyz, a1.radius/d) : [0, 0, 99999]);
		var p2 = (a2.radius < d ? tm._getPointAlong(a2.xyz, xyz, a2.radius/d) : [0, 0, 99999]);
		this.elements[ept++] = this.bonds[i] = {type:1, atoms: [a1, a2], xyz:xyz, pts:[p1,p2], sxyz:[0, 0, 0], spts:[[0,0,0],[0,0,0]], order: order, color: 0};
	}
	//alert(this.atoms.length + " atoms " + this.bonds.length + " bonds")
}

proto._getColor = function(c, f) {
	if (c == "#FFFFFF") c = "#D0D0D0"; 
	var n = parseInt("0x" + c.substring(1), 16);
	var r=(n >> 16) & 0xFF;
	var g=(n >> 8) & 0xFF;
	var b=(n) & 0xFF;
	if (r != 255)
		r += Math.floor((255 - r) * f);
	if (g != 255)
		g += Math.floor((255 - g) * f);
	if (b != 255)
		b += Math.floor((255 - b) * f);
	var s = "000000" + ((r << 16) | (g << 8) | b).toString(16);
	return "#"+s.substring(s.length - 6, s.length);
}

proto._draw = function() {
	if (!this.atoms)return;
	this._setParams();
	this._setScreenCoords();
	var es = this.elements;  
	es.sort(tm._zorder);
	var ctx=this._canvas.getContext("2d");
	ctx.fillStyle=this._color;
	ctx.fillRect(0,0,this._canvas.width, this._canvas.height);
	for (var i = es.length; --i >= 0;) {
		var e = es[i];
		if (e.type == 0)
			this._drawAtom(ctx, e);
		else
			this._drawBond(ctx, e);
	}
}

///////////////// drawing ///////////////////



proto._drawAtom = function(ctx, a) {
	ctx.beginPath();
	// context.fillStyle = "rgba(0, 0, 0, 1)";

	if (this.shadeAtoms) {
		var offset = a.srad / 4;
		var grad = ctx.createRadialGradient(a.sxyz[0]-offset,a.sxyz[1]-offset, offset, a.sxyz[0],a.sxyz[1], a.srad);
		grad.addColorStop(0.0, a.color50);
		grad.addColorStop(1.0, (a.color == "#FFFFFF" ? "#D0D0D0" : a.color));
		ctx.fillStyle = grad;  
		ctx.arc(a.sxyz[0],a.sxyz[1],a.srad,0,tm._pi2);
		ctx.fill();
	} else {
		ctx.fillStyle = a.color;
		ctx.arc(a.sxyz[0],a.sxyz[1],a.srad,0,tm._pi2);
		ctx.fill();
		ctx.strokeStyle = "#000000";
		ctx.lineWidth = 1;
		ctx.stroke();
	}
}

proto._drawBond = function(ctx, b) {
	if (b.pts[0][2] != 99999)
		this._drawBondLines(ctx, b, b.spts[0], b.atoms[0].color);       
	if (b.pts[1][2] != 99999)
		this._drawBondLines(ctx, b, b.spts[1], b.atoms[1].color); 
}

proto._drawBondLines = function(ctx, b, pt, ac) {
	var t = tm._temp;
	t.width = this.bondWidth;
	t.color = (b.color ? b.color : ac);
	if (b.order == 1) {
		t.pt1 = pt;
		t.pt2 = b.sxyz;
		this._drawLine(ctx, t);
		return;
	}
	t.pt1 = tm._temp1;
	t.pt2 = tm._temp2;
	t.pta = pt;
	t.ptb = b.sxyz;
	t.dx = t.ptb[0] - t.pta[0];
	t.dy = t.ptb[1] - t.pta[1];
	t.mag2d = Math.round(Math.sqrt(t.dx * t.dx + t.dy * t.dy));
	t.bondOrder = b.order;
	this._resetAxisCoordinates(t);
	while (t.bondOrder > 0) {
		this._drawLine(ctx, t);
		t.bondOrder--;
		this._stepAxisCoordinates(t);
	}
}

proto._drawLine = function(ctx, t) {
		ctx.beginPath();
		ctx.strokeStyle = t.color;
		ctx.moveTo(t.pt1[0], t.pt1[1]);
		ctx.lineTo(t.pt2[0], t.pt2[1]);
		ctx.lineWidth = t.width;
		ctx.stroke();
}

proto._resetAxisCoordinates = function(t) {
	var space = t.mag2d >> 3;
	if (this.multipleBondSpacing != 1) space *= this.multipleBondSpacing;
	var step = t.width + space;
	t.dxStep = Math.round(step * t.dy / t.mag2d);
	t.dyStep = Math.round(step * -t.dx / t.mag2d);
	t.pt1[0] = t.pta[0];
	t.pt1[1] = t.pta[1];
	t.pt2[0] = t.ptb[0];
	t.pt2[1] = t.ptb[1];
	var f = (t.bondOrder - 1);
	t.pt1[0] -= Math.round(t.dxStep * f / 2);
	t.pt1[1] -= Math.round(t.dyStep * f / 2);
	t.pt2[0] -= Math.round(t.dxStep * f / 2);
	t.pt2[1] -= Math.round(t.dyStep * f / 2);
}

proto._stepAxisCoordinates = function(t) {
	t.pt1[0] += t.dxStep;
	t.pt1[1] += t.dyStep;
	t.pt2[0] += t.dxStep;
	t.pt2[1] += t.dyStep;
}

})(tm.prototype)

tm._distance = function(a,b){
	var d = 0;
	for (var i = 3; --i >= 0;) {
		var dx = a[i] - b[i];
		d += dx * dx;
	}
	return Math.sqrt(d);
} 

tm._dot = function(a,b){
	var d = 0;
	for (var i = 3; --i >= 0;)
		d += a[i] * b[i];
	return d;
} 

tm._setPt = function(a,b){
	for (var i = 3; --i >= 0;)
		b[i] = a[i];
	return b;
} 

tm._scale = function(a, f){
	for (var i = 3; --i >= 0;)
		a[i] *= f;
} 

tm._length = function(a){
	var d = 0;
	for (var i = 3; --i >= 0;)
		d += a[i] * a[i];
	return Math.sqrt(d);
} 

tm._getPointAlong = function(a, b, f) {
	return [(b[0] - a[0]) * f + a[0], 
		(b[1] - a[1]) * f + a[1], 
		(b[2] - a[2]) * f + a[2]];
}
	 
tm._zorder = function(a, b) {
	var x = a.sxyz[2];
	var y = b.sxyz[2];
	return (x < y ? -1 : x > y ? 1 : 0);
}

tm.M3 = function() {
	this.m00 = 1.0;
	this.m01 = 0.0;
	this.m02 = 0.0;
	this.m10 = 0.0;
	this.m11 = 1.0;
	this.m12 = 0.0;
	this.m20 = 0.0;
	this.m21 = 0.0;
	this.m22 = 1.0;
}

tm.M3.prototype.set = function (m00, m01, m02, m10, m11, m12, m20, m21, m22) {
	this.m00 = m00;
	this.m01 = m01;
	this.m02 = m02;
	this.m10 = m10;
	this.m11 = m11;
	this.m12 = m12;
	this.m20 = m20;
	this.m21 = m21;
	this.m22 = m22;
}

tm.M3.prototype.transform = function(a, b) {
	b[0] = this.m00 * a[0] + this.m01 * a[1] + this.m02 * a[2];
	b[1] = this.m10 * a[0] + this.m11 * a[1] + this.m12 * a[2];
	b[2] = this.m20 * a[0] + this.m21 * a[1] + this.m22 * a[2];
}

tm.M3.prototype.rotX = function(angle) {
	var c = Math.cos(angle);
	var s = Math.sin(angle);
	this.m00 = 1.0;
	this.m01 = 0.0;
	this.m02 = 0.0;
	this.m10 = 0.0;
	this.m11 = c;
	this.m12 = -s;
	this.m20 = 0.0;
	this.m21 = s;
	this.m22 = c;
}

tm.M3.prototype.rotY = function(angle) {
	var c = Math.cos(angle);
	var s = Math.sin(angle);
	this.m00 = c;
	this.m01 = 0.0;
	this.m02 = s;
	this.m10 = 0.0;
	this.m11 = 1.0;
	this.m12 = 0.0;
	this.m20 = -s;
	this.m21 = 0.0;
	this.m22 = c;
}  

tm.M3.prototype.mul = function(m) {
	var m2 = this;
	this.set(m.m00 * m2.m00 + m.m01 * m2.m10 + m.m02 * m2.m20, m.m00 * m2.m01 + m.m01 * m2.m11 + m.m02 * m2.m21, m.m00 * m2.m02 + m.m01 * m2.m12 + m.m02 * m2.m22, m.m10 * m2.m00 + m.m11 * m2.m10 + m.m12 * m2.m20, m.m10 * m2.m01 + m.m11 * m2.m11 + m.m12 * m2.m21, m.m10 * m2.m02 + m.m11 * m2.m12 + m.m12 * m2.m22, m.m20 * m2.m00 + m.m21 * m2.m10 + m.m22 * m2.m20, m.m20 * m2.m01 + m.m21 * m2.m11 + m.m22 * m2.m21, m.m20 * m2.m02 + m.m21 * m2.m12 + m.m22 * m2.m22);
}

tm._pi2 = 2*Math.PI;
tm.deg_per_radian = 180 / Math.PI;
tm._z3 = [0,0,0];
tm._temp1 = [0,0,0];
tm._temp2 = [0,0,0];
tm._temp = {};
tm._m3 = new tm.M3();

})(Jmol._TMApplet);
