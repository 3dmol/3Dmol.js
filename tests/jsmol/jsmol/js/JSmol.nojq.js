___JmolDate="$Date: 2014-06-28 18:17:32 -0500 (Sat, 28 Jun 2014) $"
___fullJmolProperties="src/org/jmol/viewer/Jmol.properties"
___JmolVersion="14.2.2_2014.06.29"
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
// BH 5/8/2014 11:16:40 AM j2sPath starting with "/" fails to add idiom
// BH 1/16/2014 8:44:03 PM   Jmol.__execDelayMS = 100; // FF bug when loading into a tab that is not 
//                           immediately focused and not using jQuery for adding the applet and having  
//                           multiple applets.
// BH 12/6/2013 10:12:30 AM adding corejmoljsv.z.js
// BH 9/17/2013 10:18:40 AM  file transfer functions moved to JSmolCore 
// BH 3/5/2013 9:54:16 PM added support for a cover image: Info.coverImage, coverScript, coverTitle, deferApplet, deferUncover
// BH 1/3/2013 4:54:01 AM mouse binding should return false -- see d.bind(...), and d.bind("contextmenu") is not necessary

// JSmol.js -- Jmol pure JavaScript version
// author: Bob Hanson, hansonr@stolaf.edu	4/16/2012
// author: Takanori Nakane biochem_fan 6/12/2012

// This library requires prior inclusion of 

//  jQuery 
//	JSmoljQueryExt.js
//	JSmolCore.js
//  JSmolApplet.js
//  JSmolApi.js
//  j2sjmol.js    (Clazz and associated classes)

// these:
//
//  JSmolThree.js
//  JSmolGLmol.js
//  
//  are optional 

;(function (Jmol) {

	Jmol._coreFiles = []; // required for package.js


  ///////////////////
 // This section provides an asynchronous loading sequence
//

	Jmol.__execLog = [];
	Jmol.__execStack = [];
	Jmol.__execTimer = 0;
	Jmol.__coreSet = [];
	Jmol.__coreMore = [];
	Jmol.__execDelayMS = 100; // must be > 55 ms for FF

	Jmol.showExecLog = function() { return Jmol.__execLog.join("\n") }; 

	Jmol.__addExec = function(e) {
		var s = "JSmol load " + e[0]._id + " " + e[3];
		if (self.console)console.log(s + "...")
		Jmol.__execLog.push(s);   
		Jmol.__execStack.push(e);
	}

	Jmol.__addCoreFile = function(type, path, more) {
		if (Jmol.__coreSet.join("").indexOf(type) >= 0) return;
		Jmol.__coreSet.push(type);
		Jmol.__coreSet.sort();
		var f = Jmol.__coreSet.join("");
		Jmol._coreFiles = [path + "/core/core" + (f == "jmol" ? "" : f) + ".z.js" ];
		if (more && (more = more.split(" ")))
			for (var i = 0; i < more.length; i++)
				if (Jmol.__coreMore.join("").indexOf(more[i]) < 0)
					Jmol.__coreMore.push(path + "/core/core" + more[i] + ".z.js")
		for (var i = 0; i < Jmol.__coreMore.length; i++)
			Jmol._coreFiles.push(Jmol.__coreMore[i]);			 
	}      		

	Jmol.__nextExecution = function(trigger) {
		delete Jmol.__execTimer;
		var es = Jmol.__execStack;
		var e;
		while (es.length > 0 && (e = es[0])[4] == "done")
			es.shift();
		if (es.length == 0)
			return;
		if (!trigger) {
			setTimeout("Jmol.__nextExecution(true)",10)
			return;
		}
		e.push("done");
		var s = "JSmol exec " + e[0]._id + " " + e[3] + " " + e[2];
		if (self.console)console.log(s + " -- OK")
		Jmol.__execLog.push(s);
		e[1](e[0],e[2]);	
	};

	Jmol.__loadClazz = function(applet) {
		// problems with multiple applets?
		if (!Jmol.__clazzLoaded) {
			Jmol.__clazzLoaded = true;
			LoadClazz();
			if (applet._noMonitor)
				Clazz._LoaderProgressMonitor.showStatus = function() {}
			LoadClazz = null;

			Clazz._Loader.globalLoaded = function (file) {
			 // not really.... just nothing more yet to do yet
				Clazz._LoaderProgressMonitor.showStatus ("Application loaded.", true);
				if (!Jmol._debugCode || !Jmol.haveCore) {
					Jmol.haveCore = true;
					Jmol.__nextExecution();
				}
			};
			Clazz._Loader.packageClasspath ("java", null, true);
			return;
		}
		Jmol.__nextExecution();
	};

	Jmol.__loadClass = function(applet, javaClass) {
		Clazz._Loader.loadClass(javaClass, function() {Jmol.__nextExecution()});
	};

	Jmol._Canvas2D = function(id, Info, type, checkOnly){
		// type: Jmol or JSV
		this._uniqueId = ("" + Math.random()).substring(3);
		this._id = id;
		this._is2D = true;
		this._isJava = false;
		this._jmolType = "Jmol._Canvas2D (" + type + ")";
		switch (type) {
		case "Jmol":
			this._platform = "J.awtjs2d.Platform";
			break;
		case "JSV":
			this._isJSV = true;
			this._isLayered = true;
			this._platform = "JSV.awtjs2d.Platform";
			break;
		}
		if (checkOnly)
			return this;
		window[id] = this;
		this._createCanvas(id, Info);
		if (!Jmol._document || this._deferApplet)
			return this;
		this._init();
		return this;
	};

	Jmol._jsSetPrototype = function(proto) {
		proto._init = function() {
	 		this._setupJS();
			this._showInfo(true); 
			if (this._disableInitialConsole)
				this._showInfo(false);
		};

		proto._createCanvas = function(id, Info, glmol) {
			Jmol._setObject(this, id, Info);
			if (glmol) {
				this._GLmol = glmol;
		 		this._GLmol.applet = this;
				this._GLmol.id = this._id;
			}      
			var t = Jmol._getWrapper(this, true);
			if (this._deferApplet) {
			} else if (Jmol._document) {
				Jmol._documentWrite(t);
				this._newCanvas(false);				        
				t = "";
			} else {
				this._deferApplet = true;
				t += '<script type="text/javascript">' + id + '._cover(false)</script>';
			}
			t += Jmol._getWrapper(this, false);
			if (Info.addSelectionOptions)
				t += Jmol._getGrabberOptions(this);
			if (Jmol._debugAlert && !Jmol._document)
				alert(t);
			this._code = Jmol._documentWrite(t);
		};

		proto._cover = function (doCover) {
			if (doCover || !this._deferApplet) {
				this._displayCoverImage(doCover);
				return;
			}
			// uncovering UNMADE applet upon clicking image
			var s = (this._coverScript ? this._coverScript : "");
			this._coverScript = "";
			if (this._deferUncover)
				s += ";refresh;javascript " + this._id + "._displayCoverImage(false)";
			this._script(s, true);
			if (this._deferUncover && this._coverTitle == "activate 3D model")
				Jmol._getElement(this, "coverimage").title = "3D model is loading...";
			if (this._deferApplet)
				this._newCanvas(false);
			if (this._defaultModel)	
				Jmol._search(this, this._defaultModel);
			this._showInfo(false);
			if (!this._deferUncover)
				this._displayCoverImage(doCover);
			if (this._init)
				this._init();
		};

		proto._displayCoverImage = function(TF) {
			if (!this._coverImage || this._isCovered == TF) return;
			this._isCovered = TF;
			Jmol._getElement(this, "coverdiv").style.display = (TF ? "block" : "none");
		};

		proto._newCanvas = function(doReplace) {
			if (this._is2D)
				this._createCanvas2d(doReplace);
			else
				this._GLmol.create();
		};

		proto._createCanvas2d = function(doReplace) {
			var container = Jmol.$(this, "appletdiv");
			//if (doReplace) {
				try {
				container[0].removeChild(this._canvas);
				if (this._canvas.frontLayer)
					container[0].removeChild(this._canvas.frontLayer);
				if (this._canvas.rearLayer)
					container[0].removeChild(this._canvas.rearLayer);
				if (this._canvas.imageLayer)
					container[0].removeChild(this._canvas.imageLayer);
				Jmol._jsUnsetMouse(this._mouseInterface);
				} catch (e) {}
			//}
			var w = Math.round(container.width());
			var h = Math.round(container.height());
			var canvas = document.createElement( 'canvas' );
			canvas.applet = this;
			this._canvas = canvas;
			canvas.style.width = "100%";
			canvas.style.height = "100%";
			canvas.width = w;
			canvas.height = h; // w and h used in setScreenDimension
			canvas.id = this._id + "_canvas2d";
			container.append(canvas);
			Jmol._$(canvas.id).css({"z-index":Jmol._getZ(this, "main")});
			if (this._isLayered){
				var img = document.createElement("div");
				canvas.imageLayer = img;
				img.id = this._id + "_imagelayer";
				container.append(img);
				Jmol._$(img.id).css({zIndex:Jmol._getZ(this, "image"),position:"absolute",left:"0px",top:"0px", width:"0px", height:"0px", overflow:"hidden"});
				this._mouseInterface = this._getLayer("front", container, w, h, false);
				//this._getLayer("rear", container, w, h, true);
				//Jmol._$(canvas.id).css({background:"rgb(0,0,0,0.001)", "z-index":Jmol._z.main}); 
			} else {
				this._mouseInterface = canvas;
			}
			Jmol._jsSetMouse(this._mouseInterface);
		}
    
    proto._getLayer = function(name, container, w, h, isOpaque) {
  		var c = document.createElement("canvas");
			this._canvas[name + "Layer"] = c;
			c.style.width = "100%";
			c.style.height = "100%";
			c.id = this._id + "_" + name + "Layer";
			c.width = w;
			c.height = h; // w and h used in setScreenDimension
			container.append(c);
			c.applet = this;
			Jmol._$(c.id).css({background:(isOpaque ? "rgb(0,0,0,1)" : "rgb(0,0,0,0.001)"), "z-index": Jmol._getZ(this,name),position:"absolute",left:"0px",top:"0px",overflow:"hidden"});
			return c;	
    }
    
		proto._setupJS = function() {
			window["j2s.lib"] = {
				base : this._j2sPath + "/",
				alias : ".",
				console : this._console,
				monitorZIndex : Jmol._getZ(this, "monitorZIndex")
			};
			var isFirst = (Jmol.__execStack.length == 0);
			if (isFirst)
				Jmol.__addExec([this, Jmol.__loadClazz, null, "loadClazz"]);
			if (this._isJSV) {
				Jmol.__addCoreFile("jsv", this._j2sPath, this.__Info.preloadCore);
				if (Jmol._debugCode) {
				// no min package for that
					Jmol.__addExec([this, Jmol.__loadClass, "JSV.appletjs.JSVApplet", "load JSV"]);
					if (this._isPro)
						Jmol.__addExec([this, Jmol.__loadClass, "JSV.appletjs.JSVAppletPro", "load JSV(signed)"]);
				}
			} else {
				Jmol.__addCoreFile("jmol", this._j2sPath, this.__Info.preloadCore);
				if (!this._is2D) {
			 		Jmol.__addExec([this, Jmol.__loadClass, "J.export.JSExporter","load JSExporter"])
			//		Jmol.__addExec([this, this.__addExportHook, null, "addExportHook"])
				}			 			
				if (Jmol._debugCode)
					Jmol.__addExec([this, Jmol.__loadClass, "J.appletjs.Jmol", "load Jmol"]);
			}
			Jmol.__addExec([this, this.__startAppletJS, null, "start applet"])

			this._isSigned = true; // access all files via URL hook
			this._ready = false; 
			this._applet = null;
			this._canScript = function(script) {return true;};
			this._savedOrientations = [];
			Jmol.__execTimer && clearTimeout(Jmol.__execTimer);
			Jmol.__execTimer = setTimeout(Jmol.__nextExecution, Jmol.__execDelayMS);
		};

//		proto.__addExportHook = function(applet) {
	//		Jmol.GLmol.addExportHook(applet);
		//	Jmol.__nextExecution();
		//};

		proto.__startAppletJS = function(applet) {
			var viewerOptions =  new java.util.Hashtable ();
			Jmol._setJmolParams(viewerOptions, applet.__Info, true);
			viewerOptions.put("appletReadyCallback","Jmol._readyCallback");
			viewerOptions.put("applet", true);
			viewerOptions.put("name", applet._id);// + "_object");
			viewerOptions.put("syncId", Jmol._syncId);      
			if (applet._color) 
				viewerOptions.put("bgcolor", applet._color);
			if (!applet._is2D)  
				viewerOptions.put("script", "set multipleBondSpacing 0.35;");
			else if (applet._startupScript)
				viewerOptions.put("script", applet._startupScript)

			if (Jmol._syncedApplets.length) {
				viewerOptions.put("synccallback", "Jmol._mySyncCallback");
			}
			viewerOptions.put("signedApplet", "true");
			viewerOptions.put("platform", applet._platform);
			if (applet._is2D)
				viewerOptions.put("display",applet._id + "_canvas2d");

			// viewerOptions.put("repaintManager", "J.render");
			viewerOptions.put("documentBase", document.location.href);
			var codePath = applet._j2sPath + "/";
			if (codePath.indexOf("://") < 0) {
				var base = document.location.href.split("#")[0].split("?")[0].split("/");
				if (codePath.indexOf("/") == 0)
					base = [base[0], codePath.substring(1)];
				else
					base[base.length - 1] = codePath;
				codePath = base.join("/");
			}
			viewerOptions.put ("codePath", codePath);

			Jmol._registerApplet(applet._id, applet);
			applet._applet = (!applet._isJSV ? new J.appletjs.Jmol(viewerOptions) 
				: applet._isPro ? new JSV.appletjs.JSVAppletPro(viewerOptions) 
				: new JSV.appletjs.JSVApplet(viewerOptions));

			if (!applet._is2D)
				applet._GLmol.applet = applet;
			applet._jsSetScreenDimensions();      
			Jmol.__nextExecution();
		};

		proto._jsSetScreenDimensions = function() {
				if (!this._applet)return
				// strangely, if CTRL+/CTRL- are used repeatedly, then the
				// applet div can be not the same size as the canvas if there
				// is a border in place.
				var d = Jmol._getElement(this, (this._is2D ? "canvas2d" : "canvas"));
				this._applet.viewer.setScreenDimension(
				d.width, d.height);
// Math.floor(Jmol.$(this, "appletdiv").height()));

		};

		proto._show = function(tf) {
			Jmol.$setVisible(Jmol.$(this,"appletdiv"), tf);
			if (tf)
				Jmol._repaint(this, true);
		};

		proto._canScript = function(script) {return true};

		proto._delay = function(eval, sc, millis) {
		// does not take into account that scripts may be added after this and
		// need to be cached.
			this._delayID = setTimeout(function(){eval.resumeEval(sc,false)}, millis);		
		}

		proto._createDomNode = function(id, data) {
			id = this._id + "_" + id;
			var d = document.getElementById(id);
			if (d)
				document.body.removeChild(d);
			if (!data)
				return;
			if (data.indexOf("<?") == 0)
				data = data.substring(data.indexOf("<", 1));
			if (data.indexOf("/>") >= 0) {
				// no doubt there is a more efficient way to do this.
				// Firefox, at least, does not recognize "/>" in HTML blocks
				// that are added this way.
				var D = data.split("/>");
				for (var i = D.length - 1; --i >= 0;) {
					var s = D[i];
					var pt = s.lastIndexOf("<") + 1;
					var pt2 = pt;
					var len = s.length;
					var name = "";
					while (++pt2 < len) {
						if (" \t\n\r".indexOf(s.charAt(pt2))>= 0) {
							var name = s.substring(pt, pt2);
							D[i] = s + "></"+name+">";
							break;
						}	  	
					}
				}
				data = D.join('');
			}
			d = document.createElement("_xml")
			d.id = id;
			d.innerHTML = data;
			d.style.display = "none";
			document.body.appendChild(d);
			return d;
		}		

		proto.equals = function(a) { return this == a };
		proto.clone = function() { return this };
		proto.hashCode = function() { return parseInt(this._uniqueId) };  


		proto._processGesture = function(touches) {
			return this._applet.viewer.mouse.processTwoPointGesture(touches);
		}

		proto._processEvent = function(type, xym) {
			this._applet.viewer.processMouseEvent(type,xym[0],xym[1],xym[2],System.currentTimeMillis());
		}

		proto._resize = function() {
			var s = "__resizeTimeout_" + this._id;
			// only at end
			if (Jmol[s])
				clearTimeout(Jmol[s]);
			var me = this;
			Jmol[s] = setTimeout(function() {Jmol._repaint(me, true);Jmol[s]=null}, 100);
		}

		return proto;
	};

	Jmol._repaint = function(applet, asNewThread) {
		// asNewThread: true is from RepaintManager.repaintNow()
		// false is from Repaintmanager.requestRepaintAndWait()
		//

		// alert("_repaint " + arguments.callee.caller.caller.exName)
		if (!applet || !applet._applet)return;

		// asNewThread = false;
		var container = Jmol.$(applet, "appletdiv");
		var w = Math.round(container.width());
		var h = Math.round(container.height());
		if (applet._is2D && (applet._canvas.width != w || applet._canvas.height != h)) {
			applet._newCanvas(true);
			applet._applet.viewer.setDisplay(applet._canvas);
		}
		applet._applet.viewer.setScreenDimension(w, h);

		if (asNewThread) {
			setTimeout(function(){ applet._applet && applet._applet.viewer.updateJS(0,0)});
		} else {
			applet._applet.viewer.updateJS(0,0);
		}
		// System.out.println(applet._applet.fullName)
	}

	Jmol._getHiddenCanvas = function(applet, id, width, height, forceNew) {
		id = applet._id + "_" + id;
		var d = document.getElementById(id);
		if (d && forceNew) {
			d = null;
		}
		if (!d)
			d = document.createElement( 'canvas' );
			// for some reason both these need to be set, or maybe just d.width?
		d.width = d.style.width = width;
		d.height = d.style.height = height;
		// d.style.display = "none";
		if (d.id != id) {
			d.id = id;
		}
		return d;
	}

	Jmol._loadImage = function(platform, echoNameAndPath, bytes, fOnload, image) {
	// bytes would be from a ZIP file -- will have to reflect those back from
	// server as an image after conversion to base64
	// ah, but that's a problem, because that image would be needed to be
	// posted, but you can't post an image call.
	// so we would have to go with <image data:> which does not work in all
	// browsers. Hmm.

		var path = echoNameAndPath[1];

		if (image == null) {
			var image = new Image();
			image.onload = function() {Jmol._loadImage(platform, echoNameAndPath, null, fOnload, image)};

			if (bytes != null) {      
				bytes = J.io.Base64.getBase64(bytes).toString();      
				var filename = path.substring(url.lastIndexOf("/") + 1);                                    
				var mimetype = (filename.indexOf(".png") >= 0 ? "image/png" : filename.indexOf(".jpg") >= 0 ? "image/jpg" : "");
				 // now what?
			}
			image.src = path;
			return true; // as far as we can tell!
		}
		var width = image.width;
		var height = image.height; 
		var id = "echo_" + echoNameAndPath[0];
		var canvas = Jmol._getHiddenCanvas(platform.vwr.applet, id, width, height, true);
		canvas.imageWidth = width;
		canvas.imageHeight = height;
		canvas.id = id;
		canvas.image=image;
		Jmol._setCanvasImage(canvas, width, height);
		// return a null canvas and the error in path if there is a problem
		fOnload(canvas,path);
	};

	Jmol._setCanvasImage = function(canvas, width, height) {
		canvas.buf32 = null;
		canvas.width = width;
		canvas.height = height;
		canvas.getContext("2d").drawImage(canvas.image, 0, 0, width, height);
	};

})(Jmol);
// JmolApplet.js -- Jmol._Applet and Jmol._Image

// BH 5/8/2014 11:20:21 AM trying to fix AH nd JG problem with multiple applets
// BH 1/27/2014 8:36:43 AM adding Info.viewSet
// BH 12/13/2013 9:04:53 AM _evaluate DEPRECATED (see JSmolApi.js Jmol.evaulateVar
// BH 11/24/2013 11:41:31 AM streamlined createApplet, with added JNLP for local reading
// BH 10/11/2013 7:17:10 AM streamlined and made consistent with JSV and JSME
// BH 7/16/2012 1:50:03 PM adds server-side scripting for image
// BH 8/11/2012 11:00:01 AM adds Jmol._readyCallback for MSIE not in Quirks mode
// BH 8/12/2012 3:56:40 AM allows .min.png to be replaced by .all.png in Image file name
// BH 8/13/2012 6:16:55 PM fix for no-java message not displaying
// BH 11/18/2012 1:06:39 PM adds option ">" in database query box for quick command execution
// BH 12/17/2012 6:25:00 AM change ">" option to "!"

;(function (Jmol, document) {

	// _Applet -- the main, full-featured, Jmol object

	Jmol._Applet = function(id, Info, checkOnly){
		window[id] = this;
		this._jmolType = "Jmol._Applet" + (Info.isSigned ? " (signed)" : "");
		this._viewType = "Jmol";
		this._isJava = true;
		this._syncKeyword = "Select:";
		if (checkOnly)
			return this;
		this._isSigned = Info.isSigned;
		this._readyFunction = Info.readyFunction;
		this._ready = false;
		this._isJava = true; 
		this._isInfoVisible = false;
		this._applet = null;
		this._memoryLimit = Info.memoryLimit || 512;
		this._canScript = function(script) {return true;};
		this._savedOrientations = [];
		this._initialize = function(jarPath, jarFile) {
			var doReport = false;
			Jmol._jarFile && (jarFile = Jmol._jarFile);
			if(this._jarFile) {
				var f = this._jarFile;
				if(f.indexOf("/") >= 0) {
					alert ("This web page URL is requesting that the applet used be " + f + ". This is a possible security risk, particularly if the applet is signed, because signed applets can read and write files on your local machine or network.");
					var ok = prompt("Do you want to use applet " + f + "? ", "yes or no")
					if(ok == "yes") {
						jarPath = f.substring(0, f.lastIndexOf("/"));
						jarFile = f.substring(f.lastIndexOf("/") + 1);
					} else {
						doReport = true;
					}
				} else {
					jarFile = f;
				}
				this_isSigned = Info.isSigned = (jarFile.indexOf("Signed") >= 0);
			}
 			this._jarPath = Info.jarPath = jarPath || ".";
			this._jarFile = Info.jarFile = (typeof(jarFile) == "string" ? jarFile : (jarFile ?  "JmolAppletSigned" : "JmolApplet") + "0.jar");
			if (doReport)
				alert ("The web page URL was ignored. Continuing using " + this._jarFile + ' in directory "' + this._jarPath + '"');
			Jmol.controls == undefined || Jmol.controls._onloadResetForms();		
		}		
		this._create(id, Info);
		return this;
	}

	;(function(Applet, proto) {
	Applet._get = function(id, Info, checkOnly) {

	// note that the variable name the return is assigned to MUST match the first parameter in quotes
	// applet = Jmol.getApplet("applet", Info)

		checkOnly || (checkOnly = false);
		Info || (Info = {});
		var DefaultInfo = {
			color: "#FFFFFF", // applet object background color, as for older jmolSetBackgroundColor(s)
			width: 300,
			height: 300,
			addSelectionOptions: false,
			serverURL: "http://your.server.here/jsmol.php",
			defaultModel: "",
			script: null,
			src: null,
			readyFunction: null,
			use: "HTML5",//other options include JAVA, WEBGL, and IMAGE
			jarPath: "java",
			jarFile: "JmolApplet0.jar",
			isSigned: false,
			j2sPath: "j2s",
			coverImage: null,     // URL for image to display
			coverTitle: "",       // tip that is displayed before model starts to load
			coverCommand: "",     // Jmol command executed upon clicking image
			deferApplet: false,   // true == the model should not be loaded until the image is clicked
			deferUncover: false,  // true == the image should remain until command execution is complete 
			disableJ2SLoadMonitor: false,
			disableInitialConsole: false,
			debug: false
		};	 
		Jmol._addDefaultInfo(Info, DefaultInfo);
		Jmol._debugAlert = Info.debug;
			if (!Jmol.featureDetection.allowHTML5)Info.use = "JAVA";

		Info.serverURL && (Jmol._serverUrl = Info.serverURL);

		var javaAllowed = false;
		var applet = null;
		var List = Info.use.toUpperCase().split("#")[0].split(" ");
		for (var i = 0; i < List.length; i++) {
			switch (List[i]) {
			case "JAVA":
				javaAllowed = true;
				if (Jmol.featureDetection.supportsJava())
					applet = new Applet(id, Info, checkOnly);
				break;
			case "WEBGL":
				applet = Applet._getCanvas(id, Info, checkOnly, true);
				break;
			case "HTML5":
				applet = Applet._getCanvas(id, Info, checkOnly, false);
				break;
			case "IMAGE":
				applet = new Jmol._Image(id, Info, checkOnly);
				break;
			}
			if (applet != null)
				break;		  
		}
		if (applet == null) {
			if (checkOnly || !javaAllowed)
				applet = {_jmolType : "none" };
			else if (javaAllowed)
 		  	applet = new Applet(id, Info);
		}

		// keyed to both its string id and itself
		return (checkOnly ? applet : Jmol._registerApplet(id, applet));  
	}

	Applet._getCanvas = function(id, Info, checkOnly, webGL) {
		if (webGL && Jmol.featureDetection.supportsWebGL()) {
			Jmol._Canvas3D.prototype = Jmol.GLmol.extendApplet(Jmol._jsSetPrototype(new Applet(id, Info, true)));
			return new Jmol._Canvas3D(id, Info, "Jmol", checkOnly);
		}
		if (!webGL) {
			Jmol._Canvas2D.prototype = Jmol._jsSetPrototype(new Applet(id, Info, true));
			return new Jmol._Canvas2D(id, Info, "Jmol", checkOnly);
		}
		return null;
	};


	/*  AngelH, mar2007:
		By (re)setting these variables in the webpage before calling Jmol.getApplet(),
		a custom message can be provided (e.g. localized for user's language) when no Java is installed.
	*/
	Applet._noJavaMsg =
			"Either you do not have Java applets enabled in your web<br />browser or your browser is blocking this applet.<br />\
			Check the warning message from your browser and/or enable Java applets in<br />\
			your web browser preferences, or install the Java Runtime Environment from <a href='http://www.java.com'>www.java.com</a>";

	Applet._setCommonMethods = function(p) {
		p._showInfo = proto._showInfo;	
		p._search = proto._search;
		p._getName = proto._getName;
		p._readyCallback = proto._readyCallback;
	}

	Applet._createApplet = function(applet, Info, params) {
		applet._initialize(Info.jarPath, Info.jarFile);
		var jarFile = applet._jarFile;
		var jnlp = ""
		if (Jmol._isFile) {
			// local installations need jnlp here and should reference JmolApplet(Signed).jar, not JmolApplet(Signed)0.jar  
			jarFile = jarFile.replace(/0\.jar/,".jar");
			jnlp = " jnlp_href=\"" + jarFile.replace(/\.jar/,".jnlp") + "\"";
		}
		// size is set to 100% of containers' size, but only if resizable. 
		// Note that resizability in MSIE requires: 
		// <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
		var w = (applet._containerWidth.indexOf("px") >= 0 ? applet._containerWidth : "100%");
		var h = (applet._containerHeight.indexOf("px") >= 0 ? applet._containerHeight : "100%");
		var widthAndHeight = " style=\"width:" + w + ";height:" + h + "\" ";
		var attributes = "name='" + applet._id + "_object' id='" + applet._id + "_object' " + "\n"
				+ widthAndHeight + jnlp + "\n"
		params.codebase = applet._jarPath;
		params.codePath = params.codebase + "/";
		if (params.codePath.indexOf("://") < 0) {
			var base = document.location.href.split("#")[0].split("?")[0].split("/");
			base[base.length - 1] = params.codePath;
			params.codePath = base.join("/");
		}
		params.archive = jarFile;
		params.mayscript = 'true';
		params.java_arguments = "-Xmx" + Math.round(Info.memoryLimit || applet._memoryLimit) + "m";
		params.permissions = (applet._isSigned ? "all-permissions" : "sandbox");
		params.documentLocation = document.location.href;
		params.documentBase = document.location.href.split("#")[0].split("?")[0];

		params.jarPath = Info.jarPath;
		Jmol._syncedApplets.length && (params.synccallback = "Jmol._mySyncCallback");
		applet._startupScript && (params.script = applet._startupScript);
		var t = "\n"; 
 		for (var i in params)
			if(params[i])
		 		t += "  <param name='"+i+"' value='"+params[i]+"' />\n";
		if (Jmol.featureDetection.useIEObject || Jmol.featureDetection.useHtml4Object) {
			t = "<object " + attributes
				+ (Jmol.featureDetection.useIEObject ? 
					 " classid='clsid:8AD9C840-044E-11D1-B3E9-00805F499D93' codebase='http://java.sun.com/update/1.6.0/jinstall-6u22-windows-i586.cab'>"
				 : " type='application/x-java-applet'>")
				 + t + "<p style='background-color:yellow;" + widthAndHeight.split('"')[1] 
				+ ";text-align:center;vertical-align:middle;'>\n" + Applet._noJavaMsg + "</p></object>\n";
		} else { // use applet tag
			t = "<applet " + attributes
				+ " code='" + params.code + "' codebase='" + applet._jarPath + "' archive='" + jarFile + "' mayscript='true'>\n"
				+ t + "<table bgcolor='yellow'><tr><td align='center' valign='middle' " + widthAndHeight + ">\n"
				+ Applet._noJavaMsg + "</td></tr></table></applet>\n";
		}	
		t = Jmol._getWrapper(applet, true) + t + Jmol._getWrapper(applet, false) 
			+ (Info.addSelectionOptions ? Jmol._getGrabberOptions(applet) : "");
		if (Jmol._debugAlert)
			alert (t);
		applet._code = Jmol._documentWrite(t);
	}

	proto._create = function(id, Info){
		Jmol._setObject(this, id, Info);
		var params = {
			syncId: Jmol._syncId,
			progressbar: "true",                      
			progresscolor: "blue",
			boxbgcolor: this._color || "black",
			boxfgcolor: "white",
			boxmessage: "Downloading JmolApplet ...",
			script: (this._color ? "background \"" + this._color +"\"": ""),
			code: "JmolApplet.class"
		};

		Jmol._setJmolParams(params, Info);
		function sterilizeInline(model) {
			model = model.replace(/\r|\n|\r\n/g, (model.indexOf("|") >= 0 ? "\\/n" : "|")).replace(/'/g, "&#39;");
			if(Jmol._debugAlert)
				alert ("inline model:\n" + model);
			return model;
		}

		params.loadInline = (Info.inlineModel ? sterilizeInline(Info.inlineModel) : "");
		params.appletReadyCallback = "Jmol._readyCallback";
		if (Jmol._syncedApplets.length)
			params.synccallback = "Jmol._mySyncCallback";
		params.java_arguments = "-Xmx" + Math.round(Info.memoryLimit || this._memoryLimit) + "m";

		this._initialize(Info.jarPath, Info.jarFile);
		Applet._createApplet(this, Info, params);
	}

	proto._readyCallback = function(id, fullid, isReady, applet) {
		if (!isReady)
			return; // ignore -- page is closing
		Jmol._setDestroy(this);
		this._ready = true;
		var script = this._readyScript;
		this._applet = applet;
		if (this._defaultModel)
			Jmol._search(this, this._defaultModel, (script ? ";" + script : ""));
		else if (script)
			this._script(script);
		else if (this._src)
			this._script('load "' + this._src + '"');
		this._showInfo(true);
		this._showInfo(false);
		Jmol.Cache.setDragDrop(this);
		this._readyFunction && this._readyFunction(this);
		Jmol._setReady(this);
		var app = this._2dapplet;
		if (app && app._isEmbedded && app._ready && app.__Info.visible)
			this._show2d(true);
	}

	proto._showInfo = function(tf) {
		if(tf && this._2dapplet)
			this._2dapplet._show(false);
		Jmol.$html(Jmol.$(this, "infoheaderspan"), this._infoHeader);
		if (this._info)
			Jmol.$html(Jmol.$(this, "infodiv"), this._info);
		if ((!this._isInfoVisible) == (!tf))
			return;
		this._isInfoVisible = tf;
		// 1px does not work for MSIE
		if (this._isJava) {
			var x = (tf ? 2 : "100%");
			Jmol.$setSize(Jmol.$(this, "appletdiv"), x, x);
		}
		Jmol.$setVisible(Jmol.$(this, "infotablediv"), tf);
		Jmol.$setVisible(Jmol.$(this, "infoheaderdiv"), tf);
		this._show(!tf);
	}

	proto._show2d = function(tf) {
		this._2dapplet._show2d(tf);
		if (this._2dapplet._isEmbedded) {
			this._showInfo(false);
			this._show(!tf);
		}
	}

  proto._getAtomCorrelation = function(molData) {
    // get the first atom mapping available by loading the model structure into model 2, 
    this._loadMolData(molData, "atommap = compare({1.1} {2.1} 'MAP' 'H'); zap 2.1", true);
    var map = jmol._evaluate("atommap");
    var n = jmol._evaluate("{*}.count");
    var A = [];
    var B = [];
    // these are Jmol atom indexes. The second number will be >= n, and all must be incremented by 1.
		for (var i = 0; i < map.length; i++) {
		  var c = map[i];
		  A[c[0] + 1] = c[1] - n + 1;
		  B[c[1] - n + 1] = c[0] + 1;
		}
		return {fromJmol:A, toJmol:B}; // forward and rev.		
  }
  
	proto._show = function(tf) {
		var x = (!tf ? 2 : "100%");
		Jmol.$setSize(Jmol.$(this, "object"), x, x);
		if (!this._isJava)
			Jmol.$setVisible(Jmol.$(this, "appletdiv"), tf);
	}

	proto._clearConsole = function () {
			if (this._console == this._id + "_infodiv")
				this.info = "";
			if (!self.Clazz)return;
			Jmol._setConsoleDiv(this._console);
			Clazz.Console.clear();
		}


	proto._addScript = function(script) {      
		this._readyScript || (this.readyScript = "");
		this._readyScript && (this._readyScript += ";");
		this._readyScript += script;
		return true;
	}

	proto._script = function(script) {
		if (!this._ready)
				return this._addScript(script);
		Jmol._setConsoleDiv(this._console);
		this._applet.script(script);
	}

	proto._syncScript = function(script) {
		this._applet.syncScript(script);
	}

	proto._scriptCheck = function(script) {
		return this._ready && this._applet.scriptCheck(script);  
	}

	proto._scriptWait = function(script, noReturn) {
		var Ret = this._scriptWaitAsArray(script);
		var s = "";
		if (!noReturn)
		for(var i = Ret.length; --i >= 0; )
			for(var j = 0, jj = Ret[i].length; j < jj; j++)
				s += Ret[i][j] + "\n";
		return s;
	}

	proto._scriptEcho = function(script) {
		// returns a newline-separated list of all echos from a script
		var Ret = this._scriptWaitAsArray(script);
		var s = "";
		for(var i = Ret.length; --i >= 0; )
			for(var j = Ret[i].length; --j >= 0; )
				if(Ret[i][j][1] == "scriptEcho")
					s += Ret[i][j][3] + "\n";
		return s.replace(/ \| /g, "\n");
	}

	proto._scriptMessage = function(script) {
		// returns a newline-separated list of all messages from a script, ending with "script completed\n"
		var Ret = this._scriptWaitAsArray(script);
		var s = "";
		for(var i = Ret.length; --i >= 0; )
			for(var j = Ret[i].length; --j >= 0; )
				if(Ret[i][j][1] == "scriptStatus")
					s += Ret[i][j][3] + "\n";
		return s.replace(/ \| /g, "\n");
	}

	proto._scriptWaitOutput = function(script) {
		var ret = "";
		try {
			if(script) {
				ret += this._applet.scriptWaitOutput(script);
			}
		} catch(e) {
		}
		return ret;
	}

	proto._scriptWaitAsArray = function(script) {
		var ret = "";
		try {
			this._getStatus("scriptEcho,scriptMessage,scriptStatus,scriptError");
			if(script) {
				ret += this._applet.scriptWait(script);
				ret = Jmol._evalJSON(ret, "jmolStatus");
				if( typeof ret == "object")
					return ret;
			}
		} catch(e) {
		}
		return [[ret]];
	}

	proto._getStatus = function(strStatus) {
		return Jmol._sortMessages(this._getPropertyAsArray("jmolStatus",strStatus));
	}

	proto._getPropertyAsArray = function(sKey,sValue) {
		return Jmol._evalJSON(this._getPropertyAsJSON(sKey,sValue),sKey);
	}

	proto._getPropertyAsString = function(sKey,sValue) {
		sValue == undefined && ( sValue = "");
		return this._applet.getPropertyAsString(sKey, sValue) + "";
	}

	proto._getPropertyAsJSON = function(sKey,sValue) {
		sValue == undefined && ( sValue = "");
		try {
			return (this._applet.getPropertyAsJSON(sKey, sValue) + "");
		} catch(e) {
			return "";
		}
	}

	proto._getPropertyAsJavaObject = function(sKey,sValue) {		
		sValue == undefined && ( sValue = "");
		return this._applet.getProperty(sKey,sValue);
	}

  proto._evaluate = function(expr) {  
  	expr != null || (expr = "");
		return this._getPropertyAsArray("variableInfo", expr);

  }

	proto._evaluateDEPRECATED = function(molecularMath) {   // DEPRECATED!!!	
	// DEPRECATED!!!	
		//carries out molecular math on a model
		var result = "" + this._getPropertyAsJavaObject("evaluate", molecularMath);
		var s = result.replace(/\-*\d+/, "");
		if(s == "" && !isNaN(parseInt(result)))
			return parseInt(result);
		var s = result.replace(/\-*\d*\.\d*/, "")
		if(s == "" && !isNaN(parseFloat(result)))
			return parseFloat(result);
		return result;
	// DEPRECATED!!!	
	}

	proto._saveOrientation = function(id) {	
		return this._savedOrientations[id] = this._getPropertyAsArray("orientationInfo","info").moveTo;
	}

	proto._restoreOrientation = function(id) {
		var s = this._savedOrientations[id];
		if(!s || s == "")
			return s = s.replace(/1\.0/, "0");
		return this._scriptWait(s);
	}


	proto._restoreOrientationDelayed = function(id,delay) {
		arguments.length < 1 && ( delay = 1);
		var s = this._savedOrientations[id];
		if(!s || s == "")
			return s = s.replace(/1\.0/, delay);
		return this._scriptWait(s);
	}

	proto._resizeApplet = function(size) {
		// See _jmolGetAppletSize() for the formats accepted as size [same used by jmolApplet()]
		//  Special case: an empty value for width or height is accepted, meaning no change in that dimension.

		/*
		 * private functions
		 */
		function _getAppletSize(size, units) {
			/* Accepts single number, 2-value array, or object with width and height as mroperties, each one can be one of:
			 percent (text string ending %), decimal 0 to 1 (percent/100), number, or text string (interpreted as nr.)
			 [width, height] array of strings is returned, with units added if specified.
			 Percent is relative to container div or element (which should have explicitly set size).
			 */
			var width, height;
			if(( typeof size) == "object" && size != null) {
				width = size[0]||size.width;
				height = size[1]||size.height;
			} else {
				width = height = size;
			}
			return [_fixDim(width, units), _fixDim(height, units)];
		}

		function _fixDim(x, units) {
			var sx = "" + x;
			return (sx.length == 0 ? (units ? "" : Jmol._allowedJmolSize[2]) 
				: sx.indexOf("%") == sx.length - 1 ? sx 
				: (x = parseFloat(x)) <= 1 && x > 0 ? x * 100 + "%" 
				: (isNaN(x = Math.floor(x)) ? Jmol._allowedJmolSize[2] 
				: x < Jmol._allowedJmolSize[0] ? Jmol._allowedJmolSize[0] 
				: x > Jmol._allowedJmolSize[1] ? Jmol._allowedJmolSize[1] 
				: x)
				+ (units ? units : "")
			);
		}

		var sz = _getAppletSize(size, "px");
		var d = Jmol._getElement(this, "appletinfotablediv");
		d.style.width = sz[0];
		d.style.height = sz[1];
		this._containerWidth = sz[0];
		this._containerHeight = sz[1];
		if (this._is2D)
			Jmol._repaint(this, true);
	}

	proto._search = function(query, script){
		Jmol._search(this, query, script);
	}

	proto._searchDatabase = function(query, database, script, _jmol_searchDatabase){
		if (this._2dapplet && this._2dapplet._isEmbedded && Jmol.$(this, "2dappletdiv:visible")[0])
			return this._2dapplet._searchDatabase(query, database, script); 
		this._showInfo(false);
		if (query.indexOf("?") >= 0) {
			Jmol._getInfoFromDatabase(this, database, query.split("?")[0]);
			return;
		}
		script || (script = Jmol._getScriptForDatabase(database));
		var dm = database + query;
		this._currentView = null;
		this._searchQuery = dm;
 		this._loadFile(dm, script, dm);
	}

	proto._loadFile = function(fileName, script, chemID, _jmol_loadFile){
		this._showInfo(false);
		script || (script = "");
		this._thisJmolModel = "" + Math.random();
		this._fileName = fileName;
		if (!this._scriptLoad(fileName, script)) {
			// we load the data here instead of in Jmol in the case of
			// JSmol/Java/Sandboxed or when part of a view set 
			var me = this;      
			Jmol._loadFileData(this, fileName, 
				function(data){me.__loadModel(data, script, chemID)},
				function() {me.__loadModel(null)}
			);
		}
	}
					 
	proto._scriptLoad = function(file, script, _jmol_scriptLoad) {
		script || (script = "");
		var doscript = (this._isJava || !this._noscript);
		if (doscript)
			this._script("zap;set echo middle center;echo Retrieving data...");
		if (!this._isSigned || this._viewSet != null)
			return false;
		if (doscript)
			this._script("load async \"" + file + "\";" + script);
		else
			this._applet.viewer.openFile(file);
		this._checkDeferred("");
		return true;
	}

	proto.__loadModel = function(mol, script, chemID, _jmol__loadModel) {
		if (mol == null)
			return;
		if (this._viewSet != null) {
			script || (script = ""); 
			// first component only
			script += ";if ({*}.molecule.max > 1 || {*}.modelindex.max > 0){ delete molecule > 1 or modelindex > 0;x = getProperty('extractModel',{*});load inline @x};"
		}
		if (!script && this._noscript) {
			this._applet.viewer.loadInline(mol, '\0');
		} else {
		  this._loadMolData(mol, script, false);
		}
		if (this._viewSet != null) {
			Jmol.View.updateView(this, {chemID:chemID, data:mol});
		}      
	}

  proto._loadMolData = function(mol, script, isAppend) {
		script || (script = "");
		var name = (isAppend ? "append" : "model");
		script = 'load DATA "' + name + '"' + mol + '\nEND "'+ name +'" ;' + script;
		this._applet.scriptWait(script);
  }
  
	proto._loadModelFromView = function(view, _jmol_loadModelFromView) {
		// request from Jmol.View to update view with view.JME.data==null or needs changing
		this._currentView = view;
		var rec = view.Jmol;
		if (rec.data != null) {
			this.__loadModel(rec.data, null, view.info.chemID);
			return;
		}
		if (view.info.chemID != null) {
			Jmol._searchMol(this, view.info.chemID, null, false);
			return;
		}
		rec = view.JME;
		if (rec) {
			rec.applet._show2d(false, this);
			return;
		}
	}

	proto._updateView = function(_jmol_updateView) {
		if (this._viewSet == null || !this._applet)
			return;
		// called from model change without chemical identifier, possibly by user action and call to Jmol.updateView(applet)
		chemID = "" + this._getPropertyAsJavaObject("variableInfo","script('show chemical inchiKey')");
		if (chemID.length() < 36) // InChIKey=RZVAJINKPMORJF-BGGKNDAXNA-N
			chemID = null;
		else
			chemID = chemID.substring(36).split('\n')[0];
		Jmol.View.updateView(this, {chemID:chemID, data: "" + this._getPropertyAsJavaObject("evaluate", "extractModel", "{visible}")});
	}

	proto._atomPickedCallback = function(imodel, iatom, _jmol_atomPickedCallback) {
		// direct callback from Jmol HTML5 applet
		if (iatom < 0) {
		// TODO could be a model change? 
		} else {
			var A = [iatom + 1];
			Jmol.View.updateAtomPick(this, A);
			this._updateAtomPick(A);
		}
	}

	proto._updateAtomPick = function(A) {
			this._script(A.length == 0 ? "select none" : "select on visible and (@" + A.join(",@") + ")");
	}

  proto._isDeferred = function () {
      return !this._canvas && this._cover && this._isCovered && this._deferApplet
  }

	proto._checkDeferred = function(script) {
		if (this._isDeferred()) {
			this._coverScript = script;
			this._cover(false);
			return true;
		}
		return false;
	}

  proto._getSmiles = function() {
		return this._evaluate("{visible}.find('SMILES')");   
  }
  
  proto._getMol = function() {
		return this._evaluate("getProperty('ExtractModel',{visible})");   
  }

  proto._getMol2D = function() {
		return jmol._evaluate("script('select visible;show chemical sdf')"); // 2D equivalent
  }
  
  
})(Jmol._Applet, Jmol._Applet.prototype);

/* ****************************************


	// _Image -- an alternative to _Applet
	// commented out here, as it has found no use

	Jmol._Image = function(id, Info, checkOnly){
		this._jmolType = "image";
		if (checkOnly)
			return this;
		this._create(id, Info);
		return this;
	}

;(function (Image, proto) {

	Jmol._Applet._setCommonMethods(proto);

	proto._create = function(id, Info) {
		Jmol._setObject(this, id, Info);
		this._src || (this._src = "");
		var t = Jmol._getWrapper(this, true) 
			+ '<img id="'+id+'_image" width="' + Info.width + '" height="' + Info.height + '" src=""/>'
		 	+	Jmol._getWrapper(this, false)
			+ (Info.addSelectionOptions ? Jmol._getGrabberOptions(this) : "");
		if (Jmol._debugAlert)
			alert (t);
		this._code = Jmol._documentWrite(t);
		this._ready = false;
		if (Jmol._document)
			this._readyCallback(id, null, this._ready = true, null);
	}

	proto._canScript = function(script) {
		var slc = script.toLowerCase().replace(/[\",\']/g, '');
		var ipt = slc.length;
		return (script.indexOf("#alt:LOAD") >= 0 || slc.indexOf(";") < 0 && slc.indexOf("\n") < 0
			&& (slc.indexOf("script ") == 0 || slc.indexOf("load ") == 0)
			&& (slc.indexOf(".png") == ipt - 4 || slc.indexOf(".jpg") == ipt - 4));
	}

	proto._script = function(script) {
		var slc = script.toLowerCase().replace(/[\",\']/g, '');
		// single command only
		// "script ..." or "load ..." only
		// PNG or PNGJ or JPG only
		// automatically switches to .all.png(j) from .min.png(j)
		var ipt = slc.length;
		if (slc.indexOf(";") < 0 && slc.indexOf("\n") < 0
			&& (slc.indexOf("script ") == 0 || slc.indexOf("load ") == 0)
			&& (slc.indexOf(".png") == ipt - 4 || slc.indexOf(".pngj") == ipt - 5 || slc.indexOf(".jpg") == ipt - 4)) {
			var imageFile = script.substring(script.indexOf(" ") + 1);
			ipt = imageFile.length;
			for (var i = 0; i < ipt; i++) {
				switch (imageFile.charAt(i)) {
				case " ":
					continue;
				case '"':
					imageFile = imageFile.substring(i + 1, imageFile.indexOf('"', i + 1))
					i = ipt;
					continue;
				case "'":
					imageFile = imageFile.substring(i + 1, imageFile.indexOf("'", i + 1))
					i = ipt;
					continue;
				default:
					imageFile = imageFile.substring(i)
					i = ipt;
					continue;
				}
			}
			imageFile = imageFile.replace(/\.min\.png/,".all.png")
			document.getElementById(this._id + "_image").src = imageFile
		} else if (script.indexOf("#alt:LOAD ") >= 0) {
			imageFile = script.split("#alt:LOAD ")[1]
			if (imageFile.indexOf("??") >= 0) {
				var db = imageFile.split("??")[0];
				imageFile = prompt(imageFile.split("??")[1], "");
				if (!imageFile)
					return;
				if (!Jmol.db._DirectDatabaseCalls[imageFile.substring(0,1)])
					imageFile = db + imageFile;
			}
			this._loadFile(imageFile);
		}
	}

	proto._show = function(tf) {
		Jmol._getElement(this, "appletdiv").style.display = (tf ? "block" : "none");
	}

	proto._loadFile = function(fileName, params){
		this._showInfo(false);
		this._thisJmolModel = "" + Math.random();
		params = (params ? params : "");
		var database = "";
		if (Jmol._isDatabaseCall(fileName)) {
			database = fileName.substring(0, 1); 
			fileName = Jmol._getDirectDatabaseCall(fileName, false);
		} else if (fileName.indexOf("://") < 0) {
			var ref = document.location.href
			var pt = ref.lastIndexOf("/");
			fileName = ref.substring(0, pt + 1) + fileName;
		}

		var src = Jmol._serverUrl 
				+ "?call=getImageForFileLoad"
				+ "&file=" + escape(fileName)
				+ "&width=" + this._width
				+ "&height=" + this._height
				+ "&params=" + encodeURIComponent(params + ";frank off;");
		Jmol._getElement(this, "image").src = src;
	}

	proto._searchDatabase = function(query, database, script){
		if (query.indexOf("?") == query.length - 1) {
			Jmol._getInfoFromDatabase(this, database, query.split("?")[0]);
			return;
		}
		this._showInfo(false);
		script || (script = Jmol._getScriptForDatabase(database));
		var src = Jmol._serverUrl 
			+ "?call=getImageFromDatabase"
			+ "&database=" + database
			+ "&query=" + query
			+ "&width=" + this._width
			+ "&height=" + this._height
			+ "&script=" + encodeURIComponent(script + ";frank off;");
		Jmol._getElement(this, "image").src = src;
	}
})(Jmol._Image, Jmol._Image.prototype);

************************************ */

	Jmol.jmolSmiles = function(jmol, withStereoChemistry) {
		return jmol._getSmiles();
	}


})(Jmol, document);
// JSmolControls.js
//
// BH 5/29/2014 8:14:06 AM added default command for command input box
// BH 5/15/2014 -- removed script check prior to execution
// BH 12/3/2013 12:39:48 PM added up/down arrow key-driven command history for commandInput (changed keypress to keydown)
// BH 5/16/2013 8:14:47 AM fix for checkbox groups and default radio names
// BH 8:36 AM 7/27/2012  adds name/id for cmd button 
// BH 8/12/2012 6:51:53 AM adds function() {...} option for all controls:
//    Jmol.jmolButton(jmol, function(jmol) {...}, "xxxx")

(function(Jmol) {

	// private

	var c = Jmol.controls = {

		_hasResetForms: false,	
		_scripts: [""],
		_checkboxMasters: {},
		_checkboxItems: {},
		_actions: {},

		_buttonCount: 0,
		_checkboxCount: 0,
		_radioGroupCount: 0,
		_radioCount: 0,
		_linkCount: 0,
		_cmdCount: 0,
		_menuCount: 0,

		_previousOnloadHandler: null,	
		_control: null,
		_element: null,

		_appletCssClass: null,
		_appletCssText: "",
		_buttonCssClass: null,
		_buttonCssText: "",
		_checkboxCssClass: null,
		_checkboxCssText: "",
		_radioCssClass: null,
		_radioCssText: "",
		_linkCssClass: null,
		_linkCssText: "",
		_menuCssClass: null,
		_menuCssText: ""
	};

	c._addScript = function(appId,script) {
		var index = c._scripts.length;
		c._scripts[index] = [appId, script];
		return index;
	}

	c._getIdForControl = function(appletOrId, script) {
	//alert(appletOrId + " " + typeof appletOrId + " " + script + appletOrId._canScript)
		return (typeof appletOrId == "string" ? appletOrId 
			: !script || !appletOrId._canScript || appletOrId._canScript(script) ? appletOrId._id
			: null);
	}

	c._radio = function(appletOrId, script, labelHtml, isChecked, separatorHtml, groupName, id, title) {
		var appId = c._getIdForControl(appletOrId, script);
		if (appId == null)
			return null;
		++c._radioCount;
		groupName != undefined && groupName != null || (groupName = "jmolRadioGroup" + (c._radioGroupCount - 1));
		if (!script)
			return "";
		id != undefined && id != null || (id = "jmolRadio" + (c._radioCount - 1));
		labelHtml != undefined && labelHtml != null || (labelHtml = script.substring(0, 32));
		separatorHtml || (separatorHtml = "");
		var eospan = "</span>";
		c._actions[id] = c._addScript(appId, script);
		var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input name='"
			+ groupName + "' id='"+id+"' type='radio'"
			+ " onclick='Jmol.controls._click(this);return true;'"
			+ " onmouseover='Jmol.controls._mouseOver(this);return true;'"
			+ " onmouseout='Jmol.controls._mouseOut()' " +
		 (isChecked ? "checked='true' " : "") + c._radioCssText + " />";
		if (labelHtml.toLowerCase().indexOf("<td>")>=0) {
			t += eospan;
			eospan = "";
		}
		t += "<label for=\"" + id + "\">" + labelHtml + "</label>" +eospan + separatorHtml;
		return t;
	}

/////////// events //////////

	c._scriptExecute = function(element, scriptInfo) {
		var applet = Jmol._applets[scriptInfo[0]];
		var script = scriptInfo[1];
		if (typeof(script) == "object")
			script[0](element, script, applet);
		else if (typeof(script) == "function")
			script(applet);
		else
			Jmol.script(applet, script);
	}

	c.__checkScript = function(applet, d) {
		var ok = (d.value.indexOf("JSCONSOLE ") >= 0 || applet._scriptCheck(d.value) === "");
		d.style.color = (ok ? "black" : "red");
		return ok;
	} 

	c.__getCmd = function(dir, d) {
		if (!d._cmds || !d._cmds.length)return
		var s = d._cmds[d._cmdpt = (d._cmdpt + d._cmds.length + dir) % d._cmds.length]
		setTimeout(function(){d.value = s},10);    
		d._cmdadd = 1;
		d._cmddir = dir;
	}

	c._commandKeyPress = function(e, id, appId) {
	var keycode = (e == 13 ? 13 : window.event ? window.event.keyCode : e ? e.keyCode || e.which : 0);
	var d = document.getElementById(id);
		var applet = Jmol._applets[appId];
	switch (keycode) {
	case 13:
		var v = d.value;
		if ((c._scriptExecute(d, [appId, v]) || 1)) {
			 if (!d._cmds){
				 d._cmds = [];
				 d._cmddir = 0;
				 d._cmdpt = -1;
				 d._cmdadd = 0;      
	}
			 if (v && d._cmdadd == 0) {
					++d._cmdpt;
					d._cmds.splice(d._cmdpt, 0, v);
					d._cmdadd = 0;
					d._cmddir = 0;
			 } else {
					//d._cmdpt -= d._cmddir;
					d._cmdadd = 0;
			 }
			 d.value = "";
		}
		return false;
	case 27:
		setTimeout(function() {d.value = ""}, 20);
		return false;
	case 38: // up
		c.__getCmd(-1, d);
		break;
	case 40: // dn
		c.__getCmd(1, d);
		break;
	default:
		d._cmdadd = 0;
	}
	setTimeout(function() {c.__checkScript(applet, d)}, 20);
	return true;
 }

	c._click = function(obj, scriptIndex) {
		c._element = obj;
		if (arguments.length == 1)
			scriptIndex = c._actions[obj.id];
		c._scriptExecute(obj, c._scripts[scriptIndex]);
	}

	c._menuSelected = function(menuObject) {
		var scriptIndex = menuObject.value;
		if (scriptIndex != undefined) {
			c._scriptExecute(menuObject, c._scripts[scriptIndex]);
			return;
		}
		var len = menuObject.length;
		if (typeof len == "number")
			for (var i = 0; i < len; ++i)
				if (menuObject[i].selected) {
					c._click(menuObject[i], menuObject[i].value);
					return;
				}
		alert("?Que? menu selected bug #8734");
	}

	c._cbNotifyMaster = function(m){
		//called when a group item is checked
		var allOn = true;
		var allOff = true;
		for (var chkBox in m.chkGroup){
			if(m.chkGroup[chkBox].checked)
				allOff = false;
			else
				allOn = false;
		}
		if (allOn)m.chkMaster.checked = true;
		if (allOff)m.chkMaster.checked = false;
		if ((allOn || allOff) && c._checkboxItems[m.chkMaster.id])
			c._cbNotifyMaster(c._checkboxItems[m.chkMaster.id])
	}

	c._cbNotifyGroup = function(m, isOn){
		//called when a master item is checked
		for (var chkBox in m.chkGroup){
			var item = m.chkGroup[chkBox]
			if (item.checked != isOn) {
				item.checked = isOn;
				c._cbClick(item);
			}
			if (c._checkboxMasters[item.id])
				c._cbNotifyGroup(c._checkboxMasters[item.id], isOn)
		}
	}

	c._cbSetCheckboxGroup = function(chkMaster, chkboxes, args){
		var id = chkMaster;
		if(typeof(id)=="number")id = "jmolCheckbox" + id;
		chkMaster = document.getElementById(id);
		if (!chkMaster)alert("jmolSetCheckboxGroup: master checkbox not found: " + id);
		var m = c._checkboxMasters[id] = {};
		m.chkMaster = chkMaster;
		m.chkGroup = {};
		var i0;
		if (typeof(chkboxes)=="string") {
			chkboxes = args;
			i0 = 1;
		} else {
			i0 = 0;
		}
		for (var i = i0; i < chkboxes.length; i++){
			var id = chkboxes[i];
			if(typeof(id)=="number")id = "jmolCheckbox" + id;
			checkboxItem = document.getElementById(id);
			if (!checkboxItem)alert("jmolSetCheckboxGroup: group checkbox not found: " + id);
			m.chkGroup[id] = checkboxItem;
			c._checkboxItems[id] = m;
		}
	}

	c._cbClick = function(ckbox) {
		c._control = ckbox;
		var whenChecked = c._actions[ckbox.id][0];
		var whenUnchecked = c._actions[ckbox.id][1];
		c._click(ckbox, ckbox.checked ? whenChecked : whenUnchecked);
		if(c._checkboxMasters[ckbox.id])
			c._cbNotifyGroup(c._checkboxMasters[ckbox.id], ckbox.checked)
		if(c._checkboxItems[ckbox.id])
			c._cbNotifyMaster(c._checkboxItems[ckbox.id])
	}

	c._cbOver = function(ckbox) {
		var whenChecked = c._actions[ckbox.id][0];
		var whenUnchecked = c._actions[ckbox.id][1];
		window.status = c._scripts[ckbox.checked ? whenUnchecked : whenChecked];
	}

	c._mouseOver = function(obj, scriptIndex) {
		if (arguments.length == 1)
			scriptIndex = c._actions[obj.id];
		window.status = c._scripts[scriptIndex];
	}

	c._mouseOut = function() {
		window.status = " ";
		return true;
	}

// from JmolApplet

	c._onloadResetForms = function() {
		// must be evaluated ONLY once -- is this compatible with jQuery?
		if (c._hasResetForms)
			return;
		c._hasResetForms = true;
		c._previousOnloadHandler = window.onload;
		window.onload = function() {
			if (c._buttonCount+c._checkboxCount+c._menuCount+c._radioCount+c._radioGroupCount > 0) {
				var forms = document.forms;
				for (var i = forms.length; --i >= 0; )
					forms[i].reset();
			}
			if (c._previousOnloadHandler)
				c._previousOnloadHandler();
		}
	}

// from JmolApi

	c._getButton = function(appletOrId, script, label, id, title) {
		var appId = c._getIdForControl(appletOrId, script);
		if (appId == null)
			return "";
		//_jmolInitCheck();
		id != undefined && id != null || (id = "jmolButton" + c._buttonCount);
		label != undefined && label != null || (label = script.substring(0, 32));
		++c._buttonCount;
		c._actions[id] = c._addScript(appId, script);
		var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input type='button' name='" + id + "' id='" + id +
						"' value='" + label +
						"' onclick='Jmol.controls._click(this)' onmouseover='Jmol.controls._mouseOver(this);return true' onmouseout='Jmol.controls._mouseOut()' " +
						c._buttonCssText + " /></span>";
		if (Jmol._debugAlert)
			alert(t);
		return Jmol._documentWrite(t);
	}

	c._getCheckbox = function(appletOrId, scriptWhenChecked, scriptWhenUnchecked,
			labelHtml, isChecked, id, title) {

		var appId = c._getIdForControl(appletOrId, scriptWhenChecked);
		if (appId != null)
			appId = c._getIdForControl(appletOrId, scriptWhenUnchecked);
		if (appId == null)
			return "";

		//_jmolInitCheck();
		id != undefined && id != null || (id = "jmolCheckbox" + c._checkboxCount);
		++c._checkboxCount;
		if (scriptWhenChecked == undefined || scriptWhenChecked == null ||
				scriptWhenUnchecked == undefined || scriptWhenUnchecked == null) {
			alert("jmolCheckbox requires two scripts");
			return;
		}
		if (labelHtml == undefined || labelHtml == null) {
			alert("jmolCheckbox requires a label");
			return;
		}
		c._actions[id] = [c._addScript(appId, scriptWhenChecked),c._addScript(appId, scriptWhenUnchecked)];
		var eospan = "</span>"
		var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input type='checkbox' name='" + id + "' id='" + id +
						"' onclick='Jmol.controls._cbClick(this)" +
						"' onmouseover='Jmol.controls._cbOver(this)" +
						";return true' onmouseout='Jmol.controls._mouseOut()' " +
			(isChecked ? "checked='true' " : "")+ c._checkboxCssText + " />"
		if (labelHtml.toLowerCase().indexOf("<td>")>=0) {
			t += eospan
			eospan = "";
		}
		t += "<label for=\"" + id + "\">" + labelHtml + "</label>" +eospan;
		if (Jmol._debugAlert)
			alert(t);
		return Jmol._documentWrite(t);
	}

	c._getCommandInput = function(appletOrId, label, size, id, title, cmd0) {
		var appId = c._getIdForControl(appletOrId, "x");
		if (appId == null)
			return "";
		//_jmolInitCheck();
		id != undefined && id != null || (id = "jmolCmd" + c._cmdCount);
		label != undefined && label != null || (label = "Execute");
		size != undefined && !isNaN(size) || (size = 60);
		cmd0 != undefined || (cmd0 = "help");
		++c._cmdCount;
		var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input name='" + id + "' id='" + id +
						"' size='"+size+"' onkeydown='return Jmol.controls._commandKeyPress(event,\""+id+"\",\"" + appId + "\")' value='" + cmd0 + "'/><input " +
						" type='button' name='" + id + "Btn' id='" + id + "Btn' value = '"+label+"' onclick='Jmol.controls._commandKeyPress(13,\""+id+"\",\"" + appId + "\")' /></span>";
		if (Jmol._debugAlert)
			alert(t);
		return Jmol._documentWrite(t);
	}

	c._getLink = function(appletOrId, script, label, id, title) {
		var appId = c._getIdForControl(appletOrId, script);
		if (appId == null)
			return "";
		//_jmolInitCheck();
		id != undefined && id != null || (id = "jmolLink" + c._linkCount);
		label != undefined && label != null || (label = script.substring(0, 32));
		++c._linkCount;
		var scriptIndex = c._addScript(appId, script);
		var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><a name='" + id + "' id='" + id +
						"' href='javascript:Jmol.controls._click(null,"+scriptIndex+");' onmouseover='Jmol.controls._mouseOver(null,"+scriptIndex+");return true;' onmouseout='Jmol.controls._mouseOut()' " +
						c._linkCssText + ">" + label + "</a></span>";
		if (Jmol._debugAlert)
			alert(t);
		return Jmol._documentWrite(t);
	}

	c._getMenu = function(appletOrId, arrayOfMenuItems, size, id, title) {
		var appId = c._getIdForControl(appletOrId, null);
		var optgroup = null;
		//_jmolInitCheck();
		id != undefined && id != null || (id = "jmolMenu" + c._menuCount);
		++c._menuCount;
		var type = typeof arrayOfMenuItems;
		if (type != null && type == "object" && arrayOfMenuItems.length) {
			var len = arrayOfMenuItems.length;
			if (typeof size != "number" || size == 1)
				size = null;
			else if (size < 0)
				size = len;
			var sizeText = size ? " size='" + size + "' " : "";
			var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><select name='" + id + "' id='" + id +
							"' onChange='Jmol.controls._menuSelected(this)'" +
							sizeText + c._menuCssText + ">";
			for (var i = 0; i < len; ++i) {
				var menuItem = arrayOfMenuItems[i];
				type = typeof menuItem;
				var script = null;
				var text = null;
				var isSelected = null;
				if (type == "object" && menuItem != null) {
					script = menuItem[0];
					text = menuItem[1];
					isSelected = menuItem[2];
				} else {
					script = text = menuItem;
				}
				appId = c._getIdForControl(appletOrId, script);
				if (appId == null)
					return "";
				text == null && (text = script);
				if (script=="#optgroup") {
					t += "<optgroup label='" + text + "'>";
				} else if (script=="#optgroupEnd") {
					t += "</optgroup>";
				} else {
					var scriptIndex = c._addScript(appId, script);
					var selectedText = isSelected ? "' selected='true'>" : "'>";
					t += "<option value='" + scriptIndex + selectedText + text + "</option>";
				}
			}
			t += "</select></span>";
			if (Jmol._debugAlert)
				alert(t);
			return Jmol._documentWrite(t);
		}
	}

	c._getRadio = function(appletOrId, script, labelHtml, isChecked, separatorHtml, groupName, id, title) {
		//_jmolInitCheck();
		if (c._radioGroupCount == 0)
			++c._radioGroupCount;
		groupName || (groupName = "jmolRadioGroup" + (c._radioGroupCount - 1));
		var t = c._radio(appletOrId, script, labelHtml, isChecked, separatorHtml, groupName, (id ? id : groupName + "_" + c._radioCount), title ? title : 0);
		if (t == null)
			return "";
		if (Jmol._debugAlert)
			alert(t);
		return Jmol._documentWrite(t);
	}

	c._getRadioGroup = function(appletOrId, arrayOfRadioButtons, separatorHtml, groupName, id, title) {
		/*

			array: [radio1,radio2,radio3...]
			where radioN = ["script","label",isSelected,"id","title"]

		*/

		//_jmolInitCheck();
		var type = typeof arrayOfRadioButtons;
		if (type != "object" || type == null || ! arrayOfRadioButtons.length) {
			alert("invalid arrayOfRadioButtons");
			return;
		}
		separatorHtml != undefined && separatorHtml != null || (separatorHtml = "&#xa0; ");
		var len = arrayOfRadioButtons.length;
		++c._radioGroupCount;
		groupName || (groupName = "jmolRadioGroup" + (c._radioGroupCount - 1));
		var t = "<span id='"+(id ? id : groupName)+"'>";
		for (var i = 0; i < len; ++i) {
			if (i == len - 1)
				separatorHtml = "";
			var radio = arrayOfRadioButtons[i];
			type = typeof radio;
			var s = null;
			if (type == "object") {
				t += (s = c._radio(appletOrId, radio[0], radio[1], radio[2], separatorHtml, groupName, (radio.length > 3 ? radio[3]: (id ? id : groupName)+"_"+i), (radio.length > 4 ? radio[4] : 0), title));
			} else {
				t += (s = c._radio(appletOrId, radio, null, null, separatorHtml, groupName, (id ? id : groupName)+"_"+i, title));
			}
			if (s == null)
				return "";
		}
		t+="</span>"
		if (Jmol._debugAlert)
			alert(t);
		return Jmol._documentWrite(t);
	}


})(Jmol);
// JmolApi.js -- Jmol user functions  Bob Hanson hansonr@stolaf.edu

// BH 5/29/2014 8:14:06 AM added default command for command input box
// BH 3/10/2014 10:35:25 AM adds Jmol.saveImage(applet)
// BH 1/22/2014 7:31:59 AM Jmol._Image removed -- just never found useful to have
//    a server-side process with only a client-side image. Response time is too slow.
// BH 12/13/2013 8:39:00 AM Jmol.evaulate is DEPRECATED -- use Jmol.evaluateVar
// BH 11/25/2013 6:55:53 AM adds URL flags _USE=, _JAR=, _J2S=
// BH 9/3/2013 5:48:03 PM simplification of Jmol.getAppletHTML()
// BH 5/16/2013 9:01:41 AM checkbox group fix
// BH 1/15/2013 10:55:06 AM updated to default to HTML5 not JAVA

// This file is part of JSmol.min.js. 
// If you do not use that, then along with this file you need several other files. See JSmolCore.js for details.

// default settings are below. Generally you would do something like this:

// jmol = "jmol"
// Info = {.....your settings if not default....}
// Jmol.jmolButton(jmol,....)
// jmol = Jmol.getApplet(jmol, Info)
// Jmol.script(jmol,"....")
// Jmol.jmolLink(jmol,....)
// etc. 
// first parameter is always the applet id, either the string "jmol" or the object defined by Jmol.getApplet()
// no need for waiting to start giving script commands. You can also define a callback function as part of Info.

// see JmolCore.js for details

// BH 8/12/2012 5:15:11 PM added Jmol.getAppletHtml()

;(function (Jmol) {
	var getField = function(key) {
		key = "&" + key + "=";
		return decodeURI(("&" + document.location.search.substring(1) + key).split(key)[1].split("&")[0]);
	}
	Jmol._j2sPath = getField("_J2S");
		// allows URL-line setting of Info.j2sPath
	Jmol._jarFile = getField("_JAR");
		// allows URL-line setting of Info.jarPath and Info.jarFile
	Jmol._use = getField("_USE");
		// allows URL-line setting of Info.use
		// defaults to "HTML5"
		// looking for "_USE=xxxx" 
		// _USE=SIGNED implies JAVA, sets Info.isSigned, and adds "Signed" to applet jar name if necessary

	Jmol.getVersion = function(){return Jmol._jmolInfo.version};

	Jmol.getApplet = function(id, Info, checkOnly) {
		// requires JmolApplet.js and, if JAVA, java/JmolApplet*.jar
		// or if HTML5, then j2s/ subdirectory (core, java, JZ, J)
		/*
		var DefaultInfo = {
			color: "#FFFFFF", // applet object background color, as for older jmolSetBackgroundColor(s)
			width: 300,
			height: 300,
			addSelectionOptions: false,
			serverURL: "http://your.server.here/jsmol.php",
	 console: null,  // div for where the JavaScript console will be.
			defaultModel: "",
			script: null,
			src: null,
			readyFunction: null,
			use: "HTML5",//other options include JAVA, WEBGL//, and IMAGE (removed)
			jarPath: "java",
			jarFile: "JmolApplet0.jar",
			isSigned: false,
			j2sPath: "j2s",
			coverImage: null,     // URL for image to display
			coverTitle: "",       // tip that is displayed before model starts to load
			coverCommand: "",     // Jmol command executed upon clicking image
			deferApplet: false,   // true == the model should not be loaded until the image is clicked
			deferUncover: false,  // true == the image should remain until command execution is complete 
			disableJ2SLoadMonitor: false,
			disableInitialConsole: false,
			debug: false
		};	 

		*/
		return Jmol._Applet._get(id, Info, checkOnly);
	}

	Jmol.getJMEApplet = function(id, Info, linkedApplet, checkOnly) {
		// Java Molecular Editor
		// requires JmolJME.js and jme/ subdirectory
		/*
		var DefaultInfo = {
			width: 300,
			height: 300,
			jarPath: "jme",
			jarFile: "JME.jar",
			use: "HTML", // or JAVA
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
		*/
		return Jmol._JMEApplet._get(id, Info, linkedApplet, checkOnly);
	}

	Jmol.getJSVApplet = function(id, Info, checkOnly) {
		// JSpecView
		// requires JmolJSV.js and, if JAVA, either JSpecViewApplet.jar or JSpecViewAppletSigned.jar
		// or if HTML5, then j2s/ subdirectory (core, java, JZ, J, JSV)
		/*
		var DefaultInfo = {
			width: 500,
			height: 300,
			debug: false,
			jarPath: ".",
			jarFile: "JSpecViewApplet.jar", // or "JSpecViewAppletSigned.jar"
			uee: "HTML5", // or JAVA
			isSigned: false,
			initParams: null,
			readyFunction: null,
			script: null
		};
		*/
		return Jmol._JSVApplet._get(id, Info, checkOnly);
	}	


////////////////// scripting ///////////////////

	Jmol.loadFile = function(applet, fileName, params){
		applet._loadFile(fileName, params);
	}

	Jmol.script = function(applet, script) {
		if (applet._checkDeferred(script)) 
			return;
		applet._script(script);
	}

/**
 * returns false if cannot check, empty string if OK, or error message if not OK
 */	
	Jmol.scriptCheck = function(applet, script) {
		return applet && applet._scriptCheck && applet._ready && applet._scriptCheck(script);
	}

	Jmol.scriptWait = function(applet, script) {
		return applet._scriptWait(script);
	}

	Jmol.scriptEcho = function(applet, script) {
		return applet._scriptEcho(script);
	}

	Jmol.scriptMessage = function(applet, script) {
		return applet._scriptMessage(script);
	}

	Jmol.scriptWaitOutput = function(applet, script) {
		return applet._scriptWait(script);
	}

	Jmol.scriptWaitAsArray = function(applet, script) {
		return applet._scriptWaitAsArray(script);
	}

	Jmol.search = function(applet, query, script) {
		applet._search(query, script);
	}

////////////////// "get" methods ///////////////////


	Jmol.evaluateVar = function(applet,expr) {
		return applet._evaluate(expr);
	}

	// DEPRECATED -- use Jmol.evaluateVar
	Jmol.evaluate = function(applet,molecularMath) {
		return applet._evaluateDEPRECATED(molecularMath);
	}

	// optional Info here	
	Jmol.getAppletHtml = function(applet, Info) {
		if (Info) {
			var d = Jmol._document;
			Jmol._document = null;
			applet = Jmol.getApplet(applet, Info);
			Jmol._document = d;
		}  
		return applet._code;
	}

	Jmol.getPropertyAsArray = function(applet,sKey,sValue) {
		return applet._getPropertyAsArray(sKey,sValue);
	}

	Jmol.getPropertyAsJavaObject = function(applet,sKey,sValue) {
		return applet._getPropertyAsJavaObject(sKey,sValue);
	}

	Jmol.getPropertyAsJSON = function(applet,sKey,sValue) {
		return applet._getPropertyAsJSON(sKey,sValue);
	}

	Jmol.getPropertyAsString = function(applet,sKey,sValue) {
		return applet._getPropertyAsString(sKey,sValue);
	}

	Jmol.getStatus = function(applet,strStatus) {
		return applet._getStatus(strStatus);
	}


////////////////// general methods ///////////////////

	Jmol.resizeApplet = function(applet,size) {
		return applet._resizeApplet(size);
	}

	Jmol.restoreOrientation = function(applet,id) {
		return applet._restoreOrientation(id);
	}

	Jmol.restoreOrientationDelayed = function(applet,id,delay) {
		return applet._restoreOrientationDelayed(id,delay);
	}

	Jmol.saveOrientation = function(applet,id) {
		return applet._saveOrientation(id);
	}

	Jmol.say = function(msg) {
		alert(msg);
	}

//////////// console functions /////////////

	Jmol.clearConsole = function(applet) {
		applet._clearConsole();
	}

	Jmol.getInfo = function(applet) {
		return applet._info;
	}

	Jmol.setInfo = function(applet, info, isShown) {
		applet._info = info;
		if (arguments.length > 2)
			applet._showInfo(isShown);
	}

	Jmol.showInfo = function(applet, tf) {
		applet._showInfo(tf);
	}

	Jmol.show2d = function(applet, tf) {
		// only when JME or JSME is synced with Jmol
		applet._show2d(tf);
	}


//////////// controls and HTML /////////////


	Jmol.jmolBr = function() {
		return Jmol._documentWrite("<br />");
	}

	Jmol.jmolButton = function(appletOrId, script, label, id, title) {
		return Jmol.controls._getButton(appletOrId, script, label, id, title);
	}

	Jmol.jmolCheckbox = function(appletOrId, scriptWhenChecked, scriptWhenUnchecked,
			labelHtml, isChecked, id, title) {
		return Jmol.controls._getCheckbox(appletOrId, scriptWhenChecked, scriptWhenUnchecked,
			labelHtml, isChecked, id, title);
	}


	Jmol.jmolCommandInput = function(appletOrId, label, size, id, title, cmd0) {
		return Jmol.controls._getCommandInput(appletOrId, label, size, id, title, cmd0);
	}

	Jmol.jmolHtml = function(html) {
		return Jmol._documentWrite(html);
	}

	Jmol.jmolLink = function(appletOrId, script, label, id, title) {
		return Jmol.controls._getLink(appletOrId, script, label, id, title);
	}

	Jmol.jmolMenu = function(appletOrId, arrayOfMenuItems, size, id, title) {
		return Jmol.controls._getMenu(appletOrId, arrayOfMenuItems, size, id, title);
	}

	Jmol.jmolRadio = function(appletOrId, script, labelHtml, isChecked, separatorHtml, groupName, id, title) {
		return Jmol.controls._getRadio(appletOrId, script, labelHtml, isChecked, separatorHtml, groupName, id, title);
	}

	Jmol.jmolRadioGroup = function (appletOrId, arrayOfRadioButtons, separatorHtml, groupName, id, title) {
		return Jmol.controls._getRadioGroup(appletOrId, arrayOfRadioButtons, separatorHtml, groupName, id, title);
	}

	Jmol.setCheckboxGroup = function(chkMaster, chkBoxes) {
		// chkBoxes can be an array or any number of additional string arguments
		Jmol.controls._cbSetCheckboxGroup(chkMaster, chkBoxes, arguments);
	}

	Jmol.setDocument = function(doc) {

		// If doc is null or 0, Jmol.getApplet() will still return an Object, but the HTML will
		// put in applet._code and not written to the page. This can be nice, because then you 
		// can still refer to the applet, but place it on the page after the controls are made. 
		//
		// This really isn't necessary, though, because there is a simpler way: Just define the 
		// applet variable like this:
		//
		// jmolApplet0 = "jmolApplet0"
		//
		// and then, in the getApplet command, use
		//
		// jmolapplet0 = Jmol.getApplet(jmolApplet0,....)
		// 
		// prior to this, "jmolApplet0" will suffice, and after it, the Object will work as well
		// in any button creation 
		//		 
		//  Bob Hanson 25.04.2012

		Jmol._document = doc;
	}

	Jmol.setXHTML = function(id) {
		Jmol._isXHTML = true;
		Jmol._XhtmlElement = null;
		Jmol._XhtmlAppendChild = false;
		if (id){
			Jmol._XhtmlElement = document.getElementById(id);
			Jmol._XhtmlAppendChild = true;
		}
	}

	////////////////////////////////////////////////////////////////
	// Cascading Style Sheet Class support
	////////////////////////////////////////////////////////////////

	// BH 4/25 -- added text option. setAppletCss(null, "style=\"xxxx\"")
	// note that since you must add the style keyword, this can be used to add any attribute to these tags, not just css. 

	Jmol.setAppletCss = function(cssClass, text) {
		cssClass != null && (Jmol.controls._appletCssClass = cssClass);
		Jmol.controls._appletCssText = text ? text + " " : cssClass ? "class=\"" + cssClass + "\" " : "";
	}

	Jmol.setButtonCss = function(cssClass, text) {
		cssClass != null && (Jmol.controls._buttonCssClass = cssClass);
		Jmol.controls._buttonCssText = text ? text + " " : cssClass ? "class=\"" + cssClass + "\" " : "";
	}

	Jmol.setCheckboxCss = function(cssClass, text) {
		cssClass != null && (Jmol.controls._checkboxCssClass = cssClass);
		Jmol.controls._checkboxCssText = text ? text + " " : cssClass ? "class=\"" + cssClass + "\" " : "";
	}

	Jmol.setRadioCss = function(cssClass, text) {
		cssClass != null && (Jmol.controls._radioCssClass = cssClass);
		Jmol.controls._radioCssText = text ? text + " " : cssClass ? "class=\"" + cssClass + "\" " : "";
	}

	Jmol.setLinkCss = function(cssClass, text) {
		cssClass != null && (Jmol.controls._linkCssClass = cssClass);
		Jmol.controls._linkCssText = text ? text + " " : cssClass ? "class=\"" + cssClass + "\" " : "";
	}

	Jmol.setMenuCss = function(cssClass, text) {
		cssClass != null && (Jmol.controls._menuCssClass = cssClass);
		Jmol.controls._menuCssText = text ? text + " ": cssClass ? "class=\"" + cssClass + "\" " : "";
	}

	Jmol.setAppletSync = function(applets, commands, isJmolJSV) {
		Jmol._syncedApplets = applets;   // an array of appletIDs
		Jmol._syncedCommands = commands; // an array of commands; one or more may be null 
		Jmol._syncedReady = {};
		Jmol._isJmolJSVSync = isJmolJSV;
	}

	/*
	Jmol._grabberOptions = [
		["$", "NCI(small molecules)"],
		[":", "PubChem(small molecules)"],
		["=", "RCSB(macromolecules)"]
	];
	*/

	Jmol.setGrabberOptions = function(options) {
		Jmol._grabberOptions = options;
	}

	Jmol.setAppletHtml = function (applet, divid) {
		if (!applet._code) 
			return;
		Jmol.$html(divid, applet._code);
		if (applet._init && !applet._deferApplet)
			applet._init();
	}

	Jmol.coverApplet = function(applet, doCover) {
		if (applet._cover)
			applet._cover(doCover);
	}

	Jmol.setFileCaching = function(applet, doCache) {
		if (applet) {
			applet._cacheFiles = doCache;
		} else {
			Jmol.fileCache = (doCache ? {} : null);
		}
	}  

	Jmol.updateView = function(applet, param1, param2) {
		applet._updateView(param1, param2);
	}

	Jmol.getChemicalInfo = function(appletOrIdentifier, what, fCallback) {
		what || (what = "name");
		if (typeof applet_or_Identifier != "string") 
			appletOrIdentifier = appletOrIdentifier._getSmiles();
		return Jmol._getNCIInfo(appletOrIdentifier, what, fCallback);
	}

	Jmol.saveImage = function(app) {
		// see: https://svgopen.org/2010/papers/62-From_SVG_to_Canvas_and_Back/index.html
		// From SVG to Canvas and Back
		// Samuli Kaipiainen University of Helsinki, Department of Computer Science samuli.kaipiainen@cs.helsinki.fi
		// Matti Paksula University of Helsinki, Department of Computer Science matti.paksula@cs.helsinki.fi
		switch (app._viewType) {
		case "Jmol":
			app._script("write PNGJ \"" + app._id + ".png\"");
			break;
		case "JSV":
			app._script("write PDF");
			break;
		case "JME":
			app._script("print");
			break;
		}
	}
		
})(Jmol);
// j2sjmol.js 

// Java programming notes by Bob Hanson:
//   
//   There are a few motifs to avoid when optimizing Java code to work smoothly
//   with the J2S compiler:
//   
//   arrays: 
//   
// 1. an array with null elements cannot be typed and must be avoided.
// 2. instances of Java "instance of" involving arrays must be found and convered to calls to Clazz.isA...
// 3. new int[n][] must not be used. Use instead JU.AU.newInt2(n);
// 4. new int[] { 1, 2, 3 } has problems because it creates simply [ ] and not IntArray32
//   
//   numbers:
//   
// 1. Remember that EVERY number in JavaScript is a double -- doesn't matter if it is in IntArray32 or not. 
// 2. You cannot reliably use Java long, because doubles consume bits for the exponent which cannot be tested.
// 3. Bit 31 of an integer is unreliable, since (int) -1 is now  , not just 0zFFFFFFFF, and 
//    FFFFFFFF + 1 = 100000000, not 0. In JavaScript, 0xFFFFFFFF is 4294967295, not -1.
//    This means that writeInt(b) will fail if b is negative. What you need is instead
//    writeInt((int)(b & 0xFFFFFFFFl) so that JavaScript knocks off the high bits explicitly. 
//
//   general:
//
// 1. j2sRequireImport xxxx is needed if xxxx is a method used in a static function
// 2. URL.getContent() is not supported. Use other means based on URL.toString()
// 3. It is critical for performance to avoid any significant amount of function overloading.
//    In particular, methods such as xxx(int a, int b) and xxx(float a, int b) MUST be renamed,
//    because JavaScript only has Number, and there is absolutely no way to tell these apart.
//    It's probably bad Java programming, anyway.
// 4. Calls to super(...) can almost always be avoided. These trigger the SAEM
//    (searchAndExecuteMethod) call, and it is very destructive to performance.
//    Just find another way to do it.   

 // NOTES by Bob Hanson: 
 
 // This version of j2slib requires jQuery and works in both Chrome and MSIE locally.

 // J2S class changes:

 // BH 6/1/2014 10:58:46 AM fix for Clazz.isAP() not working
 // BH 5/26/2014 5:19:29 PM removing superConstructor call in creating Enum constants
 // BH 4/1/2014 7:55:54 PM removing all $fz references and instances where sub/super classes have same private function names
 // BH 4/1/2014 4:47:30 PM all $_X removed; this is taken care of by Google Closure Compiler
 // BH 4/1/2014 6:40:08 AM removing ClassLoader -- equals Clazz._Loader
 // BH 4/1/2014 6:40:08 AM removing ClassLoaderProgressMonitor -- equals ClazzLoaderProgressMonitor
 // BH 4/1/2014 6:17:21 AM removing Class  -- only used for "Class.forName" in Jmol, which ANT will now change to "Clazz._4Name"
 // BH 3/7/2014 9:05:06 AM Array.prototype.toString should not be aliased. -- http://sourceforge.net/p/jmol/bugs/560/ with Google Visualization

 // BH 1/30/2014 12:54:22 PM gave all field variables prefix underscore. This allows Google Closure Compiler to skip them.  
 // BH 12/3/2013 3:39:57 PM window["j2s.lib"].base implemented
 // BH 12/1/2013 5:34:21 AM removed ClazzLoaderProgressMonitor.initialize and all Clazz.event business; handled by Jmol.clearVars()
 // BH 11/30/2013 12:43:58 PM adding Clazz.arrayIs() -- avoids Number.constructor.toString() infinite recursion
 // BH 11/29/2013 6:33:51 AM adding Clazz._profiler -- reports use of SAEM
 // BH 11/10/2013 9:02:20 AM fixing fading in MSIE  
 // BH 11/3/2013 7:21:39 AM additional wrapping functions for better compressibility
 // BH 10/30/2013 8:10:58 AM added getClass().getResource() -- returning a relative string, not a URL
 // BH 10/30/2013 6:43:00 AM removed second System def and added System.$props and default System.property "line.separator" 
 // BH 6/15/2013 8:02:07 AM corrections to Class.isAS to return true if first element is null
 // BH 6/14/2013 4:41:09 PM corrections to Clazz.isAI and related methods to include check for null object
 // BH 3/17/2013 11:54:28 AM adds stackTrace for ERROR 

 // BH 3/13/2013 6:58:26 PM adds Clazz.clone(me) for BS clone 
 // BH 3/12/2013 6:30:53 AM fixes Clazz.exceptionOf for ERROR condition trapping
 // BH 3/2/2013 9:09:53 AM delete globals c$ and $fz
 // BH 3/2/2013 9:10:45 AM optimizing defineMethod using "look no further" "@" parameter designation (see "\\@" below -- removed 3/23/13)
 // BH 2/27/2013 optimizing getParamsType for common cases () and (Number)
 // BH 2/27/2013 optimizing SAEM delegation for hashCode and equals -- disallows overloading of equals(Object)
 
 // BH 2/23/2013 found String.replaceAll does not work -- solution was to never call it.
 // BH 2/9/2013 9:18:03 PM Int32Array/Float64Array fixed for MSIE9
 // BH 1/25/2013 1:55:31 AM moved package.js from j2s/java to j2s/core 
 // BH 1/17/2013 4:37:17 PM String.compareTo() added
 // BH 1/17/2013 4:52:22 PM Int32Array and Float64Array may not have .prototype.sort method
 // BH 1/16/2013 6:20:34 PM Float64Array not available in Safari 5.1
 // BH 1/14/2013 11:28:58 PM  Going to all doubles in JavaScript (Float64Array, not Float32Array)
 //   so that (new float[] {13.48f})[0] == 13.48f, effectively

 // BH 1/14/2013 12:53:41 AM  Fix for Opera 10 not loading any files
 // BH 1/13/2013 11:50:11 PM  Fix for MSIE not loading (nonbinary) files locally
 
 // BH 12/1/2012 9:52:26 AM Compiler note: Thread.start() cannot be executed within the constructor;
 
 // BH 11/24/2012 11:08:39 AM removed unneeded sections
 // BH 11/24/2012 10:23:22 AM  all XHR uses sync loading (ClazzLoader.setLoadingMode)
 // BH 11/21/2012 7:30:06 PM 	if (base != null)	map["@" + pkg] = base;  critical for multiple applets

 // BH 10/8/2012 3:27:41 PM         if (clazzName.indexOf("Array") >= 0) return "Array"; in Clazz.getClassName for function
 // BH removed Clazz.ie$plit = "\\2".split (/\\/).length == 1; unnecessary; using RegEx slows process significantly in all browsers
 // BH 10/6/12 added Int32Array, Float32Array, newArrayBH, upgraded java.lang and java.io
 // BH added Integer.bitCount in core.z.js
 // BH changed alert to Clazz.alert in java.lang.Class.js *.ClassLoader.js, java.lang.thread.js
 // BH removed toString from innerFunctionNames due to infinite recursion
 // BH note: Logger.error(null, e) does not work -- get no constructor for (String) (TypeError)
 // BH added j2s.lib.console
 // BH allowed for alias="."
 // BH removed alert def --> Clazz.alert
 // BH added wrapper at line 2856 
 // BH newArray fix at line 2205
 // BH System.getProperty fix at line 6693
 // BH added Enum .value() method at line 2183
 // BH added System.getSecurityManager() at end
 // BH added String.contains() at end
 // BH added System.gc() at end
 // BH added Clazz.exceptionOf = updated
 // BH added String.getBytes() at end
 

LoadClazz = function() {

if (!window["j2s.clazzloaded"])
	window["j2s.clazzloaded"] = false;

if (window["j2s.clazzloaded"])return;

window["j2s.clazzloaded"] = true;

window["j2s.object.native"] = true;

 // Clazz changes:

 /* http://j2s.sf.net/ *//******************************************************************************
 * Copyright (c) 2007 java2script.org and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Zhou Renjian - initial API and implementation
 *****************************************************************************/
/*******
 * @author zhou renjian
 * @create Nov 5, 2005
 *******/
 
//if (window["Clazz"] == null) {

/*
 * The following *-# are used to compress the JavaScript file into small file.
 * For more details, please read /net.sf.j2s.lib/build/build.xml
 */
/*-#
 # _x_CLASS_NAME__ -> C$N
 # _x_PKG_NAME__ -> P$N
 #
 # clazzThis -> Tz
 # objThis -> To
 # clazzHost -> Hz
 # hostThis -> Th
 # hostSuper -> Sh
 # clazzFun -> Fc
 # clazzName -> Nc
 # funName -> Nf
 # funBody -> Bf
 # objType -> oT
 #  
 # qClazzName ->Nq
 #-*/
/**
 * Class Clazz. All the methods are static in this class.
 */
/* static */
/*Class = */ Clazz = function () {};

;(function(Clazz) {

Clazz.__debuggingBH = false;
Clazz._globals = ["j2s.clazzloaded", "j2s.object.native"];
Clazz.setGlobal = function(a, v) {
	Clazz._globals.push(a);
	window[a] = v;
}

Clazz.getGlobals = function() {
	return Clazz._globals.sort().join("\n");
}

// BH Clazz.getProfile monitors exactly what is being delegated with SAEM,
// which could be a bottle-neck for function calling.

// Jmol.getProfile()

Clazz._profile = (window["j2s.doProfile"]  && self.JSON ? {} : null);

Clazz.getProfile = function() {
	var s = "";
	if (Clazz._profile) {
		var l = [];
		for (var i in Clazz._profile) {
			var n = "" + Clazz._profile[i];
			l.push("        ".substring(n.length) + n + "\t" + i);
		}
		s = l.sort().reverse().join("\r\n");
		Clazz._profile = {};
	}
	return s;
}

Clazz.addProfile = function(c, f, p) {
	var s = c.__CLASS_NAME__ + " " + f + " " + JSON.stringify(p);
	Clazz._profile[s] || (Clazz._profile[s] = 0);
	Clazz._profile[s]++;
}

NullObject = function () {};

/* protected */
Clazz._supportsNativeObject = window["j2s.object.native"];

if (Clazz._supportsNativeObject) {
	Clazz._O = function () {};
	Clazz._O.__CLASS_NAME__ = "Object";
	Clazz._O["getClass"] = function () { return Clazz._O; }; 
} else {
	Clazz._O = Object;
}

Clazz.Console = {};
Clazz.dateToString = Date.prototype.toString;


Clazz.getSignature = function(proto, name, func, isNew) {
	// BH pointing to signatures based on number of parameters
	// would only make SAEM somewhat faster; not worth it? 
	// better to just avoid SAEM altogether
	return (isNew ? proto[name] = func : proto[name]);
/*
	if (!isNew) return proto[name];
	if (proto[name])
		func.sigs = proto[name].sigs;
	proto[name] = func;
	if (!func.sigs) func.sigs = [];
	var n = (func.arguments ? func.arguments.length : 0);
	if (func.sigs[n])
		func.sigs[n] == -1; // overloaded for this number of parameters
	else                  // unique for this function
		func.sigs[n] = func;
	return func;
*/
};


Clazz.addProto = function(proto, name, func) {
	Clazz.getSignature(proto, name, func, true); // BH
};

;(function(proto) {
	Clazz.addProto(proto, "equals", function (obj) {
		return this == obj;
	});

	Clazz.addProto(proto, "hashCode", function () {
		try {
			return this.toString ().hashCode ();
		} catch (e) {
			var str = ":";
			for (var s in this) {
				str += s + ":"
			}
			return str.hashCode ();
		}
	});

	Clazz.addProto(proto, "getClass", function () {
	 return Clazz.getClass (this);
	});

	Clazz.addProto(proto, "clone", function () {
		return Clazz.clone(this);
	});

	Clazz.clone = function(me) {
		// BH allows @j2sNative access without super constructor
		var o = new me.constructor ();
		for (var i in me)
			o[i] = me[i];
		return o;
	}
/*
 * Methods for thread in Object
 */
	Clazz.addProto(proto, "finalize", function () {});
	Clazz.addProto(proto, "notify", function () {});
	Clazz.addProto(proto, "notifyAll", function () {});
	Clazz.addProto(proto, "wait", function () {});

	Clazz.addProto(proto, "to$tring", Object.prototype.toString);
	Clazz.addProto(proto, "toString", function () {
		if (this.__CLASS_NAME__ != null) {
			return "[" + this.__CLASS_NAME__ + " object]";
		} else {
			return this.to$tring.apply (this, arguments);
		}
	});

	Clazz._extendedObjectMethods = [
			"equals", "hashCode", "getClass", "clone", "finalize", "notify", "notifyAll", "wait", "to$tring", "toString"
	];

})(Clazz._O.prototype);

Clazz.extendJO = function(c, name) {  
	if (name)
		c.__CLASS_NAME__ = c.prototype.__CLASS_NAME__ = name;
	if (Clazz._supportsNativeObject) {
		for (var i = 0; i < Clazz._extendedObjectMethods.length; i++) {
			var p = Clazz._extendedObjectMethods[i];
			Clazz.getSignature(c.prototype, p, Clazz._O.prototype[p], true);
		}
	}
};

/**
 * Try to fix bug on Safari
 */
//InternalFunction = Object;

Clazz.extractClassName = function(clazzStr) {
	// [object Int32Array]
	var clazzName = clazzStr.substring (1, clazzStr.length - 1);
	return (clazzName.indexOf("Array") >= 0 ? "Array" // BH -- for Float64Array and Int32Array
		: clazzName.indexOf ("object ") >= 0 ? clazzName.substring (7) // IE
		: clazzName);
}
/**
 * Return the class name of the given class or object.
 *
 * @param clazzHost given class or object
 * @return class name
 */
/* public */
Clazz.getClassName = function (obj) {
	if (obj == null) {
		/* 
		 * null is always treated as Object.
		 * But what about "undefined"?
		 */
		return "NullObject";
	}

	if (obj instanceof Clazz.CastedNull)
		return obj.clazzName;
	switch(typeof obj) {
	case "number":
		return "Number";
	case "boolean":
		return "Boolean";
	case "string":
		/* 
		 * Always treat the constant string as String object.
		 * This will be compatiable with Java String instance.
		 */
		return "String";
	case "function":
		if (obj.__CLASS_NAME__ != null)
			return (arguments[1] ? obj.__CLASS_NAME__ : "Class"); /* user defined class name */
		var s = obj.toString();
		var idx0 = s.indexOf("function");
		if (idx0 < 0)
			return (s.charAt(0) == '[' ? Clazz.extractClassName(s) : s.replace(/[^a-zA-Z0-9]/g, ''));
		var idx1 = idx0 + 8;
		var idx2 = s.indexOf ("(", idx1);
		if (idx2 < 0)
			return "Object";
		s = s.substring (idx1, idx2);
		if (s.indexOf("Array") >= 0)
			return "Array";  // BH -- for Float64Array and Int32Array
		s = s.replace (/^\s+/, "").replace (/\s+$/, ""); // .trim ()
		return (s == "anonymous" || s == "" ? "Function" : s);
		 // BH -- for general functions, clazzName may be ""
	case "object":
		if (obj.__CLASS_NAME__ != null) // user defined class name
			return obj.__CLASS_NAME__;
		if (obj.constructor == null)
			return "Object"; // For HTML Element in IE
		if (obj.constructor.__CLASS_NAME__ == null) {
			if (obj instanceof Number)
				return "Number";
			if (obj instanceof Boolean)
				return "Boolean";
			if (obj instanceof Array)
				return "Array";
			var s = obj.toString();
			if (s.charAt(0) == '[')
				return Clazz.extractClassName(s);
		}
	}
	return Clazz.getClassName (obj.constructor, true);
};
/**
 * Return the class of the given class or object.
 *
 * @param clazzHost given class or object
 * @return class name
 */
/* public */
Clazz.getClass = function (clazzHost) {
	if (clazzHost == null) {
		/* 
		 * null is always treated as Object.
		 * But what about "undefined"?
		 */
		return Clazz._O;
	}
	if (typeof clazzHost == "function") {
		return clazzHost;
	} else {
		var clazzName = null;
		var obj = clazzHost;
		if (obj instanceof Clazz.CastedNull) {
			clazzName = obj.clazzName;
		} else {
			var objType = typeof obj;
			if (objType == "string") {
				return String;
			} else if (typeof obj == "object") {
				/* user defined class name */
				if (obj.__CLASS_NAME__ != null) {
					clazzName = obj.__CLASS_NAME__;
				} else if (obj.constructor == null) {
					return Clazz._O; // Is it safe?
				} else {
					return obj.constructor;
				}
			}
		}
		if (clazzName != null) {
			return Clazz.evalType (clazzName, true);
		} else {
			return obj.constructor;
		}
	}
};

/*
 * Be used to copy members of class
 */
/* protected */
/*-# extendsProperties -> eP #-*/
Clazz.extendsProperties = function (hostThis, hostSuper) {
	for (var o in hostSuper) {
		if (o != "b$" && o != "prototype" && o != "superClazz"
				&& o != "__CLASS_NAME__" && o != "implementz"
				&& !Clazz.checkInnerFunction (hostSuper, o)) {
			hostThis[o] = hostSuper[o];
		}
	}
};

/* private */
/*-# checkInnerFunction -> cIF #-*/
Clazz.checkInnerFunction = function (hostSuper, funName) {
	for (var k = 0; k < Clazz.innerFunctionNames.length; k++) {
		if (funName == Clazz.innerFunctionNames[k] && 
				Clazz._innerFunctions[funName] === hostSuper[funName]) {
			return true;
		}
	}
	return false;
};

/*
 * Be used to copy members of interface
 */
/* protected */
/*-# implementsProperties -> ip #-*/
Clazz.implementsProperties = function (hostThis, hostSuper) {
	for (var o in hostSuper) {
		if (o != "b$" && o != "prototype" && o != "superClazz"
				&& o != "__CLASS_NAME__" && o != "implementz") {
			if (typeof hostSuper[o] == "function") {
				/*
				 * static final member of interface may be a class, which may
				 * be function.
				 */
				if (Clazz.checkInnerFunction (hostSuper, o)) {
					continue;
				}
			}
			hostThis[o] = hostSuper[o];
			hostThis.prototype[o] = hostSuper[o];
		}
	}
};

/*-# args4InheritClass -> aIC #-*/
Clazz.args4InheritClass = function () {};
Clazz.inheritArgs = new Clazz.args4InheritClass ();

/**
 * Inherit class with "extends" keyword and also copy those static members. 
 * Example, as in Java, if NAME is a static member of ClassA, and ClassB 
 * extends ClassA then ClassB.NAME can be accessed in some ways.
 *
 * @param clazzThis child class to be extended
 * @param clazzSuper super class which is inherited from
 * @param objSuper super class instance
 */
/* protected */
/*-#
 # inheritClass -> xic 
 #
 # objSuper -> oSp
 #-*/
Clazz.inheritClass = function (clazzThis, clazzSuper, objSuper) {
	//var thisClassName = Clazz.getClassName (clazzThis);
	Clazz.extendsProperties (clazzThis, clazzSuper);
	if (Clazz.isClassUnloaded (clazzThis)) {
		// Don't change clazzThis.protoype! Keep it!
	} else if (objSuper != null) {
		// ! Unsafe reference prototype to an instance!
		// Feb 19, 2006 --josson
		// OK for this reference to an instance, as this is anonymous instance,
		// which is not referenced elsewhere.
		// March 13, 2006
		clazzThis.prototype = objSuper; 
	} else if (clazzSuper !== Number) {
		clazzThis.prototype = new clazzSuper (Clazz.inheritArgs);
	} else { // Number
		clazzThis.prototype = new Number ();
	}
	clazzThis.superClazz = clazzSuper;
	/*
	 * Is it necessary to reassign the class name?
	 * Mar 10, 2006 --josson
	 */
	//clazzThis.__CLASS_NAME__ = thisClassName;
	clazzThis.prototype.__CLASS_NAME__ = clazzThis.__CLASS_NAME__;
};

/**
 * Implementation of Java's keyword "implements".
 * As in JavaScript there are on "implements" keyword implemented, a property
 * of "implementz" is added to the class to record the interfaces the class
 * is implemented.
 * 
 * @param clazzThis the class to implement
 * @param interfacez Array of interfaces
 */
/* public */
Clazz.implementOf = function (clazzThis, interfacez) {
	if (arguments.length >= 2) {
		if (clazzThis.implementz == null) {
			clazzThis.implementz = new Array ();
		}
		var impls = clazzThis.implementz;
		if (arguments.length == 2) {
			if (typeof interfacez == "function") {
				impls[impls.length] = interfacez;
				Clazz.implementsProperties (clazzThis, interfacez);
			} else if (interfacez instanceof Array) {
				for (var i = 0; i < interfacez.length; i++) {
					impls[impls.length] = interfacez[i];
					Clazz.implementsProperties (clazzThis, interfacez[i]);
				}
			}
		} else {
			for (var i = 1; i < arguments.length; i++) {
				impls[impls.length] = arguments[i];
				Clazz.implementsProperties (clazzThis, arguments[i]);
			}
		}
	}
};

/**
 * TODO: More should be done for interface's inheritance
 */
/* public */
Clazz.extendInterface = Clazz.implementOf;

/* protected */
/*-#
 # equalsOrExtendsLevel -> eOE 
 #
 # clazzAncestor -> anc
 #-*/
Clazz.equalsOrExtendsLevel = function (clazzThis, clazzAncestor) {
	if (clazzThis === clazzAncestor) {
		return 0;
	}
	if (clazzThis.implementz != null) {
		var impls = clazzThis.implementz;
		for (var i = 0; i < impls.length; i++) {
			var level = Clazz.equalsOrExtendsLevel (impls[i], clazzAncestor);
			if (level >= 0) {
				return level + 1;
			}
		}
	}
	return -1;
};

/* protected */
/*-#
 # getInheritedLevel -> gIL 
 #
 # clazzBase -> bs
 # clazzTarget -> tg
 #-*/
Clazz.getInheritedLevel = function (clazzTarget, clazzBase) {
	if (clazzTarget === clazzBase) {
		return 0;
	}
	var isTgtStr = (typeof clazzTarget == "string");
	var isBaseStr = (typeof clazzBase == "string");
	if ((isTgtStr && ("void" == clazzTarget || "unknown" == clazzTarget)) 
			|| (isBaseStr && ("void" == clazzBase 
					|| "unknown" == clazzBase))) {
		return -1;
	}
	/*
	 * ? The following lines are confusing
	 * March 10, 2006
	 */
	if ((isTgtStr && "NullObject" == clazzTarget) 
			|| NullObject === clazzTarget) {
		if (clazzBase !== Number && clazzBase !== Boolean
				&& clazzBase !== NullObject) {
			return 0;
		}
	}
	if (isTgtStr) {
		clazzTarget = Clazz.evalType (clazzTarget);
	}
	if (isBaseStr) {
		clazzBase = Clazz.evalType (clazzBase);
	}
	if (clazzBase == null || clazzTarget == null) {
		return -1;
	}
	var level = 0;
	var zzalc = clazzTarget; // zzalc <--> clazz
	while (zzalc !== clazzBase && level < 10) {
		/* maybe clazzBase is interface */
		if (zzalc.implementz != null) {
			var impls = zzalc.implementz;
			for (var i = 0; i < impls.length; i++) {
				var implsLevel = Clazz.equalsOrExtendsLevel (impls[i], 
						clazzBase);
				if (implsLevel >= 0) {
					return level + implsLevel + 1;
				}
			}
		}

		zzalc = zzalc.superClazz;
		if (zzalc == null) {
			if (clazzBase === Object || clazzBase === Clazz._O) {
				/*
				 * getInheritedLevel(String, CharSequence) == 1
				 * getInheritedLevel(String, Object) == 1.5
				 * So if both #test(CharSequence) and #test(Object) existed,
				 * #test("hello") will correctly call #test(CharSequence)
				 * insted of #test(Object).
				 */
				return level + 1.5; // 1.5! Special!
			} else {
				return -1;
			}
		}
		level++;
	}
	return level;
};


/**
 * Implements Java's keyword "instanceof" in JavaScript's way.
 * As in JavaScript part of the object inheritance is implemented in only-
 * JavaScript way.
 *
 * @param obj the object to be tested
 * @param clazz the class to be checked
 * @return whether the object is an instance of the class
 */
/* public */
Clazz.instanceOf = function (obj, clazz) {
	if (obj == null) {
		return clazz == false; // should usually false
	}
	if (clazz == null) {
		return false;
	}
	if (obj instanceof clazz) {
		return true;
	} else {
		/*
		 * To check all the inherited interfaces.
		 */
		var clazzName = Clazz.getClassName (obj);
		return Clazz.getInheritedLevel (clazzName, clazz) >= 0;
	}
};

/**
 * Call super method of the class. 
 * The same effect as Java's expression:
 * <code> super.* () </code>
 * 
 * @param objThis host object
 * @param clazzThis class of declaring method scope. It's hard to determine 
 * which super class is right class for "super.*()" call when it's in runtime
 * environment. For example,
 * 1. ClasssA has method #run()
 * 2. ClassB extends ClassA overriding method #run() with "super.run()" call
 * 3. ClassC extends ClassB
 * 4. objC is an instance of ClassC
 * Now we have to decide which super #run() method is to be invoked. Without
 * explicit clazzThis parameter, we only know that objC.getClass() is ClassC 
 * and current method scope is #run(). We do not known we are in scope 
 * ClassA#run() or scope of ClassB#run(). if ClassB is given, Clazz can search
 * all super methods that are before ClassB and get the correct super method.
 * This is the reason why there must be an extra clazzThis parameter.
 * @param funName method name to be called
 * @param funParams Array of method parameters
 */
/* public */
Clazz.superCall = function (objThis, clazzThis, funName, funParams) {
	var fx = null;
	var i = -1;
	var clazzFun = objThis[funName];

	if (clazzFun != null) {
		if (clazzFun.claxxOwner != null) { 
			// claxxOwner is a mark for methods that is single.
			if (clazzFun.claxxOwner !== clazzThis) {
				// This is a single method, call directly!
				fx = clazzFun;
			}
		} else if (clazzFun.stacks == null && !(clazzFun.lastClaxxRef != null
					&& clazzFun.lastClaxxRef.prototype[funName] != null
					&& clazzFun.lastClaxxRef.prototype[funName].stacks != null)) { // super.toString
			fx = clazzFun;
		} else { // normal wrapped method
			var stacks = clazzFun.stacks;
			if (stacks == null) {
				stacks = clazzFun.lastClaxxRef.prototype[funName].stacks;
			}
			var length = stacks.length;
			for (i = length - 1; i >= 0; i--) {
				/*
				 * Once super call is computed precisely, there are no need 
				 * to calculate the inherited level but just an equals
				 * comparision
				 */
				//var level = Clazz.getInheritedLevel (clazzThis, stacks[i]);
				if (clazzThis === stacks[i]) { // level == 0
					if (i > 0) {
						i--;
						fx = stacks[i].prototype[funName];
					} else {
						/*
						 * Will this case be reachable?
						 * March 4, 2006
						 * Should never reach here if all things are converted
						 * by Java2Script
						 */
						fx = stacks[0].prototype[funName]["\\unknown"];







				//if (funName == "clone")alert("fx=" + fx) 
 












					}
					break;
				} else if (Clazz.getInheritedLevel (clazzThis, 
						stacks[i]) > 0) {
					fx = stacks[i].prototype[funName];

					break;
				}
			} // end of for loop
		} // end of normal wrapped method
	} // end of clazzFun != null


	if (fx != null) {
		/* there are members which are initialized out of the constructor */
		if (i == 0 && funName == "construct") {
			var ss = clazzFun.stacks;
			if (ss != null && ss[0].superClazz == null
					&& ss[0].con$truct != null) {
				ss[0].con$truct.apply (objThis, []);
			}
		}
		/*# {$no.debug.support} >>x #*/
		if (Clazz.tracingCalling) {
			var caller = arguments.callee.caller;
			if (caller === Clazz.superConstructor) {
				caller = caller.arguments.callee.caller;
			}
			Clazz.pu$hCalling (new Clazz.callingStack (caller, clazzThis));
			var ret = fx.apply (objThis, (funParams == null) ? [] : funParams);
			Clazz.p0pCalling ();
			return ret;
		}
		/*# x<< #*/


		return fx.apply (objThis, (funParams == null) ? [] : funParams);
	} else if (funName == "construct") {
		/* there are members which are initialized out of the constructor */
		/* No super constructor! */
		return ;
	}
	Clazz.alert(["j2slib","no class found",(funParams).typeString])
	throw new Clazz.MethodNotFoundException (objThis, clazzThis, funName, 
			Clazz.getParamsType (funParams).typeString);
};

/**
 * Call super constructor of the class. 
 * The same effect as Java's expression: 
 * <code> super () </code>
 */
/* public */
Clazz.superConstructor = function (objThis, clazzThis, funParams) {
	Clazz.superCall (objThis, clazzThis, "construct", funParams);
	/* If there are members which are initialized out of the constructor */
	if (clazzThis.con$truct != null) {
		clazzThis.con$truct.apply (objThis, []);
	}
};

/**
 * Class for null with a given class as to be casted.
 * This class will be used as an implementation of Java's casting way.
 * For example,
 * <code> this.call ((String) null); </code>
 */
/* protcted */
Clazz.CastedNull = function (asClazz) {
	if (asClazz != null) {
		if (asClazz instanceof String) {
			this.clazzName = asClazz;
		} else if (asClazz instanceof Function) {
			this.clazzName = Clazz.getClassName (asClazz, true);
		} else {
			this.clazzName = "" + asClazz;
		}
	} else {
		this.clazzName = "Object";
	}
	this.toString = function () {
		return null;
	};
	this.valueOf = function () {
		return null;
	};
};

/**
 * API for Java's casting null.
 * @see Clazz.CastedNull
 *
 * @param asClazz given class
 * @return an instance of class Clazz.CastedNull
 */
/* public */
Clazz.castNullAs = function (asClazz) {
	return new Clazz.CastedNull (asClazz);
};

/** 
 * MethodException will be used as a signal to notify that the method is
 * not found in the current clazz hierarchy.
 */
/* private */
Clazz.MethodException = function () {
};
/* protected */
Clazz.MethodNotFoundException = function () {
	this.toString = function () {
		return "MethodNotFoundException";
	};
};

/* private */
/*x-# getParamsType -> gPT #-x*/
Clazz.getParamsType = function (funParams) {
	// bh: optimization here for very common cases
	var n = funParams.length;
	switch (n) {
	case 0:
		var params = ["void"];
		params.typeString = "\\void";
		return params;
	case 1:
		// just so common
		var obj = funParams[0];
		if (obj != null && typeof obj == "number") {
			var params = ["Number"];
			params.typeString = "\\Number";
			return params;
		}
	}

	var params = [];
	params.hasCastedNull = false;
	if (funParams != null) {
		for (var i = 0; i < n; i++) {
			params[i] = Clazz.getClassName (funParams[i]);
			if (funParams[i] instanceof Clazz.CastedNull) {
				params.hasCastedNull = true;
			}
		}
	}
	params.typeString = "\\" + params.join ('\\');
	return params;
};




/**
 * Search the given class prototype, find the method with the same
 * method name and the same parameter signatures by the given 
 * parameters, and then run the method with the given parameters.
 *
 * @param objThis the current host object
 * @param claxxRef the current host object's class
 * @param fxName the method name
 * @param funParams the given arguments
 * @return the result of the specified method of the host object,
 * the return maybe void.
 * @throws Clazz.MethodNotFoundException if no matched method is found
 */
/* protected */
/*-# searchAndExecuteMethod -> saem #-*/


Clazz.searchAndExecuteMethod = function (objThis, claxxRef, fxName, funParams) {
	var fx = objThis[fxName];
	var params = Clazz.getParamsType (funParams);
	Clazz._profile && Clazz.addProfile(claxxRef, fxName, params);

	/*
	 * Cache last matched method
	 */
	if (fx.lastParams == params.typeString && fx.lastClaxxRef === claxxRef) {
		var methodParams = null;
		if (params.hasCastedNull) {
			methodParams = new Array ();
			for (var k = 0; k < funParams.length; k++) {
				if (funParams[k] instanceof Clazz.CastedNull) {
					/*
					 * For Clazz.CastedNull instances, the type name is
					 * already used to indentified the method in Clazz#
					 * searchMethod.
					 */
					methodParams[k] = null;
				} else {
					methodParams[k] = funParams[k];
				}
			}
		} else {
			methodParams = funParams;
		}
		if (fx.lastMethod != null) {
			return fx.lastMethod.apply (objThis, methodParams);
		} else { // missed default constructor ?
			return ;
		}
	}
	fx.lastParams = params.typeString;
	fx.lastClaxxRef = claxxRef;

	var stacks = fx.stacks;
	if (stacks == null) {
		stacks = claxxRef.prototype[fxName].stacks;
	}
	var length = stacks.length;

	/*
	 * Search the inheritance stacks to get the given class' function
	 */
	var began = false; // began to search its super classes
	for (var i = length - 1; i > -1; i--) {
		//if (Clazz.getInheritedLevel (claxxRef, stacks[i]) >= 0) {
		/*
		 * No need to calculate the inherited level as there always exist a 
		 * right claxxRef in the stacks, and the inherited level of stacks
		 * are in order.
		 */
		if (began || stacks[i] === claxxRef) {
			/*
			 * First try to search method within the same class scope
			 * with stacks[i] === claxxRef
			 */
			var clazzFun = stacks[i].prototype[fxName];

			var ret = Clazz.tryToSearchAndExecute (fxName, objThis, clazzFun, params,
					funParams/*, isSuper, clazzThis*/, fx);
			if (!(ret instanceof Clazz.MethodException)) {
				return ret;
			}
			/*
			 * As there are no such methods in current class, Clazz will try 
			 * to search its super class stacks. Here variable began indicates
			 * that super searchi is began, and there is no need checking
			 * <code>stacks[i] === claxxRef</code>
			 */
			began = true; 
		} // end of if
	} // end of for
	if ("construct" == fxName) {
		/*
		 * For non existed constructors, just return without throwing
		 * exceptions. In Java codes, extending Object can call super
		 * default Object#constructor, which is not defined in JS.
		 */
		return ;
	}
	// TODO: should be java.lang.NoSuchMethodException
	throw new Clazz.MethodNotFoundException (objThis, claxxRef, 
			fxName, params.typeString);
};


/*# {$no.debug.support} >>x #*/
Clazz.tracingCalling = false;
/*# x<< #*/

/* private */
/*-# tryToSearchAndExecute -> tsae #-*/
Clazz.tryToSearchAndExecute = function (fxName, objThis, clazzFun, params, funParams/*, 
		isSuper, clazzThis*/, fx) {
		//if (fxName != "construct")System.out.println("tsae " + fxName + " " + funParams);
	var methods = new Array ();
	//var xfparams = null;
	var generic = true;
	for (var fn in clazzFun) {
		//if (fn.indexOf ('\\') == 0) {
		if (fn.charCodeAt (0) == 92) { // 92 == '\\'.charCodeAt (0)
			var ps = fn.substring (1).split ("\\");
			if (ps.length == params.length) {
				methods[methods.length] = ps;
			}
			generic = false;
			continue;
		}
		/*
		 * When there are only one method in the class, use the funParams
		 * to identify the parameter type.
		 *
		 * AbstractCollection.remove (Object)
		 * AbstractList.remove (int)
		 * ArrayList.remove (int)
		 *
		 * Then calling #remove (Object) method on ArrayList instance will 
		 * need to search up to the AbstractCollection.remove (Object),
		 * which contains only one method.
		 */
		/*
		 * See Clazz#defineMethod --Mar 10, 2006, josson
		 */
		if (generic && fn == "funParams" && clazzFun.funParams != null) {
			//xfparams = clazzFun.funParams;
			fn = clazzFun.funParams;
			var ps = fn.substring (1).split ("\\");
			if (ps.length == params.length) {
				methods[0] = ps;
			}
			break;
		}
	}
	if (methods.length == 0) {
		//throw new Clazz.MethodException ();
		return new Clazz.MethodException ();
	}
	var method = Clazz.searchMethod (methods, params);
	if (method != null) {
		var f = null;
		if (generic) { /* Use the generic method */
			/*
			 * Will this case be reachable?
			 * March 4, 2006 josson
			 * 
			 * Reachable for calling #remove (Object) method on 
			 * ArrayList instance
			 * May 5, 2006 josson
			 */
			f = clazzFun; // call it directly
		} else {
			f = clazzFun["\\" + method];
		}
		//if (f != null) { // always not null
			var methodParams = null;
			if (params.hasCastedNull) {
				methodParams = new Array ();
				for (var k = 0; k < funParams.length; k++) {
					if (funParams[k] instanceof Clazz.CastedNull) {
						/*
						 * For Clazz.CastedNull instances, the type name is
						 * already used to indentified the method in Clazz#
						 * searchMethod.
						 */
						methodParams[k] = null;
					} else {
						methodParams[k] = funParams[k];
					}
				}
			} else {
				methodParams = funParams;
			}
			/*# {$no.debug.support} >>x #*/
			if (Clazz.tracingCalling) {
				var caller = arguments.callee.caller; // SAEM
				caller = caller.arguments.callee.caller; // Delegating
				caller = caller.arguments.callee.caller; 
				var xpushed = f.exName == "construct" 
						&& Clazz.getInheritedLevel (f.exClazz, Throwable) >= 0
						&& !Clazz._initializingException;
				if (xpushed) {
					Clazz._initializingException = true;
					// constructor is wrapped
					var xcaller = caller.arguments.callee.caller // Delegate
							.arguments.callee.caller; // last method
					var fun = xcaller.arguments.callee;
					var owner = fun.claxxReference;
					if (owner == null) {
						owner = fun.exClazz;
					}
					if (owner == null) {
						owner = fun.claxxOwner;
					}
					/*
					 * Keep the environment that Throwable instance is created
					 */
					Clazz.pu$hCalling (new Clazz.callingStack (xcaller, owner));
				}

				var noInnerWrapper = caller !== Clazz.instantialize 
						&& caller !== Clazz.superCall;
				if (noInnerWrapper) {
					var fun = caller.arguments.callee;
					var owner = fun.claxxReference;
					if (owner == null) {
						owner = fun.exClazz;
					}
					if (owner == null) {
						owner = fun.claxxOwner;
					}
					Clazz.pu$hCalling (new Clazz.callingStack (caller, owner));
				}
				fx.lastMethod = f;
				var ret = f.apply (objThis, methodParams);
				if (noInnerWrapper) {
					Clazz.p0pCalling ();
				}
				if (xpushed) {
					Clazz.p0pCalling ();
				}
				return ret;
			}
			/*# x<< #*/
			fx.lastMethod = f;
			return f.apply (objThis, methodParams);
		//}
	}
	//throw new Clazz.MethodException ();
	return new Clazz.MethodException ();
};

/*# {$no.debug.support} >>x #*/
Clazz._initializingException = false;
/*# x<< #*/

/**
 * Search the existed polymorphic methods to get the matched method with
 * the given parameter types.
 *
 * @param existedMethods Array of string which contains method parameters
 * @param paramTypes Array of string that is parameter type.
 * @return string of method parameters seperated by "\\"
 */
/* private */
/*-# 
 # searchMethod -> sM 
 #
 # roundOne -> rO
 # paramTypes -> pts
 #-*/
Clazz.searchMethod = function (roundOne, paramTypes) {
	/*
	 * Filter out all the fitted methods for the given parameters
	 */
	/*-# roundTwo -> rT #-*/
	var roundTwo = new Array ();
	for (var i = 0; i < roundOne.length; i++) {
		/*-# fittedLevel -> fL #-*/
		var fittedLevel = new Array ();
		var isFitted = true;
		for (var j = 0; j < roundOne[i].length; j++) {
			fittedLevel[j] = Clazz.getInheritedLevel (paramTypes[j], 
					roundOne[i][j]);
			if (fittedLevel[j] < 0) {
				isFitted = false;
				break;
			}
		}
		if (isFitted) {
			fittedLevel[paramTypes.length] = i; // Keep index for later use
			roundTwo[roundTwo.length] = fittedLevel;
		}
	}
	if (roundTwo.length == 0) {
		return null;
	}
	/*
	 * Find out the best method according to the inheritance.
	 */
	/*-# resultTwo -> rtT #-*/
	var resultTwo = roundTwo;
	var min = resultTwo[0];
	for (var i = 1; i < resultTwo.length; i++) {
		/*-# isVectorLesser -> vl #-*/
		var isVectorLesser = true;
		for (var j = 0; j < paramTypes.length; j++) {
			if (min[j] < resultTwo[i][j]) {
				isVectorLesser = false;;
				break;
			}
		}
		if (isVectorLesser) {
			min = resultTwo[i];
		}
	}
	var index = min[paramTypes.length]; // Get the previously stored index
	/*
	 * Return the method parameters' type string as indentifier of the
	 * choosen method.
	 */
	return roundOne[index].join ('\\');
};

/**
 * Generate delegating function for the given method name.
 *
 * @param claxxRef the specified class for the method
 * @funName method name of the specified method
 * @return the method delegate which will try to search the method
 * from the given class by the parameters
 */
/* private */
/*-# generateDelegatingMethod -> gDM #-*/

Clazz.generateDelegatingMethod = function (claxxRef, funName, fCall) {
	/*
	 * Delegating method.
	 * Each time the following expression will generate a new 
	 * function object.
	 */
	var delegating = function () {
			var r = arguments;
			return SAEM (this, r.callee.claxxReference, r.callee.methodName, r);
	};
	delegating.methodName = funName;
	delegating.claxxReference = claxxRef;
	return delegating;
};



SAEM = Clazz.searchAndExecuteMethod;

/* private */
Clazz.expExpandParameters = function ($0, $1) {
	if ($1 == 'N') {
		return "Number";
	} else if ($1 == 'B') {
		return "Boolean"
	} else if ($1 == 'S') {
		return "String";
	} else if ($1 == 'O') {
		return "Object";
	} else if ($1 == 'A') {
		return "Array"
	}
	return "Unknown";
};

/*
 * Other developers may need to extend this formatParameters method
 * to deal complicated situation.
 */
/* protected */
Clazz.formatParameters = function (funParams) {
	if (funParams == null || funParams.length == 0) {
		return "\\void";
	}
	return funParams.replace (/~([NABSO])/g, Clazz.expExpandParameters)
				.replace (/\s+/g, "").replace (/^|,/g, "\\")
				.replace (/\$/g, "org.eclipse.s");

};

/*
 * Override the existed methods which are in the same name.
 * Overriding methods is provided for the purpose that the JavaScript
 * does not need to search the whole hierarchied methods to find the
 * correct method to execute.
 * Be cautious about this method. Incorrectly using this method may
 * break the inheritance system.
 *
 * @param clazzThis host class in which the method to be defined
 * @param funName method name
 * @param funBody function object, e.g function () { ... }
 * @param funParams paramether signature, e.g ["string", "number"]
 */
/* public */
Clazz.overrideMethod = function (clazzThis, funName, funBody, funParams) {
	if (Clazz.assureInnerClass) Clazz.assureInnerClass (clazzThis, funBody);
	funBody.exName = funName;
	var fpName = Clazz.formatParameters (funParams);
	/*
	 * Replace old methods with new method. No super methods are kept.
	 */
	funBody.funParams = fpName; 
	funBody.claxxOwner = clazzThis;
	return Clazz.getSignature(clazzThis.prototype, funName, funBody, true);
};

/*
 * Define method for the class with the given method name and method
 * body and method parameter signature.
 *
 * @param clazzThis host class in which the method to be defined
 * @param funName method name
 * @param funBody function object, e.g function () { ... }
 * @param funParams paramether signature, e.g ["string", "number"]
 * @return method of the given name. The method may be funBody or a wrapper
 * of the given funBody.
 */
/* public */
Clazz.defineMethod = function (clazzThis, funName, funBody, funParams) {
	if (Clazz.assureInnerClass) Clazz.assureInnerClass (clazzThis, funBody);
	funBody.exName = funName;
	var fpName = Clazz.formatParameters (funParams);

	/*
	 * For method the first time is defined, just keep it rather than
	 * wrapping into deep hierarchies!
	 */

	// BH : signature based on nParams
	var proto = clazzThis.prototype;
	var f$ = Clazz.getSignature(proto, funName, funBody, false);
	if (f$ == null || (f$.claxxOwner === clazzThis && f$.funParams == fpName)) {
		// property "funParams" will be used as a mark of only-one method
		funBody.funParams = fpName; 
		funBody.claxxOwner = clazzThis;
		funBody.exClazz = clazzThis; // make it traceable
		//delete $fz;           // BH -- delete global variables when no longer needed
		return Clazz.getSignature(proto, funName, funBody, true);
	}
	var oldFun = null;
	var oldStacks = new Array ();
		if (f$.stacks == null) {
			/* method is not defined by Clazz.defineMethod () */
			oldFun = f$;
			if (f$.claxxOwner != null) {
				oldStacks[0] = oldFun.claxxOwner;
			}
		} else {
			oldStacks = f$.stacks;
		}
		/*
	 * Method that is already defined in super class will be overridden
	 * with a new proxy method with class hierarchy stored in a stack.
	 * That is to say, the super methods are lost in this class' proxy
	 * method. 
	 * When method are being called, methods defined in the new proxy 
	 * method will be searched through first. And if no method fitted,
	 * it will then try to search method in the super class stacks.
	 */
	/* method has not been defined yet */
	/* method is not defined by Clazz.defineMethod () */
	/* method is defined in super class */
	if (f$.stacks == null 
			|| f$.claxxReference !== clazzThis) {
			//if (funName == "hashCode") {
				//alert(f$.claxxReference + " " + clazzThis)
			//}
		/*
		 * Generate a new delegating method for the class
		 */
		f$ = Clazz.getSignature(proto, funName, Clazz.generateDelegatingMethod (clazzThis, funName, f$), true);				
		//if (funName != "construct" && Clazz.__debuggingBH)
			// System.out.println("delegating " + clazzThis.__CLASS_NAME__ + " " + funName + " " + funParams);
		/*
		 * Keep the class inheritance stacks
		 */
		var arr = new Array ();
		for (var i = 0; i < oldStacks.length; i++) {
			arr[i] = oldStacks[i];
		}
		f$.stacks = arr;
	}
	var ss = f$.stacks;

	if (ss.length == 0) {
		ss[0] = clazzThis;
	} else {
		var existed = false;
		for (var i = ss.length - 1; i >= 0; i--) {
			if (ss[i] === clazzThis) {
				existed = true;
				break;
			}
		}
		if (!existed) {
			ss[ss.length] = clazzThis;
		}
	}

	if (oldFun != null) {
		if (oldFun.claxxOwner === clazzThis) {
			f$[oldFun.funParams] = oldFun;
			oldFun.claxxOwner = null;
			// property "funParams" will be used as a mark of only-one method
			oldFun.funParams = null; // null ? safe ? // safe for " != null"
		} else if (oldFun.claxxOwner == null) {
			/*
			 * The function is not defined Clazz.defineMethod ().
			 * Try to fixup the method ...
			 * As a matter of lost method information, I just suppose
			 * the method to be fixed is with void parameter!
			 */
			f$["\\unknown"] = oldFun;
		}
	}
	funBody.exClazz = clazzThis; // make it traceable
	f$[fpName] = funBody;
	//delete $fz;           // BH -- delete global variables when no longer needed
	return f$;
};                                                

/**
 * Make constructor for the class with the given function body and parameters
 * signature.
 * 
 * @param clazzThis host class
 * @param funBody constructor body
 * @param funParams constructor parameters signature
 */
/* public */
Clazz.makeConstructor = function (clazzThis, funBody, funParams) {
	Clazz.defineMethod (clazzThis, "construct", funBody, funParams);
	if (clazzThis.con$truct != null) {
		clazzThis.con$truct.index = clazzThis.con$truct.stacks.length;
	}
	//clazzThis.con$truct = clazzThis.prototype.con$truct = null;
};

/**
 * Override constructor for the class with the given function body and
 * parameters signature.
 * 
 * @param clazzThis host class
 * @param funBody constructor body
 * @param funParams constructor parameters signature
 */
/* public */
Clazz.overrideConstructor = function (clazzThis, funBody, funParams) {
// $_k  @j2sOverrideConstructor
	Clazz.overrideMethod (clazzThis, "construct", funBody, funParams);
	if (clazzThis.con$truct != null) {
		clazzThis.con$truct.index = clazzThis.con$truct.stacks.length;
	}
	//clazzThis.con$truct = clazzThis.prototype.con$truct = null;
};

/*
 * all root packages. e.g. java.*, org.*, com.*
 */
/* protected */
Clazz.allPackage = {};

/**
 * Will be used to keep value of whether the class is defined or not.
 */
/* protected */
Clazz.allClasses = {};

Clazz.lastPackageName = null;
Clazz.lastPackage = null;

/* protected */
Clazz.unloadedClasses = new Array ();

/* public */
Clazz.isClassUnloaded = function (clzz) {
	var thisClassName = Clazz.getClassName (clzz, true);
	return Clazz.unloadedClasses[thisClassName] != null;
};

/* public */
Clazz.declarePackage = function (pkgName) {
	if (Clazz.lastPackageName == pkgName) {
		return Clazz.lastPackage;
	}
	if (pkgName != null && pkgName.length != 0) {
		var pkgFrags = pkgName.split (/\./);
		var pkg = Clazz.allPackage;
		for (var i = 0; i < pkgFrags.length; i++) {
			if (pkg[pkgFrags[i]] == null) {
				pkg[pkgFrags[i]] = { 
					__PKG_NAME__ : ((pkg.__PKG_NAME__ != null) ? 
						pkg.__PKG_NAME__ + "." + pkgFrags[i] : pkgFrags[i])
				}; 
				// pkg[pkgFrags[i]] = {};
				if (i == 0) {
					// eval ...
					Clazz.setGlobal(pkgFrags[i], pkg[pkgFrags[i]]);
				}
			}
			pkg = pkg[pkgFrags[i]]
		}
		Clazz.lastPackageName = pkgName;
		Clazz.lastPackage = pkg;
		return pkg;
	}
};

/* protected */
/*x-# evalType -> eT  #-x*/
Clazz.evalType = function (typeStr, isQualified) {
	var idx = typeStr.lastIndexOf (".");
	if (idx != -1) {
		var pkgName = typeStr.substring (0, idx);
		var pkg = Clazz.declarePackage (pkgName);
		var clazzName = typeStr.substring (idx + 1);
		return pkg[clazzName];
	} else if (isQualified) {
		return window[typeStr];
	} else if (typeStr == "number") {
		return Number;
	} else if (typeStr == "object") {
		return Clazz._O;
	} else if (typeStr == "string") {
		return String;
	} else if (typeStr == "boolean") {
		return Boolean;
	} else if (typeStr == "function") {
		return Function;
	} else if (typeStr == "void" || typeStr == "undefined"
			|| typeStr == "unknown") {
		return typeStr;
	} else if (typeStr == "NullObject") {
		return NullObject;
	} else {
		return window[typeStr];
	}
};

/**
 * Define a class or interface.
 *
 * @param qClazzName String presents the qualified name of the class
 * @param clazzFun Function of the body
 * @param clazzParent Clazz to inherit from, may be null
 * @param interfacez Clazz may implement one or many interfaces
 *   interfacez can be Clazz object or Array of Clazz objects.
 * @return Ruturn the modified Clazz object
 */
/* public */
Clazz.defineType = function (qClazzName, clazzFun, clazzParent, interfacez) {
	var cf = Clazz.unloadedClasses[qClazzName];
	if (cf != null) {
		clazzFun = cf;
	}
	var idx = qClazzName.lastIndexOf (".");
	if (idx != -1) {
		var pkgName = qClazzName.substring (0, idx);
		var pkg = Clazz.declarePackage (pkgName);
		var clazzName = qClazzName.substring (idx + 1);
		if (pkg[clazzName] != null) {
			// already defined! Should throw exception!
			return pkg[clazzName];
		}
		pkg[clazzName] = clazzFun;
	} else {
		if (window[qClazzName] != null) {
			// already defined! Should throw exception!
			return window[qClazzName];
		}
		Clazz.setGlobal(qClazzName, clazzFun);
	}
	Clazz.decorateAsType (clazzFun, qClazzName, clazzParent, interfacez);
	/*# {$no.javascript.support} >>x #*/
	var iFun = Clazz._innerFunctions;
	clazzFun.defineMethod = iFun.defineMethod;
	clazzFun.defineStaticMethod = iFun.defineStaticMethod;
	clazzFun.makeConstructor = iFun.makeConstructor;
	/*# x<< #*/
	return clazzFun;
};

Clazz.isSafari = (navigator.userAgent.indexOf ("Safari") != -1);
Clazz.isSafari4Plus = false;
if (Clazz.isSafari) {
	var ua = navigator.userAgent;
	var verIdx = ua.indexOf ("Version/");
	if (verIdx  != -1) {
		var verStr = ua.substring (verIdx + 8);
		var verNumber = parseFloat (verStr);
		Clazz.isSafari4Plus = verNumber >= 4.0;
	}
}

/* protected */
Clazz.instantialize = function (objThis, args) {
	if (args != null && args.length == 1 && args[0] != null 
			&& args[0] instanceof Clazz.args4InheritClass) {
		return ;
	}
	/*
	if (objThis.con$truct != null) {
		objThis.con$truct.apply (objThis, args);
	}
	if (objThis.construct != null) {
		objThis.construct.apply (objThis, args);
	}
	*/
	if (objThis instanceof Number) {
		objThis.valueOf = function () {
			return this;
		};
	}
	if (Clazz.isSafari4Plus) { // Fix bug of Safari 4.0+'s over-optimization
		var argsClone = new Array ();
		for (var k = 0; k < args.length; k++) {
			argsClone[k] = args[k];
		}
		args = argsClone;
	}
	var c = objThis.construct;
	if (c != null) {
		if (objThis.con$truct == null) { // no need to init fields
			c.apply (objThis, args);
		} else if (objThis.getClass ().superClazz == null) { // the base class
			objThis.con$truct.apply (objThis, []);
			c.apply (objThis, args);
		} else if ((c.claxxOwner != null 
				&& c.claxxOwner === objThis.getClass ())
				|| (c.stacks != null 
				&& c.stacks[c.stacks.length - 1] == objThis.getClass ())) {
			/*
			 * This #construct is defined by this class itself.
			 * #construct will call Clazz.superConstructor, which will
			 * call #con$truct back
			 */
			c.apply (objThis, args);
		} else { // constructor is a super constructor
			if (c.claxxOwner != null && c.claxxOwner.superClazz == null 
						&& c.claxxOwner.con$truct != null) {
				c.claxxOwner.con$truct.apply (objThis, []);
			} else if (c.stacks != null && c.stacks.length == 1
					&& c.stacks[0].superClazz == null) {
				c.stacks[0].con$truct.apply (objThis, []);
			}
			c.apply (objThis, args);
			objThis.con$truct.apply (objThis, []);
		}
	} else if (objThis.con$truct != null) {
		objThis.con$truct.apply (objThis, []);
	}
};

/**
 * Once there are other methods registered to the Function.prototype, 
 * those method names should be add to the following Array.
 */
/*
 * static final member of interface may be a class, which may
 * be function.
 */
/* protected */
/*-# innerFunctionNames -> iFN #-*/
Clazz.innerFunctionNames = [
	"equals", "hashCode", /*"toString",*/ "getName", "getClassLoader", "getResource", "getResourceAsStream" /*# {$no.javascript.support} >>x #*/, "defineMethod", "defineStaticMethod",
	"makeConstructor" /*# x<< #*/
];

/*
 * Static methods
 */
/*x-# innerFunctions -> inF #-x*/
Clazz._innerFunctions = {
	/*
	 * Similar to Object#equals
	 */
	equals : function (aFun) {
		return this === aFun;
	},

	hashCode : function () {
		return this.getName ().hashCode ();
	},

	toString : function () {
		return "class " + this.getName ();
	},

	/*
	 * Similar to Class#getName
	 */
	getName : function () {
		return Clazz.getClassName (this, true);
	},
	getClassLoader : function () {
		var clazzName = this.__CLASS_NAME__;
		var baseFolder = Clazz._Loader.getClasspathFor (clazzName);
		var x = baseFolder.lastIndexOf (clazzName.replace (/\./g, "/"));
		if (x != -1) {
			baseFolder = baseFolder.substring (0, x);
		} else {
			baseFolder = Clazz._Loader.getClasspathFor (clazzName, true);
		}
		var loader = Clazz._Loader.requireLoaderByBase (baseFolder);
		loader.getResourceAsStream = Clazz._innerFunctions.getResourceAsStream;
		loader.getResource = Clazz._innerFunctions.getResource; // BH
		return loader;
	},

	getResource : function(name) {
		return this.getResourceAsStream(name).url;
	},

	getResourceAsStream : function (name) {
		var is = null;
		if (name == null) {
			return is;
		}
		if (java.io.InputStream != null) {
			is = new java.io.InputStream ();
		} else {
			is = new Clazz._O ();
			is.__CLASS_NAME__ = "java.io.InputStream";
			is.close = NullObject; // function () {};
		}
		is.read = function () { return 0; };
		name = name.replace (/\\/g, '/');
		/*-# baseFolder -> bFr #-*/
		var baseFolder = null;
		var clazzName = this.__CLASS_NAME__;
		if (arguments.length == 2 && name.indexOf ('/') != 0) { // additional argument
			name = "/" + name;
		}
		if (name.indexOf ('/') == 0) {
			//is.url = name.substring (1);
			if (arguments.length == 2) { // additional argument
				baseFolder = arguments[1];
				if (baseFolder == null) {
					baseFolder = Clazz._Loader.binaryFolders[0];
				}
			} else if (Clazz._Loader != null) {
				baseFolder = Clazz._Loader.getClasspathFor (clazzName, true);
			}
			if (baseFolder == null || baseFolder.length == 0) {
				is.url = name.substring (1);
			} else {
				baseFolder = baseFolder.replace (/\\/g, '/');
				var length = baseFolder.length;
				var lastChar = baseFolder.charAt (length - 1);
				if (lastChar != '/') {
					baseFolder += "/";
				}
				is.url = baseFolder + name.substring (1);
			}
		} else {
			if (this.base != null) {
				baseFolder = this.base;
			} else if (Clazz._Loader != null) {
				baseFolder = Clazz._Loader.getClasspathFor (clazzName);
				var x = baseFolder.lastIndexOf (clazzName.replace (/\./g, "/"));
				if (x != -1) {
					baseFolder = baseFolder.substring (0, x);
				} else {
					//baseFolder = null;
					var y = -1;
					if (baseFolder.indexOf (".z.js") == baseFolder.length - 5
							&& (y = baseFolder.lastIndexOf ("/")) != -1) {
						baseFolder = baseFolder.substring (0, y + 1);
						var pkgs = clazzName.split (/\./);
						for (var k = 1; k < pkgs.length; k++) {
							var pkgURL = "/";
							for (var j = 0; j < k; j++) {
								pkgURL += pkgs[j] + "/";
							}
							if (pkgURL.length > baseFolder.length) {
								break;
							}
							if (baseFolder.indexOf (pkgURL) == baseFolder.length - pkgURL.length) {
								baseFolder = baseFolder.substring (0, baseFolder.length - pkgURL.length + 1);
								break;
							}
						}
					} else {
						baseFolder = Clazz._Loader.getClasspathFor (clazzName, true);
					}
				}
			} else {
				var bins = Clazz.binaryFolders;
				if (bins != null && bins.length != 0) {
					baseFolder = bins[0];
				}
			}
			if (baseFolder == null || baseFolder.length == 0) {
				baseFolder = "j2s/";
			}
			baseFolder = baseFolder.replace (/\\/g, '/');
			var length = baseFolder.length;
			var lastChar = baseFolder.charAt (length - 1);
			if (lastChar != '/') {
				baseFolder += "/";
			}
			/*
			 * FIXME: bug here for "/"
			 */
			//if (baseFolder.indexOf ('/') == 0) {
			//	baseFolder = baseFolder.substring (1);
			//}
			if (this.base != null) {
				is.url = baseFolder + name;
			} else {
				var idx = clazzName.lastIndexOf ('.');
				if (idx == -1 || this.base != null) {
					is.url = baseFolder + name;
				} else {
					is.url = baseFolder + clazzName.substring (0, idx)
							.replace (/\./g, '/') +  "/" + name;
				}
			}
		}
		return is;
	}/*# {$no.javascript.support} >>x #*/,

	/*
	 * For JavaScript programmers
	 */
	defineMethod : function (methodName, funBody, paramTypes) {
		Clazz.defineMethod (this, methodName, funBody, paramTypes);
	},

	/*
	 * For JavaScript programmers
	 */
	defineStaticMethod : function (methodName, funBody, paramTypes) {
		Clazz.defineMethod (this, methodName, funBody, paramTypes);
		this[methodName] = this.prototype[methodName];
	},

	/*
	 * For JavaScript programmers
	 */
	makeConstructor : function (funBody, paramTypes) {
		Clazz.makeConstructor (this, funBody, paramTypes);
	}
	/*# x<< #*/
};

/* private */
/*-# decorateFunction -> dF #-*/
Clazz.decorateFunction = function (clazzFun, prefix, name) {
	if (Clazz._Loader != null) {
		//alert ("decorate " + name);
		Clazz._Loader.checkInteractive ();
	}
	var qName = null;
	if (prefix == null) {
		// e.g. Clazz.declareInterface (null, "ICorePlugin", 
		//		org.eclipse.ui.IPlugin);
		qName = name;
		Clazz.setGlobal(name, clazzFun);
	} else if (prefix.__PKG_NAME__ != null) {
		// e.g. Clazz.declareInterface (org.eclipse.ui, "ICorePlugin", 
		//		org.eclipse.ui.IPlugin);
		qName = prefix.__PKG_NAME__ + "." + name;
		prefix[name] = clazzFun;
		if (prefix === java.lang) {
			Clazz.setGlobal(name, clazzFun);
		}
	} else {
		// e.g. Clazz.declareInterface (org.eclipse.ui.Plugin, "ICorePlugin", 
		//		org.eclipse.ui.IPlugin);
		qName = prefix.__CLASS_NAME__ + "." + name;
		prefix[name] = clazzFun;
		//alert("j2slib.z Clazz.decorateFunction qname=" + qName)
	}
	Clazz.extendJO(clazzFun, qName);
	var inF = Clazz.innerFunctionNames;
	for (var i = 0; i < inF.length; i++) {
		clazzFun[inF[i]] = Clazz._innerFunctions[inF[i]];
	}

	if (Clazz._Loader != null) {
		/*-# findClass -> fC #-*/
		var node = Clazz._Loader.findClass (qName);
		/*-#
		 # ClazzNode.STATUS_KNOWN -> 1
		 #-*/
		if (node != null && node.status == Clazz._Node.STATUS_KNOWN) {
			/*-# 
			 # updateNode -> uN 
			 #-*/
			window.setTimeout((function(nnn) {
				return function() {
					Clazz._Loader.updateNode (nnn);
				};
			})(node), 1);
			/*
			 * #updateNode should be delayed! Or the class itself won't
			 * be initialized completely before marking itself as loaded.
			 */
			// ClazzLoader.updateNode (node);
		}
	}
};
Clazz.currentPath= "";
/* proected */
Clazz.declareInterface = function (prefix, name, interfacez) {
	var clazzFun = function () {};
	Clazz.decorateFunction (clazzFun, prefix, name, "Clazz.declareInterface");
	if (interfacez != null) {
		Clazz.implementOf (clazzFun, interfacez);
	}
	return clazzFun;
};

/* protected */
/*-# 
 # parentClazzInstance -> pi
 # clazzParent -> pc
 #-*/
Clazz.decorateAsClass = function (clazzFun, prefix, name, clazzParent, 
		interfacez, parentClazzInstance, fromWhere) {
	var prefixName = null;
	if (prefix != null) {
		prefixName = prefix.__PKG_NAME__;
		if (prefixName == null) {
			prefixName = prefix.__CLASS_NAME__;
		}
	}
	var qName = (prefixName == null ? "" : prefixName + ".") + name;
	var cf = Clazz.unloadedClasses[qName];
	if (cf != null) {
		clazzFun = cf;
	}
	var qName = null;
	Clazz.decorateFunction (clazzFun, prefix, name, fromWhere + "...Clazz.decorateAsClass");
	if (parentClazzInstance != null) {
		Clazz.inheritClass (clazzFun, clazzParent, parentClazzInstance);
	} else if (clazzParent != null) {
		Clazz.inheritClass (clazzFun, clazzParent);
	}
	if (interfacez != null) {
		Clazz.implementOf (clazzFun, interfacez);
	}
	return clazzFun;
};

/* public */
Clazz.declareType = function (prefix, name, clazzParent, interfacez, 
		parentClazzInstance) {
	var f = function () {
		Clazz.instantialize (this, arguments);
	};
	return Clazz.decorateAsClass (f, prefix, name, clazzParent, interfacez, 
			parentClazzInstance, "Clazz.declareType");
};

/* public */
Clazz.declareAnonymous = function (prefix, name, clazzParent, interfacez, 
		parentClazzInstance) {
	var f = function () {
		Clazz.prepareCallback (this, arguments);
		Clazz.instantialize (this, arguments);
	};
	return Clazz.decorateAsClass (f, prefix, name, clazzParent, interfacez, 
			parentClazzInstance, "Clazz.declareAnonymous");
};

/* protected */
Clazz.decorateAsType = function (clazzFun, qClazzName, clazzParent, 
		interfacez, parentClazzInstance, inheritClazzFuns) {
	Clazz.extendJO(clazzFun, qClazzName);
	clazzFun.equals = Clazz._innerFunctions.equals;
	clazzFun.getName = Clazz._innerFunctions.getName;
	if (inheritClazzFuns) {
		for (var i = 0; i < Clazz.innerFunctionNames.length; i++) {
			var methodName = Clazz.innerFunctionNames[i];
			clazzFun[methodName] = Clazz._innerFunctions[methodName];
		}
	}
	if (parentClazzInstance != null) {
		Clazz.inheritClass (clazzFun, clazzParent, parentClazzInstance);
	} else if (clazzParent != null) {
		Clazz.inheritClass (clazzFun, clazzParent);
	}
	if (interfacez != null) {
		Clazz.implementOf (clazzFun, interfacez);
	}
	return clazzFun;
};

/* sgurin: native exception detection mechanism. Only NullPointerException detected and wrapped to java excepions */
/** private utility method for creating a general regexp that can be used later  
 * for detecting a certain kind of native exceptions. use with error messages like "blabla IDENTIFIER blabla"
 * @param msg String - the error message
 * @param spliterName String, must be contained once in msg
 * spliterRegex String, a string with the regexp literal for identifying the spitter in exception further error messages.
 */
Clazz._ex_reg=function (msg, spliterName, spliterRegex) {
	if(!spliterRegex) 
		spliterRegex="[^\\s]+";	
	var idx = msg.indexOf (spliterName), 
		str = msg.substring (0, idx) + spliterRegex + msg.substring(idx + spliterName.length), 
		regexp = new RegExp("^"+str+"$");
	return regexp;
};
// reproduce NullPointerException for knowing how to detect them, and create detector function Clazz._isNPEExceptionPredicate
var $$o$$ = null;
try {
	$$o$$.hello ();
} catch (e) {
	if(/Opera[\/\s](\d+\.\d+)/.test(navigator.userAgent)) {// opera throws an exception with fixed messages like "Statement on line 23: Cannot convert undefined or null to Object Backtrace: Line....long text... " 
		var idx1 = e.message.indexOf(":"), idx2 = e.message.indexOf(":", idx1+2);
		Clazz._NPEMsgFragment = e.message.substr(idx1+1, idx2-idx1-20);
		Clazz._isNPEExceptionPredicate = function(e) {
			return e.message.indexOf(Clazz._NPEMsgFragment)!=-1;
		};
	}	
	else if(navigator.userAgent.toLowerCase().indexOf("webkit")!=-1) { //webkit, google chrome prints the property name accessed. 
		Clazz._exceptionNPERegExp = Clazz._ex_reg(e.message, "hello");
		Clazz._isNPEExceptionPredicate = function(e) {
			return Clazz._exceptionNPERegExp.test(e.message);
		};
	}
	else {// ie, firefox and others print the name of the object accessed: 
		Clazz._exceptionNPERegExp = Clazz._ex_reg(e.message, "$$o$$");
		Clazz._isNPEExceptionPredicate = function(e) {
			return Clazz._exceptionNPERegExp.test(e.message);
		};
	}		
};
/**sgurin
 * Implements Java's keyword "instanceof" in JavaScript's way **for exception objects**.
 * 
 * calls Clazz.instanceOf if e is a Java exception. If not, try to detect known native 
 * exceptions, like native NullPointerExceptions and wrap it into a Java exception and 
 * call Clazz.instanceOf again. if the native exception can't be wrapped, false is returned.
 * 
 * @param obj the object to be tested
 * @param clazz the class to be checked
 * @return whether the object is an instance of the class
 * @author: sgurin
 */
Clazz.exceptionOf=function(e, clazz) {
	if(e.__CLASS_NAME__) {
		return Clazz.instanceOf(e, clazz);
	}
	if(clazz == Error) {
		if (("" + e).indexOf("Error") >= 0) {
			System.out.println(Clazz.getStackTrace());
		}
		return (("" + e).indexOf("Error") >= 0);
		// everything here is a Java Exception, not a Java Error
	}
	return (clazz == Exception || clazz == Throwable
		|| clazz == NullPointerException && Clazz._isNPEExceptionPredicate(e));
};

Clazz.getStackTrace = function(n) {
	var s = "\n";
	n || (n = 25);
	var c = arguments.callee.caller;
		for (var i = 0; i < n; i++) {
			if (!c)break;
			s += (i + " " + (c.exName ? (c.claxxOwner ? c.claxxOwner.__CLASS_NAME__ + "."  : "") + c.exName 
			: (c.toString ? c.toString().substring(0, c.toString().indexOf("{")) : "<native method>"))) + "\n";
			c = c.caller
		}
		return s;
}

/* sgurin: preserve Number.prototype.toString */
Number.prototype._numberToString=Number.prototype.toString;


Clazz.declarePackage ("java.io");
//Clazz.declarePackage ("java.lang");
Clazz.declarePackage ("java.lang.annotation"); // java.lang
Clazz.declarePackage ("java.lang.instrument"); // java.lang
Clazz.declarePackage ("java.lang.management"); // java.lang
Clazz.declarePackage ("java.lang.reflect"); // java.lang
Clazz.declarePackage ("java.lang.ref");  // java.lang.ref
java.lang.ref.reflect = java.lang.reflect;
Clazz.declarePackage ("java.util");

/*
 * Consider these interfaces are basic!
 */
Clazz.declareInterface (java.io,"Closeable");
Clazz.declareInterface (java.io,"DataInput");
Clazz.declareInterface (java.io,"DataOutput");
Clazz.declareInterface (java.io,"Externalizable");
Clazz.declareInterface (java.io,"Flushable");
Clazz.declareInterface (java.io,"Serializable");
Clazz.declareInterface (java.lang,"Iterable");
Clazz.declareInterface (java.lang,"CharSequence");
Clazz.declareInterface (java.lang,"Cloneable");
Clazz.declareInterface (java.lang,"Appendable");
Clazz.declareInterface (java.lang,"Comparable");
Clazz.declareInterface (java.lang,"Runnable");
Clazz.declareInterface (java.util,"Comparator");

java.lang.ClassLoader = {
	__CLASS_NAME__ : "ClassLoader"
};

/******************************************************************************
 * Copyright (c) 2007 java2script.org and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Zhou Renjian - initial API and implementation
 *****************************************************************************/
/*******
 * @author zhou renjian
 * @create March 10, 2006
 *******/

if (window["Clazz"] != null && window["Clazz"].unloadClass == null) {
/**
 * Once ClassExt.js is part of Class.js.
 * In order to make the Class.js as small as possible, part of its content
 * is moved into this ClassExt.js.
 *
 * See also http://j2s.sourceforge.net/j2sclazz/
 */
 
/**
 * Clazz.MethodNotFoundException is used to notify the developer about calling
 * methods with incorrect parameters.
 */
/* protected */
// Override the Clazz.MethodNotFoundException in Class.js to give details
Clazz.MethodNotFoundException = function (obj, clazz, method, params) {
	var paramStr = "";
	if (params != null) {
		paramStr = params.substring (1).replace (/\\/g, ",");
	}
	var leadingStr = "";
	if (method != null && method != "construct") {
		leadingStr = "Method";
	} else {
		leadingStr = "Constructor";
	}
	this.message = leadingStr + " " + Clazz.getClassName (clazz, true) + "." 
					+ method + "(" + paramStr + ") is not found!";
	this.toString = function () {
		return "MethodNotFoundException:" + this.message;
	}
};

/**
 * Prepare callback for instance of anonymous Class.
 * For example for the callback:
 *     this.callbacks.MyEditor.sayHello();
 *
 * @param objThis the host object for callback
 * @param args arguments object. args[0] will be classThisObj -- the "this"
 * object to be hooked
 * 
 * Attention: parameters should not be null!
 */
/* protected */
Clazz.prepareCallback = function (objThis, args) {
	var classThisObj = args[0];
	var cbName = "b$"; // "callbacks";
	if (objThis != null && classThisObj != null && classThisObj !== window) {
		var obs = new Array ();
		if (objThis[cbName] == null) {
			objThis[cbName] = obs;
		} else { // must make a copy!
			for (var s in objThis[cbName]) {
				if (s != "length") {
					obs[s] = objThis[cbName][s];
				}
			}
			objThis[cbName] = obs;
		}
		var className = Clazz.getClassName (classThisObj, true);
		//if (obs[className] == null) { /* == null make no sense! */
			//obs[className] = classThisObj;
			/*
			 * TODO: the following line is SWT-specific! Try to move it out!
			 */
			obs[className.replace (/org\.eclipse\.swt\./, "$wt.")] = classThisObj;
			var clazz = Clazz.getClass (classThisObj);
			while (clazz.superClazz != null) {
				clazz = clazz.superClazz;
				//obs[Clazz.getClassName (clazz)] = classThisObj;
				/*
				 * TODO: the following line is SWT-specific! Try to move it out!
				 */
				obs[Clazz.getClassName (clazz, true)
						.replace (/org\.eclipse\.swt\./, "$wt.")] = classThisObj;
			}
		//}
		var cbs = classThisObj[cbName];
		if (cbs != null && cbs instanceof Array) {
			for (var s in cbs) {
				if (s != "length") {
					obs[s] = cbs[s];
				}
			}
		}
	}
	// Shift the arguments
	for (var i = 0; i < args.length - 1; i++) {
		args[i] = args[i + 1];
	}
	args.length--;
	// arguments will be returned!
};

/**
 * Construct instance of the given inner class.
 *
 * @param classInner given inner class, alway with name like "*$*"
 * @param objThis this instance which can be used to call back.
 * @param finalVars final variables which the inner class may use
 * @return the constructed object
 *
 * @see Clazz#cloneFinals
 */
/* public */
Clazz.innerTypeInstance = function (clazzInner, objThis, finalVars) {
	if (clazzInner == null) {
		clazzInner = arguments.callee.caller;
	}
	var obj = null;
	if (finalVars == null && objThis.$finals == null) {
		/*if (arguments.length == 2) {
			obj = new clazzInner (objThis);
		} else */if (arguments.length == 3) {
			obj = new clazzInner (objThis);
		} else if (arguments.length == 4) {
			if (objThis.__CLASS_NAME__ == clazzInner.__CLASS_NAME__
					&& arguments[3] === Clazz.inheritArgs) {
				obj = objThis;
			} else {
				obj = new clazzInner (objThis, arguments[3]);
			}
		} else if (arguments.length == 5) {
			obj = new clazzInner (objThis, arguments[3], arguments[4]);
		} else if (arguments.length == 6) {
			obj = new clazzInner (objThis, arguments[3], arguments[4], 
					arguments[5]);
		} else if (arguments.length == 7) {
			obj = new clazzInner (objThis, arguments[3], arguments[4], 
					arguments[5], arguments[6]);
		} else if (arguments.length == 8) {
			obj = new clazzInner (objThis, arguments[3], arguments[4], 
					arguments[5], arguments[6], arguments[7]);
		} else if (arguments.length == 9) {
			obj = new clazzInner (objThis, arguments[3], arguments[4], 
					arguments[5], arguments[6], arguments[7], arguments[8]);
		} else if (arguments.length == 10) {
			obj = new clazzInner (objThis, arguments[3], arguments[4], 
					arguments[5], arguments[6], arguments[7], arguments[8],
					arguments[9]);
		} else {
			/*
			 * Should construct instance manually.
			 */
			obj = new clazzInner (objThis, Clazz.inheritArgs);
			//if (obj.construct == null) {
			//	throw new String ("No support anonymous class constructor with " 
			//			+ "more than 7 parameters.");
			//}
			var args = new Array ();
			for (var i = 3; i < arguments.length; i++) {
				args[i - 3] = arguments[i];
			}
			//obj.construct.apply (obj, args);
			Clazz.instantialize (obj, args);
		}
	} else {
		obj = new clazzInner (objThis, Clazz.inheritArgs);
		// f$ is short for the once choosen "$finals"
		if (finalVars != null && objThis.f$ == null) {
			obj.f$ = finalVars;
		} else if (finalVars == null && objThis.f$ != null) {
			obj.f$ = objThis.f$;
		} else if (finalVars != null && objThis.f$ != null) {
			var o = new Object ();
			for (var attr in objThis.f$) {
				o[attr] = objThis.f$[attr];
			}
			for (var attr in finalVars) {
				o[attr] = finalVars[attr];
			}
			obj.f$ = o;
		}

		var args = new Array ();
		for (var i = 3; i < arguments.length; i++) {
			args[i - 3] = arguments[i];
		}
		Clazz.instantialize (obj, args);
	}

	/*
	if (finalVars != null && objThis.$finals == null) {
		obj.$finals = finalVars;
	} else if (finalVars == null && objThis.$finals != null) {
		obj.$finals = objThis.$finals;
	} else if (finalVars != null && objThis.$finals != null) {
		var o = {};
		for (var attr in objThis.$finals) {
			o[attr] = objThis.$finals[attr];
		}
		for (var attr in finalVars) {
			o[attr] = finalVars[attr];
		}
		obj.$finals = o;
	}
	*/
	//Clazz.prepareCallback (obj, objThis);
	return obj;
};

/**
 * Clone variables whose modifier is "final".
 * Usage: var o = Clazz.cloneFinals ("name", name, "age", age);
 *
 * @return Object with all final variables
 */
/* protected */
Clazz.cloneFinals = function () {
	var o = {};
	var length = arguments.length / 2;
	for (var i = 0; i < length; i++) {
		o[arguments[i + i]] = arguments[i + i + 1];
	}
	return o;
};

/* public */
Clazz.isClassDefined = Clazz.isDefinedClass = function (clazzName) {
	if (clazzName != null && clazzName.length != 0) {
		if (Clazz.allClasses[clazzName]) {
			return true;
		}
		var pkgFrags = clazzName.split (/\./);
		var pkg = null;
		for (var i = 0; i < pkgFrags.length; i++) {
			if (pkg == null) {
				if (Clazz.allPackage[pkgFrags[0]] == null) {
					//error (clazzName + " / " + false);
					return false;
				}
				pkg = Clazz.allPackage[pkgFrags[0]];
			} else {
				if (pkg[pkgFrags[i]] == null) {
					//error (clazzName + " / " + false);
					return false;
				}
				pkg = pkg[pkgFrags[i]]
			}
		}
		//error (clazzName + " / " + (pkg != null));
		//return pkg != null;
		if (pkg != null) {
			Clazz.allClasses[clazzName] = true;
			return true;
		} else {
			return false;
		}
	} else {
		/* consider null or empty name as non-defined class */
		return false;
	}
};
/**
 * Define the enum constant.
 * @param classEnum enum type
 * @param enumName enum constant
 * @param enumOrdinal enum ordinal
 * @param initialParams enum constant constructor parameters
 * @return return defined enum constant
 */
/* public */
Clazz.defineEnumConstant = function (clazzEnum, enumName, enumOrdinal, initialParams, clazzEnumExt) {
	var o = null;
	if (clazzEnumExt != null) {
		o = new clazzEnumExt ();
	} else {
		o = new clazzEnum ();
	}
	// BH avoids unnecessary calls to SAEM
	o.$name = enumName;
	o.$ordinal = enumOrdinal;
	//Clazz.superConstructor (o, clazzEnum, [enumName, enumOrdinal]);
	if (initialParams != null && initialParams.length != 0) {
		o.construct.apply (o, initialParams);
	}
	clazzEnum[enumName] = o;
	clazzEnum.prototype[enumName] = o;
	if (clazzEnum["$ values"] == null) {  // BH added
		clazzEnum["$ values"] = []          // BH added
		clazzEnum.values = function() {     // BH added
			return this["$ values"];          // BH added
		};                                  // BH added
	}
	clazzEnum["$ values"].push(o);
	return o;
};

//////// (int) conversions //////////

Clazz.floatToInt = function (x) {
	return x < 0 ? Math.ceil(x) : Math.floor(x);
};

Clazz.floatToByte = Clazz.floatToShort = Clazz.floatToLong = Clazz.floatToInt;
Clazz.doubleToByte = Clazz.doubleToShort = Clazz.doubleToLong = Clazz.doubleToInt = Clazz.floatToInt;

Clazz.floatToChar = function (x) {
	return String.fromCharCode (x < 0 ? Math.ceil(x) : Math.floor(x));
};

Clazz.doubleToChar = Clazz.floatToChar;



//////// Array additions ////////////

if (self.Int32Array && self.Int32Array != Array) {
	Clazz.haveInt32 = true;
	if (!Int32Array.prototype.sort)
		Int32Array.prototype.sort = Array.prototype.sort
} else {
	Int32Array = function(n) {
		if (!n) n = 0;
		var b = new Array(n);
		b.toString = function(){return "[object Int32Array]"}
		for (var i = 0; i < n; i++)b[i] = 0
		return b;
	}
	Clazz.haveInt32 = false;
	Int32Array.prototype.sort = Array.prototype.sort
	Int32Array.prototype.int32Fake = function(){};
}

if (self.Float64Array && self.Float64Array != Array) {
	Clazz.haveFloat64 = true;
	if (!Float64Array.prototype.sort)
		Float64Array.prototype.sort = Array.prototype.sort
} else {
	Clazz.haveFloat64 = false;
	Float64Array = function(n) {
		if (!n) n = 0;
		var b = new Array(n);
		for (var i = 0; i < n; i++)b[i] = 0.0
		return b;
	};
	Float64Array.prototype.sort = Array.prototype.sort
	Float64Array.prototype.float64Fake = function() {}; // "present"
	Float64Array.prototype.toString = function() {return "[object Float64Array]"};
// Darn! Mozilla makes this a double, not a float. It's 64-bit.
// and Safari 5.1 doesn't have Float64Array 
}

/**
 * Make arrays.
 *
 * @return the created Array object
 */
/* public */
Clazz.newArray  = function () {
	if (arguments[0] instanceof Array) {
		// recursive, from newArray(n,m,value)
		// as newArray([m, value], newInt32Array)
		var args = arguments[0];
		var f = arguments[1];
	} else {
		var args = arguments;
		var f = Array;
	}
	if (args.length <= 1) return new Array(); // maybe never?
	var dim = args[0];
	if (typeof dim == "string") {
		dim = dim.charCodeAt (0); // char
	}
	var len = args.length - 1;
	var val = args[len];
	if (args.length == 2) {
		if (val == null)
			return new Array(dim);
		if (f === true && Clazz.haveInt32) return new Int32Array(dim);
		if (f === false && Clazz.haveFloat64) return new Float64Array(dim);
		if (f == Array && val == null) return new Array(dim);
		var arr = (f === true ? new Int32Array() : f === false ? new Float64Array() : new Array(dim));
		for (var i = dim; --i >= 0;)
			arr[i] = val;
		return arr;
	}
	var xargs = new Array (len);
	for (var i = 0; i < len; i++) {
		xargs[i] = args[i + 1];
	}
	var arr = new Array (dim);
	if (val == null || val >= 0 || len > 2)
		for (var i = 0; i < dim; i++) {
	 	// Call recursively!
			arr[i] = Clazz.newArray (xargs, f);
		}
	return arr;
};

Clazz.newArray32 = function(args, isInt32) {
	var dim = args[0];
	if (typeof dim == "string") {
		dim = dim.charCodeAt (0); // char
	}
	var len = args.length - 1;
	var val = args[len];
	switch (args.length) {
	case 0:
	case 1:
		alert("ERROR IN newArray32 -- args length < 2");
		return new Array(0);
	case 2:
		if (val < 0)
			return new Array(dim);
	//  if (isInt32 ? !Clazz.haveInt32 : !Clazz.haveFloat64 ) {
			// no support for Int32Array in MSIE or Float64Array in Safari 5.1
			// so we must initialize ourselves
	//    var f = (isInt32 ? new Int32Array() : new Float64Array());
	//    for (var i = dim; --i >= 0;)
	//      f[i] = 0;
	//    return f;
	//  }
		try {
			return (isInt32 ? new Int32Array(dim) : new Float64Array(dim));
		}catch (e) {
			alert(dim + " " + arguments.callee.caller.arguments.callee.caller + e)
		}
	}
	var xargs = new Array(len);
	for (var i = len; --i >= 0;) {
		xargs[i] = args[i + 1];
	}
	var arr = new Array (dim);
	for (var i = 0; i < dim; i++) {
		// Call newArray referencing this array type
		// only for the final iteration, and only if val === 0
		arr[i] = Clazz.newArray (xargs, isInt32);
	}
	return arr;
};


/**
 * Make arrays.
 *
 * @return the created Array object
 */
/* public */
Clazz.newInt32Array  = function () {
	return Clazz.newArray32(arguments, true);
}

/**
 * Make arrays.
 *
 * @return the created Array object
 */
/* public */
Clazz.newFloat64Array  = function () {
	return Clazz.newArray32(arguments, false);
}

Clazz.newFloatArray = Clazz.newDoubleArray = Clazz.newFloat64Array;
Clazz.newIntArray = Clazz.newLongArray = Clazz.newShortArray = Clazz.newByteArray = Clazz.newInt32Array;
Clazz.newCharArray = Clazz.newBooleanArray = Clazz.newArray;

//$_AI=Clazz.newIntArray;
//$_AF=Clazz.newFloatArray;
//$_AD=Clazz.newDoubleArray;
//$_AL=Clazz.newLongArray;
//$_AS=Clazz.newShortArray;
//$_AB=Clazz.newByteArray;
//$_AC=Clazz.newCharArray;
//$_Ab=Clazz.newBooleanArray;


Clazz.arrayIs = function(a, what) {
	// for some reason, Number.constructor.toString() now gives "too much recursion"
	return a.constructor && a.constructor != Number && a.constructor.toString().indexOf(what) >= 0
}

Clazz.isAS = function(a) { // just checking first parameter
	return (a && typeof a == "object" && Clazz.arrayIs(a, " Array") && (typeof a[0] == "string" || typeof a[0] == "undefined"));
}

Clazz.isASS = function(a) {
	return (a && typeof a == "object" && Clazz.isAS(a[0]));
}

Clazz.isAP = function(a) {
	return (a && Clazz.getClassName(a[0]) == "JU.P3");
}

Clazz.isAI = function(a) {
	return (a && typeof a == "object" && (Clazz.haveInt32 ? Clazz.arrayIs(a, "Int32Array") : a.int32Fake ? true : false));
}

Clazz.isAII = function(a) { // assumes non-null a[0]
	return (a && typeof a == "object" && Clazz.isAI(a[0]));
}

Clazz.isAF = function(a) {
	return (a && typeof a == "object" && (Clazz.haveFloat64 ? Clazz.arrayIs(a, "Float64Array") : a.float64Fake ? true : false));
}

Clazz.isAFF = function(a) { // assumes non-null a[0]
	return (a && typeof a == "object" && Clazz.isAF(a[0]));
}

Clazz.isAFFF = function(a) { // assumes non-null a[0]
	return (a && typeof a == "object" && Clazz.isAFF(a[0]));
}

Clazz.isAFloat = function(a) { // just checking first parameter
	return (a && typeof a == "object" && Clazz.arrayIs(a, " Array") && Clazz.instanceOf(a[0], Float));
}


/**
 * Make the RunnableCompatiability instance as a JavaScript function.
 *
 * @param jsr Instance of RunnableCompatiability
 * @return JavaScript function instance represents the method run of jsr.
 */
/* public */
Clazz.makeFunction = function (jsr) {
	return function (e) {
		if (e == null) {
			e = window.event;
		}
		if (jsr.setEvent != null) {
			jsr.setEvent (e);
		}
		jsr.run ();
		/*
		if (e != null && jsr.isReturned != null && jsr.isReturned()) {
			// Is it correct to stopPropagation here? --Feb 19, 2006
			e.cancelBubble = true;
			if (e.stopPropagation) {
				e.stopPropagation();
			}
		}
		*/
		if (jsr.returnSet == 1) {
			return jsr.returnNumber;
		} else if (jsr.returnSet == 2) {
			return jsr.returnBoolean;
		} else if (jsr.returnSet == 3) {
			return jsr.returnObject;
		}
	};
};

/* protected */
Clazz.defineStatics = function (clazz) {
	for (var j = arguments.length, i = (j - 1) / 2; --i >= 0;) {
		var val = arguments[--j]
		var name = arguments[--j];
		clazz[name] = clazz.prototype[name] = val;
	}
};

/* protected */
Clazz.prepareFields = function (clazz, fieldsFun) {
	var stacks = new Array ();
 //System.out.println(fieldsFun)
	if (clazz.con$truct != null) {
		var ss = clazz.con$truct.stacks;
		var idx = 0;//clazz.con$truct.index;
		//System.out.println("clazz index = " + idx + " sslen=" + ss.length)
		for (var i = idx; i < ss.length; i++) {
			stacks[i] = ss[i];
		}
	}
	//System.out.println(JSON.stringify(stacks))
	Clazz.addProto(clazz.prototype, "con$truct", clazz.con$truct = function () {
		var stacks = arguments.callee.stacks;
		if (stacks != null) {
			for (var i = 0; i < stacks.length; i++) {
				stacks[i].apply (this, []);
			}
		}
	});
	stacks[stacks.length] = fieldsFun;
	clazz.con$truct.stacks = stacks;
	clazz.con$truct.index = 0;
};

/*
 * Serialize those public or protected fields in class 
 * net.sf.j2s.ajax.SimpleSerializable.
 */
/* protected */
Clazz.registerSerializableFields = function (clazz) {
	var args = arguments;
	var length = args.length;
	var newArr = new Array ();
	if (clazz.declared$Fields != null) {
		for (var i = 0; i < clazz.declared$Fields.length; i++) {
			newArr[i] = clazz.declared$Fields[i];
		}
	}
	clazz.declared$Fields = newArr;

	if (length > 0 && length % 2 == 1) {
		var fs = clazz.declared$Fields;
		for (var i = 1; i <= (length - 1) / 2; i++) {
			var o = { name : args[i + i - 1], type : args[i + i] };
			var existed = false;
			for (var j = 0; j < fs.length; j++) {
				if (fs[j].name == o.name) { // reloaded classes
					fs[j].type = o.type; // update type
					existed = true;
					break;
				}
			}
			if (!existed) {
				fs[fs.length] = o;
			}
		}
	}
};

/*
 * Get the caller method for those methods that are wrapped by 
 * Clazz.searchAndExecuteMethod.
 *
 * @param args caller method's arguments
 * @return caller method, null if there is not wrapped by 
 * Clazz.searchAndExecuteMethod or is called directly.
 */
/* protected */
/*-# getMixedCallerMethod -> gMCM  #-*/
Clazz.getMixedCallerMethod = function (args) {
	var o = {};
	var argc = args.callee.caller; // Clazz.tryToSearchAndExecute
	if (argc == null) return null;
	if (argc !== Clazz.tryToSearchAndExecute) { // inherited method's apply
		argc = argc.arguments.callee.caller;
		if (argc == null) return null;
	}
	if (argc !== Clazz.tryToSearchAndExecute) return null;
	argc = argc.arguments.callee.caller; // Clazz.searchAndExecuteMethod
	if (argc == null || argc !== Clazz.searchAndExecuteMethod) return null;
	o.claxxRef = argc.arguments[1];
	o.fxName = argc.arguments[2];
	o.paramTypes = Clazz.getParamsType (argc.arguments[3]);
	argc = argc.arguments.callee.caller; // Clazz.generateDelegatingMethod
	if (argc == null) return null;
	argc = argc.arguments.callee.caller; // the private method's caller
	if (argc == null) return null;
	o.caller = argc;
	return o;
};

/* BH -- The issue here is a subclass calling its private method FOO when
 *       there is also a private method of the same name in its super class.
 *       This can ALWAYS be avoided and, one could argue, is bad 
 *       program design anyway. In Jmol, the presence of this possibility
 *       creates over 8000 references to the global $fx, which was only
 *       checked in a few rare cases. We can then also remove $fz references.
 *         
 */

/*
 * Check and return super private method.
 * In order make private methods be executed correctly, some extra javascript
 * must be inserted into the beggining of the method body of the non-private 
 * methods that with the same method signature as following:
 * <code>
 *			var $private = Clazz.checkPrivateMethod (arguments);
 *			if ($private != null) {
 *				return $private.apply (this, arguments);
 *			}
 * </code>
 * Be cautious about this. The above codes should be insert by Java2Script
 * compiler or with double checks to make sure things work correctly.
 *
 * @param args caller method's arguments
 * @return private method if there are private method fitted for the current 
 * calling environment
 */
/* public */

/*** BH
Clazz.checkPrivateMethod = function (args) {
	var m = Clazz.getMixedCallerMethod (args);
	if (m == null) return null;
	var callerFx = m.claxxRef.prototype[m.caller.exName];
	if (callerFx == null) return null; // may not be in the class hierarchies
	var ppFun = null;
	if (callerFx.claxxOwner != null) {
		ppFun = callerFx.claxxOwner.prototype[m.fxName];
	} else {
		var stacks = callerFx.stacks;
		for (var i = stacks.length - 1; i >= 0; i--) {
			var fx = stacks[i].prototype[m.caller.exName];
			if (fx === m.caller) {
				ppFun = stacks[i].prototype[m.fxName];
			} else if (fx != null) {
				for (var fn in fx) {
					if (fn.indexOf ('\\') == 0 && fx[fn] === m.caller) {
						ppFun = stacks[i].prototype[m.fxName];
						break;
					}
				}
			}
			if (ppFun != null) {
				break;
			}
		}
	}
	if (ppFun != null && ppFun.claxxOwner == null) {
		ppFun = ppFun["\\" + m.paramTypes];
	}
	if (ppFun != null && ppFun.isPrivate && ppFun !== args.callee) {
		return ppFun;
	}
	return null;
};
*****/
//$fz = null; // for private method declaration
c$ = null;
/*-# cla$$$tack -> cst  #-*/
Clazz.cStack = new Array ();
Clazz.pu$h = function (c) {
	c && Clazz.cStack.push(c);
};

Clazz.p0p = function () {
	return Clazz.cStack.pop();
};

/*# {$no.debug.support} >>x #*/
/*
 * Option to switch on/off of stack traces.
 */
/* protect */
Clazz.tracingCalling = false;

/*
 * Use to mark that the Throwable instance is created or not.
 */
/* private */
Clazz._initializingException = false;

/* private */
Clazz.callingStack = function (caller, owner) {
	this.caller = caller;
	this.owner = owner;
};
Clazz._callingStackTraces = new Array ();
Clazz.pu$hCalling = function (stack) {
	Clazz._callingStackTraces[Clazz._callingStackTraces.length] = stack;
};
Clazz.p0pCalling = function () {
	var length = Clazz._callingStackTraces.length;
	if (length > 0) {
		var stack = Clazz._callingStackTraces[length - 1];
		Clazz._callingStackTraces.length--;
		return stack;
	} else {
		return null;
	}
};
/*# x<< #*/

/**
 * The first folder is considered as the primary folder.
 * And try to be compatiable with ClazzLoader system.
 */
/* private */


/*** not used in Jmol
 * *
if (window["ClazzLoader"] != null && ClazzLoader.binaryFolders != null) {
	Clazz.binaryFolders = ClazzLoader.binaryFolders;
} else {
	Clazz.binaryFolders = ["j2s/", "", "j2slib/"];
}

Clazz.addBinaryFolder = function (bin) {
	if (bin != null) {
		var bins = Clazz.binaryFolders;
		for (var i = 0; i < bins.length; i++) {
			if (bins[i] == bin) {
				return ;
			}
		}
		bins[bins.length] = bin;
	}
};
Clazz.removeBinaryFolder = function (bin) {
	if (bin != null) {
		var bins = Clazz.binaryFolders;
		for (var i = 0; i < bins.length; i++) {
			if (bins[i] == bin) {
				for (var j = i; j < bins.length - 1; j++) {
					bins[j] = bins[j + 1];
				}
				bins.length--;
				return bin;
			}
		}
	}
	return null;
};
Clazz.setPrimaryFolder = function (bin) {
	if (bin != null) {
		Clazz.removeBinaryFolder (bin);
		var bins = Clazz.binaryFolders;
		for (var i = bins.length - 1; i >= 0; i--) {
			bins[i + 1] = bins[i];
		}
		bins[0] = bin;
	}
};

***/



/**
 * This is a simple implementation for Clazz#load. It just ignore dependencies
 * of the class. This will be fine for jar *.z.js file.
 * It will be overriden by ClazzLoader#load.
 * For more details, see ClazzLoader.js
 */
/* protected */
Clazz.load = function (musts, clazz, optionals, declaration) {
	if (declaration != null) {
		declaration ();
	}
	delete c$;           // BH -- delete global variables when no longer needed?
};

/*
 * Invade the Object prototype!
 * TODO: make sure that invading Object prototype does not affect other
 * existed library, such as Dojo, YUI, Prototype, ...
 */
java.lang.Object = Clazz._O;

Clazz._O.getName = Clazz._innerFunctions.getName;


System = {
	props : null, //new java.util.Properties (),
	$props : {},
	arraycopy : function (src, srcPos, dest, destPos, length) {
		if (src !== dest) {
			for (var i = 0; i < length; i++) {
				dest[destPos + i] = src[srcPos + i];
			}
		} else {
			var swap = [];
			for (var i = 0; i < length; i++) {
				swap[i] = src[srcPos + i];
			}
			for (var i = 0; i < length; i++) {
				dest[destPos + i] = swap[i];
			}
		}
	},
	currentTimeMillis : function () {
		return new Date ().getTime ();
	},
	gc : function() {}, // bh
	getProperties : function () {
		return System.props;
	},
	getProperty : function (key, def) {
		if (System.props)
			return System.props.getProperty (key, def);
		var v = System.$props[key];
		return (v != null ? v : arguments.length == 1 ? null : def != null ? def : key); // BH
	},
	getSecurityManager : function() { return null },  // bh
	setProperties : function (props) {
		System.props = props;
	},
	setProperty : function (key, val) {
		if (System.props == null)
			return System.$props[key] = val; // BH
		System.props.setProperty (key, val);
	}
};

System.setProperty("line.separator", navigator.userAgent.indexOf("Windows")>=0?"\r\n" : "\n") //BH
System.setProperty("os.name", navigator.userAgent) //BH

System.out = new Clazz._O ();
System.out.__CLASS_NAME__ = "java.io.PrintStream";
System.out.print = function () {};
System.out.printf = function () {};
System.out.println = function () {};

System.err = new Clazz._O ();
System.err.__CLASS_NAME__ = "java.io.PrintStream";
System.err.print = function () {};
System.err.printf = function () {};
System.err.println = function () {};

Clazz.popup = Clazz.assert = Clazz.log = Clazz.error = window.alert;

Thread = function () {};
Thread.J2S_THREAD = Thread.prototype.J2S_THREAD = new Thread ();
Thread.currentThread = Thread.prototype.currentThread = function () {
	return this.J2S_THREAD;
};

/* public */
Clazz.intCast = function (n) { // 32bit
	var b1 = (n & 0xff000000) >> 24;
	var b2 = (n & 0xff0000) >> 16;
	var b3 = (n & 0xff00) >> 8;
	var b4 = n & 0xff;
	if ((b1 & 0x80) != 0) {
		return -(((b1 & 0x7f) << 24) + (b2 << 16) + (b3 << 8) + b4 + 1);
	} else {
		return (b1 << 24) + (b2 << 16) + (b3 << 8) + b4;
	}
};

/* public */
Clazz.shortCast = function (s) { // 16bit
	var b1 = (n & 0xff00) >> 8;
	var b2 = n & 0xff;
	if ((b1 & 0x80) != 0) {
		return -(((b1 & 0x7f) << 8) + b2 + 1);
	} else {
		return (b1 << 8) + b4;
	}
};

/* public */
Clazz.byteCast = function (b) { // 8bit
	if ((b & 0x80) != 0) {
		return -((b & 0x7f) + 1);
	} else {
		return b & 0xff;
	}
};

/* public */
Clazz.charCast = function (c) { // 8bit
	return String.fromCharCode (c & 0xff).charAt (0);
};

/**
 * Warning: Unsafe conversion!
 */
/* public */
Clazz.floatCast = function (f) { // 32bit
	return f;
};

/*
 * Try to fix JavaScript's shift operator defects on long type numbers.
 */

Clazz.longMasks = [];

Clazz.longReverseMasks = [];

Clazz.longBits = [];

(function () {
	var arr = [1];
	for (var i = 1; i < 53; i++) {
		arr[i] = arr[i - 1] + arr[i - 1]; // * 2 or << 1
	}
	Clazz.longBits = arr;
	Clazz.longMasks[52] = arr[52];
	for (var i = 51; i >= 0; i--) {
		Clazz.longMasks[i] = Clazz.longMasks[i + 1] + arr[i];
	}
	Clazz.longReverseMasks[0] = arr[0];
	for (var i = 1; i < 52; i++) {
		Clazz.longReverseMasks[i] = Clazz.longReverseMasks[i - 1] + arr[i];
	}
}) ();


/* public */
Clazz.longLeftShift = function (l, o) { // 64bit
	if (o == 0) return l;
	if (o >= 64) return 0;
	if (o > 52) {
		error ("[Java2Script] Error : JavaScript does not support long shift!");
		return l;
	}
	if ((l & Clazz.longMasks[o - 1]) != 0) {
		error ("[Java2Script] Error : Such shift operator results in wrong calculation!");
		return l;
	}
	var high = l & Clazz.longMasks[52 - 32 + o];
	if (high != 0) {
		return high * Clazz.longBits[o] + (l & Clazz.longReverseMasks[32 - o]) << 0;
	} else {
		return l << o;
	}
};

/* public */
Clazz.intLeftShift = function (n, o) { // 32bit
	return (n << o) & 0xffffffff;
};

/* public */
Clazz.longRightShift = function (l, o) { // 64bit
	if ((l & Clazz.longMasks[52 - 32]) != 0) {
		return Math.round((l & Clazz.longMasks[52 - 32]) / Clazz.longBits[32 - o]) + (l & Clazz.longReverseMasks[o]) >> o;
	} else {
		return l >> o;
	}
};

/* public */
Clazz.intRightShift = function (n, o) { // 32bit
	return n >> o; // no needs for this shifting wrapper
};

/* public */
Clazz.long0RightShift = function (l, o) { // 64bit
	return l >>> o;
};

/* public */
Clazz.int0RightShift = function (n, o) { // 64bit
	return n >>> o; // no needs for this shifting wrapper
};

// Compress the common public API method in shorter name
//$_L=Clazz.load;
//$_W=Clazz.declareAnonymous;$_T=Clazz.declareType;
//$_J=Clazz.declarePackage;$_C=Clazz.decorateAsClass;
//$_Z=Clazz.instantialize;$_I=Clazz.declareInterface;$_D=Clazz.isClassDefined;
//$_H=Clazz.pu$h;$_P=Clazz.p0p;$_B=Clazz.prepareCallback;
//$_N=Clazz.innerTypeInstance;$_K=Clazz.makeConstructor;$_U=Clazz.superCall;$_R=Clazz.superConstructor;
//$_M=Clazz.defineMethod;$_V=Clazz.overrideMethod;$_S=Clazz.defineStatics;
//$_E=Clazz.defineEnumConstant;
//$_F=Clazz.cloneFinals;
//$_Y=Clazz.prepareFields;$_A=Clazz.newArray;$_O=Clazz.instanceOf;
//$_G=Clazz.inheritArgs;$_X=Clazz.checkPrivateMethod;$_Q=Clazz.makeFunction;
//$_s=Clazz.registerSerializableFields;
//$_k=Clazz.overrideConstructor;


var reflect = Clazz.declarePackage ("java.lang.reflect");
Clazz.declarePackage ("java.security");

Clazz.innerFunctionNames = Clazz.innerFunctionNames.concat (["getSuperclass",
		"isAssignableFrom", "getMethods", "getMethod", "getDeclaredMethods", 
		"getDeclaredMethod", "getConstructor", "getModifiers", "isArray", "newInstance"]);

Clazz._innerFunctions.getSuperclass = function () {
	return this.superClazz;	
};
Clazz._innerFunctions.isAssignableFrom = function (clazz) {
	return Clazz.getInheritedLevel (clazz, this) >= 0;	
};
Clazz._innerFunctions.getConstructor = function () {
	return new java.lang.reflect.Constructor (this, [], [], 
			java.lang.reflect.Modifier.PUBLIC);
};
/**
 * TODO: fix bug for polymorphic methods!
 */
Clazz._innerFunctions.getDeclaredMethods = Clazz._innerFunctions.getMethods = function () {
	var ms = new Array ();
	var p = this.prototype;
	for (var attr in p) {
		if (typeof p[attr] == "function" && p[attr].__CLASS_NAME__ == null) {
			/* there are polynormical methods. */
			ms[ms.length] = new java.lang.reflect.Method (this, attr,
					[], java.lang.Void, [], java.lang.reflect.Modifier.PUBLIC);
		}
	}
	p = this;
	for (var attr in p) {
		if (typeof p[attr] == "function" && p[attr].__CLASS_NAME__ == null) {
			ms[ms.length] = new java.lang.reflect.Method (this, attr,
					[], java.lang.Void, [], java.lang.reflect.Modifier.PUBLIC
					| java.lang.reflect.Modifier.STATIC);
		}
	}
	return ms;
};
Clazz._innerFunctions.getDeclaredMethod = Clazz._innerFunctions.getMethod = function (name, clazzes) {
	var p = this.prototype;
	for (var attr in p) {
		if (name == attr && typeof p[attr] == "function" 
				&& p[attr].__CLASS_NAME__ == null) {
			/* there are polynormical methods. */
			return new java.lang.reflect.Method (this, attr,
					[], java.lang.Void, [], java.lang.reflect.Modifier.PUBLIC);
		}
	}
	p = this;
	for (var attr in p) {
		if (name == attr && typeof p[attr] == "function" 
				&& p[attr].__CLASS_NAME__ == null) {
			return new java.lang.reflect.Method (this, attr,
					[], java.lang.Void, [], java.lang.reflect.Modifier.PUBLIC
					| java.lang.reflect.Modifier.STATIC);
		}
	}
	return null;
};
Clazz._innerFunctions.getModifiers = function () {
	return java.lang.reflect.Modifier.PUBLIC;
};
Clazz._innerFunctions.isArray = function () {
	return false;
};
Clazz._innerFunctions.newInstance = function () {
	var clz = this;
	return new clz ();
};

//Object.newInstance = Clazz._innerFunctions.newInstance;
;(function(){  // BH added wrapper here
	var inF = Clazz.innerFunctionNames;
	for (var i = 0; i < inF.length; i++) {
		Clazz._O[inF[i]] = Clazz._innerFunctions[inF[i]];
		Array[inF[i]] = Clazz._innerFunctions[inF[i]];
	}
	Array["isArray"] = function () {
		return true;
	};
})();

/* public */
Clazz._4Name = function (clazzName) {
	if (Clazz.isClassDefined (clazzName)) {
		return Clazz.evalType (clazzName);
	}
	if (Clazz._Loader != null) {
		Clazz._Loader.setLoadingMode ("xhr.sync");
		Clazz._Loader.loadClass (clazzName);
		//alert("TESTING HERE in Clazz.forName")
		return Clazz.evalType (clazzName);
	} else {
		Clazz.alert ("[Java2Script] Error: No ClassLoader!");
	}
};

/* For hotspot and unloading */

/* private */
Clazz.cleanDelegateMethod = function (m) {
	if (m == null) return;
	if (typeof m == "function" && m.lastMethod != null
			&& m.lastParams != null && m.lastClaxxRef != null) {
		m.lastMethod = null;
		m.lastParams = null;
		m.lastClaxxRef = null;
	}
};

/* public */
Clazz.unloadClass = function (qClazzName) {
	var cc = Clazz.evalType (qClazzName);
	if (cc != null) {
		Clazz.unloadedClasses[qClazzName] = cc;
		var clazzName = qClazzName;
		var pkgFrags = clazzName.split (/\./);
		var pkg = null;
		for (var i = 0; i < pkgFrags.length - 1; i++) {
			if (pkg == null) {
				pkg = Clazz.allPackage[pkgFrags[0]];
			} else {
				pkg = pkg[pkgFrags[i]]
			}
		}
		if (pkg == null) {
			Clazz.allPackage[pkgFrags[0]] = null;
			window[pkgFrags[0]] = null;
			// also try to unload inner or anonymous classes
			for (var c in window) {
				if (c.indexOf (qClazzName + "$") == 0) {
					Clazz.unloadClass (c);
					window[c] = null;
				}
			}
		} else {
			pkg[pkgFrags[pkgFrags.length - 1]] = null;
			// also try to unload inner or anonymous classes
			for (var c in pkg) {
				if (c.indexOf (pkgFrags[pkgFrags.length - 1] + "$") == 0) {
					Clazz.unloadClass (pkg.__PKG_NAME__ + "." + c);
					pkg[c] = null;
				}
			}
		}

		if (Clazz.allClasses[qClazzName] == true) {
			Clazz.allClasses[qClazzName] = false;
			// also try to unload inner or anonymous classes
			for (var c in Clazz.allClasses) {
				if (c.indexOf (qClazzName + "$") == 0) {
					Clazz.allClasses[c] = false;
				}
			}
		}

		for (var m in cc) {
			Clazz.cleanDelegateMethod (cc[m]);
		}
		for (var m in cc.prototype) {
			Clazz.cleanDelegateMethod (cc.prototype[m]);
		}

		if (Clazz._Loader != null) {
			Clazz._Loader.unloadClassExt (qClazzName);
		}

		return true;
	}
	return false;
};

}
/******************************************************************************
 * Copyright (c) 2007 java2script.org and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Zhou Renjian - initial API and implementation
 *****************************************************************************/
/*******
 * @author zhou renjian
 * @create July 10, 2006
 *******/

//if (window["ClazzNode"] == null) {
/**
 * TODO:
 * Make optimization over class dependency tree.
 */

/*
 * ClassLoader Summary
 * 
 * ClassLoader creates SCRIPT elements and setup class path and onload 
 * callback to continue class loading.
 *
 * In the onload callbacks, ClazzLoader will try to calculate the next-to-be-
 * load *.js and load it. In *.js, it will contains some codes like
 * Clazz.load (..., "$wt.widgets.Control", ...);
 * to provide information to build up the class dependency tree.
 *
 * Some known problems of different browsers:
 * 1. In IE, loading *.js through SCRIPT will first triggers onreadstatechange 
 * event, and then executes inner *.js source.
 * 2. In Firefox, loading *.js will first executes *.js source and then 
 * triggers onload event.
 * 3. In Opera, similar to IE, but trigger onload event. (TODO: More details 
 * should be studied. Currently, Opera supports no multiple-thread-loading)
 * 
 * For class dependency tree, actually, it is not a tree. It is a reference
 * net with nodes have n parents and n children. There is a root, which 
 * ClassLoader knows where to start searching and loading classes, for such
 * a net. Each node is a class. Each class may require a set of must-classes, 
 * which must be loaded before itself getting initialized, and also need a set
 * of optional classes, which also be loaded before being called.
 *
 * The class loading status will be in 6 stages.
 * 1. Unknown, the class is newly introduced by other class.
 * 2. Known, the class is already mentioned by other class.
 * 3. Loaded, *.js source is in memory, but may not be initialized yet. It 
 * requires all its must-classes be intiailized, which is in the next stage.
 * 4. Musts loaded, all must classes is already loaded and declared.
 * 5. Delcared, the class is already declared (ClazzLoader#isClassDefined).
 * 6. Optionals loaded, all optional classes is loaded and declared.
 *
 * The ClassLoader tries to load all necessary classes in order, and intialize
 * them in order. For such job, it will traverse the dependency tree, and try 
 * to next class to-be-loaded. Sometime, the class dependencies may be in one
 * or more cycles, which must be broken down so classes is loaded in correct
 * order.
 *
 * Loading order and intializing order is very important for the ClassLoader.
 * The following technical options are considered:
 * 1. SCRIPT is loading asynchronously, which means controling order must use
 * callback methods to continue.
 * 2. Multiple loading threads are later introduced, which requires the 
 * ClassLoader should use variables to record the class status.
 * 3. Different browsers have different loading orders, which means extra tests
 * should be tested to make sure loading order won't be broken.
 * 4. Java2Script simulator itself have some loading orders that must be 
 * honored, which means it should be integrated seamlessly to Clazz system.
 * 5. Packed *.z.js is introduced to avoid lots of small *.js which requires 
 * lots of HTTP connections, which means that packed *.z.js should be treated
 * specially (There will be mappings for such packed classes).
 * 6. *.js or *.css loading may fail according to network status, which means
 * another loading try should be performed, so ClazzLoader is more robust.
 * 7. SWT lazy loading is later introduced, which means that class loading
 * process may be paused and should be resumed later.
 *
 * Some known bugs:
 * <code>$_L(["$wt.graphics.Drawable","$wt.widgets.Widget"],
 *  "$wt.widgets.Control", ...</code>
 * has errors while must classes in different order such as
 * <code>$_L(["$wt.widgets.Widget", "$wt.graphics.Drawable"],
 *  "$wt.widgets.Control", ...</code>
 * has no error.
 * 
 * Other maybe bug scenarios:
 * 1. In <code>ClazzLoader.maxLoadingThreads = 1;</code> single loading thread 
 * mode, there are no errors, but in default multiple thread loading mode, 
 * there are errors.
 * 2. No errors in one browser, but has errors on other browsers (Browser 
 * script loading order differences).
 * 3. First time loading has errors, but reloading it gets no errors (Maybe 
 * HTTP connections timeout, but should not accur in local file system, or it
 * is a loading bug by using JavaScript timeout thread).
 */

/*
 * The following comments with "#" are special configurations for a much
 * smaller *.js file size.
 *
 * @see net.sf.j2s.lib/src/net/sf/j2s/lib/build/SmartJSCompressor.java
 */
/*-#
 # ClazzNode -> $CN$
 # ClazzLoader -> $CL$
 # <<< ClazzLoader = $CL$;
 #-*/
 
/*-#
 # parents -> sp
 # musts -> sm
 # xxxoptionals -> so
 # declaration -> dcl
 # optionalsLoaded -> oled
 # qClazzName ->Nq
 #-*/
 
/**
 * Static class loader class
 */
Clazz._Loader = function () {};


/**
 * Class dependency tree node
 */
/* private */
Clazz._Node = function () {
	Clazz._Loader.initNode(this);
};

;(function(Clazz, ClazzLoader, ClazzNode) {
/*-#
 # ClazzNode.STATUS_UNKNOWN = 0
 # ClazzNode.STATUS_KNOWN -> 1
 # ClazzNode.STATUS_CONTENT_LOADED -> 2
 # ClazzNode.STATUS_MUSTS_LOADED -> 3
 # ClazzNode.STATUS_DECLARED -> 4
 # ClazzNode.STATUS_OPTIONALS_LOADED -> 5
 #-*/
/*# >>x #*/

ClazzLoader.initNode = function(node) {
	node.parents = new Array ();
	node.musts = new Array ();
	node.optionals = new Array ();
	node.declaration = null;
	node.name = null; // id
	node.path = null;
	node.status = 0;
	node.random = 0.13412;
	node.optionalsLoaded = null;
}

ClazzNode.prototype.toString = function () {
	return this.name || this.path || "ClazzNode";
}

ClazzNode.STATUS_UNKNOWN = 0;
ClazzNode.STATUS_KNOWN = 1;
ClazzNode.STATUS_CONTENT_LOADED = 2;
ClazzNode.STATUS_MUSTS_LOADED = 3;
ClazzNode.STATUS_DECLARED = 4;
ClazzNode.STATUS_OPTIONALS_LOADED = 5;
/*# x<< #*/

						 
ClazzLoader.loaders = [];

ClazzLoader.requireLoaderByBase = function (base) {
	for (var i = 0; i < ClazzLoader.loaders.length; i++) {
		if (ClazzLoader.loaders[i].base == base) {
			return ClazzLoader.loaders[i];
		}
	}
	var loader = new ClazzLoader ();
	loader.base = base; 
	ClazzLoader.loaders[ClazzLoader.loaders.length] = loader;
	return loader;
};

/**
 * Class dependency tree
 */
/*-# clazzTreeRoot -> tr #-*/
ClazzLoader.clazzTreeRoot = new ClazzNode ();

/**
 * Used to keep the status whether a given *.js path is loaded or not.
 */
/* private */
/*-# loadedScripts -> ls #-*/
ClazzLoader.loadedScripts = {};

/**
 * Multiple threads are used to speed up *.js loading.
 */
/* private */
/*-# inLoadingThreads -> ilt #-*/
ClazzLoader.inLoadingThreads = 0;

/**
 * Maximum of loading threads
 */
/* protected */
ClazzLoader.maxLoadingThreads = 6;

ClazzLoader.userAgent = navigator.userAgent.toLowerCase ();
ClazzLoader.isOpera = (ClazzLoader.userAgent.indexOf ("opera") != -1);
ClazzLoader.isIE = (ClazzLoader.userAgent.indexOf ("msie") != -1) && !ClazzLoader.isOpera;
ClazzLoader.isGecko = (ClazzLoader.userAgent.indexOf ("gecko") != -1);
//ClazzLoader.isChrome = (ClazzLoader.userAgent.indexOf ("chrome") != -1);

/*
 * Opera has different loading order which will result in performance degrade!
 * So just return to single thread loading in Opera!
 *
 * FIXME: This different loading order also causes bugs in single thread!
 */
if (ClazzLoader.isOpera) {
	ClazzLoader.maxLoadingThreads = 1;
	var index = ClazzLoader.userAgent.indexOf ("opera/");
	if (index != -1) {
		var verNumber = 9.0;
		try {
			verNumber = parseFloat(ClazzLoader.userAgent.subString (index + 6));
		} catch (e) {}
		if (verNumber >= 9.6) {
			ClazzLoader.maxLoadingThreads = 6;
		}
	} 
}

/**
 * Try to be compatiable with Clazz system.
 * In original design ClazzLoader and Clazz are independent!
 *  -- zhourenjian @ December 23, 2006
 */
if (window["Clazz"] != null && Clazz.isClassDefined) {
	ClazzLoader.isClassDefined = Clazz.isClassDefined;
} else {
	/*-# definedClasses -> dC #-*/
	ClazzLoader.definedClasses = {};
	ClazzLoader.isClassDefined = function (clazzName) {
		return ClazzLoader.definedClasses[clazzName] == true;
	};
}

/*
 * binaryFolders will be used for ResourceBundle to check *.properties
 * files. The default should always be "bin/"!
 */

/***
 * 
	if (window["Clazz"] != null && Clazz.binaryFolders != null) {
	ClazzLoader.binaryFolders = Clazz.binaryFolders;
} else {
	ClazzLoader.binaryFolders = ["bin/", "", "j2slib/"];
}

ClazzLoader.addBinaryFolder = function (bin) {
	if (bin != null) {
		var bins = ClazzLoader.binaryFolders;
		for (var i = 0; i < bins.length; i++) {
			if (bins[i] == bin) {
				return;
			}
		}
		bins[bins.length] = bin;
	}
};
ClazzLoader.removeBinaryFolder = function (bin) {
	if (bin != null) {
		var bins = ClazzLoader.binaryFolders;
		for (var i = 0; i < bins.length; i++) {
			if (bins[i] == bin) {
				for (var j = i; j < bins.length - 1; j++) {
					bins[j] = bins[j + 1];
				}
				bins.length--;
				return bin;
			}
		}
	}
	return null;
};
ClazzLoader.setPrimaryFolder = function (bin) {
	if (bin != null) {
		ClazzLoader.removeBinaryFolder (bin);
		var bins = ClazzLoader.binaryFolders;
		for (var i = bins.length - 1; i >= 0; i--) {
			bins[i + 1] = bins[i];
		}
		bins[0] = bin;
	}
};

***/

/*# x<<
 # ClazzLoader.addBinaryFolder = Clazz.addBinaryFolder;
 # ClazzLoader.removeBinaryFolder = Clazz.removeBinaryFolder;
 # ClazzLoader.setPrimaryFolder = Clazz.setPrimaryFolder;
 #*/

/**
 * Indicate whether ClazzLoader is loading script synchronously or 
 * asynchronously.
 */
/* protected */
/*-# isAsynchronousLoading -> async #-*/
ClazzLoader.isAsynchronousLoading = true;

/* protected */
/*-# isUsingXMLHttpRequest -> xhr #-*/
ClazzLoader.isUsingXMLHttpRequest = false;

/* protected */
/*-# loadingTimeLag -> ltl #-*/
ClazzLoader.loadingTimeLag = -1;

/**
 * String mode:
 * asynchronous modes:
 * async(...).script, async(...).xhr, async(...).xmlhttprequest,
 * script.async(...), xhr.async(...), xmlhttprequest.async(...),
 * script
 * 
 * synchronous modes:
 * sync(...).xhr, sync(...).xmlhttprequest,
 * xhr.sync(...), xmlhttprequest.sync(...),
 * xmlhttprequest, xhr
 *                                                    
 * Integer mode:
 * Script/XHR bit: 1, 
 * 0: Script, 1: XHR
 * Asynchronous/Synchronous bit: 2
 * 0: Asynchronous, 2: Synchronous
 *
 * 0: Script & Asynchronous
 * 1: XHR & Asynchronous
 * 2: Script & Synchronous [Never]
 * 3: XHR & Synchronous
 */
/* public */
ClazzLoader.setLoadingMode = function (mode, timeLag) {
	if (mode == null) {
		if (ClazzLoader.isAsynchronousLoading && timeLag >= 0) {
			ClazzLoader.loadingTimeLag = timeLag;
		} else {
			ClazzLoader.loadingTimeLag = -1;
		}
		return;
	}
	if (typeof mode == "string") {
		mode = mode.toLowerCase ();
		if (mode.length == 0 || mode.indexOf ("script") != -1) {
			ClazzLoader.isUsingXMLHttpRequest = false;
			ClazzLoader.isAsynchronousLoading = true;
		} else {
			ClazzLoader.isUsingXMLHttpRequest = true;
			if (mode.indexOf ("async") != -1) {
				ClazzLoader.isAsynchronousLoading = true;
			} else {
				ClazzLoader.isAsynchronousLoading = false;
			}
		}
			ClazzLoader.isAsynchronousLoading = false; // BH
	/*# {$no.clazzloader.mode} >>x #*/
	} else {
		if (mode == ClazzLoader.MODE_SCRIPT) {
			ClazzLoader.isUsingXMLHttpRequest = false;
			ClazzLoader.isAsynchronousLoading = true;
		} else {
			ClazzLoader.isUsingXMLHttpRequest = true;
			if (mode == ClazzLoader.MODE_XHR_ASYNC) {
				ClazzLoader.isAsynchronousLoading = true;
			} else {
				ClazzLoader.isAsynchronousLoading = false;
			}
		}
	/*# x<< #*/
	}
	if (ClazzLoader.isAsynchronousLoading && timeLag >= 0) {
		ClazzLoader.loadingTimeLag = timeLag;
	} else {
		ClazzLoader.loadingTimeLag = -1;
	}
};

/*# {$no.clazzloader.mode} >>x #*/
ClazzLoader.MODE_SCRIPT = 0;
ClazzLoader.MODE_SCRIPT_ASYNC = 0;
ClazzLoader.MODE_XHR = 3;
ClazzLoader.MODE_XHR_SYNC = 3;
ClazzLoader.MODE_XHR_ASYNC = 1;


/*# x<< #*/
/**
 * Expand the shorten list of class names.
 * For example:
 * $wt.widgets.Shell, $.Display, $.Decorations
 * will be expanded to 
 * org.eclipse.swt.widgets.Shell, org.eclipse.swt.widgets.Display, 
 * org.eclipse.swt.widgets.Decorations ....
 * in which "$wt." stands for "org.eclipse.swt.", and "$." stands for
 * the previous class name's package.
 *
 * This method will be used to unwrap the required/optional classes list and 
 * the ignored classes list.
 */
/* private */
/*x-# unwrapArray -> uA #-x*/
ClazzLoader.unwrapArray = function (arr) {
	if (arr == null || arr.length == 0) {
		return arr;
	}
	var last = null;
	for (var i = 0; i < arr.length; i++) {
		if (arr[i] == null) {
			continue;
		}
		if (arr[i].charAt (0) == '$') {
			if (arr[i].charAt (1) == '.') {
				if (last == null) {
					continue;
				}
				var idx = last.lastIndexOf (".");
				if (idx != -1) {
					var prefix = last.substring (0, idx);
					arr[i] = prefix + arr[i].substring (1);
				}
			} else {
				arr[i] = "org.eclipse.s" + arr[i].substring (1);
			}
		}
		last = arr[i];
	}
	return arr;
};

/**
 * Used to keep to-be-loaded classes.
 */
/* private */
/*-# classQueue -> cq #-*/
ClazzLoader.classQueue = new Array ();

/* private */
/*-# classpathMap -> cm #-*/

//System.out.println("resetting classpathMap")
ClazzLoader.classpathMap = {};

/* public */
ClazzLoader.packageClasspath = function (pkg, base, index) {
	var map = ClazzLoader.classpathMap;

	/*
	 * In some situation, maybe,
	 * ClazzLoader.packageClasspath ("java", ..., true);
	 * is called after other ClazzLoader#packageClasspath, e.g.
	 * <code>
	 * ClazzLoader.packageClasspath ("org.eclipse.swt", "...", true);
	 * ClazzLoader.packageClasspath ("java", "...", true);
	 * </code>
	 * which is not recommended. But ClazzLoader should try to adjust orders
	 * which requires "java" to be declared before normal ClazzLoader
	 * #packageClasspath call before that line! And later that line
	 * should never initialize "java/package.js" again!
	 */
	var isPkgDeclared = (index == true && map["@" + pkg] != null);
	if (index && map["@java"] == null && pkg.indexOf ("java") != 0) {
		ClazzLoader.assurePackageClasspath ("java");
	}
	if (pkg instanceof Array) {
		ClazzLoader.unwrapArray (pkg);
		for (var i = 0; i < pkg.length; i++) {
			ClazzLoader.packageClasspath (pkg[i], base, index);
		}
		return;
	}
	if (pkg == "java" || pkg == "java.*") {
		// support ajax for default
		var key = "@net.sf.j2s.ajax";
		if (map[key] == null && base != null) {
			map[key] = base;
		}
		key = "@net.sf.j2s";
		if (map[key] == null && base != null) {
			map[key] = base;
		}
	} else if (pkg == "swt") { //abbrev
		pkg = "org.eclipse.swt";
	} else if (pkg == "ajax") { //abbrev
		pkg = "net.sf.j2s.ajax";
	} else if (pkg == "j2s") { //abbrev
		pkg = "net.sf.j2s";
	}
	if (pkg.lastIndexOf (".*") == pkg.length - 2) {
		pkg = pkg.substring (0, pkg.length - 2);
	}
	if (base != null) // critical for multiple applets
		map["@" + pkg] = base;
	if (index == true && window[pkg + ".registered"] != true && !isPkgDeclared) {
		ClazzLoader.pkgRefCount++;
		if (pkg == "java")pkg = "core" // JSmol -- moves java/package.js to core/package.js
		ClazzLoader.loadClass (pkg + ".package", function () {
					ClazzLoader.pkgRefCount--;
					if (ClazzLoader.pkgRefCount == 0) {
						ClazzLoader.runtimeLoaded ();
					}
				}, true);
	}
};

ClazzLoader.pkgRefCount = 0;

/**
 * Register classes to a given *.z.js path, so only a single *.z.js is loaded
 * for all those classes.
 */
/* public */
/*-# clazzes -> zs #-*/
ClazzLoader.jarClasspath = function (jar, clazzes) {
	if (!(clazzes instanceof Array))
		clazzes = [classes];
	ClazzLoader.unwrapArray (clazzes);
	for (var i = 0; i < clazzes.length; i++) {
		ClazzLoader.classpathMap["#" + clazzes[i]] = jar;
	}
	ClazzLoader.classpathMap["$" + jar] = clazzes;
};

/**
 * Usually be used in .../package.js. All given packages will be registered
 * to the same classpath of given prefix package.
 */
/* public */
ClazzLoader.registerPackages = function (prefix, pkgs) {
	//System.out.println("package " + prefix);
	ClazzLoader.checkInteractive ();
	var base = ClazzLoader.getClasspathFor (prefix + ".*", true);
	//System.out.println(base);
	for (var i = 0; i < pkgs.length; i++) {
		if (window["Clazz"] != null) {
			Clazz.declarePackage (prefix + "." + pkgs[i]);
		}
		ClazzLoader.packageClasspath (prefix + "." + pkgs[i], base);
	}
};

/**
 * Using multiple sites to load *.js in multiple threads. Using multiple
 * sites may avoid 2 HTTP 1.1 connections recommendation limit.
 * Here is a default implementation for http://archive.java2script.org.
 * In site archive.java2script.org, there are 6 sites:
 * 1. http://archive.java2script.org or http://a.java2script.org
 * 2. http://erchive.java2script.org or http://e.java2script.org
 * 3. http://irchive.java2script.org or http://i.java2script.org
 * 4. http://orchive.java2script.org or http://o.java2script.org
 * 5. http://urchive.java2script.org or http://u.java2script.org
 * 6. http://yrchive.java2script.org or http://y.java2script.org
 */
/* protected */
ClazzLoader.multipleSites = function (path) {
	var deltas = window["j2s.update.delta"];
	if (deltas != null && deltas instanceof Array && deltas.length >= 3) {
		var lastOldVersion = null;
		var lastNewVersion = null;
		for (var i = 0; i < deltas.length / 3; i++) {
			var oldVersion = deltas[i + i + i];
			if (oldVersion != "$") {
				lastOldVersion = oldVersion;
			}
			var newVersion = deltas[i + i + i + 1];
			if (newVersion != "$") {
				lastNewVersion = newVersion;
			}
			var relativePath = deltas[i + i + i + 2];
			var key = lastOldVersion + "/" + relativePath;
			var idx = path.indexOf (key);
			if (idx != -1 && idx == path.length - key.length) {
				path = path.substring (0, idx) + lastNewVersion + "/" + relativePath;
				break;
			}
		}
	}
	var length = path.length;
	if (ClazzLoader.maxLoadingThreads > 1 
			&& ((length > 15 && path.substring (0, 15) == "http://archive.")
			|| (length > 9 && path.substring (0, 9) == "http://a."))) {
		var index = path.lastIndexOf ("/");
		if (index < length - 3) {
			var arr = ['a', 'e', 'i', 'o', 'u', 'y'];
			var c1 = path.charCodeAt (index + 1);
			var c2 = path.charCodeAt (index + 2);
			var idx = (length - index) * 3 + c1 * 5 + c2 * 7; // Hash
			return path.substring (0, 7) + arr[idx % 6] + path.substring (8);
		}
	}
	return path;
};

/**
 * Return the *.js path of the given class. Maybe the class is contained
 * in a *.z.js jar file.
 * @param clazz Given class that the path is to be calculated for. May
 * be java.package, or java.lang.String
 * @param forRoot Optional argument, if true, the return path will be root
 * of the given classs' package root path.
 * @param ext Optional argument, if given, it will replace the default ".js"
 * extension.
 */
/* public */
ClazzLoader.getClasspathFor = function (clazz, forRoot, ext) {


	//System.out.println("check js path : " + arguments.callee.caller);
	//System.out.println("gcf " + clazz + " "  + forRoot + " " + ext);
	var path = ClazzLoader.classpathMap["#" + clazz];
	//System.out.println(path);
	var base = null;

	if (path != null) {
	//System.out.println("testing" + clazz + " " + forRoot + " " + ext)
		if (!forRoot && ext == null) { // return directly
			return ClazzLoader.multipleSites (path);
		} else {
			var idx = path.lastIndexOf (clazz.replace (/\./g, "/"));
			if (idx != -1) {
				base = path.substring (0, idx);
			} else {
				/*
				 * Check more: Maybe the same class' *.css is located
				 * in the same folder.
				 */
				idx = clazz.lastIndexOf (".");
				if (idx != -1) {
					idx = path.lastIndexOf (clazz.substring (0, idx)
							.replace (/\./g, "/"));
					if (idx != -1) {
						base = path.substring (0, idx);
					}
				}
			}
		}
	} else {
		/*
		path = ClazzLoader.classpathMap["@" + clazz]; // package
		if (path != null) {
			return ClazzLoader.assureBase (path) + clazz.replace (/\./g, "/") + "/";
		}
		*/
		//System.out.println("path was null")

		var idx = clazz.lastIndexOf (".");
		//*
		while (idx != -1) {
			var pkg = clazz.substring (0, idx);
			base = ClazzLoader.classpathMap["@" + pkg];
			//for (var jj in ClazzLoader.classpathMap)System.out.println("map[" + jj + "]=" + ClazzLoader.classpathMap[jj])
			//System.out.println("found " + base + " for @" + pkg)
			if (base != null) {
				break;
			}
			idx = clazz.lastIndexOf (".", idx - 2);
		}
		//*/
		/*
		if (idx != -1) {
			var pkg = clazz.substring (0, idx);
			base = ClazzLoader.classpathMap["@" + pkg];
		}
		//*/
	}



	//alert("gcf2 base=" + base);
	base = ClazzLoader.assureBase (base);
	//alert("gcf2za assurebase=" + base);
	if (forRoot) {
		return ClazzLoader.multipleSites (base);
	}
	if (clazz.lastIndexOf (".*") == clazz.length - 2) {
		var s = base + clazz.substring (0, idx + 1)
				.replace (/\./g, "/");
		//System.out.println("gcf2b " + s)
		return ClazzLoader.multipleSites (s);
	}
	if (ext == null) {
		ext = ".js";
	} else if (ext.charAt (0) != '.') {
		ext = "." + ext;
	}
	var jsPath = base + clazz.replace (/\./g, "/") + ext;

			//System.out.println("gcf2xx base=" + base + " clazz=" + clazz + " jsPath=" + jsPath)

	return ClazzLoader.multipleSites (jsPath);
};

/* private */
/*-# assureBase -> aB #-*/
ClazzLoader.assureBase = function (base) {

	if (base == null) {
		// Try to be compatiable with Clazz system.
		var bins = "binaryFolders";
		if (window["Clazz"] != null && Clazz[bins] != null
				&& Clazz[bins].length != 0) {
			base = Clazz[bins][0];
		} else if (ClazzLoader[bins] != null 
				&& ClazzLoader[bins].length != 0) {
			base = ClazzLoader[bins][0];
		} else {
			base = "j2s";
		}
	}
	if (base.lastIndexOf ("/") != base.length - 1) {
		base += "/";
	}
	return base;
};

/* Used to keep ignored classes */
/* private */
/*-# excludeClassMap -> exmap #-*/
ClazzLoader.excludeClassMap = {};

/**
 * To ignore some classes.
 */
/* public */
ClazzLoader.ignore = function () {
	var clazzes = null;
	if (arguments.length == 1) {
		if (arguments[0] instanceof Array) {
			clazzes = arguments[0];
		}
	}
	if (clazzes == null) {
		clazzes = new Array ();
		for (var i = 0; i < arguments.length; i++) {
			clazzes[clazzes.length] = arguments[i];
		}
	}
	ClazzLoader.unwrapArray (clazzes);
	for (var i = 0; i < clazzes.length; i++) {
		ClazzLoader.excludeClassMap["@" + clazzes[i]] = true;
	}
};

/* private */
/*-# isClassExcluded -> isEx #-*/
ClazzLoader.isClassExcluded = function (clazz) {
	return ClazzLoader.excludeClassMap["@" + clazz] == true;
};

/**
 * The following *.script* can be overriden to indicate the 
 * status of classes loading.
 *
 * TODO: There should be a Java interface with name like INativeLoaderStatus
 */
/* protected */
ClazzLoader.scriptLoading = function (file) {};

/* protected */
ClazzLoader.scriptLoaded = function (file) {};

/* protected */
ClazzLoader.scriptInited = function (file) {
};

/* protected */
ClazzLoader.scriptCompleted = function (file) {};

/* protected */
ClazzLoader.classUnloaded = function (clazz) {};

/* protected */
ClazzLoader.classReloaded = function (clazz) {};

/**
 * After all the classes are loaded, this method will be called.
 * Should be overriden to run *.main([]).
 */
/* protected */
ClazzLoader.globalLoaded = function () {};

/* protected */
ClazzLoader.keepOnLoading = true;

/* private */
/*-# mapPath2ClassNode -> p2node #-*/
ClazzLoader.mapPath2ClassNode = {};

/* private */
ClazzLoader.xhrOnload = function (transport, file) {

	if (transport.status >= 400 || transport.responseText == null
			|| transport.responseText.length == 0) { // error
			Clazz.alert("xhronload error" + transport.responseText);
		var fs = ClazzLoader.failedScripts;
		if (fs[file] == null) {
			// Silently take another try for bad network
			fs[file] = 1;
			ClazzLoader.loadedScripts[file] = false;
			ClazzLoader.loadScript (file, "xhrOnload 2nd try");
			return;
		} else {
			Clazz.alert ("[Java2Script] ClazzLoader.xhrOnload Error in loading " + file + "!");
		}
		ClazzLoader.tryToLoadNext (file);
	} else {
		ClazzLoader.evaluate(file, transport.responseText);
	}
};

ClazzLoader.evaluate = function(file, js) {
 		try {
			eval(js);
		} catch (e) {
			var s = "[Java2Script] Script error: " + e.message + " \n" + file + "\n\n" + js;
			alert(s)
			Clazz.alert (s);
			throw e;
		}
		ClazzLoader.scriptLoaded (file);
		ClazzLoader.tryToLoadNext (file);
}
/**
 * Empty onreadystatechange for fixing IE's memeory leak on XMLHttpRequest
 */
/* private */
/*-# emptyOnRSC -> rsc #-*/
ClazzLoader.emptyOnRSC = function () {
};

/* protected */
/*-# failedScripts -> fss #-*/
ClazzLoader.failedScripts = {};

/* protected */
/*-# failedHandles -> fhs #-*/
ClazzLoader.failedHandles = {};

/* protected */
ClazzLoader.takeAnotherTry = true;

/* private */
/*-# generateRemovingFunction -> gRF #-*/
ClazzLoader.generateRemovingFunction = function (node) {
	return function () {
		if (node.readyState != "interactive") {
			try {
				if (node.parentNode != null) {
					//alert("removing script " + node.src);
					node.parentNode.removeChild (node);
				}
			} catch (e) { }
			node = null;
		}
	};
};

/*
 * Dynamically SCRIPT elements are removed after they are parsed into memory.
 * And removed *.js may not be fetched again by "Refresh" action. And it will
 * save loading time for Java2Script applications.
 *
 * This may disturb debugging tools such as Firebug. Setting
 * window["j2s.script.debugging"] = true;
 * will ignore removing SCRIPT elements requests.
 */
/* private */
/*-# removeScriptNode -> RsN #-*/
ClazzLoader.removeScriptNode = function (n) {
	if (window["j2s.script.debugging"]) {
		return;
	}
	// lazily remove script nodes.
	window.setTimeout (ClazzLoader.generateRemovingFunction (n), 1);
};

/* private */
/*-# generatingXHROnload -> gXOd #-*/
ClazzLoader.generatingXHROnload = function (transport, file) {
	return function () {
		ClazzLoader.xhrOnload (transport, file);
		transport = null;
		file = null;
	};
};

/* private */
/*-# generatingXHRCallback -> gXcb #-*/
ClazzLoader.generatingXHRCallback = function (transport, file) {
	return function () {
		if (transport.readyState == 4) {
			if (ClazzLoader.inLoadingThreads > 0) {
				ClazzLoader.inLoadingThreads--;
			}
			var lazyFun = ClazzLoader.generatingXHROnload (transport, file);
			if (isActiveX) {
				transport.onreadystatechange = ClazzLoader.emptyOnRSC;
				// For IE, try to avoid stack overflow errors
				window.setTimeout (lazyFun,
						ClazzLoader.loadingTimeLag < 0 ? 0 : ClazzLoader.loadingTimeLag);
			} else {
				transport.onreadystatechange = null;
				if (ClazzLoader.loadingTimeLag >= 0) {
					window.setTimeout (lazyFun, ClazzLoader.loadingTimeLag);
				} else {
					ClazzLoader.xhrOnload (transport, file);
				}
			}
			transport = null;
			file = null;
		}
	};
};

/* private */
/*-# loadingNextByPath -> lNBP #-*/
ClazzLoader.loadingNextByPath = function (path) {
	if (ClazzLoader.loadingTimeLag >= 0) {
		window.setTimeout (function () {
					ClazzLoader.tryToLoadNext (path);
				}, ClazzLoader.loadingTimeLag);
	} else {
		ClazzLoader.tryToLoadNext (path);
	}
};

/* private */
/*-# ieToLoadScriptAgain -> iTLA #-*/
ClazzLoader.ieToLoadScriptAgain = function (path, local) {
	var fun = function () {
		if (!ClazzLoader.takeAnotherTry) {
			return;
		}
		// next time in "loading" state won't get waiting!
		ClazzLoader.failedScripts[path] = 0;
		ClazzLoader.loadedScripts[path] = false;
		// failed count down!
		if (ClazzLoader.inLoadingThreads > 0) {
			ClazzLoader.inLoadingThreads--;
		}
		// Take another try!
		// log ("re - loading ... " + path);
		ClazzLoader.loadScript (path, "ietoloadscriptAgain");
	};
	// consider 30 seconds available after failing!
	/*
	 * Set 1s waiting in local file system. Is it 1s enough?
	 * What about big *.z.js need more than 1s to initialize?
	 */
	var waitingTime = (local ? 500 : 15000); // 0.5s : 15s
	//alert ("waiting:" + waitingTime + " . " + path);
	return window.setTimeout (fun, waitingTime);
};

/* private */
/*-# w3cFailedLoadingTest -> wFLT #-*/
ClazzLoader.w3cFailedLoadingTest = function (script) {
	return window.setTimeout (function () {
				script.onerror ();
				script.timeoutHandle = null;
				script = null;
			}, 500); // 0.5s for loading a local file is considered enough long
};

/* private */
/*-# generatingW3CScriptOnCallback -> gWSC #-*/
ClazzLoader.generatingW3CScriptOnCallback = function (path, forError) {
	return function () {
	if (forError && Clazz.__debuggingBH)Clazz.alert("############ forError=" + forError + " path=" + path + " ####" + (forError ? "NOT" : "") + "LOADED###");

		if (ClazzLoader.isGecko && this.timeoutHandle != null) {
			window.clearTimeout (this.timeoutHandle);
			this.timeoutHandle = null;
		}
		if (ClazzLoader.inLoadingThreads > 0) {
			ClazzLoader.inLoadingThreads--;
		}
		this.onload = null;
		this.onerror = null;
		var checkOpera = false;   // BH - disabled this for Opera 10
		if (!forError && checkOpera && ClazzLoader.isOpera // Opera only
				&& !ClazzLoader.innerLoadedScripts[this.src]) {
			ClazzLoader.checkInteractive();
		}
		if (forError || (checkOpera && !ClazzLoader.innerLoadedScripts[this.src]
					&& ClazzLoader.isOpera)) { 
			// Opera will not take another try.
			var fss = ClazzLoader.failedScripts;
			if (fss[path] == null && ClazzLoader.takeAnotherTry) {
				// silently take another try for bad network
				//alert ("re loading " + path + " ... ");
				fss[path] = 1;
				if (!forError) {
					ClazzLoader.innerLoadedScripts[this.src] = false;
				}
				ClazzLoader.loadedScripts[path] = false;
				ClazzLoader.loadScript (path, "w3c script failed");
				ClazzLoader.removeScriptNode (this);
				//alert("............[Java2Script] "  + " inLoadingThreads=" + ClazzLoader.inLoadingThreads + " ClazzLoader.generatingW3CScriptOnCallback loading " + path);
				return;
			} else {
				//alert("[Java2Script] ClazzLoader.generatingW3CScriptOnCallback Error in loading " + path + "! inLoadingThreads=" + ClazzLoader.inLoadingThreads);
			}
			if (forError) {
				ClazzLoader.scriptLoaded (path);
			}
		} else {
			//System.out.println("path loaded: " + path)
			ClazzLoader.scriptLoaded (path);
		}
		ClazzLoader.loadingNextByPath (path);
		ClazzLoader.removeScriptNode (this);
	};
};

/* private */
/*-# generatingIEScriptOnCallback -> gISC #-*/
ClazzLoader.generatingIEScriptOnCallback = function (path) {
	return function () {
		var fhs = ClazzLoader.failedHandles;
		var fss = ClazzLoader.failedScripts;
		var state = "" + this.readyState;


		var local = state == "loading" 
				&& (this.src.indexOf ("file:") == 0 
				|| (window.location.protocol == "file:"
				&& this.src.indexOf ("http") != 0));

		// alert (state + "/" + this.src);
		if (state != "loaded" && state != "complete") {
			/*
			 * When no such *.js existed, IE will be
			 * stuck here without loaded event!
			 */
			/*
			if ((window.location.protocol != "file:" 
					&& this.src.indexOf ("file:") != 0)
					|| this.readyState != "loading") {
				return;
			}
			*/
			if (fss[path] == null) {
				fhs[path] = ClazzLoader.ieToLoadScriptAgain (path, local);
				return;
			}
			if (fss[path] == 1) { // above function will be executed?!
				return;
			}
		}
		if (fhs[path] != null) {
			window.clearTimeout (fhs[path]);
			fhs[path] = null;
		}
		if ((local || state == "loaded")
				&& !ClazzLoader.innerLoadedScripts[this.src]) {
			if (!local && (fss[path] == null || fss[path] == 0)
					&& ClazzLoader.takeAnotherTry) {
				// failed! count down
				if (ClazzLoader.inLoadingThreads > 0) {
					ClazzLoader.inLoadingThreads--;
				}
				// silently take another try for bad network
				fss[path] = 1;
				// log ("reloading ... " + path);
				ClazzLoader.loadedScripts[path] = false;
				ClazzLoader.loadScript (path, "local or loaded failed");
				ClazzLoader.removeScriptNode (this);
				return;
			} else {
				Clazz.alert ("[Java2Script] generatingIEScriptOnCallback Error in loading " + path + "!");
			}
		}
		if (ClazzLoader.inLoadingThreads > 0) {
			ClazzLoader.inLoadingThreads--;
		}
		ClazzLoader.scriptLoaded (path);
		// Unset onreadystatechange, leaks mem in IE
		this.onreadystatechange = null; 
		ClazzLoader.loadingNextByPath (path);
		ClazzLoader.removeScriptNode (this);
	};
};

/*
 * There is another thread trying to remove j2slib.z.js or similar SCRIPT.
 * See the end of this file.
 */


		Clazz.transport = null;

/**
 * Load *.js by adding script elements into head. Hook the onload event to
 * load the next class in dependency tree.
 */
/* protected */
/*-#
 # loadScript -> xrpt
 #
 # transport -> tt
 # isActiveX -> iX
 # ignoreOnload -> iol
 #-*/
ClazzLoader.loadScript = function (file, why) {
		Clazz.currentPath = file;
		//alert("loadScript" + file)
		//System.out.println(("call loadScript " + file.replace(/\//g,"\\") + " Z").replace(/j2s[\\\.]/g,""))
	// maybe some scripts are to be loaded without needs to know onload event.


 
	var ignoreOnload = (arguments[2] == true);
	if (ClazzLoader.loadedScripts[file] && !ignoreOnload) {
		ClazzLoader.tryToLoadNext (file);
		return;
	}
	ClazzLoader.loadedScripts[file] = true;
	/* also remove from those queue */
	var cq = ClazzLoader.classQueue;
	for (var i = 0; i < cq.length; i++) {
		if (cq[i] == file) {
			for (var j = i; j < cq.length - 1; j++) {
				cq[i] = cq[i + 1];
			}
			cq.length--;
			break;
		}
	}

	System.out.println("loading... " + file + (why ? " -- required by " + why : ""))

	if (ClazzLoader.isUsingXMLHttpRequest) {
		ClazzLoader.scriptLoading (file);
		//var transport = null;

		if (!ClazzLoader.isAsynchronousLoading) {
		 if (self.Jmol) {
			// works in MSIE locally :)
			 var data = Jmol._getFileData(file);
			} else {
				var info = {type:"GET",dataType:"text",async:false,url:file,cached:true};
				var xhr = jQuery.ajax(info);
				var data = xhr.responseText;
				if (data == null && xhr.state)
					data = "alert('error loading file " + file + " state=" + xhr.state() + "')";
			}
			ClazzLoader.evaluate(file, data); 
			return;
		}

		var isActiveX = false;
		if (!Clazz.transport) {
		if (window.XMLHttpRequest) {
			Clazz.transport = new XMLHttpRequest();
		} else {
			isActiveX = true;
			try {
				Clazz.transport = new ActiveXObject("Msxml2.XMLHTTP");
			} catch (e) {
				Clazz.transport = new ActiveXObject("Microsoft.XMLHTTP");
			}
		}

		if (Clazz.transport == null) { // should never happen in modern browsers.
			Clazz.alert ("[Java2Script] XMLHttpRequest not supported!");
			return;
		}
	}
	var transport = Clazz.transport

	//alert("transport=" + transport + " " + file)


	try{

		transport.open ("GET", file, ClazzLoader.isAsynchronousLoading);

	} catch (e) { alert(e + " " + file)}


		// transport.setRequestHeader ("User-Agent",
		// 		"Java2Script-Pacemaker/1.0 (+http://j2s.sourceforge.net)");
//		if (ClazzLoader.isAsynchronousLoading) {
			transport.onreadystatechange = ClazzLoader.generatingXHRCallback (
					transport, file);
			ClazzLoader.inLoadingThreads++;
			try {
				transport.send (null);
			} catch (e) {
				Clazz.alert ("[Java2Script] Loading file error: " + e.message);
				ClazzLoader.xhrOnload (transport, file);
				//throw e;
			}
//		} else {
//			try {
//				transport.send (null);
//			} catch (e) {
//				Clazz.alert ("[Java2Script] Loading file error: " + e.message);
//				//throw e;
//			}
//			ClazzLoader.xhrOnload (transport, file);
//		}
		return;
	}
	// Create script DOM element
	var script = document.createElement ("SCRIPT");
	script.type = "text/javascript";


	if (ClazzLoader.isChrome && ClazzLoader.reloadingClasses[file]) {
		script.src = file + "?" + Math.floor (100000 * Math.random ());
	} else {
		script.src = file;
	}

	var head = document.getElementsByTagName ("HEAD")[0];
	if (ignoreOnload) {
		head.appendChild (script);
		// ignore onload event and no call of ClazzLoader.scriptLoading
		return;
	}
	script.defer = true;
	// Alert when the script is loaded
	if (typeof (script.onreadystatechange) == "undefined" || !ClazzLoader.isIE) { // W3C
		if (ClazzLoader.isGecko && (file.indexOf ("file:") == 0 
				|| (window.location.protocol == "file:" && file.indexOf ("http") != 0))) {
			script.timeoutHandle = ClazzLoader.w3cFailedLoadingTest (script);
		}

		/*
		 * What about Safari and Google Chrome?
		 */
		/*
		 * Opera will trigger onload event even there are no *.js existed
		 */
		script.onload = ClazzLoader.generatingW3CScriptOnCallback (file, false);

		/*
		 * For Firefox/Mozilla, unexisted *.js will result in errors.
		 */
		script.onerror = ClazzLoader.generatingW3CScriptOnCallback (file, true);
		if (ClazzLoader.isOpera) {
			ClazzLoader.needOnloadCheck = true;
		}
	} else { // IE
		ClazzLoader.needOnloadCheck = true;
		script.onreadystatechange = ClazzLoader.generatingIEScriptOnCallback (file);
	}
	ClazzLoader.inLoadingThreads++;
	//alert("threads:"+ClazzLoader.inLoadingThreads);
	// Add script DOM element to document tree
	head.appendChild (script);
	ClazzLoader.scriptLoading (file);
};

/* protected */
ClazzLoader.isResourceExisted = function (id, path, base) {
	if (id != null && document.getElementById (id) != null) {
		return true;
	}
	if (path != null) {
		var key = path;
		if (base != null) {
			if (path.indexOf (base) == 0) {
				key = path.substring (base.length);
			}
		}
		if (path.lastIndexOf (".css") == path.length - 4) {
			var resLinks = document.getElementsByTagName ("LINK");
			for (var i = 0; i < resLinks.length; i++) {
				var cssPath = resLinks[i].href;
				var idx = cssPath.lastIndexOf (key);
				if (idx != -1 && idx == cssPath.length - key.length) {
					return true;
				}
			}
			if (window["css." + id] == true) {
				return true;
			}
		} else if (path.lastIndexOf (".js") == path.length - 4) {
			var resScripts = document.getElementsByTagName ("SCRIPT");
			for (var i = 0; i < resScripts.length; i++) {
				var jsPath = resScripts[i].src;
				var idx = jsPath.lastIndexOf (key);
				if (idx != -1 && idx == jsPath.length - key.length) {
					return true;
				}
			}
		}
	}
	return false;
};

/* private */
/*-# queueBe4SWT -> q4T #-*/
ClazzLoader.queueBe4SWT = [];

/* private */
/*-# lockQueueBe4SWT -> l4T #-*/
ClazzLoader.lockQueueBe4SWT = true;

/* private */
/*-# isLoadingEntryClass -> lec #-*/
ClazzLoader.isLoadingEntryClass = true;

/* private */
/*-# besidesJavaPackage -> bJP #-*/
ClazzLoader.besidesJavaPackage = false;

/**
 * After class is loaded, this method will be executed to check whether there
 * are classes in the dependency tree that need to be loaded.
 */
/* private */
/*-# tryToLoadNext -> next #-*/
ClazzLoader.tryToLoadNext = function (file) {
	/*
	 * Try to check whether current status is in SWT lazy loading mode. If
	 * yes, try to keep ClazzLoader#tryToLoadNext in queue and wait until
	 * all SWT core is loaded.
	 */
	if (ClazzLoader.lockQueueBe4SWT && ClazzLoader.pkgRefCount != 0 
			&& file.lastIndexOf ("package.js") != file.length - 10
			&& !ClazzLoader.isOpera) { // No Opera! Opera is in single thread.
		var qbs = ClazzLoader.queueBe4SWT;
		qbs[qbs.length] = file;
		return;
	}

	var node = ClazzLoader.mapPath2ClassNode["@" + file];
	if (node == null) { // maybe class tree root
		//error (" null node ?" + file);
		return;
	}
	var clazzes = ClazzLoader.classpathMap["$" + file];
	if (clazzes != null) {
		for (var i = 0; i < clazzes.length; i++) {
			var nm = clazzes[i];
			if (nm != node.name) {
				var n = ClazzLoader.findClass (nm);
				if (n != null) {
				//*
					if (n.status < ClazzNode.STATUS_CONTENT_LOADED) {
						n.status = ClazzNode.STATUS_CONTENT_LOADED;
						ClazzLoader.updateNode (n);
					}
				//*/
				//ClazzLoader.tryToLoadNext (n.path);
				} else {
					n = new ClazzNode ();
					n.name = nm;
					var pp = ClazzLoader.classpathMap["#" + nm];
					if (pp == null) {
						alert (nm);
						error ("Java2Script implementation error! Please report this bug!");
					}
					//n.path = ClazzLoader.lastScriptPath;
					n.path = pp;
					//error ("..." + node.path + "//" + node.name);
					ClazzLoader.mappingPathNameNode (n.path, nm, n);
					n.status = ClazzNode.STATUS_CONTENT_LOADED;
					ClazzLoader.addChildClassNode(ClazzLoader.clazzTreeRoot, n, -1);
					ClazzLoader.updateNode (n);
				}
			}
		}
	}
	if (node instanceof Array) {
		//log ("array of node " + node.length + ">>>>" + file);
		/*
		for (var i = 0; i < node.length; i++) {
			//alert ("array of node " + node[i].name);
		}
		*/
		for (var i = 0; i < node.length; i++) {
			if (node[i].status < ClazzNode.STATUS_CONTENT_LOADED) {
				node[i].status = ClazzNode.STATUS_CONTENT_LOADED;
				//error ("updating array : " + node[i].name + "..");
				ClazzLoader.updateNode (node[i]);
			}
		}
	} else {
		if (node.status < ClazzNode.STATUS_CONTENT_LOADED) {
			var stillLoading = false;
			var ss = document.getElementsByTagName ("SCRIPT");
			for (var i = 0; i < ss.length; i++) {
				if (ClazzLoader.isIE) {
					if (ss[i].onreadystatechange != null && ss[i].onreadystatechange.path == node.path
							&& ss[i].readyState == "interactive") {
						stillLoading = true;
						break;
					}
				} else {
					if (ss[i].onload != null && ss[i].onload.path == node.path) {
						stillLoading = true;
						break;
					}
				}
			}
			if (!stillLoading) {
				node.status = ClazzNode.STATUS_CONTENT_LOADED;
				ClazzLoader.updateNode (node);
			}
		}
	}
	/*
	 * Maybe in #optinalLoaded inside above ClazzLoader#updateNode calls, 
	 * ClazzLoader.keepOnLoading is set false (Already loaded the wanted
	 * classes), so here check to stop.
	 */
	 
	if (!ClazzLoader.keepOnLoading) {
		return;
	}

	var loadFurther = false;
	var n = ClazzLoader.findNextMustClass (ClazzLoader.clazzTreeRoot, 
			ClazzNode.STATUS_KNOWN);
	//System.out.println(file + " next ..." + n) ;
	if (n != null) {
		//log ("next ..." + n.name);
		ClazzLoader.loadClassNode (n);
		while (ClazzLoader.inLoadingThreads < ClazzLoader.maxLoadingThreads) {
			var nn = ClazzLoader.findNextMustClass (ClazzLoader.clazzTreeRoot,
					ClazzNode.STATUS_KNOWN);
			if (nn == null) break;
			ClazzLoader.loadClassNode (nn); // will increase inLoadingThreads!
		}
	} else {
		var cq = ClazzLoader.classQueue;
		if (cq.length != 0) { 
			/* queue must be loaded in order! */
			n = cq[0]; // popup class from the queue
			//alert ("load from queue");
			//alert (cq.length + ":" + cq);
			for (var i = 0; i < cq.length - 1; i++) {
				cq[i] = cq[i + 1];
			}
			cq.length--;
			//log (cq.length + ":" + cq);
			if (!ClazzLoader.loadedScripts[n.path] || cq.length != 0 
					|| !ClazzLoader.isLoadingEntryClass
					|| (n.musts != null && n.musts.length != 0)
					|| (n.optionals != null && n.optionals.length != 0)/*
					|| window["org.eclipse.swt.registered"] != null*/) {
				ClazzLoader.addChildClassNode(ClazzLoader.clazzTreeRoot, n, 1);
				ClazzLoader.loadScript (n.path, n.requiredBy);
				//alert("part1")
			} else {
				if (ClazzLoader.isLoadingEntryClass) {
					/*
					 * The first time reaching here is the time when ClassLoader
					 * is trying to load entry class. Class with #main method and
					 * is to be executed is called Entry Class.
					 *
					 * Here when loading entry class, ClassLoader should not call
					 * the next following loading script. This is because, those
					 * scripts will try to mark the class as loaded directly and
					 * then continue to call #optionalsLoaded callback method,
					 * which results in an script error!
					 */
					ClazzLoader.isLoadingEntryClass = false;
				}
				//alert ("Continue loading by SCRIPT onload event!");
				//ClazzLoader.addChildClassNode(ClazzLoader.clazzTreeRoot, n, 1);
				//ClazzLoader.loadScript (n.path);
			}
		} else { // Optionals
			n = ClazzLoader.findNextOptionalClass (ClazzNode.STATUS_KNOWN);
			//log ("options " + file);
			if (n != null) {
				//System.out.println("in optionals unknown..." + file + " " + n.name + " " + ClazzLoader.inLoadingThreads + "/" + ClazzLoader.maxLoadingThreads);
				ClazzLoader.loadClassNode (n);
				while (ClazzLoader.inLoadingThreads < ClazzLoader.maxLoadingThreads) {
					var nn = ClazzLoader.findNextOptionalClass (ClazzNode.STATUS_KNOWN);
					//log ("in second loading " + nn);
					if (nn == null) break;
				//System.out.println("in optionals unknown2..." + file + " " + nn.name);
					ClazzLoader.loadClassNode (nn); // will increase inLoadingThreads!
				}
			} else {
				//System.out.println("no optionals");
				loadFurther = true;
			}
		}
	}
		//log("test3 " + loadFurther + " " + ClazzLoader.inLoadingThreads )
	/*
	 * The following codes still need more tests, e.g. cyclic tests.
	 * And they also need optimization.
	 */
	if (loadFurther && ClazzLoader.inLoadingThreads == 0) {
		while ((n = ClazzLoader.findNextMustClass (ClazzLoader.clazzTreeRoot, ClazzNode.STATUS_CONTENT_LOADED)) != null) {
			ClazzLoader.updateNode (n);
		}
		var lastNode = null;
		while ((n = ClazzLoader.findNextOptionalClass (ClazzNode.STATUS_CONTENT_LOADED)) != null) {
			if (lastNode === n) { // Already existed cycle ?
				n.status = ClazzNode.STATUS_OPTIONALS_LOADED;
			}
			ClazzLoader.updateNode (n);
			lastNode = n;
		}
		while (true) {
			ClazzLoader.tracks = new Array ();
			if (!ClazzLoader.checkOptionalCycle (ClazzLoader.clazzTreeRoot)) {
				break;
			}
		}
		lastNode = null;
		while ((n = ClazzLoader.findNextMustClass (ClazzLoader.clazzTreeRoot, ClazzNode.STATUS_DECLARED)) != null) {
			if (lastNode === n) break;
			ClazzLoader.updateNode (n);
			lastNode = n;
		}
		lastNode = null;
		while ((n = ClazzLoader.findNextOptionalClass (ClazzNode.STATUS_DECLARED)) != null) {
			if (lastNode === n) break;
			ClazzLoader.updateNode (n);
			lastNode = n;
		}
		var dList = [];
		while ((n = ClazzLoader.findNextMustClass (ClazzLoader.clazzTreeRoot, ClazzNode.STATUS_DECLARED)) != null) {
			dList[dList.length] = n;
			n.status = ClazzNode.STATUS_OPTIONALS_LOADED;
		}
		while ((n = ClazzLoader.findNextOptionalClass (ClazzNode.STATUS_DECLARED)) != null) {
			dList[dList.length] = n;
			n.status = ClazzNode.STATUS_OPTIONALS_LOADED;
		}
		for (var i = 0; i < dList.length; i++) {
			// ClazzLoader.updateNode (dList[i]);
			ClazzLoader.destroyClassNode (dList[i]); // Same as above
		}
		for (var i = 0; i < dList.length; i++) {
			var optLoaded = dList[i].optionalsLoaded;
			if (optLoaded != null) {
				dList[i].optionalsLoaded = null;
				//window.setTimeout (optLoaded, 25);
				optLoaded ();
			}
		}
		/*
		 * It seems ClazzLoader#globalLoaded is seldom overrided.
		 */
		ClazzLoader.globalLoaded ();
		//error ("end ?");
	}
};

ClazzLoader.tracks = new Array ();

/*
 * There are classes reference cycles. Try to detect and break those cycles.
 * TODO: Reference cycles should be broken down carefully. Or there will be 
 * bugs!
 */
/* protected */
ClazzLoader.checkOptionalCycle = function (node) {
	var ts = ClazzLoader.tracks;
	var length = ts.length;
	var cycleFound = -1;
	for (var i = 0; i < ts.length; i++) {
		if (ts[i] === node && ts[i].status >= ClazzNode.STATUS_DECLARED) { 
			// Cycle is found;
			cycleFound = i;
			break;
		}
	}
	ts[ts.length] = node;
	if (cycleFound != -1) {
		/*
		for (var i = cycleFound; i < ts.length; i++) {
			//alert (ts[i].name + ":::" + ts[i].status);
		}
		//alert ("===");
		*/
		for (var i = cycleFound; i < ts.length; i++) {
			ts[i].status = ClazzNode.STATUS_OPTIONALS_LOADED;
			//ClazzLoader.updateNode (ts[i]);
			ClazzLoader.destroyClassNode (ts[i]); // Same as above
			for (var k = 0; k < ts[i].parents.length; k++) {
				//log ("updating parent ::" + ts[i].parents[k].name);
				ClazzLoader.updateNode (ts[i].parents[k]);
			}
			ts[i].parents = new Array ();
			var optLoaded = ts[i].optionalsLoaded;
			if (optLoaded != null) {
				ts[i].optionalsLoaded = null;
				//alert ("check cycle.");
				//window.setTimeout (optLoaded, 25);
				optLoaded ();
			}
		}
		ts.length = 0;
		return true;
	}
	for (var i = 0; i < node.musts.length; i++) {
		if (node.musts[i].status == ClazzNode.STATUS_DECLARED) {
			if (ClazzLoader.checkOptionalCycle (node.musts[i])) {
				return true;
			}
		}
	}
	for (var i = 0; i < node.optionals.length; i++) {
		if (node.optionals[i].status == ClazzNode.STATUS_DECLARED) {
			if (ClazzLoader.checkOptionalCycle (node.optionals[i])) {
				return true;
			}
		}
	}
	ts.length = length;
	return false;
};


/**
 * Update the dependency tree nodes recursively.
 */
/* private */
/*-# 
 # updateNode -> uN 
 #
 # isMustsOK -> mOK
 # isOptionsOK -> oOK
 #-*/
ClazzLoader.updateNode = function (node) {
	if (node.name == null 
			|| node.status >= ClazzNode.STATUS_OPTIONALS_LOADED) {
		//System.out.println("destroying node " + node.name + " " + node.path)
		ClazzLoader.destroyClassNode (node);
		return;
	}
	var isMustsOK = false;
	if (node.musts == null || node.musts.length == 0 
			|| node.declaration == null) {
		isMustsOK = true;
	} else {
		isMustsOK = true;
		var mustLength = node.musts.length;
		for (var i = mustLength - 1; i >= 0; i--) {

			/*
			 * Soheil reported a strange bug:
The problem here is that, Widget functions and field are not added to
the Control !
Control is not an instance of Widget! I got an error that says
this.checkOrientation is not a function ( in Control's constructor).

When I changed the order of Drawable and Widget in the definition of
Control's constructor, the line bellow, it works fine!
$_L(["$wt.graphics.Drawable","$wt.widgets.Widget"],"$wt.widgets.Control",
... : has the error
$_L(["$wt.widgets.Widget","$wt.graphics.Drawable"],"$wt.widgets.Control", 
... : does not have the error 
			 *
			 * In the bug fix procedure, it's known that node.musts will
			 * be changed according to the later codes:
			 * ClazzLoader.updateNode (n); // (see about 20 lines below)
			 * 
			 * As node.musts may become smaller, node.musts should be 
			 * traversed in reverse order, so all musts are checked.
			 *
			 * TODO:
			 */
			var n = node.musts[i];
			n.requiredBy = node;
			if (n.status < ClazzNode.STATUS_DECLARED) {
				if (ClazzLoader.isClassDefined (n.name)) {
					var nns = new Array (); // for optional loaded events!

		//System.out.println("destroying node2 " + n.name + " " + n.path)

					n.status = ClazzNode.STATUS_OPTIONALS_LOADED;
					// ClazzLoader.updateNode (n); // musts may be changed
					ClazzLoader.destroyClassNode (n); // Same as above
					/*
					 * For those classes within one *.js file, update
					 * them synchronously.
					 */
					if (n.declaration != null 
							&& n.declaration.clazzList != null) {
						var list = n.declaration.clazzList;
						for (var j = 0; j < list.length; j++) {
							var nn = ClazzLoader.findClass (list[j]);
							if (nn.status != ClazzNode.STATUS_OPTIONALS_LOADED
									&& nn !== n) {
								nn.status = n.status;
								nn.declaration = null;
								// ClazzLoader.updateNode (nn);
								// Same as above
		//System.out.println("destroying node3 " + nn.name + " " + nn.path)

								ClazzLoader.destroyClassNode (nn);
								if (nn.optionalsLoaded != null) {
									nns[nns.length] = nn;
								}
							}
						}
						n.declaration = null;
					}
					if (n.optionalsLoaded != null) {
						nns[nns.length] = n;
					}
					for (var j = 0; j < nns.length; j++) {
						var optLoaded = nns[j].optionalsLoaded;
						if (optLoaded != null) {
							nns[j].optionalsLoaded = null;
							//window.setTimeout (optLoaded, 25);
							optLoaded ();
						}
					}
				} else { // why not break? -Zhou Renjian @ Nov 28, 2006
					if (n.status == ClazzNode.STATUS_CONTENT_LOADED) {
						// lazy loading script doesn't work! - 2/26/2007
						ClazzLoader.updateNode (n); // musts may be changed
					}
					if (n.status < ClazzNode.STATUS_DECLARED) {
						isMustsOK = false;
					}
				}
				// fix the above strange bug
				if (node.musts.length != mustLength) {
					// length changed!
					mustLength = node.musts.length;
					i = mustLength; // -1
					isMustsOK = true;
				}
			}
		}
	}
	/*
	var scripts = document.getElementsByTagName ("SCRIPT");
	var count = 0;
	for (var i = 0; i < scripts.length; i++) {
		var s = scripts[i];
		if (s.onload != null) {
			log ("---:---" + s.src);
			count++;
		}
	}
	//alert ("There are " + count + " script loading ...");
	*/
	//System.out.println("testing " + node.name + " " + isMustsOK + " " + node.status + " " +  ClazzNode.STATUS_DECLARED + "isnull?" + (node.declaration == null))

	if (isMustsOK) {
		if (node.status < ClazzNode.STATUS_DECLARED) {
			var decl = node.declaration;
			if (decl != null) {
				if (decl.executed == false) {
					decl ();
					decl.executed = true;
				} else {
					decl ();
				}
			}
			node.status = ClazzNode.STATUS_DECLARED;
			if (ClazzLoader.definedClasses != null) {
				ClazzLoader.definedClasses[node.name] = true;
			}
			ClazzLoader.scriptInited (node.path);
					/*
					 * For those classes within one *.js file, update
					 * them synchronously.
					 */
					if (node.declaration != null 
							&& node.declaration.clazzList != null) {
						var list = node.declaration.clazzList;
						for (var j = 0; j < list.length; j++) {
							var nn = ClazzLoader.findClass (list[j]);
							if (nn.status != ClazzNode.STATUS_DECLARED
									&& nn !== node) {
			nn.status = ClazzNode.STATUS_DECLARED;
			if (ClazzLoader.definedClasses != null) {
				ClazzLoader.definedClasses[nn.name] = true;
			}
			ClazzLoader.scriptInited (nn.path);
							}
						}
					}
		}
		var level = ClazzNode.STATUS_DECLARED;
		var isOptionsOK = false;

		if (((node.optionals == null || node.optionals.length == 0) 
				&& (node.musts == null || node.musts.length == 0))
				|| (node.status > ClazzNode.STATUS_KNOWN 
				&& node.declaration == null)) {
			isOptionsOK = true;
		} else {
			isOptionsOK = true;


			for (var i = 0; i < node.musts.length; i++) {
				var n = node.musts[i];
				if (n.status < ClazzNode.STATUS_OPTIONALS_LOADED) {
					isOptionsOK = false;
					break;
				}
			}
			if (isOptionsOK) {
				for (var i = 0; i < node.optionals.length; i++) {
					var n = node.optionals[i];
					if (n.status < ClazzNode.STATUS_OPTIONALS_LOADED) {
						isOptionsOK = false;
						break;
					}
				}
			}
		}
		if (isOptionsOK) {
			level = ClazzNode.STATUS_OPTIONALS_LOADED;
			node.status = level;
			ClazzLoader.scriptCompleted (node.path);
			var optLoaded = node.optionalsLoaded;
			if (optLoaded != null) {
				node.optionalsLoaded = null;
				//window.setTimeout (optLoaded, 25);
				optLoaded ();
				if (!ClazzLoader.keepOnLoading) {
					return false;
				}
			}
			ClazzLoader.destroyClassNode (node);
					/*
					 * For those classes within one *.js file, update
					 * them synchronously.
					 */
					if (node.declaration != null 
							&& node.declaration.clazzList != null) {
						var list = node.declaration.clazzList;
						for (var j = 0; j < list.length; j++) {
							var nn = ClazzLoader.findClass (list[j]);
							if (nn.status != level && nn !== node) {
			nn.status = level;
			nn.declaration = null;
			ClazzLoader.scriptCompleted (nn.path);
			var optLoaded = nn.optionalsLoaded;
			if (optLoaded != null) {
				nn.optionalsLoaded = null;
				//window.setTimeout (optLoaded, 25);
				optLoaded ();
				if (!ClazzLoader.keepOnLoading) {
					return false;
				}
			}
			ClazzLoader.destroyClassNode (node);
							}
						}
					}
		}
		ClazzLoader.updateParents (node, level);
	}
};

/* private */
/*-# updateParents -> uP #-*/
ClazzLoader.updateParents = function (node, level) {
	if (node.parents == null || node.parents.length == 0) {
		return;
	}
	for (var i = 0; i < node.parents.length; i++) {
		var p = node.parents[i];
		if (p.status >= level) {
			continue;
		}
		ClazzLoader.updateNode (p);
	}
	if (level == ClazzNode.STATUS_OPTIONALS_LOADED) {
		node.parents = new Array ();
	}
};

/* private */
/*-# findNextMustClass -> fNM #-*/
ClazzLoader.findNextMustClass = function (node, status) {
	if (node != null) {
		/*
		if (ClazzLoader.isClassDefined (node.name)) {
			node.status = ClazzNode.STATUS_OPTIONALS_LOADED;
		}
		*/
		if (node.musts != null && node.musts.length != 0) {
			for (var i = 0; i < node.musts.length; i++) {
				var n = node.musts[i];
				if (n.status == status && (status != ClazzNode.STATUS_KNOWN 
						|| ClazzLoader.loadedScripts[n.path] != true)
						&& (status == ClazzNode.STATUS_DECLARED
						|| !ClazzLoader.isClassDefined (n.name))) {
					return n;
				} else {
					var nn = ClazzLoader.findNextMustClass (n, status);
					if (nn != null) {
						return nn;
					}
				}
			}
		}
		if (node.status == status && (status != ClazzNode.STATUS_KNOWN 
				|| ClazzLoader.loadedScripts[node.path] != true)
				&& (status == ClazzNode.STATUS_DECLARED
				|| !ClazzLoader.isClassDefined (node.name))) {
			return node;
		}
	}
	return null;
};

/*
 * Be used to record already used random numbers. And next new random
 * number should not be in the property set.
 */
/* private */
/*-# usedRandoms -> Rms #-*/
ClazzLoader.usedRandoms = {};
ClazzLoader.usedRandoms["r" + 0.13412] = 0.13412;

/* private */
/*-# findNextOptionalClass -> fNO #-*/
ClazzLoader.findNextOptionalClass = function (status) {
	var rnd = 0;
	while (true) { // try to generate a never used random number
		rnd = Math.random ();
		var s = "r" + rnd;
		if (ClazzLoader.usedRandoms[s] != rnd) {
			ClazzLoader.usedRandoms[s] = rnd;
			break;
		}
	}
	ClazzLoader.clazzTreeRoot.random = rnd;
	var node = ClazzLoader.clazzTreeRoot;
	return ClazzLoader.findNodeNextOptionalClass (node, status);
};

/* private */
/*-# findNodeNextOptionalClass -> fNNO #-*/
ClazzLoader.findNodeNextOptionalClass = function (node, status) {
	var rnd = ClazzLoader.clazzTreeRoot.random;
	// search musts first
	if (node.musts != null && node.musts.length != 0) {
		var n = ClazzLoader.searchClassArray (node.musts, rnd, status);
		if (n != null && (status != ClazzNode.STATUS_KNOWN 
				|| ClazzLoader.loadedScripts[n.path] != true)
				&& (status == ClazzNode.STATUS_DECLARED
				|| !ClazzLoader.isClassDefined (n.name))) {
			return n;
		}
	}
	// search optionals second
	if (node.optionals != null && node.optionals.length != 0) {
		var n = ClazzLoader.searchClassArray (node.optionals, rnd, status);
		if (n != null && (status != ClazzNode.STATUS_KNOWN 
				|| ClazzLoader.loadedScripts[n.path] != true)
				&& (status == ClazzNode.STATUS_DECLARED
				|| !ClazzLoader.isClassDefined (n.name))) {
			return n;
		}
	}
	// search itself
	if (node.status == status && (status != ClazzNode.STATUS_KNOWN 
			|| ClazzLoader.loadedScripts[node.path] != true)
			&& (status == ClazzNode.STATUS_DECLARED
			|| !ClazzLoader.isClassDefined (node.name))) {
		return node;
	}
	return null;
};

/* private */
ClazzLoader.searchClassArray = function (arr, rnd, status) {
	for (var i = 0; i < arr.length; i++) {
		var n = arr[i];
		if (n.status == status && (status != ClazzNode.STATUS_KNOWN 
				|| ClazzLoader.loadedScripts[n.path] != true)
				&& (status == ClazzNode.STATUS_DECLARED
				|| !ClazzLoader.isClassDefined (n.name))) {
			return n;
		} else {
			if (n.random == rnd) {
				continue;
			}
			n.random = rnd; // mark as visited!
			var nn = ClazzLoader.findNodeNextOptionalClass (n, status);
			if (nn != null) {
				return nn;
			}
		}
	}
	return null;
};

/**
 * This map variable is used to mark that *.js is correctly loaded.
 * In IE, ClazzLoader has defects to detect whether a *.js is correctly
 * loaded or not, so inner loading mark is used for detecting.
 */
/* private */
/*-# innerLoadedScripts -> ilss #-*/
ClazzLoader.innerLoadedScripts = {};

/**
 * This variable is used to keep the latest interactive SCRIPT element.
 */
/* private */
/*-# interactiveScript -> itst #-*/
ClazzLoader.interactiveScript = null;

/**
 * IE and Firefox/Mozilla are different in using <SCRIPT> tag to load *.js.
 */
/* private */
ClazzLoader.needOnloadCheck = false;

/**
 * Check the interactive status of SCRIPT elements to determine whether a
 * *.js file is correctly loaded or not.
 *
 * Only make senses for IE.
 */
/* protected */
ClazzLoader.checkInteractive = function () {
	//alert ("checking...");
	if (!ClazzLoader.needOnloadCheck) {
		return;
	}
	var is = ClazzLoader.interactiveScript;
	if (is != null && is.readyState == "interactive") { // IE
		return;
	}
	ClazzLoader.interactiveScript = null;
	var ss = document.getElementsByTagName ("SCRIPT");
	for (var i = 0; i < ss.length; i++) {
		if (ss[i].readyState == "interactive"
				&& ss[i].onreadystatechange != null) { // IE
			ClazzLoader.interactiveScript = ss[i];
			ClazzLoader.innerLoadedScripts[ss[i].src] = true;
		} else if (ClazzLoader.isOpera) { // Opera
			if (ss[i].readyState == "loaded" 
					&& ss[i].src != null && ss[i].src.length != 0) {
				ClazzLoader.innerLoadedScripts[ss[i].src] = true;
			}
		}
	}
};

/**
 * This method will be called in almost every *.js generated by Java2Script
 * compiler.
 */
/* protected */
ClazzLoader.load = function (musts, clazz, optionals, declaration) {

	ClazzLoader.checkInteractive ();

	if (clazz instanceof Array) {
		ClazzLoader.unwrapArray (clazz);
		for (var i = 0; i < clazz.length; i++) {
			ClazzLoader.load (musts, clazz[i], optionals, declaration, clazz);
		}
		return;
	}
	if (clazz.charAt (0) == '$') {
		clazz = "org.eclipse.s" + clazz.substring (1);
	}

	var node = ClazzLoader.mapPath2ClassNode["#" + clazz];

		//System.out.println("Loading " + clazz + " ... " + node);

	if (node == null) { // ClazzLoader.load called inside *.z.js?
		var n = ClazzLoader.findClass (clazz);
		if (n != null) {
			node = n;
		} else {
			node = new ClazzNode ();
		}
		node.name = clazz;
		var pp = ClazzLoader.classpathMap["#" + clazz];
		if (pp == null) { // TODO: Remove this test in final release
			//alert ("error finding classpathMap for " + clazz + " " );
			//error ("Java2Script implementation error! Please report this bug!");
			pp = "unknown"
		}
		//node.path = ClazzLoader.lastScriptPath;
		node.path = pp;
		//error ("..." + node.path + "//" + node.name);
		ClazzLoader.mappingPathNameNode (node.path, clazz, node);
		node.status = ClazzNode.STATUS_KNOWN;
		ClazzLoader.addChildClassNode(ClazzLoader.clazzTreeRoot, node, -1);
		//log (clazz);
		//alert ("[Java2Script] ClazzLoader#load is not executed correctly!");
		//*/
		/*
		if (declaration != null) {
			declaration ();
		}
		//alert ("[Java2Script] ClazzLoader#load is not executed correctly!");
		return;
		//*/
	}
	var okToInit = true;
	if (musts != null && musts.length != 0) {
		ClazzLoader.unwrapArray (musts);
		for (var i = 0; i < musts.length; i++) {
			var name = musts[i];
			if (name == null || name.length == 0) {
				continue;
			}
			//System.out.println(node.name + " must have " + name)
			if (ClazzLoader.isClassDefined (name)
					|| ClazzLoader.isClassExcluded (name)) {
			//System.out.println("which it does")
				continue;
			}
			okToInit = false;
			var n = ClazzLoader.findClass (name);
			if (n == null) {
				//System.out.println(clazz + " requires " + name)
				n = new ClazzNode ();
				n.name = musts[i];
				n.status = ClazzNode.STATUS_KNOWN;
			}
			n.requiredBy = node;
			ClazzLoader.addChildClassNode (node, n, 1);
		}
	}

	/*
	 * The following lines are commented intentionally.
	 * So lots of class is not declared until there is a must?
	 *
	 * TODO: Test the commented won't break up the dependency tree.
	 */
	/*
	if (okToInit) {
		declaration ();
		node.declaration = null;
		node.status = ClazzNode.STATUS_DECLARED;
	} else {
		node.declaration = declaration;
	}
	*/
	if (arguments.length == 5 && declaration != null) {
		declaration.status = node.status;
		declaration.clazzList = arguments[4];
	}
	node.declaration = declaration;
	if (declaration != null) {
		//System.out.println("declaration found for " + node.name + " " + node.path)
		node.status = ClazzNode.STATUS_CONTENT_LOADED;
	}

	var isOptionalsOK = true;
	if (optionals != null && optionals.length != 0) {
		ClazzLoader.unwrapArray (optionals);
		for (var i = 0; i < optionals.length; i++) {
			var name = optionals[i];
			if (name == null || name.length == 0) {
				continue;
			}
			if (ClazzLoader.isClassDefined (name) 
					|| ClazzLoader.isClassExcluded (name)) {
				continue;
			}
			isOptionalsOK = false;
			var n = ClazzLoader.findClass (name);
			if (n == null) {
				n = new ClazzNode ();
				n.name = optionals[i];
				n.status = ClazzNode.STATUS_KNOWN;
			}
			ClazzLoader.addChildClassNode (node, n, -1);
		}
	}
};

/*
 * Try to be compatiable of Clazz
 */
if (window["Clazz"] != null) {
	Clazz.load = ClazzLoader.load;
	//if (window["$_L"] != null) {
	//	$_L = Clazz.load;
	//}
}

/**
 *
 */
/* protected */
/*-# findClass -> fC #-*/
ClazzLoader.findClass = function (clazzName) {
	var rnd = 0;
	while (true) { // try to generate a never used random number
		rnd = Math.random ();
		var s = "r" + rnd;
		if (ClazzLoader.usedRandoms[s] != rnd) {
			ClazzLoader.usedRandoms[s] = rnd;
			break;
		}
	}
	ClazzLoader.clazzTreeRoot.random = rnd;
	return ClazzLoader.findClassUnderNode (clazzName, 
			ClazzLoader.clazzTreeRoot);
};

/* private */
/*-# findClassUnderNode -> fCU #-*/
ClazzLoader.findClassUnderNode = function (clazzName, node) {
	var rnd = ClazzLoader.clazzTreeRoot.random;
	if (node.name == clazzName) {
		return node;
	}
	// musts first
	for (var i = 0; i < node.musts.length; i++) {
		var n = node.musts[i];
		if (n.name == clazzName) {
			return n;
		}
		if (n.random == rnd) {
			continue;
		}
		n.random = rnd;
		var nn = ClazzLoader.findClassUnderNode (clazzName, n);
		if (nn != null) {
			return nn;
		}
	}
	// optionals last
	for (var i = 0; i < node.optionals.length; i++) {
		var n = node.optionals[i];
		if (n.name == clazzName) {
			return n;
		}
		if (n.random == rnd) {
			continue;
		}
		n.random = rnd;
		var nn = ClazzLoader.findClassUnderNode (clazzName, n);
		if (nn != null) {
			return nn;
		}
	}
	return null;
};

/**
 * Map different class to the same path! Many classes may be packed into
 * a *.z.js already.
 *
 * @path *.js path
 * @name class name
 * @node ClazzNode object
 */
/* private */
/*-# mappingPathNameNode -> mpp #-*/
ClazzLoader.mappingPathNameNode = function (path, name, node) {
	var map = ClazzLoader.mapPath2ClassNode;
	var keyPath = "@" + path;
	var v = map[keyPath];
	if (v != null) {
		if (v instanceof Array) {
			var existed = false;
			for (var i = 0; i < v.length; i++) {
				if (v[i].name == name) {
					existed = true;
					break;
				}
			}
			if (!existed) {
				v[v.length] = node;
			}
		} else {
			var arr = new Array ();
			arr[0] = v;
			arr[1] = node;
			map[keyPath] = arr;
		}
	} else {
		map[keyPath] = node;
	}
	map["#" + name] = node;
};

/* protected */
/*-# loadClassNode -> lCN #-*/
ClazzLoader.loadClassNode = function (node) {
	var name = node.name;
	if (!ClazzLoader.isClassDefined (name) 
			&& !ClazzLoader.isClassExcluded (name)) {
		var path = ClazzLoader.getClasspathFor (name/*, true*/);
		node.path = path;
		ClazzLoader.mappingPathNameNode (path, name, node);
		if (!ClazzLoader.loadedScripts[path]) {
			ClazzLoader.loadScript (path, node.requiredBy);
			return true;
		}
	}
	return false;
};


/* protected */
ClazzLoader.runtimeKeyClass = "java.lang.String";

/**
 * Queue used to store classes before key class is loaded.
 */
/* private */
ClazzLoader.queueBe4KeyClazz = new Array ();

/**
 * Return J2SLib base path from existed SCRIPT src attribute.
 */
/* private */
/*-# getJ2SLibBase -> gLB #-*/
ClazzLoader.getJ2SLibBase = function () {
	var o = window["j2s.lib"];
	if (o != null) {
	 if (o.base == null)
			o.base = "http://archive.java2script.org/";
		return o.base + (o.alias == "." ? "" : (o.alias ? o.alias : (o.version ? o.version : "1.0.0")) + "/");
	}
	var ss = document.getElementsByTagName ("SCRIPT");
	for (var i = 0; i < ss.length; i++) {
		var src = ss[i].src;
		var idx = src.indexOf ("j2slib.z.js"); // consider it as key string!
		if (idx != -1) {
			return src.substring (0, idx);
		}
		idx = src.indexOf ("j2slibcore.z.js"); // consider it as key string!
		if (idx != -1) {
			return src.substring (0, idx);
		}
		var base = ClazzLoader.classpathMap["@java"];
		if (base != null) {
			return base;
		}
		idx = src.indexOf ("java/lang/ClassLoader.js"); // may be not packed yet
		if (idx != -1) {
			return src.substring (0, idx);
		}
	}
	return null;
};

/* private static */
/*-# J2SLibBase -> JLB #-*/
ClazzLoader.J2SLibBase = null;
/*-# fastGetJ2SLibBase -> fgLB #-*/
ClazzLoader.fastGetJ2SLibBase = function () {
	if (ClazzLoader.J2SLibBase == null) {
		ClazzLoader.J2SLibBase = ClazzLoader.getJ2SLibBase ();
	}
	return ClazzLoader.J2SLibBase;
};

/*
 * Check whether given package's classpath is setup or not.
 * Only "java" and "org.eclipse.swt" are accepted in argument.
 */
/* private */
/*-# assurePackageClasspath -> acp #-*/
ClazzLoader.assurePackageClasspath = function (pkg) {
	var r = window[pkg + ".registered"];
	if (r != false && r != true && ClazzLoader.classpathMap["@" + pkg] == null) {
		window[pkg + ".registered"] = false;
		var base = ClazzLoader.fastGetJ2SLibBase ();
		if (base == null) {
			base = "j2s"; 
		}
		ClazzLoader.packageClasspath (pkg, base, true);
	}
};

/**
 * Load the given class ant its related classes.
 */
/* public */
ClazzLoader.loadClass = function (name, optionalsLoaded, forced, async) {

 	if (typeof optionalsLoaded == "boolean") {
		return Clazz.evalType (name);
	}

	/*
	 * Make sure that
	 * ClazzLoader.packageClasspath ("java", base, true); 
	 * is called before any ClazzLoader#loadClass is called.
	 */
	ClazzLoader.assurePackageClasspath ("java");
	ClazzLoader.assurePackageClasspath ("core"); // JSmol uses j2s/core, not j2s/java for package.js

	var swtPkg = "org.eclipse.swt";
	if (name.indexOf (swtPkg) == 0 || name.indexOf ("$wt") == 0) {
		ClazzLoader.assurePackageClasspath (swtPkg);
	}
	if (name.indexOf ("junit") == 0) {
		ClazzLoader.assurePackageClasspath ("junit");
	}

	/*
	 * Any ClazzLoader#loadClass calls will be queued until java.* core classes
	 * is loaded.
	 */
	ClazzLoader.keepOnLoading = true;
	if (!forced && (
			 ClazzLoader.pkgRefCount != 0 && name.lastIndexOf (".package") != name.length - 8
		|| !ClazzLoader.isClassDefined (ClazzLoader.runtimeKeyClass) && name.indexOf ("java.") != 0
	)) {
		var qbs = ClazzLoader.queueBe4KeyClazz;
		//alert("keep on loading " + name)
		qbs[qbs.length] = [name, optionalsLoaded];

		return;
	}
	if (!ClazzLoader.isClassDefined (name) 
			&& !ClazzLoader.isClassExcluded (name)) {
		var path = ClazzLoader.getClasspathFor (name/*, true*/);
		var existed = ClazzLoader.loadedScripts[path];
		var qq = ClazzLoader.classQueue;
		if (!existed) {
			for (var i = qq.length - 1; i >= 0; i--) {
				if (qq[i].path == path || qq[i].name == name) {
					existed = true;
				}
			}
		}
		//alert("@#@#@#@# " + name);
		if (!existed) {
			var n = null;
			if (Clazz.unloadedClasses[name] != null) {
				n = ClazzLoader.findClass (name);
			}
			if (n == null) {
				n = new ClazzNode ();
			}
			n.name = name;
			n.path = path;
			ClazzLoader.mappingPathNameNode (path, name, n);
			n.optionalsLoaded = optionalsLoaded;
			n.status = ClazzNode.STATUS_KNOWN;
			/*-# needBeingQueued -> nQ #-*/
			var needBeingQueued = false;
			//error (qq.length + ":" + qq);
			//error (path);
			for (var i = qq.length - 1; i >= 0; i--) {
				if (qq[i].status != ClazzNode.STATUS_OPTIONALS_LOADED) {
					needBeingQueued = true;
					break;
				}
			}
			if (path.lastIndexOf ("package.js") == path.length - 10) {//forced
				// push class to queue
				var inserted = false;
				for (var i = qq.length - 1; i >= 0; i--) {
					var name = qq[i].name;
					if (name.lastIndexOf ("package.js") == name.length - 10) {
						qq[i + 1] = n;
						inserted = true;
						break;
					}
					qq[i + 1] = qq[i];
				}
				if (!inserted) {
					qq[0] = n;
				}
			} else if (needBeingQueued) {
				qq[qq.length] = n;
			}
//alert(["-------------"]);
			if (!needBeingQueued) { // can be loaded directly
				/*-# bakEntryClassLoading -> bkECL #-*/
				var bakEntryClassLoading = false;
				if (optionalsLoaded != null) {	
					bakEntryClassLoading = ClazzLoader.isLoadingEntryClass;
					ClazzLoader.isLoadingEntryClass = true;
				}
				ClazzLoader.addChildClassNode(ClazzLoader.clazzTreeRoot, n, 1);
				//System.out.println("Clazz loading " + n.path)
				ClazzLoader.loadScript (n.path);
				if (optionalsLoaded != null) {
					ClazzLoader.isLoadingEntryClass = bakEntryClassLoading;
				}
			}
		} else if (optionalsLoaded != null) {
			var n = ClazzLoader.findClass (name);
			if (n != null) {
				if (n.optionalsLoaded == null) {
					n.optionalsLoaded = optionalsLoaded;
				} else if (optionalsLoaded != n.optionalsLoaded) {
					n.optionalsLoaded = (function (oF, nF) {
						return function () {
							oF();
							nF();
						};
					}) (n.optionalsLoaded, optionalsLoaded);
				}
			}
		}
	} else if (optionalsLoaded != null && ClazzLoader.isClassDefined (name)) {

		var nn = ClazzLoader.findClass (name);

		if (nn == null || nn.status >= ClazzNode.STATUS_OPTIONALS_LOADED) {
			/*if (nn != null) {
				ClazzLoader.destroyClassNode (nn);
			}*/
			if (async) {
				window.setTimeout (optionalsLoaded, 25);
			} else {
				optionalsLoaded ();
			}
		} // else ... should be called later
	}
};

/**
 * Load the application by the given class name and run its static main method.
 */
/* public */
ClazzLoader.loadJ2SApp = function (clazz, args, loaded) {
	if (clazz == null) {
		return;
	}
	var clazzStr = clazz;
	if (clazz.charAt (0) == '$') {
		clazzStr = "org.eclipse.s" + clazz.substring (1);
	}
	var idx = -1;
	if ((idx = clazzStr.indexOf ("@")) != -1) {
		var path = clazzStr.substring (idx + 1);
		ClazzLoader.setPrimaryFolder (path); // BH JSmol will not use this
		clazzStr = clazzStr.substring (0, idx);
		idx = clazzStr.lastIndexOf (".");
		if (idx != -1) {
			var pkgName = clazzStr.substring (0, idx);
			ClazzLoader.packageClasspath (pkgName, path);
		}
	}
	var agmts = args;
	if (agmts == null || !(agmts instanceof Array)) {
		agmts = [];
	}
	var afterLoaded = loaded;
	if (afterLoaded == null) {
		afterLoaded = (function (clazzName, argv) {
			return function () {
				Clazz.evalType (clazzName).main (argv);
			};
		}) (clazzStr, agmts);
	} else {
		afterLoaded = loaded (clazzStr, agmts);
	}
	ClazzLoader.loadClass (clazzStr, afterLoaded);
};
/**
 * Load JUnit tests by the given class name.
 */
/* public */
ClazzLoader.loadJUnit = function (clazz, args) {
	var afterLoaded = function (clazzName, argv) {
		return function () {
			ClazzLoader.loadClass ("junit.textui.TestRunner", function () {
				junit.textui.TestRunner.run (Clazz.evalType (clazzName));
			});
		};
	};
	ClazzLoader.loadJ2SApp (clazz, args, afterLoaded);
};

/* private */
ClazzLoader.runtimeLoaded = function () {
	if (ClazzLoader.pkgRefCount != 0 
			|| !ClazzLoader.isClassDefined (ClazzLoader.runtimeKeyClass)) {
		return;
	}
	var qbs = ClazzLoader.queueBe4KeyClazz;
	for (var i = 0; i < qbs.length; i++) {
		ClazzLoader.loadClass (qbs[i][0], qbs[i][1]);
	}
	ClazzLoader.queueBe4KeyClazz = [];
	/*
	 * Should not set to empty function! Some later package may need this
	 * runtimeLoaded function. For example, lazily loading SWT package may
	 * require this runtimeLoaded function. 
	 * -- zhou renjian @ Dec 17, 2006
	 */
	// ClazzLoader.runtimeLoaded = function () {};
};

/*
 * Load those key *.z.js. This *.z.js will be surely loaded before other 
 * queued *.js.
 */
/* public */
ClazzLoader.loadZJar = function (zjarPath, keyClazz) {
	var keyClass = keyClazz;
	var isArr = (keyClazz instanceof Array);
	if (isArr) {
		keyClass = keyClazz[keyClazz.length - 1];
	}
	ClazzLoader.jarClasspath (zjarPath, isArr ? keyClazz
			: [keyClazz]);
	if (keyClazz == ClazzLoader.runtimeKeyClass) {
		ClazzLoader.loadClass (keyClass, ClazzLoader.runtimeLoaded, true);
	} else {
		ClazzLoader.loadClass (keyClass, null, true);
	}
};

ClazzLoader._nodeMap = {};
ClazzLoader._allNodes = [];
/**
 * The method help constructing the multiple-binary class dependency tree.
 */
/* private */
/*-# addChildClassNode -> addCCN #-*/
ClazzLoader.addChildClassNode = function (parent, child, type) {
	var existed = false;
	var arr = null;
	if (type == 1) {
		arr = parent.musts;
		if (!child.requiredBy)child.requiredBy = parent;
		if (!parent.requires){parent.requires = [];parent.requiresMap = {};}
		if (!parent.requiresMap[child.name]) {
			parent.requiresMap[child.name] = child;
			parent.requires.push[child];
		}
	} else {
		arr = parent.optionals;
	}
	if (!ClazzLoader._nodeMap[child.name]) {
		ClazzLoader._allNodes.push(child)
		ClazzLoader._nodeMap[child.name]=child
	}
	for (var i = 0; i < arr.length; i++) {
		if (arr[i].name == child.name) {
			existed = true;
			break;
		}
	}
	if (!existed) {
		/*
		if (type != 1) { // test cyclic optionals
			existed = false;
			for (var j = 0; j < child.optionals.length; j++) {
				if (child.optionals[j].name == parent.name) {
					existed = true;
					break;
				}
			}
		}
		*/
		/* above existed tests are commented */
		//if (!existed) {
			arr[arr.length] = child;
			var swtPkg = "org.eclipse.swt";
			if (child.name.indexOf (swtPkg) == 0 
					|| child.name.indexOf ("$wt") == 0) {
				window["swt.lazy.loading.callback"] = ClazzLoader.swtLazyLoading;
				ClazzLoader.assurePackageClasspath (swtPkg);
			}
			if (ClazzLoader.isLoadingEntryClass 
					&& child.name.indexOf ("java") != 0 
					&& child.name.indexOf ("net.sf.j2s.ajax") != 0) {
				if (ClazzLoader.besidesJavaPackage) {
					ClazzLoader.isLoadingEntryClass = false;
				}
				ClazzLoader.besidesJavaPackage = true;
			}
		//}
	}
	existed = false;
	for (var i = 0; i < child.parents.length; i++) {
		if (child.parents[i].name == parent.name) {
			existed = true;
			break;
		}
	}
	if (!existed && parent.name != null && parent != ClazzLoader.clazzTreeRoot
			&& parent != child) {
		child.parents[child.parents.length] = parent;
	}
};

/*
 * Some SWT classes may already skip ClazzLoader#tryToLoadNext when it's
 * detected that SWT is in lazy loading mode. Here it will try to re-execute
 * those ClazzLoader#tryToLoadNext
 */
/* private */
ClazzLoader.swtLazyLoading = function () {
	ClazzLoader.lockQueueBe4SWT = false;
	var qbs = ClazzLoader.queueBe4SWT;
	for (var i = 0; i < qbs.length; i++) {
		ClazzLoader.tryToLoadNext (qbs[i]);
	}
	ClazzLoader.queueBe4SWT  = [];
};

/* private */
ClazzLoader.removeFromArray = function (node, arr) {
	if (arr == null || node == null) {
		return false;
	}
	/*var isPackedJS = (node.path != null
			&& node.path.indexOf (".z.js") == node.path.length - 5);
	log ("... remove " + node.path + " :: " + isPackedJS);*/
	var j = 0;
	for (var i = 0; i < arr.length; i++) {
		if (!(arr[i] === node/* || (isPackedJS && arr[i].path == node.path)*/)) {
			if (j < i) {
				arr[j] = arr[i];
			}
			j++;
		}
	}
	arr.length = j;
	return false;
};

/* private */
/*-# destroyClassNode -> dCN #-*/
ClazzLoader.destroyClassNode = function (node) {
	//log (node.name + " // " + node.path);
	var parents = node.parents;
	if (parents != null) {
		for (var k = 0; k < parents.length; k++) {
			if (!ClazzLoader.removeFromArray (node, parents[k].musts)) {
				ClazzLoader.removeFromArray (node, parents[k].optionals);
			}
		}
	}
	/*
	if (node.optionalsLoaded != null) {
		node.optionalsLoaded ();
	}
	if (!ClazzLoader.removeFromArray (node, ClazzLoader.clazzTreeRoot.musts)) {
		ClazzLoader.removeFromArray (node, ClazzLoader.clazzTreeRoot.optionals);
	}
	*/
};

/* For hotspot and unloading */

/* protected */
ClazzLoader.unloadClassExt = function (qClazzName) {
	if (ClazzLoader.definedClasses != null) {
		ClazzLoader.definedClasses[qClazzName] = false;
	}
	if (ClazzLoader.classpathMap["#" + qClazzName] != null) {
		var pp = ClazzLoader.classpathMap["#" + qClazzName];
		ClazzLoader.classpathMap["#" + qClazzName] = null;
		var arr = ClazzLoader.classpathMap["$" + pp];
		var removed = false;
		for (var i = 0; i < arr.length; i++) {
			if (arr[i] == qClazzName) {
				for (var j = i; j < arr.length - 1; j++) {
					arr[j] = arr[j + 1];
				}
				arr.length--;
				removed = true;
				break;
			}
		}
		if (removed) {
			ClazzLoader.classpathMap["$" + pp] = arr;
		}
	}
	var n = ClazzLoader.findClass (qClazzName);
	if (n != null) {
		n.status = ClazzNode.STATUS_KNOWN;
		ClazzLoader.loadedScripts[n.path] = false;
	}
	var path = ClazzLoader.getClasspathFor (qClazzName);
	ClazzLoader.loadedScripts[path] = false;
	if (ClazzLoader.innerLoadedScripts[path]) {
		ClazzLoader.innerLoadedScripts[path] = false;
	}

	ClazzLoader.classUnloaded (qClazzName);
};

/* protected */ /* Clazz#assureInnerClass */
ClazzLoader.assureInnerClass = function (clzz, fun) {
	var clzzName = clzz.__CLASS_NAME__;
	if (Clazz.unloadedClasses[clzzName]) {
		if (clzzName.indexOf ("$") != -1) return;
		var list = new Array ();
		var key = clzzName + "$";
		for (var s in Clazz.unloadedClasses) {
			if (Clazz.unloadedClasses[s] != null && s.indexOf (key) == 0) {
				list[list.length] = s;
			}
		}
		if (list.length == 0) return;
		var funStr = "" + fun;
		var idx1 = funStr.indexOf (key);
		if (idx1 == -1) return;
		var idx2 = funStr.indexOf ("\"", idx1 + key.length);
		if (idx2 == -1) return; // idx2 should never be -1;
		var anonyClazz = funStr.substring (idx1, idx2);
		if (Clazz.unloadedClasses[anonyClazz] == null) return;
		var idx3 = funStr.indexOf ("{", idx2);
		if (idx3 == -1) return;
		idx3++;
		var idx4 = funStr.indexOf ("(" + anonyClazz + ",", idx3 + 3);
		if (idx4 == -1) return; // idx3 should never be -1;
		var idx5 = funStr.lastIndexOf ("}", idx4 - 1);
		if (idx5 == -1) return;
		var innerClazzStr = funStr.substring (idx3, idx5);
		eval (innerClazzStr);
		Clazz.unloadedClasses[anonyClazz] = null;
		/*
		window.setTimeout ((function (clz, str) {
			return function () {
				eval (str);
				Clazz.unloadedClasses[clz] = null;
			};
		}) (anonyClazz, innerClazzStr), 10);
		*/
	}
};


//ClassLoader = ClazzLoader;

})(Clazz, Clazz._Loader, Clazz._Node);

//}
/******************************************************************************
 * Copyright (c) 2007 java2script.org and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Zhou Renjian - initial API and implementation
 *****************************************************************************/
/*******
 * @author zhou renjian
 * @create Jan 11, 2007
 *******/

Clazz._LoaderProgressMonitor = function() {};

;(function(clpm) {
clpm.fadeOutTimer = null;
clpm.fadeAlpha = 0;
clpm.monitorEl = null;
clpm.lastScrollTop = 0;
clpm.bindingParent = null;
clpm.DEFAULT_OPACITY = (self.Jmol && Jmol._j2sLoadMonitorOpacity ? Jmol._j2sLoadMonitorOpacity : 55);
/* private static */ clpm.clearChildren = function (el) {
	if (el == null) return;
	for (var i = el.childNodes.length - 1; i >= 0; i--) {
		var child = el.childNodes[i];
		if (child == null) continue;
		if (child.childNodes != null && child.childNodes.length != 0) {
			this.clearChildren (child);
		}
		try {
			el.removeChild (child);
		} catch (e) {};
	}
};
/* private */ clpm.setAlpha = function (alpha) {
	if (this.fadeOutTimer != null && alpha == this.DEFAULT_OPACITY) {
		window.clearTimeout (this.fadeOutTimer);
		this.fadeOutTimer = null;
	}
	this.fadeAlpha = alpha;
	var ua = navigator.userAgent.toLowerCase ();
	this.monitorEl.style.filter = "Alpha(Opacity=" + alpha + ")";
	this.monitorEl.style.opacity = alpha / 100.0;
};
/* private */ clpm.hiddingOnMouseOver = function () {
	this.style.display = "none";
};
/* private */ clpm.attached = false;
/* private */ clpm.cleanup = function () {
	var oThis = clpm;
	//if (oThis.monitorEl != null) {
	//	oThis.monitorEl.onmouseover = null;
	//}
	oThis.monitorEl = null;
	oThis.bindingParent = null;
	//Clazz.removeEvent (window, "unload", oThis.cleanup);
	//window.detachEvent ("onunload", oThis.cleanup);
	oThis.attached = false;
};
/* private */ clpm.createHandle = function () {
	var div = document.createElement ("DIV");
	div.id = "clazzloader-status";
	div.style.cssText = "position:absolute;bottom:4px;left:4px;padding:2px 8px;"
			+ "z-index:" + (window["j2s.lib"].monitorZIndex || 10000) + ";background-color:#8e0000;color:yellow;" 
			+ "font-family:Arial, sans-serif;font-size:10pt;white-space:nowrap;";
	div.onmouseover = this.hiddingOnMouseOver;
	this.monitorEl = div;
	if (this.bindingParent == null) {
		document.body.appendChild (div);
	} else {
		this.bindingParent.appendChild (div);
	}
	return div;
};
/* private */ clpm.fadeOut = function () {
	if (this.monitorEl.style.display == "none") return;
	if (this.fadeAlpha == this.DEFAULT_OPACITY) {
		this.fadeOutTimer = window.setTimeout (function () {
					clpm.fadeOut ();
				}, 750);
		this.fadeAlpha -= 5;
	} else if (this.fadeAlpha - 10 >= 0) {
		this.setAlpha (this.fadeAlpha - 10);
		this.fadeOutTimer = window.setTimeout (function () {
					clpm.fadeOut ();
				}, 40);
	} else {
		this.monitorEl.style.display = "none";
	}
};
/* private */ clpm.getFixedOffsetTop = function (){
	if (this.bindingParent != null) {
		var b = this.bindingParent;
		return b.scrollTop;
	}
	var dua = navigator.userAgent;
	var b = document.body;
	var p = b.parentNode;
	var pcHeight = p.clientHeight;
	var bcScrollTop = b.scrollTop + b.offsetTop;
	var pcScrollTop = p.scrollTop + p.offsetTop;
	if (dua.indexOf("Opera") == -1 && document.all) {
		return (pcHeight == 0) ? bcScrollTop : pcScrollTop;
	} else if (dua.indexOf("Gecko") != -1) {
		return (pcHeight == p.offsetHeight 
				&& pcHeight == p.scrollHeight) ? bcScrollTop : pcScrollTop;
	}
	return bcScrollTop;
};
/* public */
/*clpm.initialize = function (parent) {
alert("clpm.innit")
	this.bindingParent = parent;
	if (parent != null && !this.attached) {
		this.attached = true;
		//Clazz.addEvent (window, "unload", this.cleanup);
		// window.attachEvent ("onunload", this.cleanup);
	}
};
*/
/* public */
clpm.showStatus = function (msg, fading) {
	if (this.monitorEl == null) {
		this.createHandle ();
		if (!this.attached) {
			this.attached = true;
			//Clazz.addEvent (window, "unload", this.cleanup);
			// window.attachEvent ("onunload", this.cleanup);
		}
	}
	this.clearChildren (this.monitorEl);
	this.monitorEl.appendChild (document.createTextNode ("" + msg));
	if (this.monitorEl.style.display == "none") {
		this.monitorEl.style.display = "";
	}
	this.setAlpha (this.DEFAULT_OPACITY);
	var offTop = this.getFixedOffsetTop ();
	if (this.lastScrollTop != offTop) {
		this.lastScrollTop = offTop;
		this.monitorEl.style.bottom = (this.lastScrollTop + 4) + "px";
	}
	if (fading) {
		this.fadeOut();
	}
};

if (window["ClazzLoader"] != null) {
	ClazzLoader.scriptLoading = function (file) {
		clpm.showStatus ("Loading " + file + "...");
	};
	ClazzLoader.scriptLoaded = function (file) {
		clpm.showStatus (file + " loaded.", true);
	};
	ClazzLoader.globalLoaded = function (file) {
		clpm.showStatus ("Application loaded.", true);
	};
	ClazzLoader.classUnloaded = function (clazz) {
		clpm.showStatus ("Class " + clazz + " is unloaded.", true);
	};
	ClazzLoader.classReloaded = function (clazz) {
		clpm.showStatus ("Class " + clazz + " is reloaded.", true);
	};

	var ua = navigator.userAgent.toLowerCase ();
	if (ua.indexOf ("msie") != -1 && ua.indexOf ("opera") == -1) {
		ClazzLoader.setLoadingMode ("script", 5);
	}
}

})(Clazz._LoaderProgressMonitor);

//}
/******************************************************************************
 * Copyright (c) 2007 java2script.org and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Zhou Renjian - initial API and implementation
 *****************************************************************************/
/*******
 * @author zhou renjian
 * @create Nov 5, 2005
 *******/

(function(Console) {
/**
 * Setting maxTotalLines to -1 will not limit the console result
 */
/* protected */
/*-# maxTotalLines -> mtl #-*/
Console.maxTotalLines =	10000;

/* protected */
Console.setMaxTotalLines = function (lines) {
	if (lines <= 0) {
		Console.maxTotalLines = 999999; // Won't reach before browser cracks
	} else {
		Console.maxTotalLines = lines;
	}
};

/*
 * The console window will be flicking badly in some situation. Try to use
 * double buffer to avoid flicking.
 */
/* protected */
/*-# buffering -> bi #-*/
//Console.buffering = false;

/* protected */
//Console.enableBuffering = function (enabled) {
//	Console.buffering = enabled;
//};

/* protected */
/*-# maxBufferedLines -> mbl #-*/
//Console.maxBufferedLines = 20;

/* protected */
//Console.setMaxBufferedLines = function (lines) {
//	if (lines <= 0) {
//		Console.maxBufferedLines = 20;
//	} else {
//		Console.maxBufferedLines = lines;
//	}
//};

/* protected */
/*-# maxLatency -> mlc #-*/
Console.maxLatency = 40;

/* protected */
Console.setMaxLatency = function (latency) {
	if (latency <= 0) {
		Console.maxLatency = 40;
	} else {
		Console.maxLatency = latency;
	}
};

/* protected */
/*-# pinning -> pi #-*/
Console.pinning  = false;

/* protected */
Console.enablePinning = function (enabled) {
	Console.pinning = enabled;
};

/* private */
/*-# linesCount -> lc #-*/
Console.linesCount = 0;

/* private */
/*-# metLineBreak -> mbr #-*/
Console.metLineBreak = false;

/* private */
Console.splitNeedFixed = "\n".split (/\n/).length != 2; // IE

/**
 * For IE to get a correct split result.
 */
/* private */
/*-# splitIntoLineByR -> slr #-*/
Console.splitIntoLineByR = function (s) {
	var arr = new Array ();
	var i = 0;
	var last = -1;
	while (true) {
		i = s.indexOf ('\r', last + 1);
		if (i != -1) {
			arr[arr.length] = s.substring (last + 1, i);
			last = i;
			if (last + 1 == s.length) {
				arr[arr.length] = "";
				break;
			}
		} else {
			arr[arr.length] = s.substring (last + 1);
			break;
		}
	}
	return arr;
};

/**
 * For IE to get a correct split result.
 */
/* private */
/*-# splitIntoLines -> sil #-*/
Console.splitIntoLines = function (s) {
	var arr = new Array ();
	if (s == null) {
		return arr;
	}
	var i = 0;
	var last = -1;
	while (true) {
		i = s.indexOf ('\n', last + 1);
		var str = null;
		if (i != -1) {
			if (i > 0 && s.charAt (i - 1) == '\r') {
				str = s.substring (last + 1, i - 1);
			} else {
				str = s.substring (last + 1, i);
			}
			last = i;
		} else {
			str = s.substring (last + 1);
		}
		var rArr = Console.splitIntoLineByR (str);
		for (var k = 0; k < rArr.length; k++) {
			arr[arr.length] = rArr[k];
		}
		if (i == -1) {
			break;
		} else if (last + 1 == s.length) {
			arr[arr.length] = "";
			break;
		}
	}
	return arr;
};

/**
 * Cache the console until the document.body is ready and console element
 * is created.
 */
/* private */
/*-# consoleBuffer -> cB #-*/
//Console.consoleBuffer = new Array ();

/* private */
/*-# lastOutputTime -> oT #-*/
//Console.lastOutputTime = new Date ().getTime ();

/* private */
/*-# checkingTimer -> lct #-*/
//Console.checkingTimer = 0;

/* private */
/*-# loopChecking -> li  #-*/
//Console.loopChecking = function () {
//	if (Console.consoleBuffer.length == 0) {
//		return ;
//	}
//	var console = document.getElementById ("_console_");;
//	if (console == null) {
//		if (document.body == null) {
//			if (Console.checkingTimer == 0) {
//				Console.checkingTimer = window.setTimeout (
//						"Console.loopChecking ();", Console.maxLatency);
//			}
//			return ;
//		}
//	}
//	Console.consoleOutput ();
//};

/*
 * Give an extension point so external script can create and bind the console
 * themself.
 *
 * TODO: provide more template of binding console window to browser.
 */
/* protected */
Console.createConsoleWindow = function (parentEl) {
	var console = document.createElement ("DIV");
	console.style.cssText = "font-family:monospace, Arial, sans-serif;";
	document.body.appendChild (console);
	return console;
};

/* protected */
/*-#
 # consoleOutput -> cot 
 #
 #-*/
Console.consoleOutput = function (s, color) {
	var o = window["j2s.lib"];
	var console = (o && o.console);
	if (console && typeof console == "string")
		console = document.getElementById(console)
		// console)
	if (!console)
		return false; // BH this just means we have turned off all console action
	if (Console.linesCount > Console.maxTotalLines) {
		for (var i = 0; i < Console.linesCount - Console.maxTotalLines; i++) {
			if (console != null && console.childNodes.length > 0) {
				console.removeChild (console.childNodes[0]);
			}
		}
		Console.linesCount = Console.maxTotalLines;
	}

	/*-# willMeetLineBreak -> wbr #-*/
	var willMeetLineBreak = false;
	if (typeof s == "undefined") {
		s = "";
	} else if (s == null) {
		s = "null";
	} else {
		s = "" + s;
	}
	if (s.length > 0) {
		/*-# lastChar -> lc #-*/
		var lastChar = s.charAt (s.length - 1);
		if (lastChar == '\n') {
			if (s.length > 1) {
				var preLastChar = s.charAt (s.length - 2);
				if (preLastChar == '\r') {
					s = s.substring (0, s.length - 2);
				} else {
					s = s.substring (0, s.length - 1);
				}
			} else {
				s = "";
			}
			willMeetLineBreak = true;
		} else if (lastChar == '\r') {
			s = s.substring (0, s.length - 1);
			willMeetLineBreak = true;
		}
	}

	var lines = null;
	s = s.replace (/\t/g, Console.c160);
	if (Console.splitNeedFixed) { // IE
		try {
			lines = Console.splitIntoLines (s);
		} catch (e) {
			window.popup (e.message);
		}
	} else { // Mozilla/Firefox, Opera
		lines = s.split (/\r\n|\r|\n/g);
	}
	for (var i = 0; i < lines.length; i++) {
		/*-# lastLineEl -> lE #-*/
		var lastLineEl = null;
		if (Console.metLineBreak || Console.linesCount == 0 
				|| console.childNodes.length < 1) {
			lastLineEl = document.createElement ("DIV");
			console.appendChild (lastLineEl);
			lastLineEl.style.whiteSpace = "nowrap";
			Console.linesCount++;
		} else {
			try {
				lastLineEl = console.childNodes[console.childNodes.length - 1];
			} catch (e) {
				lastLineEl = document.createElement ("DIV");
				console.appendChild (lastLineEl);
				lastLineEl.style.whiteSpace = "nowrap";
				Console.linesCount++;
			}
		}
		var el = document.createElement ("SPAN");
		lastLineEl.appendChild (el);
		el.style.whiteSpace = "nowrap";
		if (color != null) {
			el.style.color = color;
		}
		if (lines[i].length == 0) {
			lines[i] = String.fromCharCode (160);
			//el.style.height = "1em";
		}
		el.appendChild (document.createTextNode (lines[i]));
		if (!Console.pinning) {
			console.scrollTop += 100;
		}

		if (i != lines.length - 1) {
			Console.metLineBreak = true;
		} else {
			Console.metLineBreak = willMeetLineBreak;
		}
	}

	var cssClazzName = console.parentNode.className;
	if (!Console.pinning && cssClazzName != null
			&& cssClazzName.indexOf ("composite") != -1) {
		console.parentNode.scrollTop = console.parentNode.scrollHeight;
	}
	Console.lastOutputTime = new Date ().getTime ();
};

/*
 * Clear all contents inside the console.
 */
/* public */
Console.clear = function () {
	try {
		Console.metLineBreak = true;
		var o = window["j2s.lib"];
		var console = o && o.console;
		if (!console || !(console = document.getElementById (console)))
			return;
		var childNodes = console.childNodes;
		for (var i = childNodes.length - 1; i >= 0; i--)
			console.removeChild (childNodes[i]);
		Console.linesCount = 0;
	} catch(e){};
};

/* public */
Clazz.alert = function (s) {
	Console.consoleOutput (s + "\r\n");
};

	Console.c160 = String.fromCharCode (160); //nbsp;
	Console.c160 += Console.c160+Console.c160+Console.c160;


	/* public */
	System.out.print = function (s) { 
		Console.consoleOutput (s);
	};

	/* public */
	System.out.println = function (s) {
		if (typeof s == "undefined") {
			s = "\r\n";
		} else if (s == null) {
			s = "null\r\n";
		} else {
			s = s + "\r\n";
		}
		Console.consoleOutput (s);
	};

	/* public */
	System.err.__CLASS_NAME__ = "java.io.PrintStream";

	/* public */
	System.err.print = function (s) { 
		Console.consoleOutput (s, "red");
	};

	/* public */
	System.err.println = function (s) {
		if (typeof s == "undefined") {
			s = "\r\n";
		} else if (s == null) {
			s = "null\r\n";
		} else {
			s = s + "\r\n";
		}
		Console.consoleOutput (s, "red");
	};



})(Clazz.Console, System);

Clazz.setConsoleDiv = function(d) {
	window["j2s.lib"] && (window["j2s.lib"].console = d);
};

Clazz.binaryFolders =  Clazz._Loader.binaryFolders = [ Clazz._Loader.getJ2SLibBase() ];

})(Clazz);

};
