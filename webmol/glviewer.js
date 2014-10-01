//a molecular viewer based on GLMol

//Adapted from the text sprite example from http://stemkoski.github.io/Three.js/index.html

WebMol.LabelCount = 0;

WebMol.Label = function(text, parameters) {

	this.id = WebMol.LabelCount++;
	this.stylespec = parameters || {};

	this.canvas = document.createElement('canvas');

	this.context = this.canvas.getContext('2d');

	this.sprite = new WebMol.Sprite();
	this.text = text;

};

WebMol.Label.prototype = {

	constructor : WebMol.Label,

	setContext : function() {
		// function for drawing rounded rectangles - for Label drawing
		var roundRect = function(ctx, x, y, w, h, r) {

			ctx.beginPath();
			ctx.moveTo(x + r, y);
			ctx.lineTo(x + w - r, y);
			ctx.quadraticCurveTo(x + w, y, x + w, y + r);
			ctx.lineTo(x + w, y + h - r);
			ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
			ctx.lineTo(x + r, y + h);
			ctx.quadraticCurveTo(x, y + h, x, y + h - r);
			ctx.lineTo(x, y + r);
			ctx.quadraticCurveTo(x, y, x + r, y);
			ctx.closePath();
			ctx.fill();

		};

		return function() {
			var fontMult = 2.0;
			this.showBackground = this.stylespec.showBackground;
			if(typeof(this.showBackground) == "undefined") this.showBackground = true; //default
			this.font = this.stylespec.font = this.stylespec.font ? this.stylespec.font
					: "Verdana";

			this.fontSize = this.stylespec.fontSize = this.stylespec.fontSize ? this.stylespec.fontSize
					: 20;
			this.fontSize *= fontMult;
			/** @type {colorlike} */
			this.fontColor = this.stylespec.fontColor = this.stylespec.fontColor ? this.stylespec.fontColor
					: {
						r : 255,
						g : 255,
						b : 255,
						a : 1.0
					};

			this.borderThickness = this.stylespec.borderThickness = this.stylespec.borderThickness ? this.stylespec.borderThickness
					: 4;

			this.borderColor = this.stylespec.borderColor = this.stylespec.borderColor ? this.stylespec.borderColor
					: {
						r : 0,
						g : 0,
						b : 0,
						a : 1.0
					};

			this.backgroundColor = this.stylespec.backgroundColor = this.stylespec.backgroundColor ? this.stylespec.backgroundColor
					: {
						r : 0,
						g : 0,
						b : 0,
						a : 1.0
					};

			this.position = this.stylespec.position = this.stylespec.position ? this.stylespec.position
					: {
						x : -10,
						y : 1,
						z : 1
					};
					
			//convert colors from 0-1.0 to 255
			if(this.backgroundColor instanceof WebMol.Color) this.backgroundColor = this.backgroundColor.scaled();
			if(this.borderColor instanceof WebMol.Color) this.borderColor = this.borderColor.scaled();
			if(this.fontColor instanceof WebMol.Color) this.fontColor = this.fontColor.scaled();
		

			// Should labels always be in front of model?
			this.inFront = this.stylespec.inFront = (this.stylespec.inFront !== undefined) ? this.stylespec.inFront
					: true;

			// clear canvas
			this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);

			var spriteAlignment = this.stylespec.alignment || WebMol.SpriteAlignment.topLeft;

			var bold = "";
			if(this.stylespec.bold)
				bold = "bold ";
			this.context.font = bold+this.fontSize + "px  " + this.font;

			var metrics = this.context.measureText(this.text);
			var textWidth = metrics.width;

			// background color
			this.context.fillStyle = "rgba(" + this.backgroundColor.r + ","
					+ this.backgroundColor.g + "," + this.backgroundColor.b
					+ "," + this.backgroundColor.a + ")";
			// border color
			this.context.strokeStyle = "rgba(" + this.borderColor.r + ","
					+ this.borderColor.g + "," + this.borderColor.b + ","
					+ this.borderColor.a + ")";

			this.context.lineWidth = this.borderThickness;
			if(this.showBackground) {
				roundRect(this.context, this.borderThickness / 2,
					this.borderThickness / 2, textWidth + 2*this.borderThickness,
					this.fontSize * 1.25 + 2*this.borderThickness, 6);
			// 1.25 is extra height factor for text below baseline: g,j,p,q.
			}
			
			// text color
			this.context.fillStyle = "rgba(" + this.fontColor.r + ","
					+ this.fontColor.g + "," + this.fontColor.b + ","
					+ this.fontColor.a + ")";

			this.context.fillText(this.text, this.borderThickness,
					this.fontSize + this.borderThickness, textWidth);

			// canvas contents will be used for a texture
			var texture = new WebMol.Texture(this.canvas);
			texture.needsUpdate = true;

			this.sprite.material = new WebMol.SpriteMaterial({
				map : texture,
				useScreenCoordinates : this.stylespec.useScreen,
				alignment : spriteAlignment,
				depthTest : !this.inFront
			});

			//TODO: figure out why the magic number 2.0 is needed
			this.sprite.scale.set(2.0*this.fontSize/fontMult, this.fontSize/fontMult, 1);
			this.sprite.position.set(this.position.x, this.position.y,
					this.position.z);
		};

	}(),

	// clean up material and texture
	dispose : function() {

		if (this.sprite.material.map !== undefined)
			this.sprite.material.map.dispose();
		if (this.sprite.material !== undefined)
			this.sprite.material.dispose();
	}

};

// a webmol unified interace to gmol
WebMol.GLViewer = (function() {
	// private class variables
	var numWorkers = 4; // number of threads for surface generation
	var maxVolume = 64000; // how much to break up surface calculations

	// private class helper functions

	// computes the bounding box around the provided atoms
	/**
	 * @param {Array.
	 *            <AtomSpec>} atomlist
	 * @return {Array}
	 */
	var getExtent = function(atomlist) {
		var xmin, ymin, zmin, xmax, ymax, zmax, xsum, ysum, zsum, cnt;

		xmin = ymin = zmin = 9999;
		xmax = ymax = zmax = -9999;
		xsum = ysum = zsum = cnt = 0;

		if (atomlist.length === 0)
			return [ [ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ] ];
		for (var i = 0; i < atomlist.length; i++) {
			var atom = atomlist[i];
			if (atom === undefined)
				continue;
			cnt++;
			xsum += atom.x;
			ysum += atom.y;
			zsum += atom.z;

			xmin = (xmin < atom.x) ? xmin : atom.x;
			ymin = (ymin < atom.y) ? ymin : atom.y;
			zmin = (zmin < atom.z) ? zmin : atom.z;
			xmax = (xmax > atom.x) ? xmax : atom.x;
			ymax = (ymax > atom.y) ? ymax : atom.y;
			zmax = (zmax > atom.z) ? zmax : atom.z;
		}

		return [ [ xmin, ymin, zmin ], [ xmax, ymax, zmax ],
				[ xsum / cnt, ysum / cnt, zsum / cnt ] ];
	};

	function GLViewer(element, callback, defaultcolors, nomouse) {
		// set variables
		var _viewer = this;
		var container = element;
		var id = container.id;

		var models = []; // atomistic molecular models
		var surfaces = [];
		var shapes = []; // Generic shapes
		var labels = [];
		var WIDTH = container.width();
		var HEIGHT = container.height();

		var spinner = $('<div class="glviewerSpinnerWrap" style = "position: absolute; width: 100%; height: 100%; display: table; z-index: 1;"><div class="glviewerSpinner" style="display: table-cell; text-align: center; vertical-align: middle; z-index:1"><img src="webmol/spinner.gif"></div></div>');
		$(element).append(spinner);
		spinner.hide();
		// set dimensions
		// $(container).width(WIDTH);
		// $(container).height(HEIGHT);

		var ASPECT = WIDTH / HEIGHT;
		var NEAR = 1, FAR = 800;
		var CAMERA_Z = 150;

		var renderer = new WebMol.Renderer({
			antialias : true
		});
		// renderer.sortObjects = false; // hopefully improve performance

		renderer.domElement.style.width = "100%";
		renderer.domElement.style.height = "100%";
		renderer.domElement.style.position = "absolute";
		renderer.domElement.style.top = "0px";
		renderer.domElement.style.zIndex = "0";
		container.append(renderer.domElement);
		renderer.setSize(WIDTH, HEIGHT);
		var camera = new WebMol.Camera(10, ASPECT, 1, 800);
		camera.position = new WebMol.Vector3(0, 0, CAMERA_Z);
		var vec = new WebMol.Vector3();
		camera.lookAt(vec);

		var raycaster = new WebMol.Raycaster(new WebMol.Vector3(0, 0, 0),
				new WebMol.Vector3(0, 0, 0));
		var projector = new WebMol.Projector();
		var mouseVector = new WebMol.Vector3(0, 0, 0);

		var scene = null;
		var rotationGroup = null; // which contains modelGroup
		var modelGroup = null;

		var bgColor = 0x000000;
		var fov = 20;
		var fogStart = 0.4;
		var slabNear = -50; // relative to the center of rotationGroup
		var slabFar = 50;

		// UI variables
		var cq = new WebMol.Quaternion(0, 0, 0, 1);
		var dq = new WebMol.Quaternion(0, 0, 0, 1);
		var isDragging = false;
		var mouseStartX = 0;
		var mouseStartY = 0;
		var touchDistanceStart = 0;
		var currentModelPos = 0;
		var cz = 0;
		var cslabNear = 0;
		var cslabFar = 0;

		var setSlabAndFog = function() {
			var center = camera.position.z - rotationGroup.position.z;
			if (center < 1)
				center = 1;
			camera.near = center + slabNear;
			if (camera.near < 1)
				camera.near = 1;
			camera.far = center + slabFar;
			if (camera.near + 1 > camera.far)
				camera.far = camera.near + 1;
			if (camera instanceof WebMol.Camera) {
				camera.fov = fov;
			} else {
				camera.right = center * Math.tan(Math.PI / 180 * fov);
				camera.left = -camera.right;
				camera.top = camera.right / ASPECT;
				camera.bottom = -camera.top;
			}
			camera.updateProjectionMatrix();
			scene.fog.near = camera.near + fogStart
					* (camera.far - camera.near);
			// if (scene.fog.near > center) scene.fog.near = center;
			scene.fog.far = camera.far;
		};

		// display scene
		var show = function() {
			if (!scene)
				return;

			// var time = new Date();
			setSlabAndFog();
			renderer.render(scene, camera);
			// console.log("rendered in " + (+new Date() - time) + "ms");
		};

		var initializeScene = function() {

			scene = new WebMol.Scene();
			scene.fog = new WebMol.Fog(bgColor, 100, 200);

			modelGroup = new WebMol.Object3D();
			rotationGroup = new WebMol.Object3D();
			rotationGroup.useQuaternion = true;
			rotationGroup.quaternion = new WebMol.Quaternion(0, 0, 0, 1);
			rotationGroup.add(modelGroup);

			scene.add(rotationGroup);

			// setup lights
			var directionalLight = new WebMol.Light(0xFFFFFF);
			directionalLight.position = new WebMol.Vector3(0.2, 0.2, 1)
					.normalize();
			directionalLight.intensity = 1.0;
			scene.add(directionalLight);
		};

		initializeScene();

		renderer.setClearColorHex(bgColor, 1.0);
		scene.fog.color = WebMol.CC.color(bgColor);

		var clickedAtom = null;
		// enable mouse support
		var glDOM = $(renderer.domElement);

		// Checks for selection intersects on mousedown
		var handleClickSelection = function(mouseX, mouseY) {

			var mouse = {
				x : mouseX,
				y : mouseY,
				z : -1.0
			};
			mouseVector.set(mouse.x, mouse.y, mouse.z);
			projector.unprojectVector(mouseVector, camera);
			mouseVector.sub(camera.position).normalize();

			raycaster.set(camera.position, mouseVector);

			var clickables = [], intersects = [];
			var i, il;

			for (i = 0, il = models.length; i < il; i++) {
				var model = models[i];

				var atoms = model.selectedAtoms({
					clickable : true
				});
				clickables = clickables.concat(atoms);

			}

			for (i = 0, il = shapes.length; i < il; i++) {

				var shape = shapes[i];
				if (shape.clickable) {
					clickables.push(shape);
				}

			}

			intersects = raycaster.intersectObjects(modelGroup, clickables);

			if (intersects.length) {
				var selected = intersects[0].clickable;
				if (selected.callback !== undefined
						&& typeof (selected.callback) === "function") {
					selected.callback(selected, _viewer);
				}
			}

			show();
		};

		var calcTouchDistance = function(ev) { // distance between first two
												// fingers
			var xdiff = ev.originalEvent.targetTouches[0].pageX
					- ev.originalEvent.targetTouches[1].pageX;
			var ydiff = ev.originalEvent.targetTouches[0].pageY
					- ev.originalEvent.targetTouches[1].pageY;
			return Math.sqrt(xdiff * xdiff + ydiff * ydiff);
		}

		if (!nomouse) {
			// user can request that the mouse handlers not be installed
			glDOM.bind('mousedown touchstart', function(ev) {
				ev.preventDefault();
				if (!scene)
					return;
				var x = ev.pageX, y = ev.pageY;
				if (ev.originalEvent.targetTouches
						&& ev.originalEvent.targetTouches[0]) {
					x = ev.originalEvent.targetTouches[0].pageX;
					y = ev.originalEvent.targetTouches[0].pageY;
				}
				if (x === undefined)
					return;
				isDragging = true;
				clickedAtom = null;
				mouseButton = ev.which;
				mouseStartX = x;
				mouseStartY = y;
				touchDistanceStart = 0;
				if (ev.originalEvent.targetTouches
						&& ev.originalEvent.targetTouches.length == 2) {
					touchDistanceStart = calcTouchDistance(ev);
				}
				cq = rotationGroup.quaternion;
				cz = rotationGroup.position.z;
				currentModelPos = modelGroup.position.clone();
				cslabNear = slabNear;
				cslabFar = slabFar;

				// handle selection
				var mouseX = (x / $(window).width()) * 2 - 1;
				var mouseY = -(y / HEIGHT) * 2 + 1;
				handleClickSelection(mouseX, mouseY, ev, container);

			});

			glDOM.bind('DOMMouseScroll mousewheel', function(ev) { // Zoom
				ev.preventDefault();
				if (!scene)
					return;
				var scaleFactor = (CAMERA_Z - rotationGroup.position.z) * 0.85;
				if (ev.originalEvent.detail) { // Webkit
					rotationGroup.position.z += scaleFactor
							* ev.originalEvent.detail / 10;
				} else if (ev.originalEvent.wheelDelta) { // Firefox
					rotationGroup.position.z -= scaleFactor
							* ev.originalEvent.wheelDelta / 400;
				}

				show();
			});

			glDOM.bind("contextmenu", function(ev) {
				ev.preventDefault();
			});
			$('body').bind('mouseup touchend', function(ev) {
				isDragging = false;
			});

			glDOM
					.bind(
							'mousemove touchmove',
							function(ev) { // touchmove
								ev.preventDefault();
								if (!scene)
									return;
								if (!isDragging)
									return;
								var mode = 0;
								var modeRadio = $('input[name=' + id
										+ '_mouseMode]:checked');
								if (modeRadio.length > 0)
									mode = parseInt(modeRadio.val());

								var x = ev.pageX, y = ev.pageY;
								if (ev.originalEvent.targetTouches
										&& ev.originalEvent.targetTouches[0]) {
									x = ev.originalEvent.targetTouches[0].pageX;
									y = ev.originalEvent.targetTouches[0].pageY;
								}
								if (x === undefined)
									return;
								var dx = (x - mouseStartX) / WIDTH;
								var dy = (y - mouseStartY) / HEIGHT;
								// check for pinch
								if (touchDistanceStart != 0
										&& ev.originalEvent.targetTouches
										&& ev.originalEvent.targetTouches.length == 2) {
									var newdist = calcTouchDistance(ev);
									// change to zoom
									mode = 2;
									dy = (touchDistanceStart - newdist) * 2
											/ (WIDTH + HEIGHT);
								} else if (ev.originalEvent.targetTouches
										&& ev.originalEvent.targetTouches.length == 3) {
									// translate
									mode = 1;
								}

								var r = Math.sqrt(dx * dx + dy * dy);
								var scaleFactor;
								if (mode == 3
										|| (mouseButton == 3 && ev.ctrlKey)) { // Slab
									slabNear = cslabNear + dx * 100;
									slabFar = cslabFar + dy * 100;
								} else if (mode == 2 || mouseButton == 3
										|| ev.shiftKey) { // Zoom
									scaleFactor = (CAMERA_Z - rotationGroup.position.z) * 0.85;
									if (scaleFactor < 80)
										scaleFactor = 80;
									rotationGroup.position.z = cz - dy
											* scaleFactor;
								} else if (mode == 1 || mouseButton == 2
										|| ev.ctrlKey) { // Translate
									scaleFactor = (CAMERA_Z - rotationGroup.position.z) * 0.85;
									if (scaleFactor < 20)
										scaleFactor = 20;
									var translationByScreen = new WebMol.Vector3(
											dx * scaleFactor,
											-dy * scaleFactor, 0);
									var q = rotationGroup.quaternion;
									var qinv = new WebMol.Quaternion(q.x, q.y,
											q.z, q.w).inverse().normalize();
									var translation = translationByScreen
											.applyQuaternion(qinv);
									modelGroup.position.x = currentModelPos.x
											+ translation.x;
									modelGroup.position.y = currentModelPos.y
											+ translation.y;
									modelGroup.position.z = currentModelPos.z
											+ translation.z;
								} else if ((mode === 0 || mouseButton == 1)
										&& r !== 0) { // Rotate
									var rs = Math.sin(r * Math.PI) / r;
									dq.x = Math.cos(r * Math.PI);
									dq.y = 0;
									dq.z = rs * dx;
									dq.w = -rs * dy;
									rotationGroup.quaternion = new WebMol.Quaternion(
											1, 0, 0, 0);
									rotationGroup.quaternion.multiply(dq);
									rotationGroup.quaternion.multiply(cq);
								}
								show();
							});
		}
		// public methods
		/**
		 * Set the background color (default white)
		 * 
		 * @function WebMol.GLViewer#setBackgroundColor
		 * @param {number}
		 *            hex Hexcode specified background color
		 * @param {number}
		 *            a Alpha level (default 1.0)
		 * 
		 * @example
		 * 
		 * //Set 'myviewer' background color to white
		 * myviewer.setBackgroundColor(0xffffff)
		 * 
		 */
		this.setBackgroundColor = function(hex, a) {
			a = a | 1.0;
			bgColor = hex;
			renderer.setClearColorHex(hex, a);
			scene.fog.color = WebMol.CC.color(hex);
			show();
		};

		/**
		 * Set viewer width
		 * 
		 * @function WebMol.GLViewer#setWidth
		 * @param {number}
		 *            w Width in pixels
		 */
		this.setWidth = function(w) {
			WIDTH = w || WIDTH;
			renderer.setSize(WIDTH, HEIGHT);
		};

		/**
		 * Set viewer height
		 * 
		 * @function WebMol.GLViewer#setHeight
		 * @param {number}
		 *            h Height in pixels
		 */
		this.setHeight = function(h) {
			HEIGHT = h || HEIGHT;
			renderer.setSize(WIDTH, HEIGHT);
		};

		/**
		 * Resize viewer according to containing HTML element's dimensions
		 * 
		 * @function WebMol.GLViewer#resize
		 */
		this.resize = function() {
			WIDTH = container.width();
			HEIGHT = container.height();
			ASPECT = WIDTH / HEIGHT;
			renderer.setSize(WIDTH, HEIGHT);
			camera.aspect = ASPECT;
			camera.updateProjectionMatrix();
			show();
		};

		$(window).resize(this.resize);

		/**
		 * Return specified model
		 * 
		 * @function WebMol.GLViewer#getModel
		 * @param {number}
		 *            [id=last model id] - Retrieve model with specified id
		 * @default Returns last model added to viewer
		 * @return {GLModel}
		 * 
		 * @example // Retrieve reference to first GLModel added var m =
		 *          glviewer.getModel(0);
		 */
		this.getModel = function(id) {
			id = id || models.length - 1;
			return models[id];
		};

		/**
		 * Rotate scene by angle degrees around axis
		 * 
		 * @function WebMol.GLViewer#rotate
		 * @param {number}
		 *            [angle] - Angle, in degrees, to rotate by.
		 * @param {string}
		 *            [angle] - Axis ("x", "y", or "z") to rotate around.
		 *            Default "y"
		 * 
		 */
		this.rotate = function(angle, axis) {
			if (typeof (axis) === "undefined") {
				axis = "y";
			}
			var i = 0, j = 0, k = 0;
			var rangle = Math.PI * angle / 180.0;
			var s = Math.sin(rangle / 2.0);
			var c = Math.cos(rangle / 2.0);
			if (axis == "x")
				i = s;
			if (axis == "y")
				j = s;
			if (axis == "z")
				k = s;

			var q = new WebMol.Quaternion(i, j, k, c).normalize();
			rotationGroup.quaternion.multiply(q);
			show();
		};

		this.getView = function() {
			if (!modelGroup)
				return [ 0, 0, 0, 0, 0, 0, 0, 1 ];
			var pos = modelGroup.position;
			var q = rotationGroup.quaternion;
			return [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y,
					q.z, q.w ];
		};

		this.setView = function(arg) {

			if (arg === undefined
					|| !(arg instanceof Array || arg.length !== 8))
				return;

			if (!modelGroup || !rotationGroup)
				return;
			modelGroup.position.x = arg[0];
			modelGroup.position.y = arg[1];
			modelGroup.position.z = arg[2];
			rotationGroup.position.z = arg[3];
			rotationGroup.quaternion.x = arg[4];
			rotationGroup.quaternion.y = arg[5];
			rotationGroup.quaternion.z = arg[6];
			rotationGroup.quaternion.w = arg[7];
			if(typeof(arg[8]) != "undefined") {
				rotationGroup.position.x = arg[8];
				rotationGroup.position.y = arg[9];
			}
			show();
		};

		// apply styles, models, etc in viewer
		/**
		 * Render current state of viewer, after adding/removing models,
		 * applying styles, etc.
		 * 
		 * @function WebMol.GLViewer#render
		 */
		this.render = function() {

			// spinner.show();
			var time1 = new Date();
			var view = this.getView();
			var i;
			for (i = 0; i < models.length; i++) {
				if (models[i]) {
					models[i].globj(modelGroup);
				}
			}

			for (i = 0; i < shapes.length; i++) {
				if (shapes[i]) {
					shapes[i].globj(modelGroup);
				}
			}

			for (i in surfaces) { // this is an array with possible holes
				if (surfaces.hasOwnProperty(i)) {
					var geo = surfaces[i].geo;
					// async surface generation can cause
					// the geometry to be webgl initialized before it is fully
					// formed; force various recalculations until full surface
					// is
					// available
					if (!surfaces[i].finished) {
						geo.verticesNeedUpdate = true;
						geo.elementsNeedUpdate = true;
						geo.normalsNeedUpdate = true;
						geo.colorsNeedUpdate = true;
						geo.buffersNeedUpdate = true;
						geo.boundingSphere = null;

						if (surfaces[i].done)
							surfaces[i].finished = true;

						// remove partially rendered surface
						if (surfaces[i].lastGL)
							modelGroup.remove(surfaces[i].lastGL);

						// create new surface
						var smesh = null;
						if(surfaces[i].mat instanceof WebMol.LineBasicMaterial) {
							//special case line meshes
							smesh = new WebMol.Line(geo, surfaces[i].mat);
						}
						else {
							smesh = new WebMol.Mesh(geo, surfaces[i].mat);
						}
						surfaces[i].lastGL = smesh;
						modelGroup.add(smesh);
					} // else final surface already there
				}
			}
			this.setView(view); // Calls show() => renderer render
			var time2 = new Date();
			spinner.hide();
			console.log("render time: " + (time2 - time1));
		};

		/**
		 * 
		 * @param {AtomSpec}
		 *            sel
		 * @return {Array.<AtomSpec>}
		 */
		function getAtomsFromSel(sel) {
			var atoms = [];
			if (typeof (sel) === "undefined")
				sel = {};

			var ms = [];
			var i;

			if (typeof sel.model === "undefined") {
				for (i = 0; i < models.length; i++) {
					if (models[i])
						ms.push(models[i]);
				}
			} else { // specific to some models
				ms = sel.model;
				if (!$.isArray(ms))
					ms = [ ms ];
			}

			for (i = 0; i < ms.length; i++) {
				atoms = atoms.concat(ms[i].selectedAtoms(sel));
			}

			return atoms;
		}

		/**
		 * 
		 * @param {AtomSpec}
		 *            atom
		 * @param {AtomSpec}
		 *            sel
		 * @return {boolean}
		 */
		function atomIsSelected(atom, sel) {
			if (typeof (sel) === "undefined")
				sel = {};

			var ms = [];
			var i;

			if (typeof sel.model === "undefined") {
				for (i = 0; i < models.length; i++) {
					if (models[i])
						ms.push(models[i]);
				}
			} else { // specific to some models
				ms = sel.model;
				if (!$.isArray(ms))
					ms = [ ms ];
			}

			for (i = 0; i < ms.length; i++) {
				if (ms[i].atomIsSelected(atom, sel))
					return true;
			}

			return false;
		}

		/**
		 * Return pdb output of selected atoms (if atoms from pdb input)
		 * 
		 * @function WebMol.GLViewer#pdbData
		 * @param {Object}
		 *            [sel] - Selection specification specifying model and atom
		 *            properties to select. Default: all atoms in viewer
		 * @return {string} PDB string of selected atoms
		 */
		this.pdbData = function(sel) {
			var atoms = getAtomsFromSel(sel);
			var ret = "";
			for (var i = 0, n = atoms.length; i < n; ++i) {
				ret += atoms[i].pdbline + "\n";
			}
			return ret;
		};

		/**
		 * Zoom current view by a constant factor
		 * 
		 * @function WebMol.GLViewer#zoom
		 * @param {number}
		 *            [factor] - Magnification factor. Values greater than 1
		 *            will zoom in, less than one will zoom out. Default 2.
		 * 
		 */
		this.zoom = function(factor) {
			var factor = factor || 2;
			var scale = (CAMERA_Z - rotationGroup.position.z) / factor;
			rotationGroup.position.z = CAMERA_Z - scale;
			show();
		};
		

		/**
		 * Zoom to center of atom selection
		 * 
		 * @function WebMol.GLViewer#zoomTo
		 * @param {Object}
		 *            [sel] - Selection specification specifying model and atom
		 *            properties to select. Default: all atoms in viewer
		 * 
		 * @example // Assuming we have created a model of a protein with
		 *          multiple chains (e.g. from a PDB file), focus on atoms in
		 *          chain B glviewer.zoomTo({chain: 'B'});
		 *  // Focus on centroid of all atoms of all models in this
		 * viewer glviewer.zoomTo(); // (equivalent to glviewer.zoomTo({}) )
		 */
		this.zoomTo = function(sel) {
			var atoms = getAtomsFromSel(sel).concat(shapes);
			var allatoms = getAtomsFromSel({}).concat(shapes);
			var tmp = getExtent(atoms);
			var alltmp = getExtent(allatoms);
			// use selection for center
			var center = new WebMol.Vector3(tmp[2][0], tmp[2][1], tmp[2][2]);
			modelGroup.position = center.multiplyScalar(-1);
			// but all for bounding box
			var x = alltmp[1][0] - alltmp[0][0], y = alltmp[1][1]
					- alltmp[0][1], z = alltmp[1][2] - alltmp[0][2];

			var maxD = Math.sqrt(x * x + y * y + z * z);
			if (maxD < 25)
				maxD = 25;

			// use full bounding box for slab/fog
			slabNear = -maxD / 1.9;
			slabFar = maxD / 3;

			// for zoom, use selection box
			x = tmp[1][0] - tmp[0][0];
			y = tmp[1][1] - tmp[0][1];
			z = tmp[1][2] - tmp[0][2];
			maxD = Math.sqrt(x * x + y * y + z * z);
			if (maxD < 25)
				maxD = 25;

			rotationGroup.position.z = -(maxD * 0.35
					/ Math.tan(Math.PI / 180.0 * camera.fov / 2) - 150);

			show();
		};

		/**
		 * Add label to viewer
		 * 
		 * @function WebMol.GLViewer#addLabel
		 * @param {string}
		 *            text - Label text
		 * @param {Object}
		 *            data - Label style specification
		 * @return {WebMol.Label}
		 * 
		 * @example
		 *  // Assuming glviewer contains a model representing a protein, label
		 * all alpha carbons with their residue name
		 *  // Select all alpha carbons (have property atom : "CA") from last
		 * model added var atoms =
		 * glviewer.getModel().selectedAtoms({atom:"CA"}); var labels = [];
		 * 
		 * for (var a in atoms) { var atom = atoms[a];
		 *  // Create label at alpha carbon's position displaying atom's residue
		 * and residue number var labelText = atom.resname + " " + atom.resi;
		 * 
		 * var l = glviewer.createLabel(labelText, {fontSize: 12, position: {x:
		 * atom.x, y: atom.y, z: atom.z});
		 * 
		 * labels.push(l); }
		 *  // Render labels glviewer.render();
		 */
		this.addLabel = function(text, data) {
			var label = new WebMol.Label(text, data);
			label.setContext();
			modelGroup.add(label.sprite);
			labels.push(label);
			show();
			return label;
		};

		/**
		 * Remove label from viewer
		 * 
		 * @function WebMol.GLViewer#removeLabel
		 * @param {WebMol.Label}
		 *            label - WebMol label
		 * 
		 * @example // Remove labels created in [addLabel example]{@link WebMol.GLViewer#addLabel}
		 * 
		 * for (var i = 0; i < labels.length; i++) {
		 * glviewer.removeLabel(label); }
		 * 
		 * glviewer.render();
		 */
		this.removeLabel = function(label) {
			labels.remove(label);
			label.dispose();
			modelGroup.remove(label.sprite);
		};

		this.removeAllLabels = function() {
			for (var i = 0; i < labels.length; i++) {
				modelGroup.remove(labels[i].sprite);
			}
			labels = [];
		};
		
		// Modify label style
		/**
		 * Modify existing label's style
		 * 
		 * @function WebMol.GLViewer#setLabelStyle
		 * @param {WebMol.Label}
		 *            label - WebMol label
		 * @param {Object}
		 *            stylespec - Label style specification
		 * @return {WebMol.Label}
		 */
		this.setLabelStyle = function(label, stylespec) {

			label.dispose();
			label.stylespec = stylespec;
			label.setContext();
			modelGroup.add(label.sprite);

			return label;

		};

		// Change label text
		/**
		 * Modify existing label's text
		 * 
		 * @function WebMol.GLViewer#setLabelText
		 * @param {WebMol.Label}
		 *            label - WebMol label
		 * @param {String}
		 *            text - Label text
		 * @return {WebMol.Label}
		 */
		this.setLabelText = function(label, text) {

			label.dispose();
			label.text = text;
			label.setContext();
			modelGroup.add(label.sprite);

			return label;

		};

		this.addShape = function(shapeSpec) {
			shapeSpec = shapeSpec || {};
			var shape = new WebMol.GLShape(shapes.length, shapeSpec);
			shapes.push(shape);

			return shape;

		};

		this.removeShape = function(shape) {
			if (!shape)
				return;
			shape.removegl(modelGroup);
			delete shapes[shape.id];
			// clear off back of model array
			while (shapes.length > 0
					&& typeof (shapes[shapes.length - 1]) === "undefined")
				shapes.pop();
		};
		
		this.removeAllShapes = function() {
			for (var i = 0; i < shapes.length; i++) {
				var shape = shapes[i];
				shape.removegl(modelGroup);
			}
			shapes = [];
		}

		this.addSphere = function(spec) {
			var s = new WebMol.GLShape(shapes.length);
			spec = spec || {};
			s.addSphere(spec);
			shapes.push(s);

			return s;
		};

		this.addArrow = function(spec) {
			var s = new WebMol.GLShape(shapes.length);
			spec = spec || {};
			s.addArrow(spec);
			shapes.push(s);

			return s;
		};
		
		this.addCylinder = function(spec) {
			var s = new WebMol.GLShape(shapes.length);
			spec = spec || {};
			s.addCylinder(spec);
			shapes.push(s);

			return s;
		};

		this.addCustom = function(spec) {
			var s = new WebMol.GLShape(shapes.length);
			spec = spec || {};
			s.addCustom(spec);
			shapes.push(s);

			return s;
		};

		this.addVolumetricData = function(data, format, spec) {
			var s = new WebMol.GLShape(shapes.length);
			spec = spec || {};
			s.addVolumetricData(data, format, spec);
			shapes.push(s);

			return s;
		};

		this.addModel = function(data, format) {

			var m = new WebMol.GLModel(models.length, defaultcolors);
			m.addMolData(data, format);
			models.push(m);

			return m;
		};

		this.removeModel = function(model) {
			if (!model)
				return;
			model.removegl(modelGroup);
			delete models[model.getID()];
			// clear off back of model array
			while (models.length > 0
					&& typeof (models[models.length - 1]) === "undefined")
				models.pop();
		};

		this.removeAllModels = function() {
			for (var i = 0; i < models.length; i++) {
				var model = models[i];
				model.removegl(modelGroup);

			}
			models = [];
		};

		this.createModelFrom = function(sel, extract) {
			var m = new WebMol.GLModel(models.length, defaultcolors);
			for (var i = 0; i < models.length; i++) {
				if (models[i]) {
					var atoms = models[i].selectedAtoms(sel);
					m.addAtoms(atoms);
					if (extract)
						models[i].removeAtoms(atoms);
				}
			}
			models.push(m);
			return m;
		};

		function applyToModels(func, sel, value1, value2) {
			for (var i = 0; i < models.length; i++) {
				if (models[i]) {
					models[i][func](sel, value1, value2);
				}
			}
		}

		this.setStyle = function(sel, style) {
			applyToModels("setStyle", sel, style, false);
		};

		this.addStyle = function(sel, style) {
			applyToModels("setStyle", sel, style, true);
		};

		this.setColorByProperty = function(sel, prop, scheme) {
			applyToModels("setColorByProperty", sel, prop, scheme);
		};

		this.setColorByElement = function(sel, colors) {
			applyToModels("setColorByElement", sel, colors);
		};

		/**
		 * 
		 * @param {Array.
		 *            <AtomSpec>} atomlist
		 * @param {Array}
		 *            extent
		 * @return {Array}
		 */
		var getAtomsWithin = function(atomlist, extent) {
			var ret = [];

			for (var i = 0; i < atomlist.length; i++) {
				var atom = atomlist[i];
				if (typeof (atom) == "undefined")
					continue;

				if (atom.x < extent[0][0] || atom.x > extent[1][0])
					continue;
				if (atom.y < extent[0][1] || atom.y > extent[1][1])
					continue;
				if (atom.z < extent[0][2] || atom.z > extent[1][2])
					continue;
				ret.push(i);
			}
			return ret;
		};

		// return volume of extent
		var volume = function(extent) {
			var w = extent[1][0] - extent[0][0];
			var h = extent[1][1] - extent[0][1];
			var d = extent[1][2] - extent[0][2];
			return w * h * d;
		}; // volume
		/*
		 * Break up bounding box/atoms into smaller pieces so we can parallelize
		 * with webworkers and also limit the size of the working memory Returns
		 * a list of bounding boxes with the corresponding atoms. These extents
		 * are expanded by 4 angstroms on each side.
		 */
		/**
		 * 
		 * @param {Array}
		 *            extent
		 * @param {Array.
		 *            <AtomSpec>} atomlist
		 * @param {Array.
		 *            <AtomSpec>} atomstoshow
		 * @return {Array}
		 */
		var carveUpExtent = function(extent, atomlist, atomstoshow) {
			var ret = [];

			var copyExtent = function(extent) {
				// copy just the dimensions
				var ret = [];
				ret[0] = [ extent[0][0], extent[0][1], extent[0][2] ];
				ret[1] = [ extent[1][0], extent[1][1], extent[1][2] ];
				return ret;
			}; // copyExtent
			var splitExtentR = function(extent) {
				// recursively split until volume is below maxVol
				if (volume(extent) < maxVolume) {
					return [ extent ];
				} else {
					// find longest edge
					var w = extent[1][0] - extent[0][0];
					var h = extent[1][1] - extent[0][1];
					var d = extent[1][2] - extent[0][2];

					var index;

					if (w > h && w > d) {
						index = 0;
					} else if (h > w && h > d) {
						index = 1;
					} else {
						index = 2;
					}

					// create two halves, splitting at index
					var a = copyExtent(extent);
					var b = copyExtent(extent);
					var mid = (extent[1][index] - extent[0][index]) / 2
							+ extent[0][index];
					a[1][index] = mid;
					b[0][index] = mid;

					var alist = splitExtentR(a);
					var blist = splitExtentR(b);
					return alist.concat(blist);
				}
			}; // splitExtentR

			// divide up extent
			var splits = splitExtentR(extent);
			// now compute atoms within expanded (this could be more efficient)
			var off = 6; // enough for water and 2*r, also depends on scale
			// factor
			for (var i = 0, n = splits.length; i < n; i++) {
				var e = copyExtent(splits[i]);
				e[0][0] -= off;
				e[0][1] -= off;
				e[0][2] -= off;
				e[1][0] += off;
				e[1][1] += off;
				e[1][2] += off;

				var atoms = getAtomsWithin(atomlist, e);
				var toshow = getAtomsWithin(atomstoshow, splits[i]);

				// ultimately, divide up by atom for best meshing
				ret.push({
					extent : splits[i],
					atoms : atoms,
					toshow : toshow
				});
			}

			return ret;
		};

		// create a mesh defined from the passed vertices and faces and material
		// Just create a single geometry chunk - broken up whether sync or not
		/**
		 * 
		 * @param {Array.
		 *            <AtomSpec>} atoms
		 * @param {{vertices:number,faces:number}}
		 *            VandF
		 * @param {WebMol.MeshLambertMaterial}
		 *            mat
		 * @return {WebMol.Mesh}
		 */
		var generateSurfaceMesh = function(atoms, VandF, mat) {

			var geo = new WebMol.Geometry(true);
			// Only one group per call to generate surface mesh (addSurface
			// should split up mesh render)
			var geoGroup = geo.updateGeoGroup(0);

			var vertexArray = geoGroup.vertexArray;
			// reconstruct vertices and faces
			var v = VandF['vertices'];
			var offset;
			var i, il;
			for (i = 0, il = v.length; i < il; i++) {
				offset = geoGroup.vertices * 3;
				vertexArray[offset] = v[i].x;
				vertexArray[offset + 1] = v[i].y;
				vertexArray[offset + 2] = v[i].z;
				geoGroup.vertices++;
			}

			var faces = VandF['faces'];
			geoGroup.faceidx = faces.length;// *3;
			geo.initTypedArrays();

			// set colors for vertices
			var colors = [];
			for (i = 0, il = atoms.length; i < il; i++) {
				var atom = atoms[i];
				if (atom) {
					if (typeof (atom.surfaceColor) != "undefined") {
						colors[i] = WebMol.CC.color(atom.surfaceColor);
					} else if (atom.color) // map from atom
						colors[i] = WebMol.CC.color(atom.color);
				}
			}

			var verts = geoGroup.vertexArray;
			var colorArray = geoGroup.colorArray;
			var normalArray = geoGroup.normalArray;
			var vA, vB, vC, norm;

			// Setup colors, faces, and normals
			for (i = 0, il = faces.length; i < il; i += 3) {

				// var a = faces[i].a, b = faces[i].b, c = faces[i].c;
				var a = faces[i], b = faces[i + 1], c = faces[i + 2];
				var A = v[a]['atomid'];
				var B = v[b]['atomid'];
				var C = v[c]['atomid'];

				var offsetA = a * 3, offsetB = b * 3, offsetC = c * 3;

				colorArray[offsetA] = colors[A].r;
				colorArray[offsetA + 1] = colors[A].g;
				colorArray[offsetA + 2] = colors[A].b;
				colorArray[offsetB] = colors[B].r;
				colorArray[offsetB + 1] = colors[B].g;
				colorArray[offsetB + 2] = colors[B].b;
				colorArray[offsetC] = colors[C].r;
				colorArray[offsetC + 1] = colors[C].g;
				colorArray[offsetC + 2] = colors[C].b;

				// setup Normals

				vA = new WebMol.Vector3(verts[offsetA], verts[offsetA + 1],
						verts[offsetA + 2]);
				vB = new WebMol.Vector3(verts[offsetB], verts[offsetB + 1],
						verts[offsetB + 2]);
				vC = new WebMol.Vector3(verts[offsetC], verts[offsetC + 1],
						verts[offsetC + 2]);

				vC.subVectors(vC, vB);
				vA.subVectors(vA, vB);
				vC.cross(vA);

				// face normal
				norm = vC;
				norm.normalize();

				normalArray[offsetA] += norm.x;
				normalArray[offsetB] += norm.x;
				normalArray[offsetC] += norm.x;
				normalArray[offsetA + 1] += norm.y;
				normalArray[offsetB + 1] += norm.y;
				normalArray[offsetC + 1] += norm.y;
				normalArray[offsetA + 2] += norm.z;
				normalArray[offsetB + 2] += norm.z;
				normalArray[offsetC + 2] += norm.z;

			}
			geoGroup.faceArray = new Uint16Array(faces);
			var mesh = new WebMol.Mesh(geo, mat);
			mesh.doubleSided = true;

			return mesh;
		};

		// do same thing as worker in main thread
		/**
		 * 
		 * @param {WebMol.SurfaceType}
		 *            type
		 * @param {Array}
		 *            expandedExtent
		 * @param {Array}
		 *            extendedAtoms
		 * @param {Array}
		 *            atomsToShow
		 * @param {Array.
		 *            <AtomSpec>} atoms
		 * @param {number}
		 *            vol
		 * @return {Object}
		 */
		var generateMeshSyncHelper = function(type, expandedExtent,
				extendedAtoms, atomsToShow, atoms, vol) {
			var time = new Date();
			var ps = new WebMol.ProteinSurface();
			ps.initparm(expandedExtent, (type === 1) ? false : true, vol);

			var time2 = new Date();
			console.log("initialize " + (time2 - time) + "ms");

			ps.fillvoxels(atoms, extendedAtoms);

			var time3 = new Date();
			console.log("fillvoxels " + (time3 - time2) + "  " + (time3 - time)
					+ "ms");

			ps.buildboundary();

			if (type == WebMol.SurfaceType.SES) {
				ps.fastdistancemap();
				ps.boundingatom(false);
				ps.fillvoxelswaals(atoms, extendedAtoms);
			}

			var time4 = new Date();
			console.log("buildboundaryetc " + (time4 - time3) + "  "
					+ (time4 - time) + "ms");

			ps.marchingcube(type);

			var time5 = new Date();
			console.log("marching cube " + (time5 - time4) + "  "
					+ (time5 - time) + "ms");

			return ps.getFacesAndVertices(atomsToShow);
		};

		/**
		 * 
		 * @param {matSpec}
		 *            style
		 * @return {WebMol.MeshLambertMaterial}
		 */
		function getMatWithStyle(style) {
			var mat = new WebMol.MeshLambertMaterial();
			mat.vertexColors = WebMol.VertexColors;

			for ( var prop in style) {
				if (prop === "color") {
					mat[prop] = WebMol.CC.color(style.color);
					delete mat.vertexColors; // ignore
				} else if (prop === "map") {
					// ignore
				} else if (style.hasOwnProperty(prop))
					mat[prop] = style[prop];
			}
			if (style.opacity !== undefined) {
				if (style.opacity === 1)
					mat.transparent = false;
				else
					mat.transparent = true;
			}

			return mat;
		}

		// get the min and max values of the specified property in the provided
		// atoms
		function getPropertyRange(atomlist, prop) {
			var min = Number.POSITIVE_INFINITY;
			var max = Number.NEGATIVE_INFINITY;

			for (var i = 0, n = atomlist.length; i < n; i++) {
				var atom = atomlist[i];
				if (atom.properties
						&& typeof (atom.properties[prop]) != "undefined") {
					var val = atom.properties[prop];
					if (val < min)
						min = val;
					if (val > max)
						max = val;
				}
			}

			if (!isFinite(min) && !isFinite(max))
				min = max = 0;
			else if (!isFinite(min))
				min = max;
			else if (!isFinite(max))
				max = min;

			return [ min, max ];
		}

		
		// Adds an explicit mesh as a surface object.
		this.addMesh = function(mesh) {
			var surfobj = {
				geo : mesh.geometry,
				mat : mesh.material,
				done : true,
				finished : false //the rendered finishes surfaces when they are done
			};
			var surfid = surfaces.length;
			surfaces[surfid] = surfobj;
			return surfid;
		}

		// add a surface
		this.addSurface = function(type, style, atomsel, allsel, focus) {
			// type 1: VDW 3: SAS 4: MS 2: SES
			// if sync is true, does all work in main thread, otherwise uses
			// workers
			// with workers, must ensure group is the actual modelgroup since
			// surface
			// will get added asynchronously
			// all atoms in atomlist are used to compute surfaces, but only the
			// surfaces
			// of atomsToShow are displayed (e.g., for showing cavities)
			// if focusSele is specified, will start rending surface around the
			// atoms specified by this selection
			var atomsToShow = getAtomsFromSel(atomsel);
			var atomlist = getAtomsFromSel(allsel);
			var focusSele = getAtomsFromSel(focus);
			var atom;

			var time = new Date();

			var mat = getMatWithStyle(style);

			var extent = getExtent(atomsToShow);

			var i, il;
			if (style['map'] && style['map']['prop']) {
				// map color space using already set atom properties
				/** @type {AtomSpec} */
				var prop = style['map']['prop'];
				/** @type {ColorScheme} */
				var scheme = style['map']['scheme'] || new WebMol.RWB();
				var range = scheme.range();
				if (!range) {
					range = getPropertyRange(atomsToShow, prop);
				}

				for (i = 0, il = atomsToShow.length; i < il; i++) {
					atom = atomsToShow[i];
					atom.surfaceColor = scheme.valueToHex(
							atom.properties[prop], range);
				}
			}

			var totalVol = volume(extent); // used to scale resolution
			var extents = carveUpExtent(extent, atomlist, atomsToShow);

			if (focusSele && focusSele.length && focusSele.length > 0) {
				var seleExtent = getExtent(focusSele);
				// sort by how close to center of seleExtent
				var sortFunc = function(a, b) {
					var distSq = function(ex, sele) {
						// distance from e (which has no center of mass) and
						// sele which does
						var e = ex.extent;
						var x = e[1][0] - e[0][0];
						var y = e[1][1] - e[0][1];
						var z = e[1][2] - e[0][2];
						var dx = (x - sele[2][0]);
						dx *= dx;
						var dy = (y - sele[2][1]);
						dy *= dy;
						var dz = (z - sele[2][2]);
						dz *= dz;

						return dx + dy + dz;
					};
					var d1 = distSq(a, seleExtent);
					var d2 = distSq(b, seleExtent);
					return d1 - d2;
				};
				extents.sort(sortFunc);
			}

			console.log("Extents " + extents.length + "  "
					+ (+new Date() - time) + "ms");

			var surfobj = {
				geo : new WebMol.Geometry(true),
				mat : mat,
				done : false,
				finished : false
			// also webgl initialized
			};
			var surfid = surfaces.length;
			surfaces[surfid] = surfobj;
			var reducedAtoms = [];
			// to reduce amount data transfered, just pass x,y,z,serial and elem
			for (i = 0, il = atomlist.length; i < il; i++) {
				atom = atomlist[i];
				reducedAtoms[i] = {
					x : atom.x,
					y : atom.y,
					z : atom.z,
					serial : i,
					elem : atom.elem
				};
			}

			var sync = !!(WebMol.syncSurface);
			if (sync) { // don't use worker, still break up for memory purposes

				// to keep the browser from locking up, call through setTimeout
				var callSyncHelper = function callSyncHelper(i) {
					if (i >= extents.length)
						return;

					var VandF = generateMeshSyncHelper(type, extents[i].extent,
							extents[i].atoms, extents[i].toshow, reducedAtoms,
							totalVol);
					var mesh = generateSurfaceMesh(atomlist, VandF, mat);
					WebMol.mergeGeos(surfobj.geo, mesh);
					_viewer.render();

					setTimeout(callSyncHelper, 1, i + 1);
				}

				setTimeout(callSyncHelper, 1, 0);

				// TODO: Asynchronously generate geometryGroups (not separate
				// meshes) and merge them into a single geometry
			} else { // use worker

				var workers = [];
				if (type < 0)
					type = 0; // negative reserved for atom data
				for (i = 0, il = numWorkers; i < il; i++) {
					// var w = new Worker('webmol/SurfaceWorker.js');
					var w = new Worker(WebMol.SurfaceWorker);
					workers.push(w);
					w.postMessage({
						'type' : -1,
						'atoms' : reducedAtoms,
						'volume' : totalVol
					});
				}
				var cnt = 0;

				var rfunction = function(event) {
					var VandF = event.data;
					var mesh = generateSurfaceMesh(atomlist, VandF, mat);
					WebMol.mergeGeos(surfobj.geo, mesh);
					_viewer.render();
					console.log("async mesh generation " + (+new Date() - time)
							+ "ms");
					cnt++;
					if (cnt == extents.length)
						surfobj.done = true;
				};

				var efunction = function(event) {
					console.log(event.message + " (" + event.filename + ":"
							+ event.lineno + ")");
				};

				for (i = 0; i < extents.length; i++) {
					var worker = workers[i % workers.length];
					worker.onmessage = rfunction;

					worker.onerror = efunction;

					worker.postMessage({
						'type' : type,
						'expandedExtent' : extents[i].extent,
						'extendedAtoms' : extents[i].atoms,
						'atomsToShow' : extents[i].toshow
					});
				}
			}

			// NOTE: This is misleading if 'async' mesh generation - returns
			// immediately
			console.log("full mesh generation " + (+new Date() - time) + "ms");

			return surfid;
		};

		// set the material to something else, must render change
		this.setSurfaceMaterialStyle = function(surf, style) {
			if (surfaces[surf]) {
				surfaces[surf].mat = getMatWithStyle(style);
				surfaces[surf].mat.side = WebMol.FrontSide;
				surfaces[surf].finished = false; // trigger redraw
			}
		};

		// given the id returned by surfid, remove surface
		this.removeSurface = function(surf) {
			if (surfaces[surf] && surfaces[surf].lastGL) {
				if (surfaces[surf].geo !== undefined)
					surfaces[surf].geo.dispose();
				if (surfaces[surf].mat !== undefined)
					surfaces[surf].mat.dispose();
				modelGroup.remove(surfaces[surf].lastGL); // remove from scene
			}
			delete surfaces[surf];
			show();
		};
		
		this.removeAllSurfaces = function() {
			for(var i = 0; i < surfaces.length; i++) {
				if (surfaces[i] && surfaces[i].lastGL) {
					if (surfaces[i].geo !== undefined)
						surfaces[i].geo.dispose();
					if (surfaces[i].mat !== undefined)
						surfaces[i].mat.dispose();
					modelGroup.remove(surfaces[i].lastGL); // remove from scene
				}
				delete surfaces[i];
			}
			show();
		};

		// return jmol moveto command to position this scene
		this.jmolMoveTo = function() {
			var pos = modelGroup.position;
			// center on same position
			var ret = "center { " + (-pos.x) + " " + (-pos.y) + " " + (-pos.z)
					+ " }; ";
			// apply rotation
			var q = rotationGroup.quaternion;
			ret += "moveto .5 quaternion { " + q.x + " " + q.y + " " + q.z
					+ " " + q.w + " };";
			// zoom is tricky.. maybe i would be best to let callee zoom on
			// selection?
			// can either do a bunch of math, or maybe zoom to the center with a
			// fixed
			// but reasonable percentage

			return ret;
		};

		this.clear = function() {
			this.removeAllSurfaces();
			this.removeAllModels();
			this.removeAllLabels();
			this.removeAllShapes();
			show();
		};

		// props is a list of objects that select certain atoms and enumerate
		// properties for those atoms
		this.mapAtomProperties = function(props) {
			var atoms = getAtomsFromSel({});
			for (var a = 0, numa = atoms.length; a < numa; a++) {
				var atom = atoms[a];
				for (var i = 0, n = props.length; i < n; i++) {
					var prop = props[i];
					if (prop.props) {
						for ( var p in prop.props) {
							if (prop.props.hasOwnProperty(p)) {
								// check the atom
								if (atomIsSelected(atom, prop)) {
									if (!atom.properties)
										atom.properties = {};
									atom.properties[p] = prop.props[p];
								}
							}
						}
					}
				}
			}
		};

		var getModelGroup = function() {
			return modelGroup;
		};

		try {
			if (typeof (callback) === "function")
				callback(this);
		} catch (e) {
			// errors in callback shouldn't invalidate the viewer
			console.log("error with glviewer callback: " + e);
		}
	}

	return GLViewer;

})();

WebMol['glmolViewer'] = WebMol.GLViewer;
