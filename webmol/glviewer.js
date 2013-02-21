//a molecular viewer based on GLMol

var WebMol = WebMol || {};

// a webmol unified interace to gmol
WebMol.glmolViewer = (function() {
	// private class variables

	var aaScale = 1; // or 2



	// private class helper functions
	
	//computes the bounding box around the provided atoms
	var getExtent = function(atomlist) {
		var xmin = ymin = zmin = 9999;
		var xmax = ymax = zmax = -9999;
		var xsum = ysum = zsum = cnt = 0;

		if (atomlist.length == 0)
			return [ [ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ] ];
		for ( var i = 0; i < atomlist.length; i++) {
			var atom = atomlist[i];
			if (atom == undefined)
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

	// The constructor
	function GLViewer(element, width, height, callback) {
		// check dependencies
		if (typeof (THREE) === "undefined") {
			// three.js not loaded, take matters into our own hands
			throw "Missing Three.js";
		}
		if (typeof (GLmol) === "undefined") {
			throw "Missing GLmol.js";
		}

		// set variables
		var container = element;
		var id = container.id;
		
		var models = []; //atomistic molecular models

		var WIDTH = container.width() * aaScale;
		var HEIGHT = container.height() * aaScale;
		//set dimensions
		$(container).width(WIDTH);
		$(container).height(HEIGHT);
		
		var ASPECT = WIDTH / HEIGHT;
		var NEAR = 1, FAR = 800;
		var CAMERA_Z = -150;
		var renderer = new THREE.WebGLRenderer({
			antialias : true
		});
		renderer.sortObjects = false; // hopefully improve performance
		// 'antialias: true' now works in Firefox too!
		// setting this.aaScale = 2 will enable antialias in older Firefox but
		// GPU load increases.
		renderer.domElement.style.width = "100%";
		renderer.domElement.style.height = "100%";
		container.append(renderer.domElement);
		renderer.setSize(WIDTH, HEIGHT);

		var camera = new THREE.PerspectiveCamera(20, ASPECT, 1, 800);
		camera.position = new TV3(0, 0, CAMERA_Z);
		camera.lookAt(new TV3(0, 0, 0));
		var perspectiveCamera = camera;
		var orthoscopicCamera = new THREE.OrthographicCamera();
		orthoscopicCamera.position.z = CAMERA_Z;
		orthoscopicCamera.lookAt(new TV3(0, 0, 0));

		var self = this;
		$(window).resize(function() { // only window can capture resize event
			self.WIDTH = self.container.width() * self.aaScale;
			self.HEIGHT = self.container.height() * self.aaScale;
			self.ASPECT = self.WIDTH / self.HEIGHT;
			self.renderer.setSize(self.WIDTH, self.HEIGHT);
			self.camera.aspect = self.ASPECT;
			self.camera.updateProjectionMatrix();
			self.show();
		});

		var scene = null;
		var rotationGroup = null; // which contains modelGroup
		var modelGroup = null;

		var bgColor = 0x000000;
		var fov = 20;
		var fogStart = 0.4;
		var slabNear = -50; // relative to the center of rotationGroup
		var slabFar = +50;

		// UI variables
		var cq = new THREE.Quaternion(1, 0, 0, 0);
		var dq = new THREE.Quaternion(1, 0, 0, 0);
		var isDragging = false;
		var mouseStartX = 0;
		var mouseStartY = 0;
		var currentModelPos = 0;
		var cz = 0;
		var cslabNear = 0;
		var cslabFar = 0;
		
		var setSlabAndFog = function() {
			var center = rotationGroup.position.z - camera.position.z;
			if (center < 1)
				center = 1;
			camera.near = center + slabNear;
			if (camera.near < 1)
				camera.near = 1;
			camera.far = center + slabFar;
			if (camera.near + 1 > camera.far)
				camera.far = camera.near + 1;
			if (camera instanceof THREE.PerspectiveCamera) {
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
		
		//display scene
		var show = function() {
			if (!scene)
				return;

			// var time = new Date();
			setSlabAndFog();
			renderer.render(scene, camera);
			// console.log("rendered in " + (+new Date() - time) + "ms");
		};

		var initializeScene = function() {
			// CHECK: Should I explicitly call scene.deallocateObject?
			scene = new THREE.Scene();
			scene.fog = new THREE.Fog(bgColor, 100, 200);

			modelGroup = new THREE.Object3D();
			rotationGroup = new THREE.Object3D();
			rotationGroup.useQuaternion = true;
			rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
			rotationGroup.add(modelGroup);

			scene.add(rotationGroup);
			
			//setup lights
			var directionalLight = new THREE.DirectionalLight(0xFFFFFF);
			directionalLight.position = new TV3(0.2, 0.2, -1).normalize();
			directionalLight.intensity = 1.2;
			scene.add(directionalLight);
			var ambientLight = new THREE.AmbientLight(0x202020);
			scene.add(ambientLight);
		};

		initializeScene();

		
		// enable mouse support
		var glDOM = $(renderer.domElement);

		// TODO: Better touch panel support.
		// Contribution is needed as I don't own any iOS or Android device
		// with
		// WebGL support.
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
			if (x == undefined)
				return;
			isDragging = true;
			mouseButton = ev.which;
			mouseStartX = x;
			mouseStartY = y;
			cq = rotationGroup.quaternion;
			cz = rotationGroup.position.z;
			currentModelPos = modelGroup.position.clone();
			cslabNear = slabNear;
			cslabFar = slabFar;
		});

		glDOM.bind('DOMMouseScroll mousewheel',
				function(ev) { // Zoom
			ev.preventDefault();
			if (!scene)
				return;
			var scaleFactor = (rotationGroup.position.z - CAMERA_Z) * 0.85;
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

		glDOM.bind('mousemove touchmove',
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
			if (x == undefined)
				return;
			var dx = (x - mouseStartX) / WIDTH;
			var dy = (y - mouseStartY) / HEIGHT;
			var r = Math.sqrt(dx * dx + dy * dy);
			if (mode == 3
					|| (mouseButton == 3 && ev.ctrlKey)) { // Slab
				slabNear = cslabNear + dx * 100;
				slabFar = cslabFar + dy * 100;
			} else if (mode == 2 || mouseButton == 3
					|| ev.shiftKey) { // Zoom
				var scaleFactor = (rotationGroup.position.z - CAMERA_Z) * 0.85;
				if (scaleFactor < 80)
					scaleFactor = 80;
				rotationGroup.position.z = cz - dy* scaleFactor;
			} else if (mode == 1 || mouseButton == 2
					|| ev.ctrlKey) { // Translate
				var scaleFactor = (rotationGroup.position.z - CAMERA_Z) * 0.85;
				if (scaleFactor < 20)
					scaleFactor = 20;
				var translationByScreen = new TV3(-dx
						* scaleFactor, -dy * scaleFactor, 0);
				var q = rotationGroup.quaternion;
				var qinv = new THREE.Quaternion(q.x, q.y, q.z,
						q.w).inverse().normalize();
				var translation = translationByScreen.applyQuaternion(qinv);
				modelGroup.position.x = currentModelPos.x + translation.x;
				modelGroup.position.y = currentModelPos.y + translation.y;
				modelGroup.position.z = currentModelPos.z + translation.z;
			} else if ((mode == 0 || mouseButton == 1) && r != 0) { // Rotate
				var rs = Math.sin(r * Math.PI) / r;
				dq.x = Math.cos(r * Math.PI);
				dq.y = 0;
				dq.z = rs * dx;
				dq.w = rs * dy;
				rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
				rotationGroup.quaternion.multiply(dq);
				rotationGroup.quaternion.multiply(cq);
			}
			show();
		});
		

		// public methods
		this.setBackgroundColor = function(hex, a) {
			a = a | 1.0;
			bgColor = hex;
			renderer.setClearColorHex(hex, a);
			scene.fog.color = new TCo(hex);
			show();
		};

		this.setWidth = function(w) {
			WIDTH = w;
			$(htmlElement).width(WIDTH);
			renderer.setSize(WIDTH, HEIGHT);
		};

		this.setHeight = function(h) {
			HEIGHT = h;
			$(htmlElement).height(h);
			renderer.setSize(WIDTH, HEIGHT);
		};


		// return specified model
		this.getModel = function(id) {
			return models[id];
		};

		this.getView = function() {
			if (!modelGroup)
				return [ 0, 0, 0, 0, 0, 0, 0, 1 ];
			var pos = modelGroup.position;
			var q = rotationGroup.quaternion;
			return [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y,
					q.z, q.w ];
		}

		this.setView = function(arg) {
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
			show();
		}
		
		// apply styles, models, etc in viewer
		this.render = function() {

			var view = this.getView();
			initializeScene();
			
			for(var i = 0; i < models.length; i++) {
				modelGroup.add(models[i].globj());
			}				
			
			this.setView(view);
		};
		
		//zoom to atom selection
		this.zoomTo = function(sel) {
			var atoms = [];
			for(var i = 0; i < models.length; i++) {
				if(models[i]) {
					atoms = atoms.concat(models[i].selectedAtoms(sel));
				}
			}
			var tmp = getExtent(atoms);
			var center = new TV3(tmp[2][0], tmp[2][1], tmp[2][2]);
			modelGroup.position = center.multiplyScalar(-1);
			var x = tmp[1][0] - tmp[0][0], y = tmp[1][1] - tmp[0][1], z = tmp[1][2]
					- tmp[0][2];

			var maxD = Math.sqrt(x * x + y * y + z * z);
			if (maxD < 25)
				maxD = 25;

			slabNear = -maxD / 1.9;
			slabFar = maxD / 3;
			

			rotationGroup.position.z = maxD * 0.35
					/ Math.tan(Math.PI / 180.0 * camera.fov / 2) - 150;
			rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
		};
		
		//given molecular data and its format (pdb, sdf or xyz)
		//create a model and add it, returning the model identifier
		this.addModel = function(data, format) {
			var m = new WebMol.GLModel(models.length);
			m.addMolData(data, format);
			models.push(m);
			return m;
		};
		
		this.removeModel = function(model) {
			delete models[m.getID()];
			//clear off back of model array
			while(model.length > 0 && typeof(models[models.length-1]) === "undefined")
				models.pop();
		}

		//add a surface
		this.addSurface = function(type, style, atomsel, allsel) {
			
		}

		this.removeSurface = function(surf) {
			
		}
	}

	return GLViewer;
})();
