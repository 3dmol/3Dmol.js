//a molecular viewer based on GLMol


/**
 * WebGL-based 3Dmol.js viewer
 * Note: The preferred method of instantiating a GLViewer is through {@link $3Dmol.createViewer} 
 * 
 * @constructor 
 * @param {Object} element HTML element within which to create viewer
 * @param {ViewerSpec} config Object containing optional configuration for the viewer
 */
$3Dmol.GLViewer = (function() {
    // private class variables
    var numWorkers = 4; // number of threads for surface generation
    var maxVolume = 64000; // how much to break up surface calculations

    // private class helper functions

    function GLViewer(element, config) { 
        // set variables
        config = config || {};
        var callback = config.callback;
        var defaultcolors = config.defaultcolors;       
        if(!defaultcolors)
            defaultcolors = $3Dmol.elementColors.defaultColors;
        var nomouse = config.nomouse;
        var bgColor = 0;
        config.backgroundColor = config.backgroundColor || "#ffffff";
       //config.disableFog= config.disableFog || false;
        if(typeof(config.backgroundColor) != undefined) {
            bgColor = $3Dmol.CC.color(config.backgroundColor).getHex();
        }
        config.backgroundAlpha = config.backgroundAlpha == undefined ? 1.0 : config.backgroundAlpha;
        

        var camerax = 0;
        if(typeof(config.camerax) != undefined) {
            camerax = parseFloat(config.camerax);
        }
        var _viewer = this;
        var container = $(element); //we expect container to be jquery
        var glDOM = null;

        var models = []; // atomistic molecular models
        var surfaces = {};
        var shapes = []; // Generic shapes
        var labels = [];
        var fixed_labels = [];
        var clickables = []; //things you can click on
        var hoverables = []; //things you can hover over
        var current_hover = null;
        var hoverDuration = 500;
        var viewer_frame = 0;
        if(config.hoverDuration != undefined) {
            hoverDuration = config.hoverDuration;
        }
        if(config.antialias === undefined) config.antialias = true;
        if(config.cartoonQuality === undefined) config.cartoonQuality = 5;
        
        //reimplement jquery getwidth/height        
        var getRect = function() {
          let div = container[0];
          let rect = div.getBoundingClientRect();
          if(rect.width == 0 && rect.height == 0 && div.style.display === 'none' ) {
            let oldpos = div.style.position;
            let oldvis = div.style.visibility;
            div.style.display = 'block';
            div.style.visibility = 'hidden';
            div.style.position = 'absolute';
            rect = div.getBoundingClientRect();
            div.style.display = 'none';
            div.style.visibility = oldvis;
            div.style.position = oldpos;
          }
          return rect;          
        };
        
        var getWidth = function() {
          return getRect().width;
        };
        
        var getHeight = function() {
          return getRect().height;
        };
        
        var WIDTH = getWidth();
        var HEIGHT = getHeight();

        var viewChangeCallback = null;
        var stateChangeCallback = null;
        
        var NEAR = 1, FAR = 800;
        var CAMERA_Z = 150;
        var fov = 20;

        var linkedViewers = [];
        var renderer = new $3Dmol.Renderer({
            antialias : config.antialias,
            preserveDrawingBuffer: true, //so we can export images
            premultipliedAlpha : false,/* more traditional compositing with background */
            id:config.id,
            row:config.row,
            col:config.col,
            rows:config.rows,
            cols:config.cols,
            canvas:config.canvas,
            containerWidth:WIDTH,
            containerHeight:HEIGHT,
        });
        renderer.domElement.style.width = "100%";
        renderer.domElement.style.height = "100%";
        renderer.domElement.style.padding = "0";
        renderer.domElement.style.position = "absolute"; //TODO: get rid of this
        renderer.domElement.style.top = "0px";
        renderer.domElement.style.left = "0px";
        renderer.domElement.style.zIndex = "0";       
        
        var row = config.row;
        var col = config.col;
        var cols = config.cols;
        var rows = config.rows;
        var viewers = config.viewers;
        var control_all = config.control_all;

        var ASPECT =renderer.getAspect(WIDTH,HEIGHT);

        var camera = new $3Dmol.Camera(fov, ASPECT, NEAR, FAR, config.orthographic);
        camera.position = new $3Dmol.Vector3(camerax, 0, CAMERA_Z);
        var lookingAt = new $3Dmol.Vector3();
        camera.lookAt(lookingAt);

        var raycaster = new $3Dmol.Raycaster(new $3Dmol.Vector3(0, 0, 0),
                new $3Dmol.Vector3(0, 0, 0));
        var projector = new $3Dmol.Projector();

        var scene = null;
        var rotationGroup = null; // which contains modelGroup
        var modelGroup = null;

        var fogStart = 0.4;
        var slabNear = -50; // relative to the center of rotationGroup
        var slabFar = 50;

        // UI variables
        var cq = new $3Dmol.Quaternion(0, 0, 0, 1);
        var dq = new $3Dmol.Quaternion(0, 0, 0, 1);
        var animated = 0;
        var animationTimers = new Set();
        var isDragging = false;
        var mouseStartX = 0;
        var mouseStartY = 0;
        var touchDistanceStart = 0;
        var currentModelPos = 0;
        var cz = 0;
        var cslabNear = 0;
        var cslabFar = 0;      
        
        var decAnim = function() {
            //decrement the number of animations currently
            animated--;
            if(animated < 0) animated = 0;
        };
        var incAnim = function() {
            animated++;
        };
        var nextSurfID = function() {
            //compute the next highest surface id directly from surfaces
            //this is necessary to support linking of model data
            var max = 0;
            for (var i in surfaces) { // this is an object with possible holes
                if(!surfaces.hasOwnProperty(i)) continue;
                var val = parseInt(i);
                if(!isNaN(val)) i = val;
                if(i > max) max = i;
            }
            return max+1;
        };

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

            camera.fov = fov;
            camera.right = center * Math.tan(Math.PI / 180 * fov);
            camera.left = -camera.right;
            camera.top = camera.right / ASPECT;
            camera.bottom = -camera.top;
            
            camera.updateProjectionMatrix();

            scene.fog.near = camera.near + fogStart * (camera.far - camera.near);
            // if (scene.fog.near > center) scene.fog.near = center;
            scene.fog.far = camera.far;
            
            if(config.disableFog){
                scene.fog.near=scene.fog.far;
            }
        };
        
        // display scene
        //if nolink is set/true, don't propagate changes to linked viewers
        var show = function(nolink) {
            renderer.setViewport();
            if (!scene)
                return;
            // var time = new Date();
            setSlabAndFog();
            renderer.render(scene, camera);
            // console.log("rendered in " + (+new Date() - time) + "ms");
            
            //have any scene change trigger a callback
            if(viewChangeCallback) viewChangeCallback(_viewer.getView());

            if(!nolink && linkedViewers.length > 0) {                
                var view = _viewer.getView();
                for(var i = 0; i < linkedViewers.length; i++) {
                    var other = linkedViewers[i];
                    other.setView(view, true);
                }
            }
        };

        var initializeScene = function() {

            scene = new $3Dmol.Scene();
            scene.fog = new $3Dmol.Fog(bgColor, 100, 200);

            modelGroup = new $3Dmol.Object3D();
            rotationGroup = new $3Dmol.Object3D();
            rotationGroup.useQuaternion = true;
            rotationGroup.quaternion = new $3Dmol.Quaternion(0, 0, 0, 1);
            rotationGroup.add(modelGroup);

            scene.add(rotationGroup);

            // setup lights
            var directionalLight = new $3Dmol.Light(0xFFFFFF);
            directionalLight.position = new $3Dmol.Vector3(0.2, 0.2, 1)
                    .normalize();
            directionalLight.intensity = 1.0;
            scene.add(directionalLight);
        };

        initializeScene();
        renderer.setClearColorHex(bgColor, config.backgroundAlpha);
        scene.fog.color = $3Dmol.CC.color(bgColor);

        var clickedAtom = null;

        // enable mouse support

        //regenerate the list of clickables
        //also updates hoverables
        var updateClickables = function() {
            clickables.splice(0,clickables.length);
            hoverables.splice(0,hoverables.length);
            
            for (let i = 0, il = models.length; i < il; i++) {
                var model = models[i];
                if(model) {
                    let atoms = model.selectedAtoms({
                        clickable : true
                    });
                    
                    let hoverable_atoms = model.selectedAtoms({
                        hoverable : true
                    });
                    // Array.prototype.push.apply(hoverables,hoverable_atoms);
                    for (let n = 0; n < hoverable_atoms.length; n++) {
                        hoverables.push(hoverable_atoms[n]);
                    }

                    // Array.prototype.push.apply(clickables, atoms); //add atoms into clickables
                    for (let m = 0; m < atoms.length; m++) {
                        clickables.push(atoms[m]);
                    }
                    
                }
            }
            for (let i = 0, il = shapes.length; i < il; i++) {

                let shape = shapes[i];
                if (shape && shape.clickable) {
                    clickables.push(shape);
                }
                if( shape && shape.hoverable){
                    hoverables.push(shape);
                }
            }
        };
        
         /**
         * Return a list of objects that intersect that at the specified viewer position.
         * 
         * @function $3Dmol.GLViewer#targetedObjects
         * @param {x} - x position in screen coordinates
         * @param {y} - y position in screen coordinates
         * @param {objects} - list of objects or selection object specifying what object to check for targeting 
        */
        let targetedObjects = this.targetedObjects = function(x,y,objects) {
            var mouse = {
                x : x,
                y : y,
                z : -1.0
            };
            if(!Array.isArray(objects)) { //assume selection object
                objects = this.selectedAtoms(objects);
            }
            if(objects.length == 0) return [];
            raycaster.setFromCamera(mouse,camera);
            return raycaster.intersectObjects(modelGroup, objects);
        };
        
        //return offset of container
        var canvasOffset = function() {
          let canvas = glDOM.get(0);
          let rect = canvas.getBoundingClientRect();
          let doc = canvas.ownerDocument;
          let docElem = doc.documentElement;
          let win = doc.defaultView;
          return {
            top: rect.top + win.pageYOffset - docElem.clientTop,
            left: rect.left + win.pageXOffset - docElem.clientLeft
          };
        };
        
        /** Convert model coordinates to screen coordinates.
         * @function $3Dmol.GLViewer#modelToScreen
         * @param {object | list} - an object or list of objects with x,y,z attributes (e.g. an atom)
         * @return {object | list} - and object or list of {x: screenX, y: screenY}
         */
        this.modelToScreen = function(coords) {
            let returnsingle = false;
            if(!Array.isArray(coords)) {
                coords = [coords];
                returnsingle = true;
            }
            
            let results = [];
            let offset = canvasOffset();
            coords.forEach(coord => {
                let t = new $3Dmol.Vector3(coord.x,coord.y,coord.z);
                t.applyMatrix4(modelGroup.matrixWorld);   
                projector.projectVector(t, camera);       
                let screenX = WIDTH*(t.x+1)/2.0+offset.left;
                let screenY = -HEIGHT*(t.y-1)/2.0+offset.top;
                results.push({x:screenX,y:screenY});
            }); 
            if(returnsingle) results = results[0];
            return results;
        };
        
        // Checks for selection intersects on mousedown
        var handleClickSelection = function(mouseX, mouseY, event) {            
            let intersects = targetedObjects(mouseX,mouseY,clickables);
            if (intersects.length) {
                var selected = intersects[0].clickable;
                if (selected.callback !== undefined) {
	                 if(typeof (selected.callback) != "function") {
	                    selected.callback = $3Dmol.makeFunction(selected.callback);
	                }
                    if(typeof (selected.callback) === "function") {
                        selected.callback(selected, _viewer, event, container);
                    }
                }
            }
        };
        
        //set current_hover to sel (which can be null), calling appropraite callbacks
        var setHover = function(selected, event) {
            if(current_hover == selected) return;
            if(current_hover) {
                if(typeof (current_hover.unhover_callback) != "function") {
                    current_hover.unhover_callback = $3Dmol.makeFunction(current_hover.unhover_callback);
                } 
                current_hover.unhover_callback(current_hover, _viewer, event, container);
            }
            current_hover=selected;
            
            if (selected && selected.hover_callback !== undefined) {
                if(typeof (selected.hover_callback) != "function") {
                    selected.hover_callback = $3Dmol.makeFunction(selected.hover_callback);
                }
                if(typeof (selected.hover_callback) === "function") {
                    selected.hover_callback(selected, _viewer, event, container);
                }
            }  
            
        };
        
        //checks for selection intersects on hover
        var handleHoverSelection = function(mouseX, mouseY){
            if(hoverables.length == 0) return;
            let intersects = targetedObjects(mouseX,mouseY,hoverables);
            if (intersects.length) {
                var selected = intersects[0].clickable;
                setHover(selected);
                current_hover=selected;
            }
            else{
                setHover(null);
            }
        };
        
        //sees if the mouse is still on the object that invoked a hover event and if not then the unhover callback is called
        var handleHoverContinue = function(mouseX,mouseY){
            let intersects = targetedObjects(mouseX,mouseY,hoverables);
            if(intersects.length == 0 || intersects[0] === undefined){
                setHover(null);
            }
            if(intersects[0]!== undefined && intersects[0].clickable !== current_hover){
                setHover(null);
            }
        };

        var calcTouchDistance = function(ev) { // distance between first two
                                                // fingers
            var xdiff = ev.originalEvent.targetTouches[0].pageX -
                    ev.originalEvent.targetTouches[1].pageX;
            var ydiff = ev.originalEvent.targetTouches[0].pageY -
                    ev.originalEvent.targetTouches[1].pageY;
            return Math.sqrt(xdiff * xdiff + ydiff * ydiff);
        };
        
        //check targetTouches as well
        var getX = function(ev) {
            var x = ev.pageX;
            if(x == undefined) x = ev.originalEvent.pageX; //firefox
            if (ev.originalEvent.targetTouches &&
                    ev.originalEvent.targetTouches[0]) {
                x = ev.originalEvent.targetTouches[0].pageX;
            }
            else if (ev.originalEvent.changedTouches &&
                    ev.originalEvent.changedTouches[0]) {
                x = ev.originalEvent.changedTouches[0].pageX;
            }            
            return x;
        };
        
        var getY = function(ev) {
            var y = ev.pageY;
            if(y == undefined) y = ev.originalEvent.pageY;
            if (ev.originalEvent.targetTouches &&
                    ev.originalEvent.targetTouches[0]) {
                y = ev.originalEvent.targetTouches[0].pageY;
            }
            else if (ev.originalEvent.changedTouches &&
                    ev.originalEvent.changedTouches[0]) {
                y = ev.originalEvent.changedTouches[0].pageY;
            }            
            return y;
        };
        
        /**
         * For a given screen (x,y) displacement return model displacement
         * @param{x} x displacement in screen coordinates
         * @param{y} y displacement in screen corodinates
         * @param{modelz} z coordinate in model coordinates to compute offset for, default is model axis 
         * @function $3Dmol.GLViewer#screenOffsetToModel
        */        
        var screenOffsetToModel = this.screenOffsetToModel = function(x,y,modelz) {
            var dx = x/WIDTH;
            var dy = y/HEIGHT;
            var zpos = (modelz === undefined ? rotationGroup.position.z : modelz); 
            var q = rotationGroup.quaternion;                        
            var t = new $3Dmol.Vector3(0,0,zpos);
            projector.projectVector(t, camera);
            t.x += dx*2;
            t.y -= dy*2;
            projector.unprojectVector(t, camera);
            t.z = 0;                            
            t.applyQuaternion(q);
            return t;
        };
        
        /**
         * Distance from screen coordinate to model coordinate assuming screen point
         * is projected to the same depth as model coordinate
         * @param{screen} xy screen coordinate
         * @param{model} xyz model coordinate
         * @function $3Dmol.GLViewer#screenToModelDistance
        */   
        this.screenToModelDistance = function(screen,model) {
            let offset = canvasOffset();
            
            //convert model to screen to get screen z                     
            let mvec = new $3Dmol.Vector3(model.x,model.y,model.z);
            mvec.applyMatrix4(modelGroup.matrixWorld);   
            let m = mvec.clone();
            projector.projectVector(mvec, camera);
            
            let t = new $3Dmol.Vector3((screen.x-offset.left)*2/WIDTH-1,(screen.y-offset.top)*2/-HEIGHT+1,mvec.z);
            projector.unprojectVector(t, camera);
            
            return t.distanceTo(m);
        };                         
        
        //for grid viewers, return true if point is in this viewer
        var isInViewer = function(x,y) {
            if(viewers != undefined && !control_all){
                var width = WIDTH/cols;
                var height = HEIGHT/rows;
                var offset = canvasOffset();
                var relx = (x - offset.left);
                var rely = (y - offset.top) ;
                    
                var r = rows-Math.floor(rely/height)-1;
                var c = Math.floor(relx/width);

                if(r != row || c != col)
                    return false;
            }
            return true;
        };

        // this event is bound to the body element, not the container,
        // so no need to put it inside initContainer()
        $('body').bind('mouseup touchend', function(ev) {
            // handle selection
            if(isDragging && scene) { //saw mousedown, haven't moved
                var x = getX(ev);
                var y = getY(ev);
                if(x == mouseStartX && y == mouseStartY) {
                    var offset = canvasOffset();
                    var mouseX = ((x - offset.left) / WIDTH) * 2 - 1;
                    var mouseY = -((y - offset.top) / HEIGHT) * 2 + 1;
                    handleClickSelection(mouseX, mouseY, ev, container);
                }
            }
            
            isDragging = false;

        });

        
        //if the user has specify zoom limits, readjust to fit within them
        //also, make sure we don't go past CAMERA_Z
        var adjustZoomToLimits = function(z) {
            //a lower limit of 0 is at CAMERA_Z
            if(config.lowerZoomLimit && config.lowerZoomLimit > 0) {
                var lower = CAMERA_Z-config.lowerZoomLimit;
                if(z > lower) z = lower;
            }
            
            if(config.upperZoomLimit && config.upperZoomLimit > 0) {
                var upper = CAMERA_Z-config.upperZoomLimit;
                if(z < upper) z = upper;
            }
            
            if(z > CAMERA_Z) {
                z = CAMERA_Z*0.999; //avoid getting stuck
            }
            return z;
        };
        
        /**
         * Set a callback to call when the view has potentially changed.
         * 
         * @function $3Dmol.GLViewer#setViewChangeCallback
        */
        this.setViewChangeCallback = function(callback) {
            if(typeof(callback) === 'function' || callback == null) 
                viewChangeCallback = callback;
        };

        /**
         * Set a callback to call when the view has potentially changed.
         * 
         * @function $3Dmol.GLViewer#setStateChangeCallback
        */
        this.setStateChangeCallback = function(callback) {
            if(typeof(callback) === 'function' || callback == null) 
                stateChangeCallback = callback;
        };

        /**
         * Return object representing internal state of 
         * the viewer appropriate for passing to setInternalState
         * 
         * @function $3Dmol.GLViewer#getInternalState
        */
        this.getInternalState = function() {
          var ret = {'models': [], 'surfaces': [], 'shapes': [], 'labels': [] };
          for (let i = 0; i < models.length; i++) {
            if (models[i]) {
              ret.models[i] = models[i].getInternalState();
            }
          }
          
          //todo: labels, shapes, surfaces
          
          return ret;
        };
        
        /**
         * Overwrite internal state of the viewer with passed  object
         * which should come from getInternalState.
         * 
         * @function $3Dmol.GLViewer#setInternalState
        */
        this.setInternalState = function(state) {
          
          //clear out current viewer
          this.clear();
          
          //set model state
          var newm = state.models;          
          for(let i = 0; i < newm.length; i++) {
            if(newm[i]) {
              models[i] = new $3Dmol.GLModel(i);
              models[i].setInternalState(newm[i]);
            }
          }
          
          //todo: labels, shapes, surfaces
          this.render();
        };
                                             
        /**
         * Set lower and upper limit stops for zoom.
         * 
         * @function $3Dmol.GLViewer#setZoomLimits
         * @param {lower} - limit on zoom in (positive number).  Default 0.
         * @param {upper} - limit on zoom out (positive number).  Default infinite.
         * @example
          $.get("data/set1_122_complex.mol2", function(moldata) {
                var m = viewer.addModel(moldata);
                viewer.setStyle({stick:{colorscheme:"Jmol"}});
                viewer.setZoomLimits(100,200);
                viewer.zoomTo();
                viewer.zoom(10); //will not zoom all the way
                viewer.render();
            });
        */
        this.setZoomLimits = function(lower, upper) {
            if(typeof(lower) !== 'undefined') config.lowerZoomLimit = lower;
            if(upper) config.upperZoomLimit = upper;
            rotationGroup.position.z = adjustZoomToLimits(rotationGroup.position.z);
            show();
        };

        /**
         * Set camera parameters (distance to the origin and field of view)
         * 
         * @function $3Dmol.GLViewer#setCameraParameters
         * @param {parameters} - new camera parameters, with possible fields 
         *                       being fov for the field of view and z for the 
         *                       distance to the origin.
         * @example
          $.get("data/set1_122_complex.mol2", function(data) {
                var m = viewer.addModel(data);
                viewer.setStyle({stick:{}});
                viewer.zoomTo();
                viewer.setCameraParameters({ fov: 10 , z: 300 });
                viewer.render();
            });
        */
        this.setCameraParameters = function(parameters) {
            if (parameters.fov !== undefined) {
                fov = parameters.fov;
                camera.fov = fov;
            }

            if (parameters.z !== undefined) {
                CAMERA_Z = parameters.z;
                camera.z = CAMERA_Z;
            }
        };

        var mouseButton;
        var _handleMouseDown = this._handleMouseDown = function(ev) {
            ev.preventDefault();
            if (!scene)
                return;
            var x = getX(ev);
            var y = getY(ev);
            if (x === undefined)
                return;
            isDragging = true;
            clickedAtom = null;
            mouseButton = ev.which;
            mouseStartX = x;
            mouseStartY = y;
            touchDistanceStart = 0;
            if (ev.originalEvent.targetTouches &&
                    ev.originalEvent.targetTouches.length == 2) {
                touchDistanceStart = calcTouchDistance(ev);
            }
            cq = rotationGroup.quaternion.clone();
            cz = rotationGroup.position.z;
            currentModelPos = modelGroup.position.clone();
            cslabNear = slabNear;
            cslabFar = slabFar;
           
        };
        
        var _handleMouseScroll  = this._handleMouseScroll = function(ev) { // Zoom
            ev.preventDefault();
            if (!scene)
                return;

            var x = getX(ev);
            var y = getY(ev);
            if (x === undefined)
                return;
            if(!isInViewer(x,y)) {
                return;
            }

            var scaleFactor = (CAMERA_Z - rotationGroup.position.z) * 0.85;
            var mult = 1.0;
            if(ev.originalEvent.ctrlKey) {
                mult = -1.0; //this is a pinch event turned into a wheel event (or they're just holding down the ctrl)
            }
            if (ev.originalEvent.detail) { 
                rotationGroup.position.z += mult * scaleFactor * ev.originalEvent.detail / 10;
            } else if (ev.originalEvent.wheelDelta) { 
                rotationGroup.position.z -= mult * scaleFactor * ev.originalEvent.wheelDelta / 400;
            }
            rotationGroup.position.z = adjustZoomToLimits(rotationGroup.position.z);            
            show();            
        };        
        /**
         * Return image URI of viewer contents (base64 encoded).
         * @function $3Dmol.GLViewer#pngURI
         * 
         */
        this.pngURI = function() {
            return this.getCanvas().toDataURL('image/png');
        };
        
        /**
         * Return underlying canvas element.         
         */
        this.getCanvas = function() {
            return glDOM.get(0);
        };
        
        /**
         * Return renderer element.
         */
        this.getRenderer = function() {
            return renderer;
        };
      /**
           * Set the duration of the hover delay
           * 
           * @function $3Dmol.GLViewer#setHoverDuration
           * @param {number}
           *            [hoverDuration] - an optional parameter that denotes
           *            the duration of the hover delay (in milliseconds) before the hover action is called
           * 
       */
        this.setHoverDuration = function(duration) {
            hoverDuration = duration;
        };
        
        var hoverTimeout;
        var _handleMouseMove = this._handleMouseMove = function(ev) { // touchmove

            clearTimeout(hoverTimeout);
            var offset = canvasOffset();
            var mouseX = ((getX(ev) - offset.left) / WIDTH) * 2 - 1;
            var mouseY = -((getY(ev) - offset.top) / HEIGHT) * 2 + 1;
            
            // hover timeout
            if(current_hover !== null) {
                handleHoverContinue(mouseX,mouseY,ev);
            }
            
            if(hoverables.length > 0) {
                hoverTimeout=setTimeout(
                        function(){
                            handleHoverSelection(mouseX,mouseY);
                        },
                    hoverDuration);
            }

            ev.preventDefault();
            if (!scene)
                return;
            if (!isDragging)
                return;
            var mode = 0;

            var x = getX(ev);
            var y = getY(ev);
            if (x === undefined)
                return;

            if(!isInViewer(x,y)) {
                return;
            }


            var dx = (x - mouseStartX) / WIDTH;
            var dy = (y - mouseStartY) / HEIGHT;
            // check for pinch
            if (touchDistanceStart != 0 &&
                    ev.originalEvent.targetTouches &&
                    ev.originalEvent.targetTouches.length == 2) {
                var newdist = calcTouchDistance(ev);
                // change to zoom
                mode = 2;
                dy = (newdist - touchDistanceStart) * 2 / (WIDTH + HEIGHT);
            } else if (ev.originalEvent.targetTouches &&
                    ev.originalEvent.targetTouches.length == 3) {
                // translate
                mode = 1;
            }
            var ratioX = renderer.getXRatio();
            var ratioY = renderer.getYRatio();
            dx*=ratioX;
            dy*=ratioY;
            var r = Math.sqrt(dx * dx + dy * dy);
            var scaleFactor;
            if (mode == 3 || (mouseButton == 3 && ev.ctrlKey)) { // Slab
                slabNear = cslabNear + dx * 100;
                slabFar = cslabFar - dy * 100;
            } else if (mode == 2 || mouseButton == 3 || ev.shiftKey) { // Zoom
                scaleFactor = (CAMERA_Z - rotationGroup.position.z) * 0.85;
                if (scaleFactor < 80)
                    scaleFactor = 80;
                rotationGroup.position.z = cz + dy * scaleFactor;
                rotationGroup.position.z = adjustZoomToLimits(rotationGroup.position.z); 
            } else if (mode == 1 || mouseButton == 2 || ev.ctrlKey) { // Translate
                var t = screenOffsetToModel(ratioX*(x-mouseStartX), ratioY*(y-mouseStartY));
                modelGroup.position.addVectors(currentModelPos,t);
                
            } else if ((mode === 0 || mouseButton == 1) && r !== 0) { // Rotate
                var rs = Math.sin(r * Math.PI) / r;
                dq.x = Math.cos(r * Math.PI);
                dq.y = 0;
                dq.z = rs * dx;
                dq.w = -rs * dy;
                rotationGroup.quaternion.set(1, 0, 0, 0);
                rotationGroup.quaternion.multiply(dq);
                rotationGroup.quaternion.multiply(cq);
            }
            show();
        };
        
        var initContainer = function(element) {
            container = element;
            WIDTH = getWidth();
            HEIGHT = getHeight();
            ASPECT = renderer.getAspect(WIDTH,HEIGHT);
            renderer.setSize(WIDTH, HEIGHT);
            container.append(renderer.domElement);
            glDOM = $(renderer.domElement);

            if (!nomouse) {
                // user can request that the mouse handlers not be installed
                glDOM.bind('mousedown touchstart', _handleMouseDown);
                glDOM.bind('DOMMouseScroll mousewheel', _handleMouseScroll);
                glDOM.bind('mousemove touchmove', _handleMouseMove);
                
                glDOM.bind("contextmenu", function(ev) {
                    ev.preventDefault();
                });
                
            }
        };
        initContainer(container);

        // public methods
        /**
         * Change the viewer's container element 
         * Also useful if the original container element was removed from the DOM.
         * 
         * @function $3Dmol.GLViewer#setContainer
         *
         * @param {Object | string} element
         *            Either HTML element or string identifier. Defaults to the element used to initialize the viewer.

         */
        this.setContainer = function(element) {
            if(typeof(element) === "string")
                element = $("#"+element);
            if(!element) {
                element = container;
            }
            initContainer(element);
            return this;
        };
        
        /**
         * Set the background color (default white)
         * 
         * @function $3Dmol.GLViewer#setBackgroundColor
         * @param {number}
         *            hex Hexcode specified background color, or standard color spec
         * @param {number}
         *            a Alpha level (default 1.0)
         * 
         * @example
         * 
         * viewer.setBackgroundColor(0x000000);


         * 
         */
        this.setBackgroundColor = function(hex, a) {
            if(typeof(a) == "undefined") {
                a = 1.0;
            }
            else if(a < 0 || a > 1.0) {
                a = 1.0;
            }
            var c = $3Dmol.CC.color(hex);
            scene.fog.color = c;
            bgColor = c.getHex();
            renderer.setClearColorHex(c.getHex(), a);
            show();
            
            return this;
        };
        
        /**
         * Set view projection scheme.  Either orthographic or perspective.  
         * Default is perspective.  Orthographic can also be enabled on viewer creation
         * by setting orthographic to true in the config object.
         * 
         * @function $3Dmol.GLViewer#setProjection
         * 
         * @example
         viewer.setViewStyle({style:"outline"});
              $.get('data/1fas.pqr', function(data){
                  viewer.addModel(data, "pqr");
                  $.get("data/1fas.cube",function(volumedata){
                      viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
                  });
                  viewer.zoomTo();

                  viewer.setProjection("orthographic");
                  viewer.render(callback);
              });
         * 
         */
        this.setProjection = function(proj) {
            camera.ortho = (proj === "orthographic");
            setSlabAndFog();            
        };
        
        /**
         * Set global view styles.  
         * @function $3Dmol.GLViewer#setViewStyle
         * 
         * @example
         *   viewer.setViewStyle({style:"outline"});
              $.get('data/1fas.pqr', function(data){
                  viewer.addModel(data, "pqr");
                  $.get("data/1fas.cube",function(volumedata){
                      viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
                  });
                  viewer.zoomTo();
                  viewer.render(callback);
              });
         * 
         */
         this.setViewStyle = function(parameters) {
            if (parameters.style === "outline") {
                var params = {};
                if(parameters.color) params.color =  $3Dmol.CC.color(parameters.color);
                if(parameters.width) params.width = parameters.width;
                renderer.enableOutline(params);
            } else {
                renderer.disableOutline();
            }           
            return this;
        };
         
        if(config.style) { //enable setting style in constructor
             this.setViewStyle(config);
        }

        /**
         * Set viewer width
         * 
         * @function $3Dmol.GLViewer#setWidth
         * @param {number}
         *            w Width in pixels
         */
        this.setWidth = function(w) {
            WIDTH = w || WIDTH;
            renderer.setSize(WIDTH, HEIGHT);
            return this;
        };

        /**
         * Set viewer height
         * 
         * @function $3Dmol.GLViewer#setHeight
         * @param {number}
         *            h Height in pixels
         */
        this.setHeight = function(h) {
            HEIGHT = h || HEIGHT;
            renderer.setSize(WIDTH, HEIGHT);
            return this;
        };

        /**
         * Resize viewer according to containing HTML element's dimensions
         * 
         * @function $3Dmol.GLViewer#resize
         */
        this.resize = function() {
            WIDTH = getWidth();
            HEIGHT = getHeight();
            ASPECT = renderer.getAspect(WIDTH,HEIGHT);
            renderer.setSize(WIDTH, HEIGHT);
            camera.aspect = ASPECT;
            camera.updateProjectionMatrix();
            show();
            return this;
        };

        $(window).resize(this.resize);
        
        if(typeof(window.ResizeObserver) !== undefined) {
            var divwatcher = new window.ResizeObserver(this.resize);
            divwatcher.observe(container[0]);
        }

        /**
         * Return specified model
         * 
         * @function $3Dmol.GLViewer#getModel
         * @param {number}
         *            [id=last model id] - Retrieve model with specified id
         * @default Returns last model added to viewer or null if there are no models
         * @return {GLModel}
         * 
         * @example // Retrieve reference to first GLModel added var m =
         *    $3Dmol.download("pdb:1UBQ",viewer,{},function(m1){
                  $3Dmol.download("pdb:1UBI", viewer,{}, function(m2) {
                    viewer.zoomTo();
                    m1.setStyle({cartoon: {color:'green'}});
                    //could use m2 here as well
                    viewer.getModel().setStyle({cartoon: {color:'blue'}});
                    viewer.render();
                })
              });     
         */
        this.getModel = function(id) {
            if(id === undefined) {
                return models.length == 0 ? null : models[models.length-1];
            } 
            if(id instanceof $3Dmol.GLModel) {
                return id;
            }
            if(!(id in models)) {
                if(models.length == 0) 
                    return null;
                else
                    return models[models.length-1]; //get last model if no (or invalid) id specified
            }
            return models[id];
        };
        
        //interpolate between two normalized quaternions (t between 0 and 1)
        //https://en.wikipedia.org/wiki/Slerp
        var slerp = function(v0, v1, t) {
            // Compute the cosine of the angle between the two vectors.
            //dot product
            if(t == 1) return v1;
            else if(t == 0) return v0;
            var dot = v0.x*v1.x+v0.y*v1.y+v0.z*v1.z+v0.w*v1.w;
            if (dot > 0.9995) {
                // If the inputs are too close for comfort, linearly interpolate
                // and normalize the result.
                var result = new $3Dmol.Quaternion(
                        v0.x+t*(v1.x-v0.x),
                        v0.y+t*(v1.y-v0.y),
                        v0.z+t*(v1.z-v0.z),
                        v0.w+t*(v1.w-v0.w));
                        
                result.normalize();
                return result;
            }

            // If the dot product is negative, the quaternions
            // have opposite handed-ness and slerp won't take
            // the shorted path. Fix by reversing one quaternion.
            if (dot < 0.0) {
                v1 = v1.clone().multiplyScalar(-1);
                dot = -dot;
            }  

            if(dot > 1) dot = 1.0;
            else if(dot < -1) dot = -1.0;

            var theta_0 = Math.acos(dot);  // theta_0 = angle between input vectors
            var theta = theta_0*t;    // theta = angle between v0 and result 

            var v2 = v1.clone();
            v2.sub(v0.clone().multiplyScalar(dot));
            v2.normalize();              // { v0, v2 } is now an orthonormal basis

            var c = Math.cos(theta);
            var s = Math.sin(theta);
            var ret = new $3Dmol.Quaternion(
                    v0.x*c+v2.x*s,
                    v0.y*c+v2.y*s,
                    v0.z*c+v2.z*s,
                    v0.w*c+v2.w*s
            );
            ret.normalize();
            return ret;
        };
        
        var spinInterval;
        /**
         * Continuously rotate a scene around the specified axis.
         *
         * Call `$3Dmol.GLViewer.spin(false)` to stop spinning.
         * 
         * @function $3Dmol.GLViewer#spin
         * @param {string}
         *            [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around.
         *            Default "y".  View relative (rather than model relative) axes are prefixed with v. 
         * @param {number}
         *            [speed] - Speed multiplier for spinning the viewer. 1 is default and a negative
         *             value reverses the direction of the spin.        
         *  
         */        
        this.spin = function(axis, speed){
            clearInterval(spinInterval);
            if(typeof axis == 'undefined')
                axis = 'y';
            if(typeof axis == "boolean"){
                if(!axis)
                    return;
                else
                    axis = 'y';
            }
            if(typeof speed != 'number'){
                speed = 1;
            }

            if(Array.isArray(axis)){
               axis = {x:axis[0],y:axis[1],z:axis[2]} ;
            }
            //out of bounds check

            var viewer = this;

            spinInterval = setInterval(
                function(){
                    viewer.rotate(1 * speed,axis);
                }, 25);            
            
        };         

        //animate motion between current position and passed position
        // can set some parameters to null
        //if fixed is true will enforce the request animation, otherwise
        //does relative updates
        //positions objects have modelggroup position, rotation group position.z,
        //and rotationgroup quaternion
        //return array includes final position, but not current 
        //the returned array includes an animate method
        var animateMotion = function(duration, fixed, mpos, rz, rot, cam) {
            var interval = 20;
            var steps = Math.ceil(duration/interval);
            if(steps < 1) steps = 1;
            incAnim();
            
            var curr = {mpos:modelGroup.position.clone(),
                    rz: rotationGroup.position.z,
                    rot: rotationGroup.quaternion.clone(),
                    cam: lookingAt.clone()};
            
            if(fixed) { //precompute path and stick to it
                steps = new Array(steps);
                var n = steps.length;
                for(var i = 0; i < n; i++) {
                    let frac = (i+1)/n;
                    let next = {mpos: curr.mpos, rz:curr.rz, rot:curr.rot};
                    if(mpos) {
                        next.mpos = mpos.clone().sub(curr.mpos).multiplyScalar(frac).add(curr.mpos);
                    }
                    if(typeof(rz) != 'undefined' && rz != null) {
                        next.rz = curr.rz+frac*(rz-curr.rz);
                    }
                    if(rot) {
                        next.rot = slerp(curr.rot,rot,frac);
                    }
                    if(cam) {
                        next.cam = cam.clone().sub(curr.cam).multiplyScalar(frac).add(curr.cam);
                    }
                    
                    steps[i] = next;
                }
                
                let step = 0;
                let callback = function() {
                    var p = steps[step];
                    step += 1;
                    if(p.mpos) {
                        modelGroup.position = p.mpos;
                    }
                    if(p.rz) {
                        rotationGroup.position.z = p.rz;
                    }
                    if(p.rot) {
                        rotationGroup.quaternion = p.rot;
                    }
                    if(p.cam) {
                        camera.lookAt(p.cam);
                    }
                    
                    if(step < steps.length) {
                        setTimeout(callback, interval);
                    } else {
                        decAnim();
                    }
                    show();
                };
                setTimeout(callback, interval);
               
            } else { //relative update
                var delta = {};
                let frac = 1.0/steps;
                if(mpos) {
                    delta.mpos = mpos.clone().sub(curr.mpos).multiplyScalar(frac);
                }
                if(typeof(rz) != 'undefined' && rz != null) {
                    delta.rz = frac*(rz-curr.rz);
                }
                if(rot) {
                    var next = slerp(curr.rot,rot,frac);
                    //comptute step delta rotation
                    delta.rot = curr.rot.clone().inverse().multiply(next);
                }
                if(cam) {
                    delta.cam = cam.clone().sub(curr.cam).multiplyScalar(frac);
                }
                let step = 0.0;
                let callback = function() {
                    step += 1;
                    if(delta.mpos) {
                        modelGroup.position.add(delta.mpos);
                    }
                    if(delta.rz) {
                        rotationGroup.position.z += delta.rz;
                    }
                    if(delta.rot) {
                        rotationGroup.quaternion.multiply(delta.rot);
                    }
                    if(delta.cam) {
                        lookingAt.add(delta.cam);
                        camera.lookAt(lookingAt);
                    }
                    
                    if(step < steps) {
                        setTimeout(callback, interval);
                    } else {
                        decAnim();
                    }
                    show();
                };
                setTimeout(callback, interval);
            }
        };
        
        /**
         * Rotate scene by angle degrees around axis
         * 
         * @function $3Dmol.GLViewer#rotate
         * @param {number}
         *            [angle] - Angle, in degrees, to rotate by.
         * @param {string}
         *            [axis] - Axis ("x", "y", "z", "vx", "vy", or "vz") to rotate around.
         *            Default "y".  View relative (rather than model relative) axes are prefixed with v.
         *            Axis can also be specified as a vector.
         * @param {number}
         *            [animationDuration] - an optional parameter that denotes
         *            the duration of the rotation animation. Default 0 (no animation)
         * @param {boolean} [fixedPath] - if true animation is constrained to 
         *      requested motion, overriding updates that happen during the animation         *            
         * @example     $3Dmol.download('cid:4000', viewer, {}, function() {
      viewer.setStyle({stick:{}});
      viewer.zoomTo();
      viewer.rotate(90,'y',1);
      viewer.render(callback);
    });

         *  
         */
        this.rotate = function(angle, axis, animationDuration, fixedPath) {
            animationDuration = animationDuration!==undefined ? animationDuration : 0;

            if (typeof (axis) === "undefined") {
                axis = "y";
            }

            if(axis == "x"){
                axis = {x:1,y:0,z:0};
            }else if(axis =="y"){
                axis = {x:0,y:1,z:0};
            }else if(axis =="z"){
                axis = {x:0,y:0,z:1};
            }
            
            //support rotating with respect to view axis, not model
            if(axis == "vx"){
                axis = {vx:1,vy:0,vz:0};
            }else if(axis =="vy"){
                axis = {vx:0,vy:1,vz:0};
            }else if(axis =="vz"){
                axis = {vx:0,vy:0,vz:1};
            }
            
            if(typeof(axis.vx) !== 'undefined') {
              var vaxis = new $3Dmol.Vector3(axis.vx,axis.vy,axis.vz);
              vaxis.applyQuaternion(rotationGroup.quaternion);
              axis = {x:vaxis.x, y:vaxis.y, z: vaxis.z};
            }
                        
            var qFromAngle = function(rangle) {
                var s = Math.sin(rangle / 2.0);
                var c = Math.cos(rangle / 2.0);
                var i = 0, j = 0, k = 0;

                i = axis.x * s;
                j = axis.y * s;
                k = axis.z * s;

                return new $3Dmol.Quaternion(i, j, k, c).normalize();
            };
            
            var rangle = Math.PI * angle / 180.0;
            var q = qFromAngle(rangle);
            
            if(animationDuration ){
                var final = new $3Dmol.Quaternion().copy(rotationGroup.quaternion).multiply(q);//final
                animateMotion(animationDuration,fixedPath,
                        modelGroup.position,
                        rotationGroup.position.z, 
                        final,
                        lookingAt);              
            } else { //not animated
                rotationGroup.quaternion.multiply(q);
                show();
            }
            return this;

        };

        this.surfacesFinished= function() {
              for(var key in surfaces){
                if(!surfaces[key][0].done){
                    return false;
                }
            }
            return true;


        };

        /** Returns an array representing the current viewpoint.
         * Translation, zoom, and rotation quaternion. 
         * @function $3Dmol.GLViewer#getView
         * @returns {Array.<number>} [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y, q.z, q.w ]
         *  */
        this.getView = function() {
            if (!modelGroup)
                return [ 0, 0, 0, 0, 0, 0, 0, 1 ];
            var pos = modelGroup.position;
            var q = rotationGroup.quaternion;
            return [ pos.x, pos.y, pos.z, rotationGroup.position.z, q.x, q.y,
                    q.z, q.w ];
        };

        /** Sets the view to the specified translation, zoom, and rotation.
         * 
         * @function $3Dmol.GLViewer#setView
         * @param {Array.<number>} arg Array formatted identically to the return value of getView */
        this.setView = function(arg, nolink) {

            if (arg === undefined ||
                    !(arg instanceof Array || arg.length !== 8))
                return this;

            if (!modelGroup || !rotationGroup)
                return this;
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

            show(nolink);
            return this;

        };

        // apply styles, models, etc in viewer
        /**
         * Render current state of viewer, after 
         * adding/removing models, applying styles, etc.
         * 
         * @function $3Dmol.GLViewer#render
         */
        this.render = function(callback, exts) {
            renderer.setViewport();
            updateClickables(); //must render for clickable styles to take effect
            var view = this.getView();
            
            if(stateChangeCallback) {
              //todo: have ability to only send delta updates
              stateChangeCallback(this.getInternalState());
            }
            
            var i, n;
            if(!exts) exts = renderer.supportedExtensions();
            for (i = 0; i < models.length; i++) {
                if (models[i]) {
                    models[i].globj(modelGroup, exts);
                }
            }

            for (i = 0; i < shapes.length; i++) {
                if (shapes[i]) { //exists
                    if ((typeof(shapes[i].frame) === 'undefined' || viewer_frame < 0 ||
                                    shapes[i].frame < 0 || shapes[i].frame == viewer_frame)) {
                        shapes[i].globj(modelGroup, exts);
                    } else { //should not be displayed in current frame
                        shapes[i].removegl(modelGroup);
                    }
                }
            }
            
            for (i = 0; i < labels.length; i++) {
                if (labels[i] && typeof(labels[i].frame) != 'undefined' && labels[i].frame >= 0) { //exists and has frame specifier
                    modelGroup.remove(labels[i].sprite);
                    if (viewer_frame < 0 || labels[i].frame == viewer_frame) {
                        modelGroup.add(labels[i].sprite);
                    }
                }
            }            
            
            for (i in surfaces) { // this is an object with possible holes
                if(!surfaces.hasOwnProperty(i)) continue;
                var surfArr = surfaces[i];
                for (n = 0; n < surfArr.length; n++) {
                    if (surfArr.hasOwnProperty(n)) {
                        var geo = surfArr[n].geo;
                        // async surface generation can cause
                        // the geometry to be webgl initialized before it is fully
                        // formed; force various recalculations until full surface
                        // is
                        // available
                        if (!surfArr[n].finished || exts.regen) {
                            geo.verticesNeedUpdate = true;
                            geo.elementsNeedUpdate = true;
                            geo.normalsNeedUpdate = true;
                            geo.colorsNeedUpdate = true;
                            geo.buffersNeedUpdate = true;
                            geo.boundingSphere = null;

                            if (surfArr[n].done)
                                surfArr[n].finished = true;

                            // remove partially rendered surface
                            if (surfArr[n].lastGL)
                                modelGroup.remove(surfArr[n].lastGL);

                            // create new surface
                            var smesh = null;

                            if(surfArr[n].mat instanceof $3Dmol.LineBasicMaterial) {
                                //special case line meshes
                                smesh = new $3Dmol.Line(geo, surfArr[n].mat);
                            }
                            else {
                                smesh = new $3Dmol.Mesh(geo, surfArr[n].mat);
                            }
                            if(surfArr[n].mat.transparent && surfArr[n].mat.opacity == 0) {
                                //don't bother with hidden surfaces
                                smesh.visible = false;
                            } else {
                                smesh.visible = true;
                            }
                            if (surfArr[n].symmetries.length > 1 || 
                            (surfArr[n].symmetries.length == 1 && 
                            !(surfArr[n].symmetries[n].isIdentity()))) {
                                var j;
                                var tmeshes = new $3Dmol.Object3D(); //transformed meshes
                                for (j = 0; j < surfArr[n].symmetries.length; j++) {
                                    var tmesh = smesh.clone();
                                    tmesh.matrix = surfArr[n].symmetries[j];
                                    tmesh.matrixAutoUpdate = false;
                                    tmeshes.add(tmesh);
                                }
                                surfArr[n].lastGL = tmeshes;
                                modelGroup.add(tmeshes);
                            }
                            else {
                                surfArr[n].lastGL = smesh;
                                modelGroup.add(smesh);
                            }
                        } // else final surface already there
                    }
                }
            }
            
            this.setView(view); // Calls show() => renderer render
            if(typeof callback ==='function'){
                callback(this);
            }
            return this;
        };

        /** @param {ATomSelectionSpec} sel
         * @return list of models specified by sel
         */
        function getModelList(sel) {
            var ms = [];
            if (typeof sel === 'undefined' || typeof sel.model === "undefined") {
                for (let i = 0; i < models.length; i++) {
                    if (models[i])
                        ms.push(models[i]);
                }
            } else { // specific to some models
                ms = sel.model;
                if (!Array.isArray(ms))
                    ms = [ ms ];
                
                for (let i = 0; i < ms.length; i++) {
                        //allow referencing models by order of creation
                    if(typeof ms[i] === 'number') {
                        var index = ms[i];
                        //support python backward indexing
                        if(index < 0) index += models.length;
                        ms[i] = models[index];
                    }             
                }
            }
            
            return ms;            
        }
        /**
         * 
         * @param {AtomSelectionSpec}
         *            sel
         * @return {AtomSpec[]}
         */
        function getAtomsFromSel(sel) {
            var atoms = [];
            if (typeof (sel) === "undefined")
                sel = {};

            var ms = getModelList(sel);

            for (var i = 0; i < ms.length; i++) {
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

            var ms = getModelList(sel);
            
            for (var i = 0; i < ms.length; i++) {
                if (ms[i].atomIsSelected(atom, sel))
                    return true;
            }

            return false;
        }


        /** return list of atoms selected by sel
         * 
         * @function $3Dmol.GLViewer#selectedAtoms
         * @param {AtomSelectionSpec} sel
         * @return {Array.<Object>}
         */
        this.selectedAtoms = function(sel) {
            return getAtomsFromSel(sel);
        };

        /**
        * Returns valid values for the specified attribute in the given selection
        * @function $3Dmol.GlViewer#getUniqueValues 
        * @param {string} attribute
        * @param {AtomSelectionSpec} sel
        * @return {Array.<Object>}
        *
        */
        this.getUniqueValues = function(attribute, sel){
            if (typeof (sel) === "undefined")
                sel = {};            
            var atoms = getAtomsFromSel(sel);
            var values = {};

            for(var atom in atoms){
                if(atoms[atom].hasOwnProperty(attribute)){
                    var value = atoms[atom][attribute];
                    values[value] = true;
                }
            }

            return Object.keys(values);
        };
        
        /**
         * Return pdb output of selected atoms (if atoms from pdb input)
         * 
         * @function $3Dmol.GLViewer#pdbData  
         * @param {Object=} [sel] - Selection specification specifying model and atom properties to select.  Default: all atoms in viewer
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
         * @function $3Dmol.GLViewer#zoom
         * @param {number}
         *            [factor] - Magnification factor. Values greater than 1
         *            will zoom in, less than one will zoom out. Default 2.
         * @param {number}
         *            [animationDuration] - an optional parameter that denotes
         *            the duration of a zoom animation
         * @param {Boolean} [fixedPath] - if true animation is constrained to 
         *      requested motion, overriding updates that happen during the animation
         * @example   
    $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo()
      viewer.zoom(2,1000);
      viewer.render();
    });
    
             */
        this.zoom = function(factor,animationDuration,fixedPath) {
            factor = factor || 2;
            animationDuration = animationDuration!==undefined ? animationDuration : 0;
            var scale = (CAMERA_Z - rotationGroup.position.z) / factor;
            var final_z = CAMERA_Z - scale;

            if(animationDuration>0){
                animateMotion(animationDuration,fixedPath,
                        modelGroup.position, 
                        adjustZoomToLimits(final_z), 
                        rotationGroup.quaternion,
                        lookingAt);
            } else { //no animation
                rotationGroup.position.z = adjustZoomToLimits(final_z);
                show();
            }
            return this;
        };
        
        /**
         * Translate current view by x,y screen coordinates
         * This pans the camera rather than translating the model.
         * 
         * @function $3Dmol.GLViewer#translate
         * @param {number} x Relative change in view coordinates of camera
         * @param {number} y Relative change in view coordinates of camera
         * @param {number}
         *            [animationDuration] - an optional parameter that denotes
         *            the duration of a zoom animation
         * @param {Boolean} [fixedPath] - if true animation is constrained to 
         *      requested motion, overriding updates that happen during the animation         *            
         * @example     $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      viewer.translate(200,50);         
      viewer.rotate(90,'z');
      viewer.render(callback);
    });
         */
        this.translate = function(x, y, animationDuration, fixedPath) {
            animationDuration = animationDuration!==undefined ? animationDuration : 0;
            var dx = x/WIDTH;
            var dy = y/HEIGHT;
            var v = new $3Dmol.Vector3(0,0,-CAMERA_Z);
            
            projector.projectVector(v, camera);
            v.x -= dx;
            v.y -= dy;
            projector.unprojectVector(v, camera);
            v.z = 0;            

            var final_position=lookingAt.clone().add(v);
            if(animationDuration>0){
                animateMotion(animationDuration,fixedPath,
                        modelGroup.position,
                        rotationGroup.position.z, 
                        rotationGroup.quaternion,
                        final_position);
            } else { //no animation
                lookingAt = final_position;
                camera.lookAt(lookingAt);
                show();
            }
            return this;
        };
        
        /**
         * Translate current models by x,y screen coordinates
         * This translates the models relative to the current view. It does
         * not change the center of rotation.
         * 
         * @function $3Dmol.GLViewer#translateScene
         * @param {number} x Relative change in x screen coordinate
         * @param {number} y Relative change in y screen coordinate
         * @param {number}
         *            [animationDuration] - an optional parameter that denotes
         *            the duration of a zoom animation
         * @param {Boolean} [fixedPath] - if true animation is constrained to 
         *      requested motion, overriding updates that happen during the animation         *            
         * @example     $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      viewer.translateScene(200,50);         
      viewer.rotate(90,'z'); // will no longer be around model center
      viewer.render(callback);
    });
         */
        this.translateScene = function(x, y, animationDuration, fixedPath) {
            animationDuration = animationDuration!==undefined ? animationDuration : 0;
            
            var t = screenOffsetToModel(x,y);
            var final_position=modelGroup.position.clone().add(t);
                
            if(animationDuration>0){
                animateMotion(animationDuration,fixedPath,
                        modelGroup.position,
                        rotationGroup.position.z, 
                        rotationGroup.quaternion,
                        lookingAt);
            } else { //no animation
                modelGroup.position = final_position;
                show();
            }
            return this;
        };

        /**
         * Adjust slab to fully enclose selection (default everything).
         * 
         * @function $3Dmol.GLViewer#fitSlab
         * @param {Object}
         *            [sel] - Selection specification specifying model and atom
         *            properties to select. Default: all atoms in viewer
         */
        this.fitSlab = function(sel) {
            sel = sel || {};
            var atoms = getAtomsFromSel(sel);
            var tmp = $3Dmol.getExtent(atoms);

            // fit to bounding box
            var x = tmp[1][0] - tmp[0][0], 
                y = tmp[1][1] - tmp[0][1], 
                z = tmp[1][2] - tmp[0][2];

            var maxD = Math.sqrt(x * x + y * y + z * z);
            if (maxD < 5)
                maxD = 5;

            // use full bounding box for slab/fog
            slabNear = -maxD / 1.9;
            slabFar = maxD / 2;

            return this;
        };        
        
        /**
         * Re-center the viewer around the provided selection (unlike zoomTo, does not zoom).
         * 
         * @function $3Dmol.GLViewer#center
         * @param {Object}
         *            [sel] - Selection specification specifying model and atom
         *            properties to select. Default: all atoms in viewer
         * @param {number}
         *            [animationDuration] - an optional parameter that denotes
         *            the duration of a zoom animation
         * @param {Boolean} [fixedPath] - if true animation is constrained to 
         *      requested motion, overriding updates that happen during the animation         *            
         * @example // if the user were to pass the animationDuration value to 
         *           // the function like so viewer.zoomTo({resn:'STI'},1000);
         *         //   the program would center on resn 'STI' over the course 
         *         //   of 1 second(1000 milleseconds).
         *  // Reposition to centroid of all atoms of all models in this
         * //viewer glviewer.center();
    $.get('data/4csv.pdb', function(data) {
      viewer.addModel(data,'pdb');
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.center();
      viewer.render(callback);
    });
         */
        this.center = function(sel,animationDuration,fixedPath){
             animationDuration=animationDuration!==undefined ? animationDuration : 0;
            var allatoms, alltmp;
            sel = sel || {};
            var atoms = getAtomsFromSel(sel);
            var tmp = $3Dmol.getExtent(atoms);

            if($3Dmol.isEmptyObject(sel)) {
                //include shapes when zooming to full scene
                //TODO: figure out a good way to specify shapes as part of a selection
                shapes.forEach((shape) => {
                    if(shape && shape.boundingSphere && shape.boundingSphere.center) {
                        var c = shape.boundingSphere.center;
                        var r = shape.boundingSphere.radius;
                        if(r > 0) {
                            //make sure full shape is visible
                            atoms.push(new $3Dmol.Vector3(c.x+r,c.y,c.z));
                            atoms.push(new $3Dmol.Vector3(c.x-r,c.y,c.z));
                            atoms.push(new $3Dmol.Vector3(c.x,c.y+r,c.z));
                            atoms.push(new $3Dmol.Vector3(c.x,c.y-r,c.z));
                            atoms.push(new $3Dmol.Vector3(c.x,c.y,c.z+r));
                            atoms.push(new $3Dmol.Vector3(c.x,c.y,c.z-r));
                        } else {
                            atoms.push(c);
                        }
                    }
                });
                tmp = $3Dmol.getExtent(atoms);
                allatoms = atoms;
                alltmp = tmp;

            }
            else {
                allatoms = getAtomsFromSel({});
                alltmp = $3Dmol.getExtent(allatoms);
            }

            // use selection for center
            var center = new $3Dmol.Vector3(tmp[2][0], tmp[2][1], tmp[2][2]);

            // but all for bounding box
            var x = alltmp[1][0] - alltmp[0][0], y = alltmp[1][1] -
                     alltmp[0][1], z = alltmp[1][2] - alltmp[0][2];

            var maxD = Math.sqrt(x * x + y * y + z * z);
            if (maxD < 5)
                maxD = 5;

            // use full bounding box for slab/fog
            slabNear = -maxD / 1.9;
            slabFar = maxD / 2;

            // for zoom, use selection box
            x = tmp[1][0] - tmp[0][0];
            y = tmp[1][1] - tmp[0][1];
            z = tmp[1][2] - tmp[0][2];
            maxD = Math.sqrt(x * x + y * y + z * z);
            if (maxD < 5)
                maxD = 5;
            
            //find the farthest atom from center to get max distance needed for view
            var maxDsq = 25;
            for (var i = 0; i < atoms.length; i++) {
                if(atoms[i]) {
                    var dsq = center.distanceToSquared(atoms[i]);
                    if(dsq > maxDsq)
                        maxDsq = dsq;
                }
            }
            
            maxD = Math.sqrt(maxDsq)*2;
            var finalpos = center.clone().multiplyScalar(-1);
            if(animationDuration>0){
                animateMotion(animationDuration,fixedPath,
                        finalpos, 
                        rotationGroup.position.z, 
                        rotationGroup.quaternion,
                        lookingAt);
            } else { //no animation 
                modelGroup.position = finalpos;
                show();
            }
            return this;
        };
        
        /**
         * Zoom to center of atom selection.  The slab will be set appropriately for
         * the selection, unless an empty selection is provided, in which case there will be no slab.
         * 
         * @function $3Dmol.GLViewer#zoomTo
         * @param {Object}
         *            [sel] - Selection specification specifying model and atom
         *            properties to select. Default: all atoms in viewer
         * @param {number}
         *            [animationDuration] - an optional parameter that denotes
         *            the duration of a zoom animation
         * @param {Boolean} [fixedPath] - if true animation is constrained to 
         *      requested motion, overriding updates that happen during the animation         *            
          * @example   
    

              $.get('data/1fas.pqr', function(data){
                  viewer.addModel(data, "pqr");
                  $.get("data/1fas.cube",function(volumedata){
                      viewer.addSurface($3Dmol.SurfaceType.VDW, {
                          opacity:0.85,
                          voldata: new $3Dmol.VolumeData(volumedata, "cube"),
                          volscheme: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'charge'))
                      },{});
                      
                  viewer.render();
                  });
                  viewer.zoomTo();
                });
         */
        this.zoomTo = function(sel, animationDuration,fixedPath) {
            animationDuration=animationDuration!==undefined ? animationDuration : 0;
            sel = sel || {};
            let atoms = getAtomsFromSel(sel);
            let atombox = $3Dmol.getExtent(atoms);
            let allbox = atombox;

            if($3Dmol.isEmptyObject(sel)) {
                //include shapes when zooming to full scene
                //TODO: figure out a good way to specify shapes as part of a selection
                let natoms = atoms && atoms.length;
                shapes.forEach((shape) => {
                if(shape && shape.boundingSphere) {
                    if(shape.boundingSphere.box) {
                        let box = shape.boundingSphere.box;
                        atoms.push(new $3Dmol.Vector3(box.min.x,box.min.y,box.min.z));
                        atoms.push(new $3Dmol.Vector3(box.max.x,box.max.y,box.max.z));                        
                    } else if(shape.boundingSphere.center) {
                        var c = shape.boundingSphere.center;
                        var r = shape.boundingSphere.radius;
                        if(r > 0) {
                            //make sure full shape is visible
                                atoms.push(new $3Dmol.Vector3(c.x+r,c.y,c.z));
                                atoms.push(new $3Dmol.Vector3(c.x-r,c.y,c.z));
                                atoms.push(new $3Dmol.Vector3(c.x,c.y+r,c.z));
                                atoms.push(new $3Dmol.Vector3(c.x,c.y-r,c.z));
                                atoms.push(new $3Dmol.Vector3(c.x,c.y,c.z+r));
                                atoms.push(new $3Dmol.Vector3(c.x,c.y,c.z-r));
                        } else {
                                atoms.push(c);
                        }
                    }
                  }
                });
                allbox = $3Dmol.getExtent(atoms);
                if(!natoms) { //if no atoms, use shapes for center
                    for(let i = 0; i < 3; i++) { //center of bounding box
                        atombox[2][i] = (allbox[0][i]+allbox[1][i])/2; 
                     }
                }
            } else { //include all atoms in slab calculation
                let allatoms = getAtomsFromSel({});
                allbox = $3Dmol.getExtent(allatoms);  
            }

            // use selection for center
            var center = new $3Dmol.Vector3(atombox[2][0], atombox[2][1], atombox[2][2]);
            
            // but all for bounding box
            var x = allbox[1][0] - allbox[0][0], y = allbox[1][1]
                    - allbox[0][1], z = allbox[1][2] - allbox[0][2];

            var maxD = Math.sqrt(x * x + y * y + z * z);
            if (maxD < 5)
                maxD = 5;

            // use full bounding box for slab/fog
            slabNear = -maxD / 1.9;
            slabFar = maxD / 2;
            
            //if we are selecting everything, do not slab
            if(Object.keys(sel).length === 0) {
                slabNear = -999999;
                slabFar = 999999;
            }

            // keep at least this much space in view
            var MAXD = config.minimumZoomToDistance || 5;
            // for zoom, use selection box
            x = atombox[1][0] - atombox[0][0];
            y = atombox[1][1] - atombox[0][1];
            z = atombox[1][2] - atombox[0][2];
            maxD = Math.sqrt(x * x + y * y + z * z);           
            if (maxD < MAXD)
                maxD = MAXD;
            
            //find the farthest atom from center to get max distance needed for view
            var maxDsq = MAXD*MAXD;
            for (var i = 0; i < atoms.length; i++) {
                if(atoms[i]) {
                    var dsq = center.distanceToSquared(atoms[i]);
                    if(dsq > maxDsq)
                        maxDsq = dsq;
                }
            }
            
            maxD = Math.sqrt(maxDsq)*2;
            var finalpos = center.clone().multiplyScalar(-1);
            var finalz =  -(maxD * 0.5
                    / Math.tan(Math.PI / 180.0 * camera.fov / 2) - CAMERA_Z);
                    
            finalz = adjustZoomToLimits(finalz);
            if(animationDuration>0){
                animateMotion(animationDuration,fixedPath,
                        finalpos,
                        finalz, 
                        rotationGroup.quaternion,
                        lookingAt);                
            } else {
                modelGroup.position = finalpos;
                rotationGroup.position.z = finalz;
                show();
            }
            return this;
        
        };

        /**
         * Set slab of view (contents outside of slab are clipped).
         * Must call render to update.
         * 
         * @function $3Dmol.GLViewer#setSlab
         * @param {number} near near clipping plane distance
         * @param {number} far far clipping plane distance
         */
        this.setSlab = function(near, far) {
            slabNear = near;
            slabFar = far;
        };
        
        /**
         * Get slab of view (contents outside of slab are clipped).
         * 
         * @function $3Dmol.GLViewer#getSlab
         * @return {Object}
         *      @property {number} near - near clipping plane distance
         *      @property {number} far - far clipping plane distance
         */
        this.getSlab = function() {
            return {near: slabNear, far: slabFar};
        };
                
        /**
         * Add label to viewer
         * 
         * @function $3Dmol.GLViewer#addLabel
         * @param {string}
         *            text - Label text
         * @param {LabelSpec}
         *            options - Label style specification
          @param {AtomSelection}
         *            sel - Set position of label to center of this selection
         * @param {boolean} noshow - if true, do not immediately display label - when adding multiple labels this is more efficient
         * @return {$3Dmol.Label}
         * 
         * @example
         *  $3Dmol.download("pdb:2EJ0",viewer,{},function(){
                  
                  viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
                  viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                           borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                           position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
                  viewer.setStyle({chain:'A'},{cross:{hidden:true}});
                  viewer.setStyle({chain:'B'},{cross:{hidden:false,
                                                      linewidth:1.0,
                                                      colorscheme:'greenCarbon'}});
                  viewer.setStyle({chain:'C'},{cross:{hidden:false,
                                                      linewidth:1.0,
                                                      radius:0.5}});
                  viewer.setStyle({chain:'D'},{cross:{hidden:false,
                                                      linewidth:10.0}});
                  viewer.setStyle({chain:'E'},{cross:{hidden:false,
                                                      linewidth:1.0,
                                                      color:'black'}});
                  
                  viewer.render();

                  
                });
            
         */
        this.addLabel = function(text, options, sel, noshow) {
            options = options || {};
            if(sel) {
                var extent = $3Dmol.getExtent(getAtomsFromSel(sel));
                options.position = {x: extent[2][0], y: extent[2][1], z: extent[2][2]};
            }
            var label = new $3Dmol.Label(text, options);
            label.setContext();
            modelGroup.add(label.sprite);
            if(options.fixed)
                fixed_labels.push(labels.length);
            labels.push(label);

            if(!noshow) show();
            return label;
        };
        


        /** Add residue labels.  This will generate one label per a
         * residue within the selected atoms.  The label will be at the
         * centroid of the atoms and styled according to the passed style.
         * The label text will be [resn][resi]
         * 
         * @function $3Dmol.GLViewer#addResLabels
         * @param {Object} sel
         * @param {Object} style
         * @param {boolean} byframe - if true, create labels for every individual frame, not just current
         * 
         * * @example  
             $3Dmol.download("mmtf:2ll5",viewer,{},function(){
                  viewer.setStyle({stick:{radius:0.15},cartoon:{}});
                  viewer.addResLabels({hetflag:false}, {font: 'Arial', fontColor:'black',showBackground:false, screenOffset: {x:0,y:0}});
                  viewer.zoomTo();
                  viewer.render();                  
                });
         */
        this.addResLabels = function(sel, style, byframe) {
            let start = labels.length;
            applyToModels("addResLabels", sel, this, style, byframe);
            show();
            return labels.slice(start);
        };

        /** Add property labels.  This will generate one label per a selected
         * atom at the atom's coordinates with the property value as the label text.
         * 
         * @function $3Dmol.GLViewer#addPropertyLabels
         * @param {string} prop - property name
         * @param {Object} sel
         * @param {Object} style
         * 
         * * @example  
             $3Dmol.download("cid:5291",viewer,{},function(){
                  viewer.setStyle({stick: {radius:.2}});
                  viewer.addPropertyLabels("index",{not:{elem:'H'}}, {fontColor:'black',font: 'sans-serif', fontSize: 28, showBackground:false,alignment:'center'});
                  viewer.zoomTo();
                  viewer.render();                  
                });
         */
        this.addPropertyLabels = function(prop, sel, style) {
            applyToModels("addPropertyLabels", prop, sel, this, style);
            show();
            return this;
        };
        
        /**
         * Remove label from viewer
         * 
         * @function $3Dmol.GLViewer#removeLabel
         * @param {$3Dmol.Label}
         *            label - $3Dmol label
         * 
         * @example // Remove labels created in 
         $3Dmol.download("pdb:2EJ0",viewer,{},function(){
                  var toremove = viewer.addLabel("Aromatic", {position: {x:-6.89, y:0.75, z:0.35}, backgroundColor: 0x800080, backgroundOpacity: 0.8});
                  viewer.addLabel("Label",{font:'sans-serif',fontSize:18,fontColor:'white',fontOpacity:1,borderThickness:1.0,
                                           borderColor:'red',borderOpacity:0.5,backgroundColor:'black',backgroundOpacity:0.5,
                                           position:{x:50.0,y:0.0,z:0.0},inFront:true,showBackground:true});
                  viewer.removeLabel(toremove);
                  viewer.render();

                  
                });

         */
        this.removeLabel = function(label) {
            //todo: don't do the linear search
            for(var i = 0; i < labels.length; i++) {
                if(labels[i] == label) {
                    labels.splice(i,1);
                    label.dispose();
                    modelGroup.remove(label.sprite);
                    break;
                }
            }
            show();
            return this;
        };



        /**
         * Remove all labels from viewer
         * 
         * @function $3Dmol.GLViewer#removeAllLabels
         *         @example             
        $3Dmol.download("pdb:1ubq",viewer,{},function(){
                          
               viewer.addResLabels();
               viewer.setStyle({},{stick:{}}); 
               viewer.render( ); //show labels

               viewer.removeAllLabels();              
               viewer.render(); //hide labels
        });
         */
        this.removeAllLabels = function() {
            for (var i = 0; i < labels.length; i++) {
              if(labels[i] && labels[i].sprite) {
                modelGroup.remove(labels[i].sprite);
              }
            }
            labels.splice(0,labels.length); //don't overwrite in case linked
            show();
            return this;
        };
        
        // Modify label style
        /**
         * Modify existing label's style
         * 
         * @function $3Dmol.GLViewer#setLabelStyle
         * @param {$3Dmol.Label}
         *            label - $3Dmol label
         * @param {Object}
         *            stylespec - Label style specification
         * @return {$3Dmol.Label}
         */
        this.setLabelStyle = function(label, stylespec) {
            modelGroup.remove(label.sprite);
            label.dispose();
            label.stylespec = stylespec;
            label.setContext();
            modelGroup.add(label.sprite);
            show();
            return label;

        };

        // Change label text
        /**
         * Modify existing label's text
         * 
         * @function $3Dmol.GLViewer#setLabelText
         * @param {$3Dmol.Label}
         *            label - $3Dmol label
         * @param {String}
         *            text - Label text
         * @return {$3Dmol.Label}
         */
        this.setLabelText = function(label, text) {
            modelGroup.remove(label.sprite);
            label.dispose();
            label.text = text;
            label.setContext();
            modelGroup.add(label.sprite);
            show();
            return label;

        };

        /**
         * Add shape object to viewer 
         * @see {@link $3Dmol.GLShape}
         * 
         * @function $3Dmol.GLViewer#addShape
         * @param {ShapeSpec} shapeSpec - style specification for label
         * @return {$3Dmol.GLShape}
         */
        this.addShape = function(shapeSpec) {
            shapeSpec = shapeSpec || {};
            var shape = new $3Dmol.GLShape(shapeSpec);
            shape.shapePosition = shapes.length;
            shapes.push(shape);

            return shape;

        };

        /**
         * Remove shape object from viewer
         *
         * @function $3Dmol.GLViewer#removeShape
         * @param {$3Dmol.GLShape} shape - Reference to shape object to remove
         */
        this.removeShape = function(shape) {
            if (!shape)
                return this;
            shape.removegl(modelGroup);
            delete shapes[shape.shapePosition];
            // clear off back of model array
            while (shapes.length > 0
                    && typeof (shapes[shapes.length - 1]) === "undefined")
                shapes.pop();
            return this;
        };
        
        /**
         * Remove all shape objects from viewer
         * @function $3Dmol.GLViewer#removeAllShapes
         */
        this.removeAllShapes = function() {
            for (var i = 0; i < shapes.length; i++) {
                var shape = shapes[i];
                if(shape) shape.removegl(modelGroup);
            }
            shapes.splice(0,shapes.length);
            return this;
        };

        //gets the center of the selection 
        var getSelectionCenter = function(spec){
            if(spec.hasOwnProperty("x") && spec.hasOwnProperty("y") && spec.hasOwnProperty("z"))
                return spec;
            var atoms = getAtomsFromSel(spec);
            if(atoms.length == 0)
                return {x:0,y:0,z:0};

            var extent = $3Dmol.getExtent(atoms);
            return {x:extent[0][0]+(extent[1][0]-extent[0][0])/2,y:extent[0][1]+(extent[1][1]-extent[0][1])/2,z:extent[0][2]+(extent[1][2]-extent[0][2])/2};
        };

        /**
         * Create and add sphere shape. This method provides a shorthand 
         * way to create a spherical shape object
         * 
         * @function $3Dmol.GLViewer#addSphere
         * @param {SphereShapeSpec} spec - Sphere shape style specification
         * @return {$3Dmol.GLShape}
         @example
         
         viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
         
         viewer.render();
         */
        this.addSphere = function(spec) {
            spec = spec || {};

            spec.center = getSelectionCenter(spec.center);

            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            s.addSphere(spec);
            shapes.push(s);
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended
            return s;
        };

        /**
         * Create and add box shape. This method provides a shorthand 
         * way to create a box shape object
         * 
         * @function $3Dmol.GLViewer#addBox
         * @param {BoxSpec} spec - Box shape style specification
         * @return {$3Dmol.GLShape}
         @example
         
         viewer.addLine({color:'red',start:{x:0,y:0,z:0},end:{x:5,y:0,z:0}});
         viewer.addLine({color:'blue',start:{x:0,y:0,z:0},end:{x:0,y:5,z:0}});
         viewer.addLine({color:'green',start:{x:0,y:0,z:0},end:{x:0,y:0,z:5}});

         viewer.addBox({center:{x:0,y:0,z:0},dimensions: {w:3,h:4,d:2},color:'magenta'});
         viewer.zoomTo();
         viewer.rotate(45, {x:1,y:1,z:1});
         viewer.render();
         */
        this.addBox = function(spec) {
            spec = spec || {};

            if(spec.corner != undefined) {
                spec.corner = getSelectionCenter(spec.corner);
            }
            if(spec.center != undefined) {
                spec.center = getSelectionCenter(spec.center);
            }

            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            s.addBox(spec);
            shapes.push(s);
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

            return s;
        };
        
        /**
         * Create and add arrow shape
         * 
         * @function $3Dmol.GLViewer#addArrow
         * @param {ArrowSpec} spec - Style specification
         * @return {$3Dmol.GLShape}
         @example
        $3Dmol.download("pdb:4DM7",viewer,{},function(){

                  viewer.setBackgroundColor(0xffffffff);
                  viewer.addArrow({
                      start: {x:-10.0, y:0.0, z:0.0},
                      end: {x:0.0, y:-10.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      clickable:true,
                      callback:function(){
                          this.color.setHex(0xFF0000FF);
                          viewer.render( );
                      }
                  });
                  viewer.render();
                });
         */
        this.addArrow = function(spec) {
            spec = spec || {};
            
            spec.start = getSelectionCenter(spec.start);
            spec.end = getSelectionCenter(spec.end);      
            
            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            s.addArrow(spec);
            shapes.push(s);
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

            return s;
        };
        
        /**
         * Create and add cylinder shape
         * 
         * @function $3Dmol.GLViewer#addCylinder
         * @param {CylinderSpec} spec - Style specification
         * @return {$3Dmol.GLShape}

          @example
         viewer.setBackgroundColor(0xffffffff);
              viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
                                  end:{x:10.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  fromCap:1,
                                  toCap:2,
                                  color:'red',
                                  hoverable:true,
                                  clickable:true,
                                  callback:function(){ this.color.setHex(0x00FFFF00);viewer.render( );},
                                  hover_callback: function(){ viewer.render( );},
                                  unhover_callback: function(){ this.color.setHex(0xFF000000);viewer.render( );}
                                 });
              viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
                                  end:{x:0.0,y:10.0,z:0.0},
                                  radius:0.5,
                                  fromCap:false,
                                  toCap:true,
                                  color:'teal'});
              viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
                                  end:{x:20.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  color:'black',
                                  fromCap:false,
                                  toCap:false});
              viewer.render();
         */
        this.addCylinder = function(spec) {
            spec = spec || {};

            spec.start = getSelectionCenter(spec.start);
            spec.end = getSelectionCenter(spec.end);

            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            if(spec.dashed)
                s.addDashedCylinder(spec);
            else
                s.addCylinder(spec);
            shapes.push(s);
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

            return s;
        };

        /**
         * Create and add Curve shape
         * 
         * @function $3Dmol.GLViewer#addCurve
         * @param {CurveSpec} spec - Style specification
         * @return {$3Dmol.GLShape}

         @example
              viewer.addCurve({points: [{x:0.0,y:0.0,z:0.0}, {x:5.0,y:3.0,z:0.0}, {x:5.0,y:7.0,z:0.0}, {x:0.0,y:10.0,z:0.0}],
                                  radius:0.5,
                                  smooth: 10,
                                  fromArrow:false,
                                  toArrow: true,
                                  color:'orange',                                  
                                  });
              viewer.addCurve({points: [{x:-1,y:0.0,z:0.0}, {x:-5.0,y:5.0,z:0.0}, {x:-2,y:10.0,z:0.0}],
                                  radius:1,
                                  fromArrow:true,
                                  toArrow: false,
                                  color:'purple',                                  
                                  });
              viewer.zoomTo();
              viewer.render();
         */
        this.addCurve = function(spec) {
            spec = spec || {};            
            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            s.addCurve(spec);
            shapes.push(s);
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

            return s;
        };


        /**
         * Create and add line shape
         * 
         * @function $3Dmol.GLViewer#addLine
         * @param {LineSpec} spec - Style specification, can specify dashed, dashLength, and gapLength
         * @return {$3Dmol.GLShape}
         @example
         $3Dmol.download("pdb:2ABJ",viewer,{},function(){
                  
                  viewer.setViewStyle({style:"outline"});
                  viewer.setStyle({chain:'A'},{sphere:{hidden:true}});
                  viewer.setStyle({chain:'D'},{sphere:{radius:3.0}});
                  viewer.setStyle({chain:'G'},{sphere:{colorscheme:'greenCarbon'}});
                  viewer.setStyle({chain:'J'},{sphere:{color:'blue'}});
                  viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
                  viewer.render();
              });

         */
        this.addLine = function(spec) {
            spec = spec || {};

            spec.start = getSelectionCenter(spec.start);
            spec.end = getSelectionCenter(spec.end);

            spec.wireframe = true;
            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            if (spec.dashed)
                s = addLineDashed(spec, s);
            else
                s.addLine(spec);
            shapes.push(s);
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

            return s;
        };
        
        /**
         * Unit Cell shape specification. 
         * @typedef UnitCellStyleSpec
         * @prop {LineStyleSpec} box - line style used to draw box
         * @prop {ArrowSpec} astyle - arrow specification of "a" axis
         * @prop {ArrowSpec} bstyle - arrow specification of "b" axis
         * @prop {ArrowSpec} cstyle - arrow specification of "c" axis
         * @prop {string} alabel - label for a axis
         * @prop {LabelSpec} alabelstyle - label style for a axis
         * @prop {string} blabel - label for b axis
         * @prop {LabelSpec} blabelstyle - label style for b axis
         * @prop {string} clabel - label for c axis
         * @prop {LabelSpec} clabelstyle - label style for c axis
         
         */
        
        /**
         * Create and add unit cell visualization.
         *
         * @function $3Dmol.GLViewer#addUnitCell
         * @param {GLModel} model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
         * @param {UnitCellStyleSpec} spec - visualization style
           @example

                $.get('data/1jpy.cif', function(data) {
                  let m = viewer.addModel(data);
                  viewer.addUnitCell(m, {box:{color:'purple'},alabel:'X',blabel:'Y',clabel:'Z',alabelstyle: {fontColor: 'black',backgroundColor:'white',inFront:true,fontSize:40},astyle:{color:'darkred', radius:5,midpos: -10}});
                  viewer.zoomTo();
                  viewer.render();
    });         
         */
        this.addUnitCell = function(model, spec) {
            model = this.getModel(model);
            spec = spec || {alabel: 'a', blabel: 'b', clabel: 'c'};
            
            spec.box = spec.box || {};
            spec.astyle = spec.astyle || {color: 'red', radius: 0.1,midpos: -1};
            spec.bstyle = spec.bstyle || {color: 'green', radius: 0.1,midpos: -1};
            spec.cstyle = spec.cstyle || {color: 'blue', radius: 0.1,midpos: -1};
            spec.alabelstyle = spec.alabelstyle || {fontColor: 'red',showBackground: false, alignment: 'center',inFront:false};
            spec.blabelstyle = spec.blabelstyle || {fontColor: 'green',showBackground: false, alignment: 'center',inFront:false};
            spec.clabelstyle = spec.clabelstyle || {fontColor: 'blue',showBackground: false, alignment: 'center',inFront:false};

            //clear any previous box
            if(model.unitCellObjects) {
                this.removeUnitCell(model);
            }
            model.unitCellObjects = {shapes:[],labels:[]};
            //calculate points
            var data = model.getCrystData();
            var matrix = null;
            if (data) {

                if (data.matrix) {
                    matrix = data.matrix;
                } else {
                    var a = data.a, b = data.b, c = data.c, alpha = data.alpha, beta = data.beta, gamma = data.gamma;
                    alpha = alpha * Math.PI/180.0;
                    beta = beta * Math.PI/180.0;
                    gamma = gamma * Math.PI/180.0;
            
                    var u, v, w;
            
                    u = Math.cos(beta);
                    v = (Math.cos(alpha) - Math.cos(beta)*Math.cos(gamma))/Math.sin(gamma);
                    w = Math.sqrt(Math.max(0, 1-u*u-v*v));
            
                    matrix = new $3Dmol.Matrix3(a, b*Math.cos(gamma), c*u, 
                                                0, b*Math.sin(gamma), c*v,
                                                0, 0,                 c*w); 
                }  
         
                var points = [  new $3Dmol.Vector3(0, 0, 0),
                                new $3Dmol.Vector3(1, 0, 0),
                                new $3Dmol.Vector3(0, 1, 0),
                                new $3Dmol.Vector3(0, 0, 1),
                                new $3Dmol.Vector3(1, 1, 0),
                                new $3Dmol.Vector3(0, 1, 1),
                                new $3Dmol.Vector3(1, 0, 1),
                                new $3Dmol.Vector3(1, 1, 1)  ];
                            
                if(data.matrix4) {
                    for (let i = 0; i < points.length; i++) {
                        if(data.size) points[i].multiplyVectors(points[i],data.size); //matrix is for unit vectors, not whole box
                        points[i] = points[i].applyMatrix4(data.matrix4);
                    }
                } else {
                    for (let i = 0; i < points.length; i++) {
                        points[i] = points[i].applyMatrix3(matrix);
                    }
                }
            
                //draw box
                if(spec.box && !spec.box.hidden) {
                    spec.box.wireframe = true;
                    var s = new $3Dmol.GLShape(spec.box);
                    s.shapePosition = shapes.length;
                
                    s.addLine({start: points[0], end: points[1]});
                    s.addLine({start: points[0], end: points[2]});
                    s.addLine({start: points[1], end: points[4]});
                    s.addLine({start: points[2], end: points[4]});
                
                    s.addLine({start: points[0], end: points[3]});
                    s.addLine({start: points[3], end: points[5]});
                    s.addLine({start: points[2], end: points[5]});
                
                    s.addLine({start: points[1], end: points[6]});
                    s.addLine({start: points[4], end: points[7]});
                    s.addLine({start: points[6], end: points[7]});
                
                    s.addLine({start: points[3], end: points[6]});
                    s.addLine({start: points[5], end: points[7]});
            
                    shapes.push(s);
                    model.unitCellObjects.shapes.push(s);
                    s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended
                }
                
                //draw arrows
                if(!spec.astyle.hidden) {
                    spec.astyle.start = points[0];
                    spec.astyle.end = points[1];
                    let arrow = this.addArrow(spec.astyle);
                    model.unitCellObjects.shapes.push(arrow);
                }

                if(!spec.bstyle.hidden) {
                    spec.bstyle.start = points[0];
                    spec.bstyle.end = points[2];
                    let arrow = this.addArrow(spec.bstyle);
                    model.unitCellObjects.shapes.push(arrow);
                }
                
                if(!spec.cstyle.hidden) {
                    spec.cstyle.start = points[0];
                    spec.cstyle.end = points[3];
                    let arrow = this.addArrow(spec.cstyle);
                    model.unitCellObjects.shapes.push(arrow);
                }
                
                if(spec.alabel) {
                    spec.alabelstyle.position = points[1];
                    let label = this.addLabel(spec.alabel, spec.alabelstyle);
                    model.unitCellObjects.labels.push(label);

                }
                if(spec.blabel) {
                    spec.blabelstyle.position = points[2];
                    let label = this.addLabel(spec.blabel, spec.blabelstyle);
                    model.unitCellObjects.labels.push(label);
                }
                if(spec.clabel) {
                    spec.clabelstyle.position = points[3];
                    let label = this.addLabel(spec.clabel, spec.clabelstyle);
                    model.unitCellObjects.labels.push(label);
                }

            }
            
        };
        
         /**
         * Remove unit cell visualization from model.
         *
         * @function $3Dmol.GLViewer#removeUnitCell
         * @param {GLModel} model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
           @example
                $.get('data/icsd_200866.cif', function(data) {
                  let m = viewer.addModel(data);
                  viewer.setStyle({sphere:{}})
                  viewer.addUnitCell();
                  viewer.zoomTo();
                  viewer.removeUnitCell();
                  viewer.render();
            });                  
         */
        this.removeUnitCell = function(model) {
            model = this.getModel(model);
            if(model.unitCellObjects) {
                let viewer = this;
                model.unitCellObjects.shapes.forEach(function(s) {viewer.removeShape(s);});
                model.unitCellObjects.labels.forEach(function(l) {viewer.removeLabel(l);});
            }
            delete model.unitCellObjects;
        };
        
         /**
         * Replicate atoms in model to form a super cell of the specified dimensions.
         * Original cell will be centered as much as possible.
         *
         * @function $3Dmol.GLViewer#replicateUnitCell
         * @param {integer} A - number of times to replicate cell in X dimension.
         * @param {integer} B - number of times to replicate cell in Y dimension.  If absent, X value is used.
         * @param {integer} C - number of times to replicate cell in Z dimension.  If absent, Y value is used.
         * @param {GLModel} model - Model with unit cell information (e.g., pdb derived).  If omitted uses most recently added model.
           @example
                $.get('data/icsd_200866.cif', function(data) {
                  let m = viewer.addModel(data);
                  viewer.setStyle({sphere:{scale:.25}})
                  viewer.addUnitCell();
                  viewer.zoomTo();
                  viewer.replicateUnitCell(3,2,1,m);
                  viewer.render();
            });                  
         */
        this.replicateUnitCell = function(A,B,C,model) {
            model = this.getModel(model);
            A = A || 3;
            B = B || A;
            C = C || B;
            let cryst = model.getCrystData();
            if(cryst) {
                const atoms = model.selectedAtoms({});
                const matrix = cryst.matrix;
                let makeoff = function(I) {
                    //alternate around zero: 1,-1,2,-2...
                    if(I%2 == 0) return -I/2;
                    else return Math.ceil(I/2);
                };

                for(let i = 0; i < A; i++) {
                    for(let j = 0; j < B; j++) {
                        for(let k = 0; k < C; k++) {
                            if(i == 0 && j == 0 && k == 0) continue; //actual unit cell
                            let offset = new $3Dmol.Vector3(makeoff(i),makeoff(j),makeoff(k));
                            offset.applyMatrix3(matrix);
                            
                            let newatoms = [];
                            for(let a = 0; a < atoms.length; a++) {
                                let newAtom = {};
                                for (let p in atoms[a]) {
                                    newAtom[p] = atoms[a][p];
                                }
                                newAtom.x += offset.x;
                                newAtom.y += offset.y;
                                newAtom.z += offset.z;
                                newatoms.push(newAtom);     
                            }
                            model.addAtoms(newatoms);
                         }
                    }
                }
            }
        };

        function addLineDashed(spec, s) {
            spec.dashLength = spec.dashLength || 0.5;
            spec.gapLength = spec.gapLength || 0.5;
            spec.start = spec.start || {};
            spec.end = spec.end || {};
            
            var p1 = new $3Dmol.Vector3(spec.start.x || 0,
                    spec.start.y || 0, spec.start.z || 0);
            var p2 = new $3Dmol.Vector3(spec.end.x,
                    spec.end.y || 0, spec.end.z || 0);
                    
            var dir = new $3Dmol.Vector3();
            var dash = new $3Dmol.Vector3();
            var gap = new $3Dmol.Vector3();
            var length, dashAmt, gapAmt;
            var temp = p1.clone();
            var drawn = 0;

            dir.subVectors(p2, p1);
            length = dir.length();
            dir.normalize();
            dash = dir.clone();
            gap = dir.clone();
            dash.multiplyScalar(spec.dashLength);
            gap.multiplyScalar(spec.gapLength);
            dashAmt = dash.length();
            gapAmt = gap.length();

            while (drawn < length) {
                if ((drawn + dashAmt) > length) { 
                    spec.start = p1;
                    spec.end = p2;
                    s.addLine(spec);
                    break;
                }
                temp.addVectors(p1, dash); 
                spec.start = p1;
                spec.end = temp;
                s.addLine(spec);
                p1 = temp.clone();
                drawn += dashAmt;

                temp.addVectors(p1, gap);
                p1 = temp.clone();   
                drawn += gapAmt;
            }
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended
                    
            return s;
        }

        

        /**
         * Add custom shape component from user supplied function
         * 
         * @function $3Dmol.GLViewer#addCustom
         * @param {CustomSpec} spec - Style specification
         * @return {$3Dmol.GLShape}
         @example
         function triangle(viewer) {
    var vertices = [];
    var normals = [];
    var colors = [];
    var r = 20;
    //triangle
    vertices.push(new $3Dmol.Vector3(0,0,0));
    vertices.push(new $3Dmol.Vector3(r,0,0));
    vertices.push(new $3Dmol.Vector3(0,r,0));
    
    normals.push(new $3Dmol.Vector3(0,0,1));
    normals.push(new $3Dmol.Vector3(0,0,1));
    normals.push(new $3Dmol.Vector3(0,0,1));
    
    colors.push({r:1,g:0,b:0});
    colors.push({r:0,g:1,b:0});
    colors.push({r:0,g:0,b:1});

    var faces = [ 0,1,2 ];
    
    var spec = {vertexArr:vertices, normalArr: normals, faceArr:faces,color:colors};
    viewer.addCustom(spec);
}
            triangle(viewer);
            viewer.render();
         */
        this.addCustom = function(spec) {
            spec = spec || {};
            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            s.addCustom(spec);
            shapes.push(s);
            s.finalize(); //finalize shape for memory efficiency, assume shape won't be extended

            return s;
        };

        /**
         * Construct isosurface from volumetric data in gaussian cube format
         * @function $3Dmol.GLViewer#addVolumetricData
         * @param {String} data - Input file contents 
         * @param {String} format - Input file format 
         * @param {IsoSurfaceSpec} or {VolumetricRenderSpec} spec - Shape style specification
         * @return {$3Dmol.GLShape}
         * 
         * @example

    
    $.get('data/bohr.cube', function(data) {
      
      viewer.addVolumetricData(data, "cube", {isoval: -0.01, color: "red", opacity: 0.95}); 
      viewer.setStyle({cartoon:{},stick:{}});
      viewer.zoomTo();
      viewer.render();
    });

                
         */
        this.addVolumetricData = function(data, format, spec) {
            spec = spec || {};
            
            var voldata = new $3Dmol.VolumeData(data, format);
            if(spec.transferfn) { //volumetric rendering
                return this.addVolumetricRender(voldata, spec);
            } else {
                return this.addIsosurface(voldata, spec);
            }
        };
        
        /**
         * Construct isosurface from volumetric data.  This is more flexible
     * than addVolumetricData, but can not be used with py3Dmol.
         * @function $3Dmol.GLViewer#addIsosurface
         * @param {$3Dmol.VolumeData} data - volumetric data
         * @param {IsoSurfaceSpec} spec - Shape style specification
         * @return {$3Dmol.GLShape}
         * 
         @example 
         $.get('../test_structs/benzene-homo.cube', function(data){
                  var voldata = new $3Dmol.VolumeData(data, "cube");
                  viewer.addIsosurface(voldata, {isoval: 0.01,
                                                 color: "blue"});
                  viewer.addIsosurface(voldata, {isoval: -0.01,
                                                 color: "red"});
                  viewer.zoomTo();
                  viewer.render();
                });
         */
        this.addIsosurface = function(data,  spec,callback) {
            spec = spec || {};
            var s = new $3Dmol.GLShape(spec);
            s.shapePosition = shapes.length;
            s.addIsosurface(data, spec, callback);
            shapes.push(s);
            return s;
        };
        
        /**
         * Create volumetric renderer for volumetricData
         * @function $3Dmol.GLViewer#addVolumetricRender
         * @param {$3Dmol.VolumeData} data - volumetric data
         * @param {VolumetricRenderSpec} spec - specification of volumetric render 
         * 
         * @return {$3Dmol.GLShape}
         * 
         */
        this.addVolumetricRender = function(data,  spec) {            
            spec = spec || {};
            var s = new $3Dmol.GLVolumetricRender(data, spec);
            s.shapePosition = shapes.length;     
            shapes.push(s);
            return s;
        };
        
        /**
         * Return true if volumetric rendering is supported (WebGL 2.0 required)
         *
         * @function $3Dmol.GLViewer#hasVolumetricRender
         * @return {boolean}
         */
        this.hasVolumetricRender = function() {
          return renderer.supportsVolumetric();  
        };
        
        /**
         * Enable/disable fog for content far from the camera
         *
         * @function $3Dmol.GLViewer#enableFog
         * @param {boolean} fog whether to enable or disable the fog
         */
        this.enableFog = function(fog){
            if (fog) {
                scene.fog = new $3Dmol.Fog(bgColor, 100, 200);
            } else {
                config.disableFog = true;
                show();
            }
        };

        /**
         * Sets the atomlists of all models in the viewer to specified frame.
         * Shapes and labels can also be displayed by frame.
         * Sets to last frame if framenum out of range
         * 
         * @function $3Dmol.GLViewer#setFrame
         * @param {number} framenum - fame index to use, starts at zero
         * @return {Promise}
         */
        this.setFrame = function (framenum) {
            viewer_frame = framenum;
            let viewer = this;
            return new Promise(function (resolve) {
                var modelMap = models.map(function (model) {
                    return model.setFrame(framenum,viewer);
                });
                Promise.all(modelMap)
                    .then(function() {resolve();});
            });
        };
        
        /**
         * Gets the current viewer frame.
         * 
         * @function $3Dmol.GLViewer#getFrame
         */
        this.getFrame = function () {
            return viewer_frame;
        };
                
        /**
         * Returns the number of frames that the model with the most frames in the viewer has
         * 
         * @function $3Dmol.GLViewer#getNumFrames
         * @return {number}
         */
        this.getNumFrames = function() {
            var mostFrames = 0;
            for (let i = 0; i < models.length; i++) {
                if (models[i].getNumFrames() > mostFrames) {
                    mostFrames = models[i].getNumFrames();
                }
            }
            for (let i = 0; i < shapes.length; i++) {
                if (shapes[i].frame && shapes[i].frame >= mostFrames) {
                    mostFrames = shapes[i].frame+1;
                }
            }         
            for (let i = 0; i < labels.length; i++) {
                if (labels[i].frame && labels[i].frame >= mostFrames) {
                    mostFrames = labels[i].frame+1;
                }
            }                     
            return mostFrames;
        };
        

        /**
         * Animate all models in viewer from their respective frames
         * @function $3Dmol.GLViewer#animate
         * @param {Object} options - can specify interval (speed of animation), loop (direction
         * of looping, 'backward', 'forward' or 'backAndForth'), step interval between frames ('step'), and reps (numer of repetitions, 0 indicates infinite loop)
         *      
         */
         
        this.animate = function(options) {
            incAnim();
            var interval = 100;
            var loop = "forward";
            var reps = 0;
            options = options || {};
            if (options.interval) {
                interval = options.interval;
            }
            if (options.loop) {
                loop = options.loop;
            }
            if (options.reps) {
                reps = options.reps;
            }
            var mostFrames = this.getNumFrames();
            var that = this;
            var currFrame = 0;
            var inc = 1;
            if (options.step) {
                inc = options.step;
                reps /= inc;
            }
            var displayCount = 0;
            var displayMax = mostFrames * reps;
            var time = new Date();
            var resolve, intervalID;
            var display = function(direction) {
                time = new Date();
                if (direction == "forward") {
                    that.setFrame(currFrame)
                    .then(function () {
                        currFrame = (currFrame + inc) % mostFrames;
                        resolve();
                    });
                }
                else if (direction == "backward") {
                    that.setFrame((mostFrames-1) - currFrame)
                    .then(function () {
                        currFrame = (currFrame + inc) % mostFrames;
                        resolve();
                    });
                }
                else { //back and forth
                    that.setFrame(currFrame)
                    .then(function () {
                        currFrame += inc;
                        inc *= (((currFrame % (mostFrames-1)) == 0) ? -1 : 1);
                        resolve();
                    });          
                }
            };            
            resolve = function() {
                that.render();
                if (++displayCount == displayMax || !that.isAnimated()) {
                    clearTimeout(intervalID);
                    animationTimers.delete(intervalID);
                    decAnim(); 
                }
                else {
                    var newInterval = interval - (new Date() - time);
                    newInterval = (newInterval>0)?newInterval:0;
                    animationTimers.delete(intervalID);
                    intervalID = setTimeout(display, newInterval, loop);
                    animationTimers.add(intervalID);
                }
            };

            intervalID = setTimeout(display, 0, loop);
            animationTimers.add(intervalID);
            return this;
        };
        
        /**
         * Stop animation of all models in viewer
         * @function $3Dmol.GLViewer#stopAnimate
         */
        this.stopAnimate = function() {
            animated = 0;
            animationTimers.forEach(function(timer) { clearTimeout(timer);});
            animationTimers = new Set();
            return this;
        };
        
        /**
         * Return true if viewer is currently being animated, false otherwise
         * @function $3Dmol.GLViewer#isAnimated
         * @return {boolean}
         */
        this.isAnimated = function() {
            return animated > 0;
        };
        

        /**
         * Create and add model to viewer, given molecular data and its format 
         * 
         * @function $3Dmol.GLViewer#addModel
         * @param {string} data - Input data
         * @param {string} format - Input format ('pdb', 'sdf', 'xyz', 'pqr', or 'mol2')
         * @param {ParserOptionsSpec} options - format dependent options. Attributes depend on the input file format.
         * @example
         

              viewer.setViewStyle({style:"outline"});
              $.get('data/1fas.pqr', function(data){
                  viewer.addModel(data, "pqr");
                  $.get("data/1fas.cube",function(volumedata){
                      viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.85,voldata: new $3Dmol.VolumeData(volumedata, "cube"), volscheme: new $3Dmol.Gradient.RWB(-10,10)},{});
                      
                  viewer.render();
                  });
                  viewer.zoomTo();
              });
         *
         * @return {$3Dmol.GLModel} 
         */
        this.addModel =  function(data, format, options) {
            if(options && !options.defaultcolors) {
                options.defaultcolors = defaultcolors;
                options.cartoonQuality = config.cartoonQuality;
            } else if(typeof(options) === 'undefined') {
                options = {defaultcolors:defaultcolors, cartoonQuality:config.cartoonQuality};
            }
            var m = new $3Dmol.GLModel(models.length, options);
            m.addMolData(data, format, options);
            models.push(m);

            return m;
        };
        
        /**
         * Given multimodel file and its format, add atom data to the viewer as separate models
         * and return list of these models
         * 
         * @function $3Dmol.GLViewer#addModels
         * @param {string} data - Input data
         * @param {string} format - Input format (see {@link FileFormats})
         * @return {Array<$3Dmol.GLModel>}
         */
        this.addModels = function(data, format, options) {
            options = options || {};
            options.multimodel = true;
            options.frames = true;

            var modelatoms = $3Dmol.GLModel.parseMolData(data, format, options);

            for (var i = 0; i < modelatoms.length; i++) {
                var newModel = new $3Dmol.GLModel(models.length, defaultcolors);
                newModel.setAtomDefaults(modelatoms[i]);
                newModel.addFrame(modelatoms[i]);
                newModel.setFrame(0);
                if(modelatoms.modelData)
                    newModel.setModelData(modelatoms.modelData[i]);
                newModel.setDontDuplicateAtoms(!options.duplicateAssemblyAtoms);
                models.push(newModel);
            }
            
            return models;
        };
        
        /**
         * Create and add model to viewer. Given multimodel file and its format, 
         * different atomlists are stored in model's frame
         * property and model's atoms are set to the 0th frame
         * 
         * @function $3Dmol.GLViewer#addModelsAsFrames
         * @param {string} data - Input data
         * @param {string} format - Input format (see {@link FileFormats})
         * @return {$3Dmol.GLModel}
         * 
         * @example
                $.get('../test_structs/multiple2.xyz', function(data){
                  viewer.addModelsAsFrames(data, "xyz");
                  viewer.animate({loop: "forward",reps: 1});
                  viewer.setStyle({stick:{colorscheme:'magentaCarbon'}});
                  viewer.zoomTo();
                  viewer.render();
              });
         */
        this.addModelsAsFrames = function(data, format, options) {
            options = options || {};
            options.multimodel = true;
            options.frames = true;
            var m = new $3Dmol.GLModel(models.length, defaultcolors);
            m.addMolData(data, format, options);
            models.push(m);

            return m;
        };
        
        /**
         * Create and add model to viewer. Given multimodel file and its format,
         * all atoms are added to one model
         * 
         * @function $3Dmol.GLViewer#addAsOneMolecule
         * @param {string} data - Input data
         * @param {string} format - Input format (see {@link FileFormats})
         * @return {$3Dmol.GLModel}
         @example
          

              $.get('../test_structs/multiple.sdf', function(data){
                  viewer.addAsOneMolecule(data, "sdf");
                  viewer.zoomTo();
                  viewer.render();
              });
         */
        this.addAsOneMolecule = function(data, format, options) {
            options = options || {};
            options.multimodel = true;
            options.onemol = true;
            var m = new $3Dmol.GLModel(models.length, defaultcolors);
            m.addMolData(data, format, options);
            models.push(m);
            
            return m;
        };
        

        /**
         * Delete specified model from viewer
         * 
         * @function $3Dmol.GLViewer#removeModel
         * @param {$3Dmol.GLModel} model
         */
        this.removeModel = function(model) {
            model = this.getModel(model);
            if (!model)
                return;
            model.removegl(modelGroup);
            delete models[model.getID()];
            // clear off back of model array
            while (models.length > 0
                    && typeof (models[models.length - 1]) === "undefined")
                models.pop();
            return this;
        };

        /** 
         * Delete all existing models
         * @function $3Dmol.GLViewer#removeAllModels
         */
        this.removeAllModels = function() {
            for (var i = 0; i < models.length; i++) {
                var model = models[i];
                if(model) model.removegl(modelGroup);

            }
            models.splice(0,models.length); //don't simply overwrite array in case linked
            return this;
        };

        /**
         * Export one or all of the loaded models into ChemDoodle compatible JSON.
         * @function $3Dmol.GLViewer#exportJSON
         * @param {boolean} includeStyles - Whether or not to include style information.
         * @param {number} modelID - Optional parameter for which model to export. If left out, export all of them.
         * @return {string}
         */
        this.exportJSON = function(includeStyles, modelID) {
            var object = {};
            if (modelID === undefined) {
                object.m = models.map(function(model) {
                    return model.toCDObject(includeStyles);
                });
            } else {
                object.m = [ models[modelID].toCDObject() ];
            }
            return JSON.stringify(object);
        };

        /** return a VRML string representation of the scene.  Include VRML header information
         * @function $3Dmol.GLViewer#exportVRML
         * @return VRML
         */
        this.exportVRML = function() {
            var savedmodelGroup = modelGroup;
            applyToModels("removegl",modelGroup); //cleanup
            modelGroup = new $3Dmol.Object3D();
            //rendering with plain mesh
            this.render(null, {supportsImposters: false, supportsAIA: false, regen: true});
            var ret = '#VRML V2.0 utf8\n' + modelGroup.vrml() + '\n';
            applyToModels("removegl",modelGroup); //cleanup
            modelGroup = savedmodelGroup;
            return ret;
        };
        
        /**
         * Create a new model from atoms specified by sel.
         * If extract, removes selected atoms from existing models 
         * 
         * @function $3Dmol.GLViewer#createModelFrom
         * @param {Object} sel - Atom selection specification
         * @param {boolean=} extract - If true, remove selected atoms from existing models
         * @return {$3Dmol.GLModel}
         */
        this.createModelFrom = function(sel, extract) {
            var m = new $3Dmol.GLModel(models.length, defaultcolors);
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

        function applyToModels(func, sel, value1, value2, value3, value4, value5) {
            
            //apply func to all models that are selected by sel with value1 and 2
            var ms = getModelList(sel);
            for (var i = 0; i < ms.length; i++) {
                ms[i][func](sel, value1, value2, value3, value4, value5);
            }
        }

        /**
         * Set style properties to all selected atoms
         * 
         * @function $3Dmol.GLViewer#setStyle
         * @param {AtomSelectionSpec} sel - Atom selection specification
         * @param {AtomStyleSpec} style - Style spec to apply to specified atoms
         * 
         * @example
         viewer.setBackgroundColor(0xffffffff);
       $3Dmol.download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
        m.setStyle({chain:'A'},{'cartoon':{color:'spectrum'}});
        m.setStyle({chain:'C'},{'cartoon':{style:'trace',color:'blue'}});
        m.setStyle({chain:'E'},{'cartoon':{tubes:true,arrows:true,color:'green',opacity:0.75}});
        m.setStyle({chain:'B'},{'cartoon':{color:'red',opacity:0.5}});
        m.setStyle({chain:'D'},{'cartoon':{style:'trace',color:'grey',opacity:0.75}});
        m.setStyle({chain:'F'},{'cartoon':{arrows:true,color:'white'}});
       // viewer.addStyle({chain:'B'},{line:{}});
       viewer.zoomTo();
       viewer.render();
    });
         */
        this.setStyle = function(sel, style) {
            if(typeof(style) === 'undefined') {
                //if a single argument is provided, assume it is a style and select all
                style = sel;
                sel = {};
            }
            
            applyToModels("setStyle", sel, style, false);
            return this;
        };

        /**
         * Add style properties to all selected atoms
         * 
         * @function $3Dmol.GLViewer#addStyle
         * @param {AtomSelectionSpec} sel - Atom selection specification
         * @param {AtomStyleSpec} style - style spec to add to specified atoms
         @example
         
       $3Dmol.download('pdb:5IRE',viewer,{doAssembly: false},function(m) {
       viewer.setStyle({cartoon:{}});
       //keep cartoon style, but show thick sticks for chain A
       viewer.addStyle({chain:'A'},{stick:{radius:.5,colorscheme:"magentaCarbon"}});
       viewer.zoomTo();
       viewer.render();
       });
         */
        this.addStyle = function(sel, style) {
            if(typeof(style) === 'undefined') {
                //if a single argument is provided, assume it is a style and select all
                style = sel;
                sel = {};
            }
            applyToModels("setStyle", sel, style, true);
            return this;
        };


        /**
         * Set click-handling properties to all selected atomsthis.
         * 
         * @function $3Dmol.GLViewer#setClickable
         * @param {AtomSelectionSpec} sel - atom selection to apply clickable settings to
         * @param {boolean} clickable - whether click-handling is enabled for the selection
         * @param {function} callback - function called when an atom in the selection is clicked
         * 
         * @example       
            $3Dmol.download("cid:307900",viewer,{},function(){
                              
                   viewer.setStyle({},{sphere:{}});                
                   viewer.setClickable({},true,function(atom,viewer,event,container) {
                       viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8});
                   });
                   viewer.render();
        });
         */
        this.setClickable = function(sel, clickable, callback) {
            applyToModels("setClickable", sel, clickable, callback);
            return this;
        };
        /** Set hoverable and callback of selected atoms
         * 
         * @function $3Dmol.GLViewer#setHoverable
         * @param {AtomSelectionSpec} sel - atom selection to apply hoverable settings to
         * @param {boolean} hoverable - whether hover-handling is enabled for the selection
         * @param {function} hover_callback - function called when an atom in the selection is hovered over
         * @param {function} unhover_callback - function called when the mouse moves out of the hover area  
        @example             
        $3Dmol.download("pdb:1ubq",viewer,{},function(){
                          
               viewer.setHoverable({},true,function(atom,viewer,event,container) {
                   if(!atom.label) {
                    atom.label = viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                   }
               },
               function(atom) { 
                   if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                   }
                }
               );
               viewer.setStyle({},{stick:{}});               
               viewer.render();
        });

         */
        this.setHoverable = function(sel,hoverable,hover_callback,unhover_callback){
            applyToModels("setHoverable", sel,hoverable, hover_callback,unhover_callback);
            return this;
        };
        
        /**
         * If  atoms have dx, dy, dz properties (in some xyz files), vibrate populates each model's frame property based on parameters.
         * Models can then be animated
         * 
         * @function $3Dmol.GLViewer#vibrate
         * @param {number} numFrames - number of frames to be created, default to 10
         * @param {number} amplitude - amplitude of distortion, default to 1 (full)
         * @param {boolean} bothWays - if true, extend both in positive and negative directions by numFrames
         * @param {ArrowSpec} arrowSpec - specification for drawing animated arrows. If color isn't specified, atom color (sphere, stick, line preference) is used.
         */
        this.vibrate = function(numFrames, amplitude, bothways, arrowSpec) {
            applyToModels("vibrate", numFrames, amplitude, bothways, this, arrowSpec);
            return this;
        };
        
        /**
         * @function $3Dmol.GLViewer#setColorByProperty
         * @param {AtomSelectionSpec} sel
         * @param {type} prop
         * @param {type} scheme
         */
        this.setColorByProperty = function(sel, prop, scheme, range) {
            applyToModels("setColorByProperty", sel, prop, scheme, range);
            return this;
        };

        /**
         * @function $3Dmol.GLViewer#setColorByElement
         * @param {AtomSelectionSpec} sel
         * @param {type} colors
         */
        this.setColorByElement = function(sel, colors) {
            applyToModels("setColorByElement", sel, colors);
            return this;
        };

        /**
         * 
         * @param {AtomSpec[]} atomlist
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
                ret.push(atom);
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
         * @param {AtomSpec[]} atomlist
         * @param {AtomSpec[]} atomstoshow
         * @return {Array}
         */
        var carveUpExtent = function(extent, atomlist, atomstoshow) {
            var ret = [];

            var index2atomlist = {}; //map from atom.index to position in atomlist
            for(var i = 0, n = atomlist.length; i < n; i++) {
                index2atomlist[atomlist[i].index] = i;
            }
            
            var atomsToListIndex = function(atoms) {
            //return a list of indices into atomlist
                var ret = [];
                for(var i = 0, n = atoms.length; i < n; i++) {
                    if(atoms[i].index in index2atomlist)
                        ret.push(index2atomlist[atoms[i].index]);
                }
                return ret;
            };
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
            for (let i = 0, n = splits.length; i < n; i++) {
                let e = copyExtent(splits[i]);
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
                    atoms : atomsToListIndex(atoms),
                    toshow : atomsToListIndex(toshow)
                });
            }

            return ret;
        };

        // create a mesh defined from the passed vertices and faces and material
        // Just create a single geometry chunk - broken up whether sync or not
        /**
         * 
         * @param {AtomSpec[]} atoms
         * @param {{vertices:number,faces:number}}
         *            VandF
         * @param {$3Dmol.MeshLambertMaterial}
         *            mat
         * @return {$3Dmol.Mesh}
         */
        var generateSurfaceMesh = function(atoms, VandF, mat) {
            var geo = new $3Dmol.Geometry(true);
            // Only one group per call to generate surface mesh (addSurface
            // should split up mesh render)
            var geoGroup = geo.updateGeoGroup(0);

            // set colors for vertices
            var colors = [];
            for (let i = 0, il = atoms.length; i < il; i++) {
                var atom = atoms[i];
                if (atom) {
                    if (typeof (atom.surfaceColor) != "undefined") {
                        colors[i] = atom.surfaceColor;
                    } else if (atom.color) // map from atom
                        colors[i] = $3Dmol.CC.color(atom.color);
                }
            }
            
            var vertexArray = geoGroup.vertexArray;

            // reconstruct vertices and faces
            var v = VandF.vertices;
            for (let i = 0, il = v.length; i < il; i++) {
                let offset = geoGroup.vertices * 3;
                vertexArray[offset] = v[i].x;
                vertexArray[offset + 1] = v[i].y;
                vertexArray[offset + 2] = v[i].z;
                geoGroup.vertices++;                
            }

            //set colorArray of there are per-atom colors
            var colorArray = geoGroup.colorArray;
            
            if(mat.voldata && mat.volscheme) {
                //convert volumetric data into colors
                var scheme = mat.volscheme;
                var voldata = mat.voldata;
                var range = scheme.range() || [-1,1];
                for (let i = 0, il = v.length; i < il; i++) {
                    let val = voldata.getVal(v[i].x,v[i].y,v[i].z);
                    let col =  $3Dmol.CC.color(scheme.valueToHex(val, range));
                    let offset = i * 3;
                    colorArray[offset] = col.r;
                    colorArray[offset + 1] = col.g;
                    colorArray[offset + 2] = col.b;
                }
            }
            else if(colors.length > 0) { //have atom colors
                for (let i = 0, il = v.length; i < il; i++) {
                    let A = v[i].atomid;
                    let offsetA = i * 3;

                    colorArray[offsetA] = colors[A].r;
                    colorArray[offsetA + 1] = colors[A].g;
                    colorArray[offsetA + 2] = colors[A].b;
                }
            }
            
            var faces = VandF.faces;
            geoGroup.faceidx = faces.length;// *3;
            geo.initTypedArrays();

            var verts = geoGroup.vertexArray;
            var normalArray = geoGroup.normalArray;
            var vA, vB, vC, norm;

            // Setup colors, faces, and normals
            for (let i = 0, il = faces.length; i < il; i += 3) {

                // var a = faces[i].a, b = faces[i].b, c = faces[i].c;
                var a = faces[i], b = faces[i + 1], c = faces[i + 2];
                var offsetA = a * 3, offsetB = b * 3, offsetC = c * 3;

                // setup Normals
                // todo - calculate normals in parallel code
                vA = new $3Dmol.Vector3(verts[offsetA], verts[offsetA + 1],
                        verts[offsetA + 2]);
                vB = new $3Dmol.Vector3(verts[offsetB], verts[offsetB + 1],
                        verts[offsetB + 2]);
                vC = new $3Dmol.Vector3(verts[offsetC], verts[offsetC + 1],
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
            var mesh = new $3Dmol.Mesh(geo, mat);
            mesh.doubleSided = true;        
            return mesh;
        };

        // do same thing as worker in main thread
        /**
         * 
         * @param {$3Dmol.SurfaceType}
         *            type
         * @param {Array}
         *            expandedExtent
         * @param {Array}
         *            extendedAtoms
         * @param {Array}
         *            atomsToShow
         * @param {AtomSpec[]} atoms
         * @param {number}
         *            vol
         * @return {Object}
         */
        var generateMeshSyncHelper = function(type, expandedExtent,
                extendedAtoms, atomsToShow, atoms, vol) {
//            var time = new Date();
            var ps = new $3Dmol.ProteinSurface();
            ps.initparm(expandedExtent, (type === 1) ? false : true, vol);

//            var time2 = new Date();
            //console.log("initialize " + (time2 - time) + "ms");

            ps.fillvoxels(atoms, extendedAtoms);

//            var time3 = new Date();
            //console.log("fillvoxels " + (time3 - time2) + "  " + (time3 - time) + "ms");

            ps.buildboundary();

            if (type == $3Dmol.SurfaceType.SES || type == $3Dmol.SurfaceType.MS) {
                ps.fastdistancemap();
                ps.boundingatom(false);
                ps.fillvoxelswaals(atoms, extendedAtoms);
            }

//            var time4 = new Date();
            //console.log("buildboundaryetc " + (time4 - time3) + "  " + (time4 - time) + "ms");

            ps.marchingcube(type);

//            var time5 = new Date();
            //console.log("marching cube " + (time5 - time4) + "  "+ (time5 - time) + "ms");

            return ps.getFacesAndVertices(atomsToShow);
        };

        /**
         * 
         * @param {matSpec}
         *            style
         * @return {$3Dmol.MeshLambertMaterial}
         */
        function getMatWithStyle(style) {
            var mat = new $3Dmol.MeshLambertMaterial();
            mat.vertexColors = $3Dmol.VertexColors;

            for ( var prop in style) {
                if (prop === "color" || prop === "map") {
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

        
        /**
         * Adds an explicit mesh as a surface object.
         * @function $3Dmol.GLViewer#addMesh
         * @param {$3Dmol.Mesh}
         *            mesh
         * @param {Object}
         *            style
         * @returns {number} surfid
         */
        this.addMesh = function(mesh) {
            var surfobj = {
                geo : mesh.geometry,
                mat : mesh.material,
                done : true,
                finished : false //the rendered finishes surfaces when they are done
            };
            var surfid = nextSurfID();
            surfaces[surfid] = surfobj;
            return surfid;
        };

        //return a shallow copy of list l, e.g., for atoms so we can
        //ignore superficial changes (ie surfacecolor, position) that happen
        //while we're surface building
        var shallowCopy = function(l) {
            var ret = [];
            let length = l.length;
            for(let i = 0; i < length; i++) {
                ret[i] = $3Dmol.extend({},l[i]);
            }
            return ret;
        };

        var surfaceTypeMap={
            "VDW":$3Dmol.SurfaceType.VDW,
            "MS":$3Dmol.SurfaceType.MS,
            "SAS":$3Dmol.SurfaceType.SAS,
            "SES":$3Dmol.SurfaceType.SES
        };

        
        /**
         * Add surface representation to atoms
         * @function $3Dmol.GLViewer#addSurface
         * @param {$3Dmol.SurfaceType|string} type - Surface type (VDW, MS, SAS, or SES)
         * @param {SurfaceStyleSpec} style - optional style specification for surface material (e.g. for different coloring scheme, etc)
         * @param {AtomSelectionSpec} atomsel - Show surface for atoms in this selection
         * @param {AtomSelectionSpec} allsel - Use atoms in this selection to calculate surface; may be larger group than 'atomsel' 
         * @param {AtomSelectionSpec} focus - Optionally begin rendering surface specified atoms
         * @param {function} surfacecallback - function to be called after setting the surface
         * @return {Promise} promise - Returns a promise that ultimately resovles to the surfid.  Returns surfid immediately if surfacecallback is specified.  Returned promise has a surfid field for immediate access.
         */
        this.addSurface = function(type, style, atomsel, allsel, focus, surfacecallback) {
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
            
            //surfacecallback gets called when done
            var surfid = nextSurfID();
            var mat = null;
            if(typeof type =="string"){
                if(surfaceTypeMap[type]!== undefined)
                    type = surfaceTypeMap[type];
                else{
                    console.log("Surface type : " + type + " is not recognized");
                } 
            }
            else if(type===undefined){
                type = $3Dmol.SurfaceType.VDW; //default
            }
            // atoms specified by this selection
            var atomlist = null, focusSele = null;
            //TODO: currently generating a shallow copy to avoid problems when atoms are chagned
            //during surface generation - come up with a better solution
            var atomsToShow = shallowCopy(getAtomsFromSel(atomsel));
            if(!allsel) {
                atomlist = atomsToShow;
            }
            else {
                atomlist = shallowCopy(getAtomsFromSel(allsel));
            }
            
            $3Dmol.adjustVolumeStyle(style);
            var symmetries = false;
            var n;
            for (n = 0; n < models.length; n++) { 
                if(models[n]) {
                    var symMatrices = models[n].getSymmetries();
                    if (symMatrices.length > 1 || (symMatrices.length == 1 && !(symMatrices[0].isIdentity()))) {
                        symmetries = true;
                        break;
                    }
                }
            }

            var addSurfaceHelper = function addSurfaceHelper(surfobj, atomlist, atomsToShow) {
                //function returns promise with surfid resolved
                if(!focus) {
                    focusSele = atomsToShow;
                } else {
                    focusSele = shallowCopy(getAtomsFromSel(focus));
                }

                var atom;
//                var time = new Date();
                var extent = $3Dmol.getExtent(atomsToShow, true);
                if (style.map && style.map.prop) {
                    // map color space using already set atom properties
                    /** @type {AtomSpec} */
                    var prop = style.map.prop;
                    /** @type {Gradient} */
                    var scheme = style.map.scheme || style.map.gradient || new $3Dmol.Gradient.RWB();
                    var range = scheme.range();
                    if (!range) {
                        range = $3Dmol.getPropertyRange(atomsToShow, prop);
                    }
                    style.colorscheme = {prop: prop, gradient: scheme};

                }
                
                //cache surface color on each atom
                for (let i = 0, il = atomlist.length; i < il; i++) {
                    atom = atomlist[i];
                    atom.surfaceColor = $3Dmol.getColorFromStyle(atom, style);
                }                

                var totalVol = volume(extent); // used to scale resolution
                var extents = carveUpExtent(extent, atomlist, atomsToShow);

                if (focusSele && focusSele.length && focusSele.length > 0) {
                    var seleExtent = $3Dmol.getExtent(focusSele, true);
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

                //console.log("Extents " + extents.length + "  "+ (+new Date() - time) + "ms");


                var reducedAtoms = [];
                // to reduce amount data transfered, just pass x,y,z,serial and elem
                for (let i = 0, il = atomlist.length; i < il; i++) {
                    atom = atomlist[i];
                    reducedAtoms[i] = {
                        x : atom.x,
                        y : atom.y,
                        z : atom.z,
                        serial : i,
                        elem : atom.elem
                    };
                }

                var sync = !!($3Dmol.syncSurface);
                if (sync) { // don't use worker, still break up for memory purposes

                    // to keep the browser from locking up, call through setTimeout
                    var callSyncHelper = function callSyncHelper(i) {
                        return new Promise(function(resolve) { 
                            var VandF = generateMeshSyncHelper(type, extents[i].extent,
                                    extents[i].atoms, extents[i].toshow, reducedAtoms,
                                    totalVol);
                            //complicated surfaces sometimes have > 2^16 vertices
                            var VandFs = $3Dmol.splitMesh({vertexArr:VandF.vertices, faceArr:VandF.faces});
                            for(var vi=0,vl=VandFs.length;vi<vl;vi++){
                                VandF={vertices:VandFs[vi].vertexArr,
                                        faces:VandFs[vi].faceArr};                            
                                var mesh = generateSurfaceMesh(atomlist, VandF, mat);
                                $3Dmol.mergeGeos(surfobj.geo, mesh);
                            }
                            _viewer.render();
                            resolve();
                        });
                    };
                    var promises = [];
                    for (let i = 0; i < extents.length; i++) {
                        promises.push(callSyncHelper(i));
                    }
                    return Promise.all(promises)
                    .then(function() {
                        surfobj.done = true;
                        return Promise.resolve(surfid);
                    });

                    // TODO: Asynchronously generate geometryGroups (not separate
                    // meshes) and merge them into a single geometry
                } else { // use worker

                   var workers = [];
                    if (type < 0)
                        type = 0; // negative reserved for atom data
                    for (let i = 0, il = numWorkers; i < il; i++) {
                        var w = new Worker($3Dmol.SurfaceWorker);
                        workers.push(w);
                        w.postMessage({
                            'type' : -1,
                            'atoms' : reducedAtoms,
                            'volume' : totalVol
                        });
                    }
                    
                    return new Promise(function(resolve,reject) {
                        var cnt = 0;

                        var rfunction = function(event) {
                            var VandFs = $3Dmol.splitMesh({vertexArr:event.data.vertices,
                                                           faceArr:event.data.faces});
                            for(var i=0,vl=VandFs.length;i<vl;i++){
                                var VandF={vertices:VandFs[i].vertexArr,
                                           faces:VandFs[i].faceArr};
                                var mesh = generateSurfaceMesh(atomlist, VandF, mat);
                                $3Dmol.mergeGeos(surfobj.geo, mesh);
                            }
                            _viewer.render();

                        //    console.log("async mesh generation " + (+new Date() - time) + "ms");
                            cnt++;
                            if (cnt == extents.length) {
                                surfobj.done = true;
                                resolve(surfid); //caller of helper will resolve callback if present
                            }
                        };

                        var efunction = function(event) {
                            console.log(event.message + " (" + event.filename + ":" + event.lineno + ")");
                            reject(event);
                        };

                        for (let i = 0; i < extents.length; i++) {
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
                    });
                }                
            };
            
            style = style || {};
            mat = getMatWithStyle(style);
            var surfobj = [];
            var promise = null;
            if (symmetries) { //do preprocessing
                var modelsAtomList = {};
                var modelsAtomsToShow = {};
                for (n = 0; n < models.length; n++) {
                    modelsAtomList[n] = [];
                    modelsAtomsToShow[n] = [];
                }
                for (n = 0; n < atomlist.length; n++) {
                    modelsAtomList[atomlist[n].model].push(atomlist[n]);
                }
                for (n = 0; n < atomsToShow.length; n++) {
                    modelsAtomsToShow[atomsToShow[n].model].push(atomsToShow[n]);
                }
                var promises = [];
                for (n = 0; n < models.length; n++) {
                    if(modelsAtomsToShow[n].length > 0) {
                        surfobj.push({
                            geo : new $3Dmol.Geometry(true),
                            mat : mat,
                            done : false,
                            finished : false,
                            symmetries : models[n].getSymmetries()
                        // also webgl initialized
                        });
                        promises.push(addSurfaceHelper(surfobj[n], modelsAtomList[n], modelsAtomsToShow[n]));
                    }
                }
                promise = Promise.all(promises);
            }
            else {
                surfobj.push({
                    geo : new $3Dmol.Geometry(true),
                    mat : mat,
                    done : false,
                    finished : false,
                    symmetries : [new $3Dmol.Matrix4()]
                });
                promise = addSurfaceHelper(surfobj[surfobj.length-1], atomlist, atomsToShow);
            }
            surfaces[surfid] = surfobj;
            promise.surfid = surfid;
            if(surfacecallback && typeof(surfacecallback) == "function") {
                promise.then(function(surfid) {
                    surfacecallback(surfid);
                });
                return surfid;
            }
            else return promise;
        };

        /**
         * Set the surface material to something else, must render change
        *  @function $3Dmol.GLViewer#setSurfaceMaterialStyle
         * @param {number} surf - Surface ID to apply changes to
         * @param {SurfaceStyleSpec} style - new material style specification
         @example
         $.get("data/9002806.cif",function(data){
            viewer.addModel(data);
            viewer.setStyle({stick:{}});
            let surf = viewer.addSurface("SAS");
            surf.then(function() {
                viewer.setSurfaceMaterialStyle(surf.surfid, {color:'blue',opacity:0.5});
                viewer.render();
                });
           });
         */ 
        this.setSurfaceMaterialStyle = function(surf, style) {
            $3Dmol.adjustVolumeStyle(style);
            if (surfaces[surf]) {
                var surfArr = surfaces[surf];
                for (var i = 0; i < surfArr.length; i++) {
                    var mat = surfArr[i].mat = getMatWithStyle(style);
                    surfArr[i].mat.side = $3Dmol.FrontSide;
                    if(style.color) {
                        surfArr[i].mat.color = style.color;
                        surfArr[i].geo.colorsNeedUpdate = true;
                        const c = $3Dmol.CC.color(style.color);
                        surfArr[i].geo.setColors(function() { return c;});
                    }
                    else if(mat.voldata && mat.volscheme) {
                        //convert volumetric data into colors
                        const scheme = mat.volscheme;
                        const voldata = mat.voldata;
                        const cc = $3Dmol.CC;
                        const range = scheme.range() || [-1,1];
                        surfArr[i].geo.setColors(function(x,y,z) {
                            let val = voldata.getVal(x,y,z);
                            let col =  cc.color(scheme.valueToHex(val, range));
                            return col;
                        });
                    }
                    surfArr[i].finished = false; // trigger redraw
                }
            }
            return this;
        };

        /**
         * Remove surface with given ID
         * @function $3Dmol.GLViewer#removeSurface
         * @param {number} surf - surface id
         */
        this.removeSurface = function(surf) {
            var surfArr = surfaces[surf];
            for (var i = 0; i < surfArr.length; i++) {
                if (surfArr[i] && surfArr[i].lastGL) {
                    if (surfArr[i].geo !== undefined)
                        surfArr[i].geo.dispose();
                    if (surfArr[i].mat !== undefined)
                        surfArr[i].mat.dispose();
                    modelGroup.remove(surfArr[i].lastGL); // remove from scene
                }
            }
            delete surfaces[surf];
            show();
            return this;
        };
        
        /** Remove all surfaces.
         * @function $3Dmol.GLViewer#removeAllSurfaces */
        this.removeAllSurfaces = function() {
            for (var n in  surfaces) {
                if(!surfaces.hasOwnProperty(n)) continue;
                var surfArr = surfaces[n];
                for(var i = 0; i < surfArr.length; i++) {
                    if (surfArr[i] && surfArr[i].lastGL) {
                        if (surfArr[i].geo !== undefined)
                            surfArr[i].geo.dispose();
                        if (surfArr[i].mat !== undefined)
                            surfArr[i].mat.dispose();
                        modelGroup.remove(surfArr[i].lastGL); // remove from scene
                    }
                }
                delete surfaces[n];
            }
            show();
            return this;
        };

        /** return Jmol moveto command to position this scene */
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

        /** Clear scene of all objects 
         * @function $3Dmol.GLViewer#clear
         * */
        this.clear = function() {
            this.removeAllSurfaces();
            this.removeAllModels();
            this.removeAllLabels();
            this.removeAllShapes();
            show();
            return this;
        };

        // props is a list of objects that select certain atoms and enumerate
        // properties for those atoms
        /**
         * Add specified properties to all atoms matching input argument
         * @function $3Dmol.GLViewer#mapAtomProperties
         * @param {Object} props, either array of atom selectors with associated props, or function that takes atom and sets its properties
         * @param {AtomSelectionSpec} sel  - subset of atoms to work on - model selection must be specified here
             @example 
             $.get('../test_structs/b.sdf', function(data){
                      viewer.addModel(data,'sdf');
                      let props = [];
                      //make the atom index a property x
                      for(let i = 0; i < 8; i++) {
                        props.push({index:i,props:{'x':i}});
                      }
                      viewer.mapAtomProperties(props);
                      viewer.setStyle({sphere:{colorscheme:{gradient:'roygb',prop:'x',min:0,max:8}}});
                      viewer.zoomTo();
                      viewer.render();
                    });         
         */
        this.mapAtomProperties = function(props, sel) {
            sel = sel || {};
            var atoms = getAtomsFromSel(sel);
            
            if(typeof(props) == "function") {
                for (let a = 0, numa = atoms.length; a < numa; a++) {
                    let atom = atoms[a];
                    props(atom);
                }
            }
            else {
                for (let a = 0, numa = atoms.length; a < numa; a++) {
                    var atom = atoms[a];
                    for (let i = 0, n = props.length; i < n; i++) {
                        let prop = props[i];
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
            }
            return this;
        };


        /**
         * Synchronize this view matrix of this viewer to the passed viewer.
         * When the viewpoint of this viewer changes, the other viewer will
         * be set to this viewer's view.
         * @function $3Dmol.GLViewer#linkViewer
         * @param {$3Dmol.GLViewer} otherview 
         */
        this.linkViewer = function(otherviewer) {
           linkedViewers.push(otherviewer);
           return this;
        };
        

        try {
            if (typeof (callback) === "function")
                callback(this);
        } catch (e) {
            // errors in callback shouldn't invalidate the viewer
            console.log("error with glviewer callback: " + e);
        }

        /**
         * Return the z distance between the model and the camera
         * @function $3Dmol.GLViewer#getPerceivedDistance
         * @return {number} distance
         */
        this.getPerceivedDistance = function() {
            return CAMERA_Z - rotationGroup.position.z;
        };

        /**
         * Set the distance between the model and the camera
         * Essentially zooming. Useful while stereo rendering.
         * @function $3Dmol.GLViewer#setPerceivedDistance
         */
        this.setPerceivedDistance = function(dist) {
            rotationGroup.position.z = CAMERA_Z - dist;
        };

        /**
         * Used for setting an approx value of eyeSeparation. Created for calling by StereoViewer object
         * @function $3Dmol.GLViewer#setAutoEyeSeparation
         * @return {number} camera x position
         */
        this.setAutoEyeSeparation = function(isright) {
            var dist = this.getPerceivedDistance();
            if (isright || camera.position.x > 0) //setting a value of dist*tan(5)
                camera.position.x = dist*Math.tan(Math.PI / 180.0 * 5.0);
            else
                camera.position.x = -dist*Math.tan(Math.PI / 180.0 * 5.0);
            camera.lookAt(new $3Dmol.Vector3(0,0,rotationGroup.position.z));
            return camera.position.x;
        };
        
        /**
         * Set the default cartoon quality for newly created models.  Default is 5.
         * Current models are not affected.
         * @number quality, higher results in higher resolution renders
         * @function $3Dmol.GLViewer#setDefaultCartoonQuality
         */
        this.setDefaultCartoonQuality = function(val) {
            config.cartoonQuality = val;
        };

    }

    return GLViewer;

})();

$3Dmol.glmolViewer = $3Dmol.GLViewer;

