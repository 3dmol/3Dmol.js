//a molecular viewer based on GLMol

var WebMol = WebMol || {};

//Adapted from the text sprite example from http://stemkoski.github.io/Three.js/index.html

// function for drawing rounded rectangles
var roundRect = function(ctx, x, y, w, h, r) {
    //ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
    ctx.beginPath();
    ctx.moveTo(x+r, y);
    ctx.lineTo(x+w-r, y);
    ctx.quadraticCurveTo(x+w, y, x+w, y+r);
    ctx.lineTo(x+w, y+h-r);
    ctx.quadraticCurveTo(x+w, y+h, x+w-r, y+h);
    ctx.lineTo(x+r, y+h);
    ctx.quadraticCurveTo(x, y+h, x, y+h-r);
    ctx.lineTo(x, y+r);
    ctx.quadraticCurveTo(x, y, x+r, y);
    ctx.closePath();
    ctx.fill();
    //ctx.stroke();   
};

WebMol.LabelCount = 0;

WebMol.Label = function(message, parameters) {
        
    this.id = WebMol.LabelCount++;    
    this.stylespec = parameters || {};

    var canvas = document.createElement('canvas');
    
    this.context = canvas.getContext('2d');

    this.sprite = null;
    this.text = message;
    
};

//TODO: How should this class be structured?
WebMol.Label.prototype = {
    
    constructor : WebMol.Label,
    
    setContext : function() {

        var font = this.stylespec.font ? 
            this.stylespec.font : "Arial";
    
        var fontsize = this.stylespec.fontsize ? 
            this.stylespec.fontsize : 20;
            
        var fontcolor = this.stylespec.fontcolor ?
            this.stylespec.fontcolor : { r:255, g:255, b:255, a:1.0};
    
        var borderThickness = this.stylespec.borderThickness ? 
            this.stylespec.borderThickness : 4;
    
        var borderColor = this.stylespec.borderColor ?
            this.stylespec.borderColor : { r:0, g:0, b:0, a:1.0 };
    
        var backgroundColor = this.stylespec.backgroundColor ?
            this.stylespec.backgroundColor : { r:0, g:0, b:0, a:1.0 };
            
        var position = this.stylespec.position ?
            this.stylespec.position : { x:-10, y:1, z:1 };
        
        //Should labels always be in front of model? 
        var inFront = (this.stylespec.inFront !== undefined) ?
            this.stylespec.inFront : true;
            
        var spriteAlignment = WebMol.SpriteAlignment.topLeft;
    
        // get size data (height depends only on font size)
        
        var metrics = this.context.measureText(this.text);
        var textWidth = metrics.width;
        var textHeight = textWidth * 100/260;
        
        textWidth = textWidth / textHeight * fontsize;
        textHeight = fontsize;
        
        
        this.context.font = textHeight + "px " + font;
        
        // background color
        this.context.fillStyle   = "rgba(" + backgroundColor.r + "," + backgroundColor.g + ","
                                                                  + backgroundColor.b + "," + backgroundColor.a + ")";
        // border color
        this.context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + ","
                                                                  + borderColor.b + "," + borderColor.a + ")";
    
        this.context.lineWidth = borderThickness;
        roundRect(this.context, borderThickness/2, borderThickness/2, textWidth + borderThickness, textHeight * 1.4 + borderThickness, 6);
        // 1.4 is extra height factor for text below baseline: g,j,p,q.
    
        // text color
        this.context.fillStyle = "rgba(" + fontcolor.r + "," + fontcolor.g + ","
                                                                + fontcolor.b + "," + fontcolor.a + ")";
    
        this.context.fillText(this.text, borderThickness, textHeight + borderThickness, textWidth);
        
        // canvas contents will be used for a texture
        var texture = new WebMol.Texture(this.context.canvas);
        texture.needsUpdate = true;
    
        var spriteMaterial = new WebMol.SpriteMaterial( 
                { map: texture, useScreenCoordinates: false, alignment: spriteAlignment, depthTest: !inFront } );
                
        this.sprite = new WebMol.Sprite( spriteMaterial );

        //this.sprite.scale.set(textWidth / textHeight * fontsize, fontsize, 1);
        this.sprite.scale.set(textWidth, textHeight, 1);
        this.sprite.position.set(position.x, position.y, position.z);
    }
    
};

// a webmol unified interace to gmol
WebMol.glmolViewer = (function() {
    // private class variables
    var numWorkers = 4; // number of threads for surface generation
    var maxVolume = 64000; // how much to break up surface calculations

    // private class helper functions

    // computes the bounding box around the provided atoms
    var getExtent = function(atomlist) {
        var xmin = ymin = zmin = 9999;
        var xmax = ymax = zmax = -9999;
        var xsum = ysum = zsum = cnt = 0;

        if (atomlist.length === 0)
            return [ [ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ] ];
        for ( var i = 0; i < atomlist.length; i++) {
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

    // The constructor
    function GLViewer(element, callback, defaultcolors) {

        // set variables
        var container = element;
        var id = container.id;

        var models = []; // atomistic molecular models
        var surfaces = [];

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
        renderer.domElement.style["zIndex"] = "0";
        container.append(renderer.domElement);
        renderer.setSize(WIDTH, HEIGHT);

        var camera = new WebMol.Camera(20, ASPECT, 1, 800);
        camera.position = new TV3(0, 0, CAMERA_Z);
        camera.lookAt(new TV3(0, 0, 0));

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
            // CHECK: Should I explicitly call scene.deallocateObject?
            scene = new WebMol.Scene();
            //scene = new WebMol.Scene();
            scene.fog = new WebMol.Fog(bgColor, 100, 200);

            modelGroup = new WebMol.Object3D();
            rotationGroup = new WebMol.Object3D();
            rotationGroup.useQuaternion = true;
            rotationGroup.quaternion = new WebMol.Quaternion(0, 0, 0, 1);
            rotationGroup.add(modelGroup);

            scene.add(rotationGroup);

            // setup lights
            var directionalLight = new WebMol.Light(0xFFFFFF);
            directionalLight.position = new TV3(0.2, 0.2, 1).normalize();
            directionalLight.intensity = 1.0;
            scene.add(directionalLight);
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
            if (x === undefined)
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

        glDOM.bind('mousemove touchmove', function(ev) { // touchmove
            ev.preventDefault();
            if (!scene)
                return;
            if (!isDragging)
                return;
            var mode = 0;
            var modeRadio = $('input[name=' + id + '_mouseMode]:checked');
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
            if (mode == 3 || (mouseButton == 3 && ev.ctrlKey)) { // Slab
                slabNear = cslabNear + dx * 100;
                slabFar = cslabFar + dy * 100;
            } else if (mode == 2 || mouseButton == 3 || ev.shiftKey) { // Zoom
                var scaleFactor = (CAMERA_Z - rotationGroup.position.z) * 0.85;
                if (scaleFactor < 80)
                    scaleFactor = 80;
                rotationGroup.position.z = cz - dy * scaleFactor;
            } else if (mode == 1 || mouseButton == 2 || ev.ctrlKey) { // Translate
                var scaleFactor = (CAMERA_Z - rotationGroup.position.z) * 0.85;
                if (scaleFactor < 20)
                    scaleFactor = 20;
                var translationByScreen = new TV3(dx * scaleFactor, -dy
                        * scaleFactor, 0);
                var q = rotationGroup.quaternion;
                var qinv = new WebMol.Quaternion(q.x, q.y, q.z, q.w).inverse()
                        .normalize();
                var translation = translationByScreen.applyQuaternion(qinv);
                modelGroup.position.x = currentModelPos.x + translation.x;
                modelGroup.position.y = currentModelPos.y + translation.y;
                modelGroup.position.z = currentModelPos.z + translation.z;
            } else if ((mode == 0 || mouseButton == 1) && r != 0) { // Rotate
                var rs = Math.sin(r * Math.PI) / r;
                dq.x = Math.cos(r * Math.PI);
                dq.y = 0;
                dq.z = rs * dx;
                dq.w = -rs * dy;
                rotationGroup.quaternion = new WebMol.Quaternion(1, 0, 0, 0);
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
            scene.fog.color = WebMol.CC.color(hex);
            show();
        };

        this.setWidth = function(w) {
            WIDTH = w;
            renderer.setSize(WIDTH, HEIGHT);
        };

        this.setHeight = function(h) {
            HEIGHT = h;
            renderer.setSize(WIDTH, HEIGHT);
        };

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
        };

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
        };

        // apply styles, models, etc in viewer
        this.render = function() {

            //spinner.show();
            var time1 = new Date();
            var view = this.getView();

            for ( var i = 0; i < models.length; i++) {
                if (models[i]) {
                    models[i].globj(modelGroup);
                }
            }

            for ( var i in surfaces) { // this is an array with possible holes
                if (surfaces.hasOwnProperty(i)) {
                    var geo = surfaces[i].geo;
                    // async surface generation can cause
                    // the geometry to be webgl initialized before it is fully
                    // formed; force various recalculations until full surface is
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
                        
                        initBuffers(geo, true);
                        // create new surface
                        var smesh = new WebMol.Mesh(geo, surfaces[i].mat);
                        surfaces[i].lastGL = smesh;
                        modelGroup.add(smesh);
                    } // else final surface already there
                }
            }
            this.setView(view);  //Calls show() => renderer render
            var time2 = new Date();
            spinner.hide();
            console.log("render time: " + (time2 - time1));
        };

        function getAtomsFromSel(sel) {
            var atoms = [];
            if (typeof (sel) === "undefined")
                sel = {};

            var ms = [];
            if (typeof sel.model === "undefined") {
                for ( var i = 0; i < models.length; i++) {
                    if (models[i])
                        ms.push(models[i]);
                }
            } else { // specific to some models
                var ms = sel.model;
                if (!$.isArray(ms))
                    ms = [ ms ];
            }

            for ( var i = 0; i < ms.length; i++) {
                atoms = atoms.concat(ms[i].selectedAtoms(sel));
            }
            return atoms;
        };
        
        function atomIsSelected(atom,sel) {
            if (typeof (sel) === "undefined")
                sel = {};

            var ms = [];
            if (typeof sel.model === "undefined") {
                for ( var i = 0; i < models.length; i++) {
                    if (models[i])
                        ms.push(models[i]);
                }
            } else { // specific to some models
                var ms = sel.model;
                if (!$.isArray(ms))
                    ms = [ ms ];
            }

            for ( var i = 0; i < ms.length; i++) {
                if(ms[i].atomIsSelected(atom, sel))
                    return true;
            }
            return false;
        };

        // return pdb output of selected atoms
        // currently only works if input was pdb
        this.pdbData = function(sel) {
            var atoms = getAtomsFromSel(sel);
            var ret = "";
            for ( var i = 0, n = atoms.length; i < n; ++i) {
                ret += atoms[i].pdbline + "\n";
            }
            return ret;
        };

        // zoom to atom selection
        this.zoomTo = function(sel) {
            var atoms = getAtomsFromSel(sel);
            var allatoms = getAtomsFromSel({});
            var tmp = getExtent(atoms);
            var alltmp = getExtent(allatoms);
            // use selection for center
            var center = new TV3(tmp[2][0], tmp[2][1], tmp[2][2]);
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
        
        // add a label to viewer
        this.addLabel = function(text, data) {
            var label = new WebMol.Label(text, data); 
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;
        };
        
        this.removeLabel = function(label) {
            label.sprite.material.map.dispose(); // Texture clean up
            label.sprite.material.dispose(); // SpriteMaterial clean up
            modelGroup.remove(label.sprite);                       
        };
        
        //Modify label style
        this.setLabelStyle = function(label, stylespec) {    
            this.removeLabel(label);

            label.stylespec = stylespec;
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;
            
        };
        
        //Change label text
        this.setLabelText = function(label, text) {
            this.removeLabel(label);
            
            label.text = text;
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;

        };

        // given molecular data and its format (pdb, sdf, xyz or mol2)
        // create a model and add it, returning the model identifier
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
            for (var i = 0; i < models.length; i++){
                var model = models[i];
                model.removegl(modelGroup);
                
            }
            models = [];
        };

        // create a new model out of sel,
        // if extract is true, removes sel form this model
        // updates bond indices appropriately
        this.createModelFrom = function(sel, extract) {
            var m = new WebMol.GLModel(models.length, defaultcolors);
            for ( var i = 0; i < models.length; i++) {
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
            for ( var i = 0; i < models.length; i++) {
                if (models[i]) {
                    models[i][func](sel, value1, value2);
                }
            }
        }

        // apply sel to all models and apply style
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

        var getAtomsWithin = function(atomlist, extent) {
            var ret = [];

            for ( var i = 0; i < atomlist.length; i++) {
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
                    var index = 0;
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
            var ret = [];
            // now compute atoms within expanded (this could be more efficient)
            var off = 6; // enough for water and 2*r, also depends on scale
            // factor
            for ( var i = 0, n = splits.length; i < n; i++) {
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
        var generateSurfaceMesh = function(atoms, VandF, mat) {
            var geo = new WebMol.Geometry();
            geo.geometryChunks = [];
            geo.geometryChunks.push( new geometryChunk() );
            
            var geoGroup = geo.geometryChunks[0];
            
            // reconstruct vertices and faces
            geo.vertices = [];
            var v = VandF.vertices;
            
            for ( var i = 0; i < v.length; i++) {
                
                geoGroup.vertexArr.push(v[i].x), geoGroup.vertexArr.push(v[i].y), geoGroup.vertexArr.push(v[i].z);

                geoGroup.colorArr.push(0.0), geoGroup.colorArr.push(0.0), geoGroup.colorArr.push(0.0);            
                geoGroup.normalArr.push(0.0), geoGroup.normalArr.push(0.0), geoGroup.normalArr.push(0.0);
                
                geoGroup.vertices++;
            }

            var faces = VandF.faces;

            // set colors for vertices
            var colors = [];
            for ( var i = 0; i < atoms.length; i++) {
                var atom = atoms[i];
                if (atom) {
                    if (typeof (atom.surfaceColor) != "undefined") {
                        colors[i] = WebMol.CC.color(atom.surfaceColor);
                    } else if (atom.color) // map from atom
                        colors[i] = WebMol.CC.color(atom.color);
                }
            }
            var verts = geoGroup.vertexArr;
            var vA, vB, vC, norm;
            //Setup colors, faces, and normals
            for ( var i = 0; i < faces.length; i++) {
                
                var a = faces[i].a, b = faces[i].b, c = faces[i].c;
                var A = v[a].atomid;
                var B = v[b].atomid;
                var C = v[c].atomid;
                
                var offsetA = a * 3, offsetB = b * 3, offsetC = c * 3;

                geoGroup.faceArr.push(faces[i].a), geoGroup.faceArr.push(faces[i].b), geoGroup.faceArr.push(faces[i].c);
                
                geoGroup.colorArr[offsetA] = colors[A].r, geoGroup.colorArr[offsetA+1] = colors[A].g,
                         geoGroup.colorArr[offsetA+2] = colors[A].b;
                geoGroup.colorArr[offsetB] = colors[B].r, geoGroup.colorArr[offsetB+1] = colors[B].g,
                         geoGroup.colorArr[offsetB+2] = colors[B].b;
                geoGroup.colorArr[offsetC] = colors[C].r, geoGroup.colorArr[offsetC+1] = colors[C].g,
                         geoGroup.colorArr[offsetC+2] = colors[C].b;
                 
                //setup Normals
                
                vA = new TV3(verts[offsetA], verts[offsetA+1], verts[offsetA+2]);
                vB = new TV3(verts[offsetB], verts[offsetB+1], verts[offsetB+2]);
                vC = new TV3(verts[offsetC], verts[offsetC+1], verts[offsetC+2]);
                
                vC.subVectors(vC, vB);
                vA.subVectors(vA, vB);
                vC.cross(vA);

                //face normal
                norm = vC;
                norm.normalize();
                
                geoGroup.normalArr[offsetA] += norm.x, geoGroup.normalArr[offsetB] += norm.x, geoGroup.normalArr[offsetC] += norm.x;
                geoGroup.normalArr[offsetA+1] += norm.y, geoGroup.normalArr[offsetB+1] += norm.y, geoGroup.normalArr[offsetC+1] += norm.y;
                geoGroup.normalArr[offsetA+2] += norm.z, geoGroup.normalArr[offsetB+2] += norm.z, geoGroup.normalArr[offsetC+2] += norm.z;
                
            }

            var mesh = new WebMol.Mesh(geo, mat);
            mesh.doubleSided = true;

            return mesh;
        };

        // do same thing as worker in main thread
        var generateMeshSyncHelper = function(type, expandedExtent,
                extendedAtoms, atomsToShow, atoms, vol) {
            var time = new Date();
            var ps = new ProteinSurface();
            ps.initparm(expandedExtent, (type == 1) ? false : true, vol);

            var time2 = new Date();
            console.log("initialize " + (time2 - time) + "ms");

            ps.fillvoxels(atoms, extendedAtoms);

            var time3 = new Date();
            console.log("fillvoxels " + (time3 - time2) + "  " + (time3 - time)
                    + "ms");

            ps.buildboundary();

            if (type == 4 || type == 2)
                ps.fastdistancemap();
            if (type == 2) {
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
            ps.laplaciansmooth(1);
            return ps.getFacesAndVertices(atomsToShow);
        };

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
            if ( style.opacity !== undefined) {
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

            for ( var i = 0, n = atomlist.length; i < n; i++) {
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

            var time = new Date();

            var mat = getMatWithStyle(style);

            var extent = getExtent(atomsToShow);

            if (style.map && style.map.prop) {
                // map color space using already set atom properties
                var prop = style.map.prop;
                var scheme = style.map.scheme || new WebMol.RWB();
                var range = scheme.range();
                if (!range) {
                    range = getPropertyRange(atomsToShow, prop);
                }

                for ( var i = 0, n = atomsToShow.length; i < n; i++) {
                    var atom = atomsToShow[i];
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
                geo : new WebMol.Geometry(),
                mat : mat,
                done : false,
                finished : false
            // also webgl initialized
            };
            var surfid = surfaces.length;
            surfaces[surfid] = surfobj;
            var reducedAtoms = [];
            // to reduce amount data transfered, just pass x,y,z,serial and elem
            for ( var i = 0, n = atomlist.length; i < n; i++) {
                var atom = atomlist[i];
                reducedAtoms[i] = {
                    x : atom.x,
                    y : atom.y,
                    z : atom.z,
                    serial : i,
                    elem : atom.elem
                };
            }

            var sync = false;
            var view = this; //export render function to worker
            if (sync) { // don't use worker, still break up for memory purposes

                for ( var i = 0; i < extents.length; i++) {
                    var VandF = generateMeshSyncHelper(type, extents[i].extent,
                            extents[i].atoms, extents[i].toshow, reducedAtoms,
                            totalVol);
                    var mesh = generateSurfaceMesh(atomlist, VandF, mat);
                    mergeGeos(surfobj.geo, mesh);
                    view.render();
                }
            } else { // use worker
                
                var workers = [];
                if (type < 0)
                    type = 0; // negative reserved for atom data
                for ( var i = 0; i < numWorkers; i++) {
                    var w = new Worker('webmol/SurfaceWorker.js');
                    workers.push(w);
                    w.postMessage({
                        type : -1,
                        atoms : reducedAtoms,
                        volume : totalVol
                    });
                }
                var cnt = 0;
                for ( var i = 0; i < extents.length; i++) {
                    var worker = workers[i % workers.length];
                    worker.onmessage = function(event) {
                        var VandF = event.data;
                        var mesh = generateSurfaceMesh(atomlist, VandF, mat);
                        mergeGeos(surfobj.geo, mesh);
                        view.render();
                        console.log("async mesh generation "
                                + (+new Date() - time) + "ms");
                        cnt++;
                        if (cnt == extents.length)
                            surfobj.done = true;
                    };

                    worker.onerror = function(event) {
                        console.log(event.message + " (" + event.filename + ":"
                                + event.lineno + ")");
                    };

                    worker.postMessage({
                        type : type,
                        expandedExtent : extents[i].extent,
                        extendedAtoms : extents[i].atoms,
                        atomsToShow : extents[i].toshow,
                    });
                }
            }

            //NOTE: This is misleading if 'async' mesh generation - returns immediately
            console.log("full mesh generation " + (+new Date() - time) + "ms");

            return surfid;
        };

        // set the material to something else, must render change
        this.setSurfaceMaterialStyle = function(surf, style) {
            if (surfaces[surf]) {
                surfaces[surf].mat = getMatWithStyle(style);
                surfaces[surf].finished = false; //trigger redraw
            }
        };

        // given the id returned by surfid, remove surface
        this.removeSurface = function(surf) {
            if (surfaces[surf] && surfaces[surf].lastGL) {
                if (surfaces[surf].geo !== undefined) surfaces[surf].geo.dispose();             
                if (surfaces[surf].mat !== undefined) surfaces[surf].mat.dispose();
                modelGroup.remove(surfaces[surf].lastGL); // remove from scene
            }
            delete surfaces[surf];
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

            surfaces = [];
            //models = [];
            this.removeAllModels();
            show();
        };

        // props is a list of objects that select certain atoms and enumerate
        // properties for those atoms
        this.mapAtomProperties = function(props) {
            var atoms = getAtomsFromSel({});
            for(var a = 0, numa = atoms.length; a < numa; a++) {
                var atom = atoms[a];
                for ( var i = 0, n = props.length; i < n; i++) {
                    var prop = props[i];
                    if (prop.props) {
                        for ( var p in prop.props) {
                            if (prop.props.hasOwnProperty(p)) {
                                // check the atom
                                if(atomIsSelected(atom, prop)) {
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
        
        this.getModelGroup = function() {
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
