//auto-initialization
//Create embedded viewer from HTML attributes if true
//viewer and callback are used by the testing harness
$3Dmol.autoload=function(viewer,callback){
    var i, dataname, type;
    if ($(".viewer_3Dmoljs")[0] !== undefined)
        $3Dmol.autoinit = true;

    if ($3Dmol.autoinit) {
        viewer =(viewer!= undefined) ? viewer :null;
        
        $3Dmol.viewers = {};
        var nviewers = 0;
        $(".viewer_3Dmoljs").each( function() {
            var viewerdiv = $(this);
            var datauri = [];
            var datatypes = [];
            var uri = '';
            var showUI = false;
            if(viewerdiv.css('position') == 'static') {
                //slight hack - canvas needs this element to be positioned
                viewerdiv.css('position','relative');
            }

            if(viewerdiv.data('ui'))
                showUI = true;

            type = null;
            if (viewerdiv.data("pdb")) {
                datauri.push("https://files.rcsb.org/view/" + viewerdiv.data("pdb") + ".pdb");
                datatypes.push("pdb");
            } else if(viewerdiv.data("cid")) {
                //this doesn't actually work since pubchem does have CORS enabled
                datatypes.push("sdf");
                datauri.push("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + viewerdiv.data("cid") + 
                "/SDF?record_type=3d");
            }
            else if (viewerdiv.data("href") || viewerdiv.data("url")){
                uri = viewerdiv.data("href");
                datauri.push(uri);
                type = uri.substr(uri.lastIndexOf('.')+1);                
                datatypes.push(type);

                var molName = uri.substring(uri.lastIndexOf('/') + 1, uri.lastIndexOf('.'));
                if(molName == '/') 
                    molName = uri.substring(uri.lastIndexOf('/') + 1);

                viewerdiv.data(datatypes[datatypes.length -1 ], molName);
            }
            
            var divdata=viewerdiv.data();
            for(i in divdata){
                if((i.substring(0,3) ==="pdb" && (i !== "pdb"))){
                    datauri.push("https://files.rcsb.org/view/" +divdata[i]+".pdb");
                    datatypes.push('pdb');

                }else if(i.substring(0,4) ==="href" && (i!=="href")){
                    uri = divdata[i];
                    datauri.push(uri);
                    datatypes.push(uri.substr(uri.lastIndexOf('.')+1));
                }else if(i.substring(0,3)==="cid" && (i!=="cid")){
                    datauri.push("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + divdata[i] +  "/SDF?record_type=3d");
                    datatypes.push('sdf');
                }
            }
            var options = {};
            if(viewerdiv.data("options"))
                options = $3Dmol.specStringToObject(viewerdiv.data("options"));
                
            //note that data tags must be lowercase
            var bgcolor = $3Dmol.CC.color(viewerdiv.data("backgroundcolor"));
            var bgalpha = viewerdiv.data("backgroundalpha");
            bgalpha = bgalpha == undefined ? 1.0 : parseFloat(bgalpha);
            var style = {line:{}};
            if(viewerdiv.data("style")) style = $3Dmol.specStringToObject(viewerdiv.data("style"));
            var select = {};
            if(viewerdiv.data("select")) select = $3Dmol.specStringToObject(viewerdiv.data("select"));
            var selectstylelist = [];
            var surfaces = [];
            var labels = [];
            var zoomto = {};
            var spin = null;
            var d = viewerdiv.data();
            
            //let users specify individual but matching select/style/surface tags, eg.
            //data-select1 data-style1
            var stylere = /style(.+)/;
            var surfre = /surface(.*)/;
            var reslabre = /labelres(.*)/;
            var keys = [];
            for(dataname in d) {
                if(Object.prototype.hasOwnProperty.call(d,dataname)) {
                    keys.push(dataname);
                }
            }
            keys.sort();
            for(i = 0; i < keys.length; i++) {
                dataname = keys[i];
                var m = stylere.exec(dataname);
                var selname, newsel, styleobj;
                if(m) {
                    selname = "select"+m[1];
                    newsel = $3Dmol.specStringToObject(d[selname]);
                    styleobj = $3Dmol.specStringToObject(d[dataname]);
                    selectstylelist.push([newsel,styleobj]);
                }         
                m = surfre.exec(dataname);
                if(m) {
                    selname = "select"+m[1];
                    newsel = $3Dmol.specStringToObject(d[selname]);
                    styleobj = $3Dmol.specStringToObject(d[dataname]);
                    surfaces.push([newsel,styleobj]);
                }
                m = reslabre.exec(dataname);
                if(m) {
                    selname = "select"+m[1];
                    newsel = $3Dmol.specStringToObject(d[selname]);
                    styleobj = $3Dmol.specStringToObject(d[dataname]);
                    labels.push([newsel,styleobj]);
                }
                if(dataname == "zoomto") {
                    zoomto = $3Dmol.specStringToObject(d[dataname]);
                }
                if(dataname == "spin") {
                    spin = $3Dmol.specStringToObject(d[dataname]);
                }
            }
            
            //apply all the selections/styles parsed out above to the passed viewer
            var applyStyles = function(glviewer) {
                glviewer.setStyle(select,style);

                if(showUI){
                    glviewer.ui.loadSelectionStyle(select, style);
                }

                for(i = 0; i < selectstylelist.length; i++) {
                    let sel = selectstylelist[i][0] || {};
                    let sty = selectstylelist[i][1] || {"line":{}};
                    glviewer.setStyle(sel, sty);
                    if(showUI){
                        glviewer.ui.loadSelectionStyle(sel, sty);
                    }
                }
                for(i = 0; i < surfaces.length; i++) {
                    let sel = surfaces[i][0] || {};
                    let sty = surfaces[i][1] || {};
                    let viewer = glviewer;

                    if(showUI){
                        //seemingly unnecessary capturing of values due to jshint
                        //not understanding let?
                        let doload = function($3D, viewer, sel, sty) {
                            viewer.addSurface($3D.SurfaceType.VDW, sty, sel, sel).then((surfid)=>{
                                viewer.ui.loadSurface("VDW", sel, sty, surfid);
                            });
                        };
                        doload($3Dmol, viewer, sel, sty);                 
                    }
                    else {
                        glviewer.addSurface($3Dmol.SurfaceType.VDW, sty, sel, sel);
                    }
                    
                }
                for(i = 0; i < labels.length; i++) {
                    let sel = labels[i][0] || {};
                    let sty = labels[i][1] || {};
                    glviewer.addResLabels(sel, sty);
                }               
                
                glviewer.zoomTo(zoomto);
                if(spin) {
                    glviewer.spin(spin.axis,spin.speed);
                }
                glviewer.render();             
            };
            
            var glviewer = viewer;
            try {
                if(glviewer==null) {
                    var config = viewerdiv.data('config') || {};
                    config.defaultcolors = config.defaultcolors || $3Dmol.rasmolElementColors;
                    if(config.backgroundColor === undefined) config.backgroundColor = bgcolor;
                    if(config.backgroundAlpha === undefined) config.backgroundAlpha = bgalpha;                     
                    config.ui = showUI;
                    glviewer = $3Dmol.viewers[this.id || nviewers++] = $3Dmol.createViewer(viewerdiv, config);

                } else {
                    glviewer.setBackgroundColor(bgcolor, bgalpha);
                    if(showUI) glviewer.ui.initiateUI();
                } 
            } catch ( error ) {
                console.log(error);
                //for autoload, provide a useful error message
                viewerdiv.text("WebGL appears to be disabled.");
            }           
            
            if (datauri.length!=0) {
                //load multiple data elements in serial
                i = 0;
                var process = function(moldata) {
                    //add moldata to viewer and load next model
                    uri = datauri[i]; //this is where the moldata came from
                    var type = viewerdiv.data("type") || viewerdiv.data("datatype") || datatypes[i]; 
                    glviewer.addModel(moldata, type, options);
                    if(showUI){
                        var modelName = viewerdiv.data(datatypes[i]);

                        glviewer.ui.setModelTitle(modelName);
                    }
                    i += 1;
                    if(i < datauri.length) {
                        $.get(datauri[i], process, 'text');
                    }
                    else {
                        // or finalize if this is the last model
                        applyStyles(glviewer);
                        if(viewerdiv.data("callback")) {
                            //evaluate javascript in the string, if it resolves to a function,
                            //call it with the viewer
                            /*jshint -W061 */ var runres = eval(viewerdiv.data("callback"));  
                            if(typeof(runres) == 'function') {
                                runres(glviewer);
                            }
                        }
                        if(callback) callback(glviewer);
                    }
                };
                $.get(datauri[0], process, 'text');
         
            }           
            else {
                
                if (viewerdiv.data("element")) {
                    var moldata = $("#" + viewerdiv.data("element")).val() || "";
                    type = viewerdiv.data("type") || viewerdiv.data("datatype");

                    if (!type){

                        console.log("Warning: No type specified for embedded viewer with moldata from " + viewerdiv.data("element") +
                                    "\n assuming type 'pdb'");
                        type = 'pdb';
                    }

                    glviewer.addModel(moldata, type, options);        
                }

                applyStyles(glviewer);
                if(viewerdiv.data("callback")) {
                    //evaluate javascript in the string, if it resolves to a function,
                    //call it with the viewer
                    /*jshint -W061 */ var runres = eval(viewerdiv.data("callback"));
                    if(typeof(runres) == 'function') {
                        runres(glviewer);
                    }
                }                
                if (callback) 
                    callback(glviewer);
            }            
        });                      
    }};
$(document).ready(function() {
    $3Dmol.autoload();
    
});
    
 
