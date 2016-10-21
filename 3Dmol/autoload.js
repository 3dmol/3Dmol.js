//auto-initialization
//Create embedded viewer from HTML attributes if true
$3Dmol.autoload=function(callback,viewer){
    if ($(".viewer_3Dmoljs")[0] !== undefined)
        $3Dmol.autoinit = true;
        
    if ($3Dmol.autoinit) {
        viewer =(viewer!= undefined) ?
            viewer :null;
        var callback = (callback!=undefined) ? 
                callback : null;
        $3Dmol.viewers = {};
        var nviewers = 0;
        $(".viewer_3Dmoljs").each( function() {
            var viewerdiv = $(this);
            var datauri = null;
            if(viewerdiv.css('position') == 'static') {
                //slight hack - canvas needs this element to be positioned
                viewerdiv.css('position','relative');
            }
            var type = null;
            if (viewerdiv.data("pdb")) {
                datauri = "http://files.rcsb.org/view/" + viewerdiv.data("pdb") + ".pdb";
                type = "pdb";
            } else if(viewerdiv.data("cid")) {
                //this doesn't actually work since pubchem does have CORS enabled
                type = "sdf";
                datauri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + viewerdiv.data("cid") + 
                "/SDF?record_type=3d";
            }
            else if (viewerdiv.data("href"))
                datauri = viewerdiv.data("href");
                
            var options = {}
            if(viewerdiv.data("options"))
                options = $3Dmol.specStringToObject(viewerdiv.data("options"));
                
            var bgcolor = $3Dmol.CC.color(viewerdiv.data("backgroundcolor"));
            var style = {line:{}};
            if(viewerdiv.data("style")) style = $3Dmol.specStringToObject(viewerdiv.data("style"));
            var select = {};
            if(viewerdiv.data("select")) select = $3Dmol.specStringToObject(viewerdiv.data("select"));
            var selectstylelist = [];
            var surfaces = [];
            var labels = [];
            var d = viewerdiv.data();
            
            //let users specify individual but matching select/style tags, eg.
            //data-select1 data-style1
            var stylere = /style(.+)/;
            var surfre = /surface(.*)/;
            var reslabre = /labelres(.*)/;
            var keys = [];
            for(var dataname in d) {
                if(d.hasOwnProperty(dataname)) {
                    keys.push(dataname);
                }
            }
            keys.sort();
            for(var i = 0; i < keys.length; i++) {
                var dataname = keys[i];
                var m = stylere.exec(dataname);
                if(m) {
                    var selname = "select"+m[1];
                    var newsel = $3Dmol.specStringToObject(d[selname]);
                    var styleobj = $3Dmol.specStringToObject(d[dataname]);
                    selectstylelist.push([newsel,styleobj]);
                }         
                m = surfre.exec(dataname);
                if(m) {
                    var selname = "select"+m[1];
                    var newsel = $3Dmol.specStringToObject(d[selname]);
                    var styleobj = $3Dmol.specStringToObject(d[dataname]);
                    surfaces.push([newsel,styleobj]);
                }
                m = reslabre.exec(dataname);
                if(m) {
                    var selname = "select"+m[1];
                    var newsel = $3Dmol.specStringToObject(d[selname]);
                    var styleobj = $3Dmol.specStringToObject(d[dataname]);
                    labels.push([newsel,styleobj]);
                }
            }
            
            //apply all the selections/styles parsed out above to the passed viewer
            var applyStyles = function(glviewer) {
                glviewer.setStyle(select,style);
                for(var i = 0; i < selectstylelist.length; i++) {
                    var sel = selectstylelist[i][0] || {};
                    var sty = selectstylelist[i][1] || {"line":{}}
                    glviewer.setStyle(sel, sty);
                }
                for(var i = 0; i < surfaces.length; i++) {
                    var sel = surfaces[i][0] || {};
                    var sty = surfaces[i][1] || {}
                    glviewer.addSurface($3Dmol.SurfaceType.VDW, sty, sel, sel);
                }
                for(var i = 0; i < labels.length; i++) {
                    var sel = labels[i][0] || {};
                    var sty = labels[i][1] || {}
                    glviewer.addResLabels(sel, sty);
                }               
                
                glviewer.zoomTo();
                glviewer.render();             
            }
            
            
            var glviewer = viewer;
            try {
                if(glviewer==null)
                    glviewer = $3Dmol.viewers[this.id || nviewers++] = $3Dmol.createViewer(viewerdiv, {defaultcolors: $3Dmol.rasmolElementColors});
                glviewer.setBackgroundColor(bgcolor);                            
            } catch ( error ) {
                //for autoload, provide a useful error message
                window.location = "http://get.webgl.org";                    
            }           
            
            if (datauri) {  
                
                type = viewerdiv.data("type") || viewerdiv.data("datatype") || type;
                if(!type) {
                    type = datauri.substr(datauri.lastIndexOf('.')+1);
                }
                                
                $.get(datauri, function(ret) {
                    glviewer.addModel(ret, type, options);
                    applyStyles(glviewer);       
                    if (callback) 
                        callback(glviewer);
                }, 'text');         
            }           
            else {
                
                if (viewerdiv.data("element")) {
                    var moldata = $("#" + viewerdiv.data("element")).val() || "";
                    var type = viewerdiv.data("type") || viewerdiv.data("datatype");

                    if (!type){

                        console.log("Warning: No type specified for embedded viewer with moldata from " + viewerdiv.data("element") +
                                    "\n assuming type 'pdb'")

                        type = 'pdb';
                    }

                    glviewer.addModel(moldata, type, options);        
                }

                applyStyles(glviewer);                
                if (callback) 
                    callback(glviewer);
            }
            
        });              
    }}
$(document).ready(function() {
    $3Dmol.autoload();
    
});
    
 