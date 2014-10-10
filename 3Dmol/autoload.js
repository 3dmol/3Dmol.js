//auto-initialization
//Create embedded viewer from HTML attributes if true

$(document).ready(function() {

    if ($(".viewer_3Dmoljs")[0] !== undefined)
        $3Dmol.autoinit = true;
        
    if ($3Dmol.autoinit) { 
        $3Dmol.viewers = {};
        var nviewers = 0;
        $(".viewer_3Dmoljs").each( function() {
            var viewerdiv = $(this);
            var datauri = null;
            
        
            var callback = (typeof(window[viewerdiv.data("callback")]) === 'function') ? 
                    window[viewerdiv.data("callback")] : null;
            
            if (viewerdiv.data("pdb"))
                datauri = "http://www.pdb.org/pdb/files/" + viewerdiv.data("pdb") + ".pdb";
            else if (viewerdiv.data("href"))
                datauri = viewerdiv.data("href");
                
            var bgcolor = Number(viewerdiv.data("backgroundcolor")) || 0x000000;
            var style = {line:{}};
            if(viewerdiv.data("style")) style = $3Dmol.specStringToObject(viewerdiv.data("style"));
            var select = {};
            if(viewerdiv.data("select")) select = $3Dmol.specStringToObject(viewerdiv.data("select"));
            var selectstylelist = [];
            var surfaces = []
            var d = viewerdiv.data();
            
            //let users specify individual but matching select/style tags, eg.
            //data-select1 data-style1
            var stylere = /style(.+)/;
            var surfre = /surface(.*)/;
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
            }
            
            var glviewer = $3Dmol.viewers[this.id || nviewers++] = $3Dmol.createViewer(viewerdiv, {defaultcolors: $3Dmol.rasmolElementColors, callback: function(viewer) {            
                viewer.setBackgroundColor(bgcolor);            
            }});
            
            
            if (datauri) {  
                
                var type = viewerdiv.data("datatype") || "pdb";
                 
                $.get(datauri, function(ret) {
                    glviewer.addModel(ret, type);
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
                    // Allowing us to fire callback after viewer has added model
                    if (callback) 
                        callback(glviewer);                    
                    
                    glviewer.zoomTo();
                    glviewer.render();          
                    
                }, 'text');
           
            }
            
            else {
                
                if (viewerdiv.data("element")) {
                    var moldata = $("#" + viewerdiv.data("element")).val() || "";
                    var type = viewerdiv.data("datatype");

                    if (!type){

                        console.log("Warning: No type specified for embedded viewer with moldata from " + viewerdiv.data("element") +
                                    "\n assuming type 'pdb'")

                        type = 'pdb';
                    }

                    glviewer.addModel(moldata, type);
                	glviewer.setStyle(select, style);
                    for(var i = 0; i < selectstylelist.length; i++) {
                    	var sel = selectstylelist[i][0] || {};
                    	var sty = selectstylelist[i][1] || {"line":{}}
                    	glviewer.setStyle(sel, sty);
                    }                
                }


                if (callback) 
                    callback(glviewer);
                
                glviewer.zoomTo();
                glviewer.render();
            }
            
        });              
    }
});
    
