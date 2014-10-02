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
            
            if(typeof($(viewerdiv).css('position')) === "undefined") {
            	$(viewerdiv).css('position','relative'); //to overlap spinner need positioned container
            }
            var callback = (typeof(window[viewerdiv.data("callback")]) === 'function') ? 
                    window[viewerdiv.data("callback")] : null;
            
            if (viewerdiv.data("pdb"))
                datauri = "http://www.pdb.org/pdb/files/" + viewerdiv.data("pdb") + ".pdb";
            else if (viewerdiv.data("href"))
                datauri = viewerdiv.data("href");
                
            var bgcolor = Number(viewerdiv.data("backgroundcolor")) || 0x000000;
            var style = viewerdiv.data("style") || {line:{}};
            var select = viewerdiv.data("select") || {};
            
            var glviewer = $3Dmol.viewers[this.id || nviewers++] = $3Dmol.createViewer(viewerdiv, {defaultcolors: $3Dmol.rasmolElementColors, callback: function(viewer) {            
                viewer.setBackgroundColor(bgcolor);            
            }});
            
            
            if (datauri) {  
                
                var type = viewerdiv.data("datatype") || "pdb";
                 
                $.get(datauri, function(ret) {
                    glviewer.addModel(ret, type);
                    glviewer.setStyle(select, style);
                    
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
                }


                if (callback) 
                    callback(glviewer);
                
                glviewer.zoomTo();
                glviewer.render();
            }
            
        });              
    }
});
    
