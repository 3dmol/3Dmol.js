//auto-initialization
//Create embedded viewer from HTML attributes if true

$(document).ready(function() {

    if ($(".webmoljs_viewer")[0] !== undefined)
        WebMol.autoinit = true;
        
    if (WebMol.autoinit) { 
        WebMol.viewers = {};
        var nviewers = 0;
        $(".webmoljs_viewer").each( function() {
            var viewerdiv = $(this);
            var datauri = null;
            
            var callback = (typeof(window[viewerdiv.data("callback")]) === 'function') ? 
                    window[viewerdiv.data("callback")] : null;
            
            if (viewerdiv.data("pdb"))
                datauri = "http://www.pdb.org/pdb/files/" + viewerdiv.data("pdb") + ".pdb";
            else if (viewerdiv.data("href"))
                datauri = viewerdiv.data("href");
                
            var bgcolor = Number(viewerdiv.data("backgroundcolor")) || 0x000000;
            var style = viewerdiv.data("style") || {line:{}};
            
            var glviewer = WebMol.viewers[this.id || nviewers++] = WebMol.createViewer(viewerdiv, {defaultcolors: WebMol.rasmolElementColors, callback: function(viewer) {            
                viewer.setBackgroundColor(bgcolor);            
            }});
            
            if (datauri) {  
                
                var type = viewerdiv.data("datatype") || "pdb";
                 
                $.get(datauri, function(ret) {
                    glviewer.addModel(ret, type);
                    glviewer.setStyle({}, style);
                    
                    // Allowing us to fire callback after viewer has added model
                    if (callback) 
                        callback(glviewer);                    
                    
                    glviewer.zoomTo();
                    glviewer.render();          
                    
                }, 'text');
           
            }
            
            else {
                console.log("Executing callback on viewer " + this.id + ". No input data specified");
                if (callback) 
                    callback(glviewer);
                
                glviewer.zoomTo();
                glviewer.render();
            }
            
        });              
    }
});
    