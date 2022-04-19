// auto-initialization
import $ from "jquery"
import { elementColors, CC } from "../colors";
import SurfaceType from "../enum/SurfaceType";
import createViewer from "../util/createViewer";
import specStringToObject from "../util/specStringToObject";
import Viewers from "../singletons/Viewers";

// eslint-disable-next-line import/no-mutable-exports
export let autoinit = false;
// eslint-disable-next-line import/no-mutable-exports
export let processing_autoinit = false;

// Create embedded viewer from HTML attributes if true
// viewer and callback are used by the testing harness
const autoload=function(viewer,callback){
    let i; let dataname; let type;
    if ($(".viewer_3Dmoljs")[0] !== undefined)
        autoinit = true;

    if (autoinit) {
        processing_autoinit = true;
        viewer =(viewer!== undefined) ? viewer :null;
        
        let nviewers = 0;
        $(".viewer_3Dmoljs").each(function() {
            // @ts-ignore
            const viewerdiv = $(this);
            const datauri = [];
            const datatypes = [];
            let uri = '';
            let showUI = false;
            if(viewerdiv.css('position') === 'static') {
                // slight hack - canvas needs this element to be positioned
                viewerdiv.css('position','relative');
            }

            if(viewerdiv.data('ui'))
                showUI = true;

            type = null;
            if (viewerdiv.data("pdb")) {
                datauri.push(`https://files.rcsb.org/view/${  viewerdiv.data("pdb")  }.pdb`);
                datatypes.push("pdb");
            } else if(viewerdiv.data("cid")) {
                // this doesn't actually work since pubchem does have CORS enabled
                datatypes.push("sdf");
                datauri.push(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${  viewerdiv.data("cid")  
                }/SDF?record_type=3d`);
            }
            else if (viewerdiv.data("href") || viewerdiv.data("url")){
                uri = viewerdiv.data("href");
                datauri.push(uri);
                type = uri.substr(uri.lastIndexOf('.')+1);                
                datatypes.push(type);

                let molName = uri.substring(uri.lastIndexOf('/') + 1, uri.lastIndexOf('.'));
                if(molName === '/') 
                    molName = uri.substring(uri.lastIndexOf('/') + 1);

                viewerdiv.data(datatypes[datatypes.length -1 ], molName);
            }
            
            const divdata=viewerdiv.data();
            for(i in divdata){
                if((i.substring(0,3) ==="pdb" && (i !== "pdb"))){
                    datauri.push(`https://files.rcsb.org/view/${ divdata[i]}.pdb`);
                    datatypes.push('pdb');

                }else if(i.substring(0,4) ==="href" && (i!=="href")){
                    uri = divdata[i];
                    datauri.push(uri);
                    datatypes.push(uri.substr(uri.lastIndexOf('.')+1));
                }else if(i.substring(0,3)==="cid" && (i!=="cid")){
                    datauri.push(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${  divdata[i]   }/SDF?record_type=3d`);
                    datatypes.push('sdf');
                }
            }
            let options = {};
            if(viewerdiv.data("options"))
                options = specStringToObject(viewerdiv.data("options"));
                
            // note that data tags must be lowercase
            const bgcolor = CC.color(viewerdiv.data("backgroundcolor"));
            let bgalpha = viewerdiv.data("backgroundalpha");
            bgalpha = bgalpha === undefined ? 1.0 : parseFloat(bgalpha);
            let style = {line:{}};
            if(viewerdiv.data("style")) style = specStringToObject(viewerdiv.data("style"));
            let select = {};
            if(viewerdiv.data("select")) select = specStringToObject(viewerdiv.data("select"));
            const selectstylelist = [];
            const surfaces = [];
            const labels = [];
            let zoomto = {};
            let spin = null;
            const d = viewerdiv.data();
            
            // let users specify individual but matching select/style/surface tags, eg.
            // data-select1 data-style1
            const stylere = /style(.+)/;
            const surfre = /surface(.*)/;
            const reslabre = /labelres(.*)/;
            const keys = [];
            for(dataname in d) {
                if(Object.prototype.hasOwnProperty.call(d,dataname)) {
                    keys.push(dataname);
                }
            }
            keys.sort();
            for(i = 0; i < keys.length; i++) {
                dataname = keys[i];
                let m = stylere.exec(dataname);
                let selname; let newsel; let styleobj;
                if(m) {
                    selname = `select${m[1]}`;
                    newsel = specStringToObject(d[selname]);
                    styleobj = specStringToObject(d[dataname]);
                    selectstylelist.push([newsel,styleobj]);
                }         
                m = surfre.exec(dataname);
                if(m) {
                    selname = `select${m[1]}`;
                    newsel = specStringToObject(d[selname]);
                    styleobj = specStringToObject(d[dataname]);
                    surfaces.push([newsel,styleobj]);
                }
                m = reslabre.exec(dataname);
                if(m) {
                    selname = `select${m[1]}`;
                    newsel = specStringToObject(d[selname]);
                    styleobj = specStringToObject(d[dataname]);
                    labels.push([newsel,styleobj]);
                }
                if(dataname === "zoomto") {
                    zoomto = specStringToObject(d[dataname]);
                }
                if(dataname === "spin") {
                    spin = specStringToObject(d[dataname]);
                }
            }
            
            // apply all the selections/styles parsed out above to the passed viewer
            const applyStyles = function(glviewer) {
                glviewer.setStyle(select,style);

                if(showUI){
                    glviewer.ui.loadSelectionStyle(select, style);
                }

                for(i = 0; i < selectstylelist.length; i++) {
                    const sel = selectstylelist[i][0] || {};
                    const sty = selectstylelist[i][1] || {"line":{}};
                    glviewer.setStyle(sel, sty);
                    if(showUI){
                        glviewer.ui.loadSelectionStyle(sel, sty);
                    }
                }
                for(i = 0; i < surfaces.length; i++) {
                    const sel = surfaces[i][0] || {};
                    const sty = surfaces[i][1] || {};
                    const viewer = glviewer;

                    if(showUI){
                        // seemingly unnecessary capturing of values due to jshint
                        // not understanding let?
                        const doload = function(viewer, sel, sty) {
                            viewer.addSurface(SurfaceType.VDW, sty, sel, sel).then((surfid)=>{
                                viewer.ui.loadSurface("VDW", sel, sty, surfid);
                            });
                        };
                        doload(viewer, sel, sty);                 
                    }
                    else {
                        glviewer.addSurface(SurfaceType.VDW, sty, sel, sel);
                    }
                    
                }
                for(i = 0; i < labels.length; i++) {
                    const sel = labels[i][0] || {};
                    const sty = labels[i][1] || {};
                    glviewer.addResLabels(sel, sty);
                }               
                
                glviewer.render();             
                glviewer.zoomTo(zoomto);
                
                if(spin) {
                    glviewer.spin(spin.axis,spin.speed);
                }
            };
            
            let glviewer = viewer;
            try {
                const config = viewerdiv.data('config') || {};
                config.defaultcolors = config.defaultcolors || elementColors.rasmol;
                if(config.backgroundColor === undefined) config.backgroundColor = bgcolor;
                if(config.backgroundAlpha === undefined) config.backgroundAlpha = bgalpha;
                config.ui = showUI;		    
                if(glviewer==null) {
                    // @ts-ignore
                    glviewer = Viewers[this.id || nviewers++] = createViewer(viewerdiv, config);
                } else {
                    glviewer.setBackgroundColor(bgcolor, bgalpha);
		            glviewer.setConfig(config);
                    if(showUI) glviewer.ui.initiateUI();
                } 
            } catch ( error ) {
                console.log(error);
                // for autoload, provide a useful error message
                viewerdiv.text("WebGL appears to be disabled.");
            }           
            
            if (datauri.length!==0) {
                // load multiple data elements in serial
                i = 0;
                const process = function(moldata) {
                    // add moldata to viewer and load next model
                    uri = datauri[i]; // this is where the moldata came from
                    const type = viewerdiv.data("type") || viewerdiv.data("datatype") || datatypes[i]; 
                    glviewer.addModel(moldata, type, options);
                    if(showUI){
                        const modelName = viewerdiv.data(datatypes[i]);
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
                            // evaluate javascript in the string, if it resolves to a function,
                            // call it with the viewer
                            /* jshint -W061 */ const runres = eval(viewerdiv.data("callback"));  
                            if(typeof(runres) === 'function') {
                                runres(glviewer);
                            }
                        }
                        const processing_autoinit = false;
                        if(callback) callback(glviewer);
                    }
                };
                $.get(datauri[0], process, 'text');
         
            }           
            else {
                
                if (viewerdiv.data("element")) {
                    const moldata = $(`#${  viewerdiv.data("element")}`).val() || "";
                    type = viewerdiv.data("type") || viewerdiv.data("datatype");

                    if (!type){

                        console.log(`Warning: No type specified for embedded viewer with moldata from ${  viewerdiv.data("element") 
                                    }\n assuming type 'pdb'`);
                        type = 'pdb';
                    }

                    glviewer.addModel(moldata, type, options);        
                }

                applyStyles(glviewer);
                if(viewerdiv.data("callback")) {
                    // evaluate javascript in the string, if it resolves to a function,
                    // call it with the viewer
                    /* jshint -W061 */ const runres = eval(viewerdiv.data("callback"));
                    if(typeof(runres) === 'function') {
                        runres(glviewer);
                    }
                }
                const processing_autoinit = false;                
                if (callback) 
                    callback(glviewer);
            }            
        });                      
    }};
    
document.addEventListener('DOMContentLoaded', () => {
    autoload();    
});
    
 
