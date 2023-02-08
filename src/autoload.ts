//auto-initialization

import { createViewer } from "./GLViewer";
import { SurfaceType } from "./ProteinSurface4";
import { get, specStringToObject } from "./utilities";
import { CC } from "./colors";

export var autoinit = false;
export var processing_autoinit = false;

/**
 * Contains a dictionary of embedded viewers created from HTML elements
 * with a the viewer_3Dmoljs css class indexed by their id (or numerically
 * if they do not have an id).
*/
export var viewers: any = {};

//Create embedded viewer from HTML attributes if true
//viewer and callback are used by the testing harness
export function autoload(viewer?, callback?) {
    var i, dataname, type;
    if (document.querySelector(".viewer_3Dmoljs") != null)
        autoinit = true;

    if (autoinit) {
        processing_autoinit = true;
        viewer = (viewer != undefined) ? viewer : null;
        var nviewers = 0;
        document.querySelectorAll<HTMLInputElement>(".viewer_3Dmoljs").forEach(viewerdiv => {
            var datauri = [];
            var datatypes = [];
            var uri = '';

            if (viewerdiv.style.position == 'static') {
                //slight hack - canvas needs this element to be positioned
                viewerdiv.style.position = 'relative';
            }

            var UI:any = null;

            type = null;
            if (viewerdiv.dataset.pdb) {
                datauri.push("https://files.rcsb.org/view/" + viewerdiv.dataset.pdb + ".pdb");
                datatypes.push("pdb");
            } else if (viewerdiv.dataset.cid) {
                //this doesn't actually work since pubchem does have CORS enabled
                datatypes.push("sdf");
                datauri.push("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + viewerdiv.dataset.cid +
                    "/SDF?record_type=3d");
            }
            else if (viewerdiv.dataset.href || viewerdiv.dataset.url) {
                if (viewerdiv.dataset.href)
                    uri = viewerdiv.dataset.href;
                else
                    uri = viewerdiv.dataset.url;
                datauri.push(uri);
                type = uri.substring(uri.lastIndexOf('.') + 1);
                datatypes.push(type);

                var molName = uri.substring(uri.lastIndexOf('/') + 1, uri.lastIndexOf('.'));
                if (molName == '/')
                    molName = uri.substring(uri.lastIndexOf('/') + 1);

                viewerdiv.dataset[datatypes[datatypes.length - 1]] = molName;
            }

            var divdata = viewerdiv.dataset;
            for (i in divdata) {
                if ((i.substring(0, 3) === "pdb" && (i !== "pdb"))) {
                    datauri.push("https://files.rcsb.org/view/" + divdata[i] + ".pdb");
                    datatypes.push('pdb');

                } else if (i.substring(0, 4) === "href" && (i !== "href")) {
                    uri = divdata[i];
                    datauri.push(uri);
                    datatypes.push(uri.substr(uri.lastIndexOf('.') + 1));
                } else if (i.substring(0, 3) === "cid" && (i !== "cid")) {
                    datauri.push("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + divdata[i] + "/SDF?record_type=3d");
                    datatypes.push('sdf');
                }
            }
            var options = {};
            if (viewerdiv.dataset.options)
                options = specStringToObject(viewerdiv.dataset.options);

            //note that data tags must be lowercase
            var bgcolor = CC.color(viewerdiv.dataset.backgroundcolor);
            var bgalpha: string | number = viewerdiv.dataset.backgroundalpha;
            bgalpha = (bgalpha == undefined) ? 1.0 : parseFloat(bgalpha);
            var style = { line: {} };
            if (viewerdiv.dataset.style) style = specStringToObject(viewerdiv.dataset.style);
            var select = {};
            if (viewerdiv.dataset.select) select = specStringToObject(viewerdiv.dataset.select);
            var selectstylelist = [];
            var surfaces = [];
            var labels = [];
            var zoomto = {};
            var spin = null;
            var d = viewerdiv.dataset;

            //let users specify individual but matching select/style/surface tags, eg.
            //data-select1 data-style1
            var stylere = /style(.+)/;
            var surfre = /surface(.*)/;
            var reslabre = /labelres(.*)/;
            var keys = [];
            for (dataname in d) {
                if (Object.prototype.hasOwnProperty.call(d, dataname)) {
                    keys.push(dataname);
                }
            }
            keys.sort();
            for (i = 0; i < keys.length; i++) {
                dataname = keys[i];
                var m = stylere.exec(dataname);
                var selname, newsel, styleobj;
                if (m) {
                    selname = "select" + m[1];
                    newsel = specStringToObject(d[selname]);
                    styleobj = specStringToObject(d[dataname]);
                    selectstylelist.push([newsel, styleobj]);
                }
                m = surfre.exec(dataname);
                if (m) {
                    selname = "select" + m[1];
                    newsel = specStringToObject(d[selname]);
                    styleobj = specStringToObject(d[dataname]);
                    surfaces.push([newsel, styleobj]);
                }
                m = reslabre.exec(dataname);
                if (m) {
                    selname = "select" + m[1];
                    newsel = specStringToObject(d[selname]);
                    styleobj = specStringToObject(d[dataname]);
                    labels.push([newsel, styleobj]);
                }
                if (dataname == "zoomto") {
                    zoomto = specStringToObject(d[dataname]);
                }
                if (dataname == "spin") {
                    spin = specStringToObject(d[dataname]);
                }
            }

            //apply all the selections/styles parsed out above to the passed viewer
            var applyStyles = function (glviewer) {
                glviewer.setStyle(select, style);

                if (UI) {
                    UI.createSelectionAndStyle(select, style);
                }

                for (i = 0; i < selectstylelist.length; i++) {
                    let sel = selectstylelist[i][0] || {};
                    let sty = selectstylelist[i][1] || { "line": {} };
                    glviewer.setStyle(sel, sty);
                    if (UI) {
                        UI.createSelectionAndStyle(select, style);
                    }
                }
                for (i = 0; i < surfaces.length; i++) {
                    let sel = surfaces[i][0] || {};
                    let sty = surfaces[i][1] || {};
                    let viewer = glviewer;

                    if (UI) {
                        viewer.addSurface(SurfaceType.VDW, sty, sel, sel).then((surfid) => {
                            UI.loadSurface("VDW", sel, sty, surfid);
                        });
                    }
                    else {
                        glviewer.addSurface(SurfaceType.VDW, sty, sel, sel);
                    }

                }
                for (i = 0; i < labels.length; i++) {
                    let sel = labels[i][0] || {};
                    let sty = labels[i][1] || {};
                    glviewer.addResLabels(sel, sty);
                }

                glviewer.render();
                glviewer.zoomTo(zoomto);

                if (spin) {
                    glviewer.spin(spin.axis, spin.speed);
                }
            };

            var glviewer = viewer;
            try {
                var config: any = specStringToObject(viewerdiv.dataset.config) || {};
                if (config.backgroundColor === undefined) config.backgroundColor = bgcolor;
                if (config.backgroundAlpha === undefined) config.backgroundAlpha = bgalpha;
                if (glviewer == null) {
                    glviewer = viewers[viewerdiv.id || nviewers++] = createViewer(viewerdiv, config);
                } else {
                    glviewer.setBackgroundColor(bgcolor, bgalpha);
                    glviewer.setConfig(config);
                    if (UI) 
                        UI.initiateUI();
                }

                if(viewerdiv.dataset.ui && $3Dmol.StateManager) {
                    UI = new $3Dmol.StateManager(glviewer); // Creates the UI state management tool
                }
            } catch (error) {
                console.log(error);
                //for autoload, provide a useful error message
                viewerdiv.textContent = "WebGL appears to be disabled.";
            }

            if (datauri.length != 0) {
                //load multiple data elements in serial
                i = 0;
                var process = function (moldata) {
                    //add moldata to viewer and load next model
                    uri = datauri[i]; //this is where the moldata came from
                    var type = viewerdiv.dataset.type || viewerdiv.dataset.datatype || datatypes[i];
                    glviewer.addModel(moldata, type, options);
                    if (UI) {
                        var modelName = viewerdiv.dataset[datatypes[i]];
                        UI.setModelTitle(modelName);
                    }
                    i += 1;
                    if (i < datauri.length) {
                        get(datauri[i]).then(process);
                    }
                    else {
                        // or finalize if this is the last model
                        applyStyles(glviewer);
                        if (viewerdiv.dataset.callback) {
                            //evaluate javascript in the string, if it resolves to a function,
                            //call it with the viewer
                            /*jshint -W061 */ var runres = eval(viewerdiv.dataset.callback);
                            if (typeof (runres) == 'function') {
                                runres(glviewer);
                            }
                        }
                        processing_autoinit = false;
                        if (callback) callback(glviewer);
                    }
                };
                get(datauri[0]).then(process);

            }
            else {

                if (viewerdiv.dataset.element) {
                    var moldataid = "#" + viewerdiv.dataset.element;
                    var molelem = document.querySelector(moldataid);
                    var moldata = molelem ? molelem.textContent : "";
                    type = viewerdiv.dataset.type || viewerdiv.dataset.datatype;
                    glviewer.addModel(moldata, type, options);
                }

                applyStyles(glviewer);
                if (viewerdiv.dataset.callback) {
                    //evaluate javascript in the string, if it resolves to a function,
                    //call it with the viewer
                    /*jshint -W061 */ var runres = eval(viewerdiv.dataset.callback);
                    if (typeof (runres) == 'function') {
                        runres(glviewer);
                    }
                }
                processing_autoinit = false;
                if (callback)
                    callback(glviewer);
            }
        });
    }
};



document.onreadystatechange = () => {
    if (document.readyState === "complete") {
        autoload();
    }
};
