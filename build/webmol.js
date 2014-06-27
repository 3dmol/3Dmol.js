
//This defines the WebMol object which is used to create viewers
//and configure system-wide settings

//TODO: automated build system.  (doc gen system, too?)

//the presence of jquery is assumed
/** 
 * WebMol global namespace
 * 
 * @namespace
 * 
 * @property JmolElementColors - Jmol style atom colors (default color scheme)
 * @property rasmolElementColors - Rasmol style atom colors
 */
var WebMol = {
    
    /**
     * Create and initialize an appropriate viewer at supplied HTML element using specification in config
     * @function WebMol.createViewer
     * @param {Object | string} element Either HTML element or string identifier
     * @param {Object} config Viewer specification
     * @returns {WebMol.GLViewer} GLViewer
     * 
     * @example
     * // Assume there exists an HTML div with id "gldiv"
     * var element = $("#gldiv");
     * 
     * // Viewer config - properties 'defaultcolors' and 'callback'
     * var config = {defaultcolors: WebMol.rasmolElementColors,
     *               callback : function(viewer) {
     *                            //'data' is a string containing molecule data in pdb format  
     *                            viewer.addModel(data, "pdb");
     *                            viewer.zoomTo();
     *                            viewer.render();
     *                          }  
     *                        
     *               };
     * 
     * // Create GLViewer within 'gldiv' and execute callback
     * var myviewer = WebMol.createViewer(element, config);
     *      
     */
    createViewer : function(element, config)
    {
        if($.type(element) === "string")
            element = $("#"+element);
        if(!element) return;

        config = config || {};
        if(!config.order)
            config.order = ["glmol","jmol"];
        if(!config.defaultcolors)
            config.defaultcolors = WebMol.defaultElementColors;

        //try to create the appropriate viewer
        for(var i = 0; i < config.order.length; i++) {
            var kind = config.order[i];
            var fname =kind+"Viewer";

            if(typeof(this[fname]) === "function")
            {
                try {
                    return new this[fname](element, config.callback, config.defaultcolors);
                }
                catch(e) {
                    console.log("error with "+kind+":"+e);
                }
            }
        }
        alert("Unable to instantiate webmol viewer: "+config.order);
        return null;
    },

    /**
     * Load a PDB/PubChem structure into existing viewer. Automatically calls 'zoomTo' and 'render' on viewer after loading model
     * 
     * @function WebMol.download
     * @param {string} query String specifying pdb or pubchem id; must be prefaced with "pdb: " or "cid: ", respectively
     * @param {Object} viewer Add new model to existing viewer
     * @example
     * var myviewer = WebMol.createViewer(gldiv);
     * 
     * // GLModel 'm' created and loaded into glviewer for PDB id 2POR
     * var m = WebMol.download('pdb: 2POR', myviewer);
     * 
     * @returns {WebMol.GLModel} GLModel
     */    
    download : function(query, viewer) {
        var baseURL = '';
        var type = "";
        var m = null;
        if (query.substr(0, 4) === 'pdb:') {
            type = "pdb";
            query = query.substr(4).toUpperCase();
            if (!query.match(/^[1-9][A-Za-z0-9]{3}$/)) {
               alert("Wrong PDB ID"); return;
            }
            uri = "http://www.pdb.org/pdb/files/" + query + ".pdb";
        } else if (query.substr(0, 4) == 'cid:') {
            type = "sdf";
            query = query.substr(4);
            if (!query.match(/^[1-9]+$/)) {
               alert("Wrong Compound ID"); return;
            }
            uri = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + query + 
              "/SDF?record_type=3d";
        }

       $.get(uri, function(ret) {
          viewer.addModel(ret, type);
          viewer.zoomTo();
          viewer.render();                            
       });
       
       return m;
    }
       

};

WebMol.SurfaceType = {
    VDW : 1,
    MS : 2,
    SAS : 3,
    SES  : 4
};

// in an attempt to reduce memory overhead, cache all WebMol.Colors
//this makes things a little faster
WebMol.CC = {
    cache : {},
    color : function(hex) {
        if(typeof(this.cache[hex]) !== "undefined") {
            return this.cache[hex];
        }
        else {
            var c = new WebMol.Color(hex);
            this.cache[hex] = c;
            return c;
        }
    }
};

WebMol.Color = function( color ){
    
    if ( arguments.length > 1) {
            this.r = arguments[0] || 0.0;
            this.g = arguments[1] || 0.0;
            this.b = arguments[2] || 0.0;

            return this;
    }
    
    return this.set(color);
                
};

WebMol.Color.prototype = {
    
    constructor: WebMol.Color,
    
    r: 0.0, g: 0.0, b: 0.0,
    
    set : function(val) {
        
            if (val instanceof WebMol.Color) 
                return val.clone();

            else if (typeof val === 'number')
                this.setHex(val);
    },
    
    setHex: function(hex) {
        
            hex = Math.floor(hex);

            this.r = (hex >> 16 & 255) / 255;
            this.g = (hex >> 8 & 255) / 255;
            this.b = (hex & 255) / 255;                                                                                     
        
            return this;
    },
    
    clone : function() {
            return new WebMol.Color(this.r, this.g, this.b);
    },
        
    copy : function(color) {
        this.r = color.r;
        this.g = color.g;
        this.b = color.b;
        
        return this;
    }
    
};

//Miscellaneous functions and classes - to be incorporated into WebMol proper

WebMol.mergeGeos = function(geometry, mesh) {
    
    var meshGeo = mesh.geometry;
    
    if (meshGeo === undefined) 
        return;
    
    geometry.geometryGroups.push( meshGeo.geometryGroups[0] );
    
};

WebMol.multiLineString = function(f) {
    return f.toString()
            .replace(/^[^\/]+\/\*!?/, '')
            .replace(/\*\/[^\/]+$/, '');
            
};

//Synchronized (i.e. not threaded) surface gen? Used mainly for debugging
WebMol.syncSurface = false;
//WebMol constants (replaces needed THREE constants)

//material constants

// sides
WebMol.FrontSide = 0;
WebMol.BackSide = 1;
WebMol.DoubleSide = 2;

// blending modes
WebMol.NoBlending = 0;
WebMol.NormalBlending = 1;
WebMol.AdditiveBlending = 2;
WebMol.SubtractiveBlending = 3;
WebMol.MultiplyBlending = 4;
WebMol.CustomBlending = 5;

// shading
WebMol.NoShading = 0;
WebMol.FlatShading = 1;
WebMol.SmoothShading = 2;

// colors
WebMol.NoColors = 0;
WebMol.FaceColors = 1;
WebMol.VertexColors = 2;

//Texture constants
//TODO: Which of these do I need (since I only use textures to display label sprites) ?
WebMol.MultiplyOperation = 0;
WebMol.MixOperation = 1;
WebMol.AddOperation = 2;

// mapping modes

WebMol.UVMapping = function() {};

// wrapping modes
WebMol.ClampToEdgeWrapping = 1001;

//Filters
WebMol.LinearFilter = 1006;
WebMol.LinearMipMapLinearFilter = 1008;

//Data types
WebMol.UnsignedByteType = 1009;

//Pixel formats
WebMol.RGBAFormat = 1021;


WebMol = WebMol || {};
//Encapsulate marching cube algorithm for isosurface generation
// (currently used by protein surface rendering and generic volumetric data reading)
WebMol.MarchingCube = (function() {
    
    //Marching cube algorithm - assume data has been pre-treated so isovalue is 0 
    // (i.e. select points greater than 0)
    //origin -  vector of origin of volumetric data (default is (0,0,0))
    // nX, nY, nZ - specifies number of voxels in each dimension
    // scale - cube diagonal unit vector scale (webmol vector) (specifying distance between data points); diagonal of cube
    // - default is 1 - assumes unit cube (1,1,1) diag)
    // fulltable - if true, use full marching cubes and tritables - else use trimmed table (e.g. surf render)
    // voxel - if true, draws with a blocky voxel style (default false)
    // verts, faces - vertex and face arrays to fill up
    
    //to match with protein surface...
    var ISDONE = 2;
    var my = {};
    
    my.march = function(data, verts, faces, spec) {

        var fulltable = !!(spec.fulltable);
        var origin = (spec.hasOwnProperty('origin') && spec.origin.hasOwnProperty('x')) ? spec.origin : {x:0, y:0, z:0};
        var voxel = !!(spec.voxel);
        
        var nX = spec.nX || 0;
        var nY = spec.nY || 0;
        var nZ = spec.nZ || 0;
        
        var scale = spec.scale || 1.0;
        
        var unitCube = new WebMol.Vector3(1,1,1).multiplyScalar(scale);
        
        //keep track of calculated vertices to avoid repeats
        var vertnums = new Int32Array(nX*nY*nZ);
        
        var i, il;
        
        for (i = 0, il = vertnums.length; i < il; ++i)
            vertnums[i] = -1;

        // create (or retrieve) a vertex at the appropriate point for
        // the edge (p1,p2)
        
        var getVertex = function(i, j, k, code, p1, p2) {
            var pt = new WebMol.Vector3();
            pt.copy(origin);
            var val1 = !!(code & (1 << p1));
            var val2 = !!(code & (1 << p2));
             
            // p1 if they are the same or if !val1
            var p = p1;
            if (!val1 && val2)
                p = p2;
            
            // adjust i,j,k by p
            if (p & 1)
                k++;
            if (p & 2)
                j++;
            if (p & 4)
                i++;
    
            pt.x += unitCube.x*i;
            pt.y += unitCube.y*j;
            pt.z += unitCube.z*k;
    
            var index = ((nY * i) + j) * nZ + k;
            
            //Have to add option to do voxels
            if (!voxel) {
            
                if (vertnums[index] < 0) // not created yet
                {
                    vertnums[index] = verts.length;
                    verts.push( pt );
                }
                return vertnums[index];
            
            }
            
            else {
                verts.push(pt);
                return verts.length - 1;
            }
            
        };
            
        var intersects = new Int32Array(12);
        
        var etable = (fulltable) ? edgeTable2 : edgeTable;
        var tritable = (fulltable) ? triTable2 : triTable;
                
        //Run marching cubes algorithm
        for (i = 0; i < nX-1; ++i) {
            
            for (var j = 0; j < nY-1; ++j){
                
                for (var k = 0; k < nZ-1; ++k){
                    
                    var code = 0;
                    
                    for (var p = 0; p < 8; ++p) {
                        var index = ((nY * (i + ((p & 4) >> 2))) + j + ((p & 2) >> 1)) *
                                        nZ + k + (p & 1);

                        //TODO: Need to fix vpBits in protein surface for this to work
                        var val = !!(data[index] & ISDONE);
                        //var val = !!(data[index] > 0);   
                        
                        code |= val << p;                        
                    }
                    
                    if (code === 0 || code === 255)
                        continue;
                    
                    var ecode = etable[code];
                    
                    if (ecode === 0)
                        continue;
                        
                    var ttable = tritable[code];                        
                    
                    if (ecode & 1)
                        intersects[0] = getVertex(i, j, k, code, 0, 1);
                    if (ecode & 2)
                        intersects[1] = getVertex(i, j, k, code, 1, 3);
                    if (ecode & 4)
                        intersects[2] = getVertex(i, j, k, code, 3, 2);
                    if (ecode & 8)
                        intersects[3] = getVertex(i, j, k, code, 2, 0);
                    if (ecode & 16)
                        intersects[4] = getVertex(i, j, k, code, 4, 5);
                    if (ecode & 32)
                        intersects[5] = getVertex(i, j, k, code, 5, 7);
                    if (ecode & 64)
                        intersects[6] = getVertex(i, j, k, code, 7, 6);
                    if (ecode & 128)
                        intersects[7] = getVertex(i, j, k, code, 6, 4);
                    if (ecode & 256)
                        intersects[8] = getVertex(i, j, k, code, 0, 4);
                    if (ecode & 512)
                        intersects[9] = getVertex(i, j, k, code, 1, 5);
                    if (ecode & 1024)
                        intersects[10] = getVertex(i, j, k, code, 3, 7);
                    if (ecode & 2048)
                        intersects[11] = getVertex(i, j, k, code, 2, 6);       
                        
                    for (var t = 0; t < ttable.length; t += 3) {
                        
                        var a = intersects[ttable[t]],
                            b = intersects[ttable[t+1]],
                            c = intersects[ttable[t+2]];         
                                           
                        if (voxel && t >= 3) {
                            verts.push(verts[a]); a = verts.length - 1;
                            verts.push(verts[b]); b = verts.length - 1;
                            verts.push(verts[c]); c = verts.length - 1;
                        }

                        
                        faces.push(a); faces.push(b); faces.push(c);                               
                    }              
                    
                }
                
            }
            
        }
             
        
    };

    my.laplacianSmooth = function(numiter, verts, faces) {
            var tps = new Array(verts.length);
            var i, il, j, jl, k, kl;
            for (i = 0, il = verts.length; i < il; i++)
                    tps[i] = {
                        x : 0,
                        y : 0,
                        z : 0
                    };
            var vertdeg = new Array(20);
            var flagvert;
            for (i = 0; i < 20; i++)
                    vertdeg[i] = new Array(verts.length);
            for (i = 0, il = verts.length; i < il; i++)
                    vertdeg[0][i] = 0;
            for (i = 0, il = faces.length / 3; i < il; i++) {
                var aoffset = i*3, boffset = i*3 + 1, coffset = i*3 + 2;
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[aoffset]]; j < jl; j++) {
                    if (faces[boffset] == vertdeg[j + 1][faces[aoffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[aoffset]]++;
                    vertdeg[vertdeg[0][faces[aoffset]]][faces[aoffset]] = faces[boffset];
                }
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[aoffset]]; j < jl; j++) {
                    if (faces[coffset] == vertdeg[j + 1][faces[aoffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[aoffset]]++;
                    vertdeg[vertdeg[0][faces[aoffset]]][faces[aoffset]] = faces[coffset];
                }
                // b
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[boffset]]; j < jl; j++) {
                    if (faces[aoffset] == vertdeg[j + 1][faces[boffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[boffset]]++;
                    vertdeg[vertdeg[0][faces[boffset]]][faces[boffset]] = faces[aoffset];
                }
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[boffset]]; j < jl; j++) {
                    if (faces[coffset] == vertdeg[j + 1][faces[boffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[boffset]]++;
                    vertdeg[vertdeg[0][faces[boffset]]][faces[boffset]] = faces[coffset];
                }
                // c
                flagvert = true;
                for (j = 0; j < vertdeg[0][faces[coffset]]; j++) {
                    if (faces[aoffset] == vertdeg[j + 1][faces[coffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[coffset]]++;
                    vertdeg[vertdeg[0][faces[coffset]]][faces[coffset]] = faces[aoffset];
                }
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[coffset]]; j < jl; j++) {
                    if (faces[boffset] == vertdeg[j + 1][faces[coffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[coffset]]++;
                    vertdeg[vertdeg[0][faces[coffset]]][faces[coffset]] = faces[boffset];
                }
            }

            var wt = 1.00;
            var wt2 = 0.50;
            var ssign;
            var scaleFactor = 1;
            var outwt = 0.75 / (scaleFactor + 3.5); // area-preserving
            for (k = 0; k < numiter; k++) {
                    for (i = 0, il = verts.length; i < il; i++) {
                            if (vertdeg[0][i] < 3) {
                                    tps[i].x = verts[i].x;
                                    tps[i].y = verts[i].y;
                                    tps[i].z = verts[i].z;
                            } else if (vertdeg[0][i] == 3 || vertdeg[0][i] == 4) {
                                    tps[i].x = 0;
                                    tps[i].y = 0;
                                    tps[i].z = 0;
                                    for (j = 0, jl = vertdeg[0][i]; j < jl; j++) {
                                            tps[i].x += verts[vertdeg[j + 1][i]].x;
                                            tps[i].y += verts[vertdeg[j + 1][i]].y;
                                            tps[i].z += verts[vertdeg[j + 1][i]].z;
                                    }
                                    tps[i].x += wt2 * verts[i].x;
                                    tps[i].y += wt2 * verts[i].y;
                                    tps[i].z += wt2 * verts[i].z;
                                    tps[i].x /= wt2 + vertdeg[0][i];
                                    tps[i].y /= wt2 + vertdeg[0][i];
                                    tps[i].z /= wt2 + vertdeg[0][i];
                            } else {
                                    tps[i].x = 0;
                                    tps[i].y = 0;
                                    tps[i].z = 0;
                                    for (j = 0, jl = vertdeg[0][i]; j < jl; j++) {
                                            tps[i].x += verts[vertdeg[j + 1][i]].x;
                                            tps[i].y += verts[vertdeg[j + 1][i]].y;
                                            tps[i].z += verts[vertdeg[j + 1][i]].z;
                                    }
                                    tps[i].x += wt * verts[i].x;
                                    tps[i].y += wt * verts[i].y;
                                    tps[i].z += wt * verts[i].z;
                                    tps[i].x /= wt + vertdeg[0][i];
                                    tps[i].y /= wt + vertdeg[0][i];
                                    tps[i].z /= wt + vertdeg[0][i];
                            }
                    }
                    for (i = 0, il = verts.length; i < il; i++) {
                            verts[i].x = tps[i].x;
                            verts[i].y = tps[i].y;
                            verts[i].z = tps[i].z;
                    }
                    /*
                     * computenorm(); for (var i = 0; i < vertnumber; i++) { if
                     * (verts[i].inout) ssign = 1; else ssign = -1; verts[i].x += ssign *
                     * outwt * verts[i].pn.x; verts[i].y += ssign * outwt *
                     * verts[i].pn.y; verts[i].z += ssign * outwt * verts[i].pn.z; }
                     */
            }
    };


    /*
     * These tables are based off those by Paul Bourke and Geoffrey Heller:
     * http://paulbourke.net/geometry/polygonise/
     * http://paulbourke.net/geometry/polygonise/table2.txt
     * 
     * However, they have been substantially modified to reflect a more 
     * sensible corner numbering scheme and the discrete nature of our voxel data
     * (resulting in fewer faces).
     */
    my.edgeTable = [ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
            0xb00, 0x0, 0x0, 0x0, 0x700, 0x0, 0xd00, 0xe00, 0xf00, 0x0, 0x0, 0x0,
            0x8a, 0x0, 0x15, 0x0, 0x86, 0x0, 0x0, 0x0, 0x28c, 0x0, 0x813, 0xf19,
            0xe10, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x0, 0x0, 0x126, 0x0, 0x0, 0x15, 0x1c,
            0x0, 0xf23, 0x419, 0xd20, 0x0, 0xa8, 0xa2, 0xaa, 0x0, 0x285, 0x9ab,
            0x8a2, 0x0, 0x2af, 0x125, 0xac, 0xfaa, 0xea3, 0xda9, 0xca0, 0x0, 0x0,
            0x0, 0x0, 0x0, 0x45, 0x0, 0x384, 0x0, 0x0, 0x0, 0x700, 0x8a, 0x83,
            0x648, 0x780, 0x0, 0x51, 0x0, 0x81a, 0x54, 0x55, 0x54, 0x56, 0x0, 0x51,
            0x0, 0xe5c, 0x14a, 0x451, 0x759, 0x650, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x45,
            0x0, 0x1f6, 0x0, 0x0, 0x15, 0xdfc, 0x8a, 0x7f3, 0x4f9, 0x5f0, 0xb00,
            0x68, 0x921, 0x6a, 0x348, 0x245, 0x16f, 0x66, 0xb00, 0xe6f, 0xd65,
            0xc6c, 0x76a, 0x663, 0x569, 0x460, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
            0xf46, 0x0, 0x0, 0x45, 0x24c, 0x2a, 0x823, 0x29, 0xb40, 0x0, 0x0, 0x0,
            0x6ba, 0x0, 0x8f5, 0xfff, 0xef6, 0x0, 0xff, 0x2f5, 0x2fc, 0x9ea, 0x8f3,
            0xbf9, 0xaf0, 0x0, 0x0, 0x51, 0x152, 0x0, 0xf55, 0x45f, 0xd56, 0x54,
            0x357, 0x55, 0x154, 0x852, 0xb53, 0x59, 0x950, 0x700, 0x2c8, 0xc2,
            0x48a, 0xfc4, 0xec5, 0xdcf, 0xcc6, 0x2c4, 0x2cf, 0xc5, 0xcc, 0xbca,
            0xac3, 0x9c9, 0x8c0, 0x0, 0x0, 0x0, 0x0, 0xa8, 0x1a4, 0xa8, 0x7a6,
            0xa2, 0xa2, 0x2a4, 0xbac, 0xaa, 0xa3, 0x2a8, 0x3a0, 0xd00, 0xc18,
            0xd00, 0xe3a, 0x34, 0x35, 0x73f, 0x636, 0x924, 0x83f, 0xb35, 0xa3c,
            0x12a, 0x33, 0x339, 0x230, 0xe00, 0xe00, 0xc12, 0xd9a, 0x684, 0x795,
            0x49f, 0x596, 0x92, 0xb9f, 0x815, 0x99c, 0x9a, 0x393, 0x99, 0x190,
            0xf00, 0xe08, 0xd01, 0xc0a, 0x704, 0x605, 0x50f, 0x406, 0xb02, 0xa0f,
            0x905, 0x80c, 0x30a, 0x203, 0x109, 0x0 ];
    
    var edgeTable = new Uint32Array(my.edgeTable);
    
    var triTable = my.triTable = [ [], [], [], [], [], [], [], [ 11, 9, 8 ], [], [], [],
            [ 8, 10, 9 ], [], [ 10, 8, 11 ], [ 9, 11, 10 ],
            [ 8, 10, 9, 8, 11, 10 ], [], [], [], [ 1, 7, 3 ], [], [ 4, 2, 0 ], [],
            [ 2, 1, 7 ], [], [], [], [ 2, 7, 3, 2, 9, 7 ], [],
            [ 1, 4, 11, 1, 0, 4 ], [ 3, 8, 0, 11, 9, 4, 11, 10, 9 ],
            [ 4, 11, 9, 11, 10, 9 ], [], [], [], [ 5, 3, 1 ], [], [], [],
            [ 2, 5, 8, 2, 1, 5 ], [], [], [ 2, 4, 0 ], [ 3, 2, 4 ], [],
            [ 0, 9, 1, 8, 10, 5, 8, 11, 10 ], [ 3, 4, 0, 3, 10, 4 ],
            [ 5, 8, 10, 8, 11, 10 ], [], [ 3, 5, 7 ], [ 7, 1, 5 ],
            [ 1, 7, 3, 1, 5, 7 ], [], [ 9, 2, 0, 9, 7, 2 ],
            [ 0, 3, 8, 1, 7, 11, 1, 5, 7 ], [ 11, 1, 7, 1, 5, 7 ], [],
            [ 9, 1, 0, 5, 3, 2, 5, 7, 3 ], [ 8, 2, 5, 8, 0, 2 ],
            [ 2, 5, 3, 5, 7, 3 ], [ 3, 9, 1, 3, 8, 9, 7, 11, 10, 7, 10, 5 ],
            [ 9, 1, 0, 10, 7, 11, 10, 5, 7 ], [ 3, 8, 0, 7, 10, 5, 7, 11, 10 ],
            [ 11, 5, 7, 11, 10, 5 ], [], [], [], [], [], [ 0, 6, 2 ], [],
            [ 7, 2, 9, 7, 9, 8 ], [], [], [], [ 8, 10, 9 ], [ 7, 1, 3 ],
            [ 7, 1, 0 ], [ 6, 9, 3, 6, 10, 9 ], [ 7, 10, 8, 10, 9, 8 ], [],
            [ 6, 0, 4 ], [], [ 11, 1, 4, 11, 3, 1 ], [ 2, 4, 6 ],
            [ 2, 0, 4, 2, 4, 6 ], [ 2, 4, 6 ], [ 1, 4, 2, 4, 6, 2 ], [],
            [ 6, 0, 4 ], [], [ 2, 11, 3, 6, 9, 4, 6, 10, 9 ], [ 8, 6, 1, 8, 1, 3 ],
            [ 10, 0, 6, 0, 4, 6 ], [ 8, 0, 3, 9, 6, 10, 9, 4, 6 ],
            [ 10, 4, 6, 10, 9, 4 ], [], [], [], [ 5, 3, 1 ], [], [ 0, 6, 2 ], [],
            [ 7, 4, 8, 5, 2, 1, 5, 6, 2 ], [], [], [ 2, 4, 0 ],
            [ 7, 4, 8, 2, 11, 3, 10, 5, 6 ], [ 7, 1, 3 ],
            [ 5, 6, 10, 0, 9, 1, 8, 7, 4 ], [ 5, 6, 10, 7, 0, 3, 7, 4, 0 ],
            [ 10, 5, 6, 4, 8, 7 ], [ 9, 11, 8 ], [ 3, 5, 6 ],
            [ 0, 5, 11, 0, 11, 8 ], [ 6, 3, 5, 3, 1, 5 ], [ 3, 9, 6, 3, 8, 9 ],
            [ 9, 6, 0, 6, 2, 0 ], [ 0, 3, 8, 2, 5, 6, 2, 1, 5 ],
            [ 1, 6, 2, 1, 5, 6 ], [ 9, 11, 8 ], [ 1, 0, 9, 6, 10, 5, 11, 3, 2 ],
            [ 6, 10, 5, 2, 8, 0, 2, 11, 8 ], [ 3, 2, 11, 10, 5, 6 ],
            [ 10, 5, 6, 9, 3, 8, 9, 1, 3 ], [ 0, 9, 1, 5, 6, 10 ],
            [ 8, 0, 3, 10, 5, 6 ], [ 10, 5, 6 ], [], [], [], [], [], [], [],
            [ 1, 10, 2, 9, 11, 6, 9, 8, 11 ], [], [], [ 6, 0, 2 ],
            [ 3, 6, 9, 3, 2, 6 ], [ 3, 5, 1 ], [ 0, 5, 1, 0, 11, 5 ], [ 0, 3, 5 ],
            [ 6, 9, 11, 9, 8, 11 ], [], [], [], [ 4, 5, 9, 7, 1, 10, 7, 3, 1 ], [],
            [ 11, 6, 7, 2, 4, 5, 2, 0, 4 ],
            [ 11, 6, 7, 8, 0, 3, 1, 10, 2, 9, 4, 5 ],
            [ 6, 7, 11, 1, 10, 2, 9, 4, 5 ], [],
            [ 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2 ], [ 9, 4, 5, 0, 6, 7, 0, 2, 6 ],
            [ 4, 5, 9, 6, 3, 2, 6, 7, 3 ], [ 6, 7, 11, 5, 3, 8, 5, 1, 3 ],
            [ 6, 7, 11, 4, 1, 0, 4, 5, 1 ], [ 4, 5, 9, 3, 8, 0, 11, 6, 7 ],
            [ 9, 4, 5, 7, 11, 6 ], [], [], [ 0, 6, 4 ], [ 8, 6, 4, 8, 1, 6 ], [],
            [ 0, 10, 2, 0, 9, 10, 4, 8, 11, 4, 11, 6 ],
            [ 10, 2, 1, 6, 0, 3, 6, 4, 0 ], [ 10, 2, 1, 11, 4, 8, 11, 6, 4 ],
            [ 4, 2, 6 ], [ 1, 0, 9, 2, 4, 8, 2, 6, 4 ], [ 2, 4, 0, 2, 6, 4 ],
            [ 8, 2, 4, 2, 6, 4 ], [ 11, 4, 1, 11, 6, 4 ],
            [ 0, 9, 1, 4, 11, 6, 4, 8, 11 ], [ 3, 6, 0, 6, 4, 0 ],
            [ 8, 6, 4, 8, 11, 6 ], [ 10, 8, 9 ], [ 6, 3, 9, 6, 7, 3 ], [ 6, 7, 1 ],
            [ 10, 7, 1, 7, 3, 1 ], [ 7, 11, 6, 8, 10, 2, 8, 9, 10 ],
            [ 11, 6, 7, 10, 0, 9, 10, 2, 0 ], [ 2, 1, 10, 7, 11, 6, 8, 0, 3 ],
            [ 1, 10, 2, 6, 7, 11 ], [ 7, 2, 6, 7, 9, 2 ],
            [ 1, 0, 9, 3, 6, 7, 3, 2, 6 ], [ 7, 0, 6, 0, 2, 6 ],
            [ 2, 7, 3, 2, 6, 7 ], [ 7, 11, 6, 3, 9, 1, 3, 8, 9 ],
            [ 9, 1, 0, 11, 6, 7 ], [ 0, 3, 8, 11, 6, 7 ], [ 11, 6, 7 ], [], [], [],
            [], [ 5, 3, 7 ], [ 8, 5, 2, 8, 7, 5 ], [ 5, 3, 7 ],
            [ 1, 10, 2, 5, 8, 7, 5, 9, 8 ], [ 1, 7, 5 ], [ 1, 7, 5 ],
            [ 9, 2, 7, 9, 7, 5 ], [ 11, 3, 2, 8, 5, 9, 8, 7, 5 ],
            [ 1, 3, 7, 1, 7, 5 ], [ 0, 7, 1, 7, 5, 1 ], [ 9, 3, 5, 3, 7, 5 ],
            [ 9, 7, 5, 9, 8, 7 ], [ 8, 10, 11 ], [ 3, 4, 10, 3, 10, 11 ],
            [ 8, 10, 11 ], [ 5, 9, 4, 1, 11, 3, 1, 10, 11 ], [ 2, 4, 5 ],
            [ 5, 2, 4, 2, 0, 4 ], [ 0, 3, 8, 5, 9, 4, 10, 2, 1 ],
            [ 2, 1, 10, 9, 4, 5 ], [ 2, 8, 5, 2, 11, 8 ],
            [ 3, 2, 11, 1, 4, 5, 1, 0, 4 ], [ 9, 4, 5, 8, 2, 11, 8, 0, 2 ],
            [ 11, 3, 2, 9, 4, 5 ], [ 8, 5, 3, 5, 1, 3 ], [ 5, 0, 4, 5, 1, 0 ],
            [ 3, 8, 0, 4, 5, 9 ], [ 9, 4, 5 ], [ 11, 9, 10 ], [ 11, 9, 10 ],
            [ 1, 11, 4, 1, 10, 11 ], [ 8, 7, 4, 11, 1, 10, 11, 3, 1 ],
            [ 2, 7, 9, 2, 9, 10 ], [ 4, 8, 7, 0, 10, 2, 0, 9, 10 ],
            [ 2, 1, 10, 0, 7, 4, 0, 3, 7 ], [ 10, 2, 1, 8, 7, 4 ], [ 1, 7, 4 ],
            [ 3, 2, 11, 4, 8, 7, 9, 1, 0 ], [ 11, 4, 2, 4, 0, 2 ],
            [ 2, 11, 3, 7, 4, 8 ], [ 4, 1, 7, 1, 3, 7 ], [ 1, 0, 9, 8, 7, 4 ],
            [ 3, 4, 0, 3, 7, 4 ], [ 8, 7, 4 ], [ 8, 9, 10, 8, 10, 11 ],
            [ 3, 9, 11, 9, 10, 11 ], [ 0, 10, 8, 10, 11, 8 ],
            [ 10, 3, 1, 10, 11, 3 ], [ 2, 8, 10, 8, 9, 10 ], [ 9, 2, 0, 9, 10, 2 ],
            [ 8, 0, 3, 1, 10, 2 ], [ 10, 2, 1 ], [ 1, 11, 9, 11, 8, 9 ],
            [ 11, 3, 2, 0, 9, 1 ], [ 11, 0, 2, 11, 8, 0 ], [ 11, 3, 2 ],
            [ 8, 1, 3, 8, 9, 1 ], [ 9, 1, 0 ], [ 8, 0, 3 ], [] ];
     
    var edgeTable2 = [ 0x0, 0x109, 0x203, 0x30a, 0x80c, 0x905, 0xa0f,
            0xb06, 0x406, 0x50f, 0x605, 0x70c, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190,
            0x99, 0x393, 0x29a, 0x99c, 0x895, 0xb9f, 0xa96, 0x596, 0x49f, 0x795,
            0x69c, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0xa3c,
            0xb35, 0x83f, 0x936, 0x636, 0x73f, 0x435, 0x53c, 0xe3a, 0xf33, 0xc39,
            0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa, 0xbac, 0xaa5, 0x9af, 0x8a6, 0x7a6,
            0x6af, 0x5a5, 0x4ac, 0xfaa, 0xea3, 0xda9, 0xca0, 0x8c0, 0x9c9, 0xac3,
            0xbca, 0xcc, 0x1c5, 0x2cf, 0x3c6, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x4ca,
            0x5c3, 0x6c9, 0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0x15c, 0x55, 0x35f,
            0x256, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x55a, 0x453, 0x759, 0x650, 0xaf0,
            0xbf9, 0x8f3, 0x9fa, 0x2fc, 0x3f5, 0xff, 0x1f6, 0xef6, 0xfff, 0xcf5,
            0xdfc, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0x36c,
            0x265, 0x16f, 0x66, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x76a, 0x663, 0x569,
            0x460, 0x460, 0x569, 0x663, 0x76a, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x66,
            0x16f, 0x265, 0x36c, 0x86a, 0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3,
            0x6fa, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x1f6, 0xff, 0x3f5, 0x2fc, 0x9fa,
            0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0xe5c, 0xf55, 0xc5f,
            0xd56, 0x256, 0x35f, 0x55, 0x15c, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0,
            0x6c9, 0x5c3, 0x4ca, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0x3c6, 0x2cf, 0x1c5,
            0xcc, 0xbca, 0xac3, 0x9c9, 0x8c0, 0xca0, 0xda9, 0xea3, 0xfaa, 0x4ac,
            0x5a5, 0x6af, 0x7a6, 0x8a6, 0x9af, 0xaa5, 0xbac, 0xaa, 0x1a3, 0x2a9,
            0x3a0, 0xd30, 0xc39, 0xf33, 0xe3a, 0x53c, 0x435, 0x73f, 0x636, 0x936,
            0x83f, 0xb35, 0xa3c, 0x13a, 0x33, 0x339, 0x230, 0xe90, 0xf99, 0xc93,
            0xd9a, 0x69c, 0x795, 0x49f, 0x596, 0xa96, 0xb9f, 0x895, 0x99c, 0x29a,
            0x393, 0x99, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0x70c, 0x605, 0x50f,
            0x406, 0xb06, 0xa0f, 0x905, 0x80c, 0x30a, 0x203, 0x109, 0x0 ];
     
    var triTable2 = [ [], [ 8, 3, 0 ], [ 9, 0, 1 ], [ 8, 3, 1, 8, 1, 9 ],
            [ 11, 2, 3 ], [ 11, 2, 0, 11, 0, 8 ], [ 11, 2, 3, 0, 1, 9 ],
            [ 2, 1, 11, 1, 9, 11, 11, 9, 8 ], [ 10, 1, 2 ], [ 8, 3, 0, 1, 2, 10 ],
            [ 9, 0, 2, 9, 2, 10 ], [ 3, 2, 8, 2, 10, 8, 8, 10, 9 ],
            [ 10, 1, 3, 10, 3, 11 ], [ 1, 0, 10, 0, 8, 10, 10, 8, 11 ],
            [ 0, 3, 9, 3, 11, 9, 9, 11, 10 ], [ 8, 10, 9, 8, 11, 10 ], [ 8, 4, 7 ],
            [ 3, 0, 4, 3, 4, 7 ], [ 1, 9, 0, 8, 4, 7 ],
            [ 9, 4, 1, 4, 7, 1, 1, 7, 3 ], [ 2, 3, 11, 7, 8, 4 ],
            [ 7, 11, 4, 11, 2, 4, 4, 2, 0 ], [ 3, 11, 2, 4, 7, 8, 9, 0, 1 ],
            [ 2, 7, 11, 2, 1, 7, 1, 4, 7, 1, 9, 4 ], [ 10, 1, 2, 8, 4, 7 ],
            [ 2, 10, 1, 0, 4, 7, 0, 7, 3 ], [ 4, 7, 8, 0, 2, 10, 0, 10, 9 ],
            [ 2, 7, 3, 2, 9, 7, 7, 9, 4, 2, 10, 9 ],
            [ 8, 4, 7, 11, 10, 1, 11, 1, 3 ],
            [ 11, 4, 7, 1, 4, 11, 1, 11, 10, 1, 0, 4 ],
            [ 3, 8, 0, 7, 11, 4, 11, 9, 4, 11, 10, 9 ],
            [ 7, 11, 4, 4, 11, 9, 11, 10, 9 ], [ 9, 5, 4 ], [ 3, 0, 8, 4, 9, 5 ],
            [ 5, 4, 0, 5, 0, 1 ], [ 4, 8, 5, 8, 3, 5, 5, 3, 1 ],
            [ 11, 2, 3, 9, 5, 4 ], [ 9, 5, 4, 8, 11, 2, 8, 2, 0 ],
            [ 3, 11, 2, 1, 5, 4, 1, 4, 0 ],
            [ 8, 5, 4, 2, 5, 8, 2, 8, 11, 2, 1, 5 ], [ 2, 10, 1, 9, 5, 4 ],
            [ 0, 8, 3, 5, 4, 9, 10, 1, 2 ], [ 10, 5, 2, 5, 4, 2, 2, 4, 0 ],
            [ 3, 4, 8, 3, 2, 4, 2, 5, 4, 2, 10, 5 ],
            [ 5, 4, 9, 1, 3, 11, 1, 11, 10 ],
            [ 0, 9, 1, 4, 8, 5, 8, 10, 5, 8, 11, 10 ],
            [ 3, 4, 0, 3, 10, 4, 4, 10, 5, 3, 11, 10 ],
            [ 4, 8, 5, 5, 8, 10, 8, 11, 10 ], [ 9, 5, 7, 9, 7, 8 ],
            [ 0, 9, 3, 9, 5, 3, 3, 5, 7 ], [ 8, 0, 7, 0, 1, 7, 7, 1, 5 ],
            [ 1, 7, 3, 1, 5, 7 ], [ 11, 2, 3, 8, 9, 5, 8, 5, 7 ],
            [ 9, 2, 0, 9, 7, 2, 2, 7, 11, 9, 5, 7 ],
            [ 0, 3, 8, 2, 1, 11, 1, 7, 11, 1, 5, 7 ],
            [ 2, 1, 11, 11, 1, 7, 1, 5, 7 ], [ 1, 2, 10, 5, 7, 8, 5, 8, 9 ],
            [ 9, 1, 0, 10, 5, 2, 5, 3, 2, 5, 7, 3 ],
            [ 5, 2, 10, 8, 2, 5, 8, 5, 7, 8, 0, 2 ],
            [ 10, 5, 2, 2, 5, 3, 5, 7, 3 ],
            [ 3, 9, 1, 3, 8, 9, 7, 11, 10, 7, 10, 5 ],
            [ 9, 1, 0, 10, 7, 11, 10, 5, 7 ], [ 3, 8, 0, 7, 10, 5, 7, 11, 10 ],
            [ 11, 5, 7, 11, 10, 5 ], [ 11, 7, 6 ], [ 0, 8, 3, 11, 7, 6 ],
            [ 9, 0, 1, 11, 7, 6 ], [ 7, 6, 11, 3, 1, 9, 3, 9, 8 ],
            [ 2, 3, 7, 2, 7, 6 ], [ 8, 7, 0, 7, 6, 0, 0, 6, 2 ],
            [ 1, 9, 0, 3, 7, 6, 3, 6, 2 ], [ 7, 6, 2, 7, 2, 9, 2, 1, 9, 7, 9, 8 ],
            [ 1, 2, 10, 6, 11, 7 ], [ 2, 10, 1, 7, 6, 11, 8, 3, 0 ],
            [ 11, 7, 6, 10, 9, 0, 10, 0, 2 ],
            [ 7, 6, 11, 3, 2, 8, 8, 2, 10, 8, 10, 9 ],
            [ 6, 10, 7, 10, 1, 7, 7, 1, 3 ],
            [ 6, 10, 1, 6, 1, 7, 7, 1, 0, 7, 0, 8 ],
            [ 9, 0, 3, 6, 9, 3, 6, 10, 9, 6, 3, 7 ],
            [ 6, 10, 7, 7, 10, 8, 10, 9, 8 ], [ 8, 4, 6, 8, 6, 11 ],
            [ 11, 3, 6, 3, 0, 6, 6, 0, 4 ], [ 0, 1, 9, 4, 6, 11, 4, 11, 8 ],
            [ 1, 9, 4, 11, 1, 4, 11, 3, 1, 11, 4, 6 ],
            [ 3, 8, 2, 8, 4, 2, 2, 4, 6 ], [ 2, 0, 4, 2, 4, 6 ],
            [ 1, 9, 0, 3, 8, 2, 2, 8, 4, 2, 4, 6 ], [ 9, 4, 1, 1, 4, 2, 4, 6, 2 ],
            [ 10, 1, 2, 11, 8, 4, 11, 4, 6 ],
            [ 10, 1, 2, 11, 3, 6, 6, 3, 0, 6, 0, 4 ],
            [ 0, 2, 10, 0, 10, 9, 4, 11, 8, 4, 6, 11 ],
            [ 2, 11, 3, 6, 9, 4, 6, 10, 9 ],
            [ 8, 4, 6, 8, 6, 1, 6, 10, 1, 8, 1, 3 ],
            [ 1, 0, 10, 10, 0, 6, 0, 4, 6 ], [ 8, 0, 3, 9, 6, 10, 9, 4, 6 ],
            [ 10, 4, 6, 10, 9, 4 ], [ 9, 5, 4, 7, 6, 11 ],
            [ 4, 9, 5, 3, 0, 8, 11, 7, 6 ], [ 6, 11, 7, 4, 0, 1, 4, 1, 5 ],
            [ 6, 11, 7, 4, 8, 5, 5, 8, 3, 5, 3, 1 ], [ 4, 9, 5, 6, 2, 3, 6, 3, 7 ],
            [ 9, 5, 4, 8, 7, 0, 0, 7, 6, 0, 6, 2 ],
            [ 4, 0, 1, 4, 1, 5, 6, 3, 7, 6, 2, 3 ], [ 7, 4, 8, 5, 2, 1, 5, 6, 2 ],
            [ 6, 11, 7, 1, 2, 10, 9, 5, 4 ],
            [ 11, 7, 6, 8, 3, 0, 1, 2, 10, 9, 5, 4 ],
            [ 11, 7, 6, 10, 5, 2, 2, 5, 4, 2, 4, 0 ],
            [ 7, 4, 8, 2, 11, 3, 10, 5, 6 ],
            [ 4, 9, 5, 6, 10, 7, 7, 10, 1, 7, 1, 3 ],
            [ 5, 6, 10, 0, 9, 1, 8, 7, 4 ], [ 5, 6, 10, 7, 0, 3, 7, 4, 0 ],
            [ 10, 5, 6, 4, 8, 7 ], [ 5, 6, 9, 6, 11, 9, 9, 11, 8 ],
            [ 0, 9, 5, 0, 5, 3, 3, 5, 6, 3, 6, 11 ],
            [ 0, 1, 5, 0, 5, 11, 5, 6, 11, 0, 11, 8 ],
            [ 11, 3, 6, 6, 3, 5, 3, 1, 5 ], [ 9, 5, 6, 3, 9, 6, 3, 8, 9, 3, 6, 2 ],
            [ 5, 6, 9, 9, 6, 0, 6, 2, 0 ], [ 0, 3, 8, 2, 5, 6, 2, 1, 5 ],
            [ 1, 6, 2, 1, 5, 6 ], [ 1, 2, 10, 5, 6, 9, 9, 6, 11, 9, 11, 8 ],
            [ 1, 0, 9, 6, 10, 5, 11, 3, 2 ], [ 6, 10, 5, 2, 8, 0, 2, 11, 8 ],
            [ 3, 2, 11, 10, 5, 6 ], [ 10, 5, 6, 9, 3, 8, 9, 1, 3 ],
            [ 0, 9, 1, 5, 6, 10 ], [ 8, 0, 3, 10, 5, 6 ], [ 10, 5, 6 ],
            [ 10, 6, 5 ], [ 8, 3, 0, 10, 6, 5 ], [ 0, 1, 9, 5, 10, 6 ],
            [ 10, 6, 5, 9, 8, 3, 9, 3, 1 ], [ 3, 11, 2, 10, 6, 5 ],
            [ 6, 5, 10, 2, 0, 8, 2, 8, 11 ], [ 1, 9, 0, 6, 5, 10, 11, 2, 3 ],
            [ 1, 10, 2, 5, 9, 6, 9, 11, 6, 9, 8, 11 ], [ 1, 2, 6, 1, 6, 5 ],
            [ 0, 8, 3, 2, 6, 5, 2, 5, 1 ], [ 5, 9, 6, 9, 0, 6, 6, 0, 2 ],
            [ 9, 6, 5, 3, 6, 9, 3, 9, 8, 3, 2, 6 ], [ 11, 6, 3, 6, 5, 3, 3, 5, 1 ],
            [ 0, 5, 1, 0, 11, 5, 5, 11, 6, 0, 8, 11 ],
            [ 0, 5, 9, 0, 3, 5, 3, 6, 5, 3, 11, 6 ],
            [ 5, 9, 6, 6, 9, 11, 9, 8, 11 ], [ 10, 6, 5, 4, 7, 8 ],
            [ 5, 10, 6, 7, 3, 0, 7, 0, 4 ], [ 5, 10, 6, 0, 1, 9, 8, 4, 7 ],
            [ 4, 5, 9, 6, 7, 10, 7, 1, 10, 7, 3, 1 ],
            [ 7, 8, 4, 2, 3, 11, 10, 6, 5 ],
            [ 11, 6, 7, 10, 2, 5, 2, 4, 5, 2, 0, 4 ],
            [ 11, 6, 7, 8, 0, 3, 1, 10, 2, 9, 4, 5 ],
            [ 6, 7, 11, 1, 10, 2, 9, 4, 5 ], [ 7, 8, 4, 5, 1, 2, 5, 2, 6 ],
            [ 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2 ],
            [ 9, 4, 5, 8, 0, 7, 0, 6, 7, 0, 2, 6 ], [ 4, 5, 9, 6, 3, 2, 6, 7, 3 ],
            [ 6, 7, 11, 4, 5, 8, 5, 3, 8, 5, 1, 3 ],
            [ 6, 7, 11, 4, 1, 0, 4, 5, 1 ], [ 4, 5, 9, 3, 8, 0, 11, 6, 7 ],
            [ 9, 4, 5, 7, 11, 6 ], [ 10, 6, 4, 10, 4, 9 ],
            [ 8, 3, 0, 9, 10, 6, 9, 6, 4 ], [ 1, 10, 0, 10, 6, 0, 0, 6, 4 ],
            [ 8, 6, 4, 8, 1, 6, 6, 1, 10, 8, 3, 1 ],
            [ 2, 3, 11, 6, 4, 9, 6, 9, 10 ],
            [ 0, 10, 2, 0, 9, 10, 4, 8, 11, 4, 11, 6 ],
            [ 10, 2, 1, 11, 6, 3, 6, 0, 3, 6, 4, 0 ],
            [ 10, 2, 1, 11, 4, 8, 11, 6, 4 ], [ 9, 1, 4, 1, 2, 4, 4, 2, 6 ],
            [ 1, 0, 9, 3, 2, 8, 2, 4, 8, 2, 6, 4 ], [ 2, 4, 0, 2, 6, 4 ],
            [ 3, 2, 8, 8, 2, 4, 2, 6, 4 ],
            [ 1, 4, 9, 11, 4, 1, 11, 1, 3, 11, 6, 4 ],
            [ 0, 9, 1, 4, 11, 6, 4, 8, 11 ], [ 11, 6, 3, 3, 6, 0, 6, 4, 0 ],
            [ 8, 6, 4, 8, 11, 6 ], [ 6, 7, 10, 7, 8, 10, 10, 8, 9 ],
            [ 9, 3, 0, 6, 3, 9, 6, 9, 10, 6, 7, 3 ],
            [ 6, 1, 10, 6, 7, 1, 7, 0, 1, 7, 8, 0 ],
            [ 6, 7, 10, 10, 7, 1, 7, 3, 1 ],
            [ 7, 11, 6, 3, 8, 2, 8, 10, 2, 8, 9, 10 ],
            [ 11, 6, 7, 10, 0, 9, 10, 2, 0 ], [ 2, 1, 10, 7, 11, 6, 8, 0, 3 ],
            [ 1, 10, 2, 6, 7, 11 ], [ 7, 2, 6, 7, 9, 2, 2, 9, 1, 7, 8, 9 ],
            [ 1, 0, 9, 3, 6, 7, 3, 2, 6 ], [ 8, 0, 7, 7, 0, 6, 0, 2, 6 ],
            [ 2, 7, 3, 2, 6, 7 ], [ 7, 11, 6, 3, 9, 1, 3, 8, 9 ],
            [ 9, 1, 0, 11, 6, 7 ], [ 0, 3, 8, 11, 6, 7 ], [ 11, 6, 7 ],
            [ 11, 7, 5, 11, 5, 10 ], [ 3, 0, 8, 7, 5, 10, 7, 10, 11 ],
            [ 9, 0, 1, 10, 11, 7, 10, 7, 5 ],
            [ 3, 1, 9, 3, 9, 8, 7, 10, 11, 7, 5, 10 ],
            [ 10, 2, 5, 2, 3, 5, 5, 3, 7 ],
            [ 5, 10, 2, 8, 5, 2, 8, 7, 5, 8, 2, 0 ],
            [ 9, 0, 1, 10, 2, 5, 5, 2, 3, 5, 3, 7 ],
            [ 1, 10, 2, 5, 8, 7, 5, 9, 8 ], [ 2, 11, 1, 11, 7, 1, 1, 7, 5 ],
            [ 0, 8, 3, 2, 11, 1, 1, 11, 7, 1, 7, 5 ],
            [ 9, 0, 2, 9, 2, 7, 2, 11, 7, 9, 7, 5 ],
            [ 11, 3, 2, 8, 5, 9, 8, 7, 5 ], [ 1, 3, 7, 1, 7, 5 ],
            [ 8, 7, 0, 0, 7, 1, 7, 5, 1 ], [ 0, 3, 9, 9, 3, 5, 3, 7, 5 ],
            [ 9, 7, 5, 9, 8, 7 ], [ 4, 5, 8, 5, 10, 8, 8, 10, 11 ],
            [ 3, 0, 4, 3, 4, 10, 4, 5, 10, 3, 10, 11 ],
            [ 0, 1, 9, 4, 5, 8, 8, 5, 10, 8, 10, 11 ],
            [ 5, 9, 4, 1, 11, 3, 1, 10, 11 ],
            [ 3, 8, 4, 3, 4, 2, 2, 4, 5, 2, 5, 10 ],
            [ 10, 2, 5, 5, 2, 4, 2, 0, 4 ], [ 0, 3, 8, 5, 9, 4, 10, 2, 1 ],
            [ 2, 1, 10, 9, 4, 5 ], [ 8, 4, 5, 2, 8, 5, 2, 11, 8, 2, 5, 1 ],
            [ 3, 2, 11, 1, 4, 5, 1, 0, 4 ], [ 9, 4, 5, 8, 2, 11, 8, 0, 2 ],
            [ 11, 3, 2, 9, 4, 5 ], [ 4, 5, 8, 8, 5, 3, 5, 1, 3 ],
            [ 5, 0, 4, 5, 1, 0 ], [ 3, 8, 0, 4, 5, 9 ], [ 9, 4, 5 ],
            [ 7, 4, 11, 4, 9, 11, 11, 9, 10 ],
            [ 3, 0, 8, 7, 4, 11, 11, 4, 9, 11, 9, 10 ],
            [ 11, 7, 4, 1, 11, 4, 1, 10, 11, 1, 4, 0 ],
            [ 8, 7, 4, 11, 1, 10, 11, 3, 1 ],
            [ 2, 3, 7, 2, 7, 9, 7, 4, 9, 2, 9, 10 ],
            [ 4, 8, 7, 0, 10, 2, 0, 9, 10 ], [ 2, 1, 10, 0, 7, 4, 0, 3, 7 ],
            [ 10, 2, 1, 8, 7, 4 ], [ 2, 11, 7, 2, 7, 1, 1, 7, 4, 1, 4, 9 ],
            [ 3, 2, 11, 4, 8, 7, 9, 1, 0 ], [ 7, 4, 11, 11, 4, 2, 4, 0, 2 ],
            [ 2, 11, 3, 7, 4, 8 ], [ 9, 1, 4, 4, 1, 7, 1, 3, 7 ],
            [ 1, 0, 9, 8, 7, 4 ], [ 3, 4, 0, 3, 7, 4 ], [ 8, 7, 4 ],
            [ 8, 9, 10, 8, 10, 11 ], [ 0, 9, 3, 3, 9, 11, 9, 10, 11 ],
            [ 1, 10, 0, 0, 10, 8, 10, 11, 8 ], [ 10, 3, 1, 10, 11, 3 ],
            [ 3, 8, 2, 2, 8, 10, 8, 9, 10 ], [ 9, 2, 0, 9, 10, 2 ],
            [ 8, 0, 3, 1, 10, 2 ], [ 10, 2, 1 ], [ 2, 11, 1, 1, 11, 9, 11, 8, 9 ],
            [ 11, 3, 2, 0, 9, 1 ], [ 11, 0, 2, 11, 8, 0 ], [ 11, 3, 2 ],
            [ 8, 1, 3, 8, 9, 1 ], [ 9, 1, 0 ], [ 8, 0, 3 ], [] ];
            
            return my;
})();



/*  ProteinSurface.js by biochem_fan

Ported and modified for Javascript based on EDTSurf,
  whose license is as follows.

Permission to use, copy, modify, and distribute this program for any
purpose, with or without fee, is hereby granted, provided that this
copyright notice and the reference information appear in all copies or
substantial portions of the Software. It is provided "as is" without
express or implied warranty. 

Reference:
http://zhanglab.ccmb.med.umich.edu/EDTSurf/
D. Xu, Y. Zhang (2009) Generating Triangulated Macromolecular Surfaces
by Euclidean Distance Transform. PLoS ONE 4(12): e8140.

=======

TODO: Improved performance on Firefox
      Reduce memory consumption
      Refactor!
 */


// dkoes
// Surface calculations.  This must be safe to use within a web worker.
if (typeof console === 'undefined') {
    // this should only be true inside of a webworker
    console = {
        log : function() {
        }
    };
}

WebMol.ProteinSurface = function() {

    // constants for vpbits bitmasks
    var INOUT = 1;
    var ISDONE = 2;
    var ISBOUND = 4;

    var ptranx = 0, ptrany = 0, ptranz = 0;
    var probeRadius = 1.4;
    var defaultScaleFactor = 2;
    var scaleFactor = defaultScaleFactor; // 2 is .5A grid; if this is made user configurable,
                            // also have to adjust offset used to find non-shown
                            // atoms
    var pHeight = 0, pWidth = 0, pLength = 0;
    var cutRadius = 0;
    var vpBits = null; // uint8 array of bitmasks
    var vpDistance = null; // floatarray of _squared_ distances
    var vpAtomID = null; // intarray
    var vertnumber = 0, facenumber = 0;
    var pminx = 0, pminy = 0, pminz = 0, pmaxx = 0, pmaxy = 0, pmaxz = 0;

    var vdwRadii = {
            "H" : 1.2,
            "Li" : 1.82,
            "Na" : 2.27,
            "K" : 2.75,
            "C" : 1.7,
            "N" : 1.55,
            "O" : 1.52,
            "F" : 1.47,
            "P" : 1.80,
            "S" : 1.80,
            "CL" : 1.75,
            "BR" : 1.85,
            "SE" : 1.90,
            "ZN" : 1.39,
            "CU" : 1.4,
            "NI" : 1.63,
            "X" : 2
        };

    var getVDWIndex = function(atom) {
        if(!atom.elem || typeof(vdwRadii[atom.elem]) == "undefined") {
            return "X";
        }
        return atom.elem;
    };
    
    var depty = {}, widxz = {};
    var faces, verts;
    var nb = [ new Int32Array([ 1, 0, 0 ]), new Int32Array([ -1, 0, 0 ]), 
               new Int32Array([ 0, 1, 0 ]), new Int32Array([ 0, -1, 0 ]),
               new Int32Array([ 0, 0, 1 ]), 
               new Int32Array([ 0, 0, -1 ]), 
               new Int32Array([ 1, 1, 0 ]), 
               new Int32Array([ 1, -1, 0 ]), 
               new Int32Array([ -1, 1, 0 ]),
               new Int32Array([ -1, -1, 0 ]), 
               new Int32Array([ 1, 0, 1 ]), 
               new Int32Array([ 1, 0, -1 ]), 
               new Int32Array([ -1, 0, 1 ]),
               new Int32Array([ -1, 0, -1 ]), 
               new Int32Array([ 0, 1, 1 ]), 
               new Int32Array([ 0, 1, -1 ]), 
               new Int32Array([ 0, -1, 1 ]),
               new Int32Array([ 0, -1, -1 ]), 
               new Int32Array([ 1, 1, 1 ]), 
               new Int32Array([ 1, 1, -1 ]), 
               new Int32Array([ 1, -1, 1 ]),
               new Int32Array([ -1, 1, 1 ]), 
               new Int32Array([ 1, -1, -1 ]), 
               new Int32Array([ -1, -1, 1 ]), 
               new Int32Array([ -1, 1, -1 ]),
               new Int32Array([ -1, -1, -1 ]) ];

    var origextent;

    var inOrigExtent = function(x, y, z) {
        if (x < origextent[0][0] || x > origextent[1][0])
            return false;
        if (y < origextent[0][1] || y > origextent[1][1])
            return false;
        if (z < origextent[0][2] || z > origextent[1][2])
            return false;
        return true;
    };

    this.getFacesAndVertices = function(atomlist) {
        var atomsToShow = {};
        var i, il;
        for (i = 0, il = atomlist.length; i < il; i++)
            atomsToShow[atomlist[i]] = true;
        var vertices = verts;
        for (i = 0, il = vertices.length; i < il; i++) {
            vertices[i].x = vertices[i].x / scaleFactor - ptranx;
            vertices[i].y = vertices[i].y / scaleFactor - ptrany;
            vertices[i].z = vertices[i].z / scaleFactor - ptranz;
        }

        var finalfaces = [];
        for (i = 0, il = faces.length; i < il; i += 3) {
            //var f = faces[i];
            var fa = faces[i], fb = faces[i+1], fc = faces[i+2];
            var a = vertices[fa].atomid, b = vertices[fb].atomid, c = vertices[fc].atomid;

            // must be a unique face for each atom
            var which = a;
            if (b < which)
                which = b;
            if (c < which)
                which = c;
            if (!atomsToShow[which]) {
                continue;
            }
            var av = vertices[faces[i]];
            var bv = vertices[faces[i+1]];
            var cv = vertices[faces[i+2]];

            if (fa !== fb && fb !== fc && fa !== fc){
                finalfaces.push(fa); 
                finalfaces.push(fb); 
                finalfaces.push(fc); 
            }
               
        }

        //try to help the garbage collector
        vpBits = null; // uint8 array of bitmasks
        vpDistance = null; // floatarray
        vpAtomID = null; // intarray
        
        return {
            vertices : vertices,
            faces : finalfaces
        };
    };


    this.initparm = function(extent, btype, volume) {
        if(volume > 1000000) //heuristical decrease resolution to avoid large memory consumption
            scaleFactor = defaultScaleFactor/2;
        
        var margin = (1 / scaleFactor) * 5.5; // need margin to avoid
                                                // boundary/round off effects
        origextent = extent;
        pminx = extent[0][0]; pmaxx = extent[1][0];
        pminy = extent[0][1]; pmaxy = extent[1][1];
        pminz = extent[0][2]; pmaxz = extent[1][2];

        if (!btype) {
            pminx -= margin;
            pminy -= margin;
            pminz -= margin;
            pmaxx += margin;
            pmaxy += margin;
            pmaxz += margin;
        } else {
            pminx -= probeRadius + margin;
            pminy -= probeRadius + margin;
            pminz -= probeRadius + margin;
            pmaxx += probeRadius + margin;
            pmaxy += probeRadius + margin;
            pmaxz += probeRadius + margin;
        }

        pminx = Math.floor(pminx * scaleFactor) / scaleFactor;
        pminy = Math.floor(pminy * scaleFactor) / scaleFactor;
        pminz = Math.floor(pminz * scaleFactor) / scaleFactor;
        pmaxx = Math.ceil(pmaxx * scaleFactor) / scaleFactor;
        pmaxy = Math.ceil(pmaxy * scaleFactor) / scaleFactor;
        pmaxz = Math.ceil(pmaxz * scaleFactor) / scaleFactor;

        ptranx = -pminx;
        ptrany = -pminy;
        ptranz = -pminz;

        pLength = Math.ceil(scaleFactor * (pmaxx - pminx)) + 1;
        pWidth = Math.ceil(scaleFactor * (pmaxy - pminy)) + 1;
        pHeight = Math.ceil(scaleFactor * (pmaxz - pminz)) + 1;

        this.boundingatom(btype);
        cutRadius = probeRadius * scaleFactor;

        vpBits = new Uint8Array(pLength * pWidth * pHeight);
        vpDistance = new Float64Array(pLength * pWidth * pHeight); // float 32
        // doesn't
        // play
        // nicely
        // with
        // native
        // floats
        vpAtomID = new Int32Array(pLength * pWidth * pHeight);
        console.log("Box size: ", pLength, pWidth, pHeight, vpBits.length);
    };

    this.boundingatom = function(btype) {
        var tradius = [];
        var txz, tdept, sradius, idx;
        flagradius = btype;

        for ( var i in vdwRadii) {
            if(!vdwRadii.hasOwnProperty(i))
                continue;
            var r = vdwRadii[i];
            if (!btype)
                tradius[i] = r * scaleFactor + 0.5;
            else
                tradius[i] = (r + probeRadius) * scaleFactor + 0.5;

            sradius = tradius[i] * tradius[i];
            widxz[i] = Math.floor(tradius[i]) + 1;
            depty[i] = new Int32Array(widxz[i] * widxz[i]);
            indx = 0;
            for (j = 0; j < widxz[i]; j++) {
                for (k = 0; k < widxz[i]; k++) {
                    txz = j * j + k * k;
                    if (txz > sradius)
                        depty[i][indx] = -1; // outside
                    else {
                        tdept = Math.sqrt(sradius - txz);
                        depty[i][indx] = Math.floor(tdept);
                    }
                    indx++;
                }
            }
        }
    };

    this.fillvoxels = function(atoms, atomlist) { // (int seqinit,int
        // seqterm,bool
        // atomtype,atom*
        // proseq,bool bcolor)
        var i, il;
        for (i = 0, il = vpBits.length; i < il; i++) {
            vpBits[i] = 0;
            vpDistance[i] = -1.0;
            vpAtomID[i] = -1;
        }

        for (i in atomlist) {
            var atom = atoms[atomlist[i]];
            if (atom === undefined)
                continue;
            this.fillAtom(atom, atoms);
        }

        for (i = 0, il = vpBits.length; i < il; i++)
            if (vpBits[i] & INOUT)
                vpBits[i] |= ISDONE;

    };


    this.fillAtom = function(atom, atoms) {
        var cx, cy, cz, ox, oy, oz, mi, mj, mk, i, j, k, si, sj, sk;
        var ii, jj, kk, n;
        cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
        cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
        cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

        var at = getVDWIndex(atom);
        var nind = 0;
        var cnt = 0;
        var pWH = pWidth*pHeight;
        
        for (i = 0, n = widxz[at]; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (depty[at][nind] != -1) {
                    for (ii = -1; ii < 2; ii++) {
                        for (jj = -1; jj < 2; jj++) {
                            for (kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    mi = ii * i;
                                    mk = kk * j;
                                    for (k = 0; k <= depty[at][nind]; k++) {
                                        mj = k * jj;
                                        si = cx + mi;
                                        sj = cy + mj;
                                        sk = cz + mk;
                                        if (si < 0 || sj < 0 || 
                                                sk < 0 ||
                                                si >= pLength || 
                                                sj >= pWidth || 
                                                sk >= pHeight)
                                            continue;
                                        var index = si * pWH + sj * pHeight + sk;

                                        if (!(vpBits[index] & INOUT)) {
                                            vpBits[index] |= INOUT;
                                            vpAtomID[index] = atom.serial;
                                        } else {
                                            var atom2 = atoms[vpAtomID[index]];
                                            ox = Math.floor(0.5 + scaleFactor *
                                                    (atom2.x + ptranx));
                                            oy = Math.floor(0.5 + scaleFactor *
                                                    (atom2.y + ptrany));
                                            oz = Math.floor(0.5 + scaleFactor *
                                                    (atom2.z + ptranz));
                                            if (mi * mi + mj * mj + mk * mk < ox *
                                                    ox + oy * oy + oz * oz)
                                                vpAtomID[index] = atom.serial;
                                        }

                                    }// k
                                }// if
                            }// kk
                        }// jj
                    }// ii
                }// if
                nind++;
            }// j
        }// i
    };

    this.fillvoxelswaals = function(atoms, atomlist) {
        var i, il;
        for (i = 0, il = vpBits.length; i < il; i++)
            vpBits[i] &= ~ISDONE; // not isdone

        for (i in atomlist) {
            var atom = atoms[atomlist[i]];
            if (atom === undefined)
                continue;

            this.fillAtomWaals(atom, atoms);
        }
    };

    this.fillAtomWaals = function(atom, atoms) {
        var cx, cy, cz, ox, oy, oz, nind = 0;
        var mi, mj, mk, si, sj, sk, i, j, k, ii, jj, kk, n;
        cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
        cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
        cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

        var at = getVDWIndex(atom);
        var pWH = pWidth*pHeight;
        for (i = 0, n = widxz[at]; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (depty[at][nind] != -1) {
                    for (ii = -1; ii < 2; ii++) {
                        for (jj = -1; jj < 2; jj++) {
                            for (kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    mi = ii * i;
                                    mk = kk * j;
                                    for (k = 0; k <= depty[at][nind]; k++) {
                                        mj = k * jj;
                                        si = cx + mi;
                                        sj = cy + mj;
                                        sk = cz + mk;
                                        if (si < 0 || sj < 0 || 
                                                sk < 0 || 
                                                si >= pLength || 
                                                sj >= pWidth || 
                                                sk >= pHeight)
                                            continue;
                                        var index = si * pWH + sj * pHeight + sk;
                                        if (!(vpBits[index] & ISDONE)) {
                                            vpBits[index] |= ISDONE;
                                            vpAtomID[index] = atom.serial;
                                        } else {
                                            var atom2 = atoms[vpAtomID[index]];
                                            ox = Math.floor(0.5 + scaleFactor * (atom2.x + ptranx));
                                            oy = Math.floor(0.5 + scaleFactor * (atom2.y + ptrany));
                                            oz = Math.floor(0.5 + scaleFactor * (atom2.z + ptranz));
                                            if (mi * mi + mj * mj + mk * mk < ox *
                                                    ox + oy * oy + oz * oz)
                                                vpAtomID[index] = atom.serial;
                                        }

                                    }// k
                                }// if
                            }// kk
                        }// jj
                    }// ii
                }// if
                nind++;
            }// j
        }// i
    };

    this.buildboundary = function() {
        var pWH = pWidth*pHeight;
        for (i = 0; i < pLength; i++) {
            for (j = 0; j < pHeight; j++) {
                for (k = 0; k < pWidth; k++) {
                    var index = i * pWH + k * pHeight + j;
                    if (vpBits[index] & INOUT) {
                        var flagbound = false;
                        var ii = 0;
                        while (ii < 26) {
                            var ti = i + nb[ii][0], tj = j + nb[ii][2], tk = k +
                                    nb[ii][1];
                            if (ti > -1 && 
                                ti < pLength && 
                                tk > -1 && 
                                tk < pWidth && 
                                tj > -1 && 
                                tj < pHeight && 
                                !(vpBits[ti * pWH + tk * pHeight + tj] & INOUT)) {
                                vpBits[index] |= ISBOUND;
                                break;
                            } else
                                ii++;
                        }
                    }
                }
            }
        }
    };

    // a little class for 3d array, should really generalize this and
    // use throughout...
    var PointGrid = function(length, width, height) {
        // the standard says this is zero initialized
        var data = new Int32Array(length * width * height * 3);

        // set position x,y,z to pt, which has ix,iy,and iz
        this.set = function(x, y, z, pt) {
            var index = ((((x * width) + y) * height) + z) * 3;
            data[index] = pt.ix;
            data[index + 1] = pt.iy;
            data[index + 2] = pt.iz;
        };

        // return point at x,y,z
        this.get = function(x, y, z) {
            var index = ((((x * width) + y) * height) + z) * 3;
            return {
                ix : data[index],
                iy : data[index + 1],
                iz : data[index + 2]
            };
        };
    };

    this.fastdistancemap = function() {
        var eliminate = 0;
        var certificate;
        var i, j, k, n;

        var boundPoint = new PointGrid(pLength, pWidth, pHeight);
        var pWH = pWidth*pHeight;
        var cutRSq = cutRadius*cutRadius;
        
        var inarray = [];
        var outarray = [];
        
        var index;
        
        for (i = 0; i < pLength; i++) {
            for (j = 0; j < pWidth; j++) {
                for (k = 0; k < pHeight; k++) {
                    index = i * pWH + j * pHeight + k;
                    vpBits[index] &= ~ISDONE; // isdone = false
                    if (vpBits[index] & INOUT) {
                        if (vpBits[index] & ISBOUND) {
                            var triple = {
                                ix : i,
                                iy : j,
                                iz : k
                            };
                            boundPoint.set(i, j, k, triple);
                            inarray.push(triple);
                            vpDistance[index] = 0;
                            vpBits[index] |= ISDONE;
                            vpBits[index] &= ~ISBOUND;
                        } 
                    }
                }
            }
        }

        do {
            outarray = this.fastoneshell(inarray, boundPoint);
            inarray = [];
            for (i = 0, n = outarray.length; i < n; i++) {
                index = pWH * outarray[i].ix + pHeight *
                    outarray[i].iy + outarray[i].iz;
                vpBits[index] &= ~ISBOUND;
                if (vpDistance[index] <= 1.0404 * cutRSq) {
                    inarray.push({
                        ix : outarray[i].ix,
                        iy : outarray[i].iy,
                        iz : outarray[i].iz
                    });
                }
            }
        } while (inarray.length !== 0);

        inarray = [];
        outarray = [];
        boundPoint = null;
        
        var cutsf = scaleFactor - 0.5;
        if (cutsf < 0)
            cutsf = 0;
        var cutoff = cutRSq - 0.50 / (0.1 + cutsf);
        for (i = 0; i < pLength; i++) {
            for (j = 0; j < pWidth; j++) {
                for (k = 0; k < pHeight; k++) {
                    index = i * pWH + j * pHeight + k;
                    vpBits[index] &= ~ISBOUND;
                    // ses solid
                    if (vpBits[index] & INOUT) {
                        if (!(vpBits[index] & ISDONE) ||
                                ((vpBits[index] & ISDONE) && vpDistance[index] >= cutoff)) {
                            vpBits[index] |= ISBOUND;
                        }
                    }
                }
            }
        }

    };

    this.fastoneshell = function(inarray, boundPoint) { // (int* innum,int
        // *allocout,voxel2
        // ***boundPoint, int*
        // outnum, int *elimi)
        var tx, ty, tz;
        var dx, dy, dz;
        var i, j, n;
        var square;
        var bp, index;
        var outarray = [];
        if (inarray.length === 0)
            return outarray;

        tnv = {
            ix : -1,
            iy : -1,
            iz : -1
        };
        var pWH = pWidth*pHeight;
        for ( i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (j = 0; j < 6; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];
                
                if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth &&
                        tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + pHeight * tnv.iy + tnv.iz;
                    
                    if ((vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
    
                        boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        vpDistance[index] = square;
                        vpBits[index] |= ISDONE;
                        vpBits[index] |= ISBOUND;
    
                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
    
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
    
                            vpDistance[index] = square;
                            if (!(vpBits[index] & ISBOUND)) {
                                vpBits[index] |= ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        // console.log("part1", positout);

        for (i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (j = 6; j < 18; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];

                if(tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth &&
                        tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + pHeight * tnv.iy + tnv.iz;
                    
                    if ((vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
                        boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);
    
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        vpDistance[index] = square;
                        vpBits[index] |= ISDONE;
                        vpBits[index] |= ISBOUND;
    
                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
                            vpDistance[index] = square;
                            if (!(vpBits[index] & ISBOUND)) {
                                vpBits[index] |= ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        // console.log("part2", positout);

        for (i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (j = 18; j < 26; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];

                if (tnv.ix < pLength && tnv.ix > -1 && tnv.iy < pWidth &&
                        tnv.iy > -1 && tnv.iz < pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + pHeight * tnv.iy + tnv.iz;

                    if ((vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
                        boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);

                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        vpDistance[index] = square;
                        vpBits[index] |= ISDONE;
                        vpBits[index] |= ISBOUND;

                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((vpBits[index] & INOUT)  && (vpBits[index] & ISDONE)) {
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

                            vpDistance[index] = square;
                            if (!(vpBits[index] & ISBOUND)) {
                                vpBits[index] |= ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        // console.log("part3", positout);
        return outarray;
    };

    this.marchingcubeinit = function(stype) {
        for ( var i = 0, lim = vpBits.length; i < lim; i++) {
            if (stype == 1) {// vdw
                vpBits[i] &= ~ISBOUND;
            } else if (stype == 4) { // ses
                vpBits[i] &= ~ISDONE;
                if (vpBits[i] & ISBOUND)
                    vpBits[i] |= ISDONE;
                vpBits[i] &= ~ISBOUND;
            } else if (stype == 2) {// after vdw
                if ((vpBits[i] & ISBOUND) && (vpBits[i] & ISDONE))
                    vpBits[i] &= ~ISBOUND;
                else if ((vpBits[i] & ISBOUND) && !(vpBits[i] & ISDONE))
                    vpBits[i] |= ISDONE;
            } else if (stype == 3) { // sas
                vpBits[i] &= ~ISBOUND;
            }
        }
    };

    // this code allows me to empirically prune the marching cubes code tables
    // to more efficiently handle discrete data
    var counter = function() {
        var data = Array(256);
        for ( var i = 0; i < 256; i++)
            data[i] = [];

        this.incrementUsed = function(i, j) {
            if (typeof data[i][j] === 'undefined')
                data[i][j] = {
                    used : 0,
                    unused : 0
                };
            data[i][j].used++;
        };

        this.incrementUnused = function(i, j) {
            if (typeof data[i][j] === 'undefined')
                data[i][j] = {
                    used : 0,
                    unused : 0
                };
            data[i][j].unused++;

        };

        var redoTable = function(triTable) {
            var str = "[";
            for ( var i = 0; i < triTable.length; i++) {
                var code = 0;
                var table = triTable[i];
                for ( var j = 0; j < table.length; j++) {
                    code |= (1 << (table[j]));
                }
                str += "0x" + code.toString(16) + ", ";
            }
            str += "]";
            console.log(str);
        };

        this.print = function() {

            var table = MarchingCube.triTable;
            var str;
            var newtable = [];
            for ( var i = 0; i < table.length; i++) {
                var newarr = [];
                for ( var j = 0; j < table[i].length; j += 3) {
                    var k = j / 3;
                    if (typeof data[i][k] === 'undefined' || !data[i][k].unused) {
                        newarr.push(table[i][j]);
                        newarr.push(table[i][j + 1]);
                        newarr.push(table[i][j + 2]);
                    }
                    if (typeof data[i][k] === 'undefined')
                        console.log("undef " + i + "," + k);
                }
                newtable.push(newarr);
            }
            console.log(JSON.stringify(newtable));
            redoTable(newtable);
        };
    };
    
    this.marchingcube = function(stype) {
        this.marchingcubeinit(stype);
        verts = []; faces = [];   
        WebMol.MarchingCube.march(vpBits, verts, faces, {
            smooth : 1,
            nX : pLength,
            nY : pWidth,
            nZ : pHeight        
        });      


        var pWH = pWidth*pHeight;
        for (var i = 0, vlen = verts.length; i < vlen; i++) {
            verts[i].atomid = vpAtomID[verts[i].x * pWH + pHeight *
                    verts[i].y + verts[i].z];
        }  

        WebMol.MarchingCube.laplacianSmooth(1, verts, faces);

    };


};
/*
* math-like functionality
* quaternion, vector, matrix
*/

//Math functions
var WebMol = WebMol || {};

WebMol.Math = {

    clamp : function(x, min, max) {
        return Math.min( Math.max( x, min ), max );
    },

    degToRad : function() {
       var degreeToRadiansFactor = Math.PI / 180;
       
       return function(deg) {
           return deg * degreeToRadiansFactor;
       };
    
    }()
    
};


// Quaternion

WebMol.Quaternion = function(x, y, z, w) {

    this.x = x || 0;
    this.y = y || 0;
    this.z = z || 0;
    this.w = (w !== undefined) ? w : 1;

};

WebMol.Quaternion.prototype = {

    constructor : WebMol.Quaternion,

    set : function(x, y, z, w) {
        
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;

        return this;
    },

    copy : function(q) {
        
        this.x = q.x;
        this.y = q.y;
        this.z = q.z;
        this.w = q.w;

        return this;
    },

    conjugate : function() {
        
        this.x *= -1;
        this.y *= -1;
        this.z *= -1;

        return this;
    },

    inverse : function() {
        
        return this.conjugate().normalize();
    },

    length : function() {
        
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w);
    },

    normalize : function() {
        
        var l = this.length();

        if (l === 0) {
            this.x = 0;
            this.y = 0;
            this.z = 0;
            this.w = 1;
        } else {
            l = 1 / l;

            this.x *= l;
            this.y *= l;
            this.z *= l;
            this.w *= l;
        }

        return this;

    },

    multiply : function(q) {
        
        return this.multiplyQuaternions(this, q);
    },

    multiplyQuaternions : function(a, b) {

        var qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
        var qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;

        this.x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
        this.y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
        this.z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
        this.w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;

    }
};

//A 2 Vector
WebMol.Vector2 = function(x, y) {
    
    this.x = x || 0.0;
    this.y = y || 0.0;
};

WebMol.Vector2.prototype = {
    
    constructor : WebMol.Vector2,
   
    set : function(x, y) {
       
        this.x = x;
        this.y = y;
       
        return this;
    },
    
    subVectors : function(a, b) {
        
      this.x = a.x - b.x;
      this.y = a.y - b.y;
      
      return this;
    },
   
    copy : function(v) {
       
        this.x = v.x;
        this.y = v.y;
       
        return this;
    },
   
    clone : function() {
        
        return new WebMol.Vector2(this.x, this.y);
    }    
   
};

//A 3 Vector

WebMol.Vector3 = function(x, y, z) {
    this.x = x || 0.0;
    this.y = y || 0.0;
    this.z = z || 0.0;
};

WebMol.Vector3.prototype =  {
    
    constructor : WebMol.Vector3,
    
    set : function(x, y, z) {
        
        this.x = x;
        this.y = y;
        this.z = z;
        
        return this;
    },
    
    copy : function(v) {
        
        this.x = v.x;
        this.y = v.y;
        this.z = v.z;
        
        return this;  
    },
    
    add : function(v) {
        
        this.x += v.x;
        this.y += v.y;
        this.z += v.z;  
        
        return this;
    },
    
    addVectors : function(a, b) {
        
        this.x = a.x + b.x;
        this.y = a.y + b.y;
        this.z = a.z + b.z;
        
        return this;
    },
    
    sub : function(v) {
        
        this.x -= v.x;
        this.y -= v.y;
        this.z -= v.z;
        
        return this;
    },
    
    subVectors : function(a, b) {
        
        this.x = a.x - b.x;
        this.y = a.y - b.y;
        this.z = a.z - b.z;
        
        return this;
    },
    
    multiplyScalar : function(s) {
        
        this.x *= s;
        this.y *= s;
        this.z *= s;
        
        return this;
    },
    
    divideScalar : function(s) {
        
        if (s !== 0) {
            this.x /= s;
            this.y /= s;
            this.z /= s;
        }
        
        else {
            this.x = 0;
            this.y = 0;
            this.z = 0;
        }
        
        return this;
    },
    

    distanceTo: function(v) {
        return Math.sqrt(this.distanceToSquared(v));
    },

    distanceToSquared: function(v) {
        var dx = this.x - v.x;
        var dy = this.y - v.y;
        var dz = this.z - v.z;

        return dx * dx + dy * dy + dz * dz;
    },
    
    applyMatrix4 : function(m) {
    
        var x = this.x, y = this.y, z = this.z;
        
        var e = m.elements;
        
        this.x = e[0]*x + e[4]*y + e[8]*z + e[12];
        this.y = e[1]*x + e[5]*y + e[9]*z + e[13];
        this.z = e[2]*x + e[6]*y + e[10]*z + e[14];
        
        return this;
    },
    
    applyProjection : function(m) {
        
        //input: WebMol.Matrix4 projection matrix
        
        var x = this.x, y = this.y, z = this.z;
        
        var e = m.elements;
        var d = ( e[3]*x + e[7]*y + e[11]*z + e[15]);
        
        this.x = (e[0]*x + e[4]*y + e[8]*z + e[12]) / d;
        this.y = (e[1]*x + e[5]*y + e[9]*z + e[13]) / d;
        this.z = (e[2]*x + e[6]*y + e[10]*z + e[14]) / d;
        
        return this;
    },
    
    applyQuaternion : function(q) { 
        
        var x = this.x;
        var y = this.y;
        var z = this.z;
        
        var qx = q.x;
        var qy = q.y;
        var qz = q.z;
        var qw = q.w;
        
        // calculate quaternion * vector
        
        var ix = qw * x + qy * z - qz * y;
        var iy = qw * y + qz * x - qx * z;
        var iz = qw * z + qx * y - qy * x;
        var iw = -qw * x - qy * y - qz * z;
        
        // calculate result * inverse quaternion
        
        this.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
        this.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
        this.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;
        
        return this;
    },
    
    negate : function() {
        
        return this.multiplyScalar(-1);
    },
    
    dot : function(v) {
        
        return this.x * v.x + this.y * v.y + this.z * v.z;
    },
    
    length : function() {
        
        return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    },
    
    lengthSq : function() {
    
        return (this.x*this.x + this.y*this.y + this.z*this.z);
    },
    
    normalize : function() {
        
        return this.divideScalar( this.length() );
    },
    
    cross : function (v) {
        
        var x = this.x, y = this.y, z = this.z;
        
        this.x = y * v.z - z * v.y;
        this.y = z * v.x - x * v.z;
        this.z = x * v.y - y * v.x;
        
        return this;
    },
    
    crossVectors : function(a, b) {
        
        this.x = a.y * b.z - a.z * b.y;
        this.y = a.z * b.x - a.x * b.z;
        this.z = a.x * b.y - a.y * b.x;
        
        return this;
    },
    
    getPositionFromMatrix : function(m) {
        
        this.x = m.elements[12];
        this.y = m.elements[13];
        this.z = m.elements[14];
        
        return this;
    },

    setEulerFromRotationMatrix : function (m, order) {

        // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

        var te = m.elements;
        var m11 = te[0], m12 = te[4], m13 = te[8];
        var m21 = te[1], m22 = te[5], m23 = te[9];
        var m31 = te[2], m32 = te[6], m33 = te[10];

        if ( order === undefined || order === 'XYZ' ) {

            this.y = Math.asin( WebMol.Math.clamp( m13, -1, 1 ) );

            if ( Math.abs( m13 ) < 0.99999 ) {

                this.x = Math.atan2( - m23, m33 );
                this.z = Math.atan2( - m12, m11 );

            } else {

                this.x = Math.atan2( m32, m22 );
                this.z = 0;

            }
        }
        
        else {
            console.error("Error with vector's setEulerFromRotationMatrix: Unknown order: " + order);
        }
        
        return this;

    },
    
    clone : function() {
        return new WebMol.Vector3(this.x, this.y, this.z);
    }
    
};

//Matrices

//Matrix3

WebMol.Matrix3 = function(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
    
    this.elements = new Float32Array(9);
    
    this.set(
        (n11 !== undefined) ? n11 : 1, n12 || 0, n13 || 0,
        n21 || 0, (n22 !== undefined) ? n22 : 1, n23 || 0,
        n31 || 0, n32 || 0, (n33 !== undefined) ? n33 : 1
    );
    
};

WebMol.Matrix3.prototype = {
    
    constructor : WebMol.Matrix3,    
   
    set : function(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
        var te = this.elements;
        
        te[0] = n11; te[3] = n12; te[6] = n13;
        te[1] = n21; te[4] = n22; te[7] = n23;
        te[2] = n31; te[5] = n32; te[8] = n33;
        
        return this;
    },
    
    identity : function() {   
        this.set(
            1,0,0,
            0,1,0,
            0,0,1
        );
        
        return this;
    },
   
    copy : function(m) {
        var me = m.elements;
       
        this.set(
            me[0], me[3], me[6],
            me[1], me[4], me[7],
            me[2], me[5], me[8]
        );
    },
    
    multiplyScalar: function ( s ) {
        var te = this.elements;

        te[0] *= s; te[3] *= s; te[6] *= s;
        te[1] *= s; te[4] *= s; te[7] *= s;
        te[2] *= s; te[5] *= s; te[8] *= s;

        return this;
    },

    getInverse: function ( matrix, throwOnInvertible ) {
        // input: Matrix4

        var me = matrix.elements;
        var te = this.elements;

        te[ 0 ] =   me[10] * me[5] - me[6] * me[9];
        te[ 1 ] = - me[10] * me[1] + me[2] * me[9];
        te[ 2 ] =   me[6] * me[1] - me[2] * me[5];
        te[ 3 ] = - me[10] * me[4] + me[6] * me[8];
        te[ 4 ] =   me[10] * me[0] - me[2] * me[8];
        te[ 5 ] = - me[6] * me[0] + me[2] * me[4];
        te[ 6 ] =   me[9] * me[4] - me[5] * me[8];
        te[ 7 ] = - me[9] * me[0] + me[1] * me[8];
        te[ 8 ] =   me[5] * me[0] - me[1] * me[4];

        var det = me[ 0 ] * te[ 0 ] + me[ 1 ] * te[ 3 ] + me[ 2 ] * te[ 6 ];

        // no inverse

        if ( det === 0 ) {

            var msg = "Matrix3.getInverse(): can't invert matrix, determinant is 0";

            if ( throwOnInvertible || false ) {

                throw new Error( msg ); 

            } else {

                console.warn( msg );

            }

            this.identity();

            return this;

        }

        this.multiplyScalar( 1.0 / det );

        return this;
    },
    
    transpose: function () {
        var tmp, m = this.elements;

        tmp = m[1]; m[1] = m[3]; m[3] = tmp;
        tmp = m[2]; m[2] = m[6]; m[6] = tmp;
        tmp = m[5]; m[5] = m[7]; m[7] = tmp;

        return this;
    },
    
    clone: function () {
        var te = this.elements;

        return new WebMol.Matrix3(

            te[0], te[3], te[6],
            te[1], te[4], te[7],
            te[2], te[5], te[8]

        );
    }
   
};

//Matrix 4

WebMol.Matrix4 = function(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44) {

    var te = this.elements = new Float32Array( 16 );
    
    te[0] = ( n11 !== undefined ) ? n11 : 1; te[4] = n12 || 0; te[8] = n13 || 0; te[12] = n14 || 0;
    te[1] = n21 || 0; te[5] = ( n22 !== undefined ) ? n22 : 1; te[9] = n23 || 0; te[13] = n24 || 0;
    te[2] = n31 || 0; te[6] = n32 || 0; te[10] = ( n33 !== undefined ) ? n33 : 1; te[14] = n34 || 0;
    te[3] = n41 || 0; te[7] = n42 || 0; te[11] = n43 || 0; te[15] = ( n44 !== undefined ) ? n44 : 1;

};

WebMol.Matrix4.prototype = {

    constructor : WebMol.Matrix4,

    set: function ( n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44 ) {
        var te = this.elements;

        te[0] = n11; te[4] = n12; te[8] = n13; te[12] = n14;
        te[1] = n21; te[5] = n22; te[9] = n23; te[13] = n24;
        te[2] = n31; te[6] = n32; te[10] = n33; te[14] = n34;
        te[3] = n41; te[7] = n42; te[11] = n43; te[15] = n44;

        return this;
    },

    identity: function () {
        this.set(

            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1

        );

        return this;
    },

    copy: function ( m ) {
        var me = m.elements;

        this.set(

            me[0], me[4], me[8], me[12],
            me[1], me[5], me[9], me[13],
            me[2], me[6], me[10], me[14],
            me[3], me[7], me[11], me[15]

        );

        return this;
    },

    setRotationFromEuler: function ( v, order ) {

        var te = this.elements;

        var x = v.x, y = v.y, z = v.z;
        var a = Math.cos( x ), b = Math.sin( x );
        var c = Math.cos( y ), d = Math.sin( y );
        var e = Math.cos( z ), f = Math.sin( z );

        if ( order === undefined || order === 'XYZ' ) {

            var ae = a * e, af = a * f, be = b * e, bf = b * f;

            te[0] = c * e;
            te[4] = - c * f;
            te[8] = d;

            te[1] = af + be * d;
            te[5] = ae - bf * d;
            te[9] = - b * c;

            te[2] = bf - ae * d;
            te[6] = be + af * d;
            te[10] = a * c;

        } 
        
        else
            console.error("Error with matrix4 setRotationFromEuler. Order: " + order);

        return this;

    },

    setRotationFromQuaternion: function ( q ) {
        var te = this.elements;

        var x = q.x, y = q.y, z = q.z, w = q.w;
        var x2 = x + x, y2 = y + y, z2 = z + z;
        var xx = x * x2, xy = x * y2, xz = x * z2;
        var yy = y * y2, yz = y * z2, zz = z * z2;
        var wx = w * x2, wy = w * y2, wz = w * z2;

        te[0] = 1 - ( yy + zz );
        te[4] = xy - wz;
        te[8] = xz + wy;

        te[1] = xy + wz;
        te[5] = 1 - ( xx + zz );
        te[9] = yz - wx;

        te[2] = xz - wy;
        te[6] = yz + wx;
        te[10] = 1 - ( xx + yy );

        return this;
    },

    lookAt: function() {
        var x = new WebMol.Vector3();
        var y = new WebMol.Vector3();
        var z = new WebMol.Vector3();

        return function ( eye, target, up ) {

            var te = this.elements;

            z.subVectors( eye, target ).normalize();

            if ( z.length() === 0 ) {

                z.z = 1;

            }

            x.crossVectors( up, z ).normalize();

            if ( x.length() === 0 ) {

                z.x += 0.0001;
                x.crossVectors( up, z ).normalize();

            }

            y.crossVectors( z, x );


            te[0] = x.x; te[4] = y.x; te[8] = z.x;
            te[1] = x.y; te[5] = y.y; te[9] = z.y;
            te[2] = x.z; te[6] = y.z; te[10] = z.z;

            return this;
        };

    }(),

    multiplyMatrices: function ( a, b ) {
        var ae = a.elements;
        var be = b.elements;
        var te = this.elements;

        var a11 = ae[0], a12 = ae[4], a13 = ae[8], a14 = ae[12];
        var a21 = ae[1], a22 = ae[5], a23 = ae[9], a24 = ae[13];
        var a31 = ae[2], a32 = ae[6], a33 = ae[10], a34 = ae[14];
        var a41 = ae[3], a42 = ae[7], a43 = ae[11], a44 = ae[15];

        var b11 = be[0], b12 = be[4], b13 = be[8], b14 = be[12];
        var b21 = be[1], b22 = be[5], b23 = be[9], b24 = be[13];
        var b31 = be[2], b32 = be[6], b33 = be[10], b34 = be[14];
        var b41 = be[3], b42 = be[7], b43 = be[11], b44 = be[15];

        te[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
        te[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
        te[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
        te[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

        te[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
        te[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
        te[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
        te[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

        te[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
        te[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
        te[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
        te[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

        te[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
        te[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
        te[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
        te[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

        return this;
    },
    
    multiplyScalar: function ( s ) {
        var te = this.elements;
    
        te[0] *= s; te[4] *= s; te[8] *= s; te[12] *= s;
        te[1] *= s; te[5] *= s; te[9] *= s; te[13] *= s;
        te[2] *= s; te[6] *= s; te[10] *= s; te[14] *= s;
        te[3] *= s; te[7] *= s; te[11] *= s; te[15] *= s;
    
        return this;
    },
    
    transpose: function () {
        var te = this.elements;
        var tmp;

        tmp = te[1]; te[1] = te[4]; te[4] = tmp;
        tmp = te[2]; te[2] = te[8]; te[8] = tmp;
        tmp = te[6]; te[6] = te[9]; te[9] = tmp;

        tmp = te[3]; te[3] = te[12]; te[12] = tmp;
        tmp = te[7]; te[7] = te[13]; te[13] = tmp;
        tmp = te[11]; te[11] = te[14]; te[14] = tmp;

        return this;
    },

    getPosition: function() {
        var v1 = new WebMol.Vector3();

        return function () {

            console.warn( 'DEPRECATED: Matrix4\'s .getPosition() has been removed. Use Vector3.getPositionFromMatrix( matrix ) instead.' );

            var te = this.elements;
            return v1.set( te[12], te[13], te[14] );
        };

    }(),

    setPosition: function ( v ) {
        var te = this.elements;

        te[12] = v.x;
        te[13] = v.y;
        te[14] = v.z;

        return this;
    },

    getInverse: function ( m, throwOnInvertible ) {
        // based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
        var te = this.elements;
        var me = m.elements;

        var n11 = me[0], n12 = me[4], n13 = me[8], n14 = me[12];
        var n21 = me[1], n22 = me[5], n23 = me[9], n24 = me[13];
        var n31 = me[2], n32 = me[6], n33 = me[10], n34 = me[14];
        var n41 = me[3], n42 = me[7], n43 = me[11], n44 = me[15];

        te[0] = n23*n34*n42 - n24*n33*n42 + n24*n32*n43 - n22*n34*n43 - n23*n32*n44 + n22*n33*n44;
        te[4] = n14*n33*n42 - n13*n34*n42 - n14*n32*n43 + n12*n34*n43 + n13*n32*n44 - n12*n33*n44;
        te[8] = n13*n24*n42 - n14*n23*n42 + n14*n22*n43 - n12*n24*n43 - n13*n22*n44 + n12*n23*n44;
        te[12] = n14*n23*n32 - n13*n24*n32 - n14*n22*n33 + n12*n24*n33 + n13*n22*n34 - n12*n23*n34;
        te[1] = n24*n33*n41 - n23*n34*n41 - n24*n31*n43 + n21*n34*n43 + n23*n31*n44 - n21*n33*n44;
        te[5] = n13*n34*n41 - n14*n33*n41 + n14*n31*n43 - n11*n34*n43 - n13*n31*n44 + n11*n33*n44;
        te[9] = n14*n23*n41 - n13*n24*n41 - n14*n21*n43 + n11*n24*n43 + n13*n21*n44 - n11*n23*n44;
        te[13] = n13*n24*n31 - n14*n23*n31 + n14*n21*n33 - n11*n24*n33 - n13*n21*n34 + n11*n23*n34;
        te[2] = n22*n34*n41 - n24*n32*n41 + n24*n31*n42 - n21*n34*n42 - n22*n31*n44 + n21*n32*n44;
        te[6] = n14*n32*n41 - n12*n34*n41 - n14*n31*n42 + n11*n34*n42 + n12*n31*n44 - n11*n32*n44;
        te[10] = n12*n24*n41 - n14*n22*n41 + n14*n21*n42 - n11*n24*n42 - n12*n21*n44 + n11*n22*n44;
        te[14] = n14*n22*n31 - n12*n24*n31 - n14*n21*n32 + n11*n24*n32 + n12*n21*n34 - n11*n22*n34;
        te[3] = n23*n32*n41 - n22*n33*n41 - n23*n31*n42 + n21*n33*n42 + n22*n31*n43 - n21*n32*n43;
        te[7] = n12*n33*n41 - n13*n32*n41 + n13*n31*n42 - n11*n33*n42 - n12*n31*n43 + n11*n32*n43;
        te[11] = n13*n22*n41 - n12*n23*n41 - n13*n21*n42 + n11*n23*n42 + n12*n21*n43 - n11*n22*n43;
        te[15] = n12*n23*n31 - n13*n22*n31 + n13*n21*n32 - n11*n23*n32 - n12*n21*n33 + n11*n22*n33;

        var det = me[ 0 ] * te[ 0 ] + me[ 1 ] * te[ 4 ] + me[ 2 ] * te[ 8 ] + me[ 3 ] * te[ 12 ];

        if ( det === 0 ) {

            var msg = "Matrix4.getInverse(): can't invert matrix, determinant is 0";

            if ( throwOnInvertible || false ) {

                throw new Error( msg ); 

            } else {

                console.warn( msg );

            }

            this.identity();

            return this;
        }

        this.multiplyScalar( 1 / det );

        return this;
    },

    compose: function() {
        var mRotation = new WebMol.Matrix4(),
            mScale = new WebMol.Matrix4();
        
        return function ( translation, rotation, scale ) {

            var te = this.elements;

            mRotation.identity();
            mRotation.setRotationFromQuaternion( rotation );

            mScale.makeScale( scale.x, scale.y, scale.z );

            this.multiplyMatrices( mRotation, mScale );

            te[12] = translation.x;
            te[13] = translation.y;
            te[14] = translation.z;

            return this;

        };
    }(),

    decompose: function() {
        var x = new WebMol.Vector3(),
            y = new WebMol.Vector3(),
            z = new WebMol.Vector3(),
            matrix = new WebMol.Matrix4();

        return function ( translation, rotation, scale ) {

            var te = this.elements;

            // grab the axis vectors
            x.set( te[0], te[1], te[2] );
            y.set( te[4], te[5], te[6] );
            z.set( te[8], te[9], te[10] );

            translation = ( translation instanceof WebMol.Vector3 ) ? translation : new WebMol.Vector3();
            rotation = ( rotation instanceof WebMol.Quaternion ) ? rotation : new WebMol.Quaternion();
            scale = ( scale instanceof Webmol.Vector3 ) ? scale : new WebMol.Vector3();

            scale.x = x.length();
            scale.y = y.length();
            scale.z = z.length();

            translation.x = te[12];
            translation.y = te[13];
            translation.z = te[14];

            // scale the rotation part

            matrix.copy( this );

            matrix.elements[0] /= scale.x;
            matrix.elements[1] /= scale.x;
            matrix.elements[2] /= scale.x;

            matrix.elements[4] /= scale.y;
            matrix.elements[5] /= scale.y;
            matrix.elements[6] /= scale.y;

            matrix.elements[8] /= scale.z;
            matrix.elements[9] /= scale.z;
            matrix.elements[10] /= scale.z;

            rotation.setFromRotationMatrix( matrix );

            return [ translation, rotation, scale ];

        };
    }(),

    scale: function ( v ) {
        var te = this.elements;
        var x = v.x, y = v.y, z = v.z;

        te[0] *= x; te[4] *= y; te[8] *= z;
        te[1] *= x; te[5] *= y; te[9] *= z;
        te[2] *= x; te[6] *= y; te[10] *= z;
        te[3] *= x; te[7] *= y; te[11] *= z;

        return this;
    },
    
    getMaxScaleOnAxis : function() {
        
        var te = this.elements;
        
        var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
        var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
        var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];
        
        return Math.sqrt(Math.max(scaleXSq, Math.max(scaleYSq, scaleZSq)));
        
    },

    makeFrustum: function ( left, right, bottom, top, near, far ) {
        var te = this.elements;
        var x = 2 * near / ( right - left );
        var y = 2 * near / ( top - bottom );

        var a = ( right + left ) / ( right - left );
        var b = ( top + bottom ) / ( top - bottom );
        var c = - ( far + near ) / ( far - near );
        var d = - 2 * far * near / ( far - near );

        te[0] = x;  te[4] = 0;  te[8] = a;  te[12] = 0;
        te[1] = 0;  te[5] = y;  te[9] = b;  te[13] = 0;
        te[2] = 0;  te[6] = 0;  te[10] = c; te[14] = d;
        te[3] = 0;  te[7] = 0;  te[11] = - 1;   te[15] = 0;

        return this;
    },

    makePerspective: function ( fov, aspect, near, far ) {
        var ymax = near * Math.tan( WebMol.Math.degToRad( fov * 0.5 ) );
        var ymin = - ymax;
        var xmin = ymin * aspect;
        var xmax = ymax * aspect;

        return this.makeFrustum( xmin, xmax, ymin, ymax, near, far );
    },
    

    clone: function () {
        var te = this.elements;

        return new WebMol.Matrix4(

            te[0], te[4], te[8], te[12],
            te[1], te[5], te[9], te[13],
            te[2], te[6], te[10], te[14],
            te[3], te[7], te[11], te[15]

        );
    }
    
};

WebMol.Ray = function(origin, direction) {
    
    this.origin = (origin !== undefined) ? 
        origin : new WebMol.Vector3();
        
    this.direction = (direction !== undefined) ?
        direction : new WebMol.Vector3();
      
};

//TODO: Remove methods we don't need (intersectPlane ??)
WebMol.Ray.prototype = {
    
    constructor : WebMol.Ray,
     
    set : function(origin, direction){
        
        this.origin.copy(origin);
        this.direction.copy(direction);
        
        return this;
    
    },
    
    copy : function(ray) {
        
        this.origin.copy(ray.origin);
        this.direction.copy(ray.direction);
        
        return this;
        
    },
    
    at : function(t, optionalTarget) {
        
        var result = optionalTarget || new WebMol.Vector3();
        
        return result.copy(this.direction).multiplyScalar(t).add(this.origin);
        
    },
    
    recast : function() {
        
        var v1 = new WebMol.Vector3();
        
        return function(t) {
            this.origin.copy(this.at(t, v1));
            
            return this;
        };
        
    }(),
    
    closestPointToPoint : function(point, optionalTarget) {
        
        var result = optionalTarget || new WebMol.Vector3();
        result.subVectors(point, this.origin);
        var directionDistance = result.dot(this.direction);
        
        //returns a point on this ray
        return result.copy(this.direction).multiplyScalar(directionDistance).add(this.origin);
        
    },
    
    distanceToPoint : function() {
        
        var v1 = new WebMol.Vector3();
        
        return function(point) {
            var directionDistance = v1.subVectors(point, this.origin).dot(this.direction);
            v1.copy(this.direction).multiplyScalar(directionDistance).add(this.origin);
            
            return v1.distanceTo(point);
        };
        
    }(),
    
    isIntersectionCylinder : function() {
        
    },
    
    isIntersectionSphere : function(sphere) {
       
       return (this.distanceToPoint(sphere.center) <= sphere.radius);
          
    },
    
    isIntersectionPlane : function(plane) {
        
        var denominator = plane.normal.dot(this.direction);
        
        //plane and ray are not perpendicular
        if (denominator !== 0) 
            return true;
        
        if (plane.distanceToPoint(this.origin) === 0) 
            return true;
        
        return false;
        
    },
    
    distanceToPlane : function(plane) {
       
       var denominator = plane.normal.dot(this.direction);
       if (denominator === 0) {
           
           //line is coplanar
       if (plane.distanceToPoint(this.origin) === 0)
           return 0;
       
       //ray is parallel
           return undefined;
       }
       
       var t = - (this.origin.dot(plane.normal) + plane.constant) / denominator;
       
       return t;
       
    },
    
    intersectPlane : function(plane, optionalTarget) {
       
       var t = this.distanceToPlane(plane);
       
       if (t === undefined)
           return undefined;
       
       return this.at(t, optionalTarget);
       
    },
    
    
    
    applyMatrix4 : function(matrix4) {
       
       this.direction.add(this.origin).applyMatrix4(matrix4);
       this.origin.applyMatrix4(matrix4);
       this.direction.sub(this.origin);
       
       return this;
       
    },
    
    equals : function(ray) {
       
       return ray.origin.equals(this.origin) && ray.direction.equals(this.direction);
       
    },
    
    clone : function() {
    
       return new WebMol.Ray().copy(this);
    
    }
 
     
};

//Intersection sphere and box shapes.  


//Intersection sphere for sphere, stick render
WebMol.Sphere = function(center, radius) {

    this.center = (center !== undefined) ? 
        center : new WebMol.Vector3();
        
    this.radius = (radius !== undefined) ?
        radius : 0;
        
};

WebMol.Sphere.prototype = {
    
    constructor : WebMol.Sphere,
    
    set : function(center, radius) {
        
        this.center.copy(center);
        this.radius = radius;
        
        return this;
        
    },
    
    copy : function(sphere) {
        
        this.center.copy(sphere.center);
        this.radius = sphere.radius;
        
        return this;
        
    },
    
    applyMatrix4 : function(matrix) {
        
        this.center.applyMatrix4(matrix);
        this.radius = this.radius * matrix.getMaxScaleOnAxis();
        
        return this;
        
    },
    
    translate : function(offset) {
        
        this.center.add(offset);
        
        return this;
        
    },
    
    equals : function(sphere) {
        
        return sphere.center.equals(this.center) && (sphere.radius === this.radius);
        
    },
       
    clone : function() {
        
        return new WebMol.Sphere().copy(this);
        
    }

};


//Bounding cylinder for stick render  
WebMol.Cylinder = function(c1, c2, radius) {

    this.c1 = (c1 !== undefined) ?
        c1 : new WebMol.Vector3();

    this.c2 = (c2 !== undefined) ?
        c2 : new WebMol.Vector3();
        
    this.direction = new WebMol.Vector3().subVectors(this.c2, this.c1).normalize();

    this.radius = (radius !== undefined) ?
        radius : 0;
    
};

WebMol.Cylinder.prototype = {

    constructor : WebMol.Cylinder,

    copy : function(cylinder) {

        this.c1.copy(cylinder.c1);
        this.c2.copy(cylinder.c2);
        this.direction.copy(cylinder.direction);
        this.radius = cylinder.radius;

        return this;

    },
    
    lengthSq : function() {
    
        var vector = new WebMol.Vector3();
        
        return function(){
            return vector.subVectors(this.c2, this.c1).lengthSq();
        };
        
    }(),

    applyMatrix4 : function(matrix) {
        
        this.direction.add(this.c1).applyMatrix4(matrix);
        this.c1.applyMatrix4(matrix);
        this.c2.applyMatrix4(matrix);
        this.direction.sub(this.c1).normalize();
        this.radius = this.radius * matrix.getMaxScaleOnAxis();

        return this;

    }

};


//plane specified by three points
WebMol.Triangle = function(a, b, c){
   
    this.a = (a !== undefined) ?
        a : new WebMol.Vector3();

    this.b = (b !== undefined) ?
        b : new WebMol.Vector3();
    
    this.c = (c !== undefined) ?
        c : new WebMol.Vector3();   
  
};

WebMol.Triangle.prototype = {

    constructor : WebMol.Triangle,
    
    copy : function(triangle) {
        
        this.a.copy(triangle.a);
        this.b.copy(triangle.b);
        this.c.copy(triangle.c);
        
        return this;
        
    },
    
    applyMatrix4 : function(matrix) {
        
        this.a.applyMatrix4(matrix);
        this.b.applyMatrix4(matrix);
        this.c.applyMatrix4(matrix);
        
        return this;
        
    },
    
    getNormal : function() {
        
        var v1 = new WebMol.Vector3();
        
        return function() {
            
            var norm = this.a.clone();
            norm.sub(this.b);
            v1.subVectors(this.c, this.b);
            
            norm.cross(v1);
            norm.normalize();
            
            return norm;
            
        };
        
    }()

};


/* core Object3D
 * Base class for Scene, Camera, Geometry
 * Geometry class
 */

var WebMol = WebMol || {};

//Event Handling
WebMol.EventDispatcher = function() {
  
    var listeners = {};
    
    this.addEventListener = function(type, listener) {
        if (listeners[type] === undefined)
            listeners[type] = [];
        
        if (listeners[type].indexOf(listener) === -1)
            listeners[type].push(listener);
    };  
    
    this.removeEventListener = function(type, listener) {
        
        var index = listeners[type].indexOf(listener);
        
        if (index !== -1)
            listeners[type].splice(index, 1);
              
    };
    
    this.dispatchEvent = function(event) {
        
        var listenerArray = listeners[event.type];
        
        if (listenerArray !== undefined) {
            event.target = this;
            
            for (var i = 0, l = listenerArray.length; i < l; i++)
                listenerArray[i].call(this, event);
                
        }
            
    };
    
};


//Object3D base constructor function
WebMol.Object3D = function() {
    
    this.id = WebMol.Object3DIDCount++;
    
    this.name = "";
    
    this.parent = undefined;
    this.children = [];
    
    this.position = new WebMol.Vector3();
    this.rotation = new WebMol.Vector3();
    this.matrix = new WebMol.Matrix4();
    this.matrixWorld = new WebMol.Matrix4();
    this.quaternion = new WebMol.Quaternion();
    this.eulerOrder = 'XYZ';
    
    this.up = new WebMol.Vector3(0, 1, 0);
    this.scale = new WebMol.Vector3(1, 1, 1);
    
    this.matrixAutoUpdate = true;
    this.matrixWorldNeedsUpdate = true;
    this.rotationAutoUpdate = true;
    this.useQuaternion = false;
    
    this.visible = true;
    
};

WebMol.Object3D.prototype = {
    
    constructor : WebMol.Object3D,
    
    lookAt : function(vector) {
        
        this.matrix.lookAt(vector, this.position, this.up);
        
        if (this.rotationAutoUpdate) {
            
            if (this.useQuaternion === true) 
                this.quaternion.copy(this.matrix.decompose()[1]);
            else
                this.rotation.setEulerFromRotationMatrix(this.matrix, this.eulerOrder);
        }
    },
    
    //add child object
    add : function(object) {
        if (object === this){
            console.error("Can't add WebMol.Object3D to itself");
            return;
        }
        
        object.parent = this;
        this.children.push(object);
        
        //add to the scene (i.e. follow up this instance's parents until reach the top)
        
        var scene = this;
        
        while (scene.parent !== undefined)
            scene = scene.parent;
            
        if (scene !== undefined && scene instanceof WebMol.Scene) 
            scene.__addObject(object);
        
    },
    
    remove : function(object) {
        
        var index = this.children.indexOf(object);
        
        if (index !== -1) {
            
            object.parent = undefined;
            this.children.splice(index, 1);
            
            //Remove from scene
            
            var scene = this;
            
            while (scene.parent !== undefined)
                scene = scene.parent;
                
            if (scene !== undefined && scene instanceof WebMol.Scene)
                scene.__removeObject(object);
                
        }
    },
    
    updateMatrix : function() {
        
        this.matrix.setPosition(this.position);
        
        if (this.useQuaternion === false) 
            this.matrix.setRotationFromEuler(this.rotation, this.eulerOrder);
        else
            this.matrix.setRotationFromQuaternion(this.quaternion);
        
        //TODO: Do I need this??
        if (this.scale.x !== 1 || this.scale.y !== 1 || this.scale.z !== 1)
            this.matrix.scale(this.scale);
            
        this.matrixWorldNeedsUpdate = true;
        
    },
    
    updateMatrixWorld : function(force) {
        
        if (this.matrixAutoUpdate === true) 
            this.updateMatrix();
        
        if (this.matrixWorldNeedsUpdate === true || force === true) {
            
            if (this.parent === undefined)
                this.matrixWorld.copy(this.matrix);
            else
                this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
                
        }
        
        this.matrixWorldNeedsUpdate = false;
        
        //Update matrices of all children
        for (var i in this.children) {
            this.children[i].updateMatrixWorld(true);
        }
    },
    
    clone : function(object) {
        
        if (object === undefined)
            object = new WebMol.Object3D();
            
        object.name = this.name;
        
        object.up.copy(this.up);
        object.position.copy(this.position);
        object.rotation.copy(this.rotation);
        object.eulerOrder = this.eulerOrder;
        object.scale.copy(this.scale);
        
        object.rotationAutoUpdate = this.rotationAutoUpdate;
        object.matrix.copy(this.matrix);
        object.matrixWorld.copy(this.matrixWorld);
        object.quaternion.copy(this.quaternion);
        
        object.useQuaternion = this.useQuaternion;
        
        object.visible = this.visible;
        
        for (var i in this.children) {
            var child = this.children[i];
            object.add(child.clone());
        }
        
        return object;
        
    }
    
};

WebMol.Object3DIDCount = 0;

//Geometry class
//TODO: What can I remove - how can I optimize ?
WebMol.Geometry = (function() {

    //return truncated typed array, including its buffer
    // type == 0 => Uint16Array; type == 1 => Float32Array
    //TODO: Should integrate this directly into geometryGroup's truncateArrayBuffers method
    var truncateArrayBuffer = function(arr, type, end) {
        
        if (arr === null || arr === undefined) {
            return (type === 0) ? new Uint16Array() : new Float32Array();
        }
        
        if (type === 0)
            return new Uint16Array(arr.buffer.slice(arr.byteOffset, end*2));
        else if (type === 1) 
            return new Float32Array(arr.buffer.slice(arr.byteOffset, end*4));
    };
    
    
    var geometryGroup = function(id) {
        this.id = id || 0;
        this.__vertexArray = null;
        this.__colorArray = null;
        this.__normalArray = null;
        this.__faceArray = null;
        this.__lineArray = null;
        this.vertices = 0;
        this.faceidx = 0;
        this.lineidx = 0;
    };
    
    geometryGroup.prototype.getCentroid = function() {
        
        var centroid = new WebMol.Vector3();
        var offset, x, y, z;
        
        for (var i = 0; i < this.vertices; ++i) {
            offset = i*3;
            
            x = this.__vertexArray[offset]; y = this.__vertexArray[offset+1]; z = this.__vertexArray[offset+2];
            
            centroid.x += x; centroid.y += y; centroid.z += z;
        }
        
        //divideScalar checks for 0 denom
        centroid.divideScalar(this.vertices);
        
        return centroid;
    };
    
    //setup normals - vertex and face array must exist
    geometryGroup.prototype.setNormals = function() {        
        
        var faces = this.__faceArray;
        var verts = this.__vertexArray;
        var norms = this.__normalArray;
        
        if (! this.vertices || ! this.faceidx) 
            return;
        
        //vertex indices
        var a, b, c, d,
        //and actual vertices
        vA, vB, vC, norm;
            
        for (var i = 0; i < faces.length / 3; ++i) {
            
            a = faces[i * 3] * 3;
            b = faces[i * 3 + 1] * 3;
            c = faces[i * 3 + 2] * 3;
            
            vA = new WebMol.Vector3(verts[a], verts[a+1], verts[a+2]);
            vB = new WebMol.Vector3(verts[b], verts[b+1], verts[b+2]);
            vC = new WebMol.Vector3(verts[c], verts[c+1], verts[c+2]);
            
            vA.subVectors(vA, vB);
            vC.subVectors(vC, vB);
            vC.cross(vA);
            
            //face normal
            norm = vC;
            norm.normalize();
            
            norms[a] += norm.x; norms[b] += norm.x; norms[c] += norm.x;
            norms[a + 1] += norm.y; norms[b + 1] += norm.y; norms[c + 1] += norm.y;
            norms[a + 2] += norm.z; norms[b + 2] += norm.z; norms[c + 2] += norm.z;
            
        }             
                
    };
    
    //sets line index array from face arr
    //Note - assumes all faces are triangles (i.e. there will
    //be an extra diagonal for four-sided faces - user should 
    //specify linearr for custom shape generation to show wireframe squares
    //as rectangles rather than two triangles)
    geometryGroup.prototype.setLineIndices = function() {
        
        if (! this.faceidx)
            return;
                    
        var faceArr = this.__faceArray, lineArr = this.__lineArray = new Uint16Array(this.faceidx*2);      
        this.lineidx = this.faceidx*2;         
        var faceoffset;
            
        for (var i = 0; i < this.faceidx / 3; ++i) {
            
            faceoffset = i*3; lineoffset = faceoffset*2;          
            var a = faceArr[faceoffset], b = faceArr[faceoffset+1], c = faceArr[faceoffset+2];
            
            lineArr[lineoffset] = a; lineArr[lineoffset+1] = b;
            lineArr[lineoffset+2] = a; lineArr[lineoffset+3] = c;
            lineArr[lineoffset+4] = b; lineArr[lineoffset+5] = c;
            
        }
    };
    
    geometryGroup.prototype.truncateArrayBuffers = function(mesh) {
        
        mesh = (mesh === true) ? true : false;
        
        var vertexArr = this.__vertexArray,
            colorArr = this.__colorArray,
            normalArr = this.__normalArray,
            faceArr = this.__faceArray,
            lineArr = this.__lineArray;
                       
        this.__vertexArray = truncateArrayBuffer(vertexArr, 1, this.vertices*3);
        this.__colorArray = truncateArrayBuffer(colorArr, 1, this.vertices*3);
        
        if (mesh) {
            this.__normalArray = truncateArrayBuffer(normalArr, 1, this.vertices*3);
            this.__faceArray = truncateArrayBuffer(faceArr, 0, this.faceidx);
            this.__lineArray = truncateArrayBuffer(lineArr, 0, this.lineidx);
        }
        else {
            this.__normalArray = truncateArrayBuffer(normalArr, 1, 0);
            this.__faceArray = truncateArrayBuffer(faceArr, 0, 0);
            this.__lineArray = truncateArrayBuffer(lineArr, 0, 0);            
        }
        
        this.__inittedArrays = true;        
        
    };
    
    var addGroup = function(geo) {
        var ret = new geometryGroup(geo.geometryGroups.length);
        geo.geometryGroups.push(ret);
        geo.groups = geo.geometryGroups.length;
        
        ret.__vertexArray = new Float32Array(65535*3);
        ret.__colorArray = new Float32Array(65535*3);
        
        //TODO: instantiating uint arrays according to max number of vertices
        // is dangerous, since there exists the possibility that there will be 
        // more face or line indices than vertex points - but so far that doesn't
        // seem to be the case for any of the renders 
        if (geo.mesh) {
            ret.__normalArray = new Float32Array(65535*3);
            ret.__faceArray = new Uint16Array(65535*6);
            ret.__lineArray = new Uint16Array(65535*6);
        }
        
        
        return ret;
    };
        
    Geometry = function(mesh) {
        
        WebMol.EventDispatcher.call(this);
        
        this.id = WebMol.GeometryIDCount++;
    
        this.name = '';
    
        this.hasTangents = false;
    
        this.dynamic = true; // the intermediate typed arrays will be deleted when set to false
        this.mesh = (mesh === true) ? true : false; // Does this geometry represent a mesh (i.e. do we need Face/Line index buffers?)
        // update flags
    
        this.verticesNeedUpdate = false;
        this.elementsNeedUpdate = false;
        this.normalsNeedUpdate = false;
        this.colorsNeedUpdate = false;
    
        this.buffersNeedUpdate = false;
        
        this.geometryGroups = [];
        this.groups = 0;
        
    };
    
    
    Geometry.prototype = {
        
        constructor : Geometry,

        //Get geometry group to accomodate addVertices new vertices - create 
        // new group if necessary       
        updateGeoGroup : function(addVertices) {
        
            addVertices = addVertices || 0;
            
            var retGroup = this.groups > 0 ? this.geometryGroups[ this.groups - 1 ] : null;
            
            if (!retGroup || retGroup.vertices + addVertices > 65535) 
                retGroup = addGroup(this);
                
            return retGroup;
            
        },
        
        addGeoGroup : function() {
            return addGroup(this);  
        },
        
        setUpNormals : function(three) {
            
            three = three || false;
            
            for ( var g in this.geometryGroups ) {
            
                var geoGroup = this.geometryGroups[g];            
                
                geoGroup.setNormals(three);
                
            }  
                      
        },
        
        setUpWireframe : function() {
            for (var g in this.geometryGroups ) {
                var geoGroup = this.geometryGroups[g];
                
                geoGroup.setLineIndices();
            }
        },
        
        //After vertices, colors, etc are collected in regular or typed arrays,
        // either create typed arrays from regular arrays if they don't already exist,
        // or shorten last typed array
        initTypedArrays : function() {
                
            for (var g in this.geometryGroups) {
                
                var group = this.geometryGroups[g];
                
                if (group.__inittedArrays === true)
                    continue;
                
                group.truncateArrayBuffers(this.mesh);
            }
            
        
        },
        
        dispose : function() {
            this.dispatchEvent( {type : 'dispose'} );
        }
    };

    
    return Geometry;
    
})();

Object.defineProperty(WebMol.Geometry.prototype, "vertices", {
    
    get : function() {
        var vertices = 0;
        for (var g in this.geometryGroups)
            vertices += this.geometryGroups[g].vertices;
            
        return vertices;
    } 
        
});

WebMol.GeometryIDCount = 0;


//Raycaster

WebMol.Raycaster = (function() {
    
    var Raycaster = function(origin, direction, far, near) {
        
        this.ray = new WebMol.Ray(origin, direction);
        
        if (this.ray.direction.lengthSq() > 0) 
            this.ray.direction.normalize();
        
        this.near = near || 0;
        this.far = far || Infinity;
    
    };
    
    var sphere = new WebMol.Sphere();
    var cylinder = new WebMol.Cylinder();
    var triangle = new WebMol.Triangle();
    var w_0 = new WebMol.Vector3(); // for cylinders, cylinder.c1 - ray.origin
    var v1 = new WebMol.Vector3(); // all purpose local vector
    var v2 = new WebMol.Vector3();
    var v3 = new WebMol.Vector3();
    //var facePlane = new WebMol.Plane();
    var localRay = new WebMol.Ray();
    var intersectPoint = new WebMol.Vector3();
    var matrixPosition = new WebMol.Vector3();
    
    var inverseMatrix = new WebMol.Matrix4();
        
    var descSort = function(a, b) {
        return a.distance - b.distance;
    };

    // [-1, 1]
    var clamp = function(x) {
        return Math.min(Math.max(x, -1), 1);
    };
    
    //object is a Sphere or (Bounding) Box
    var intersectObject = function(group, clickable, raycaster, intersects) {
        
        matrixPosition.getPositionFromMatrix(group.matrixWorld);
        
        if ((clickable.clickable !== true) || (clickable.intersectionShape === undefined))
            return intersects;
        
        var intersectionShape = clickable.intersectionShape;
        var precision = raycaster.linePrecision;
        precision *= group.matrixWorld.getMaxScaleOnAxis();
        var precisionSq = precision*precision;

        //Check for intersection with clickable's bounding sphere, if it exists
        if (clickable.boundingSphere !== undefined && clickable.boundingSphere instanceof WebMol.Sphere) {
            sphere.copy(clickable.boundingSphere);
            sphere.applyMatrix4(group.matrixWorld);
            
            if (!raycaster.ray.isIntersectionSphere(sphere)) {
                return intersects;
            }
        }
        
        //Iterate through intersection objects
        var i, il,
            norm, normProj, cylProj, rayProj,
            distance, closestDistSq, denom, discriminant,
            s, t, s_c, t_c;
        //triangle faces
        for (i = 0, il = intersectionShape.triangle.length; i < il; i++) {
            
            if (intersectionShape.triangle[i] instanceof WebMol.Triangle) {
                
                triangle.copy(intersectionShape.triangle[i]);
                triangle.applyMatrix4(group.matrixWorld);
                
                norm = triangle.getNormal();
                
                normProj = raycaster.ray.direction.dot(norm);
                
                //face culling
                if (normProj >= 0)
                    continue;
                
                w_0.subVectors(triangle.a, raycaster.ray.origin);
                
                distance = (norm.dot(w_0)) / normProj;
                
                if (distance < 0)
                    continue;
                    
                //intersects with plane, check if P inside triangle
                v1.copy(raycaster.ray.direction).multiplyScalar(distance).add(raycaster.ray.origin);
                v1.sub(triangle.a); // from pt a to intersection point P
                
                v2.copy(triangle.b).sub(triangle.a); // from pt a to b
                v3.copy(triangle.c).sub(triangle.a); // from pt a to c
                var b_dot_c = v2.dot(v3);
                var b_sq = v2.lengthSq();
                var c_sq = v3.lengthSq();
                
                // P = A + s(v2) + t(v3), inside trianle if 0 <= s, t <=1  and (s + t) <=0
                
                t = ( b_sq*v1.dot(v3) - b_dot_c*v1.dot(v2) ) / ( b_sq*c_sq - b_dot_c*b_dot_c );
                
                if (t < 0 || t > 1)
                    continue;
                
                s = ( v1.dot(v2) - t*b_dot_c ) / b_sq;
                
                if ( (s < 0 || s > 1) || s + t > 1)
                    continue;
                    
                else
                    intersects.push({clickable : clickable,
                                     distance : distance});  
            }
        }
        
        //cylinders
        for (i = 0, il = intersectionShape.cylinder.length; i < il; i++) {
            
            if (intersectionShape.cylinder[i] instanceof WebMol.Cylinder){
                
                cylinder.copy(intersectionShape.cylinder[i]);
                cylinder.applyMatrix4(group.matrixWorld);
                
                w_0.subVectors(cylinder.c1, raycaster.ray.origin); 
                
                cylProj = w_0.dot(cylinder.direction); // Dela
                rayProj = w_0.dot(raycaster.ray.direction); // Epsilon
                
                normProj = clamp(raycaster.ray.direction.dot(cylinder.direction)); // Beta
                
                denom = 1 - normProj*normProj;
                
                if (denom === 0.0)
                    continue;
                
                s_c = (normProj*rayProj - cylProj) / denom;
                t_c = (rayProj - normProj*cylProj) / denom;
                
                v1.copy(cylinder.direction).multiplyScalar(s_c).add(cylinder.c1);  // Q_c
                v2.copy(raycaster.ray.direction).multiplyScalar(t_c).add(raycaster.ray.origin); // P_c
                
                closestDistSq = v3.subVectors(v1, v2).lengthSq();
                var radiusSq = cylinder.radius*cylinder.radius;
                
                //Smoothing?
                //if (closestDistSq > radiusSq) radiusSq += precisionSq;
                
                // closest distance between ray and cylinder axis not greater than cylinder radius;
                // might intersect this cylinder between atom and bond midpoint
                if (closestDistSq <= radiusSq){

                    //Find points where ray intersects sides of cylinder
                    discriminant = (normProj*cylProj - rayProj)*(normProj*cylProj - rayProj) - 
                            denom*(w_0.lengthSq() - cylProj*cylProj - radiusSq);
                    
                    // ray tangent to cylinder?
                    if (discriminant <= 0)
                        t = distance = Math.sqrt(closestDistSq);
                    else
                        t = distance = ( (rayProj - normProj*cylProj) - Math.sqrt(discriminant) ) / denom; 
                    
                    //find closest intersection point; make sure it's between atom's position and cylinder midpoint
                    
                    s = normProj*t - cylProj;
                    
                    //does not intersect cylinder between atom and midpoint,
                    // or intersects cylinder behind camera
                    if (s < 0 || s*s > cylinder.lengthSq() || t < 0)
                        continue;
                    
                    else
                        intersects.push({clickable : clickable,
                                         distance : distance});
                    
                }
                    
                
            }
            
        }
         
        //lines
        for (i = 0, il = intersectionShape.line.length; i < il; i += 2) {
            
            v1.copy(intersectionShape.line[i]);
            v1.applyMatrix4(group.matrixWorld);
            v2.copy(intersectionShape.line[i+1]);
            v2.applyMatrix4(group.matrixWorld);
            
            v3.subVectors(v2, v1);
            var bondLengthSq = v3.lengthSq();
            v3.normalize();
            
            w_0.subVectors(v1, raycaster.ray.origin);
            
            lineProj = w_0.dot(v3);
            rayProj = w_0.dot(raycaster.ray.direction);
            
            normProj = clamp(raycaster.ray.direction.dot(v3));
            
            denom = 1 - normProj*normProj;
            
            if (denom === 0.0)
                continue;
            
            s_c = (normProj*rayProj - lineProj) / denom;
            t_c = (rayProj - normProj*lineProj) / denom;
            
            v1.add(v3.multiplyScalar(s_c)); // Q_c
            v2.copy(raycaster.ray.direction).multiplyScalar(t_c).add(raycaster.ray.origin); // P_c
            
            closestDistSq = v3.subVectors(v2, v1).lengthSq();
            
            if (closestDistSq < precisionSq && s_c*s_c < bondLengthSq)
                intersects.push({clickable : clickable,
                                 distance : t_c
                                });
            
        }

        for (i = 0, il = intersectionShape.sphere.length; i < il; i++) {
            //sphere
            if (intersectionShape.sphere[i] instanceof WebMol.Sphere) {
                
                sphere.copy(intersectionShape.sphere[i]);
                sphere.applyMatrix4(group.matrixWorld);
                
                if (raycaster.ray.isIntersectionSphere(sphere)) {
                    
                    v1.subVectors(sphere.center, raycaster.ray.origin);
                    
                    //distance from ray origin to point on the ray normal to sphere's center
                    //must be less than sphere's radius (since ray intersects sphere)
                    var distanceToCenter = v1.dot(raycaster.ray.direction);
                    
                    discriminant = distanceToCenter*distanceToCenter - (v1.lengthSq() - sphere.radius*sphere.radius);
                    
                    //Don't select if sphere center behind camera
                    if (distanceToCenter < 0) 
                        return intersects;
                    
                    //ray tangent to sphere?
                    if (discriminant <= 0)
                        distance = distanceToCenter;
                    
                    //This is reversed if sphere is closer than ray origin.  Do we have 
                    //to worry about handling that case?
                    else 
                        distance = distanceToCenter - Math.sqrt(discriminant);
    
                    intersects.push({clickable : clickable, 
                                     distance : distance});
                    return intersects;
                }
            }        
       }
        
    };   
       
    Raycaster.prototype.precision = 0.0001;
    Raycaster.prototype.linePrecision = 0.2;
    
    Raycaster.prototype.set = function(origin, direction) {
        
        this.ray.set(origin, direction);
          
    };
    
    Raycaster.prototype.intersectObjects = function(group, objects) {
        
        var intersects = [];
        
        for (var i = 0, l = objects.length; i < l; i++)            
            intersectObject(group, objects[i], this, intersects);
            
        intersects.sort(descSort);
        
        return intersects;
        
    };
    
    return Raycaster;
    
})();


//WebMol Projecion 
//TODO: can probably strip this down a lot (only used for selection handling)
WebMol.Projector = function () {

    var _viewMatrix = new WebMol.Matrix4(),
    _viewProjectionMatrix = new WebMol.Matrix4();

    this.projectVector = function ( vector, camera ) {

        camera.matrixWorldInverse.getInverse( camera.matrixWorld );

        _viewProjectionMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

        return vector.applyProjection( _viewProjectionMatrix );

    };

    this.unprojectVector = function ( vector, camera ) {

        camera.projectionMatrixInverse.getInverse(camera.projectionMatrix);

        _viewProjectionMatrix.multiplyMatrices(camera.matrixWorld, camera.projectionMatrixInverse);

        return vector.applyProjection( _viewProjectionMatrix );

    };

};/*
 * Simplified Perspective Camera
 */


WebMol.Camera = function(fov, aspect, near, far) {
    
    WebMol.Object3D.call(this);
    
    this.fov = fov !== undefined ? fov : 50;
    this.aspect = aspect !== undefined ? aspect : 1;
    this.near = near !== undefined ? near : 0.1;
    this.far = far !== undefined ? far : 2000;

    this.projectionMatrix = new WebMol.Matrix4();
    this.projectionMatrixInverse = new WebMol.Matrix4();
    this.matrixWorldInverse = new WebMol.Matrix4();
    
    this.updateProjectionMatrix();
        
};

//Inherit Object3D's prototyped methods
WebMol.Camera.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Camera.prototype.lookAt = function(vector){
    
    //Why is the parameter order switched (compared to Object3D)?
    this.matrix.lookAt(this.position, vector, this.up);
    
    if (this.rotationAutoUpdate) {    
        
        if (this.useQuaternion === false) 
            this.rotation.setEulerFromRotationMatrix( this.matrix, this.eulerOrder );
        else
            this.quaternion.copy( this.matrix.decompose()[ 1 ] );    
            
    }
    
};

WebMol.Camera.prototype.updateProjectionMatrix = function () {

    this.projectionMatrix.makePerspective( this.fov, this.aspect, this.near, this.far );

};


//Render plugins go here

/**
 * Sprite render plugin
 */

WebMol.SpritePlugin = function () {

    var _gl, _renderer, _precision, _sprite = {};

    this.init = function ( renderer ) {

        _gl = renderer.context;
        _renderer = renderer;

        _precision = renderer.getPrecision();

        _sprite.vertices = new Float32Array( 8 + 8 );
        _sprite.faces    = new Uint16Array( 6 );

        var i = 0;

        _sprite.vertices[ i++ ] = -1; _sprite.vertices[ i++ ] = -1; // vertex 0
        _sprite.vertices[ i++ ] = 0;  _sprite.vertices[ i++ ] = 0;  // uv 0

        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = -1; // vertex 1
        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = 0;  // uv 1

        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = 1;  // vertex 2
        _sprite.vertices[ i++ ] = 1;  _sprite.vertices[ i++ ] = 1;  // uv 2

        _sprite.vertices[ i++ ] = -1; _sprite.vertices[ i++ ] = 1;  // vertex 3
        _sprite.vertices[ i++ ] = 0;  _sprite.vertices[ i++ ] = 1;  // uv 3

        i = 0;

        _sprite.faces[ i++ ] = 0; _sprite.faces[ i++ ] = 1; _sprite.faces[ i++ ] = 2;
        _sprite.faces[ i++ ] = 0; _sprite.faces[ i++ ] = 2; _sprite.faces[ i++ ] = 3;

        _sprite.vertexBuffer  = _gl.createBuffer();
        _sprite.elementBuffer = _gl.createBuffer();

        _gl.bindBuffer( _gl.ARRAY_BUFFER, _sprite.vertexBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, _sprite.vertices, _gl.STATIC_DRAW );

        _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, _sprite.elementBuffer );
        _gl.bufferData( _gl.ELEMENT_ARRAY_BUFFER, _sprite.faces, _gl.STATIC_DRAW );

        _sprite.program = createProgram( WebMol.ShaderLib.sprite, _precision );

        _sprite.attributes = {};
        _sprite.uniforms = {};

        _sprite.attributes.position           = _gl.getAttribLocation ( _sprite.program, "position" );
        _sprite.attributes.uv                 = _gl.getAttribLocation ( _sprite.program, "uv" );

        _sprite.uniforms.uvOffset             = _gl.getUniformLocation( _sprite.program, "uvOffset" );
        _sprite.uniforms.uvScale              = _gl.getUniformLocation( _sprite.program, "uvScale" );

        _sprite.uniforms.rotation             = _gl.getUniformLocation( _sprite.program, "rotation" );
        _sprite.uniforms.scale                = _gl.getUniformLocation( _sprite.program, "scale" );
        _sprite.uniforms.alignment            = _gl.getUniformLocation( _sprite.program, "alignment" );

        _sprite.uniforms.color                = _gl.getUniformLocation( _sprite.program, "color" );
        _sprite.uniforms.map                  = _gl.getUniformLocation( _sprite.program, "map" );
        _sprite.uniforms.opacity              = _gl.getUniformLocation( _sprite.program, "opacity" );

        _sprite.uniforms.useScreenCoordinates = _gl.getUniformLocation( _sprite.program, "useScreenCoordinates" );
        _sprite.uniforms.sizeAttenuation      = _gl.getUniformLocation( _sprite.program, "sizeAttenuation" );
        _sprite.uniforms.screenPosition       = _gl.getUniformLocation( _sprite.program, "screenPosition" );
        _sprite.uniforms.modelViewMatrix      = _gl.getUniformLocation( _sprite.program, "modelViewMatrix" );
        _sprite.uniforms.projectionMatrix     = _gl.getUniformLocation( _sprite.program, "projectionMatrix" );

        _sprite.uniforms.fogType              = _gl.getUniformLocation( _sprite.program, "fogType" );
        _sprite.uniforms.fogDensity           = _gl.getUniformLocation( _sprite.program, "fogDensity" );
        _sprite.uniforms.fogNear              = _gl.getUniformLocation( _sprite.program, "fogNear" );
        _sprite.uniforms.fogFar               = _gl.getUniformLocation( _sprite.program, "fogFar" );
        _sprite.uniforms.fogColor             = _gl.getUniformLocation( _sprite.program, "fogColor" );

        _sprite.uniforms.alphaTest            = _gl.getUniformLocation( _sprite.program, "alphaTest" );

    };

    this.render = function ( scene, camera, viewportWidth, viewportHeight ) {

        var sprites = scene.__webglSprites,
            nSprites = sprites.length;

        if ( ! nSprites ) return;

        var attributes = _sprite.attributes,
            uniforms = _sprite.uniforms;

        var invAspect = viewportHeight / viewportWidth;

        var halfViewportWidth = viewportWidth * 0.5,
            halfViewportHeight = viewportHeight * 0.5;

        // setup gl

        _gl.useProgram( _sprite.program );

        _gl.enableVertexAttribArray( attributes.position );
        _gl.enableVertexAttribArray( attributes.uv );

        _gl.disable( _gl.CULL_FACE );
        _gl.enable( _gl.BLEND );

        _gl.bindBuffer( _gl.ARRAY_BUFFER, _sprite.vertexBuffer );
        _gl.vertexAttribPointer( attributes.position, 2, _gl.FLOAT, false, 2 * 8, 0 );
        _gl.vertexAttribPointer( attributes.uv, 2, _gl.FLOAT, false, 2 * 8, 8 );

        _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, _sprite.elementBuffer );

        _gl.uniformMatrix4fv( uniforms.projectionMatrix, false, camera.projectionMatrix.elements );

        _gl.activeTexture( _gl.TEXTURE0 );
        _gl.uniform1i( uniforms.map, 0 );

        var oldFogType = 0;
        var sceneFogType = 0;
        var fog = scene.fog;

        if ( fog ) {

            _gl.uniform3f( uniforms.fogColor, fog.color.r, fog.color.g, fog.color.b );

            _gl.uniform1f( uniforms.fogNear, fog.near );
            _gl.uniform1f( uniforms.fogFar, fog.far );

            _gl.uniform1i( uniforms.fogType, 1 );
            oldFogType = 1;
            sceneFogType = 1;


        } 
        
        else {

            _gl.uniform1i( uniforms.fogType, 0 );
            oldFogType = 0;
            sceneFogType = 0;

        }


        // update positions and sort

        var i, sprite, material, screenPosition, size, fogType, scale = [];

        for( i = 0; i < nSprites; i ++ ) {

            sprite = sprites[ i ];
            material = sprite.material;

            if ( ! sprite.visible || material.opacity === 0 ) continue;

            if ( ! material.useScreenCoordinates ) {

                sprite._modelViewMatrix.multiplyMatrices( camera.matrixWorldInverse, sprite.matrixWorld );
                sprite.z = - sprite._modelViewMatrix.elements[ 14 ];

            } else {

                sprite.z = - sprite.position.z;

            }

        }

        sprites.sort( painterSortStable );

        // render all sprites

        for( i = 0; i < nSprites; i ++ ) {

            sprite = sprites[ i ];
            material = sprite.material;

            if ( ! sprite.visible || material.opacity === 0 ) continue;

            if ( material.map && material.map.image && material.map.image.width ) {

                _gl.uniform1f( uniforms.alphaTest, material.alphaTest );

                if ( material.useScreenCoordinates === true ) {

                    _gl.uniform1i( uniforms.useScreenCoordinates, 1 );
                    _gl.uniform3f(
                        uniforms.screenPosition,
                        ( ( sprite.position.x * _renderer.devicePixelRatio ) - halfViewportWidth  ) / halfViewportWidth,
                        ( halfViewportHeight - ( sprite.position.y * _renderer.devicePixelRatio ) ) / halfViewportHeight,
                        Math.max( 0, Math.min( 1, sprite.position.z ) )
                    );

                    scale[ 0 ] = _renderer.devicePixelRatio;
                    scale[ 1 ] = _renderer.devicePixelRatio;

                } else {

                    _gl.uniform1i( uniforms.useScreenCoordinates, 0 );
                    _gl.uniform1i( uniforms.sizeAttenuation, material.sizeAttenuation ? 1 : 0 );
                    _gl.uniformMatrix4fv( uniforms.modelViewMatrix, false, sprite._modelViewMatrix.elements );

                    scale[ 0 ] = 1;
                    scale[ 1 ] = 1;

                }

                if ( scene.fog && material.fog ) {

                    fogType = sceneFogType;

                } else {

                    fogType = 0;

                }

                if ( oldFogType !== fogType ) {

                    _gl.uniform1i( uniforms.fogType, fogType );
                    oldFogType = fogType;

                }

                size = 1 / ( material.scaleByViewport ? viewportHeight : 1 );

                scale[ 0 ] *= size * invAspect * sprite.scale.x;
                scale[ 1 ] *= size * sprite.scale.y;

                _gl.uniform2f( uniforms.uvScale, material.uvScale.x, material.uvScale.y );
                _gl.uniform2f( uniforms.uvOffset, material.uvOffset.x, material.uvOffset.y );
                _gl.uniform2f( uniforms.alignment, material.alignment.x, material.alignment.y );

                _gl.uniform1f( uniforms.opacity, material.opacity );
                _gl.uniform3f( uniforms.color, material.color.r, material.color.g, material.color.b );

                _gl.uniform1f( uniforms.rotation, sprite.rotation );
                _gl.uniform2fv( uniforms.scale, scale );

                //_renderer.setBlending( material.blending, material.blendEquation, material.blendSrc, material.blendDst );
                _renderer.setDepthTest( material.depthTest );
                _renderer.setDepthWrite( material.depthWrite );
                _renderer.setTexture( material.map, 0 );

                _gl.drawElements( _gl.TRIANGLES, 6, _gl.UNSIGNED_SHORT, 0 );

            }

        }

        // restore gl

        _gl.enable( _gl.CULL_FACE );

    };

    function createProgram ( shader, precision ) {

        var program = _gl.createProgram();

        var fragmentShader = _gl.createShader( _gl.FRAGMENT_SHADER );
        var vertexShader = _gl.createShader( _gl.VERTEX_SHADER );

        var prefix = "precision " + precision + " float;\n";

        _gl.shaderSource( fragmentShader, prefix + shader.fragmentShader );
        _gl.shaderSource( vertexShader, prefix + shader.vertexShader );

        _gl.compileShader( fragmentShader );
        _gl.compileShader( vertexShader );
        
        if ( ! _gl.getShaderParameter(fragmentShader, _gl.COMPILE_STATUS) || ! _gl.getShaderParameter(vertexShader,_gl.COMPILE_STATUS) ) {

                console.error(_gl.getShaderInfoLog(fragmentShader));
                console.error("could not initialize shader");
                return null;
        }

        _gl.attachShader( program, fragmentShader );
        _gl.attachShader( program, vertexShader );

        _gl.linkProgram( program );

        if (! _gl.getProgramParameter(program, _gl.LINK_STATUS) )
                console.error("Could not initialize shader");

        return program;

    }

    function painterSortStable ( a, b ) {

        if ( a.z !== b.z ) {

            return b.z - a.z;

        } else {

            return b.id - a.id;

        }

    }

};
/* 
 * WebMol Lighting
 */

//TODO: Strip down this class - do I really use all of these instance variables?
WebMol.Light = function(hex, intensity) {
    
    WebMol.Object3D.call(this);
    
    this.color = new WebMol.Color(hex);
    this.position = new WebMol.Vector3( 0, 1, 0 );
    this.target = new WebMol.Object3D();

    this.intensity = ( intensity !== undefined ) ? intensity : 1;

    this.castShadow = false;
    this.onlyShadow = false;
    
};

WebMol.Light.prototype = Object.create(WebMol.Object3D.prototype);
/* 
 * Line and Mesh material types
 */

//Material base class

WebMol.Material = function () {

    WebMol.EventDispatcher.call( this );

    this.id = WebMol.MaterialIdCount ++;

    this.name = '';
    
    //TODO: Which of these instance variables can I remove??
    this.side = WebMol.FrontSide;

    this.opacity = 1;
    this.transparent = false;

    this.blending = WebMol.NormalBlending;

    this.depthTest = true;
    this.depthWrite = true;

    this.polygonOffset = false;
    this.polygonOffsetFactor = 0;
    this.polygonOffsetUnits = 0;

    this.alphaTest = 0;

    this.visible = true;

    this.needsUpdate = true;

};


WebMol.Material.prototype.setValues = function ( values ) {

    if ( values === undefined ) return;

    for ( var key in values ) {

        var newValue = values[ key ];

        if ( newValue === undefined ) {

            console.warn( 'WebMol.Material: \'' + key + '\' parameter is undefined.' );
            continue;

        }

        if ( key in this ) {

            var currentValue = this[ key ];

            if ( currentValue instanceof WebMol.Color && newValue instanceof WebMol.Color ) {

                currentValue.copy( newValue );

            } else if ( currentValue instanceof WebMol.Color ) {

                currentValue.set( newValue );

            } else if ( currentValue instanceof WebMol.Vector3 && newValue instanceof WebMol.Vector3 ) {

                currentValue.copy( newValue );

            } else {

                this[ key ] = newValue;

            }

        }

    }

};
//TODO: might want to look into blending equations
WebMol.Material.prototype.clone = function ( material ) {

    if ( material === undefined ) material = new WebMol.Material();

    material.name = this.name;

    material.side = this.side;

    material.opacity = this.opacity;
    material.transparent = this.transparent;

    material.blending = this.blending;

    material.depthTest = this.depthTest;
    material.depthWrite = this.depthWrite;

    material.polygonOffset = this.polygonOffset;
    material.polygonOffsetFactor = this.polygonOffsetFactor;
    material.polygonOffsetUnits = this.polygonOffsetUnits;

    material.alphaTest = this.alphaTest;

    material.overdraw = this.overdraw;

    material.visible = this.visible;

    return material;

};

WebMol.Material.prototype.dispose = function () {

    this.dispatchEvent( { type: 'dispose' } );

};

WebMol.MaterialIdCount = 0;

//Line basic material

WebMol.LineBasicMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    
    this.linewidth = 1;
    this.linecap = 'round';
    this.linejoin = 'round';
    
    this.vertexColors = false;
    
    this.fog = true;
    
    this.setValues(parameters);
    
};

WebMol.LineBasicMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.LineBasicMaterial.prototype.clone = function() {
  
    var material = new WebMol.LineBasicMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy();
    
};

//Mesh Lambert material

WebMol.MeshLambertMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    this.ambient = new WebMol.Color(0xfffff);
    this.emissive = new WebMol.Color(0x000000);
    
    //TODO: Which of these instance variables do I really need?
    this.wrapAround = false;
    this.wrapRGB = new WebMol.Vector3(1,1,1);
    
    this.map = null;
    
    this.lightMap = null;
    
    this.specularMap = null;
    
    this.envMap = null;
    this.reflectivity = 1;
    this.refractionRatio = 0.98;
    
    this.fog = true;
    
    this.wireframe = false;
    this.wireframeLinewidth = 1;
    this.wireframeLinecap = 'round';
    this.wireframeLinejoin = 'round';
    
    this.shading = WebMol.SmoothShading;
    
    this.vertexColors = WebMol.NoColors;
    
    this.skinning = false;
    
    this.setValues(parameters);
    
};

WebMol.MeshLambertMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.MeshLambertMaterial.prototype.clone = function() {
  
    var material = new WebMol.MeshLambertMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy(this.color);
    material.ambient.copy(this.ambient);
    material.emissive.copy(this.emissive);
    
    material.wrapAround = this.wrapAround;
    material.wrapRGB.copy(this.wrapRGB);
    
    material.map = this.map;
    
    material.lightMap = this.lightMap;
    
    material.specularMap = this.specularMap;
    
    material.envMap = this.envMap;
    material.combine = this.combine;
    material.reflectivity = this.reflectivity;
    material.refractionRatio = this.refractionRatio;
    
    material.fog = this.fog;
    
    material.shading = this.shading;
    
    material.vertexColors = this.vertexColors;
    
    material.skinning = this.skinning;
    material.morphTargets = this.morphTargets;
    material.morphNormals = this.morphNormals;
    
    return material;
    
};


//Sprite material
WebMol.SpriteMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    this.map = new WebMol.Texture();
    
    this.useScreenCoordinates = true;
    this.depthTest = !this.useScreenCoordinates;
    this.sizeAttenuation = !this.useScreenCoordinates;
    this.scaleByViewPort = !this.sizeAttenuation;
    this.alignment = WebMol.SpriteAlignment.center.clone();
    
    this.fog = false; // use scene fog
    
    this.uvOffset = new WebMol.Vector2(0, 0);
    this.uvScale = new WebMol.Vector2(1, 1);
    
    this.setValues(parameters);
    
    parameters = parameters || {};
    
    if (parameters.depthTest === undefined)
        this.depthTest = !this.useScreenCoordinates;
    if (parameters.sizeAttenuation === undefined)
        this.sizeAttenuation = !this.useScreenCoordinates;
    if (parameters.scaleByViewPort === undefined)
        this.scaleByViewPort = !this.sizeAttenuation;
    
};

WebMol.SpriteMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.SpriteMaterial.prototype.clone = function() {
    
    var material = new WebMol.SpriteMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy(this.color);
    material.map = this.map;
    
    material.useScreenCoordinates = useScreenCoordinates;
    material.sizeAttenuation = this.sizeAttenuation;
    material.scaleByViewport = this.scaleByViewPort;
    material.alignment.copy(this.alignment);
    
    material.uvOffset.copy(this.uvOffset);
    
    return material;
    
};

//Alignment for Sprites

WebMol.SpriteAlignment = {};
WebMol.SpriteAlignment.topLeft = new WebMol.Vector2(1, -1);
WebMol.SpriteAlignment.topCenter = new WebMol.Vector2(0, -1);
WebMol.SpriteAlignment.topRight = new WebMol.Vector2(-1, -1);
WebMol.SpriteAlignment.centerLeft = new WebMol.Vector2(1, 0);
WebMol.SpriteAlignment.center = new WebMol.Vector2(0, 0);
WebMol.SpriteAlignment.centerRight = new WebMol.Vector2(-1, 0);
WebMol.SpriteAlignment.bottomLeft = new WebMol.Vector2(1, 1);
WebMol.SpriteAlignment.bottomCenter = new WebMol.Vector2(0, 1);
WebMol.SpriteAlignment.bottomRight = new WebMol.Vector2(-1, 1);


//Texture
//We really only create textures from 2d rendering contexts (to display text labels)

WebMol.Texture = function(image) {

    WebMol.EventDispatcher.call(this);
    
    this.id = WebMol.TextureIdCount++;
    
    this.name = "";
    
    this.image = image;
    this.mipmaps = [];
    
    this.mapping = new WebMol.UVMapping();
    
    this.wrapS = WebMol.ClampToEdgeWrapping;
    this.wrapT = WebMol.ClampToEdgeWrapping;
    
    this.magFilter = WebMol.LinearFilter;
    this.minFilter = WebMol.LinearMipMapLinearFilter;
    
    this.anisotropy = 1;
    
    this.format = WebMol.RGBAFormat;
    this.type = WebMol.UnsignedByteType;
    
    this.offset = new WebMol.Vector2(0, 0);
    this.repeat = new WebMol.Vector2(1, 1);
    
    this.generateMipmaps = true;
    this.premultiplyAlpha = false;
    this.flipY = true;
    this.unpackAlignment = 4;
    
    this.needsUpdate = false;
    this.onUpdate = null;
    
};

WebMol.Texture.prototype = {

    constructor : WebMol.Texture,
    
    clone : function(texture) {
        
        if (texture === undefined)
            texture = new WebMol.Texture();
        
        texture.image = this.image;
        texture.mipmaps = this.mipmaps.slice(0);
        
        texture.mapping = this.mapping;
        
        texture.wrapS = this.wrapS;
        texture.wrapT = this.wrapT;
        
        texture.magFilter = this.magFilter;
        texture.minFilter = this.minFilter;
        
        texture.anisotropy = this.anisotropy;
        
        texture.format = this.format;
        texture.type = this.type;
        
        texture.offset.copy(this.offset);
        texture.repeat.copy(this.repeat);
        
        texture.generateMipmaps = this.generateMipmaps;
        texture.premultiplyAlpha = this.premultiplyAlpha;
        texture.flipY = this.flipY;
        texture.unpackAlignment = this.unpackAlignment;
        
        return texture;
        
    },
    
    dispose : function() {
        
        this.dispatchEvent( {type: 'dispose'});
        
    }    
    
};

WebMol.TextureIdCount = 0;


/* 
 * WebMol Mesh and Line objects
 */


//Line Object

WebMol.Line = function (geometry, material, type) {

    WebMol.Object3D.call(this);

    this.geometry = geometry;
        //TODO: update material and type to webgl
    this.material = (material !== undefined) ? material : new WebMol.LineBasicMaterial( { color: Math.random() * 0xffffff } );
    this.type = (type !== undefined) ? type : WebMol.LineStrip;

};

WebMol.LineStrip = 0;
WebMol.LinePieces = 1;

WebMol.Line.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Line.prototype.clone = function (object) {

    if (object === undefined) object = new WebMol.Line(this.geometry, this.material, this.type);

    WebMol.Object3D.prototype.clone.call(this, object);

    return object;

};


//Mesh Object

WebMol.Mesh = function(geometry, material) {

    WebMol.Object3D.call(this);

    this.geometry = geometry;
    this.material = (material !== undefined) ? material : new WebMol.MeshBasicMaterial( { color: Math.random() * 0xffffff, wireframe: true } );

};

WebMol.Mesh.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Mesh.prototype.clone = function (object) {

    if (object === undefined) object = new WebMol.Mesh(this.geometry, this.material);

    WebMol.Object3D.prototype.clone.call(this, object);

    return object;

};


//Sprite object

WebMol.Sprite = function(material) {
    
    WebMol.Object3D.call(this);
    
    this.material = (material !== undefined) ? material : new WebMol.SpriteMaterial();

    this.rotation3d = this.rotation;
    this.rotation = 0;
    
};

WebMol.Sprite.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Sprite.prototype.updateMatrix = function() {
    
    this.matrix.setPosition(this.position);
    
    this.rotation3d.set(0, 0, this.rotation);
    this.matrix.setRotationFromEuler(this.rotation3d);
    
    if (this.scale.x !== 1 || this.scale.y !== 1)
        this.matrix.scale(this.scale);
    
    this.matrixWorldNeedsUpdate = true;
    
};

WebMol.Sprite.prototype.clone = function(object) {
    
    if (object === undefined)
        object = new WebMol.Sprite(this.material);
    
    WebMol.Object3D.prototype.clone.call(this, object);
    
    return object;
    
};
/**
Simplified webGL renderer 
 */

WebMol.Renderer = function ( parameters ) {
    
    parameters = parameters || {};
    
    var _canvas = parameters.canvas !== undefined ? parameters.canvas : document.createElement( 'canvas' ),

    _precision = parameters.precision !== undefined ? parameters.precision : 'highp',

    _alpha = parameters.alpha !== undefined ? parameters.alpha : true,
    _premultipliedAlpha = parameters.premultipliedAlpha !== undefined ? parameters.premultipliedAlpha : true,
    _antialias = parameters.antialias !== undefined ? parameters.antialias : false,
    _stencil = parameters.stencil !== undefined ? parameters.stencil : true,
    _preserveDrawingBuffer = parameters.preserveDrawingBuffer !== undefined ? parameters.preserveDrawingBuffer : false,

    _clearColor = parameters.clearColor !== undefined ? new WebMol.Color( parameters.clearColor ) : new WebMol.Color( 0x000000 ),
    _clearAlpha = parameters.clearAlpha !== undefined ? parameters.clearAlpha : 0;
    
    this.domElement = _canvas;
    this.context = null;
    this.devicePixelRatio = parameters.devicePixelRatio !== undefined ? 
    parameters.devicePixelRatio : (self.devicePixelRatio !== undefined) ? 
                                   self.devicePixelRatio : 1;

    // clearing

    this.autoClear = true;
    this.autoClearColor = true;
    this.autoClearDepth = true;
    this.autoClearStencil = true;

    // scene graph

    this.sortObjects = true;

    this.autoUpdateObjects = true;
    this.autoUpdateScene = true;
    
    this.renderPluginsPost = [];
    
    // info

    this.info = {

        memory: {
    
        programs: 0,
        geometries: 0,
        textures: 0
    
        },
    
        render: {
    
        calls: 0,
        vertices: 0,
        faces: 0,
        points: 0
    
        }

    };

    // internal properties

    var _this = this,

    _programs = [],
    _programs_counter = 0,

    // internal state cache

    _currentProgram = null,
    _currentFramebuffer = null,
    _currentMaterialId = -1,
    _currentGeometryGroupHash = null,
    _currentCamera = null,
    _geometryGroupCounter = 0,

    _usedTextureUnits = 0,

    // GL state cache

    _oldDoubleSided = -1,
    _oldFlipSided = -1,

    _oldBlending = -1,

    _oldBlendEquation = -1,
    _oldBlendSrc = -1,
    _oldBlendDst = -1,

    _oldDepthTest = -1,
    _oldDepthWrite = -1,

    _oldPolygonOffset = null,
    _oldPolygonOffsetFactor = null,
    _oldPolygonOffsetUnits = null,

    _oldLineWidth = null,

    _viewportX = 0,
    _viewportY = 0,
    _viewportWidth = 0,
    _viewportHeight = 0,
    _currentWidth = 0,
    _currentHeight = 0,

    _enabledAttributes = {},

     // camera matrices cache

    _projScreenMatrix = new WebMol.Matrix4(),

    _vector3 = new WebMol.Vector3(),

    // light arrays cache

    _direction = new WebMol.Vector3(),

    _lightsNeedUpdate = true,

    _lights = {

            ambient: [0,0,0],
            directional: { length: 0, colors: [], positions: [] },
            point: { length: 0, colors: [], positions: [], distances: [] },
            spot: { length: 0, colors: [], positions: [], distances: [], directions: [], anglesCos: [], exponents: [] },
            hemi: { length: 0, skyColors: [], groundColors: [], positions: [] }

    };

    // initialize

    var _gl;

    initGL();

    setDefaultGLState();

    this.context = _gl;    

    // API

    this.getContext = function () {

            return _gl;

    };

    this.getPrecision = function () {

            return _precision;

    };
    
    this.setClearColorHex = function ( hex, alpha ) {

            _clearColor.setHex( hex );
            _clearAlpha = alpha;

            _gl.clearColor( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha );

    };

    this.setSize = function ( width, height ) {

            _canvas.width = width * this.devicePixelRatio;
            _canvas.height = height * this.devicePixelRatio;

            _canvas.style.width = width + 'px';
            _canvas.style.height = height + 'px';

            this.setViewport( 0, 0, _canvas.width, _canvas.height );

    };

    this.setViewport = function ( x, y, width, height ) {

            _viewportX = x !== undefined ? x : 0;
            _viewportY = y !== undefined ? y : 0;

            _viewportWidth = width !== undefined ? width : _canvas.width;
            _viewportHeight = height !== undefined ? height : _canvas.height;

            _gl.viewport( _viewportX, _viewportY, _viewportWidth, _viewportHeight );

    };

    this.clear = function ( color, depth, stencil ) {

            var bits = 0;

            if ( color === undefined || color ) bits |= _gl.COLOR_BUFFER_BIT;
            if ( depth === undefined || depth ) bits |= _gl.DEPTH_BUFFER_BIT;
            if ( stencil === undefined || stencil ) bits |= _gl.STENCIL_BUFFER_BIT;

            _gl.clear( bits );

    };

    this.clearTarget = function ( renderTarget, color, depth, stencil ) {

            this.setRenderTarget( renderTarget );
            this.clear( color, depth, stencil );

    };

    this.setMaterialFaces = function ( material ) {

            var doubleSided = material.side === WebMol.DoubleSide;
            var flipSided = material.side === WebMol.BackSide;

            if ( _oldDoubleSided !== doubleSided ) {

                if ( doubleSided ) {

                    _gl.disable( _gl.CULL_FACE );

                } else {

                    _gl.enable( _gl.CULL_FACE );

                }

                _oldDoubleSided = doubleSided;

            }

            if ( _oldFlipSided !== flipSided ) {

                if ( flipSided ) {

                    _gl.frontFace( _gl.CW );

                } else {

                    _gl.frontFace( _gl.CCW );

                }

                _oldFlipSided = flipSided;

            }    

    };
    
    this.setDepthTest = function ( depthTest ) {

            if ( _oldDepthTest !== depthTest ) {

                if ( depthTest ) {

                    _gl.enable( _gl.DEPTH_TEST );

                } else {

                    _gl.disable( _gl.DEPTH_TEST );

                }

                _oldDepthTest = depthTest;

            }

    };

    this.setDepthWrite = function ( depthWrite ) {

            if ( _oldDepthWrite !== depthWrite ) {

                    _gl.depthMask( depthWrite );
                    _oldDepthWrite = depthWrite;

            }

    };

    this.setBlending = function( blending ) {

            if (blending === WebMol.NoBlending) 
                    _gl.disable( _gl.BLEND );

            else {
                    _gl.enable( _gl.BLEND );
                    _gl.blendEquationSeparate( _gl.FUNC_ADD, _gl.FUNC_ADD );
                    _gl.blendFuncSeparate( _gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA, _gl.ONE, _gl.ONE_MINUS_SRC_ALPHA );
            }

            _oldBlending = blending;
    };
    
    // Plugins
    
    this.addPostPlugin = function(plugin) {

        plugin.init(this);
        this.renderPluginsPost.push(plugin);

    };

    // Sorting

    function numericalSort ( a, b ) {

            return b[ 0 ] - a[ 0 ];

    }

    function enableAttribute( attribute ) {

        if ( ! _enabledAttributes[ attribute ] ) {

            _gl.enableVertexAttribArray( attribute );
            _enabledAttributes[ attribute ] = true;

        }

    }

    function disableAttributes() {

        for ( var attribute in _enabledAttributes ) {

            if ( _enabledAttributes[ attribute ] ) {

                _gl.disableVertexAttribArray( attribute );
                _enabledAttributes[ attribute ] = false;

            }

        }

    } 

    function setPolygonOffset ( polygonOffset, factor, units) {

        if ( _oldPolygonOffset !== polygonOffset ) {

            if (polygonOffset)
                _gl.enable( _gl.POLYGON_OFFSET_FILL );
            else
                _gl.disable( _gl.POLYGON_OFFSET_FILL );
        }
    }

    function setLineWidth ( width ) {

        if ( width !== _oldLineWidth ) {
            _gl.lineWidth(width);
            _oldLineWidth = width;
        }

    }
    
    var onGeometryDispose = function(event) {
        
        var geometry = event.target;
        geometry.removeEventListener('dispose', onGeometryDispose);
        
        deallocateGeometry(geometry);
        
        _this.info.memory.geometries--;
        
    };
    
    var onTextureDispose = function(event) {

        var texture = event.target;

        texture.removeEventListener('dispose', onTextureDispose);

        deallocateTexture(texture);

        _this.info.memory.textures--;


    };
    
    var onMaterialDispose = function(event) {
        
        var material = event.target;
        material.removeEventListener('dispose', onMaterialDispose);
        
        deallocateMaterial(material);
        
    };
    
    var deallocateGeometry = function(geometry) {
        
        geometry.__webglInit = undefined;
        
        if (geometry.__webglVertexBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglVertexBuffer);
        
        if (geometry.__webglColorBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglColorBuffer);
        
        if (geometry.geometryGroups !== undefined) {
            
            for (var g = 0, gl = geometry.groups; g < gl; g++) {  
                
                var geometryGroup = geometry.geometryGroups[g];

                if (geometryGroup.__webglVertexBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglVertexBuffer);

                if (geometryGroup.__webglColorBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglColorBuffer);
                
                if (geometryGroup.__webglNormalBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglNormalBuffer);  
                
                if (geometryGroup.__webglFaceBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglFaceBuffer);
                    
                if (geometryGroup.__webglLineBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglLineBuffer);
                    
            }
        }
    };
    
    var deallocateMaterial = function (material) {

        var program = material.program;

        if ( program === undefined ) return;

        material.program = undefined;

        // only deallocate GL program if this was the last use of shared program
        // assumed there is only single copy of any program in the _programs list
        // (that's how it's constructed)

        var i, il, programInfo;
        var deleteProgram = false;

        for ( i = 0, il = _programs.length; i < il; i ++ ) {

            programInfo = _programs[ i ];

            if ( programInfo.program === program ) {

                programInfo.usedTimes --;

                if ( programInfo.usedTimes === 0 ) {

                    deleteProgram = true;

                }

                break;

            }

        }

        if ( deleteProgram === true ) {

            // avoid using array.splice, this is costlier than creating new array from scratch

            var newPrograms = [];

            for ( i = 0, il = _programs.length; i < il; i ++ ) {

                programInfo = _programs[ i ];

                if ( programInfo.program !== program ) {

                    newPrograms.push( programInfo );

                }

            }

            _programs = newPrograms;

            _gl.deleteProgram( program );

            _this.info.memory.programs --;

        }

    };
    
    var deallocateTexture = function(texture) {

        if (texture.image && texture.image.__webglTextureCube) {

            // cube texture

            _gl.deleteTexture(texture.image.__webglTextureCube);

        } 
        
        else {

            // 2D texture

            if ( ! texture.__webglInit ) return;

            texture.__webglInit = false;
            _gl.deleteTexture( texture.__webglTexture );

        }

    };

    //Compile and return shader
    function getShader (type, str) {

        var shader;

        if (type === "fragment")
            shader = _gl.createShader( _gl.FRAGMENT_SHADER );
        else if (type === "vertex")
            shader = _gl.createShader( _gl.VERTEX_SHADER );

        _gl.shaderSource(shader, str);
        _gl.compileShader(shader);

        if ( ! _gl.getShaderParameter(shader, _gl.COMPILE_STATUS) ) {

            console.error(_gl.getShaderInfoLog(shader));
            console.error("could not initialize shader");
            return null;

        }

        return shader;

    } 


    //Compile appropriate shaders (if necessary) from source code and attach to gl program.
    function buildProgram(fragmentShader, vertexShader, uniforms, parameters) {

        var p, pl, d, program, code;
        var chunks = [];

        chunks.push(fragmentShader);
        chunks.push(vertexShader);
        
        for (p in parameters) {
            chunks.push(p);
            chunks.push(parameters[p]);
        }
        
        code = chunks.join();

        //check if program has already been compiled

        for (p = 0, pl = _programs.length; p < pl; p++) {

            var programInfo = _programs[p];

            if (programInfo.code === code) {

                programInfo.usedTimes++;

                return programInfo.program;
            }
        }

        //Set up new program and compile shaders

        program = _gl.createProgram();
        
        //set up precision
        var precision = _precision;
        var prefix = "precision " + precision + " float;";
        
        var prefix_vertex = [
                             prefix
                            ].join("\n");
                            
        var prefix_fragment = [                               
                               parameters.wireframe ? "#define WIREFRAME 1" : "",
                               prefix
                              ].join("\n");
        
        var glFragmentShader = getShader("fragment", prefix_fragment + fragmentShader);
        var glVertexShader = getShader("vertex", prefix_vertex + vertexShader);

        _gl.attachShader(program, glVertexShader);
        _gl.attachShader(program, glFragmentShader);

        _gl.linkProgram(program);

        if (! _gl.getProgramParameter(program, _gl.LINK_STATUS) )
            console.error("Could not initialize shader");

        //gather and cache uniform variables and attributes

        program.uniforms = {};
        program.attributes = {};

        var identifiers, u, a, i;

        //uniform vars
        identifiers = 
            [ 'viewMatrix', 'modelViewMatrix', 'projectionMatrix', 'normalMatrix', 'modelMatrix', 'cameraPosition' ];

        //custom uniform vars
        for (u in uniforms) 
            identifiers.push(u);

        for (i = 0; i < identifiers.length; i++) {

            var uniformVar = identifiers[i];
            program.uniforms[uniformVar] = _gl.getUniformLocation(program, uniformVar);

        }

        //attributes
        identifiers = 
            [ 'position', 'normal', 'color', 'lineDistance' ];

        /*
        for (a in attributes)
                identifiers.push(a);
        */

        for (i = 0; i < identifiers.length; i++) {

            var attributeVar = identifiers[i];
            program.attributes[attributeVar] = _gl.getAttribLocation(program, attributeVar);
        }

        program.id = _programs_counter++;
        _programs.push( {program: program, code: code, usedTimes: 1} );
        _this.info.memory.programs = _programs.length;

        return program;
    }

    //TODO: need to set up shader attributes and uniforms as attributes on material object after attaching prgm
    //We need to attach appropriate uniform variables to material after shaders have been chosen
    this.initMaterial = function ( material, lights, fog, object ) {

        material.addEventListener('dispose', onMaterialDispose);

        var u, a, identifiers, i, parameters, maxLightCount, maxBones, maxShadows, shaderID;

        if (material instanceof WebMol.LineBasicMaterial)
            shaderID = "basic";
        else if (material instanceof WebMol.MeshLambertMaterial)
            shaderID = "lambert";

        if (shaderID) {

            var shader = WebMol.ShaderLib[shaderID];
            material.shaderType = shaderID;
            material.vertexShader = shader.vertexShader;
            material.fragmentShader = shader.fragmentShader;
            material.uniforms = WebMol.ShaderUtils.clone(shader.uniforms);
            //TODO: set material uniforms to shader uniform variables

        }
        
        parameters = {
            wireframe: material.wireframe
        };

        material.program = buildProgram(material.fragmentShader, material.vertexShader, material.uniforms, parameters);

    };

    function setProgram( camera, lights, fog, material, object ) {

        if ( material.needsUpdate ) {

            if (material.program)
                deallocateMaterial(material);

                _this.initMaterial( material, lights, fog, object );
                material.needsUpdate = false;
        }

        var refreshMaterial = false;

        //p_uniforms: uniformVarName => uniformLocation
        //m_uniforms: uniformVarName => uniformJsVal
        var program = material.program,
            p_uniforms = program.uniforms,
            m_uniforms = material.uniforms;

        if (program != _currentProgram) {        
            _gl.useProgram(program);
            _currentProgram = program;

            refreshMaterial = true;
        }

        if (material.id != _currentMaterialId) {
            _currentMaterialId = material.id;
            refreshMaterial = true;
        }

        if (camera != _currentCamera) {    
            _currentCamera = camera;
            refreshMaterial = true;
        }

        //Send projection matrix to uniform variable in shader
        if (refreshMaterial) {

            //Load projection, model-view matrices for perspective
            _gl.uniformMatrix4fv(p_uniforms.projectionMatrix, false, camera.projectionMatrix.elements);
            _gl.uniformMatrix4fv(p_uniforms.modelViewMatrix, false, object._modelViewMatrix.elements);

            //Set up correct fog uniform vals
            m_uniforms.fogColor.value = fog.color;
            m_uniforms.fogNear.value = fog.near;
            m_uniforms.fogFar.value = fog.far;

            //Set up lights for lambert shader
            if (material.shaderType === "lambert") {

                //load view and normal matrices for directional and object lighting
                _gl.uniformMatrix4fv(p_uniforms.viewMatrix, false, camera.matrixWorldInverse.elements);
                _gl.uniformMatrix3fv(p_uniforms.normalMatrix, false, object._normalMatrix.elements);
                //_gl.uniformMatrix4fv(p_uniforms.modelMatrix, false, object.matrixWorld.elements);

                if (_lightsNeedUpdate) {
                    setupLights(program, lights);
                    _lightsNeedUpdate = false;
                }

                //Set up correct light uniform var vals
                m_uniforms.ambientLightColor.value = _lights.ambient;
                m_uniforms.directionalLightColor.value = _lights.directional.colors;
                m_uniforms.directionalLightDirection.value = _lights.directional.positions;
                m_uniforms.ambient.value = material.ambient;
                m_uniforms.emissive.value = material.emissive;

            }

            //opacity, diffuse, emissive, etc
            m_uniforms.opacity.value = material.opacity;
            m_uniforms.diffuse.value = material.color;

            //Load any other material specific uniform variables to gl shaders
            loadMaterialUniforms(p_uniforms, m_uniforms);

        }

        return program;

    }

    function loadMaterialUniforms(p_uniforms, m_uniforms) {
        var uniformVar, type, uniformVal, uniformLoc;

        for (uniformVar in m_uniforms) {
            if (! p_uniforms[uniformVar])
                continue;

            type = m_uniforms[uniformVar].type;
            uniformVal = m_uniforms[uniformVar].value;
            uniformLoc = p_uniforms[uniformVar];

            //single float
            if (type === 'f')
                _gl.uniform1f(uniformLoc, uniformVal);
            //array of floats
            else if (type === 'fv')
                _gl.uniform3fv(uniformLoc, uniformVal);
            //color - r,g,b floats
            else if (type === 'c')
                _gl.uniform3f(uniformLoc, uniformVal.r, uniformVal.g, uniformVal.b);

        }

    }

    this.renderBuffer = function ( camera, lights, fog, material, geometryGroup, object ) {

        if ( ! material.visible )
            return;

        var program, attributes, linewidth, primitives, a, attribute, i, il;

        //Sets up proper vertex and fragment shaders and attaches them to webGL program
        //Also sets appropriate uniform variables 
        program = setProgram(camera, lights, fog, material, object);

        attributes = program.attributes;

        var updateBuffers = false,
            wireframeBit = material.wireframe ? 1 : 0,
            geometryGroupHash = (geometryGroup.id * 0xffffff) + (program.id * 2) + wireframeBit;

        if (geometryGroupHash !== _currentGeometryGroupHash) {
            _currentGeometryGroupHash = geometryGroupHash;
            updateBuffers = true;
        }

        //rebind shader attributes to appropriate (and already initialized) gl buffers
        if (updateBuffers) {

            disableAttributes();

            // Vertices
            if (attributes.position >= 0) {            
                _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglVertexBuffer );
                enableAttribute( attributes.position );
                _gl.vertexAttribPointer( attributes.position, 3, _gl.FLOAT, false, 0, 0 );    
            }

            // Colors
            if (attributes.color >= 0) {
                _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglColorBuffer);
                enableAttribute( attributes.color );
                _gl.vertexAttribPointer( attributes.color, 3, _gl.FLOAT, false, 0, 0 );
            }

            // Normals (lambert shader only)
            if (attributes.normal >=0) {
                _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglNormalBuffer );
                enableAttribute( attributes.normal );
                _gl.vertexAttribPointer( attributes.normal, 3, _gl.FLOAT, false, 0, 0 );
            }

        }

        //Render
        var faceCount, lineCount;
        //lambert shaders - draw triangles
        //TODO: make sure geometryGroup's face count is setup correctly
        if (object instanceof WebMol.Mesh) {
            
            if (material.wireframe) {
                lineCount = geometryGroup.lineidx;
                setLineWidth(material.wireframeLinewidth);
                
                if (updateBuffers)
                    _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglLineBuffer );
                
                _gl.drawElements( _gl.LINES, lineCount, _gl.UNSIGNED_SHORT, 0 );
            }
            
            else {
                faceCount = geometryGroup.faceidx;

                if (updateBuffers)
                    _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglFaceBuffer );
                
                _gl.drawElements( _gl.TRIANGLES, faceCount, _gl.UNSIGNED_SHORT, 0 );
                
            }


            _this.info.render.calls++;
            _this.info.render.vertices += faceCount;
            _this.info.render.faces += faceCount / 3;
        }

        //basic shaders - draw lines
        else if (object instanceof WebMol.Line) {
            lineCount = geometryGroup.vertices;

            setLineWidth(material.linewidth);
            _gl.drawArrays( _gl.LINES, 0, lineCount );

            _this.info.render.calls++;
        }

    };

    //rendering
    function renderObjects ( renderList, reverse, materialType, camera, lights, fog, useBlending, overrideMaterial)  {

        var webglObject, object, buffer, material, start, end, delta;

        //Forward or backward render

        if (reverse) {
            start = renderList.length - 1;
            end = -1;
            delta = -1;
        }

        else {
            start = 0;
            end = renderList.length;
            delta = 1;
        }

        for (var i = start; i !== end; i += delta) {

            webglObject = renderList[i];

            if (webglObject.render) {

                object = webglObject.object;
                buffer = webglObject.buffer;
                material = webglObject[materialType];

                if ( ! material )
                    continue;

                if (useBlending)
                    _this.setBlending(material.blending);

                _this.setDepthTest(material.depthTest);
                _this.setDepthWrite(material.depthWrite);
                setPolygonOffset(material.polygonOffset, material.polygonOffsetFactor, material.polygonOffsetUnits);

                _this.setMaterialFaces(material);

                _this.renderBuffer(camera, lights, fog, material, buffer, object);
            }
        }

    }
    
    this.render = function ( scene, camera, renderTarget, forceClear ) {

        if ( camera instanceof WebMol.Camera === false )  {

            console.error( 'WebMol.Renderer.render: camera is not an instance of WebMol.Camera.' );
            return;

        }

        var i, il,

        webglObject, object,
        renderList,

        lights = scene.__lights,
        fog = scene.fog;

        // reset caching for this frame

        _currentMaterialId = -1;
        _lightsNeedUpdate = true;

        // update scene graph

        if ( this.autoUpdateScene ) scene.updateMatrixWorld();

        // update camera matrices
        //Pretty sure camera's parent is always going to be undefined for our purposes...
        if ( camera.parent === undefined ) camera.updateMatrixWorld();

        camera.matrixWorldInverse.getInverse( camera.matrixWorld );

        _projScreenMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

        // update WebGL objects

        if ( this.autoUpdateObjects ) this.initWebGLObjects( scene );


        _this.info.render.calls = 0;
        _this.info.render.vertices = 0;
        _this.info.render.faces = 0;
        _this.info.render.points = 0;

        _currentWidth = _viewportWidth;
        _currentHeight = _viewportHeight;

        if ( this.autoClear || forceClear ) {

            this.clear( this.autoClearColor, this.autoClearDepth, this.autoClearStencil );

        }

        // set matrices for regular objects (frustum culled)

        renderList = scene.__webglObjects;

        for ( i = 0, il = renderList.length; i < il; i ++ ) {

            webglObject = renderList[ i ];
            object = webglObject.object;

            webglObject.render = false;

            if ( object.visible ) {        
                setupMatrices( object, camera );
                unrollBufferMaterial( webglObject );
                webglObject.render = true;
            }
        }

        // set matrices for immediate objects

        var material = null;

        // opaque pass (front-to-back order)

        this.setBlending( WebMol.NoBlending );

        renderObjects( scene.__webglObjects, true, "opaque", camera, lights, fog, false, material );
        
        //prime depth buffer
        renderObjects( scene.__webglObjects, true, "blank", camera, lights, fog, true, material );

        // transparent pass (back-to-front order)

        renderObjects( scene.__webglObjects, false, "transparent", camera, lights, fog, true, material );

        // Render plugins (e.g. sprites), and reset state
        
        renderPlugins(this.renderPluginsPost, scene, camera);

        // Ensure depth buffer writing is enabled so it can be cleared on next render

        this.setDepthTest( true );
        this.setDepthWrite( true );

        //_gl.finish();

    };
    
    function renderPlugins(plugins, scene, camera) {
        
        //Reset state once regardless
        //This should also fix cartoon render bug (after transparent surface render)
        
        _currentGeometryGroupHash = -1;
        _currentProgram = null;
        _currentCamera = null;
        _oldBlending = -1;
        _oldDepthWrite = -1;
        _oldDepthTest = -1;
        _oldDoubleSided = -1;
        _currentMaterialId = -1;
        _oldFlipSided = -1;
        
        
        if (!plugins.length)
            return;
        
        for (var i = 0, il = plugins.length; i < il; i++) {
            
            _lightsNeedUpdate = true;
            
            plugins[i].render(scene, camera, _currentWidth, _currentHeight);
            
            //Reset state after plugin render
            _currentGeometryGroupHash = -1;
            _currentProgram = null;
            _currentCamera = null;
            _oldBlending = -1;
            _oldDepthWrite = -1;
            _oldDepthTest = -1;
            _oldDoubleSided = -1;
            _currentMaterialId = -1;
            _oldFlipSided = -1;       
                
        }  
        
    }

    this.initWebGLObjects = function ( scene ) {

        if ( !scene.__webglObjects ) {

            scene.__webglObjects = [];
            scene.__webglObjectsImmediate = [];
            scene.__webglSprites = [];
            scene.__webglFlares = [];

        }

        //Add objects; this sets up buffers for each geometryGroup
        if (scene.__objectsAdded.length) {
            
            while(scene.__objectsAdded.length){
                addObject(scene.__objectsAdded[0], scene);
                scene.__objectsAdded.splice(0, 1);
            }
            
            //Force buffer update during render
            //Hackish fix for initial cartoon-render-then-transparent-surface bug
            _currentGeometryGroupHash = -1;
            
        }

        while (scene.__objectsRemoved.length) {

            removeObject(scene.__objectsRemoved[ 0 ], scene);
            scene.__objectsRemoved.splice(0, 1);

        }

        // update must be called after objects adding / removal
        //This sends typed arrays to GL buffers for each geometryGroup
        for ( var o = 0, ol = scene.__webglObjects.length; o < ol; o ++ ) {

            updateObject(scene.__webglObjects[ o ].object);

        }

    };
    
    // Objects adding

    function addObject (object, scene) {

        var g, gl, geometry, material, geometryGroup;

        if ( !object.__webglInit ) {

            object.__webglInit = true;

            object._modelViewMatrix = new WebMol.Matrix4();
            object._normalMatrix = new WebMol.Matrix3();

            if (object.geometry !== undefined && object.geometry.__webglInit === undefined) {

                object.geometry.__webglInit = true;
                object.geometry.addEventListener('dispose', onGeometryDispose);

            }
            
            if (object instanceof WebMol.Mesh || object instanceof WebMol.Line) {
                geometry = object.geometry;
                material = object.material;           
    
                for (g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
    
                    geometryGroup = geometry.geometryGroups[ g ];
                    
                    geometryGroup.id = _geometryGroupCounter++;

                    // initialise VBO on the first access

                    if ( !geometryGroup.__webglVertexBuffer ) {
                            
                        if (object instanceof WebMol.Mesh) {
                            createMeshBuffers(geometryGroup);
                            geometry.elementsNeedUpdate = true;
                            geometry.normalsNeedUpdate = true;
                        }
                            
                        else if (object instanceof WebMol.Line)
                            createLineBuffers(geometryGroup);

                        geometry.verticesNeedUpdate = true;
                        geometry.colorsNeedUpdate = true;

                    }
                        
                }
                
            }
        
        }
        
        if ( ! object.__webglActive ) {
            
            if (object instanceof WebMol.Mesh || object instanceof WebMol.Line) {
                
                geometry = object.geometry;

                for ( g = 0, gl = geometry.geometryGroups.length; g < gl; g++ ) {
                    geometryGroup = geometry.geometryGroups[g];

                    addBuffer(scene.__webglObjects, geometryGroup, object);
                }
                
            }
            
            //Sprite
            else if (object instanceof WebMol.Sprite) 
                scene.__webglSprites.push(object);
         
                     
            object.__webglActive = true;
            
        }

    }

    function updateObject ( object ) {

        var geometry = object.geometry, material = object.material,
                geometryGroup, customAttributesDirty;
        
        if ( object instanceof WebMol.Mesh || object instanceof WebMol.Line ) {
            
            for (var g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
                
                geometryGroup = geometry.geometryGroups[ g ];

                if ( geometry.verticesNeedUpdate || geometry.elementsNeedUpdate || geometry.colorsNeedUpdate || geometry.normalsNeedUpdate) {
                    setBuffers( geometryGroup, _gl.DYNAMIC_DRAW );
                }
            }
            
            geometry.verticesNeedUpdate = false;
            geometry.elementsNeedUpdate = false;
            geometry.normalsNeedUpdate = false;
            geometry.colorsNeedUpdate = false;

            geometry.buffersNeedUpdate = false;

        }

    }
    
    function removeObject( object, scene ) {

        if (object instanceof WebMol.Mesh || object instanceof WebMol.Line )
            removeInstances(scene.__webglObjects, object);

        else if (object instanceof WebMol.Sprite)
            removeInstancesDirect(scene.__webglSprites, object);
            
        object.__webglActive = false;

    }

    function removeInstances( objList, object ) {

        for (var o = objList.length - 1; o >= 0; --o) {

            if (objList[o].object === object) 
                objList.splice(o, 1);

        }
    }

    function removeInstancesDirect( objList, object ) {

        for (var o = objList.length - 1; o >= 0; --o) {

            if (objList[o] === object) 
                objList.splice(o, 1);

        }
    }

    function unrollBufferMaterial( globject ) {

        var object = globject.object;
        var material = object.material;

        if ( material.transparent) {                    
            globject.opaque = null;
            globject.transparent = material;
            var blankMaterial = material.clone();
            blankMaterial.opacity = 0;
            globject.blank = blankMaterial;
        }

        else {
            globject.opaque = material;
            globject.transparent = null;
        }

    }

    function setBuffers( geometryGroup, hint, line ) {

        var vertexArray = geometryGroup.__vertexArray;
        var colorArray = geometryGroup.__colorArray;
         
        //vertex buffers
        _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglVertexBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, vertexArray, hint );        

        //color buffers
        _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglColorBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, colorArray, hint );    
              
        
        //normal buffers
        if (geometryGroup.__normalArray !== undefined && geometryGroup.__webglNormalBuffer !== undefined) {
            var normalArray = geometryGroup.__normalArray;
            _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglNormalBuffer );
            _gl.bufferData( _gl.ARRAY_BUFFER, normalArray, hint );       
             
        }
        
        //face (index) buffers
        if (geometryGroup.__faceArray !== undefined && geometryGroup.__webglFaceBuffer !== undefined) {
            var faceArray = geometryGroup.__faceArray;
            _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglFaceBuffer );
            _gl.bufferData( _gl.ELEMENT_ARRAY_BUFFER, faceArray, hint );  
                      
        }
        
        //line (index) buffers (for wireframe)
        if (geometryGroup.__lineArray !== undefined && geometryGroup.__webglLineBuffer !== undefined) {
            var lineArray = geometryGroup.__lineArray;            
            _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglLineBuffer );
            _gl.bufferData( _gl.ELEMENT_ARRAY_BUFFER, lineArray, hint );
        }

    }
    
    //Creates appropriate gl buffers for geometry chunk
    //TODO: do we need line buffer for mesh objects?
    //Also, can we integrate this with createLineBuffers?
    function createMeshBuffers ( geometryGroup ) {

        geometryGroup.__webglVertexBuffer = _gl.createBuffer();
        geometryGroup.__webglNormalBuffer = _gl.createBuffer();
        geometryGroup.__webglColorBuffer = _gl.createBuffer();

        geometryGroup.__webglFaceBuffer = _gl.createBuffer();
        geometryGroup.__webglLineBuffer = _gl.createBuffer();

        _this.info.memory.geometries++;
    }
    
    function createLineBuffers ( geometry ) {
        
        geometry.__webglVertexBuffer = _gl.createBuffer();
        geometry.__webglColorBuffer = _gl.createBuffer();
        
        _this.info.memory.geometries++;
    }

    function addBuffer (objlist, buffer, object) {

        objlist.push(
            {
                buffer: buffer,
                object: object,
                opaque: null,
                transparent: null
            }
        );

    }

    function setupMatrices (object, camera) {

        object._modelViewMatrix.multiplyMatrices( camera.matrixWorldInverse, object.matrixWorld );

        object._normalMatrix.getInverse( object._modelViewMatrix );
        object._normalMatrix.transpose();

    }

    function isPowerOfTwo ( value ) {

        return ( value & ( value - 1 ) ) === 0;

    }
    
    // Fallback filters for non-power-of-2 textures

    function filterFallback ( f ) {

        return _gl.LINEAR;

    }

    function setTextureParameters ( textureType, texture, isImagePowerOfTwo ) {

        if ( isImagePowerOfTwo ) {

            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, paramToGL( texture.wrapS ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, paramToGL( texture.wrapT ) );

            _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, paramToGL( texture.magFilter ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, paramToGL( texture.minFilter ) );

        } else {

            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE );
            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE );

            _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, filterFallback( texture.magFilter ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, filterFallback( texture.minFilter ) );

        }

    }
    
    this.setTexture = function (texture, slot) {

        if (texture.needsUpdate) {

            if ( !texture.__webglInit ) {

                texture.__webglInit = true;

                texture.addEventListener('dispose', onTextureDispose);

                texture.__webglTexture = _gl.createTexture();

                _this.info.memory.textures++;

            }

            _gl.activeTexture(_gl.TEXTURE0 + slot);
            _gl.bindTexture(_gl.TEXTURE_2D, texture.__webglTexture);

            _gl.pixelStorei(_gl.UNPACK_FLIP_Y_WEBGL, texture.flipY);
            _gl.pixelStorei(_gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, texture.premultiplyAlpha);
            _gl.pixelStorei(_gl.UNPACK_ALIGNMENT, texture.unpackAlignment);

            var image = texture.image,
            isImagePowerOfTwo = isPowerOfTwo(image.width) && isPowerOfTwo(image.height),
            glFormat = paramToGL(texture.format),
            glType = paramToGL(texture.type);

            setTextureParameters(_gl.TEXTURE_2D, texture, isImagePowerOfTwo);

            var mipmap, mipmaps = texture.mipmaps;

            // regular Texture (image, video, canvas)

            // use manually created mipmaps if available
            // if there are no manual mipmaps
            // set 0 level mipmap and then use GL to generate other mipmap levels

            if ( mipmaps.length > 0 && isImagePowerOfTwo ) {

                for ( var i = 0, il = mipmaps.length; i < il; i ++ ) {
                    mipmap = mipmaps[ i ];
                    _gl.texImage2D( _gl.TEXTURE_2D, i, glFormat, glFormat, glType, mipmap );
                }
                
                texture.generateMipmaps = false;
            } 
            
            else 
                _gl.texImage2D( _gl.TEXTURE_2D, 0, glFormat, glFormat, glType, texture.image );

            
            if ( texture.generateMipmaps && isImagePowerOfTwo ) _gl.generateMipmap( _gl.TEXTURE_2D );

            texture.needsUpdate = false;

            if ( texture.onUpdate ) texture.onUpdate();

        } else {

            _gl.activeTexture( _gl.TEXTURE0 + slot );
            _gl.bindTexture( _gl.TEXTURE_2D, texture.__webglTexture );

        }

    };
    
    // Map constants to WebGL constants

    function paramToGL ( p ) {

        if ( p === WebMol.UnsignedByteType ) return _gl.UNSIGNED_BYTE;
        if ( p === WebMol.RGBAFormat ) return _gl.RGBA;

        return 0;

    }
    
    function setupLights ( program, lights ) {
        var l, ll, light, n,
        r = 0, g = 0, b = 0,
        color,
        position,
        intensity,
        distance,
        
        zlights = _lights,
        
        dirColors = zlights.directional.colors,
        dirPositions = zlights.directional.positions,
        
        dirCount = 0,
        dirLength = 0,
        dirOffset = 0;
        
        for ( l = 0, ll = lights.length; l < ll; l++) {
            
            light = lights[l];
            
            color = light.color;
            intensity = light.intensity;
            distance = light.distance;
            
            if (light instanceof WebMol.Light) {
                
                dirCount++;
                
                _direction.getPositionFromMatrix(light.matrixWorld);
                _vector3.getPositionFromMatrix(light.target.matrixWorld);
                _direction.sub(_vector3);
                _direction.normalize();
                
                if (_direction.x === 0 && _direction.y === 0 && _direction.z === 0)
                    continue;
                
                dirPositions[dirOffset] = _direction.x;
                dirPositions[dirOffset + 1] = _direction.y;
                dirPositions[dirOffset + 2] = _direction.z;

                dirColors[dirOffset] = color.r * intensity;
                dirColors[dirOffset + 1] = color.g * intensity;
                dirColors[dirOffset + 2] = color.b * intensity;
                
                dirOffset += 3;
                
                dirLength++;    
            }
        
        }

        zlights.ambient[0] = r;
        zlights.ambient[1] = g;
        zlights.ambient[2] = b;
        zlights.directional.length = dirLength;
    }

    function initGL () {

        try {

            if ( ! ( _gl = _canvas.getContext( 'experimental-webgl', { alpha: _alpha, premultipliedAlpha: _premultipliedAlpha, antialias: _antialias, stencil: _stencil, preserveDrawingBuffer: _preserveDrawingBuffer } ) ) ) {

                throw 'Error creating WebGL context.';

            }

        } catch ( error ) {

            console.error( error );

        }

    }

    function setDefaultGLState () {

        _gl.clearColor( 0, 0, 0, 1 );
        _gl.clearDepth( 1 );
        _gl.clearStencil( 0 );

        _gl.enable( _gl.DEPTH_TEST );
        _gl.depthFunc( _gl.LEQUAL );

        _gl.frontFace( _gl.CCW );
        _gl.cullFace( _gl.BACK );
        _gl.enable( _gl.CULL_FACE );

        _gl.enable( _gl.BLEND );
        _gl.blendEquation( _gl.FUNC_ADD );
        _gl.blendFunc( _gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA );

        _gl.clearColor( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha );

    }
    
    this.addPostPlugin(new WebMol.SpritePlugin());
        
};

/*
 * Scene class
 */

WebMol.Scene = function() {
    
    WebMol.Object3D.call(this);
    
    this.fog = null;
    
    //May not need...
    this.overrideMaterial = null;
    
    this.matrixAutoUpdate = false;
    
    this.__objects = [];
    this.__lights = [];
    
    this.__objectsAdded = [];
    this.__objectsRemoved = [];
    
};

WebMol.Scene.prototype = Object.create(WebMol.Object3D.prototype);

WebMol.Scene.prototype.__addObject = function(object) {
    
    var i;
    
    //Directional Lighting
    if (object instanceof WebMol.Light) {
        
        if (this.__lights.indexOf(object) === -1)
            this.__lights.push(object);
        
        //TODO: Do I need this??
        if (object.target && object.target.parent === undefined)
            this.add(object.target);
            
    }
    
    //Rotation group
    else {
        
        if (this.__objects.indexOf(object) === -1) {
            
            this.__objects.push(object);
            this.__objectsAdded.push(object);
            
            //Check if previously removed
            
            i = this.__objectsRemoved.indexOf(object);
            
            if (i !== -1)
                this.__objectsRemoved.splice(i, 1);
                
        }
    }
    
    //Add object's children
    
    for (i= 0, il = object.children.length; i < il; i++) 
        this.__addObject(object.children[i]);
    
};

WebMol.Scene.prototype.__removeObject = function(object) {
    
    var i, il;
    
    if (object instanceof WebMol.Light) {
        
        i = this.__lights.indexOf(object);
        
        if (i !== -1)
            this.__lights.splice(i, 1);
            
    }
    
    //Object3D
    else {
        
        i = this.__objects.indexOf(object);
        
        if (i !== -1) {
            
            this.__objects.splice(i, 1);
            this.__objectsRemoved.push(object);
            
            //Check if previously added
            
            var ai = this.__objectsAdded.indexOf(object);
            
            if (ai !== -1) 
                this.__objectsAdded.splice(i, 1);
                
        }
    
    }
    
    //Remove object's children
    for (i = 0, il = object.children.length; i < il; i++)
        this.__removeObject(object.children[i]);
    
};


/*
 * Fog Class
 */


WebMol.Fog = function ( hex, near, far ) {

    this.name = '';

    this.color = new WebMol.Color( hex );

    this.near = ( near !== undefined ) ? near : 1;
    this.far = ( far !== undefined ) ? far : 1000;

};

WebMol.Fog.prototype.clone = function () {

    return new WebMol.Fog( this.color.getHex(), this.near, this.far );

};/* 
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

WebMol.ShaderUtils = {
    
    clone: function ( uniforms_src ) {
        
        var u, p, parameter, parameter_src, uniforms_clone = {};
        
        for (u in uniforms_src) {
            uniforms_clone[u] = {};
            uniforms_clone[u].type = uniforms_src[u].type;
            
            var srcValue = uniforms_src[u].value;
            
            if (srcValue instanceof WebMol.Color)
                uniforms_clone[u].value = srcValue.clone();
            else if (typeof srcValue === "number")
                uniforms_clone[u].value = srcValue;
            else if (srcValue instanceof Array) 
                uniforms_clone[u].value = [];
            else
                console.error("Error copying shader uniforms from ShaderLib: unknown type for uniform");
            
        }
        
        return uniforms_clone;
    }
};

WebMol.ShaderLib = { 
    basic : {
        fragmentShader : [
"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 diffuse;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vColor;",

"void main() {",
    
"    gl_FragColor = vec4( diffuse, opacity );",
"    gl_FragColor = gl_FragColor * vec4( vColor, opacity );",
    
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",    
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
    
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"
                                                     
].join("\n"),
        
        vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",

"attribute vec3 position;",
"attribute vec3 color;",

"varying vec3 vColor;",

"void main() {",

"    vColor = color;",
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
"    gl_Position = projectionMatrix * mvPosition;",

"}"
        
].join("\n"),
    
        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            diffuse: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogColor: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000}
        }

    },
    
    
    lambert : { 
        fragmentShader : [

"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vLightFront;",
"varying vec3 vColor;",

"void main() {",
    
"    gl_FragColor = vec4( vec3 ( 1.0 ), opacity );",
    
"    #ifndef WIREFRAME",
"    gl_FragColor.xyz *= vLightFront;",
"    #endif",
    
"    gl_FragColor = gl_FragColor * vec4( vColor, opacity );",
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",
    
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
    
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"


].join("\n"),
       
       vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 ambient;",
"uniform vec3 diffuse;",
"uniform vec3 emissive;",
"uniform vec3 ambientLightColor;",
"uniform vec3 directionalLightColor[ 1 ];",
"uniform vec3 directionalLightDirection[ 1 ];",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"varying vec3 vColor;",
"varying vec3 vLightFront;",

"void main() {",
    
"    vColor = color;",
    
"    vec3 objectNormal = normal;",  
"    vec3 transformedNormal = normalMatrix * objectNormal;",    
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
    
"    vLightFront = vec3( 0.0 );",
    
"    transformedNormal = normalize( transformedNormal );",
    
"    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );",
"    vec3 dirVector = normalize( lDirection.xyz );",
"    float dotProduct = dot( transformedNormal, dirVector );",
"    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );",
    
"    vLightFront += directionalLightColor[ 0 ] * directionalLightWeighting;",
"    vLightFront = vLightFront * diffuse + ambient * ambientLightColor + emissive;",
    
"    gl_Position = projectionMatrix * mvPosition;",
"}"
           
].join("\n"),

        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            diffuse: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogColor: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000},           
            ambient: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            emissive: { type: 'c', value: new WebMol.Color(1.0, 1.0, 1.0) },
            ambientLightColor: { type: 'fv', value: [] },
            directionalLightColor: { type: 'fv', value: [] },
            directionalLightDirection: { type: 'fv', value: [] }
        }

    },
    
    
    sprite: {
        
        fragmentShader: [
                                                         
"uniform vec3 color;",
"uniform sampler2D map;",
"uniform float opacity;",

"uniform int fogType;",
"uniform vec3 fogColor;",
"uniform float fogDensity;",
"uniform float fogNear;",
"uniform float fogFar;",
"uniform float alphaTest;",

"varying vec2 vUV;",

"void main() {",
    
"    vec4 texture = texture2D(map, vUV);",
    
"    if (texture.a < alphaTest) discard;",
    
"    gl_FragColor = vec4(color * texture.xyz, texture.a * opacity);",
    
"    if (fogType > 0) {",
        
"        float depth = gl_FragCoord.z / gl_FragCoord.w;",
"        float fogFactor = 0.0;",
        
"        if (fogType == 1) {",
"            fogFactor = smoothstep(fogNear, fogFar, depth);",
"        }",
        
"        else {",
"            const float LOG2 = 1.442695;",
"            float fogFactor = exp2(- fogDensity * fogDensity * depth * depth * LOG2);",
"            fogFactor = 1.0 - clamp(fogFactor, 0.0, 1.0);",
"        }",
        
"        gl_FragColor = mix(gl_FragColor, vec4(fogColor, gl_FragColor.w), fogFactor);",
        
"    }",
"}"                                              
            
].join("\n"),
        
        vertexShader: [

"uniform int useScreenCoordinates;",
"uniform int sizeAttenuation;",     
"uniform vec3 screenPosition;",
"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform float rotation;",
"uniform vec2 scale;",
"uniform vec2 alignment;",
"uniform vec2 uvOffset;",
"uniform vec2 uvScale;",

"attribute vec2 position;",
"attribute vec2 uv;",

"varying vec2 vUV;",

"void main() {",
    
"    vUV = uvOffset + uv * uvScale;",
    
"    vec2 alignedPosition = position + alignment;",
    
"    vec2 rotatedPosition;",
"    rotatedPosition.x = ( cos(rotation) * alignedPosition.x - sin(rotation) * alignedPosition.y ) * scale.x;",
"    rotatedPosition.y = ( sin(rotation) * alignedPosition.x + cos(rotation) * alignedPosition.y ) * scale.y;",
    
"    vec4 finalPosition;",
    
"    if(useScreenCoordinates != 0) {",
"        finalPosition = vec4(screenPosition.xy + rotatedPosition, screenPosition.z, 1.0);",
"    }",
    
"    else {",
"        finalPosition = projectionMatrix * modelViewMatrix * vec4(0.0, 0.0, 0.0, 1.0);",
"        finalPosition.xy += rotatedPosition * (sizeAttenuation == 1 ? 1.0 : finalPosition.z);",
"    }",
    
"    gl_Position = finalPosition;",
    
"}"
       
].join("\n"),

        uniforms : {
            
        }
        
    }
    
};
//Hackish way to create webworker (independent of WebMol namespace) within minified file
WebMol.workerString = function(){

    self.onmessage = function(oEvent) {
    	var obj = oEvent.data;
    	var type = obj.type;
    	if (type < 0) // sending atom data, initialize
    	{
    		self.atomData = obj.atoms;
    		self.volume = obj.volume;
    		self.ps = new ProteinSurface();
    	} else {
    		var ps = self.ps;
    		ps.initparm(obj.expandedExtent, (type == 1) ? false : true, self.volume);
    		ps.fillvoxels(self.atomData, obj.extendedAtoms);
    		ps.buildboundary();
    		if (type === 4 || type === 2) {
    			ps.fastdistancemap();
                ps.boundingatom(false);
                ps.fillvoxelswaals(self.atomData, obj.extendedAtoms);	
            }		
    		ps.marchingcube(type);
    		var VandF = ps.getFacesAndVertices(obj.atomsToShow);
    		self.postMessage(VandF);
    	}
    };
    
    var ISDONE = 2;
    
    var Vector3 = function(x, y, z) {
        this.x = x || 0.0;
        this.y = y || 0.0;
        this.z = z || 0.0;
    };    
    
    Vector3.prototype =  {
        
        constructor : Vector3,
        
        copy : function(v) {
            
            this.x = v.x;
            this.y = v.y;
            this.z = v.z;
            
            return this;  
        },
        
        multiplyScalar : function(s) {
            
            this.x *= s;
            this.y *= s;
            this.z *= s;
            
            return this;
        }
    }; 
    
}.toString().replace(/(^.*?\{|\}$)/g, "");

WebMol.workerString += ",ISDONE=2";
WebMol.workerString += ",ProteinSurface=" + WebMol.ProteinSurface.toString().replace(/WebMol.MarchingCube./g, "");
WebMol.workerString += ",march=" + WebMol.MarchingCube.march.toString().replace(/WebMol./g, "");
WebMol.workerString += ",laplacianSmooth=" + WebMol.MarchingCube.laplacianSmooth.toString();

WebMol.workerString += ",edgeTable=new Uint32Array([" + WebMol.MarchingCube.edgeTable.toString() + "])";

WebMol.workerString += ",triTable=[";

for (var i = 0, il = WebMol.MarchingCube.triTable.length; i < il - 1; i++)
    WebMol.workerString += "[" + WebMol.MarchingCube.triTable[i].toString() + "],";

WebMol.workerString += "[]]";


WebMol.SurfaceWorker = window.URL.createObjectURL(new Blob([WebMol.workerString]));

 
//auto-initialization
//Create embedded viewer from HTML attributes if true

$(document).ready(function() {

    if ($("#webmoljs_viewer")[0] !== undefined)
        WebMol.autoinit = true;
        
    if (WebMol.autoinit) { 
  
        var viewerdiv = WebMol.viewerdiv = $("#webmoljs_viewer");
        
        var datauri = null;
        
        if (viewerdiv.data("pdb"))
            datauri = "http://www.pdb.org/pdb/files/" + viewerdiv.data("pdb") + ".pdb";
        else if (viewerdiv.data("href"))
            datauri = viewerdiv.data("href");
            
        var bgcolor = Number(viewerdiv.data("backgroundcolor")) || 0x000000;
        var style = viewerdiv.data("style") || {line:{}};
        
        WebMol.glviewer = WebMol.createViewer("webmoljs_viewer", {defaultcolors: WebMol.rasmolElementColors, callback: function(viewer) {            
            viewer.setBackgroundColor(bgcolor);            
        }});
        
        if (datauri) {  
            
            var type = viewerdiv.data("datatype") || "pdb";
             
            $.get(datauri, function(ret) {
                WebMol.glviewer.addModel(ret, type);
                WebMol.glviewer.setStyle({}, style);
                WebMol.glviewer.zoomTo();
                WebMol.glviewer.render();                                           
            }, 'text');
       
        }
               
    }
});
    
var WebMol = WebMol || {};
WebMol.defaultElementColor = 0xff1493;

WebMol.JmolElementColors = {
		H: 0xFFFFFF,
		He: 0xD9FFFF,
		HE: 0xD9FFFF,
		Li: 0xCC80FF,
		LI: 0xCC80FF,
		B: 0xFFB5B5,
		C: 0x909090,
		N: 0x3050F8,
		O: 0xFF0D0D,
		F: 0x90E050,
		Na: 0xAB5CF2,
		NA: 0xAB5CF2,
		Mg: 0x8AFF00,
		MG: 0x8AFF00,
		Al: 0xBFA6A6,
		AL: 0xBFA6A6,
		Si: 0xF0C8A0,
		SI: 0xF0C8A0,
		P: 0xFF8000,
		S: 0xFFFF30,
		Cl: 0x1FF01F,
		CL: 0x1FF01F,
		Ca: 0x3DFF00,
		CA: 0x3DFF00,
		Ti: 0xBFC2C7,
		TI: 0xBFC2C7,
		Cr: 0x8A99C7,
		CR: 0x8A99C7,
		Mn: 0x9C7AC7,
		MN: 0x9C7AC7,
		Fe: 0xE06633,
		FE: 0xE06633,
		Ni: 0x50D050,
		NI: 0x50D050,
		Cu: 0xC88033,
		CU: 0xC88033,
		Zn: 0x7D80B0,
		ZN: 0x7D80B0,
		Br: 0xA62929,
		BR: 0xA62929,
		Ag: 0xC0C0C0,
		AG: 0xC0C0C0,
		I: 0x940094,
		Ba: 0x00C900,
		BA: 0x00C900,
		Au: 0xFFD123,
		AU: 0xFFD123
};

WebMol.rasmolElementColors = {
		H: 0xFFFFFF,
		He: 0xFFC0CB,
		HE: 0xFFC0CB,
		Li: 0xB22222,
		LI: 0xB22222,
		B: 0x00FF00,
		C: 0xC8C8C8,
		N: 0x8F8FFF,
		O: 0xF00000,
		F: 0xDAA520,
		Na: 0x0000FF,
		NA: 0x0000FF,
		Mg: 0x228B22,
		MG: 0x228B22,
		Al: 0x808090,
		AL: 0x808090,
		Si: 0xDAA520,
		SI: 0xDAA520,
		P: 0xFFA500,
		S: 0xFFC832,
		Cl: 0x00FF00,
		CL: 0x00FF00,
		Ca: 0x808090,
		CA: 0x808090,
		Ti: 0x808090,
		TI: 0x808090,
		Cr: 0x808090,
		CR: 0x808090,
		Mn: 0x808090,
		MN: 0x808090,
		Fe: 0xFFA500,
		FE: 0xFFA500,
		Ni: 0xA52A2A,
		NI: 0xA52A2A,
		Cu: 0xA52A2A,
		CU: 0xA52A2A,
		Zn: 0xA52A2A,
		ZN: 0xA52A2A,
		Br: 0xA52A2A,
		BR: 0xA52A2A,
		Ag: 0x808090,
		AG: 0x808090,
		I: 0xA020F0,
		Ba: 0xFFA500,
		BA: 0xFFA500,
		Au: 0xDAA520,
		AU: 0xDAA520	
};

WebMol.defaultElementColors = WebMol.JmolElementColors;//color scheme mappings
var WebMol = WebMol || {};

//red to white to blue, for charges
WebMol.RWB = function(min, max) {
	
	//map value to hex color, range is provided
	this.valueToHex = function(val, range) {
		var lo, hi;
		if(range) {
			lo = range[0];
			hi = range[1];
		}
		else {
			lo = min;
			hi = max;
		}
	
		if(val === undefined)
			return 0xffffff;
		
		if(val < lo) val = lo;
		if(val > hi) val = hi;
		
		var middle = (hi+lo)/2;
		var scale, color;
		
		//scale bottom from red to white
		if(val <= middle) {
			scale = Math.floor(255*Math.sqrt((val-lo)/(middle-lo)));
			color = 0xff0000 + 0x100*scale + scale;
			return color;
		}
		else { //form white to blue
			scale = Math.floor(255*Math.sqrt((1-(val-middle)/(hi-middle))));
			color =  0x10000*scale+0x100*scale+0xff;
			return color;
		}
	};
	
	this.jmolID = function() {
		return "rwb";
	};

	//return range used for color mapping, null if none set
	this.range = function() {
		if(typeof(min) != "undefined" && typeof(max) != "undefined") {
			return [min,max];
		}
		return null;
	};

};

//rainbow gradient, but without purple to match jmol
WebMol.ROYGB = function(min, max) {
	
	//map value to hex color, range is provided
	this.valueToHex = function(val, range) {
		var lo, hi;
		if(range) {
			lo = range[0];
			hi = range[1];
		}
		else {
			lo = min;
			hi = max;
		}
	
		if(typeof(val) == "undefined")
			return 0xffffff;
		
		if(val < lo) val = lo;
		if(val > hi) val = hi;
		
		var mid = (lo+hi)/2;
		var q1 = (lo+mid)/2;
		var q3 = (mid+hi)/2;
		
		var scale, color;
		
		if(val < q1) { //scale green up, red up, blue down
			scale = Math.floor(255*Math.sqrt((val-lo)/(q1-lo)));
			color = 0xff0000 + 0x100*scale + 0;
			return color;
		}
		else if(val < mid) { //scale red down, green up, blue down
			scale = Math.floor(255*Math.sqrt((1-(val-q1)/(mid-q1))));
			color =  0x010000*scale+0xff00+0x0;
			return color;
		}
		else if(val < q3) { //scale blue up, red down, green up
			scale = Math.floor(255*Math.sqrt((val-mid)/(q3-mid)));
			color = 0x000000 + 0xff00 + 0x1*scale;
			return color;
		}
		else { //scale green down, blue up, red down
			scale = Math.floor(255*Math.sqrt((1-(val-q3)/(hi-q3))));
			color =  0x000000+0x0100*scale+0xff;
			return color;
		}		
	};
	
	this.jmolID = function() {
		return "roygb";
	};

	//return range used for color mapping, null if none set
	this.range = function() {
		if(typeof(min) != "undefined" && typeof(max) != "undefined") {
			return [min,max];
		}
		return null;
	};

};

//rainbow gradient with constant saturation, all the way to purple!
WebMol.Sinebow = function(min, max) {
	
	//map value to hex color, range is provided
	this.valueToHex = function(val, range) {
		var lo, hi;
		if(range) {
			lo = range[0];
			hi = range[1];
		}
		else {
			lo = min;
			hi = max;
		}
	
		if(typeof(val) == "undefined")
			return 0xffffff;
		
		if(val < lo) val = lo;
		if(val > hi) val = hi;
		
		var scale = (val-lo)/(hi-lo);
		var h = (5*scale/6.0+0.5);
		var r = Math.sin(Math.PI*h);
		r *= r*255;
		var g = Math.sin(Math.PI*(h+1/3.0));
		g *= g*255;
		var b = Math.sin(Math.PI*(h+2/3.0));
		b *= b*255;
		
		return 0x10000*Math.floor(r)+0x100*Math.floor(b)+0x1*Math.floor(g);
	};
	
	this.jmolID = function() {
		return "sinebow";
	};

	//return range used for color mapping, null if none set
	this.range = function() {
		if(typeof(min) != "undefined" && typeof(max) != "undefined") {
			return [min,max];
		}
		return null;
	};

};
//glcartoon.js
//This contains all the routines for rendering a cartoon given a set
//of atoms with assigned secondary structure
//TODO: secondary structure calculation

//TODO: generate normals directly in drawStrip and drawThinStrip

var WebMol = WebMol || {};

WebMol.drawCartoon = (function() {

    var axisDIV = 5; // 3 still gives acceptable quality
    var strandDIV = 6;
    var nucleicAcidStrandDIV = 4;
    var tubeDIV = 8;
    var coilWidth = 0.3;
    var helixSheetWidth = 1.3;
    var nucleicAcidWidth = 0.8;
    var thickness = 0.4; 

    // helper functions

    // Catmull-Rom subdivision
    var subdivide = function(_points, DIV) { // points as Vector3
        var ret = [];
        var points = _points;
        points = new Array(); // Smoothing test
        points.push(_points[0]);
        for ( var i = 1, lim = _points.length - 1; i < lim; i++) {
            var p1 = _points[i], p2 = _points[i + 1];
            if (p1.smoothen)
                points.push(new WebMol.Vector3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2,
                        (p1.z + p2.z) / 2));
            else
                points.push(p1);
        }
        points.push(_points[_points.length - 1]);

        
        for ( var i = -1, size = points.length; i <= size - 3; i++) {
            var p0 = points[(i === -1) ? 0 : i];
            var p1 = points[i + 1], p2 = points[i + 2];
            var p3 = points[(i === size - 3) ? size - 1 : i + 3];
            var v0 = new WebMol.Vector3().subVectors(p2, p0).multiplyScalar(0.5);
            var v1 = new WebMol.Vector3().subVectors(p3, p1).multiplyScalar(0.5);

            for ( var j = 0; j < DIV; j++) {
                var t = 1.0 / DIV * j;
                var x = p1.x + t * v0.x + t * t
                        * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x) + t * t * t
                        * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
                var y = p1.y + t * v0.y + t * t
                        * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y) + t * t * t
                        * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
                var z = p1.z + t * v0.z + t * t
                        * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z) + t * t * t
                        * (2 * p1.z - 2 * p2.z + v0.z + v1.z);
                        
                var pt = new WebMol.Vector3(x, y, z);
                
                var atomIndex = Math.floor( (ret.length+2) / DIV);
                
                if (_points[atomIndex] !== undefined && _points[atomIndex].atom !== undefined)
                    pt.atom = _points[atomIndex].atom;
                    
                ret.push(pt);
            }
        }
        ret.push(points[points.length - 1]);
        return ret;
    };

    var drawThinStrip = function(group, p1, p2, colors, div) {
    
        var geo = new WebMol.Geometry(true);       
        var offset, vertoffset;
        var color;
        
        for ( var i = 0, lim = p1.length; i < lim; i++) {
            
            color = WebMol.CC.color(colors[Math.round((i - 1) / div)]);
           
            geoGroup = geo.updateGeoGroup(2);
            offset = geoGroup.vertices, vertoffset = offset*3;
            
            geoGroup.__vertexArray[vertoffset] = p1[i].x;
            geoGroup.__vertexArray[vertoffset+1] = p1[i].y;
            geoGroup.__vertexArray[vertoffset+2] = p1[i].z;
            
            geoGroup.__vertexArray[vertoffset+3] = p2[i].x;
            geoGroup.__vertexArray[vertoffset+4] = p2[i].y;
            geoGroup.__vertexArray[vertoffset+5] = p2[i].z;
            
            for (var j = 0; j < 6; ++j) {                
                geoGroup.__colorArray[vertoffset+3*j] = color.r; geoGroup.__colorArray[vertoffset+1+3*j] = color.g; geoGroup.__colorArray[vertoffset+2+3*j] = color.b;                
            }            
           
            if (i > 0) {
                var faces = [offset, offset + 1, offset - 1, offset - 2];
                var faceoffset = geoGroup.faceidx;
                
                geoGroup.__faceArray[faceoffset] = faces[0]; geoGroup.__faceArray[faceoffset+1] = faces[1]; geoGroup.__faceArray[faceoffset+2] = faces[3];
                geoGroup.__faceArray[faceoffset+3] = faces[1]; geoGroup.__faceArray[faceoffset+4] = faces[2]; geoGroup.__faceArray[faceoffset+5] = faces[3];

                geoGroup.faceidx += 6;
            }
            
            geoGroup.vertices += 2;
        }
        
        geo.initTypedArrays();
        geo.setUpNormals();
        
        var material = new WebMol.MeshLambertMaterial();
        material.vertexColors = WebMol.FaceColors;
                material.side = WebMol.DoubleSide;
        var mesh = new WebMol.Mesh(geo, material);
        group.add(mesh);
    };

    var drawStrip = function(group, p1, p2, colors, div, thickness) {
        if ((p1.length) < 2)
            return;
        div = div || axisDIV;
        p1 = subdivide(p1, div);
        p2 = subdivide(p2, div);
        if (!thickness)
            return drawThinStrip(group, p1, p2, colors, div);

        var geo = new WebMol.Geometry(true);
        
        //var vs = geo.vertices, fs = geo.faces;
                var vs = [], fs = [];
        var axis, p1v, p2v, a1v, a2v;
        
        var faces = [ [ 0, 2, -6, -8 ], [ -4, -2, 6, 4 ], [ 7, -1, -5, 3 ],
                [ -3, 5, 1, -7 ] ];
                
        var offset, vertoffset, faceoffset;
        var color;
        var currentAtom, lastAtom;
        
        for ( var i = 0, lim = p1.length; i < lim; i++) {
        
            color = WebMol.CC.color(colors[Math.round((i - 1) / div)]);
            
            vs.push(p1v = p1[i]); // 0
            vs.push(p1v); // 1
            vs.push(p2v = p2[i]); // 2
            vs.push(p2v); // 3
            if (i < lim - 1) {
                var toNext = p1[i + 1].clone().sub(p1[i]);
                var toSide = p2[i].clone().sub(p1[i]);
                axis = toSide.cross(toNext).normalize().multiplyScalar(
                        thickness);
            }
            vs.push(a1v = p1[i].clone().add(axis)); // 4
            vs.push(a1v); // 5
            vs.push(a2v = p2[i].clone().add(axis)); // 6
            vs.push(a2v); // 7
            
            if (p1v.atom !== undefined)
                currentAtom = p1v.atom;
            
            var geoGroup = geo.updateGeoGroup(8);
            offset = geoGroup.vertices, vertoffset = offset*3;
            
            geoGroup.__vertexArray[vertoffset] = p1v.x, geoGroup.__vertexArray[vertoffset+1] = p1v.y, geoGroup.__vertexArray[vertoffset+2] = p1v.z;
            geoGroup.__vertexArray[vertoffset+3] = p1v.x, geoGroup.__vertexArray[vertoffset+4] = p1v.y, geoGroup.__vertexArray[vertoffset+5] = p1v.z;
            geoGroup.__vertexArray[vertoffset+6] = p2v.x, geoGroup.__vertexArray[vertoffset+7] = p2v.y, geoGroup.__vertexArray[vertoffset+8] = p2v.z;
            geoGroup.__vertexArray[vertoffset+9] = p2v.x, geoGroup.__vertexArray[vertoffset+10] = p2v.y, geoGroup.__vertexArray[vertoffset+11] = p2v.z;
            geoGroup.__vertexArray[vertoffset+12] = a1v.x, geoGroup.__vertexArray[vertoffset+13] = a1v.y, geoGroup.__vertexArray[vertoffset+14] = a1v.z;
            geoGroup.__vertexArray[vertoffset+15] = a1v.x, geoGroup.__vertexArray[vertoffset+16] = a1v.y, geoGroup.__vertexArray[vertoffset+17] = a1v.z;
            geoGroup.__vertexArray[vertoffset+18] = a2v.x, geoGroup.__vertexArray[vertoffset+19] = a2v.y, geoGroup.__vertexArray[vertoffset+20] = a2v.z;
            geoGroup.__vertexArray[vertoffset+21] = a2v.x, geoGroup.__vertexArray[vertoffset+22] = a2v.y, geoGroup.__vertexArray[vertoffset+23] = a2v.z;
            
            for (var j = 0; j < 8; ++j) {                
                geoGroup.__colorArray[vertoffset+3*j] = color.r; geoGroup.__colorArray[vertoffset+1+3*j] = color.g; geoGroup.__colorArray[vertoffset+2+3*j] = color.b;                
            }
            
            if (i > 0) {
             
                //both points have distinct atoms
                var diffAtoms = ((lastAtom !== undefined && currentAtom !== undefined) && lastAtom.serial !== currentAtom.serial);
                
                for ( var j = 0; j < 4; j++ ) {
                
                    var face = [offset + faces[j][0], offset
                        + faces[j][1], offset + faces[j][2], offset
                        + faces[j][3]];
                    
                    faceoffset = geoGroup.faceidx;    
                    
                    geoGroup.__faceArray[faceoffset] = face[0], geoGroup.__faceArray[faceoffset+1] = face[1], geoGroup.__faceArray[faceoffset+2] = face[3];             
                    geoGroup.__faceArray[faceoffset+3] = face[1], geoGroup.__faceArray[faceoffset+4] = face[2], geoGroup.__faceArray[faceoffset+5] = face[3];
                    
                    geoGroup.faceidx += 6;
                    
                    if (currentAtom.clickable || lastAtom.clickable) {
                        
                        var p1a = vs[face[3]].clone(), p1b = vs[face[0]].clone(),
                            p2a = vs[face[2]].clone(), p2b = vs[face[1]].clone();
                        
                        p1a.atom = vs[face[3]].atom || null; //should be same
                        p2a.atom = vs[face[2]].atom || null; 
                        
                        
                        p1b.atom = vs[face[0]].atom || null; //should be same                      
                        p2b.atom = vs[face[1]].atom || null; 
                            
                        var face1, face2, face3;
                        
                        if (diffAtoms) {
                            var m1 = p1a.clone().add(p1b).multiplyScalar(0.5);
                            var m2 = p2a.clone().add(p2b).multiplyScalar(0.5);
                            var m = p1a.clone().add(p2b).multiplyScalar(0.5);
                            
                            if (j % 2 === 0)
                            {
                                if (lastAtom.clickable) {
                                    face1 = new WebMol.Triangle(m1, m, p1a);
                                    face2 = new WebMol.Triangle(m2, p2a, m);
                                    face3 = new WebMol.Triangle(m, p2a, p1a);
                                    lastAtom.intersectionShape.triangle.push(face1);
                                    lastAtom.intersectionShape.triangle.push(face2);
                                    lastAtom.intersectionShape.triangle.push(face3);
                                }
                                
                                if (currentAtom.clickable) {
                                    face1 = new WebMol.Triangle(p1b, p2b, m);
                                    face2 = new WebMol.Triangle(p2b, m2, m);
                                    face3 = new WebMol.Triangle(p1b, m, m1);
                                    currentAtom.intersectionShape.triangle.push(face1);
                                    currentAtom.intersectionShape.triangle.push(face2);
                                    currentAtom.intersectionShape.triangle.push(face3);
                                }
                            }
                            else {
                                if (currentAtom.clickable) {
                                    face1 = new WebMol.Triangle(m1, m, p1a);
                                    face2 = new WebMol.Triangle(m2, p2a, m);
                                    face3 = new WebMol.Triangle(m, p2a, p1a);
                                    currentAtom.intersectionShape.triangle.push(face1);
                                    currentAtom.intersectionShape.triangle.push(face2);
                                    currentAtom.intersectionShape.triangle.push(face3);
                                }
                                
                                if (lastAtom.clickable) {
                                    face1 = new WebMol.Triangle(p1b, p2b, m);
                                    face2 = new WebMol.Triangle(p2b, m2, m);
                                    face3 = new WebMol.Triangle(p1b, m, m1);
                                    lastAtom.intersectionShape.triangle.push(face1);
                                    lastAtom.intersectionShape.triangle.push(face2);
                                    lastAtom.intersectionShape.triangle.push(face3);
                                }                          
                            }
                            
                        }
                        
                        //face for single atom
                        else if (currentAtom.clickable) {
                            face1 = new WebMol.Triangle(p1b, p2b, p1a);
                            face2 = new WebMol.Triangle(p2b, p2a, p1a);
                            currentAtom.intersectionShape.triangle.push(face1);
                            currentAtom.intersectionShape.triangle.push(face2);
                        }
                        
                    }
                    
                }
            }
            
            geoGroup.vertices += 8;
            lastAtom = currentAtom;
        }
        

        var vsize = vs.length - 8; // Cap
        
        geoGroup = geo.updateGeoGroup(8);
        offset = geoGroup.vertices, vertoffset = offset*3, faceoffset = geoGroup.faceidx;
        
        for ( var i = 0; i < 4; i++) {
            vs.push(vs[i * 2]);
            vs.push(vs[vsize + i * 2]);
            
            var v1 = vs[i * 2], v2 = vs[vsize + i * 2];
            
            geoGroup.__vertexArray[vertoffset+6*i] = v1.x, geoGroup.__vertexArray[vertoffset+1+6*i] = v1.y, geoGroup.__vertexArray[vertoffset+2+6*i] = v1.z;
            geoGroup.__vertexArray[vertoffset+3+6*i] = v2.x, geoGroup.__vertexArray[vertoffset+4+6*i] = v2.y, geoGroup.__vertexArray[vertoffset+5+6*i] = v2.z;
            
            geoGroup.__colorArray[vertoffset+6*i] = color.r, geoGroup.__colorArray[vertoffset+1+6*i] = color.g, geoGroup.__colorArray[vertoffset+2+6*i] = color.b;
            geoGroup.__colorArray[vertoffset+3+6*i] = color.r, geoGroup.__colorArray[vertoffset+4+6*i] = color.g, geoGroup.__colorArray[vertoffset+5+6*i] = color.b;

        }
        
        vsize += 8;
                
        var face1 = [offset, offset + 2, offset + 6, offset + 4];
        var face2 = [offset + 1, offset + 5, offset + 7, offset + 3];
        
        geoGroup.__faceArray[faceoffset] = face1[0], geoGroup.__faceArray[faceoffset+1] = face1[1], geoGroup.__faceArray[faceoffset+2] = face1[3];
        geoGroup.__faceArray[faceoffset+3] = face1[1], geoGroup.__faceArray[faceoffset+4] = face1[2], geoGroup.__faceArray[faceoffset+5] = face1[3];
        geoGroup.__faceArray[faceoffset+6] = face2[0], geoGroup.__faceArray[faceoffset+7] = face2[1], geoGroup.__faceArray[faceoffset+8] = face2[3];
        geoGroup.__faceArray[faceoffset+9] = face2[1], geoGroup.__faceArray[faceoffset+10] = face2[2], geoGroup.__faceArray[faceoffset+11] = face2[3];
        
        geoGroup.faceidx += 12;
        geoGroup.vertices += 8;
        
        //TODO: Add intersection planes for caps
        
        geo.initTypedArrays();
        geo.setUpNormals();
        
        var material = new WebMol.MeshLambertMaterial();
        material.vertexColors = WebMol.FaceColors;
        material.side = WebMol.DoubleSide;
        var mesh = new WebMol.Mesh(geo, material);
        group.add(mesh);
        
    };

    //TODO: Need to update this (will we ever use this?)
    var drawSmoothCurve = function(group, _points, width, colors, div) {
        if (_points.length == 0)
            return;

        div = (div == undefined) ? 5 : div;

        var geo = new WebMol.Geometry();
        var points = subdivide(_points, div);
                /*
        for ( var i = 0; i < points.length; i++) {
            geo.vertices.push(points[i]);
            geo.colors.push(WebMol.color(colors[(i == 0) ? 0 : Math.round((i - 1)
                    / div)]));
        }
                */
        var lineMaterial = new WebMol.LineBasicMaterial({
            linewidth : width
        });
        lineMaterial.vertexColors = true;
        var line = new WebMol.Line(geo, lineMaterial);
        line.type = WebMol.LineStrip;
        group.add(line);
    };

    var drawStrand = function(group, atomlist, num, div, fill, coilWidth,
            helixSheetWidth, doNotSmoothen, thickness) {
        num = num || strandDIV;
        div = div || axisDIV;
        doNotSmoothen == (doNotSmoothen == undefined) ? false : doNotSmoothen;
        var points = [];
        for ( var k = 0; k < num; k++)
            points[k] = [];
        var colors = [];
        var currentChain, currentReschain, currentResi, currentCA;
        var prevCO = null, ss = null, ssborder = false;

        for ( var i in atomlist) {
            var atom = atomlist[i];
            if (atom == undefined)
                continue;

            if ((atom.atom == 'O' || atom.atom == 'CA') && !atom.hetflag) {
                if (atom.atom == 'CA') {
                    if (currentChain != atom.chain 
                            || currentResi + 1 != atom.resi || currentReschain != atom.reschain) {
                        for ( var j = 0; !thickness && j < num; j++)
                            drawSmoothCurve(group, points[j], 1, colors, div);
                        if (fill)
                            drawStrip(group, points[0], points[num - 1],
                                    colors, div, thickness);
                        var points = [];
                        for ( var k = 0; k < num; k++)
                            points[k] = [];
                        colors = [];
                        prevCO = null;
                        ss = null;
                        ssborder = false;
                    }
                    currentCA = new WebMol.Vector3(atom.x, atom.y, atom.z);
                    currentAtom = atom;
                    currentChain = atom.chain;
                    currentReschain = atom.reschain;
                    currentResi = atom.resi;
                    ss = atom.ss;
                    ssborder = atom.ssbegin || atom.ssend;
                    var atomcolor = atom.color;
                    if(typeof(atom.style.cartoon.color) != "undefined") {
                        atomcolor = atom.style.cartoon.color;
                    }
                    colors.push(atomcolor);
                    
                    if (atom.clickable === true && (atom.intersectionShape === undefined || atom.intersectionShape.triangle === undefined)) 
                        atom.intersectionShape = {sphere : null, cylinder : [], line : [], triangle : []};
                    
                } else { // O
                    var O = new WebMol.Vector3(atom.x, atom.y, atom.z);
                    O.sub(currentCA);
                    O.normalize(); // can be omitted for performance
                    O.multiplyScalar((ss == 'c') ? coilWidth : helixSheetWidth);
                    if (prevCO != undefined && O.dot(prevCO) < 0)
                        O.negate();
                    prevCO = O;
                    for ( var j = 0; j < num; j++) {
                        var delta = -1 + 2 / (num - 1) * j;
                        var v = new WebMol.Vector3(currentCA.x + prevCO.x * delta,
                                currentCA.y + prevCO.y * delta, currentCA.z
                                        + prevCO.z * delta);
                        v.atom = currentAtom;
                        if (!doNotSmoothen && ss == 's')
                            v.smoothen = true;
                        points[j].push(v);
                    }
                }
            }
        }
        for ( var j = 0; !thickness && j < num; j++)
            drawSmoothCurve(group, points[j], 1, colors, div);
        if (fill)
            drawStrip(group, points[0], points[num - 1], colors, div, thickness);
    };

    // actual function call
    var drawCartoon = function(group, atomlist) {
        
        drawStrand(group, atomlist, 2, undefined, true, coilWidth, helixSheetWidth,
                false, thickness);
    };

    return drawCartoon;
})();
// A model is a collection of related atoms.  Bonds are only allowed between
//atoms in the same model.  An atom is uniquely specified by its model id and
//its serial number.
//A glmodel knows how to apply the styles on each atom to create a gl object

var WebMol = WebMol || {};

WebMol.GLModel = (function() {

    // class variables go here
    var defaultAtomStyle = {
        line : {},
    };

    var Nucleotides = [ '  G', '  A', '  T', '  C', '  U', ' DG', ' DA', ' DT',
            ' DC', ' DU' ];

    var defaultlineWidth = 1.0;

    // Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
    var vdwRadii = {
        "H" : 1.2,
        "Li" : 1.82,
        "Na" : 2.27,
        "K" : 2.75,
        "C" : 1.7,
        "N" : 1.55,
        "O" : 1.52,
        "F" : 1.47,
        "P" : 1.80,
        "S" : 1.80,
        "CL" : 1.75,
        "BR" : 1.85,
        "SE" : 1.90,
        "ZN" : 1.39,
        "CU" : 1.4,
        "NI" : 1.63
    };

    // class functions

    // return true if a and b represent the same style
    var sameObj = function(a,b) {
        if(a && b)
            return JSON.stringify(a) == JSON.stringify(b);
        else
            return a == b;
    };
    

    // return true if atom1 and atom2 are probably bonded to each other
    // based on distance alone
    var areConnected = function(atom1, atom2) {
        var maxsq = 3.6;

        var xdiff = atom1.x - atom2.x;
        xdiff *= xdiff;
        if (xdiff > maxsq)
            return false;
        var ydiff = atom1.y - atom2.y;
        ydiff *= ydiff;
        if (ydiff > maxsq)
            return false;
        var zdiff = atom1.z - atom2.z;
        zdiff *= zdiff;
        if (zdiff > maxsq)
            return false;

        var distSquared = xdiff + ydiff + zdiff;

        if (isNaN(distSquared))
            return false;
        if (distSquared < 0.5)
            return false; // maybe duplicate position.

        if (distSquared > 1.3
                && (atom1.elem == 'H' || atom2.elem == 'H' || atom1.elem == 'D' || atom2.elem == 'D'))
            return false;
        if (distSquared < 3.6 && (atom1.elem == 'S' || atom2.elem == 'S'))
            return true;
        if (distSquared > 2.78)
            return false;
        return true;
    };

    var assignBonds = function(atomsarray) {
        // assign bonds - yuck, can't count on connect records
        var atoms = atomsarray.slice(0);
        for ( var i = 0; i < atomsarray.length; i++)
        {
            //Don't reindex if atoms are already indexed 
            if (!atomsarray[i].index)
                atomsarray[i].index = i;
        }
        
        atoms.sort(function(a, b) {
            return a.z - b.z;
        });
        for ( var i = 0; i < atoms.length; i++) {
            var ai = atoms[i];

            for ( var j = i + 1; j < atoms.length; j++) {
                var aj = atoms[j];
                if (aj.z - ai.z > 1.9) // can't be connected
                    break;
                if (areConnected(ai, aj)) {
                    if (ai.bonds.indexOf(aj.index) == -1) {
                        // only add if not already there
                        ai.bonds.push(aj.index);
                        ai.bondOrder.push(1);
                        aj.bonds.push(ai.index);
                        aj.bondOrder.push(1);
                    }
                }
            }
        }
    };
    
    // this is optimized for proteins where it is assumed connected
    // atoms are on the same or next residue
    var assignPDBBonds = function(atomsarray) {
        // assign bonds - yuck, can't count on connect records
        var protatoms = [];
        var hetatoms = [];        
        
        for ( var i = 0; i < atomsarray.length; i++) {
            var atom = atomsarray[i];
            atom.index = i;
            if(atom.hetflag)
                hetatoms.push(atom);
            else
                protatoms.push(atom);
        }

        assignBonds(hetatoms);
        
        // sort by resid
        protatoms.sort(function(a, b) {
            if(a.chain != b.chain)
                return a.chain < b.chain ? -1 : 1;
            return a.resi - b.resi;
        });
        
        //for identifying connected residues
        var currentResi = -1;
        var reschain = -1;
        var lastResConnected;
        
        for ( var i = 0; i < protatoms.length; i++) {
            var ai = protatoms[i];
            
            if (ai.resi !== currentResi) {
                currentResi = ai.resi;
                if (!lastResConnected)
                    reschain++;
                    
                lastResConnected = false;
            }
            
            ai.reschain = reschain;

            for ( var j = i + 1; j < protatoms.length; j++) {
                var aj = protatoms[j];
                if(aj.chain != ai.chain)
                    break;
                if (aj.resi - ai.resi > 1) // can't be connected
                    break;
                if (areConnected(ai, aj)) {
                    if (ai.bonds.indexOf(aj.index) === -1) {
                        // only add if not already there
                        ai.bonds.push(aj.index);
                        ai.bondOrder.push(1);
                        aj.bonds.push(ai.index);
                        aj.bondOrder.push(1);
                    }
                    
                    if (ai.resi !== aj.resi) 
                        lastResConnected = true;                   
                        
                }
            }
        }
        
    };

    // return distance between donor-acceptor, if not valid pair, return inf
    var hbondDistance = function(a1, a2, maxlength) {
        if(a1.chain == a2.chain) { // ignore if residues too close
            if(Math.abs(a1.resi-a2.resi) < 4)
                return Number.POSITIVE_INFINITY;
        }
        if ((a1.atom === "O" && a2.atom === "N")
                || (a1.atom === "N" && a2.atom === "O")) {
            var xdiff = a1.x - a2.x;
            if (xdiff > maxlength)
                return Number.POSITIVE_INFINITY;
            var ydiff = a1.y - a2.y;
            if (ydiff > maxlength)
                return Number.POSITIVE_INFINITY;
            var zdiff = a1.z - a2.z;
            if (zdiff > maxlength)
                return Number.POSITIVE_INFINITY;
            
            var dist = Math.sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff);
            if(dist < maxlength)
                return dist;
        }
        return Number.POSITIVE_INFINITY;
    };

    // this will identify all hydrogen bonds between backbone
    // atoms; assume atom names are correct, only identifies
    // single closest hbond
    var assignBackboneHBonds = function(atomsarray) {
        var maxlength = 3.5; // ver generous hbond distance
        var atoms = [];
        for ( var i = 0; i < atomsarray.length; i++) {
            atomsarray[i].index = i;
            // only consider 'N' and 'O'
            var atom = atomsarray[i];
            if (!atom.hetflag && (atom.atom === "N" || atom.atom === "O")) {
                atoms.push(atom);
                atom.hbondOther = null;
                atom.hbondDistance = Number.POSITIVE_INFINITY;                
            }
        }

        atoms.sort(function(a, b) {
            return a.z - b.z;
        });
        for ( var i = 0; i < atoms.length; i++) {
            var ai = atoms[i];

            for ( var j = i + 1; j < atoms.length; j++) {
                var aj = atoms[j];
                if (aj.z - ai.z > maxlength) // can't be connected
                    break;
                var dist = hbondDistance(ai,aj,maxlength);
                if (dist < ai.hbondDistance) {
                    ai.hbondOther = aj;
                    ai.hbondDistance = dist;
                }
                if(dist < aj.hbondDistance) {
                    aj.hbondOther = ai;
                    aj.hbondDistance = dist;
                }
            }
        }
    };

    var computeSecondaryStructure = function(atomsarray) {
        assignBackboneHBonds(atomsarray);
        
        // compute, per residue, what the secondary structure is
        var chres = {}; // lookup by chain and resid
        for ( var i = 0; i < atomsarray.length; i++) {
            var atom = atomsarray[i];
            
            if(typeof(chres[atom.chain]) === "undefined")
                chres[atom.chain] = [];
            
            if(isFinite(atom.hbondDistance)) {
                var other = atom.hbondOther;
                if(Math.abs(other.resi - atom.resi) === 4) { 
                    // helix
                    chres[atom.chain][atom.resi] = 'h';
                }
                else { // otherwise assume sheet
                    chres[atom.chain][atom.resi] = 's';
                }
            }
        }
        
        // plug gaps and remove singletons
        for(c in chres) {
            for(var r = 1; r < chres[c].length-1; r++) {
                var valbefore = chres[c][r-1];
                var valafter = chres[c][r+1];
                var val = chres[c][r];
                if(valbefore == valafter && val != valbefore) {
                    chres[c][r] = valbefore;
                }
            }
            for(var r = 0; r < chres[c].length; r++) {
                var val = chres[c][r];
                if(val == 'h' || val == 's') {
                    if(chres[c][r-1] != val && chres[c][r+1] != val)
                        delete chres[c][r];
                }
            }
        }
        
        // assign to all atoms in residue, keep track of start
        var curres = null;
        for ( var i = 0; i < atomsarray.length; i++) {
            var atom = atomsarray[i];
            var val = chres[atom.chain][atom.resi];
            if(typeof(val) == "undefined")
                continue;
            atom.ss = val;
            if(chres[atom.chain][atom.resi-1] != val)
                atom.ssbegin = true;
            if(chres[atom.chain][atom.resi+1] != val)
                atom.ssend = true;
        }
    };
    
    var parseCube = function(atoms, str) {
        var lines = str.replace(/^\s+/, "").split(/[\n\r]+/);
        
        if (lines.length < 6)
            return;
            
        var lineArr = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");       
          
        var natoms = Math.abs(parseFloat(lineArr[0]));        
        
        lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        
        //might have to convert from bohr units to angstroms
        var convFactor = (parseFloat(lineArr[0]) > 0) ? 0.529177 : 1;
        
        //Extract atom portion; send to new GLModel...
        lines = lines.splice(6, natoms);
       
        var start = atoms.length;
        var end = start + lines.length;
        
        for (var i = start; i < end; ++i) {
            var atom = {};
            atom.serial = i;
            var line = lines[i - start];
            var tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
                    " ");
            
            if (tokens[0] == 6) 
                atom.elem = "C";
                        
            else if (tokens[0] == 1) 
                atom.elem = "H";
            
            else if (tokens[0] == 8)
                atom.elem = "O";
                
            else if (tokens[0] == 17)
                atom.elem = "CL";
                
            atom.x = parseFloat(tokens[2]) * convFactor;
            atom.y = parseFloat(tokens[3]) * convFactor;
            atom.z = parseFloat(tokens[4]) * convFactor;
            
            atom.hetflag = true;
            atom.singleBonds = true;
            atom.bonds = [];
            atom.bondOrder = [];
            atom.properties = {};
            atoms.push(atom);
            
        }   
        
        assignBonds(atoms);
        
        return true; 
    };
        
    // read an XYZ file from str and put the result in atoms
    var parseXYZ = function(atoms, str) {

        var lines = str.split("\n");
        if (lines.length < 3)
            return;
        var atomCount = parseInt(lines[0].substr(0, 3));
        if (isNaN(atomCount) || atomCount <= 0)
            return;
        if (lines.length < atomCount + 2)
            return;
        var offset = 2;
        var start = atoms.length;
        var end = start + atomCount;
        for ( var i = start; i < end; i++) {
            var line = lines[offset++];
            var tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
                    " ");
            var atom = {};
            atom.serial = i;
            atom.atom = atom.elem = tokens[0];
            atom.x = parseFloat(tokens[1]);
            atom.y = parseFloat(tokens[2]);
            atom.z = parseFloat(tokens[3]);
            atom.hetflag = true;
            atom.bonds = [];
            atom.bondOrder = [];
            atom.singleBonds = true;
            atom.properties = {};
            atoms[i] = atom;
        }
        assignBonds(atoms);

        return true;
    };

    // put atoms specified in sdf fromat in str into atoms
    // adds to atoms, does not replace
    var parseSDF = function(atoms, str) {

        var lines = str.split("\n");
        if (lines.length < 4)
            return;
        var atomCount = parseInt(lines[3].substr(0, 3));
        if (isNaN(atomCount) || atomCount <= 0)
            return;
        var bondCount = parseInt(lines[3].substr(3, 3));
        var offset = 4;
        if (lines.length < 4 + atomCount + bondCount)
            return;
        var start = atoms.length;
        var end = start + atomCount;
        for ( var i = start; i < end; i++) {
            var line = lines[offset];
            offset++;
            var atom = {};
            atom.serial = i;
            atom.x = parseFloat(line.substr(0, 10));
            atom.y = parseFloat(line.substr(10, 10));
            atom.z = parseFloat(line.substr(20, 10));
            atom.hetflag = true;
            atom.singleBonds = true; //atom only makes single bonds ?
            atom.atom = atom.elem = line.substr(31, 3).replace(/ /g, "");
            atom.bonds = [];
            atom.bondOrder = [];
            atom.properties = {};
            atoms[i] = atom;
        }
        for (i = 0; i < bondCount; i++) {
            var line = lines[offset];
            offset++;
            var from = parseInt(line.substr(0, 3)) - 1 + start;
            var to = parseInt(line.substr(3, 3)) - 1 + start;
            var order = parseInt(line.substr(6, 3));
            if (order > 1)
                atoms[from].singleBonds = false, atoms[to].singleBonds = false;
            atoms[from].bonds.push(to);
            atoms[from].bondOrder.push(order);
            atoms[to].bonds.push(from);
            atoms[to].bondOrder.push(order);
        }

        return true;
    };

    // parse SYBYL mol2 file from string - assumed to only contain one molecule
    // tag
    // TODO: Figure out how to handle multi molecule files (for SDF, too)
    var parseMOL2 = function(atoms, str, keepH) {
        
        var noH = !keepH; // again, suppress H's by default
        
        // Note: these regex's work, though they don't match '<TRIPOS>'
        // correctly - something to do with angle brackets
        var mol_pos = str.search(/@<TRIPOS>MOLECULE/);
        var atom_pos = str.search(/@<TRIPOS>ATOM/);
        
        // Assuming both Molecule and Atom sections exist
        if (mol_pos == -1 || atom_pos == -1)
            return;
        
        // serial is atom's index in file; index is atoms index in 'atoms'
        var serialToIndex = [];
        

        // assert (mol_pos < atom_pos), "Unexpected formatting of mol2 file
        // (expected 'molecule' section before 'atom' section)";
        

        var lines = str.substr(mol_pos, str.length).split("\n");
        var tokens = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        var natoms = parseInt(tokens[0]);
        if (tokens.length > 1)
            var nbonds = parseInt(tokens[1]); // NOTE: this may or may not be
                                                // set
        else
            var nbonds = 0;
        
        var offset = 4;
        // Continue until 'Atom' section
        for (var i = 3; ;i++)
        {
            if (lines[i] == "@<TRIPOS>ATOM")
            {
                offset = i+1;
                break;
            }
        }
        
        var start = atoms.length;
        var end = start + natoms;
        
        // Process ATOMS
        for (var i = start; i < end; i++)
        {
            var line = lines[offset++];
            var tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
                    " ");
            var atom = {};
            
            // 'index' is this atom's index in 'atoms'; 'serial' is this atom's
            // serial id in mol2 file
            var index = i;
            var serial = parseInt(tokens[0]);
            atom.serial = serial;
            // atom.serial = i;
            
            atom.x = parseFloat(tokens[2]);
            atom.y = parseFloat(tokens[3]);
            atom.z = parseFloat(tokens[4]);
            atom.atom = atom.elem = tokens[5].split('.')[0];
            
            atom.singleBonds = true;
            
            // TODO: Add capability to ignore H's

            if (atom.elem == 'H' && noH)
                continue;
                
            atom.bonds = [];
            atom.bondOrder = [];
            atom.properties = {};
            
            serialToIndex[serial] = index;
            atoms.push(atom);
        }
        
        // Process BONDS
        var bonds_found = false;
        while (offset < lines.length)
        {
            if (lines[offset++] == "@<TRIPOS>BOND")
            {
                bonds_found = true;
                break;            
            }        
        }
        
        if (bonds_found && nbonds)
        {
            for (var i = 0; i < nbonds; i++)
            {
                var line = lines[offset++];

                var tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
                            " ");
                var from = parseInt(tokens[1]);
                fromAtom = atoms[serialToIndex[from]];
                var to = parseInt(tokens[2]);
                toAtom = atoms[serialToIndex[to]];              
                    
                // Won't be able to read aromatic bonds correctly...
                var order = parseInt(tokens[3]);
                if (isNaN(order))
                    order = 1;
                
                if (order > 1)
                    fromAtom.singleBonds = false, toAtom.singleBonds = false;
                
                if (fromAtom != undefined && toAtom != undefined){
                    fromAtom.bonds.push(serialToIndex[to]);
                    fromAtom.bondOrder.push(order);
                    toAtom.bonds.push(serialToIndex[from]);
                    toAtom.bondOrder.push(order);
                }    

                
                /*
                 * atoms[from].bonds.push(to);
                 * atoms[from].bondOrder.push(order);
                 * atoms[to].bonds.push(from); atoms[to].bondOrder.push(order);
                 */

            }
        }
        
        return true;
        
    };
    
    // parse pdb file from str and create atoms
    //if computeStruct is true will always perform secondary structure analysis,
    //otherwise only do analysis of SHEET/HELIX comments are missing
    var parsePDB = function(atoms, str, keepH, computeStruct) {

        var atoms_cnt = 0;
        var noH = !keepH; // suppress hydrogens by default
        var start = atoms.length;
        var protein = {
            sheet : [],
            helix : []
        }; // get secondary structure straight from pdb

        var hasStruct = false;
        var serialToIndex = []; // map from pdb serial to index in atoms
        lines = str.split("\n");
        for ( var i = 0; i < lines.length; i++) {
            line = lines[i].replace(/^\s*/, ''); // remove indent
            var recordName = line.substr(0, 6);
            if (recordName == 'ATOM  ' || recordName == 'HETATM') {
                var atom, resn, chain, resi, icode, x, y, z, hetflag, elem, serial, altLoc, b;
                altLoc = line.substr(16, 1);
                if (altLoc != ' ' && altLoc != 'A')
                    continue; // FIXME: ad hoc
                serial = parseInt(line.substr(6, 5));
                atom = line.substr(12, 4).replace(/ /g, "");
                resn = line.substr(17, 3);
                chain = line.substr(21, 1);
                resi = parseInt(line.substr(22, 4));
                icode = line.substr(26, 1);
                x = parseFloat(line.substr(30, 8));
                y = parseFloat(line.substr(38, 8));
                z = parseFloat(line.substr(46, 8));
                b = parseFloat(line.substr(60, 8));
                elem = line.substr(76, 2).replace(/ /g, "");
                if (elem == '') { // for some incorrect PDB files
                    elem = line.substr(12, 2).replace(/ /g, "");
                }
                if((elem == 'H' || elem == 'HH' || elem == 'HD') && noH)
                    continue;
                if (line[0] == 'H')
                    hetflag = true;
                else
                    hetflag = false;
                serialToIndex[serial] = atoms.length;
                atoms.push({
                    'resn' : resn,
                    'x' : x,
                    'y' : y,
                    'z' : z,
                    'elem' : elem,
                    'hetflag' : hetflag,
                    'chain' : chain,
                    'resi' : resi,
                    'icode' : icode,
                    'rescode': resi + (icode != ' ' ? "^"+icode: ""), // combo
                                                                        // resi
                                                                        // and
                                                                        // icode
                    'serial' : serial,
                    'atom' : atom,
                    'bonds' : [],
                    'ss' : 'c',
                    'singleBonds' : true,
                    'bonds' : [],
                    'bondOrder' : [],
                    'properties' : {},
                    'b' : b,
                    'pdbline' : line
                });
            } else if (recordName == 'SHEET ') {
            	hasStruct = true;
                var startChain = line.substr(21, 1);
                var startResi = parseInt(line.substr(22, 4));
                var endChain = line.substr(32, 1);
                var endResi = parseInt(line.substr(33, 4));
                protein.sheet
                        .push([ startChain, startResi, endChain, endResi ]);
            } else if (recordName == 'CONECT') {
                // MEMO: We don't have to parse SSBOND, LINK because both are
                // also
                // described in CONECT. But what about 2JYT???
                var from = parseInt(line.substr(6, 5));
                var fromAtom = atoms[serialToIndex[from]];
                for ( var j = 0; j < 4; j++) {
                    var to = parseInt(line.substr([ 11, 16, 21, 26 ][j], 5));
                    var toAtom = atoms[serialToIndex[to]];
                    if (fromAtom != undefined && toAtom != undefined) {
                        fromAtom.bonds.push(serialToIndex[to]);
                        fromAtom.bondOrder.push(1);
                    }
                }
            } else if (recordName == 'HELIX ') {
            	hasStruct = true;
                var startChain = line.substr(19, 1);
                var startResi = parseInt(line.substr(21, 4));
                var endChain = line.substr(31, 1);
                var endResi = parseInt(line.substr(33, 4));
                protein.helix
                        .push([ startChain, startResi, endChain, endResi ]);
            }

        }

        var starttime = (new Date()).getTime();
        // assign bonds - yuck, can't count on connect records
        assignPDBBonds(atoms);
        console.log("bond connecting " + ((new Date()).getTime() - starttime));
        
        
        if(computeStruct || !hasStruct) {
            var starttime = (new Date()).getTime();
        	computeSecondaryStructure(atoms);
        	console.log("secondary structure " + ((new Date()).getTime() - starttime));
        }
        
        // Assign secondary structures from pdb file
        for (i = start; i < atoms.length; i++) {
            atom = atoms[i];
            if (atom == undefined)
                continue;

            var found = false;
            // MEMO: Can start chain and end chain differ?
            for (j = 0; j < protein.sheet.length; j++) {
                if (atom.chain != protein.sheet[j][0])
                    continue;
                if (atom.resi < protein.sheet[j][1])
                    continue;
                if (atom.resi > protein.sheet[j][3])
                    continue;
                atom.ss = 's';
                if (atom.resi == protein.sheet[j][1])
                    atom.ssbegin = true;
                if (atom.resi == protein.sheet[j][3])
                    atom.ssend = true;
            }
            for (j = 0; j < protein.helix.length; j++) {
                if (atom.chain != protein.helix[j][0])
                    continue;
                if (atom.resi < protein.helix[j][1])
                    continue;
                if (atom.resi > protein.helix[j][3])
                    continue;
                atom.ss = 'h';
                if (atom.resi == protein.helix[j][1])
                    atom.ssbegin = true;
                else if (atom.resi == protein.helix[j][3])
                    atom.ssend = true;
            }
        }
        return true;
    };

    /**
     * GLModel represents a group of related atoms
     * @constructor WebMol.GLModel
     * @param {number} mid 
     * @param {Object} defaultcolors Object defining default atom colors as atom => color property value pairs
     * @see WebMol.download
     */
    function GLModel(mid, defaultcolors) {
        // private variables
        var atoms = [];
        var id = mid;
        var molObj = null;
        var renderedMolObj = null;
        var lastStyle = null; // cache previous styles to avoid recomputation
        var lastColors = null;
        
        var defaultColor = WebMol.defaultElementColor;

        if (defaultcolors)
            ElementColors = defaultcolors;
        else
            ElementColors = WebMol.defaultElementColors;

        // drawing functions must be associated with model object since
        // geometries can't span multiple canvases

        // sphere drawing
        var defaultSphereRadius = 1.5;

        // return proper radius for atom given style
        var getRadiusFromStyle = function(atom, style) {
            var r = defaultSphereRadius;
            if (typeof (style.radius) != "undefined")
                r = style.radius;
            else if (vdwRadii[atom.elem])
                r = vdwRadii[atom.elem];

            if (typeof (style.scale) != "undefined")
                r *= style.scale;
            return r;
        };

        // memoize capped cylinder for given radius
        var cylVertexCache = {
            
            //Ortho normal vectors for cylinder radius/ sphere cap equator
            // Direction is j basis (0,1,0)
            basisVectors : function() {
                
                var ret = {vertices : [], norms : []};
                
                var nvecs = [];
                
                nvecs[0] = new WebMol.Vector3(-1,0,0);
                nvecs[4] = new WebMol.Vector3(0,0,1);
                nvecs[8] = new WebMol.Vector3(1,0,0);
                nvecs[12] = new WebMol.Vector3(0,0,-1);
    
                // now quarter positions
                nvecs[2] = nvecs[0].clone().add(nvecs[4]).normalize();
                nvecs[6] = nvecs[4].clone().add(nvecs[8]).normalize();
                nvecs[10] = nvecs[8].clone().add(nvecs[12]).normalize();
                nvecs[14] = nvecs[12].clone().add(nvecs[0]).normalize();
    
                // eights
                nvecs[1] = nvecs[0].clone().add(nvecs[2]).normalize();
                nvecs[3] = nvecs[2].clone().add(nvecs[4]).normalize();
                nvecs[5] = nvecs[4].clone().add(nvecs[6]).normalize();
                nvecs[7] = nvecs[6].clone().add(nvecs[8]).normalize();
                nvecs[9] = nvecs[8].clone().add(nvecs[10]).normalize();
                nvecs[11] = nvecs[10].clone().add(nvecs[12]).normalize();
                nvecs[13] = nvecs[12].clone().add(nvecs[14]).normalize();
                nvecs[15] = nvecs[14].clone().add(nvecs[0]).normalize(); 
                
                /*
                nvecs[0] = new WebMol.Vector3(-1,0,0);
                nvecs[1] = new WebMol.Vector3(0,0,1);
                nvecs[2] = new WebMol.Vector3(1,0,0);
                nvecs[3] = new WebMol.Vector3(0,0,-1);
                */
                return nvecs;
                                        
            }(),
            
            cache : {},
            
            getVerticesForRadius : function(radius) {
                
                if (this.cache[radius] !== undefined)
                    return this.cache[radius];
                
                var dir = new WebMol.Vector3(0,1,0);    
                var w = this.basisVectors.length;
                var nvecs = [], norms = [];
                
                for (var i = 0; i < w; i++) {
                    //bottom
                    nvecs.push(this.basisVectors[i].clone().multiplyScalar(radius));
                    //top
                    nvecs.push(this.basisVectors[i].clone().multiplyScalar(radius));
                    
                    //NOTE: this normal is used for constructing sphere caps - 
                    // cylinder normals taken care of in drawCylinder
                    var n = this.basisVectors[i].clone().normalize();
                    norms.push(n);
                    norms.push(n);
                }

                //norms[0]   
                
                var verticesRows = [];
                
                //Require that heightSegments is even and >= 2
                //Equator points at h/2 (theta = pi/2)
                //(repeated) polar points at 0 and h (theta = 0 and pi)
                var heightSegments = 10, widthSegments = w; // 16 or however many basis vectors for cylinder
                
                if (heightSegments % 2 !== 0 || !heightSegments) {
                    console.error("heightSegments must be even");
                    
                    return null;
                }        
                
                var phiStart = 0;
                var phiLength = Math.PI * 2;

                var thetaStart = 0;
                var thetaLength = Math.PI;

                var x, y;
                var polar = false, equator = false;
                
                for (y = 0; y <= heightSegments; y++) {        
                    
                    polar = (y === 0 || y === heightSegments) ? true : false;
                    equator = (y === heightSegments/2) ? true : false;                 
                    
                    var verticesRow = [], toRow = [];
                    
                    for (x = 0; x <= widthSegments; x++) {
                        
                        // Two vertices rows for equator pointing to previously constructed cyl points
                        if (equator) {
                            var xi = (x < widthSegments) ? 2*x : 0;
                            toRow.push(xi+1), verticesRow.push(xi);
                            
                            continue;
                        }
                        
                        var u = x / widthSegments;
                        var v = y / heightSegments;
                        
                        //Only push first polar point
                        
                        if (!polar || x === 0) {
                            
                            if (x < widthSegments) {
                                var vertex = new WebMol.Vector3();
                                vertex.x = -radius * Math.cos(phiStart + u * phiLength)
                                        * Math.sin(thetaStart + v * thetaLength);
                                vertex.y = radius
                                        * Math.cos(thetaStart + v * thetaLength);
                                vertex.z = radius * Math.sin(phiStart + u * phiLength)
                                        * Math.sin(thetaStart + v * thetaLength);

                                if (Math.abs(vertex.x) < 1e-5) vertex.x = 0;
                                if (Math.abs(vertex.y) < 1e-5) vertex.y = 0;
                                if (Math.abs(vertex.z) < 1e-5) vertex.z = 0;

                                var n = new WebMol.Vector3(vertex.x, vertex.y, vertex.z);
                                n.normalize();

                                nvecs.push(vertex);
                                norms.push(n);   
                                
                                verticesRow.push(nvecs.length - 1);
                            }
                            
                            //last point is just the first point for this row
                            else {
                                verticesRow.push(nvecs.length - widthSegments);
                            }
                                                       
                        }
                        
                        // x > 0; index to already added point
                        else if (polar) 
                            verticesRow.push(nvecs.length - 1);
                        
                    }
                    
                    //extra equator row
                    if (equator)
                        verticesRows.push(toRow);
                    
                    verticesRows.push(verticesRow);
                    
                }         
                
                var obj = {
                    vertices : nvecs,
                    normals : norms,
                    verticesRows : verticesRows,
                    w : widthSegments,
                    h : heightSegments
                };  
                
                this.cache[radius] = obj;
                
                return obj;
                
            }
        };

        // construct vertices around origin for given radius, memoize results
        var sphereVertexCache = {
            cache : {},
            getVerticesForRadius : function(radius) {

                if (typeof (this.cache[radius]) != "undefined")
                    return this.cache[radius];

                var obj = {
                    vertices : [],
                    verticesRows : [],
                    normals : []
                };
                // scale quality with radius heuristically
                var widthSegments = 12;
                var heightSegments = 10;
                if (radius < 1) {
                    widthSegments = 8;
                    heightSegments = 6;
                }

                var phiStart = 0;
                var phiLength = Math.PI * 2;

                var thetaStart = 0;
                var thetaLength = Math.PI;

                var x, y, vertices = [], uvs = [];

                for (y = 0; y <= heightSegments; y++) {

                    var verticesRow = [];
                    for (x = 0; x <= widthSegments; x++) {

                        var u = x / widthSegments;
                        var v = y / heightSegments;

                        var vertex = {};
                        vertex.x = -radius * Math.cos(phiStart + u * phiLength)
                                * Math.sin(thetaStart + v * thetaLength);
                        vertex.y = radius
                                * Math.cos(thetaStart + v * thetaLength);
                        vertex.z = radius * Math.sin(phiStart + u * phiLength)
                                * Math.sin(thetaStart + v * thetaLength);

                        var n = new WebMol.Vector3(vertex.x, vertex.y, vertex.z);
                        n.normalize();

                        obj.vertices.push(vertex);
                        obj.normals.push(n);

                        verticesRow.push(obj.vertices.length - 1);

                    }

                    obj.verticesRows.push(verticesRow);

                }

                this.cache[radius] = obj;
                return obj;
            }
        };
        
        // cross drawing
        var drawAtomCross = function(atom, geos) {
            if (!atom.style.cross)
                return;
            var style = atom.style.cross;
            if (style.hidden)
                return;
            var linewidth = (style.linewidth || defaultlineWidth);
            if (!geos[linewidth])
                geos[linewidth] = new WebMol.Geometry();
                
            var geoGroup = geos[linewidth].updateGeoGroup(6);
            
            var delta = getRadiusFromStyle(atom, style);

            var points = [ [ delta, 0, 0 ], [ -delta, 0, 0 ], [ 0, delta, 0 ],
                    [ 0, -delta, 0 ], [ 0, 0, delta ], [ 0, 0, -delta ] ];

            var clickable = atom.clickable;
            if (clickable && atom.intersectionShape === undefined)
                atom.intersectionShape = {sphere : [], cylinder : [], line : []};
            
            var c = WebMol.CC.color(atom.color);
            
            for ( var j = 0; j < 6; j++) {
                
                var offset = geoGroup.vertices*3;
                
                geoGroup.vertices++;
                geoGroup.__vertexArray[offset] = atom.x + points[j][0];
                geoGroup.__vertexArray[offset+1] = atom.y + points[j][1];
                geoGroup.__vertexArray[offset+2] = atom.z + points[j][2];
                geoGroup.__colorArray[offset] = c.r;
                geoGroup.__colorArray[offset+1] = c.g;
                geoGroup.__colorArray[offset+2] = c.b;
                
                if (clickable){
                    var point = new WebMol.Vector3(points[j][0], points[j][1], points[j][2]);
                    
                    //decrease cross size for selection to prevent misselection from atom overlap
                    point.multiplyScalar(0.1);
                    point.set(point.x+atom.x, point.y+atom.y, point.z+atom.z);
                    atom.intersectionShape.line.push(point);
                }

            }
                        
        };

        // bonds - both atoms must match bond style
        // standardize on only drawing for lowest to highest
        var drawBondLines = function(atom, atoms, geos) {
            if (!atom.style.line)
                return;
            var style = atom.style.line;
            if (style.hidden)
                return;

            // have a separate geometry for each linewidth
            var linewidth = (style.linewidth || defaultlineWidth);

            if (!geos[linewidth])
                geos[linewidth] = new WebMol.Geometry();
            var geoGroup = geos[linewidth].updateGeoGroup(2*atom.bonds.length);
            

            for ( var i = 0; i < atom.bonds.length; i++) {
                
                var j = atom.bonds[i]; // our neighbor
                // TODO: handle bond orders
                var atom2 = atoms[j];
                if (!atom2.style.line)
                    continue; // don't sweat the details

                var p1 = new WebMol.Vector3(atom.x, atom.y, atom.z);
                var p2 = new WebMol.Vector3(atom2.x, atom2.y, atom2.z);
                var mp = p1.clone().add(p2).multiplyScalar(0.5);

                if (atom.clickable){
                    if (atom.intersectionShape === undefined)
                        atom.intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};
                    atom.intersectionShape.line.push(p1);
                    atom.intersectionShape.line.push(mp);
                }

                var c1 = WebMol.CC.color(atom.color);
                var offset = geoGroup.vertices*3;
                geoGroup.vertices += 2;

                geoGroup.__vertexArray[offset] = p1.x, geoGroup.__vertexArray[offset+1] = p1.y, geoGroup.__vertexArray[offset+2] = p1.z;
                geoGroup.__colorArray[offset] = c1.r, geoGroup.__colorArray[offset+1] = c1.g, geoGroup.__colorArray[offset+2] = c1.b;
                geoGroup.__vertexArray[offset+3] = mp.x, geoGroup.__vertexArray[offset+4] = mp.y, geoGroup.__vertexArray[offset+5] = mp.z;
                geoGroup.__colorArray[offset+3] = c1.r, geoGroup.__colorArray[offset+4] = c1.g, geoGroup.__colorArray[offset+5] = c1.b;

            }

        };

        // bonds as cylinders
        var defaultStickRadius = .25;

        //sphere drawing
        //See also: drawCylinder
        var drawAtomSphere = function(atom, geo) {
            
            if (!atom.style.sphere)
                return;
            var style = atom.style.sphere;
            if (style.hidden)
                return;
                                                                 
            var color = atom.color;
            if (typeof (style.color) != "undefined")
                color = style.color;
            var C = WebMol.CC.color(color);

            var x, y;
            var radius = getRadiusFromStyle(atom, style);
            
            if ((atom.clickable === true) && (atom.intersectionShape !== undefined)) {
                var center = new WebMol.Vector3(atom.x, atom.y, atom.z);
                atom.intersectionShape.sphere.push(new WebMol.Sphere(center, radius));
            }
            
            var vobj = sphereVertexCache.getVerticesForRadius(radius);                
                        
            var vertices = vobj.vertices;
            var normals = vobj.normals;
            
            var geoGroup = geo.updateGeoGroup(vertices.length);
            var start = geoGroup.vertices;
            
            for (var i = 0, il = vertices.length; i < il; ++i) {
                var offset = 3*(start + i);   
                var v = vertices[i];
                
                geoGroup.__vertexArray[offset] = (v.x + atom.x);
                geoGroup.__vertexArray[offset+1] = (v.y + atom.y);
                geoGroup.__vertexArray[offset+2] = (v.z + atom.z);
                
                geoGroup.__colorArray[offset] = C.r;
                geoGroup.__colorArray[offset+1] = C.g;
                geoGroup.__colorArray[offset+2] = C.b;
               
            }
            
            geoGroup.vertices += vertices.length;
            
            var verticesRows = vobj.verticesRows;
            var h = verticesRows.length - 1;
            
            //var color = [C.r, C.g, C.b];
            for (y = 0; y < h; y++) {
                var w = verticesRows[y].length - 1;
                for (x = 0; x < w; x++) {
                    
                    var faceoffset = geoGroup.faceidx;
                    
                    var v1 = verticesRows[y][x + 1] + start, v1offset = v1 * 3;
                    var v2 = verticesRows[y][x] + start, v2offset = v2 * 3;
                    var v3 = verticesRows[y + 1][x] + start, v3offset = v3 * 3;
                    var v4 = verticesRows[y + 1][x + 1] + start, v4offset = v4 * 3;

                    var n1 = normals[v1 - start];
                    var n2 = normals[v2 - start];
                    var n3 = normals[v3 - start];
                    var n4 = normals[v4 - start];
                    var face, norm;
                    if (Math.abs(vertices[v1 - start].y) === radius) {
                        //face = [v1, v3, v4];
                        //norm = [n1, n3, n4];
                        
                        geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v3offset] = n3.x, geoGroup.__normalArray[v4offset] = n4.x;
                        geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v3offset+1] = n3.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                        geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v3offset+2] = n3.z, geoGroup.__normalArray[v4offset+2] = n4.z;

                        geoGroup.__faceArray[faceoffset] = v1; 
                        geoGroup.__faceArray[faceoffset+1] = v3;
                        geoGroup.__faceArray[faceoffset+2] = v4;
                        
                        geoGroup.faceidx += 3;
                        
                    } else if (Math.abs(vertices[v3 - start].y) === radius) {
                        //face = [v1, v2, v3];            
                        //norm = [n1, n2, n3];
                        
                        geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v3offset] = n3.x;
                        geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v3offset+1] = n3.y;
                        geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v3offset+2] = n3.z;

                        geoGroup.__faceArray[faceoffset] = v1;
                        geoGroup.__faceArray[faceoffset+1] = v2;
                        geoGroup.__faceArray[faceoffset+2] = v3;
                        
                        geoGroup.faceidx += 3;
                        
                    } else {
                        //face = [v1, v2, v3, v4];
                        //norm = [n1, n2, n3, n4];
                        
                        geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v4offset] = n4.x;
                        geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                        geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v4offset+2] = n4.z;
                        
                        geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v3offset] = n3.x, geoGroup.__normalArray[v4offset] = n4.x;
                        geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v3offset+1] = n3.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                        geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v3offset+2] = n3.z, geoGroup.__normalArray[v4offset+2] = n4.z;
                        
                        geoGroup.__faceArray[faceoffset] = v1;
                        geoGroup.__faceArray[faceoffset+1] = v2;
                        geoGroup.__faceArray[faceoffset+2] = v4;
                        
                        geoGroup.__faceArray[faceoffset+3] = v2;
                        geoGroup.__faceArray[faceoffset+4] = v3;
                        geoGroup.__faceArray[faceoffset+5] = v4;
                        
                        geoGroup.faceidx += 6;
                    }

                }
            }

        };
        
        // Rotation matrix around z and x axis - 
        // according to y basis vector
        // TODO: Try to optimize this (square roots?)
        var getRotationMatrix = function() {
                      
            var d = new WebMol.Vector3();
            //var rot = new Float32Array(9);
           
            return function(dir) {
               
                d.copy(dir);
                
                var dx = d.x, dy = d.y, dz = d.z;
                
                var dxy = Math.sqrt(dx*dx + dy*dy);
                var dxz, dyz;
                
                var sinA, cosA, sinB, cosB, sinC, cosC;
                
                // about z axis - Phi
                if (dxy < 0.0001)
                    sinA = 0, cosA = 1;
                else 
                    sinA = -dx / dxy, cosA = dy / dxy;
                
                /* 
                Rphi.set( cosA, -sinA,  0, 0,
                          sinA,  cosA,  0, 0,
                          0,     0,     1, 0,
                          0,     0,     0, 1 );   
                */
               
                //recast dy in terms of new axes - z is the same
                
                dy = -sinA*dx + cosA*dy;
                dyz = Math.sqrt(dy*dy + dz*dz);    
                 
                // about new x axis - Theta
                
                if (dyz < 0.0001)
                    sinB = 0, cosB = 1;
                else 
                    sinB =  dz / dyz, cosB = dy / dyz;                                     

                /*        
                Rtheta.set( 1,  0,     0,    0,
                            0,  cosB, -sinB, 0,
                            0,  sinB,  cosB, 0,
                            0,  0,     0,    1 );

                rot.multiplyMatrices(Rphi, Rtheta);
                */
               

                                   
                rot = new Float32Array(9);
                rot[0] = cosA, rot[1] = sinA, rot[2] = 0;
                rot[3] = -sinA*cosB, rot[4] = cosA*cosB, rot[5] = sinB;
                rot[6] = sinA*sinB, rot[7] = -cosA*sinB, rot[8] = cosB;
                return rot;
            
            };
            
        }();
        
        // creates a cylinder
        // TODO: create it ourselves in the hopes of getting a speed up
        var drawnC = 0;
        var drawCylinder = function(geo, from, to, radius, color, fromCap, toCap) {
            if (!from || !to)
                return;
            drawnC++;
            // vertices
            var drawcaps = fromCap || toCap;
            //drawcaps = false;
            
            var dir = to.clone();
            dir.sub(from);
            
            var e = getRotationMatrix(dir);
            //get orthonormal vectors from cache
            //TODO: Will have orient with model view matrix according to direction
            var vobj = cylVertexCache.getVerticesForRadius(radius);
            
            //w (n) corresponds to the number of orthonormal vectors for cylinder (default 16)
            var n = vobj.w, h = vobj.h;
            var w = n;           
            // get orthonormal vector
            var n_verts = (drawcaps) ? h*n + 2 : 2*n;
            
            var geoGroup = geo.updateGeoGroup(n_verts);
            
            var vertices = vobj.vertices, normals = vobj.normals, verticesRows = vobj.verticesRows;
            var toRow = verticesRows[h/2], fromRow = verticesRows[h/2 + 1];
            
            var start = geoGroup.vertices;
            var offset, faceoffset;
            
            // add vertices, opposing vertices paired together
            for ( var i = 0; i < n; ++i) {
                
                var vi = 2*i;
                
                var x = e[0]*vertices[vi].x + e[3]*vertices[vi].y + e[6]*vertices[vi].z,
                    y = e[1]*vertices[vi].x + e[4]*vertices[vi].y + e[7]*vertices[vi].z,
                    z =                       e[5]*vertices[vi].y + e[8]*vertices[vi].z;
                              
                //var xn = x/radius, yn = y/radius, zn = z/radius;
                
                offset = 3*(start + vi), faceoffset = geoGroup.faceidx;
                
                //from
                geoGroup.__vertexArray[offset] = x + from.x;
                geoGroup.__vertexArray[offset+1] = y + from.y;
                geoGroup.__vertexArray[offset+2] = z + from.z;             
                //to
                geoGroup.__vertexArray[offset+3] = x + to.x;
                geoGroup.__vertexArray[offset+4] = y + to.y;
                geoGroup.__vertexArray[offset+5] = z + to.z;
                
                //normals
                geoGroup.__normalArray[offset] = x, geoGroup.__normalArray[offset+3] = x;
                geoGroup.__normalArray[offset+1] = y, geoGroup.__normalArray[offset+4] = y;
                geoGroup.__normalArray[offset+2] = z, geoGroup.__normalArray[offset+5] = z;
                
                //colors               
                geoGroup.__colorArray[offset] = color.r, geoGroup.__colorArray[offset+3] = color.r;
                geoGroup.__colorArray[offset+1] = color.g, geoGroup.__colorArray[offset+4] = color.g;
                geoGroup.__colorArray[offset+2] = color.b, geoGroup.__colorArray[offset+5] = color.b;  
                
                //faces
                // 0 - 2 - 1
                geoGroup.__faceArray[faceoffset] = fromRow[i] + start,
                geoGroup.__faceArray[faceoffset+1] = fromRow[i+1] + start,
                geoGroup.__faceArray[faceoffset+2] = toRow[i] + start;
                // 1 - 2 - 3
                geoGroup.__faceArray[faceoffset+3] = toRow[i] + start,
                geoGroup.__faceArray[faceoffset+4] = fromRow[i+1] + start,
                geoGroup.__faceArray[faceoffset+5] = toRow[i+1] + start;
                
                geoGroup.faceidx += 6;
                
            }
            
         
            //SPHERE CAPS         

            if (drawcaps) {

                // h - sphere rows, verticesRows.length - 2
                var ystart = (toCap) ? 0 : h/2;
                var yend = (fromCap) ? h + 1 : h/2+1;
                
                var v1, v2, v3, v4, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4,
                nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4, nz1, nz2, nz3, nz4;
                
                for (var y = ystart; y < yend; y++) {
                    if (y === h/2)
                        continue;
                    // n number of points for each level (verticesRows[i].length - 1)
                    var cap = (y <= h/2) ? to : from;
                    
                    for (var x = 0; x < n; x++) {
                                           
                        faceoffset = geoGroup.faceidx;
                        
                        v1 = verticesRows[y][x + 1], v1offset = (v1 + start) * 3;
                        v2 = verticesRows[y][x], v2offset = (v2 + start) * 3;
                        v3 = verticesRows[y + 1][x], v3offset = (v3 + start) * 3;
                        v4 = verticesRows[y + 1][x + 1], v4offset = (v4 + start) * 3;
                        
                        //rotate sphere vectors
                        x1 = e[0]*vertices[v1].x + e[3]*vertices[v1].y + e[6]*vertices[v1].z;
                        x2 = e[0]*vertices[v2].x + e[3]*vertices[v2].y + e[6]*vertices[v2].z;
                        x3 = e[0]*vertices[v3].x + e[3]*vertices[v3].y + e[6]*vertices[v3].z;
                        x4 = e[0]*vertices[v4].x + e[3]*vertices[v4].y + e[6]*vertices[v4].z;
                        
                        y1 = e[1]*vertices[v1].x + e[4]*vertices[v1].y + e[7]*vertices[v1].z;
                        y2 = e[1]*vertices[v2].x + e[4]*vertices[v2].y + e[7]*vertices[v2].z;
                        y3 = e[1]*vertices[v3].x + e[4]*vertices[v3].y + e[7]*vertices[v3].z;
                        y4 = e[1]*vertices[v4].x + e[4]*vertices[v4].y + e[7]*vertices[v4].z;

                        z1 =                       e[5]*vertices[v1].y + e[8]*vertices[v1].z;
                        z2 =                       e[5]*vertices[v2].y + e[8]*vertices[v2].z;
                        z3 =                       e[5]*vertices[v3].y + e[8]*vertices[v3].z;
                        z4 =                       e[5]*vertices[v4].y + e[8]*vertices[v4].z;
                        
                        geoGroup.__vertexArray[v1offset] = x1 + cap.x, geoGroup.__vertexArray[v2offset] = x2 + cap.x,
                        geoGroup.__vertexArray[v3offset] = x3 + cap.x, geoGroup.__vertexArray[v4offset] = x4 + cap.x;
    
                        geoGroup.__vertexArray[v1offset+1] = y1 + cap.y, geoGroup.__vertexArray[v2offset+1] = y2 + cap.y,
                        geoGroup.__vertexArray[v3offset+1] = y3 + cap.y, geoGroup.__vertexArray[v4offset+1] = y4 + cap.y;
    
                        geoGroup.__vertexArray[v1offset+2] = z1 + cap.z, geoGroup.__vertexArray[v2offset+2] = z2 + cap.z,
                        geoGroup.__vertexArray[v3offset+2] = z3 + cap.z, geoGroup.__vertexArray[v4offset+2] = z4 + cap.z;
    
                        geoGroup.__colorArray[v1offset] = color.r, geoGroup.__colorArray[v2offset] = color.r,
                        geoGroup.__colorArray[v3offset] = color.r, geoGroup.__colorArray[v4offset] = color.r;
    
                        geoGroup.__colorArray[v1offset+1] = color.g, geoGroup.__colorArray[v2offset+1] = color.g,
                        geoGroup.__colorArray[v3offset+1] = color.g, geoGroup.__colorArray[v4offset+1] = color.g;
    
                        geoGroup.__colorArray[v1offset+2] = color.b, geoGroup.__colorArray[v2offset+2] = color.b,
                        geoGroup.__colorArray[v3offset+2] = color.b, geoGroup.__colorArray[v4offset+2] = color.b;                                      
                             
                        
                                             
                        nx1 = e[0]*normals[v1].x + e[3]*normals[v1].y + e[6]*normals[v1].z;
                        nx2 = e[0]*normals[v2].x + e[3]*normals[v2].y + e[6]*normals[v2].z;
                        nx3 = e[0]*normals[v3].x + e[3]*normals[v3].y + e[6]*normals[v3].z;
                        nx4 = e[0]*normals[v4].x + e[3]*normals[v4].y + e[6]*normals[v4].z;
                        
                        ny1 = e[1]*normals[v1].x + e[4]*normals[v1].y + e[7]*normals[v1].z;
                        ny2 = e[1]*normals[v2].x + e[4]*normals[v2].y + e[7]*normals[v2].z;
                        ny3 = e[1]*normals[v3].x + e[4]*normals[v3].y + e[7]*normals[v3].z;
                        ny4 = e[1]*normals[v4].x + e[4]*normals[v4].y + e[7]*normals[v4].z;

                        nz1 =                       e[5]*normals[v1].y + e[8]*normals[v1].z;
                        nz2 =                       e[5]*normals[v2].y + e[8]*normals[v2].z;
                        nz3 =                       e[5]*normals[v3].y + e[8]*normals[v3].z;
                        nz4 =                       e[5]*normals[v4].y + e[8]*normals[v4].z;    
                                         
                        //if (Math.abs(vobj.sphereVertices[v1].y) === radius) {
                        if (y === 0) {
                            //face = [v1, v3, v4];
                            //norm = [n1, n3, n4];
                            
                            geoGroup.__normalArray[v1offset] = nx1, geoGroup.__normalArray[v3offset] = nx3, geoGroup.__normalArray[v4offset] = nx4;
                            geoGroup.__normalArray[v1offset+1] = ny1, geoGroup.__normalArray[v3offset+1] = ny3, geoGroup.__normalArray[v4offset+1] = ny4;
                            geoGroup.__normalArray[v1offset+2] = nz1, geoGroup.__normalArray[v3offset+2] = nz3, geoGroup.__normalArray[v4offset+2] = nz4;
    
                            geoGroup.__faceArray[faceoffset] = v1 + start; 
                            geoGroup.__faceArray[faceoffset+1] = v3 + start;
                            geoGroup.__faceArray[faceoffset+2] = v4 + start;
                            
                            geoGroup.faceidx += 3;
                            
                        } 
                        
                        //else if (Math.abs(vobj.sphereVertices[v3].y) === radius) {
                        else if (y === yend - 1) {
                            //face = [v1, v2, v3];            
                            //norm = [n1, n2, n3];
                            
                            geoGroup.__normalArray[v1offset] = nx1, geoGroup.__normalArray[v2offset] = nx2, geoGroup.__normalArray[v3offset] = nx3;
                            geoGroup.__normalArray[v1offset+1] = ny1, geoGroup.__normalArray[v2offset+1] = ny2, geoGroup.__normalArray[v3offset+1] = ny3;
                            geoGroup.__normalArray[v1offset+2] = nz1, geoGroup.__normalArray[v2offset+2] = nz2, geoGroup.__normalArray[v3offset+2] = nz3;
                            
                            geoGroup.__faceArray[faceoffset] = v1 + start;
                            geoGroup.__faceArray[faceoffset+1] = v2 + start;
                            geoGroup.__faceArray[faceoffset+2] = v3 + start;
                            
                            geoGroup.faceidx += 3;
                            
                        } 
                        
                        else {
                            //face = [v1, v2, v3, v4];
                            //norm = [n1, n2, n3, n4];
                            
                            geoGroup.__normalArray[v1offset] = nx1, geoGroup.__normalArray[v2offset] = nx2, geoGroup.__normalArray[v4offset] = nx4;
                            geoGroup.__normalArray[v1offset+1] = ny1, geoGroup.__normalArray[v2offset+1] = ny2, geoGroup.__normalArray[v4offset+1] = ny4;
                            geoGroup.__normalArray[v1offset+2] = nz1, geoGroup.__normalArray[v2offset+2] = nz2, geoGroup.__normalArray[v4offset+2] = nz4;
                            
                            geoGroup.__normalArray[v2offset] = nx2, geoGroup.__normalArray[v3offset] = nx3, geoGroup.__normalArray[v4offset] = nx4;
                            geoGroup.__normalArray[v2offset+1] = ny2, geoGroup.__normalArray[v3offset+1] = ny3, geoGroup.__normalArray[v4offset+1] = ny4;
                            geoGroup.__normalArray[v2offset+2] = nz2, geoGroup.__normalArray[v3offset+2] = nz3, geoGroup.__normalArray[v4offset+2] = nz4;
                            
                            geoGroup.__faceArray[faceoffset] = v1 + start;
                            geoGroup.__faceArray[faceoffset+1] = v2 + start;
                            geoGroup.__faceArray[faceoffset+2] = v4 + start;
                            
                            geoGroup.__faceArray[faceoffset+3] = v2 + start;
                            geoGroup.__faceArray[faceoffset+4] = v3 + start;
                            geoGroup.__faceArray[faceoffset+5] = v4 + start;
                            
                            geoGroup.faceidx += 6;
                        }
                    
                    }
                }
                                           
            }
            
            geoGroup.vertices += n_verts;
           
        };
        
        // draws cylinders and small spheres (at bond radius)
        var drawBondSticks = function(atom, atoms, geo) {
            if (!atom.style.stick)
                return;
            var style = atom.style.stick;
            if (style.hidden)
                return;

            var bondR = style.radius || defaultStickRadius;
            var fromCap = false, toCap = false;

            var c1 = atom.color;
            if (typeof (style.color) != "undefined") {
                c1 = style.color;
            }
            var C1 = WebMol.CC.color(c1);
            var mp, mp1, mp2;
            
            if (!atom.capDrawn && atom.bonds.length < 4)
                fromCap = true;              
                
            for (var i = 0; i < atom.bonds.length; i++) {
                var j = atom.bonds[i]; // our neighbor
                var atom2 = atoms[j]; //parsePDB, etc should only add defined bonds
                
                if (atom.serial < atom2.serial) {// only draw if less, this
                    // lets us combine
                    // cylinders of the same
                    // color
                    // TODO: handle bond orders
                    if (!atom2.style.stick)
                        continue; // don't sweat the details                     

                    var p1 = new WebMol.Vector3(atom.x, atom.y, atom.z);
                    var p2 = new WebMol.Vector3(atom2.x, atom2.y, atom2.z);

                    var c2 = atom2.color;
                    if (typeof (style.color) != "undefined") {
                        c2 = style.color;
                    }
                    var C2 = WebMol.CC.color(c2);

                    // draw cylinders
                    if (atom.bondOrder[i] === 1) {

                        if (!atom2.capDrawn && atom2.bonds.length < 4)
                            toCap = true;       
                                                
                        if (c1 != c2) {
                            mp = new WebMol.Vector3().addVectors(p1, p2)
                                    .multiplyScalar(0.5);
                            drawCylinder(geo, p1, mp, bondR, C1, fromCap, false);
                            drawCylinder(geo, mp, p2, bondR, C2, false, toCap);
                        } else {
                            drawCylinder(geo, p1, p2, bondR, C1, fromCap, toCap);
                        }
                        
                        if (atom.clickable || atom2.clickable) {
                            mp = new WebMol.Vector3().addVectors(p1, p2).multiplyScalar(0.5);
                            if (atom.clickable){
                                var cylinder1 = new WebMol.Cylinder(p1.clone(), mp.clone(), bondR);
                                var sphere1 = new WebMol.Sphere(p1.clone(), bondR);
                                atom.intersectionShape.cylinder.push(cylinder1);   
                                atom.intersectionShape.sphere.push(sphere1);                             
                            }
                            if (atom2.clickable){
                                var cylinder2 = new WebMol.Cylinder(p2.clone(), mp.clone(), bondR);
                                var sphere2 = new WebMol.Sphere(p2.clone(), bondR);
                                atom2.intersectionShape.cylinder.push(cylinder2);
                                atom2.intersectionShape.sphere.push(sphere2);
                            }

                        }
                        
                    } else if (atom.bondOrder[i] > 1) {
                        fromCap = false, toCap = false;
                        var dir = p2.clone();
                        var v = null;
                        dir.sub(p1);

                        if (atom.bonds.length === 1) {
                            if (atom2.bonds.length === 1) {
                                v = dir.clone();
                                if (Math.abs(v.x) > .0001)
                                    v.y += 1;
                                else
                                    v.x += 1;
                            } else {
                                var i2 = (i + 1) % atom2.bonds.length;
                                var j2 = atom2.bonds[i2];
                                var atom3 = atoms[j2];
                                var p3 = new WebMol.Vector3(atom3.x, atom3.y, atom3.z);

                                var dir2 = p3.clone();
                                dir2.sub(p1);

                                v = dir2.clone();
                                v.cross(dir);
                            }
                        } else {
                            // get vector 2 different neighboring atom
                            var i2 = (i + 1) % atom.bonds.length;
                            var j2 = atom.bonds[i2];
                            var atom3 = atoms[j2];
                            var p3 = new WebMol.Vector3(atom3.x, atom3.y, atom3.z);

                            var dir2 = p3.clone();
                            dir2.sub(p1);

                            v = dir2.clone();
                            v.cross(dir);
                        }

                        if (atom.bondOrder[i] == 2) {
                            var r = bondR / 2.5;
                            v.cross(dir);
                            v.normalize();
                            v.multiplyScalar(r * 1.5);

                            var p1a = p1.clone();
                            p1a.add(v);
                            var p1b = p1.clone();
                            p1b.sub(v);

                            var p2a = p1a.clone();
                            p2a.add(dir);
                            var p2b = p1b.clone();
                            p2b.add(dir);
                                                                 
                            if (c1 != c2) {
                                mp = new WebMol.Vector3().addVectors(p1a, p2a)
                                        .multiplyScalar(0.5);
                                mp2 = new WebMol.Vector3().addVectors(p1b, p2b)
                                        .multiplyScalar(0.5);
                                drawCylinder(geo, p1a, mp, r, C1, fromCap, false);
                                drawCylinder(geo, mp, p2a, r, C2, false, toCap);
                                drawCylinder(geo, p1b, mp2, r, C1, fromCap, false);
                                drawCylinder(geo, mp2, p2b, r, C2, false, toCap);
                            } else {
                                drawCylinder(geo, p1a, p2a, r, C1, fromCap, toCap);
                                drawCylinder(geo, p1b, p2b, r, C1, fromCap, toCap);
                            }
                            if (atom.clickable || atom2.clickable){
                                mp = new WebMol.Vector3().addVectors(p1a, p2a)
                                               .multiplyScalar(0.5);
                                mp2 = new WebMol.Vector3().addVectors(p1b, p2b)
                                                .multiplyScalar(0.5);
                                if (atom.clickable) {
                                    var cylinder1a = new WebMol.Cylinder(p1a.clone(), mp.clone(), r);
                                    var cylinder1b = new WebMol.Cylinder(p1b.clone(), mp2.clone(), r);
                                    atom.intersectionShape.cylinder.push(cylinder1a);
                                    atom.intersectionShape.cylinder.push(cylinder1b);
                                }
                                if (atom2.clickable) {
                                    var cylinder2a = new WebMol.Cylinder(p2a.clone(), mp.clone(), r);
                                    var cylinder2b = new WebMol.Cylinder(p2b.clone(), mp2.clone(), r);
                                    atom2.intersectionShape.cylinder.push(cylinder2a);
                                    atom2.intersectionShape.cylinder.push(cylinder2b);                               
                                }
                            }
                        } else if (atom.bondOrder[i] == 3) {
                            var r = bondR / 4;
                            v.cross(dir);
                            v.normalize();
                            v.multiplyScalar(r * 3);

                            var p1a = p1.clone();
                            p1a.add(v);
                            var p1b = p1.clone();
                            p1b.sub(v);

                            var p2a = p1a.clone();
                            p2a.add(dir);
                            var p2b = p1b.clone();
                            p2b.add(dir);

                            if (c1 != c2) {
                                mp = new WebMol.Vector3().addVectors(p1a, p2a)
                                        .multiplyScalar(0.5);
                                mp2 = new WebMol.Vector3().addVectors(p1b, p2b)
                                        .multiplyScalar(0.5);
                                mp3 = new WebMol.Vector3().addVectors(p1, p2)
                                        .multiplyScalar(0.5);
                                drawCylinder(geo, p1a, mp, r, C1, fromCap, false);
                                drawCylinder(geo, mp, p2a, r, C2, false, toCap);
                                drawCylinder(geo, p1, mp3, r, C1, fromCap, false);
                                drawCylinder(geo, mp3, p2, r, C2, false, toCap);
                                drawCylinder(geo, p1b, mp2, r, C1, fromCap, false);
                                drawCylinder(geo, mp2, p2b, r, C2, false, toCap);
                            } else {
                                drawCylinder(geo, p1a, p2a, r, C1, fromCap, toCap);
                                drawCylinder(geo, p1, p2, r, C1, fromCap, toCap);
                                drawCylinder(geo, p1b, p2b, r, C1, fromCap, toCap);

                            }
                            if (atom.clickable || atom2.clickable) {
                                mp = new WebMol.Vector3().addVectors(p1a, p2a)
                                        .multiplyScalar(0.5);
                                mp2 = new WebMol.Vector3().addVectors(p1b, p2b)
                                        .multiplyScalar(0.5);
                                mp3 = new WebMol.Vector3().addVectors(p1, p2)
                                        .multiplyScalar(0.5);
                                
                                if (atom.clickable) {
                                    var cylinder1a = new WebMol.Cylinder(p1a.clone(), mp.clone(), r);
                                    var cylinder1b = new WebMol.Cylinder(p1b.clone(), mp2.clone(), r);
                                    var cylinder1c = new WebMol.Cylinder(p1.clone(), mp3.clone(), r);
                                    atom.intersectionShape.cylinder.push(cylinder1a);
                                    atom.intersectionShape.cylinder.push(cylinder1b);
                                    atom.intersectionShape.cylinder.push(cylinder1c);
                                } 
                                if (atom2.clickable) {                               
                                    var cylinder2a = new WebMol.Cylinder(p2a.clone(), mp.clone(), r);
                                    var cylinder2b = new WebMol.Cylinder(p2b.clone(), mp2.clone(), r);
                                    var cylinder2c = new WebMol.Cylinder(p2.clone(), mp3.clone(), r);
                                    atom2.intersectionShape.cylinder.push(cylinder2a);
                                    atom2.intersectionShape.cylinder.push(cylinder2b);
                                    atom2.intersectionShape.cylinder.push(cylinder2c);                                
                                }
                            }
                        }
                    }
                    if (toCap || atom2.bonds.length > 3)
                        atom2.capDrawn = true;
                    if (fromCap || atom.bonds.length > 3)
                        atom.capDrawn = true;                        
                }                   
                                 
            }            

            // draw non bonded heteroatoms as spheres   
            var drawSphere = (!atom.singleBonds && atom.bonds.length === 1) || (!atom.bonds.length);
            if (drawSphere) {
                var savedstyle = atom.style;
                atom.style = {
                    sphere : {
                        radius : bondR,
                        color : c1
                    }
                };
                drawAtomSphere(atom, geo);                
                atom.style = savedstyle;
            }
            
        };

        // go through all the atoms and regenerate their geometries
        // we try to have one geometry for each style since this is much much
        // faster
        // at some point we should optimize this to avoid unnecessary
        // recalculation
        var createMolObj = function(atoms) {

            console.log("creating for "+id);
            var ret = new WebMol.Object3D();
            var cartoonAtoms = [];
            var lineGeometries = {};
            var crossGeometries = {};
            var sphereGeometry = new WebMol.Geometry(true);                                                         
            var stickGeometry = new WebMol.Geometry(true);
            
            for ( var i = 0; i < atoms.length; i++) {
                var atom = atoms[i];
                // recreate gl info for each atom as necessary
                // set up appropriate intersection spheres for clickable atoms
                if (atom && atom.style) {
                    if (atom.clickable && atom.intersectionShape === undefined)
                        atom.intersectionShape = {sphere: [], cylinder: [], line: [], triangle : []};                    
                    drawAtomSphere(atom, sphereGeometry);
                    drawAtomCross(atom, crossGeometries);
                    drawBondLines(atom, atoms, lineGeometries);
                    drawBondSticks(atom, atoms, stickGeometry);
                    if (typeof (atom.style.cartoon) !== "undefined"
                            && !atom.style.cartoon.hidden) {
                        cartoonAtoms.push(atom);
                    }

                }
            }
            // create cartoon if needed - this is a whole model analysis
            if (cartoonAtoms.length > 0) {
                WebMol.drawCartoon(ret, cartoonAtoms, false);
                
                for (var i = 0; i < ret.children.length; i++){
                    var geo = ret.children[i].geometry;
                }
            }

            // add sphere geometry
            if (sphereGeometry.vertices > 0) {
                var sphereMaterial = new WebMol.MeshLambertMaterial({
                    ambient : 0x000000,
                    vertexColors : true,
                    reflectivity : 0
                });
                
                //Initialize buffers in geometry                
                sphereGeometry.initTypedArrays();
                
                var sphere = new WebMol.Mesh(sphereGeometry, sphereMaterial);
                console
                        .log("sphere geometry "
                                + sphereGeometry.vertices.length);

                ret.add(sphere);
            }

            // add stick geometry
            if (stickGeometry.vertices > 0) {
                var cylinderMaterial = new WebMol.MeshLambertMaterial({
                    vertexColors : true,
                    ambient : 0x000000,
                    reflectivity : 0
                });

                //Initialize buffers in geometry                
                stickGeometry.initTypedArrays();
                
                if (cylinderMaterial.wireframe)
                    stickGeometry.setUpWireframe();
                
                var sticks = new WebMol.Mesh(stickGeometry, cylinderMaterial);
                ret.add(sticks);
            }

            // add any line geometries, distinguished by line width
            for ( var i in lineGeometries) {
                if (lineGeometries.hasOwnProperty(i)) {
                    var linewidth = i;
                    var lineMaterial = new WebMol.LineBasicMaterial({
                        linewidth : linewidth,
                        vertexColors : true
                    });
                    
                    lineGeometries[i].initTypedArrays();
                    
                    var line = new WebMol.Line(lineGeometries[i], lineMaterial,
                            WebMol.LinePieces);

                    ret.add(line);
                }
            }

            // add any cross geometries
            for ( var i in crossGeometries) {
                if (crossGeometries.hasOwnProperty(i)) {
                    var linewidth = i;
                    var lineMaterial = new WebMol.LineBasicMaterial({
                        linewidth : linewidth,
                        vertexColors : true
                    });

                    crossGeometries[i].initTypedArrays();
                    
                    var line = new WebMol.Line(crossGeometries[i], lineMaterial,
                            WebMol.LinePieces);

                    ret.add(line);
                }
            }

            return ret;
        };

        this.getID = function() {
            return id;
        };

        // set default style and colors for atoms
        var setAtomDefaults = function(atoms, id) {
            for ( var i = 0; i < atoms.length; i++) {
                var atom = atoms[i];
                if (atom) {
                    atom.style = atom.style || defaultAtomStyle;
                    atom.color = atom.color || ElementColors[atom.elem]
                            || defaultColor;
                    atom.model = id;
                    if (atom.clickable)
                        atom.intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};
                }
            }
        };

        // add atoms to this model from molecular data string
        this.addMolData = function(data, format) {
            
            if (!data)
                console.error("Erorr with addMolData: No input data specified");
            
            switch (format) {
            case "xyz":
                parseXYZ(atoms, data);
                break;
            case "pdb":
                parsePDB(atoms, data, false, true);
                break;
            case "sdf":
                parseSDF(atoms, data);
                break;
            case "mol2":
                parseMOL2(atoms, data);
                break;
            case "cube":
                parseCube(atoms, data);
                break;
            }
            setAtomDefaults(atoms, id);
        };
        
        // given a selection specification, return true if atom is selected
        this.atomIsSelected = function(atom, sel) {
            if (typeof (sel) === "undefined")
                return true; // undef gets all
            var invert = !!sel.invert;
            var ret = true;
            for ( var key in sel) {
                if (sel.hasOwnProperty(key) && key != "props" && key != "invert" && key != "model") {
                    // if something is in sel, atom must have it
                    if (typeof (atom[key]) === "undefined") {
                        ret = false;
                        break;
                    }
                    var isokay = false;
                    if ($.isArray(sel[key])) {
                        // can be any of the listed values
                        var valarr = sel[key];
                        for ( var i = 0; i < valarr.length; i++) {
                            if (atom[key] == valarr[i]) {
                                isokay = true;
                                break;
                            }
                        }
                        if (!isokay) {
                            ret = false;
                            break;
                        }
                    } else { // single match
                        var val = sel[key];
                        if (atom[key] != val) {
                            ret = false;
                            break;
                        }
                    }
                }
            }
            
            return invert ? !ret : ret;
        };

        // return list of atoms selected by sel, this is specific to glmodel
        this.selectedAtoms = function(sel) {
            var ret = [];
            for ( var i = 0; i < atoms.length; i++) {
                var atom = atoms[i];
                if (atom) {
                    if (this.atomIsSelected(atom, sel))
                        ret.push(atom);
                }
            }
            return ret;
        };
        
        // copy new atoms into this model, adjust bonds appropriately
        this.addAtoms = function(newatoms) {
            molObj = null;
            var start = atoms.length;
            var indexmap = [];
            // mapping from old index to new index
            for(var i = 0; i < newatoms.length; i++) {
                indexmap[newatoms[i].index] = start+i;
            }
            
            // copy and push newatoms onto atoms
            for(var i = 0; i < newatoms.length; i++) {
                var olda = newatoms[i];
                var nindex = indexmap[olda.index];
                var a = $.extend(false, {}, olda);
                a.index = nindex;
                a.bonds = [];
                a.bondOrder = [];
                // copy over all bonds contained in selection,
                // updating indices appropriately
                for(var j = 0; j < olda.bonds.length; j++) {
                    var neigh = indexmap[olda.bonds[j]];
                    if(typeof(neigh) != "undefined") {
                        a.bonds.push(neigh);
                        a.bondOrder.push(olda.bondOrder[j]);
                    }                
                }
                atoms.push(a);
            }
        };

        // remove badatoms from model
        this.removeAtoms = function(badatoms) {
            molObj = null;
            // make map of all baddies
            var baddies = [];
            for(var i = 0; i < badatoms.length; i++) {
                baddies[badatoms[i].index] = true;
            }
            
            // create list of good atoms
            var newatoms = [];
            for(var i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                if(!baddies[a.index])
                    newatoms.push(a);
            }
            
            // clear it all out
            atoms = [];
            // and add back in to get updated bonds
            this.addAtoms(newatoms);
        };
        
        
        // style the select atoms with style
        this.setStyle = function(sel, style, add) {
            
            if(!add && molObj != null && sameObj(style, lastStyle))
                return; // no need to recompute
            
            if(add) lastStyle = null; // todo: compute merged style
            else lastStyle = style;
            
            var atoms = this.selectedAtoms(sel);
            if(atoms.length > 0)
                molObj = null; // force rebuild
            // do a copy to enforce style changes through this function
            var mystyle = $.extend(true, {}, style);

            // somethings we only calculate if there is a change in a certain
            // style, although these checks will only catch cases where both
            // are either null or undefined
            for ( var i = 0; i < atoms.length; i++) {
                
                if (atoms[i].clickable) 
                    atoms[i].intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};                    
                
                atoms[i].capDrawn = false; //reset for proper stick render
               
                if(!add) atoms[i].style = {};
                for(var s in mystyle) {
                    if(mystyle.hasOwnProperty(s)) {
                        atoms[i].style[s] = mystyle[s];
                    }
                }
            }
        };
        
        // given a mapping from element to color, set atom colors
        this.setColorByElement = function(sel, colors) {
            
            if(molObj !== null && sameObj(colors,lastColors))
                return; // don't recompute
            lastColors = colors;
            var atoms = this.selectedAtoms(sel);
            if(atoms.length > 0)
                molObj = null; // force rebuild
            for ( var i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                if(typeof(colors[a.elem]) !== "undefined") {
                    a.color = colors[a.elem];
                }
            }
        };
        
        this.setColorByProperty = function(sel, prop, scheme) {
            var atoms = this.selectedAtoms(sel);
            lastColors = null; // don't bother memoizing
            if(atoms.length > 0)
                molObj = null; // force rebuild
            var min =  Number.POSITIVE_INFINITY;
            var max =  Number.NEGATIVE_INFINITY;
            
            // compute the range
            for ( var i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                if(a.properties && typeof(a.properties[prop]) != undefined) {
                    var p = parseFloat(a.properties[prop]);
                    if(p < min) min = p;
                    if(p > max) max = p;
                }                    
            }
            // now apply colors using scheme
            for ( var i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                if(a.properties && typeof(a.properties[prop]) != undefined) {
                    var c = scheme.valueToHex(parseFloat(a.properties[prop]), [min,max]);
                    a.color = c;
                }                    
            }
        };


        // manage the globj for this model in the possed modelGroup -
        // if it has to be regenerated, remove and add

        this.globj = function(group) {
            var time = new Date();
            if(molObj === null) { // have to regenerate
                molObj = createMolObj(atoms);
                var time2 = new Date();
                console.log("object creation time: " + (time2 - time));
                if(renderedMolObj) { // previously rendered, remove
                    group.remove(renderedMolObj);
                    renderedMolObj = null;
                }
                renderedMolObj = molObj.clone();
                group.add(renderedMolObj);
            }
        };
        
        // remove any rendered object from the scene
        this.removegl = function(group) {
            if(renderedMolObj) {
                //dispose of geos and materials
                if (renderedMolObj.geometry !== undefined) renderedMolObj.geometry.dispose();             
                if (renderedMolObj.material !== undefined) renderedMolObj.material.dispose();
                group.remove(renderedMolObj);
                renderedMolObj = null;
            }
            molObj = null;
        };

    };
    
    GLModel.prototype.testMethod = function() {
          
    };

    return GLModel;
    
})();
//A GLShape is a collection of user specified shapes. Includes
// build in sphere and arrow shapes, as well as custom user specified shapes

WebMol.GLShape = (function() {
    
    //Marching cube, to match with protein surface generation
    var ISDONE = 2;

    var updateColor = function(geo, color) {
        
        var C = color || WebMol.CC.color(color);
        geo.colorsNeedUpdate = true;
        
        for (var g in geo.geometryGroups) {
            
            var geoGroup = geo.geometryGroups[g];
            var colorArr = geoGroup.__colorArray;
            
            for (var i = 0, il = geoGroup.vertices; i < il; ++i) {
                colorArr[i*3] = C.r, colorArr[i*3+1] = C.g, colorArr[i*3+2] = C.b;                       
            }
        }
        
    };
    
    //Preset component builders
    
    //Sphere component 
    var sphereVertexCache = {
        cache : {},
        getVerticesForRadius : function(radius) {

            if (typeof (this.cache[radius]) !== "undefined")
                return this.cache[radius];

            var obj = {
                vertices : [],
                verticesRows : [],
                normals : []
            };
            // scale quality with radius heuristically
            var widthSegments = 16;
            var heightSegments = 10;
            if (radius < 1) {
                widthSegments = 8;
                heightSegments = 6;
            }

            var phiStart = 0;
            var phiLength = Math.PI * 2;

            var thetaStart = 0;
            var thetaLength = Math.PI;

            var x, y, vertices = [], uvs = [];

            for (y = 0; y <= heightSegments; y++) {

                var verticesRow = [];
                for (x = 0; x <= widthSegments; x++) {

                    var u = x / widthSegments;
                    var v = y / heightSegments;

                    var vertex = {};
                    vertex.x = -radius * Math.cos(phiStart + u * phiLength)
                            * Math.sin(thetaStart + v * thetaLength);
                    vertex.y = radius
                            * Math.cos(thetaStart + v * thetaLength);
                    vertex.z = radius * Math.sin(phiStart + u * phiLength)
                            * Math.sin(thetaStart + v * thetaLength);

                    var n = new WebMol.Vector3(vertex.x, vertex.y, vertex.z);
                    n.normalize();

                    obj.vertices.push(vertex);
                    obj.normals.push(n);

                    verticesRow.push(obj.vertices.length - 1);

                }

                obj.verticesRows.push(verticesRow);

            }

            this.cache[radius] = obj;
            return obj;
        }
        
    }; 
    
    var drawSphere = function(shape, geoGroup, spec) {
        
        var pos = spec.center, radius = spec.radius;        
        
        var center = new WebMol.Vector3(pos.x, pos.y, pos.z);
        shape.intersectionShape.sphere.push( new WebMol.Sphere(center, radius) );                                                                  

        var x, y;
        var vobj = sphereVertexCache.getVerticesForRadius(radius);                
                    
        var vertices = vobj.vertices;
        var normals = vobj.normals;
        
        var start = geoGroup.vertices;
        
        for (var i = 0, il = vertices.length; i < il; ++i) {
            var offset = 3*(start + i);   
            var v = vertices[i];
            
            geoGroup.__vertexArray[offset] = (v.x + pos.x);
            geoGroup.__vertexArray[offset+1] = (v.y + pos.y);
            geoGroup.__vertexArray[offset+2] = (v.z + pos.z);            
           
        }
        
        geoGroup.vertices += vertices.length;
        
        var verticesRows = vobj.verticesRows;
        var h = verticesRows.length - 1;
        
        //var color = [C.r, C.g, C.b];
        for (y = 0; y < h; y++) {
            var w = verticesRows[y].length - 1;
            for (x = 0; x < w; x++) {
                
                var faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;
                
                var v1 = verticesRows[y][x + 1] + start, v1offset = v1 * 3;
                var v2 = verticesRows[y][x] + start, v2offset = v2 * 3;
                var v3 = verticesRows[y + 1][x] + start, v3offset = v3 * 3;
                var v4 = verticesRows[y + 1][x + 1] + start, v4offset = v4 * 3;

                var n1 = normals[v1 - start];
                var n2 = normals[v2 - start];
                var n3 = normals[v3 - start];
                var n4 = normals[v4 - start];
                var face, norm;
                if (Math.abs(vertices[v1 - start].y) === radius) {
                    //face = [v1, v3, v4];
                    //norm = [n1, n3, n4];
                    
                    geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v3offset] = n3.x, geoGroup.__normalArray[v4offset] = n4.x;
                    geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v3offset+1] = n3.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                    geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v3offset+2] = n3.z, geoGroup.__normalArray[v4offset+2] = n4.z;

                    geoGroup.__faceArray[faceoffset] = v1; 
                    geoGroup.__faceArray[faceoffset+1] = v3;
                    geoGroup.__faceArray[faceoffset+2] = v4;
                    
                    geoGroup.__lineArray[lineoffset] = v1, geoGroup.__lineArray[lineoffset+1] = v3, geoGroup.__lineArray[lineoffset+2] = v1;
                    geoGroup.__lineArray[lineoffset+3] = v4, geoGroup.__lineArray[lineoffset+4] = v3, geoGroup.__lineArray[lineoffset+5] = v4;
                    
                    geoGroup.faceidx += 3;
                    geoGroup.lineidx += 6;
                    
                } else if (Math.abs(vertices[v3 - start].y) === radius) {
                    //face = [v1, v2, v3];            
                    //norm = [n1, n2, n3];
                    
                    geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v3offset] = n3.x;
                    geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v3offset+1] = n3.y;
                    geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v3offset+2] = n3.z;

                    geoGroup.__faceArray[faceoffset] = v1;
                    geoGroup.__faceArray[faceoffset+1] = v2;
                    geoGroup.__faceArray[faceoffset+2] = v3;
                    
                    geoGroup.__lineArray[lineoffset] = v1, geoGroup.__lineArray[lineoffset+1] = v2, geoGroup.__lineArray[lineoffset+2] = v1;
                    geoGroup.__lineArray[lineoffset+3] = v3, geoGroup.__lineArray[lineoffset+4] = v2, geoGroup.__lineArray[lineoffset+5] = v3;
                    
                    geoGroup.faceidx += 3;
                    geoGroup.lineidx += 6;
                    
                } else {
                    //face = [v1, v2, v3, v4];
                    //norm = [n1, n2, n3, n4];
                    
                    geoGroup.__normalArray[v1offset] = n1.x, geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v4offset] = n4.x;
                    geoGroup.__normalArray[v1offset+1] = n1.y, geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                    geoGroup.__normalArray[v1offset+2] = n1.z, geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v4offset+2] = n4.z;
                    
                    geoGroup.__normalArray[v2offset] = n2.x, geoGroup.__normalArray[v3offset] = n3.x, geoGroup.__normalArray[v4offset] = n4.x;
                    geoGroup.__normalArray[v2offset+1] = n2.y, geoGroup.__normalArray[v3offset+1] = n3.y, geoGroup.__normalArray[v4offset+1] = n4.y;
                    geoGroup.__normalArray[v2offset+2] = n2.z, geoGroup.__normalArray[v3offset+2] = n3.z, geoGroup.__normalArray[v4offset+2] = n4.z;
                    
                    geoGroup.__faceArray[faceoffset] = v1;
                    geoGroup.__faceArray[faceoffset+1] = v2;
                    geoGroup.__faceArray[faceoffset+2] = v4;
                    
                    geoGroup.__faceArray[faceoffset+3] = v2;
                    geoGroup.__faceArray[faceoffset+4] = v3;
                    geoGroup.__faceArray[faceoffset+5] = v4;
                    
                    geoGroup.__lineArray[lineoffset] = v1, geoGroup.__lineArray[lineoffset+1] = v2,
                    geoGroup.__lineArray[lineoffset+2] = v1, geoGroup.__lineArray[lineoffset+3] = v4;
                    
                    geoGroup.__lineArray[lineoffset+4] = v2, geoGroup.__lineArray[lineoffset+5] = v3,
                    geoGroup.__lineArray[lineoffset+6] = v3, geoGroup.__lineArray[lineoffset+7] = v4;
                    
                    geoGroup.faceidx += 6;
                    geoGroup.lineidx += 8;
                    
                }

            }
        }

    };
    
    var drawArrow = function(shape, geoGroup, spec) {
        
        var from = spec.start, end = spec.end, radius = spec.radius, radiusRatio = spec.radiusRatio, mid = spec.mid;

        if (!from || !end)
            return;
        
        // vertices
        
        var dir = end.clone();
        dir.sub(from).multiplyScalar(mid);
        var to = from.clone().add(dir);
        var negDir = dir.clone().negate();
        
        shape.intersectionShape.cylinder.push( new WebMol.Cylinder(from.clone(), to.clone(), radius) );
        shape.intersectionShape.sphere.push( new WebMol.Sphere(from.clone(), radius) );
        
        // get orthonormal vector
        var nvecs = [];
        nvecs[0] = dir.clone();
        if (Math.abs(nvecs[0].x) > .0001)
            nvecs[0].y += 1;
        else
            nvecs[0].x += 1;
        nvecs[0].cross(dir);
        nvecs[0].normalize();

        nvecs[0] = nvecs[0];
        // another orth vector
        nvecs[4] = nvecs[0].clone();
        nvecs[4].crossVectors(nvecs[0], dir);
        nvecs[4].normalize();
        nvecs[8] = nvecs[0].clone().negate();
        nvecs[12] = nvecs[4].clone().negate();

        // now quarter positions
        nvecs[2] = nvecs[0].clone().add(nvecs[4]).normalize();
        nvecs[6] = nvecs[4].clone().add(nvecs[8]).normalize();
        nvecs[10] = nvecs[8].clone().add(nvecs[12]).normalize();
        nvecs[14] = nvecs[12].clone().add(nvecs[0]).normalize();

        // eights
        nvecs[1] = nvecs[0].clone().add(nvecs[2]).normalize();
        nvecs[3] = nvecs[2].clone().add(nvecs[4]).normalize();
        nvecs[5] = nvecs[4].clone().add(nvecs[6]).normalize();
        nvecs[7] = nvecs[6].clone().add(nvecs[8]).normalize();
        nvecs[9] = nvecs[8].clone().add(nvecs[10]).normalize();
        nvecs[11] = nvecs[10].clone().add(nvecs[12]).normalize();
        nvecs[13] = nvecs[12].clone().add(nvecs[14]).normalize();
        nvecs[15] = nvecs[14].clone().add(nvecs[0]).normalize();

        //var start = geo.vertices.length;
        var start = geoGroup.vertices;
        // add vertices, opposing vertices paired together
        for ( var i = 0, n = nvecs.length; i < n; ++i) {
            var offset = 3*(start + 3*i);
            var bottom = nvecs[i].clone().multiplyScalar(radius).add(from);
            var top = nvecs[i].clone().multiplyScalar(radius).add(to);
            var conebase = nvecs[i].clone().multiplyScalar(radius*radiusRatio).add(to);

            geoGroup.__vertexArray[offset] = bottom.x;
            geoGroup.__vertexArray[offset+1] = bottom.y;
            geoGroup.__vertexArray[offset+2] = bottom.z;             
            
            geoGroup.__vertexArray[offset+3] = top.x;
            geoGroup.__vertexArray[offset+4] = top.y;
            geoGroup.__vertexArray[offset+5] = top.z; 
            
            geoGroup.__vertexArray[offset+6] = conebase.x;
            geoGroup.__vertexArray[offset+7] = conebase.y;
            geoGroup.__vertexArray[offset+8] = conebase.z;
            
            if (i > 0) {
                var prev_x = geoGroup.__vertexArray[offset-3];
                var prev_y = geoGroup.__vertexArray[offset-2];
                var prev_z = geoGroup.__vertexArray[offset-1];
                
                var c = new WebMol.Vector3(prev_x, prev_y, prev_z);
                var b = end.clone(), b2 = to.clone();
                var a = new WebMol.Vector3(conebase.x, conebase.y, conebase.z);
                
                shape.intersectionShape.triangle.push( new WebMol.Triangle(a, b, c) );
                shape.intersectionShape.triangle.push( new WebMol.Triangle(c.clone(), b2, a.clone()) );
            }
        }
        
        geoGroup.vertices += 48;
        offset = geoGroup.vertices*3;
        
        //caps
        geoGroup.__vertexArray[offset] = from.x;
        geoGroup.__vertexArray[offset+1] = from.y;
        geoGroup.__vertexArray[offset+2] = from.z;
        
        geoGroup.__vertexArray[offset+3] = to.x;
        geoGroup.__vertexArray[offset+4] = to.y;
        geoGroup.__vertexArray[offset+5] = to.z;
        
        geoGroup.__vertexArray[offset+6] = end.x;
        geoGroup.__vertexArray[offset+7] = end.y;
        geoGroup.__vertexArray[offset+8] = end.z;
        
        geoGroup.vertices += 3;
        
        // now faces
        var face, norm, offset, faceoffset, lineoffset;
        var n_vertices = 0;
        var fromi = geoGroup.vertices - 3, toi = geoGroup.vertices - 2, endi = geoGroup.vertices - 1;
        var fromoffset = fromi*3, tooffset = toi*3, endoffset = endi*3;
        for (var i = 0, n = nvecs.length - 1; i < n; ++i) {
        
            var ti = start + 3 * i, offset = ti * 3;
            faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;
            
            var t1 = ti, t1offset = t1 * 3;
            var t2 = ti + 1, t2offset = t2 * 3;
            var t2b = ti + 2, t2boffset = t2b * 3;
            var t3 = ti + 4, t3offset = t3 * 3;
            var t3b = ti + 5, t3boffset = t3b * 3;
            var t4 = ti + 3, t4offset = t4 * 3;
            
            //face = [t1, t2, t4], [t2, t3, t4];    
            //face = [t1, t2, t3, t4];
                
            norm = [ nvecs[i], nvecs[i], nvecs[i + 1], nvecs[i + 1]];
            var n1, n2, n3, n4;
            n1 = n2 = nvecs[i];
            n3 = n4 = nvecs[i + 1];
            
            geoGroup.__normalArray[t1offset] = n1.x, geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t4offset] = n4.x;
            geoGroup.__normalArray[t1offset+1] = n1.y, geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t4offset+1] = n4.y;
            geoGroup.__normalArray[t1offset+2] = n1.z, geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t4offset+2] = n4.z;
            
            geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t3offset] = n3.x, geoGroup.__normalArray[t4offset] = n4.x;
            geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t3offset+1] = n3.y, geoGroup.__normalArray[t4offset+1] = n4.y;
            geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t3offset+2] = n3.z, geoGroup.__normalArray[t4offset+2] = n4.z;
            
            geoGroup.__normalArray[t2boffset] = n2.x, geoGroup.__normalArray[t3boffset] = n3.x;
            geoGroup.__normalArray[t2boffset+1] = n2.y, geoGroup.__normalArray[t3boffset+1] = n3.y;
            geoGroup.__normalArray[t2boffset+2] = n2.z, geoGroup.__normalArray[t3boffset+2] = n3.z;
            
            //sides
            geoGroup.__faceArray[faceoffset] = t1, geoGroup.__faceArray[faceoffset+1] = t2, geoGroup.__faceArray[faceoffset+2] = t4;
            geoGroup.__faceArray[faceoffset+3] = t2, geoGroup.__faceArray[faceoffset+4] = t3, geoGroup.__faceArray[faceoffset+5] = t4;
            //caps
            geoGroup.__faceArray[faceoffset+6] = t1, geoGroup.__faceArray[faceoffset+7] = t4, geoGroup.__faceArray[faceoffset+8] = fromi;
            geoGroup.__faceArray[faceoffset+9] = t2b, geoGroup.__faceArray[faceoffset+10] = toi, geoGroup.__faceArray[faceoffset+11] = t3b;
            //arrowhead
            geoGroup.__faceArray[faceoffset+12] = t2b, geoGroup.__faceArray[faceoffset+13] = endi, geoGroup.__faceArray[faceoffset+14] = t3b;
            
            //sides
            geoGroup.__lineArray[lineoffset] = t1, geoGroup.__lineArray[lineoffset+1] = t2;
            geoGroup.__lineArray[lineoffset+2] = t1, geoGroup.__lineArray[lineoffset+3] = t4;           
            //geoGroup.__lineArray[lineoffset+4] = t2, geoGroup.__lineArray[lineoffset+5] = t3;
            geoGroup.__lineArray[lineoffset+4] = t3, geoGroup.__lineArray[lineoffset+5] = t4;
            //caps
            geoGroup.__lineArray[lineoffset+6] = t1, geoGroup.__lineArray[lineoffset+7] = t4;
            //geoGroup.__lineArray[lineoffset+10] = t1, geoGroup.__lineArray[lineoffset+11] = fromi;           
            //geoGroup.__lineArray[lineoffset+12] = t4, geoGroup.__lineArray[lineoffset+13] = fromi;
            
            geoGroup.__lineArray[lineoffset+8] = t2b, geoGroup.__lineArray[lineoffset+9] = t2; //toi   
            geoGroup.__lineArray[lineoffset+10] = t2b, geoGroup.__lineArray[lineoffset+11] = t3b;
            geoGroup.__lineArray[lineoffset+12] = t3, geoGroup.__lineArray[lineoffset+13] = t3b; //toi
            //arrowhead
            geoGroup.__lineArray[lineoffset+14] = t2b, geoGroup.__lineArray[lineoffset+15] = endi;
            geoGroup.__lineArray[lineoffset+16] = t2b, geoGroup.__lineArray[lineoffset+17] = t3b;
            geoGroup.__lineArray[lineoffset+18] = endi, geoGroup.__lineArray[lineoffset+19] = t3b;
                             
            geoGroup.faceidx += 15;
            geoGroup.lineidx += 20;
            
        }
        // final face

        face = [start + 45, start + 46, start + 1, start, start + 47, start + 2];
        norm = [ nvecs[15], nvecs[15], nvecs[0], nvecs[0] ];
        
        faceoffset = geoGroup.faceidx, lineoffset = geoGroup.lineidx;
        
        var t1 = face[0], t1offset = t1 * 3;
        var t2 = face[1], t2offset = t2 * 3;
        var t2b = face[4], t2boffset = t2b * 3;
        var t3 = face[2], t3offset = t3 * 3;
        var t3b = face[5], t3boffset = t3b * 3;
        var t4 = face[3], t4offset = t4 * 3;        
        var n1, n2, n3, n4;
        
        n1 = n2 = nvecs[15];
        n3 = n4 = nvecs[0];
                
        geoGroup.__normalArray[t1offset] = n1.x, geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t4offset] = n4.x;
        geoGroup.__normalArray[t1offset+1] = n1.y, geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t4offset+1] = n4.y;
        geoGroup.__normalArray[t1offset+2] = n1.z, geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t4offset+2] = n4.z;
        
        geoGroup.__normalArray[t2offset] = n2.x, geoGroup.__normalArray[t3offset] = n3.x, geoGroup.__normalArray[t4offset] = n4.x;
        geoGroup.__normalArray[t2offset+1] = n2.y, geoGroup.__normalArray[t3offset+1] = n3.y, geoGroup.__normalArray[t4offset+1] = n4.y;
        geoGroup.__normalArray[t2offset+2] = n2.z, geoGroup.__normalArray[t3offset+2] = n3.z, geoGroup.__normalArray[t4offset+2] = n4.z;
        
        geoGroup.__normalArray[t2boffset] = n2.x, geoGroup.__normalArray[t3boffset] = n3.x;
        geoGroup.__normalArray[t2boffset+1] = n2.y, geoGroup.__normalArray[t3boffset+1] = n3.y;
        geoGroup.__normalArray[t2boffset+2] = n2.z, geoGroup.__normalArray[t3boffset+2] = n3.z;
        
        //Cap normals
        dir.normalize(), negDir.normalize();
        geoGroup.__normalArray[fromoffset] = negDir.x, geoGroup.__normalArray[tooffset] = geoGroup.__normalArray[endoffset] = dir.x;
        geoGroup.__normalArray[fromoffset+1] = negDir.y, geoGroup.__normalArray[tooffset+1] = geoGroup.__normalArray[endoffset+1] = dir.y;
        geoGroup.__normalArray[fromoffset+2] = negDir.z, geoGroup.__normalArray[tooffset+2] = geoGroup.__normalArray[endoffset+2] = dir.z;
        
        //Final side
        geoGroup.__faceArray[faceoffset] = t1, geoGroup.__faceArray[faceoffset+1] = t2, geoGroup.__faceArray[faceoffset+2] = t4;
        geoGroup.__faceArray[faceoffset+3] = t2, geoGroup.__faceArray[faceoffset+4] = t3, geoGroup.__faceArray[faceoffset+5] = t4;
        //final caps
        geoGroup.__faceArray[faceoffset+6] = t1, geoGroup.__faceArray[faceoffset+7] = t4, geoGroup.__faceArray[faceoffset+8] = fromi;
        geoGroup.__faceArray[faceoffset+9] = t2b, geoGroup.__faceArray[faceoffset+10] = toi, geoGroup.__faceArray[faceoffset+11] = t3b;
        //final arrowhead
        geoGroup.__faceArray[faceoffset+12] = t2b, geoGroup.__faceArray[faceoffset+13] = endi, geoGroup.__faceArray[faceoffset+14] = t3b;
        
        //sides
        geoGroup.__lineArray[lineoffset] = t1, geoGroup.__lineArray[lineoffset+1] = t2;
        geoGroup.__lineArray[lineoffset+2] = t1, geoGroup.__lineArray[lineoffset+3] = t4;           
        //geoGroup.__lineArray[lineoffset+4] = t2, geoGroup.__lineArray[lineoffset+5] = t3;
        geoGroup.__lineArray[lineoffset+4] = t3, geoGroup.__lineArray[lineoffset+5] = t4;
        //caps
        geoGroup.__lineArray[lineoffset+6] = t1, geoGroup.__lineArray[lineoffset+7] = t4;
        //geoGroup.__lineArray[lineoffset+10] = t1, geoGroup.__lineArray[lineoffset+11] = fromi;           
        //geoGroup.__lineArray[lineoffset+12] = t4, geoGroup.__lineArray[lineoffset+13] = fromi;

        geoGroup.__lineArray[lineoffset+8] = t2b, geoGroup.__lineArray[lineoffset+9] = t2; //toi        
        geoGroup.__lineArray[lineoffset+10] = t2b, geoGroup.__lineArray[lineoffset+11] = t3b;
        geoGroup.__lineArray[lineoffset+12] = t3, geoGroup.__lineArray[lineoffset+13] = t3b; //toi
        //arrowhead
        geoGroup.__lineArray[lineoffset+14] = t2b, geoGroup.__lineArray[lineoffset+15] = endi;
        geoGroup.__lineArray[lineoffset+16] = t2b, geoGroup.__lineArray[lineoffset+17] = t3b;
        geoGroup.__lineArray[lineoffset+18] = endi, geoGroup.__lineArray[lineoffset+19] = t3b; 
        
        geoGroup.faceidx += 15;        
        geoGroup.lineidx += 20;    

        
    };
    
    //handles custom shape generation from user supplied arrays
    //May need to generate normal and/or line indices
    var drawCustom = function(shape, geoGroup, customSpec) {
        
        var vertexArr = customSpec.vertexArr, normalArr = customSpec.normalArr, faceArr = customSpec.faceArr, lineArr = customSpec.lineArr;        
        
        if (vertexArr.length === 0 || faceArr.length === 0) {
            console.warn("Error adding custom shape component: No vertices and/or face indices supplied!");
        }
        
        geoGroup.vertices = vertexArr.length, geoGroup.faceidx = faceArr.length;
        
        var offset, v, a, b, c;
        
        for (var i = 0; i < geoGroup.vertices; ++i) {            
            offset = i*3;
            v = vertexArr[i];    
            geoGroup.__vertexArray[offset] = v.x, geoGroup.__vertexArray[offset+1] = v.y, geoGroup.__vertexArray[offset+2] = v.z;           
        }
        
        for (var i = 0; i < geoGroup.faceidx / 3; ++i) {
            offset = i*3;
            a = faceArr[offset], b = faceArr[offset+1], c = faceArr[offset+2];
            var vA = new WebMol.Vector3(), vB = new WebMol.Vector3(), vC = new WebMol.Vector3();
            shape.intersectionShape.triangle.push( new WebMol.Triangle( vA.copy(vertexArr[a]), vB.copy(vertexArr[b]), vC.copy(vertexArr[c]) ) );
        }
                  
        geoGroup.__faceArray = new Uint16Array(faceArr);
     
        geoGroup.truncateArrayBuffers(true);
        
        if (normalArr.length < geoGroup.vertices)
            geoGroup.setNormals();
        else {
            
            geoGroup.__normalArray = new Float32Array(geoGroup.vertices*3);
            var n;
            for (var i = 0; i < geoGroup.vertices; ++i) {
                offset = i*3;
                n = normalArr[i];
                geoGroup.__normalArray[offset] = n.x, geoGroup.__normalArray[offset+1] = n.y, geoGroup.__normalArray[offset+2] = n.z;
            }
        }
            
        
        if (! lineArr.length)
            geoGroup.setLineIndices(); 
        else
            geoGroup.__lineArray = new Uint16Array(lineArr);
        
        geoGroup.lineidx = geoGroup.__lineArray.length;
             
    };
    
    //Read a cube file - generate model and possibly shape(s)
    var parseCube = function(shape, geoGroup, str, isoval, voxel) {
        
        var lines = str.replace(/^\s+/, "").split(/[\n\r]+/);
        
        if (lines.length < 6)
            return;
            
        var lineArr = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");       
          
        var natoms = Math.abs(parseFloat(lineArr[0]));        
        var origin = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3]));
        
        lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        
        //might have to convert from bohr units to angstroms
        var convFactor = (parseFloat(lineArr[0]) > 0) ? 0.529177 : 1;
        
        origin.multiplyScalar(convFactor);
        
        var nX = Math.abs(lineArr[0]);
        var xVec = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3])).multiplyScalar(convFactor);
        
        lineArr = lines[4].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        
        var nY = Math.abs(lineArr[0]);
        var yVec = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3])).multiplyScalar(convFactor);
        
        lineArr = lines[5].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
        
        var nZ = Math.abs(lineArr[0]);
        var zVec = new WebMol.Vector3(parseFloat(lineArr[1]), parseFloat(lineArr[2]), parseFloat(lineArr[3])).multiplyScalar(convFactor);
        
        //lines.splice(6, natoms).join("\n");
        
        lines = new Float32Array(lines.splice(natoms+7).join(" ").replace(/^\s+/, "").split(/[\s\r]+/));        
        
        var vertnums = new Int16Array(nX*nY*nZ);        
        for (var i = 0; i < vertnums.length; ++i)
            vertnums[i] = -1;

        var bitdata = new Uint8Array(nX*nY*nZ);
        
        for (var i = 0; i < lines.length; ++i) {
            var val = (isoval >= 0) ? lines[i] - isoval : isoval - lines[i];
            
            if (val > 0)
                bitdata[i] |= ISDONE;
            
        }
        
        var verts = [], faces = [];
        
        WebMol.MarchingCube.march(bitdata, verts, faces, {
            fulltable : true,
            voxel : voxel,
            scale : xVec.length(),
            origin : origin,
            nX : nX,
            nY : nY,
            nZ : nZ        
        });
        
        if (!voxel)
            WebMol.MarchingCube.laplacianSmooth(10, verts, faces);        
        
        drawCustom(shape, geoGroup, {vertexArr:verts, 
                                     faceArr:faces,
                                     normalArr:[],
                                     lineArr:[]});
        
    };
    
    //Update a bounding sphere's position and radius
    //from list of centroids and new points
    var updateBoundingFromPoints = function(sphere, components, points) {       
           
        sphere.center.set(0,0,0);
            
        if (components.length > 0) {
                            
            for (var i = 0, il = components.length; i < il; ++i) {
                var centroid = components[i].centroid;
                sphere.center.add(centroid);            
            }                
            
            sphere.center.divideScalar(components.length);           
        }
       
        var maxRadiusSq = sphere.radius*sphere.radius;
        
        for (var i = 0, il = points.length / 3; i < il; i++) {              
            var x = points[i*3], y = points[i*3 + 1], z = points[i*3 + 2];
            var radiusSq = sphere.center.distanceToSquared({x:x, y:y, z:z});                
            maxRadiusSq = Math.max(maxRadiusSq, radiusSq);
        }
        
        sphere.radius = Math.sqrt(maxRadiusSq);                        

    };
    
    var updateFromStyle = function(shape, stylespec) {
        shape.color = stylespec.color || new WebMol.Color();
        shape.wireframe = stylespec.wireframe ? true : false;
        shape.alpha = stylespec.alpha ? WebMol.Math.clamp(stylespec.alpha, 0.0, 1.0) : 1.0;
        shape.side = (stylespec.side !== undefined) ? stylespec.side : WebMol.DoubleSide;            

        //Click handling
        shape.clickable = stylespec.clickable ? true : false;
        shape.callback = typeof(stylespec.callback) === "function" ? stylespec.callback : null;
    };
    
    /**
     * Custom renderable shape
     * 
     * @constructor WebMol.GLShape
     * 
     * @param {Number} sid - Unique identifier
     * @param {Object} stylespec
     * @returns {WebMol.GLShape}
     */
    var GLShape = function(sid, stylespec) {
        
        stylespec = stylespec || {};
        WebMol.ShapeIDCount++;
        this.id = sid;
               
        this.boundingSphere = new WebMol.Sphere();
        this.intersectionShape = {sphere: [], cylinder: [], line: [], triangle: []};
        
        updateFromStyle(this, stylespec);
        
        //Keep track of shape components and their centroids
        var components = [];
        var shapeObj = null;
        var renderedShapeObj = null;
        
        var geo = new WebMol.Geometry(true);
                
        this.updateStyle = function(newspec) {
            
            for (var prop in newspec) {
                stylespec[prop] = newspec[prop];
            }    
            
            updateFromStyle(this, stylespec);
        };
        
        this.addCustom = function(customSpec) {
            
            customSpec.vertexArr = customSpec.vertexArr || [];
            customSpec.faceArr = customSpec.faceArr || [];
            customSpec.normalArr = customSpec.normalArr || [];
            customSpec.lineArr = customSpec.lineArr || [];
            
            //Force creation of new geometryGroup for each added component
            var geoGroup = geo.addGeoGroup();
            drawCustom(this, geoGroup, customSpec);            
            geoGroup.truncateArrayBuffers(true);
    
            for (var i = 0; i < geoGroup.__colorArray.length / 3; ++i) {
                geoGroup.__colorArray[i*3] = this.color.r;
                geoGroup.__colorArray[i*3 + 1] = this.color.g;
                geoGroup.__colorArray[i*3 + 2] = this.color.b;      
            }            
            
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup,
                centroid : geoGroup.getCentroid() 
            });
            
            updateBoundingFromPoints( this.boundingSphere, components, geoGroup.__vertexArray );            
        };
        
        //TODO: Refactor so 'drawSphere' method automatically updates bounding sphere as vertices are added
        this.addSphere = function(sphereSpec) {
          
            sphereSpec.center = sphereSpec.center || {x: 0, y: 0, z: 0};
            sphereSpec.radius = sphereSpec.radius ? WebMol.Math.clamp(sphereSpec.radius, 0, Infinity) : 1.5;
            
            var geoGroup = geo.addGeoGroup();
            drawSphere(this, geoGroup, sphereSpec);
            geoGroup.truncateArrayBuffers(true);            
            
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup, //has to be last group added
                centroid : new WebMol.Vector3(sphereSpec.center.x, sphereSpec.center.y, sphereSpec.center.z)
            });
            
            updateBoundingFromPoints( this.boundingSphere, components, geoGroup.__vertexArray );
        };        
        
        this.addArrow = function(arrowSpec) {
            
            arrowSpec.start = arrowSpec.start || {};
            arrowSpec.end = arrowSpec.end || {};
            
            arrowSpec.start = new WebMol.Vector3(arrowSpec.start.x || 0, arrowSpec.start.y || 0, arrowSpec.start.z || 0);
            arrowSpec.end = new WebMol.Vector3(arrowSpec.end.x || 3, arrowSpec.end.y || 0, arrowSpec.end.z || 0);
            arrowSpec.radius = arrowSpec.radius || 0.25;
            
            arrowSpec.radiusRatio = arrowSpec.radiusRatio || 1.618034;
            arrowSpec.mid = ( 0 < arrowSpec.mid && arrowSpec.mid < 1) ? arrowSpec.mid : 0.618034;
            
            var geoGroup = geo.addGeoGroup();
            
            drawArrow(this, geoGroup, arrowSpec);
            geoGroup.truncateArrayBuffers(true);
            
            var centroid = new WebMol.Vector3();
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup,
                centroid : centroid.addVectors(arrowSpec.start, arrowSpec.end).multiplyScalar(0.5)
            });
            
            updateBoundingFromPoints( this.boundingSphere, components, geoGroup.__vertexArray );
                            
        };
        
        this.addVolumetricData = function(data, fmt, volSpec) {
                            
            //str, fmt, isoval, vxl
            var isoval = (volSpec.isoval !== undefined && typeof(volSpec.isoval) === "number") ? volSpec.isoval : 0.0;
            var vxl = (volSpec.voxel) ? true : false;
            
            var geoGroup = geo.addGeoGroup();
            
            //TODO: Initialize geometry group here (parseCube currently calls addCustom...)
            switch(fmt) {
                case "cube":
                    parseCube(this, geoGroup, data, isoval, vxl);
                    break;
            }              
            
            components.push({
                id : geoGroup.id,
                geoGroup : geoGroup,
                centroid : geoGroup.getCentroid()    
            });
            
            this.updateStyle(volSpec);
            
            updateBoundingFromPoints( this.boundingSphere, components, geoGroup.__vertexArray );
            
        };
    
        this.globj = function(group) {
            
            geo.initTypedArrays();
            
            updateColor(geo, this.color);
            
            shapeObj = new WebMol.Object3D();
            var material = new WebMol.MeshLambertMaterial({
                wireframe : this.wireframe,
                vertexColors : true,
                ambient : 0x000000,
                reflectivity : 0,
                side : this.side,
                transparent : (this.alpha < 1) ? true : false,
                opacity : this.alpha
            });
            
            var mesh = new WebMol.Mesh(geo, material);
            shapeObj.add(mesh);
                 
            if (renderedShapeObj) {
                group.remove(renderedShapeObj);
                renderedShapeObj = null;
            }
            renderedShapeObj = shapeObj.clone();
            group.add(renderedShapeObj);
            
        };
    
    };

    Object.defineProperty(GLShape.prototype, "position", {
        
        get : function() {
            return this.boundingSphere.center;
        } 
            
    });    
    
    Object.defineProperty(GLShape.prototype, "x", {
        
        get : function() {
            return this.boundingSphere.center.x;
        } 
            
    });  
 
    Object.defineProperty(GLShape.prototype, "y", {
        
        get : function() {
            return this.boundingSphere.center.y;
        } 
            
    });   

    Object.defineProperty(GLShape.prototype, "z", {
        
        get : function() {
            return this.boundingSphere.center.z;
        } 
            
    });  
  
    return GLShape;
    
}());

WebMol.ShapeIDCount = 0;//a molecular viewer based on GLMol

var WebMol = WebMol || {};

//Adapted from the text sprite example from http://stemkoski.github.io/Three.js/index.html

// function for drawing rounded rectangles
var roundRect = function(ctx, x, y, w, h, r) {

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
 
};

WebMol.LabelCount = 0;

/**
 * Renderable labels
 * @constructor WebMol.Label
 * @param {string} tag Label text
 * @param {Object} parameters Label style and font specifications
 */
WebMol.Label = function(text, parameters) {
        
    this.id = WebMol.LabelCount++;    
    this.stylespec = parameters || {};

    this.canvas = document.createElement('canvas');
    
    this.context = this.canvas.getContext('2d');

    this.sprite = new WebMol.Sprite();
    this.text = text;
    
};


WebMol.Label.prototype = {
    
    constructor : WebMol.Label,
    
    setContext : function() {
        
        //Update stylespec
        /** Label text font style
         * @public
         * @type {string} */
        this.font = this.stylespec.font = 
            this.stylespec.font ? this.stylespec.font : "Arial";
    
         
        /** Label text font pt size
         * @type {number} */
        this.fontSize = this.stylespec.fontSize =
            this.stylespec.fontSize ? this.stylespec.fontSize : 20;
            
        /** Label font color - specify with an object with r, g, b, and a (alpha) values
         * @type {Object | WebMol.Color} */
        this.fontColor = this.stylespec.fontColor =
            this.stylespec.fontColor ? this.stylespec.fontColor : { r:255, g:255, b:255, a:1.0};
    
        this.borderThickness = this.stylespec.borderThickness =
            this.stylespec.borderThickness ? this.stylespec.borderThickness : 4;
    
        this.borderColor = this.stylespec.borderColor =
            this.stylespec.borderColor ? this.stylespec.borderColor : { r:0, g:0, b:0, a:1.0 };
    
        this.backgroundColor = this.stylespec.backgroundColor =
            this.stylespec.backgroundColor ? this.stylespec.backgroundColor : { r:0, g:0, b:0, a:1.0 };
            
        this.position = this.stylespec.position =
            this.stylespec.position ? this.stylespec.position : { x:-10, y:1, z:1 };
        
        //Should labels always be in front of model? 
        this.inFront = this.stylespec.inFront = 
            (this.stylespec.inFront !== undefined) ? this.stylespec.inFront : true;
            
        //clear canvas
        this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);
            
        var spriteAlignment = WebMol.SpriteAlignment.topLeft;
        
        this.context.font = this.fontSize + "pt " + this.font;
        
        var metrics = this.context.measureText(this.text);           
        var textWidth = metrics.width;
        
        // background color
        this.context.fillStyle   = "rgba(" + this.backgroundColor.r + "," + this.backgroundColor.g + ","
                                                                  + this.backgroundColor.b + "," + this.backgroundColor.a + ")";
        // border color
        this.context.strokeStyle = "rgba(" + this.borderColor.r + "," + this.borderColor.g + ","
                                                                  + this.borderColor.b + "," + this.borderColor.a + ")";
    
        this.context.lineWidth = this.borderThickness;
        roundRect(this.context, this.borderThickness/2, this.borderThickness/2, textWidth + this.borderThickness, this.fontSize * 1.4 + this.borderThickness, 6);
        // 1.4 is extra height factor for text below baseline: g,j,p,q.
    
        // text color
        this.context.fillStyle = "rgba(" + this.fontColor.r + "," + this.fontColor.g + ","
                                                                + this.fontColor.b + "," + this.fontColor.a + ")";
    
        this.context.fillText(this.text, this.borderThickness, this.fontSize + this.borderThickness, textWidth);
        
        // canvas contents will be used for a texture
        var texture = new WebMol.Texture(this.canvas);
        texture.needsUpdate = true;
        
        this.sprite.material = new WebMol.SpriteMaterial( 
                { map: texture, useScreenCoordinates: false, alignment: spriteAlignment, depthTest: !this.inFront } );
                    
        this.sprite.scale.set(2 * this.fontSize, this.fontSize, 1);
        this.sprite.position.set(this.position.x, this.position.y, this.position.z);
        
    },
    
    //clean up material and texture
    dispose : function() {
        
        if (this.sprite.material.map !== undefined)
            this.sprite.material.map.dispose();
        if (this.sprite.material !== undefined)
            this.sprite.material.dispose();        
    }
    
};

// a webmol unified interace to gmol
WebMol.GLViewer = (function() {
    // private class variables
    var numWorkers = 4; // number of threads for surface generation
    var maxVolume = 64000; // how much to break up surface calculations

    // private class helper functions

    // computes the bounding box around the provided atoms
    var getExtent = function(atomlist) {
        var xmin, ymin, zmin,
            xmax, ymax, zmax,
            xsum, ysum, zsum, cnt;
        
        xmin = ymin = zmin = 9999;
        xmax = ymax = zmax = -9999;
        xsum = ysum = zsum = cnt = 0;

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
    /**
     * WebGL WebMol viewer
     * Note: The preferred method of instantiating a GLViewer is through {@link WebMol.createViewer} 
     * 
     * @constructor WebMol.GLViewer
     * @param {Object} element HTML element within which to create viewer
     * @param {Function} callback - Callback function to be immediately executed on this viewer
     * @param {Object} defaultcolors - Object defining default atom colors as atom => color property value pairs for all models within this viewer
     */
    function GLViewer(element, callback, defaultcolors) {

        // set variables
        var _viewer = this;
        var container = element;
        var id = container.id;

        var models = []; // atomistic molecular models
        var surfaces = [];
        var shapes = []; //Generic shapes

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
        camera.position = new WebMol.Vector3(0, 0, CAMERA_Z);
        camera.lookAt(new WebMol.Vector3(0, 0, 0));
        
        var raycaster = new WebMol.Raycaster(new WebMol.Vector3(0,0,0), new WebMol.Vector3(0,0,0));
        var projector = new WebMol.Projector();
        var mouseVector = new WebMol.Vector3(0, 0, 0);

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
            
            //var time = new Date();
            setSlabAndFog();
            renderer.render(scene, camera);
            //console.log("rendered in " + (+new Date() - time) + "ms");
        };

        var initializeScene = function() {
            
            scene = new WebMol.Scene();
            scene.fog = new WebMol.Fog(bgColor, 100, 200);

            modelGroup = new WebMol.Object3D();
            rotationGroup = new WebMol.Object3D();
            rotationGroup.useQuaternion = true;
            rotationGroup.quaternion = new WebMol.Quaternion(0, 0, 0, 1);
            rotationGroup.add(modelGroup);

            scene.add(rotationGroup);

            // setup lights
            var directionalLight = new WebMol.Light(0xFFFFFF);
            directionalLight.position = new WebMol.Vector3(0.2, 0.2, 1).normalize();
            directionalLight.intensity = 1.0;
            scene.add(directionalLight);
        };

        initializeScene();
        
        renderer.setClearColorHex(bgColor, 1.0);
        scene.fog.color = WebMol.CC.color(bgColor);
        
        var clickedAtom = null;
        // enable mouse support
        var glDOM = $(renderer.domElement);

        //Checks for selection intersects on mousedown
        var handleClickSelection = function(mouseX, mouseY) {
            
            var mouse = {x : mouseX, y : mouseY, z : -1.0};
            mouseVector.set(mouse.x, mouse.y, mouse.z);
            projector.unprojectVector(mouseVector, camera);
            mouseVector.sub(camera.position).normalize();
             
            raycaster.set(camera.position, mouseVector);

            var clickables = [], intersects = [];
            
            for (var i in models) {
                var model = models[i];
                
                var atoms = model.selectedAtoms({clickable: true});
                clickables = clickables.concat(atoms);

            }
            
            for (var i in shapes) {
                
                var shape = shapes[i];
                if (shape.clickable) {
                    clickables.push(shape);
                }

            }
            
            intersects = raycaster.intersectObjects(modelGroup, clickables);
            
            if (intersects.length) {
                var selected = intersects[0].clickable;
                if (selected.callback !== undefined && typeof(selected.callback) === "function") {
                    selected.callback(selected, _viewer);
                }
            }
            
            show();        
        }; 
        
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
            clickedAtom = null;
            mouseButton = ev.which;
            mouseStartX = x;
            mouseStartY = y;
            cq = rotationGroup.quaternion;
            cz = rotationGroup.position.z;
            currentModelPos = modelGroup.position.clone();
            cslabNear = slabNear;
            cslabFar = slabFar;
            
            //handle selection
            var mouseX = (x / $(window).width())*2 - 1;
            var mouseY = -(y / HEIGHT)*2 + 1;
            handleClickSelection(mouseX, mouseY, ev, container);
            
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
                var translationByScreen = new WebMol.Vector3(dx * scaleFactor, -dy
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
        /**
         * Set the background color (default white)
         * 
         * @function WebMol.GLViewer#setBackgroundColor
         * @param {number} hex Hexcode specified background color
         * @param {number} a Alpha level (default 1.0)
         * 
         * @example
         * 
         * //Set 'myviewer' background color to white
         * myviewer.setBackgroundColor(0xffffff)
         * 
         */
        this.setBackgroundColor = function(hex, a) {
            a = a | 1.0;
            bgColor = hex;
            renderer.setClearColorHex(hex, a);
            scene.fog.color = WebMol.CC.color(hex);
            show();
        };
        
        /**
         * Set viewer width
         * 
         * @function WebMol.GLViewer#setWidth
         * @param {number} w Width in pixels
         */
        this.setWidth = function(w) {
            WIDTH = w || WIDTH;
            renderer.setSize(WIDTH, HEIGHT);
        };
        
        /**
         * Set viewer height
         * 
         * @function WebMol.GLViewer#setHeight
         * @param {number} h Height in pixels
         */
        this.setHeight = function(h) {
            HEIGHT = h || HEIGHT;
            renderer.setSize(WIDTH, HEIGHT);
        };
        
        /**
         * Resize viewer according to containing HTML element's dimensions
         * 
         * @function WebMol.GLViewer#resize
         */
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

        /**
         * Return specified model
         * 
         * @function WebMol.GLViewer#getModel
         * @param {number} [id=last model id] - Retrieve model with specified id
         * @default Returns last model added to viewer
         * @returns {GLModel}
         * 
         * @example
         * // Retrieve reference to first GLModel added
         * var m = glviewer.getModel(0);
         */
        this.getModel = function(id) {
            id = id || models.length - 1;
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
            
            if (arg === undefined || !(arg instanceof Array || arg.length !== 8))
                return;                       
                
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
        /**
         * Render current state of viewer, after 
         * adding/removing models, applying styles, etc.
         * 
         * @function WebMol.GLViewer#render
         */
        this.render = function() {

            //spinner.show();
            var time1 = new Date();
            var view = this.getView();

            for ( var i = 0; i < models.length; i++) {
                if (models[i]) {
                    models[i].globj(modelGroup);
                }
            }
            
            for ( var i = 0; i < shapes.length; i++ ) {
                if (shapes[i]) {
                    shapes[i].globj(modelGroup);
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

        /**
         * Return pdb output of selected atoms (if atoms from pdb input)
         * 
         * @function WebMol.GLViewer#pdbData  
         * @param {Object} [sel] - Selection specification specifying model and atom properties to select.  Default: all atoms in viewer
         * @returns {string} PDB string of selected atoms
         */
        this.pdbData = function(sel) {
            var atoms = getAtomsFromSel(sel);
            var ret = "";
            for ( var i = 0, n = atoms.length; i < n; ++i) {
                ret += atoms[i].pdbline + "\n";
            }
            return ret;
        };

        /**
         * Zoom to center of atom selection
         * 
         * @function WebMol.GLViewer#zoomTo
         * @param {Object} [sel] - Selection specification specifying model and atom properties to select. Default: all atoms in viewer
         * 
         * @example
         * // Assuming we have created a model of a protein with multiple chains (e.g. from a PDB file), focus on atoms in chain B
         * glviewer.zoomTo({chain: 'B'});
         * 
         * // Focus on centroid of all atoms of all models in this viewer
         * glviewer.zoomTo();  // (equivalent to glviewer.zoomTo({}) )
         */
        this.zoomTo = function(sel) {
            var atoms = getAtomsFromSel(sel).concat(shapes);
            var allatoms = getAtomsFromSel({}).concat(shapes);
            var tmp = getExtent(atoms);
            var alltmp = getExtent(allatoms);
            // use selection for center
            var center = new WebMol.Vector3(tmp[2][0], tmp[2][1], tmp[2][2]);
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
        
        /**
         * Add label to viewer
         * 
         * @function WebMol.GLViewer#addLabel
         * @param {string} text - Label text
         * @param {Object} data - Label style specification
         * @returns {WebMol.Label}
         * 
         * @example
         * 
         * // Assuming glviewer contains a model representing a protein, label all alpha carbons with their residue name
         * 
         * // Select all alpha carbons (have property atom : "CA") from last model added
         * var atoms = glviewer.getModel().selectedAtoms({atom:"CA"});
         * var labels = [];
         * 
         * for (var a in atoms) {
         *     var atom = atoms[a];
         * 
         *     // Create label at alpha carbon's position displaying atom's residue and residue number
         *     var labelText = atom.resname + " " + atom.resi;
         *      
         *     var l = glviewer.createLabel(labelText, {fontSize: 12, position: {x: atom.x, y: atom.y, z: atom.z});
         * 
         *     labels.push(l);
         * }
         * 
         * // Render labels
         * glviewer.render();
         */
        this.addLabel = function(text, data) {
            var label = new WebMol.Label(text, data); 
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;
        };
        
        /**
         * Remove label from viewer
         * 
         * @function WebMol.GLViewer#removeLabel
         * @param {WebMol.Label} label - WebMol label
         * 
         * @example
         * // Remove labels created in [addLabel example]{@link WebMol.GLViewer#addLabel}
         * 
         * for (var i = 0; i < labels.length; i++) {
         *     glviewer.removeLabel(label);
         * }
         * 
         * glviewer.render();
         */
        this.removeLabel = function(label) {

            label.dispose();
            modelGroup.remove(label.sprite);                       
        };
        
        //Modify label style
        /**
         * Modify existing label's style
         * 
         * @function WebMol.GLViewer#setLabelStyle
         * @param {WebMol.Label} label - WebMol label
         * @param {Object} stylespec - Label style specification
         * @returns {WebMol.Label}
         */
        this.setLabelStyle = function(label, stylespec) {   
             
            label.dispose();
            label.stylespec = stylespec;
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;
            
        };
        
        //Change label text
        /**
         * Modify existing label's text
         * 
         * @function WebMol.GLViewer#setLabelText
         * @param {WebMol.Label} label - WebMol label
         * @param {String} text - Label text
         * @returns {WebMol.Label}
         */
        this.setLabelText = function(label, text) {
         
            label.dispose();
            label.text = text;
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;

        };
        
        /**
         * Add shape object to viewer 
         * @see {@link WebMol.GLShape}
         * 
         * @function WebMol.GLViewer#addShape
         * @param {Object} shapeSpec - style specification for label
         * @returns {WebMol.GLShape}
         */
        this.addShape = function(shapeSpec) {
            shapeSpec = shapeSpec || {};
            var shape = new WebMol.GLShape(shapes.length, shapeSpec);
            shapes.push(shape);
            
            return shape;
              
        };
        
        /**
         * Create and add sphere shape. This method provides a shorthand 
         * way to create a spherical shape object
         * 
         * @param {Object} spec - Sphere shape style specification
         * @returns {WebMol.GLShape}
         */
        this.addSphere = function(spec) {
            var s = new WebMol.GLShape(shapes.length);
            spec = spec || {};
            s.addSphere(spec);      
            shapes.push(s);
            
            return s;
        };
        
        /**
         * Create and add arrow shape
         * 
         * @param {Object} spec - Style specification
         * @returns {WebMol.GLShape}
         */
        this.addArrow = function(spec) {            
            var s = new WebMol.GLShape(shapes.length);            
            spec = spec || {};
            s.addArrow(spec);
            shapes.push(s);
            
            return s;
        };
        
        /**
         * Add custom shape component from user supplied function
         * 
         * @param {Object} spec - Style specification
         * @returns {WebMol.GLShape}
         */
        this.addCustom = function(spec) {   
            var s = new WebMol.GLShape(shapes.length);                         
            spec = spec || {};
            s.addCustom(spec);     
            shapes.push(s);
            
            return s;                       
        };
        
        /**
         * Construct isosurface from volumetric data in gaussian cube format
         * 
         * @param {String} data - Input file contents 
         * @param {String} format - Input file format (currently only supports "cube")
         * @param {Object} spec - Shape style specification
         * @returns {WebMol.GLShape}
         */
        this.addVolumetricData = function(data, format, spec) {
            var s = new WebMol.GLShape(shapes.length);
            spec = spec || {};            
            s.addVolumetricData(data, format, spec);   
            shapes.push(s);
            
            return s;       
        };

        /**
         * Create and add model to viewer, given molecular data and its format 
         * (pdb, sdf, xyz, or mol2)
         * 
         * @param {String} data - Input data
         * @param {String} format - Input format ('pdb', 'sdf', 'xyz', or 'mol2')
         * @returns {WebMol.GLModel}
         */
        this.addModel = function(data, format) {
           
            var m = new WebMol.GLModel(models.length, defaultcolors);
            m.addMolData(data, format);
            models.push(m);
            
            return m;
        };

        /**
         * Delete specified model from viewer
         * 
         * @param {WebMol.GLModel} model
         */
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

        /** 
         * Delete all existing models
         */
        this.removeAllModels = function() {
            for (var i = 0; i < models.length; i++) {
                var model = models[i];
                model.removegl(modelGroup);
                
            }
            models = [];
        };

        /**
         * Create a new model from atoms specified by sel.
         * If extract, removes selected atoms from existing models 
         * @param {Object} sel - Atom selection specification
         * @param {Boolean} extract - If true, remove selected atoms from existing models
         * @returns {WebMol.GLModel}
         */
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

        /**
         * Set style properties to all selected atoms
         * 
         * @param {Object} sel - Atom selection specification
         * @param {Object} style - Style spec to apply to specified atoms
         */
        this.setStyle = function(sel, style) {
            applyToModels("setStyle", sel, style, false);
        };
        
        /**
         * Add style properties to all selected atoms
         * 
         * @param {Object} sel - Atom selection specification
         * @param {Object} style - style spec to add to specified atoms
         */
        this.addStyle = function(sel, style) {
            applyToModels("setStyle", sel, style, true);
        };

        /**
         * 
         * @param {type} sel
         * @param {type} prop
         * @param {type} scheme
         * @returns {undefined}
         */
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
                    
                    var index;
                    
                    if (w > h && w > d) {
                        index = 0;
                    } 
                    else if (h > w && h > d) {
                        index = 1;
                    } 
                    else {
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
        
            var geo = new WebMol.Geometry(true);  
            //Only one group per call to generate surface mesh (addSurface should split up mesh render)     
            var geoGroup = geo.updateGeoGroup(0);
            
            // reconstruct vertices and faces
            var v = VandF.vertices;
            var offset;
            for ( var i = 0; i < v.length; i++) {            
                offset = geoGroup.vertices*3;
                geoGroup.__vertexArray[offset] = v[i].x; geoGroup.__vertexArray[offset+1] = v[i].y, geoGroup.__vertexArray[offset+2] =v[i].z;                
                geoGroup.vertices++;
            }
                       
            var faces = VandF.faces;
            geoGroup.faceidx = faces.length;//*3;
            geo.initTypedArrays();

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
            
            var verts = geoGroup.__vertexArray;
            var vA, vB, vC, norm;
            var faceoffset;
            
            //Setup colors, faces, and normals
            for ( var i = 0; i < faces.length; i+=3) {
                
                faceoffset = i;
                //var a = faces[i].a, b = faces[i].b, c = faces[i].c;
                var a = faces[i], b = faces[i+1], c = faces[i+2];
                var A = v[a].atomid;
                var B = v[b].atomid;
                var C = v[c].atomid;
                
                var offsetA = a * 3, offsetB = b * 3, offsetC = c * 3;

                geoGroup.__faceArray[faceoffset] = faces[i].a, geoGroup.__faceArray[faceoffset+1] = faces[i].b, 
                    geoGroup.__faceArray[faceoffset+2] = faces[i].c;
                
                geoGroup.__colorArray[offsetA] = colors[A].r, geoGroup.__colorArray[offsetA+1] = colors[A].g,
                         geoGroup.__colorArray[offsetA+2] = colors[A].b;
                geoGroup.__colorArray[offsetB] = colors[B].r, geoGroup.__colorArray[offsetB+1] = colors[B].g,
                         geoGroup.__colorArray[offsetB+2] = colors[B].b;
                geoGroup.__colorArray[offsetC] = colors[C].r, geoGroup.__colorArray[offsetC+1] = colors[C].g,
                         geoGroup.__colorArray[offsetC+2] = colors[C].b;
                 
                //setup Normals
                
                vA = new WebMol.Vector3(verts[offsetA], verts[offsetA+1], verts[offsetA+2]);
                vB = new WebMol.Vector3(verts[offsetB], verts[offsetB+1], verts[offsetB+2]);
                vC = new WebMol.Vector3(verts[offsetC], verts[offsetC+1], verts[offsetC+2]);
                
                vC.subVectors(vC, vB);
                vA.subVectors(vA, vB);
                vC.cross(vA);

                //face normal
                norm = vC;
                norm.normalize();
                
                geoGroup.__normalArray[offsetA] += norm.x, geoGroup.__normalArray[offsetB] += norm.x, geoGroup.__normalArray[offsetC] += norm.x;
                geoGroup.__normalArray[offsetA+1] += norm.y, geoGroup.__normalArray[offsetB+1] += norm.y, geoGroup.__normalArray[offsetC+1] += norm.y;
                geoGroup.__normalArray[offsetA+2] += norm.z, geoGroup.__normalArray[offsetB+2] += norm.z, geoGroup.__normalArray[offsetC+2] += norm.z;
                
            }
            geoGroup.__faceArray = new Uint16Array(faces);
            var mesh = new WebMol.Mesh(geo, mat);
            mesh.doubleSided = true;

            return mesh;
        };

        // do same thing as worker in main thread
        var generateMeshSyncHelper = function(type, expandedExtent,
                extendedAtoms, atomsToShow, atoms, vol) {
            var time = new Date();
            var ps = new WebMol.ProteinSurface();
            ps.initparm(expandedExtent, (type === 1) ? false : true, vol);

            var time2 = new Date();
            console.log("initialize " + (time2 - time) + "ms");

            ps.fillvoxels(atoms, extendedAtoms);

            var time3 = new Date();
            console.log("fillvoxels " + (time3 - time2) + "  " + (time3 - time)
                    + "ms");

            ps.buildboundary();


            if (type == WebMol.SurfaceType.SES) {
                ps.fastdistancemap();
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
                geo : new WebMol.Geometry(true),
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

            var sync = !!(WebMol.syncSurface);
            var view = this; //export render function to worker
            if (sync) { // don't use worker, still break up for memory purposes

                for ( var i = 0; i < extents.length; i++) {
                    //console.profile();
                    var VandF = generateMeshSyncHelper(type, extents[i].extent,
                            extents[i].atoms, extents[i].toshow, reducedAtoms,
                            totalVol);
                    var mesh = generateSurfaceMesh(atomlist, VandF, mat);
                    WebMol.mergeGeos(surfobj.geo, mesh);
                    view.render();
                    //console.profileEnd();
                }
            //TODO: Asynchronously generate geometryGroups (not separate meshes) and merge them into a single geometry
            } else { // use worker
                
                var workers = [];
                if (type < 0)
                    type = 0; // negative reserved for atom data
                for ( var i = 0; i < numWorkers; i++) {
                    //var w = new Worker('webmol/SurfaceWorker.js');
                    var w = new Worker(WebMol.SurfaceWorker);
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
                        WebMol.mergeGeos(surfobj.geo, mesh);
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
                surfaces[surf].mat.side = WebMol.FrontSide;
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
        
        getModelGroup = function() {
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

WebMol.glmolViewer = WebMol.GLViewer;
var WebMol = WebMol || {};

//properties for mapping
WebMol.partialCharges = [
{ resn: "ALA", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "ALA", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "ALA", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "ALA", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "ALA", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "ARG", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "ARG", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "ARG", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "ARG", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "ARG", atom: "CD", props: {partialCharge: 0.10}},
	{ resn: "ARG", atom: "NE", props: {partialCharge: -0.10}},
	{ resn: "ARG", atom: "CZ", props: {partialCharge: 0.50}},
	{ resn: "ARG", atom: "NH1", props: {partialCharge: 0.25}},
	{ resn: "ARG", atom: "NH2", props: {partialCharge: 0.25}},
	{ resn: "ARG", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "ARG", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "ASN", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "ASN", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "ASN", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "ASN", atom: "CG", props: {partialCharge: 0.55}},
	{ resn: "ASN", atom: "OD1", props: {partialCharge: -0.55}},
	{ resn: "ASN", atom: "ND2", props: {partialCharge: 0.00}},
	{ resn: "ASN", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "ASN", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "ASP", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "ASP", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "ASP", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "ASP", atom: "CG", props: {partialCharge: 0.14}},
	{ resn: "ASP", atom: "OD1", props: {partialCharge: -0.57}},
	{ resn: "ASP", atom: "OD2", props: {partialCharge: -0.57}},
	{ resn: "ASP", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "ASP", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "CYS", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "CYS", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "CYS", atom: "CB", props: {partialCharge: 0.19}},
	{ resn: "CYS", atom: "SG", props: {partialCharge: -0.19}},
	{ resn: "CYS", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "CYS", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "GLN", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "GLN", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "GLN", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "GLN", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "GLN", atom: "CD", props: {partialCharge: 0.55}},
	{ resn: "GLN", atom: "OE1", props: {partialCharge: -0.55}},
	{ resn: "GLN", atom: "NE2", props: {partialCharge: 0.00}},
	{ resn: "GLN", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "GLN", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "GLU", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "GLU", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "GLU", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "GLU", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "GLU", atom: "CD", props: {partialCharge: 0.14}},
	{ resn: "GLU", atom: "OE1", props: {partialCharge: -0.57}},
	{ resn: "GLU", atom: "OE2", props: {partialCharge: -0.57}},
	{ resn: "GLU", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "GLU", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "GLY", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "GLY", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "GLY", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "GLY", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "HIS", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "HIS", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "HIS", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "HIS", atom: "CG", props: {partialCharge: 0.10}},
	{ resn: "HIS", atom: "ND1", props: {partialCharge: -0.10}},
	{ resn: "HIS", atom: "CD2", props: {partialCharge: 0.10}},
	{ resn: "HIS", atom: "NE2", props: {partialCharge: -0.40}},
	{ resn: "HIS", atom: "CE1", props: {partialCharge: 0.30}},
	{ resn: "HIS", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "HIS", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "ILE", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "ILE", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "ILE", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "ILE", atom: "CG2", props: {partialCharge: 0.00}},
	{ resn: "ILE", atom: "CG1", props: {partialCharge: 0.00}},
	{ resn: "ILE", atom: "CD", props: {partialCharge: 0.00}},
	{ resn: "ILE", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "ILE", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "LEU", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "LEU", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "LEU", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "LEU", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "LEU", atom: "CD1", props: {partialCharge: 0.00}},
	{ resn: "LEU", atom: "CD2", props: {partialCharge: 0.00}},
	{ resn: "LEU", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "LEU", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "LYS", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "LYS", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "LYS", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "LYS", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "LYS", atom: "CD", props: {partialCharge: 0.00}},
	{ resn: "LYS", atom: "CE", props: {partialCharge: 0.25}},
	{ resn: "LYS", atom: "NZ", props: {partialCharge: 0.75}},
	{ resn: "LYS", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "LYS", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "MET", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "MET", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "MET", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "MET", atom: "CG", props: {partialCharge: 0.06}},
	{ resn: "MET", atom: "SD", props: {partialCharge: -0.12}},
	{ resn: "MET", atom: "CE", props: {partialCharge: 0.06}},
	{ resn: "MET", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "MET", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "PHE", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "PHE", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "PHE", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "PHE", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "PHE", atom: "CD1", props: {partialCharge: 0.00}},
	{ resn: "PHE", atom: "CD2", props: {partialCharge: 0.00}},
	{ resn: "PHE", atom: "CE1", props: {partialCharge: 0.00}},
	{ resn: "PHE", atom: "CE2", props: {partialCharge: 0.00}},
	{ resn: "PHE", atom: "CZ", props: {partialCharge: 0.00}},
	{ resn: "PHE", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "PHE", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "PRO", atom: "N", props: {partialCharge: -0.25}},
	{ resn: "PRO", atom: "CD", props: {partialCharge: 0.10}},
	{ resn: "PRO", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "PRO", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "PRO", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "PRO", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "PRO", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "SER", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "SER", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "SER", atom: "CB", props: {partialCharge: 0.25}},
	{ resn: "SER", atom: "OG", props: {partialCharge: -0.25}},
	{ resn: "SER", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "SER", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "THR", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "THR", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "THR", atom: "CB", props: {partialCharge: 0.25}},
	{ resn: "THR", atom: "OG1", props: {partialCharge: -0.25}},
	{ resn: "THR", atom: "CG2", props: {partialCharge: 0.00}},
	{ resn: "THR", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "THR", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "TRP", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "TRP", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "TRP", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "TRP", atom: "CG", props: {partialCharge: -0.03}},
	{ resn: "TRP", atom: "CD2", props: {partialCharge: 0.10}},
	{ resn: "TRP", atom: "CE2", props: {partialCharge: -0.04}},
	{ resn: "TRP", atom: "CE3", props: {partialCharge: -0.03}},
	{ resn: "TRP", atom: "CD1", props: {partialCharge: 0.06}},
	{ resn: "TRP", atom: "NE1", props: {partialCharge: -0.06}},
	{ resn: "TRP", atom: "CZ2", props: {partialCharge: 0.00}},
	{ resn: "TRP", atom: "CZ3", props: {partialCharge: 0.00}},
	{ resn: "TRP", atom: "CH2", props: {partialCharge: 0.00}},
	{ resn: "TRP", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "TRP", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "TYR", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "TYR", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "TYR", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "TYR", atom: "CG", props: {partialCharge: 0.00}},
	{ resn: "TYR", atom: "CD1", props: {partialCharge: 0.00}},
	{ resn: "TYR", atom: "CE1", props: {partialCharge: 0.00}},
	{ resn: "TYR", atom: "CD2", props: {partialCharge: 0.00}},
	{ resn: "TYR", atom: "CE2", props: {partialCharge: 0.00}},
	{ resn: "TYR", atom: "CZ", props: {partialCharge: 0.25}},
	{ resn: "TYR", atom: "OH", props: {partialCharge: -0.25}},
	{ resn: "TYR", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "TYR", atom: "O", props: {partialCharge: -0.55}},
	{ resn: "VAL", atom: "N", props: {partialCharge: -0.15}},
	{ resn: "VAL", atom: "CA", props: {partialCharge: 0.10}},
	{ resn: "VAL", atom: "CB", props: {partialCharge: 0.00}},
	{ resn: "VAL", atom: "CG1", props: {partialCharge: 0.00}},
	{ resn: "VAL", atom: "CG2", props: {partialCharge: 0.00}},
	{ resn: "VAL", atom: "C", props: {partialCharge: 0.60}},
	{ resn: "VAL", atom: "O", props: {partialCharge: -0.55}}
]; 
	
