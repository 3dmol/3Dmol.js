//a molecular viewer based on GLMol

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

WebMol.Label = function(message, parameters) {
        
    this.id = WebMol.LabelCount++;    
    this.stylespec = parameters || {};

    this.canvas = document.createElement('canvas');
    
    this.context = this.canvas.getContext('2d');

    this.sprite = new WebMol.Sprite();
    this.text = message;
    
};

//TODO: How should this class be structured?
WebMol.Label.prototype = {
    
    constructor : WebMol.Label,
    
    setContext : function() {
        
        //Update stylespec
        this.font = this.stylespec.font = 
            this.stylespec.font ? this.stylespec.font : "Arial";
    
        this.fontSize = this.stylespec.fontSize =
            this.stylespec.fontSize ? this.stylespec.fontSize : 20;
            
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
    
    var zSort = function(a, b) {
        return a.z < b.z;
    };
    
    //Read a cube file - generate model and possibly shape(s)
    var parseCube = function(str, viewer) {
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
        
        //Extract atom portion; send to new GLModel...
        var atomStr = lines.splice(6, natoms).join("\n");
        atomStr = convFactor + "\n" + atomStr;
        
        lines = lines.splice(7).join(" ").replace(/^\s+/, "").split(/[\s\r]+/);
        lines = new Float32Array(lines);
        
        var isoval = 0.005;
        
        var cubepts = [
            origin.clone(), origin.clone().add(zVec),
            origin.clone().add(yVec), origin.clone().add(yVec).add(zVec),
            
            origin.clone().add(xVec), origin.clone().add(xVec).add(zVec),
            origin.clone().add(xVec).add(yVec), origin.clone().add(xVec).add(yVec).add(zVec)
        ];        
        
        //voxel values for current position
        var grid = new Float32Array(8);
                
        
        var p1 = new WebMol.Vector3(), p2 = new WebMol.Vector3();

        var vertnums = new Int16Array(nX*nY*nZ*12);
        
        for (var i = 0; i < vertnums.length; ++i)
            vertnums[i] = -1;

        //TODO: Need a good way to compute hard vertex normals for non-smoothed voxel (to get facetted look)
        smooth = true;
        
        //Also TODO:  vertnums must be signed (to initialize at -1) -> but this means we can't have more than
        // 32,768 vertices per geoGroup (rather than 65,536) - should probably enforce (or else use Int32Array for vertnums...)
        for (var neg = 0; neg < 2; ++neg) {
            
            if (neg)
                isoval = -isoval;
            
            var verts = [], faces = [], norms = [];
            
            for (var i = 0; i < nX - 1; ++i) {
                for (var j = 0; j < nY - 1; ++j) {
                    for (var k = 0; k < nZ - 1; ++k) {
                        
                        
                        //unpack voxels for this cube
                        
                        offset = (i*nY*nZ) + (j*nZ) + k;
                       
                        for (var p = 0; p < 8; p++) {
                            var index = ((nY * (i + ((p & 4) >> 2))) + j + ((p & 2) >> 1))
                                            * nZ + k + (p & 1);
                            grid[p] = parseFloat(lines[index]) * convFactor;
                        }
                        
                        var bit = 0;
                        
                        if ((grid[0] > isoval && isoval >= 0) || (grid[0] < isoval && isoval < 0))
                            bit |= 1;
                        
                        if ((grid[1] > isoval && isoval >= 0) || (grid[1] < isoval && isoval < 0))
                            bit |= 2;
                            
                        if ((grid[2] > isoval && isoval >= 0) || (grid[2] < isoval && isoval < 0))
                            bit |= 4;
                        
                        if ((grid[3] > isoval && isoval >= 0) || (grid[3] < isoval && isoval < 0))
                            bit |= 8;
                        
                        if ((grid[4] > isoval && isoval >= 0) || (grid[4] < isoval && isoval < 0))
                            bit |= 16;
                        
                        if ((grid[5] > isoval && isoval >= 0) || (grid[5] < isoval && isoval < 0))
                            bit |= 32;
                        
                        if ((grid[6] > isoval && isoval >= 0) || (grid[6] < isoval && isoval < 0))
                            bit |= 64;
                            
                        if ((grid[7] > isoval && isoval >= 0) || (grid[7] < isoval && isoval < 0))
                            bit |= 128;
                        
                        if (bit == 0 || bit == 255) 
                            continue;
                            
                        var edgeIdx = MarchingCube.edgeTable2[bit];
                        var triangles = MarchingCube.triTable2[bit];
                        
                        //Not on isosurface
                        if (edgeIdx == 0)
                            continue;
                            
                        //check edges
                        
                        var xV = xVec.clone().multiplyScalar(i);
                        var yV = yVec.clone().multiplyScalar(j);
                        var zV = zVec.clone().multiplyScalar(k);
                        
                        //Cube points
                        var cube = [null, null, null, null,
                                    null, null, null, null,];
                        
                        for (var c = 0; c < 8; c++) {
                            cube[c] = cubepts[c].clone().add(xV).add(yV).add(zV);
                        }
                        
                        var intersects = [null, null, null, null,
                                          null, null, null, null,
                                          null, null, null, null];
                                    
                        var v1, v2, idx;
                        var index = offset*12;      
                        //0 to 1
                        if (edgeIdx & 1) {
                            p1.addVectors(cubepts[0], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[1], xV).add(yV).add(zV);
                            v1 = grid[0];
                            v2 = grid[1];
                            idx = index+0;
                            intersects[0] = linearInterpolate(i,j,k,cube,grid,0,1,verts,vertnums,bit,isoval,smooth);
                        }
                        //1 to 3
                        if (edgeIdx & 2) {
                            p1.addVectors(cubepts[1], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[3], xV).add(yV).add(zV);
                            v1 = grid[1];
                            v2 = grid[3];
                            idx = index+1;
                            intersects[1] = linearInterpolate(i,j,k,cube,grid,1,3,verts,vertnums,bit,isoval,smooth);
                        }
                        //3 to 2
                        if (edgeIdx & 4) {
                            p1.addVectors(cubepts[3], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[2], xV).add(yV).add(zV);
                            v1 = grid[3];
                            v2 = grid[2];
                            idx = index+2;
                            intersects[2] = linearInterpolate(i,j,k,cube,grid,3,2,verts,vertnums,bit,isoval,smooth);
                        }
                        //2 to 0
                        if (edgeIdx & 8) {
                            p1.addVectors(cubepts[2], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[0], xV).add(yV).add(zV);
                            v1 = grid[2];
                            v2 = grid[0];
                            idx = index+3;
                            intersects[3] = linearInterpolate(i,j,k,cube,grid,2,0,verts,vertnums,bit,isoval,smooth);
                        }     
                        //4 to 5
                        if (edgeIdx & 16) {
                            p1.addVectors(cubepts[4], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[5], xV).add(yV).add(zV);
                            v1 = grid[4];
                            v2 = grid[5];
                            idx = index+4;
                            intersects[4] = linearInterpolate(i,j,k,cube,grid,4,5,verts,vertnums,bit,isoval,smooth);
                        }    
                        //5 to 7
                        if (edgeIdx & 32) {
                            p1.addVectors(cubepts[5], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[7], xV).add(yV).add(zV);
                            v1 = grid[5];
                            v2 = grid[7];
                            idx = index+5;
                            intersects[5] = linearInterpolate(i,j,k,cube,grid,5,7,verts,vertnums,bit,isoval,smooth);
                        }
                        //7 to 6
                        if (edgeIdx & 64) {
                            p1.addVectors(cubepts[7], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[6], xV).add(yV).add(zV);
                            v1 = grid[7];
                            v2 = grid[6];
                            idx = index+6;
                            intersects[6] = linearInterpolate(i,j,k,cube,grid,7,6,verts,vertnums,bit,isoval,smooth);
                        }
                        //6 to 4
                        if (edgeIdx & 128) {
                            p1.addVectors(cubepts[6], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[4], xV).add(yV).add(zV);
                            v1 = grid[6];
                            v2 = grid[4];
                            idx = index+7;
                            intersects[7] = linearInterpolate(i,j,k,cube,grid,6,4,verts,vertnums,bit,isoval,smooth);
                        }
                        //0 to 4
                        if (edgeIdx & 256) {
                            p1.addVectors(cubepts[0], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[4], xV).add(yV).add(zV);
                            v1 = grid[0];
                            v2 = grid[4];
                            idx = index+8;
                            intersects[8] = linearInterpolate(i,j,k,cube,grid,0,4,verts,vertnums,bit,isoval,smooth);
                        }
                        //1 to 5
                        if (edgeIdx & 512) {
                            p1.addVectors(cubepts[1], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[5], xV).add(yV).add(zV);
                            v1 = grid[1];
                            v2 = grid[5];
                            idx = index+9;
                            intersects[9] = linearInterpolate(i,j,k,cube,grid,1,5,verts,vertnums,bit,isoval,smooth);
                        }  
                        //3 to 7
                        if (edgeIdx & 1024) {
                            p1.addVectors(cubepts[3], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[7], xV).add(yV).add(zV);
                            v1 = grid[3];
                            v2 = grid[7];
                            idx = index+10;
                            intersects[10] = linearInterpolate(i,j,k,cube,grid,3,7,verts,vertnums,bit,isoval,smooth);
                        }
                        //2 to 6
                        if (edgeIdx & 2048) {
                            p1.addVectors(cubepts[2], xV).add(yV).add(zV);                        
                            p2.addVectors(cubepts[6], xV).add(yV).add(zV);
                            v1 = grid[2];
                            v2 = grid[6];
                            idx = index+11;
                            intersects[11] = linearInterpolate(i,j,k,cube,grid,2,6,verts,vertnums,bit,isoval,smooth);
                        }
                        
                        //add Vectors
                        
                        for (var itri = 0; itri < triangles.length / 3; ++itri) {
                            var trioffset = itri*3;
                            
                            var a = intersects[triangles[trioffset]];                        
                            var b = intersects[triangles[trioffset + 2]], c = intersects[triangles[trioffset + 1]];
                            
                            var vA = verts[a], vB = verts[b], vC = verts[c];
                            //var normA = norms[a], normB = norms[b], normC = norms[c];
                            //normA.subVectors(vA, vB);
                            //normC.subVectors(vC, vB);
                            
                            //normA.cross(normC).normalize();
                            //norms[b].copy(normA);
                            //norms[c].copy(normA);
                            
                            if (! smooth && itri > 0) {
                                faces.push(verts.length);
                                verts.push(vA);
                                faces.push(verts.length);
                                verts.push(vB);
                                faces.push(verts.length);
                                verts.push(vC);
                                
                            }
                            else {
                                //faces.push(verts.length);
                                faces.push(a);
                                //faces.push(verts.length);
                                faces.push(b);
                                //faces.push(verts.length);
                                faces.push(c);                           
                            }
    
                            
                        }
                            
                    }
    
                }
    
            }
            
    
            if (smooth) 
                laplacianSmooth(1, verts, faces);
                
            var color = neg ? new WebMol.Color(1,0,0) : new WebMol.Color(0,0,1);
            
            var shape = viewer.addShape({
                wireframe : false,
                color : color,
                alpha : 0.8,
                side : WebMol.FrontSide
            });
            
            
            viewer.addCustom(shape, {vertexArr:verts, 
                                     faceArr:faces,
                                     normalArr:[]});
                      
        }

        return atomStr;  
                
    };
    
    var linearInterpolate = function() {       
        
        var nY = 55, nZ = 40;
        
        return function(i, j, k, cube, grid, p1, p2, verts, vertnums, code, isoval, smooth) {
            
            var pt = new WebMol.Vector3();
            
            var v1 = grid[p1], v2 = grid[p2];
            var pt1 = cube[p1], pt2 = cube[p2];
            
            if (smooth) {
                
                var val1 = !!(code & (1 << p1));
                var val2 = !!(code & (1 << p2));
                
                var p = p1;
                if (!val1 && val2)
                    p = p2;
                    
                if (p & 1)
                    k++;
                if (p & 2)
                    j++;
                if (p & 4)
                    i++;
                    
                var index = (i*nY*nZ) + (j*nZ) + k;
                
                if (Math.abs(isoval-v1) < 0.000001)
                    pt = pt1.clone();
                else if (Math.abs(isoval-v2) < 0.000001) 
                    pt = pt2.clone();
                else if (Math.abs(v1 - v2) < 0.000001)
                    pt = pt1.clone().add(pt2).multiplyScalar(0.5);
                    
                else {            
                    pt.subVectors(pt2,pt1);
                    var scale = (isoval-v1)/(v2-v1);                   
                    pt.multiplyScalar(scale).add(pt1);            
                }     
                   
                if (vertnums[index] < 1) {
                    vertnums[index] = verts.length;
                    verts.push(pt);          
                    //norms.push(new WebMol.Vector3());  
                }          
                
                return vertnums[index];                      
            }

            else {
                
                pt.addVectors(pt1, pt2).multiplyScalar(0.5); 
                
                verts.push(pt);      
                
                return verts.length - 1;              
            }

        };

       
    }();
    

    
    laplacianSmooth = function(numiter, verts, faces) {
            var tps = new Array(verts.length);
            for ( var i = 0; i < verts.length; i++)
                    tps[i] = {
                            x : 0,
                            y : 0,
                            z : 0
                    };
            var vertdeg = new Array(20);
            var flagvert;
            for ( var i = 0; i < 20; i++)
                    vertdeg[i] = new Array(verts.length);
            for ( var i = 0; i < verts.length; i++)
                    vertdeg[0][i] = 0;
            for ( var i = 0; i < faces.length / 3; i++) {
                    var aoffset = i*3, boffset = i*3 + 1, coffset = i*3 + 2;
                    flagvert = true;
                    for ( var j = 0; j < vertdeg[0][faces[aoffset]]; j++) {
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
                    for ( var j = 0; j < vertdeg[0][faces[aoffset]]; j++) {
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
                    for (j = 0; j < vertdeg[0][faces[boffset]]; j++) {
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
                    for (j = 0; j < vertdeg[0][faces[boffset]]; j++) {
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
                    for (j = 0; j < vertdeg[0][faces[coffset]]; j++) {
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
            for ( var k = 0; k < numiter; k++) {
                    for ( var i = 0; i < verts.length; i++) {
                            if (vertdeg[0][i] < 3) {
                                    tps[i].x = verts[i].x;
                                    tps[i].y = verts[i].y;
                                    tps[i].z = verts[i].z;
                            } else if (vertdeg[0][i] == 3 || vertdeg[0][i] == 4) {
                                    tps[i].x = 0;
                                    tps[i].y = 0;
                                    tps[i].z = 0;
                                    for (j = 0; j < vertdeg[0][i]; j++) {
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
                                    for ( var j = 0; j < vertdeg[0][i]; j++) {
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
                    for ( var i = 0; i < verts.length; i++) {
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
    
    // The constructor
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
        camera.position = new TV3(0, 0, CAMERA_Z);
        camera.lookAt(new TV3(0, 0, 0));
        
        var raycaster = new WebMol.Raycaster(new TV3(0,0,0), new TV3(0,0,0));
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
            var atoms = getAtomsFromSel(sel).concat(shapes);
            var allatoms = getAtomsFromSel({}).concat(shapes);
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

            label.dispose();
            modelGroup.remove(label.sprite);                       
        };
        
        //Modify label style
        this.setLabelStyle = function(label, stylespec) {   
             
            label.dispose();
            label.stylespec = stylespec;
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;
            
        };
        
        //Change label text
        this.setLabelText = function(label, text) {
         
            label.dispose();
            label.text = text;
            label.setContext();
            modelGroup.add(label.sprite);
            
            return label;

        };
        
        //Add generic GLShape to viewer
        this.addShape = function(shapeSpec) {
            shapeSpec = shapeSpec || {};
            var shape = new WebMol.GLShape(shapeSpec);
            shapes.push(shape);
            
            return shape;
              
        };
        
        this.addSphere = function(shape, spec) {
            spec = spec || {};
            shape.addSphere(spec);      
        };
        
        this.addArrow = function(shape, spec) {
            spec = spec || {};
            shape.addArrow(spec);
        };
        
        //Add custom shape component from user supplied function
        this.addCustom = function(shape, spec) {                            
            spec = spec || {};
            shape.addCustom(spec);                            
        };

        // given molecular data and its format (pdb, sdf, xyz or mol2)
        // create a model and add it, returning the model identifier
        this.addModel = function(data, format) {
            var m = new WebMol.GLModel(models.length, defaultcolors);
            if (format === "cube") {
                data = parseCube(data, this);
            }
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
            geoGroup.faceidx = faces.length*3;
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
            for ( var i = 0; i < faces.length; i++) {
                
                faceoffset = i*3;
                var a = faces[i].a, b = faces[i].b, c = faces[i].c;
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
                
                vA = new TV3(verts[offsetA], verts[offsetA+1], verts[offsetA+2]);
                vB = new TV3(verts[offsetB], verts[offsetB+1], verts[offsetB+2]);
                vC = new TV3(verts[offsetC], verts[offsetC+1], verts[offsetC+2]);
                
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
            //TODO: Asynchronously generate geometryGroups (not separate meshes) and merge them into a single geometry
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
        
        //TODO: Probably want to keep modelgroup hidden in deployment
        //  used for debugging
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
