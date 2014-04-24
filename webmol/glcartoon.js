//glcartoon.js
//This contains all the routines for rendering a cartoon given a set
//of atoms with assigned secondary structure
//TODO: secondary structure calculation

var WebMol = WebMol || {};

WebMol.drawCartoon = (function() {

    var axisDIV = 5; // 3 still gives acceptable quality
    var strandDIV = 6;
    var nucleicAcidStrandDIV = 4;
    var tubeDIV = 8;
    var coilWidth = 0.3;
    var helixSheetWidth = 1.3;
    var nucleicAcidWidth = 0.8;
    var thickness = 0; 

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
                points.push(new TV3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2,
                        (p1.z + p2.z) / 2));
            else
                points.push(p1);
        }
        points.push(_points[_points.length - 1]);

        
        for ( var i = -1, size = points.length; i <= size - 3; i++) {
            var p0 = points[(i == -1) ? 0 : i];
            var p1 = points[i + 1], p2 = points[i + 2];
            var p3 = points[(i == size - 3) ? size - 1 : i + 3];
            var v0 = new TV3().subVectors(p2, p0).multiplyScalar(0.5);
            var v1 = new TV3().subVectors(p3, p1).multiplyScalar(0.5);

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
                        
                var pt = new TV3(x, y, z);
                
                var atomIndex = Math.floor( (ret.length+2) / DIV);
                
                if (_points[atomIndex].atom !== undefined)
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
        
        for ( var i = 0, lim = p1.length; i < lim; i++) {
           
            geoGroup = geo.updateGeoGroup(2);
            offset = geoGroup.vertices, vertoffset = offset*3;
            
            geoGroup.__vertexArray[vertoffset] = p1[i].x;
            geoGroup.__vertexArray[vertoffset+1] = p1[i].y;
            geoGroup.__vertexArray[vertoffset+2] = p1[i].z;
            
            geoGroup.__vertexArray[vertoffset+3] = p2[i].x;
            geoGroup.__vertexArray[vertoffset+4] = p2[i].y;
            geoGroup.__vertexArray[vertoffset+5] = p2[i].z;
           
            if (i > 0) {
                var faces = [offset, offset + 1, offset - 1, offset - 2];
                var faceoffset = geoGroup.faceidx;
                
                geoGroup.__faceArray[faceoffset] = faces[0]; geoGroup.__faceArray[faceoffset+1] = faces[1]; geoGroup.__faceArray[faceoffset+2] = faces[3];
                geoGroup.__faceArray[faceoffset+3] = faces[1]; geoGroup.__faceArray[faceoffset+4] = faces[2]; geoGroup.__faceArray[faceoffset+5] = faces[3];

                geoGroup.faceidx += 6;
            }
            
            geoGroup.vertices += 2;
        }
        
        setUpNormals(geo);
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
        
        var faces = [ [ 0, 2, -6, -8 ], [ -4, -2, 6, 4 ], [ 7, 3, -5, -1 ],
                [ -3, -7, 1, 5 ] ];
                
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
        
        //geo.computeFaceNormals();
        //geo.computeVertexNormals(false);
        geo.initTypedArrays();
        setUpNormals(geo);
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
        var currentChain, currentResi, currentCA;
        var prevCO = null, ss = null, ssborder = false;

        for ( var i in atomlist) {
            var atom = atomlist[i];
            if (atom == undefined)
                continue;

            if ((atom.atom == 'O' || atom.atom == 'CA') && !atom.hetflag) {
                if (atom.atom == 'CA') {
                    if (currentChain != atom.chain
                            || currentResi + 1 != atom.resi) {
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
                    currentCA = new TV3(atom.x, atom.y, atom.z);
                    currentAtom = atom;
                    currentChain = atom.chain;
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
                    var O = new TV3(atom.x, atom.y, atom.z);
                    O.sub(currentCA);
                    O.normalize(); // can be omitted for performance
                    O.multiplyScalar((ss == 'c') ? coilWidth : helixSheetWidth);
                    if (prevCO != undefined && O.dot(prevCO) < 0)
                        O.negate();
                    prevCO = O;
                    for ( var j = 0; j < num; j++) {
                        var delta = -1 + 2 / (num - 1) * j;
                        var v = new TV3(currentCA.x + prevCO.x * delta,
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
