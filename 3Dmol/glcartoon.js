//glcartoon.js
//This contains all the routines for rendering a cartoon given a set
//of atoms with assigned secondary structure

//TODO: generate normals directly in drawStrip and drawThinStrip

var $3Dmol = $3Dmol || {};


/**@typedef CartoonStyleSpec
 * @prop {string} color - solid color, may specify as 'spectrum'
 * @prop {string} style - style of cartoon rendering (currently just default and trace)
 */

/**
 * @ignore
 * @param {$3Dmol.Object3D} group
 * @param {AtomSpec} atomlist
 * @param {$3Dmol.Gradient} gradientscheme
 */
$3Dmol.drawCartoon = (function() {

    var axisDIV = 5; // 3 still gives acceptable quality
    var strandDIV = 6;
    var nucleicAcidStrandDIV = 4;
    var tubeDIV = 8;
    var coilWidth = 0.3;
    var helixSheetWidth = 1.3;
    var nucleicAcidWidth = 0.8;
    var defaultThickness = 0.4; 

    // helper functions

    // Catmull-Rom subdivision
    var subdivide = function(_points, DIV) { // points as Vector3
        var ret = [];
        var points = _points;
        points = []; // Smoothing test
        points.push(_points[0]);
        
        var i, lim, size;
        var p0, p1, p2, p3, v0, v1;
        
        for (i = 1, lim = _points.length - 1; i < lim; i++) {
            p1 = _points[i]; p2 = _points[i + 1];
            if (p1.smoothen)
                points.push(new $3Dmol.Vector3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2,
                        (p1.z + p2.z) / 2));
            else
                points.push(p1);
        }
        points.push(_points[_points.length - 1]);

        
        for (i = -1, size = points.length; i <= size - 3; i++) {
            p0 = points[(i === -1) ? 0 : i];
            p1 = points[i + 1]; p2 = points[i + 2];
            p3 = points[(i === size - 3) ? size - 1 : i + 3];
            v0 = new $3Dmol.Vector3().subVectors(p2, p0).multiplyScalar(0.5);
            v1 = new $3Dmol.Vector3().subVectors(p3, p1).multiplyScalar(0.5);

            for ( var j = 0; j < DIV; j++) {
                var t = 1.0 / DIV * j;
                var x = p1.x + t * v0.x + t * t * 
                        (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x) + t * t * t *
                        (2 * p1.x - 2 * p2.x + v0.x + v1.x);
                var y = p1.y + t * v0.y + t * t * 
                        (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y) + t * t * t *
                        (2 * p1.y - 2 * p2.y + v0.y + v1.y);
                var z = p1.z + t * v0.z + t * t * 
                        (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z) + t * t * t * 
                        (2 * p1.z - 2 * p2.z + v0.z + v1.z);
                        
                var pt = new $3Dmol.Vector3(x, y, z);
                
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
    
        var geo = new $3Dmol.Geometry(true);       
        var offset, vertoffset;
        var color;

        
        for ( var i = 0, lim = p1.length; i < lim; i++) {
            
            color = $3Dmol.CC.color(colors[Math.round((i - 1) / div)]);
           
            geoGroup = geo.updateGeoGroup(2);
            var vertexArray = geoGroup.vertexArray;
            var colorArray = geoGroup.colorArray;
            var faceArray = geoGroup.faceArray;
            offset = geoGroup.vertices; vertoffset = offset*3;
            
            vertexArray[vertoffset] = p1[i].x;
            vertexArray[vertoffset+1] = p1[i].y;
            vertexArray[vertoffset+2] = p1[i].z;
            
            vertexArray[vertoffset+3] = p2[i].x;
            vertexArray[vertoffset+4] = p2[i].y;
            vertexArray[vertoffset+5] = p2[i].z;
            
            for (var j = 0; j < 6; ++j) {                
                colorArray[vertoffset+3*j] = color.r; colorArray[vertoffset+1+3*j] = color.g; colorArray[vertoffset+2+3*j] = color.b;
            }            
           
            if (i > 0) {
                var faces = [offset, offset + 1, offset - 1, offset - 2];
                var faceoffset = geoGroup.faceidx;
                
                faceArray[faceoffset] = faces[0]; faceArray[faceoffset+1] = faces[1]; faceArray[faceoffset+2] = faces[3];
                faceArray[faceoffset+3] = faces[1]; faceArray[faceoffset+4] = faces[2]; faceArray[faceoffset+5] = faces[3];

                geoGroup.faceidx += 6;
            }
            
            geoGroup.vertices += 2;
        }
        
        geo.initTypedArrays();
        geo.setUpNormals();
        
        var material = new $3Dmol.MeshLambertMaterial();
        material.vertexColors = $3Dmol.FaceColors;
                material.side = $3Dmol.DoubleSide;
        var mesh = new $3Dmol.Mesh(geo, material);
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

        var geo = new $3Dmol.Geometry(true);
        
        //var vs = geo.vertices, fs = geo.faces;
                var vs = [], fs = [];
        var axis, p1v, p2v, a1v, a2v;
        
        var faces = [ [ 0, 2, -6, -8 ], [ -4, -2, 6, 4 ], [ 7, -1, -5, 3 ],
                [ -3, 5, 1, -7 ] ];
                
        var offset, vertoffset, faceoffset;
        var color;
        var currentAtom, lastAtom;
        var i, lim, j;
        var face1, face2, face3;
        var geoGroup;
        
        for (i = 0, lim = p1.length; i < lim; i++) {
        
            color = $3Dmol.CC.color(colors[Math.round((i - 1) / div)]);
            
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
            
            geoGroup = geo.updateGeoGroup(8);
            var vertexArray = geoGroup.vertexArray;
            var colorArray = geoGroup.colorArray;
            var faceArray = geoGroup.faceArray;
            offset = geoGroup.vertices; vertoffset = offset*3;
            
            vertexArray[vertoffset] = p1v.x; vertexArray[vertoffset+1] = p1v.y; vertexArray[vertoffset+2] = p1v.z;
            vertexArray[vertoffset+3] = p1v.x; vertexArray[vertoffset+4] = p1v.y; vertexArray[vertoffset+5] = p1v.z;
            vertexArray[vertoffset+6] = p2v.x; vertexArray[vertoffset+7] = p2v.y; vertexArray[vertoffset+8] = p2v.z;
            vertexArray[vertoffset+9] = p2v.x; vertexArray[vertoffset+10] = p2v.y; vertexArray[vertoffset+11] = p2v.z;
            vertexArray[vertoffset+12] = a1v.x; vertexArray[vertoffset+13] = a1v.y; vertexArray[vertoffset+14] = a1v.z;
            vertexArray[vertoffset+15] = a1v.x; vertexArray[vertoffset+16] = a1v.y; vertexArray[vertoffset+17] = a1v.z;
            vertexArray[vertoffset+18] = a2v.x; vertexArray[vertoffset+19] = a2v.y; vertexArray[vertoffset+20] = a2v.z;
            vertexArray[vertoffset+21] = a2v.x; vertexArray[vertoffset+22] = a2v.y; vertexArray[vertoffset+23] = a2v.z;
            
            for (j = 0; j < 8; ++j) {                
                colorArray[vertoffset+3*j] = color.r; colorArray[vertoffset+1+3*j] = color.g; colorArray[vertoffset+2+3*j] = color.b;
            }
            
            if (i > 0) {
             
                //both points have distinct atoms
                var diffAtoms = ((lastAtom !== undefined && currentAtom !== undefined) && lastAtom.serial !== currentAtom.serial);
                
                for (j = 0; j < 4; j++ ) {
                
                    var face = [offset + faces[j][0], offset + faces[j][1], offset + faces[j][2], offset + faces[j][3]];
                    
                    faceoffset = geoGroup.faceidx;    
                    
                    faceArray[faceoffset] = face[0]; faceArray[faceoffset+1] = face[1]; faceArray[faceoffset+2] = face[3];
                    faceArray[faceoffset+3] = face[1]; faceArray[faceoffset+4] = face[2]; faceArray[faceoffset+5] = face[3];
                    
                    geoGroup.faceidx += 6;
                    
                    if (currentAtom.clickable || lastAtom.clickable) {
                        
                        var p1a = vs[face[3]].clone(), p1b = vs[face[0]].clone(),
                            p2a = vs[face[2]].clone(), p2b = vs[face[1]].clone();
                        
                        p1a.atom = vs[face[3]].atom || null; //should be same
                        p2a.atom = vs[face[2]].atom || null; 
                        
                        
                        p1b.atom = vs[face[0]].atom || null; //should be same                      
                        p2b.atom = vs[face[1]].atom || null; 
                        
                        if (diffAtoms) {
                            var m1 = p1a.clone().add(p1b).multiplyScalar(0.5);
                            var m2 = p2a.clone().add(p2b).multiplyScalar(0.5);
                            var m = p1a.clone().add(p2b).multiplyScalar(0.5);
                            
                            if (j % 2 === 0)
                            {
                                if (lastAtom.clickable) {
                                    face1 = new $3Dmol.Triangle(m1, m, p1a);
                                    face2 = new $3Dmol.Triangle(m2, p2a, m);
                                    face3 = new $3Dmol.Triangle(m, p2a, p1a);
                                    lastAtom.intersectionShape.triangle.push(face1);
                                    lastAtom.intersectionShape.triangle.push(face2);
                                    lastAtom.intersectionShape.triangle.push(face3);
                                }
                                
                                if (currentAtom.clickable) {
                                    face1 = new $3Dmol.Triangle(p1b, p2b, m);
                                    face2 = new $3Dmol.Triangle(p2b, m2, m);
                                    face3 = new $3Dmol.Triangle(p1b, m, m1);
                                    currentAtom.intersectionShape.triangle.push(face1);
                                    currentAtom.intersectionShape.triangle.push(face2);
                                    currentAtom.intersectionShape.triangle.push(face3);
                                }
                            }
                            else {
                                if (currentAtom.clickable) {
                                    face1 = new $3Dmol.Triangle(m1, m, p1a);
                                    face2 = new $3Dmol.Triangle(m2, p2a, m);
                                    face3 = new $3Dmol.Triangle(m, p2a, p1a);
                                    currentAtom.intersectionShape.triangle.push(face1);
                                    currentAtom.intersectionShape.triangle.push(face2);
                                    currentAtom.intersectionShape.triangle.push(face3);
                                }
                                
                                if (lastAtom.clickable) {
                                    face1 = new $3Dmol.Triangle(p1b, p2b, m);
                                    face2 = new $3Dmol.Triangle(p2b, m2, m);
                                    face3 = new $3Dmol.Triangle(p1b, m, m1);
                                    lastAtom.intersectionShape.triangle.push(face1);
                                    lastAtom.intersectionShape.triangle.push(face2);
                                    lastAtom.intersectionShape.triangle.push(face3);
                                }                          
                            }
                            
                        }
                        
                        //face for single atom
                        else if (currentAtom.clickable) {
                            face1 = new $3Dmol.Triangle(p1b, p2b, p1a);
                            face2 = new $3Dmol.Triangle(p2b, p2a, p1a);
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
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;
        var faceArray = geoGroup.faceArray;
        offset = geoGroup.vertices; vertoffset = offset*3; faceoffset = geoGroup.faceidx;
        
        for (i = 0; i < 4; i++) {
            vs.push(vs[i * 2]);
            vs.push(vs[vsize + i * 2]);
            
            var v1 = vs[i * 2], v2 = vs[vsize + i * 2];
            
            vertexArray[vertoffset+6*i] = v1.x; vertexArray[vertoffset+1+6*i] = v1.y; vertexArray[vertoffset+2+6*i] = v1.z;
            vertexArray[vertoffset+3+6*i] = v2.x; vertexArray[vertoffset+4+6*i] = v2.y; vertexArray[vertoffset+5+6*i] = v2.z;
            
            colorArray[vertoffset+6*i] = color.r; colorArray[vertoffset+1+6*i] = color.g; colorArray[vertoffset+2+6*i] = color.b;
            colorArray[vertoffset+3+6*i] = color.r; colorArray[vertoffset+4+6*i] = color.g; colorArray[vertoffset+5+6*i] = color.b;

        }
        
        vsize += 8;
                
        face1 = [offset, offset + 2, offset + 6, offset + 4];
        face2 = [offset + 1, offset + 5, offset + 7, offset + 3];
        
        faceArray[faceoffset] = face1[0]; faceArray[faceoffset+1] = face1[1]; faceArray[faceoffset+2] = face1[3];
        faceArray[faceoffset+3] = face1[1]; faceArray[faceoffset+4] = face1[2]; faceArray[faceoffset+5] = face1[3];
        faceArray[faceoffset+6] = face2[0]; faceArray[faceoffset+7] = face2[1]; faceArray[faceoffset+8] = face2[3];
        faceArray[faceoffset+9] = face2[1]; faceArray[faceoffset+10] = face2[2]; faceArray[faceoffset+11] = face2[3];
        
        geoGroup.faceidx += 12;
        geoGroup.vertices += 8;
        
        //TODO: Add intersection planes for caps
        
        geo.initTypedArrays();
        geo.setUpNormals();
        
        var material = new $3Dmol.MeshLambertMaterial();
        material.vertexColors = $3Dmol.FaceColors;
        material.side = $3Dmol.DoubleSide;
        var mesh = new $3Dmol.Mesh(geo, material);
        group.add(mesh);
        
    };

    //TODO: Need to update this (will we ever use this?)
    var drawSmoothCurve = function(group, _points, width, colors, div) {
        if (_points.length === 0)
            return;

        div = (div === undefined) ? 5 : div;

        var geo = new $3Dmol.Geometry();
        var points = subdivide(_points, div);
                /*
        for ( var i = 0; i < points.length; i++) {
            geo.vertices.push(points[i]);
            geo.colors.push($3Dmol.color(colors[(i == 0) ? 0 : Math.round((i - 1)
                    / div)]));
        }
                */
        var lineMaterial = new $3Dmol.LineBasicMaterial({
            linewidth : width
        });
        lineMaterial.vertexColors = true;
        var line = new $3Dmol.Line(geo, lineMaterial);
        line.type = $3Dmol.LineStrip;
        group.add(line);
    };

    var drawStrand = function(group, atomList, num, div, fill, coilWidth, helixSheetWidth, doNotSmoothen, gradientScheme, geo)
    {
        num = num || strandDIV;
        div = div || axisDIV;
        doNotSmoothen = !!(doNotSmoothen);

        var cartoonAtoms = ["CA", "O", "P", "O5'", "C4'", "OP2", "N1", "N3"];
        var cartoon, curr, next, currColor, nextColor, thickness, i;
        var backbonePt, prevOrientPt, baseStartPt, baseEndPt;
        var traceGeo = null;
        var colors = [];
        var points = [];
        for (var i = 0; i < num; i++)
            points[i] = [];

        for (i in atomList)
        {
            next = atomList[i];
            
            if (next === undefined || $.inArray(next.atom, cartoonAtoms) === -1 || next.hetflag)
                continue; // skip array holes, heteroatoms, and atoms not involved in cartoon drawing

            // determine cartoon style
            cartoon = next.style.cartoon;
            if (cartoon.style === "trace") // draw cylinders connecting consecutive 'backbone' atoms
            {
                /* "trace" style just draws cylinders between consecutive 'backbone' atoms,
                    such as alpha carbon for polypeptides and phosphorus for DNA. */

                if (!traceGeo) traceGeo = new $3Dmol.Geometry(true);

                if (next.atom === "CA" || next.atom === "P")
                {
                    // determine cylinder color
                    if (gradientScheme)
                        nextColor = gradientScheme.valueToHex(next.resi, gradientScheme.range());
                    else
                        nextColor = $3Dmol.getColorFromStyle(next, cartoon).getHex();
                    colors.push(nextColor);

                    // determine cylinder thickness
                    if ($.isNumeric(cartoon.thickness))
                        thickness = cartoon.thickness;
                    else
                        thickness = defaultThickness;

                    /* do not draw connections between different chains, but ignore
                       differences in reschain to properly support CA-only files */
                    if (curr && curr.chain === next.chain && curr.resi + 1 === next.resi)
                    {
                        // if both atoms are same color, draw single cylinder
                        if (nextColor == currColor)
                        {
                            var color = $3Dmol.CC.color(nextColor);
                            $3Dmol.GLDraw.drawCylinder(traceGeo, curr, next, thickness, color, true, true);
                        }

                        else // otherwise draw cylinders for each color (split down the middle)
                        {
                            var midpoint = new $3Dmol.Vector3().addVectors(curr, next).multiplyScalar(0.5);
                            var color1 = $3Dmol.CC.color(currColor);
                            var color2 = $3Dmol.CC.color(nextColor);
                            $3Dmol.GLDraw.drawCylinder(traceGeo, curr, midpoint, thickness, color1, true, false);
                            $3Dmol.GLDraw.drawCylinder(traceGeo, midpoint, next, thickness, color2, false, true);
                        } // note that an atom object can be duck-typed as a 3-vector in this case
                    }

                    curr = next;
                    currColor = nextColor;
                }
            }

            else // draw default-style cartoons based on secondary structure
            {
                // draw bases for dna ladder
                if (baseEndPt && (next.atom === "P" || next.atom === "O5'"))
                {
                    //console.log(next.chain + " " + next.resi + " " + next.resn);
                    if (curr && curr.chain === next.chain) // use the midpoint
                        baseStartPt = new $3Dmol.Vector3().addVectors(curr, next).multiplyScalar(0.5);
                    else if (next.atom === "P") // when starting a new chain
                        baseStartPt = new $3Dmol.Vector3(next.x, next.y, next.z);
                    else
                        baseStartPt = new $3Dmol.Vector3(curr.x, curr.y, curr.z);
                    var startFix = baseStartPt.clone().sub(baseEndPt).multiplyScalar(0.04);
                    baseStartPt.add(startFix);
                    if (!currColor) currColor = next.color;
                    $3Dmol.GLDraw.drawCylinder(geo, baseStartPt, baseEndPt, 0.4, $3Dmol.CC.color(currColor), false, true);
                    baseStartPt = null;
                    baseEndPt = null;
                }

                if (next.atom === "CA" || next.atom === "P" || next.atom === "O5'" && next.resi === 1)
                {

                    if (curr != undefined &&
                        (curr.chain != next.chain || curr.resi + 1 != next.resi || curr.reschain != next.reschain))
                    {
                        // reached end of a chain of connected residues, so draw accumulated points
                        for (i = 0; !thickness && i < num; i++)
                            drawSmoothCurve(group, points[i], 1, colors, div);
                        if (fill)
                            drawStrip(group, points[0], points[num - 1], colors, div, thickness);

                        // forget features of previous chain (points, colors, location of most recent oxygen)
                        points = [];
                        for (i = 0; i < num; i++)
                            points[i] = [];
                        colors = [];
                        //backbonePt = null;
                        //prevOrientPt = null;
                    }

                    // determine segment color and thickness
                    if (gradientScheme)
                        nextColor = gradientScheme.valueToHex(next.resi, gradientScheme.range());
                    else
                        nextColor = $3Dmol.getColorFromStyle(next, cartoon).getHex();
                    colors.push(nextColor);
                    if ($.isNumeric(cartoon.thickness))
                        thickness = cartoon.thickness;
                    else
                        thickness = defaultThickness;

                    curr = next; // advance current-backbone-atom pointer
                    backbonePt = new $3Dmol.Vector3(curr.x, curr.y, curr.z);
                    backbonePt.atom = curr.atom;
                    currColor = nextColor; // used for dna bases, which lag due to midpoint calculation

                    // click handling
                    if (next.clickable === true &&
                        (next.intersectionShape === undefined || next.intersectionShape.triangle === undefined)) 
                        next.intersectionShape = {sphere : null, cylinder : [], line : [], triangle : []};
                }

                else if (next.atom === "O" || next.atom === "OP2" || next.atom === "C4'" && next.resi === 1)
                {
                    var orientPt = new $3Dmol.Vector3(next.x, next.y, next.z);
                    orientPt.sub(backbonePt);
                    orientPt.normalize();

                    var widthScalar;
                    if (curr.ss === "c")
                    {
                        if (curr.atom === "P")
                            widthScalar = nucleicAcidWidth;
                        else
                            widthScalar = coilWidth;
                    } else
                        widthScalar = helixSheetWidth;
                    orientPt.multiplyScalar(widthScalar);

                    if (prevOrientPt && orientPt.dot(prevOrientPt) < 0)
                        orientPt.negate();
                    prevOrientPt = orientPt;

                    for (i = 0; i < num; i++)
                    {
                        var delta = -1 + (2 * i) / (num - 1); // produces num increments from -1 to 1
                        var v = new $3Dmol.Vector3(backbonePt.x + delta * orientPt.x,
                                                   backbonePt.y + delta * orientPt.y,
                                                   backbonePt.z + delta * orientPt.z);
                        v.atom = curr;
                        if (!doNotSmoothen && curr.ss === "s")
                            v.smoothen = true;
                        points[i].push(v);
                    }
                }

                // dna base end point is different for pyramidines and purines
                else if ((next.atom === "N1" && (next.resn === " DG" || next.resn ===  " DA") ||
                         next.atom === "N3" && (next.resn === " DC" || next.resn === " DT")))
                {
                    baseEndPt = new $3Dmol.Vector3(next.x, next.y, next.z);
                    baseEndPt.resi = next.resi;
                    baseEndPt.chain = next.chain;
                }
            }  
        }

        // for default style, draw the last chain
        for (i = 0; !thickness && i < num; i++)
            drawSmoothCurve(group, points[i], 1, colors, div);
        if (fill)
            drawStrip(group, points[0], points[num - 1], colors, div, thickness);

        if (traceGeo) // generate mesh for trace geometry
        {
            var traceMaterial = new $3Dmol.MeshLambertMaterial();
            traceMaterial.vertexColors = $3Dmol.FaceColors;
            traceMaterial.side = $3Dmol.DoubleSide;
            var traceMesh = new $3Dmol.Mesh(traceGeo, traceMaterial);
            group.add(traceMesh);
        }
    };

    // actual function call
    var drawCartoon = function(group, atomlist, geo, gradientscheme) {
        
        drawStrand(group, atomlist, 2, undefined, true, coilWidth, helixSheetWidth,
                false, gradientscheme, geo);
    };

    return drawCartoon;
})();
