// A model is a collection of related atoms.  Bonds are only allowed between
//atoms in the same model.  An atom is uniquely specified by its model id and
//its serial number.
//A glmodel knows how to apply the styles on each atom to create a gl object

var $3Dmol = $3Dmol || {};

/**
 * GLModel represents a group of related atoms
 * @constructor 
 * @param {number=} mid 
 * @param {Object=} defaultcolors Object defining default atom colors as atom => color property value pairs
 * @see $3Dmol.download
 */
$3Dmol.GLModel = (function() {

    // class variables go here
    var defaultAtomStyle = {
        line : {}
    };

    var Nucleotides = [ '  G', '  A', '  T', '  C', '  U', ' DG', ' DA', ' DT',
            ' DC', ' DU' ];

    var defaultlineWidth = 1.0;

    // Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
    var vdwRadii = {
        "H" : 1.2,
        "Li" : 1.82,
        "LI" : 1.82,
        "Na" : 2.27,
        "NA" : 2.27,
        "K" : 2.75,
        "C" : 1.7,
        "N" : 1.55,
        "O" : 1.52,
        "F" : 1.47,
        "P" : 1.80,
        "S" : 1.80,
        "CL" : 1.75,
        "Cl" : 1.75,
        "BR" : 1.85,
        "Br" : 1.85,
        "SE" : 1.90,
        "Se" : 1.90,
        "ZN" : 1.39,
        "Zn" : 1.39,
        "CU" : 1.4,
        "Cu" : 1.4,
        "NI" : 1.63,
        "Ni" : 1.63
    };

    var validAtomSpecs = [
        "resn", // Parent residue name
        "x", // Atom's x coordinate
        "y", // Atom's y coordinate
        "z", // Atom's z coordinate
        "color", // Atom's color, as hex code
        "surfaceColor", // Hex code for color to be used for surface patch over this atom
        "elem", // Element abbreviation (e.g. 'H', 'Ca', etc)
        "hetflag", // Set to true if atom is a heteroatom
        "chain", // Chain this atom belongs to, if specified in input file (e.g 'A' for chain A)
        "resi", // Residue number 
        "icode",
        "rescode",
        "serial", // Atom's serial id numbermodels
        "atom", // Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon)
        "bonds", // Array of atom ids this atom is bonded to
        "ss", // Secondary structure identifier (for cartoon render; e.g. 'h' for helix)
        "singleBonds", // true if this atom forms only single bonds or no bonds at all
        "bondOrder", // Array of this atom's bond orders, corresponding to bonds identfied by 'bonds'
        "properties", // Optional mapping of additional properties
        "b", // Atom b factor data
        "pdbline", // If applicable, this atom's record entry from the input PDB file (used to output new PDB from models)
        "clickable", // Set this flag to true to enable click selection handling for this atom
        "callback", // Callback click handler function to be executed on this atom and its parent viewer
        "invert" // for selection, inverts the meaning of the selection
    ];

    var validAtomSelectionSpecs = validAtomSpecs.concat([  // valid atom specs are ok too
        "model", // a single model or list of models from which atoms should be selected
        "bonds", // overloaded to select number of bonds, e.g. {bonds: 0} will select all nonbonded atoms
        "predicate", // user supplied function that gets passed an {AtomSpec} and should return true if the atom should be selected
        "invert", // if set, inverts the meaning of the selection
        "byres", // if set, expands the selection to include all atoms of any residue that has any atom selected
        "expand", // expands the selection to include all atoms within a given distance from the selection
        "within" // intersects the selection with the set of atoms within a given distance from another selection
    ]);

    var validAtomStyleSpecs = [
        "line", // draw bonds as lines
        "cross", // draw atoms as crossed lines (aka stars)
        "stick", // draw bonds as capped cylinders
        "sphere", // draw atoms as spheres
        "cartoon" // draw cartoon representation of secondary structure
    ];

    // class functions

    // return true if a and b represent the same style
    var sameObj = function(a,b) {
        if(a && b)
            return JSON.stringify(a) == JSON.stringify(b);
        else
            return a == b;
    };    

   
    function GLModel(mid, defaultcolors) {
        // private variables
        var atoms = [];
        var frames = [];
        var id = mid;
        var hidden = false;
        var molObj = null;
        var renderedMolObj = null;
        var lastColors = null;
        var modelData = {};
        var idMatrix = new $3Dmol.Matrix4();
        var dontDuplicateAtoms;
        var defaultColor = $3Dmol.elementColors.defaultColor;
        
        var ElementColors = (defaultcolors) ? defaultcolors : $3Dmol.elementColors.defaultColors;


        // drawing functions must be associated with model object since
        // geometries can't span multiple canvases

        // sphere drawing
        var defaultSphereRadius = 1.5;

        // return proper radius for atom given style
        /** 
         * 
         * @param {AtomSpec} atom
         * @param {atomstyle} style
         * @return {number} 
         * 
         */
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

        // cross drawing
        /** @typedef CrossStyleSpec
         * @prop {boolean} hidden - do not show 
         * @prop {number} linewidth 
         * @prop {number} radius 
         * @prop {number} scale - scale radius by specified amount
         * @prop {ColorschemeSpec} colorscheme - element based coloring
         * @prop {ColorSpec} color - fixed coloring, overrides colorscheme
         */
        
        /**
         * 
         * @param {AtomSpec} atom
         * @param {$3Dmol.Geometry[]} geos
         */
        var drawAtomCross = function(atom, geos) {
            if (!atom.style.cross)
                return;
            var style = atom.style.cross;
            if (style.hidden)
                return;
            var linewidth = (style.linewidth || defaultlineWidth);
            if (!geos[linewidth])
                geos[linewidth] = new $3Dmol.Geometry();
                
            var geoGroup = geos[linewidth].updateGeoGroup(6);
            
            var delta = getRadiusFromStyle(atom, style);

            var points = [ [ delta, 0, 0 ], [ -delta, 0, 0 ], [ 0, delta, 0 ],
                    [ 0, -delta, 0 ], [ 0, 0, delta ], [ 0, 0, -delta ] ];

            var clickable = atom.clickable;
            if (clickable && atom.intersectionShape === undefined)
                atom.intersectionShape = {sphere : [], cylinder : [], line : []};
            
            var c = $3Dmol.getColorFromStyle(atom, style);
            
            var vertexArray = geoGroup.vertexArray;
            var colorArray = geoGroup.colorArray;
            
            for ( var j = 0; j < 6; j++) {
                
                var offset = geoGroup.vertices*3;
                
                geoGroup.vertices++;
                vertexArray[offset] = atom.x + points[j][0];
                vertexArray[offset+1] = atom.y + points[j][1];
                vertexArray[offset+2] = atom.z + points[j][2];
                colorArray[offset] = c.r;
                colorArray[offset+1] = c.g;
                colorArray[offset+2] = c.b;
                
                if (clickable){
                    var point = new $3Dmol.Vector3(points[j][0], points[j][1], points[j][2]);
                    
                    //decrease cross size for selection to prevent misselection from atom overlap
                    point.multiplyScalar(0.1);
                    point.set(point.x+atom.x, point.y+atom.y, point.z+atom.z);
                    atom.intersectionShape.line.push(point);
                }

            }
                        
        };

        //from atom, return a normalized vector v that is orthogonal and along which
        //it is appropraite to draw multiple bonds
        var getSideBondV = function(atom, atom2, i) {

            var p1 = new $3Dmol.Vector3(atom.x, atom.y, atom.z);
            var p2 = new $3Dmol.Vector3(atom2.x, atom2.y, atom2.z);

            var dir = p2.clone();
            var v = null;
            dir.sub(p1);

            var p1a, p1b, p2a, p2b;
            var i2, j2, atom3, p3, dir2;
            if (atom.bonds.length === 1) {
                if (atom2.bonds.length === 1) {
                    v = dir.clone();
                    if (Math.abs(v.x) > 0.0001)
                        v.y += 1;
                    else
                        v.x += 1;
                } else {
                    i2 = (i + 1) % atom2.bonds.length;
                    j2 = atom2.bonds[i2];
                    atom3 = atoms[j2];
                    p3 = new $3Dmol.Vector3(atom3.x, atom3.y, atom3.z);

                    dir2 = p3.clone();
                    dir2.sub(p1);

                    v = dir2.clone();
                    v.cross(dir);
                }
            } else {
                // get vector 2 different neighboring atom
                i2 = (i + 1) % atom.bonds.length;
                j2 = atom.bonds[i2];
                atom3 = atoms[j2];
                p3 = new $3Dmol.Vector3(atom3.x, atom3.y, atom3.z);

                dir2 = p3.clone();
                dir2.sub(p1);

                v = dir2.clone();
                v.cross(dir);
            }

            // especially for C#C (triple bond) dir and dir2
            // may be opposites resulting in a zero v
            if (v.lengthSq() < 0.01) {
                v = dir.clone();
                if (Math.abs(v.x) > 0.0001)
                    v.y += 1;
                else
                    v.x += 1;
            }

            v.cross(dir);
            v.normalize();
            
            return v;
            
            //v.multiplyScalar(r * 1.5);

        }
        
        var getTripleBondPoints = function() {
            
            v.cross(dir);
            v.normalize();
            v.multiplyScalar(r * 3);

            p1a = p1.clone();
            p1a.add(v);
            p1b = p1.clone();
            p1b.sub(v);

            p2a = p1a.clone();
            p2a.add(dir);
            p2b = p1b.clone();
            p2b.add(dir);
        }
        
        var addLine = function(vertexArray, colorArray, offset, p1, p2, c1) {
            //make line from p1 to p2, does not incremeant counts
            vertexArray[offset] = p1.x; vertexArray[offset+1] = p1.y; vertexArray[offset+2] = p1.z;
            colorArray[offset] = c1.r; colorArray[offset+1] = c1.g; colorArray[offset+2] = c1.b;
            vertexArray[offset+3] = p2.x; vertexArray[offset+4] = p2.y; vertexArray[offset+5] = p2.z;
            colorArray[offset+3] = c1.r; colorArray[offset+4] = c1.g; colorArray[offset+5] = c1.b;            
        }
        
        /**@typedef LineStyleSpec
         * @prop {boolean} hidden - do not show line
         * @prop {number} linewidth 
         * @prop {ColorschemeSpec} colorscheme - element based coloring
         * @prop {ColorSpec} color - fixed coloring, overrides colorscheme
         */
        
        // bonds - both atoms must match bond style
        // standardize on only drawing for lowest to highest
        /**
         * 
         * @param {AtomSpec}
         *            atom
         * @param {AtomSpec[]} atoms
         * @param {$3Dmol.Geometry[]} geos
         */
        var drawBondLines = function(atom, atoms, geos) {
            if (!atom.style.line)
                return;
            var style = atom.style.line;
            if (style.hidden)
                return;

            // have a separate geometry for each linewidth
            var linewidth = (style.linewidth || defaultlineWidth);

            if (!geos[linewidth])
                geos[linewidth] = new $3Dmol.Geometry();
            /** @type {geometryGroup} */
            var geoGroup = geos[linewidth].updateGeoGroup(6*atom.bonds.length); //reserve enough space even for triple bonds
            
            var vertexArray = geoGroup.vertexArray;
            var colorArray = geoGroup.colorArray;
            
            for ( var i = 0; i < atom.bonds.length; i++) {
                var j = atom.bonds[i]; // our neighbor
                
                var atom2 = atoms[j];
                if (!atom2.style.line)
                    continue; // don't sweat the details

                if (atom.serial >= atom2.serial) // only draw if less, this way we can do multi bonds correctly
                    continue;
                var p1 = new $3Dmol.Vector3(atom.x, atom.y, atom.z);
                var p2 = new $3Dmol.Vector3(atom2.x, atom2.y, atom2.z);                
                var mp = p1.clone().add(p2).multiplyScalar(0.5);
                var singleBond = false;               
                
                if (atom.clickable){
                    if (atom.intersectionShape === undefined)
                        atom.intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};
                    atom.intersectionShape.line.push(p1);
                    atom.intersectionShape.line.push(p2);
                }

                var c1 = $3Dmol.getColorFromStyle(atom, atom.style.line);
                var c2 = $3Dmol.getColorFromStyle(atom2, atom2.style.line);
               
                if(atom.bondStyles && atom.bondStyles[i]) {
                    var bstyle = atom.bondStyles[i];
                    if(!bstyle.iswire) {
                        continue;
                    }
                    if(bstyle.radius) bondR = bstyle.radius;
                    if(bstyle.singleBond) singleBond = true;
                    if(typeof(bstyle.color1) != "undefined") {
                        c1 = $3Dmol.CC.color(bstyle.color1);
                    }
                    if(typeof(bstyle.color2) != "undefined") {
                        c2 = $3Dmol.CC.color(bstyle.color2);
                    }
                }

                var offset = geoGroup.vertices*3;
                
                if(atom.bondOrder[i] > 1 && atom.bondOrder[i] < 4 && !singleBond) {
                    var v = getSideBondV(atom, atom2, i);
                    var dir = p2.clone();
                    dir.sub(p1);
                    
                    if(atom.bondOrder[i] == 2) { //double
                        
                        v.multiplyScalar(.1);
                           p1a = p1.clone();
                        p1a.add(v);
                        p1b = p1.clone();
                        p1b.sub(v);

                        p2a = p1a.clone();
                        p2a.add(dir);
                        p2b = p1b.clone();
                        p2b.add(dir);
                        
                        if(c1 == c2) {
                            geoGroup.vertices += 4;
                            addLine(vertexArray, colorArray, offset, p1a, p2a, c1);                            
                            addLine(vertexArray, colorArray, offset+6, p1b, p2b, c1);                            
                        }
                        else {
                            geoGroup.vertices += 8;
                            dir.multiplyScalar(0.5);
                            var mpa = p1a.clone();
                            mpa.add(dir);
                            var mpb = p1b.clone();
                            mpb.add(dir);
                            
                            addLine(vertexArray, colorArray, offset, p1a, mpa, c1);                            
                            addLine(vertexArray, colorArray, offset+6, mpa, p2a, c2);                            
                            addLine(vertexArray, colorArray, offset+12, p1b, mpb, c1); 
                            addLine(vertexArray, colorArray, offset+18, mpb, p2b, c2); 
                        }
                    }
                    else if(atom.bondOrder[i] == 3) { //triple
                        
                        v.multiplyScalar(.1);
                           p1a = p1.clone();
                        p1a.add(v);
                        p1b = p1.clone();
                        p1b.sub(v);

                        p2a = p1a.clone();
                        p2a.add(dir);
                        p2b = p1b.clone();
                        p2b.add(dir);
                        
                        if(c1 == c2) {
                            geoGroup.vertices += 6;
                            addLine(vertexArray, colorArray, offset, p1, p2, c1);                            
                            addLine(vertexArray, colorArray, offset+6, p1a, p2a, c1);                            
                            addLine(vertexArray, colorArray, offset+12, p1b, p2b, c1);                            
                        }
                        else {
                            geoGroup.vertices += 12;
                            dir.multiplyScalar(0.5);
                            var mpa = p1a.clone();
                            mpa.add(dir);
                            var mpb = p1b.clone();
                            mpb.add(dir);

                            addLine(vertexArray, colorArray, offset, p1, mp, c1);                            
                            addLine(vertexArray, colorArray, offset+6, mp, p2, c2);
                            addLine(vertexArray, colorArray, offset+12, p1a, mpa, c1);                            
                            addLine(vertexArray, colorArray, offset+18, mpa, p2a, c2);                            
                            addLine(vertexArray, colorArray, offset+24, p1b, mpb, c1); 
                            addLine(vertexArray, colorArray, offset+30, mpb, p2b, c2); 
                        }
                    }
                }
                else { //single bond                                    
                    if(c1 == c2) {
                        geoGroup.vertices += 2;
                        addLine(vertexArray, colorArray, offset, p1, p2, c1);
                    } else {
                        geoGroup.vertices += 4;
                        addLine(vertexArray, colorArray, offset, p1, mp, c1);
                        addLine(vertexArray, colorArray, offset+6, mp, p2, c2);                        
                    }
                    
                }
            }

        };

        // bonds as cylinders
        var defaultStickRadius = 0.25;

        /**@typedef SphereStyleSpec
         * @prop {boolean} hidden - do not show atom
         * @prop {number} radius - override van der waals radius
         * @prop {number} scale - scale radius by specified amount
         * @prop {ColorschemeSpec} colorscheme - element based coloring
         * @prop {ColorSpec} color - fixed coloring, overrides colorscheme
         */
        
        //sphere drawing
        //See also: drawCylinder
        /** 
         * 
         * @param {AtomSpec} atom
         * @param {$3Dmol.Geometry} geo
         */
        var drawAtomSphere = function(atom, geo) {
            
            if (!atom.style.sphere)
                return;
            var style = atom.style.sphere;
            if (style.hidden)
                return;
                                                                 
            var C = $3Dmol.getColorFromStyle(atom, style);
            
            var x, y;
            var radius = getRadiusFromStyle(atom, style);
            
            if ((atom.clickable === true) && (atom.intersectionShape !== undefined)) {
                var center = new $3Dmol.Vector3(atom.x, atom.y, atom.z);
                atom.intersectionShape.sphere.push(new $3Dmol.Sphere(center, radius));
            }
            
            $3Dmol.GLDraw.drawSphere(geo, atom, radius, C);    
            
        };
        
        //dkoes - test code for sphere imposters
        var drawAtomImposter = function(atom, geo) {
            
            if (!atom.style.spherei)
                return;
            var style = atom.style.spherei;
            if (style.hidden)
                return;
            
            var radius = getRadiusFromStyle(atom, style);
            var C = $3Dmol.getColorFromStyle(atom, style);
            
            //create flat square                       
            
            var geoGroup = geo.updateGeoGroup(4);
            var startv =  geoGroup.vertices;
            var start = startv*3;
            var vertexArray = geoGroup.vertexArray;
            var colorArray = geoGroup.colorArray;
            
            //use center point for each vertex
            for(var i = 0; i < 4; i++) {
                vertexArray[start+3*i] = atom.x;
                vertexArray[start+3*i+1] = atom.y ;
                vertexArray[start+3*i+2] = atom.z;                           
            }
            

            //same colors for all 4 vertices
            var normalArray = geoGroup.normalArray;
            var colorArray = geoGroup.colorArray;
            for(var i = 0; i < 4; i++) {
                colorArray[start+3*i] = C.r;
                colorArray[start+3*i+1] = C.g;
                colorArray[start+3*i+2] = C.b;
                
            }
            
            normalArray[start+0] = -radius;
            normalArray[start+1] = radius;
            normalArray[start+2] = 0;
            
            normalArray[start+3] = -radius;
            normalArray[start+4] = -radius;
            normalArray[start+5] = 0;
            
            normalArray[start+6] = radius;
            normalArray[start+7] = -radius;
            normalArray[start+8] = 0;
            
            normalArray[start+9] = radius;
            normalArray[start+10] = radius;
            normalArray[start+11] = 0;
            
            geoGroup.vertices += 4;
            
            //two faces
            var faceArray = geoGroup.faceArray;
            var faceoffset = geoGroup.faceidx; //not number faces, but index
            faceArray[faceoffset+0] = startv;
            faceArray[faceoffset+1] = startv+1;
            faceArray[faceoffset+2] = startv+2;
            faceArray[faceoffset+3] = startv+2;
            faceArray[faceoffset+4] = startv+3;
            faceArray[faceoffset+5] = startv;
            geoGroup.faceidx += 6;
            
        };
                
           
        /**@typedef StickStyleSpec
         * @prop {boolean} hidden - do not show 
         * @prop {number} radius 
         * @prop {boolean} singleBonds - draw all bonds as single bonds if set
         * @prop {ColorschemeSpec} colorscheme - element based coloring
         * @prop {ColorSpec} color - fixed coloring, overrides colorscheme
         */
        
        // draws cylinders and small spheres (at bond radius)
        var drawBondSticks = function(atom, atoms, geo) {
            if (!atom.style.stick)
                return;
            var style = atom.style.stick;
            if (style.hidden)
                return;

            var atomBondR = style.radius || defaultStickRadius;
            var bondR = atomBondR;
            var atomSingleBond = style.singleBonds || false;
            var fromCap = 0, toCap = 0;

            var C1 = $3Dmol.getColorFromStyle(atom, style);

            var mp, mp1, mp2;
            
            if (!atom.capDrawn && atom.bonds.length < 4)
                fromCap = 2;              
                
            for (var i = 0; i < atom.bonds.length; i++) {
                var j = atom.bonds[i]; // our neighbor
                var atom2 = atoms[j]; //parsePDB, etc should only add defined bonds
                
                if (atom.serial < atom2.serial) {// only draw if less, this
                    // lets us combine
                    // cylinders of the same
                    // color
                    var style2 = atom2.style;
                    if (!style2.stick)
                        continue; // don't sweat the details                     
                   
                    var C2 = $3Dmol.getColorFromStyle(atom2, style2.stick);
                    
                    //support bond specific styles
                    bondR = atomBondR;                    
                    var singleBond = atomSingleBond;
                    if(atom.bondStyles && atom.bondStyles[i]) {
                        var bstyle = atom.bondStyles[i];
                        if(bstyle.iswire) {
                            continue;
                        }
                        if(bstyle.radius) bondR = bstyle.radius;
                        if(bstyle.singleBond) singleBond = true;
                        if(typeof(bstyle.color1) != "undefined") {
                            C1 = $3Dmol.CC.color(bstyle.color1);
                        }
                        if(typeof(bstyle.color2) != "undefined") {
                            C2 = $3Dmol.CC.color(bstyle.color2);
                        }
                    }
                    var p1 = new $3Dmol.Vector3(atom.x, atom.y, atom.z);
                    var p2 = new $3Dmol.Vector3(atom2.x, atom2.y, atom2.z);

                    // draw cylinders
                    if (atom.bondOrder[i] === 1 || singleBond) {

                        if (!atom2.capDrawn && atom2.bonds.length < 4)
                            toCap = 2;       
                                                
                        if (C1 != C2) {
                            mp = new $3Dmol.Vector3().addVectors(p1, p2)
                                    .multiplyScalar(0.5);
                            $3Dmol.GLDraw.drawCylinder(geo, p1, mp, bondR, C1, fromCap, 0);
                            $3Dmol.GLDraw.drawCylinder(geo, mp, p2, bondR, C2, 0, toCap);
                        } else {
                            $3Dmol.GLDraw.drawCylinder(geo, p1, p2, bondR, C1, fromCap, toCap);
                        }
                        
                        if (atom.clickable || atom2.clickable) {
                            mp = new $3Dmol.Vector3().addVectors(p1, p2).multiplyScalar(0.5);
                            if (atom.clickable){
                                var cylinder1 = new $3Dmol.Cylinder(p1 , mp , bondR);
                                var sphere1 = new $3Dmol.Sphere(p1 , bondR);
                                atom.intersectionShape.cylinder.push(cylinder1);   
                                atom.intersectionShape.sphere.push(sphere1);                             
                            }
                            if (atom2.clickable){
                                var cylinder2 = new $3Dmol.Cylinder(p2 , mp , bondR);
                                var sphere2 = new $3Dmol.Sphere(p2 , bondR);
                                atom2.intersectionShape.cylinder.push(cylinder2);
                                atom2.intersectionShape.sphere.push(sphere2);
                            }

                        }
                        
                    } 
                    
                    else if (atom.bondOrder[i] > 1) {
                        var mfromCap = 0; mtoCap = 0; //multi bond caps
                        
                        if(bondR != atomBondR) {
                            //assume jmol style multiple bonds - the radius doesn't fit within atom sphere
                            mfromCap = 2;
                            mtoCap = 2;
                        }
                        
                        var dir = p2.clone();
                        var v = null;
                        dir.sub(p1);
                        
                        var r, p1a, p1b, p2a, p2b;
                        var v = getSideBondV(atom, atom2, i);
                        
                        if (atom.bondOrder[i] == 2) {
                            var r = bondR/2.5;
                            var v = getSideBondV(atom, atom2, i);
                            
                            v.multiplyScalar(r*1.5);
                            p1a = p1.clone();
                            p1a.add(v);
                            p1b = p1.clone();
                            p1b.sub(v);

                            p2a = p1a.clone();
                            p2a.add(dir);
                            p2b = p1b.clone();
                            p2b.add(dir);

                                                                 
                            if (C1 != C2) {
                                mp = new $3Dmol.Vector3().addVectors(p1a, p2a)
                                        .multiplyScalar(0.5);
                                mp2 = new $3Dmol.Vector3().addVectors(p1b, p2b)
                                        .multiplyScalar(0.5);
                                $3Dmol.GLDraw.drawCylinder(geo, p1a, mp, r, C1, mfromCap, 0);
                                $3Dmol.GLDraw.drawCylinder(geo, mp, p2a, r, C2, 0, mtoCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1b, mp2, r, C1, mfromCap, 0);
                                $3Dmol.GLDraw.drawCylinder(geo, mp2, p2b, r, C2, 0, mtoCap);
                            } else {
                                $3Dmol.GLDraw.drawCylinder(geo, p1a, p2a, r, C1, mfromCap, mtoCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1b, p2b, r, C1, mfromCap, mtoCap);
                            }
                            if (atom.clickable || atom2.clickable){
                                mp = new $3Dmol.Vector3().addVectors(p1a, p2a)
                                               .multiplyScalar(0.5);
                                mp2 = new $3Dmol.Vector3().addVectors(p1b, p2b)
                                                .multiplyScalar(0.5);
                                if (atom.clickable) {
                                    cylinder1a = new $3Dmol.Cylinder(p1a , mp , r);
                                    cylinder1b = new $3Dmol.Cylinder(p1b , mp2 , r);
                                    atom.intersectionShape.cylinder.push(cylinder1a);
                                    atom.intersectionShape.cylinder.push(cylinder1b);
                                }
                                if (atom2.clickable) {
                                    cylinder2a = new $3Dmol.Cylinder(p2a , mp , r);
                                    cylinder2b = new $3Dmol.Cylinder(p2b , mp2 , r);
                                    atom2.intersectionShape.cylinder.push(cylinder2a);
                                    atom2.intersectionShape.cylinder.push(cylinder2b);                               
                                }
                            }
                        } 
                        else if (atom.bondOrder[i] == 3) {
                            r = bondR / 4;
                            v.cross(dir);
                            v.normalize();
                            v.multiplyScalar(r * 3);

                            p1a = p1.clone();
                            p1a.add(v);
                            p1b = p1.clone();
                            p1b.sub(v);

                            p2a = p1a.clone();
                            p2a.add(dir);
                            p2b = p1b.clone();
                            p2b.add(dir);

                            if (C1 != C2) {
                                mp = new $3Dmol.Vector3().addVectors(p1a, p2a)
                                        .multiplyScalar(0.5);
                                mp2 = new $3Dmol.Vector3().addVectors(p1b, p2b)
                                        .multiplyScalar(0.5);
                                mp3 = new $3Dmol.Vector3().addVectors(p1, p2)
                                        .multiplyScalar(0.5);
                                $3Dmol.GLDraw.drawCylinder(geo, p1a, mp, r, C1, mfromCap, 0);
                                $3Dmol.GLDraw.drawCylinder(geo, mp, p2a, r, C2, 0, mtoCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1, mp3, r, C1, fromCap, 0);
                                $3Dmol.GLDraw.drawCylinder(geo, mp3, p2, r, C2, 0, toCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1b, mp2, r, C1, mfromCap, 0);
                                $3Dmol.GLDraw.drawCylinder(geo, mp2, p2b, r, C2, 0, mtoCap);
                            } else {
                                $3Dmol.GLDraw.drawCylinder(geo, p1a, p2a, r, C1, mfromCap, mtoCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1, p2, r, C1, fromCap, toCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1b, p2b, r, C1, mfromCap, mtoCap);

                            }
                            if (atom.clickable || atom2.clickable) {
                                mp = new $3Dmol.Vector3().addVectors(p1a, p2a)
                                        .multiplyScalar(0.5);
                                mp2 = new $3Dmol.Vector3().addVectors(p1b, p2b)
                                        .multiplyScalar(0.5);
                                mp3 = new $3Dmol.Vector3().addVectors(p1, p2)
                                        .multiplyScalar(0.5);
                                                                
                                if (atom.clickable) {
                                    cylinder1a = new $3Dmol.Cylinder(p1a.clone(), mp.clone(), r);
                                    cylinder1b = new $3Dmol.Cylinder(p1b.clone(), mp2.clone(), r);
                                    cylinder1c = new $3Dmol.Cylinder(p1.clone(), mp3.clone(), r);
                                    atom.intersectionShape.cylinder.push(cylinder1a);
                                    atom.intersectionShape.cylinder.push(cylinder1b);
                                    atom.intersectionShape.cylinder.push(cylinder1c);
                                } 
                                if (atom2.clickable) {                               
                                    cylinder2a = new $3Dmol.Cylinder(p2a.clone(), mp.clone(), r);
                                    cylinder2b = new $3Dmol.Cylinder(p2b.clone(), mp2.clone(), r);
                                    cylinder2c = new $3Dmol.Cylinder(p2.clone(), mp3.clone(), r);
                                    atom2.intersectionShape.cylinder.push(cylinder2a);
                                    atom2.intersectionShape.cylinder.push(cylinder2b);
                                    atom2.intersectionShape.cylinder.push(cylinder2c);                                
                                }
                            }
                        }
                    }
                     
                }                   
                                 
            }            

            // draw non bonded heteroatoms as spheres
            var drawSphere = false;
            var numsinglebonds = 0;
            var differentradii = false;
            //also, if any bonds were drawn as multiples, need sphere
            for(var i = 0; i < atom.bonds.length; i++) {
                var singleBond = atomSingleBond;
                if(atom.bondStyles && atom.bondStyles[i]) {
                    var bstyle = atom.bondStyles[i];
                    if(bstyle.singleBond) singleBond = true;
                    if(bstyle.radius && bstyle.radius != atomBondR) {
                        differentradii = true;
                    }
                }
                if(singleBond || atom.bondOrder[i] == 1) {
                    numsinglebonds++;
                }
            }
            
            if(differentradii) { //jmol style double/triple bonds - no sphere
                if(numsinglebonds > 0) drawSphere = true; //unless needed as a cap
            }
            else if(numsinglebonds == 0 && atom.bonds.length > 0) {
                drawSphere = true;
            }
           
            if (drawSphere) {
                var savedstyle = atom.style;
                bondR = atomBondR;
                //do not use bond style as this can be variable, particularly
                //with jmol export of double/triple bonds
                $3Dmol.GLDraw.drawSphere(geo, atom, bondR, C1);    
            }
            
        };

        // go through all the atoms and regenerate their geometries
        // we try to have one geometry for each style since this is much much
        // faster
        // at some point we should optimize this to avoid unnecessary
        // recalculation
        /** param {AtomSpec[]} atoms */
        var createMolObj = function(atoms) {

            var ret = new $3Dmol.Object3D();
            var cartoonAtoms = [];
            var lineGeometries = {};
            var crossGeometries = {};
            var sphereGeometry = new $3Dmol.Geometry(true);                                                         
            var imposterGeometry = new $3Dmol.Geometry(true);                                                         
            var stickGeometry = new $3Dmol.Geometry(true);
            var i, j, n, testOpacities;
            var opacities = {};
            var range = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
            for (i = 0, n = atoms.length; i < n; i++) {
                var atom = atoms[i];
                // recreate gl info for each atom as necessary
                // set up appropriate intersection spheres for clickable atoms

                if (atom && atom.style) {

                    if (atom.clickable && atom.intersectionShape === undefined)
                        atom.intersectionShape = {sphere: [], cylinder: [], line: [], triangle : []};

                    testOpacities = {line:undefined, cross:undefined, stick:undefined, sphere:undefined};
                    for (j in testOpacities)
                    {
                        if (atom.style[j])
                        {
                            if (atom.style[j].opacity)
                                testOpacities[j] = parseFloat(atom.style[j].opacity);
                            else
                                testOpacities[j] = 1;

                        } else testOpacities[j] = undefined;

                        if (opacities[j])
                        {
                            if (testOpacities[j] != undefined && opacities[j] != testOpacities[j])
                            {
                                console.log("Warning: " + j + " opacity is ambiguous");
                                opacities[j] = 1;
                            }

                        } else opacities[j] = testOpacities[j];
                    }
                    
                    drawAtomSphere(atom, sphereGeometry);
                    drawAtomImposter(atom, imposterGeometry);
                    drawAtomCross(atom, crossGeometries);
                    drawBondLines(atom, atoms, lineGeometries);
                    drawBondSticks(atom, atoms, stickGeometry);

                    if (typeof (atom.style.cartoon) !== "undefined" && !atom.style.cartoon.hidden) {
                        //gradient color scheme range
                        if (atom.style.cartoon.color === "spectrum" && typeof(atom.resi) === "number" && !atom.hetflag) {                            
                            if (atom.resi < range[0])
                                range[0] = atom.resi;
                            if (atom.resi > range[1])
                                range[1] = atom.resi;
                        }

                        cartoonAtoms.push(atom);
                    }                   
                }
            }
            // create cartoon if needed - this is a whole model analysis
            if (cartoonAtoms.length > 0) {
                var gradientscheme = null;
                //TODO: Should have an option to choose gradient type
                if (range[0] < range[1])
                    gradientscheme = new $3Dmol.Gradient.Sinebow(range[0], range[1]);

                $3Dmol.drawCartoon(ret, cartoonAtoms, gradientscheme);

            }

            // add sphere geometry
            if (sphereGeometry.vertices > 0) {
                var sphereMaterial = new $3Dmol.MeshLambertMaterial({
                    ambient : 0x000000,
                    vertexColors : true,
                    reflectivity : 0,
                });
                if (opacities.sphere < 1 && opacities.sphere >= 0)
                {
                    sphereMaterial.transparent = true;
                    sphereMaterial.opacity = opacities.sphere;
                }
                
                //Initialize buffers in geometry                
                sphereGeometry.initTypedArrays();   
				
                var sphere = new $3Dmol.Mesh(sphereGeometry, sphereMaterial);
                ret.add(sphere);
            }

            // add imposter geometry
            if (imposterGeometry.vertices > 0) {
                var imposterMaterial = new $3Dmol.ImposterMaterial({
                    ambient : 0x000000,
                    vertexColors : true,
                    reflectivity : 0
                });
                
                //Initialize buffers in geometry                
                imposterGeometry.initTypedArrays();
                
                var spherei = new $3Dmol.Mesh(imposterGeometry, imposterMaterial);
                console.log("spherei geometry " + imposterGeometry.vertices.length);

                ret.add(spherei);
            }
            
            // add stick geometry
            if (stickGeometry.vertices > 0) {
                var cylinderMaterial = new $3Dmol.MeshLambertMaterial({
                    vertexColors : true,
                    ambient : 0x000000,
                    reflectivity : 0
                });
                if (opacities.stick < 1 && opacities.stick >= 0)
                {
                    cylinderMaterial.transparent = true;
                    cylinderMaterial.opacity = opacities.stick;
                }

                //Initialize buffers in geometry                
                stickGeometry.initTypedArrays();
                
                if (cylinderMaterial.wireframe)
                    stickGeometry.setUpWireframe();
            
                var sticks = new $3Dmol.Mesh(stickGeometry, cylinderMaterial);
                ret.add(sticks);
            }
            
            //var linewidth;
            // add any line geometries, distinguished by line width
            for (i in lineGeometries) {
                if (lineGeometries.hasOwnProperty(i)) {
                    var linewidth = i;
                    var lineMaterial = new $3Dmol.LineBasicMaterial({
                        linewidth : linewidth,
                        vertexColors : true
                    });
                    if (opacities.line < 1 && opacities.line >= 0)
                    {
                        lineMaterial.transparent = true;
                        lineMaterial.opacity = opacities.line;
                    }
                    
                    lineGeometries[i].initTypedArrays();
                    
                    var line = new $3Dmol.Line(lineGeometries[i], lineMaterial,
                            $3Dmol.LinePieces);

                    ret.add(line);
                }
            }

            // add any cross geometries
            for (i in crossGeometries) {
                if (crossGeometries.hasOwnProperty(i)) {
                    var linewidth = i;
                    var crossMaterial = new $3Dmol.LineBasicMaterial({
                        linewidth : linewidth,
                        vertexColors : true
                    });
                    if (opacities.cross < 1 && opacities.cross >= 0)
                    {
                        crossMaterial.transparent = true;
                        crossMaterial.opacity = opacities.cross;
                    }

                    crossGeometries[i].initTypedArrays();
                    
                    var cross = new $3Dmol.Line(crossGeometries[i], crossMaterial,
                            $3Dmol.LinePieces);

                    ret.add(cross);
                }
            }
            
            //for BIOMT assembly
            if (dontDuplicateAtoms && modelData.symmetries && modelData.symmetries.length > 0) {
                var finalRet = new $3Dmol.Object3D();
                var t;
                for (t = 0; t < modelData.symmetries.length; t++) {
                    var transformedRet = new $3Dmol.Object3D();
                    transformedRet = ret.clone();
                    transformedRet.matrix.copy(modelData.symmetries[t]);
                    transformedRet.matrixAutoUpdate = false;
                    finalRet.add(transformedRet);
                }
                return finalRet;
            }

            return ret;
        };
        
        
        this.getCrystData = function() {
            if (modelData.cryst) {
                return modelData.cryst;
            }
            else {
                return null;
            }
        }
        
        /**
         * Returns list of rotational/translational matrices if there is BIOMT data
         * Otherwise returns a list of just the ID matrix
         *
         * @function $3Dmol.GlModel#getSymmetries
         * @return {Array<$3Dmol.Matrix4>}
         *
         */
        this.getSymmetries = function() {
            
            if (typeof(modelData.symmetries) == 'undefined') {
                modelData.symmetries = [idMatrix];
            }
            return modelData.symmetries; 
        };
        
        /**
         * Sets symmetries based on specified matrices in list
         *
         * @function $3Dmol.GlModel#setSymmetries
         * @param {Array<$3Dmol.Matrix4>} list
         *
         */
        this.setSymmetries = function(list) {
            if (typeof(list) == "undefined") { //delete sym data
                modelData.symmetries = [idMatrix];
            }
            else {
                modelData.symmetries = list;
            }
        };

        /**
         * Returns model id number
         * 
         * @function $3Dmol.GLModel#getID
         * @return {number} Model ID
         */
        this.getID = function() {
            return id;
        };
        
        /**
         * Returns model's frames property, a list of atom lists
         * 
         * @function $3Dmol.GLModel#getFrames
         * @return {Array.<Object>}
         */
        this.getFrames = function() {
            return frames;
        };
        
        /**
         * Sets model's atomlist to specified frame
         * Sets to last frame if framenum out of range
         * 
         * @function $3Dmol.GLModel#setFrame
         * @param {number} framenum - model's atoms are set to this index in frames list
         */
        this.setFrame = function(framenum) {
            if (frames.length == 0) {
                return;
            }
            if (framenum >= 0 && framenum < frames.length) {
                atoms = frames[framenum];
            }
            else {
                atoms = frames[frames.length-1];
            }
            molObj = null;
        };
        
        /**
         * Add atoms as frames of model
         * 
         * @function $3Dmol.GLModel#addFrame
         * @param {AtomSpec} atom - atoms to be added
         */
        this.addFrame = function(atoms) {
            frames.push(atoms);
        };

        
        /**
         * If model atoms have dx, dy, dz properties (in some xyz files), vibrate populates the model's frame property based on parameters.
         * Model can then be animated
         * 
         * @function $3Dmol.GLModel#vibrate
         * @param {number} numFrames - number of frames to be created, default to 10
         * @param {number} amplitude - amplitude of distortion, default to 1 (full)
         */
        this.vibrate = function(numFrames, amplitude) {
            var vectors = [];
            var amplitude = amplitude || 1;
            var numFrames = numFrames || 10; 
            for (var i = 0; i < atoms.length; i++) {
                var vector = new $3Dmol.Vector3(atoms[i].dx, atoms[i].dy, atoms[i].dz);
                vectors.push(vector);
            }
            numFrames--;
            for (var i = 1; i <= numFrames; i++) {
                var newAtoms = [];
                for (var j = 0; j < atoms.length; j++) {
                    var newVector = new $3Dmol.Vector3(atoms[j].dx, atoms[j].dy, atoms[j].dz);
                    var starting = new $3Dmol.Vector3(atoms[j].x, atoms[j].y, atoms[j].z);
                    newVector.sub(starting);
                    newVector.multiplyScalar((i*amplitude)/numFrames);
                    starting.add(newVector);
                    var newAtom = {};
                    for (var k in atoms[j]) {
                        newAtom[k] = atoms[j][k];
                    }
                    newAtom.x = starting.x;
                    newAtom.y = starting.y;
                    newAtom.z = starting.z;
                    newAtoms.push(newAtom);
                }
                frames.push(newAtoms);
            }
            frames.unshift(atoms); //add 1st frame
        };
        
        // set default style and colors for atoms
        this.setAtomDefaults = function(atoms) {
            for ( var i = 0; i < atoms.length; i++) {
                var atom = atoms[i];
                if (atom) {
                    atom.style = atom.style || defaultAtomStyle;
                    atom.color = atom.color || ElementColors[atom.elem] || defaultColor;
                    atom.model = id;
                    if (atom.clickable)
                        atom.intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};
                }
            }
        };

        /** add atoms to this model from molecular data string
         * 
         * @function $3Dmol.GLModel#addMolData
         * @param {string|ArrayBuffer} data - atom structure file input data string, for gzipped input use ArrayBuffer
         * @param {string} format - input file string format (e.g 'pdb', 'sdf', 'sdf.gz', etc.)
         * @param {Object} options - format dependent options (e.g. 'options.keepH' to keep hydrogens)
         */
        this.addMolData = function(data, format, options) {
            options = options || {};
            var parsedAtoms = $3Dmol.GLModel.parseMolData(data, format, options);
            dontDuplicateAtoms = !options.duplicateAssemblyAtoms;
            var mData = parsedAtoms.modelData;
            if (mData) {
                if (Array.isArray(mData)) {
                    modelData = mData[0];
                } else {
                    modelData = mData;
                }
            }

            if (frames.length == 0) { //first call
                for (var i = 0; i < parsedAtoms.length; i++) {
                    if (parsedAtoms[i].length != 0)
                        frames.push(parsedAtoms[i]);
                }
                if(frames[0])
                    atoms = frames[0];
            }
            
            else { //subsequent calls
                if (options.frames) { //add to new frame
                    for (var i = 0; i < parsedAtoms.length; i++) {
                        frames.push(parsedAtoms[i]);
                    }
                }
                else { //add atoms to current frame
                    for (var i = 0; i < parsedAtoms.length; i++) {
                        this.addAtoms(parsedAtoms[i]); 
                    }
                }
            }
            
            for (var i = 0; i < frames.length; i++) {
                this.setAtomDefaults(frames[i], id);
            }

        };

        this.setDontDuplicateAtoms = function(dup) {
            dontDuplicateAtoms = dup;
        }

        this.setModelData = function(mData) {
            modelData = mData;
        }
        
        /** given a selection specification, return true if atom is selected
         * 
         * @function $3Dmol.GLModel#atomIsSelected
         * @param {AtomSpec} atom
         * @param {AtomSelectionSpec} sel
         * @return {boolean}
         */
        this.atomIsSelected = function(atom, sel) {
            if (typeof (sel) === "undefined")
                return true; // undef gets all
            var invert = !!sel.invert;
            var ret = true;
            for ( var key in sel) {
                if(key === 'predicate') { //a user supplied function for evaluating atoms
                    if(!sel['predicate'](atom) ) {
                        ret = false;
                        break;
                    }
                }
                else if(key == "properties" && atom[key]) {
                    for (var propkey in sel.properties) {
                        if(typeof(atom.properties[propkey]) === 'undefined') {
                            ret = false;
                            break
                        }
                        if(atom.properties[propkey] != sel.properties[propkey]) {
                            ret = false;
                            break;
                        }
                    }
                }
                else if (sel.hasOwnProperty(key) && key != "props" && key != "invert" && key != "model" && key != "byres" && key != "expand" && key != "within") {

                    // if something is in sel, atom must have it                    
                    if (typeof (atom[key]) === "undefined") {
                        ret = false;
                        break;
                    }
                    var isokay = false;
                    if(key === "bonds") {
                        //special case counting number of bonds, for selecting nonbonded mostly
                        var val = sel[key];
                        if(val != atom.bonds.length) {
                            ret = false;
                            break;
                        }
                    }
                    else if ($.isArray(sel[key])) {
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


        /** return list of atoms selected by sel, this is specific to glmodel
         * 
         * @function $3Dmol.GLModel#selectedAtoms
         * @param {AtomSelectionSpec} sel
         * @return {Array.<Object>}
         */
        this.selectedAtoms = function(sel, from) {
            var ret = [];
            sel = sel || {};
            if (!from) from = atoms;
            var aLength = from.length;
            for ( var i = 0; i < aLength; i++) {
                var atom = from[i];
                if (atom) {
                    if (this.atomIsSelected(atom, sel))
                        ret.push(atom);
                }
            }

            // expand selection by some distance
            if (sel.hasOwnProperty("expand")) {

                // get atoms in expanded bounding box

                var expand = expandAtomList(ret, parseFloat(sel.expand));
                var retlen = ret.length;
                for (var i = 0; i < expand.length; i++) {
                    for (var j = 0; j < retlen; j++) {

                        var dist = squaredDistance(expand[i], ret[j]);
                        var thresh = Math.pow(sel.expand, 2);
                        if (dist < thresh && dist > 0) {
                            ret.push(expand[i]);
                        }
                    }
                }
            }

            // selection within distance of sub-selection
            if (sel.hasOwnProperty("within") && sel.within.hasOwnProperty("sel") && sel.within.hasOwnProperty("distance")) {

                // get atoms in second selection
                var sel2 = this.selectedAtoms(sel.within.sel, atoms);
                var within = [];
                for (var i = 0; i < sel2.length; i++) {
                    for (var j = 0; j < ret.length; j++) {

                        var dist = squaredDistance(sel2[i], ret[j]);
                        var thresh = Math.pow(parseFloat(sel.within.distance), 2);
                        if (dist < thresh && dist > 0) {
                            within.push(ret[j]);
                        }
                    }
                }
                ret = within;
            }

            // byres selection flag
            if (sel.hasOwnProperty("byres")) {

                // Keep track of visited residues, visited atoms, and atom stack
                var vResis = {};
                var vAtoms = [];
                var stack = [];

                for (var i = 0; i < ret.length; i++) {
                    
                    // Check if atom is part of a residue, and that the residue hasn't been traversed yet
                    var atom = ret[i];
                    var c = atom.chain;
                    var r = atom.resi;
                    if (vResis[c] === undefined) vResis[c] = {};
                    if (atom.hasOwnProperty("resi") && vResis[c][r] === undefined) {

                        // Perform a depth-first search of atoms with the same resi
                        vResis[c][r] = true;
                        stack.push(atom);
                        while(stack.length > 0) {
                            atom = stack.pop();
                            c = atom.chain;
                            r = atom.resi;
                            if (vAtoms[atom.index] === undefined) {
                                vAtoms[atom.index] = true;
                                for (var j = 0; j < atom.bonds.length; j++) {
                                    var atom2 = atoms[atom.bonds[j]];
                                    if (vAtoms[atom2.index] === undefined && atom2.hasOwnProperty("resi") && atom2.chain == c && atom2.resi == r) {
                                        stack.push(atom2);
                                        ret.push(atom2);
                                    }
                                }
                            }
                        }
                    }   
                }
            }

            return ret;
        };

        var squaredDistance = function(atom1, atom2) {
            var xd = atom2.x - atom1.x;
            var yd = atom2.y - atom1.y;
            var zd = atom2.z - atom1.z;
            return (Math.pow(xd, 2) + Math.pow(yd, 2) + Math.pow(zd, 2));
        };

        /** returns a list of atoms in the expanded bounding box, but not in the current one
         *
         *  Bounding box:
         *
         *    [ [ xmin, ymin, zmin ],
         *      [ xmax, ymax, zmax ],
         *      [ xctr, yctr, zctr ] ]
         *
         **/
        var expandAtomList = function(atomList, amt) {

            if (amt <= 0) return atomList;

            var pb = $3Dmol.getExtent(atomList); // previous bounding box
            var nb = [[], [], []]; // expanded bounding box

            for (var i = 0; i < 3; i++)
            {
                nb[0][i] = pb[0][i]-amt;
                nb[1][i] = pb[1][i]+amt;
                nb[2][i] = pb[2][i];
            }

            // look in added box "shell" for new atoms
            var expand = [];
            for (var i = 0; i < atoms.length; i++) {

                var x = atoms[i].x;
                var y = atoms[i].y;
                var z = atoms[i].z;

                if (x >= nb[0][0] && x <= nb[1][0] && y >= nb[0][1] && y <= nb[1][1] && z >= nb[0][2] && z <= nb[1][2]) {
                    if (!(x >= pb[0][0] && x <= pb[1][0] && y >= pb[0][1] && y <= pb[1][1] && z >= pb[0][2] && z <= pb[1][2])) {
                        expand.push(atoms[i]);
                    }
                }
            }
            return expand;
        };
        
        /** Add list of new atoms to model.  Adjusts bonds appropriately.
         * 
         * @function $3Dmol.GLModel#addAtoms
         * @param {type} newatoms
         */        
        this.addAtoms = function(newatoms) {
            molObj = null;
            var start = atoms.length;
            var indexmap = [];
            // mapping from old index to new index
            var i;
            for(i = 0; i < newatoms.length; i++) {
                if(typeof(newatoms[i].index) == "undefined")
                    newatoms[i].index = i;
                if(typeof(newatoms[i].serial) == "undefined")
                    newatoms[i].serial = i;
                indexmap[newatoms[i].index] = start+i;
            }
            
            // copy and push newatoms onto atoms
            for(i = 0; i < newatoms.length; i++) {
                var olda = newatoms[i];
                var nindex = indexmap[olda.index];
                var a = $.extend(false, {}, olda);
                a.index = nindex;
                a.bonds = [];
                a.bondOrder = [];
                a.model = id;
                a.style = a.style || defaultAtomStyle;
                if(typeof(a.color) == "undefined")
                    a.color = ElementColors[a.elem] || defaultColor;                
                // copy over all bonds contained in selection,
                // updating indices appropriately
                var nbonds = olda.bonds ? olda.bonds.length : 0;
                for(var j = 0; j < nbonds; j++) {
                    var neigh = indexmap[olda.bonds[j]];
                    if(typeof(neigh) != "undefined") {
                        a.bonds.push(neigh);
                        a.bondOrder.push(olda.bondOrder ? olda.bondOrder[j] : 1);
                    }                
                }
                atoms.push(a);
            }
        };

        /** Remove specified atoms from model
         * 
         * @function $3Dmol.GLModel#removeAtoms
         * @param {type} badatoms
         * @return {removeAtoms}
         */
        this.removeAtoms = function(badatoms) {
            molObj = null;
            // make map of all baddies
            var baddies = [];
            var i;
            for(i = 0; i < badatoms.length; i++) {
                baddies[badatoms[i].index] = true;
            }
            
            // create list of good atoms
            var newatoms = [];
            for(i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                if(!baddies[a.index])
                    newatoms.push(a);
            }
            
            // clear it all out
            atoms = [];
            // and add back in to get updated bonds
            this.addAtoms(newatoms);
        };
        
        
        /** Set atom style of selected atoms
         * 
         * @function $3Dmol.GLModel#setStyle
         * @param {AtomSelectionSpec} sel
         * @param {AtomStyleSpec} style
         * @param {boolean} add - if true, add to current style, don't replace
         */
        this.setStyle = function(sel, style, add) {
            
            if(typeof(style) === 'undefined' && typeof(add) == 'undefined') {
                //if a single argument is provided, assume it is a style and select all
                style = sel;
                sel = {};
            }
            
            // report to console if this is not a valid selector
            var s;
            for (s in sel) {
                if(validAtomSelectionSpecs.indexOf(s) === -1) {
                    console.log('Unknown selector ' + s);
                }
            }
            // report to console if this is not a valid style
            for (s in style) {
                if(validAtomStyleSpecs.indexOf(s) === -1) {
                    console.log('Unknown style ' + s);
                }
            }

            var changedAtoms = false;
            // somethings we only calculate if there is a change in a certain
            // style, although these checks will only catch cases where both
            // are either null or undefined
            
            var setStyleHelper = function(atomArr) {
                var selected = that.selectedAtoms(sel, atomArr);
                for (var i = 0; i < atomArr.length; i++) {
                    if (atomArr[i]) atomArr[i].capDrawn = false; //reset for proper stick render
                }
            
                for ( var i = 0; i < selected.length; i++) {                
                    changedAtoms = true;
                    if (selected[i].clickable) 
                        selected[i].intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};                    
                   

                    if(!add) selected[i].style = {};
                    for(s in style) {
                        if(style.hasOwnProperty(s)) {
                            selected[i].style[s]=selected[i].style[s]||{}; //create distinct object for each atom
                            for(var prop in style[s]){
                                selected[i].style[s][prop]=style[s][prop];
                            }
                        }
                    }
                }
            }
            
            var that = this;
            setStyleHelper(atoms);
            for (var i = 0; i < frames.length; i++) {
                setStyleHelper(frames[i]);
            }
            
            if (changedAtoms)
                molObj = null; //force rebuild
            
        };

        /** Set clickable and callback of selected atoms
         * 
         * @function $3Dmol.GLModel#setClickable
         * @param {AtomSelectionSpec} sel - atom selection to apply clickable settings to
         * @param {boolean} clickable - whether click-handling is enabled for the selection
         * @param {function} callback - function called when an atom in the selection is clicked
         */
        this.setClickable = function(sel, clickable, callback) {           

            // report to console if this is not a valid selector
            var s;
            for (s in sel) {
                if (validAtomSelectionSpecs.indexOf(s) === -1) {
                    console.log('Unknown selector ' + s);
                }
            }

            // make sure clickable is a boolean
            clickable = !!clickable;

            // report to console if callback is not a valid function
            if (callback && typeof callback != "function") {
                console.log("Callback is not a function");
                return;
            }

            var i;
            var selected = this.selectedAtoms(sel, atoms);
            var len = selected.length;
            for (i = 0; i < len; i++) {                

                selected[i].intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};
                selected[i].clickable = clickable;
                if (callback) selected[i].callback = callback;

            }

            if (len > 0) molObj = null; // force rebuild to get correct intersection shapes         
        };
        
        /** given a mapping from element to color, set atom colors
         * 
         * @function $3Dmol.GLModel#setColorByElement
         * @param {type} sel
         * @param {type} colors
         */
        this.setColorByElement = function(sel, colors) {
            
            if(molObj !== null && sameObj(colors,lastColors))
                return; // don't recompute
            lastColors = colors;
            var atoms = this.selectedAtoms(sel, atoms);
            if(atoms.length > 0)
                molObj = null; // force rebuild
            for ( var i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                if(typeof(colors[a.elem]) !== "undefined") {
                    a.color = colors[a.elem];
                }
            }
        };
        
        /**
         * @function $3Dmol.GLModel.setColorByProperty
         * @param {type} sel
         * @param {type} prop
         * @param {type} scheme
         */
        this.setColorByProperty = function(sel, prop, scheme, range) {
            var atoms = this.selectedAtoms(sel, atoms);
            lastColors = null; // don't bother memoizing
            if(atoms.length > 0)
                molObj = null; // force rebuild
            var min =  Number.POSITIVE_INFINITY;
            var max =  Number.NEGATIVE_INFINITY;
            var i, a;
            
            if(!range) { //no explicit range, get from scheme
                range = scheme.range();
            }
            
            if(!range) { //no range in scheme, compute the range for this model
                range = $3Dmol.getPropertyRange(atoms, prop);
            }
            // now apply colors using scheme
            for (i = 0; i < atoms.length; i++) {
                a = atoms[i];
                var val = $3Dmol.getAtomProperty(a, prop);
                if(val != null) {
                    a.color = scheme.valueToHex(parseFloat(a.properties[prop]), [range[0],range[1]]);
                }                    
            }
        };
        
        /**
         * @function $3Dmol.GLModel.setColorByFunction
         * @param {type} sel
         * @param {type} func
         */
        this.setColorByFunction = function(sel, colorfun) {
            var atoms = this.selectedAtoms(sel, atoms);
            lastColors = null; // don't bother memoizing
            if(atoms.length > 0)
                molObj = null; // force rebuild
            
            // now apply colorfun
            for (i = 0; i < atoms.length; i++) {
                a = atoms[i];
                a.color = colorfun(a);
            }
        };

        /** Convert the model into an object in the format of a ChemDoodle JSON model.
         *
         * @function $3Dmol.GLModel#toCDObject
         * @param {boolean} whether or not to include style information. Defaults to false.
         * @return {Object}
         */
        this.toCDObject = function(includeStyles) {
            var out = { a:[], b:[] };
            if (includeStyles) {
                out.s = [];
            }
            for (var i = 0; i < atoms.length; i++) {
                var atomJSON = {};
                var atom = atoms[i];
                atomJSON.x = atom.x;
                atomJSON.y = atom.y;
                atomJSON.z = atom.z;
                if (atom.elem != "C") {
                    atomJSON.l = atom.elem;
                }
                if (includeStyles) {
                    var s = 0;
                    while (s < out.s.length &&
                          (JSON.stringify(atom.style) !== JSON.stringify(out.s[s]))) {
                        s++;
                    }
                    if (s === out.s.length) {
                        out.s.push(atom.style);
                    }
                    if (s !== 0) {
                        atomJSON.s = s;
                    }
                }
                
                out.a.push(atomJSON);

                for (var b = 0; b < atom.bonds.length; b++) {
                    var firstAtom = i;
                    var secondAtom = atom.bonds[b];
                    if (firstAtom >= secondAtom)
                        continue;
                    var bond = {
                        b: firstAtom,
                        e: secondAtom
                    };
                    var bondOrder =  atom.bondOrder[b];
                    if (bondOrder != 1) {
                        bond.o = bondOrder;
                    }
                    out.b.push(bond);
                }
            }
            return out;
        }


        /** manage the globj for this model in the possed modelGroup - if it has to be regenerated, remove and add
         * 
         * @function $3Dmol.GLModel#globj
         * @param {$3Dmol.Object3D} group
         */
        this.globj = function(group) {
            var time = new Date();
            if(molObj === null) { // have to regenerate
                molObj = createMolObj(atoms);
                var time2 = new Date();
                //console.log("object creation time: " + (time2 - time));
                if(renderedMolObj) { // previously rendered, remove
                    group.remove(renderedMolObj);
                    renderedMolObj = null;
                }
                renderedMolObj = molObj.clone();
                if(hidden) {
                    renderedMolObj.setVisible(false);
                    molObj.setVisible(false);
                }
                group.add(renderedMolObj);              
            }
        };
        
        /** Remove any renderable mol object from scene
         * 
         * @function $3Dmol.GLModel#removegl
         * @param {$3Dmol.Object3D} group
         */
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
        
        /** Don't show this model is future renderings.  Keep all styles and state
         * so it can be efficiencly shown again.
         * 
         * @function $3Dmol.GLModel#hide
         */
        this.hide = function() {
            hidden = true;
            if(renderedMolObj) renderedMolObj.setVisible(false);
            if(molObj) molObj.setVisible(false);
        }
        
        this.show = function() {
            hidden = false;
            if(renderedMolObj) renderedMolObj.setVisible(true);
            if(molObj) molObj.setVisible(true);
        }
        
        /** Create labels for residues of selected atoms.
         * Will create a single label at the center of mass of all atoms
         * with the same chain,resn, and resi.
         * @function $3Dmol.GLModel#addResLabels
         * 
         * @param {AtomSelectionSpec} sel
         * @param {$3Dmol.GLViewer} viewer
         */
        this.addResLabels = function(sel, viewer, style) {
            var atoms = this.selectedAtoms(sel, atoms);
            var bylabel = {}
            //collect by chain:resn:resi
            for(var i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                var c = a.chain;
                var resn = a.resn;
                var resi = a.resi;
                var label =  resn + '' + resi;
                if(!bylabel[c]) bylabel[c] = {};
                if(!bylabel[c][label]) bylabel[c][label] = []
                bylabel[c][label].push(a);
            }
            
            var mystyle = $.extend(true, {}, style);
            //now compute centers of mass
            for(var c in bylabel) {
                if(bylabel.hasOwnProperty(c)) {
                    var labels = bylabel[c];
                    for(var label in labels) {
                        if(labels.hasOwnProperty(label)) {
                            var atoms = labels[label];
                            var sum = new $3Dmol.Vector3(0,0,0);
                            for(var i = 0; i < atoms.length; i++) {
                                var a = atoms[i];
                                sum.x += a.x;
                                sum.y += a.y;
                                sum.z += a.z;
                            }
                            sum.divideScalar(atoms.length);
                            mystyle.position = sum;
                            viewer.addLabel(label, mystyle);
                        }                        
                    }
                }
            }
        }

    }

    GLModel.parseMolData = function(data, format, options) {
        format = format || "";

        if (!data)
            return []; //leave an empty model

        if(/\.gz$/.test(format)) {
            //unzip gzipped files
            format = format.replace(/\.gz$/,'');
            try {
                data = pako.inflate(data, {to: 'string'});
            } catch(err) {
                console.log(err);
            }
        }

        if(typeof($3Dmol.Parsers[format]) == "undefined") {
            //let someone provide a file name and get format from extension
            format = format.split('.').pop();
            if(typeof($3Dmol.Parsers[format]) == "undefined") {
                console.log("Unknown format: "+format);
                //try to guess correct format from data contents
                if(data.match(/^@<TRIPOS>MOLECULE/gm)) {
                    format = "mol2";
                } else if(data.match(/^HETATM/gm) || data.match(/^ATOM/gm)) {
                    format = "pdb";
                } else if(data.match(/^.*\n.*\n.\s*(\d+)\s+(\d+)/gm)){
                    format = "sdf"; //could look at line 3
                } else {
                    format = "xyz";
                }
                console.log("Best guess: "+format);
            }
        }
        var parse = $3Dmol.Parsers[format];
        var parsedAtoms = parse(data, options);

        return parsedAtoms;
    };


    // set default style and colors for atoms
    GLModel.setAtomDefaults = function(atoms, id) {
        for ( var i = 0; i < atoms.length; i++) {
            var atom = atoms[i];
            if (atom) {
                atom.style = atom.style || defaultAtomStyle;
                atom.color = atom.color || ElementColors[atom.elem] || defaultColor;
                atom.model = id;
                if (atom.clickable)
                    atom.intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};
            }
        }
    };

    return GLModel;
    
})();
