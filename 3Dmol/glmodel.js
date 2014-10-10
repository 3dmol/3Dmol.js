// A model is a collection of related atoms.  Bonds are only allowed between
//atoms in the same model.  An atom is uniquely specified by its model id and
//its serial number.
//A glmodel knows how to apply the styles on each atom to create a gl object

var $3Dmol = $3Dmol || {};

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

        if (distSquared > 1.3 && (atom1.elem == 'H' || atom2.elem == 'H' || atom1.elem == 'D' || atom2.elem == 'D'))
            return false;
        if (distSquared < 3.6 && (atom1.elem == 'S' || atom2.elem == 'S'))
            return true;
        if (distSquared > 2.78)
            return false;
        return true;
    };

    /** @param {Array.<AtomSpec>} atomsarray */
    var assignBonds = function(atomsarray) {
        // assign bonds - yuck, can't count on connect records
        var atoms = atomsarray.slice(0);
        var i, j, n;
        for (i = 0, n = atomsarray.length; i < n; i++)
        {
            //Don't reindex if atoms are already indexed 
            if (!atomsarray[i].index)
                atomsarray[i].index = i;
        }
        
        atoms.sort(function(a, b) {
            return a.z - b.z;
        });
        for (i = 0, n = atoms.length; i < n; i++) {
            var ai = atoms[i];

            for (j = i + 1; j < n; j++) {
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
    /** @param {Array.<AtomSpec>} atomsarray */
    var assignPDBBonds = function(atomsarray) {
        // assign bonds - yuck, can't count on connect records
        var protatoms = [];
        var hetatoms = [];        
        var i, n;
        for (i = 0, n = atomsarray.length; i < n; i++) {
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
        
        for (i = 0, n = protatoms.length; i < n; i++) {
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

    // this will identify all hydrogen bonds between backbone
    // atoms; assume atom names are correct, only identifies
    // single closest hbond
    var assignBackboneHBonds = function(atomsarray) {
	var maxlength = 3.2;
	var maxlengthSq = 10.24;
        var atoms = [];
        var i, j, n;
        for (i = 0, n = atomsarray.length; i < n; i++) {
            atomsarray[i].index = i;
            // only consider 'N' and 'O'
            var atom = atomsarray[i];
            if (!atom.hetflag && (atom.atom === "N" || atom.atom === "O")) {
                atoms.push(atom);
                atom.hbondOther = null;
                atom.hbondDistanceSq = Number.POSITIVE_INFINITY;                
            }
        }

        atoms.sort(function(a, b) {
            return a.z - b.z;
        });
        for (i = 0, n = atoms.length; i < n; i++) {
            var ai = atoms[i];

            for (j = i + 1; j < n; j++) {
                var aj = atoms[j];
		var zdiff = aj.z - ai.z;
                if (zdiff > maxlength) // can't be connected
                    break;
		if (aj.atom == ai.atom)
		    continue; //can't be connected, but later might be	
		var ydiff = Math.abs(aj.y - ai.y);
		if( ydiff > maxlength)
		    continue;
		var xdiff = Math.abs(aj.x - ai.x);
		if(xdiff > maxlength)
		    continue;
                var dist = xdiff*xdiff+ydiff*ydiff+zdiff*zdiff;
		if (dist >  maxlengthSq)
		    continue;

		if(aj.chain == ai.chain && Math.abs(aj.resi - ai.resi) < 4)
		    continue; //ignore bonds between too close residues
		//select closest hbond
                if (dist < ai.hbondDistanceSq) {
                    ai.hbondOther = aj;
                    ai.hbondDistanceSq = dist;
                }
                if(dist < aj.hbondDistanceSq) {
                    aj.hbondOther = ai;
                    aj.hbondDistanceSq = dist;
                }
            }
        }
    };

    var computeSecondaryStructure = function(atomsarray) {
        assignBackboneHBonds(atomsarray);
        
        // compute, per residue, what the secondary structure is
        var chres = {}; // lookup by chain and resid
        var i, il, c, r;
        var atom, val;
        
        for (i = 0, il = atomsarray.length; i < il; i++) {
            atom = atomsarray[i];
            
            if (typeof(chres[atom.chain]) === "undefined")
                chres[atom.chain] = [];
            
            if (isFinite(atom.hbondDistanceSq)) {
                var other = atom.hbondOther;
                if (Math.abs(other.resi - atom.resi) === 4) { 
                    // helix
                    chres[atom.chain][atom.resi] = 'h';
                }
                else { // otherwise assume sheet
                    chres[atom.chain][atom.resi] = 's';
                }
            }
        }
        
        // plug gaps and remove singletons
        for (c in chres) {
            for (r = 1; r < chres[c].length-1; r++) {
                var valbefore = chres[c][r-1];
                var valafter = chres[c][r+1];
                val = chres[c][r];
                if(valbefore == valafter && val != valbefore) {
                    chres[c][r] = valbefore;
                }
            }
            for (r = 0; r < chres[c].length; r++) {
                val = chres[c][r];
                if (val == 'h' || val == 's') {
                    if (chres[c][r-1] != val && chres[c][r+1] != val)
                        delete chres[c][r];
                }
            }
        }
        
        // assign to all atoms in residue, keep track of start
        var curres = null;
        for (i = 0, il = atomsarray.length; i < il; i++) {
            atom = atomsarray[i];
            val = chres[atom.chain][atom.resi];
            if(typeof(val) == "undefined")
                continue;
            atom.ss = val;
            if(chres[atom.chain][atom.resi-1] != val)
                atom.ssbegin = true;
            if(chres[atom.chain][atom.resi+1] != val)
                atom.ssend = true;
        }
    };
    
    /**
     * @param {Array.<AtomSpec>} atoms
     * @param {string} str
     */
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
            atom.bonds = [];
            atom.bondOrder = [];
            atom.properties = {};
            atoms.push(atom);
            
        }   
        
        assignBonds(atoms);
        
        return true; 
    };
        
    // read an XYZ file from str and put the result in atoms
    /**
     * @param {Array.<AtomSpec>} atoms
     * @param {string} str
     */
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
            atom.properties = {};
            atoms[i] = atom;
        }
        assignBonds(atoms);

        return true;
    };

    // put atoms specified in sdf fromat in str into atoms
    // adds to atoms, does not replace
    /** 
     * @param {Array.<AtomSpec>} atoms
     * @param {string} str
     */
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
        var i, line;
        for (i = start; i < end; i++) {
            line = lines[offset];
            offset++;
            var atom = {};
            atom.serial = i;
            atom.x = parseFloat(line.substr(0, 10));
            atom.y = parseFloat(line.substr(10, 10));
            atom.z = parseFloat(line.substr(20, 10));
            atom.hetflag = true;
            atom.atom = atom.elem = line.substr(31, 3).replace(/ /g, "");
            atom.bonds = [];
            atom.bondOrder = [];
            atom.properties = {};
            atoms[i] = atom;
        }
        
        for (i = 0; i < bondCount; i++) {
            line = lines[offset];
            offset++;
            var from = parseInt(line.substr(0, 3)) - 1 + start;
            var to = parseInt(line.substr(3, 3)) - 1 + start;
            var order = parseInt(line.substr(6, 3));            
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
    /**
     * @param {Array.<AtomSpec>} atoms
     * @param {string} str
     * @param {boolean=} keepH
     */
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
        var nbonds = 0;
        
        if (tokens.length > 1)
            nbonds = parseInt(tokens[1]); 
        
        var offset = 4;
        var i;
        // Continue until 'Atom' section
        for (i = 3; i < lines.length; i++)
        {
            if (lines[i] == "@<TRIPOS>ATOM")
            {
                offset = i+1;
                break;
            }
        }
        
        var start = atoms.length;
        var end = start + natoms;
        var line;
        // Process ATOMS
        for (i = start; i < end; i++) {
            line = lines[offset++];
            tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
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
            for (i = 0; i < nbonds; i++)
            {
                line = lines[offset++];

                tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
                            " ");
                var from = parseInt(tokens[1]);
                fromAtom = atoms[serialToIndex[from]];
                var to = parseInt(tokens[2]);
                toAtom = atoms[serialToIndex[to]];              
                    
                // Won't be able to read aromatic bonds correctly...
                var order = parseInt(tokens[3]);
                if (isNaN(order))
                    order = 1;
                
                if (fromAtom !== undefined && toAtom !== undefined){
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
    /**
     * @param {Array.<AtomSpec>} atoms
     * @param {string} str
     * @param {keepH=} boolean
     * @param {computeStruct=} boolean
     */
    var parsePDB = function(atoms, str, keepH, computeStruct) {

        var atoms_cnt = 0;
        var noH = !keepH; // suppress hydrogens by default
        var start = atoms.length;
        var atom;
        var protein = {
            sheet : [],
            helix : []
        }; // get secondary structure straight from pdb

        var hasStruct = false;
        var serialToIndex = []; // map from pdb serial to index in atoms
        var lines = str.split("\n");
        var i, j, k, line;
        for (i = 0; i < lines.length; i++) {
            line = lines[i].replace(/^\s*/, ''); // remove indent
            var recordName = line.substr(0, 6);
            var startChain, startResi, endChain, endResi;
            if (recordName == 'ATOM  ' || recordName == 'HETATM') {
                var resn, chain, resi, icode, x, y, z, hetflag, elem, serial, altLoc, b;
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
                if (elem === '') { // for some incorrect PDB files
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
                    'bondOrder' : [],
                    'properties' : {},
                    'b' : b,
                    'pdbline' : line
                });
            } else if (recordName == 'SHEET ') {
            	hasStruct = true;
                startChain = line.substr(21, 1);
                startResi = parseInt(line.substr(22, 4));
                endChain = line.substr(32, 1);
                endResi = parseInt(line.substr(33, 4));
                protein.sheet
                        .push([ startChain, startResi, endChain, endResi ]);
            } else if (recordName == 'CONECT') {
                // MEMO: We don't have to parse SSBOND, LINK because both are
                // also
                // described in CONECT. But what about 2JYT???
                var from = parseInt(line.substr(6, 5));
                var fromAtom = atoms[serialToIndex[from]];
                for (j = 0; j < 4; j++) {
                    var to = parseInt(line.substr([ 11, 16, 21, 26 ][j], 5));
                    var toAtom = atoms[serialToIndex[to]];
                    if (fromAtom !== undefined && toAtom !== undefined) {
                        fromAtom.bonds.push(serialToIndex[to]);
                        fromAtom.bondOrder.push(1);
                    }
                }
            } else if (recordName == 'HELIX ') {
            	hasStruct = true;
                startChain = line.substr(19, 1);
                startResi = parseInt(line.substr(21, 4));
                endChain = line.substr(31, 1);
                endResi = parseInt(line.substr(33, 4));
                protein.helix
                        .push([ startChain, startResi, endChain, endResi ]);
            }

        }

        var starttime = (new Date()).getTime();
        // assign bonds - yuck, can't count on connect records
        assignPDBBonds(atoms);
        console.log("bond connecting " + ((new Date()).getTime() - starttime));
        
        
        if(computeStruct || !hasStruct) {
            starttime = (new Date()).getTime();
        	computeSecondaryStructure(atoms);
        	console.log("secondary structure " + ((new Date()).getTime() - starttime));
        }
        
        // Assign secondary structures from pdb file
        for (i = start; i < atoms.length; i++) {
            atom = atoms[i];
            if (atom === undefined)
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

    function GLModel(mid, defaultcolors) {
        // private variables
        var atoms = [];
        var id = mid;
        var molObj = null;
        var renderedMolObj = null;
        var lastStyle = null; // cache previous styles to avoid recomputation
        var lastColors = null;
        
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

  
        /** Return proper color for atom given style
         * @param {AtomSpec} atom
         * @param {AtomStyle} style
         * @return {$3Dmol.Color}
         */
        var getColorFromStyle = function(atom, style) {
            var color = atom.color;
            if (typeof (style.color) != "undefined")
                color = style.color;
            if(typeof(style.colorscheme) != "undefined" &&
            		typeof($3Dmol.elementColors[style.colorscheme]) != "undefined") {
            	var scheme = $3Dmol.elementColors[style.colorscheme];
            	if(typeof(scheme[atom.elem]) != "undefined") {
            		color = scheme[atom.elem];
            	}
            }
            var C = $3Dmol.CC.color(color);
            return C;
        }

        // cross drawing
        /**
         * 
         * @param {AtomSpec} atom
         * @param {Object.<numlike,$3Dmol.Geometry>} geos
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
            
            var c = getColorFromStyle(atom, style);
            
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
			
			v.multiplyScalar(r * 1.5);

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
		
        // bonds - both atoms must match bond style
        // standardize on only drawing for lowest to highest
        /**
		 * 
		 * @param {AtomSpec}
		 *            atom
		 * @param {Array.
		 *            <AtomSpec>} atoms
		 * @param {Object.
		 *            <numlike, $3Dmol.Geometry>} geos
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
            var geoGroup = geos[linewidth].updateGeoGroup(2*atom.bonds.length);
            
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

                var c1 = getColorFromStyle(atom, atom.style.line);
                var c2 = getColorFromStyle(atom2, atom2.style.line);
               
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
                                                                 
            var C = getColorFromStyle(atom, style);
            
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
            var C = getColorFromStyle(atom, style);
            
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
        	normalArray[start+1] = -radius;
        	normalArray[start+2] = 0;
        	
        	normalArray[start+3] = -radius;
        	normalArray[start+4] = radius;
        	normalArray[start+5] = 0;
        	
        	normalArray[start+6] = radius;
        	normalArray[start+7] = radius;
        	normalArray[start+8] = 0;
        	
        	normalArray[start+9] = radius;
        	normalArray[start+10] = -radius;
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
            var fromCap = false, toCap = false;

            var C1 = getColorFromStyle(atom, style);

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
                	var style2 = atom2.style;
                    if (!style2.stick)
                        continue; // don't sweat the details                     
                   
                    var C2 = getColorFromStyle(atom2, style2.stick);
                    
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
                            toCap = true;       
                                                
                        if (C1 != C2) {
                            mp = new $3Dmol.Vector3().addVectors(p1, p2)
                                    .multiplyScalar(0.5);
                            $3Dmol.GLDraw.drawCylinder(geo, p1, mp, bondR, C1, fromCap, false);
                            $3Dmol.GLDraw.drawCylinder(geo, mp, p2, bondR, C2, false, toCap);
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
                        var mfromCap = false; mtoCap = false; //multi bond caps
                        
                        if(bondR != atomBondR) {
                        	//assume jmol style multiple bonds - the radius doesn't fit within atom sphere
                        	mfromCap = true;
                        	mtoCap = true;
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
                                $3Dmol.GLDraw.drawCylinder(geo, p1a, mp, r, C1, mfromCap, false);
                                $3Dmol.GLDraw.drawCylinder(geo, mp, p2a, r, C2, false, mtoCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1b, mp2, r, C1, mfromCap, false);
                                $3Dmol.GLDraw.drawCylinder(geo, mp2, p2b, r, C2, false, mtoCap);
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
                                $3Dmol.GLDraw.drawCylinder(geo, p1a, mp, r, C1, mfromCap, false);
                                $3Dmol.GLDraw.drawCylinder(geo, mp, p2a, r, C2, false, mtoCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1, mp3, r, C1, fromCap, false);
                                $3Dmol.GLDraw.drawCylinder(geo, mp3, p2, r, C2, false, toCap);
                                $3Dmol.GLDraw.drawCylinder(geo, p1b, mp2, r, C1, mfromCap, false);
                                $3Dmol.GLDraw.drawCylinder(geo, mp2, p2b, r, C2, false, mtoCap);
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
        /** @type {Array.<AtomSpec>} atoms */
        var createMolObj = function(atoms) {

            console.log("creating for "+id);
            var ret = new $3Dmol.Object3D();
            var cartoonAtoms = [];
            var lineGeometries = {};
            var crossGeometries = {};
            var sphereGeometry = new $3Dmol.Geometry(true);                                                         
            var imposterGeometry = new $3Dmol.Geometry(true);                                                         
            var stickGeometry = new $3Dmol.Geometry(true);
            var i, n;
            var range = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
            for (i = 0, n = atoms.length; i < n; i++) {
                var atom = atoms[i];
                // recreate gl info for each atom as necessary
                // set up appropriate intersection spheres for clickable atoms
                if (atom && atom.style) {
                    if (atom.clickable && atom.intersectionShape === undefined)
                        atom.intersectionShape = {sphere: [], cylinder: [], line: [], triangle : []};                    
                    drawAtomSphere(atom, sphereGeometry);
                    drawAtomImposter(atom, imposterGeometry);
                    drawAtomCross(atom, crossGeometries);
                    drawBondLines(atom, atoms, lineGeometries);
                    drawBondSticks(atom, atoms, stickGeometry);
                    if (typeof (atom.style.cartoon) !== "undefined" && !atom.style.cartoon.hidden) {
                        //gradient color scheme range
                        if (atom.style.cartoon.color === 'spectrum' && typeof(atom.resi) === "number") {                            
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
                //TODO: Should have an option to choose color scheme
                if (range[0] < range[1])
                    gradientscheme = new $3Dmol.Sinebow(range[0], range[1]);

                $3Dmol.drawCartoon(ret, cartoonAtoms, gradientscheme);
                
                for (i = 0; i < ret.children.length; i++){
                    var geo = ret.children[i].geometry;
                }
            }

            // add sphere geometry
            if (sphereGeometry.vertices > 0) {
                var sphereMaterial = new $3Dmol.MeshLambertMaterial({
                    ambient : 0x000000,
                    vertexColors : true,
                    reflectivity : 0,
                });
                
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
                console
                        .log("spherei geometry " + imposterGeometry.vertices.length);

                ret.add(spherei);
            }
            
            // add stick geometry
            if (stickGeometry.vertices > 0) {
                var cylinderMaterial = new $3Dmol.MeshLambertMaterial({
                    vertexColors : true,
                    ambient : 0x000000,
                    reflectivity : 0
                });

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

                    crossGeometries[i].initTypedArrays();
                    
                    var cross = new $3Dmol.Line(crossGeometries[i], crossMaterial,
                            $3Dmol.LinePieces);

                    ret.add(cross);
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
                    atom.color = atom.color || ElementColors[atom.elem] || defaultColor;
                    atom.model = id;
                    if (atom.clickable)
                        atom.intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};
                }
            }
        };

        // add atoms to this model from molecular data string
        this.addMolData = function(data, format, options) {
            options = options || {}; 
            if (!data)
                console.error("Erorr with addMolData: No input data specified");
            
            switch (format) {
            case "xyz":
                parseXYZ(atoms, data);
                break;
            case "pdb":
                parsePDB(atoms, data, options.keepH, true, options.computeStruct);
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
            var i;
            for(i = 0; i < newatoms.length; i++) {
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
        
        
        // style the select atoms with style
        this.setStyle = function(sel, style, add) {
            
            if(!add && molObj !== null && sameObj(style, lastStyle))
                return; // no need to recompute
            
            if(add) lastStyle = null; // todo: compute merged style
            else lastStyle = style;

            // do a copy to enforce style changes through this function
            var mystyle = $.extend(true, {}, style);
            var changedAtoms = false;
            // somethings we only calculate if there is a change in a certain
            // style, although these checks will only catch cases where both
            // are either null or undefined
            for ( var i = 0; i < atoms.length; i++) {
                atoms[i].capDrawn = false; //reset for proper stick render
                
                if (this.atomIsSelected(atoms[i], sel)) {
                    changedAtoms = true;
                    if (atoms[i].clickable) 
                        atoms[i].intersectionShape = {sphere : [], cylinder : [], line : [], triangle : []};                    
    
                   
                    if(!add) atoms[i].style = {};
                    for(var s in mystyle) {
                        if(mystyle.hasOwnProperty(s)) {
                            atoms[i].style[s] = mystyle[s];
                        }
                    }
                }
            }
            
            if (changedAtoms)
                molObj = null; //force rebuild
            
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
            var i, a;
            // compute the range            
            for (i = 0; i < atoms.length; i++) {
                a = atoms[i];
                if(a.properties && typeof(a.properties[prop]) !== undefined) {
                    var p = parseFloat(a.properties[prop]);
                    if(p < min) min = p;
                    if(p > max) max = p;
                }                    
            }
            // now apply colors using scheme
            for (i = 0; i < atoms.length; i++) {
                a = atoms[i];
                if(a.properties && typeof(a.properties[prop]) !== undefined) {
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

    }

    return GLModel;
    
})();
