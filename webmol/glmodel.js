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
	var maxlength = 3.5;
	var maxlengthSq = 12.25;
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
            atom.singleBonds = true;
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
            atom.singleBonds = true; //atom only makes single bonds ?
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
            if (order > 1) {
                atoms[from].singleBonds = false; atoms[to].singleBonds = false;
            }                
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
                
                if (order > 1) {
                    fromAtom.singleBonds = false; toAtom.singleBonds = false;   
                }                   
                
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
                    'singleBonds' : true,
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
        
        var defaultColor = WebMol.defaultElementColor;
        
        var ElementColors = (defaultcolors) ? defaultcolors : WebMol.defaultElementColors;


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
                
                if (this.cache[radius] !== undefined) return this.cache[radius];
                
                var dir = new WebMol.Vector3(0,1,0);    
                var w = this.basisVectors.length;
                var nvecs = [], norms = [];
                var n;
                
                for (var i = 0; i < w; i++) {
                    //bottom
                    nvecs.push(this.basisVectors[i].clone().multiplyScalar(radius));
                    //top
                    nvecs.push(this.basisVectors[i].clone().multiplyScalar(radius));
                    
                    //NOTE: this normal is used for constructing sphere caps - 
                    // cylinder normals taken care of in drawCylinder
                    n = this.basisVectors[i].clone().normalize();
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
                            toRow.push(xi+1); verticesRow.push(xi);
                            
                            continue;
                        }
                        
                        var u = x / widthSegments;
                        var v = y / heightSegments;
                        
                        //Only push first polar point
                        
                        if (!polar || x === 0) {
                            
                            if (x < widthSegments) {
                                var vertex = new WebMol.Vector3();
                                vertex.x = -radius * Math.cos(phiStart + u * phiLength) * Math.sin(thetaStart + v * thetaLength);
                                vertex.y = radius * Math.cos(thetaStart + v * thetaLength);
                                vertex.z = radius * Math.sin(phiStart + u * phiLength) * Math.sin(thetaStart + v * thetaLength);

                                if (Math.abs(vertex.x) < 1e-5) vertex.x = 0;
                                if (Math.abs(vertex.y) < 1e-5) vertex.y = 0;
                                if (Math.abs(vertex.z) < 1e-5) vertex.z = 0;

                                n = new WebMol.Vector3(vertex.x, vertex.y, vertex.z);
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
                        vertex.x = -radius * Math.cos(phiStart + u * phiLength) * Math.sin(thetaStart + v * thetaLength);
                        vertex.y = radius * Math.cos(thetaStart + v * thetaLength);
                        vertex.z = radius * Math.sin(phiStart + u * phiLength) * Math.sin(thetaStart + v * thetaLength);

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
        /**
         * 
         * @param {AtomSpec} atom
         * @param {Object.<numlike,WebMol.Geometry>} geos
         */
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
        /**
         * 
         * @param {AtomSpec} atom
         * @param {Array.<AtomSpec>} atoms
         * @param {Object.<numlike, WebMol.Geometry>} geos
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
                geos[linewidth] = new WebMol.Geometry();
            /** @type {geometryGroup} */
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

                geoGroup.__vertexArray[offset] = p1.x; geoGroup.__vertexArray[offset+1] = p1.y; geoGroup.__vertexArray[offset+2] = p1.z;
                geoGroup.__colorArray[offset] = c1.r; geoGroup.__colorArray[offset+1] = c1.g; geoGroup.__colorArray[offset+2] = c1.b;
                geoGroup.__vertexArray[offset+3] = mp.x; geoGroup.__vertexArray[offset+4] = mp.y; geoGroup.__vertexArray[offset+5] = mp.z;
                geoGroup.__colorArray[offset+3] = c1.r; geoGroup.__colorArray[offset+4] = c1.g; geoGroup.__colorArray[offset+5] = c1.b;

            }

        };

        // bonds as cylinders
        var defaultStickRadius = 0.25;

        //sphere drawing
        //See also: drawCylinder
        /** 
         * 
         * @param {AtomSpec} atom
         * @param {WebMol.Geometry} geo
         */
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
                        
                        geoGroup.__normalArray[v1offset] = n1.x; geoGroup.__normalArray[v3offset] = n3.x; geoGroup.__normalArray[v4offset] = n4.x;
                        geoGroup.__normalArray[v1offset+1] = n1.y; geoGroup.__normalArray[v3offset+1] = n3.y; geoGroup.__normalArray[v4offset+1] = n4.y;
                        geoGroup.__normalArray[v1offset+2] = n1.z; geoGroup.__normalArray[v3offset+2] = n3.z; geoGroup.__normalArray[v4offset+2] = n4.z;

                        geoGroup.__faceArray[faceoffset] = v1; 
                        geoGroup.__faceArray[faceoffset+1] = v3;
                        geoGroup.__faceArray[faceoffset+2] = v4;
                        
                        geoGroup.faceidx += 3;
                        
                    } else if (Math.abs(vertices[v3 - start].y) === radius) {
                        //face = [v1, v2, v3];            
                        //norm = [n1, n2, n3];
                        
                        geoGroup.__normalArray[v1offset] = n1.x; geoGroup.__normalArray[v2offset] = n2.x; geoGroup.__normalArray[v3offset] = n3.x;
                        geoGroup.__normalArray[v1offset+1] = n1.y; geoGroup.__normalArray[v2offset+1] = n2.y; geoGroup.__normalArray[v3offset+1] = n3.y;
                        geoGroup.__normalArray[v1offset+2] = n1.z; geoGroup.__normalArray[v2offset+2] = n2.z; geoGroup.__normalArray[v3offset+2] = n3.z;

                        geoGroup.__faceArray[faceoffset] = v1;
                        geoGroup.__faceArray[faceoffset+1] = v2;
                        geoGroup.__faceArray[faceoffset+2] = v3;
                        
                        geoGroup.faceidx += 3;
                        
                    } else {
                        //face = [v1, v2, v3, v4];
                        //norm = [n1, n2, n3, n4];
                        
                        geoGroup.__normalArray[v1offset] = n1.x; geoGroup.__normalArray[v2offset] = n2.x; geoGroup.__normalArray[v4offset] = n4.x;
                        geoGroup.__normalArray[v1offset+1] = n1.y; geoGroup.__normalArray[v2offset+1] = n2.y; geoGroup.__normalArray[v4offset+1] = n4.y;
                        geoGroup.__normalArray[v1offset+2] = n1.z; geoGroup.__normalArray[v2offset+2] = n2.z; geoGroup.__normalArray[v4offset+2] = n4.z;
                        
                        geoGroup.__normalArray[v2offset] = n2.x; geoGroup.__normalArray[v3offset] = n3.x; geoGroup.__normalArray[v4offset] = n4.x;
                        geoGroup.__normalArray[v2offset+1] = n2.y; geoGroup.__normalArray[v3offset+1] = n3.y; geoGroup.__normalArray[v4offset+1] = n4.y;
                        geoGroup.__normalArray[v2offset+2] = n2.z; geoGroup.__normalArray[v3offset+2] = n3.z; geoGroup.__normalArray[v4offset+2] = n4.z;
                        
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
               
                d.set(dir[0], dir[1], dir[2]);
                
                var dx = d.x, dy = d.y, dz = d.z;
                
                var dxy = Math.sqrt(dx*dx + dy*dy);
                var dxz, dyz;
                
                var sinA, cosA, sinB, cosB, sinC, cosC;
                
                // about z axis - Phi
                if (dxy < 0.0001) {
                   sinA = 0; cosA = 1; 
                }
                    
                else {
                    sinA = -dx / dxy; cosA = dy / dxy; 
                }                   
                  
                //recast dy in terms of new axes - z is the same
                
                dy = -sinA*dx + cosA*dy;
                dyz = Math.sqrt(dy*dy + dz*dz);    
                 
                // about new x axis - Theta
                
                if (dyz < 0.0001) {
                    sinB = 0; cosB = 1;
                }
                    
                else {
                    sinB =  dz / dyz; cosB = dy / dyz; 
                }                                                       
               
                rot = new Float32Array(9);
                rot[0] = cosA; rot[1] = sinA; rot[2] = 0;
                rot[3] = -sinA*cosB; rot[4] = cosA*cosB; rot[5] = sinB;
                rot[6] = sinA*sinB; rot[7] = -cosA*sinB; rot[8] = cosB;
                
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
            
            /** @type {Array.<number>} */
            var dir = [to.x, to.y, to.z];
            dir[0] -= from.x; dir[1] -= from.y; dir[2] -= from.z;
            
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
            var i, x, y, z;
            
            // add vertices, opposing vertices paired together
            for (i = 0; i < n; ++i) {
                
                var vi = 2*i;
                
                x = e[0]*vertices[vi].x + e[3]*vertices[vi].y + e[6]*vertices[vi].z;
                y = e[1]*vertices[vi].x + e[4]*vertices[vi].y + e[7]*vertices[vi].z;
                z =                       e[5]*vertices[vi].y + e[8]*vertices[vi].z;
                              
                //var xn = x/radius, yn = y/radius, zn = z/radius;
                
                offset = 3*(start + vi); faceoffset = geoGroup.faceidx;
                
                //from
                geoGroup.__vertexArray[offset] = x + from.x;
                geoGroup.__vertexArray[offset+1] = y + from.y;
                geoGroup.__vertexArray[offset+2] = z + from.z;             
                //to
                geoGroup.__vertexArray[offset+3] = x + to.x;
                geoGroup.__vertexArray[offset+4] = y + to.y;
                geoGroup.__vertexArray[offset+5] = z + to.z;
                
                //normals
                geoGroup.__normalArray[offset] = x; geoGroup.__normalArray[offset+3] = x;
                geoGroup.__normalArray[offset+1] = y; geoGroup.__normalArray[offset+4] = y;
                geoGroup.__normalArray[offset+2] = z; geoGroup.__normalArray[offset+5] = z;
                
                //colors               
                geoGroup.__colorArray[offset] = color.r; geoGroup.__colorArray[offset+3] = color.r;
                geoGroup.__colorArray[offset+1] = color.g; geoGroup.__colorArray[offset+4] = color.g;
                geoGroup.__colorArray[offset+2] = color.b; geoGroup.__colorArray[offset+5] = color.b;  
                
                //faces
                // 0 - 2 - 1
                geoGroup.__faceArray[faceoffset] = fromRow[i] + start;
                geoGroup.__faceArray[faceoffset+1] = fromRow[i+1] + start;
                geoGroup.__faceArray[faceoffset+2] = toRow[i] + start;
                // 1 - 2 - 3
                geoGroup.__faceArray[faceoffset+3] = toRow[i] + start;
                geoGroup.__faceArray[faceoffset+4] = fromRow[i+1] + start;
                geoGroup.__faceArray[faceoffset+5] = toRow[i+1] + start;
                
                geoGroup.faceidx += 6;
                
            }
            
         
            //SPHERE CAPS         

            if (drawcaps) {

                // h - sphere rows, verticesRows.length - 2
                var ystart = (toCap) ? 0 : h/2;
                var yend = (fromCap) ? h + 1 : h/2+1;
                
                var v1, v2, v3, v4, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4,
                nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4, nz1, nz2, nz3, nz4,
                v1offset, v2offset, v3offset, v4offset;              
                
                for (y = ystart; y < yend; y++) {
                    if (y === h/2)
                        continue;
                    // n number of points for each level (verticesRows[i].length - 1)
                    var cap = (y <= h/2) ? to : from;
                    
                    for (x = 0; x < n; x++) {
                                           
                        faceoffset = geoGroup.faceidx;
                        
                        v1 = verticesRows[y][x + 1]; v1offset = (v1 + start) * 3;
                        v2 = verticesRows[y][x]; v2offset = (v2 + start) * 3;
                        v3 = verticesRows[y + 1][x]; v3offset = (v3 + start) * 3;
                        v4 = verticesRows[y + 1][x + 1]; v4offset = (v4 + start) * 3;
                        
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
                        
                        geoGroup.__vertexArray[v1offset] = x1 + cap.x; 
                        geoGroup.__vertexArray[v2offset] = x2 + cap.x;
                        geoGroup.__vertexArray[v3offset] = x3 + cap.x; 
                        geoGroup.__vertexArray[v4offset] = x4 + cap.x;
    
                        geoGroup.__vertexArray[v1offset+1] = y1 + cap.y; 
                        geoGroup.__vertexArray[v2offset+1] = y2 + cap.y;
                        geoGroup.__vertexArray[v3offset+1] = y3 + cap.y; 
                        geoGroup.__vertexArray[v4offset+1] = y4 + cap.y;
    
                        geoGroup.__vertexArray[v1offset+2] = z1 + cap.z; 
                        geoGroup.__vertexArray[v2offset+2] = z2 + cap.z;
                        geoGroup.__vertexArray[v3offset+2] = z3 + cap.z; 
                        geoGroup.__vertexArray[v4offset+2] = z4 + cap.z;
    
                        geoGroup.__colorArray[v1offset] = color.r; geoGroup.__colorArray[v2offset] = color.r;
                        geoGroup.__colorArray[v3offset] = color.r; geoGroup.__colorArray[v4offset] = color.r;
    
                        geoGroup.__colorArray[v1offset+1] = color.g; geoGroup.__colorArray[v2offset+1] = color.g;
                        geoGroup.__colorArray[v3offset+1] = color.g; geoGroup.__colorArray[v4offset+1] = color.g;
    
                        geoGroup.__colorArray[v1offset+2] = color.b; geoGroup.__colorArray[v2offset+2] = color.b;
                        geoGroup.__colorArray[v3offset+2] = color.b; geoGroup.__colorArray[v4offset+2] = color.b;                                      
                             
                        
                                             
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
                            
                            geoGroup.__normalArray[v1offset] = nx1; geoGroup.__normalArray[v3offset] = nx3; geoGroup.__normalArray[v4offset] = nx4;
                            geoGroup.__normalArray[v1offset+1] = ny1; geoGroup.__normalArray[v3offset+1] = ny3; geoGroup.__normalArray[v4offset+1] = ny4;
                            geoGroup.__normalArray[v1offset+2] = nz1; geoGroup.__normalArray[v3offset+2] = nz3; geoGroup.__normalArray[v4offset+2] = nz4;
    
                            geoGroup.__faceArray[faceoffset] = v1 + start; 
                            geoGroup.__faceArray[faceoffset+1] = v3 + start;
                            geoGroup.__faceArray[faceoffset+2] = v4 + start;
                            
                            geoGroup.faceidx += 3;
                            
                        } 
                        
                        //else if (Math.abs(vobj.sphereVertices[v3].y) === radius) {
                        else if (y === yend - 1) {
                            //face = [v1, v2, v3];            
                            //norm = [n1, n2, n3];
                            
                            geoGroup.__normalArray[v1offset] = nx1; geoGroup.__normalArray[v2offset] = nx2; geoGroup.__normalArray[v3offset] = nx3;
                            geoGroup.__normalArray[v1offset+1] = ny1; geoGroup.__normalArray[v2offset+1] = ny2; geoGroup.__normalArray[v3offset+1] = ny3;
                            geoGroup.__normalArray[v1offset+2] = nz1; geoGroup.__normalArray[v2offset+2] = nz2; geoGroup.__normalArray[v3offset+2] = nz3;
                            
                            geoGroup.__faceArray[faceoffset] = v1 + start;
                            geoGroup.__faceArray[faceoffset+1] = v2 + start;
                            geoGroup.__faceArray[faceoffset+2] = v3 + start;
                            
                            geoGroup.faceidx += 3;
                            
                        } 
                        
                        else {
                            //face = [v1, v2, v3, v4];
                            //norm = [n1, n2, n3, n4];
                            
                            geoGroup.__normalArray[v1offset] = nx1; geoGroup.__normalArray[v2offset] = nx2; geoGroup.__normalArray[v4offset] = nx4;
                            geoGroup.__normalArray[v1offset+1] = ny1; geoGroup.__normalArray[v2offset+1] = ny2; geoGroup.__normalArray[v4offset+1] = ny4;
                            geoGroup.__normalArray[v1offset+2] = nz1; geoGroup.__normalArray[v2offset+2] = nz2; geoGroup.__normalArray[v4offset+2] = nz4;
                            
                            geoGroup.__normalArray[v2offset] = nx2; geoGroup.__normalArray[v3offset] = nx3; geoGroup.__normalArray[v4offset] = nx4;
                            geoGroup.__normalArray[v2offset+1] = ny2; geoGroup.__normalArray[v3offset+1] = ny3; geoGroup.__normalArray[v4offset+1] = ny4;
                            geoGroup.__normalArray[v2offset+2] = nz2; geoGroup.__normalArray[v3offset+2] = nz3; geoGroup.__normalArray[v4offset+2] = nz4;
                            
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
                                var cylinder1 = new WebMol.Cylinder(p1 , mp , bondR);
                                var sphere1 = new WebMol.Sphere(p1 , bondR);
                                atom.intersectionShape.cylinder.push(cylinder1);   
                                atom.intersectionShape.sphere.push(sphere1);                             
                            }
                            if (atom2.clickable){
                                var cylinder2 = new WebMol.Cylinder(p2 , mp , bondR);
                                var sphere2 = new WebMol.Sphere(p2 , bondR);
                                atom2.intersectionShape.cylinder.push(cylinder2);
                                atom2.intersectionShape.sphere.push(sphere2);
                            }

                        }
                        
                    } 
                    
                    else if (atom.bondOrder[i] > 1) {
                        fromCap = false; toCap = false;
                        var dir = p2.clone();
                        var v = null;
                        dir.sub(p1);
                        
                        var r, p1a, p1b, p2a, p2b;
                        var cylinder1a, cylinder1b, cylinder1c;
                        var i2, j2, atom3, p3, dir2;
                        if (atom.bonds.length === 1) {
                            if (atom2.bonds.length === 1) {
                                v = dir.clone();
                                if (Math.abs(v.x) > 0.0001)
                                    v.y += 1;
                                else
                                    v.x += 1;
                            } 
                            else {
                                i2 = (i + 1) % atom2.bonds.length;
                                j2 = atom2.bonds[i2];
                                atom3 = atoms[j2];
                                p3 = new WebMol.Vector3(atom3.x, atom3.y, atom3.z);

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
                            p3 = new WebMol.Vector3(atom3.x, atom3.y, atom3.z);

                            dir2 = p3.clone();
                            dir2.sub(p1);

                            v = dir2.clone();
                            v.cross(dir);
                        }
                        
                        if (atom.bondOrder[i] == 2) {
                            r = bondR / 2.5;
                            v.cross(dir);
                            v.normalize();
                            v.multiplyScalar(r * 1.5);

                            p1a = p1.clone();
                            p1a.add(v);
                            p1b = p1.clone();
                            p1b.sub(v);

                            p2a = p1a.clone();
                            p2a.add(dir);
                            p2b = p1b.clone();
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
                                    cylinder1a = new WebMol.Cylinder(p1a , mp , r);
                                    cylinder1b = new WebMol.Cylinder(p1b , mp2 , r);
                                    atom.intersectionShape.cylinder.push(cylinder1a);
                                    atom.intersectionShape.cylinder.push(cylinder1b);
                                }
                                if (atom2.clickable) {
                                    cylinder2a = new WebMol.Cylinder(p2a , mp , r);
                                    cylinder2b = new WebMol.Cylinder(p2b , mp2 , r);
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
                                    cylinder1a = new WebMol.Cylinder(p1a.clone(), mp.clone(), r);
                                    cylinder1b = new WebMol.Cylinder(p1b.clone(), mp2.clone(), r);
                                    cylinder1c = new WebMol.Cylinder(p1.clone(), mp3.clone(), r);
                                    atom.intersectionShape.cylinder.push(cylinder1a);
                                    atom.intersectionShape.cylinder.push(cylinder1b);
                                    atom.intersectionShape.cylinder.push(cylinder1c);
                                } 
                                if (atom2.clickable) {                               
                                    cylinder2a = new WebMol.Cylinder(p2a.clone(), mp.clone(), r);
                                    cylinder2b = new WebMol.Cylinder(p2b.clone(), mp2.clone(), r);
                                    cylinder2c = new WebMol.Cylinder(p2.clone(), mp3.clone(), r);
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
        /** @type {Array.<AtomSpec>} atoms */
        var createMolObj = function(atoms) {

            console.log("creating for "+id);
            var ret = new WebMol.Object3D();
            var cartoonAtoms = [];
            var lineGeometries = {};
            var crossGeometries = {};
            var sphereGeometry = new WebMol.Geometry(true);                                                         
            var stickGeometry = new WebMol.Geometry(true);
            var i, n;
            
            for (i = 0, n = atoms.length; i < n; i++) {
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
                    if (typeof (atom.style.cartoon) !== "undefined" && !atom.style.cartoon.hidden) {
                        cartoonAtoms.push(atom);
                    }
                    

                }
            }
            // create cartoon if needed - this is a whole model analysis
            if (cartoonAtoms.length > 0) {
                WebMol.drawCartoon(ret, cartoonAtoms, false);
                
                for (i = 0; i < ret.children.length; i++){
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
                        .log("sphere geometry " + sphereGeometry.vertices.length);

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
            
            //var linewidth;
            // add any line geometries, distinguished by line width
            for (i in lineGeometries) {
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
            for (i in crossGeometries) {
                if (crossGeometries.hasOwnProperty(i)) {
                    var linewidth = i;
                    var crossMaterial = new WebMol.LineBasicMaterial({
                        linewidth : linewidth,
                        vertexColors : true
                    });

                    crossGeometries[i].initTypedArrays();
                    
                    var cross = new WebMol.Line(crossGeometries[i], crossMaterial,
                            WebMol.LinePieces);

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
                parsePDB(atoms, data, options.keepH, options.computeStruct);
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
