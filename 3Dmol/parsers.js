/**
 * $3Dmol.Parsers stores functions for parsing molecular data. They all take a string of molecular data
 * and options. The default behavior is to only read the first model in the case of multimodel files, and
 * all parsers return a list of atom list(s)
 * 
 * $3Dmol.Parsers.<ext> corresponds to the parsers for files with extension ext
 */
$3Dmol.Parsers = (function() {
    var parsers = {};

    /**
     * @param {AtomSpec[]}
     *            atomsarray
     */
    var assignBonds = function(atoms) {
        // assign bonds - yuck, can't count on connect records

        for (var i = 0, n = atoms.length; i < n; i++) {
            // Don't reindex if atoms are already indexed
            if (!atoms[i].index)
                atoms[i].index = i;
        }

        var grid = {};
        var MAX_BOND_LENGTH = 4.95; // (largest bond length, Cs) 2.25 * 2 * 1.1 (fudge factor)

        for (var index = 0; index < atoms.length; index++) {
            var atom = atoms[index];
            var x = Math.floor(atom.x / MAX_BOND_LENGTH);
            var y = Math.floor(atom.y / MAX_BOND_LENGTH);
            var z = Math.floor(atom.z / MAX_BOND_LENGTH);
            if (!grid[x]) {
                grid[x] = {};
            }
            if (!grid[x][y]) {
                grid[x][y] = {};
            }
            if (!grid[x][y][z]) {
                grid[x][y][z] = [];
            }

            grid[x][y][z].push(atom);
        }

        var findConnections = function(points, otherPoints) {
            for (var i = 0; i < points.length; i++) {
                var atom1 = points[i];
                for (var j = 0; j < otherPoints.length; j++) {
                    var atom2 = otherPoints[j];

                    if (areConnected(atom1, atom2)) {
                        //gracefully handle one-sided bonds
                        var a2i = atom1.bonds.indexOf(atom2.index);
                        var a1i = atom2.bonds.indexOf(atom1.index);
                        if (a2i == -1 && a1i == -1) {
                            atom1.bonds.push(atom2.index);
                            atom1.bondOrder.push(1);
                            atom2.bonds.push(atom1.index);
                            atom2.bondOrder.push(1);
                        } else if (a2i == -1) {
                            atom1.bonds.push(atom2.index);
                            atom1.bondOrder.push(atom2.bondOrder[a1i])
                        } else if (a1i == -1) {
                            atom2.bonds.push(atom1.index);
                            atom2.bondOrder.push(atom1.bondOrder[a2i])                            
                        }
                            
                    }
                }
            }
        }


        /*const*/ var OFFSETS = [
            {x: 0, y: 0, z: 1},
            {x: 0, y: 1, z:-1},
            {x: 0, y: 1, z: 0},
            {x: 0, y: 1, z: 1},
            {x: 1, y:-1, z:-1},
            {x: 1, y:-1, z: 0},
            {x: 1, y:-1, z: 1},
            {x: 1, y: 0, z:-1},
            {x: 1, y: 0, z: 0},
            {x: 1, y: 0, z: 1},
            {x: 1, y: 1, z:-1},
            {x: 1, y: 1, z: 0},
            {x: 1, y: 1, z: 1}
        ];
        for (var x in grid) {
            x = parseInt(x);
            for (var y in grid[x]) {
                y = parseInt(y);
                for (var z in grid[x][y]) {
                    z = parseInt(z);
                    var points = grid[x][y][z];

                    for (var i = 0; i < points.length; i++) {
                        var atom1 = points[i];
                        for (var j = i + 1; j < points.length; j++) {
                            var atom2 = points[j];
                            if (areConnected(atom1, atom2)) {
                                if (atom1.bonds.indexOf(atom2.index) == -1) {
                                    atom1.bonds.push(atom2.index);
                                    atom1.bondOrder.push(1);
                                    atom2.bonds.push(atom1.index);
                                    atom2.bondOrder.push(1);
                                }
                            }
                        }
                    }

                    for (var o = 0; o < OFFSETS.length; o++) {
                        var offset = OFFSETS[o];
                        if (!grid[x+offset.x]
                            || !grid[x+offset.x][y+offset.y]
                            || !grid[x+offset.x][y+offset.y][z+offset.z]) continue;

                        var otherPoints = grid[x + offset.x][y + offset.y][z + offset.z];
                        findConnections(points, otherPoints);
                    }
                }
            }
        }
    };

    // this is optimized for proteins where it is assumed connected
    // atoms are on the same or next residue
    /**
     * @param {AtomSpec[]}
     *            atomsarray
     */
    var assignPDBBonds = function(atomsarray) {
        // assign bonds - yuck, can't count on connect records
        var protatoms = [];
        var hetatoms = [];
        var i, n;
        for (i = 0, n = atomsarray.length; i < n; i++) {
            var atom = atomsarray[i];
            atom.index = i;
            if (atom.hetflag)
                hetatoms.push(atom);
            else
                protatoms.push(atom);
        }

        assignBonds(hetatoms);

        // sort by resid
        protatoms.sort(function(a, b) {
            if (a.chain != b.chain)
                return a.chain < b.chain ? -1 : 1;
            return a.resi - b.resi;
        });

        // for identifying connected residues
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

            for (var j = i + 1; j < protatoms.length; j++) {
                var aj = protatoms[j];
                if (aj.chain != ai.chain)
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
                    continue; // can't be connected, but later might be
                var ydiff = Math.abs(aj.y - ai.y);
                if (ydiff > maxlength)
                    continue;
                var xdiff = Math.abs(aj.x - ai.x);
                if (xdiff > maxlength)
                    continue;
                var dist = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
                if (dist > maxlengthSq)
                    continue;

                if (aj.chain == ai.chain && Math.abs(aj.resi - ai.resi) < 4)
                    continue; // ignore bonds between too close residues
                // select closest hbond
                if (dist < ai.hbondDistanceSq) {
                    ai.hbondOther = aj;
                    ai.hbondDistanceSq = dist;
                }
                if (dist < aj.hbondDistanceSq) {
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
        var i, il, c, r; // i: used in for loop, il: length of atomsarray
        var atom, val;

        //identify helices first
        for (i = 0, il = atomsarray.length; i < il; i++) {
            atom = atomsarray[i];

            if (typeof (chres[atom.chain]) === "undefined")
                chres[atom.chain] = [];
            
            if (isFinite(atom.hbondDistanceSq)) {
                var other = atom.hbondOther;
                if (typeof (chres[other.chain]) === "undefined")
                    chres[other.chain] = [];
                
                if (Math.abs(other.resi - atom.resi) === 4) {
                    // helix
                    chres[atom.chain][atom.resi] = 'h';
                } 
            }
        }
        
        // plug gaps in helices
        for (c in chres) {
            for (r = 1; r < chres[c].length - 1; r++) {
                var valbefore = chres[c][r - 1];
                var valafter = chres[c][r + 1];
                val = chres[c][r];
                if (valbefore == 'h' && valbefore == valafter && val != valbefore) {
                    chres[c][r] = valbefore;
                }
            }
        }
        
        //now potential sheets - but only if mate not part of helix
        for (i = 0, il = atomsarray.length; i < il; i++) {
            atom = atomsarray[i];

            if (isFinite(atom.hbondDistanceSq) && chres[atom.chain][atom.resi] != 'h' && atom.ss != 'h') {
                chres[atom.chain][atom.resi] = 'maybesheet';
            }
        }
        
        
        //sheets must bond to other sheets
        for (i = 0, il = atomsarray.length; i < il; i++) {
            atom = atomsarray[i];

            if (isFinite(atom.hbondDistanceSq) && chres[atom.chain][atom.resi] == 'maybesheet') {
                var other = atom.hbondOther;
                var otherval = chres[other.chain][other.resi];
                if (otherval == 'maybesheet' || otherval == 's') {
                    // true sheet
                    chres[atom.chain][atom.resi] = 's';
                    chres[other.chain][other.resi] = 's';
                }
            }
        }
        
        // plug gaps in sheets and remove singletons
        for (c in chres) {
            for (r = 1; r < chres[c].length - 1; r++) {
                var valbefore = chres[c][r - 1];
                var valafter = chres[c][r + 1];
                val = chres[c][r];
                if (valbefore == 's' && valbefore == valafter && val != valbefore) {
                    chres[c][r] = valbefore;
                }
            }
            for (r = 0; r < chres[c].length; r++) {
                val = chres[c][r];
                if (val == 'h' || val == 's') {
                    if (chres[c][r - 1] != val && chres[c][r + 1] != val)
                        delete chres[c][r];
                }
            }
        }

        
        // assign to all atoms in residue, keep track of start
        var curres = null;
        for (i = 0, il = atomsarray.length; i < il; i++) {
            atom = atomsarray[i];
            val = chres[atom.chain][atom.resi];
            if (typeof (val) == "undefined" || val == 'maybesheet')
                continue;
            atom.ss = val;
            if (chres[atom.chain][atom.resi - 1] != val)
                atom.ssbegin = true;
            if (chres[atom.chain][atom.resi + 1] != val)
                atom.ssend = true;
        }
    };
    
    
    //make sure bonds are actually two way
    var validateBonds = function(atomsarray, serialToIndex) {
        for (var i = 0, n = atomsarray.length; i < n; i++) {
            var atom = atomsarray[i];
            for(var b = 0; b < atom.bonds.length; b++) {
                var a2i = atom.bonds[b];
                var atom2 = atomsarray[a2i];
                var atomi = serialToIndex[atom.serial];
                if(atom2 && atomi) {
                    var a1i = atom2.bonds.indexOf(atomi);
                    if(a1i < 0) {
                        atom2.bonds.push(atomi);
                        atom2.bondOrder.push(atom.bondOrder[b]);
                    }
                }
            }
        }
    };
        

    /**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options
     */
    parsers.vasp = parsers.VASP = function (str, options) {
      var atoms = [[]];
      var lattice = {};

      var lines = str.replace(/^\s+/, "").split(/[\n\r]/);

      if (lines.length < 3){
        return atoms;
      }

      if (lines[1].match(/\d+/)) {
        lattice.length = parseFloat(lines[1]);
      } else {
        console.log("Warning: second line of the vasp structure file must be a number");
        return atoms;
      }

      if (lattice.length<0) {
        console.log("Warning: Vasp implementation for negative lattice lengths is not yet available");
        return atoms;
      }

      lattice.xVec = new Float32Array(lines[2].replace(/^\s+/, "").split(/\s+/));
      lattice.yVec = new Float32Array(lines[3].replace(/^\s+/, "").split(/\s+/));
      lattice.zVec = new Float32Array(lines[4].replace(/^\s+/, "").split(/\s+/));

      var matrix = new $3Dmol.Matrix4(lattice.xVec[0], lattice.xVec[1], lattice.xVec[2], 0, 
                                      lattice.yVec[0], lattice.yVec[1], lattice.yVec[2], 0,
                                      lattice.zVec[0], lattice.zVec[1], lattice.zVec[2], 0,
                                      0,                             0,               0, 1);
      
	  matrix.multiplyScalar(lattice.length)

      var modelData = atoms.modelData = [{symmetries:[], cryst:{matrix:matrix}}];
	  
      var atomSymbols=lines[5].replace(/\s+/, "").replace(/\s+$/,"").split(/\s+/);
      var atomSpeciesNumber=new Int16Array(lines[6].replace(/^\s+/, "").split(/\s+/));
      var vaspMode=lines[7].replace(/\s+/, "");


      if (vaspMode.match(/C/)) {
        vaspMode = "cartesian";
      }else if (vaspMode.match(/D/)){
        vaspMode="direct";
      } else {
        console.log("Warning: Unknown vasp mode in POSCAR file: mode must be either C(artesian) or D(irect)");
        return atoms;
      }

      if (atomSymbols.length != atomSpeciesNumber.length) {
        console.log("Warning: declaration of atomary species wrong:");
        console.log(atomSymbols);
        console.log(atomSpeciesNumber);
        return atoms;
      }

      lines.splice(0,8);

      var atomCounter = 0;

      for (var i = 0, len = atomSymbols.length; i < len; i++) {
        var atomSymbol = atomSymbols[i];
       for (var j = 0, atomLen = atomSpeciesNumber[i]; j < atomLen; j++) {

        var coords = new Float32Array(lines[atomCounter + j].replace(/^\s+/, "").split(/\s+/));

        var atom={};
        atom.elem = atomSymbol;
        if (vaspMode == "cartesian") {
          atom.x = lattice.length*coords[0];
          atom.y = lattice.length*coords[1];
          atom.z = lattice.length*coords[2];
        } else {
          atom.x = lattice.length*(coords[0]*lattice.xVec[0] + coords[1]*lattice.yVec[0] + coords[2]*lattice.zVec[0]);
          atom.y = lattice.length*(coords[0]*lattice.xVec[1] + coords[1]*lattice.yVec[1] + coords[2]*lattice.zVec[1]);
          atom.z = lattice.length*(coords[0]*lattice.xVec[2] + coords[1]*lattice.yVec[2] + coords[2]*lattice.zVec[2]);
        }

        atom.bonds=[];

        atoms[0].push(atom);
       }
        atomCounter += atomSpeciesNumber[i];
      }

      return atoms;


    }

    /**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options
     */
    parsers.cube = parsers.CUBE = function(str, options) {
        var atoms = [[]];
        var lines = str.replace(/^\s+/, "").split(/\n\r|\r+/);

        if (lines.length < 6)
            return atoms;

        var lineArr = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(
                " ");

        var natoms = Math.abs(parseFloat(lineArr[0]));

        lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");

        // might have to convert from bohr units to angstroms
        var convFactor = (parseFloat(lineArr[0]) > 0) ? 0.529177 : 1;

        // Extract atom portion; send to new GLModel...
        lines = lines.splice(6, natoms);

        var start = atoms[atoms.length-1].length;
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
                atom.elem = "Cl";

            atom.x = parseFloat(tokens[2]) * convFactor;
            atom.y = parseFloat(tokens[3]) * convFactor;
            atom.z = parseFloat(tokens[4]) * convFactor;

            atom.hetflag = true;
            atom.bonds = [];
            atom.bondOrder = [];
            atom.properties = {};
            atoms[atoms.length-1].push(atom);

        }
        for (var i = 0; i < atoms.length; i++)
            assignBonds(atoms[i]);

        return atoms;
    };

    // read an XYZ file from str and return result
    /**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options
     */
    parsers.xyz = parsers.XYZ = function(str, options) {
        
        var atoms = [[]];
        var lines = str.split(/\r?\n|\r/);
        while (lines.length > 0) {
            if (lines.length < 3)
                break;
            var atomCount = parseInt(lines[0]);
            if (isNaN(atomCount) || atomCount <= 0)
                break;
            if (lines.length < atomCount + 2)
                break;
            var offset = 2;
            var start = atoms[atoms.length-1].length;
            var end = start + atomCount;
            for (var i = start; i < end; i++) {
                var line = lines[offset++];
                var tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
                        " ");
                var atom = {};
                atom.serial = i;
                var elem = tokens[0];
                atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();
                atom.x = parseFloat(tokens[1]);
                atom.y = parseFloat(tokens[2]);
                atom.z = parseFloat(tokens[3]);
                atom.hetflag = true;
                atom.bonds = [];
                atom.bondOrder = [];
                atom.properties = {};
                atoms[atoms.length-1][i] = atom;
                if (tokens.length >= 7) {
                    atom.dx = parseFloat(tokens[4]);
                    atom.dy = parseFloat(tokens[5]);
                    atom.dz = parseFloat(tokens[6]);
                }
            }

            if (options.multimodel) {
                atoms.push([]);
                lines.splice(0, offset);
            }
            else {
                break;
            }
        }
        
        for (var i = 0; i < atoms.length; i++) {
            assignBonds(atoms[i]);
        }
        
        if (options.onemol) {
            var temp = atoms;
            atoms = [];
            atoms.push(temp[0]);
            for (var i = 1; i < temp.length; i++) {
                var offset = atoms[0].length;
                for (var j = 0; j < temp[i].length; j++) {
                    var a = temp[i][j];
                    for (var k = 0; k < a.bonds.length; k++) {
                        a.bonds[k] = a.bonds[k] + offset;
                    }
                    a.index = atoms[0].length;
                    a.serial = atoms[0].length;
                    atoms[0].push(a);
                }
            }
        }

         return atoms;
    };

    // put atoms specified in sdf fromat in str into atoms
    // adds to atoms, does not replace
    /**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options
     */
    parsers.sdf = parsers.SDF = function(str, options) {

        var atoms = [[]];
        var noH = false;
        if (typeof options.keepH !== "undefined")
            noH = !options.keepH;
        var lines = str.split(/\r?\n|\r/);
        
        while(lines.length > 0) { 
            if (lines.length < 4)
                break;
            var atomCount = parseInt(lines[3].substr(0, 3));
            if (isNaN(atomCount) || atomCount <= 0)
                break;
            var bondCount = parseInt(lines[3].substr(3, 3));
            var offset = 4;
            if (lines.length < 4 + atomCount + bondCount)
                break;

            // serial is atom's index in file; index is atoms index in 'atoms'
            var serialToIndex = [];
            var start = atoms[atoms.length-1].length;
            var end = start + atomCount;
            var i, line;
            for (i = start; i < end; i++,offset++) {
                line = lines[offset];
                var atom = {};
                var elem = line.substr(31, 3).replace(/ /g, "");
                atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();

                if (atom.elem != 'H' || !noH) {
                    atom.serial = i;
                    serialToIndex[i] = atoms[atoms.length-1].length;
                    atom.x = parseFloat(line.substr(0, 10));
                    atom.y = parseFloat(line.substr(10, 10));
                    atom.z = parseFloat(line.substr(20, 10));
                    atom.hetflag = true;
                    atom.bonds = [];
                    atom.bondOrder = [];
                    atom.properties = {};
                    atom.index = atoms[atoms.length-1].length;
                    atoms[atoms.length-1].push(atom);
                }
            }

            for (i = 0; i < bondCount; i++,offset++) {
                line = lines[offset];
                var from = serialToIndex[parseInt(line.substr(0, 3)) - 1 + start];
                var to = serialToIndex[parseInt(line.substr(3, 3)) - 1 + start];
                var order = parseInt(line.substr(6, 3));
                if (typeof (from) != 'undefined' && typeof (to) != 'undefined') {
                    atoms[atoms.length-1][from].bonds.push(to);
                    atoms[atoms.length-1][from].bondOrder.push(order);
                    atoms[atoms.length-1][to].bonds.push(from);
                    atoms[atoms.length-1][to].bondOrder.push(order);
                }
            }
            if (options.multimodel) {
                if (!options.onemol)
                    atoms.push([]);
                while (lines[offset] != "$$$$")
                    offset++
                lines.splice(0, ++offset);
            }
            else {
                break;
            }
        }

        return atoms;
    };

    // This parses the ChemDoodle json file format. Although this is registered
    // for the json file extension, other chemical json file formats exist that
    // this can not parse. Check which one you have and do not assume that
    // .json can be parsed
    parsers.cdjson = parsers.json = function(str, options) {
        var atoms = [[]];
        if (typeof str === "string") { // Str is usually automatically parsed by JQuery
            str = JSON.parse(str);
        }
        var molecules = str.m;
        var atomsInFile = molecules[0].a; // Assumes there is at least one
        var bondsInFile = molecules[0].b; // molecule and ignores any more
                                          // Ignores any shapes
        var styles = molecules[0].s;
        var parseStyle = options !== undefined && options.parseStyle !== undefined ? options.parseStyle : styles !== undefined;
        
        var offset = atoms[atoms.length-1].length; // When adding atoms their index will be
                                   // Offset by the number of existing atoms
        
        for (var i = 0; i < atomsInFile.length; i++) {
            var currentAtom = atomsInFile[i];
            var atom = {};
            atom.id = currentAtom.i; // Probably won't exist. Doesn't seem to
                                     // break anything.
            atom.x = currentAtom.x;
            atom.y = currentAtom.y;
            atom.z = currentAtom.z || 0; // Default value if file is 2D

            atom.bonds = [];
            atom.bondOrder = [];
            
            var elem = currentAtom.l || 'C';
            atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();

            atom.serial = atoms[atoms.length-1].length;
            if (parseStyle) {
                atom.style = styles[currentAtom.s || 0];
            }
            atoms[atoms.length-1].push(atom);
        }
        for (var i = 0; i < bondsInFile.length; i++) {
            var currentBond = bondsInFile[i];
            var beginIndex = currentBond.b + offset;
            var endIndex = currentBond.e + offset;
            var bondOrder = currentBond.o || 1;
            
            var firstAtom = atoms[atoms.length-1][beginIndex];
            var secondAtom = atoms[atoms.length-1][endIndex];

            firstAtom.bonds.push(endIndex);
            firstAtom.bondOrder.push(bondOrder);
            secondAtom.bonds.push(beginIndex);
            secondAtom.bondOrder.push(bondOrder);
        }
        return atoms;
    }

    // puts atoms specified in mmCIF fromat in str into atoms
    /**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options
     */
    parsers.mcif = parsers.cif = function(str, options) {
        var atoms = [];
        var noAssembly = !options.doAssembly; // don't assemble by default
        var copyMatrix = !options.duplicateAssemblyAtoms;
        var modelData = atoms.modelData = [];

        // Used to handle quotes correctly
        function splitRespectingQuotes(string, separator) {
            var sections = [];
            var sectionStart = 0;
            var sectionEnd = 0;
            while (sectionEnd < string.length) {
                while (string.substr(sectionEnd, separator.length) !== separator
                        && sectionEnd < string.length) {
                    // currently does not support escaping quotes
                    if (string[sectionEnd] === "'") {
                        sectionEnd++;
                        while (sectionEnd < string.length
                                && string[sectionEnd] !== "'") {
                            sectionEnd++;
                        }
                    } else if (string[sectionEnd] === '"') {
                        sectionEnd++;
                        while (sectionEnd < string.length
                                && string[sectionEnd] !== '"') {
                            sectionEnd++;
                        }
                    }
                    sectionEnd++;

                }
                sections.push(string.substr(sectionStart, sectionEnd
                        - sectionStart));
                sectionStart = sectionEnd = sectionEnd + separator.length;
            }
            return sections;
        }


        var lines = str.split(/\r?\n|\r/);
        // Filter text to remove comments, trailing spaces, and empty lines
        var linesFiltered = [];
        var trimDisabled = false;
        for (var lineNum = 0; lineNum < lines.length; lineNum++) {
            // first remove comments
            // incorrect if #'s are allowed in strings
            // comments might only be allowed at beginning of line, not sure
            var line = lines[lineNum].split('#')[0];

            // inside data blocks, the string must be left verbatim
            // datablocks are started with a ';' at the beginning of a line
            // and ended with a ';' on its own line.
            if (trimDisabled) {
                if (line[0] === ';') {
                    trimDisabled = false;
                }
            } else {
                if (line[0] === ';') {
                    trimDisabled = true;
                }
            }

            if (trimDisabled || line !== "") {
                if (!trimDisabled) {
                    line = line.trim();
                    if (line[0] === '_') {
                        // Replace dot separating category from data item with underscore. Dots aren't guarenteed, to makes
                        // files consistent.
                        var dot = line.split(/\s/)[0].indexOf('.');
                        if (dot > -1) {
                            line[dot] = '_';
                            line = line.substr(0,dot) + '_' + line.substr(dot + 1)
                        }
                    }
                }
                linesFiltered.push(line);
            }
        }

        var lineNum = 0;
        while (lineNum < linesFiltered.length) {
            while (! linesFiltered[lineNum].startsWith("data_") ||
                   linesFiltered[lineNum] === "data_global") {
                lineNum++;
            }
            lineNum++;

            // Process the lines and puts all of the data into an object.
            var mmCIF = {};
            while (lineNum < linesFiltered.length &&
                   ! linesFiltered[lineNum].startsWith("data_")) {
                if (linesFiltered[lineNum][0] === undefined) {
                    lineNum++;
                } else if (linesFiltered[lineNum][0] === '_') {
                    var dataItemName = (linesFiltered[lineNum].split(/\s/)[0]).toLowerCase();
                    var dataItem = (mmCIF[dataItemName] = mmCIF[dataItemName] || []);

                    // if nothing left on the line go to the next one
                    var restOfLine = linesFiltered[lineNum]
                        .substr(linesFiltered[lineNum].indexOf(dataItemName)
                                + dataItemName.length);
                    if (restOfLine === "") {
                        lineNum++;
                        if (linesFiltered[lineNum][0] === ';') {
                            var dataBlock = linesFiltered[lineNum].substr(1);
                            lineNum++;
                            while (linesFiltered[lineNum] !== ';') {
                                dataBlock = dataBlock + '\n'
                                            + linesFiltered[lineNum];
                                lineNum++;
                            }
                            dataItem.push(dataBlock);
                        } else {
                            dataItem.push(linesFiltered[lineNum]);
                        }
                    } else {
                        dataItem.push(restOfLine.trim());
                    }
                    lineNum++;
                } else if (linesFiltered[lineNum].substr(0, 5) === "loop_") {
                    lineNum++;
                    var dataItems = [];
                    while (linesFiltered[lineNum] === ""
                           || linesFiltered[lineNum][0] === '_') {
                        if (linesFiltered[lineNum] !== "") {
                            var dataItemName = (linesFiltered[lineNum].split(/\s/)[0]).toLowerCase();
                            var dataItem = (mmCIF[dataItemName] = mmCIF[dataItemName] || []);
                            dataItems.push(dataItem);
                        }
                        lineNum++;
                    }

                    var currentDataItem = 0;
                    while (lineNum < linesFiltered.length
                           && linesFiltered[lineNum][0] !== '_'
                           && !linesFiltered[lineNum].startsWith("loop_")
                           && !linesFiltered[lineNum].startsWith("data_")) {
                        var line = splitRespectingQuotes(linesFiltered[lineNum], " ");
                        for (var field = 0; field < line.length; field++) {
                            if (line[field] !== "") {
                                dataItems[currentDataItem].push(line[field]);
                                currentDataItem = (currentDataItem + 1) % dataItems.length;
                            }
                        }
                        lineNum++;
                    }
                } else {
                    lineNum++;
                }
            }

            modelData.push({symmetries:[]});

            // Pulls atom information out of the data
            atoms.push([]);
            var currentIndex = 0;
            var atomCount = mmCIF._atom_site_id !== undefined ? mmCIF._atom_site_id.length
                : mmCIF._atom_site_label.length;
            function sqr(n) {
                return n*n;
            }
            var conversionMatrix;
            if (mmCIF._cell_length_a !== undefined) {
                var a = parseFloat(mmCIF._cell_length_a);
                var b = parseFloat(mmCIF._cell_length_b);
                var c = parseFloat(mmCIF._cell_length_c);
                var alpha_deg = parseFloat(mmCIF._cell_angle_alpha) || 90;
                var beta_deg = parseFloat(mmCIF._cell_angle_beta) || 90;
                var gamma_deg = parseFloat(mmCIF._cell_angle_gamma) || 90;
                var alpha = alpha_deg * Math.PI / 180;
                var beta = beta_deg * Math.PI / 180;
                var gamma = gamma_deg * Math.PI / 180;
                var cos_alpha = Math.cos(alpha);
                var cos_beta = Math.cos(beta);
                var cos_gamma = Math.cos(gamma);
                var sin_gamma = Math.sin(gamma);
                conversionMatrix = [
                    [a, b*cos_gamma, c*cos_beta],
                    [0, b*sin_gamma, c*(cos_alpha-cos_beta*cos_gamma)/sin_gamma],
                    [0, 0, c*Math.sqrt(1-sqr(cos_alpha)-sqr(cos_beta)-sqr(cos_gamma)+2*cos_alpha*cos_beta*cos_gamma)/sin_gamma]
                ];
                modelData[modelData.length-1].cryst = {'a' : a, 'b' : b, 'c' : c, 'alpha' : alpha_deg, 'beta' : beta_deg, 'gamma' : gamma_deg};
            }
            function fractionalToCartesian(a, b, c) {
                var x = conversionMatrix[0][0]*a + conversionMatrix[0][1]*b + conversionMatrix[0][2]*c;
                var y = conversionMatrix[1][0]*a + conversionMatrix[1][1]*b + conversionMatrix[1][2]*c;
                var z = conversionMatrix[2][0]*a + conversionMatrix[2][1]*b + conversionMatrix[2][2]*c;
                return {x:x, y:y, z:z};
            }
            for (var i = 0; i < atomCount; i++) {
                if (mmCIF._atom_site_group_pdb !== undefined && mmCIF._atom_site_group_pdb[i] === "TER")
                    continue;
                var atom = {};
                if (mmCIF._atom_site_cartn_x !== undefined) {
                    atom.x = parseFloat(mmCIF._atom_site_cartn_x[i]);
                    atom.y = parseFloat(mmCIF._atom_site_cartn_y[i]);
                    atom.z = parseFloat(mmCIF._atom_site_cartn_z[i]);
                }
                else {
                    var coords = fractionalToCartesian(
                        parseFloat(mmCIF._atom_site_fract_x[i]),
                        parseFloat(mmCIF._atom_site_fract_y[i]),
                        parseFloat(mmCIF._atom_site_fract_z[i]));
                    atom.x = coords.x;
                    atom.y = coords.y;
                    atom.z = coords.z;
                }
                atom.chain = mmCIF._atom_site_auth_asym_id ? mmCIF._atom_site_auth_asym_id[i] : undefined;
                atom.resi = mmCIF._atom_site_auth_seq_id ? parseInt(mmCIF._atom_site_auth_seq_id[i]) : undefined;
                atom.resn = mmCIF._atom_site_auth_comp_id ? mmCIF._atom_site_auth_comp_id[i].trim() : undefined;
                atom.atom = mmCIF._atom_site_auth_atom_id ? mmCIF._atom_site_auth_atom_id[i].replace(/"/gm,'')  : undefined; //"primed" names are in quotes
                atom.hetflag = !mmCIF._atom_site_group_pdb || mmCIF._atom_site_group_pdb[i] === "HETA" || mmCIF._atom_site_group_pdb[i] === "HETATM";
                var elem = mmCIF._atom_site_type_symbol[i];
                atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();
                atom.bonds = [];
                atom.ss = 'c';
                atom.serial = i;
                atom.bondOrder = [];
                atom.properties = {};
                atoms[atoms.length-1].push(atom);
            }

            if (mmCIF._pdbx_struct_oper_list_id !== undefined && !noAssembly) {
                for (var i = 0; i < mmCIF._pdbx_struct_oper_list_id.length; i++) {
                    var matrix11 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[1][1]'][i]);
                    var matrix12 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[1][2]'][i]);
                    var matrix13 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[1][3]'][i]);
                    var vector1 = parseFloat(mmCIF['_pdbx_struct_oper_list_vector[1]'][i]);
                    var matrix21 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[2][1]'][i]);
                    var matrix22 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[2][2]'][i]);
                    var matrix23 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[2][3]'][i]);
                    var vector2 = parseFloat(mmCIF['_pdbx_struct_oper_list_vector[2]'][i]);
	            var matrix31 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[3][1]'][i]);
                    var matrix32 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[3][2]'][i]);
                    var matrix33 = parseFloat(mmCIF['_pdbx_struct_oper_list_matrix[3][3]'][i]);
                    var vector3 = parseFloat(mmCIF['_pdbx_struct_oper_list_vector[3]'][i]);

                    var matrix = new $3Dmol.Matrix4(matrix11, matrix12, matrix13, vector1,
                                                    matrix21, matrix22, matrix23, vector2,
                                                    matrix31, matrix32, matrix33, vector3);
                    modelData[modelData.length-1].symmetries.push(matrix);
                }
                for (var i = 0; i < atoms.length; i++) {
                    processSymmetries(modelData[modelData.length-1].symmetries, copyMatrix, atoms[i]);
                }
            }
            function parseTerm(term){
                var negative = term.match('-');
                term = term.replace(/[-xyz]/g, "");
                var fractionParts = term.split('/');

                var numerator, denominator;
                if (fractionParts[1] === undefined) {
                    denominator = 1;
                }
                else {
                    denominator = parseInt(fractionParts[1]);
                }
                if (fractionParts[0] === "") {
                    numerator = 1;
                }
                else {
                    numerator = parseInt(fractionParts[0]);
                }
                return numerator / denominator * (negative ? -1 : 1);
            }
            if (mmCIF._symmetry_equiv_pos_as_xyz !== undefined) {
                for (var sym = 0; sym < mmCIF._symmetry_equiv_pos_as_xyz.length; sym++) {
                    var transform = mmCIF._symmetry_equiv_pos_as_xyz[sym].replace(/["' ]/g,"");
                    var componentStrings = transform.split(',').map(
                        function(val){
                            return val.replace(/-/g,"+-");
                        });
                    var matrix = new $3Dmol.Matrix4(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1);
                    for (var coord = 0; coord < 3; coord++) {
                        var terms = componentStrings[coord].split('+');
                        var constant = 0, xTerm = 0, yTerm = 0, zTerm = 0;
                        for (var t = 0; t < terms.length; t++) {
                            var term = terms[t];
                            if (term === "")
                                continue;
                            var coefficient = parseTerm(term);
                            if (term.match('x')) {
                                matrix.elements[coord + 0] = coefficient;
                            }
                            else if (term.match('y')) {
                                matrix.elements[coord + 4] = coefficient;
                            }
                            else if (term.match('z')) {
                                matrix.elements[coord + 8] = coefficient;
                            }
                            else {
                                matrix.elements[coord + 12] = coefficient;
                            }
                        }
                    }
                    var conversionMatrix4 = new $3Dmol.Matrix4(
                        conversionMatrix[0][0], conversionMatrix[0][1], conversionMatrix[0][2], 0,
                        conversionMatrix[1][0], conversionMatrix[1][1], conversionMatrix[1][2], 0,
                        conversionMatrix[2][0], conversionMatrix[2][1], conversionMatrix[2][2], 0);
                    var conversionInverse = (new $3Dmol.Matrix4()).getInverse(conversionMatrix4, true);
                    matrix = (new $3Dmol.Matrix4()).multiplyMatrices(matrix, conversionInverse);
                    matrix = (new $3Dmol.Matrix4()).multiplyMatrices(conversionMatrix4, matrix);
                    modelData[modelData.length-1].symmetries.push(matrix);
                }
            }
        }
        for (var i = 0; i < atoms.length; i++) {
            assignBonds(atoms[i]);
            computeSecondaryStructure(atoms[i]);
            processSymmetries(modelData[i].symmetries, copyMatrix, atoms[i]);
        }

        return atoms;
    }

    // parse SYBYL mol2 file from string - assumed to only contain one molecule
    // tag
    /**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options
     */
    parsers.mol2 = parsers.MOL2 = function(str, options) {

        var atoms = [[]];
        var noH = false;
        if (typeof options.keepH !== "undefined")
            noH = !options.keepH;

        // assert (mol_pos < atom_pos), "Unexpected formatting of mol2 file
        // (expected 'molecule' section before 'atom' section)";

        var lines = str.substr(mol_pos, str.length).split(/\r?\n|\r/);
        
        while(lines.length > 0) { 
        
            // Note: these regex's work, though they don't match '<TRIPOS>'
            // correctly - something to do with angle brackets
            var mol_pos = str.search(/@<TRIPOS>MOLECULE/);
            var atom_pos = str.search(/@<TRIPOS>ATOM/);

            // Assuming both Molecule and Atom sections exist
            if (mol_pos == -1 || atom_pos == -1)
                break;
        
            // serial is atom's index in file; index is atoms index in 'atoms'
            var serialToIndex = []; 
            var tokens = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(
                    " ");
            var natoms = parseInt(tokens[0]);
            var nbonds = 0;

            if (tokens.length > 1)
                nbonds = parseInt(tokens[1]);

            var offset = 4;
            var i;
            // Continue until 'Atom' section
            for (i = 3; i < lines.length; i++) {
                if (lines[i] == "@<TRIPOS>ATOM") {
                    offset = i + 1;
                    break;
                }
            }
        
            var start = atoms[atoms.length-1].length;
            var end = start + natoms;
            var line;
            // Process ATOMS
            for (i = start; i < end; i++) {
                line = lines[offset++];
                tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
                var atom = {};
                // get element
                var elem = tokens[5].split('.')[0];
                atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();
                if (atom.elem == 'H' && noH) {
                    // ignore
                } else {
                    // 'index' is this atom's index in 'atoms'; 'serial' is this
                    // atom's
                    // serial id in mol2 file
                    var index = atoms[atoms.length-1].length;
                    var serial = parseInt(tokens[0]);
                    atom.serial = serial;
                    // atom.serial = i;

                    atom.x = parseFloat(tokens[2]);
                    atom.y = parseFloat(tokens[3]);
                    atom.z = parseFloat(tokens[4]);
                    atom.atom = tokens[5];
                    var charge = parseFloat(tokens[8]);
                    
                    atom.index = index;
                    atom.bonds = [];
                    atom.bondOrder = [];
                    atom.properties = {
                        'charge' : charge,
                        'partialCharge' : charge
                    };
                    serialToIndex[serial] = index;

                    atoms[atoms.length-1].push(atom);
                }
            }

            // Process BONDS
            var bonds_found = false;
            while (offset < lines.length) {
                if (lines[offset++] == "@<TRIPOS>BOND") {
                    bonds_found = true;
                    break;
                }
            }

            if (bonds_found && nbonds) {
                for (i = 0; i < nbonds; i++) {
                    line = lines[offset++];

                    tokens = line.replace(/^\s+/, "").replace(/\s+/g, " ").split(
                            " ");
                    var from = parseInt(tokens[1]);
                    fromAtom = atoms[atoms.length-1][serialToIndex[from]];
                    var to = parseInt(tokens[2]);
                    toAtom = atoms[atoms.length-1][serialToIndex[to]];

                    // Won't be able to read aromatic bonds correctly...
                    var order = parseInt(tokens[3]);
                    if (isNaN(order))
                        order = 1;

                    if (fromAtom !== undefined && toAtom !== undefined) {
                        fromAtom.bonds.push(serialToIndex[to]);
                        fromAtom.bondOrder.push(order);
                        toAtom.bonds.push(serialToIndex[from]);
                        toAtom.bondOrder.push(order);
                    }

                }
            }
            if (options.multimodel) {
                if (!options.onemol)
                    atoms.push([])
                lines.splice(0, offset);
                str = lines.join("\n"); //update for str.search
                continue;
            }
            else {
                break;
            }
        }
        return atoms;

    };

    //Covalent radii
    var bondTable = {
            H :0.37,                                                                                                                                He:0.32,
            Li:1.34,Be:0.90,                                                                                B :0.82,C :0.77,N :0.75,O :0.73,F :0.71,Ne:0.69,
            Na:1.54,Mg:1.30,                                                                                Al:1.18,Si:1.11,P :1.06,S :1.02,Cl:0.99,Ar:0.97,
            K :1.96,Ca:1.74,Sc:1.44,Ti:1.56,V :1.25,/* Cr */Mn:1.39,Fe:1.25,Co:1.26,Ni:1.21,Cu:1.38,Zn:1.31,Ga:1.26,Ge:1.22,/* As */Se:1.16,Br:1.14,Kr:1.10,
            Rb:2.11,Sr:1.92,Y :1.62,Zr:1.48,Nb:1.37,Mo:1.45,Tc:1.56,Ru:1.26,Rh:1.35,Pd:1.31,Ag:1.53,Cd:1.48,In:1.44,Sn:1.41,Sb:1.38,Te:1.35,I :1.33,Xe:1.30,
            Cs:2.25,Ba:1.98,Lu:1.60,Hf:1.50,Ta:1.38,W :1.46,Re:1.59,Os:1.44,Ir:1.37,Pt:1.28,Au:1.44,Hg:1.49,Tl:1.48,Pb:1.47,Bi:1.46,/* Po *//* At */Rn:1.45,

            // None of the bottom row or any of the Lanthanides have bond lengths
    }
    var bondLength = function(elem) {
        return bondTable[elem] || 1.6;
    }
    // return true if atom1 and atom2 are probably bonded to each other
    // based on distance alone
    var areConnected = function(atom1, atom2) {
        var maxsq = bondLength(atom1.elem) + bondLength(atom2.elem);
        maxsq += 0.25;// fudge factor, especially important for md frames, also see 1i3d
        maxsq *= maxsq;

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
        else if (distSquared < 0.5)
            return false; // maybe duplicate position.
        else if (distSquared > maxsq)
            return false;
        else if(atom1.altLoc != atom2.altLoc && atom1.altLoc != ' ' && atom2.altLoc != ' ')
            return false; // don't connect across alternate locations
        else
            return true;
    };

    //adds symmetry info to either duplicate and rotate/translate biological unit later or add extra atoms now
    var processSymmetries = function(copyMatrices, copyMatrix, atoms) {
        var end = atoms.length;
        var offset = end;
        var t, l, n; // Used in for loops
        if (!copyMatrix) { // do full assembly
            for (t = 0; t < copyMatrices.length; t++) {
                if (!copyMatrices[t].isIdentity()) {
                    var xyz = new $3Dmol.Vector3();
                    for (n = 0; n < end; n++) {
                        var bondsArr = [];
                        for (l = 0; l < atoms[n].bonds.length; l++) {
                            bondsArr.push(atoms[n].bonds[l] + offset);
                        }
                        xyz.set(atoms[n].x, atoms[n].y, atoms[n].z);
                        xyz.applyMatrix4(copyMatrices[t]);
                        var newAtom = {};
                        for (var i in atoms[n]) {
                            newAtom[i] = atoms[n][i];
                        }
                        newAtom.x = xyz.x;
                        newAtom.y = xyz.y;
                        newAtom.z = xyz.z;
                        newAtom.bonds = bondsArr;
                        atoms.push(newAtom);
                    }
                    offset = atoms.length;
                }
            }
        }
        else if(copyMatrices.length > 1) {
            for (t = 0; t < atoms.length; t++) {
                var symmetries = [];
                for (l = 0; l < copyMatrices.length; l++) {
                    if (!copyMatrices[l].isIdentity()) {
                        var newXYZ = new $3Dmol.Vector3();
                        newXYZ.set(atoms[t].x, atoms[t].y, atoms[t].z);
                        newXYZ.applyMatrix4(copyMatrices[l]);
                        symmetries.push(newXYZ);
                    }
                }
                atoms[t].symmetries = symmetries;
            }
        }
    }

    var isEmpty = function( obj ) {
        var name;
        for ( name in obj ) {
            return false;
        }
        return true;
    };

    //return one model worth of pdb, returns atoms, modelData, and remaining lines
    var getSinglePDB = function(lines, options, sslookup) {
        var atoms = [];
        var noH = !options.keepH; // suppress hydrogens by default
        var ignoreStruct = !!options.noSecondaryStructure; 
        var computeStruct = !options.noComputeSecondaryStructure;
        var noAssembly = !options.doAssembly; // don't assemble by default
        var copyMatrix = !options.duplicateAssemblyAtoms; //default true
        var selAltLoc = options.altLoc ? options.altLoc : 'A'; //default alternate location to select if present
        var modelData  = {symmetries:[]};
        var atom;
        var remainingLines = [];

        var hasStruct = false;
        var serialToIndex = []; // map from pdb serial to index in atoms
        var i, j, k, line;
        var seenbonds = {}; //sometimes connect records are duplicated as an unofficial means of relaying bond orders
        
        for (i = 0; i < lines.length; i++) {
            line = lines[i].replace(/^\s*/, ''); // remove indent
            var recordName = line.substr(0, 6);
            var startChain, startResi, endChain, endResi;
            
            if(recordName.indexOf("END") == 0) {
                remainingLines = lines.slice(i+1);
                if(recordName == "END") { //as opposed to ENDMDL
                    //reset secondary structure
                    for (var prop in sslookup) {
                        if (sslookup.hasOwnProperty(prop)) {
                            delete sslookup[prop];
                        }
                    }
                }
                break;
            }
            else if (recordName == 'ATOM  ' || recordName == 'HETATM') {
                var resn, chain, resi, icode, x, y, z, hetflag, elem, serial, altLoc, b;
                altLoc = line.substr(16, 1);
                if (altLoc != ' ' && altLoc != selAltLoc && selAltLoc != '*')
                    continue; 
                serial = parseInt(line.substr(6, 5));
                atom = line.substr(12, 4).replace(/ /g, "");
                resn = line.substr(17, 3).replace(/ /g, "");
                chain = line.substr(21, 1);
                resi = parseInt(line.substr(22, 4));
                icode = line.substr(26, 1);
                x = parseFloat(line.substr(30, 8));
                y = parseFloat(line.substr(38, 8));
                z = parseFloat(line.substr(46, 8));
                b = parseFloat(line.substr(60, 8));
                elem = line.substr(76, 2).replace(/ /g, "");
                if (elem === '' || typeof(bondTable[elem]) === 'undefined') { // for some incorrect PDB files
                    elem = line.substr(12, 2).replace(/ /g, "");
                    if(elem.length > 0 && elem[0] == 'H' && elem != 'Hg') {
                        elem = 'H'; //workaround weird hydrogen names from MD, note mercury must use lowercase
                    }
                    if(elem.length > 1) {
                        elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();   
                        if(typeof(bondTable[elem]) === 'undefined') {
                            //not a known element, probably should just use first letter
                            elem = elem[0];
                        } else if(line[0] == 'A' && elem == 'Ca') { //alpha carbon, not calcium
                            elem = "C";
                        }
                    }
                } else {
                    elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();                    
                }

                if(elem == 'H' && noH)
                    continue;
                if (recordName[0] == 'H')
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
                    'altLoc' : altLoc,
                    'chain' : chain,
                    'resi' : resi,
                    'icode' : icode,
                    'rescode' : resi + (icode != ' ' ? "^" + icode : ""), // combo
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
                if(!(startChain in sslookup)) {
                    sslookup[startChain] = {};
                }
                //mark start and end with additional character
                sslookup[startChain][startResi] = 's1';
                for(var res = startResi+1; res < endResi; res++) {
                    sslookup[startChain][res] = 's';
                }
                sslookup[startChain][endResi] = 's2';

            } else if (recordName == 'CONECT') {
                // MEMO: We don't have to parse SSBOND, LINK because both are
                // also
                // described in CONECT. But what about 2JYT???
                var from = parseInt(line.substr(6, 5));
                var fromindex = serialToIndex[from];
                var fromAtom = atoms[fromindex];
                var coffsets = [ 11, 16, 21, 26 ];
                for (j = 0; j < 4; j++) {
                    var to = parseInt(line.substr(coffsets[j], 5));
                    var toindex = serialToIndex[to];
                    var toAtom = atoms[toindex];
                    if (fromAtom !== undefined && toAtom !== undefined) {
                        // duplicated conect records indicate bond order
                        if(!seenbonds[ [fromindex,toindex] ]) {
                            seenbonds[ [fromindex,toindex] ] = 1;
                            if (fromAtom.bonds.length == 0 || fromAtom.bonds[fromAtom.bonds.length - 1] != toindex) {
                                fromAtom.bonds.push(toindex);
                                fromAtom.bondOrder.push(1);
                            }
                        } else { //update bond order
                            seenbonds[ [fromindex,toindex] ] += 1;
                            
                            for(var bi = 0; bi < fromAtom.bonds.length; bi++) {
                                if(fromAtom.bonds[bi] == toindex) {
                                    var newbo = seenbonds[ [fromindex,toindex] ];
                                    if(newbo >= 4) { //aromatic
                                        fromAtom.bondOrder[bi] = 1;
                                    } else {
                                        fromAtom.bondOrder[bi] = newbo;
                                    }
                                }
                            }
                        }
                    }
                }
            } else if (recordName == 'HELIX ') {
                hasStruct = true;
                startChain = line.substr(19, 1);
                startResi = parseInt(line.substr(21, 4));
                endChain = line.substr(31, 1);
                endResi = parseInt(line.substr(33, 4));
                if(!(startChain in sslookup)) {
                    sslookup[startChain] = {};
                }
                sslookup[startChain][startResi] = 'h1';
                for(var res = startResi+1; res < endResi; res++) {
                    sslookup[startChain][res] = 'h';
                }
                sslookup[startChain][endResi] = 'h2';

            } else if ((!noAssembly) && (recordName == 'REMARK')
                    && (line.substr(13, 5) == 'BIOMT')) {
                var n;
                var matrix = new $3Dmol.Matrix4(); 
                for (n = 1; n <= 3; n++) {
                    line = lines[i].replace(/^\s*/, '');
                    if (parseInt(line.substr(18, 1)) == n) { // check for all
                                                                // three lines
                                                                // by matching #
                                                                // @ end of
                                                                // "BIOMT" to n
                        matrix.elements[(n - 1)] = parseFloat(line.substr(23,
                                10));
                        matrix.elements[(n - 1) + 4] = parseFloat(line.substr(
                                33, 10));
                        matrix.elements[(n - 1) + 8] = parseFloat(line.substr(
                                43, 10));
                        matrix.elements[(n - 1) + 12] = parseFloat(line
                                .substr(53));
                        i++;
                    } else {
                        while (line.substr(13, 5) == 'BIOMT') {
                            i++;
                            line = lines[i].replace(/^\s*/, '');
                        }
                    }
                }
                matrix.elements[3] = 0;
                matrix.elements[7] = 0;
                matrix.elements[11] = 0;
                matrix.elements[15] = 1;
                modelData.symmetries.push(matrix);
                i--; // set i back
            } else if (recordName == 'CRYST1') {
                var a, b, c, alpha, beta, gamma;
                a = parseFloat(line.substr(7, 8));
                b = parseFloat(line.substr(16, 8));
                c = parseFloat(line.substr(25, 8));
                alpha = parseFloat(line.substr(34, 6));
                beta = parseFloat(line.substr(41, 6));
                gamma = parseFloat(line.substr(48, 6));
                modelData.cryst = {'a' : a, 'b' : b, 'c' : c, 'alpha' : alpha, 'beta' : beta, 'gamma' : gamma};
            } else if (recordName == 'ANISOU') {
                var serial = parseInt(line.substr(6, 5));
                var anisouAtomIndex = serialToIndex[serial];
                var anisouAtom = atoms[anisouAtomIndex];
                if(anisouAtom) {
                    var vals = line.substr(30).trim().split(/\s+/);
                    var uMat = {u11:parseInt(vals[0]), u22:parseInt(vals[1]), u33:parseInt(vals[2]), 
                        u12:parseInt(vals[3]), u13:parseInt(vals[4]), u23:parseInt(vals[5])};
    
                    anisouAtom["uMat"] = uMat;
                }
            }
        }

        var starttime = (new Date()).getTime();
        
        //fix any "one-way" bonds in CONECT records
        validateBonds(atoms, serialToIndex);
        // assign bonds - yuck, can't count on connect records
        assignPDBBonds(atoms);
       // console.log("bond connecting " + ((new Date()).getTime() -starttime));

        if (!noAssembly)
            processSymmetries(modelData.symmetries, copyMatrix, atoms);

        if (computeStruct  && !ignoreStruct) {
            starttime = (new Date()).getTime();
            computeSecondaryStructure(atoms);
           // console.log("secondary structure " + ((new Date()).getTime() - starttime));
        }
        starttime = (new Date()).getTime();

        // Assign secondary structures from pdb file
        if(!isEmpty(sslookup)) {
            for (i = 0; i < atoms.length; i++) {
                atom = atoms[i];
                if (atom === undefined)
                    continue;
                if(atom.chain in sslookup &&
                    atom.resi in sslookup[atom.chain]) {
                    var code = sslookup[atom.chain][atom.resi];
                    atom.ss = code[0];
                    if(code.length > 1) {
                        if(code[1] == '1') atom.ssbegin = true;
                        else if(code[1] == '2') atom.ssend = true;
                    }
                }
            }
        }
    //console.log("assign structure " + ((new Date()).getTime() - starttime));
        
        return [atoms,modelData,remainingLines];
    };


    // parse pdb file from str and create atoms
    // if computeStruct is true will always perform secondary structure
    // analysis,
    // otherwise only do analysis of SHEET/HELIX comments are missing
    /**
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options - keepH (do not strip hydrogens), noSecondaryStructure
     *            (do not compute ss), altLoc (which alternate location to select, if present; '*' to load all)
     */
    parsers.pdb = parsers.PDB = parsers.pdbqt = parsers.PDBQT = function(str, options) {

        var atoms = []; //a separate list for each model
        var sslookup = {}; //stores SHEET and HELIX info, which is shared across models
        atoms.modelData = [];
        var lines = str.split(/\r?\n|\r/);
        while(lines.length > 0) {
            var pdbinfo = getSinglePDB(lines, options, sslookup);
            var modelatoms = pdbinfo[0];
            var modelData = pdbinfo[1];
            lines = pdbinfo[2];
            
            if(modelatoms.length == 0) {
                continue; //happens when there are blank lines
            }
            if(options.multimodel && options.onemol && atoms.length > 0) {
                //merge into existing atoms
                var inc = atoms[0].length;
                for(var i = 0; i < modelatoms.length; i++) {
                    //renumber
                    var atom = modelatoms[i];
                    atom.index = i;
                    for(var b = 0; b < atom.bonds.length; b++) {
                        atom.bonds[b] += inc;
                    }
                    atoms[0].push(atom);
                }
            } else  {
                atoms.modelData.push(modelData);
                atoms.push(modelatoms);
            }
            
            if(!options.multimodel) {
                break;
            }
        }
        
        return atoms;
    };

    /**
     * Parse a pqr file from str and create atoms. A pqr file is assumed to be a
     * whitespace delimited PDB with charge and radius fields.
     *
     * @param {string}
     *            str
     * @param {ParserOptionsSpec}
     *            options - noSecondaryStructure (do not compute ss)
     */
    parsers.pqr = parsers.PQR = function(str, options) {

        var atoms = [[]];
        var atoms_cnt = 0;
        var start = atoms[atoms.length-1].length;
        var atom;
        var computeStruct = !options.noSecondaryStructure;
        var noAssembly = !options.doAssembly; // don't assemble by default
        var copyMatrix = !options.duplicateAssemblyAtoms; //default true
        var modelData = atoms.modelData = [{symmetries:[]}];
        
        var serialToIndex = []; // map from pdb serial to index in atoms
        var lines = str.split(/\r?\n|\r/);
        var i, j, k, line;
        for (i = 0; i < lines.length; i++) {
            line = lines[i].replace(/^\s*/, ''); // remove indent
            var recordName = line.substr(0, 6);
            var startChain, startResi, endChain, endResi;
            
            if (recordName.indexOf("END") == 0) {
                if (options.multimodel) {
                    if (!options.onemol)
                        atoms.push([]);
                    continue;
                }
                else {
                    break;
                }
            }
            else if (recordName == 'ATOM  ' || recordName == 'HETATM') {
                // I would have liked to split based solely on whitespace, but
                // it seems that there is no guarantee that all the fields will
                // be filled out (e.g. the chain) so this doesn't work
                var hetflag;
                var serial = parseInt(line.substr(6, 5));
                var atom = line.substr(12, 4).replace(/ /g, "");
                var resn = line.substr(17, 3);
                var chain = line.substr(21, 1);
                var resi = parseInt(line.substr(22, 4));
                // however let's split the coordinates, charge and radius by
                // whitespace
                // to support extra precision
                var vals = line.substr(30).trim().split(/\s+/);
                var x = parseFloat(vals[0]);
                var y = parseFloat(vals[1]);
                var z = parseFloat(vals[2]);
                var charge = parseFloat(vals[3]);
                var radius = parseFloat(vals[4]);

                var elem = atom[0];
                if (atom.length > 1 && atom[1].toUpperCase() != atom[1]) {
                    // slight hack - identify two character elements by the
                    // second character in the atom name being lowercase
                    elem = atom.substr(0, 2);
                }

                if (line[0] == 'H')
                    hetflag = true;
                else
                    hetflag = false;
                serialToIndex[serial] = atoms[atoms.length-1].length;
                atoms[atoms.length-1].push({
                    'resn' : resn,
                    'x' : x,
                    'y' : y,
                    'z' : z,
                    'elem' : elem,
                    'hetflag' : hetflag,
                    'chain' : chain,
                    'resi' : resi,
                    'serial' : serial,
                    'atom' : atom,
                    'bonds' : [],
                    'ss' : 'c',
                    'bondOrder' : [],
                    'properties' : {
                        'charge' : charge,
                        'partialCharge' : charge,
                        'radius' : radius
                    },
                    'pdbline' : line
                });
            } else if (recordName == 'CONECT') {
                // MEMO: We don't have to parse SSBOND, LINK because both are
                // also
                // described in CONECT. But what about 2JYT???
                var from = parseInt(line.substr(6, 5));
                var fromAtom = atoms[atoms.length-1][serialToIndex[from]];
                for (j = 0; j < 4; j++) {
                    var to = parseInt(line.substr([ 11, 16, 21, 26 ][j], 5));
                    var toAtom = atoms[atoms.length-1][serialToIndex[to]];
                    if (fromAtom !== undefined && toAtom !== undefined) {
                        fromAtom.bonds.push(serialToIndex[to]);
                        fromAtom.bondOrder.push(1);
                    }
                }
            }
        }

        // assign bonds - yuck, can't count on connect records
        for (var i = 0; i < atoms.length; i++) {
            assignPDBBonds(atoms[i]);
            if (computeStruct)
                computeSecondaryStructure(atoms[i]);
        }
        
        return atoms;
    };
    
    var fromCharCode = function( charCodeArray ){
        return String.fromCharCode.apply( null, charCodeArray ).replace(/\0/g, '');
    };
    
    var convertSS = function(val) {
      //convert mmtf code to 3dmol code
        if(val == 2) return 'h';
        if(val == 3) return 's';
        return 'c';
    };

    
    //mmtf shoul be passed as a binary UInt8Array buffer
    parsers.mmtf = parsers.MMTF = function(bindata, options) {
        
        var noH = !options.keepH; // suppress hydrogens by default
        var selAltLoc = options.altLoc ? options.altLoc : 'A'; //default alternate location to select if present
        var ignoreStruct = !!options.noSecondaryStructure; 
        var computeStruct = !options.noComputeSecondaryStructure;
        //extract symmetries - only take first assembly, apply to all models (ignoring changes for now)
        var noAssembly = !options.doAssembly; // don't assemble by default
        var copyMatrix = !options.duplicateAssemblyAtoms; //default true
        var assemblyIndex = options.assemblyIndex ? options.assemblyIndex : 0; 
        
        var mmtfData = MMTF.decode( bindata );
        var atoms = [[]];
        var modelData = atoms.modelData = [];
        
        // setup index counters
        var modelIndex = 0;
        var chainIndex = 0;
        var groupIndex = 0;
        var atomIndex = 0;

        // setup optional fields
        var chainNameList = mmtfData.chainNameList;
        var secStructList = mmtfData.secStructList;
        var insCodeList = mmtfData.insCodeList;
        var sequenceIndexList = mmtfData.sequenceIndexList;
        var bFactorList = mmtfData.bFactorList;
        var altLocList = mmtfData.altLocList;
        var occupancyList = mmtfData.occupancyList;
        var bondAtomList = mmtfData.bondAtomList;
        var bondOrderList = mmtfData.bondOrderList;
        
        var numModels = mmtfData.numModels;
        if (numModels == 0) return atoms;
        if (!options.multimodel) numModels = 1; //first only
        // hoisted loop variables
        var i, j, k, kl, m, n;
        

        
        var symmetries = [];
        if(!noAssembly && mmtfData.bioAssemblyList && mmtfData.bioAssemblyList.length > 0) {
            var transforms = mmtfData.bioAssemblyList[assemblyIndex].transformList;
            for(i = 0, n = transforms.length; i < n; i++) {
                var matrix = new $3Dmol.Matrix4(transforms[i].matrix);
                matrix.transpose();
                symmetries.push(matrix);
            }
        }

        var bondAtomListStart = 0; //for current model
        //loop over models, 
        for (m = 0; m < numModels; m++ ) {
            var modelChainCount = mmtfData.chainsPerModel[m];
            var matoms = atoms[atoms.length-1];
            var serialToIndex = []; // map to matoms index, needed for noh

            modelData.push({symmetries:symmetries});
            for( i = 0; i < modelChainCount; ++i ){

                var chainGroupCount = mmtfData.groupsPerChain[ chainIndex ];
                var chainId = fromCharCode(
                    mmtfData.chainIdList.subarray( chainIndex * 4, chainIndex * 4 + 4 )
                );
                if(mmtfData.chainNameList) {
                    chainId = fromCharCode(
                            mmtfData.chainNameList.subarray( chainIndex * 4, chainIndex * 4 + 4 )
                    );
                }

                var startGroup = groupIndex;
                var prevSS = '';
                for( j = 0; j < chainGroupCount; ++j ){ //over residues (groups)

                    var groupData = mmtfData.groupList[ mmtfData.groupTypeList[ groupIndex ] ];
                    var groupAtomCount = groupData.atomNameList.length;
                    var secStruct = 0;
                    var secStructBegin = false;
                    var secStructEnd = false;
                    
                    if( secStructList ){
                        secStruct = secStructList[ groupIndex ];
                        var sscode = convertSS(secStruct)
                        if(groupIndex  == 0 || sscode != prevSS) {
                            secStructBegin = true;
                        }
                        prevSS = sscode;
                        var nextgroup = groupIndex+1;
                        if(nextgroup >= secStructList.length || convertSS(secStructList[nextgroup] != sscode)) {
                            secStructEnd = true;
                        }
                    }
                    var insCode = null;
                    if( mmtfData.insCodeList ){
                        insCode = String.fromCharCode( insCodeList[ groupIndex ] );
                    }
                    var sequenceIndex = null;
                    if( sequenceIndexList ){
                        sequenceIndex = sequenceIndexList[ groupIndex ];
                    }

                    var groupId = mmtfData.groupIdList[ groupIndex ];
                    var groupName = groupData.groupName;
                    var startAtom = atomIndex;
                    
                    for( k = 0; k < groupAtomCount; ++k ){

                        var element = groupData.elementList[ k ];
                        if(noH && element == 'H') {
                            atomIndex += 1;
                            continue;
                        }
                        
                        var bFactor = '';
                        if( bFactorList ){
                            bFactor = bFactorList[ atomIndex ];
                        }
                        var altLoc = '';
                        if( altLocList && altLocList[ atomIndex ]){ //not zero
                            altLoc = String.fromCharCode( altLocList[ atomIndex ] );
                        }
                        var occupancy = '';
                        if( occupancyList ){
                            occupancy = occupancyList[ atomIndex ];
                        }

                        if (altLoc != '' && altLoc != selAltLoc && selAltLoc != '*') {
                            atomIndex += 1;
                            continue; 
                        }
                        
                        var atomId = mmtfData.atomIdList[ atomIndex ];
                        var atomName = groupData.atomNameList[ k ];
                        var atomCharge = 0;
                        if(groupData.atomChargeList) atomCharge = groupData.atomChargeList[ k ];
                        var xCoord = mmtfData.xCoordList[ atomIndex ];
                        var yCoord = mmtfData.yCoordList[ atomIndex ];
                        var zCoord = mmtfData.zCoordList[ atomIndex ];
                            
                        serialToIndex[atomIndex] = matoms.length;
                        matoms.push({
                            'resn' : groupName,
                            'x' : xCoord,
                            'y' : yCoord,
                            'z' : zCoord,
                            'elem' : element,
                            'hetflag' : secStruct < 0,
                            'chain' : chainId,
                            'resi' : groupId,
                            'icode' : altLoc,
                            'rescode' : groupId + (altLoc != ' ' ? "^" + altLoc : ""), // combo
                            // resi
                            // and
                            // icode
                            'serial' : atomId,
                            'altLoc' : altLoc,
                            'index' : atomIndex,
                            'atom' : atomName,
                            'bonds' : [],
                            'ss' : convertSS(secStruct),
                            'ssbegin' : secStructBegin,
                            'ssend' : secStructEnd,
                            'bondOrder' : [],
                            'properties' : {charge: atomCharge, occupancy:occupancy},
                            'b' : bFactor,
                        });

                        atomIndex += 1;
                    }
                    
                    // intra group bonds
                    var groupBondAtomList = groupData.bondAtomList;
                    for( k = 0, kl = groupData.bondOrderList.length; k < kl; ++k ){
                        var atomIndex1 = startAtom + groupBondAtomList[ k * 2 ];
                        var atomIndex2 = startAtom + groupBondAtomList[ k * 2 + 1 ];
                        var bondOrder = groupData.bondOrderList[ k ];
                        
                        //I assume bonds are only recorded once
                        var i1 = serialToIndex[atomIndex1];
                        var i2 = serialToIndex[atomIndex2];
                        var a1 = matoms[i1];
                        var a2 = matoms[i2];
                        if(a1 && a2) {
                            a1.bonds.push(i2)
                            a1.bondOrder.push(bondOrder);
                            a2.bonds.push(i1);
                            a2.bondOrder.push(bondOrder);         
                        }
                    }
                    
                    groupIndex += 1;
                }
                
                //reset for bonds
                groupIndex = startGroup;
                for( j = 0; j < chainGroupCount; ++j ){ //over residues (groups)
                    
                    groupIndex += 1;

                }

                chainIndex += 1;
            }

            
            // inter group bonds
            if( bondAtomList ){
                for( k = bondAtomListStart, kl = bondAtomList.length; k < kl; k += 2 ){
                     var atomIndex1 = bondAtomList[ k ];
                     var atomIndex2 = bondAtomList[ k + 1 ];
                     var bondOrder = bondOrderList ? bondOrderList[ k / 2 ] : 1;
                     
                     if(atomIndex1 >= atomIndex) {
                         bondAtomListStart = k;
                         break; //on next model
                     }
                     //I assume bonds are only recorded once
                     var i1 = serialToIndex[atomIndex1];
                     var i2 = serialToIndex[atomIndex2];
                     var a1 = matoms[i1];
                     var a2 = matoms[i2];
                     if(a1 && a2) {
                         a1.bonds.push(i2)
                         a1.bondOrder.push(bondOrder);
                         a2.bonds.push(i1);
                         a2.bondOrder.push(bondOrder);   
                     }
                }
            }
            
            if (options.multimodel) {
                if (!options.onemol) atoms.push([]);
            }

            if(!noAssembly) {
                for (var n = 0; n < atoms.length; n++) {        
                        processSymmetries(modelData[modelIndex].symmetries, copyMatrix, atoms[n]);
                }
            }
            modelIndex += 1;
        } 
                
        
        if (computeStruct  && !ignoreStruct) {
            computeSecondaryStructure(atoms);
        }       
        
        return atoms;
    };
    
    /**
     * Parse a prmtop file from str and create atoms
     */
    parsers.prmtop = parsers.PRMTOP = function(str, options) {
	var atoms = [];
	var count = 0;
        var lines = str.split(/\r?\n|\r/);
	if(lines.length > 0 && lines[0].includes("VERSION")){
	    var sectionList = lines.filter(function (line){	//store the relevant section lists
		return line.includes("POINTERS") || line.includes("ATOM_NAME") ||
		line.includes("CHARGE") || line.includes("RADII") || line.includes("BONDS_INC_HYDROGEN") ||
		line.includes("BONDS_WITHOUT_HYDROGEN");
	    });
	    var index = getIndex("POINTERS");
	    if (index == -1)
		return [];
	    var col = getColEleSize(index);
	    var atomCount = parseInt(lines[index+1].slice(0,col[1]));
            if (isNaN(atomCount) || atomCount <= 0)
                return [];
	    index = getIndex("ATOM_NAME");
	    if (index == -1)
		return [];
	    col = getColEleSize(index);
	    var noOfCol = col[0];
	    for (var i = 0; i < atomCount/col[0]; i++){
		if (i == parseInt(atomCount/col[0]))
		    noOfCol = atomCount % col[0]; 
		for(var j=0; j < noOfCol; j++){
		    var atom = {};
		    var properties = {"charge":"", "radii":""};
	    	    atom.serial = count;
		    atom.x = 0;
		    atom.y = 0;
		    atom.z = 0;
		    atom.atom = lines[index+1].slice(col[1]*j, col[1]*(j+1));
		    atom.elem = lines[index+1].slice(col[1]*j, col[1]*j+1);
		    atom.properties = properties;
		    atom.bonds = [];
		    atom.bondOrder = [];
		    atoms.push(atom);
		    count++;
		}
		index++;
	    }
	    index = getIndex("CHARGE");
	    if (index != -1){
	        col = getColEleSize(index);
	        count = 0;
		noOfCol = col[0];
	        for (i = 0; i < atomCount/col[0]; i++){
		    if (i == parseInt(atomCount/col[0]))
			noOfCol = atomCount % col[0];
		    for(j = 0; j < noOfCol; j++){
		       atoms[count].properties["charge"] = lines[index+1].slice(col[1]*j, col[1]*(j+1));	
		        count++;
		    }
		    index++;
	        }
	    }
	    index = getIndex("RADII");
	    if (index != -1){
		col = getColEleSize(index);
		count = 0;
		noOfCol = col[0];
		for (i = 0; i < atomCount/col[0]; i++){
		    if (i == parseInt(atomCount/col[0]))
			noOfCol = atomCount % col[0];
		    for(j = 0; j < noOfCol; j++){
			atoms[count].properties.radii = lines[index+1].slice(col[1]*j, col[1]*(j+1));
			count++;
		    }
		    index++;
		}
	    }
	    index = getIndex("BONDS_WITHOUT_HYDROGEN");
	    if (index != -1){
		col = getColEleSize(index);
		count = 0;
		noOfCol = col[0];
		var atomIndex;
		for (i = 0; i < atomCount/col[0]; i++){
		    if (i == parseInt(atomCount/col[0]))
			noOfCol = atomCount % col[0];	
		    for (j = 0; j < noOfCol; j++){
			if (count%3 == 0){
			    atomIndex = parseInt(lines[index+1].slice(col[1]*j, col[1]*(j+1))/3 + 1);
			}
			if (count%3 == 1){
			    atoms[atomIndex].bonds.push(parseInt(lines[index+1].slice(col[1]*j, col[1]*(j+1))/3 + 1));
			}
		    count++;
		    }
		index++;
		}
	    }
	    index = getIndex("BONDS_INC_HYDROGEN");
	    if (index != -1){
		col = getColEleSize(index);
		count = 0;
		noOfCol = col[0];
		var atomIndex;
		for (i = 0; i < atomCount/col[0]; i++){
		    if (i == parseInt(atomCount/col[0]))
			noOfCol = atomCount % col[0];	
		    for (j = 0; j < noOfCol; j++){
			if (count%3 == 0){
			    atomIndex = parseInt(lines[index+1].slice(col[1]*j, col[1]*(j+1))/3 + 1);
			}
			if (count%3 == 1){
			    atoms[atomIndex].bonds.push(parseInt(lines[index+1].slice(col[1]*j, col[1]*(j+1))/3 + 1));
			}
		    count++;
		    }
		index++;
		}
	    }
	}
	else{
	    return [];
	}
	function getIndex(section){
	    var index = lines.indexOf(sectionList.filter(function (line){
		return line.includes(section);
	    })[0]);	//returns the index of the line containing FLAG POINTERS
	    if (Number.isInteger(index) && index > 0){
		while(!lines[index].includes("FORMAT"))  //doing this so as to take comments into consideration
		    index++;
	    	return index;
	    }
	    else{
		return -1;
	    }
	}
	function getColEleSize(i){
	    var numberOfCol = lines[i].match(/\((\d*)\S*/); // stores the number of columns
	    var elementSize = lines[i].match(/[a-zA-Z](\d*)\)\s*/);
	    if(elementSize == null){
		elementSize = lines[i].match(/[a-zA-Z](\d*)\.\d*\)\s*/); //stores the element size
	    }
	    return [numberOfCol[1], elementSize[1]];	
	}       
        return [atoms];
    };

    /**
     * Parse a gro file from str and create atoms
     */
    parsers.gro = parsers.GRO = function(str, options) {
	var atoms = [];
        var lines = str.split(/\r?\n|\r/);
        while (lines.length > 0) {
            if (lines.length < 3)
                break;
            var atomCount = parseInt(lines[1]);
            if (isNaN(atomCount) || atomCount <= 0)
                break;
            if (lines.length < atomCount + 3)
                break;
	    atoms.push([]);
            var offset = 2;
            var start = atoms[atoms.length-1].length;
            var end = start + atomCount;
            for (var i = start; i < end; i++) {
                var line = lines[offset++];
                var atom = {};
                atom.serial = i;
                atom.atom = line.slice(10,15).trim();
		if(atom.atom.charCodeAt(1) >= 97 && atom.atom.charCodeAt(1) <= 122)
		    atom.elem = atom.atom.slice(0,2);
		else
		    atom.elem = atom.atom[0];
                atom.x = parseFloat(line.slice(20,28));
                atom.y = parseFloat(line.slice(28,36));
                atom.z = parseFloat(line.slice(36,44));
		atom.resi = line.slice(5,10);
                atom.bonds = [];
                atom.bondOrder = [];
                atom.properties = {};
		if (line.length > 44){
                    atom.dx = parseFloat(line.slice(44,52));
                    atom.dy = parseFloat(line.slice(52,60));
                    atom.dz = parseFloat(line.slice(60,68));
		}
                atoms[atoms.length-1][i] = atom;
            }
	    lines.splice(0, ++offset);
        }
	for (var i=0; i<atoms.length; i++){
	    assignBonds(atoms[i]);
	}
        return atoms;
    }

    /**
     * Parse a lammps trajectory file from str and create atoms
     */
    parsers.lammpstrj = parsers.LAMMPSTRJ = function(str, options){
 	var atoms = [];
        var dic = {'id':'serial','type':'atom','element':'elem','q':'charge','radius':'radius',
                         'x':'x','xu':'x','xs':'x','xsu':'x',
                         'y':'y','yu':'y','ys':'y','ysu':'y',
                         'z':'z','zu':'z','zs':'z','zsu':'z'};
        var lines = str.split(/\r?\n|\r/);
        var offset = 0;
        var atomCount = 0;
        var start = 0;
        while (start<lines.length-9){
            for (var j=start; j<lines.length; j++){
               if (lines[j].match(/ITEM: NUMBER OF ATOMS/))
                    atomCount = parseInt(lines[j+1]);
               if (lines[j].match(/ITEM: ATOMS/)){
                    offset = j+1;
                    break;
                }
            }
            var types = lines[offset-1].replace('ITEM: ATOMS ','').split(' ');
            atoms.push([]);
            for (var j=offset; j<offset+atomCount; j++){
                var atom = {};
                var properties = {};
                var tokens = lines[j].split(' ');
                for (var k=0; k<tokens.length; k++){
                    var prop = dic[types[k]];
                    if (prop != undefined){
                        if (prop == 'serial')
                            atom[prop] = parseInt(tokens[k]);
                        else if (prop == 'x' || prop == 'y' || prop === 'z')
                            atom[prop] = parseFloat(tokens[k]);
                        else if (prop == 'charge' || prop == 'radius')
                            properties[prop] = parseFloat(tokens[k]);
                        else
                            atom[prop] = tokens[k];
                    }
                    atom.properties = properties;
                    atom.bonds = [];
                    atom.bondOrder = [];
                }
                atoms[atoms.length-1][j-offset] = atom;   
            }
            start = offset+atomCount-1;
        }
        if (options.assignbonds){
 	    for (var i=0; i<atoms.length; i++)
	        assignBonds(atoms[i]);          
        }
        return atoms;       
    }
    return parsers;
})();
