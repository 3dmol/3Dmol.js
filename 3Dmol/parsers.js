/**
 * $3Dmol.Parsers stores functions for parsing molecular data. The all take an
 * atom list (which gets filled out) and a string.
 * 
 * $3Dmol.Parsers.<ext> corresponds to the parsers for files with extension ext
 */
$3Dmol.Parsers = (function() {
    var parsers = {};

    /**
     * @param {AtomSpec[]}
     *            atomsarray
     */
    var assignBonds = function(atomsarray) {
        // assign bonds - yuck, can't count on connect records
        var atoms = atomsarray.slice(0);
        var i, j, n;
        for (i = 0, n = atomsarray.length; i < n; i++) {
            // Don't reindex if atoms are already indexed
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
                if (aj.z - ai.z > 4.725) // can't be connected
                    break;
                else if (Math.abs(aj.x - ai.x) > 4.725 || Math.abs(aj.y - ai.y) > 4.725) { // can't be connected either
                    continue;
                }
                else if (areConnected(ai, aj)) {
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
        var i, il, c, r;
        var atom, val;

        for (i = 0, il = atomsarray.length; i < il; i++) {
            atom = atomsarray[i];

            if (typeof (chres[atom.chain]) === "undefined")
                chres[atom.chain] = [];

            if (isFinite(atom.hbondDistanceSq)) {
                var other = atom.hbondOther;
                if (Math.abs(other.resi - atom.resi) === 4) {
                    // helix
                    chres[atom.chain][atom.resi] = 'h';
                } else { // otherwise assume sheet
                    chres[atom.chain][atom.resi] = 's';
                }
            }
        }

        // plug gaps and remove singletons
        for (c in chres) {
            for (r = 1; r < chres[c].length - 1; r++) {
                var valbefore = chres[c][r - 1];
                var valafter = chres[c][r + 1];
                val = chres[c][r];
                if (valbefore == valafter && val != valbefore) {
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
            if (typeof (val) == "undefined")
                continue;
            atom.ss = val;
            if (chres[atom.chain][atom.resi - 1] != val)
                atom.ssbegin = true;
            if (chres[atom.chain][atom.resi + 1] != val)
                atom.ssend = true;
        }
    };

    /**
     * @param {AtomSpec[]}
     *            atoms
     * @param {string}
     *            str
     */
    parsers.cube = parsers.CUBE = function(atoms, str, options) {
        var lines = str.replace(/^\s+/, "").split(/[\n\r]+/);

        if (lines.length < 6)
            return;

        var lineArr = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(
                " ");

        var natoms = Math.abs(parseFloat(lineArr[0]));

        lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");

        // might have to convert from bohr units to angstroms
        var convFactor = (parseFloat(lineArr[0]) > 0) ? 0.529177 : 1;

        // Extract atom portion; send to new GLModel...
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
                atom.elem = "Cl";

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
     * @param {AtomSpec[]}
     *            atoms
     * @param {string}
     *            str
     */
    parsers.xyz = parsers.XYZ = function(atoms, str, options) {

        var lines = str.split(/\r?\n/);
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
            atoms[i] = atom;
        }
        assignBonds(atoms);

        return true;
    };

    // put atoms specified in sdf fromat in str into atoms
    // adds to atoms, does not replace
    /**
     * @param {AtomSpec[]}
     *            atoms
     * @param {string}
     *            str
     */
    parsers.sdf = parsers.SDF = function(atoms, str, options) {

        var noH = false;
        if (typeof options.keepH !== "undefined")
            noH = !options.keepH;
        var lines = str.split(/\r?\n/);
        if (lines.length < 4)
            return;
        var atomCount = parseInt(lines[3].substr(0, 3));
        if (isNaN(atomCount) || atomCount <= 0)
            return;
        var bondCount = parseInt(lines[3].substr(3, 3));
        var offset = 4;
        if (lines.length < 4 + atomCount + bondCount)
            return;

        // serial is atom's index in file; index is atoms index in 'atoms'
        var serialToIndex = [];
        var start = atoms.length;
        var end = start + atomCount;
        var i, line;
        for (i = start; i < end; i++) {
            line = lines[offset];
            offset++;
            var atom = {};
            var elem = line.substr(31, 3).replace(/ /g, "");
            atom.atom = atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();

            if (atom.elem != 'H' || !noH) {
                atom.serial = i;
                serialToIndex[i] = atoms.length;
                atom.x = parseFloat(line.substr(0, 10));
                atom.y = parseFloat(line.substr(10, 10));
                atom.z = parseFloat(line.substr(20, 10));
                atom.hetflag = true;
                atom.bonds = [];
                atom.bondOrder = [];
                atom.properties = {};
                atoms.push(atom);
            }
        }

        for (i = 0; i < bondCount; i++) {
            line = lines[offset];
            offset++;
            var from = serialToIndex[parseInt(line.substr(0, 3)) - 1 + start];
            var to = serialToIndex[parseInt(line.substr(3, 3)) - 1 + start];
            var order = parseInt(line.substr(6, 3));
            if (typeof (from) != 'undefined' && typeof (to) != 'undefined') {
                atoms[from].bonds.push(to);
                atoms[from].bondOrder.push(order);
                atoms[to].bonds.push(from);
                atoms[to].bondOrder.push(order);
            }
        }

        return true;
    };

    // This parses the ChemDoodle json file format. Although this is registered
    // for the json file extension, other chemical json file formats exist that
    // this can not parse. Check which one you have and do not assume that
    // .json can be parsed
    parsers.cdj = parsers.jso = // Hack because the file format is truncated
                                // at the moment
    parsers.cdjson = parsers.json = function(atoms, str, options, modelData) {
        if (typeof str === "string") { // Str is usually automatically parsed by JQuery
            str = JSON.parse(str);
        }
        var molecules = str.m;
        var atomsInFile = molecules[0].a; // Assumes there is at least one
        var bondsInFile = molecules[0].b; // molecule and ignores any more
                                          // Ignores any shapes
        var offset = atoms.length; // When adding atoms their index will be
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
            atoms.push(atom);
        }
        for (var i = 0; i < bondsInFile.length; i++) {
            var currentBond = bondsInFile[i];
            var beginIndex = currentBond.b + offset;
            var endIndex = currentBond.e + offset;
            var bondOrder = currentBond.o || 1;
            
            var firstAtom = atoms[beginIndex];
            var secondAtom = atoms[endIndex];

            firstAtom.bonds.push(endIndex);
            firstAtom.bondOrder.push(bondOrder);
            secondAtom.bonds.push(beginIndex);
            secondAtom.bondOrder.push(bondOrder);
        }
    }

    // puts atoms specified in mmCIF fromat in str into atoms
    /**
     * @param {AtomSpec[]}
     *            atoms
     * @param {string}
     *            str
     */
    parsers.mcif = parsers.cif = function(atoms, str, options, modelData) {
    
        var noAssembly = !options.doAssembly; // don't assemble by default
        var copyMatrix = !options.duplicateAssemblyAtoms;
        modelData.symmetries = [];

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


        var lines = str.split(/\r?\n/);
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

        // Process the lines and puts all of the data into an object.
        var mmCIF = {};
        var lineNum = 0;
        while (lineNum < linesFiltered.length) {
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
                        && linesFiltered[lineNum].substr(0, 5) !== "loop_") {
                    var line = splitRespectingQuotes(linesFiltered[lineNum],
                            " ");
                    for (var field = 0; field < line.length; field++) {
                        if (line[field] !== "") {
                            dataItems[currentDataItem].push(line[field]);
                            currentDataItem = (currentDataItem + 1)
                                    % dataItems.length;
                        }
                    }
                    lineNum++;
                }
            } else {
                lineNum++;
            }
        }

        // Pulls atom information out of the data
        var atomsPreBonds = [];
        var currentIndex = 0;
        var atomCount = mmCIF._atom_site_id !== undefined ? mmCIF._atom_site_id.length
                        : mmCIF._atom_site_label.length;
        function sqr(n) {
            return n*n;
        }
        var cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma, conversionMatrix;
        if (mmCIF._cell_length_a !== undefined) {
            var a = cell_a = parseFloat(mmCIF._cell_length_a);
            var b = cell_b = parseFloat(mmCIF._cell_length_b);
            var c = cell_c = parseFloat(mmCIF._cell_length_c);
            var alpha = cell_alpha = parseFloat(mmCIF._cell_angle_alpha) * Math.PI / 180 || Math.PI / 2;
            var beta = cell_beta = parseFloat(mmCIF._cell_angle_beta) * Math.PI / 180 || Math.PI / 2;
            var gamma = cell_gamma = parseFloat(mmCIF._cell_angle_gamma) * Math.PI / 180 || Math.PI / 2;
            var cos_alpha = Math.cos(alpha);
            var cos_beta = Math.cos(beta);
            var cos_gamma = Math.cos(gamma);
            var sin_gamma = Math.sin(gamma);
            conversionMatrix = [
                [a, b*cos_gamma, c*cos_beta],
                [0, b*sin_gamma, c*(cos_alpha-cos_beta*cos_gamma)/sin_gamma],
                [0, 0, c*Math.sqrt(1-sqr(cos_alpha)-sqr(cos_beta)-sqr(cos_gamma)+2*cos_alpha*cos_beta*cos_gamma)/sin_gamma]
            ];
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
            atom.hetflag = mmCIF._atom_site_group_pdb ? mmCIF._atom_site_group_pdb[i] === "HETA" : true;
            var elem = mmCIF._atom_site_type_symbol[i];
            atom.elem = elem[0].toUpperCase() + elem.substr(1).toLowerCase();
            atom.bonds = [];
            atom.ss = 'c';
            atom.serial = i;
            atom.bondOrder = [];
            atom.properties = {};
            atom.index = currentIndex++;
            atomsPreBonds[atom.index] = atom;
        }

        // create a hash table of the atoms using label and sequence as keys
        var atomHashTable = {};
        for (var i = 0; i < atomCount; i++) {
            var label_alt = (mmCIF._atom_site_label_alt_id || [])[i];
            if (label_alt === undefined) {
                label_alt = '.';
            }
            var label_asym = (mmCIF._atom_site_label_asym_id || [])[i];
            if (label_asym === undefined) {
                label_asym = '.';
            }
            var label_atom = (mmCIF._atom_site_label_atom_id || [])[i];
            if (label_atom === undefined) {
                label_atom = '.';
            }
            var label_seq = (mmCIF._atom_site_label_seq_id || [])[i];
            if (label_seq === undefined) {
                label_seq = '.';
            }

            if (atomHashTable[label_alt] === undefined) {
                atomHashTable[label_alt] = {};
            }
            if (atomHashTable[label_alt][label_asym] === undefined) {
                atomHashTable[label_alt][label_asym] = {};
            }
            if (atomHashTable[label_alt][label_asym][label_atom] === undefined) {
                atomHashTable[label_alt][label_asym][label_atom] = {};
            }

            atomHashTable[label_alt][label_asym][label_atom][label_seq] = i;
        }

        if (false && mmCIF._struct_conn && mmCIF._struct_conn_id) {
            for (var i = 0; i < mmCIF._struct_conn_id.length; i++) {
                var offset = atoms.length;

                var alt = (mmCIF._struct_conn_ptnr1_label_alt_id || [])[i];
                if (alt === undefined) {
                    alt = ".";
                }
                var asym = (mmCIF._struct_conn_ptnr1_label_asym_id || [])[i];
                if (asym === undefined) {
                    asym = ".";
                }
                var atom = (mmCIF._struct_conn_ptnr1_label_atom_id || [])[i];
                if (atom === undefined) {
                    atom = ".";
                }
                var seq = (mmCIF._struct_conn_ptnr1_label_seq_id || [])[i];
                if (seq === undefined) {
                    seq = ".";
                }

                var id1 = atomHashTable[alt][asym][atom][seq];
                // if (atomsPreBonds[id1] === undefined) continue;
                var index1 = atomsPreBonds[id1].index;

                var alt = (mmCIF._struct_conn_ptnr2_label_alt_id || [])[i];
                if (alt === undefined) {
                    alt = ".";
                }
                var asym = (mmCIF._struct_conn_ptnr2_label_asym_id || [])[i];
                if (asym === undefined) {
                    asym = ".";
                }
                var atom = (mmCIF._struct_conn_ptnr2_label_atom_id || [])[i];
                if (atom === undefined) {
                    atom = ".";
                }
                var seq = (mmCIF._struct_conn_ptnr2_label_seq_id || [])[i];
                if (seq === undefined) {
                    seq = ".";
                }

                var id2 = atomHashTable[alt][asym][atom][seq];
                if (atomsPreBonds[id2] === undefined)
                    continue;
                var index2 = atomsPreBonds[id2].index;

                atomsPreBonds[id1].bonds.push(index2 + offset);
                atomsPreBonds[id1].bondOrder.push(1);
                atomsPreBonds[id2].bonds.push(index1 + offset);
                atomsPreBonds[id2].bondOrder.push(1);
                console.log("connected " + index1 + " and " + index2);
            }
        }

        // atoms = atoms.concat(atomsPreBonds);
        for (var i = 0; i < atomsPreBonds.length; i++) {
            delete atomsPreBonds[i].index;
            atoms.push(atomsPreBonds[i]);
        }

        assignBonds(atoms);
        computeSecondaryStructure(atoms);
        
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

                var matrix = new $3Dmol.Matrix4(matrix11, matrix12, matrix13,
                        vector1, matrix21, matrix22, matrix23, vector2,
                        matrix31, matrix32, matrix33, vector3);
                modelData.symmetries.push(matrix);
            }
            processSymmetries("mcif", modelData.symmetries, copyMatrix, atoms);
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
                var transform = mmCIF._symmetry_equiv_pos_as_xyz[sym];
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
                modelData.symmetries.push(matrix);
            }
            processSymmetries("mcif", modelData.symmetries, copyMatrix, atoms);
        }
    }

    // parse SYBYL mol2 file from string - assumed to only contain one molecule
    // tag
    // TODO: Figure out how to handle multi molecule files (for SDF, too)
    /**
     * @param {AtomSpec[]}
     *            atoms
     * @param {string}
     *            str
     * @param {Object}
     *            options - keepH (do not strip hydrogens)
     */
    parsers.mol2 = parsers.MOL2 = function(atoms, str, options) {

        var noH = false;
        if (typeof options.keepH !== "undefined")
            noH = !options.keepH;

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

        var lines = str.substr(mol_pos, str.length).split(/\r?\n/);
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

        var start = atoms.length;
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
                var index = atoms.length;
                var serial = parseInt(tokens[0]);
                atom.serial = serial;
                // atom.serial = i;

                atom.x = parseFloat(tokens[2]);
                atom.y = parseFloat(tokens[3]);
                atom.z = parseFloat(tokens[4]);
                atom.atom = tokens[5];
                var charge = parseFloat(tokens[8]);

                atom.bonds = [];
                atom.bondOrder = [];
                atom.properties = {
                    'charge' : charge,
                    'partialCharge' : charge
                };
                serialToIndex[serial] = index;

                atoms.push(atom);
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
                fromAtom = atoms[serialToIndex[from]];
                var to = parseInt(tokens[2]);
                toAtom = atoms[serialToIndex[to]];

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

        return true;

    };

	var bondTable = {
		H :0.37,                                                                                                                                He:0.32,
		Li:1.34,Be:0.90,                                                                                B :0.82,C :0.77,N :0.75,O :0.73,F :0.71,Ne:0.69,
		Na:1.54,Mg:1.30,                                                                                Al:1.18,Si:1.11,P :1.06,S :1.02,Cl:0.99,Ar:0.97,
		K :1.96,Ca:1.74,Sc:1.44,Ti:1.56,V :1.25,/* Cr */Mn:1.39,Fe:1.25,Co:1.26,Ni:1.21,Cu:1.38,Zn:1.31,Ga:1.26,Ge:1.22,/* As */Se:1.16,Br:1.14,Kr:1.10,
		Rb:2.11,Sr:1.92,Y :1.62,Zr:1.48,Nb:1.37,Mo:1.45,Tc:1.56,Ru:1.26,Rh:1.35,Pd:1.31,Ag:1.53,Cd:1.48,In:1.44,Sn:1.41,Sb:1.38,Te:1.35,I :1.33,Xe:1.30,
		Cs:2.25,Ba:1.98,Lu:1.60,Hf:1.50,Ta:1.38,W :1.46,Re:1.59,Os:1.44,Ir:1.37,Pt:1.28,Au:1.44,Hg:1.49,Tl:1.48,Pb:1.47,Bi:1.46,/* Po *//* At */Rn:1.45,

		// None of the boottom row or any of the Lanthanides have bond lengths
	}
    var bondLength = function(elem) {
        return bondTable[elem] || 1.6;
    }
    // return true if atom1 and atom2 are probably bonded to each other
    // based on distance alone
    var areConnected = function(atom1, atom2) {
        var maxsq = bondLength(atom1.elem) + bondLength(atom2.elem);
        maxsq *= maxsq;
        maxsq *= 1.1; // fudge factor, especially important for md frames

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
        else
            return true;
    };

    //adds symmetry info to either duplicate and rotate/translate biological unit later or add extra atoms now
    var processSymmetries = function(format, copyMatrices, copyMatrix, atoms) {
        var end = atoms.length;
        var offset = end;
        var t, l, n;
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
                        if (format == "pdb") {
                            atoms.push({
                                'resn' : atoms[n].resn,
                                'x' : xyz.x,
                                'y' : xyz.y,
                                'z' : xyz.z,
                                'elem' : atoms[n].elem,
                                'hetflag' : atoms[n].hetflag,
                                'chain' : atoms[n].chain,
                                'resi' : atoms[n].resi,
                                'icode' : atoms[n].icode,
                                'rescode' : atoms[n].rescode,
                                'serial' : atoms[n].serial,
                                'atom' : atoms[n].atom,
                                'bonds' : bondsArr,
                                'ss' : atoms[n].ss,
                                'bondOrder' : atoms[n].bondOrder,
                                'properties' : atoms[n].properties,
                                'b' : atoms[n].b,
                                'pdbline' : atoms[n].pdbline,
                            });
                        }
                        else if (format == "mcif") {
                            atoms.push({
                                'resn' : atoms[n].resn,
                                'x' : xyz.x,
                                'y' : xyz.y,
                                'z' : xyz.z,
                                'elem' : atoms[n].elem,
                                'hetflag' : atoms[n].hetflag,
                                'chain' : atoms[n].chain,
                                'resi' : atoms[n].resi,
                                'serial' : atoms[n].serial,
                                'atom' : atoms[n].atom,
                                'bonds' : bondsArr,
                                'ss' : atoms[n].ss,
                                'bondOrder' : atoms[n].bondOrder,
                                'properties' : atoms[n].properties,
                                'hbondDistanceSq' : atoms[n].hbondDistanceSq,
                                'hbondOther' : atoms[n].hbondOther,
                                'ssbegin' : atoms[n].ssbegin,
                                'id' : atoms[n].id,
                                'index' : atoms[n].index
                            });
                        }
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


    // parse pdb file from str and create atoms
    // if computeStruct is true will always perform secondary structure
    // analysis,
    // otherwise only do analysis of SHEET/HELIX comments are missing
    /**
     * @param {AtomSpec[]}
     *            atoms
     * @param {string}
     *            str
     * @param {Object}
     *            options - keepH (do not strip hydrogens), noSecondaryStructure
     *            (do not compute ss)
     */
    parsers.pdb = parsers.PDB = parsers.pdbqt = parsers.PDBQT = function(atoms,
            str, options, modelData) {

        var atoms_cnt = 0;
        var noH = !options.keepH; // suppress hydrogens by default
        var computeStruct = !options.noSecondaryStructure;
        var noAssembly = !options.doAssembly; // don't assemble by default
        var copyMatrix = !options.duplicateAssemblyAtoms; //default true
        modelData.symmetries = [];
    
        var start = atoms.length;
        var atom;
        var protein = {
            sheet : [],
            helix : []
        }; // get secondary structure straight from pdb

        var hasStruct = false;
        var serialToIndex = []; // map from pdb serial to index in atoms
        var lines = str.split(/\r?\n/);
        var i, j, k, line;
        for (i = 0; i < lines.length; i++) {
            line = lines[i].replace(/^\s*/, ''); // remove indent
            var recordName = line.substr(0, 6);
            var startChain, startResi, endChain, endResi;
            
            if(recordName.indexOf("END") == 0) {
                break;
            }
            else if (recordName == 'ATOM  ' || recordName == 'HETATM') {
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
                        // minimal cleanup here - pymol likes to output
                        // duplicated conect records
                        var toindex = serialToIndex[to];
                        if (fromAtom.bonds[fromAtom.bonds.length - 1] != toindex) {
                            fromAtom.bonds.push(toindex);
                            fromAtom.bondOrder.push(1);
                        }
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

            }
        }

        var starttime = (new Date()).getTime();
        // assign bonds - yuck, can't count on connect records
        assignPDBBonds(atoms);
        // console.log("bond connecting " + ((new Date()).getTime() -
        // starttime));
        
        if (!noAssembly) {
            processSymmetries("pdb", modelData.symmetries, copyMatrix, atoms);
        }

        if (computeStruct || !hasStruct) {
            starttime = (new Date()).getTime();
            computeSecondaryStructure(atoms);
            // console.log("secondary structure " + ((new Date()).getTime() -
            // starttime));
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

    /**
     * Parse a pqr file from str and create atoms. A pqr file is assumed to be a
     * whitespace delimited PDB with charge and radius fields.
     * 
     * 
     * @param {AtomSpec[]}
     *            atoms
     * @param {string}
     *            str
     * @param {Object}
     *            options - noSecondaryStructure (do not compute ss)
     */
    parsers.pqr = parsers.PQR = function(atoms, str, options) {

        var atoms_cnt = 0;
        var start = atoms.length;
        var atom;
        var computeStruct = !options.noSecondaryStructure;

        var serialToIndex = []; // map from pdb serial to index in atoms
        var lines = str.split(/\r?\n/);
        var i, j, k, line;
        for (i = 0; i < lines.length; i++) {
            line = lines[i].replace(/^\s*/, ''); // remove indent
            var recordName = line.substr(0, 6);
            var startChain, startResi, endChain, endResi;
            if (recordName == 'ATOM  ' || recordName == 'HETATM') {
                // I would have liked to split based solely on whitespace, but
                // it seems that there is no guarantee that all the fields will
                // be filled out (e.g. the chain) so this doesn't work
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
                var fromAtom = atoms[serialToIndex[from]];
                for (j = 0; j < 4; j++) {
                    var to = parseInt(line.substr([ 11, 16, 21, 26 ][j], 5));
                    var toAtom = atoms[serialToIndex[to]];
                    if (fromAtom !== undefined && toAtom !== undefined) {
                        fromAtom.bonds.push(serialToIndex[to]);
                        fromAtom.bondOrder.push(1);
                    }
                }
            }
        }

        // assign bonds - yuck, can't count on connect records
        assignPDBBonds(atoms);
        if (computeStruct)
            computeSecondaryStructure(atoms);

        return true;
    };

    return parsers;
})();
