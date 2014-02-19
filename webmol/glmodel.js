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
		
		for ( var i = 0; i < protatoms.length; i++) {
			var ai = protatoms[i];

			for ( var j = i + 1; j < protatoms.length; j++) {
				var aj = protatoms[j];
				if(aj.chain != ai.chain)
					break;
				if (aj.resi - ai.resi > 1) // can't be connected
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
		var atoms = []
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
	var parsePDB = function(atoms, str, keepH) {

		var atoms_cnt = 0;
		var noH = !keepH; // suppres hydrogens by default
		var start = atoms.length;
		var protein = {
			sheet : [],
			helix : []
		}; // get secondary structure straight from pdb

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
					'bonds' : [],
					'bondOrder' : [],
					'properties' : {},
					'b' : b,
					'pdbline' : line
				});
			} else if (recordName == 'SHEET ') {
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
		
		var starttime = (new Date()).getTime();
		//computeSecondaryStructure(atoms);
		console.log("secondary structure " + ((new Date()).getTime() - starttime));

		
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

		// construct vertices around orgin for given radius, memoize results
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

						var n = new THREE.Vector3(vertex.x, vertex.y, vertex.z);
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

		var drawAtomSphere = function(atom, geo) {
			if (!atom.style.sphere)
				return;
			var style = atom.style.sphere;
			if (style.hidden)
				return;
				
			
			var geoGroup = geo.geometryChunks[geo.geometryChunks.length - 1];
													 			
			var color = atom.color;
			if (typeof (style.color) != "undefined")
				color = style.color;
			var C = WebMol.CC.color(color);

			var x, y;
			var radius = getRadiusFromStyle(atom, style);
			var vobj = sphereVertexCache.getVerticesForRadius(radius);
			var start = geo.vertices.length;

			// now add vertices and create faces at appropriate location
			var vertices = vobj.vertices;
			var normals = vobj.normals;

			var verticesRows = vobj.verticesRows;
			var h = verticesRows.length - 1;
			
			var face, normal;
			var n_vertices = 0;
			var color = [C.r, C.g, C.b];
			for (y = 0; y < h; y++) {
				var w = verticesRows[y].length - 1;
				for (x = 0; x < w; x++) {

					var v1 = verticesRows[y][x + 1] + start;
					var v2 = verticesRows[y][x] + start;
					var v3 = verticesRows[y + 1][x] + start;
					var v4 = verticesRows[y + 1][x + 1] + start;

					var n1 = normals[v1 - start];
					var n2 = normals[v2 - start];
					var n3 = normals[v3 - start];
					var n4 = normals[v4 - start];
					var face, norm;
					if (Math.abs(vertices[v1 - start].y) === radius) {
						//geo.faces.push(new THREE.Face3(v1, v3, v4, [ n1, n3, n4 ], C));
						n_vertices = 3;
						var vertex1 = new vertex(vertices[v1].x + atom.x, vertices[v1].y + atom.y, vertices[v1].z + atom.z);
						var vertex3 = new vertex(vertices[v3].x + atom.x, vertices[v3].y + atom.y, vertices[v3].z + atom.z);
						var vertex4 = new vertex(vertices[v4].x + atom.x, vertices[v4].y + atom.y, vertices[v4].z + atom.z);
						face = [vertex1, vertex3, vertex4];
						norm = [n1, n3, n4];
					} else if (Math.abs(vertices[v3 - start].y) === radius) {
						//geo.faces.push(new THREE.Face3(v1, v2, v3, [ n1, n2, n3 ], C));
						n_vertices = 3;
						var vertex1 = new vertex(vertices[v1].x + atom.x, vertices[v1].y + atom.y, vertices[v1].z + atom.z);
						var vertex2 = new vertex(vertices[v2].x + atom.x, vertices[v2].y + atom.y, vertices[v2].z + atom.z);
						var vertex3 = new vertex(vertices[v3].x + atom.x, vertices[v3].y + atom.y, vertices[v3].z + atom.z);
						face = [vertex1, vertex2, vertex3];
						//face = [geo.vertices[v1], geo.vertices[v2], geo.vertices[v3]];
						norm = [n1, n2, n3];
					} else {
						//geo.faces.push(new THREE.Face4(v1, v2, v3, v4, [ n1, n2, n3, n4 ], C));
						n_vertices = 4;
						var vertex1 = new vertex(vertices[v1].x + atom.x, vertices[v1].y + atom.y, vertices[v1].z + atom.z);
						var vertex2 = new vertex(vertices[v2].x + atom.x, vertices[v2].y + atom.y, vertices[v2].z + atom.z);
						var vertex3 = new vertex(vertices[v3].x + atom.x, vertices[v3].y + atom.y, vertices[v3].z + atom.z);
						var vertex4 = new vertex(vertices[v4].x + atom.x, vertices[v4].y + atom.y, vertices[v4].z + atom.z);
						face = [vertex1, vertex2, vertex3, vertex4];
						//face = [geo.vertices[v1], geo.vertices[v2], geo.vertices[v3], geo.vertices[v4]];
						norm = [n1, n2, n3, n4];
					}
					
					//make new group if necessary
					if ( geoGroup.vertices + n_vertices > 65535 ) {
						geo.geometryChunks.push( new geometryChunk() );
						geoGroup = geo.geometryChunks[ geo.geometryChunks.length - 1];
					}
					
					geoGroup.vertices += n_vertices;
					var offset = geoGroup.vertexArr.length;
					
					//populate non-typed vertexArr and face arrays with vertex coordinates
					populateGroup(geoGroup, face, norm, color, offset);
				}
			}

		};
		var vertex = function(x, y, z) {
			this.x = x;
			this.y = y;
			this.z = z;
		};
		//Fill up geometry chunk's non typed arrays for a single face
		var populateGroup = function(group, face, norm, color, index) {
			
			//one triangle
			if (face.length === 3) {
				var v1 = face[ 0 ];
				var v2 = face[ 1 ];
				var v3 = face[ 2 ];
				var verts = [v1.x, v1.y, v1.z,
							 v2.x, v2.y, v2.z,
							 v3.x, v3.y, v3.z];
							 
				var n1 = norm[0];
				var n2 = norm[1];
				var n3 = norm[2];
				var norms = [n1.x, n1.y, n1.z,
							 n2.x, n2.y, n2.z,
							 n3.x, n3.y, n3.z];
				var r = color[0];
				var g = color[1];
				var b = color[2];
							 
				//group.vertexArr = group.vertexArr.concat(verts);
				for (var i in verts) {
					group.vertexArr.push(verts[i]);
				}

				for (var i in norms) {
					group.normalArr.push(norms[i]);
				}
				
				group.normalArr = group.normalArr.concat(norms);
				
				group.faceArr.push(index); 
				group.faceArr.push(index + 1);
				group.faceArr.push(index + 2);
				
				//TODO: do I really need these arrays?
				group.lineArr.push(index);
				group.lineArr.push(index + 1);
				group.lineArr.push(index);
				group.lineArr.push(index + 2);
				group.lineArr.push(index + 1);
				group.lineArr.push(index + 2);
				
				//add rgb color values for each vertex - faces are always same color
				//group.colorArr = group.colorArr.concat(color);
				//group.colorArr = group.colorArr.concat(color);
				//group.colorArr = group.colorArr.concat(color);
				group.colorArr.push(r);
				group.colorArr.push(g);
				group.colorArr.push(b);
				
				group.colorArr.push(r);
				group.colorArr.push(g);
				group.colorArr.push(b);
				
				group.colorArr.push(r);
				group.colorArr.push(g);
				group.colorArr.push(b);

			}
			
			else if (face.length === 4) {
				v1 = face[ 0 ];
				v2 = face[ 1 ];
				v3 = face[ 2 ];
				v4 = face[ 3 ];
				
				var verts = [v1.x, v1.y, v1.z,
							 v2.x, v2.y, v2.z,
							 v3.x, v3.y, v3.z,
							 v4.x, v4.y, v4.z];
							 
				var n1 = norm[0];
				var n2 = norm[1];
				var n3 = norm[2];
				var n4 = norm[3];
				var norms = [n1.x, n1.y, n1.z,
							 n2.x, n2.y, n2.z,
							 n3.x, n3.y, n3.z,
							 n4.x, n4.y, n4.z];
				
				var r = color[0];
				var g = color[1];
				var b = color[2];

				//group.vertexArr = group.vertexArr.concat(verts);
				for (var i in verts) {
					group.vertexArr.push(verts[i]);
				}
				
				//group.normalArr = group.normalArr.concat(norms);
				for (var i in norms) {
					group.normalArr.push(norms[i]);
				}
				
				group.faceArr.push(index);
				group.faceArr.push(index + 1);
				group.faceArr.push(index + 3);
				
				group.faceArr.push(index + 1);
				group.faceArr.push(index + 2);
				group.faceArr.push(index + 3);
				
				group.lineArr.push(index);
				group.lineArr.push(index + 1);
				group.lineArr.push(index);
				group.lineArr.push(index + 3);
				group.lineArr.push(index + 1);
				group.lineArr.push(index + 2);
				group.lineArr.push(index + 2);
				group.lineArr.push(index + 3);
				
				//group.colorArr = group.colorArr.concat(color);
				//group.colorArr = group.colorArr.concat(color);
				//group.colorArr = group.colorArr.concat(color);
				//group.colorArr = group.colorArr.concat(color);
				group.colorArr.push(r);
				group.colorArr.push(g);
				group.colorArr.push(b);
				
				group.colorArr.push(r);
				group.colorArr.push(g);
				group.colorArr.push(b);
				
				group.colorArr.push(r);
				group.colorArr.push(g);
				group.colorArr.push(b);		

				group.colorArr.push(r);
				group.colorArr.push(g);
				group.colorArr.push(b);			
			};
		};
		
		//Initialize typed array buffers for completed geometry
		var initBuffers = function(geometry) {
			
			for (var i = 0; i < geometry.geometryChunks.length; ++i){
				group = geometry.geometryChunks[i];
				group.__vertexArray = new Float32Array(group.vertexArr);
				group.__colorArray = new Float32Array(group.colorArr);
				group.__normalArray = new Float32Array(group.normalArr);
				group.__faceArray = new Uint16Array(group.faceArr);
				group.__lineArray = new Uint16Array(group.lineArr);
				
				//Doesn't free memory directly, but should break references for gc 
				delete group.vertexArr;
				delete group.colorArr;
				delete group.normalArr;
				delete group.faceArr;
				delete group.lineArr;
			}
		};
		
		//represents individual renderable geometry group
		var geometryChunk = function() {
			this.vertexArr = [];
			this.colorArr = [];
			this.normalArr = [];
			this.faceArr = [];
			this.lineArr = [];
			this.vertices = 0;
		};
		// cross drawing
		var drawAtomCross = function(atom, geos) {
			if (!atom.style.cross)
				return;
			var style = atom.style.cross;
			if (style.hidden)
				return;
			var linewidth = (atom.style.cross.lineWidth || defaultlineWidth);
			if (!geos[linewidth])
				geos[linewidth] = new THREE.Geometry();
			var geo = geos[linewidth];

			var delta = getRadiusFromStyle(atom, style);

			var points = [ [ delta, 0, 0 ], [ -delta, 0, 0 ], [ 0, delta, 0 ],
					[ 0, -delta, 0 ], [ 0, 0, delta ], [ 0, 0, -delta ] ];

			var c = WebMol.CC.color(atom.color);
			for ( var j = 0; j < 6; j++) {
				geo.vertices.push(new TV3(atom.x + points[j][0], atom.y
						+ points[j][1], atom.z + points[j][2]));
				geo.colors.push(c);
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
			var linewidth = (atom.style.line.lineWidth || defaultlineWidth);

			if (!geos[linewidth])
				geos[linewidth] = new THREE.Geometry();
			var geo = geos[linewidth];

			for ( var i = 0; i < atom.bonds.length; i++) {
				var j = atom.bonds[i]; // our neighbor
				if (i < j) {// only draw if less
					// TODO: handle bond orders
					var atom2 = atoms[j];
					if (!atom2.style.line)
						continue; // don't sweat the details
					var vs = geo.vertices, cs = geo.colors;
					var p1 = new TV3(atom.x, atom.y, atom.z);
					var p2 = new TV3(atom2.x, atom2.y, atom2.z);
					var mp = p1.clone().add(p2).multiplyScalar(0.5);

					var c1 = WebMol.CC.color(atom.color), c2 = WebMol.CC
							.color(atom2.color);

					if (typeof (style.color) != "undefined") {
						c1 = c2 = WebMol.CC.color(style.color);
					}
					vs.push(p1);
					cs.push(c1);
					vs.push(mp);
					cs.push(c1);
					vs.push(mp);
					cs.push(c2);
					vs.push(p2);
					cs.push(c2);
				}
			}

		};

		// bonds as cylinders
		var defaultStickRadius = .25;
		var cylinderQuality = 12;
		var cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1,
				cylinderQuality, 1, true);
		cylinderGeometry.faceUvs = []; // null these out to make merging faster
		cylinderGeometry.faceVertexUvs = [ [] ];

		// creates a cylinder
		// TODO: create it ourselves in the hopes of getting a speed up
		var drawnC = 0;
		var drawCylinder = function(geo, from, to, radius, color) {
			if (!from || !to)
				return;
			drawnC++;
			// vertices

			var dir = to.clone();
			dir.sub(from);

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

			var start = geo.vertices.length;
			// add vertices, opposing vertices paired together
			for ( var i = 0, n = nvecs.length; i < n; ++i) {
				var bottom = nvecs[i].clone().multiplyScalar(radius).add(from);
				var top = nvecs[i].clone().multiplyScalar(radius).add(to);

				geo.vertices.push(bottom);
				geo.vertices.push(top);
			}

			// now faces
			var face;
			for ( var i = 0, n = nvecs.length - 1; i < n; ++i) {
				var ti = start + 2 * i;
				face = new THREE.Face4(ti, ti + 1, ti + 3, ti + 2);
				face.color = color;
				face.vertexNormals = [ nvecs[i], nvecs[i], nvecs[i + 1],
						nvecs[i + 1] ];
				geo.faces.push(face);
			}
			// final face

			face = new THREE.Face4(start + 30, start + 31, start + 1, start);
			face.color = color;
			face.vertexNormals = [ nvecs[15], nvecs[15], nvecs[0], nvecs[0] ];
			geo.faces.push(face);

		};

		// draws cylinders and small spheres (at bond radius)
		var drawBondSticks = function(atom, atoms, geo) {
			if (!atom.style.stick)
				return;
			var style = atom.style.stick;
			if (style.hidden)
				return;

			var bondR = style.radius || defaultStickRadius;

			var c1 = atom.color;
			if (typeof (style.color) != "undefined") {
				c1 = style.color;
			}
			var C1 = WebMol.CC.color(c1);

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

					var p1 = new TV3(atom.x, atom.y, atom.z);
					var p2 = new TV3(atom2.x, atom2.y, atom2.z);

					var c2 = atom2.color;
					if (typeof (style.color) != "undefined") {
						c2 = style.color;
					}
					var C2 = WebMol.CC.color(c2);

					// draw cylinders
					if (atom.bondOrder[i] == 1) {
						if (c1 != c2) {
							var mp = new TV3().addVectors(p1, p2)
									.multiplyScalar(0.5);
							drawCylinder(geo, p1, mp, bondR, C1);
							drawCylinder(geo, mp, p2, bondR, C2);
						} else {
							drawCylinder(geo, p1, p2, bondR, C1);
						}
					} else if (atom.bondOrder[i] > 1) {
						var dir = p2.clone();
						var v = null;
						dir.sub(p1);

						if (atom.bonds.length == 1) {
							if (atom2.bonds.length == 1) {
								v = dir.clone();
								if (Math.abs(v.x) > .0001)
									v.y += 1;
								else
									v.x += 1;
							} else {
								var i2 = (i + 1) % atom2.bonds.length;
								var j2 = atom2.bonds[i2];
								var atom3 = atoms[j2];
								var p3 = new TV3(atom3.x, atom3.y, atom3.z);

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
							var p3 = new TV3(atom3.x, atom3.y, atom3.z);

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
								var mp = new TV3().addVectors(p1a, p2a)
										.multiplyScalar(0.5);
								var mp2 = new TV3().addVectors(p1b, p2b)
										.multiplyScalar(0.5);
								drawCylinder(geo, p1a, mp, r, C1);
								drawCylinder(geo, mp, p2a, r, C2);
								drawCylinder(geo, p1b, mp2, r, C1);
								drawCylinder(geo, mp2, p2b, r, C2);
							} else {
								drawCylinder(geo, p1a, p2a, r, C1);
								drawCylinder(geo, p1b, p2b, r, C1);
							}
						} else if (atom.bondOrder[i] == 3) {
							var r = bondR / 4;
							v.cross(dir);
							v.normalize();
							v.multiplyScalar(r * 3);

							var p1a = p1.clone();
							p1a.add(v);
							var p1b = p1.clone()
							p1b.sub(v);

							var p2a = p1a.clone();
							p2a.add(dir);
							var p2b = p1b.clone();
							p2b.add(dir);

							if (c1 != c2) {
								var mp = new TV3().addVectors(p1a, p2a)
										.multiplyScalar(0.5);
								var mp2 = new TV3().addVectors(p1b, p2b)
										.multiplyScalar(0.5);
								var mp3 = new TV3().addVectors(p1, p2)
										.multiplyScalar(0.5);
								drawCylinder(geo, p1a, mp, r, C1);
								drawCylinder(geo, mp, p2a, r, C2);
								drawCylinder(geo, p1, mp3, r, C1);
								drawCylinder(geo, mp3, p2, r, C2);
								drawCylinder(geo, p1b, mp2, r, C1);
								drawCylinder(geo, mp2, p2b, r, C2);
							} else {
								drawCylinder(geo, p1a, p2a, r, C1);
								drawCylinder(geo, p1, p2, r, C1);
								drawCylinder(geo, p1b, p2b, r, C1);

							}
						}
					}
				}
			}

			// for junctions draw sphere; merge the sphere geometries was really
			// really slow
			var savedstyle = atom.style;
			atom.style = {
				sphere : {
					radius : bondR,
					color : c1
				}
			};
			drawAtomSphere(atom, geo);
			atom.style = savedstyle;

		};

		// go through all the atoms and regenerate their geometries
		// we try to have one geometry for each style since this is much much
		// faster
		// at some point we should optimize this to avoid unnecessary
		// recalculation
		var createMolObj = function(atoms) {


			console.log("creating for "+id);
			var ret = new THREE.Object3D();
			var cartoonAtoms = [];
			var lineGeometries = {};
			var crossGeometries = {};
			var sphereGeometry = new THREE.Geometry();
			
			sphereGeometry.geometryChunks = [];
			sphereGeometry.geometryChunks.push( new geometryChunk() );	
									   			  	
			var stickGeometry = new THREE.Geometry();
			for ( var i = 0; i < atoms.length; i++) {
				var atom = atoms[i];
				// recreate gl info for each atom as necessary
				if (atom && atom.style) {
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
			}

			// add sphere geometry
			if (sphereGeometry.geometryChunks[0].vertices) {
				var sphereMaterial = new THREE.MeshLambertMaterial({
					ambient : 0x000000,
					vertexColors : true,
					reflectivity : 0
				});
				
				//Initialize buffers in geometry
				initBuffers(sphereGeometry);
				
				var sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
				console
						.log("sphere geometry "
								+ sphereGeometry.vertices.length);

				ret.add(sphere);
			}

			// add stick geometry
			if (stickGeometry.vertices && stickGeometry.vertices.length > 0) {
				var cylinderMaterial = new THREE.MeshLambertMaterial({
					vertexColors : true,
					ambient : 0x000000,
					reflectivity : 0
				});

				var sticks = new THREE.Mesh(stickGeometry, cylinderMaterial);
				ret.add(sticks);
			}

			// add any line geometries, distinguished by line width
			for ( var i in lineGeometries) {
				if (lineGeometries.hasOwnProperty(i)) {
					var linewidth = i;
					var lineMaterial = new THREE.LineBasicMaterial({
						linewidth : linewidth,
						vertexColors : true,
					});

					var line = new THREE.Line(lineGeometries[i], lineMaterial,
							THREE.LinePieces);

					ret.add(line);
				}
			}

			// add any cross geometries
			for ( var i in crossGeometries) {
				if (crossGeometries.hasOwnProperty(i)) {
					var linewidth = i;
					var lineMaterial = new THREE.LineBasicMaterial({
						linewidth : linewidth,
						vertexColors : true,
					});

					var line = new THREE.Line(crossGeometries[i], lineMaterial,
							THREE.LinePieces);

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
				}
			}
		};

		// add atoms to this model from molecular data string
		this.addMolData = function(data, format) {

			switch (format) {
			case "xyz":
				parseXYZ(atoms, data);
				break;
			case "pdb":
				parsePDB(atoms, data);
				break;
			case "sdf":
				parseSDF(atoms, data);
				break;
			case "mol2":
				parseMOL2(atoms, data);
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
			
			if(molObj != null && sameObj(colors,lastColors))
				return; // don't recompute
			lastColors = colors;
			var atoms = this.selectedAtoms(sel);
			if(atoms.length > 0)
				molObj = null; // force rebuild
			for ( var i = 0; i < atoms.length; i++) {
				var a = atoms[i];
				if(typeof(colors[a.elem]) != "undefined") {
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
			if(molObj == null) { // have to regenerate
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
				group.remove(renderedMolObj);
				renderedMolObj = null;
			}
			molObj = null;
		};

	};

	return GLModel;
})();
