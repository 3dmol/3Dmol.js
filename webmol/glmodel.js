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

	var defaultlineWidth = 1.5;

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

	// given a selection specification, return true if atom is selected
	var atomIsSelected = function(atom, sel) {

		if (typeof (sel) === "undefined")
			return true; // undef gets all
		for ( var key in sel) {
			if (sel.hasOwnProperty(key) && key != "props") {
				// if something is in sel, atom must have it
				if (typeof (atom[key]) === "undefined")
					return false;
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
					if (!isokay)
						return false;
				} else { // single match
					if (atom[key] != sel[key])
						return false;
				}
			}
		}
		return true;
	}
	// return true if atom1 and atom2 are probably bonded to each other
	// based on distance alone
	var areConnected = function(atom1, atom2) {
		var max = 3.42;

		var xdiff = atom1.x - atom2.x;
		if (xdiff > max)
			return false;
		var ydiff = atom1.y - atom2.y;
		if (ydiff > max)
			return false;
		var zdiff = atom1.z - atom2.z;
		if (zdiff > max)
			return false;

		var distSquared = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;

		// if (atom1.altLoc != atom2.altLoc) return false;
		if (isNaN(distSquared))
			return 0;
		if (distSquared < 0.5)
			return 0; // maybe duplicate position.

		if (distSquared > 1.3
				&& (atom1.elem == 'H' || atom2.elem == 'H' || atom1.elem == 'D' || atom2.elem == 'D'))
			return 0;
		if (distSquared < 3.6 && (atom1.elem == 'S' || atom2.elem == 'S'))
			return 1;
		if (distSquared > 2.78)
			return 0;
		return 1;
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
		for ( var i = start; i < end; i++)
			// n^2 behavior.. should at least sort by z
			for ( var j = i + 1; j < end; j++)
				if (areConnected(atoms[i], atoms[j])) {
					atoms[i].bonds.push(j);
					atoms[i].bondOrder.push(1);
					atoms[j].bonds.push(i);
					atoms[j].bondOrder.push(1);
				}

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

	// parse pdb file from str and create atoms
	var parsePDB = function(atoms, str) {

		var atoms_cnt = 0;
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
					if (isNaN(to))
						continue;
					var toAtom = atoms[serialToIndex[to]];
					if (fromAtom != undefined) {
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
		for ( var i = start; i < atoms.length; i++) {
			// n^2 behavior.. should at least sort by z
			for ( var j = i + 1; j < atoms.length; j++) {
				if (areConnected(atoms[i], atoms[j])) {
					if (atoms[i].bonds.indexOf(j) == -1) {
						// only add if not already there
						atoms[i].bonds.push(j);
						atoms[i].bondOrder.push(1);
						atoms[j].bonds.push(i);
						atoms[j].bondOrder.push(1);
					}
				}
			}
		}
		console.log("bond connecting " + ((new Date()).getTime() - starttime));
		// Assign secondary structures
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
		}

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
			for ( var i = 0, n = vertices.length; i < n; i++) {
				var v = vertices[i];
				var vertex = new THREE.Vector3(v.x + atom.x, v.y + atom.y, v.z
						+ atom.z);
				geo.vertices.push(vertex);
				geo.normals.push(normals[i]);
			}
			var verticesRows = vobj.verticesRows;
			var h = verticesRows.length - 1;
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

					if (Math.abs(vertices[v1 - start].y) === radius) {
						geo.faces.push(new THREE.Face3(v1, v3, v4,
								[ n1, n3, n4 ], C));
					} else if (Math.abs(vertices[v3 - start].y) === radius) {
						geo.faces.push(new THREE.Face3(v1, v2, v3,
								[ n1, n2, n3 ], C));
					} else {
						geo.faces.push(new THREE.Face4(v1, v2, v3, v4, [ n1,
								n2, n3, n4 ], C));
					}
				}
			}

		};

		// cross drawing
		var drawAtomCross = function(atom, geos) {
			if (!atom.style.cross)
				return;
			var style = atom.style.cross;
			if (style.hidden)
				return;
			var linewidth = (atom.style.lineWidth || defaultlineWidth)
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
			var linewidth = (atom.style.lineWidth || defaultlineWidth)
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

					var c1 = WebMol.CC.color(atom.color), c2 = WebMol.CC.color(atom2.color);

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
			drawnC ++;
			// vertices
			
			var dir = to.clone();
			dir.sub(from);
			
			//get orthonormal vector
			var nvecs = [];
			nvecs[0] = dir.clone();
			if(Math.abs(nvecs[0].x) > .0001) nvecs[0].y += 1;
			else nvecs[0].x += 1;
			nvecs[0].cross(dir);
			nvecs[0].normalize();
			
			nvecs[0] = nvecs[0];
			//another orth vector
			nvecs[4] = nvecs[0].clone();
			nvecs[4].crossVectors(nvecs[0],dir);
			nvecs[4].normalize();
			nvecs[8] = nvecs[0].clone().negate();
			nvecs[12] = nvecs[4].clone().negate();
			
			//now quarter positions
			nvecs[2] = nvecs[0].clone().add(nvecs[4]).normalize();
			nvecs[6] = nvecs[4].clone().add(nvecs[8]).normalize();
			nvecs[10] = nvecs[8].clone().add(nvecs[12]).normalize();
			nvecs[14] = nvecs[12].clone().add(nvecs[0]).normalize();
			
			//eights
			nvecs[1] = nvecs[0].clone().add(nvecs[2]).normalize();
			nvecs[3] = nvecs[2].clone().add(nvecs[4]).normalize();
			nvecs[5] = nvecs[4].clone().add(nvecs[6]).normalize();
			nvecs[7] = nvecs[6].clone().add(nvecs[8]).normalize();
			nvecs[9] = nvecs[8].clone().add(nvecs[10]).normalize();
			nvecs[11] = nvecs[10].clone().add(nvecs[12]).normalize();
			nvecs[13] = nvecs[12].clone().add(nvecs[14]).normalize();
			nvecs[15] = nvecs[14].clone().add(nvecs[0]).normalize();
			
			var start = geo.vertices.length;
			//add vertices, opposing vertices paired together
			for(var i = 0, n = nvecs.length; i < n; ++i) {
				var bottom = nvecs[i].clone().multiplyScalar(radius).add(from);
				var top = nvecs[i].clone().multiplyScalar(radius).add(to);
				
				geo.vertices.push(bottom);
				geo.vertices.push(top);				
			}
			
			//now faces
			var face;
			for(var i = 0, n = nvecs.length-1; i < n; ++i) {
				var ti = start+2*i;
				face = new THREE.Face4(ti,ti+1,ti+3,ti+2);
				face.color = color;
				face.vertexNormals = [ nvecs[i], nvecs[i], nvecs[i+1], nvecs[i+1]];
				geo.faces.push(face);
			}
			//final face

			face = new THREE.Face4(start+30,start+31,start+1,start);
			face.color = color;
			face.vertexNormals = [ nvecs[15], nvecs[15], nvecs[0], nvecs[0]];
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

			for ( var i = 0; i < atom.bonds.length; i++) {
				var j = atom.bonds[i]; // our neighbor
				var atom2 = atoms[j];
				if (atom.serial < atom2.serial) {// only draw if less, this lets use combine cylinders of the same color
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
					if (c1 != c2) {
						var mp = new TV3().addVectors(p1, p2).multiplyScalar(0.5);
						drawCylinder(geo, p1, mp, bondR, C1);
						drawCylinder(geo, mp, p2, bondR, C2);
					} else {
						drawCylinder(geo, p1, p2, bondR, C1);
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
			var ret = new THREE.Object3D();
			var cartoonAtoms = [];
			var lineGeometries = {};
			var crossGeometries = {};
			var sphereGeometry = new THREE.Geometry();
			var stickGeometry = new THREE.Geometry();
			for ( var i = 0; i < atoms.length; i++) {
				var atom = atoms[i];
				// recreate gl info for each atom as necessary
				if (atom && atom.style) {
					drawAtomSphere(atom, sphereGeometry);
					drawAtomCross(atom, crossGeometries);
					drawBondLines(atom, atoms, lineGeometries);
					drawBondSticks(atom, atoms, stickGeometry);
					if (typeof (atom.style.cartoon) != "undefined"
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
			if (sphereGeometry.vertices && sphereGeometry.vertices.length > 0) {
				var sphereMaterial = new THREE.MeshLambertMaterial({
					ambient : 0x000000,
					vertexColors : true,
					reflectivity : 0
				});
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
				console.log("stick geometry " + stickGeometry.vertices.length);
				console.log("drawC " + drawnC);
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
					atom.globj = null;
				}
			}
		}

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
			}
			setAtomDefaults(atoms, id);
		};

		// return list of atoms selected by sel, this is specific to glmodel
		this.selectedAtoms = function(sel) {
			var ret = [];
			for ( var i = 0; i < atoms.length; i++) {
				var atom = atoms[i];
				if (atom) {
					if (atomIsSelected(atom, sel))
						ret.push(atom);
				}
			}
			return ret;
		}

		// style the select atoms with style
		this.setStyle = function(style, sel) {
			var atoms = this.selectedAtoms(sel);
			// do a copy to enforce style changes through this function
			var mystyle = $.extend(true, {}, style);

			// somethings we only calculate if there is a change in a certain
			// style, although these checks will only catch cases where both
			// are either null or undefined
			for ( var i = 0; i < atoms.length; i++) {

				atoms[i].style = mystyle;
				atoms[i].globj = null; // need to recalculate

				if (atoms[i].style.line != mystyle.line
						|| atoms[i].style.stick != mystyle.stick) {
					var bonds = atoms[i].bonds;
					for ( var i = 0; i < bonds.length; i++) {
						var atomj = atoms[bonds[i]];
						atomj.globj = null;
					}
				}
			}

		};
		
		//return pdb output of selected atoms
		//currently only works if input was pdb
		this.pdbData = function(style, sel) {
			var atoms = this.selectedAtoms(sel);
			var ret = "";
			for( var i = 0, n = atoms.length; i < n; ++i) {
				ret += atoms[i].pdbline;
			}
			return ret;
		}

		// return 3d data for this model, this is specific to glmodel
		this.globj = function() {
			var time = new Date();
			molObj = createMolObj(atoms);
			var time2 = new Date();
			console.log("object creation time: " + (time2 - time));
			return molObj;
		};

	}
	;

	return GLModel;
})();