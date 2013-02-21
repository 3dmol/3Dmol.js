// A model is a collection of related atoms.  Bonds are only allowed between
//atoms in the same model.  An atom is uniquely specified by its model id and
//its serial number.
//A glmodel knows how to apply the styles on each atom to create a gl object

var WebMol = WebMol || {};

WebMol.GLModel = (function() {
	// class variables go here
	var defaultAtomStyle = {
		sphere : null,
		icosahedron : null,
		stick : null,
		line : {},
		cross : null,
		cartoon : null
	};

	var Nucleotides = [ '  G', '  A', '  T', '  C', '  U', ' DG', ' DA', ' DT',
			' DC', ' DU' ];
	var ElementColors = {
		"H" : 0xCCCCCC,
		"C" : 0xAAAAAA,
		"O" : 0xCC0000,
		"N" : 0x0000CC,
		"S" : 0xCCCC00,
		"P" : 0x6622CC,
		"F" : 0x00CC00,
		"CL" : 0x00CC00,
		"BR" : 0x882200,
		"I" : 0x6600AA,
		"FE" : 0xCC6600,
		"CA" : 0x8888AA
	};
	var defaultColor = 0xCCCCCC;
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

		if(typeof(sel) === "undefined") return true; //undef gets all
		for ( var key in sel) {
			if (sel.hasOwnProperty(key)) {
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

		var distSquared = (atom1.x - atom2.x) * (atom1.x - atom2.x)
				+ (atom1.y - atom2.y) * (atom1.y - atom2.y)
				+ (atom1.z - atom2.z) * (atom1.z - atom2.z);

		// if (atom1.altLoc != atom2.altLoc) return false;
		if (isNaN(distSquared))
			return 0;
		if (distSquared < 0.5)
			return 0; // maybe duplicate position.

		if (distSquared > 1.3
				&& (atom1.elem == 'H' || atom2.elem == 'H' || atom1.elem == 'D' || atom2.elem == 'D'))
			return 0;
		if (distSquared < 3.42 && (atom1.elem == 'S' || atom2.elem == 'S'))
			return 1;
		if (distSquared > 2.78)
			return 0;
		return 1;
	};

	// set default style and colors for atoms
	var setAtomDefaults = function(atoms, id) {
		for ( var i = 0; i < atoms.length; i++) {
			var atom = atoms[i];
			atom.style = atom.style || defaultAtomStyle;
			atom.color = atom.color || ElementColors[atom.elem] || defaultColor;
			atom.model = id;
			atom.globj = null;
		}
	}

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
	
	//put atoms specified in sdf fromat in str into atoms
	//adds to atoms, does not replace
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
		var end = start+atomCount;
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
			atoms[i] = atom;
		}
		for (i = 0; i < bondCount; i++) {
			var line = lines[offset];
			offset++;
			var from = parseInt(line.substr(0, 3))-1+start;
			var to = parseInt(line.substr(3, 3))-1+start;
			var order = parseInt(line.substr(6, 3));
			atoms[from].bonds.push(to);
			atoms[from].bondOrder.push(order);
			atoms[to].bonds.push(from);
			atoms[to].bondOrder.push(order);
		}

		return true;
	};

	// sphere drawing
	var defaultSphereRadius = 1.5;
	var sphereQuality = 16; // 16;
	var sphereGeometry = new THREE.SphereGeometry(1, sphereQuality,
			sphereQuality);
	
	var drawAtomSphere = function(atom) {
		if (!atom.style.sphere)
			return;
		var style = atom.style.sphere;
		if (style.hidden)
			return;
		
		var color = atom.color;
		if(typeof(style.color) != "undefined")
			color = style.color;
		var sphereMaterial = new THREE.MeshLambertMaterial({
			color : color
		});
		var sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
		var r = defaultSphereRadius;
		if (typeof (style.fixedRadius) != "undefined")
			r = style.fixedRadius;
		else if (vdwRadii[atom.elem])
			r = vdwRadii[atom.elem];

		if (typeof (style.scale) != "undefined")
			r *= style.scale;

		sphere.scale.x = sphere.scale.y = sphere.scale.z = r;
		sphere.position.x = atom.x;
		sphere.position.y = atom.y;
		sphere.position.z = atom.z;
		
		atom.globj = atom.globj || new THREE.Object3D();
		atom.globj.add(sphere);
	}

	// icosohedral drawing
	var icosahedronGeo = new THREE.IcosahedronGeometry(1);
	var drawAtomIcosahedron = function(atom) {
		if (!atom.style.icosahedron)
			return;
		var style = atom.style.icosahedron;
		if (style.hidden)
			return;

		var color = atom.color;
		if(typeof(style.color) != "undefined")
			color = style.color;
		var mat = new THREE.MeshLambertMaterial({
			color : color
		});
		var sphere = new THREE.Mesh(icosahedronGeo, mat);
		
		var r = defaultSphereRadius;
		if (typeof (style.fixedRadius) != "undefined")
			r = style.fixedRadius;
		else if (vdwRadii[atom.elem])
			r = vdwRadii[atom.elem];

		if (typeof (style.scale) != "undefined")
			r *= style.scale;

		sphere.scale.x = sphere.scale.y = sphere.scale.z = r;
		sphere.position.x = atom.x;
		sphere.position.y = atom.y;
		sphere.position.z = atom.z;
		
		atom.globj = atom.globj || new THREE.Object3D();
		atom.globj.add(sphere);
	};
	
	//cross drawing
	var drawAtomCross = function(atom) {
		if (!atom.style.cross)
			return;
		var style = atom.style.cross;
		if (style.hidden)
			return;
		var geo = new THREE.Geometry();
		var delta = style.delta || 0.5;
		var points = [ [ delta, 0, 0 ], [ -delta, 0, 0 ], [ 0, delta, 0 ],
				[ 0, -delta, 0 ], [ 0, 0, delta ], [ 0, 0, -delta ] ];

		var c = new TCo(atom.color);
		for ( var j = 0; j < 6; j++) {
			geo.vertices.push(new TV3(atom.x + points[j][0], atom.y
						+ points[j][1], atom.z + points[j][2]));
				geo.colors.push(c);
		}
		
		var lineMaterial = new THREE.LineBasicMaterial({
			linewidth : (style.linewidth || defaultlineWidth)
		});
		lineMaterial.vertexColors = true;
		var line = new THREE.Line(geo, lineMaterial, THREE.LinePieces);
		atom.globj = atom.globj || new THREE.Object3D();
		atom.globj.add(line);
	};
	
	//bonds - both atoms must match bond style
	//standardize on only drawing for lowest to highest
	
	var drawBondLines = function(atom, atoms) {
		if (!atom.style.line)
			return;
		var style = atom.style.line;
		if (style.hidden)
			return;
		
		var geo = new THREE.Geometry();

		for(var i = 0; i < atom.bonds.length; i++) {
			var j = atom.bonds[i]; //our neighbor
			if(i < j) {//only draw if less
				//TODO: handle bond orders
				var atom2 = atoms[j];
				if(!atom2.style.line)
					continue; //don't sweat the details
				var vs = geo.vertices, cs = geo.colors;
				var p1 = new TV3(atom.x, atom.y, atom.z);
				var p2 = new TV3(atom2.x, atom2.y, atom2.z);
				var mp = p1.clone().add(p2).multiplyScalar(0.5);

				var c1 = new TCo(atom.color), c2 = new TCo(atom2.color);
				vs.push(p1);
				cs.push(c1);
				vs.push(mp);
				cs.push(c1);
				vs.push(p2);
				cs.push(c2);
				vs.push(mp);
				cs.push(c2);	
			}
		}
		
		var lineMaterial = new THREE.LineBasicMaterial({
			linewidth : (style.lineWidth || defaultlineWidth)
		});
		
		if(!style.color)
			lineMaterial.vertexColors = true;
		else
			lineMaterial.color = new TCo(style.color);

		var line = new THREE.Line(geo, lineMaterial);
		line.type = THREE.LinePieces;
		atom.globj = atom.globj || new THREE.Object3D();
		atom.globj.add(line);		
	};
	
	//bonds as cylinders
	var defaultStickRadius = .25;
	var cylinderQuality = 16;
	var cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1,
			cylinderQuality, 1, false);
	cylinderGeometry.faceUvs = [];
	var faceVertexUvs = [];
	
	var drawCylinder = function(obj, from, to, radius, color) {
		if (!from || !to)
			return;

		var midpoint = new TV3().addVectors(from, to).multiplyScalar(0.5);
		var color = new TCo(color);

		var cylinderMaterial = new THREE.MeshLambertMaterial({
			color : color.getHex()
		});
		var cylinder = new THREE.Mesh(cylinderGeometry, cylinderMaterial);
		cylinder.position = midpoint;
		cylinder.lookAt(from);
		cylinder.updateMatrix();
		cylinder.matrixAutoUpdate = false;
		var m = new THREE.Matrix4().makeScale(radius, radius, from
				.distanceTo(to));
		m.rotateX(Math.PI / 2);
		cylinder.matrix.multiply(m);
		obj.add(cylinder);
	};
	
	var drawBondSticks = function(atom, atoms) {
		if (!atom.style.stick)
			return;
		var style = atom.style.stick;
		if (style.hidden)
			return;
		
		var geo = new THREE.Geometry();
		var bondR = style.radius || defaultStickRadius;
		
		for(var i = 0; i < atom.bonds.length; i++) {
			var j = atom.bonds[i]; //our neighbor
			if(i < j) {//only draw if less
				//TODO: handle bond orders
				var atom2 = atoms[j];
				if(!atom2.style.stick)
					continue; //don't sweat the details

				var p1 = new TV3(atom.x, atom.y, atom.z);
				var p2 = new TV3(atom2.x, atom2.y, atom2.z);
				var mp = new TV3().addVectors(p1,p2).multiplyScalar(0.5);

				var c1 = atom.color;
				var c2 = atom2.color;
				if(typeof(style.color) != "undefined") {
					c1 = style.color;
					c2 = style.color;
				}
				
				atom.globj = atom.globj || new THREE.Object3D();
				drawCylinder(atom.globj, p1, mp, bondR, c1);
				drawCylinder(atom.globj, p2, mp, bondR, c2);				
			}
		}
		
		//for junctions draw sphere
		var savedstyle = atom.style;
		atom.style = {sphere: {fixedRadius: bondR, color: style.color}};
		drawAtomSphere(atom);
		atom.style = savedstyle;
	};
	
	// go through all the atoms and regenerate their geometries
	// at some point we should optimize this to avoid unnecessary
	// recalculation
	var createMolObj = function(atoms) {
		var ret = new THREE.Object3D();

		for ( var i = 0; i < atoms.length; i++) {
			var atom = atoms[i];
			//recreate gl info for each atom as necessary
			if (atom && atom.style && atom.globj == null) {
				drawAtomSphere(atom);
				drawAtomIcosahedron(atom);
				drawAtomCross(atom);
				drawBondLines(atom, atoms);
				drawBondSticks(atom, atoms);
			}
			if(atom.globj)
				ret.add(atom.globj);
		}
		
		//TODO: create cartoons as needed
		return ret;
	};

	function GLModel(mid) {

		// private variables
		var atoms = [];
		var id = mid;
		var molObj = null;
		var redrawCartoon = true; // set when cartoon style changed

		this.getID = function() {
			return id;
		};

		// add atoms to this model from molecular data string
		this.addMolData = function(data, format) {

			switch (format) {
			case "xyz":
				parseXYZ(atoms, data);
				break;
			case "pdb":
				break;
			case "sdf":
				parseSDF(atoms, data);
				break;
			}
			setAtomDefaults(atoms, id);
		};

		// return list of atoms selected by sel
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
			//do a copy to enforce style changes through this function
			var mystyle = $.extend(true, {}, style);
			
			for(var i = 0; i < atoms.length; i++) {
				if(atoms[i].style.cartoon != mystyle.cartoon)
					redrawCartoon = true;
				atoms[i].style = mystyle;
				atoms[i].globj = null; //need to recalculate
			}
				
		};

		// return 3d data for this model
		this.globj = function() {
			molObj = createMolObj(atoms);
			return molObj;
		};


	}
	;

	return GLModel;
})();