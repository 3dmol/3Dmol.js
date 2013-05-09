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
	var thickness = 0.4;

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
				ret.push(new TV3(x, y, z));
			}
		}
		ret.push(points[points.length - 1]);
		return ret;
	};

	var drawThinStrip = function(group, p1, p2, colors, div) {
		var geo = new THREE.Geometry();
		for ( var i = 0, lim = p1.length; i < lim; i++) {
			geo.vertices.push(p1[i]); // 2i
			geo.vertices.push(p2[i]); // 2i + 1
		}
		for ( var i = 1, lim = p1.length; i < lim; i++) {
			var f = new THREE.Face4(2 * i, 2 * i + 1, 2 * i - 1, 2 * i - 2);
			f.color = WebMol.color(colors[Math.round((i - 1) / div)]);
			geo.faces.push(f);
		}
		geo.computeFaceNormals();
		geo.computeVertexNormals(false);
		var material = new THREE.MeshLambertMaterial();
		material.vertexColors = THREE.FaceColors;
		var mesh = new THREE.Mesh(geo, material);
		mesh.doubleSided = true;
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

		var geo = new THREE.Geometry();
		var vs = geo.vertices, fs = geo.faces;
		var axis, p1v, p2v, a1v, a2v;
		for ( var i = 0, lim = p1.length; i < lim; i++) {
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
		}
		var faces = [ [ 0, 2, -6, -8 ], [ -4, -2, 6, 4 ], [ 7, 3, -5, -1 ],
				[ -3, -7, 1, 5 ] ];
		for ( var i = 1, lim = p1.length; i < lim; i++) {
			var offset = 8 * i, color = WebMol.CC.color(colors[Math
					.round((i - 1) / div)]);
			for ( var j = 0; j < 4; j++) {
				var f = new THREE.Face4(offset + faces[j][0], offset
						+ faces[j][1], offset + faces[j][2], offset
						+ faces[j][3], undefined, color);
				fs.push(f);
			}
		}
		var vsize = vs.length - 8; // Cap
		for ( var i = 0; i < 4; i++) {
			vs.push(vs[i * 2]);
			vs.push(vs[vsize + i * 2])
		}
		;
		vsize += 8;
		fs.push(new THREE.Face4(vsize, vsize + 2, vsize + 6, vsize + 4,
				undefined, fs[0].color));
		fs.push(new THREE.Face4(vsize + 1, vsize + 5, vsize + 7, vsize + 3,
				undefined, fs[fs.length - 3].color));
		geo.computeFaceNormals();
		geo.computeVertexNormals(false);
		var material = new THREE.MeshLambertMaterial();
		material.vertexColors = THREE.FaceColors;
		material.side = THREE.DoubleSide;
		var mesh = new THREE.Mesh(geo, material);
		group.add(mesh);
	};

	var drawSmoothCurve = function(group, _points, width, colors, div) {
		if (_points.length == 0)
			return;

		div = (div == undefined) ? 5 : div;

		var geo = new THREE.Geometry();
		var points = subdivide(_points, div);

		for ( var i = 0; i < points.length; i++) {
			geo.vertices.push(points[i]);
			geo.colors.push(WebMol.color(colors[(i == 0) ? 0 : Math.round((i - 1)
					/ div)]));
		}
		var lineMaterial = new THREE.LineBasicMaterial({
			linewidth : width
		});
		lineMaterial.vertexColors = true;
		var line = new THREE.Line(geo, lineMaterial);
		line.type = THREE.LineStrip;
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
				if (atom.atom == 'CA') {
					currentCA = new TV3(atom.x, atom.y, atom.z);
					currentChain = atom.chain;
					currentResi = atom.resi;
					ss = atom.ss;
					ssborder = atom.ssbegin || atom.ssend;
					var atomcolor = atom.color;
					if(typeof(atom.style.cartoon.color) != "undefined") {
						atomcolor = atom.style.cartoon.color;
					}
					colors.push(atomcolor);
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
	}

	return drawCartoon;
})();
