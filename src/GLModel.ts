// A model is a collection of related atoms.  Bonds are only allowed between
//atoms in the same model.  An atom is uniquely specified by its model id and
//its serial number.
//A glmodel knows how to apply the styles on each atom to create a gl object
import { Geometry, Material, StickImposterMaterial } from "./WebGL";
import { Sphere, Cylinder } from "./WebGL/shapes";
import { Vector3, Matrix4, conversionMatrix3, Matrix3, XYZ } from "./WebGL/math";
import { Color, CC, ColorschemeSpec, ColorSpec } from "./colors";
import { InstancedMaterial, SphereImposterMaterial, MeshLambertMaterial, Object3D, Mesh, LineBasicMaterial, Line, LineStyle } from "./WebGL";
import { GLDraw } from "./GLDraw"
import { CartoonStyleSpec, drawCartoon } from "./glcartoon";
import { elementColors } from "./colors";
import { get, deepCopy, extend, getExtent, getAtomProperty, makeFunction, getPropertyRange, specStringToObject, getbin, getColorFromStyle } from "./utilities";
import { Gradient } from "./Gradient";
import { Parsers } from "./parsers";
import { NetCDFReader } from "netcdfjs"
import { inflate } from "pako"
import { AtomSelectionSpec, AtomSpec } from "./specs";
import { GLViewer } from "GLViewer";
import { ArrowSpec } from "GLShape";
import { ParserOptionsSpec } from "./parsers/ParserOptionsSpec";
import { LabelSpec } from "Label";
import { assignBonds } from "./parsers/utils/assignBonds";

/**
 * GLModel represents a group of related atoms
 * @class
 */
export class GLModel {

    // class variables go here
    static defaultAtomStyle: AtomStyleSpec = {
        line: {}
    };

    static defaultlineWidth = 1.0;

    // Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
    // https://en.wikipedia.org/wiki/Van_der_Waals_radius
    static vdwRadii = {
        "H": 1.2,
        "He": 1.4,
        "Li": 1.82,
        "Be": 1.53,
        "B": 1.92,
        "C": 1.7,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "Ne": 1.54,
        "Na": 2.27,
        "Mg": 1.73,
        "Al": 1.84,
        "Si": 2.1,
        "P": 1.8,
        "S": 1.8,
        "Cl": 1.75,
        "Ar": 1.88,
        "K": 2.75,
        "Ca": 2.31,
        "Ni": 1.63,
        "Cu": 1.4,
        "Zn": 1.39,
        "Ga": 1.87,
        "Ge": 2.11,
        "As": 1.85,
        "Se": 1.9,
        "Br": 1.85,
        "Kr": 2.02,
        "Rb": 3.03,
        "Sr": 2.49,
        "Pd": 1.63,
        "Ag": 1.72,
        "Cd": 1.58,
        "In": 1.93,
        "Sn": 2.17,
        "Sb": 2.06,
        "Te": 2.06,
        "I": 1.98,
        "Xe": 2.16,
        "Cs": 3.43,
        "Ba": 2.68,
        "Pt": 1.75,
        "Au": 1.66,
        "Hg": 1.55,
        "Tl": 1.96,
        "Pb": 2.02,
        "Bi": 2.07,
        "Po": 1.97,
        "At": 2.02,
        "Rn": 2.20,
        "Fr": 3.48,
        "Ra": 2.83,
        "U": 1.86
    };

    // class functions
    // return true if a and b represent the same style
    static sameObj(a, b) {
        if (a && b)
            return JSON.stringify(a) == JSON.stringify(b);
        else
            return a == b;
    };

    public unitCellObjects: any;

    // private variables
    private atoms: AtomSpec[] = [];
    private frames: any = [];
    private box: any = null;
    private atomdfs: any = null; //depth first search over connected components
    private id = 0;
    private hidden: any = false;
    private molObj: any = null;
    private renderedMolObj: any = null;
    private lastColors: any = null;
    private modelData: any = {};
    private modelDatas: any = null; //if there is different modelData per frame
    private idMatrix = new Matrix4();
    private dontDuplicateAtoms = true;
    private defaultColor = elementColors.defaultColor;

    private options: any;
    private ElementColors: any;


    private readonly defaultSphereRadius: number;
    private readonly defaultCartoonQuality: number;
    // bonds as cylinders
    private readonly defaultStickRadius = 0.25;

    constructor(mid, options?) {

        this.options = options || {};
        this.ElementColors = (this.options.defaultcolors) ? this.options.defaultcolors : elementColors.defaultColors;

        this.defaultSphereRadius = (this.options.defaultSphereRadius) ? this.options.defaultSphereRadius : 1.5;
        this.defaultCartoonQuality = (this.options.cartoonQuality) ? this.options.cartoonQuality : 10;
        this.id = mid;
    }
    // return proper radius for atom given style
    /**
     *
     * @param {AtomSpec} atom
     * @param {atomstyle} style
     * @return {number}
     *
     */
    private getRadiusFromStyle(atom:AtomSpec, style:SphereStyleSpec|ClickSphereStyleSpec|CrossStyleSpec) {
        var r = this.defaultSphereRadius;
        if (typeof (style.radius) != "undefined")
            r = style.radius;
        else if (GLModel.vdwRadii[atom.elem])
            r = GLModel.vdwRadii[atom.elem];
        else if (atom.elem.length > 1) { //see if adjusting case helps
            let e: string = atom.elem;
            e = e[0].toUpperCase() + e[1].toLowerCase();
            if (GLModel.vdwRadii[e])
                r = GLModel.vdwRadii[e];
        }

        if (typeof (style.scale) != "undefined")
            r *= style.scale;
        return r;
    };

    // cross drawing

    /**
     *
     * @param {AtomSpec} atom
     * @param {Record<number, Geometry>} geos
     */
    private drawAtomCross(atom:AtomSpec, geos: Record<number, Geometry>) {
        if (!atom.style.cross)
            return;
        var style = atom.style.cross;
        if (style.hidden)
            return;
        var linewidth = (style.linewidth || GLModel.defaultlineWidth);
        if (!geos[linewidth])
            geos[linewidth] = new Geometry();

        var geoGroup = geos[linewidth].updateGeoGroup(6);

        var delta = this.getRadiusFromStyle(atom, style);

        var points = [[delta, 0, 0], [-delta, 0, 0], [0, delta, 0],
        [0, -delta, 0], [0, 0, delta], [0, 0, -delta]];

        var clickable = atom.clickable || atom.hoverable;
        if (clickable && atom.intersectionShape === undefined)
            atom.intersectionShape = { sphere: [], cylinder: [], line: [] };

        var c = getColorFromStyle(atom, style);

        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;

        for (var j = 0; j < 6; j++) {

            var offset = geoGroup.vertices * 3;

            geoGroup.vertices++;
            vertexArray[offset] = atom.x + points[j][0];
            vertexArray[offset + 1] = atom.y + points[j][1];
            vertexArray[offset + 2] = atom.z + points[j][2];
            colorArray[offset] = c.r;
            colorArray[offset + 1] = c.g;
            colorArray[offset + 2] = c.b;

            if (clickable) {
                var point = new Vector3(points[j][0], points[j][1], points[j][2]);

                //decrease cross size for selection to prevent misselection from atom overlap
                point.multiplyScalar(0.1);
                point.set(point.x + atom.x, point.y + atom.y, point.z + atom.z);
                atom.intersectionShape.line.push(point);
            }

        }

    };

    private getGoodCross(atom:AtomSpec, atom2:AtomSpec, p1, dir) {
        // get vector 2 different neighboring atom
        //find most divergent neighbor
        var bestv = null;
        var bestlen = -1;
        for (var j = 0, n = atom.bonds.length; j < n; j++) {
            if (atom.bonds[j] != atom2.index) {
                let j2 = atom.bonds[j];
                let atom3 = this.atoms[j2];
                let p3 = new Vector3(atom3.x, atom3.y, atom3.z);

                let dir2 = p3.clone();
                dir2.sub(p1);

                let v = dir2.clone();
                v.cross(dir);
                var l = v.lengthSq();
                if (l > bestlen) {
                    bestlen = l;
                    bestv = v;
                    if (bestlen > 0.1) {
                        return bestv;
                    }
                }
            }
        }
        return bestv;
    };

    //from atom, return a normalized vector v that is orthogonal and along which
    //it is appropraite to draw multiple bonds
    private getSideBondV(atom:AtomSpec, atom2:AtomSpec, i:number) {

        var i2, j2, atom3, p3, dir2;
        var p1 = new Vector3(atom.x, atom.y, atom.z);
        var p2 = new Vector3(atom2.x, atom2.y, atom2.z);
        var dir = p2.clone();
        var v = null;
        dir.sub(p1);

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
                atom3 = this.atoms[j2];
                p3 = new Vector3(atom3.x, atom3.y, atom3.z);

                dir2 = p3.clone();
                dir2.sub(p1);

                v = dir2.clone();
                v.cross(dir);
            }
        } else {
            v = this.getGoodCross(atom, atom2, p1, dir);

            if (v.lengthSq() < 0.01) {
                var v2 = this.getGoodCross(atom2, atom, p1, dir);
                if (v2 != null) v = v2; //can be null if no other neighbors
            }
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

    };

    private addLine(vertexArray, colorArray, offset, p1: Vector3, p2: Vector3, c1: Color) {
        //make line from p1 to p2, does not incremeant counts
        vertexArray[offset] = p1.x; vertexArray[offset + 1] = p1.y; vertexArray[offset + 2] = p1.z;
        colorArray[offset] = c1.r; colorArray[offset + 1] = c1.g; colorArray[offset + 2] = c1.b;
        vertexArray[offset + 3] = p2.x; vertexArray[offset + 4] = p2.y; vertexArray[offset + 5] = p2.z;
        colorArray[offset + 3] = c1.r; colorArray[offset + 4] = c1.g; colorArray[offset + 5] = c1.b;
    };


    // bonds - both atoms must match bond style
    // standardize on only drawing for lowest to highest
    /**
     *
     * @param {AtomSpec}
     *            atom
     * @param {AtomSpec[]} atoms
     * @param {Record<number,Geometry>} geos
     */
    private drawBondLines(atom:AtomSpec, atoms:AtomSpec[], geos: Record<number,Geometry>) {
        if (!atom.style.line)
            return;
        var style = atom.style.line;
        if (style.hidden)
            return;
        var p1a, p1b, p2a, p2b;
        // have a separate geometry for each linewidth
        var linewidth = (style.linewidth || GLModel.defaultlineWidth);

        if (!geos[linewidth])
            geos[linewidth] = new Geometry();
        /** @type {geometryGroup} */
        var geoGroup = geos[linewidth].updateGeoGroup(6 * atom.bonds.length); //reserve enough space even for triple bonds

        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;

        for (var i = 0; i < atom.bonds.length; i++) {
            var j = atom.bonds[i]; // our neighbor

            var atom2 = atoms[j];
            if (!atom2.style.line)
                continue; // don't sweat the details

            if (atom.index >= atom2.index) // only draw if less, this way we can do multi bonds correctly
                continue;
            var p1 = new Vector3(atom.x, atom.y, atom.z);
            var p2 = new Vector3(atom2.x, atom2.y, atom2.z);
            var mp = p1.clone().add(p2).multiplyScalar(0.5);
            var singleBond = false;

            var atomneedsi = atom.clickable || atom.hoverable;
            var atom2needsi = atom2.clickable || atom2.hoverable;

            if (atomneedsi || atom2needsi) {
                if (atomneedsi) {
                    if (atom.intersectionShape === undefined)
                        atom.intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };
                    atom.intersectionShape.line.push(p1);
                    atom.intersectionShape.line.push(mp);
                }
                if (atom2needsi) {
                    if (atom2.intersectionShape === undefined)
                        atom2.intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };
                    atom2.intersectionShape.line.push(mp);
                    atom2.intersectionShape.line.push(p2);
                }
            }
            var c1 = getColorFromStyle(atom, atom.style.line);
            var c2 = getColorFromStyle(atom2, atom2.style.line);

            if (atom.bondStyles && atom.bondStyles[i]) {
                var bstyle = atom.bondStyles[i];
                if (!bstyle.iswire) {
                    continue;
                }
                if (bstyle.singleBond) singleBond = true;
                if (typeof (bstyle.color1) != "undefined") {
                    c1 = CC.color(bstyle.color1) as Color;
                }
                if (typeof (bstyle.color2) != "undefined") {
                    c2 = CC.color(bstyle.color2) as Color;
                }
            }

            var offset = geoGroup.vertices * 3;
            var mpa, mpb;

            if (atom.bondOrder[i] > 1 && atom.bondOrder[i] < 4 && !singleBond) {
                var v = this.getSideBondV(atom, atom2, i);
                var dir = p2.clone();
                dir.sub(p1);

                if (atom.bondOrder[i] == 2) { //double

                    v.multiplyScalar(0.1);
                    p1a = p1.clone();
                    p1a.add(v);
                    p1b = p1.clone();
                    p1b.sub(v);

                    p2a = p1a.clone();
                    p2a.add(dir);
                    p2b = p1b.clone();
                    p2b.add(dir);

                    if (c1 == c2) {
                        geoGroup.vertices += 4;
                        this.addLine(vertexArray, colorArray, offset, p1a, p2a, c1);
                        this.addLine(vertexArray, colorArray, offset + 6, p1b, p2b, c1);
                    }
                    else {
                        geoGroup.vertices += 8;
                        dir.multiplyScalar(0.5);
                        mpa = p1a.clone();
                        mpa.add(dir);
                        mpb = p1b.clone();
                        mpb.add(dir);

                        this.addLine(vertexArray, colorArray, offset, p1a, mpa, c1);
                        this.addLine(vertexArray, colorArray, offset + 6, mpa, p2a, c2);
                        this.addLine(vertexArray, colorArray, offset + 12, p1b, mpb, c1);
                        this.addLine(vertexArray, colorArray, offset + 18, mpb, p2b, c2);
                    }
                }
                else if (atom.bondOrder[i] == 3) { //triple

                    v.multiplyScalar(0.1);
                    p1a = p1.clone();
                    p1a.add(v);
                    p1b = p1.clone();
                    p1b.sub(v);

                    p2a = p1a.clone();
                    p2a.add(dir);
                    p2b = p1b.clone();
                    p2b.add(dir);

                    if (c1 == c2) {
                        geoGroup.vertices += 6;
                        this.addLine(vertexArray, colorArray, offset, p1, p2, c1);
                        this.addLine(vertexArray, colorArray, offset + 6, p1a, p2a, c1);
                        this.addLine(vertexArray, colorArray, offset + 12, p1b, p2b, c1);
                    }
                    else {
                        geoGroup.vertices += 12;
                        dir.multiplyScalar(0.5);
                        mpa = p1a.clone();
                        mpa.add(dir);
                        mpb = p1b.clone();
                        mpb.add(dir);

                        this.addLine(vertexArray, colorArray, offset, p1, mp, c1);
                        this.addLine(vertexArray, colorArray, offset + 6, mp, p2, c2);
                        this.addLine(vertexArray, colorArray, offset + 12, p1a, mpa, c1);
                        this.addLine(vertexArray, colorArray, offset + 18, mpa, p2a, c2);
                        this.addLine(vertexArray, colorArray, offset + 24, p1b, mpb, c1);
                        this.addLine(vertexArray, colorArray, offset + 30, mpb, p2b, c2);
                    }
                }
            }
            else { //single bond
                if (c1 == c2) {
                    geoGroup.vertices += 2;
                    this.addLine(vertexArray, colorArray, offset, p1, p2, c1);
                } else {
                    geoGroup.vertices += 4;
                    this.addLine(vertexArray, colorArray, offset, p1, mp, c1);
                    this.addLine(vertexArray, colorArray, offset + 6, mp, p2, c2);
                }

            }
        }

    };

    //sphere drawing
    //See also: drawCylinder
    /**
     *
     * @param {AtomSpec} atom
     * @param {Geometry} geo
     */
    private drawAtomSphere(atom: AtomSpec, geo: Geometry) {

        if (!atom.style.sphere)
            return;
        var style = atom.style.sphere;
        if (style.hidden)
            return;

        var C = getColorFromStyle(atom, style);

        var radius = this.getRadiusFromStyle(atom, style);

        if ((atom.clickable === true || atom.hoverable) && (atom.intersectionShape !== undefined)) {
            var center = new Vector3(atom.x, atom.y, atom.z);
            atom.intersectionShape.sphere.push(new Sphere(center, radius));
        }

        GLDraw.drawSphere(geo, atom, radius, C);
    };

    /** Register atom shaped click handlers */
    private drawAtomClickSphere(atom: AtomSpec) {

        if (!atom.style.clicksphere)
            return;
        var style = atom.style.clicksphere;
        if (style.hidden)
            return;

        var radius = this.getRadiusFromStyle(atom, style);

        if ((atom.clickable === true || atom.hoverable) && (atom.intersectionShape !== undefined)) {
            var center = new Vector3(atom.x, atom.y, atom.z);
            atom.intersectionShape.sphere.push(new Sphere(center, radius));
        }
    };

    private drawAtomInstanced(atom: AtomSpec, geo: Geometry) {

        if (!atom.style.sphere)
            return;
        var style = atom.style.sphere;
        if (style.hidden)
            return;

        var radius = this.getRadiusFromStyle(atom, style);
        var C = getColorFromStyle(atom, style);

        var geoGroup = geo.updateGeoGroup(1);
        var startv = geoGroup.vertices;
        var start = startv * 3;
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;
        var radiusArray = geoGroup.radiusArray;

        vertexArray[start] = atom.x;
        vertexArray[start + 1] = atom.y;
        vertexArray[start + 2] = atom.z;

        colorArray[start] = C.r;
        colorArray[start + 1] = C.g;
        colorArray[start + 2] = C.b;

        radiusArray[startv] = radius;

        if ((atom.clickable === true || atom.hoverable) && (atom.intersectionShape !== undefined)) {
            var center = new Vector3(atom.x, atom.y, atom.z);
            atom.intersectionShape.sphere.push(new Sphere(center, radius));
        }

        geoGroup.vertices += 1;

    };

    private drawSphereImposter(geo: Geometry, center: XYZ, radius: number, C: Color) {
        //create flat square
        var geoGroup = geo.updateGeoGroup(4);
        var i;
        var startv = geoGroup.vertices;
        var start = startv * 3;
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;

        //use center point for each vertex
        for (i = 0; i < 4; i++) {
            vertexArray[start + 3 * i] = center.x;
            vertexArray[start + 3 * i + 1] = center.y;
            vertexArray[start + 3 * i + 2] = center.z;
        }


        //same colors for all 4 vertices
        var normalArray = geoGroup.normalArray;
        for (i = 0; i < 4; i++) {
            colorArray[start + 3 * i] = C.r;
            colorArray[start + 3 * i + 1] = C.g;
            colorArray[start + 3 * i + 2] = C.b;
        }

        normalArray[start + 0] = -radius;
        normalArray[start + 1] = radius;
        normalArray[start + 2] = 0;

        normalArray[start + 3] = -radius;
        normalArray[start + 4] = -radius;
        normalArray[start + 5] = 0;

        normalArray[start + 6] = radius;
        normalArray[start + 7] = -radius;
        normalArray[start + 8] = 0;

        normalArray[start + 9] = radius;
        normalArray[start + 10] = radius;
        normalArray[start + 11] = 0;

        geoGroup.vertices += 4;

        //two faces
        var faceArray = geoGroup.faceArray;
        var faceoffset = geoGroup.faceidx; //not number faces, but index
        faceArray[faceoffset + 0] = startv;
        faceArray[faceoffset + 1] = startv + 1;
        faceArray[faceoffset + 2] = startv + 2;
        faceArray[faceoffset + 3] = startv + 2;
        faceArray[faceoffset + 4] = startv + 3;
        faceArray[faceoffset + 5] = startv;
        geoGroup.faceidx += 6;
    };

    //dkoes -  code for sphere imposters
    private drawAtomImposter(atom: AtomSpec, geo: Geometry) {

        if (!atom.style.sphere)
            return;
        var style = atom.style.sphere;
        if (style.hidden)
            return;

        var radius = this.getRadiusFromStyle(atom, style);
        var C = getColorFromStyle(atom, style);

        if ((atom.clickable === true || atom.hoverable) && (atom.intersectionShape !== undefined)) {
            var center = new Vector3(atom.x, atom.y, atom.z);
            atom.intersectionShape.sphere.push(new Sphere(center, radius));
        }

        this.drawSphereImposter(geo, atom as XYZ, radius, C);
    };


    static drawStickImposter(geo: Geometry, from: XYZ, to: XYZ, radius: number, color: Color) {
        //we need the four corners - two have from coord, two have to coord, the normal
        //is the opposing point, from which we can get the normal and length
        //also need the radius
        var geoGroup = geo.updateGeoGroup(4);
        var startv = geoGroup.vertices;
        var start = startv * 3;
        var vertexArray = geoGroup.vertexArray;
        var colorArray = geoGroup.colorArray;
        var radiusArray = geoGroup.radiusArray;
        var normalArray = geoGroup.normalArray;
        //encode extra bits of information in the color
        var r = color.r;
        var g = color.g;
        var b = color.b;

        var negateColor = function (c) {
            //set sign bit
            var n = -c;
            if (n == 0) n = -0.0001;
            return n;
        };

        /* for sticks, always draw caps, but we could in theory set caps in color */

        //4 vertices, distinguish between p1 and p2 with neg blue
        var pos = start;
        for (var i = 0; i < 4; i++) {
            vertexArray[pos] = from.x;
            normalArray[pos] = to.x;
            colorArray[pos] = r;
            pos++;
            vertexArray[pos] = from.y;
            normalArray[pos] = to.y;
            colorArray[pos] = g;
            pos++;
            vertexArray[pos] = from.z;
            normalArray[pos] = to.z;
            if (i < 2)
                colorArray[pos] = b;
            else
                colorArray[pos] = negateColor(b);
            pos++;
        }

        geoGroup.vertices += 4;

        radiusArray[startv] = -radius;
        radiusArray[startv + 1] = radius;
        radiusArray[startv + 2] = -radius;
        radiusArray[startv + 3] = radius;

        //two faces
        var faceArray = geoGroup.faceArray;
        var faceoffset = geoGroup.faceidx; //not number faces, but index
        faceArray[faceoffset + 0] = startv;
        faceArray[faceoffset + 1] = startv + 1;
        faceArray[faceoffset + 2] = startv + 2;
        faceArray[faceoffset + 3] = startv + 2;
        faceArray[faceoffset + 4] = startv + 3;
        faceArray[faceoffset + 5] = startv;
        geoGroup.faceidx += 6;
    };

    // draws cylinders and small spheres (at bond radius)
    private drawBondSticks(atom: AtomSpec, atoms: AtomSpec[], geo: Geometry) {
        if (!atom.style.stick)
            return;
        var style = atom.style.stick;
        if (style.hidden)
            return;

        var atomBondR = style.radius || this.defaultStickRadius;
        var bondR = atomBondR;
        var atomSingleBond = style.singleBonds || false;
        var fromCap = 0, toCap = 0;
        var atomneedsi, atom2needsi, i, singleBond, bstyle;
        var cylinder1a, cylinder1b, cylinder1c, cylinder2a, cylinder2b, cylinder2c;

        var C1 = getColorFromStyle(atom, style);

        var mp, mp2, mp3;

        if (!atom.capDrawn && atom.bonds.length < 4)
            fromCap = 2;

        var drawCyl = GLDraw.drawCylinder; //mesh cylinder
        if (geo.imposter)
            drawCyl = GLModel.drawStickImposter;


        for (i = 0; i < atom.bonds.length; i++) {
            var j = atom.bonds[i]; // our neighbor
            var atom2 = atoms[j]; //parsePDB, etc should only add defined bonds
            mp = mp2 = mp3 = null;
            if (atom.index < atom2.index) {// only draw if less, this
                // lets us combine
                // cylinders of the same
                // color
                var style2 = atom2.style;
                if (!style2.stick || style2.stick.hidden)
                    continue; // don't sweat the details

                var C2 = getColorFromStyle(atom2, style2.stick);

                //support bond specific styles
                bondR = atomBondR;
                singleBond = atomSingleBond;
                if (atom.bondStyles && atom.bondStyles[i]) {
                    bstyle = atom.bondStyles[i];
                    if (bstyle.iswire) {
                        continue;
                    }
                    if (bstyle.radius) bondR = bstyle.radius;
                    if (bstyle.singleBond) singleBond = true;
                    if (typeof (bstyle.color1) != "undefined") {
                        C1 = CC.color(bstyle.color1) as Color;
                    }
                    if (typeof (bstyle.color2) != "undefined") {
                        C2 = CC.color(bstyle.color2) as Color;
                    }
                }
                var p1 = new Vector3(atom.x, atom.y, atom.z);
                var p2 = new Vector3(atom2.x, atom2.y, atom2.z);

                // draw cylinders
                if (atom.bondOrder[i] <= 1 || singleBond || atom.bondOrder[i] > 3) { //TODO: aromatics at 4

                    if(atom.bondOrder[i] < 1) bondR *= atom.bondOrder[i];
                    if (!atom2.capDrawn && atom2.bonds.length < 4)
                        toCap = 2;

                    if (C1 != C2) {
                        mp = new Vector3().addVectors(p1, p2)
                            .multiplyScalar(0.5);
                        drawCyl(geo, p1, mp, bondR, C1, fromCap, 0);
                        drawCyl(geo, mp, p2, bondR, C2, 0, toCap);
                    } else {
                        drawCyl(geo, p1, p2, bondR, C1, fromCap, toCap);
                    }


                    atomneedsi = atom.clickable || atom.hoverable;
                    atom2needsi = atom2.clickable || atom2.hoverable;

                    if (atomneedsi || atom2needsi) {
                        if (!mp) mp = new Vector3().addVectors(p1, p2).multiplyScalar(0.5);
                        if (atomneedsi) {
                            var cylinder1 = new Cylinder(p1, mp, bondR);
                            var sphere1 = new Sphere(p1, bondR);
                            atom.intersectionShape.cylinder.push(cylinder1);
                            atom.intersectionShape.sphere.push(sphere1);
                        }
                        if (atom2needsi) {
                            var cylinder2 = new Cylinder(p2, mp, bondR);
                            var sphere2 = new Sphere(p2, bondR);
                            atom2.intersectionShape.cylinder.push(cylinder2);
                            atom2.intersectionShape.sphere.push(sphere2);
                        }
                    }
                }
                else if (atom.bondOrder[i] > 1) {
                    //multi bond caps
                    var mfromCap = 0;
                    var mtoCap = 0;

                    if (bondR != atomBondR) {
                        //assume jmol style multiple bonds - the radius doesn't fit within atom sphere
                        mfromCap = 2;
                        mtoCap = 2;
                    }

                    var dir = p2.clone();
                    var v = null;
                    dir.sub(p1);

                    var r, p1a, p1b, p2a, p2b;
                    v = this.getSideBondV(atom, atom2, i);

                    if (atom.bondOrder[i] == 2) {
                        r = bondR / 2.5;

                        v.multiplyScalar(r * 1.5);
                        p1a = p1.clone();
                        p1a.add(v);
                        p1b = p1.clone();
                        p1b.sub(v);

                        p2a = p1a.clone();
                        p2a.add(dir);
                        p2b = p1b.clone();
                        p2b.add(dir);

                        if (C1 != C2) {
                            mp = new Vector3().addVectors(p1a, p2a)
                                .multiplyScalar(0.5);
                            mp2 = new Vector3().addVectors(p1b, p2b)
                                .multiplyScalar(0.5);
                            drawCyl(geo, p1a, mp, r, C1, mfromCap, 0);
                            drawCyl(geo, mp, p2a, r, C2, 0, mtoCap);
                            drawCyl(geo, p1b, mp2, r, C1, mfromCap, 0);
                            drawCyl(geo, mp2, p2b, r, C2, 0, mtoCap);
                        } else {
                            drawCyl(geo, p1a, p2a, r, C1, mfromCap, mtoCap);
                            drawCyl(geo, p1b, p2b, r, C1, mfromCap, mtoCap);
                        }

                        atomneedsi = atom.clickable || atom.hoverable;
                        atom2needsi = atom2.clickable || atom2.hoverable;

                        if (atomneedsi || atom2needsi) {
                            if (!mp) mp = new Vector3().addVectors(p1a, p2a)
                                .multiplyScalar(0.5);
                            if (!mp2) mp2 = new Vector3().addVectors(p1b, p2b)
                                .multiplyScalar(0.5);
                            if (atomneedsi) {
                                cylinder1a = new Cylinder(p1a, mp, r);
                                cylinder1b = new Cylinder(p1b, mp2, r);
                                atom.intersectionShape.cylinder.push(cylinder1a);
                                atom.intersectionShape.cylinder.push(cylinder1b);
                            }
                            if (atom2needsi) {
                                cylinder2a = new Cylinder(p2a, mp, r);
                                cylinder2b = new Cylinder(p2b, mp2, r);
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
                            mp = new Vector3().addVectors(p1a, p2a)
                                .multiplyScalar(0.5);
                            mp2 = new Vector3().addVectors(p1b, p2b)
                                .multiplyScalar(0.5);
                            mp3 = new Vector3().addVectors(p1, p2)
                                .multiplyScalar(0.5);
                            drawCyl(geo, p1a, mp, r, C1, mfromCap, 0);
                            drawCyl(geo, mp, p2a, r, C2, 0, mtoCap);
                            drawCyl(geo, p1, mp3, r, C1, fromCap, 0);
                            drawCyl(geo, mp3, p2, r, C2, 0, toCap);
                            drawCyl(geo, p1b, mp2, r, C1, mfromCap, 0);
                            drawCyl(geo, mp2, p2b, r, C2, 0, mtoCap);
                        } else {
                            drawCyl(geo, p1a, p2a, r, C1, mfromCap, mtoCap);
                            drawCyl(geo, p1, p2, r, C1, fromCap, toCap);
                            drawCyl(geo, p1b, p2b, r, C1, mfromCap, mtoCap);
                        }

                        atomneedsi = atom.clickable || atom.hoverable;
                        atom2needsi = atom2.clickable || atom2.hoverable;

                        if (atomneedsi || atom2needsi) {
                            if (!mp) mp = new Vector3().addVectors(p1a, p2a)
                                .multiplyScalar(0.5);
                            if (!mp2) mp2 = new Vector3().addVectors(p1b, p2b)
                                .multiplyScalar(0.5);
                            if (!mp3) mp3 = new Vector3().addVectors(p1, p2)
                                .multiplyScalar(0.5);

                            if (atomneedsi) {
                                cylinder1a = new Cylinder(p1a.clone(), mp.clone(), r);
                                cylinder1b = new Cylinder(p1b.clone(), mp2.clone(), r);
                                cylinder1c = new Cylinder(p1.clone(), mp3.clone(), r);
                                atom.intersectionShape.cylinder.push(cylinder1a);
                                atom.intersectionShape.cylinder.push(cylinder1b);
                                atom.intersectionShape.cylinder.push(cylinder1c);
                            }
                            if (atom2needsi) {
                                cylinder2a = new Cylinder(p2a.clone(), mp.clone(), r);
                                cylinder2b = new Cylinder(p2b.clone(), mp2.clone(), r);
                                cylinder2c = new Cylinder(p2.clone(), mp3.clone(), r);
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
        for (i = 0; i < atom.bonds.length; i++) {
            singleBond = atomSingleBond;
            if (atom.bondStyles && atom.bondStyles[i]) {
                bstyle = atom.bondStyles[i];
                if (bstyle.singleBond) singleBond = true;
                if (bstyle.radius && bstyle.radius != atomBondR) {
                    differentradii = true;
                }
            }
            if (singleBond || atom.bondOrder[i] == 1) {
                numsinglebonds++;
            }
        }

        if (differentradii) { //jmol style double/triple bonds - no sphere
            if (numsinglebonds > 0) drawSphere = true; //unless needed as a cap
        }
        else if (numsinglebonds == 0 && (atom.bonds.length > 0 || style.showNonBonded)) {
            drawSphere = true;
        }

        if (drawSphere) {
            bondR = atomBondR;
            //do not use bond style as this can be variable, particularly
            //with jmol export of double/triple bonds
            if (geo.imposter) {
                this.drawSphereImposter(geo.sphereGeometry, atom as XYZ, bondR, C1);
            }
            else {
                GLDraw.drawSphere(geo, atom, bondR, C1);
            }
        }

    };



    // go through all the atoms and regenerate their geometries
    // we try to have one geometry for each style since this is much much
    // faster
    // at some point we should optimize this to avoid unnecessary
    // recalculation
    /** param {AtomSpec[]} atoms */
    private createMolObj(atoms:AtomSpec[], options?) {

        options = options || {};

        var ret = new Object3D();
        var cartoonAtoms = [];
        var lineGeometries: Record<number, Geometry> = {};
        var crossGeometries:Record<number, Geometry> = {};

        var drawSphereFunc = this.drawAtomSphere;
        var sphereGeometry: Geometry = null;
        var stickGeometry: Geometry = null;
        if (options.supportsImposters) {
            drawSphereFunc = this.drawAtomImposter;
            sphereGeometry = new Geometry(true);
            sphereGeometry.imposter = true;
            stickGeometry = new Geometry(true, true);
            stickGeometry.imposter = true;
            stickGeometry.sphereGeometry = new Geometry(true); //for caps
            stickGeometry.sphereGeometry.imposter = true;
            stickGeometry.drawnCaps = {};
        }
        else if (options.supportsAIA) {
            drawSphereFunc = this.drawAtomInstanced;
            sphereGeometry = new Geometry(false, true, true);
            sphereGeometry.instanced = true;
            stickGeometry = new Geometry(true); //don't actually have instanced sticks
        } else {
            sphereGeometry = new Geometry(true);
            stickGeometry = new Geometry(true);
        }

        var i, j, n, testOpacities;
        var opacities: any = {};
        var range = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
        for (i = 0, n = atoms.length; i < n; i++) {
            var atom = atoms[i];
            // recreate gl info for each atom as necessary
            // set up appropriate intersection spheres for clickable atoms

            if (atom && atom.style) {

                if ((atom.clickable || atom.hoverable) && atom.intersectionShape === undefined)
                    atom.intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };

                testOpacities = { line: undefined, cross: undefined, stick: undefined, sphere: undefined };
                for (j in testOpacities) {
                    if (atom.style[j]) {
                        if (atom.style[j].opacity)
                            testOpacities[j] = parseFloat(atom.style[j].opacity);
                        else
                            testOpacities[j] = 1;

                    } else testOpacities[j] = undefined;

                    if (opacities[j]) {
                        if (testOpacities[j] != undefined && opacities[j] != testOpacities[j]) {
                            console.log("Warning: " + j + " opacity is ambiguous");
                            opacities[j] = 1;
                        }

                    } else opacities[j] = testOpacities[j];
                }

                drawSphereFunc.call(this, atom, sphereGeometry);
                this.drawAtomClickSphere(atom);
                this.drawAtomCross(atom, crossGeometries);
                this.drawBondLines(atom, atoms, lineGeometries);
                this.drawBondSticks(atom, atoms, stickGeometry);

                if (typeof (atom.style.cartoon) !== "undefined" && !atom.style.cartoon.hidden) {
                    //gradient color scheme range
                    if (atom.style.cartoon.color === "spectrum" && typeof (atom.resi) === "number" && !atom.hetflag) {
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
            drawCartoon(ret, cartoonAtoms, range, this.defaultCartoonQuality);
        }

        // add sphere geometry
        if (sphereGeometry && sphereGeometry.vertices > 0) {
            //Initialize buffers in geometry
            sphereGeometry.initTypedArrays();
            var sphereMaterial = null;
            var sphere = null;

            //create appropriate material
            if (sphereGeometry.imposter) {
                sphereMaterial = new SphereImposterMaterial({
                    ambient: 0x000000,
                    vertexColors: true,
                    reflectivity: 0
                });
            }
            else if (sphereGeometry.instanced) {
                sphere = new Geometry(true);
                GLDraw.drawSphere(sphere, { x: 0, y: 0, z: 0 }, 1, new Color(0.5, 0.5, 0.5));
                sphere.initTypedArrays();
                sphereMaterial = new InstancedMaterial({
                    sphereMaterial: new MeshLambertMaterial({
                        ambient: 0x000000,
                        vertexColors: true,
                        reflectivity: 0,
                    }),
                    sphere: sphere
                });
            }
            else { //regular mesh
                sphereMaterial = new MeshLambertMaterial({
                    ambient: 0x000000,
                    vertexColors: true,
                    reflectivity: 0,
                });
            }
            if (opacities.sphere < 1 && opacities.sphere >= 0) {
                sphereMaterial.transparent = true;
                sphereMaterial.opacity = opacities.sphere;
            }

            sphere = new Mesh(sphereGeometry, sphereMaterial);
            ret.add(sphere);
        }

        // add stick geometry
        if (stickGeometry.vertices > 0) {

            var stickMaterial = null;
            var ballMaterial = null;
            var balls = stickGeometry.sphereGeometry;
            if (!balls || typeof (balls.vertices) === 'undefined' || balls.vertices == 0) balls = null; //no balls

            //Initialize buffers in geometry
            stickGeometry.initTypedArrays();
            if (balls) balls.initTypedArrays();

            //create material
            var matvals = { ambient: 0x000000, vertexColors: true, reflectivity: 0 };

            if (stickGeometry.imposter) {
                stickMaterial = new StickImposterMaterial(matvals);
                ballMaterial = new SphereImposterMaterial(matvals);
            } else {
                stickMaterial = new MeshLambertMaterial(matvals);
                ballMaterial = new MeshLambertMaterial(matvals);

                if (stickMaterial.wireframe) {
                    stickGeometry.setUpWireframe();
                    if (balls) balls.setUpWireframe();
                }
            }

            if (opacities.stick < 1 && opacities.stick >= 0) {
                stickMaterial.transparent = true;
                stickMaterial.opacity = opacities.stick;
                ballMaterial.transparent = true;
                ballMaterial.opacity = opacities.stick;
            }
            var sticks = new Mesh(stickGeometry, stickMaterial);
            ret.add(sticks);

            if (balls) {
                var stickspheres = new Mesh(balls, ballMaterial);
                ret.add(stickspheres);
            }
        }

        //var linewidth;
        // add any line geometries, distinguished by line width
        var linewidth;
        for (i in lineGeometries) {
            if (lineGeometries.hasOwnProperty(i)) {
                linewidth = i;
                var lineMaterial = new LineBasicMaterial({
                    linewidth: linewidth,
                    vertexColors: true
                });
                if (opacities.line < 1 && opacities.line >= 0) {
                    lineMaterial.transparent = true;
                    lineMaterial.opacity = opacities.line;
                }

                lineGeometries[i].initTypedArrays();

                var line = new Line(lineGeometries[i], lineMaterial as Material, LineStyle.LinePieces);

                ret.add(line);
            }
        }

        // add any cross geometries
        for (i in crossGeometries) {
            if (crossGeometries.hasOwnProperty(i)) {
                linewidth = i;
                var crossMaterial = new LineBasicMaterial({
                    linewidth: linewidth,
                    vertexColors: true
                });
                if (opacities.cross < 1 && opacities.cross >= 0) {
                    crossMaterial.transparent = true;
                    crossMaterial.opacity = opacities.cross;
                }

                crossGeometries[i].initTypedArrays();

                var cross = new Line(crossGeometries[i], crossMaterial as Material, LineStyle.LinePieces);

                ret.add(cross);
            }
        }


        //for BIOMT assembly
        if (this.dontDuplicateAtoms && this.modelData.symmetries && this.modelData.symmetries.length > 0) {
            var finalRet = new Object3D();
            var t;
            for (t = 0; t < this.modelData.symmetries.length; t++) {
                var transformedRet = new Object3D();
                transformedRet = ret.clone();
                transformedRet.matrix.copy(this.modelData.symmetries[t]);
                transformedRet.matrixAutoUpdate = false;
                finalRet.add(transformedRet);
            }
            return finalRet;
        }

        return ret;
    };



    /**
     * Return object representing internal state of
     * the model appropriate for passing to setInternalState
     *
    */
    public getInternalState() {
        return {
            'atoms': this.atoms,
            'frames': this.frames
        };
    };

    /**
     * Overwrite the internal model state with the passed state.
     *
    */
    public setInternalState(state) {
        this.atoms = state.atoms;
        this.frames = state.frames;
        this.molObj = null;
    };

    /**
     * Returns crystallographic information if present.
     *
     *
     */
    public getCrystData() {
        if (this.modelData.cryst) {
            // add the matrix if it is missing
            if (!this.modelData.cryst.matrix) {
                const cryst = this.modelData.cryst;
                this.modelData.cryst.matrix = conversionMatrix3(
                    cryst.a, cryst.b, cryst.c,
                    cryst.alpha, cryst.beta, cryst.gamma
                );
            }
            return this.modelData.cryst;
        } else {
            return null;
        }
    };

    /**
     * Set crystallographic information using three angles and three lengths
     *
     * @param {number} a - length of unit cell side
     * @param {number} b - length of unit cell side
     * @param {number} c - length of unit cell side
     * @param {number} alpha - unit cell angle in degrees (default 90)
     * @param {number} beta - unit cell angle in degrees (default 90)
     * @param {number} gamma - unit cell angle in degrees (default 90)

     */
    public setCrystData(a?: number, b?:number, c?:number, alpha?:number, beta?:number, gamma?:number) {
        //I am assuming these
        a = a || 1.0;
        b = b || 1.0;
        c = c || 1.0;
        alpha = alpha || 90;
        beta = beta || 90;
        gamma = gamma || 90;

        const matrix = conversionMatrix3(a, b, c, alpha, beta, gamma);
        this.modelData.cryst = {
            'a': a, 'b': b, 'c': c,
            'alpha': alpha, 'beta': beta, 'gamma': gamma,
            'matrix': matrix
        };
    };

    /**
     * Set the crystallographic matrix to the given matrix.
     *
     * This function removes `a`, `b`, `c`, `alpha`, `beta`, `gamma` from
     * the crystal data.
     *
     * @param {Matrix3} matrix - unit cell matrix
     */
    public setCrystMatrix(matrix: Matrix3) {
        matrix = matrix || new Matrix3(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        );

        this.modelData.cryst = {
            'matrix': matrix
        };
    };

    /**
     * Returns list of rotational/translational matrices if there is BIOMT data
     * Otherwise returns a list of just the ID matrix
     *
     * @return {Array<Matrix4>}
     *
     */
    public getSymmetries() {

        if (typeof (this.modelData.symmetries) == 'undefined') {
            this.modelData.symmetries = [this.idMatrix];
        }
        return this.modelData.symmetries;
    };

    /**
     * Sets symmetries based on specified matrices in list
     *
     * @param {Array<Matrix4>} list
     *
     */
    public setSymmetries(list) {
        if (typeof (list) == "undefined") { //delete sym data
            this.modelData.symmetries = [this.idMatrix];
        }
        else {
            this.modelData.symmetries = list;
        }
    };

    /**
     * Returns model id number
     *
     * @return {number} Model ID
     */
    public getID() {
        return this.id;
    };

    /**
     * Returns model's frames property, a list of atom lists
     *
     * @return {number}
     */
    public getNumFrames() {
        return (this.frames.numFrames != undefined) ? this.frames.numFrames : this.frames.length;
    };

    private adjustCoord(x1:number, x2:number, margin:number, adjust:number) {
        //return new value of x2 that isn't more than margin away
        var dist = x2 - x1;
        if (dist < -margin) {
            return x2 + adjust;
        } else if (dist > margin) {
            return x2 - adjust;
        }
        return x2;
    };
    //go over current atoms in depth first order and ensure that connected
    //attoms aren't split across the box
    private adjustCoordinatesToBox() {
        if (!this.box) return;
        if (!this.atomdfs) return;
        var bx = this.box[0];
        var by = this.box[1];
        var bz = this.box[2];
        var mx = bx * 0.9;
        var my = by * 0.9;
        var mz = bz * 0.9;

        for (var c = 0; c < this.atomdfs.length; c++) {
            //for each connected component
            var component = this.atomdfs[c];
            for (var i = 1; i < component.length; i++) {
                //compare each atom to its previous and prevent coordinates from wrapping
                var atom = this.atoms[component[i][0]];
                var prev = this.atoms[component[i][1]];
                atom.x = this.adjustCoord(prev.x, atom.x, mx, bx);
                atom.y = this.adjustCoord(prev.y, atom.y, my, by);
                atom.z = this.adjustCoord(prev.z, atom.z, mz, bz);
            }
        }
    };

    /**
     * Sets model's atomlist to specified frame
     * Sets to last frame if framenum out of range
     *
     * @param {number} framenum - model's atoms are set to this index in frames list
     * @return {Promise}
     */
    public setFrame(framenum: number, viewer?: GLViewer) { //viewer only passed internally for unit cell
        var numFrames = this.getNumFrames();
        let model = this;
        return new Promise<void>(function (resolve, reject) {
            if (numFrames == 0) {
                //return;
                resolve();
            }
            if (framenum < 0 || framenum >= numFrames) {
                framenum = numFrames - 1;
            }
            if (model.frames.url != undefined) {
                var url = model.frames.url;
                getbin(url + "/traj/frame/" + framenum + "/" + model.frames.path, undefined, 'POST', undefined).then(function (buffer) {
                    var values = new Float32Array(buffer, 44);
                    var count = 0;
                    for (var i = 0; i < model.atoms.length; i++) {
                        model.atoms[i].x = values[count++];
                        model.atoms[i].y = values[count++];
                        model.atoms[i].z = values[count++];
                    }
                    //if a box was provided, check to see if we need to wrap connected components
                    if (model.box && model.atomdfs) {
                        model.adjustCoordinatesToBox();
                    }
                    resolve();
                }).catch(reject);
            }
            else {
                model.atoms = model.frames[framenum];
                resolve();
            }
            model.molObj = null;
            if (model.modelDatas && framenum < model.modelDatas.length) {
                model.modelData = model.modelDatas[framenum];
                if (model.unitCellObjects && viewer) {
                    viewer.removeUnitCell(model);
                    viewer.addUnitCell(model);
                }
            }
        });
    };

    /**
     * Add atoms as frames of model
     *
     * @param {AtomSpec[]} atoms - atoms to be added
     */
    public addFrame(atoms: AtomSpec[]) {
        this.frames.push(atoms);
    };


    /**
     * If model atoms have dx, dy, dz properties (in some xyz files), vibrate populates the model's frame property based on parameters.
     * Model can then be animated
     *
     * @param {number} numFrames - number of frames to be created, default to 10
     * @param {number} amplitude - amplitude of distortion, default to 1 (full)
     * @param {boolean} bothWays - if true, extend both in positive and negative directions by numFrames
     * @param {GLViewer} viewer - required if arrowSpec is provided
     * @param {ArrowSpec} arrowSpec - specification for drawing animated arrows. If color isn't specified, atom color (sphere, stick, line preference) is used.
     *@example

      $3Dmol.download("pdb:4UAA",viewer,{},function(){
        viewer.setStyle({},{stick:{}});
        viewer.vibrate(10, 1);
        viewer.animate({loop: "forward",reps: 1});

        viewer.zoomTo();
              viewer.render();
          });
     */
    public vibrate(numFrames:number=10, amplitude:number=1, bothWays:boolean=false, viewer?:GLViewer, arrowSpec?:ArrowSpec) {
        var start = 0;
        var end = numFrames;
        if (bothWays) {
            start = -numFrames;
            end = numFrames;
        }

        //to enable multiple setting of vibrate with bothWays, must record original position
        if (this.frames !== undefined && this.frames.origIndex !== undefined) {
            this.setFrame(this.frames.origIndex);
        } else {
            this.setFrame(0);
        }

        if (start < end) this.frames = []; //clear
        if (bothWays) this.frames.origIndex = numFrames;

        for (var i = start; i < end; i++) {
            var newAtoms = [];
            var currframe = this.frames.length;
            if (i == 0 && !arrowSpec) { //still need to calculate if drawing arrows
                this.frames.push(this.atoms);
                continue;
            }
            for (var j = 0; j < this.atoms.length; j++) {
                var dx = getAtomProperty(this.atoms[j], 'dx');
                var dy = getAtomProperty(this.atoms[j], 'dy');
                var dz = getAtomProperty(this.atoms[j], 'dz');
                var newVector = new Vector3(dx, dy, dz);
                var starting = new Vector3(this.atoms[j].x, this.atoms[j].y, this.atoms[j].z);
                var mult = (i * amplitude) / numFrames;
                newVector.multiplyScalar(mult);
                starting.add(newVector);
                var newAtom: any = {};
                for (var k in this.atoms[j]) {
                    newAtom[k] = this.atoms[j][k];
                }
                newAtom.x = starting.x;
                newAtom.y = starting.y;
                newAtom.z = starting.z;
                newAtoms.push(newAtom);
                if (viewer && arrowSpec) {
                    var spec = extend({}, arrowSpec);
                    var arrowend = new Vector3(dx, dy, dz);
                    arrowend.multiplyScalar(amplitude);
                    arrowend.add(starting);

                    spec.start = starting;
                    spec.end = arrowend;
                    spec.frame = currframe;
                    if (!spec.color) {
                        var s = newAtom.style.sphere;
                        if (!s) s = newAtom.style.stick;
                        if (!s) s = newAtom.style.line;
                        spec.color = getColorFromStyle(newAtom, s);
                    }
                    viewer.addArrow(spec);
                }
            }
            this.frames.push(newAtoms);
        }
    };

    // set default style and colors for atoms
    public setAtomDefaults(atoms: AtomSpec[]) {
        for (let i = 0; i < atoms.length; i++) {
            let atom = atoms[i];
            if (atom) {
                atom.style = atom.style || deepCopy(GLModel.defaultAtomStyle);
                atom.color = atom.color || this.ElementColors[atom.elem] || this.defaultColor;
                atom.model = this.id;
                if (atom.clickable || atom.hoverable)
                    atom.intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };
            }
        }
    };

    /** add atoms to this model from molecular data string
     *
     * @param {string|ArrayBuffer} data - atom structure file input data string, for gzipped input use ArrayBuffer
     * @param {string} format - input file string format (e.g 'pdb', 'sdf', 'sdf.gz', etc.)
     * @param {ParserOptionsSpec} options - format dependent options. Attributes depend on the input format
     */
    public addMolData(data: string|ArrayBuffer, format:string, options: ParserOptionsSpec={}) {
        var parsedAtoms = GLModel.parseMolData(data, format, options);
        this.dontDuplicateAtoms = !options.duplicateAssemblyAtoms;
        var mData = parsedAtoms.modelData;
        if (mData) {
            if (Array.isArray(mData)) {
                this.modelData = mData[0];
                if (options.frames) {
                    this.modelDatas = mData;
                }
            } else {
                this.modelData = mData;
            }
        }

        if (parsedAtoms.box) {
            this.box = parsedAtoms.box;
        } else {
            this.box = null;
        }

        if (this.frames.length == 0) { //first call
            for (let i = 0; i < parsedAtoms.length; i++) {
                if (parsedAtoms[i].length != 0)
                    this.frames.push(parsedAtoms[i]);
            }
            if (this.frames[0])
                this.atoms = this.frames[0];
        }

        else { //subsequent calls
            if (options.frames) { //add to new frame
                for (let i = 0; i < parsedAtoms.length; i++) {
                    this.frames.push(parsedAtoms[i]);
                }
            }
            else { //add atoms to current frame
                for (var i = 0; i < parsedAtoms.length; i++) {
                    this.addAtoms(parsedAtoms[i]);
                }
            }
        }

        for (let i = 0; i < this.frames.length; i++) {
            this.setAtomDefaults(this.frames[i]);
        }

        if (options.vibrate && options.vibrate.frames && options.vibrate.amplitude) {
            //fill in vibrational modes
            this.vibrate(options.vibrate.frames, options.vibrate.amplitude);
        }

        if (options.style) {
            this.setStyle({}, options.style);
        }
    };

    public setDontDuplicateAtoms(dup:boolean) {
        this.dontDuplicateAtoms = dup;
    };

    public setModelData(mData) {
        this.modelData = mData;
    };

    //return true if atom value matches property val
    private propertyMatches(atomval, val) {
        if (atomval == val) {
            return true;
        } else if (typeof (val) == 'string' && typeof (atomval) == 'number') {
            //support numerical integer ranges, e.g. resi: 3-7
            var match = val.match(/(-?\d+)\s*-\s*(-?\d+)/);
            if (match) {
                var lo = parseInt(match[1]);
                var hi = parseInt(match[2]);
                if (match && atomval >= lo && atomval <= hi) {
                    return true;
                }
            }
        }
        return false;
    };

    // make a deep copy of a selection object and create caches of expensive
    // selections.  We create a copy so caches are not attached to user
    // supplied objects where the user might change them invalidating the cache.
    // This does not support arbitrary
    // javascript objects, but support enough for eveything that is
    // used in selections: number, string, boolean, functions; as well
    // as arrays and nested objects with values of the aformentioned
    // types.
    private static deepCopyAndCache(selobject, model) {
        if (typeof selobject != 'object' || selobject == null) return selobject;
        if (selobject.__cache_created) return selobject; //already done
        const copy: any = {};
        for (const key in selobject) {
            const item = selobject[key];
            if (Array.isArray(item)) {
                // handle array separatly from other typeof == "object"
                // elements
                copy[key] = [];
                for (let i = 0; i < item.length; i++) {
                    copy[key].push(GLModel.deepCopyAndCache(item[i], model));
                }
            } else if (typeof item === "object" && key != "properties" && key != "model") {
                copy[key] = GLModel.deepCopyAndCache(item, model);
            } else {
                copy[key] = item;
            }

            //create caches of expensive selection types - the cache
            //stores the atoms matching the selection type
            if (key == "and" || key == "or") {
                // create a list of sets of matching atoms indexes for
                // each sub-selection
                const results = [];
                for (const subSelection of copy[key]) {
                    const set = new Set();
                    for (const match of model.selectedAtoms(subSelection)) {
                        set.add(match.index);
                    }
                    results.push(set);
                }

                if (key == "and") {
                    // get the intersection of two sets
                    const intersect = function (first, other) {
                        const result = new Set();
                        for (const elem of other) {
                            if (first.has(elem)) {
                                result.add(elem);
                            }
                        }
                        return result;
                    };

                    let intersection = new Set(results[0]);
                    for (const set of results.splice(1)) {
                        intersection = intersect(intersection, set);
                    }
                    copy[key].__cached_results = intersection;

                } else if (key == "or") {
                    const union = new Set();
                    for (const set of results) {
                        for (const elem of set) {
                            union.add(elem);
                        }
                    }

                    copy[key].__cached_results = union;
                }
            }

        }
        copy.__cache_created = true;
        return copy;
    };

    /** given a selection specification, return true if atom is selected.
     * Does not support context-aware selectors like expand/within/byres.
     *
     * @param {AtomSpec} atom
     * @param {AtomSelectionSpec} sel
     * @return {boolean}
     */
    public atomIsSelected(atom:AtomSpec, sel?:AtomSelectionSpec) {
        if (typeof (sel) === "undefined")
            return true; // undef gets all
        var invert = !!sel.invert;
        var ret = true;
        for (var key in sel) {
            if (key == "and" || key == "or" || key == "not") {  //boolean operators
                if (key == "not") {
                    if (this.atomIsSelected(atom, sel[key])) {
                        ret = false;
                        break;
                    }
                } else { //"or" and "and"
                    // these selections are expensive so when called via
                    //selectedAtoms shoudl be cached - but if atomIsSelected
                    //is called directly create the cache
                    if (sel[key].__cached_results === undefined) {
                        sel = GLModel.deepCopyAndCache(sel, this);
                    }

                    ret = sel[key].__cached_results.has(atom.index);
                    if (!ret) {
                        break;
                    }
                }

            } else if (key === 'predicate') { //a user supplied function for evaluating atoms
                if (!sel.predicate(atom)) {
                    ret = false;
                    break;
                }
            }
            else if (key == "properties" && atom[key]) {
                for (var propkey in sel.properties) {
                    if (propkey.startsWith("__cache")) continue;
                    if (typeof (atom.properties[propkey]) === 'undefined') {
                        ret = false;
                        break;
                    }
                    if (atom.properties[propkey] != sel.properties[propkey]) {
                        ret = false;
                        break;
                    }
                }
            }
            else if (sel.hasOwnProperty(key) && key != "props" &&
                key != "invert" && key != "model" && key != "byres" &&
                key != "expand" && key != "within" && key != "and" &&
                key != "or" && key != "not" && !key.startsWith('__cache')) {

                // if something is in sel, atom must have it
                if (typeof (atom[key]) === "undefined") {
                    ret = false;
                    break;
                }
                var isokay = false;
                if (key === "bonds") {
                    //special case counting number of bonds, for selecting nonbonded mostly
                    var val = sel[key];
                    if (val != atom.bonds.length) {
                        ret = false;
                        break;
                    }
                }
                else if (Array.isArray(sel[key])) {
                    // can be any of the listed values
                    var valarr = sel[key];
                    var atomval = atom[key];
                    for (let i = 0; i < valarr.length; i++) {
                        if (this.propertyMatches(atomval, valarr[i])) {
                            isokay = true;
                            break;
                        }
                    }
                    if (!isokay) {
                        ret = false;
                        break;
                    }
                } else { // single match
                    let val = sel[key];
                    if (!this.propertyMatches(atom[key], val)) {
                        ret = false;
                        break;
                    }
                }
            }
        }

        return invert ? !ret : ret;
    };


    private static squaredDistance(atom1:XYZ|AtomSpec, atom2:XYZ|AtomSpec) {
        var xd = atom2.x - atom1.x;
        var yd = atom2.y - atom1.y;
        var zd = atom2.z - atom1.z;
        return xd * xd + yd * yd + zd * zd;
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
    private expandAtomList(atomList: AtomSpec[], amt:number) {

        if (amt <= 0) return atomList;

        var pb = getExtent(atomList, undefined); // previous bounding box
        var nb = [[], [], []]; // expanded bounding box

        for (var i = 0; i < 3; i++) {
            nb[0][i] = pb[0][i] - amt;
            nb[1][i] = pb[1][i] + amt;
            nb[2][i] = pb[2][i];
        }

        // look in added box "shell" for new atoms
        var expand = [];
        for (let i = 0; i < this.atoms.length; i++) {

            var x = this.atoms[i].x;
            var y = this.atoms[i].y;
            var z = this.atoms[i].z;

            if (x >= nb[0][0] && x <= nb[1][0] && y >= nb[0][1] && y <= nb[1][1] && z >= nb[0][2] && z <= nb[1][2]) {
                if (!(x >= pb[0][0] && x <= pb[1][0] && y >= pb[0][1] && y <= pb[1][1] && z >= pb[0][2] && z <= pb[1][2])) {
                    expand.push(this.atoms[i]);
                }
            }
        }
        return expand;
    };



    private static getFloat(val: string|number): number {
        if(typeof(val) === 'number')
            return val;
        else
            return parseFloat(val);
    }
    /** return list of atoms selected by sel, this is specific to glmodel
     *
     * @param {AtomSelectionSpec} sel
     * @return {Object[]}
     * @example
     $3Dmol.download("pdb:4wwy",viewer,{},function(){
              var atoms = viewer.selectedAtoms({chain:'A'});
              for(var i = 0, n = atoms.length; i < n; i++) {
                 atoms[i].b = 0.0;
              }
              viewer.setStyle({cartoon:{colorscheme:{prop:'b',gradient: 'roygb',min:0,max:30}}});
              viewer.render();
          });
     */
    public selectedAtoms(sel: AtomSelectionSpec, from?: AtomSpec[]): AtomSpec[] {
        var ret = [];

        // make a copy of the selection to allow caching results without
        // the possibility for the user to change the selection and this
        // code not noticing the changes
        sel = GLModel.deepCopyAndCache(sel || {}, this);

        if (!from) from = this.atoms;
        var aLength = from.length;
        for (var i = 0; i < aLength; i++) {
            var atom = from[i];
            if (atom) {
                if (this.atomIsSelected(atom, sel))
                    ret.push(atom);
            }
        }

        // expand selection by some distance
        if (sel.hasOwnProperty("expand")) {
            // get atoms in expanded bounding box
            const exdist:number = GLModel.getFloat(sel.expand);
            let expand = this.expandAtomList(ret, exdist);
            let retlen = ret.length;
            const thresh = exdist * exdist;
            for (let i = 0; i < expand.length; i++) {
                for (let j = 0; j < retlen; j++) {

                    var dist = GLModel.squaredDistance(expand[i], ret[j]);
                    if (dist < thresh && dist > 0) {
                        ret.push(expand[i]);
                    }
                }
            }
        }

        // selection within distance of sub-selection
        if (sel.hasOwnProperty("within") && sel.within.hasOwnProperty("sel") &&
            sel.within.hasOwnProperty("distance")) {

            // get atoms in second selection
            var sel2 = this.selectedAtoms(sel.within.sel, this.atoms);
            var within = {};
            const dist = GLModel.getFloat(sel.within.distance);
            const thresh = dist * dist;
            for (let i = 0; i < sel2.length; i++) {
                for (let j = 0; j < ret.length; j++) {

                    let dist = GLModel.squaredDistance(sel2[i], ret[j]);
                    if (dist < thresh && dist > 0) {
                        within[j] = 1;
                    }
                }
            }
            var newret = [];
            if (sel.within.invert) {
                for (let j = 0; j < ret.length; j++) {
                    if (!within[j]) newret.push(ret[j]);
                }
            } else {
                for (let j in within) {
                    newret.push(ret[j]);
                }
            }
            ret = newret;
        }

        // byres selection flag
        if (sel.hasOwnProperty("byres")) {

            // Keep track of visited residues, visited atoms, and atom stack
            var vResis = {};
            var vAtoms = [];
            var stack = [];

            for (let i = 0; i < ret.length; i++) {

                // Check if atom is part of a residue, and that the residue hasn't been traversed yet
                let atom = ret[i];
                var c = atom.chain;
                var r = atom.resi;
                if (vResis[c] === undefined) vResis[c] = {};
                if (atom.hasOwnProperty("resi") && vResis[c][r] === undefined) {

                    // Perform a depth-first search of atoms with the same resi
                    vResis[c][r] = true;
                    stack.push(atom);
                    while (stack.length > 0) {
                        atom = stack.pop();
                        c = atom.chain;
                        r = atom.resi;
                        if (vAtoms[atom.index] === undefined) {
                            vAtoms[atom.index] = true;
                            for (var j = 0; j < atom.bonds.length; j++) {
                                var atom2 = this.atoms[atom.bonds[j]];
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


    /** Add list of new atoms to model.  Adjusts bonds appropriately.
     *
     * @param {AtomSpec[]} newatoms
     * @example
     * var atoms = [{elem: 'C', x: 0, y: 0, z: 0, bonds: [1,2], bondOrder: [1,2]}, {elem: 'O', x: -1.5, y: 0, z: 0, bonds: [0]},{elem: 'O', x: 1.5, y: 0, z: 0, bonds: [0], bondOrder: [2]}];

        viewer.setBackgroundColor(0xffffffff);
        var m = viewer.addModel();
        m.addAtoms(atoms);
        m.setStyle({},{stick:{}});
        viewer.zoomTo();
        viewer.render();
     */
    public addAtoms(newatoms: AtomSpec[]) {
        this.molObj = null;
        var start = this.atoms.length;
        var indexmap = [];
        // mapping from old index to new index
        var i;
        for (i = 0; i < newatoms.length; i++) {
            if (typeof (newatoms[i].index) == "undefined")
                newatoms[i].index = i;
            if (typeof (newatoms[i].serial) == "undefined")
                newatoms[i].serial = i;
            indexmap[newatoms[i].index] = start + i;
        }

        // copy and push newatoms onto atoms
        for (i = 0; i < newatoms.length; i++) {
            var olda = newatoms[i];
            var nindex = indexmap[olda.index];
            var a = extend({}, olda);
            a.index = nindex;
            a.bonds = [];
            a.bondOrder = [];
            a.model = this.id;
            a.style = a.style || deepCopy(GLModel.defaultAtomStyle);
            if (typeof (a.color) == "undefined")
                a.color = this.ElementColors[a.elem] || this.defaultColor;
            // copy over all bonds contained in selection,
            // updating indices appropriately
            var nbonds = olda.bonds ? olda.bonds.length : 0;
            for (var j = 0; j < nbonds; j++) {
                var neigh = indexmap[olda.bonds[j]];
                if (typeof (neigh) != "undefined") {
                    a.bonds.push(neigh);
                    a.bondOrder.push(olda.bondOrder ? olda.bondOrder[j] : 1);
                }
            }
            this.atoms.push(a);
        }
    };

    /** Assign bonds based on atomic coordinates.
     *  This currently uses a primitive distance-based algorithm that does not
     * consider valence constraints and will only create single bonds.
     */
    public assignBonds() {
        assignBonds(this.atoms);
    }

    /** Remove specified atoms from model
     *
     * @param {AtomSpec[]} badatoms - list of atoms
     */
    public removeAtoms(badatoms: AtomSpec[]) {
        this.molObj = null;
        // make map of all baddies
        var baddies = [];
        var i;
        for (i = 0; i < badatoms.length; i++) {
            baddies[badatoms[i].index] = true;
        }

        // create list of good atoms
        var newatoms = [];
        for (i = 0; i < this.atoms.length; i++) {
            var a = this.atoms[i];
            if (!baddies[a.index])
                newatoms.push(a);
        }

        // clear it all out
        this.atoms = [];
        // and add back in to get updated bonds
        this.addAtoms(newatoms);
    };


    /** Set atom style of selected atoms
     *
     * @param {AtomSelectionSpec} sel
     * @param {AtomStyleSpec} style
     * @param {boolean} add - if true, add to current style, don't replace
     @example
    $3Dmol.download("pdb:4UB9",viewer,{},function(){
              viewer.setBackgroundColor(0xffffffff);

              viewer.setStyle({chain:'A'},{line:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.setStyle({chain:'B'},{line:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.setStyle({chain:'C'},{cross:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.setStyle({chain:'D'},{cross:{colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.setStyle({chain:'E'},{cross:{radius:2.0,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.setStyle({chain:'F'},{stick:{hidden:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.RWB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.setStyle({chain:'G'},{stick:{radius:0.8,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.setStyle({chain:'H'},{stick:{singleBonds:true,colorscheme:{prop:'b',gradient: new $3Dmol.Gradient.ROYGB($3Dmol.getPropertyRange(viewer.selectedAtoms(),'b'))}}});
              viewer.render();
          });
     */
    public setStyle(sel:AtomSelectionSpec|AtomStyleSpec|string, style?:AtomStyleSpec|string, add?) {

        if (typeof (style) === 'undefined' && typeof (add) == 'undefined') {
            //if a single argument is provided, assume it is a style and select all
            style = sel as AtomStyleSpec|string;
            sel = {};
        }

        //if type is just a string, promote it to an object
        if (typeof (style) === 'string') {
            style = specStringToObject(style);
        }

        var changedAtoms = false;
        // somethings we only calculate if there is a change in a certain
        // style, although these checks will only catch cases where both
        // are either null or undefined
        var that = this;
        var setStyleHelper = function (atomArr) {
            var selected = that.selectedAtoms(sel as AtomSelectionSpec, atomArr);
            for (let i = 0; i < atomArr.length; i++) {
                if (atomArr[i]) atomArr[i].capDrawn = false; //reset for proper stick render
            }

            for (let i = 0; i < selected.length; i++) {
                changedAtoms = true;
                if (selected[i].clickable || selected[i].hoverable)
                    selected[i].intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };


                if (!add) selected[i].style = {};
                for (let s in style as AtomStyleSpec) {
                    if (style.hasOwnProperty(s)) {
                        selected[i].style[s] = selected[i].style[s] || {}; //create distinct object for each atom
                        Object.assign(selected[i].style[s], style[s]);
                    }
                }
            }
        };

        setStyleHelper(this.atoms);
        for (var i = 0; i < this.frames.length; i++) {
            if (this.frames[i] !== this.atoms) setStyleHelper(this.frames[i]);
        }

        if (changedAtoms)
            this.molObj = null; //force rebuild

    };

    /** Set clickable and callback of selected atoms
     *
     * @param {AtomSelectionSpec} sel - atom selection to apply clickable settings to
     * @param {boolean} clickable - whether click-handling is enabled for the selection
     * @param {function} callback - function called when an atom in the selection is clicked

     */
    public setClickable(sel:AtomSelectionSpec, clickable:boolean, callback) {

        // make sure clickable is a boolean
        clickable = !!clickable;
        callback = makeFunction(callback);
        if (callback === null) {
            console.log("Callback is not a function");
            return;
        }

        var selected = this.selectedAtoms(sel, this.atoms);
        var len = selected.length;
        for (let i = 0; i < len; i++) {

            selected[i].intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };
            selected[i].clickable = clickable;
            if (callback) selected[i].callback = callback;

        }

        if (len > 0) this.molObj = null; // force rebuild to get correct intersection shapes
    };

    /** Set hoverable and callback of selected atoms
    *
    * @param {AtomSelectionSpec} sel - atom selection to apply hoverable settings to
    * @param {boolean} hoverable - whether hover-handling is enabled for the selection
    * @param {function} hover_callback - function called when an atom in the selection is hovered over
    * @param {function} unhover_callback - function called when the mouse moves out of the hover area
    */
    public setHoverable(sel:AtomSelectionSpec, hoverable:boolean, hover_callback, unhover_callback) {

        // make sure hoverable is a boolean
        hoverable = !!hoverable;
        hover_callback = makeFunction(hover_callback);
        unhover_callback = makeFunction(unhover_callback);

        // report to console if hover_callback is not a valid function
        if (hover_callback === null) {
            console.log("Hover_callback is not a function");
            return;
        }
        // report to console if unhover_callback is not a valid function
        if (unhover_callback === null) {
            console.log("Unhover_callback is not a function");
            return;
        }

        var selected = this.selectedAtoms(sel, this.atoms);
        var len = selected.length;
        for (let i = 0; i < len; i++) {

            selected[i].intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };
            selected[i].hoverable = hoverable;
            if (hover_callback) selected[i].hover_callback = hover_callback;
            if (unhover_callback) selected[i].unhover_callback = unhover_callback;

        }

        if (len > 0) this.molObj = null; // force rebuild to get correct intersection shapes
    };

    /** enable context menu of selected atoms
     *
     * @param {AtomSelectionSpec} sel - atom selection to apply hoverable settings to
     * @param {boolean} contextMenuEnabled - whether contextMenu-handling is enabled for the selection
     */
    public enableContextMenu(sel:AtomSelectionSpec, contextMenuEnabled) {
        // make sure contextMenuEnabled is a boolean
        contextMenuEnabled = !!contextMenuEnabled;

        var i;
        var selected = this.selectedAtoms(sel, this.atoms);
        var len = selected.length;
        for (i = 0; i < len; i++) {

            selected[i].intersectionShape = { sphere: [], cylinder: [], line: [], triangle: [] };
            selected[i].contextMenuEnabled = contextMenuEnabled;
        }

        if (len > 0) this.molObj = null; // force rebuild to get correct intersection shapes
    };

    /** given a mapping from element to color, set atom colors
     *
     * @param {AtomSelectionSpec} sel
     * @param {object} colors
     */
    public setColorByElement(sel: AtomSelectionSpec, colors) {

        if (this.molObj !== null && GLModel.sameObj(colors, this.lastColors))
            return; // don't recompute
        this.lastColors = colors;
        var atoms = this.selectedAtoms(sel, atoms);
        if (atoms.length > 0)
            this.molObj = null; // force rebuild
        for (var i = 0; i < atoms.length; i++) {
            var a = atoms[i];
            if (typeof (colors[a.elem]) !== "undefined") {
                a.color = colors[a.elem];
            }
        }
    };

    /**
     * @param {AtomSelectionSpec} sel
     * @param {string} prop
     * @param {Gradient|string} scheme
     */
    public setColorByProperty(sel: AtomSelectionSpec, prop: string, scheme: Gradient|string, range?) {
        var i, a;
        var atoms = this.selectedAtoms(sel, atoms);
        this.lastColors = null; // don't bother memoizing
        if (atoms.length > 0)
            this.molObj = null; // force rebuild

        if (typeof scheme === 'string' && typeof (Gradient.builtinGradients[scheme]) != "undefined") {
            scheme = new Gradient.builtinGradients[scheme]();
        }
        scheme = scheme as Gradient;
        if (!range) { //no explicit range, get from scheme
            range = scheme.range();
        }

        if (!range) { //no range in scheme, compute the range for this model
            range = getPropertyRange(atoms, prop);
        }
        // now apply colors using scheme
        for (i = 0; i < atoms.length; i++) {
            a = atoms[i];
            var val = getAtomProperty(a, prop);
            if (val != null) {
                a.color = scheme.valueToHex(parseFloat(a.properties[prop]), range);
            }
        }
    };

    /**
     * @deprecated use setStyle and colorfunc attribute
     * @param {AtomSelectionSpec} sel - selection object
     * @param {function} func - function to be used to set the color
     @example
      $3Dmol.download("pdb:4UAA",viewer,{},function(){
              viewer.setBackgroundColor(0xffffffff);
              var colorAsSnake = function(atom) {
                return atom.resi % 2 ? 'white': 'green'
              };

              viewer.setStyle( {}, { cartoon: {colorfunc: colorAsSnake }});

              viewer.render();
          });

     */
    public setColorByFunction(sel: AtomSelectionSpec, colorfun) {
        var atoms = this.selectedAtoms(sel, atoms);
        if (typeof (colorfun) !== 'function')
            return;
        this.lastColors = null; // don't bother memoizing
        if (atoms.length > 0)
            this.molObj = null; // force rebuild

        // now apply colorfun
        for (let i = 0; i < atoms.length; i++) {
            let a = atoms[i];
            a.color = colorfun(a);
        }
    };

    /** Convert the model into an object in the format of a ChemDoodle JSON model.
     *
     * @param {boolean} whether or not to include style information. Defaults to false.
     * @return {Object}
     */
    public toCDObject(includeStyles: boolean = false) {
        var out: any = { a: [], b: [] };
        if (includeStyles) {
            out.s = [];
        }
        for (let i = 0; i < this.atoms.length; i++) {
            let atomJSON: any = {};
            let atom = this.atoms[i];
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

            for (let b = 0; b < atom.bonds.length; b++) {
                let firstAtom = i;
                let secondAtom = atom.bonds[b];
                if (firstAtom >= secondAtom)
                    continue;
                let bond: any = {
                    b: firstAtom,
                    e: secondAtom
                };
                let bondOrder = atom.bondOrder[b];
                if (bondOrder != 1) {
                    bond.o = bondOrder;
                }
                out.b.push(bond);
            }
        }
        return out;
    };


    /** manage the globj for this model in the possed modelGroup - if it has to be regenerated, remove and add
     *
     * @param {Object3D} group
     * @param Object options
     */
    public globj(group, options) {
        if (this.molObj === null || options.regen) { // have to regenerate
            this.molObj = this.createMolObj(this.atoms, options);
            if (this.renderedMolObj) { // previously rendered, remove
                group.remove(this.renderedMolObj);
                this.renderedMolObj = null;
            }
            this.renderedMolObj = this.molObj.clone();
            if (this.hidden) {
                this.renderedMolObj.setVisible(false);
                this.molObj.setVisible(false);
            }
            group.add(this.renderedMolObj);
        }
    };

    /** return a VRML string representation of the model.  Does not include VRML header information
     * @return VRML
     */
    public exportVRML() {
        //todo: export spheres and cylinder objects instead of all mesh
        var tmpobj = this.createMolObj(this.atoms, { supportsImposters: false, supportsAIA: false });
        return tmpobj.vrml();
    };

    /** Remove any renderable mol object from scene
     *
     * @param {Object3D} group
     */
    public removegl(group) {
        if (this.renderedMolObj) {
            //dispose of geos and materials
            if (this.renderedMolObj.geometry !== undefined) this.renderedMolObj.geometry.dispose();
            if (this.renderedMolObj.material !== undefined) this.renderedMolObj.material.dispose();
            group.remove(this.renderedMolObj);
            this.renderedMolObj = null;
        }
        this.molObj = null;
    };

    /**
     * Don't show this model in future renderings. Keep all styles and state
     * so it can be efficiencly shown again.
     *
     * * @see GLModel#show

     * @example
        $3Dmol.download("pdb:3ucr",viewer,{},function(){
        viewer.setStyle({},{stick:{}});
        viewer.getModel().hide();
        viewer.render();
        });
     */
    public hide() {
        this.hidden = true;
        if (this.renderedMolObj) this.renderedMolObj.setVisible(false);
        if (this.molObj) this.molObj.setVisible(false);
    };

    /**
     * Unhide a hidden model
     * @see GLModel#hide
     * @example
        $3Dmol.download("pdb:3ucr",viewer,{},function(){
        viewer.setStyle({},{stick:{}});
        viewer.getModel().hide();
        viewer.render(  )
        viewer.getModel().show()
        viewer.render();
        });
     */
    public show() {
        this.hidden = false;
        if (this.renderedMolObj) this.renderedMolObj.setVisible(true);
        if (this.molObj) this.molObj.setVisible(true);
    };


    /** Create labels for atoms that show the value of the passed property.
     *
     * @param {String} prop - property name
     * @param {AtomSelectionSpec} sel
     * @param {GLViewer} viewer
     * @param {LabelSpec} options
     */
    public addPropertyLabels(prop: string, sel: AtomSelectionSpec, viewer: GLViewer, style: LabelSpec) {
        var atoms = this.selectedAtoms(sel, atoms);
        var mystyle = deepCopy(style);
        for (var i = 0; i < atoms.length; i++) {
            var a = atoms[i];
            var label = null;
            if (typeof (a[prop]) != 'undefined') {
                label = String(a[prop]);
            } else if (typeof (a.properties[prop]) != 'undefined') {
                label = String(a.properties[prop]);
            }

            if (label != null) {
                mystyle.position = a;
                viewer.addLabel(label, mystyle);
            }
        }
    };


    /** Create labels for residues of selected atoms.
     * Will create a single label at the center of mass of all atoms
     * with the same chain,resn, and resi.
     *
     * @param {AtomSelectionSpec} sel
     * @param {GLViewer} viewer
     * @param {LabelSpec} options
     * @param {boolean} byframe - if true, create labels for every individual frame, not just current; frames must be loaded already
     */
    public addResLabels(sel: AtomSelectionSpec, viewer: GLViewer, style: LabelSpec, byframe:boolean=false) {

        var created_labels = [];
        var helper = function (model, framenum?) {
            var atoms = model.selectedAtoms(sel, atoms);
            var bylabel = {};
            //collect by chain:resn:resi
            for (var i = 0; i < atoms.length; i++) {
                var a = atoms[i];
                var c = a.chain;
                var resn = a.resn;
                var resi = a.resi;
                var label = resn + '' + resi;
                if (!bylabel[c]) bylabel[c] = {};
                if (!bylabel[c][label]) bylabel[c][label] = [];
                bylabel[c][label].push(a);
            }

            var mystyle = deepCopy(style);
            //now compute centers of mass
            for (let c in bylabel) {
                if (bylabel.hasOwnProperty(c)) {
                    var labels = bylabel[c];
                    for (let label in labels) {
                        if (labels.hasOwnProperty(label)) {
                            let atoms = labels[label];
                            let sum = new Vector3(0, 0, 0);
                            for (let i = 0; i < atoms.length; i++) {
                                let a = atoms[i];
                                sum.x += a.x;
                                sum.y += a.y;
                                sum.z += a.z;
                            }
                            sum.divideScalar(atoms.length);
                            mystyle.position = sum;
                            mystyle.frame = framenum;
                            let l = viewer.addLabel(label, mystyle, undefined, true);
                            created_labels.push(l);
                        }
                    }
                }
            }
        };

        if (byframe) {
            var n = this.getNumFrames();
            let savedatoms = this.atoms;
            for (let i = 0; i < n; i++) {
                if (this.frames[i]) {
                    this.atoms = this.frames[i];
                    helper(this, i);
                }
            }
            this.atoms = savedatoms;
        } else {
            helper(this);
        }
        return created_labels;
    };


    //recurse over the current atoms to establish a depth first order
    private setupDFS() {
        this.atomdfs = [];
        var self = this;
        var visited = new Int8Array(this.atoms.length);
        visited.fill(0);

        var search = function (i, prev, component) {
            //add i to component and recursive explore connected atoms
            component.push([i, prev]);
            var atom = self.atoms[i];
            visited[i] = 1;
            for (var b = 0; b < atom.bonds.length; b++) {
                var nexti = atom.bonds[b];
                if (self.atoms[nexti] && !visited[nexti]) {
                    search(nexti, i, component);
                }
            }
        };

        for (var i = 0; i < this.atoms.length; i++) {
            var atom = this.atoms[i];
            if (atom && !visited[i]) {
                var component = [];
                search(i, -1, component);
                this.atomdfs.push(component);
            }
        }
    };

    /**
    * Set coordinates from remote trajectory file.
    * @param {string} url - contains the url where mdsrv has been hosted
    * @param {string} path - contains the path of the file (<root>/filename)
    * @return {Promise}
    */
    public setCoordinatesFromURL(url:string, path:string) {
        this.frames = [];
        var self = this;
        if (this.box) this.setupDFS();
        if (!url.startsWith('http'))
            url = 'http://' + url;
        return get(url + "/traj/numframes/" + path, function (numFrames) {
            if (!isNaN(parseInt(numFrames))) {
                self.frames.push(self.atoms);
                self.frames.numFrames = numFrames;
                self.frames.url = url;
                self.frames.path = path;
                return self.setFrame(0);
            }
        });
    };


    /**
    * Set coordinates for the atoms from provided trajectory file.
    * @param {string|ArrayBuffer} str - contains the data of the file
    * @param {string} format - contains the format of the file (mdcrd, inpcrd, pdb, netcdf, or array).  Arrays should be TxNx3 where T is the number of timesteps and N the number of atoms.
      @example
         let m = viewer.addModel()  //create an empty model
         m.addAtoms([{x:0,y:0,z:0,elem:'C'},{x:2,y:0,z:0,elem:'C'}]) //provide a list of dictionaries representing the atoms
         viewer.setStyle({'sphere':{}})
         m.setCoordinates([[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [2.8888888359069824, 0.0, 0.0]], [[0.0, 0.0, 0.0], [3.777777671813965, 0.0, 0.0]], [[0.0, 0.0, 0.0], [4.666666507720947, 0.0, 0.0]], [[0.0, 0.0, 0.0], [5.55555534362793, 0.0, 0.0]], [[0.0, 0.0, 0.0], [6.44444465637207, 0.0, 0.0]], [[0.0, 0.0, 0.0], [7.333333492279053, 0.0, 0.0]], [[0.0, 0.0, 0.0], [8.222222328186035, 0.0, 0.0]], [[0.0, 0.0, 0.0], [9.11111068725586, 0.0, 0.0]], [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]],'array');
         viewer.animate({loop: "forward",reps: 1});
         viewer.zoomTo();
         viewer.zoom(0.5);
         viewer.render();
    */

    public setCoordinates(str:string|ArrayBuffer, format:string) {
        format = format || "";
        if (!str)
            return []; // leave an empty model

        if (/\.gz$/.test(format)) {
            // unzip gzipped files
            format = format.replace(/\.gz$/, '');
            try {
                str = inflate(str, {
                    to: 'string'
                });
            } catch (err) {
                console.log(err);
            }
        }
        var supportedFormats = { "mdcrd": "", "inpcrd": "", "pdb": "", "netcdf": "", "array": "" };
        if (supportedFormats.hasOwnProperty(format)) {
            this.frames = [];
            var atomCount = this.atoms.length;
            var values = GLModel.parseCrd(str, format);
            var count = 0;
            while (count < values.length) {
                var temp = [];
                for (var i = 0; i < atomCount; i++) {
                    var newAtom = {};
                    for (var k in this.atoms[i]) {
                        newAtom[k] = this.atoms[i][k];
                    }
                    temp[i] = newAtom;
                    temp[i].x = values[count++];
                    temp[i].y = values[count++];
                    temp[i].z = values[count++];
                }

                this.frames.push(temp);
            }
            this.atoms = this.frames[0];
            return this.frames;
        }
        return [];
    };

    /**
     * add atomSpecs to validAtomSelectionSpecs
     * @deprecated
     * @param {Array} customAtomSpecs - array of strings that can be used as atomSelectionSpecs
     * this is to prevent the 'Unknown Selector x' message on the console for the strings passed.
     * These messages are no longer generated as, in theory, typescript will catch problems at compile time.
     * In practice, there may still be issues at run-time but we don't check for them...
     *
     * What we should do is use something like https://github.com/woutervh-/typescript-is to do runtime
     * type checking, but it currently doesn't work with our types...
     */

    public addAtomSpecs(customAtomSpecs) {

    };

    static parseCrd(data, format:string) {
        var values = []; // this will contain the all the float values in the
        // file.
        var counter = 0;
        if (format == "pdb") {
            var index = data.indexOf("\nATOM");
            while (index != -1) {
                while (data.slice(index, index + 5) == "\nATOM" ||
                    data.slice(index, index + 7) == "\nHETATM") {
                    values[counter++] = parseFloat(data.slice(index + 31,
                        index + 39));
                    values[counter++] = parseFloat(data.slice(index + 39,
                        index + 47));
                    values[counter++] = parseFloat(data.slice(index + 47,
                        index + 55));
                    index = data.indexOf("\n", index + 54);
                    if (data.slice(index, index + 4) == "\nTER")
                        index = data.indexOf("\n", index + 5);
                }
                index = data.indexOf("\nATOM", index);
            }

        } else if (format == "netcdf") {
            var reader = new NetCDFReader(data);
            values = [].concat.apply([], reader.getDataVariable('coordinates'));

        } else if (format == "array" || Array.isArray(data)) {
            return data.flat(2);
        } else {
            let index = data.indexOf("\n"); // remove the first line containing title
            if (format == 'inpcrd') {
                index = data.indexOf("\n", index + 1); //remove second line w/#atoms
            }

            data = data.slice(index + 1);
            values = data.match(/\S+/g).map(parseFloat);
        }
        return values;
    };

    static parseMolData(data?, format: string="", options?) {
        if (!data)
            return []; //leave an empty model

        if (/\.gz$/.test(format)) {
            //unzip gzipped files
            format = format.replace(/\.gz$/, '');
            try {
                data = inflate(data, { to: 'string' });
            } catch (err) {
                console.log(err);
            }
        }

        if (typeof (Parsers[format]) == "undefined") {
            // let someone provide a file name and get format from extension
            format = format.split('.').pop();
            if (typeof (Parsers[format]) == "undefined") {
                console.log("Unknown format: " + format);
                // try to guess correct format from data contents
                if (data instanceof Uint8Array) {
                    format = "mmtf"; //currently only supported binary format?
                } else if (data.match(/^@<TRIPOS>MOLECULE/gm)) {
                    format = "mol2";
                } else if (data.match(/^data_/gm) && data.match(/^loop_/gm)) {
                    format = "cif";
                } else if (data.match(/^HETATM/gm) || data.match(/^ATOM/gm)) {
                    format = "pdb";
                } else if (data.match(/ITEM: TIMESTEP/gm)) {
                    format = "lammpstrj";
                } else if (data.match(/^.*\n.*\n.\s*(\d+)\s+(\d+)/gm)) {
                    format = "sdf"; // could look at line 3
                } else if (data.match(/^%VERSION\s+VERSION_STAMP/gm)) {
                    format = "prmtop";
                } else {
                    format = "xyz";
                }
                console.log("Best guess: " + format);
            }
        }
        var parse = Parsers[format];
        var parsedAtoms = parse(data, options);

        return parsedAtoms;
    };

}

/** Atom style specification */
export interface AtomStyleSpec {
    /** draw bonds as lines */
    line?: LineStyleSpec;
    /** draw atoms as crossed lines (aka stars) */
    cross?: CrossStyleSpec;
    /** draw bonds as capped cylinders */
    stick?: StickStyleSpec;
    /** draw atoms as spheres */
    sphere?: SphereStyleSpec;
    /** draw cartoon representation of secondary structure */
    cartoon?: CartoonStyleSpec;
    /** invisible style for click handling only */
    clicksphere?: ClickSphereStyleSpec;
};

/** Line style specification
 */
export interface LineStyleSpec {
    /** do not show line */
    hidden?: boolean;
    /** *deprecated due to vanishing browser support*  */
    linewidth?: number;
    /** colorscheme to use on atoms */
    colorscheme?: ColorschemeSpec;
    /** fixed coloring, overrides colorscheme */
    color?: ColorSpec;
    /** opacity (zero to one), must be the same for all atoms in a model */
    opacity?: number;
    /** wireframe style */
    wireframe?: boolean;

}

/** Cross style specification
 */
export interface CrossStyleSpec {
    /** do not show line */
    hidden?: boolean;
    /** *deprecated due to vanishing browser support*  */
    linewidth?: number;
    /** radius of cross */
    radius?: number;
    /** scale VDW radius by specified amount */
    scale?: number;
    /** colorscheme to use on atoms */
    colorscheme?: ColorschemeSpec;
    /** fixed coloring, overrides colorscheme */
    color?: ColorSpec;
    /** opacity (zero to one), must be the same for all atoms in a model */
    opacity?: number;
}

/** Stick (cylinder) style specification
 */
export interface StickStyleSpec {
    /** do not show sticks */
    hidden?: boolean;
    /** radius of stick */
    radius?: number;
    /** draw all bonds as single bonds */
    singleBonds?: boolean;
    /** colorscheme to use on atoms */
    colorscheme?: ColorschemeSpec;
    /** fixed coloring, overrides colorscheme */
    color?: ColorSpec;
    /** opacity (zero to one), must be the same for all atoms in a model */
    opacity?: number;
    /** display nonbonded atoms as spheres */
    showNonBonded?: boolean;
}


/** Sphere (spacefill) style specification
 */
export interface SphereStyleSpec {
    /** do not show sticks */
    hidden?: boolean;
    /** fixed radius of sphere */
    radius?: number;
    /** scale VDW radius by specified amount */
    scale?: number;
    /** colorscheme to use on atoms */
    colorscheme?: ColorschemeSpec;
    /** fixed coloring, overrides colorscheme */
    color?: ColorSpec;
    /** opacity (zero to one), must be the same for all atoms in a model */
    opacity?: number;
}

/** Invisible click sphere style specification.  This lets you set
 * larger (or smaller) click targets on atoms then the default radii or
 * have clickable atoms even if they aren't being rendered visibly.
 */
export interface ClickSphereStyleSpec {
    /** do not show sticks */
    hidden?: boolean;
    /** fixed radius of sphere */
    radius?: number;
    /** scale VDW radius by specified amount */
    scale?: number;
}

/** Style for individual bond. */
export interface BondStyle {
    iswire?: boolean;
    /**  */
    singleBond?: boolean;
    /**  */
    radius?: number;
    /**  */
    color1?: ColorSpec;
    /**  */
    color2?: ColorSpec;
}
