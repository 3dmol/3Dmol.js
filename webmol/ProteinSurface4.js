/*  ProteinSurface.js by biochem_fan

Ported and modified for Javascript based on EDTSurf,
  whose license is as follows.

Permission to use, copy, modify, and distribute this program for any
purpose, with or without fee, is hereby granted, provided that this
copyright notice and the reference information appear in all copies or
substantial portions of the Software. It is provided "as is" without
express or implied warranty. 

Reference:
http://zhanglab.ccmb.med.umich.edu/EDTSurf/
D. Xu, Y. Zhang (2009) Generating Triangulated Macromolecular Surfaces
by Euclidean Distance Transform. PLoS ONE 4(12): e8140.

=======

TODO: Improved performance on Firefox
      Reduce memory consumption
      Refactor!
 */


// dkoes
// Surface calculations.  This must be safe to use within a web worker.
if (typeof console === 'undefined') {
    // this should only be true inside of a webworker
    console = {
        log : function() {
        }
    };
}

var ProteinSurface = function() {
    // constants for vpbits bitmasks
    var INOUT = 1;
    var ISDONE = 2;
    var ISBOUND = 4;

    var ptranx = 0, ptrany = 0, ptranz = 0;
    var probeRadius = 1.4;
    var defaultScaleFactor = 2;
    var scaleFactor = defaultScaleFactor; // 2 is .5A grid; if this is made user configurable,
                                                    // also have to adjust offset used to find non-shown
                                                    // atoms
    var nZ = 0, nY = 0, nX = 0;
    var cutRadius = 0;
    var vpBits = null; // uint8 array of bitmasks
    var vpDistance = null; // floatarray
    var vpAtomID = null; // intarray
    var vertnumber = 0, facenumber = 0;
    var pminx = 0, pminy = 0, pminz = 0, pmaxx = 0, pmaxy = 0, pmaxz = 0;

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
                    "NI" : 1.63,
                    "X" : 2
            };

    var getVDWIndex = function(atom) {
        if(!atom.elem || typeof(vdwRadii[atom.elem]) === "undefined") {
            return "X";
        }
        return atom.elem;
    };

    var depty = {}, widxz = {};
    var faces, verts;
    var nb = [ [ 1, 0, 0 ], [ -1, 0, 0 ], [ 0, 1, 0 ], [ 0, -1, 0 ],
                    [ 0, 0, 1 ], [ 0, 0, -1 ], [ 1, 1, 0 ], [ 1, -1, 0 ], [ -1, 1, 0 ],
                    [ -1, -1, 0 ], [ 1, 0, 1 ], [ 1, 0, -1 ], [ -1, 0, 1 ],
                    [ -1, 0, -1 ], [ 0, 1, 1 ], [ 0, 1, -1 ], [ 0, -1, 1 ],
                    [ 0, -1, -1 ], [ 1, 1, 1 ], [ 1, 1, -1 ], [ 1, -1, 1 ],
                    [ -1, 1, 1 ], [ 1, -1, -1 ], [ -1, -1, 1 ], [ -1, 1, -1 ],
                    [ -1, -1, -1 ] ];

    var origextent;

    var inOrigExtent = function(x, y, z) {
        if (x < origextent[0][0] || x > origextent[1][0])
            return false;
        if (y < origextent[0][1] || y > origextent[1][1])
            return false;
        if (z < origextent[0][2] || z > origextent[1][2])
            return false;
        return true;
    };

    this.getFacesAndVertices = function(atomlist) {
        var atomsToShow = new Object();
        for ( var i = 0, lim = atomlist.length; i < lim; i++)
            atomsToShow[atomlist[i]] = true;
        var vertices = verts;
        for (i = 0; i < vertices.length; i++) {
            vertices[i].x = vertices[i].x / scaleFactor - ptranx;
            vertices[i].y = vertices[i].y / scaleFactor - ptrany;
            vertices[i].z = vertices[i].z / scaleFactor - ptranz;
        }

        var finalfaces = [];
        for ( var i = 0; i < faces.length; i += 3) {
            //var f = faces[i];
            var fa = faces[i], fb = faces[i+1], fc = faces[i+2];
            var a = vertices[fa].atomid, b = vertices[fb].atomid, c = vertices[fc].atomid;

            // must be a unique face for each atom
            var which = a;
            if (b < which)
                which = b;
            if (c < which)
                which = c;
            if (!atomsToShow[which]) {
                continue;
            }
            var av = vertices[faces[i]];
            var bv = vertices[faces[i+1]];
            var cv = vertices[faces[i+2]];

            if (fa !== fb && fb !== fc && fa !== fc)
                finalfaces.push(fa), finalfaces.push(fb), finalfaces.push(fc);
        }

        //try to help the garbage collector
        vpBits = null; // uint8 array of bitmasks
        vpDistance = null; // floatarray
        vpAtomID = null; // intarray
        
        return {
            vertices : vertices,
            faces : finalfaces
        };
    };

    //btype is true if not type 1 (VDW) ...
    //TODO: Is this right (probe radius seems like an SAS surface...)?
    this.initparm = function(extent, btype, volume) {
        if(volume > 1000000) //heuristical decrease resolution to avoid large memory consumption
            scaleFactor = defaultScaleFactor/2;

        var margin = (1 / scaleFactor) * 5.5; // need margin to avoid
                                                                                        // boundary/round off effects
        origextent = extent;
        pminx = extent[0][0], pmaxx = extent[1][0];
        pminy = extent[0][1], pmaxy = extent[1][1];
        pminz = extent[0][2], pmaxz = extent[1][2];

        if (!btype) {
            pminx -= margin;
            pminy -= margin;
            pminz -= margin;
            pmaxx += margin;
            pmaxy += margin;
            pmaxz += margin;
        } else {
            pminx -= probeRadius + margin;
            pminy -= probeRadius + margin;
            pminz -= probeRadius + margin;
            pmaxx += probeRadius + margin;
            pmaxy += probeRadius + margin;
            pmaxz += probeRadius + margin;
        }

        pminx = Math.floor(pminx * scaleFactor) / scaleFactor;
        pminy = Math.floor(pminy * scaleFactor) / scaleFactor;
        pminz = Math.floor(pminz * scaleFactor) / scaleFactor;
        pmaxx = Math.ceil(pmaxx * scaleFactor) / scaleFactor;
        pmaxy = Math.ceil(pmaxy * scaleFactor) / scaleFactor;
        pmaxz = Math.ceil(pmaxz * scaleFactor) / scaleFactor;

        ptranx = -pminx;
        ptrany = -pminy;
        ptranz = -pminz;

        nX = Math.ceil(scaleFactor * (pmaxx - pminx)) + 1;
        nY = Math.ceil(scaleFactor * (pmaxy - pminy)) + 1;
        nZ = Math.ceil(scaleFactor * (pmaxz - pminz)) + 1;

        this.boundingatom(btype);
        cutRadius = probeRadius * scaleFactor;

        vpBits = new Uint8Array(nX * nY * nZ);
        vpDistance = new Float64Array(nX * nY * nZ); // float 32
        // doesn't
        // play
        // nicely
        // with
        // native
        // floats
        vpAtomID = new Int32Array(nX * nY * nZ);
        console.log("Box size: ", nX, nY, nZ, vpBits.length);
    };

    this.boundingatom = function(btype) {
        var tradius = [];
        var txz, tdept, sradius;

        for ( var i in vdwRadii) {
            if(!vdwRadii.hasOwnProperty(i))
                continue;
            var r = vdwRadii[i];
            if (!btype)
                tradius[i] = r * scaleFactor + 0.5;
            else
                tradius[i] = (r + probeRadius) * scaleFactor + 0.5;

            sradius = tradius[i] * tradius[i];
            widxz[i] = Math.floor(tradius[i]) + 1;
            depty[i] = new Array(widxz[i] * widxz[i]);
            var indx = 0;
            for (var j = 0; j < widxz[i]; j++) {
                for (var k = 0; k < widxz[i]; k++) {
                    txz = j * j + k * k;
                    if (txz > sradius)
                        depty[i][indx] = -1; // outside
                    else {
                        tdept = Math.sqrt(sradius - txz);
                        depty[i][indx] = Math.floor(tdept);
                    }
                    indx++;
                }
            }
        }
    };

    this.fillvoxels = function(atoms, atomlist) { // (int seqinit,int
        // seqterm,bool
        // atomtype,atom*
        // proseq,bool bcolor)
        var i;
        for ( var i = 0, lim = vpBits.length; i < lim; i++) {
            vpBits[i] = 0;
            vpDistance[i] = -1.0;
            vpAtomID[i] = -1;
        }

        for (i in atomlist) {
            var atom = atoms[atomlist[i]];
            if (atom === undefined)
                continue;
            this.fillAtom(atom, atoms);
        }

        for (i = 0, lim = vpBits.length; i < lim; i++)
            if (vpBits[i] & INOUT)
                vpBits[i] |= ISDONE;

    };


    this.fillAtom = function(atom, atoms) {
        var cx, cy, cz, ox, oy, oz;
        cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
        cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
        cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

        var at = getVDWIndex(atom);
        var nind = 0;
        var cnt = 0;

        for (var i = 0, n = widxz[at]; i < n; i++) {
            for (var j = 0; j < n; j++) {
                if (depty[at][nind] !== -1) {
                    for (var ii = -1; ii < 2; ii++) {
                        for (var jj = -1; jj < 2; jj++) {
                            for (var kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    var mi = ii * i;
                                    var mk = kk * j;
                                    for (var k = 0; k <= depty[at][nind]; k++) {
                                        var mj = k * jj;
                                        var si = cx + mi;
                                        var sj = cy + mj;
                                        var sk = cz + mk;
                                        if (si < 0 || sj < 0 || sk < 0
                                        || si >= nX
                                        || sj >= nY
                                        || sk >= nZ)
                                            continue;
                                        var index = si * nY * nZ + sj
                                                        * nZ + sk;

                                        if (!(vpBits[index] & INOUT)) {
                                            vpBits[index] |= INOUT;
                                            vpAtomID[index] = atom.serial;
                                        } 
                                        else {
                                            var atom2 = atoms[vpAtomID[index]];
                                            ox = Math.floor(0.5 + scaleFactor
                                                            * (atom2.x + ptranx));
                                            oy = Math.floor(0.5 + scaleFactor
                                                            * (atom2.y + ptrany));
                                            oz = Math.floor(0.5 + scaleFactor
                                                            * (atom2.z + ptranz));
                                            if (mi * mi + mj * mj + mk * mk < ox
                                                            * ox + oy * oy + oz * oz)
                                                vpAtomID[index] = atom.serial;
                                        }

                                    }// k
                                }// if
                            }// kk
                        }// jj
                    }// ii
                }// if
                nind++;
            }// j
        }// i
    };

    this.fillvoxelswaals = function(atoms, atomlist) {
        for ( var i = 0, lim = vpBits.length; i < lim; i++)
            vpBits[i] &= ~ISDONE; // not isdone

        for (i in atomlist) {
            var atom = atoms[atomlist[i]];
            if (atom == undefined)
                continue;

            this.fillAtomWaals(atom, atoms);
        }
    };

    this.fillAtomWaals = function(atom, atoms) {
        var cx, cy, cz, ox, oy, oz, nind = 0;
        cx = Math.floor(0.5 + scaleFactor * (atom.x + ptranx));
        cy = Math.floor(0.5 + scaleFactor * (atom.y + ptrany));
        cz = Math.floor(0.5 + scaleFactor * (atom.z + ptranz));

        var at = getVDWIndex(atom);
        var i, j, k, n, ii, jj, kk;
        for (i = 0, n = widxz[at]; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (depty[at][nind] !== -1) {
                    for (ii = -1; ii < 2; ii++) {
                        for (jj = -1; jj < 2; jj++) {
                            for (kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    var mi = ii * i;
                                    var mk = kk * j;
                                    for (k = 0; k <= depty[at][nind]; k++) {
                                        var mj = k * jj;
                                        var si = cx + mi;
                                        var sj = cy + mj;
                                        var sk = cz + mk;
                                        if (si < 0 || sj < 0 || sk < 0
                                                        || si >= nX
                                                        || sj >= nY
                                                        || sk >= nZ)
                                                continue;
                                        var index = si * nY * nZ + sj
                                                        * nZ + sk;
                                        if (!(vpBits[index] & ISDONE)) {
                                                vpBits[index] |= ISDONE;
                                                vpAtomID[index] = atom.serial;
                                        } 
                                        
                                        else {
                                                var atom2 = atoms[vpAtomID[index]];
                                                ox = Math.floor(0.5 + scaleFactor
                                                                * (atom2.x + ptranx));
                                                oy = Math.floor(0.5 + scaleFactor
                                                                * (atom2.y + ptrany));
                                                oz = Math.floor(0.5 + scaleFactor
                                                                * (atom2.z + ptranz));
                                                if (mi * mi + mj * mj + mk * mk < ox
                                                                * ox + oy * oy + oz * oz)
                                                        vpAtomID[index] = atom.serial;
                                        }

                                    }// k
                                }// if
                            }// kk
                        }// jj
                    }// ii
                }// if
                nind++;
            }// j
        }// i
    };
    
    this.buildboundary = function() {
        for (var i = 0; i < nX; i++) {
            for (var j = 0; j < nY; j++) {
                for (k = 0; k < nZ; k++) {
                    var index = i * nY * nZ + j * nZ + k;
                    if (vpBits[index] & INOUT) {
                        var flagbound = false;
                        var ii = 0;
                        while (!flagbound && ii < 26) {
                            var ti = i + nb[ii][0], tj = j + nb[ii][1], tk = k
                                     + nb[ii][2];
                            if (ti > -1
                                && ti < nX
                                && tj > -1
                                && tj < nY
                                && tk > -1
                                && tk < nZ                                                            
                                && !(vpBits[ti * nY * nZ + tj
                                                * nZ + tk] & INOUT)) {
                                vpBits[index] |= ISBOUND;
                                flagbound = true;
                            } 
                            else
                                ii++;
                        }
                    }
                }//k
            }//j
        }//i
    };

    // a little class for 3d array, should really generalize this and
    // use throughout...
    var PointGrid = function(length, width, height) {
        // the standard says this is zero initialized
        var data = new Int32Array(length * width * height * 3);

        // set position x,y,z to pt, which has ix,iy,and iz
        this.set = function(x, y, z, pt) {
            var index = ((((x * width) + y) * height) + z) * 3;
            data[index] = pt.ix;
            data[index + 1] = pt.iy;
            data[index + 2] = pt.iz;
        };

        // return point at x,y,z
        this.get = function(x, y, z) {
            var index = ((((x * width) + y) * height) + z) * 3;
            return {
                ix : data[index],
                iy : data[index + 1],
                iz : data[index + 2]
            };
        };
    };
    
    var inarray, outarray;
    this.fastdistancemap = function() {
        var positin = 0, positout = 0, eliminate = 0;
        var certificate;
        var i, j, k;
        var totalsurfacevox = 0;
        var totalinnervox = 0;

        var boundPoint = new PointGrid(nX, nY, nZ);

        for (i = 0; i < nX; i++) {
            for (j = 0; j < nY; j++) {
                for (k = 0; k < nZ; k++) {
                    var index = i * nY * nZ + j * nZ + k;
                    vpBits[index] &= ~ISDONE; // isdone = false
                    if (vpBits[index] & INOUT) {
                        if (vpBits[index] & ISBOUND) {
                            totalsurfacevox++;
                            boundPoint.set(i, j, k, {
                                ix : i,
                                iy : j,
                                iz : k
                            });
                            vpDistance[index] = 0;
                            vpBits[index] |= ISDONE;
                        } 
                        else 
                            totalinnervox++;
                        
                    }
                }
            }
        }

        inarray = new Array();
        outarray = new Array();
        var positin = 0;

        for (i = 0; i < nX; i++) {
            for (j = 0; j < nY; j++) {
                for (k = 0; k < nZ; k++) {
                    var index = i * nY * nZ + j * nZ + k;
                    if (vpBits[index] & ISBOUND) {
                        inarray.push({
                            ix : i,
                            iy : j,
                            iz : k
                        });
                        positin++;
                        vpBits[index] &= ~ISBOUND;
                    }
                }
            }
        }

        do {
            positout = this.fastoneshell(positin, boundPoint);
            positin = 0;
            inarray = [];
            for (i = 0; i < positout; i++) {
                var index = nY * nZ * outarray[i].ix + nZ
                            * outarray[i].iy + outarray[i].iz;
                vpBits[index] &= ~ISBOUND;
                if (vpDistance[index] <= 1.02 * cutRadius) {
                    inarray.push({
                        ix : outarray[i].ix,
                        iy : outarray[i].iy,
                        iz : outarray[i].iz
                    });
                    // inarray[positin].ix=outarray[i].ix;
                    // inarray[positin].iy=outarray[i].iy;
                    // inarray[positin].iz=outarray[i].iz;
                    positin++;
                }
            }
        } while (positin !== 0);

        var cutsf = scaleFactor - 0.5;
        if (cutsf < 0)
            cutsf = 0;
        for (i = 0; i < nX; i++) {
            for (j = 0; j < nY; j++) {
                for (k = 0; k < nZ; k++) {
                    var index = i * nY * nZ + j * nZ + k;
                    vpBits[index] &= ~ISBOUND;
                    // ses solid
                    if (vpBits[index] & INOUT) {
                        if (!(vpBits[index] & ISDONE)
                            || ((vpBits[index] & ISDONE) && vpDistance[index] >= cutRadius
                                 - 0.50 / (0.1 + cutsf))) {
                            vpBits[index] |= ISBOUND;
                        }
                    }
                }
            }
        }
        inarray = [];
        outarray = [];
        boundPoint = null;
    };

    this.fastoneshell = function(number, boundPoint) { // (int* innum,int
        // *allocout,voxel2
        // ***boundPoint, int*
        // outnum, int *elimi)
        var positout = 0;
        var tx, ty, tz;
        var dx, dy, dz;
        var square;
        if (number === 0)
                return 0;
        outarray = [];

        var tnv = {
                ix : -1,
                iy : -1,
                iz : -1
        };
        for ( var i = 0; i < number; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            var bp = boundPoint.get(tx, ty, tz);

            for ( var j = 0; j < 6; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];
                var index = tnv.ix * nY * nZ + nZ * tnv.iy
                            + tnv.iz;
                if (tnv.ix < nX && tnv.ix > -1 && tnv.iy < nY
                    && tnv.iy > -1 && tnv.iz < nZ && tnv.iz > -1
                    && (vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {

                    boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);
                    dx = tnv.ix - bp.ix;
                    dy = tnv.iy - bp.iy;
                    dz = tnv.iz - bp.iz;
                    var square = dx * dx + dy * dy + dz * dz;
                    vpDistance[index] = Math.sqrt(square);
                    vpBits[index] |= ISDONE;
                    vpBits[index] |= ISBOUND;

                    outarray.push({
                        ix : tnv.ix,
                        iy : tnv.iy,
                        iz : tnv.iz
                    });
                    positout++;
                    
                } 
                
                else if (tnv.ix < nX && tnv.ix > -1 && tnv.iy < nY
                         && tnv.iy > -1 && tnv.iz < nZ && tnv.iz > -1
                         && (vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {

                    dx = tnv.ix - bp.ix;
                    dy = tnv.iy - bp.iy;
                    dz = tnv.iz - bp.iz;
                    square = dx * dx + dy * dy + dz * dz;
                    square = Math.sqrt(square);
                    if (square < vpDistance[index]) {
                        boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

                        vpDistance[index] = square;
                        if (!(vpBits[index] & ISBOUND)) {
                            vpBits[index] |= ISBOUND;
                            outarray.push({
                                ix : tnv.ix,
                                iy : tnv.iy,
                                iz : tnv.iz
                            });
                            positout++;
                        }
                    }
                }
            }
        }

        // console.log("part1", positout);

        for (i = 0; i < number; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            var bp = boundPoint.get(tx, ty, tz);

            for (j = 6; j < 18; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];
                var index = tnv.ix * nY * nZ + nZ * tnv.iy
                            + tnv.iz;

                if (tnv.ix < nX && tnv.ix > -1 && tnv.iy < nY
                    && tnv.iy > -1 && tnv.iz < nZ && tnv.iz > -1
                    && (vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {

                    boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);

                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        vpDistance[index] = Math.sqrt(square);
                        vpBits[index] |= ISDONE;
                        vpBits[index] |= ISBOUND;

                        outarray.push({
                                ix : tnv.ix,
                                iy : tnv.iy,
                                iz : tnv.iz
                        });
                        positout++;
                } 
                
                else if (tnv.ix < nX && tnv.ix > -1 && tnv.iy < nY
                                && tnv.iy > -1 && tnv.iz < nZ && tnv.iz > -1
                                && (vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
                    dx = tnv.ix - bp.ix;
                    dy = tnv.iy - bp.iy;
                    dz = tnv.iz - bp.iz;
                    square = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    if (square < vpDistance[index]) {
                        boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
                        vpDistance[index] = square;
                        if (!(vpBits[index] & ISBOUND)) {
                            vpBits[index] |= ISBOUND;
                            outarray.push({
                                ix : tnv.ix,
                                iy : tnv.iy,
                                iz : tnv.iz
                            });
                            positout++;
                        }
                    }
                }
            }
        }

        // console.log("part2", positout);

        for (i = 0; i < number; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            var bp = boundPoint.get(tx, ty, tz);

            for (j = 18; j < 26; j++) {
                tnv.ix = tx + nb[j][0];
                tnv.iy = ty + nb[j][1];
                tnv.iz = tz + nb[j][2];
                var index = tnv.ix * nY * nZ + nZ * tnv.iy
                            + tnv.iz;

                if (tnv.ix < nX && tnv.ix > -1 && tnv.iy < nY
                    && tnv.iy > -1 && tnv.iz < nZ && tnv.iz > -1
                    && (vpBits[index] & INOUT) && !(vpBits[index] & ISDONE)) {
                        
                    boundPoint.set(tnv.ix, tnv.iy, tz + nb[j][2], bp);

                    dx = tnv.ix - bp.ix;
                    dy = tnv.iy - bp.iy;
                    dz = tnv.iz - bp.iz;
                    square = dx * dx + dy * dy + dz * dz;
                    vpDistance[index] = Math.sqrt(square);
                    vpBits[index] |= ISDONE;
                    vpBits[index] |= ISBOUND;

                    outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                    });
                    positout++;
                    
                } 
                
                else if (tnv.ix < nX && tnv.ix > -1 && tnv.iy < nY
                            && tnv.iy > -1 && tnv.iz < nZ && tnv.iz > -1
                            && (vpBits[index] & INOUT) && (vpBits[index] & ISDONE)) {
                        
                    dx = tnv.ix - bp.ix;
                    dy = tnv.iy - bp.iy;
                    dz = tnv.iz - bp.iz;
                    square = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    
                    if (square < vpDistance[index]) {
                        boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

                        vpDistance[index] = square;
                        if (!(vpBits[index] & ISBOUND)) {
                            vpBits[index] |= ISBOUND;
                            outarray.push({
                                ix : tnv.ix,
                                iy : tnv.iy,
                                iz : tnv.iz
                            });
                            positout++;
                        }
                    }
                }
            }
        }

        // console.log("part3", positout);
        return positout;
    };

    //TODO: Is this right? Concerned about consistency of surface type numbering
    this.marchingcubeinit = function(stype) {
        for ( var i = 0, lim = vpBits.length; i < lim; i++) {
            if (stype === 1) {// vdw
                vpBits[i] &= ~ISBOUND;
            } else if (stype === 2) { // ses
                vpBits[i] &= ~ISDONE;
                if (vpBits[i] & ISBOUND)
                    vpBits[i] |= ISDONE;
                vpBits[i] &= ~ISBOUND;
            } else if (stype === 2) {// after vdw
                    if ((vpBits[i] & ISBOUND) && (vpBits[i] & ISDONE))
                            vpBits[i] &= ~ISBOUND;
                    else if ((vpBits[i] & ISBOUND) && !(vpBits[i] & ISDONE))
                            vpBits[i] |= ISDONE;
            } else if (stype === 3) { // sas
                    vpBits[i] &= ~ISBOUND;
            }
        }
    };

    // this code allows me to empirically prune the marching cubes code tables
    // to more efficiently handle discrete data
    var counter = function() {
        var data = Array(256);
        for ( var i = 0; i < 256; i++)
            data[i] = [];

        this.incrementUsed = function(i, j) {
            if (typeof data[i][j] === 'undefined')
                data[i][j] = {
                    used : 0,
                    unused : 0
                };
            data[i][j].used++;
        };

        this.incrementUnused = function(i, j) {
            if (typeof data[i][j] === 'undefined')
                data[i][j] = {
                    used : 0,
                    unused : 0
                };
            data[i][j].unused++;

        };

        var redoTable = function(triTable) {
            var str = "[";
            for ( var i = 0; i < triTable.length; i++) {
                var code = 0;
                var table = triTable[i];
                for ( var j = 0; j < table.length; j++) {
                        code |= (1 << (table[j]));
                }
                str += "0x" + code.toString(16) + ", ";
            }
            str += "]";
            console.log(str);
        };

        this.print = function() {

            var table = MarchingCube.triTable;
            var str;
            var newtable = [];
            for ( var i = 0; i < table.length; i++) {
                var newarr = [];
                for ( var j = 0; j < table[i].length; j += 3) {
                    var k = j / 3;
                    if (typeof data[i][k] === 'undefined' || !data[i][k].unused) {
                        newarr.push(table[i][j]);
                        newarr.push(table[i][j + 1]);
                        newarr.push(table[i][j + 2]);
                    }
                    if (typeof data[i][k] === 'undefined')
                        console.log("undef " + i + "," + k);
                }
                newtable.push(newarr);
            }
            console.log(JSON.stringify(newtable));
            redoTable(newtable);
        };
    };

    this.marchingcube = function(stype) {
        
        this.marchingcubeinit(stype);
        verts = [], faces = [];   
        WebMol.MarchingCube.march(vpBits, verts, faces, {
            smooth : 1,
            nX : nX,
            nY : nY,
            nZ : nZ        
        });	     

        for (i = 0, vlen = verts.length; i < vlen; i++) {
            verts[i].atomid = vpAtomID[verts[i].x * nY * nZ + nZ
                    * verts[i].y + verts[i].z];
        }  

        WebMol.MarchingCube.laplacianSmooth(1, verts, faces);

    };


};
