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


 */

import {MarchingCube} from "./marchingcube"

// a little class for 3d array, should really generalize this and
// use throughout...
export class PointGrid  {

    data: Int32Array;
    width: number;
    height: number;

    constructor(length, width, height) {
        // the standard says this is zero initialized
        this.data = new Int32Array(length * width * height * 3);
        this.width = width;
        this.height = height;
    }

    // set position x,y,z to pt, which has ix,iy,and iz
    set(x:number, y:number, z:number, pt) {
        let index = ((((x * this.width) + y) * this.height) + z) * 3;
        this.data[index] = pt.ix;
        this.data[index + 1] = pt.iy;
        this.data[index + 2] = pt.iz;
    };

    // return point at x,y,z
    get(x:number, y:number, z:number) {
        let index = ((((x * this.width) + y) * this.height) + z) * 3;
        return {
            ix : this.data[index],
            iy : this.data[index + 1],
            iz : this.data[index + 2]
        };
    };
};

/*
 * @type Class
*/
export class ProteinSurface {

    // constants for vpbits bitmasks
    readonly INOUT = 1;
    readonly ISDONE = 2;
    readonly ISBOUND = 4;

    ptranx:number = 0;
    ptrany:number = 0;
    ptranz:number = 0;
    probeRadius:number = 1.4;
    defaultScaleFactor:number = 2;
    scaleFactor:number = this.defaultScaleFactor; // 2 is .5A grid; if this is made user configurable,
                            // also have to adjust offset used to find non-shown
                            // atoms
    pHeight:number = 0;
    pWidth:number = 0;
    pLength:number = 0;
    cutRadius:number = 0;
    vpBits: any = null; // uint8 array of bitmasks
    vpDistance: any = null; // floatarray of _squared_ distances
    vpAtomID: any = null; // intarray
    
    pminx:number = 0;
    pminy:number = 0; 
    pminz:number = 0;
    pmaxx:number = 0;
    pmaxy:number = 0;
    pmaxz:number = 0;

    depty = {};
    widxz = {};
    faces: number[] = [];
    verts = [];

    readonly vdwRadii = {
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
    
    private getVDWIndex(atom:any) {
        if(!atom.elem || typeof(this.vdwRadii[atom.elem]) == "undefined") {
            return "X";
        }
        return atom.elem;
    };
    

    readonly nb = [ new Int32Array([ 1, 0, 0 ]), new Int32Array([ -1, 0, 0 ]), 
               new Int32Array([ 0, 1, 0 ]), new Int32Array([ 0, -1, 0 ]),
               new Int32Array([ 0, 0, 1 ]), 
               new Int32Array([ 0, 0, -1 ]), 
               new Int32Array([ 1, 1, 0 ]), 
               new Int32Array([ 1, -1, 0 ]), 
               new Int32Array([ -1, 1, 0 ]),
               new Int32Array([ -1, -1, 0 ]), 
               new Int32Array([ 1, 0, 1 ]), 
               new Int32Array([ 1, 0, -1 ]), 
               new Int32Array([ -1, 0, 1 ]),
               new Int32Array([ -1, 0, -1 ]), 
               new Int32Array([ 0, 1, 1 ]), 
               new Int32Array([ 0, 1, -1 ]), 
               new Int32Array([ 0, -1, 1 ]),
               new Int32Array([ 0, -1, -1 ]), 
               new Int32Array([ 1, 1, 1 ]), 
               new Int32Array([ 1, 1, -1 ]), 
               new Int32Array([ 1, -1, 1 ]),
               new Int32Array([ -1, 1, 1 ]), 
               new Int32Array([ 1, -1, -1 ]), 
               new Int32Array([ -1, -1, 1 ]), 
               new Int32Array([ -1, 1, -1 ]),
               new Int32Array([ -1, -1, -1 ]) ];


    public getFacesAndVertices(atomlist: any[]) {
        let atomsToShow = {};
        for (let i = 0, il = atomlist.length; i < il; i++)
            atomsToShow[atomlist[i]] = true;
        let vertices = this.verts;
        for (let i = 0, il = vertices.length; i < il; i++) {
            vertices[i].x = vertices[i].x / this.scaleFactor - this.ptranx;
            vertices[i].y = vertices[i].y / this.scaleFactor - this.ptrany;
            vertices[i].z = vertices[i].z / this.scaleFactor - this.ptranz;
        }

        let finalfaces = [];
        for (let i = 0, il = this.faces.length; i < il; i += 3) {
            //let f = faces[i];
            let fa = this.faces[i], fb = this.faces[i+1], fc = this.faces[i+2];
            let a = vertices[fa].atomid, b = vertices[fb].atomid, c = vertices[fc].atomid;

            // must be a unique face for each atom
            let which = a;
            if (b < which)
                which = b;
            if (c < which)
                which = c;
            if (!atomsToShow[which]) {
                continue;
            }

            if (fa !== fb && fb !== fc && fa !== fc){
                finalfaces.push(fa); 
                finalfaces.push(fb); 
                finalfaces.push(fc); 
            }
               
        }

        //try to help the garbage collector
        this.vpBits = null; // uint8 array of bitmasks
        this.vpDistance = null; // floatarray
        this.vpAtomID = null; // intarray
        
        return {
            'vertices' : vertices,
            'faces' : finalfaces
        };
    };


    public initparm (extent: number[][], btype, volume) {
        if(volume > 1000000) //heuristical decrease resolution to avoid large memory consumption
            this.scaleFactor = this.defaultScaleFactor/2;
        
        let margin = (1 / this.scaleFactor) * 5.5; // need margin to avoid
                                                // boundary/round off effects
        this.pminx = extent[0][0]; this.pmaxx = extent[1][0];
        this.pminy = extent[0][1]; this.pmaxy = extent[1][1];
        this.pminz = extent[0][2]; this.pmaxz = extent[1][2];

        if (!btype) {
            this.pminx -= margin;
            this.pminy -= margin;
            this.pminz -= margin;
            this.pmaxx += margin;
            this.pmaxy += margin;
            this.pmaxz += margin;
        } else {
            this.pminx -= this.probeRadius + margin;
            this.pminy -= this.probeRadius + margin;
            this.pminz -= this.probeRadius + margin;
            this.pmaxx += this.probeRadius + margin;
            this.pmaxy += this.probeRadius + margin;
            this.pmaxz += this.probeRadius + margin;
        }

        this.pminx = Math.floor(this.pminx * this.scaleFactor) / this.scaleFactor;
        this.pminy = Math.floor(this.pminy * this.scaleFactor) / this.scaleFactor;
        this.pminz = Math.floor(this.pminz * this.scaleFactor) / this.scaleFactor;
        this.pmaxx = Math.ceil(this.pmaxx * this.scaleFactor) / this.scaleFactor;
        this.pmaxy = Math.ceil(this.pmaxy * this.scaleFactor) / this.scaleFactor;
        this.pmaxz = Math.ceil(this.pmaxz * this.scaleFactor) / this.scaleFactor;

        this.ptranx = -this.pminx;
        this.ptrany = -this.pminy;
        this.ptranz = -this.pminz;

        this.pLength = Math.ceil(this.scaleFactor * (this.pmaxx - this.pminx)) + 1;
        this.pWidth = Math.ceil(this.scaleFactor * (this.pmaxy - this.pminy)) + 1;
        this.pHeight = Math.ceil(this.scaleFactor * (this.pmaxz - this.pminz)) + 1;

        this.boundingatom(btype);
        this.cutRadius = this.probeRadius * this.scaleFactor;

        this.vpBits = new Uint8Array(this.pLength * this.pWidth * this.pHeight);
        this.vpDistance = new Float64Array(this.pLength * this.pWidth * this.pHeight); // float 32
        // doesn't
        // play
        // nicely
        // with
        // native
        // floats
        this.vpAtomID = new Int32Array(this.pLength * this.pWidth * this.pHeight);
    };

    public boundingatom(btype) {
        let tradius = {};

        for ( const i in this.vdwRadii) {
            let r = this.vdwRadii[i];
            if (!btype)
                tradius[i] = r * this.scaleFactor + 0.5;
            else
                tradius[i] = (r + this.probeRadius) * this.scaleFactor + 0.5;

            let sradius = tradius[i] * tradius[i];
            this.widxz[i] = Math.floor(tradius[i]) + 1;
            this.depty[i] = new Int32Array(this.widxz[i] * this.widxz[i]);
            let indx = 0;
            for (let j = 0; j < this.widxz[i]; j++) {
                for (let k = 0; k < this.widxz[i]; k++) {
                    let txz = j * j + k * k;
                    if (txz > sradius)
                        this.depty[i][indx] = -1; // outside
                    else {
                        let tdept = Math.sqrt(sradius - txz);
                        this.depty[i][indx] = Math.floor(tdept);
                    }
                    indx++;
                }
            }
        }
    };

    public fillvoxels(atoms, atomlist) { // (int seqinit,int
        // seqterm,bool
        // atomtype,atom*
        // proseq,bool bcolor)
        for (let i = 0, il = this.vpBits.length; i < il; i++) {
            this.vpBits[i] = 0;
            this.vpDistance[i] = -1.0;
            this.vpAtomID[i] = -1;
        }

        for (let i in atomlist) {
            let atom = atoms[atomlist[i]];
            if (atom === undefined)
                continue;
            this.fillAtom(atom, atoms);
        }

        for (let i = 0, il = this.vpBits.length; i < il; i++)
            if (this.vpBits[i] & this.INOUT)
                this.vpBits[i] |= this.ISDONE;

    };


    public fillAtom(atom, atoms) {

        let cx = Math.floor(0.5 + this.scaleFactor * (atom.x + this.ptranx));
        let cy = Math.floor(0.5 + this.scaleFactor * (atom.y + this.ptrany));
        let cz = Math.floor(0.5 + this.scaleFactor * (atom.z + this.ptranz));

        let at = this.getVDWIndex(atom);
        let nind = 0;
        let pWH = this.pWidth*this.pHeight;
        
        for (let i = 0, n = this.widxz[at]; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (this.depty[at][nind] != -1) {
                    for (let ii = -1; ii < 2; ii++) {
                        for (let jj = -1; jj < 2; jj++) {
                            for (let kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    let mi = ii * i;
                                    let mk = kk * j;
                                    for (let k = 0; k <= this.depty[at][nind]; k++) {
                                        let mj = k * jj;
                                        let si = cx + mi;
                                        let sj = cy + mj;
                                        let sk = cz + mk;
                                        if (si < 0 || sj < 0 || 
                                                sk < 0 ||
                                                si >= this.pLength || 
                                                sj >= this.pWidth || 
                                                sk >= this.pHeight)
                                            continue;
                                        let index = si * pWH + sj * this.pHeight + sk;

                                        if (!(this.vpBits[index] & this.INOUT)) {
                                            this.vpBits[index] |= this.INOUT;
                                            this.vpAtomID[index] = atom.serial;
                                        } else {
                                            let atom2 = atoms[this.vpAtomID[index]];
                                            if(atom2.serial != atom.serial) {
                                                let ox = cx + mi - Math.floor(0.5 + this.scaleFactor *
                                                        (atom2.x + this.ptranx));
                                                let oy = cy + mj - Math.floor(0.5 + this.scaleFactor *
                                                        (atom2.y + this.ptrany));
                                                let oz = cz + mk - Math.floor(0.5 + this.scaleFactor *
                                                        (atom2.z + this.ptranz));
                                                if (mi * mi + mj * mj + mk * mk < ox *
                                                        ox + oy * oy + oz * oz)
                                                    this.vpAtomID[index] = atom.serial;
                                            }
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

    public fillvoxelswaals(atoms, atomlist) {
        for (let i = 0, il = this.vpBits.length; i < il; i++)
            this.vpBits[i] &= ~this.ISDONE; // not isdone

        for (let i in atomlist) {
            let atom = atoms[atomlist[i]];
            if (atom === undefined)
                continue;

            this.fillAtomWaals(atom, atoms);
        }
    };

    public fillAtomWaals(atom, atoms) {
        let nind = 0;
        let cx = Math.floor(0.5 + this.scaleFactor * (atom.x + this.ptranx));
        let cy = Math.floor(0.5 + this.scaleFactor * (atom.y + this.ptrany));
        let cz = Math.floor(0.5 + this.scaleFactor * (atom.z + this.ptranz));

        let at = this.getVDWIndex(atom);
        let pWH = this.pWidth*this.pHeight;
        for (let i = 0, n = this.widxz[at]; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (this.depty[at][nind] != -1) {
                    for (let ii = -1; ii < 2; ii++) {
                        for (let jj = -1; jj < 2; jj++) {
                            for (let kk = -1; kk < 2; kk++) {
                                if (ii !== 0 && jj !== 0 && kk !== 0) {
                                    let mi = ii * i;
                                    let mk = kk * j;
                                    for (let k = 0; k <= this.depty[at][nind]; k++) {
                                        let mj = k * jj;
                                        let si = cx + mi;
                                        let sj = cy + mj;
                                        let sk = cz + mk;
                                        if (si < 0 || sj < 0 || 
                                                sk < 0 || 
                                                si >= this.pLength || 
                                                sj >= this.pWidth || 
                                                sk >= this.pHeight)
                                            continue;
                                        let index = si * pWH + sj * this.pHeight + sk;
                                        if (!(this.vpBits[index] & this.ISDONE)) {
                                            this.vpBits[index] |= this.ISDONE;
                                            this.vpAtomID[index] = atom.serial;
                                        }  else {
                                            let atom2 = atoms[this.vpAtomID[index]];
                                            if(atom2.serial != atom.serial) {
                                                let ox = cx + mi - Math.floor(0.5 + this.scaleFactor *
                                                        (atom2.x + this.ptranx));
                                                let oy = cy + mj - Math.floor(0.5 + this.scaleFactor *
                                                        (atom2.y + this.ptrany));
                                                let oz = cz + mk - Math.floor(0.5 + this.scaleFactor *
                                                        (atom2.z + this.ptranz));
                                                if (mi * mi + mj * mj + mk * mk < ox *
                                                        ox + oy * oy + oz * oz)
                                                    this.vpAtomID[index] = atom.serial;
                                            }
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

    public buildboundary() {
        let pWH = this.pWidth*this.pHeight;
        for (let i = 0; i < this.pLength; i++) {
            for (let j = 0; j < this.pHeight; j++) {
                for (let k = 0; k < this.pWidth; k++) {
                    let index = i * pWH + k * this.pHeight + j;
                    if (this.vpBits[index] & this.INOUT) {
                        let ii = 0;
                        while (ii < 26) {
                            let ti = i + this.nb[ii][0], tj = j + this.nb[ii][2], tk = k +
                                    this.nb[ii][1];
                            if (ti > -1 && 
                                ti < this.pLength && 
                                tk > -1 && 
                                tk < this.pWidth && 
                                tj > -1 && 
                                tj < this.pHeight && 
                                !(this.vpBits[ti * pWH + tk * this.pHeight + tj] & this.INOUT)) {
                                this.vpBits[index] |= this.ISBOUND;
                                break;
                            } else
                                ii++;
                        }
                    }
                }
            }
        }
    };

    public fastdistancemap() {
        let boundPoint = new PointGrid(this.pLength, this.pWidth, this.pHeight);
        let pWH = this.pWidth*this.pHeight;
        let cutRSq = this.cutRadius*this.cutRadius;
        
        let inarray = [];
        let outarray = [];
        
        let index;
        
        for (let i = 0; i < this.pLength; i++) {
            for (let j = 0; j < this.pWidth; j++) {
                for (let k = 0; k < this.pHeight; k++) {
                    index = i * pWH + j * this.pHeight + k;
                    this.vpBits[index] &= ~this.ISDONE; // isdone = false
                    if (this.vpBits[index] & this.INOUT) {
                        if (this.vpBits[index] & this.ISBOUND) {
                            let triple = {
                                ix : i,
                                iy : j,
                                iz : k
                            };
                            boundPoint.set(i, j, k, triple);
                            inarray.push(triple);
                            this.vpDistance[index] = 0;
                            this.vpBits[index] |= this.ISDONE;
                            this.vpBits[index] &= ~this.ISBOUND;
                        } 
                    }
                }
            }
        }

        do {
            outarray = this.fastoneshell(inarray, boundPoint);
            inarray = [];
            for (let i = 0, n = outarray.length; i < n; i++) {
                index = pWH * outarray[i].ix + this.pHeight *
                    outarray[i].iy + outarray[i].iz;
                this.vpBits[index] &= ~this.ISBOUND;
                if (this.vpDistance[index] <= 1.0404 * cutRSq) {
                    inarray.push({
                        ix : outarray[i].ix,
                        iy : outarray[i].iy,
                        iz : outarray[i].iz
                    });
                }
            }
        } while (inarray.length !== 0);

        inarray = [];
        outarray = [];
        boundPoint = null;
        
        let cutsf = this.scaleFactor - 0.5;
        if (cutsf < 0)
            cutsf = 0;
        let cutoff = cutRSq - 0.50 / (0.1 + cutsf);
        for (let i = 0; i < this.pLength; i++) {
            for (let j = 0; j < this.pWidth; j++) {
                for (let k = 0; k < this.pHeight; k++) {
                    index = i * pWH + j * this.pHeight + k;
                    this.vpBits[index] &= ~this.ISBOUND;
                    // ses solid
                    if (this.vpBits[index] & this.INOUT) {
                        if (!(this.vpBits[index] & this.ISDONE) ||
                                ((this.vpBits[index] & this.ISDONE) && this.vpDistance[index] >= cutoff)) {
                            this.vpBits[index] |= this.ISBOUND;
                        }
                    }
                }
            }
        }

    };

    public fastoneshell(inarray, boundPoint) { // (int* innum,int
        // *allocout,voxel2
        // ***boundPoint, int*
        // outnum, int *elimi)
        let tx, ty, tz;
        let dx, dy, dz;
        let square;
        let bp, index;
        let outarray = [];
        if (inarray.length === 0)
            return outarray;

        let tnv = {
            ix : -1,
            iy : -1,
            iz : -1
        };
        let pWH = this.pWidth*this.pHeight;
        for (let i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (let j = 0; j < 6; j++) {
                tnv.ix = tx + this.nb[j][0];
                tnv.iy = ty + this.nb[j][1];
                tnv.iz = tz + this.nb[j][2];
                
                if (tnv.ix < this.pLength && tnv.ix > -1 && tnv.iy < this.pWidth &&
                        tnv.iy > -1 && tnv.iz < this.pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + this.pHeight * tnv.iy + tnv.iz;
                    
                    if ((this.vpBits[index] & this.INOUT) && !(this.vpBits[index] & this.ISDONE)) {
    
                        boundPoint.set(tnv.ix, tnv.iy, tz + this.nb[j][2], bp);
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        this.vpDistance[index] = square;
                        this.vpBits[index] |= this.ISDONE;
                        this.vpBits[index] |= this.ISBOUND;
    
                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((this.vpBits[index] & this.INOUT) && (this.vpBits[index] & this.ISDONE)) {
    
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < this.vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
    
                            this.vpDistance[index] = square;
                            if (!(this.vpBits[index] & this.ISBOUND)) {
                                this.vpBits[index] |= this.ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        for (let i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (let j = 6; j < 18; j++) {
                tnv.ix = tx + this.nb[j][0];
                tnv.iy = ty + this.nb[j][1];
                tnv.iz = tz + this.nb[j][2];

                if(tnv.ix < this.pLength && tnv.ix > -1 && tnv.iy < this.pWidth &&
                        tnv.iy > -1 && tnv.iz < this.pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + this.pHeight * tnv.iy + tnv.iz;
                    
                    if ((this.vpBits[index] & this.INOUT) && !(this.vpBits[index] & this.ISDONE)) {
                        boundPoint.set(tnv.ix, tnv.iy, tz + this.nb[j][2], bp);
    
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        this.vpDistance[index] = square;
                        this.vpBits[index] |= this.ISDONE;
                        this.vpBits[index] |= this.ISBOUND;
    
                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((this.vpBits[index] & this.INOUT) && (this.vpBits[index] & this.ISDONE)) {
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < this.vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);
                            this.vpDistance[index] = square;
                            if (!(this.vpBits[index] & this.ISBOUND)) {
                                this.vpBits[index] |= this.ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        for (let i = 0, n = inarray.length; i < n; i++) {
            tx = inarray[i].ix;
            ty = inarray[i].iy;
            tz = inarray[i].iz;
            bp = boundPoint.get(tx, ty, tz);

            for (let j = 18; j < 26; j++) {
                tnv.ix = tx + this.nb[j][0];
                tnv.iy = ty + this.nb[j][1];
                tnv.iz = tz + this.nb[j][2];

                if (tnv.ix < this.pLength && tnv.ix > -1 && tnv.iy < this.pWidth &&
                        tnv.iy > -1 && tnv.iz < this.pHeight && tnv.iz > -1) {
                    index = tnv.ix * pWH + this.pHeight * tnv.iy + tnv.iz;

                    if ((this.vpBits[index] & this.INOUT) && !(this.vpBits[index] & this.ISDONE)) {
                        boundPoint.set(tnv.ix, tnv.iy, tz + this.nb[j][2], bp);

                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        this.vpDistance[index] = square;
                        this.vpBits[index] |= this.ISDONE;
                        this.vpBits[index] |= this.ISBOUND;

                        outarray.push({
                            ix : tnv.ix,
                            iy : tnv.iy,
                            iz : tnv.iz
                        });
                    } else if ((this.vpBits[index] & this.INOUT)  && (this.vpBits[index] & this.ISDONE)) {
                        dx = tnv.ix - bp.ix;
                        dy = tnv.iy - bp.iy;
                        dz = tnv.iz - bp.iz;
                        square = dx * dx + dy * dy + dz * dz;
                        if (square < this.vpDistance[index]) {
                            boundPoint.set(tnv.ix, tnv.iy, tnv.iz, bp);

                            this.vpDistance[index] = square;
                            if (!(this.vpBits[index] & this.ISBOUND)) {
                                this.vpBits[index] |= this.ISBOUND;
                                outarray.push({
                                    ix : tnv.ix,
                                    iy : tnv.iy,
                                    iz : tnv.iz
                                });
                            }
                        }
                    }
                }
            }
        }

        return outarray;
    };

    public marchingcubeinit(stype) {
        for ( let i = 0, lim = this.vpBits.length; i < lim; i++) {
            if (stype == 1) {// vdw
                this.vpBits[i] &= ~this.ISBOUND;
            } else if (stype == 4) { // ses
                this.vpBits[i] &= ~this.ISDONE;
                if (this.vpBits[i] & this.ISBOUND)
                    this.vpBits[i] |= this.ISDONE;
                this.vpBits[i] &= ~this.ISBOUND;
            } else if (stype == 2) {// after vdw
                if ((this.vpBits[i] & this.ISBOUND) && (this.vpBits[i] & this.ISDONE))
                    this.vpBits[i] &= ~this.ISBOUND;
                else if ((this.vpBits[i] & this.ISBOUND) && !(this.vpBits[i] & this.ISDONE))
                    this.vpBits[i] |= this.ISDONE;
            } else if (stype == 3) { // sas
                this.vpBits[i] &= ~this.ISBOUND;
            }
        }
    };
    
    public marchingcube(stype) {
        this.marchingcubeinit(stype);
        this.verts = []; this.faces = [];   
        MarchingCube.march(this.vpBits, this.verts, this.faces, {
            smooth : 1,
            nX : this.pLength,
            nY : this.pWidth,
            nZ : this.pHeight        
        });      

        let pWH = this.pWidth*this.pHeight;
        for (let i = 0, vlen = this.verts.length; i < vlen; i++) {
            this.verts[i].atomid = this.vpAtomID[this.verts[i].x * pWH + this.pHeight *
                    this.verts[i].y + this.verts[i].z];
        }  

        MarchingCube.laplacianSmooth(1, this.verts, this.faces);

    };


};
