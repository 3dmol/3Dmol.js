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
import {Vector3} from "./WebGL/math"


/**
 * Surface types
 * @readonly
 * @enum {number}
 * @property VDW - van der Waals surface
 * @property MS - Molecular surface
 * @property SAS - Solvent accessible surface
 * @property SES - Solvent exposed surface
 */
export enum SurfaceType {
    VDW = 1,
    MS = 2,
    SAS = 3,
    SES = 4
};

/** 
 * Render surface synchronously if true
 * @param {boolean} [$3Dmol.SyncSurface=false]
 * @type {boolean} */
export var syncSurface = false;

export function setSyncSurface(val:boolean) {
    syncSurface = val;
}

// Internet Explorer refuses to allow webworkers in data blobs.  I can find
// no way of checking for this feature directly, so must do a browser check
if (window.navigator.userAgent.indexOf('MSIE ') >= 0 ||
    window.navigator.userAgent.indexOf('Trident/') >= 0) {
    syncSurface = true; // can't use webworkers
}


//Encapsulate marching cube algorithm for isosurface generation
// (currently used by protein surface rendering and generic volumetric data reading)
export class MarchingCubeInitializer {
    
    //Marching cube algorithm - assume data has been pre-treated so isovalue is 0 
    // (i.e. select points greater than 0)
    //origin -  vector of origin of volumetric data (default is (0,0,0))
    // nX, nY, nZ - specifies number of voxels in each dimension
    // scale - cube diagonal unit vector scale (3Dmol vector) (specifying distance between data points); diagonal of cube
    // - default is 1 - assumes unit cube (1,1,1) diag)
    // fulltable - if true, use full marching cubes and tritables - else use trimmed table (e.g. surf render)
    // voxel - if true, draws with a blocky voxel style (default false)
    // verts, faces - vertex and face arrays to fill up
    
    //to match with protein surface...
    ISDONE: number = 2;

    constructor() {

    }

    march(data, verts, faces, spec) {

        let fulltable = !!(spec.fulltable);
        let origin = (spec.hasOwnProperty('origin') && spec.origin.hasOwnProperty('x')) ? spec.origin : {x:0, y:0, z:0};
        let voxel = !!(spec.voxel);
        let transform = spec.matrix; //if this is set, it overrides origin and unitCube
        
        let nX = spec.nX || 0;
        let nY = spec.nY || 0;
        let nZ = spec.nZ || 0;
        
        let scale = spec.scale || 1.0;
        let unitCube = null;
        if(spec.unitCube) {
            unitCube = spec.unitCube;
        } else {
            unitCube = {x:scale,y:scale,z:scale};
        }
        
        //keep track of calculated vertices to avoid repeats
        let vertnums = new Int32Array(nX*nY*nZ);
        
        let i, il;
        
        for (i = 0, il = vertnums.length; i < il; ++i)
            vertnums[i] = -1;

        // create (or retrieve) a vertex at the appropriate point for
        // the edge (p1,p2)
        
        let getVertex = function(i, j, k, code, p1, p2) {
            let pt = {x:0,y:0,z:0};
            let val1 = !!(code & (1 << p1));
            let val2 = !!(code & (1 << p2));
             
            // p1 if they are the same or if !val1
            let p = p1;
            if (!val1 && val2)
                p = p2;
            
            // adjust i,j,k by p
            if (p & 1)
                k++;
            if (p & 2)
                j++;
            if (p & 4)
                i++;

            if(transform) {
                let vpt = new Vector3(i,j,k);
                vpt = vpt.applyMatrix4(transform);
                pt = {x: vpt.x, y: vpt.y, z: vpt.z}; //remove vector gunk
            } else {
                pt.x = origin.x+unitCube.x*i;
                pt.y = origin.y+unitCube.y*j;
                pt.z = origin.z+unitCube.z*k;
            }
            
            let index = ((nY * i) + j) * nZ + k;
            //Have to add option to do voxels
            if (!voxel) {
            
                if (vertnums[index] < 0) // not created yet
                {
                    vertnums[index] = verts.length;
                    verts.push( pt );
                }
                return vertnums[index];
            
            }
            
            else {
                verts.push(pt);
                return verts.length - 1;
            }
            
        };
            
        let intersects = new Int32Array(12);
        
        let etable = (fulltable) ? this.edgeTable2 : this.edgeTable;
        let tritable = (fulltable) ? this.triTable2 : this.triTable;
                
        //Run marching cubes algorithm
        for (i = 0; i < nX-1; ++i) {
            
            for (let j = 0; j < nY-1; ++j){
                
                for (let k = 0; k < nZ-1; ++k){
                    
                    let code = 0;
                    
                    for (let p = 0; p < 8; ++p) {
                        let index = ((nY * (i + ((p & 4) >> 2))) + j + ((p & 2) >> 1)) *
                                        nZ + k + (p & 1);

                        //TODO: Need to fix vpBits in protein surface for this to work
                        let val:any = !!(data[index] & this.ISDONE);
                        //let val = !!(data[index] > 0);   
                        
                        code |= val << p;                        
                    }
                    
                    if (code === 0 || code === 255)
                        continue;
                    
                    let ecode = etable[code];
                    
                    if (ecode === 0)
                        continue;
                        
                    let ttable = tritable[code];                        
                    
                    if (ecode & 1)
                        intersects[0] = getVertex(i, j, k, code, 0, 1);
                    if (ecode & 2)
                        intersects[1] = getVertex(i, j, k, code, 1, 3);
                    if (ecode & 4)
                        intersects[2] = getVertex(i, j, k, code, 3, 2);
                    if (ecode & 8)
                        intersects[3] = getVertex(i, j, k, code, 2, 0);
                    if (ecode & 16)
                        intersects[4] = getVertex(i, j, k, code, 4, 5);
                    if (ecode & 32)
                        intersects[5] = getVertex(i, j, k, code, 5, 7);
                    if (ecode & 64)
                        intersects[6] = getVertex(i, j, k, code, 7, 6);
                    if (ecode & 128)
                        intersects[7] = getVertex(i, j, k, code, 6, 4);
                    if (ecode & 256)
                        intersects[8] = getVertex(i, j, k, code, 0, 4);
                    if (ecode & 512)
                        intersects[9] = getVertex(i, j, k, code, 1, 5);
                    if (ecode & 1024)
                        intersects[10] = getVertex(i, j, k, code, 3, 7);
                    if (ecode & 2048)
                        intersects[11] = getVertex(i, j, k, code, 2, 6);       

                    for (let t = 0; t < ttable.length; t += 3) {
                        
                        let a = intersects[ttable[t]],
                            b = intersects[ttable[t+1]],
                            c = intersects[ttable[t+2]];         
                                           
                        if (voxel && t >= 3) {
                            verts.push(verts[a]); a = verts.length - 1;
                            verts.push(verts[b]); b = verts.length - 1;
                            verts.push(verts[c]); c = verts.length - 1;
                        }

                        
                        faces.push(a); faces.push(b); faces.push(c);                               
                    }              
                    
                }
                
            }
            
        }
             
        
    };

    laplacianSmooth (numiter, verts, faces) {
            let tps = new Array(verts.length);
            let i, il, j, jl, k;
            for (i = 0, il = verts.length; i < il; i++)
                    tps[i] = {
                        x : 0,
                        y : 0,
                        z : 0
                    };
            let vertdeg = new Array(20);
            let flagvert;
            for (i = 0; i < 20; i++)
                    vertdeg[i] = new Array(verts.length);
            for (i = 0, il = verts.length; i < il; i++)
                    vertdeg[0][i] = 0;
            for (i = 0, il = faces.length / 3; i < il; i++) {
                let aoffset = i*3, boffset = i*3 + 1, coffset = i*3 + 2;
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[aoffset]]; j < jl; j++) {
                    if (faces[boffset] == vertdeg[j + 1][faces[aoffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[aoffset]]++;
                    vertdeg[vertdeg[0][faces[aoffset]]][faces[aoffset]] = faces[boffset];
                }
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[aoffset]]; j < jl; j++) {
                    if (faces[coffset] == vertdeg[j + 1][faces[aoffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[aoffset]]++;
                    vertdeg[vertdeg[0][faces[aoffset]]][faces[aoffset]] = faces[coffset];
                }
                // b
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[boffset]]; j < jl; j++) {
                    if (faces[aoffset] == vertdeg[j + 1][faces[boffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[boffset]]++;
                    vertdeg[vertdeg[0][faces[boffset]]][faces[boffset]] = faces[aoffset];
                }
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[boffset]]; j < jl; j++) {
                    if (faces[coffset] == vertdeg[j + 1][faces[boffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[boffset]]++;
                    vertdeg[vertdeg[0][faces[boffset]]][faces[boffset]] = faces[coffset];
                }
                // c
                flagvert = true;
                for (j = 0; j < vertdeg[0][faces[coffset]]; j++) {
                    if (faces[aoffset] == vertdeg[j + 1][faces[coffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[coffset]]++;
                    vertdeg[vertdeg[0][faces[coffset]]][faces[coffset]] = faces[aoffset];
                }
                flagvert = true;
                for (j = 0, jl = vertdeg[0][faces[coffset]]; j < jl; j++) {
                    if (faces[boffset] == vertdeg[j + 1][faces[coffset]]) {
                        flagvert = false;
                        break;
                    }
                }
                if (flagvert) {
                    vertdeg[0][faces[coffset]]++;
                    vertdeg[vertdeg[0][faces[coffset]]][faces[coffset]] = faces[boffset];
                }
            }

            let wt = 1.00;
            let wt2 = 0.50;

            for (k = 0; k < numiter; k++) {
                    for (i = 0, il = verts.length; i < il; i++) {
                            if (vertdeg[0][i] < 3) {
                                    tps[i].x = verts[i].x;
                                    tps[i].y = verts[i].y;
                                    tps[i].z = verts[i].z;
                            } else if (vertdeg[0][i] == 3 || vertdeg[0][i] == 4) {
                                    tps[i].x = 0;
                                    tps[i].y = 0;
                                    tps[i].z = 0;
                                    for (j = 0, jl = vertdeg[0][i]; j < jl; j++) {
                                            tps[i].x += verts[vertdeg[j + 1][i]].x;
                                            tps[i].y += verts[vertdeg[j + 1][i]].y;
                                            tps[i].z += verts[vertdeg[j + 1][i]].z;
                                    }
                                    tps[i].x += wt2 * verts[i].x;
                                    tps[i].y += wt2 * verts[i].y;
                                    tps[i].z += wt2 * verts[i].z;
                                    tps[i].x /= wt2 + vertdeg[0][i];
                                    tps[i].y /= wt2 + vertdeg[0][i];
                                    tps[i].z /= wt2 + vertdeg[0][i];
                            } else {
                                    tps[i].x = 0;
                                    tps[i].y = 0;
                                    tps[i].z = 0;
                                    for (j = 0, jl = vertdeg[0][i]; j < jl; j++) {
                                            tps[i].x += verts[vertdeg[j + 1][i]].x;
                                            tps[i].y += verts[vertdeg[j + 1][i]].y;
                                            tps[i].z += verts[vertdeg[j + 1][i]].z;
                                    }
                                    tps[i].x += wt * verts[i].x;
                                    tps[i].y += wt * verts[i].y;
                                    tps[i].z += wt * verts[i].z;
                                    tps[i].x /= wt + vertdeg[0][i];
                                    tps[i].y /= wt + vertdeg[0][i];
                                    tps[i].z /= wt + vertdeg[0][i];
                            }
                    }
                    for (i = 0, il = verts.length; i < il; i++) {
                            verts[i].x = tps[i].x;
                            verts[i].y = tps[i].y;
                            verts[i].z = tps[i].z;
                    }
                    /*
                     * computenorm(); for (let i = 0; i < vertnumber; i++) { if
                     * (verts[i].inout) ssign = 1; else ssign = -1; verts[i].x += ssign *
                     * outwt * verts[i].pn.x; verts[i].y += ssign * outwt *
                     * verts[i].pn.y; verts[i].z += ssign * outwt * verts[i].pn.z; }
                     */
            }
    };


    /*
     * These tables are based off those by Paul Bourke and Geoffrey Heller:
     * http://paulbourke.net/geometry/polygonise/
     * http://paulbourke.net/geometry/polygonise/table2.txt
     * 
     * However, they have been substantially modified to reflect a more 
     * sensible corner numbering scheme and the discrete nature of our voxel data
     * (resulting in fewer faces).
     */
    edgeTable: Uint32Array = new Uint32Array([ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
            0xb00, 0x0, 0x0, 0x0, 0x700, 0x0, 0xd00, 0xe00, 0xf00, 0x0, 0x0, 0x0,
            0x8a, 0x0, 0x15, 0x0, 0x86, 0x0, 0x0, 0x0, 0x28c, 0x0, 0x813, 0xf19,
            0xe10, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x0, 0x0, 0x126, 0x0, 0x0, 0x15, 0x1c,
            0x0, 0xf23, 0x419, 0xd20, 0x0, 0xa8, 0xa2, 0xaa, 0x0, 0x285, 0x9ab,
            0x8a2, 0x0, 0x2af, 0x125, 0xac, 0xfaa, 0xea3, 0xda9, 0xca0, 0x0, 0x0,
            0x0, 0x0, 0x0, 0x45, 0x0, 0x384, 0x0, 0x0, 0x0, 0x700, 0x8a, 0x83,
            0x648, 0x780, 0x0, 0x51, 0x0, 0x81a, 0x54, 0x55, 0x54, 0x56, 0x0, 0x51,
            0x0, 0xe5c, 0x14a, 0x451, 0x759, 0x650, 0x0, 0x0, 0x0, 0x2a, 0x0, 0x45,
            0x0, 0x1f6, 0x0, 0x0, 0x15, 0xdfc, 0x8a, 0x7f3, 0x4f9, 0x5f0, 0xb00,
            0x68, 0x921, 0x6a, 0x348, 0x245, 0x16f, 0x66, 0xb00, 0xe6f, 0xd65,
            0xc6c, 0x76a, 0x663, 0x569, 0x460, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
            0xf46, 0x0, 0x0, 0x45, 0x24c, 0x2a, 0x823, 0x29, 0xb40, 0x0, 0x0, 0x0,
            0x6ba, 0x0, 0x8f5, 0xfff, 0xef6, 0x0, 0xff, 0x2f5, 0x2fc, 0x9ea, 0x8f3,
            0xbf9, 0xaf0, 0x0, 0x0, 0x51, 0x152, 0x0, 0xf55, 0x45f, 0xd56, 0x54,
            0x357, 0x55, 0x154, 0x852, 0xb53, 0x59, 0x950, 0x700, 0x2c8, 0xc2,
            0x48a, 0xfc4, 0xec5, 0xdcf, 0xcc6, 0x2c4, 0x2cf, 0xc5, 0xcc, 0xbca,
            0xac3, 0x9c9, 0x8c0, 0x0, 0x0, 0x0, 0x0, 0xa8, 0x1a4, 0xa8, 0x7a6,
            0xa2, 0xa2, 0x2a4, 0xbac, 0xaa, 0xa3, 0x2a8, 0x3a0, 0xd00, 0xc18,
            0xd00, 0xe3a, 0x34, 0x35, 0x73f, 0x636, 0x924, 0x83f, 0xb35, 0xa3c,
            0x12a, 0x33, 0x339, 0x230, 0xe00, 0xe00, 0xc12, 0xd9a, 0x684, 0x795,
            0x49f, 0x596, 0x92, 0xb9f, 0x815, 0x99c, 0x9a, 0x393, 0x99, 0x190,
            0xf00, 0xe08, 0xd01, 0xc0a, 0x704, 0x605, 0x50f, 0x406, 0xb02, 0xa0f,
            0x905, 0x80c, 0x30a, 0x203, 0x109, 0x0 ]);
    
    
    triTable =  [ [], [], [], [], [], [], [], [ 11, 9, 8 ], [], [], [],
            [ 8, 10, 9 ], [], [ 10, 8, 11 ], [ 9, 11, 10 ],
            [ 8, 10, 9, 8, 11, 10 ], [], [], [], [ 1, 7, 3 ], [], [ 4, 2, 0 ], [],
            [ 2, 1, 7 ], [], [], [], [ 2, 7, 3, 2, 9, 7 ], [],
            [ 1, 4, 11, 1, 0, 4 ], [ 3, 8, 0, 11, 9, 4, 11, 10, 9 ],
            [ 4, 11, 9, 11, 10, 9 ], [], [], [], [ 5, 3, 1 ], [], [], [],
            [ 2, 5, 8, 2, 1, 5 ], [], [], [ 2, 4, 0 ], [ 3, 2, 4 ], [],
            [ 0, 9, 1, 8, 10, 5, 8, 11, 10 ], [ 3, 4, 0, 3, 10, 4 ],
            [ 5, 8, 10, 8, 11, 10 ], [], [ 3, 5, 7 ], [ 7, 1, 5 ],
            [ 1, 7, 3, 1, 5, 7 ], [], [ 9, 2, 0, 9, 7, 2 ],
            [ 0, 3, 8, 1, 7, 11, 1, 5, 7 ], [ 11, 1, 7, 1, 5, 7 ], [],
            [ 9, 1, 0, 5, 3, 2, 5, 7, 3 ], [ 8, 2, 5, 8, 0, 2 ],
            [ 2, 5, 3, 5, 7, 3 ], [ 3, 9, 1, 3, 8, 9, 7, 11, 10, 7, 10, 5 ],
            [ 9, 1, 0, 10, 7, 11, 10, 5, 7 ], [ 3, 8, 0, 7, 10, 5, 7, 11, 10 ],
            [ 11, 5, 7, 11, 10, 5 ], [], [], [], [], [], [ 0, 6, 2 ], [],
            [ 7, 2, 9, 7, 9, 8 ], [], [], [], [ 8, 10, 9 ], [ 7, 1, 3 ],
            [ 7, 1, 0 ], [ 6, 9, 3, 6, 10, 9 ], [ 7, 10, 8, 10, 9, 8 ], [],
            [ 6, 0, 4 ], [], [ 11, 1, 4, 11, 3, 1 ], [ 2, 4, 6 ],
            [ 2, 0, 4, 2, 4, 6 ], [ 2, 4, 6 ], [ 1, 4, 2, 4, 6, 2 ], [],
            [ 6, 0, 4 ], [], [ 2, 11, 3, 6, 9, 4, 6, 10, 9 ], [ 8, 6, 1, 8, 1, 3 ],
            [ 10, 0, 6, 0, 4, 6 ], [ 8, 0, 3, 9, 6, 10, 9, 4, 6 ],
            [ 10, 4, 6, 10, 9, 4 ], [], [], [], [ 5, 3, 1 ], [], [ 0, 6, 2 ], [],
            [ 7, 4, 8, 5, 2, 1, 5, 6, 2 ], [], [], [ 2, 4, 0 ],
            [ 7, 4, 8, 2, 11, 3, 10, 5, 6 ], [ 7, 1, 3 ],
            [ 5, 6, 10, 0, 9, 1, 8, 7, 4 ], [ 5, 6, 10, 7, 0, 3, 7, 4, 0 ],
            [ 10, 5, 6, 4, 8, 7 ], [ 9, 11, 8 ], [ 3, 5, 6 ],
            [ 0, 5, 11, 0, 11, 8 ], [ 6, 3, 5, 3, 1, 5 ], [ 3, 9, 6, 3, 8, 9 ],
            [ 9, 6, 0, 6, 2, 0 ], [ 0, 3, 8, 2, 5, 6, 2, 1, 5 ],
            [ 1, 6, 2, 1, 5, 6 ], [ 9, 11, 8 ], [ 1, 0, 9, 6, 10, 5, 11, 3, 2 ],
            [ 6, 10, 5, 2, 8, 0, 2, 11, 8 ], [ 3, 2, 11, 10, 5, 6 ],
            [ 10, 5, 6, 9, 3, 8, 9, 1, 3 ], [ 0, 9, 1, 5, 6, 10 ],
            [ 8, 0, 3, 10, 5, 6 ], [ 10, 5, 6 ], [], [], [], [], [], [], [],
            [ 1, 10, 2, 9, 11, 6, 9, 8, 11 ], [], [], [ 6, 0, 2 ],
            [ 3, 6, 9, 3, 2, 6 ], [ 3, 5, 1 ], [ 0, 5, 1, 0, 11, 5 ], [ 0, 3, 5 ],
            [ 6, 9, 11, 9, 8, 11 ], [], [], [], [ 4, 5, 9, 7, 1, 10, 7, 3, 1 ], [],
            [ 11, 6, 7, 2, 4, 5, 2, 0, 4 ],
            [ 11, 6, 7, 8, 0, 3, 1, 10, 2, 9, 4, 5 ],
            [ 6, 7, 11, 1, 10, 2, 9, 4, 5 ], [],
            [ 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2 ], [ 9, 4, 5, 0, 6, 7, 0, 2, 6 ],
            [ 4, 5, 9, 6, 3, 2, 6, 7, 3 ], [ 6, 7, 11, 5, 3, 8, 5, 1, 3 ],
            [ 6, 7, 11, 4, 1, 0, 4, 5, 1 ], [ 4, 5, 9, 3, 8, 0, 11, 6, 7 ],
            [ 9, 4, 5, 7, 11, 6 ], [], [], [ 0, 6, 4 ], [ 8, 6, 4, 8, 1, 6 ], [],
            [ 0, 10, 2, 0, 9, 10, 4, 8, 11, 4, 11, 6 ],
            [ 10, 2, 1, 6, 0, 3, 6, 4, 0 ], [ 10, 2, 1, 11, 4, 8, 11, 6, 4 ],
            [ 4, 2, 6 ], [ 1, 0, 9, 2, 4, 8, 2, 6, 4 ], [ 2, 4, 0, 2, 6, 4 ],
            [ 8, 2, 4, 2, 6, 4 ], [ 11, 4, 1, 11, 6, 4 ],
            [ 0, 9, 1, 4, 11, 6, 4, 8, 11 ], [ 3, 6, 0, 6, 4, 0 ],
            [ 8, 6, 4, 8, 11, 6 ], [ 10, 8, 9 ], [ 6, 3, 9, 6, 7, 3 ], [ 6, 7, 1 ],
            [ 10, 7, 1, 7, 3, 1 ], [ 7, 11, 6, 8, 10, 2, 8, 9, 10 ],
            [ 11, 6, 7, 10, 0, 9, 10, 2, 0 ], [ 2, 1, 10, 7, 11, 6, 8, 0, 3 ],
            [ 1, 10, 2, 6, 7, 11 ], [ 7, 2, 6, 7, 9, 2 ],
            [ 1, 0, 9, 3, 6, 7, 3, 2, 6 ], [ 7, 0, 6, 0, 2, 6 ],
            [ 2, 7, 3, 2, 6, 7 ], [ 7, 11, 6, 3, 9, 1, 3, 8, 9 ],
            [ 9, 1, 0, 11, 6, 7 ], [ 0, 3, 8, 11, 6, 7 ], [ 11, 6, 7 ], [], [], [],
            [], [ 5, 3, 7 ], [ 8, 5, 2, 8, 7, 5 ], [ 5, 3, 7 ],
            [ 1, 10, 2, 5, 8, 7, 5, 9, 8 ], [ 1, 7, 5 ], [ 1, 7, 5 ],
            [ 9, 2, 7, 9, 7, 5 ], [ 11, 3, 2, 8, 5, 9, 8, 7, 5 ],
            [ 1, 3, 7, 1, 7, 5 ], [ 0, 7, 1, 7, 5, 1 ], [ 9, 3, 5, 3, 7, 5 ],
            [ 9, 7, 5, 9, 8, 7 ], [ 8, 10, 11 ], [ 3, 4, 10, 3, 10, 11 ],
            [ 8, 10, 11 ], [ 5, 9, 4, 1, 11, 3, 1, 10, 11 ], [ 2, 4, 5 ],
            [ 5, 2, 4, 2, 0, 4 ], [ 0, 3, 8, 5, 9, 4, 10, 2, 1 ],
            [ 2, 1, 10, 9, 4, 5 ], [ 2, 8, 5, 2, 11, 8 ],
            [ 3, 2, 11, 1, 4, 5, 1, 0, 4 ], [ 9, 4, 5, 8, 2, 11, 8, 0, 2 ],
            [ 11, 3, 2, 9, 4, 5 ], [ 8, 5, 3, 5, 1, 3 ], [ 5, 0, 4, 5, 1, 0 ],
            [ 3, 8, 0, 4, 5, 9 ], [ 9, 4, 5 ], [ 11, 9, 10 ], [ 11, 9, 10 ],
            [ 1, 11, 4, 1, 10, 11 ], [ 8, 7, 4, 11, 1, 10, 11, 3, 1 ],
            [ 2, 7, 9, 2, 9, 10 ], [ 4, 8, 7, 0, 10, 2, 0, 9, 10 ],
            [ 2, 1, 10, 0, 7, 4, 0, 3, 7 ], [ 10, 2, 1, 8, 7, 4 ], [ 1, 7, 4 ],
            [ 3, 2, 11, 4, 8, 7, 9, 1, 0 ], [ 11, 4, 2, 4, 0, 2 ],
            [ 2, 11, 3, 7, 4, 8 ], [ 4, 1, 7, 1, 3, 7 ], [ 1, 0, 9, 8, 7, 4 ],
            [ 3, 4, 0, 3, 7, 4 ], [ 8, 7, 4 ], [ 8, 9, 10, 8, 10, 11 ],
            [ 3, 9, 11, 9, 10, 11 ], [ 0, 10, 8, 10, 11, 8 ],
            [ 10, 3, 1, 10, 11, 3 ], [ 2, 8, 10, 8, 9, 10 ], [ 9, 2, 0, 9, 10, 2 ],
            [ 8, 0, 3, 1, 10, 2 ], [ 10, 2, 1 ], [ 1, 11, 9, 11, 8, 9 ],
            [ 11, 3, 2, 0, 9, 1 ], [ 11, 0, 2, 11, 8, 0 ], [ 11, 3, 2 ],
            [ 8, 1, 3, 8, 9, 1 ], [ 9, 1, 0 ], [ 8, 0, 3 ], [] ];
     
    edgeTable2 = [ 0x0, 0x109, 0x203, 0x30a, 0x80c, 0x905, 0xa0f,
            0xb06, 0x406, 0x50f, 0x605, 0x70c, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190,
            0x99, 0x393, 0x29a, 0x99c, 0x895, 0xb9f, 0xa96, 0x596, 0x49f, 0x795,
            0x69c, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0xa3c,
            0xb35, 0x83f, 0x936, 0x636, 0x73f, 0x435, 0x53c, 0xe3a, 0xf33, 0xc39,
            0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa, 0xbac, 0xaa5, 0x9af, 0x8a6, 0x7a6,
            0x6af, 0x5a5, 0x4ac, 0xfaa, 0xea3, 0xda9, 0xca0, 0x8c0, 0x9c9, 0xac3,
            0xbca, 0xcc, 0x1c5, 0x2cf, 0x3c6, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x4ca,
            0x5c3, 0x6c9, 0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0x15c, 0x55, 0x35f,
            0x256, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x55a, 0x453, 0x759, 0x650, 0xaf0,
            0xbf9, 0x8f3, 0x9fa, 0x2fc, 0x3f5, 0xff, 0x1f6, 0xef6, 0xfff, 0xcf5,
            0xdfc, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0x36c,
            0x265, 0x16f, 0x66, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x76a, 0x663, 0x569,
            0x460, 0x460, 0x569, 0x663, 0x76a, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x66,
            0x16f, 0x265, 0x36c, 0x86a, 0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3,
            0x6fa, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x1f6, 0xff, 0x3f5, 0x2fc, 0x9fa,
            0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0xe5c, 0xf55, 0xc5f,
            0xd56, 0x256, 0x35f, 0x55, 0x15c, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0,
            0x6c9, 0x5c3, 0x4ca, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0x3c6, 0x2cf, 0x1c5,
            0xcc, 0xbca, 0xac3, 0x9c9, 0x8c0, 0xca0, 0xda9, 0xea3, 0xfaa, 0x4ac,
            0x5a5, 0x6af, 0x7a6, 0x8a6, 0x9af, 0xaa5, 0xbac, 0xaa, 0x1a3, 0x2a9,
            0x3a0, 0xd30, 0xc39, 0xf33, 0xe3a, 0x53c, 0x435, 0x73f, 0x636, 0x936,
            0x83f, 0xb35, 0xa3c, 0x13a, 0x33, 0x339, 0x230, 0xe90, 0xf99, 0xc93,
            0xd9a, 0x69c, 0x795, 0x49f, 0x596, 0xa96, 0xb9f, 0x895, 0x99c, 0x29a,
            0x393, 0x99, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0x70c, 0x605, 0x50f,
            0x406, 0xb06, 0xa0f, 0x905, 0x80c, 0x30a, 0x203, 0x109, 0x0 ];
     
    triTable2 = [ [], [ 8, 3, 0 ], [ 9, 0, 1 ], [ 8, 3, 1, 8, 1, 9 ],
            [ 11, 2, 3 ], [ 11, 2, 0, 11, 0, 8 ], [ 11, 2, 3, 0, 1, 9 ],
            [ 2, 1, 11, 1, 9, 11, 11, 9, 8 ], [ 10, 1, 2 ], [ 8, 3, 0, 1, 2, 10 ],
            [ 9, 0, 2, 9, 2, 10 ], [ 3, 2, 8, 2, 10, 8, 8, 10, 9 ],
            [ 10, 1, 3, 10, 3, 11 ], [ 1, 0, 10, 0, 8, 10, 10, 8, 11 ],
            [ 0, 3, 9, 3, 11, 9, 9, 11, 10 ], [ 8, 10, 9, 8, 11, 10 ], [ 8, 4, 7 ],
            [ 3, 0, 4, 3, 4, 7 ], [ 1, 9, 0, 8, 4, 7 ],
            [ 9, 4, 1, 4, 7, 1, 1, 7, 3 ], [ 2, 3, 11, 7, 8, 4 ],
            [ 7, 11, 4, 11, 2, 4, 4, 2, 0 ], [ 3, 11, 2, 4, 7, 8, 9, 0, 1 ],
            [ 2, 7, 11, 2, 1, 7, 1, 4, 7, 1, 9, 4 ], [ 10, 1, 2, 8, 4, 7 ],
            [ 2, 10, 1, 0, 4, 7, 0, 7, 3 ], [ 4, 7, 8, 0, 2, 10, 0, 10, 9 ],
            [ 2, 7, 3, 2, 9, 7, 7, 9, 4, 2, 10, 9 ],
            [ 8, 4, 7, 11, 10, 1, 11, 1, 3 ],
            [ 11, 4, 7, 1, 4, 11, 1, 11, 10, 1, 0, 4 ],
            [ 3, 8, 0, 7, 11, 4, 11, 9, 4, 11, 10, 9 ],
            [ 7, 11, 4, 4, 11, 9, 11, 10, 9 ], [ 9, 5, 4 ], [ 3, 0, 8, 4, 9, 5 ],
            [ 5, 4, 0, 5, 0, 1 ], [ 4, 8, 5, 8, 3, 5, 5, 3, 1 ],
            [ 11, 2, 3, 9, 5, 4 ], [ 9, 5, 4, 8, 11, 2, 8, 2, 0 ],
            [ 3, 11, 2, 1, 5, 4, 1, 4, 0 ],
            [ 8, 5, 4, 2, 5, 8, 2, 8, 11, 2, 1, 5 ], [ 2, 10, 1, 9, 5, 4 ],
            [ 0, 8, 3, 5, 4, 9, 10, 1, 2 ], [ 10, 5, 2, 5, 4, 2, 2, 4, 0 ],
            [ 3, 4, 8, 3, 2, 4, 2, 5, 4, 2, 10, 5 ],
            [ 5, 4, 9, 1, 3, 11, 1, 11, 10 ],
            [ 0, 9, 1, 4, 8, 5, 8, 10, 5, 8, 11, 10 ],
            [ 3, 4, 0, 3, 10, 4, 4, 10, 5, 3, 11, 10 ],
            [ 4, 8, 5, 5, 8, 10, 8, 11, 10 ], [ 9, 5, 7, 9, 7, 8 ],
            [ 0, 9, 3, 9, 5, 3, 3, 5, 7 ], [ 8, 0, 7, 0, 1, 7, 7, 1, 5 ],
            [ 1, 7, 3, 1, 5, 7 ], [ 11, 2, 3, 8, 9, 5, 8, 5, 7 ],
            [ 9, 2, 0, 9, 7, 2, 2, 7, 11, 9, 5, 7 ],
            [ 0, 3, 8, 2, 1, 11, 1, 7, 11, 1, 5, 7 ],
            [ 2, 1, 11, 11, 1, 7, 1, 5, 7 ], [ 1, 2, 10, 5, 7, 8, 5, 8, 9 ],
            [ 9, 1, 0, 10, 5, 2, 5, 3, 2, 5, 7, 3 ],
            [ 5, 2, 10, 8, 2, 5, 8, 5, 7, 8, 0, 2 ],
            [ 10, 5, 2, 2, 5, 3, 5, 7, 3 ],
            [ 3, 9, 1, 3, 8, 9, 7, 11, 10, 7, 10, 5 ],
            [ 9, 1, 0, 10, 7, 11, 10, 5, 7 ], [ 3, 8, 0, 7, 10, 5, 7, 11, 10 ],
            [ 11, 5, 7, 11, 10, 5 ], [ 11, 7, 6 ], [ 0, 8, 3, 11, 7, 6 ],
            [ 9, 0, 1, 11, 7, 6 ], [ 7, 6, 11, 3, 1, 9, 3, 9, 8 ],
            [ 2, 3, 7, 2, 7, 6 ], [ 8, 7, 0, 7, 6, 0, 0, 6, 2 ],
            [ 1, 9, 0, 3, 7, 6, 3, 6, 2 ], [ 7, 6, 2, 7, 2, 9, 2, 1, 9, 7, 9, 8 ],
            [ 1, 2, 10, 6, 11, 7 ], [ 2, 10, 1, 7, 6, 11, 8, 3, 0 ],
            [ 11, 7, 6, 10, 9, 0, 10, 0, 2 ],
            [ 7, 6, 11, 3, 2, 8, 8, 2, 10, 8, 10, 9 ],
            [ 6, 10, 7, 10, 1, 7, 7, 1, 3 ],
            [ 6, 10, 1, 6, 1, 7, 7, 1, 0, 7, 0, 8 ],
            [ 9, 0, 3, 6, 9, 3, 6, 10, 9, 6, 3, 7 ],
            [ 6, 10, 7, 7, 10, 8, 10, 9, 8 ], [ 8, 4, 6, 8, 6, 11 ],
            [ 11, 3, 6, 3, 0, 6, 6, 0, 4 ], [ 0, 1, 9, 4, 6, 11, 4, 11, 8 ],
            [ 1, 9, 4, 11, 1, 4, 11, 3, 1, 11, 4, 6 ],
            [ 3, 8, 2, 8, 4, 2, 2, 4, 6 ], [ 2, 0, 4, 2, 4, 6 ],
            [ 1, 9, 0, 3, 8, 2, 2, 8, 4, 2, 4, 6 ], [ 9, 4, 1, 1, 4, 2, 4, 6, 2 ],
            [ 10, 1, 2, 11, 8, 4, 11, 4, 6 ],
            [ 10, 1, 2, 11, 3, 6, 6, 3, 0, 6, 0, 4 ],
            [ 0, 2, 10, 0, 10, 9, 4, 11, 8, 4, 6, 11 ],
            [ 2, 11, 3, 6, 9, 4, 6, 10, 9 ],
            [ 8, 4, 6, 8, 6, 1, 6, 10, 1, 8, 1, 3 ],
            [ 1, 0, 10, 10, 0, 6, 0, 4, 6 ], [ 8, 0, 3, 9, 6, 10, 9, 4, 6 ],
            [ 10, 4, 6, 10, 9, 4 ], [ 9, 5, 4, 7, 6, 11 ],
            [ 4, 9, 5, 3, 0, 8, 11, 7, 6 ], [ 6, 11, 7, 4, 0, 1, 4, 1, 5 ],
            [ 6, 11, 7, 4, 8, 5, 5, 8, 3, 5, 3, 1 ], [ 4, 9, 5, 6, 2, 3, 6, 3, 7 ],
            [ 9, 5, 4, 8, 7, 0, 0, 7, 6, 0, 6, 2 ],
            [ 4, 0, 1, 4, 1, 5, 6, 3, 7, 6, 2, 3 ], [ 7, 4, 8, 5, 2, 1, 5, 6, 2 ],
            [ 6, 11, 7, 1, 2, 10, 9, 5, 4 ],
            [ 11, 7, 6, 8, 3, 0, 1, 2, 10, 9, 5, 4 ],
            [ 11, 7, 6, 10, 5, 2, 2, 5, 4, 2, 4, 0 ],
            [ 7, 4, 8, 2, 11, 3, 10, 5, 6 ],
            [ 4, 9, 5, 6, 10, 7, 7, 10, 1, 7, 1, 3 ],
            [ 5, 6, 10, 0, 9, 1, 8, 7, 4 ], [ 5, 6, 10, 7, 0, 3, 7, 4, 0 ],
            [ 10, 5, 6, 4, 8, 7 ], [ 5, 6, 9, 6, 11, 9, 9, 11, 8 ],
            [ 0, 9, 5, 0, 5, 3, 3, 5, 6, 3, 6, 11 ],
            [ 0, 1, 5, 0, 5, 11, 5, 6, 11, 0, 11, 8 ],
            [ 11, 3, 6, 6, 3, 5, 3, 1, 5 ], [ 9, 5, 6, 3, 9, 6, 3, 8, 9, 3, 6, 2 ],
            [ 5, 6, 9, 9, 6, 0, 6, 2, 0 ], [ 0, 3, 8, 2, 5, 6, 2, 1, 5 ],
            [ 1, 6, 2, 1, 5, 6 ], [ 1, 2, 10, 5, 6, 9, 9, 6, 11, 9, 11, 8 ],
            [ 1, 0, 9, 6, 10, 5, 11, 3, 2 ], [ 6, 10, 5, 2, 8, 0, 2, 11, 8 ],
            [ 3, 2, 11, 10, 5, 6 ], [ 10, 5, 6, 9, 3, 8, 9, 1, 3 ],
            [ 0, 9, 1, 5, 6, 10 ], [ 8, 0, 3, 10, 5, 6 ], [ 10, 5, 6 ],
            [ 10, 6, 5 ], [ 8, 3, 0, 10, 6, 5 ], [ 0, 1, 9, 5, 10, 6 ],
            [ 10, 6, 5, 9, 8, 3, 9, 3, 1 ], [ 3, 11, 2, 10, 6, 5 ],
            [ 6, 5, 10, 2, 0, 8, 2, 8, 11 ], [ 1, 9, 0, 6, 5, 10, 11, 2, 3 ],
            [ 1, 10, 2, 5, 9, 6, 9, 11, 6, 9, 8, 11 ], [ 1, 2, 6, 1, 6, 5 ],
            [ 0, 8, 3, 2, 6, 5, 2, 5, 1 ], [ 5, 9, 6, 9, 0, 6, 6, 0, 2 ],
            [ 9, 6, 5, 3, 6, 9, 3, 9, 8, 3, 2, 6 ], [ 11, 6, 3, 6, 5, 3, 3, 5, 1 ],
            [ 0, 5, 1, 0, 11, 5, 5, 11, 6, 0, 8, 11 ],
            [ 0, 5, 9, 0, 3, 5, 3, 6, 5, 3, 11, 6 ],
            [ 5, 9, 6, 6, 9, 11, 9, 8, 11 ], [ 10, 6, 5, 4, 7, 8 ],
            [ 5, 10, 6, 7, 3, 0, 7, 0, 4 ], [ 5, 10, 6, 0, 1, 9, 8, 4, 7 ],
            [ 4, 5, 9, 6, 7, 10, 7, 1, 10, 7, 3, 1 ],
            [ 7, 8, 4, 2, 3, 11, 10, 6, 5 ],
            [ 11, 6, 7, 10, 2, 5, 2, 4, 5, 2, 0, 4 ],
            [ 11, 6, 7, 8, 0, 3, 1, 10, 2, 9, 4, 5 ],
            [ 6, 7, 11, 1, 10, 2, 9, 4, 5 ], [ 7, 8, 4, 5, 1, 2, 5, 2, 6 ],
            [ 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2 ],
            [ 9, 4, 5, 8, 0, 7, 0, 6, 7, 0, 2, 6 ], [ 4, 5, 9, 6, 3, 2, 6, 7, 3 ],
            [ 6, 7, 11, 4, 5, 8, 5, 3, 8, 5, 1, 3 ],
            [ 6, 7, 11, 4, 1, 0, 4, 5, 1 ], [ 4, 5, 9, 3, 8, 0, 11, 6, 7 ],
            [ 9, 4, 5, 7, 11, 6 ], [ 10, 6, 4, 10, 4, 9 ],
            [ 8, 3, 0, 9, 10, 6, 9, 6, 4 ], [ 1, 10, 0, 10, 6, 0, 0, 6, 4 ],
            [ 8, 6, 4, 8, 1, 6, 6, 1, 10, 8, 3, 1 ],
            [ 2, 3, 11, 6, 4, 9, 6, 9, 10 ],
            [ 0, 10, 2, 0, 9, 10, 4, 8, 11, 4, 11, 6 ],
            [ 10, 2, 1, 11, 6, 3, 6, 0, 3, 6, 4, 0 ],
            [ 10, 2, 1, 11, 4, 8, 11, 6, 4 ], [ 9, 1, 4, 1, 2, 4, 4, 2, 6 ],
            [ 1, 0, 9, 3, 2, 8, 2, 4, 8, 2, 6, 4 ], [ 2, 4, 0, 2, 6, 4 ],
            [ 3, 2, 8, 8, 2, 4, 2, 6, 4 ],
            [ 1, 4, 9, 11, 4, 1, 11, 1, 3, 11, 6, 4 ],
            [ 0, 9, 1, 4, 11, 6, 4, 8, 11 ], [ 11, 6, 3, 3, 6, 0, 6, 4, 0 ],
            [ 8, 6, 4, 8, 11, 6 ], [ 6, 7, 10, 7, 8, 10, 10, 8, 9 ],
            [ 9, 3, 0, 6, 3, 9, 6, 9, 10, 6, 7, 3 ],
            [ 6, 1, 10, 6, 7, 1, 7, 0, 1, 7, 8, 0 ],
            [ 6, 7, 10, 10, 7, 1, 7, 3, 1 ],
            [ 7, 11, 6, 3, 8, 2, 8, 10, 2, 8, 9, 10 ],
            [ 11, 6, 7, 10, 0, 9, 10, 2, 0 ], [ 2, 1, 10, 7, 11, 6, 8, 0, 3 ],
            [ 1, 10, 2, 6, 7, 11 ], [ 7, 2, 6, 7, 9, 2, 2, 9, 1, 7, 8, 9 ],
            [ 1, 0, 9, 3, 6, 7, 3, 2, 6 ], [ 8, 0, 7, 7, 0, 6, 0, 2, 6 ],
            [ 2, 7, 3, 2, 6, 7 ], [ 7, 11, 6, 3, 9, 1, 3, 8, 9 ],
            [ 9, 1, 0, 11, 6, 7 ], [ 0, 3, 8, 11, 6, 7 ], [ 11, 6, 7 ],
            [ 11, 7, 5, 11, 5, 10 ], [ 3, 0, 8, 7, 5, 10, 7, 10, 11 ],
            [ 9, 0, 1, 10, 11, 7, 10, 7, 5 ],
            [ 3, 1, 9, 3, 9, 8, 7, 10, 11, 7, 5, 10 ],
            [ 10, 2, 5, 2, 3, 5, 5, 3, 7 ],
            [ 5, 10, 2, 8, 5, 2, 8, 7, 5, 8, 2, 0 ],
            [ 9, 0, 1, 10, 2, 5, 5, 2, 3, 5, 3, 7 ],
            [ 1, 10, 2, 5, 8, 7, 5, 9, 8 ], [ 2, 11, 1, 11, 7, 1, 1, 7, 5 ],
            [ 0, 8, 3, 2, 11, 1, 1, 11, 7, 1, 7, 5 ],
            [ 9, 0, 2, 9, 2, 7, 2, 11, 7, 9, 7, 5 ],
            [ 11, 3, 2, 8, 5, 9, 8, 7, 5 ], [ 1, 3, 7, 1, 7, 5 ],
            [ 8, 7, 0, 0, 7, 1, 7, 5, 1 ], [ 0, 3, 9, 9, 3, 5, 3, 7, 5 ],
            [ 9, 7, 5, 9, 8, 7 ], [ 4, 5, 8, 5, 10, 8, 8, 10, 11 ],
            [ 3, 0, 4, 3, 4, 10, 4, 5, 10, 3, 10, 11 ],
            [ 0, 1, 9, 4, 5, 8, 8, 5, 10, 8, 10, 11 ],
            [ 5, 9, 4, 1, 11, 3, 1, 10, 11 ],
            [ 3, 8, 4, 3, 4, 2, 2, 4, 5, 2, 5, 10 ],
            [ 10, 2, 5, 5, 2, 4, 2, 0, 4 ], [ 0, 3, 8, 5, 9, 4, 10, 2, 1 ],
            [ 2, 1, 10, 9, 4, 5 ], [ 8, 4, 5, 2, 8, 5, 2, 11, 8, 2, 5, 1 ],
            [ 3, 2, 11, 1, 4, 5, 1, 0, 4 ], [ 9, 4, 5, 8, 2, 11, 8, 0, 2 ],
            [ 11, 3, 2, 9, 4, 5 ], [ 4, 5, 8, 8, 5, 3, 5, 1, 3 ],
            [ 5, 0, 4, 5, 1, 0 ], [ 3, 8, 0, 4, 5, 9 ], [ 9, 4, 5 ],
            [ 7, 4, 11, 4, 9, 11, 11, 9, 10 ],
            [ 3, 0, 8, 7, 4, 11, 11, 4, 9, 11, 9, 10 ],
            [ 11, 7, 4, 1, 11, 4, 1, 10, 11, 1, 4, 0 ],
            [ 8, 7, 4, 11, 1, 10, 11, 3, 1 ],
            [ 2, 3, 7, 2, 7, 9, 7, 4, 9, 2, 9, 10 ],
            [ 4, 8, 7, 0, 10, 2, 0, 9, 10 ], [ 2, 1, 10, 0, 7, 4, 0, 3, 7 ],
            [ 10, 2, 1, 8, 7, 4 ], [ 2, 11, 7, 2, 7, 1, 1, 7, 4, 1, 4, 9 ],
            [ 3, 2, 11, 4, 8, 7, 9, 1, 0 ], [ 7, 4, 11, 11, 4, 2, 4, 0, 2 ],
            [ 2, 11, 3, 7, 4, 8 ], [ 9, 1, 4, 4, 1, 7, 1, 3, 7 ],
            [ 1, 0, 9, 8, 7, 4 ], [ 3, 4, 0, 3, 7, 4 ], [ 8, 7, 4 ],
            [ 8, 9, 10, 8, 10, 11 ], [ 0, 9, 3, 3, 9, 11, 9, 10, 11 ],
            [ 1, 10, 0, 0, 10, 8, 10, 11, 8 ], [ 10, 3, 1, 10, 11, 3 ],
            [ 3, 8, 2, 2, 8, 10, 8, 9, 10 ], [ 9, 2, 0, 9, 10, 2 ],
            [ 8, 0, 3, 1, 10, 2 ], [ 10, 2, 1 ], [ 2, 11, 1, 1, 11, 9, 11, 8, 9 ],
            [ 11, 3, 2, 0, 9, 1 ], [ 11, 0, 2, 11, 8, 0 ], [ 11, 3, 2 ],
            [ 8, 1, 3, 8, 9, 1 ], [ 9, 1, 0 ], [ 8, 0, 3 ], [] ];
            
};

export let MarchingCube = new MarchingCubeInitializer();
//each webworker needs its own marching cube object


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
