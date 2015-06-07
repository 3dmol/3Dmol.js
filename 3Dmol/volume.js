/**
 * $3Dmol.VolumeData stores volumetric data. This includes file parsing
 * functionality.
 * 
 * @class
 */
$3Dmol.VolumeData = function(str, format) {

    this.unit = {
        x : 1,
        y : 1,
        z : 1
    }; // scale of each voxel
    this.origin = {
        x : 0,
        y : 0,
        z : 0
    }; // origin (bottom "left", not center)
    this.size = {
        x : 0,
        y : 0,
        z : 0
    }; // number of voxels in each direction
    this.data = new Float32Array([]); // actual floating point data, arranged
                                        // x->y->z

    format = format.toLowerCase();
    if (this[format]) {
        this[format](str);
    }
};

/**
 * @function $3Dmol.VolumeData.getVal
 * @param x,y,z - the coordinates
 * @returns - value closest to provided coordinate; zero if coordinate invalid
 */
$3Dmol.VolumeData.prototype.getVal = function(x,y,z) {
    x -= this.origin.x;
    y -= this.origin.y;
    z -= this.origin.z;
    
    x /= this.unit.x;
    y /= this.unit.y;
    z /= this.unit.z;
    
    x = Math.round(x);
    y = Math.round(y);
    z = Math.round(z);
    
    if(x < 0 || x >= this.size.x) return 0;
    if(y < 0 || y >= this.size.y) return 0;
    if(z < 0 || z >= this.size.z) return 0;
    
    return this.data[x*this.size.y*this.size.z + y*this.size.z + z];
};

// parse cube data
$3Dmol.VolumeData.prototype.cube = function(str) {
    var lines = str.replace(/^\s+/, "").split(/[\n\r]+/);

    if (lines.length < 6)
        return;

    var lineArr = lines[2].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");

    var atomsnum = parseFloat(lineArr[0]); //includes sign, which indicates presence of oribital line in header
    var natoms = Math.abs(atomsnum);

    var origin = this.origin = new $3Dmol.Vector3(parseFloat(lineArr[1]),
            parseFloat(lineArr[2]), parseFloat(lineArr[3]));

    lineArr = lines[3].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");

    // might have to convert from bohr units to angstroms
    // there is a great deal of confusion here:
    // n>0 means angstroms: http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm
    // n<0 means angstroms: http://paulbourke.net/dataformats/cube/
    // always assume bohr: openbabel source code
    // always assume angstrom: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
    // we are going to go with n<0 means angstrom - note this is just the first n
    var convFactor = (lineArr[0] > 0) ? 0.529177 : 1;
    origin.multiplyScalar(convFactor);

    var nX = Math.abs(lineArr[0]);
    var xVec = new $3Dmol.Vector3(parseFloat(lineArr[1]),
            parseFloat(lineArr[2]), parseFloat(lineArr[3]))
            .multiplyScalar(convFactor);

    lineArr = lines[4].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
    var nY = Math.abs(lineArr[0]);
    var yVec = new $3Dmol.Vector3(parseFloat(lineArr[1]),
            parseFloat(lineArr[2]), parseFloat(lineArr[3]))
            .multiplyScalar(convFactor);

    lineArr = lines[5].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");
    var nZ = Math.abs(lineArr[0]);
    var zVec = new $3Dmol.Vector3(parseFloat(lineArr[1]),
            parseFloat(lineArr[2]), parseFloat(lineArr[3]))
            .multiplyScalar(convFactor);

    this.size = {x:nX, y:nY, z:nZ};
    this.unit = new $3Dmol.Vector3(xVec.x, yVec.y, zVec.z);
    
    if (xVec.y != 0 || xVec.z != 0 || yVec.x != 0 || yVec.z != 0 || zVec.x != 0
            || zVec.y != 0)
        console
                .log("Warning: Cube file is not axis aligned.  This isn't going to look right.");   
    var headerlines = 6;
    if(atomsnum < 0) headerlines++; //see: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html
    var raw = lines.splice(natoms + headerlines).join(" ");
    raw = raw.replace(/^\s+/,'');
    raw = raw.split(/[\s\r]+/);
    this.data = new Float32Array(raw);

};
