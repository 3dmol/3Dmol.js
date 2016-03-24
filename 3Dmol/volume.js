/**
 * $3Dmol.VolumeData stores volumetric data. This includes file parsing
 * functionality.
 * 
 * @class
 * @param {string} str - volumetric data
 * @param {string} format - format of supplied data (cube)
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
 * @param {number} x,y,z - the coordinates
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

//parse cp4 files
$3Dmol.VolumeData.prototype.ccp4 = function(bin) {

    // http://www.ccp4.ac.uk/html/maplib.html#description
        
        var header = {};

        var intView = new Int32Array( bin, 0, 56 );
        var floatView = new Float32Array( bin, 0, 56 );

        for (var i = 0; i < intView.length; i++) {
            console.log(intView[i]);
        }

           //# of columns (fastest changing),rows, sections (slowest changing)
        header.NX = intView[ 0 ];
        header.NY = intView[ 1 ];
        header.NZ = intView[ 2 ];

        // mode
        //  0 image : signed 8-bit bytes range -128 to 127
        //  1 image : 16-bit halfwords
        //  2 image : 32-bit reals
        //  3 transform : complex 16-bit integers
        //  4 transform : complex 32-bit reals
        //  6 image : unsigned 16-bit range 0 to 65535
        // 16 image: unsigned char * 3 (for rgb data, non-standard)
        header.MODE = intView[ 3 ];

        //position of first column, first row, and first section (voxel grid units)
        header.NCSTART = intView[ 4 ];
        header.NRSTART = intView[ 5 ];
        header.NSSTART = intView[ 6 ];

        //intervals per unit cell repeat along X,Y Z
        header.MX = Math.abs(intView[ 7 ]);
        header.MY = Math.abs(intView[ 8 ]);
        header.MZ = Math.abs(intView[ 9 ]);

       //Map lengths along X,Y,Z in Ångstroms
        header.xlen =  floatView[ 10 ]; 
        header.ylen = floatView[ 11 ];  
        header.zlen = floatView[ 12 ];  

        //single voxel lengths
        xVox = header.xlen / header.NX; 
        yVox = header.ylen / header.NY; 
        zVox = header.zlen / header.NZ;

        // cell angle
        header.alpha = floatView[ 13 ];
        header.beta  = floatView[ 14 ];
        header.gamma = floatView[ 15 ];

        //relationship of X,Y,Z axes to columns, rows, sections
        header.MAPC = intView[ 16 ];
        switch (header.MAPC){
            case 1: header.NXSTART = header.NCSTART; header.NYSTART = 0; header.NZSTART = 0; break;
            case 2: header.NYSTART = header.NCSTART; header.NXSTART = 0; header.NZSTART = 0; break;
            case 3: header.NZSTART = header.NCSTART; header.NXSTART = 0; header.NYSTART = 0; break;
        }
        header.MAPR = intView[ 17 ];
        switch (header.MAPR){
            case 1: header.NXSTART = header.NRSTART; header.NYSTART = 0; header.NZSTART = 0; break;
            case 2: header.NYSTART = header.NRSTART; header.NXSTART = 0; header.NZSTART = 0; break;
            case 3: header.NZSTART = header.NRSTART; header.NXSTART = 0; header.NYSTART = 0; break;
        }
        header.MAPS = intView[ 18 ];
         switch (header.MAPS){
            case 1: header.NXSTART = header.NSSTART; header.NYSTART = 0; header.NZSTART = 0; break;
            case 2: header.NYSTART = header.NSSTART; header.NXSTART = 0; header.NZSTART = 0; break;
            case 3: header.NZSTART = header.NSSTART; header.NXSTART = 0; header.NYSTART = 0; break;
        }

        //Origin Position
        this.origin = new $3Dmol.Vector3(header.NXSTART, header.NYSTART, header.NZSTART);

        //Minimum, maximum, average density
        header.DMIN  = intView[ 19 ];
        header.DMAX  = intView[ 20 ];
        header.DMEAN = intView[ 21 ];

        // space group number 0 or 1 (default=0)
        header.ISPG = intView[ 22 ];

        // number of bytes used for symmetry data (0 or 80)
        header.NSYMBT = intView[23];

        // machine stamp
        header.ARMS = floatView[54];

         this.size = {x:header.MX, y:header.MY, z:header.MZ};
         this.unit = new $3Dmol.Vector3(xVox, yVox, zVox ); 

        this.data = new Float32Array(bin, 
        256 * 4 + header.NSYMBT ,
        header.NX * header.NY * header.NZ
        );
          
};
