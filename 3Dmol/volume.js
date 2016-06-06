/**
 * $3Dmol.VolumeData stores volumetric data. This includes file parsing
 * functionality.
 * 
 * @class
 * @param {string} str - volumetric data
 * @param {string} format - format of supplied data (cube)
 * @param {Object} options - normalize (zero mean, unit variance), negate
 */
$3Dmol.VolumeData = function(str, format, options) {

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

    this.matrix = null; //if set must transform data
    format = format.toLowerCase();
    
    if(/\.gz$/.test(format)) {
        //unzip gzipped files
        format = format.replace(/\.gz$/,'');
        try {
            str = pako.inflate(str);
        } catch(err) {
            console.log(err);
        }
    }
    
    if (this[format]) {
        this[format](str);
    }
    
    if(options) {
        if(options.negate) {
            for(var i = 0, n = this.data.length; i < n; i++) {
                this.data[i] = -this.data[i];
            }
        }
        if(options.normalize) {
            var total = 0.0;
            for(var i = 0, n = this.data.length; i < n; i++) {
                total += this.data[i];
            }
            var mean = total/this.data.length;
            console.log("computed mean: "+mean);
            total = 0;
            for(var i = 0, n = this.data.length; i < n; i++) {
                var diff = this.data[i]-mean;
                total += diff*diff; //variance is ave of squared difference with mean
            }
            var variance = total/this.data.length;
            console.log("Computed variance: "+variance);
            //now normalize
            for(var i = 0, n = this.data.length; i < n; i++) {
                this.data[i] = (this.data[i]-mean)/variance;
            }
        }
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

$3Dmol.VolumeData.prototype.getCoordinates = function(index){
    
    var x = index/(this.size.y*this.size.z);
    var y = index % (this.size.y*this.size.z);
    var z = index % this.size.z;

    x *= this.unit.x;
    y *= this.unit.y;
    z *= this.unit.z;

    x += this.origin.x;
    y += this.origin.y;
    z += this.origin.z;
    //console.log("getCoordinates : "+x+" , "+y+" , "+z);
    //console.log("val : "+data.getVal(x,y,z));
    //console.log("actual : "+data.data[index]);
    return {x:x,y:y,z:z};
    /*
    Y = (int)(index / Width)
    X = index - (Y * Width)
    */
}

/*
 * parse vasp data
 * Essentially this parser converts the CHGCAR data into
 * cube data. It has been adapted from 'chg2cube.pl' found in
 * http://theory.cm.utexas.edu/vtsttools/
 */
$3Dmol.VolumeData.prototype.vasp = function(str) {

    var lines = str.replace(/^\s+/, "").split(/[\n\r]/);

    var atomicData = $3Dmol.Parsers.vasp(str)[0];
    var natoms = atomicData.length;

    if (natoms == 0) {
      console.log("No good formating of CHG or CHGCAR file, not atomic information provided in the file.");
      this.data = [];
      return;
    }



    // Assume atomic units 
    var unittype = "bohr/hartree";
    var l_units = 1.889725992;
    var e_units = 0.036749309;

    // copied from $3Dmol.Parsers.vasp
    var convFactor = parseFloat(lines[1]);
    // This is how Vasp reads in the basis We need the l_units in order to
    // compute the volume of the cell. Afterwards to obtain the axis for the
    // voxels we have to remove this unit and divide by the number of voxels in
    // each dimension
    var v;
    v=lines[2].replace(/^\s+/, "").split(/\s+/);
    var xVec = new $3Dmol.Vector3(parseFloat(v[0]),parseFloat(v[1]),parseFloat(v[2])).multiplyScalar(convFactor*l_units);
    v=lines[3].replace(/^\s+/, "").split(/\s+/);
    var yVec = new $3Dmol.Vector3(parseFloat(v[0]),parseFloat(v[1]),parseFloat(v[2])).multiplyScalar(convFactor*l_units);
    v=lines[4].replace(/^\s+/, "").split(/\s+/);
    var zVec = new $3Dmol.Vector3(parseFloat(v[0]),parseFloat(v[1]),parseFloat(v[2])).multiplyScalar(convFactor*l_units);

    // correct volume for non-orthognal box (expansion by minors)
    var vol = xVec.x*(yVec.y*zVec.z - zVec.y*yVec.z) - yVec.x*(xVec.y*zVec.z - zVec.y*xVec.z) + zVec.x*(xVec.y*yVec.z - yVec.y*xVec.z);

    vol = Math.abs(vol)/(Math.pow(l_units,3));
    var vol_scale = 1.0/(vol); //This Only for CHGCAR files

    // We splice the structure information
    // 2 (header) + 3 (vectors) + 2 (atoms) + 1 (vaspMode) + natoms (coords) + 1 (blank line) 
    lines.splice(0,2+3+2+1+natoms+1);


    var lineArr = lines[0].replace(/^\s+/, "").replace(/\s+/g, " ").split(" ");

    var nX = Math.abs(lineArr[0]);
    var nY = Math.abs(lineArr[1]);
    var nZ = Math.abs(lineArr[2]);

    var origin = this.origin = new $3Dmol.Vector3(0,0,0);

    this.size = {x:nX, y:nY, z:nZ};
    this.unit = new $3Dmol.Vector3(xVec.x, yVec.y, zVec.z);

    // resize the vectors accordingly
    xVec = xVec.multiplyScalar(1/(l_units*nX));
    yVec = yVec.multiplyScalar(1/(l_units*nY));
    zVec = zVec.multiplyScalar(1/(l_units*nZ));

    if (xVec.y != 0 || xVec.z != 0 || yVec.x != 0 || yVec.z != 0 || zVec.x != 0
            || zVec.y != 0) {
        //need a transformation matrix
        this.matrix =  new $3Dmol.Matrix4(xVec.x, yVec.x, zVec.x, 0, xVec.y, yVec.y, zVec.y, 0, xVec.z, yVec.z, zVec.z, 0, 0,0,0,1);
        //include translation in matrix
        this.matrix = this.matrix.multiplyMatrices(this.matrix, 
                new $3Dmol.Matrix4().makeTranslation(origin.x, origin.y, origin.z));
        //all translation and scaling done by matrix, so reset origin and unit
        this.origin = new $3Dmol.Vector3(0,0,0);
        this.unit = new $3Dmol.Vector3(1,1,1);
    }


    lines.splice(0,1); //Remove the dimension line 
    var raw = lines.join(" ");

    raw = raw.replace(/^\s+/,'');
    raw = raw.split(/[\s\r]+/);
    raw.splice(nX*nY*nZ+1);

    var preConvertedData = new Float32Array(raw); //We still have to format it to get the density

    for (var i = 0; i< preConvertedData.length; i++){
      preConvertedData[i] = preConvertedData[i]*vol_scale*e_units;
    }

    this.data = preConvertedData;

    //console.log(xVec);
    //console.log(yVec);
    //console.log(zVec);
    //console.log(this.unit);
    //console.log(this.origin);
    //console.log(this.matrix);
    //console.log(this.data);

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
            || zVec.y != 0) {
        //need a transformation matrix
        this.matrix =  new $3Dmol.Matrix4(xVec.x, yVec.x, zVec.x, 0, xVec.y, yVec.y, zVec.y, 0, xVec.z, yVec.z, zVec.z, 0, 0,0,0,1);
        //include translation in matrix
        this.matrix = this.matrix.multiplyMatrices(this.matrix, 
                new $3Dmol.Matrix4().makeTranslation(origin.x, origin.y, origin.z));
        //all translation and scaling done by matrix, so reset origin and unit
        this.origin = new $3Dmol.Vector3(0,0,0);
        this.unit = new $3Dmol.Vector3(1,1,1);
    }
    
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
    //code from ngl: https://github.com/arose/ngl/blob/master/js/ngl/parser.js
    var header = {};
    bin = new Int8Array(bin);
    var intView = new Int32Array( bin.buffer, 0, 56 );
    var floatView = new Float32Array( bin.buffer, 0, 56 );
    var dv = new DataView( bin.buffer );
    

    // 53  MAP         Character string 'MAP ' to identify file type
    header.MAP = String.fromCharCode(
        dv.getUint8( 52 * 4 ), dv.getUint8( 52 * 4 + 1 ),
        dv.getUint8( 52 * 4 + 2 ), dv.getUint8( 52 * 4 + 3 )
    );

    // 54  MACHST      Machine stamp indicating machine type which wrote file
    //                 17 and 17 for big-endian or 68 and 65 for little-endian
    header.MACHST = [ dv.getUint8( 53 * 4 ), dv.getUint8( 53 * 4 + 1 ) ];

    // swap byte order when big endian
    if( header.MACHST[ 0 ] === 17 && header.MACHST[ 1 ] === 17 ){
        var n = bin.byteLength;
        for( var i = 0; i < n; i+=4 ){
            dv.setFloat32( i, dv.getFloat32( i ), true );
        }
    }

    header.NX = intView[ 0 ];  // NC - columns (fastest changing)
    header.NY = intView[ 1 ];  // NR - rows
    header.NZ = intView[ 2 ];  // NS - sections (slowest changing)

    // mode
    //  0 image : signed 8-bit bytes range -128 to 127
    //  1 image : 16-bit halfwords
    //  2 image : 32-bit reals
    //  3 transform : complex 16-bit integers
    //  4 transform : complex 32-bit reals
    //  6 image : unsigned 16-bit range 0 to 65535
    // 16 image: unsigned char * 3 (for rgb data, non-standard)
    //
    // Note: Mode 2 is the normal mode used in the CCP4 programs.
    //       Other modes than 2 and 0 may NOT WORK
    header.MODE = intView[ 3 ];

    // start
    header.NXSTART = intView[ 4 ];  // NCSTART - first column
    header.NYSTART = intView[ 5 ];  // NRSTART - first row
    header.NZSTART = intView[ 6 ];  // NSSTART - first section

    // intervals
    header.MX = intView[ 7 ];  // intervals along x
    header.MY = intView[ 8 ];  // intervals along y
    header.MZ = intView[ 9 ];  // intervals along z

    // cell length (Angstroms in CCP4)
    header.xlen = floatView[ 10 ];
    header.ylen = floatView[ 11 ];
    header.zlen = floatView[ 12 ];

    // cell angle (Degrees)
    header.alpha = floatView[ 13 ];
    header.beta  = floatView[ 14 ];
    header.gamma = floatView[ 15 ];

    // axis correspondence (1,2,3 for X,Y,Z)
    header.MAPC = intView[ 16 ];  // column
    header.MAPR = intView[ 17 ];  // row
    header.MAPS = intView[ 18 ];  // section

    // density statistics
    header.DMIN  = floatView[ 19 ];
    header.DMAX  = floatView[ 20 ];
    header.DMEAN = floatView[ 21 ];

    // space group number 0 or 1 (default=0)
    header.ISPG = intView[ 22 ];

    // number of bytes used for symmetry data (0 or 80)
    header.NSYMBT = intView[ 23 ];

    // Flag for skew transformation, =0 none, =1 if foll
    header.LSKFLG = intView[ 24 ];

    // 26-34  SKWMAT  Skew matrix S (in order S11, S12, S13, S21 etc) if
    //                LSKFLG .ne. 0.
    // 35-37  SKWTRN  Skew translation t if LSKFLG != 0.
    //                Skew transformation is from standard orthogonal
    //                coordinate frame (as used for atoms) to orthogonal
    //                map frame, as Xo(map) = S * (Xo(atoms) - t)

    // 38      future use       (some of these are used by the MSUBSX routines
    //  .          "              in MAPBRICK, MAPCONT and FRODO)
    //  .          "   (all set to zero by default)
    //  .          "
    // 52          "

    // 50-52 origin in X,Y,Z used for transforms
    header.originX = floatView[ 49 ];
    header.originY = floatView[ 50 ];
    header.originZ = floatView[ 51 ];

    // 53  MAP         Character string 'MAP ' to identify file type
    // => see top of this parser

    // 54  MACHST      Machine stamp indicating machine type which wrote file
    // => see top of this parser

    // Rms deviation of map from mean density
    header.ARMS = floatView[ 54 ];

    // 56      NLABL           Number of labels being used
    // 57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)
    console.log("Map has min,mean,average,rmsddv: "+header.DMIN+","+header.DMAX+","+header.DMEAN+","+header.ARMS);

    //create transformation matrix, code mostly copied from ngl
    var h = header;
    var basisX = [
          h.xlen,
          0,
          0
      ];

      var basisY = [
          h.ylen * Math.cos( Math.PI / 180.0 * h.gamma ),
          h.ylen * Math.sin( Math.PI / 180.0 * h.gamma ),
          0
      ];

      var basisZ = [
          h.zlen * Math.cos( Math.PI / 180.0 * h.beta ),
          h.zlen * (
                  Math.cos( Math.PI / 180.0 * h.alpha )
                  - Math.cos( Math.PI / 180.0 * h.gamma )
                  * Math.cos( Math.PI / 180.0 * h.beta )
              ) / Math.sin( Math.PI / 180.0 * h.gamma ),
          0
      ];
      basisZ[ 2 ] = Math.sqrt(
          h.zlen * h.zlen * Math.sin( Math.PI / 180.0 * h.beta ) *
          Math.sin( Math.PI / 180.0 * h.beta ) - basisZ[ 1 ] * basisZ[ 1 ]
      );

      var basis = [ 0, basisX, basisY, basisZ ];
      var nxyz = [ 0, h.MX, h.MY, h.MZ ];
      var mapcrs = [ 0, h.MAPC, h.MAPR, h.MAPS ];

      this.matrix = new $3Dmol.Matrix4();

      this.matrix.set(

          basis[ mapcrs[1] ][0] / nxyz[ mapcrs[1] ],
          basis[ mapcrs[2] ][0] / nxyz[ mapcrs[2] ],
          basis[ mapcrs[3] ][0] / nxyz[ mapcrs[3] ],
          0,

          basis[ mapcrs[1] ][1] / nxyz[ mapcrs[1] ],
          basis[ mapcrs[2] ][1] / nxyz[ mapcrs[2] ],
          basis[ mapcrs[3] ][1] / nxyz[ mapcrs[3] ],
          0,

          basis[ mapcrs[1] ][2] / nxyz[ mapcrs[1] ],
          basis[ mapcrs[2] ][2] / nxyz[ mapcrs[2] ],
          basis[ mapcrs[3] ][2] / nxyz[ mapcrs[3] ],
          0,

          0, 0, 0, 1

      );
      //include translation in matrix
      this.matrix = this.matrix.multiplyMatrices(this.matrix, 
              new $3Dmol.Matrix4().makeTranslation(
                      h.NXSTART + h.originX,
                      h.NYSTART + h.originY,
                      h.NZSTART + h.originZ));
      //all translation and scaling done by matrix, so reset origin and unit
      this.origin = new $3Dmol.Vector3(0,0,0);
      this.unit = new $3Dmol.Vector3(1,1,1); 
      this.size = {x:header.NX, y:header.NY, z:header.NZ};
      var data = new Float32Array(bin.buffer, 1024 + header.NSYMBT);
      //data must by (slowest changing) x,y,z (fastest changing)

      var NX = header.NX, NY = header.NY, NZ = header.NZ;
      this.data = new Float32Array(NX*NY*NZ);
      for(var i = 0; i < NX; i++) {
          for(var j = 0; j < NY; j++) {
              for(var k = 0; k < NZ; k++) {
                  //should I be concerned that I'm not using mapc?
                  this.data[((i*NY)+j)*NZ+k] = data[((k*NY)+j)*NX+i];
              }
          }
      }

};
