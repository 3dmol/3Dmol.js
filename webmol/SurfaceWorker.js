
//Hackish way to create webworker (independent of WebMol namespace) within minified file
WebMol.workerString = function(){

    self.onmessage = function(oEvent) {
    	var obj = oEvent.data;
    	var type = obj.type;
    	if (type < 0) // sending atom data, initialize
    	{
    		self.atomData = obj.atoms;
    		self.volume = obj.volume;
    		self.ps = new ProteinSurface();
    	} else {
    		var ps = self.ps;
    		ps.initparm(obj.expandedExtent, (type == 1) ? false : true, self.volume);
    		ps.fillvoxels(self.atomData, obj.extendedAtoms);
    		ps.buildboundary();
    		if (type === 4 || type === 2) {
    			ps.fastdistancemap();
                ps.boundingatom(false);
                ps.fillvoxelswaals(self.atomData, obj.extendedAtoms);	
            }		
    		ps.marchingcube(type);
    		var VandF = ps.getFacesAndVertices(obj.atomsToShow);
    		self.postMessage(VandF);
    	}
    };
    
    var ISDONE = 2;
    
    var Vector3 = function(x, y, z) {
        this.x = x || 0.0;
        this.y = y || 0.0;
        this.z = z || 0.0;
    };    
    
    Vector3.prototype =  {
        
        constructor : Vector3,
        
        copy : function(v) {
            
            this.x = v.x;
            this.y = v.y;
            this.z = v.z;
            
            return this;  
        },
        
        multiplyScalar : function(s) {
            
            this.x *= s;
            this.y *= s;
            this.z *= s;
            
            return this;
        }
    }; 
    
}.toString().replace(/(^.*?\{|\}$)/g, "");

WebMol.workerString += ",ISDONE=2";
WebMol.workerString += ",ProteinSurface=" + WebMol.ProteinSurface.toString().replace(/WebMol.MarchingCube./g, "");
WebMol.workerString += ",march=" + WebMol.MarchingCube.march.toString().replace(/WebMol./g, "");
WebMol.workerString += ",laplacianSmooth=" + WebMol.MarchingCube.laplacianSmooth.toString();

WebMol.workerString += ",edgeTable=new Uint32Array([" + WebMol.MarchingCube.edgeTable.toString() + "])";

WebMol.workerString += ",triTable=[";

for (var i = 0, il = WebMol.MarchingCube.triTable.length; i < il - 1; i++)
    WebMol.workerString += "[" + WebMol.MarchingCube.triTable[i].toString() + "],";

WebMol.workerString += "[]]";


WebMol.SurfaceWorker = window.URL.createObjectURL(new Blob([WebMol.workerString]));

 
