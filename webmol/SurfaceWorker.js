//dkoes - calls protein surface as a webworker, returns the faces and vertices

var WebMol = WebMol || {};

if (WebMol.Vector3 === undefined) {
    
    WebMol.Vector3 = function(x, y, z) {
        this.x = x || 0.0;
        this.y = y || 0.0;
        this.z = z || 0.0;
    };    

    WebMol.Vector3.prototype =  {
        
        constructor : WebMol.Vector3,
        
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
}

importScripts("marchingcube.js");
importScripts("ProteinSurface4.js");
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

		if (type == 4 || type == 2) {
			ps.fastdistancemap();
            ps.boundingatom(false);
            ps.fillvoxelswaals(self.atomData, obj.extendedAtoms);	
        }		

		ps.marchingcube(type);
		var VandF = ps.getFacesAndVertices(obj.atomsToShow);
		self.postMessage(VandF);
	}
};
