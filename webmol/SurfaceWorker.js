//dkoes - calls protein surface as a webworker, returns the faces and vertices

importScripts("MarchingCubeData.js");
importScripts("ProteinSurface4.js");
self.onmessage = function(oEvent) {
	var obj = oEvent.data;
	var type = obj.type;

	if (type < 0) // sending atom data, initialize
	{
		self.atomData = obj.atoms;
		self.ps = new ProteinSurface();
	} else {
		var ps = self.ps;
		ps.initparm(obj.expandedExtent, (type == 1) ? false : true);
		ps.fillvoxels(self.atomData, obj.extendedAtoms);
		ps.buildboundary();

		if (type == 4 || type == 2)
			ps.fastdistancemap();
		if (type == 2) {
			ps.boundingatom(false);
			ps.fillvoxelswaals(self.atomData, obj.extendedAtoms);
		}

		ps.marchingcube(type);
		ps.laplaciansmooth(1);
		var VandF = ps.getFacesAndVertices(obj.atomsToShow);
		self.postMessage(VandF);
	}
};
