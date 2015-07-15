//Hackish way to create webworker (independent of $3Dmol namespace) within minified file
$3Dmol.workerString = function(){

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
    
}.toString().replace(/(^.*?\{|\}$)/g, "");

$3Dmol.workerString += "; var ProteinSurface=" + $3Dmol.ProteinSurface.toString().replace(/\$3Dmol.MarchingCube./g, "MarchingCube.");
$3Dmol.workerString += ",MarchingCube=("+$3Dmol.MarchingCubeInitializer.toString() +")();";

$3Dmol.SurfaceWorker = window.URL.createObjectURL(new Blob([$3Dmol.workerString],{type: 'text/javascript'}));

$3Dmol['workerString'] = $3Dmol.workerString;
$3Dmol['SurfaceWorker'] = $3Dmol.SurfaceWorker;
