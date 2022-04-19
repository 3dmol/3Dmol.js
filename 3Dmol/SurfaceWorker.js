// @ts-nocheck
import { MarchingCubeInitializer } from "./MarchingCube";
import ProteinSurface from "./ProteinSurface4";

// Hackish way to create webworker (independent of $3Dmol namespace) within minified file
let workerStringBody = (function(){

    self.onmessage = (function(oEvent) {
        const obj = oEvent.data;
        const {type} = obj;
        if (type < 0) // sending atom data, initialize
        {
            self.atomData = obj.atoms;
            self.volume = obj.volume;
            self.ps = new ProteinSurface();  // jshint ignore:line
        } else {
            const {ps} = self;
            ps.initparm(obj.expandedExtent, type !== 1, self.volume);
            ps.fillvoxels(self.atomData, obj.extendedAtoms);
            ps.buildboundary();
            if (type === 4 || type === 2) {
                ps.fastdistancemap();
                ps.boundingatom(false);
                ps.fillvoxelswaals(self.atomData, obj.extendedAtoms);    
            }        
            ps.marchingcube(type);
            const VandF = ps.getFacesAndVertices(obj.atomsToShow);
            self.postMessage(VandF);
        }
    });
    
}).toString().replace(/(^.*?\{|\}$)/g, "");

// NOTE: variable replacement is simplified
// (See: http://stackoverflow.com/questions/1661197/what-characters-are-valid-for-javascript-variable-names)
workerStringBody += `; var ProteinSurface=${  ProteinSurface.toString().replace(/[a-zA-Z_$]{1}[0-9a-zA-Z_$]*.MarchingCube./g, "MarchingCube.")}`;
workerStringBody += `,MarchingCube=(${MarchingCubeInitializer.toString() })();`;
export const workerString = `${workerStringBody}`;
export const SurfaceWorker = (window.URL && window.URL.createObjectURL) ? window.URL.createObjectURL(new Blob([workerString],{type: 'text/javascript'})) : {postMessage(){}};
