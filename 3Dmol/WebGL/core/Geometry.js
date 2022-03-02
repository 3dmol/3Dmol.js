// @ts-check
import { EventDispatcher } from "./EventDispatcher";

Object.defineProperty(Geometry.prototype, "vertices", {
    
    /** @this {Geometry} */
    get : function() {
       
    } 
        
});

export var GeometryIDCount = 0;

export class Geometry extends EventDispatcher{
    constructor(mesh, radii, offset) {
        super();
        this.id = GeometryIDCount++;

        this.name = '';

        this.hasTangents = false;

        this.dynamic = true; // the intermediate typed arrays will be deleted when set to false
        this.mesh = (mesh === true) ? true : false; // Does this geometry represent a mesh (i.e. do we need Face/Line index buffers?)
        this.radii = radii || false;
        this.offset = offset || false; //offset buffer used for instancing
        // update flags

        this.verticesNeedUpdate = false;
        this.elementsNeedUpdate = false;
        this.normalsNeedUpdate = false;
        this.colorsNeedUpdate = false;

        this.buffersNeedUpdate = false;

        this.geometryGroups = [];
        this.groups = 0;

    }

    get vertices() {
        var vertices = 0;
        for (var g = 0; g < this.groups; g++)
            vertices += this.geometryGroups[g].vertices;
            
        return vertices;
    }

    updateGeoGroup(addVertices) {

        addVertices = addVertices || 0;

        var retGroup = this.groups > 0 ? this.geometryGroups[this.groups - 1] : null;

        if (!retGroup || retGroup.vertices + addVertices > retGroup.vertexArray.length / 3)
            retGroup = addGroup(this);

        return retGroup;

    }

    vrml(indent, material) {
        var ret = '';
        var len = this.geometryGroups.length;
        for (var g = 0; g < len; g++) {
            var geoGroup = this.geometryGroups[g];
            ret += geoGroup.vrml(indent, material) + ',\n';
        }
        return ret;
    }

    addGeoGroup() {
        return addGroup(this);
    }

    setUpNormals(three) {

        three = three || false;

        for (var g = 0; g < this.groups; g++) {

            var geoGroup = this.geometryGroups[g];

            geoGroup.setNormals(three);

        }

    }

    setColors(setcolor) {
        var len = this.geometryGroups.length;
        for (var g = 0; g < len; g++) {

            var geoGroup = this.geometryGroups[g];
            geoGroup.setColors(setcolor);

        }
    }

    setUpWireframe() {
        for (var g = 0; g < this.groups; g++) {
            var geoGroup = this.geometryGroups[g];

            geoGroup.setLineIndices();
        }
    }

    initTypedArrays() {

        for (var g = 0; g < this.groups; g++) {

            var group = this.geometryGroups[g];

            if (group.__inittedArrays === true)
                continue;

            //do not actually reallocate smaller memory here because
            //of the performance hit - if you know your geometry is small,
            //truncate manually with the second parameter true
            group.truncateArrayBuffers(this.mesh, false);
        }


    }

    dispose() {
        this.dispatchEvent({ type: 'dispose' });
    }
}