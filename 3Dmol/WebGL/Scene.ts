import { Color } from './core/Color';
import { Object3D } from './core/Object3D';
import { Light } from './Light';
/*
 * Scene class
 */
/** @constructor */
export class Scene extends Object3D {
  fog: null;
  overrideMaterial: null;
  __objects: Object3D[];
  __lights: Light[];
  __objectsAdded: Object3D[];
  __objectsRemoved: Object3D[];
  constructor() {
    super();

    this.fog = null;

    //May not need...
    this.overrideMaterial = null;

    this.matrixAutoUpdate = false;

    this.__objects = [];
    this.__lights = [];

    this.__objectsAdded = [];
    this.__objectsRemoved = [];

  };

  __addObject(object) {

    //Directional Lighting
    if (object instanceof Light) {

      if (this.__lights.indexOf(object) === -1)
        this.__lights.push(object);

      //TODO: Do I need this??
      if (object.target && object.target.parent === undefined)
        this.add(object.target);

    }

    //Rotation group
    else {

      if (this.__objects.indexOf(object) === -1) {

        this.__objects.push(object);
        this.__objectsAdded.push(object);

        //Check if previously removed

        var idx = this.__objectsRemoved.indexOf(object);

        if (idx !== -1)
          this.__objectsRemoved.splice(idx, 1);

      }
    }

    //Add object's children

    for (var i = 0; i < object.children.length; i++)
      this.__addObject(object.children[i]);

  };

  __removeObject(object) {

    var idx;
    if (object instanceof Light) {

      idx = this.__lights.indexOf(object);

      if (idx !== -1)
        this.__lights.splice(idx, 1);

    }

    //Object3D
    else {

      idx = this.__objects.indexOf(object);

      if (idx !== -1) {

        this.__objects.splice(idx, 1);
        this.__objectsRemoved.push(object);

        //Check if previously added

        var ai = this.__objectsAdded.indexOf(object);

        if (ai !== -1)
          this.__objectsAdded.splice(idx, 1);

      }

    }

    //Remove object's children
    for (var i = 0; i < object.children.length; i++)
      this.__removeObject(object.children[i]);

  };
}


