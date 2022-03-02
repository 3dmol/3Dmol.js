// @ts-check

import { Color, Object3D } from "./core";
import { Vector3 } from "./math";


export class Light {
    constructor(hex, intensity) {
        
        Object3D.call(this);
        
        this.color = new Color(hex);
        this.position = new Vector3( 0, 1, 0 );
        this.target = new Object3D();

        this.intensity = ( intensity !== undefined ) ? intensity : 1;

        this.castShadow = false;
        this.onlyShadow = false;
        
    }
};

Light.prototype = Object.create(Object3D.prototype);
