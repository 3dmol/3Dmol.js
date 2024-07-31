import { Matrix4 } from "WebGL/math";

export const uniforms = {
    total_strength: { type: 'f', value: 1.0 },
    radius: { type: 'f', value: 5},    
    projinv: { type: 'mat4', value: [] as Matrix4[] },

};
