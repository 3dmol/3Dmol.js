import {Color} from "../core/Color";

declare type Shader = {
    fragmentShader: string;
    vertexShader: string;
    uniforms: Record<string, Color | number | []>
}

declare module "*.vert" {
    const value: string;
    export default value;
}

declare module "*.frag" {
    const value: string;
    export default value;
}