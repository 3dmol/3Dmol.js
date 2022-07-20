import {Color} from "../core/Color";

declare type Shader = {
    fragmentShader: string;
    vertexShader: string;
    uniforms: Record<string, any>
}