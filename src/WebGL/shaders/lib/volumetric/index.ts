import { uniforms } from "./uniforms";
import { Shader } from '../../ShaderType';
import fragmentShader from "./volumetric.frag";
import vertexShader from "./volumetric.vert";

export const volumetric: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}