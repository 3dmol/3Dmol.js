import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './instanced.frag';
import vertexShader from './instanced.vert';

export const instanced: Shader = {
    fragmentShader,
    vertexShader,
    uniforms,
}