import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './lambert.frag';
import vertexShader from './lambert.vert';

export const lambert: Shader = {
    fragmentShader,
    vertexShader,
    uniforms,
}