import { uniforms } from './uniforms';
import { Shader } from '../../ShaderType';
import fragmentShader from './sphereimposter.frag';
import vertexShader from './sphereimposter.vert';

export const sphereimposter: Shader = {
    vertexShader,
    fragmentShader,
    uniforms
}