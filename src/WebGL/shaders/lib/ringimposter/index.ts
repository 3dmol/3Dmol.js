import { uniforms } from './uniforms';
import { Shader } from '../../ShaderType';
import fragmentShader from './ringimposter.frag';
import vertexShader from './ringimposter.vert';

export const ringimposter: Shader = {
    vertexShader,
    fragmentShader,
    uniforms
}
