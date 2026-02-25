import { uniforms } from './uniforms';
import { Shader } from '../../ShaderType';
import fragmentShader from './ringimposter.frag';
import vertexShader from './ringimposter.vert';

export const ringimposter: Shader = {
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    uniforms
}
