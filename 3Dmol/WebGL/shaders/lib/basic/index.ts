import { uniforms } from './uniforms';
import { Shader } from '../../ShaderType';
import fragmentShader from './basic.frag';
import vertexShader from './basic.vert';

export const basic: Shader = {
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    uniforms
}