import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './screen.frag';
import vertexShader from './screen.vert';

export const screen: Shader = {
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    uniforms
}