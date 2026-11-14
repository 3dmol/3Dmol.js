import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './screen.frag';
import vertexShader from './screen.vert';

export const screen: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}