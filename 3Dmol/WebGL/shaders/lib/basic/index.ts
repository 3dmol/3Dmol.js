import { uniforms } from './uniforms';
import { Shader } from '../../shaders';
import fragmentShader from './basic.frag';
import vertexShader from './basic.vert';

export const basic: Shader = {
    vertexShader,
    fragmentShader,
    uniforms
}