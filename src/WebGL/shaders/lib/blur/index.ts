import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './blur.frag';
import vertexShader from './blur.vert';

export const blur: Shader = {
    fragmentShader: fragmentShader,
    vertexShader: vertexShader,
    uniforms
}