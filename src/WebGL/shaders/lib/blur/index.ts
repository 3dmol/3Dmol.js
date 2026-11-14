import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './blur.frag';
import vertexShader from '../screen/screen.vert';

export const blur: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}