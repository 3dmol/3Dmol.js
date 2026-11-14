import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './ssao.frag';
import vertexShader from '../screen/screen.vert';

export const ssao: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}