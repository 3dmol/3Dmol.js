import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './ssao.frag';
import vertexShader from './ssao.vert';

export const ssao: Shader = {
    fragmentShader: fragmentShader,
    vertexShader: vertexShader,
    uniforms
}