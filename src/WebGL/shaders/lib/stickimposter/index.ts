import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import stickimposterFragmentShaderStart from './stickimposterFragmentShader.partial.frag';
import stickimposterFragmentShaderEnd from './stickimposter.partial.frag';
import vertexShader from './stickimposter.vert';

const fragmentShader = [stickimposterFragmentShaderStart, stickimposterFragmentShaderEnd].join('\n');

export const stickimposter: Shader = {
    fragmentShader: fragmentShader,
    vertexShader: vertexShader,
    uniforms
}