import stickimposterFragmentShaderStart from '../../utils/stickimposterFragmentShader.partial.frag';
import stickimposterFragmentShaderEnd from './stickimposter.partial.frag';
import vertexShader from './stickimposter.vert';
import { uniforms } from './uniforms';

const fragmentShader = [stickimposterFragmentShaderStart, stickimposterFragmentShaderEnd].join('\n');

export const stickimposter = {
    fragmentShader,
    vertexShader,
    uniforms
}