import stickimposterFragmentShaderStart from '../../utils/stickimposterFragmentShader.partial.frag';
import vertexShader from './stickimposteroutline.vert';
import { uniforms } from './uniforms';

const fragmentShader = stickimposterFragmentShaderStart + 'gl_FragColor = vec4(color,1.0);}';

export const stickimposteroutline = {
    fragmentShader,
    vertexShader,
    uniforms
}