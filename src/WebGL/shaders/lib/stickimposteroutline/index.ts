import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import stickimposterFragmentShaderStart from "../stickimposter/stickimposterFragmentShader.partial.frag";
import vertexShader from "./stickimposteroutline.vert";

const fragmentShader = stickimposterFragmentShaderStart + 'gl_FragColor = vec4(color,1.0);}';

export const stickimposteroutline: Shader = {
    fragmentShader: fragmentShader,
    vertexShader: vertexShader,
    uniforms
}