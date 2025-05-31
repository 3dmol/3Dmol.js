import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import stickimposterFragmentShaderStart from "../stickimposter/stickimposterFragmentShader.partial.frag";
import stickimposterFragmentShaderEnd from './stickimposteroutline.partial.frag';

import vertexShader from "./stickimposteroutline.vert";

const fragmentShader = stickimposterFragmentShaderStart + stickimposterFragmentShaderEnd;

export const stickimposteroutline: Shader = {
    fragmentShader: fragmentShader,
    vertexShader: vertexShader,
    uniforms
}