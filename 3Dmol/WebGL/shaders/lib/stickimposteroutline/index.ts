import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import stickimposterFragmentShaderStart from "../../utils/stickimposterFragmentShader.partial.frag";
import vertexShader from "./stickimposteroutline.vert";

//import fs from 'fs';
//const stickimposterFragmentShaderStart = fs.readFileSync(__dirname + '/stickimposterFragmentShader.partial.frag', 'utf8');
//const vertexShader = fs.readFileSync(__dirname + '/stickimposteroutline.vert', 'utf8');

const fragmentShader = stickimposterFragmentShaderStart + 'gl_FragColor = vec4(color,1.0);}';

export const stickimposteroutline: Shader = {
    fragmentShader: fragmentShader.replace('#define GLSLIFY 1', ''),
    vertexShader: vertexShader.replace('#define GLSLIFY 1', ''),
    uniforms
}