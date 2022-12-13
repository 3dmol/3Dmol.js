import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import stickimposterFragmentShaderStart from '../../utils/stickimposterFragmentShader.partial.frag';
import stickimposterFragmentShaderEnd from './stickimposter.partial.frag';
import vertexShader from './stickimposter.vert';

//import fs from "fs"

//const stickimposterFragmentShaderStart = fs.readFileSync(__dirname + "../../utils/stickimposterFragmentShader.partial.frag", "utf8");
//const stickimposterFragmentShaderEnd = fs.readFileSync(__dirname + "./stickimposter.partial.frag", "utf8");
//const vertexShader = fs.readFileSync(__dirname + "./stickimposter.vert", "utf8");

const fragmentShader = [stickimposterFragmentShaderStart, stickimposterFragmentShaderEnd].join('\n');

export const stickimposter: Shader = {
    fragmentShader: fragmentShader.replace('#define GLSLIFY 1', ''),
    vertexShader: vertexShader.replace('#define GLSLIFY 1', ''),
    uniforms
}