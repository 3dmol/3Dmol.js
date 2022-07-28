import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './lambert.frag';
import vertexShader from './lambert.vert';
//import fs from 'fs';
//
//const fragmentShader = fs.readFileSync(__dirname + '/lambert.frag', 'utf8');
//const vertexShader = fs.readFileSync(__dirname + '/lambert.vert', 'utf8');

export const lambert: Shader = {
    fragmentShader: fragmentShader.replace('#define GLSLIFY 1', ''),
    vertexShader: vertexShader.replace('#define GLSLIFY 1', ''),
    uniforms,
}