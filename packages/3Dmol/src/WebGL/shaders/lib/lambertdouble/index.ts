import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './lambertdouble.frag';
import vertexShader from './lambertdouble.vert';
//import fs from 'fs';
//
//const fragmentShader = fs.readFileSync(__dirname + '/lambertdouble.frag', 'utf8');
//const vertexShader = fs.readFileSync(__dirname + '/lambertdouble.vert', 'utf8');

export const lambertdouble: Shader = {
    fragmentShader: fragmentShader.replace('#define GLSLIFY 1', ''),
    vertexShader: vertexShader.replace('#define GLSLIFY 1', ''),
    uniforms
}