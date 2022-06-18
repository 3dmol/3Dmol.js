import { Shader } from '../../shaders';
import { uniforms } from './uniforms';
import fragmentShader from './screen.frag';
import vertexShader from './screen.vert';

//import fs from 'fs';
//
//const fragmentShader = fs.readFileSync(__dirname + '/screen.frag', 'utf8');
//const vertexShader = fs.readFileSync(__dirname + '/screen.vert', 'utf8');

export const screen: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}