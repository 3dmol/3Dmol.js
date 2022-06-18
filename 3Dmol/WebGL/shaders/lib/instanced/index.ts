import { Shader } from '../../shaders';
import { uniforms } from './uniforms';
import fragmentShader from './instanced.frag';
import vertexShader from './instanced.vert';

//import fs from "fs"
//
//const fragmentShader = fs.readFileSync(__dirname + "./instanced.frag", "utf8")
//const vertexShader = fs.readFileSync(__dirname + "./instanced.vert", "utf8") 

export const instanced: Shader = {
    fragmentShader,
    vertexShader,
    uniforms,
}