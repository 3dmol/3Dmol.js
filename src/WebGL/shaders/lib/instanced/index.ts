import { Shader } from '../../ShaderType';
import { uniforms } from './uniforms';
import fragmentShader from './instanced.frag';
import vertexShader from './instanced.vert';

//import fs from "fs"
//
//const fragmentShader = fs.readFileSync(__dirname + "./instanced.frag", "utf8")
//const vertexShader = fs.readFileSync(__dirname + "./instanced.vert", "utf8") 

export const instanced: Shader = {
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    uniforms,
}