import { uniforms } from './uniforms';
import { Shader } from '../../ShaderType';
import fragmentShader from './sphereimposter.frag';
import vertexShader from './sphereimposter.vert';

//import fs from "fs"
//
//const fragmentShader = fs.readFileSync(__dirname + "./sphereimposter.frag", "utf8");
//const vertexShader = fs.readFileSync(__dirname + "./sphereimposter.vert", "utf8");


export const sphereimposter: Shader = {
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    uniforms
}