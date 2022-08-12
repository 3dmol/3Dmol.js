import { Shader } from '../../ShaderType';
import { uniforms } from "./uniforms";
import fragmentShader from "./sphereimposteroutline.frag";
import vertexShader from "./sphereimposteroutline.vert";
//import fs from 'fs';
//const fragmentShader = fs.readFileSync(__dirname + "/sphereimposteroutline.frag", "utf8");
//const vertexShader = fs.readFileSync(__dirname + "/sphereimposteroutline.vert", "utf8");

export const sphereimposteroutline: Shader = {
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    uniforms
}