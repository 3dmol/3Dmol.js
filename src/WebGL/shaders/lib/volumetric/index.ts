import { uniforms } from "./uniforms";
import { Shader } from '../../ShaderType';
import fragmentShader from "./volumetric.frag";
import vertexShader from "./volumetric.vert";
//import fs from "fs";
//
//const fragmentShader = fs.readFileSync(__dirname + "/volumetric.frag", "utf8");
//const vertexShader = fs.readFileSync(__dirname + "/volumetric.vert", "utf8");

export const volumetric: Shader = {
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    uniforms
}