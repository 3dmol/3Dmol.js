import { Shader } from '../../ShaderType';
import { uniforms } from "./uniforms";
import fragmentShader from "./outline.frag";
import vertexShader from "./outline.vert";

//import fs from "fs";
//
//const fragmentShader = fs.readFileSync(__dirname + "/outline.frag", "utf8");
//const vertexShader = fs.readFileSync(__dirname + "/outline.vert", "utf8");

export const outline: Shader = {
  fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
  vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
  uniforms,
};
