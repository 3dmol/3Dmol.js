import { Shader } from "../../shaders";
import { uniforms } from "./uniforms";
import fragmentShader from "./sphereimposteroutline.frag";
import vertexShader from "./sphereimposteroutline.vert";
//import fs from 'fs';
//const fragmentShader = fs.readFileSync(__dirname + "/sphereimposteroutline.frag", "utf8");
//const vertexShader = fs.readFileSync(__dirname + "/sphereimposteroutline.vert", "utf8");

export const sphereimposteroutline: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}