import { Shader } from "../../shaders";
import { uniforms } from "./uniforms";
import fragmentShader from "./outline.frag";
import vertexShader from "./outline.vert";

//import fs from "fs";
//
//const fragmentShader = fs.readFileSync(__dirname + "/outline.frag", "utf8");
//const vertexShader = fs.readFileSync(__dirname + "/outline.vert", "utf8");

export const outline: Shader = {
  fragmentShader,
  vertexShader,
  uniforms,
};
