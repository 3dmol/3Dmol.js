import { Shader } from '../../ShaderType';
import { uniforms } from "./uniforms"
import fragmentShader from "./screenaa.frag"
import vertexShader from "./screenaa.vert"


//import fs from "fs"
//
//const fragmentShader = fs.readFileSync(__dirname + "/screenaa.frag", "utf8")
//const vertexShader = fs.readFileSync(__dirname + "/screenaa.vert", "utf8")

export const screenaa: Shader = {
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    uniforms
}
