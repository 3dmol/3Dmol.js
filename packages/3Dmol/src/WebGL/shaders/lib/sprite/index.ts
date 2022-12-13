import { Shader } from '../../ShaderType';
import { uniforms } from "./uniforms"
//import { fragmentShader, vertexShader } from "./shaders"
import fragmentShader from "./sprite.frag"
import vertexShader from "./sprite.vert"


//const fragmentShader = fs.readFileSync(path.resove(__dirname, "./sprite.frag"), "utf8")
//const vertexShader = fs.readFileSync(path.resove(__dirname, "./sprite.vert"), "utf8")

export const sprite: Shader = {
    fragmentShader: fragmentShader.replace("#define GLSLIFY 1", ""),
    vertexShader: vertexShader.replace("#define GLSLIFY 1", ""),
    uniforms
}