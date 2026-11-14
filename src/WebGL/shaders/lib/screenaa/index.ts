import { Shader } from '../../ShaderType';
import { uniforms } from "./uniforms"
import fragmentShader from "./screenaa.frag"
import vertexShader from "../screen/screen.vert"

export const screenaa: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}
