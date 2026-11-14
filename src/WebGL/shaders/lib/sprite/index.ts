import { Shader } from '../../ShaderType';
import { uniforms } from "./uniforms"
import fragmentShader from "./sprite.frag"
import vertexShader from "./sprite.vert"

export const sprite: Shader = {
    fragmentShader,
    vertexShader,
    uniforms
}