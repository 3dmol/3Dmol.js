import { Shader } from '../../ShaderType';
import { uniforms } from "./uniforms";
import fragmentShader from "./outline.frag";
import vertexShader from "./outline.vert";

export const outline: Shader = {
  fragmentShader,
  vertexShader,
  uniforms,
};
