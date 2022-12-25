import { Color } from "../../../colors";

export type Uniform = Record<any, { type?: unknown; value?: unknown }>

export function clone(
  uniforms_src: Uniform
) {
  let uniforms_clone: Uniform = {};

  for (const u in uniforms_src) {
    uniforms_clone[u] = {};
    uniforms_clone[u].type = uniforms_src[u].type;

    var srcValue = uniforms_src[u].value;

    if (srcValue instanceof Color) uniforms_clone[u].value = srcValue.clone();
    else if (typeof srcValue === "number") uniforms_clone[u].value = srcValue;
    else if (srcValue instanceof Array) uniforms_clone[u].value = [];
    else
      console.error(
        "Error copying shader uniforms from ShaderLib: unknown type for uniform"
      );
  }

  return uniforms_clone;
}
