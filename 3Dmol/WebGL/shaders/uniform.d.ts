import type { Color } from 'three';

declare type Uniform = {
    fragmentShader: string;
    vertexShader: string;
    uniforms: Record<string, Color | number | []>
}