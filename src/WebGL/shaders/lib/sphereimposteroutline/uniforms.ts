import { Color } from "../../../../colors";

export const uniforms = {
    opacity: { type: 'f', value: 1.0 },
    outlineColor: { type: 'c', value: new Color(0.0, 0.0, 0.0) },
    fogColor: { type: 'c', value: new Color(1.0, 1.0, 1.0) },
    fogNear: { type: 'f', value: 1.0 },
    fogFar: { type: 'f', value: 2000},
    outlineWidth: { type: 'f', value: 0.1 },
    outlinePushback: { type: 'f', value: 1.0 },
}