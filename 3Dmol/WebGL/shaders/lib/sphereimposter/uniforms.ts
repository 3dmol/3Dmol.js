import { Color } from '../../../core/Color';
import { Uniforms } from '../../shaders';

export const uniforms: Uniforms = {
    opacity: { type: 'f', value: 1.0 },
    fogColor: { type: 'c', value: new Color(1.0, 1.0, 1.0) },
    fogNear: { type: 'f', value: 1.0 },
    fogFar: { type: 'f', value: 2000},
    directionalLightColor: { type: 'fv', value: [] },
    directionalLightDirection: { type: 'fv', value: [] }
}