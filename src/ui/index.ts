import { UI } from './ui';
import { StateManager } from './state';
import { LearningSession } from './learning';

export { UI, StateManager, LearningSession };
export type { LearningConfig, LearningRole, LearningAdapter } from './learning';
export type { UIConfig, UIControlOptions, UIModelSource } from './ui';
export type { UISelection, UIStyle, UISurface, UIShape, UIVolume, UIDistance, UIViewSettings } from './state';
export type { TrajectorySettings } from './trajectory';

// Keep the optional browser bundle compatible with data-ui and existing embeds.
if (typeof window !== 'undefined' && window.$3Dmol) {
    Object.assign(window.$3Dmol, { UI, StateManager, LearningSession });
}
