import { UI } from './ui';
import { StateManager } from './state';

export { UI, StateManager };
export type { UIConfig, UIControlOptions, UIModelSource } from './ui';
export type { UISelection, UIStyle, UISurface, UIViewSettings } from './state';
export type { TrajectorySettings } from './trajectory';

// Keep the optional browser bundle compatible with data-ui and existing embeds.
if (typeof window !== 'undefined' && window.$3Dmol) {
    Object.assign(window.$3Dmol, { UI, StateManager });
}
