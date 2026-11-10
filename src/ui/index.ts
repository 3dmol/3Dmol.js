import { UI } from './ui';
import { StateManager } from './state';

export { UI, StateManager };
export type { UIConfig } from './ui';

// Keep the optional browser bundle compatible with data-ui and existing embeds.
if (typeof window !== 'undefined' && window.$3Dmol) {
    Object.assign(window.$3Dmol, { UI, StateManager });
}
