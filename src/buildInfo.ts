/** Information identifying the source used to build this library. */
export interface BuildInfo {
    /** Version from package.json. */
    readonly version: string;
    /** Git tag/commit description, with -dirty for tracked changes; null without Git metadata. */
    readonly gitDescribe: string | null;
    /** Full Git commit hash, or null without Git metadata. */
    readonly commit: string | null;
    /** Whether tracked files had uncommitted changes; null without Git metadata. */
    readonly dirty: boolean | null;
}

declare const __3DMOL_BUILD_INFO__: BuildInfo;

/**
 * Build metadata, available as $3Dmol.buildInfo or an ES module named export.
 * Direct source imports without the webpack definitions report an unknown version.
 */
export const buildInfo: Readonly<BuildInfo> = Object.freeze(
    typeof __3DMOL_BUILD_INFO__ !== 'undefined' ? __3DMOL_BUILD_INFO__ : {
        version: 'unknown', gitDescribe: null, commit: null, dirty: null
    }
);
