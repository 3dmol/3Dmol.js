const { existsSync } = require('node:fs');
const path = require('node:path');
const childProcess = require('node:child_process');

// Read metadata only from this checkout, never from a source archive's parent.
function getBuildInfo(directory, version) {
    const fallback = { version, gitDescribe: null, commit: null, dirty: null };
    // .git can be a directory or a file pointing to a linked worktree.
    if (!existsSync(path.join(directory, '.git'))) return fallback;

    const git = args => childProcess.execFileSync('git', args, {
        cwd: directory,
        encoding: 'utf8',
        stdio: ['ignore', 'pipe', 'ignore'],
        timeout: 5000,
        env: { ...process.env, GIT_OPTIONAL_LOCKS: '0' }
    }).trim();

    try {
        return {
            version,
            gitDescribe: git(['describe', '--tags', '--always', '--dirty']),
            commit: git(['rev-parse', '--verify', 'HEAD']),
            dirty: git(['status', '--porcelain', '--untracked-files=no']) !== ''
        };
    } catch {
        // Git is optional when building an npm package or source archive.
        return fallback;
    }
}

module.exports = { getBuildInfo };
