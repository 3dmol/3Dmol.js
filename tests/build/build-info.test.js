/** @jest-environment node */
'use strict';

const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const childProcess = require('node:child_process');
const { getBuildInfo } = require('../../scripts/build-info.cjs');

let directory;
const version = '1.2.3';
const fallback = { version, gitDescribe: null, commit: null, dirty: null };
const git = (...args) => childProcess.execFileSync('git', [
    '-c', 'user.name=Build Test', '-c', 'user.email=build@example.invalid',
    '-c', 'commit.gpgsign=false', '-c', 'core.hooksPath=/dev/null', ...args
], { cwd: directory, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }).trim();

beforeEach(() => {
    directory = fs.mkdtempSync(path.join(os.tmpdir(), '3dmol-build-info-'));
    git('init', '--quiet');
    fs.writeFileSync(path.join(directory, 'tracked.txt'), 'initial\n');
    git('add', 'tracked.txt');
    git('commit', '--quiet', '-m', 'Initial commit');
});

afterEach(() => {
    jest.restoreAllMocks();
    fs.rmSync(directory, { recursive: true, force: true });
});

test('reports an exact lightweight tag and full commit without changing the package version', () => {
    git('tag', 'v1.0.0');
    expect(getBuildInfo(directory, version)).toEqual({
        version, gitDescribe: 'v1.0.0', commit: git('rev-parse', 'HEAD'), dirty: false
    });
});

test('describes a detached commit after an annotated tag', () => {
    git('-c', 'tag.gpgsign=false', 'tag', '-a', 'v1.0.0', '-m', 'Release');
    fs.appendFileSync(path.join(directory, 'tracked.txt'), 'next\n');
    git('commit', '--quiet', '-am', 'Next commit');
    git('checkout', '--quiet', '--detach', 'HEAD');
    const info = getBuildInfo(directory, version);
    expect(info.gitDescribe).toMatch(/^v1\.0\.0-1-g[0-9a-f]+$/);
    expect(info.commit).toBe(git('rev-parse', 'HEAD'));
    expect(info.dirty).toBe(false);
});

test.each([false, true])('marks tracked changes as dirty (staged: %s)', staged => {
    git('tag', 'v1.0.0');
    fs.appendFileSync(path.join(directory, 'tracked.txt'), 'changed\n');
    if (staged) git('add', 'tracked.txt');
    expect(getBuildInfo(directory, version)).toMatchObject({
        gitDescribe: 'v1.0.0-dirty', dirty: true
    });
});

test('uses a commit abbreviation without tags and ignores untracked files', () => {
    fs.writeFileSync(path.join(directory, 'untracked.txt'), 'untracked\n');
    const info = getBuildInfo(directory, version);
    expect(info.gitDescribe).toMatch(/^[0-9a-f]+$/);
    expect(info.commit.startsWith(info.gitDescribe)).toBe(true);
    expect(info.dirty).toBe(false);
});

test('supports linked worktrees whose .git is a file', () => {
    const worktree = path.join(directory, 'linked');
    git('worktree', 'add', '--quiet', '--detach', worktree, 'HEAD');
    expect(fs.statSync(path.join(worktree, '.git')).isFile()).toBe(true);
    expect(getBuildInfo(worktree, version)).toEqual(getBuildInfo(directory, version));
});

test('does not mistake a source archive for its enclosing Git checkout', () => {
    const archive = path.join(directory, 'archive');
    fs.mkdirSync(archive);
    expect(getBuildInfo(archive, version)).toEqual(fallback);
});

test('falls back when Git is not installed or cannot read the checkout', () => {
    jest.spyOn(childProcess, 'execFileSync').mockImplementation(() => {
        throw new Error('Git unavailable');
    });
    expect(getBuildInfo(directory, version)).toEqual(fallback);
});
