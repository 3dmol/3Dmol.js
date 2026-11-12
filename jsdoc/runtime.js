'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');
const webpack = require('webpack');

// This function is serialized into the generated browser bundle. Keep its
// dependencies limited to browser globals and the supplied search entries.
// Coverage counters reference the Node process and cannot travel with it.
/* istanbul ignore next */
function bootstrap(searchIndex) {
    const node = document.querySelector('script[data-island-props]');
    if (node) {
        const props = JSON.parse(node.textContent);
        const searchIndexUrl = 'data:application/json;charset=utf-8,'
            + encodeURIComponent(JSON.stringify(searchIndex));
        for (const value of Object.values(props)) {
            if (value && typeof value === 'object' && value.searchIndexUrl) value.searchIndexUrl = searchIndexUrl;
        }
        node.textContent = JSON.stringify(props);
    }

    // The mobile drawer generates its Home link from basePath directly,
    // whereas its other links use explicit page slugs. Directory URLs do not
    // open index.html under file://, so normalize just those root links.
    const root = document.body?.getAttribute('data-doc-root');
    if (root == null) return;
    const rootUrl = new URL(root + '/', document.baseURI);
    const rootPath = rootUrl.pathname.replace(/\/$/, '');
    function normalize(node) {
        const links = [];
        if (node.matches?.('a[href]')) links.push(node);
        links.push(...(node.querySelectorAll?.('a[href]') || []));
        for (const link of links) {
            const url = new URL(link.getAttribute('href'), document.baseURI);
            if (url.protocol === rootUrl.protocol && url.host === rootUrl.host
                && url.pathname.replace(/\/$/, '') === rootPath) {
                link.setAttribute('href', root + '/index.html' + url.search + url.hash);
            }
        }
    }
    normalize(document.body);
    new MutationObserver(records => {
        for (const record of records) {
            if (record.type === 'attributes') normalize(record.target);
            else for (const node of record.addedNodes) normalize(node);
        }
    }).observe(document.body, { subtree: true, childList: true, attributes: true, attributeFilter: ['href'] });
}

/**
 * Bundle the theme's emitted island entries as a classic browser script.
 * Browsers block ES module imports from file:// pages, even when the modules
 * are adjacent to the HTML. A single IIFE preserves the upstream controls
 * without changing the installed theme or requiring browser security flags.
 * Each upstream entry hydrates only its own data-island elements.
 */
async function prepareRuntime(destination) {
    const directory = path.resolve(destination, '_islands');
    const entries = (await fs.readdir(directory))
        .filter(name => name.endsWith('.js') && !name.startsWith('chunk-'))
        .sort()
        .map(name => path.join(directory, name));
    if (entries.length === 0) throw new Error('The documentation theme emitted no island entry scripts.');

    const homepage = await fs.readFile(path.resolve(destination, 'index.html'), 'utf8');
    const propsText = homepage.match(/<script\b[^>]*\bdata-island-props\b[^>]*>([\s\S]*?)<\/script>/)?.[1];
    if (!propsText) throw new Error('The documentation homepage has no island props.');
    const props = JSON.parse(propsText);
    const indexUrl = Object.values(props).find(value => value?.searchIndexUrl)?.searchIndexUrl;
    if (!indexUrl) throw new Error('The documentation homepage has no search index URL.');
    const indexFile = path.resolve(destination, '_assets', path.basename(new URL(indexUrl, 'https://docs.invalid/').pathname));
    const searchIndex = JSON.parse(await fs.readFile(indexFile, 'utf8')).map(entry => ({
        ...entry,
        slug: entry.slug.replace(/^([^#?]*)/, route => (route ? route.replace(/\/$/, '') + '/' : '') + 'index.html')
    }));
    const entryPath = path.resolve(destination, '_assets/.runtime-entry.cjs');
    await fs.writeFile(entryPath, `(${bootstrap.toString()})(${JSON.stringify(searchIndex)});\n`
        + entries.map(entry => `require(${JSON.stringify(entry)});`).join('\n'));

    const compiler = webpack({
        mode: 'none',
        target: ['web', 'es2022'],
        entry: entryPath,
        devtool: false,
        output: {
            path: path.resolve(destination, '_assets'),
            filename: 'runtime-[contenthash:12].js',
            iife: true
        },
        optimization: {
            // The emitted islands are already minified. Concatenation shares
            // their common module without running a second minification pass.
            concatenateModules: true,
            minimize: false
        },
        performance: false,
        stats: 'errors-warnings'
    });

    let stats;
    try {
        stats = await new Promise((resolve, reject) => {
            compiler.run((error, stats) => {
                compiler.close(closeError => {
                    if (error || closeError) reject(error || closeError);
                    else if (stats.hasErrors()) reject(new Error(stats.toString({ all: false, errors: true })));
                    else resolve(stats);
                });
            });
        });
    } finally {
        await fs.rm(entryPath, { force: true });
    }
    const assets = stats.toJson({ all: false, assets: true }).assets;
    const runtime = assets.find(asset => /^runtime-[a-f0-9]+\.js$/.test(asset.name));
    if (!runtime) throw new Error('The documentation runtime bundle was not emitted.');
    return { scriptPath: `_assets/${runtime.name}` };
}

module.exports = { prepareRuntime };
