'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');
const { JSDOM } = require('jsdom');

async function outputFiles(directory, prefix = '') {
    const result = [];
    for (const entry of await fs.readdir(directory, { withFileTypes: true })) {
        const file = path.posix.join(prefix, entry.name);
        if (entry.isDirectory()) result.push(...await outputFiles(path.join(directory, entry.name), file));
        else result.push(file);
    }
    return result;
}

/** Make the generated site work both as local files and under any HTTP path. */
async function makePortable(destination, { basePath = '/doc/', scriptPath }) {
    const files = new Set(await outputFiles(destination));
    const base = '/' + basePath.replace(/^\/+|\/+$/g, '');
    const explicitPage = value => {
        const [, pathname, suffix] = value.match(/^([^?#]*)(.*)$/);
        const index = path.posix.join(pathname, 'index.html');
        return (files.has(index) ? index : pathname) + suffix;
    };

    for (const file of files) {
        if (!file.endsWith('.html')) continue;
        const directory = path.posix.dirname(file);
        const relative = target => path.posix.relative(directory, target) || '.';
        const root = relative('.');
        const resolve = value => {
            const [, pathname, suffix] = value.match(/^([^?#]*)(.*)$/);
            if (pathname === base || pathname.startsWith(base + '/')) {
                return relative(explicitPage(pathname.slice(base.length).replace(/^\//, ''))) + suffix;
            }
            // Authored examples use sibling build/tests/tutorials directories.
            if (/^\/(build|tests|tutorials)\//.test(pathname)) {
                return relative('..' + pathname) + suffix;
            }
            return value;
        };
        const transformProps = (value, key) => {
            if (Array.isArray(value)) return value.map(item => transformProps(item));
            if (value && typeof value === 'object') {
                return Object.fromEntries(Object.entries(value).map(([name, item]) => [name, transformProps(item, name)]));
            }
            if (typeof value !== 'string') return value;
            if (key === 'basePath') return root;
            if (key === 'slug' || key === 'currentSlug') return explicitPage(value);
            return resolve(value);
        };

        const dom = new JSDOM(await fs.readFile(path.join(destination, file), 'utf8'));
        const document = dom.window.document;
        document.body.setAttribute('data-doc-root', root);
        // The header is static; leave the sidebar's DOM and props intact for
        // hydration and the mobile drawer, hiding its desktop copy with CSS.
        const menu = document.querySelector('[data-island="sidebar"] > nav > ul');
        const header = document.querySelector('header > div');
        if (menu && header) {
            const nav = document.createElement('nav');
            nav.className = 'docs-top-menu';
            nav.setAttribute('aria-label', 'Main navigation');
            for (const item of menu.querySelectorAll('a')) {
                const link = item.cloneNode(true);
                link.textContent = item.textContent.trim();
                link.removeAttribute('class');
                nav.append(link);
            }
            header.insertBefore(nav, header.lastElementChild);
            document.body.classList.add('docs-has-top-menu');
        }
        for (const node of document.querySelectorAll('[href], [src], [poster], [action], [data-pdb-url]')) {
            for (const attribute of ['href', 'src', 'poster', 'action', 'data-pdb-url']) {
                if (node.hasAttribute(attribute)) node.setAttribute(attribute, resolve(node.getAttribute(attribute)));
            }
        }
        for (const node of document.querySelectorAll('script[data-island-props]')) {
            node.textContent = JSON.stringify(transformProps(JSON.parse(node.textContent))).replace(/</g, '\\u003c');
        }
        for (const node of document.querySelectorAll('script[type="module"]')) {
            if (!node.textContent.includes("querySelectorAll('[data-island]')")) continue;
            node.removeAttribute('type');
            node.textContent = '';
            node.setAttribute('src', relative(scriptPath));
            node.setAttribute('defer', '');
        }
        const refresh = document.querySelector('meta[http-equiv="refresh"]');
        if (refresh && document.title === 'Documentation moved') {
            refresh.setAttribute('content', refresh.getAttribute('content').replace(/^(\d+;url=)(.*)$/, (_, delay, url) => delay + resolve(url)));
            // Only our generated redirect scripts contain these URL literals.
            // Leave authored example JavaScript and escaped code unchanged.
            for (const script of document.querySelectorAll('script:not([src])')) {
                script.textContent = script.textContent.replace(/"(?:[^"\\]|\\.)*"/g,
                    literal => JSON.stringify(resolve(JSON.parse(literal))).replace(/</g, '\\u003c'));
            }
        }
        await fs.writeFile(path.join(destination, file), dom.serialize());
        dom.window.close();
    }
}

module.exports = { makePortable };
