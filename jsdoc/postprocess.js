'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');
const yaml = require('js-yaml');

async function outputFiles(directory, prefix = '') {
    const files = [];
    for (const entry of await fs.readdir(directory, { withFileTypes: true })) {
        const relative = path.posix.join(prefix, entry.name);
        if (entry.isDirectory()) files.push(...await outputFiles(path.join(directory, entry.name), relative));
        else if (/\.(html|md)$/.test(entry.name)) files.push(relative);
    }
    return files;
}

function attribute(tag, name) {
    return tag.match(new RegExp(`\\s${name}\\s*=\\s*(["'])([\\s\\S]*?)\\1`))?.[2];
}

function replaceAttribute(tag, name, transform) {
    return tag.replace(new RegExp(`(\\s${name}\\s*=\\s*)(["'])([\\s\\S]*?)\\2`, 'g'),
        (match, prefix, quote, value) => `${prefix}${quote}${transform(value)}${quote}`);
}

function sourcePath(url, basePath) {
    return url.startsWith(`${basePath}source/`) || url.startsWith('/source/');
}

function removeSourceLines(text, basePath) {
    text = text.replace(/<a\b[^>]*>([\s\S]*?)<\/a>/gi, (anchor) => {
        const end = anchor.indexOf('>') + 1;
        const opening = anchor.slice(0, end);
        if (!sourcePath(attribute(opening, 'href') || '', basePath)) return anchor;
        return replaceAttribute(opening, 'href', value => value.replace(/#L\d+$/, ''))
            + anchor.slice(end).replace(/:\d+(?=\s*(?:<\/[^>]+>\s*)*$)/, '');
    });
    return text.replace(/<(SourceLink|MemberMeta|MemberCard)\b[^>]*>/g, (tag, component) => {
        const href = component === 'SourceLink' ? 'href' : 'sourceHref';
        const label = component === 'SourceLink' ? 'label' : 'sourceLabel';
        if (!sourcePath(attribute(tag, href) || '', basePath)) return tag;
        tag = replaceAttribute(tag, href, value => value.replace(/#L\d+$/, ''));
        return replaceAttribute(tag, label, value => value.replace(/:\d+$/, ''));
    });
}

function repairLinks(text, aliases, basePath) {
    const resolve = value => {
        if (value.startsWith('../')) {
            const url = new URL(value, `https://docs.invalid${basePath}`);
            return url.pathname + url.search + url.hash;
        }
        const relative = value.replace(/^\.\//, '');
        const file = relative.split(/[?#]/, 1)[0];
        return aliases.has(file) ? basePath + relative : value;
    };
    // Keep executable code and comments intact, while still fixing external
    // script src attributes. Only opening tags need their URLs adjusted.
    text = text.replace(/<!--([\s\S]*?)-->|<(script|style)\b[^>]*>[\s\S]*?<\/\2\s*>|<[A-Za-z][^>]*>/gi, tag => {
        if (tag.startsWith('<!--')) return tag;
        const end = tag.indexOf('>') + 1;
        let opening = replaceAttribute(tag.slice(0, end), 'href', resolve);
        opening = replaceAttribute(opening, 'src', resolve);
        opening = replaceAttribute(opening, 'data-pdb-url', resolve);
        return opening + tag.slice(end);
    });
    // Generated Markdown retains ordinary tutorial links alongside components.
    return text.replace(/(\]\()([^\s)]+)(\))/g,
        (match, open, href, close) => `${open}${resolve(href)}${close}`);
}

function escapeHtml(value) {
    return value.replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;');
}

function outsideCodeFences(text, transform) {
    const result = [];
    let pending = '';
    let fence;
    for (const line of text.split(/(?<=\n)/)) {
        const marker = line.match(/^ {0,3}(`{3,}|~{3,})([^\r\n]*)(?:\r?\n)?$/);
        if (fence) {
            result.push(line);
            if (marker && marker[1][0] === fence[0] && marker[1].length >= fence.length
                && !marker[2].trim()) fence = undefined;
        } else if (marker && (marker[1][0] !== '`' || !marker[2].includes('`'))) {
            result.push(transform(pending), line);
            pending = '';
            fence = marker[1];
        } else {
            pending += line;
        }
    }
    result.push(transform(pending));
    return result.join('');
}

function redirectPage(target, hashRoutes = new Map()) {
    const escaped = escapeHtml(target);
    const literal = JSON.stringify(target).replace(/</g, '\\u003c');
    const routes = JSON.stringify(Object.fromEntries(hashRoutes)).replace(/</g, '\\u003c');
    // Legacy JSDoc anchors used camelCase and a leading dot/tilde for scope.
    // v5 uses lowercase headings. Unusual hand-authored anchors are best effort.
    return `<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Documentation moved</title><link rel="canonical" href="${escaped}"><meta http-equiv="refresh" content="0;url=${escaped}"><script>var hash=location.hash.slice(1);try{hash=decodeURIComponent(hash)}catch{}var routes=${routes};if(Object.prototype.hasOwnProperty.call(routes,hash)){location.replace(routes[hash]+location.search)}else{hash=hash.replace(/^[.~]/,'').replace(/_/g,'').toLowerCase();location.replace(${literal}+location.search+(hash?'#'+encodeURIComponent(hash):''))}</script></head><body><a href="${escaped}">Continue to the documentation</a></body></html>\n`;
}

/** Finish v5 output while retaining the fork's source links and old public URLs. */
async function finish(destination, { basePath = '/doc/', legacyLinks = new Map() } = {}) {
    basePath = `/${basePath.replace(/^\/+|\/+$/g, '')}/`.replace(/^\/\//, '/');
    const files = await outputFiles(destination);
    const contents = new Map(await Promise.all(files.map(async file => [
        file, await fs.readFile(path.join(destination, file), 'utf8')
    ])));
    const aliases = new Map();
    const hashRoutes = new Map();
    const tutorials = new Map();
    const validAlias = file => /^[^/\\?#]+\.html$/.test(file) && file !== 'index.html';
    const addAlias = (file, route) => {
        // Alias destinations are generated routes; legacy names must stay in
        // the output root and must never replace the actual homepage.
        if (validAlias(file)) aliases.set(file, route);
    };

    for (const [file, text] of contents) {
        if (!file.endsWith('index.md')) continue;
        const frontmatter = text.match(/^---\r?\n([\s\S]*?)\r?\n---(?:\r?\n|$)/);
        const metadata = frontmatter ? yaml.load(frontmatter[1]) : {};
        const route = file.slice(0, -'index.md'.length);
        const legacy = legacyLinks.get(metadata?.longname);
        if (legacy) {
            const file = legacy.split(/[?#]/, 1)[0];
            const hashIndex = legacy.indexOf('#');
            if (hashIndex < 0) {
                addAlias(file, route);
            } else if (validAlias(file)) {
                // Some globals and typedefs now have their own pages. They
                // must not replace their former owner's default destination.
                let hash = legacy.slice(hashIndex + 1);
                try { hash = decodeURIComponent(hash); } catch {}
                if (!hashRoutes.has(file)) hashRoutes.set(file, new Map());
                hashRoutes.get(file).set(hash, basePath + route);
            }
        }
        if (route === 'global/') addAlias('global.html', route);
        const tutorial = route.match(/^tutorials\/([^/]+)\/$/);
        if (tutorial) {
            tutorials.set(tutorial[1], route);
            addAlias(`tutorial-${tutorial[1]}.html`, route);
        }
    }
    // Preserve original tutorial spellings where v5 normalized underscores.
    for (const text of contents.values()) {
        for (const match of text.matchAll(/\btutorial-([A-Za-z0-9_-]+)\.html\b/g)) {
            const route = tutorials.get(match[1].toLowerCase().replace(/_/g, '-'));
            if (route) addAlias(match[0], route);
        }
    }
    for (const [file, text] of contents) {
        const transform = content => repairLinks(removeSourceLines(content, basePath), aliases, basePath);
        const updated = file.endsWith('.md') ? outsideCodeFences(text, transform) : transform(text);
        if (updated !== text) await fs.writeFile(path.join(destination, file), updated);
    }
    for (const [file, route] of aliases) {
        await fs.writeFile(path.join(destination, file), redirectPage(basePath + route, hashRoutes.get(file)));
    }
    return { aliases };
}

module.exports = { finish };
