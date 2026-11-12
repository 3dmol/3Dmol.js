'use strict';

const { createHash } = require('node:crypto');
const { readdir, readFile, writeFile } = require('node:fs/promises');
const { extname, join } = require('node:path');
const { JSDOM } = require('jsdom');

// clean-jsdoc-theme v5 lowers prose HTML to Markdown, discarding the elements
// that power our embedded viewers and interactive tutorials. Protect authored
// HTML before that conversion, then restore it after the upstream publisher.
// Code examples remain ordinary text and still receive the theme's highlighting.
const EMBEDDED_HTML = [
    'script', 'style', 'iframe', 'div', 'button', 'input', 'form', 'textarea',
    'select', 'canvas', 'video', 'audio', 'object', 'embed', 'svg', 'details'
].join(',');
const TOKEN_PREFIX = 'THREEDMOLHTML';
const TOKEN_PATTERN = `${TOKEN_PREFIX}[a-f0-9]+END`;

function escapeHtml(text) {
    return text.replace(/[&<>"']/g, character => ({
        '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    })[character]);
}

/**
 * Prepare a JSDoc publish call without changing its options or tutorial tree.
 * Pass the returned opts/tutorials to the upstream publisher, await it, then
 * await restore(destination). URL rebasing belongs after restoration.
 */
function prepare(opts, tutorials) {
    const fragments = new Map();

    function protect(html) {
        if (typeof html !== 'string' || !html) return html;
        const document = JSDOM.fragment(html);
        const elements = [...document.querySelectorAll(EMBEDDED_HTML)].filter(element =>
            !element.closest('pre, code') && !element.parentElement?.closest(EMBEDDED_HTML)
        );
        if (elements.length === 0) return html;

        for (const element of elements) {
            const original = element.outerHTML;
            const hash = createHash('sha256').update(`${fragments.size}\0${original}`).digest('hex');
            const token = `${TOKEN_PREFIX}${hash}END`;
            const searchable = element.cloneNode(true);
            searchable.querySelectorAll('script, style').forEach(node => node.remove());
            const text = element.matches('script, style') ? '' : searchable.textContent.trim();
            fragments.set(token, { html: original, text });
            const placeholder = document.ownerDocument.createElement('p');
            placeholder.textContent = token;
            element.replaceWith(placeholder);
        }

        const container = document.ownerDocument.createElement('div');
        container.append(document);
        return container.innerHTML;
    }

    function prepareTutorial(node) {
        if (!node) return node;
        const result = {
            name: node.name,
            title: node.title,
            content: node.content,
            type: node.type,
            children: (node.children || []).map(prepareTutorial)
        };
        if (typeof node.content !== 'string' || !node.content) return result;
        // Tutorial.parse() is JSDoc's own Markdown-to-HTML parser. Using it
        // ensures fenced HTML is escaped before we inspect actual elements.
        const html = typeof node.parse === 'function' ? node.parse() : node.content;
        if (node.type === 2 && typeof node.parse !== 'function') {
            throw new Error('Cannot preserve HTML in a Markdown tutorial without Tutorial.parse().');
        }
        const protectedHtml = protect(html);
        if (protectedHtml !== html) {
            result.content = protectedHtml;
            result.type = 1; // jsdoc/tutorial.TYPES.HTML
        }
        return result;
    }

    const preparedOpts = { ...opts, readme: protect(opts.readme) };
    const preparedTutorials = prepareTutorial(tutorials);

    async function restore(destination) {
        if (fragments.size === 0) return;
        const restored = new Set();
        const tokenPattern = new RegExp(TOKEN_PATTERN, 'g');
        const paragraphPattern = new RegExp(`<p(?:\\s[^>]*)?>\\s*(${TOKEN_PATTERN})\\s*</p>`, 'g');

        function replaceText(text, replacement) {
            return text.replace(tokenPattern, token => {
                const fragment = fragments.get(token);
                return fragment ? replacement(fragment, token) : token;
            });
        }

        function cleanSearch(value) {
            if (typeof value === 'string') return replaceText(value, fragment => fragment.text);
            if (Array.isArray(value)) return value.map(cleanSearch);
            if (value && typeof value === 'object') {
                return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, cleanSearch(item)]));
            }
            return value;
        }

        async function visit(directory) {
            for (const entry of await readdir(directory, { withFileTypes: true })) {
                const file = join(directory, entry.name);
                if (entry.isDirectory()) {
                    await visit(file);
                    continue;
                }
                const extension = extname(entry.name);
                if (!['.html', '.md', '.json', '.txt'].includes(extension)) continue;
                const original = await readFile(file, 'utf8');
                if (!original.includes(TOKEN_PREFIX)) continue;
                let content;
                if (extension === '.html') {
                    // Restore only complete content paragraphs. The renderer
                    // also copies prose into metadata; those occurrences need
                    // plain text, never executable HTML in an attribute.
                    content = original.replace(paragraphPattern, (paragraph, token) => {
                        const fragment = fragments.get(token);
                        if (!fragment) return paragraph;
                        restored.add(token);
                        return fragment.html;
                    });
                    content = replaceText(content, fragment => escapeHtml(fragment.text));
                } else if (extension === '.json') {
                    content = JSON.stringify(cleanSearch(JSON.parse(original)));
                } else {
                    // Preserve the author's HTML in companion Markdown and
                    // optional llms-full.txt as well as in rendered pages.
                    content = replaceText(original, fragment => fragment.html);
                }
                await writeFile(file, content);
            }
        }

        await visit(destination);
        if (restored.size !== fragments.size) {
            throw new Error(`The documentation theme did not render ${fragments.size - restored.size} protected HTML fragment(s).`);
        }
    }

    return { opts: preparedOpts, tutorials: preparedTutorials, restore };
}

module.exports = { prepare };
