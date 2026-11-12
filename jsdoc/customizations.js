'use strict';

const { taffy } = require('@jsdoc/salty');

/**
 * Preserve the member/method handling from the former 3Dmol theme fork.
 *
 * better-docs emits converted TypeScript declarations after the original
 * declarations. Keep the last public definition, at the first definition's
 * position, just as the fork's container template did. Constants share the
 * member group because the old publisher changed their kind before rendering.
 *
 * The old publisher assigned IDs after sorting. At this earlier stage JSDoc
 * normally only supplies longname, which identifies the same symbol without
 * depending on a theme's anchor format. Include owner and scope so symbols on
 * different pages, and static/instance symbols, cannot overwrite one another.
 *
 * Return a new array without changing doclets, including source line metadata.
 */
function prepareDoclets(doclets, { sort = true } = {}) {
    let ordered = doclets;
    if (sort) {
        // Use the same public collection API and ordering as JSDoc 4 and the
        // former publisher, including its handling of missing sort fields.
        const collection = taffy(doclets);
        collection.sort('longname, version, since');
        ordered = collection().get();
    }

    const result = [];
    const positions = new Map();

    for (const doclet of ordered) {
        const group = doclet.kind === 'function' ? 'function'
            : ['member', 'constant'].includes(doclet.kind) ? 'member' : null;

        if (!group) {
            result.push(doclet);
            continue;
        }
        if (doclet.access === 'private') continue;

        const identity = doclet.id || doclet.longname;
        if (!identity) {
            // Do not collapse unrelated declarations that lack an identity.
            result.push(doclet);
            continue;
        }

        const key = JSON.stringify([
            group, doclet.memberof ?? null, doclet.scope ?? null, identity
        ]);
        if (positions.has(key)) {
            result[positions.get(key)] = doclet;
        } else {
            positions.set(key, result.length);
            result.push(doclet);
        }
    }

    return result;
}

module.exports = { prepareDoclets };
