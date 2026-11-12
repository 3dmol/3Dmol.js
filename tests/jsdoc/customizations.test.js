/** @jest-environment node */
'use strict';

const { prepareDoclets } = require('../../jsdoc/customizations');

const member = (name, extra = {}) => ({
    kind: 'member',
    name,
    longname: `Viewer#${name}`,
    memberof: 'Viewer',
    scope: 'instance',
    ...extra
});

describe('JSDoc theme customizations', () => {
    test('retains the converted TypeScript member and method documentation', () => {
        const originalMember = member('background', { description: 'Original member' });
        const convertedMember = member('background', {
            description: 'Converted member', type: { names: ['ColorSpec'] }
        });
        const originalMethod = member('render', {
            kind: 'function', description: 'Original method'
        });
        const convertedMethod = member('render', {
            kind: 'function', description: 'Converted method',
            params: [{ name: 'callback', type: { names: ['function'] } }]
        });

        expect(prepareDoclets([
            originalMember, originalMethod, convertedMember, convertedMethod
        ])).toEqual([convertedMember, convertedMethod]);
    });

    test('keeps the last nonprivate definition and omits private-only symbols', () => {
        const publicMember = member('background');
        const publicMethod = member('render', { kind: 'function' });
        const privateCopies = [publicMember, publicMethod].map(doclet => ({
            ...doclet, access: 'private', description: 'Private duplicate'
        }));

        expect(prepareDoclets([
            publicMember, publicMethod, ...privateCopies,
            member('secret', { access: 'private' }),
            member('secretMethod', { kind: 'function', access: 'private' }),
            member('secretConstant', { kind: 'constant', access: 'private' })
        ])).toEqual([publicMember, publicMethod]);
    });

    test('deduplicates constants and members in the same rendering group', () => {
        const constant = member('BACKGROUND', { kind: 'constant' });
        const converted = member('BACKGROUND', { description: 'Converted constant' });
        expect(prepareDoclets([constant, converted])).toEqual([converted]);
        expect(prepareDoclets([converted, constant])).toEqual([constant]);
    });

    test('keeps functions and members with matching IDs independent', () => {
        const value = member('render', { id: 'render' });
        const method = member('render', { id: 'render', kind: 'function' });
        expect(prepareDoclets([value, method])).toEqual([value, method]);
    });

    test('preserves matching IDs on different owners and scopes, including globals', () => {
        const entries = [
            member('value', { id: 'value' }),
            member('value', { id: 'value', scope: 'static', longname: 'Viewer.value' }),
            member('value', { id: 'value', memberof: 'Model', longname: 'Model#value' }),
            member('value', {
                id: 'value', memberof: undefined, scope: 'global', longname: 'value'
            })
        ];
        const convertedGlobal = { ...entries[3], description: 'Converted global' };
        const result = prepareDoclets([...entries, convertedGlobal]);
        expect(result).toHaveLength(4);
        expect(result).toEqual(expect.arrayContaining([
            ...entries.slice(0, 3), convertedGlobal
        ]));
        expect(result).not.toContain(entries[3]);
    });

    test('uses full longnames when the publisher has not assigned IDs', () => {
        const entries = [
            member('value'),
            member('value', { scope: 'static', longname: 'Viewer.value' }),
            member('value', { memberof: 'Model', longname: 'Model#value' })
        ];
        const converted = { ...entries[0], description: 'Converted instance member' };
        const result = prepareDoclets([...entries, converted]);
        expect(result).toHaveLength(3);
        expect(result).toEqual(expect.arrayContaining([...entries.slice(1), converted]));
    });

    test('does not filter or deduplicate classes, typedefs, events, or other kinds', () => {
        const entries = ['class', 'typedef', 'event', 'namespace'].flatMap(kind => [
            { kind, longname: 'Example', access: 'private', description: 'First' },
            { kind, longname: 'Example', description: 'Second' }
        ]);
        expect(prepareDoclets(entries)).toEqual(entries);
    });

    test('keeps declarations without an ID or longname separate', () => {
        const entries = [
            { kind: 'member', description: 'First unnamed declaration' },
            { kind: 'member', description: 'Second unnamed declaration' }
        ];
        expect(prepareDoclets(entries)).toEqual(entries);
    });

    test('selects the last definition after the old version/since sort', () => {
        const latest = member('value', { version: '2', since: '3' });
        const previous = member('value', { version: '2', since: '1' });
        const oldest = member('value', { version: '1', since: '9' });
        expect(prepareDoclets([latest, previous, oldest])).toEqual([latest]);
    });

    test('retains the first position while replacing the definition when sorting is off', () => {
        const first = member('z');
        const middle = member('a');
        const last = { ...first, description: 'Replacement' };
        expect(prepareDoclets([first, middle, last], { sort: false })).toEqual([last, middle]);
    });

    test('leaves the input array and source metadata unchanged', () => {
        const entries = Object.freeze([
            Object.freeze(member('z', { meta: Object.freeze({ filename: 'viewer.ts', lineno: 17 }) })),
            Object.freeze(member('a', { meta: Object.freeze({ filename: 'viewer.ts', lineno: 42 }) }))
        ]);
        const before = JSON.stringify(entries);
        const result = prepareDoclets(entries);
        expect(result).toEqual([entries[1], entries[0]]);
        expect(JSON.stringify(entries)).toBe(before);
        expect(result[0]).toBe(entries[1]);
        expect(result[1].meta.lineno).toBe(17);
    });
});
