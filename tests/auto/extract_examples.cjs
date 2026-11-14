'use strict';

// TypeScript is already a build dependency. Use its syntax tree to associate
// examples with real declarations, not comment-looking text inside strings.
const ts = require('typescript');
const vm = require('node:vm');

function sourceFile(text, filename = 'example.js') {
    return ts.createSourceFile(filename, text, ts.ScriptTarget.Latest, true);
}

function literalRanges(source) {
    const ranges = [];
    function visit(node) {
        if (ts.isStringLiteralLike(node) || ts.isTemplateExpression(node)
            || node.kind === ts.SyntaxKind.RegularExpressionLiteral) {
            ranges.push([node.getStart(source), node.end]);
            return;
        }
        ts.forEachChild(node, visit);
    }
    visit(source);
    return ranges;
}

function inside(ranges, position) {
    return ranges.some(([start, end]) => position >= start && position < end);
}

function commentTags(comment) {
    const content = comment.slice(3, -2);
    // Mask one conventional margin star without shifting source positions.
    const searchable = content.replace(/^[ \t]*\*(?!\/)[ ]?/gm, margin => ' '.repeat(margin.length));
    const pattern = /^[ \t]*(?:\*[ \t]*)?@(\w+)\b([^\r\n]*)/gm;
    const tags = [];
    let codeStart;
    let captionEnd = -1;
    for (const tag of searchable.matchAll(pattern)) {
        if (tag.index < captionEnd) continue;
        if (codeStart !== undefined) {
            // Parse only the prefix through this candidate. An unterminated
            // string/template still covers a tag-looking line inside its body;
            // later captions cannot affect the classification of this boundary.
            const prefix = searchable.slice(codeStart, tag.index + tag[0].length);
            if (inside(literalRanges(sourceFile(prefix)), tag.index - codeStart)) continue;
        }
        tags.push(tag);
        codeStart = undefined;
        if (tag[1] === 'example') {
            codeStart = tag.index + tag[0].length - tag[2].length;
            const caption = searchable.slice(codeStart).match(/^\s*<caption>[\s\S]*?<\/caption>/);
            if (caption) codeStart = captionEnd = codeStart + caption[0].length;
        }
    }
    return { content, tags };
}

function normalizeMargins(code) {
    const lines = code.split('\n');
    const body = lines.slice(1).filter(line => !/^[ \t]*\*?[ \t]*$/.test(line));
    const decorated = body.length && body.every(line => /^[ \t]*\*/.test(line));
    const result = [];
    const literals = literalRanges(sourceFile(code));
    let offset = 0;
    for (const line of lines) {
        const lineStart = offset;
        offset += line.length + 1;
        const margin = line.match(/^[ \t]*\*( ?)(.*)$/);
        if (!margin || (!decorated && inside(literals, lineStart + line.indexOf('*')))) {
            result.push(line); continue;
        }
        const scanner = ts.createScanner(ts.ScriptTarget.Latest, true, ts.LanguageVariant.Standard, result.join('\n'));
        let last, token;
        while ((token = scanner.scan()) !== ts.SyntaxKind.EndOfFileToken) last = token;
        const canEndExpression = [ts.SyntaxKind.Identifier, ts.SyntaxKind.NumericLiteral,
            ts.SyntaxKind.StringLiteral, ts.SyntaxKind.NoSubstitutionTemplateLiteral,
            ts.SyntaxKind.CloseParenToken, ts.SyntaxKind.CloseBracketToken,
            ts.SyntaxKind.TrueKeyword, ts.SyntaxKind.FalseKeyword, ts.SyntaxKind.NullKeyword].includes(last);
        // Mixed/unstarred examples can continue an expression with * or **.
        // Do not mistake those operators for a documentation margin. Fully
        // decorated examples instead retain the second star on such lines.
        const continuation = !decorated && canEndExpression && !margin[2].trim().startsWith('//')
            && !/^(?:const|let|var|return|throw|if|else|for|while|function|class|try|catch)\b/.test(margin[2].trim());
        result.push(continuation && margin[2].trim() ? line : line.replace(/^[ \t]*\*( ?)/, ''));
    }
    return result.join('\n');
}

function extractExamples(text, filename) {
    const source = sourceFile(text, filename);
    const result = [];
    const visited = new Set();
    function visit(node) {
        for (const doc of node.jsDoc || []) {
            if (visited.has(doc.pos)) continue;
            visited.add(doc.pos);
            const { content, tags } = commentTags(text.slice(doc.pos, doc.end));
            const examples = tags.filter(tag => tag[1] === 'example');
            if (!examples.length) continue;
            const typedef = tags.find(tag => tag[1] === 'typedef');
            const typeName = typedef?.[2].trim().replace(/^\{[^}]*\}\s*/, '').split(/\s+/)[0];
            const declaration = ts.isVariableStatement(node) ? node.declarationList.declarations[0] : node;
            const name = typeName || declaration.name?.getText(source).replace(/^['"]|['"]$/g, '');
            if (!name) throw new Error(`${filename}: @example has no named declaration`);
            for (const [index, tag] of examples.entries()) {
                const start = tag.index + tag[0].length - tag[2].length;
                const next = tags[tags.indexOf(tag) + 1];
                let code = normalizeMargins(content.slice(start, next ? next.index : content.length));
                const caption = code.match(/^\s*<caption>([\s\S]*?)<\/caption>/);
                if (caption) code = code.slice(caption[0].length);
                validateExample(code, `${filename}: ${name} example ${index + 1}`);
                result.push({ name, ordinal: index + 1, code });
            }
        }
        ts.forEachChild(node, visit);
    }
    visit(source);
    return result;
}

function validateExample(code, location) {
    const source = sourceFile(code);
    if (!source.statements.length) throw new Error(`${location}: empty @example`);
    if (source.parseDiagnostics.length) {
        const message = ts.flattenDiagnosticMessageText(source.parseDiagnostics[0].messageText, ' ');
        throw new Error(`${location}: invalid JavaScript: ${message}`);
    }
    function visit(node, inFunction = false) {
        if (!inFunction && (ts.isAwaitExpression(node) || ts.isReturnStatement(node)
            || (ts.isForOfStatement(node) && node.awaitModifier))) {
            throw new Error(`${location}: wrap await/return in an async function or promise callback`);
        }
        if (ts.isImportDeclaration(node) || ts.isExportDeclaration(node)) {
            throw new Error(`${location}: examples must be browser scripts, not module declarations`);
        }
        ts.forEachChild(node, child => visit(child, inFunction || ts.isFunctionLike(node)));
    }
    visit(source);
    try { new vm.Script(code); }
    catch (error) { throw new Error(`${location}: invalid browser JavaScript: ${error.message}`); }
}

function prepareExample(code) {
    const ranges = literalRanges(sourceFile(code));
    const divs = [], datas = [], prescripts = [];
    const directives = { div: divs, data: datas, script: prescripts };
    // Find harness directives only in actual block comments. Leave ordinary
    // comments intact: removing them can merge tokens or change line breaks/ASI.
    const scanner = ts.createScanner(ts.ScriptTarget.Latest, false, ts.LanguageVariant.Standard, code);
    const removals = [];
    let token;
    while ((token = scanner.scan()) !== ts.SyntaxKind.EndOfFileToken) {
        const start = scanner.getTokenPos();
        const literal = ranges.find(([from, to]) => start >= from && start < to);
        if (literal) {
            scanner.setTextPos(literal[1]);
            continue;
        }
        if (token === ts.SyntaxKind.MultiLineCommentTrivia) {
            const match = scanner.getTokenText().match(/^\/\*\s*@(div|data|script)\s*([\s\S]*?)\*\/$/);
            if (match) {
                directives[match[1]].push(match[2]);
                removals.push([start, scanner.getTextPos()]);
            }
        }
    }
    for (const [start, end] of removals.reverse()) {
        code = code.slice(0, start) + code.slice(start, end).replace(/[^\r\n]/g, ' ') + code.slice(end);
    }
    const prepared = instrumentRenders(code);
    code = prepared.script;
    if (!prepared.hasRenderCallback && !divs.length) code += '\nviewer.render(callback);\n';
    return { script: code, divs, datas, prescripts: prescripts.map(script => instrumentRenders(script).script) };
}

function instrumentRenders(code) {
    const source = sourceFile(code);
    const insertions = [];
    let hasRenderCallback = false;
    function visit(node) {
        if (ts.isCallExpression(node) && ts.isPropertyAccessExpression(node.expression)
            && ts.isIdentifier(node.expression.expression) && node.expression.expression.text === 'viewer'
            && node.expression.name.text === 'render') {
            // Existing examples deliberately use render( ) for intermediate
            // draws and event handlers. Keep that opt-out from completion.
            if (!node.arguments.length && node.getText(source) === 'viewer.render()') {
                insertions.push(node.arguments.pos);
                hasRenderCallback = true;
            } else if (node.arguments.some(arg => /\bcallback\b/.test(arg.getText(source)))) {
                hasRenderCallback = true;
            }
        }
        ts.forEachChild(node, visit);
    }
    visit(source);
    for (const position of insertions.sort((a, b) => b - a)) code = code.slice(0, position) + 'callback' + code.slice(position);
    return { script: code, hasRenderCallback };
}

function processFiles(files) {
    return files.map(file => ({
        filename: file.filename,
        examples: (file.kind === 'builtin' ? [{ code: file.text }] : extractExamples(file.text, file.filename))
            .map(example => ({ ...example, ...prepareExample(example.code) }))
    }));
}

if (require.main === module) {
    let input = '';
    process.stdin.setEncoding('utf8');
    process.stdin.on('data', chunk => { input += chunk; });
    process.stdin.on('end', () => {
        try { process.stdout.write(JSON.stringify(processFiles(JSON.parse(input)))); }
        catch (error) { console.error(error.stack); process.exitCode = 1; }
    });
}
module.exports = { extractExamples, prepareExample, processFiles };
