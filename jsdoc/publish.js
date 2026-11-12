// Keep project-specific documentation behavior outside the installed theme.
const { taffy } = require('@jsdoc/salty');
const upstream = require('clean-jsdoc-theme/publish');
const { prepareDoclets } = require('./customizations');
const { prepare } = require('./preserve-html');
const { finish } = require('./postprocess');
const { prepareRuntime } = require('./runtime');
const { makePortable } = require('./portable');
const { prepareTutorialExamples } = require('./tutorial-examples');

exports.publish = async function (data, opts, tutorials) {
    const helper = require('jsdoc/util/templateHelper');
    const doclets = prepareDoclets(helper.prune(data)().get());
    const legacyLinks = new Map();
    for (const doclet of doclets) {
        if (!doclet.longname) continue;
        const url = helper.createLink(doclet);
        helper.registerLink(doclet.longname, url);
        legacyLinks.set(doclet.longname, url);
    }

    const prepared = prepare(opts, tutorials);
    await upstream.publish(taffy(doclets), prepared.opts, prepared.tutorials);
    await prepared.restore(opts.destination);
    await finish(opts.destination, { basePath: opts.basePath, legacyLinks });
    const { scriptPath } = await prepareRuntime(opts.destination);
    await prepareTutorialExamples(opts.destination);
    await makePortable(opts.destination, { basePath: opts.basePath, scriptPath });
};
