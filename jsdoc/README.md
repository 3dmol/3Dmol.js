# Documentation theme customizations

`jsdoc.conf.json` uses this directory as a small adapter around the published
`clean-jsdoc-theme` package. The theme itself is installed unchanged from npm.

The adapter preserves the two changes from the former 3dmol theme fork:

- Keep the last public member or method definition emitted by the TypeScript
  converter, without duplicate entries. Constants share the member group.
- Link to source files without line fragments or line-number captions, because
  conversion changes the documented source positions. Source pages remain available.

`classic.css` restores the previous theme's typography, light palette, member
colors, and compact spacing while retaining v5's responsive layout and controls.
The original Work Sans, Open Sans, and Inconsolata font files and their licenses
live in `fonts/`; JSDoc copies them into the output so the fonts also load when
viewing local documentation without an internet connection. The palette and
stylesheet are configured in `jsdoc.conf.json`.

The configured menu links appear in the desktop header. `portable.js` adds
plain links during generation; CSS hides the sidebar copy while preserving
the theme's existing mobile menu and hydration.

Version 5 converts prose HTML to Markdown, so `preserve-html.js` protects and
restores the authored viewers, scripts, styles, buttons and iframes in the
homepage and tutorials. Escaped code examples stay code. `postprocess.js`
repairs relative links for the new nested pages and emits redirects for the
previous API and tutorial URLs. Legacy member anchors are normalized to the
new lowercase headings, and symbols moved from the global page redirect to
their individual pages.

`portable.js` converts the theme's deployment paths to relative links with
explicit HTML filenames. The generated `doc/index.html` can be opened directly
or served over HTTP at any path. `runtime.js` bundles the theme controls into a
classic script and includes the search index, so they also work under `file://`
without relaxing browser security. Embedded examples may still need network
access for their external scripts and molecular data. The code tutorial's
three "Try it" examples are self-contained: `tutorial-examples.js` packages
the built viewer and the 1YCR molecule into one local script. Repeated clicks
reuse the viewer and replace the previous scene.

After updating the theme, validate both the adapter and generated pages:

```bash
npm run doc
npx jest tests/jsdoc --runInBand
node tests/jsdoc/browser.cjs
node tests/jsdoc/local-browser.cjs
```

The browser checks can take a different generated output directory as their
first argument. Both click all three tutorial buttons and inspect the rendered
scenes. The local browser check covers `file://`, an HTTP site root, and a nested
HTTP path. They use local fixtures and builds, and do not send analytics.
