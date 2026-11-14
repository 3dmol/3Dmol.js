# Render test notes for 3Dmol.js contributors

## Overview
Automated "render tests" come from every `@example` in JSDoc under `src/` and from
manual JavaScript scenarios in `tests/auto/tests`. Source discovery is recursive.
The generated HTML and JavaScript harness files are build outputs; edit the source
examples instead.

GitHub Actions executes these tests in two ways after pushes:
1. glcheck runs each test, exporting the canvas image, and comparing it to a reference image (located in tests/glcheck/reference-images)
2. jest runs the tests with coverage checking enabled.

The test files are pre-processed by python code to produce the actual test cases that will get executed by the test harnesses:
1. glcheck: `glcheck/generate_glcheck_render_tests.py` creates tests/glcheck/render-tests/*
2. jest: `jest/generate_jest_render_tests.py` creates tests/jest/render.test.js

glcheck uses Puppeteer to render the canvas, whereas jest uses jsdom.
This can lead to slight differences in test behavior.

## Running test suites
1. glcheck: `npm run glcheck`
2. jest coverage: `npm run cover`

## Running individual tests:
To run test `MYTEST`:
1. glcheck
```
npm run generate:glcheck
npx glcheck --config ./tests/glcheck/glcheck.config.json --only tests/glcheck/render-tests/MYTEST.html
```

2. jest
```
npm run generate:jest
npx jest --coverage --testNamePattern "MYTEST"
```

3. browser
If you set up a way to serve tests/auto/generate_test.cgi, you can view the render test results in the browser:
eg: https://3dmol.org/tests/auto/generate_test.cgi?test=testclick

Run `npm run generate:tests` before publishing the browser tests. It creates
`tests/auto/build.js` for the comparison harness and `tests/auto/examples.json`
for the single-example CGI. Deploy both generated files with the test pages;
they are build artifacts and are not committed. Regenerate them after editing
examples. With the manifest present, the CGI needs only Python 3, so Node and
TypeScript need not be available to the web server. A development checkout
without the manifest falls back to live extraction and requires those build
dependencies on the CGI process's PATH. Extraction errors produce an error page
and a traceback in the server log instead of a blank response.

## Multiple examples and deterministic scenes

Every `@example` in one JSDoc block is extracted independently. Keep the original
example first: its existing test name and reference image stay valid. Additional
examples receive ordinal suffixes (`_2`, `_3`, ...). Optional HTML captions appear
in documentation but are removed from the JavaScript. Only JSDoc margin stars are
removed; multiplication, exponentiation, strings, templates, and regular
expressions are preserved. The extraction regression tests run with the ordinary
Jest suite and exercise the browser, Jest, and glcheck generators.

Each scenario must initialize its own scene using the supplied `viewer`. Prefer
inline data or local fixtures to live downloads. End with one `viewer.render()`
after all required asynchronous work, using an async function or a completion
handler when necessary. Keep animations finite. Do not catch and ignore errors:
a completion flag alone cannot establish that the scene is correct.
The harness rewrites the exact call `viewer.render()` to signal completion;
use `viewer.render( )` for intermediate draws that must not finish the test.
The generators require Node and the installed TypeScript dependency in addition
to Python. Run `npm ci` before generating tests in a fresh checkout.

The manual embed scenarios can include `@data`, `@div`, and `@script` block-comment
directives to set up inline data, an embedded-viewer container, and a callback.
`embed_style_order.js` uses these to test inline XYZ loading, configuration
precedence, and lexicographic style suffixes without downloading a molecule.
The generated glcheck pages use the locally installed jQuery dependency.

For a new scenario, generate the glcheck pages and run glcheck with `--only`
pointing to that page and `--capture-missing true`. Inspect the new PNG under
`tests/glcheck/reference-images`, then rerun with capture disabled to verify it
matches. Commit the new image with the example; leave existing reference images
unchanged unless their expected behavior was deliberately changed and reviewed.

Run `npm run test:examples:browser` for the Python README grid export and real
viewer/UI lifecycle assertions. These checks run as part of `npm test`, use
local bundles and inline molecule data, and need Python 3 and Chromium through
Puppeteer. Existing surface cancellation and empty-selection unit tests provide
assertions that a single rendered image cannot establish.
