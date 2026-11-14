# Contributing Guidelines to 3Dmol.js

Use Node.js 24 and npm 10 or newer for development. Run
`nvm install` and `nvm use` from the repository root to use `.nvmrc`.

The Jest and coverage commands pass `--no-sparkplug` directly to Node to
avoid the intermittent Node 24 V8 garbage-collector crash tracked in
[nodejs/node#62393](https://github.com/nodejs/node/issues/62393). Jest workers
inherit the flag; it is not allowed in `NODE_OPTIONS`. Remove this workaround
when the supported Node version includes the fix.

The security overrides in `package.json` update the modern `glob` and
`minimatch` dependencies used by Jest, coverage, and glcheck. They deliberately
leave legacy versions alone because those consumers use different APIs.
Remove the overrides when the parent packages adopt patched versions. When
changing them, check test discovery, coverage generation, and glcheck rendering
as well as the regular tests and build.

## Bug Report
 
 1. Ensure that the bug has not already been reported.
 2. If the bug is new, follow the bug report template to submit an issue.
 3. Please close your issue when you think it has been adequately resolved.

## Feature Request
 
 1. Please go through the 3Dmol.js [js documentation](https://3dmol.org/doc/index.html) to check whether the feature is already present.
 2. Ensure that the feature request is new and is not already proposed.
 3. If the feature request is new, follow the feature request template to submit an issue.
 4. Please close your issue when you think it has been adequately resolved.

## Contributions
 
 1. If you have found a patch that solves a particular issue, please submit a pull request.
 2. If you are willing to contribute to 3Dmol.js, please reach out to the community members or start by solving/ creating an issue.

 ## Doubts/Questions/Help 
 
 1. Please search the internet or go through already solved issues for a solution.
 2. If a solution is not found, feel free to reach out to the community members/ post an issue on github.
 
 Any other issues contact Professor David Koes (dkoes@pitt.edu).
 

## Documentation and examples

Keep API descriptions beside their implementation. The documentation build reads
TypeScript/JSDoc comments and tutorials; do not edit generated `doc/` pages or
maintain a second copy of method signatures in an agent guide.

| Content | Edit here |
| --- | --- |
| Viewer creation, loading, and scene operations | `src/utilities.ts`, `src/GLViewer.ts` |
| Atom selection, model representations, and parser options | `src/specs.ts`, `src/GLModel.ts`, `src/parsers/ParserOptionsSpec.ts` |
| Labels, shapes, colors, gradients, and volumetric data | Corresponding classes and specifications in `src/` |
| Optional controls and their state | `src/ui/` (separate `3Dmol.ui` bundle) |
| Getting started and hosted/embedded viewers | `doc.md`, `tutorials/` |
| Python-specific behavior | `py3Dmol/README.md`, `py3Dmol/py3Dmol/__init__.py` |
| Generation, theme behavior, and agent discovery | `jsdoc.conf.json`, `jsdoc/`, root `llms.txt` |

For public methods, document defaults, units, selection scope, return values,
side effects, and when work finishes. Distinguish drawing from waiting for a
load, frame change, or surface calculation. Link to existing specification types
instead of repeating their option lists. Identify the required bundle for APIs
outside the core viewer. Generated TypeScript declarations describe call shapes;
comments explain runtime behavior and limitations.

`npm install` runs `prepare`, which builds the development and production bundles.
After dependencies are installed, these commands have distinct purposes:

| Command | Purpose |
| --- | --- |
| `npm run build:dev` / `npm run build:prod` | Build bundles and TypeScript declarations |
| `npm run doc` | Recreate `doc/`, including HTML, Markdown, and the generated agent index; requires built bundles for the tutorial demos |
| `npm run generate:tests` | Extract API examples and manual render tests into the browser harness |
| `npm run jest -- --runInBand` | Run unit tests, excluding the generated render suite |
| `npm test` | Generate tests, run unit tests and the Python/UI browser examples, then run glcheck image comparisons |
| `npm run test:examples:browser` | Execute the Python README grid example and verify viewer/UI lifecycle in Chromium |
| `npm run cover` | Generate and run Jest render tests as well as unit tests, collecting coverage |

The root `llms.txt` is a discovery pointer. The theme generates `doc/llms.txt`
from its page inventory, linking directly to companion Markdown. Retrieve relevant
pages on demand; the build deliberately does not concatenate the whole reference.
Use documentation and declarations from the same checkout or package version.
See [the documentation adapter](jsdoc/README.md) for output validation.

Existing `@example` blocks are render-test inputs: preserve their code when making
prose-only corrections. The extractor walks `.ts` and `.js` files under `src/`
recursively and reads JSDoc attached to named declarations using the installed
TypeScript parser. It also includes `tests/auto/tests/*.js`. Python 3 and Node.js
with the development dependencies installed are required to generate tests.
`npm run generate:tests` also writes `tests/auto/examples.json` for
`generate_test.cgi`. Deploy this generated manifest with the browser harness so
the web server can serve examples using only Python; regenerate it after example
changes. See [render test notes](tests/auto/RENDER_TEST_NOTES.md) for CGI setup.

Each `@example` in a documentation block becomes an independent render test.
The first keeps its historical name, such as `GLViewer_ts_addStyle`; subsequent
examples use `_2`, `_3`, and so on. Append new examples to retain existing names
and reference images. An optional `<caption>...</caption>` describes the example
in the reference and is excluded from executable code. Examples must be standalone
JavaScript, with their own data and setup; they do not share variables with earlier
examples. Duplicate generated names fail generation rather than overwriting tests.

Markdown fences are not automatically render tests. The exported-grid example in
`py3Dmol/README.md` has an explicit browser harness that executes the marked source
block, so its documentation and test code stay together. See
[render test notes](tests/auto/RENDER_TEST_NOTES.md) for authoring and baseline review.

A new render example should use a local fixture or inline molecule, a deterministic
camera/style, and a single final `viewer.render()` after asynchronous work finishes.
The harness supplies `viewer` and rewrites the render call to signal completion.
Avoid indefinite animation and external services where local data is sufficient.
Inspect any new glcheck reference image; Jest uses mocked WebGL and cannot verify
visual correctness. See [render test notes](tests/auto/RENDER_TEST_NOTES.md).
