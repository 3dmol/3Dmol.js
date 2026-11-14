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
 
