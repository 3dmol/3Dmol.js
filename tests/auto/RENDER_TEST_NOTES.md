# Render test notes for 3Dmol.js contributors

## Overview
Automated "render tests" for the 3Dmol.js viewer live in tests/auto/tests.

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
