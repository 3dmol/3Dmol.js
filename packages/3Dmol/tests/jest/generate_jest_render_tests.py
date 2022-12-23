#!/usr/bin/env python3

#outputs all tests to /tests/jests/render.test.js
#These use a mock WebGL instance so can't actually test correctness.

import sys, re

from os import mkdir, symlink
from os.path import join as pathjoin, exists, abspath
from posixpath import split
from shutil import copytree, rmtree

suffix = ''
if len(sys.argv) > 1:
    suffix = sys.argv[1]

curDir = abspath(split(abspath(__file__))[0]) 
projectDir = abspath(pathjoin(curDir, "..", ".."))# hardcoded directory paths :-(
testSrcDir = pathjoin(projectDir, 'tests', "auto")
buildSrcDir = pathjoin(projectDir, 'build')
dataSrcDir = pathjoin(testSrcDir, 'data')
testStructsSrcDir = pathjoin(projectDir, "tests", "test_structs")
testAssetsDir = buildSrcDir
dataAssetDir = pathjoin(testSrcDir, "data")
structAssetDir = pathjoin(projectDir, "tests","test_structs")

sys.path.append(testSrcDir) # this line lets you import python files from the /tests/auto directory
import generate_tests


def preproc_src(string: str) -> str:
    #for the most part generate_tests does this, but it is necessary here
    #to handle code in divs as we aren't using the post-processed versions
    string = string.replace("../test_structs/","3Dmol/tests/test_structs/")
    if 'mdsrv' not in string: #test93
        string = string.replace("data/","3Dmol/tests/auto/data/")
    string = string.replace('viewer.render()','viewer.render(callback)')
    return string


testsys = generate_tests.TestSystem(projectDir)


output = open(pathjoin(projectDir,'tests','jest','render.test.js'),'wt')
output.write('''/**
 * @jest-environment jest-environment-jsdom
 */
util = require('util')
global.TextEncoder = util.TextEncoder
global.TextDecoder = util.TextDecoder

jest.setTimeout(30000);

//we only use this to check the
WebGLRenderingContext.prototype.getParameter = function(p) { return "WebGL 2.0 (OpenGL ES 3.0 Chromium)"; };
WebGLRenderingContext.prototype.texImage3D = function() {};
''');
       

for file in testsys.files:
    for example in file.examples:
        prescript = '\n'.join([f"//prescript\n{preproc_src(pre)}\n" for pre in example.prescripts])

        output.write(f"test('{example.name}', (jestdone) => {{\n")

        output.write(f'''    const $3Dmol = require('../../build/3Dmol.js');
    $3Dmol.setSyncSurface(true);
    const $ = $3Dmol.$;
    console.log('{example.name}');
    ''')
        output.write(prescript);
        

        output.write(f"""document.body.innerHTML = `
            <html lang="en">
                <head>
                    <title>{example.name}</title>

                </head>
                <body style="margin:0;padding:0">
                    {"".join([preproc_src(data) for data in example.datas])}
                    {"".join([preproc_src(div) for div in example.divs])}
                    {"<div id='gldiv' style='width: 400px; height: 400px; position: relative; '></div>" if not example.divs else ""}
                </body>
            </html>
            `;
            """)

        output.write(f'''var callback = function(viewer) {{ 
                let v = viewer || $3Dmol.viewers[Object.keys($3Dmol.viewers)[0]];
                if(v.surfacesFinished() && !v.isAnimated() && !$3Dmol.processing_autoinit) {{
                    jestdone();
                }} else {{
                    setTimeout(callback, 1000, v);
                }}
            }};
            ''')

        if not example.divs and not re.search(r'create(Stereo)?Viewer', example.script):
            output.write("var viewer=$3Dmol.createViewer($('#gldiv'));\n")

        output.write(preproc_src(example.script))

        if example.script.find("callback") == -1:
            output.write('$3Dmol.autoload();\n')
            if prescript.find("callback") == -1:
                output.write('callback();\n') #don't really get functioning events..

        output.write('\n});\n\n\n')
