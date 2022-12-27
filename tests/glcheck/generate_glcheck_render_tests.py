#!/usr/bin/env python3

#outputs all tests to /tests/glcheck/render-tests/<TestName>.html

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
generationTargetDir = pathjoin(curDir, 'render-tests')
dataSrcDir = pathjoin(testSrcDir, 'data')
testStructsSrcDir = pathjoin(projectDir, "tests", "test_structs")
testAssetsDir = pathjoin(generationTargetDir, "assets")
dataAssetDir = pathjoin(generationTargetDir, "data")
structAssetDir = pathjoin(generationTargetDir, "structs")

sys.path.append(testSrcDir) # this line lets you import python files from the /tests/auto directory
import generate_tests

def replace_refs_to_test_structs(srccode: str) -> str:
    search_term = "../test_structs"
    replace_term = "./structs"
    return srccode.replace(search_term, replace_term, -1)

def preproc_src(string: str) -> str:
    string = replace_refs_to_test_structs(string)
    #for the most part generate_tests does this, but it is necessary here
    #to handle code in divs as we aren't using the post-processed versions
    string = string.replace('viewer.render()','viewer.render(callback)')
    return string


testsys = generate_tests.TestSystem(projectDir)


if exists(generationTargetDir):
    rmtree(generationTargetDir)

mkdir(generationTargetDir)

symlink(buildSrcDir, testAssetsDir,target_is_directory=True)
symlink(dataSrcDir, dataAssetDir,target_is_directory=True)
symlink(testStructsSrcDir, structAssetDir,target_is_directory=True)

additionalCode = '';
if suffix == '-prof': #coverage doesn't play nicely with webworker
    additionalCode = '$3Dmol.setSyncSurface(true);'

for file in testsys.files:
    for example in file.examples:
        prescript = '\n'.join([f"<script>{preproc_src(pre)}</script>" for pre in example.prescripts])
        with open(pathjoin(generationTargetDir, f"{example.name}.html"), "w") as f:
            f.write(
f"""
<html lang="en">
    <head>
        <title>{example.name}</title>
        <script src="./assets/3Dmol{suffix}.js"></script>  
        <script src="https://code.jquery.com/jquery-3.6.3.min.js" integrity="sha256-pvPw+upLPUjgMXY0G+8O0xUf+/Im1MZjXxxgOcBQBXU=" crossorigin="anonymous"></script>

        <script>{additionalCode}</script>       
    </head>
    <body style="margin:0;padding:0">
        {prescript}
        {"".join([preproc_src(data) for data in example.datas])}
        {"".join([preproc_src(div) for div in example.divs])}
        {"<div id='gldiv' style='width: 400px; height: 400px; position: relative; '></div>" if not example.divs else ""}
        {"<script>var viewer=$3Dmol.createViewer($('#gldiv'));</script>" if not example.divs and not re.search(r'create(Stereo)?Viewer', example.script) else ""}
        <script>
            var callback = function(viewer) {{ 
                let v = viewer || $3Dmol.viewers[Object.keys($3Dmol.viewers)[0]];
                if(v.surfacesFinished() && !v.isAnimated() && !$3Dmol.processing_autoinit) {{
                    window.glcheck_renderDone = true;
                }} else {{
                    setTimeout(callback, 1000, v);
                }}
            }};
            {preproc_src(example.script)}
            {
                '$(window).on("load",function() {{callback();}});'
                if example.script.find("callback") == -1 and prescript.find("callback") == -1 
                else ""
            }
        </script>        
    </body>
</html>
"""
            )
