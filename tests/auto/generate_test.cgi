#!/usr/bin/env python3

#output a single python test that is specified as a cgi argument

print("Content-Type: text/html")     # HTML is following
print()                               # blank line, end of headers

import generate_tests
import sys
import traceback
from urllib.parse import parse_qs

testsys = generate_tests.TestSystem('../..') # hardcoded directory paths :-(
testinfo = {}

for file in testsys.files:
    for example in file.examples:
        testinfo[example.name] = example

#get name of user specified test
try:
    query_string = sys.stdin.read()
    params = parse_qs(query_string)
    test = params.get("test", ["test49"])[0]
except Exception as e:
    print("<pre>")
    traceback.print_exc()
    print("</pre>")
    sys.exit(1)

if test not in testinfo:
    print("""<html ><head>
    <title>Error</title></head>
    
    <body>
    <h1>%s does not exist</h1>
    </body>
</html>""" % test)
else:
    ex = testinfo[test]
    print("""<html><head>
        <title>%s</title>
        <script src="../../build/3Dmol.js"></script>        
        <script src="https://code.jquery.com/jquery-3.6.3.min.js"></script>      
        </head>
        
        <body style="margin:0;padding:0">
    """ % test)
    
    if ex.prescripts:
        for prescript in ex.prescripts:
            print("<script>\n%s</script>\n" % prescript)
        
    if ex.datas:
        for d in ex.datas:
            print(d)
    
    if ex.divs:
        for div in ex.divs:
            print(div)
    else:
        print("<div id='gldiv' style='width: 100vw; height: 100vh; position: relative;'></div>")
        print('<script>var viewer=$3Dmol.createViewer($("#gldiv"));</script>')
        
    print("""<script>
        var callback = function() {};
        %s
        </script>        
        </body>
    </html>""" % (ex.script))

