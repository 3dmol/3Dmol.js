#!/usr/bin/env python3

"""Display one render example, using build-time extraction when available."""

from html import escape
import json
import os
from pathlib import Path
import sys
import traceback
from urllib.parse import parse_qs

HERE = Path(__file__).resolve().parent


def load_examples():
    manifest = HERE / "examples.json"
    if manifest.is_file():
        data = json.loads(manifest.read_text(encoding="utf-8"))
        if data.get("version") != 1:
            raise ValueError("Unsupported render-example manifest version")
        return data["examples"]

    # Preserve live extraction for development checkouts with Node installed.
    # Published sites should deploy examples.json so CGI requests need only Python.
    import generate_tests
    testsys = generate_tests.TestSystem(HERE.parent.parent)
    return {
        example.name: {field: getattr(example, field)
                       for field in ("script", "divs", "datas", "prescripts")}
        for file in testsys.files for example in file.examples
    }


def error_page(status, title, message):
    print("Status: " + status)
    print("Content-Type: text/html; charset=utf-8\n")
    print("<!doctype html><html><head><title>%s</title></head>"
          "<body><h1>%s</h1><p>%s</p></body></html>" %
          (escape(title), escape(title), escape(message)))


def main():
    params = parse_qs(os.environ.get("QUERY_STRING", ""))
    test = params.get("test", ["test49"])[0]
    try:
        examples = load_examples()
    except Exception:
        traceback.print_exc(file=sys.stderr)
        error_page("500 Internal Server Error", "Unable to load render examples",
                   "Run npm run generate:tests in the build checkout, then deploy "
                   "tests/auto/examples.json alongside generate_test.cgi. "
                   "Live extraction without this file requires Node and the installed "
                   "TypeScript dependency in the CGI environment.")
        return

    if test not in examples:
        error_page("404 Not Found", "Unknown render example", test + " does not exist")
        return

    ex = examples[test]
    print("Content-Type: text/html; charset=utf-8\n")
    print("""<html><head>
        <title>%s</title>
        <script src="../../build/3Dmol.js"></script>
        <script src="https://code.jquery.com/jquery-3.6.3.min.js"></script>
        </head>
        <body style="margin:0;padding:0">
    """ % escape(test))

    for prescript in ex["prescripts"]:
        print("<script>\n%s</script>\n" % prescript)
    for data in ex["datas"]:
        print(data)
    if ex["divs"]:
        for div in ex["divs"]:
            print(div)
    else:
        print("<div id='gldiv' style='width: 100vw; height: 100vh; position: relative;'></div>")
        print('<script>var viewer=$3Dmol.createViewer($("#gldiv"));</script>')

    print("""<script>
        var callback = function() {};
        %s
        </script>
        </body>
    </html>""" % ex["script"])


if __name__ == "__main__":
    main()
