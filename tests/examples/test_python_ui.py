"""Constructor/export checks for py3Dmol's optional browser controls.

Run with Python's standard library; neither an installed py3Dmol distribution
nor IPython is required. Browser behavior is covered by browser.cjs.
"""
from html.parser import HTMLParser
import importlib.util
import json
from pathlib import Path
import re
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "py3Dmol_under_test", ROOT / "py3Dmol" / "py3Dmol" / "__init__.py"
)
py3Dmol = importlib.util.module_from_spec(spec)
with patch("importlib.metadata.version", return_value="0+example-test"), \
        patch.dict(sys.modules, {"IPython": None}):
    spec.loader.exec_module(py3Dmol)


class ScriptTags(HTMLParser):
    def __init__(self):
        super().__init__()
        self.starts = 0
        self.ends = 0

    def handle_starttag(self, tag, attrs):
        if tag == "script":
            self.starts += 1

    def handle_endtag(self, tag):
        if tag == "script":
            self.ends += 1


def controls_config(html):
    """Read the serialized public StateManager constructor options."""
    match = re.search(r"new\s+\$3Dmol\.StateManager\([^,]+,\s*", html)
    if not match:
        return None
    return json.JSONDecoder().raw_decode(html[match.end():])[0]


def bundle_urls(html):
    """Read script requests without depending on a particular release number."""
    return [json.JSONDecoder().raw_decode(html[match.end():])[0]
            for match in re.finditer(r"\bloadScriptAsync\(\s*", html)]


class PythonUIOptions(unittest.TestCase):
    def test_plain_views_use_release_bundle_without_controls(self):
        for options in ({}, {"learning": False}, {"ui_js": "unused-ui.js"}):
            with self.subTest(options=options):
                html = py3Dmol.view(**options).write_html()
                urls = bundle_urls(html)
                self.assertEqual(len(urls), 1)
                self.assertRegex(urls[0],
                    r"^https://cdn\.jsdelivr\.net/npm/3dmol@\d+\.\d+\.\d+"
                    r"(?:-[0-9A-Za-z.-]+)?(?:\+[0-9A-Za-z.-]+)?/build/3Dmol-min\.js$")
                self.assertIsNone(controls_config(html))
                self.assertNotIn("unused-ui.js", html)

    def test_controls_are_independently_enabled_and_use_matching_release_bundles(self):
        # npm and Python can be released separately. Check matching bundles,
        # not equality with this checkout's package.json version.
        core_url, = bundle_urls(py3Dmol.view().write_html())
        ui_url = core_url.replace("/3Dmol-min.js", "/3Dmol.ui-min.js")
        for options in ({"edit": True}, {"trajectory": True},
                        {"learning": True}, {"learning": {}},
                        {"edit": True, "trajectory": True, "learning": True}):
            with self.subTest(options=options):
                html = py3Dmol.view(**options).write_html()
                self.assertEqual(bundle_urls(html), [core_url, ui_url])
                config = controls_config(html)
                self.assertIsNotNone(config)
                self.assertEqual(config["edit"], options.get("edit", False))
                self.assertEqual(config["trajectory"], options.get("trajectory", False))
                if "learning" in options:
                    self.assertEqual(config["learning"], options["learning"])

    def test_learning_configuration_survives_html_serialization(self):
        config = {
            "mode": "hidden",
            "session": 'Room "A" \\ lecture\n</script><script>unexpected()</script>',
            "server": "http://127.0.0.1:5072",
            "hotkey": "Ctrl+Alt+K",
            "clientURL": "https://example.test/socket.io.js",
        }
        html = py3Dmol.view(learning=config).write_html()
        self.assertEqual(controls_config(html)["learning"], config)
        tags = ScriptTags()
        tags.feed(html)
        self.assertEqual((tags.starts, tags.ends), (1, 1))

    def test_invalid_learning_configuration_fails_before_display(self):
        for value in (None, 0, 1, "hidden", [], ["button"]):
            with self.subTest(learning=value):
                with self.assertRaises((TypeError, ValueError)):
                    py3Dmol.view(learning=value)
        for mode in ("buttons", "HIDDEN", "", 1):
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError):
                    py3Dmol.view(learning={"mode": mode})
        for session in (None, "", " \t\n", 123):
            with self.subTest(session=session):
                with self.assertRaises((TypeError, ValueError)):
                    py3Dmol.view(learning={"mode": "hidden", "session": session})

    def test_grid_controls_require_one_cell(self):
        for options in ({"edit": True}, {"trajectory": True},
                        {"learning": True}, {"learning": {}}):
            with self.subTest(options=options):
                with self.assertRaises(ValueError):
                    py3Dmol.view(viewergrid=(1, 2), **options)
                self.assertIsNotNone(controls_config(
                    py3Dmol.view(viewergrid=(1, 1), **options).write_html()
                ))
        self.assertIsNone(controls_config(
            py3Dmol.view(viewergrid=(2, 2), ui_js="unused-ui.js").write_html()
        ))

    def test_ui_bundle_is_inferred_from_conventional_core_urls(self):
        for core, ui in (
            ("https://example.test/build/3Dmol-min.js?v=2#bundle",
             "https://example.test/build/3Dmol.ui-min.js?v=2#bundle"),
            ("/assets/3Dmol.js?version=local", "/assets/3Dmol.ui.js?version=local"),
            ("3Dmol.js", "3Dmol.ui.js"),
        ):
            with self.subTest(core=core):
                html = py3Dmol.view(js=core, edit=True).write_html()
                self.assertIn(core, html)
                self.assertIn(ui, html)

    def test_custom_core_names_require_an_explicit_ui_bundle_only_for_controls(self):
        core = "https://example.test/custom-core.js"
        ui = "https://example.test/custom-controls.js"
        self.assertIn(core, py3Dmol.view(js=core).write_html())
        with self.assertRaises(ValueError):
            py3Dmol.view(js=core, edit=True)
        html = py3Dmol.view(js=core, ui_js=ui, edit=True).write_html()
        self.assertIn(core, html)
        self.assertIn(ui, html)
        self.assertTrue(controls_config(html)["edit"])


if __name__ == "__main__":
    unittest.main()
