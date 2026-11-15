"""Execute a marked README example using this checkout's Python wrapper."""
import json
from pathlib import Path
import re
import sys
from types import SimpleNamespace
from unittest.mock import patch

root = Path(__file__).resolve().parents[2]
readme = root / "py3Dmol" / "README.md"
example = sys.argv[1] if len(sys.argv) > 1 else "exported-grid"
match = re.search(
    r"<!-- example: " + re.escape(example) + r" -->\s*```python\n(.*?)\n```",
    readme.read_text(encoding="utf-8"), re.DOTALL,
)
if not match and example != "query-update":
    raise RuntimeError(f"The {example} README example is missing")

# Exercise the local source without requiring an installed Python distribution.
# Stub package metadata and make optional IPython unavailable; all exported
# HTML/JavaScript still comes from the real wrapper, with no notebook dependency.
sys.path.insert(0, str(root / "py3Dmol"))
namespace = {}
with patch("importlib.metadata.version", return_value="0+example-test"), \
        patch.dict(sys.modules, {"IPython": None}):
    if example == "query-update":
        # Additional regression: notebook updates must wait for query downloads.
        import py3Dmol
        view = py3Dmol.view(query="url:https://example.test/atoms.xyz", edit=True)
        namespace = {"py3Dmol": py3Dmol, "view": view, "html": view.write_html()}
    else:
        exec(compile(match.group(1), str(readme), "exec"), namespace)

result = {"html": namespace["html"], "id": namespace["view"].uniqueid}
if "--update" in sys.argv:
    # Capture the real update() output without installing/running a notebook.
    wrapper = namespace["py3Dmol"]
    display = SimpleNamespace(publish_display_data=lambda data, metadata: data["text/html"])
    with patch.object(wrapper, "_has_IPython", True), \
            patch.object(wrapper, "IPython", SimpleNamespace(display=display), create=True):
        namespace["view"].addLabel("Python update", {"position": {"x": 0, "y": 0, "z": 0}})
        result["update"] = namespace["view"].update()

print(json.dumps(result))
