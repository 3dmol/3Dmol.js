"""Execute the exported-grid README example using this checkout's Python wrapper."""
import json
from pathlib import Path
import re
import sys
from unittest.mock import patch

root = Path(__file__).resolve().parents[2]
readme = root / "py3Dmol" / "README.md"
match = re.search(
    r"<!-- example: exported-grid -->\s*```python\n(.*?)\n```",
    readme.read_text(encoding="utf-8"), re.DOTALL,
)
if not match:
    raise RuntimeError("The exported-grid README example is missing")

# Exercise the local source without requiring an installed Python distribution.
# Stub package metadata and make optional IPython unavailable; all exported
# HTML/JavaScript still comes from the real wrapper, with no notebook dependency.
sys.path.insert(0, str(root / "py3Dmol"))
namespace = {}
with patch("importlib.metadata.version", return_value="0+example-test"), \
        patch.dict(sys.modules, {"IPython": None}):
    exec(compile(match.group(1), str(readme), "exec"), namespace)

print(json.dumps({"html": namespace["html"], "id": namespace["view"].uniqueid}))
