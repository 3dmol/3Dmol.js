py3Dmol
=======

A simple [IPython/Jupyter](http://jupyter.org/) widget to
embed an interactive [3Dmol.js](http://3dmol.org) viewer in a notebook.

The widget emits HTML and JavaScript. Once displayed, it remains interactive
without a running IPython kernel, but communication is one-way: Python queues
commands for the browser and cannot read JavaScript return values or receive
viewer events. Exported pages still need JavaScript and WebGL enabled, plus access
to the configured library URL and any remotely loaded molecular data.

If you experience problems, please file 
an [issue](https://github.com/3dmol/3Dmol.js/issues).


[An example notebook](http://nbviewer.jupyter.org/github/3dmol/3Dmol.js/blob/master/py3Dmol/examples.ipynb)

Installation
------------

From PyPI:

    pip install py3Dmol


Usage
-----

Open a notebook

    jupyter notebook

and issue

```Python
import py3Dmol
view = py3Dmol.view(query='pdb:1ubq')
view.setStyle({'cartoon':{'color':'spectrum'}})
view
```

API
---

Use the [$3Dmol.GLViewer reference](https://3dmol.org/doc/GLViewer.html) for
method names and argument objects. The wrapper translates positional arguments
into JavaScript using JSON serialization; use Python dictionaries for selections
and styles, and Python `True`, `False`, and `None` for their JSON equivalents.
Python functions cannot be passed as JavaScript callbacks.

Most forwarded methods queue a command and return the Python `view` for chaining,
not the JavaScript result. For example, a forwarded `addModel` call does not return
a model object, and a forwarded `selectedAtoms` call does not return atoms to
Python. `getModel()` is a wrapper-specific exception: it returns a model proxy
whose methods queue calls to the JavaScript [GLModel](https://3dmol.org/doc/GLModel.html).
Promises and browser-side values do not cross back into Python.

Display and export helpers
--------------------------

| Operation | Python helper |
| --- | --- |
| Display a new viewer | Evaluate the view as the last notebook expression, or call `show()` |
| Apply queued commands to the displayed viewer | `update()` |
| Return embeddable HTML as a Python string | `write_html()` |
| Write an HTML page to a path or writable file object | `write_html(f)`; `fullpage=True` fills the page when writing a file |
| Display a still or animated image in notebook output | `png()` or `apng(nframes)` after the viewer has been displayed |

Repeated `show()` calls create new displays; subsequent `update()` calls target the
most recently displayed instance. Image helpers publish browser output rather
than returning image bytes to Python. Notebook output must allow the embedded
JavaScript to run; if the viewer fails to load, inspect the browser console.

Use `viewergrid=(rows, columns)` to create a grid. Forwarded commands apply to
every cell unless `viewer=(row, column)` selects one, with zero-based indices.
`getModel` requires an explicit `viewer` when a grid is enabled. Constructor
option `linked=True` links mouse controls across the grid.

The `js` constructor argument selects the 3Dmol.js bundle URL. Match that bundle
to the API documentation you use; the hosted development reference can include
features absent from the pinned default bundle. The library is loaded once per
page, so changing `js` on a later view does not replace an already loaded library.


Export a grid with independent styles
------------------------------------

This example puts the same inline molecule in two cells, then styles each cell's
model separately. `addModel` and `zoomTo` target both cells; `getModel` selects one
cell. `write_html()` returns an HTML fragment that can be saved or embedded in a
page without an active notebook kernel. The fragment still loads the JavaScript
bundle from the configured URL.

<!-- example: exported-grid -->
```python
import py3Dmol

xyz = "2\nCarbon monoxide\nC 0 0 0\nO 1.2 0 0\n"
view = py3Dmol.view(width=640, height=320, viewergrid=(1, 2), linked=False)
view.addModel(xyz, "xyz")
view.getModel(viewer=(0, 0)).setStyle({}, {"sphere": {"radius": 0.45, "color": "red"}})
view.getModel(viewer=(0, 1)).setStyle({}, {"stick": {"radius": 0.15, "color": "blue"}})
view.zoomTo()
html = view.write_html()
```

Run `npm run test:examples:browser` from the repository root to exercise this exact
example with the local Python wrapper and JavaScript bundle, including both cells'
rendered colors. The check uses Python's standard library and Chromium via the
existing Puppeteer dependency; it does not require Jupyter or network downloads.


License
-------

MIT
