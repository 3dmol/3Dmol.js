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
to the API documentation you use. The library is loaded once per page, so changing 
`js` on a later view does not replace an already loaded library.


Editing, trajectory, and classroom controls
-------------------------------------------

These optional constructor arguments enable browser controls in notebooks and
exported HTML. All three features are disabled by default.

| Argument | Behavior |
| --- | --- |
| `edit=True` | Show the edit button for models, styles, surfaces, volumes, and measurements |
| `trajectory=True` | Enable trajectory controls for models loaded with frames, such as through `addModelsAsFrames` |
| `learning=True` | Show the classroom button and session-name dialog |
| `learning={...}` | Configure the classroom with the options below |
| `ui_js="..."` | Override the optional UI bundle URL |

All views default to the same pinned 3Dmol.js release on jsDelivr, with or without
controls. The UI bundle loads automatically when any controls are requested;
its URL is the sibling `3Dmol.ui.js` or `3Dmol.ui-min.js` of the core `js` URL.
To use another release or a hosted build, set `js` explicitly and use matching
core and UI builds. For a custom core filename, supply `ui_js` explicitly.
Use the same bundle URLs for every view on a page: the first core and UI bundles
loaded are reused by later views. Restart the notebook output page to switch
bundles. The examples below use the default released bundles.

Controls apply to one viewer. Grids with more than one cell share a container and
canvas, so they cannot use these options; create separate `view` objects for
independent controls or classrooms. A `viewergrid=(1, 1)` is supported.

<!-- example: viewer-controls -->
```python
import py3Dmol

frames = ("2\nFrame 0\nC 0 0 0\nO 1.2 0 0\n"
          "2\nFrame 1\nC 0 0 0\nO 1.5 0 0\n")
view = py3Dmol.view(width=640, height=400, edit=True, trajectory=True,
                   learning=True)
view.addModelsAsFrames(frames, "xyz")
view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.3}})
view.zoomTo()
html = view.write_html()
```

Display `view` in a notebook, or save `html` to embed the same controls in a page.
Use the classroom button to create a session or join an existing one. Pass a
dictionary such as `learning={"mode": "button", "session": "my-class"}` to
prefill the session name. Browser edits and classroom events are not returned to
Python; `write_html()` exports the commands queued in Python, not subsequent
browser edits.

For hidden activation, set `mode` to `"hidden"` and provide a nonempty, unique
`session` identifier. Give every participant the same identifier and server.
The classroom button stays hidden while waiting. Focus the instructor's viewer
and press **Ctrl+Alt+I** to start; other active viewers with that configuration
automatically join as students, including viewers opened later.

<!-- example: hidden-classroom -->
```python
import py3Dmol

xyz = "2\nCarbon monoxide\nC 0 0 0\nO 1.2 0 0\n"
view = py3Dmol.view(width=640, height=400, data=xyz, format="xyz",
                   style={"sphere": {"scale": 0.3}, "stick": {}},
                   learning={"mode": "hidden", "session": "py3dmol-bonding-class-a"})
html = view.write_html()
```

Replace the example session identifier with one unique to your class and viewer.
The instructor gets the usual editing, trajectory, and classroom controls.
Student editing and trajectory controls are disabled while following; students
can still rotate, zoom, and answer atom queries. A session ends when the
instructor leaves or closes their viewer.

The `learning` dictionary accepts the JSON options from
[`LearningConfig`](https://3dmol.org/doc/LearningConfig.html): `mode`, `session`,
`server`, `hotkey`, and `clientURL`. `server` defaults to `"https://3dmol.csb.pitt.edu"`;
set it to your session server's origin, including protocol and port when needed,
for example `"http://127.0.0.1:5072"`. `hotkey` changes the hidden activation
shortcut, and `clientURL` overrides the Socket.IO browser-client script. Python
callbacks, including `socketFactory` and adapter hooks, are not supported.

See the [learning-environment guide](https://github.com/3dmol/3Dmol.js/blob/master/tutorials/learning_environment.md)
for shared scene behavior, classroom controls, and server setup. Hidden mode
changes how controls appear; it does not authenticate instructors.


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

Run `npm run test:examples:browser` from the repository root to exercise the grid
and both controls examples above with the local Python wrapper and JavaScript
bundles. The checks use Python's standard library and Chromium via the existing
Puppeteer dependency; they do not require Jupyter or network downloads. Run
`npm run test:python` for constructor validation and HTML serialization checks.


License
-------

MIT
