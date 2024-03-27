py3Dmol
=======

A simple [IPython/Jupyter](http://jupyter.org/) widget to
embed an interactive [3Dmol.js](http://3dmol.org) viewer in a notebook.

The widget is completely static, which means the viewer doesn't need a running
IPython kernel to be useful and web pages and presentations generated from
the notebook will work as expected.  However, this also means there is only
one-way communication between the notebook and the viewer.

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

The returned view object has the exact same API as [$3Dmol.GLViewer](https://3dmol.org/doc/GLViewer.html)
with the exception that functions return None.


License
-------

MIT
