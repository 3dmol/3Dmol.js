## Overview

<style>
.page-title {visibility: hidden; height: 0px; width: 0px;} //hack to get rid of Index
</style>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<div  style="float: right; height: 250px; width: 250px; position: relative;" class='viewer_3Dmoljs' data-pdb='1UBQ' data-backgroundcolor='0xffffff' data-style='{"cartoon":{"color": "spectrum"}}'></div>
<script>
  setInterval(function() {
    if($3Dmol.viewers) if($3Dmol.viewers[0]) {
      var view = $3Dmol.viewers[0];
      view.rotate(1);
    }
  }, 50);
</script>

3Dmol.js is an object-oriented JavaScript library for molecular visualization in browsers with WebGL support. With 3Dmol.js, you can add beautifully rendered molecular visualizations to your web applications. 

Features include:

- Molecular structure and trajectory formats listed in {@link FileFormats}, including PDB, mmCIF, BinaryCIF, SDF, MOL2, and XYZ
- Volumetric data for isosurfaces and volume rendering (see {@link VolumeData})
- Parallelized molecular surface computation
- Sphere, stick, line, cross, cartoon, and surface styles
- Atom property based selection and styling
- Labels
- Clickable interactivity with molecular data
- Geometric shapes including spheres and arrows

## Getting Started

Molecular data can be shared and visualized without writing any HTML using only a declarative URL specification and our hosted viewer (see [Viewing Molecules with the 3Dmol Server](tutorial-url.html)).

Viewers can be embedded in an HTML document using a script and a sized container (see [Embedding a 3Dmol Viewer](tutorial-embeddable.html)). For application code, start with [Using the JavaScript API](tutorial-code.html). For notebooks, use the [py3Dmol wrapper](https://github.com/3dmol/3Dmol.js/tree/master/py3Dmol).

### Finding the right API

| Task | Reference |
| --- | --- |
| Create and configure a viewer | {@link createViewer}, {@link ViewerSpec}, {@link GLViewer} |
| Load molecular data or trajectories | {@link GLViewer#addModel}, {@link GLViewer#addModelsAsFrames}, {@link ParserOptionsSpec}, {@link FileFormats} |
| Select atoms and set representations | {@link AtomSelectionSpec}, {@link AtomStyleSpec} |
| Add surfaces, shapes, and labels | {@link SurfaceStyleSpec}, {@link GLShape}, {@link LabelSpec} |
| Load volumetric data | {@link VolumeData}, {@link IsoSurfaceSpec}, {@link VolumetricRendererSpec} |

The [API index for agents](llms.txt) links to Markdown generated from the same source comments as this reference. It is a navigation entry point, not a separately maintained API description. Use the documentation and TypeScript declarations from the version of 3Dmol.js you are running; published development documentation may include APIs absent from an older release.

#### Mouse Controls

| Movement    | Mouse Input                                          | Touch Input          |
| ----------- | ---------------------------------------------------- | -------------------- |
| Rotation    | Primary Mouse Button                                 | Single touch         |
| Translation | Middle Mouse Button or Ctrl+Primary                  | Triple touch         |
| Zoom        | Scroll Wheel or Second Mouse Button or Shift+Primary | Pinch (double touch) |
| Slab        | Ctrl+Second                                          | Not Available        |

## Using 3Dmol.js

3Dmol.js provides a full-featured API for manipulating and styling molecular data.

### Including it in your web project

#### Method 1: Using the latest version

The library is available as a single minified JavaScript:

```{@lang xml}
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
```

An un-minified file is also provided for debugging purposes:

```{@lang xml}
<script src="https://3Dmol.org/build/3Dmol.js"></script>
```


#### Method 2: Using a CDN

The files hosted by 3Dmol.org track development and may change without a version change in your page. For reproducible deployments, pin a release on a CDN, as below, or serve a specific build yourself. `3Dmol-min.js` is already minified; use that file rather than applying another CDN minifier such as `3Dmol-min.min.js`.

```{@lang html}
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"></script>
```

A pinned release does not automatically receive later fixes or features. Update the version deliberately and validate your application against that release.


#### Method 3: Using an ES6 style import

If you are using a bundler such as [Webpack](https://webpack.js.org/) or [Parcel](https://parceljs.org/) and your code is in [TypeScript](https://www.typescriptlang.org/) or using ES6 imports (and after maybe using [Babel](https://babeljs.io/) to transpile it), then you will want to import it like this:

First add the library to your dependencies:

```{@lang bash}
npm install 3dmol
# or yarn add 3dmol
```

Next, you need to import 3dmol.

```{@lang javascript}
import("3dmol/build/3Dmol.js").then( ($3Dmol) => {
console.log($3Dmol);
//can do things with $3Dmol here
});
```

Alternatively, a bundler can use a static import. The core viewer does not require a separately loaded jQuery global.

```{@lang javascript}
import * as $3Dmol from '3dmol/build/3Dmol.js';
```

The `3Dmol.js` bundle used above also populates `window.$3Dmol`. Native ES module builds are available as `build/3Dmol.es6.js` and `build/3Dmol.es6-min.js`; these export the core API without assigning that global. The optional UI is a separate `3Dmol.ui` bundle. Load browser-dependent code on the client, after a DOM is available; importing it during server-side rendering is not supported.

For TypeScript projects, importing the package name `3dmol` lets the compiler resolve its bundled `build/types/index.d.ts` declarations through the package metadata.

### Using the source code

The full source distribution is available [from github](https://github.com/3dmol/3Dmol.js).

```{@lang bash}
git clone https://github.com/3dmol/3Dmol.js
cd 3Dmol.js
npm install
```

Use Node.js 24 and npm 10 or newer as specified in the repository. With lifecycle scripts enabled, `npm install` builds the development and production bundles through `prepare`. `npm run build` also generates the example tests and documentation. The development bundle is `build/3Dmol.js`; see [CONTRIBUTING.md](https://github.com/3dmol/3Dmol.js/blob/master/CONTRIBUTING.md) for the development workflow.

Since 3Dmol.js is licensed under the permissive BSD open-source license, you are free
to copy this code and use it in any project, as long as the code is properly acknowledged.

A {@link GLViewer} manages one view. A viewer grid can share a canvas between multiple viewers. The viewer object
includes methods for setting viewer properties and for creating and manipulating molecular models, surfaces
and custom geometric shapes within the view.

A {@link GLModel} object represents molecular data. Each model stores its own
atoms and rendering data. {@link AtomSpec} describes an individual atom;
{@link AtomSelectionSpec} identifies atoms to operate on, and {@link AtomStyleSpec}
defines their representations. These objects serve different purposes.

After changing models or styles, call {@link GLViewer#render} to update the scene.
Use {@link GLViewer#zoomTo} to fit the camera to the selected atoms. Loading data
and computing surfaces can be asynchronous; follow the completion contract on
the method before using its result.

An example of a viewer that manipulates the styles of the embedded objects is shown below.  View the source code for the implementation details.

<iframe width=800, height=800 src="../tests/example.html"></iframe>

### Citing 3Dmol.js

If this software is useful in your work, please use the following citation:

> Nicholas Rego and David Koes
> 3Dmol.js: molecular visualization with WebGL
> Bioinformatics (2015) 31 (8): 1322-1324 [doi:10.1093/bioinformatics/btu829](http://doi.org/10.1093/bioinformatics/btu829)

## FAQ

### What are your future plans for 3Dmol.js?

Due to limited resources, our focus is on developing 3Dmol.js as a molecular viewer library, rather than a full featured cheminformatics/bioinformatics toolkit.

The hosted viewer and optional UI provide scene controls; molecular editing and chemical analysis are separate concerns from the core visualization API.

We hope others use 3Dmol.js as a building block for such web applications and look forward to collaborating with web developers to deliver the visualization functionality needed to enable them. Of course, if additional resources become available we may expand the scope of our efforts. However, the goal is to keep the core of 3Dmol.js as small as possible and focused on visualization.

### Are you planning on supporting additional file formats?

Time permitting, we will add support for additional formats as they are requested. We also hope to integrate with [other cheminformatics libraries](https://sourceforge.net/projects/jsmol/) to automatically provide additional formats and analysis.

### Does 3Dmol.js work with touch devices?

Touch controls are supported in browsers with WebGL enabled. The controls are listed above. Available graphics memory and browser limits can affect the size of structures and surfaces a device can display.

## Contact

Please address any questions or concerns to [dkoes@pitt.edu](mailto:dkoes@pitt.edu).
You may also [submit an issue](https://github.com/3dmol/3Dmol.js/issues/new/choose) on github.

## Funding

3DMol.js is funded through R35GM140753 from the National Institute of General Medical Sciences. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institute of General Medical Sciences or the National Institutes of Health.

<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-55629183-1', 'auto');
  ga('send', 'pageview');

</script>
