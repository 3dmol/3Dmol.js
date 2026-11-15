<script src="../build/3Dmol-min.js"></script> 
<script src="../build/3Dmol.ui-min.js"></script> 
A viewer can be automatically instantiated by simply assigning the class `viewer_3Dmoljs` to a `div`.
The viewer will be styled according to the containing `div`, so be sure to set a height and width.
The code below is all that is needed to create the displayed viewer.

```{@lang xml} 
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>     
<script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>     

<div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='2POR' data-backgroundcolor='0xffffff' data-style='stick' data-ui='true'></div>       
```

<div style="height: 500px; width: 500px; position: relative;" class='viewer_3Dmoljs' data-pdb='2POR' data-backgroundcolor='0xffffff' data-style='stick' data-ui='true'></div>       

The core script creates viewers when the document finishes loading. Load the separate UI script after the core script when using the controls shown above.

Use lowercase HTML `data-` attributes to configure the viewer. Selection and style values accept JSON or the shorthand described in [the URL tutorial](tutorial-url.html). Supported attributes include:

 - **`data-pdb`** The value describes a PDB ID to be loaded into the viewer.
 - **`data-cid`** The value describes a PubChem compound id to be loaded into the viewer.
 - **`data-href`** The value is a URL to a molecular data file (`data-url` is an alias). The data server must permit cross-origin requests if it is on another origin.
 - **`data-element`** The value is the id of an HTML element on the current page that has molecular data embedded in it.
 - **`data-type`** The file format passed to {@link GLViewer#addModel}. Without it, PDB and PubChem sources use their known formats, URLs use their extensions, and embedded data uses content detection. Specify the type explicitly for ambiguous content or URLs without a usable extension. See {@link FileFormats}.
 - **`data-backgroundcolor`** The background color (default black).
 - **`data-backgroundalpha`** The background alpha (default opaque: 1.0).
 - **`data-options`** {@link ParserOptionsSpec} options passed to each loaded model.
 - **`data-config`** {@link ViewerSpec} options passed to the viewer. Explicit `backgroundColor` and `backgroundAlpha` values here take precedence over the separate background attributes.
 - **`data-select`** An {@link AtomSelectionSpec}; omitted selects all atoms.
 - **`data-style`** An {@link AtomStyleSpec}. An explicit style replaces the selected atoms' representations. If omitted, line style is added while retaining styles supplied by the molecular file.
 - **`data-surface`** A {@link SurfaceStyleSpec} for a Van der Waals surface of the selected atoms.
 - **`data-labelres`** A {@link LabelSpec} for one label per selected residue.
 - **`data-spin`** If set will spin the model using {@link GLViewer#spin}. Can specify axis and speed (e.g. `data-spin='axis:z;speed:10'`)
 - **`data-callback`** A JavaScript function called after the models are loaded and initial styles and camera fitting are applied. Receives the {@link GLViewer}. Asynchronous surface computation may still be running.
 - **`data-ui`** Set to `true` to enable both the edit dialog (hamburger button) and trajectory controls. Requires the separate `3Dmol.ui` bundle. Set to `false` to disable both by default.
 - **`data-ui-edit`** Set to `true` or `false` to override whether the edit dialog is available.
 - **`data-ui-trajectory`** Set to `true` or `false` to override whether trajectory controls are available. The controls only appear when the viewer has multiple trajectory frames.
 - **`data-learning`** Set to `true` or `button` to add the classroom button, or `hidden` to start an instructor session with a hotkey. Requires the separate `3Dmol.ui` bundle. Omitted or `false` disables classroom support.
 - **`data-learning-session`** The classroom session identifier. Required for hidden mode; viewers using the same server and identifier can join the same classroom, including across pages.
 - **`data-learning-server`** The session server's origin, including protocol and optional port. Defaults to `https://3dmol.org`.
 - **`data-learning-hotkey`** The hidden-mode instructor shortcut. Defaults to `Ctrl+Alt+I` while the viewer is focused.
 - **`data-zoomto`** An {@link AtomSelectionSpec} to fit in the view; omitted fits all atoms.

Classroom support and editing controls are independent. For example,
`data-learning="true" data-ui="true"` enables both; hidden classroom viewers can
omit `data-ui` and receive editing controls when activated as the instructor.
See [embedded classrooms and the JavaScript UI API](tutorial-learning_environment.html#embedded-viewers)
for examples, session discovery, server setup, and cleanup.

 
Trajectory controls can be enabled independently of the edit dialog. This example
loads the repository's [aluminum trajectory](https://3dmol.csb.pitt.edu/tests/auto/data/Al.lammpstrj),
which contains five frames of 256 atoms. Use the play button or frame slider to
explore the trajectory; only the trajectory controls are enabled.

```{@lang xml}
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>

<div id="trajectory-example"
     style="height: 400px; width: 400px; max-width: 100%; position: relative;"
     class="viewer_3Dmoljs"
     data-href="https://3dmol.org/tests/auto/data/Al.lammpstrj"
     data-type="lammpstrj" data-options="multimodel:true;frames:true"
     data-backgroundcolor="white" data-style="sphere:radius=0.05,color=royalblue"
     data-ui-trajectory="true" data-ui-edit="false"
     data-callback="function(viewer) { viewer.zoom(4); viewer.render(); }"></div>
```

<div id="trajectory-example"
     style="height: 400px; width: 400px; max-width: 100%; position: relative;"
     class="viewer_3Dmoljs"
     data-href="https://3dmol.org/tests/auto/data/Al.lammpstrj"
     data-type="lammpstrj" data-options="multimodel:true;frames:true"
     data-backgroundcolor="white" data-style="sphere:radius=0.05,color=royalblue"
     data-ui-trajectory="true" data-ui-edit="false"
     data-callback="function(viewer) { viewer.zoom(4); viewer.render(); }"></div>

The callback zooms in on this file's fractional coordinates. The separate UI
bundle provides the playback controls; `data-ui-edit="false"` disables the edit dialog.

Alternatively, use `data-ui="true" data-ui-edit="false"`. The trajectory bar has a
minimize button that collapses it to a playback icon in the upper left. When the
edit dialog is disabled, that icon uses the hamburger button's position. It stays
hidden when no trajectory is loaded. Minimizing does not stop playback.

The same options are available in JavaScript:

```js
const controls = new $3Dmol.StateManager(viewer, { edit: false, trajectory: true });
controls.setControls({ edit: true, trajectory: false });
controls.ui.setTrajectoryMinimized(true);
```

Omitting `edit` and `trajectory` when constructing a `StateManager` enables both.
`setControls` changes only the options provided.

Multiple selections, styles, residue labels, and surfaces can be provided by appending a matching suffix
to each attribute. `data-select1` supplies the selection for `data-style1`,
`data-surface1`, and `data-labelres1`. Missing selections mean all atoms. Suffixed
styles are applied after the unsuffixed style, in lexicographic attribute order
(so `style10` precedes `style2`); later styles replace earlier representations on
overlapping atoms.
 
For example:
 
```{@lang xml} 
         <div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='1YCR' 
         data-select1='chain:A' data-style1='cartoon:color=spectrum' data-surface1='opacity:.7;color:white' data-select2='chain:B' data-style2='stick'></div>       
```

<div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='1YCR'   data-select1='chain:A' data-style1='cartoon:color=spectrum' data-surface1='opacity:.7;color:white' data-select2='chain:B' data-style2='stick'></div>  
 
 Once created, the viewer is available in {@link viewers} (`$3Dmol.viewers`) under the container's ID, or a numeric key if the container has no ID. The presence of a viewer in this registry does not mean a remote file has finished loading; use `data-callback` before working with its atoms.

To insert viewers after the initial document load, use {@link createViewer} and the [JavaScript API](tutorial-code.html). This gives your application ownership of the instance and its loading and cleanup lifecycle.
 

