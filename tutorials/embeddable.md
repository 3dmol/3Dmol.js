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

The contents of the viewer can be set and manipulated through the use of `data-` tags, as shown above.  Supported tags include:
 - **data-pdb** The value describes a PDB ID to be loaded into the viewer.
 - **data-cid** The value describes a PubChem compound id to be loaded into the viewer.
 - **data-href** The value is a URL to a molecular data file.
 - **data-element** The value is the id of an HTML element on the current page that has molecular data embedded in it.
 - **data-type** The value is the file format (default pdb; examples include pdb, sdf, xyz, mol2, cube, and lammpstrj).
 - **data-backgroundcolor** The background color (default black).
 - **data-backgroundalpha** The background alpha (default opaque: 1.0).
 - **data-options** Options for parser.
 - **data-select** The value is an {@link AtomSpec} selection specification.
 - **data-style** The value is a style specification.
 - **data-surface** A surface style specification.
 - **data-labelres** A residue label style specification.
 - **data-spin** If set will spin the model using {@link GLViewer#spin}. Can specify axis and speed (e.g. data-spin='axis:z;speed:10')
 - **data-callback** JavaScript function to call upon viewer creation.  Receives {@link GLViewer} object as argument.
 - **data-ui** Set to `true` to enable both the edit dialog (hamburger button) and trajectory controls. Requires the separate `3Dmol.ui` bundle. Set to `false` to disable both by default.
 - **data-ui-edit** Set to `true` or `false` to override whether the edit dialog is available.
 - **data-ui-trajectory** Set to `true` or `false` to override whether trajectory controls are available. The controls only appear when the viewer has multiple trajectory frames.
  - **data-zoomto** An {@link AtomSpec} selection specification to zoom to.

 
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

Multiple selections, styles, residue labels, and surfaces can be provided by appending a suffix
 each tag.  
 
For example:
 
```{@lang xml} 
         <div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='1YCR' 
         data-select1='chain:A' data-style1='cartoon:color=spectrum' data-surface1='opacity:.7;color:white' data-select2='chain:B' data-style2='stick'></div>       
```

<div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='1YCR'   data-select1='chain:A' data-style1='cartoon:color=spectrum' data-surface1='opacity:.7;color:white' data-select2='chain:B' data-style2='stick'></div>  
 
 Once created, the 3Dmol viewer can be accessed using the id of the container div in `$3Dmol.viewers`.
 
