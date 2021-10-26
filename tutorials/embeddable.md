<script src="../build/3Dmol-min.js"></script> 
A viewer can be automatically instantiated by simply assigning the class `viewer_3Dmoljs` to a `div`.
The viewer will be styled according to the containing `div`, so be sure to set a height and width.
The code below is all that is needed to create the displayed viewer.

```
{@lang xml} <script src="https://3Dmol.org/build/3Dmol-min.js" async></script>     
         <div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='2POR' data-backgroundcolor='0xffffff' data-style='stick' ></div>       
```

<div style="height: 500px; width: 500px; position: relative;" class='viewer_3Dmoljs' data-pdb='2POR' data-backgroundcolor='0xffffff' data-style='stick' data-ui="true"></div>       

The contents of the viewer can be set and manipulated through the use of `data-` tags, as shown above.  Supported tags include:
 - **data-pdb** The value describes a PDB ID to be loaded into the viewer.
 - **data-cid** The value describes a PubChem compound id to be loaded into the viewer.
 - **data-href** The value is a URL to a molecular data file.
 - **data-element** The value is the id of an HTML element on the current page that has molecular data embedded in it.
 - **data-type** The value is the file format (default pdb; can be pdb, sdf, xyz, mol2, or cube).  
 - **data-backgroundcolor** The background color (default black).
 - **data-backgroundalpha** The background alpha (default opaque: 1.0).
 - **data-select** The value is an {@link AtomSpec} selection specification.
 - **data-style** The value is a style specification.
 - **data-surface** A surface style specification.
 - **data-labelres** A residue label style specification.
 - **data-zoomto** An {@link AtomSpec} selection specification to zoom to.
 - **data-spin** If set will spin the model using {@link $3Dmol.GLViewer#spin}. Can specify axis and speed (e.g. data-spin='axis:z;speed:10')
 - **data-ui** If set will show the UI for the viewer.
 
 Multiple selections, styles, residue labels, and surfaces can be provided by appending a suffix
 each tag.  For example.
 
 ```
{@lang xml} 
         <div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='1YCR' data-backgroundcolor='0xffffff' 
         data-select1='chain:A' data-style1='cartoon:color=spectrum' data-surface1='opacity:.7;color:white' data-select2='chain:B' data-style2='stick'></div>       
```

    <div style="height: 400px; width: 400px; position: relative;" class='viewer_3Dmoljs' data-pdb='1YCR' data-backgroundcolor='0xffffff' 
         data-select1='chain:A' data-style1='cartoon:color=spectrum' data-surface1='opacity:.7;color:white' data-select2='chain:B' data-style2='stick'></div>  
 
 Once created, the 3Dmol viewer can be accessed using the id of the container div in `$3Dmol.viewers`.




 
