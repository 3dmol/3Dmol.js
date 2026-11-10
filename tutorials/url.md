
You can view and share molecular structures with the [3Dmol.js viewer](../viewer.html), without writing HTML. Open the viewer with a structure and styles in its URL, then use the controls to edit the scene. The address bar follows the current settings, so copying the URL shares the edited view.

A 3Dmol viewer URL takes the form `https://3Dmol.org/viewer.html?[query string]`, where the `query string` specifies a structure (i.e. a PDB ID) and specific 3Dmol styles to apply.  Click the URL below for an example.

<a href="../viewer.html?pdb=1YCR&select=chain:A&style=cartoon;stick:radius~0.1&surface=opacity:0.8;colorscheme:whiteCarbon&select=chain:B&style=cartoon;line&select=resi:19,23,26;chain:B&style=cartoon;stick&labelres=backgroundOpacity:0.8;fontSize:14" style="word-wrap: break-word;">https://3Dmol.org/viewer.html?pdb=1YCR&select=chain:A&style=cartoon;stick:radius~0.1&surface=opacity:0.8;colorscheme:whiteCarbon&select=chain:B&style=cartoon;line&select=resi:19,23,26;chain:B&style=cartoon;stick&labelres=backgroundOpacity:0.8;fontSize:14</a>

### Sharing the current view

The query-string syntax below remains supported. After editing the scene, the viewer writes a complete snapshot into the URL fragment, using `#state=` followed by `encodeURIComponent(JSON.stringify(snapshot))`. When a URL includes both this fragment and query-string commands, the fragment describes the scene to restore.

Opening that link restores:

- Loaded models, their source and format, parser options, names, and visibility.
- Named selections, selection criteria, visibility, and every representation and its settings. A selection with Representation set to None is also preserved.
- Residue labels and their styling, including parameters supplied in the original URL. Static labels retain the frame at which they were created.
- Surfaces, their type and colors, the atoms shown and the atoms used for calculation, and the frame at which each surface was created.
- The camera position, clipping planes, background, projection, and trajectory frame, speed, loop, and play/pause state.

For local files, the file contents are included in the fragment; they are not uploaded to a server. Such links can be large. Use a hosted molecular-data URL for large files, and make sure it is accessible to the people receiving the link. The fragment carries the scene settings while remote structure URLs continue to refer to their original data sources.

Browser Back and Forward restore the corresponding saved scene. Moving or collapsing controls does not change the molecular settings in a link.


#### Mouse Controls ####

| Movement    | Mouse Input                                          | Touch Input          |
| ----------- | ---------------------------------------------------- | -------------------- |
| Rotation    | Primary Mouse Button                                 | Single touch         |
| Translation | Middle Mouse Button or Ctrl+Primary                  | Triple touch         |
| Zoom        | Scroll Wheel or Second Mouse Button or Shift+Primary | Pinch (double touch) |
| Slab        | Ctrl+Second                                          | Not Available        |


### Building a Query String ###

The URL query can be composed of three types of specifiers: **structure identifiers**, **atom selectors**, and **style specifications**. The viewer loads the specified structures, then applies selectors and style specifications in order.
Specifiers are read in the order they are added, and are separated by an '&' character.


#### Specifying a Structure ####

Each **structure identifier** is formatted as `pdb=[PDB ID]`, `cid=[PubChem CID]`, or `url=[URL]`
to fetch molecules from the PDB, the PubChem database, or an arbitrary URL.

The `cif=[PDB ID]` alias loads an mmCIF structure. The legacy `mmtf=[PDB ID]` alias is also accepted and resolves to the corresponding mmCIF structure from the PDB. To read an actual MMTF file, use `url=[URL]&type=mmtf`.

Repeat structure parameters to load multiple models, for example `?pdb=1YCR&url=molecule.sdf&type=sdf`. A `type` parameter applies to the most recent structure source. Models are numbered from zero in source order; use `select=model:0` to style the first model, or a JSON selection such as `{"model":[0,1]}` for several models. A selection without a model restriction can include atoms from every loaded model.

Multiframe structures load as trajectories by default. Add `frames=false` after a source to load only its first model instead, for example `?url=motion.xyz&type=xyz&frames=false&style=stick`. As with `type`, this setting applies to the most recent source.

Use `frame=N` to choose a zero-based trajectory frame. Any integer is accepted and wraps modulo the number of frames: `frame=-1` selects the last frame and `frame=-2` the second-to-last. For a trajectory with ten frames, `frame=10` selects frame zero and `frame=12` selects frame two.


Let's try viewing a structure of [green fluorescent protein](https://www.rcsb.org/structure/4KW4) (GFP).

For this example, we'll use the crystal structure with PDB ID 4KW4.

So, the structure identifier portion of the url is `pdb=4KW4`, and the url to view the structure on the 3Dmol.js server is:

[`https://3Dmol.org/viewer.html?pdb=4KW4`](../viewer.html?pdb=4KW4)

Arbitrary URLs can also be provided.

[`https://3Dmol.org/viewer.html?url=https://3dmol.org/tests/test_structs/benzene.sdf&type=sdf`](../viewer.html?url=https://3dmol.org/tests/test_structs/benzene.sdf&type=sdf)

By default, the file format will be inferred from an extension (e.g., `.sdf`) present in the URL. The file format may be
manually specified using the type parameter, as shown above.

The viewer accepts the formats supported by the core parser, including `cube`, `vasp`, and the `mcif` alias. If a URL has no extension or an unknown extension, the core attempts to recognize the format from the downloaded contents. Use `type` when that detection would be ambiguous.

Relative structure URLs are resolved against the viewer page, so a self-hosted viewer can load `?url=molecule.sdf&style=stick`. Compressed files such as `molecule.sdf.gz` retain their compression suffix when a format is inferred.

Percent-encode each parameter value that contains URL delimiters. For example, use `encodeURIComponent(structureURL)` when a data URL contains its own `&`, `=`, or `#`. Do not encode the whole query as one value. Legacy literal `+` characters are preserved; `%2B` also encodes a plus sign. Encode spaces as `%20` and a hex color's `#` as `%23`.

#### Selecting atoms ####

Whenever the viewer encounters a **style specification**, it applies the style to all of the currently selected atoms in the viewer. By default, all atoms are selected.

You can choose to apply styles to select groups of atoms by adding **atom selectors** to the URL before the **style specifications**.

Atoms can be selected based upon properties defined in the {@link AtomSelectionSpec}. They are formatted with atom properties in a semicolon-separated list and the atom property values following a colon. All specified properties must hold for an atom to be selected. For example, to select all tryptophan residues on chain B you would use `select=resn:TRP;chain:B&style=stick`:

[`https://3Dmol.org/viewer.html?pdb=1YCR&select=resn:TRP;chain:B&style=stick`](../viewer.html?pdb=1YCR&select=resn:TRP;chain:B&style=stick)

Some atom properties, such as residue ids and names, can be specified in a comma separated list.  In this case, *any* property may match to select the atom. For example, `select=resi:19,23,26;chain:B` selects the three residues 19, 23, and 26 on chain B:

[`https://3Dmol.org/viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&style=stick`](../viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&style=stick)

Selections and styles are processed in the order they are specified in the URL.  These directives can be chained together to produce complex scenes.

Use `select=all` to return to all atoms. Residue numbers also accept ranges, such as `select=resi:19-26;chain:B`. For nested selections such as `and`, `or`, `not`, or `within`, pass a JSON object as the parameter value and percent-encode it. The same JSON form is available for style specifications when the shorthand cannot express a nested object.

#### Adding style ####

After specifying a structure and (optionally) selecting a set of atoms, styles can be applied to atoms with **style specification** in the URL.  The three modes of styling available are `style`, `surface`, and `labelres`.

##### Molecular Styles #####

The available molecular styles are **line** (default), **cross**, **cartoon**, **stick**, and **sphere**. The viewer initially applies the **line** style to all atoms in the structure. Each `style` directive replaces the representations on its selected atoms; atoms outside that selection keep their existing representations.

To render GFP with a **cartoon** representation, enter

[`https://3Dmol.org/viewer.html?pdb=4KW4&style=cartoon`](../viewer.html?pdb=4KW4&style=cartoon)

You can apply multiple styles at once to the selected atoms using a semi-colon separated list:

[`https://3Dmol.org/viewer.html?pdb=4KW4&style=cartoon;stick`](../viewer.html?pdb=4KW4&style=cartoon;stick)

Characteristics of each style, such as line width and color, may also be set.  These key-value pairs take the form of a comma separated list after a colon.  Since `=` has special meaning within a URL, `~` is used to associate the key-value pairs:

[`https://3Dmol.org/viewer.html?pdb=4KW4&style=cartoon:color~spectrum;stick:radius~0.25,colorscheme~greenCarbon&select=bonds:0&style=sphere:radius~0.5`](../viewer.html?pdb=4KW4&style=cartoon:color~spectrum;stick:radius~0.25,colorscheme~greenCarbon&select=bonds:0&style=sphere:radius~0.5)


{@link AtomStyleSpec} provides for a list of possible atom style specification options.

The `addstyle` directive adds representations without clearing existing ones. Legacy `colorbyelement` and `colorbyproperty` directives are also accepted. Repeated selectors and style directives retain their order when the scene is opened and edited.

##### Adding a surface #####

A surface style specification draws the Van der Waals surface of the currently selected atoms.  It can specify the opacity and color of the surface.

[`https://3Dmol.org/viewer.html?pdb=4KW4&select=resn:HOH;invert:1&style=cartoon&surface=opacity:0.8;colorscheme:whiteCarbon`](../viewer.html?pdb=4KW4&select=resn:HOH;invert:1&style=cartoon&surface=opacity:0.8;colorscheme:whiteCarbon)

##### Labeling Residues #####

A single label per unique residue can be generated for an atom selection. The label is positioned at the center of the selected atoms within each residue and displays the residue name and number. The label can be styled according to the {@link LabelSpec}. For example, `labelres=backgroundOpacity:0.8;fontSize:14` sets the font size to 14 and makes the background slightly transparent.

[`https://3Dmol.org/viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&labelres=backgroundOpacity:0.8;fontSize:14`](../viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&labelres=backgroundOpacity:0.8;fontSize:14)

