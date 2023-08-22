
You can view a PDB structure immediately by visiting the 3Dmol.js server.  Simply type `https://3Dmol.org/viewer.html` with an appropriately formatted [URL query string](http://en.wikipedia.org/wiki/Query_string) into your browser.

A 3Dmol viewer URL takes the form `https://3Dmol.org/viewer.html?[query string]`, where the `query string` specifies a structure (i.e. a PDB ID) and specific 3Dmol styles to apply.  Click the URL below for an example.

<a href="../viewer.html?pdb=1YCR&select=chain:A&style=cartoon;stick:radius~0.1&surface=opacity:0.8;colorscheme:whiteCarbon&select=chain:B&style=cartoon;line&select=resi:19,23,26;chain:B&style=cartoon;stick&labelres=backgroundOpacity:0.8;fontSize:14" style="word-wrap: break-word;">https://3Dmol.org/viewer.html?pdb=1YCR&select=chain:A&style=cartoon;stick:radius~0.1&surface=opacity:0.8;colorscheme:whiteCarbon&select=chain:B&style=cartoon;line&select=resi:19,23,26;chain:B&style=stick&labelres=backgroundOpacity:0.8;fontSize:14</a>


#### Mouse Controls ####

| Movement    | Mouse Input                                          | Touch Input          |
| ----------- | ---------------------------------------------------- | -------------------- |
| Rotation    | Primary Mouse Button                                 | Single touch         |
| Translation | Middle Mouse Button or Ctrl+Primary                  | Triple touch         |
| Zoom        | Scroll Wheel or Second Mouse Button or Shift+Primary | Pinch (double touch) |
| Slab        | Ctrl+Second                                          | Not Available        |


### Building a Query String ###

The URL query can be composed of three types of specifiers: a single **structure identifier** that is optionally followed by a number of alternating **atom selectors** and **style specifications**. The viewer will load the specified structure, and then apply the selectors and style specifications.
Specifiers are read in the order they are added, and are separated by an '&' character.


#### Specifying a Structure ####

The **structure identifier** portion of the URL is a single selector formatted as `pdb=[PDB ID]`, `cid=[PubChem CID]`, or `url=[URL]`
to fetch molecules form the the PDB, the PubChem database, or any arbitrary URL.


Let's try viewing a structure of  [green fluorescent protein](http://www.rcsb.org/pdb/explore/explore.do?structureId=4KW4) (GFP).  

For this example, we'll use the crystal structure with PDB ID 4KW4.

So, the structure identifier portion of the url is `pdb=4KW4`, and the url to view the structure on the 3Dmol.js server is:

[`https://3Dmol.org/viewer.html?pdb=4KW4`](../viewer.html?pdb=4KW4)

Arbitrary URLs can also be provided.

[`https://3Dmol.org/viewer.html?url=https://3dmol.org/tests/test_structs/benzene.sdf&type=sdf`](../viewer.html?url=https://3dmol.org/tests/test_structs/benzene.sdf&type=sdf)

Be default, the file format will be inferred by any extension (e.g., `.sdf`) present in the URL.  The file format may be
manually specified using the type parameter, as shown above.


</br>

#### Selecting atoms ####

Whenever the viewer encounters a **style specification**, it applies the style to all of the currently selected atoms in the viewer. By default, all atoms are selected.

You can choose to apply styles to select groups of atoms by adding **atom selectors** to the URL before the **style specifications**.

Atoms can be selected based upon properties defined in the {@link AtomSpec}.  They are formated with atom properties in a semi-colon separated list and the atom property values following a colon.  All specified properties must hold for an atom to be selected.  For example, to select all tryptophans on chain B you would use `select=resn:TRP;chain:B&style=stick`:

[`https://3Dmol.org/viewer.html?pdb=1YCR&select=resn:TRP;chain:B&style=stick`](../viewer.html?pdb=1YCR&select=resn:TRP;chain:B&style=stick)

Some atom properties, such as residue ids and names, can be specified in a comma separated list.  In this case, *any* property may match to select the atom. For example, `select=resi:19,23,26;chain:B` selects the three residues 19, 23, and 26 on chain B:

[`https://3Dmol.org/viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&style=stick`](../viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&style=stick)

Selections and styles are processed in the order they are specified in the URL.  These directives can be chained together to produce complex scenes.

#### Adding style ####

After specifying a structure and (optionally) selecting a set of atoms, styles can be applied to atoms with **style specification** in the URL.  The three modes of styling available are `style`, `surface`, and `labelres`.

##### Molecular Styles #####

The available molecular styles are **line** (default), **cross**, **cartoon**, **stick**, and **sphere**.  Note that the viewer initially applies the **line** style to all atoms in the structure

To render GFP with a **cartoon** representation, enter

[`https://3Dmol.org/viewer.html?pdb=4KW4&style=cartoon`](http://3Dmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon)

You can apply multiple styles at once to the selected atoms using a semi-colon separated list:

[`https://3Dmol.org/viewer.html?pdb=4KW4&style=cartoon;stick`](../viewer.html?pdb=4KW4&style=cartoon;stick)

Characteristics of each style, such as line width and color, may also be set.  These key-value pairs take the form of a comma separated list after a colon.  Since `=` has special meaning within a URL, `~` is used to associate the key-value pairs:

[`https://3Dmol.org/viewer.html?pdb=4KW4&style=cartoon:color~spectrum;stick:radius~0.25,colorscheme~greenCarbon&select=bonds:0&style=sphere:radius~0.5`](../viewer.html?pdb=4KW4&style=cartoon:color~spectrum;stick:radius~0.25,colorscheme~greenCarbon&select=bonds:0&style=sphere:radius~0.5)


{@link AtomStyleSpec} provides for a list of possible atom style specification options.

##### Adding a surface #####

A surface style specification draws the Van der Waals surface of the currently selected atoms.  It can specify the opacity and color of the surface.

[`https://3Dmol.org/viewer.html?pdb=4KW4&select=resn:HOH;invert:1&style=cartoon&surface=opacity:0.8;colorscheme:whiteCarbon`](../viewer.html?pdb=4KW4&select=resn:HOH;invert:1&style=cartoon&surface=opacity:0.8;colorscheme:whiteCarbon)

##### Labeling Residues #####

A single label per a unique residue can be generated for an atom selection.  The label is positioned at the center of the selected atoms within each residue and displays the residue name and number.  The label can be styled according to the {@link LabelSpec}.  For example, `labelres=backgroundOpacity:0.8;fontSize:14` sets the font size to 14 and makes the background slightly transparent.

[`https://3Dmol.org/viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&labelres=backgroundOpacity:0.8;fontSize:14`](../viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&labelres=backgroundOpacity:0.8;fontSize:14)

