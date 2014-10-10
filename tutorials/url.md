

You can view a PDB structure immediately by visiting the 3Dmol.js server.  Simply type `3Dmol.csb.pitt.edu/viewer.html` with an appropriately formatted [URL query string](http://en.wikipedia.org/wiki/Query_string) into your browser.

A 3Dmol viewer URL takes the form `3Dmol.csb.pitt.edu/viewer.html?[query string]`, where the `query string` specifies a structure (i.e. a PDB ID) and specific 3Dmol styles to apply.  Click the URL below for an example.

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=1YCR&select=chain:A&style=cartoon;stick:radius~0.1&surface=opacity:0.8;colorscheme:whiteCarbon&select=chain:B&style=cartoon;line&select=resi:19,23,26;chain:B&style=stick`](../viewer.html?pdb=1YCR&select=chain:A&style=cartoon;stick:radius~0.1&surface=opacity:0.8;colorscheme:whiteCarbon&select=chain:B&style=cartoon;line&select=resi:19,23,26;chain:B&style=cartoon;stick)


### Building a Query String ###

The URL query can be composed of three types of specifiers: a single **structure identifier** that is optionally followed by a number of alternating **atom selectors** and **style specifications**. The viewer will load the specified structure, and then apply the selectors and style specifications.
Specifiers are read in the order they are added, and are separated by an '&' character.



#### Specifying a Structure ####

The **structure identifier** portion of the URL is a single selector formatted as `pdb=[PDB ID]`, `cid=[PubChem CID]`, or `url=[URL]`
to fetch molecules form the the PDB, the PubChem database, or any arbitrary URL.


Let's try viewing a structure of  [green fluorescent protein](http://www.rcsb.org/pdb/explore/explore.do?structureId=4KW4) (GFP).  

For this example, we'll use the crystal structure with PDB ID 4KW4.

So, the structure identifier portion of the url is `pdb=4KW4`, and the url to view the structure on the 3Dmol.js server is:

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=4KW4`](../viewer.html?pdb=4KW4)

Arbitrary URLs can also be provided.

[`http://3Dmol.csb.pitt.edu/viewer.html?url=http://3dmol.csb.pitt.edu/tests/test_structs/benzene.sdf&type=sdf`](../viewer.html?url=http://3dmol.csb.pitt.edu/tests/test_structs/benzene.sdf&type=sdf)

Be default, the file format will be inferred by any extension (e.g., `.sdf`) present in the URL.  The file format may be
manually specified using the type parameter, as shown above.



</br>
#### Selecting atoms and specifying style ####

##### Adding style #####

After specifying a structure, we can apply styles to specific atoms by adding a **style specification** to the URL, formatted as `style=[style spec]`


Currently, the available styles are **line** (default), **cross**, **cartoon**, **stick**, and **sphere**.  Note that the viewer initially applies the **line** style to all atoms in the structure

To render GFP with a **cartoon** representation, enter

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon`](http://3Dmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon`)

You can apply multiple styles at once to the selected atoms using a semi-colon separated list:

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon;stick`](../viewer.html?pdb=4KW4&style=cartoon;stick`)

Characteristics of each style, such as line width and color, may also be set.  These key-value pairs take the form of a comma separated list after a colon.  Since `=` has special meaning within a URL, `~` is used to associate the key-value pairs:

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon:color~spectrum;stick:radius~0.25,colorscheme~greenCarbon&select=bonds:0&style=sphere:radius~0.5`](../viewer.html?pdb=4KW4&style=cartoon:color~spectrum;stick:radius~0.25,colorscheme~greenCarbon&select=bonds:0&style=sphere:radius~0.5)


{@link AtomStyleSpec} provides for a list of possible atom style specification options.

##### Adding a surface #####

A surface style specification draws the Van der Waals surface of the currently selected atoms.  It can specify the opacity and color of the surface.

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon&surface=opacity:0.8;colorscheme:whiteCarbon`](../viewer.html?pdb=4KW4&style=cartoon&surface=opacity:0.8;colorscheme:whiteCarbon)


##### Selecting atoms #####

Whenever the viewer encounters a **style specification**, it applies the style to all of the currently selected atoms in the viewer. By default, all atoms are selected.

You can choose to apply styles to select groups of atoms by adding **atom selectors** to the URL before the **style specification**.

Atoms can be selected based upon properties defined in the {@link AtomSpec}.  They are formated like styles, with atom properties in a semi-colon separated list and the atom property values following a colon.  All specified properties must hold for an atom to be selected.  For example, to select all tryptophans on chain B:

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=1YCR&select=resn:TRP;chain:B&style=stick`](../viewer.html?pdb=1YCR&select=resn:TRP;chain:B&style=stick)

Some atom properties, such as residue ids and names, can be specified in a comma separated list.  In this case, *any* property may match to select the atom. For example, to select the three residues 19, 23, and 26:

[`http://3Dmol.csb.pitt.edu/viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&style=stick`](../viewer.html?pdb=1YCR&select=resi:19,23,26;chain:B&style=stick)

Selections and styles are processed in the order they are specified in the URL.  These directives can be chained together to produce complex scenes:


 

