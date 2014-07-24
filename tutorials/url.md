

You can view a PDB structure immediately by visiting the WebMol.js server.  Simply type `webmol.csb.pitt.edu/viewer.html` with an appropriately formatted [URL query string](http://en.wikipedia.org/wiki/Query_string) into your browser.

A WebMol viewer URL takes the form `webmol.csb.pitt.edu/viewer.html?[query string]`, where the `query string` specifies a structure (i.e. a PDB ID) and specific WebMol styles to apply.


### Building a Query String ###

The URL query can be composed of three types of specifiers: a single **structure identifier** that is optionally followed by a number of alternating **atom selectors** and **style specifications**. The viewer will load the specified structure, and then alternately select a group of atoms and apply a specified viewing style.


Specifiers are read in the order they are added, and are separated by an '&' character.

The full url specification, then, is:

`webmol.csb.pitt.edu/viewer.html?[structure identifier[[atom selector]style specification]`

Let's work through a simple example:


#### Specifying a Structure ####

The **structure identifier** portion of the URL is a single selector formatted as `pdb=[PDB ID]`.  

Let's try viewing a structure of  [green fluorescent protein](http://www.rcsb.org/pdb/explore/explore.do?structureId=4KW4) (GFP).  

For this example, we'll use the crystal structure with PDB ID 4KW4.

So, the structure identifier portion of the url is `pdb=4KW4`, and the url to view the structure on the WebMol.js server is:

`webmol.csb.pitt.edu/viewer.html?pdb=4KW4`

Try copying and pasting this URL into your browser to view GFP!

You can click and drag to rotate the structure, and right click and drag (or use your mouse wheel) to zoom.

Of course, you can change the pdb ID in the url to view any PDB structure you wish.

Next, let's experiment with some different styles.

</br>
#### Selecting atoms and specifying style ####

##### Adding style #####

After specifying a structure, we can apply styles to specific atoms by adding a **style specification** to the URL, formatted as `style=[style spec]`

Currently, the available styles are **line** (default), **cross**, **cartoon**, **stick**, and **sphere**.  Note that the viewer initially applies the **line** style to all atoms in the structure

To render GFP with a **cartoon** representation, enter

`webmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon`

Besides specifying the style type, you can also tune various characteristics of the style, such as line width and color, by appending comma separated key-value pairs, formatted as:

`style=[stylespec],[key~val],[key~val],...`  

Check out {@link AtomStyleSpec} for a list of possible atom style specification key value pairs.

For example, in order to change the color of the cartoon representation to blue, change the style specification to `style=cartoon,color~blue`:

`webmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon,color~blue`

##### Selecting atoms #####

Whenever the viewer encounters a **style specification**, it applies the style to all of the currently selected atoms in the viewer. By default, all atoms are selected.

You can choose to apply styles to select groups of atoms by adding **atom selectors** to the URL before the **style specification**.

Atoms can be selected based upon properties defined in the {@link AtomSpec}, are formatted as key-value pairs `select=[key~value,key~value]`.

GFP has a modified residue called a *chromophore* that is nestled within the protein's beta barrel structure.  The PDB entry for GFP names this residue 'CRO'.

We can select the chromophore atoms by selecting atoms that have the *resn* property set to 'CRO': `select=resn~CRO`.  

Any subsequent style specifications will only be applied to these atoms, until another atom specification is supplied in the URL.

So, to view the chromophore as a stick structure within our cartoon representation, use:

`webmol.csb.pitt.edu/viewer.html?pdb=4KW4&style=cartoon&select=resn~CRO&style=stick`

