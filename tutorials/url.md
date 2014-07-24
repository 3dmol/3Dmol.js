
## Viewing Molecules with the WebMol.js Server ##

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

So, the structure identifier portion of the url is `pdb=4K4W`, and the url to view the structure on the WebMol.js server is:

`webmol.csb.pitt.edu/viewer.html?pdb=4KW4`

Try copying and pasting this URL into your browser to view GFP!

You can click and drag to rotate the structure, and right click and drag (or use your mouse wheel) to zoom.

Of course, you can change the pdb ID in the url to view any PDB structure you wish.

Next, let's experiment with some different styles.

#### Selecting atoms and specifying style ####

After specifying a structure, we can apply styles to specific atoms by adding a **style specification** to the URL, formatted as `style=[style spec]`



