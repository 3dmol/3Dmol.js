#3Dmol.js

##Overview    

3Dmol.js is an object-oriented, webGL based JavaScript library for online molecular visualization - No Java required!
With 3Dmol.js, you can add beautifully rendered molecular visualizations to your web applications.  Features include
 * support for pdb, sdf, mol2, xyz, and cube formats
 * parallelized molecular surface computation
 * sphere, stick, line, cross, cartoon, and surface styles
 * atom property based selection and styling
 * labels
 * clickable interactivity with molecular data
 * geometric shapes including spheres and arrows

## Getting Started ##

Molecular data can be shared and visualized without writing any HTML
 using only a declarative URL specification and our hosted viewer (see [Viewing Molecules with the 3Dmol Server](http://3dmol.csb.pitt.edu/doc/tutorial-url.html)).

Viewers can be quickly embedded in any HTML document using just two lines of source code (see [Embedding a 3Dmol Viewer](http://3dmol.csb.pitt.edu/doc/tutorial-embeddable.html)).

####Mouse Controls####

| Movement | Mouse Input | Touch Input |
| -------- | ----------- | ----------- |
| Rotation | Primary Mouse Button | Single touch |
| Translation | Middle Mouse Button or Ctrl+Primary | Triple touch |
| Zoom | Scroll Wheel or Second Mouse Button or Shift+Primary | Pinch (double touch) |
| Slab | Ctrl+Second | Not Available |

##Developing with 3Dmol.js##

3Dmol.js provides a full-featured API for manipulating and styling molecular data.

###Getting the source code###
The library is available as a single minified JavaScript file (includes jQuery):

```xml
<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script> 
```

An un-minified file is also provided for debugging purposes.

```xml
<script src="http://3Dmol.csb.pitt.edu/build/3Dmol.js"></script> 
```

There is also an unminified version provided that is missing jQuery for use when compiling your own minified libraries.
```xml
<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-nojquery.js"></script>
```

The files hosted by 3Dmol.csb.pitt.edu closely track the development version and so
will change frequently.  If you desire more stability you may copy the files into your
own project or, alternatively, we host release snapshots on the [jsDeliver](http://www.jsdelivr.com)
content delivery network (CDN).

```xml
<script src="http://cdn.jsdelivr.net/3dmol.js/latest/3Dmol-min.js"></script>
<script src="http://cdn.jsdelivr.net/3dmol.js/latest/3Dmol-nojquery-min.js"></script>
```

Using the CDN will likely provide the best network performance to your users, but features will
lag behind the development branch as we only plan to tag new releases a few times a year.

The full source distribution is available [from github](https://github.com/dkoes/3Dmol.js).

```
git clone https://github.com/dkoes/3Dmol.js.git
``` 

Since 3Dmol.js is licensed under the permissive BSD open-source license, you are free
to copy this code and use it any any project, as long the code is properly acknowledged.

###Using the source code###

Every 3Dmol.js viewer canvas corresponds to a [$3Dmol.GLViewer](http://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html) object. The viewer object
includes methods for setting viewer properties and for creating and manipulating molecular models, surfaces
and custom geometric shapes within the view.

A [$3Dmol.GLModel](http://3dmol.csb.pitt.edu/doc/$3Dmol.GLModel.html) object represents molecular data.  Each model object stores its own
rendering data and provides a convient way to reference a defined part of the scene.

Models are manipulated and styled using [AtomSpec](http://3dmol.csb.pitt.edu/doc/types.html#AtomSpec) JavaScript objects. 

An example of a viewer that manipulates the styles of the embedded objects is shown below.  View the source code for the implementation details: http://3dmol.csb.pitt.edu/doc/example.html


##FAQ

**What are your future plans for 3Dmol.js?**

Due to limited resources, our focus is on developing 3Dmol.js as a molecular viewer library,
rather than a full featured cheminformatics/bioinformatics toolkit. 
For example, adding support for WebGL 2.0 has a higher priority than adding
editing functionality to the hosted viewer.
We hope others use 3Dmol.js as a building block for such web applications and look 
forward to collaborating with web developers to deliver the visualization functionality 
needed to enable them.  Of course, if additional resources become available we may
expand the scope of our efforts.  However, the goal is to keep the core of 
3Dmol.js as small as possible and focused on visualization.

**Are you planning on supporting additional file formats?**

Time permitting, we will add support for additional formats as they are requested.  We also hope
to integrate with [other cheminformatics libraries](http://sourceforge.net/projects/jsmol/) to
automatically provide additional formats and analyses.

**Does 3Dmol.js work with touch devices?**

Yes, as long as they support WebGL.  For example, it runs great in Safari on an iPad running iOS 8.

##Contact

Please address any questions or concerns to [dkoes@pitt.edu](mailto:dkoes+3dmol@pitt.edu).  
You may also [submit an issue](https://github.com/dkoes/3Dmol.js/issues) on github.

##Funding

3DMol.js is funded  is funded through R01GM108340 from the National Institute of General Medical Sciences. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institute of General Medical Sciences or the National Institutes of Health. 

