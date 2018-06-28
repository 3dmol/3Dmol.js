<style>
.page-title {visibility: hidden; height: 0px; width: 0px;} //hack to get rid of Index
</style>
<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script> 
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-55629183-1', 'auto');
  ga('send', 'pageview');

</script>

#3Dmol.js

<div  style="float: right; height: 250px; width: 250px; position: relative;" class='viewer_3Dmoljs' data-pdb='1UBQ' data-backgroundcolor='0xffffff' data-style='{"cartoon":{"color": "spectrum"}}'></div>  
<script>
setInterval(function() {
 if($3Dmol.viewers) if($3Dmol.viewers[0]) {
    var view = $3Dmol.viewers[0];
    view.rotate(1);
 }
}, 50);
</script>


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
 using only a declarative URL specification and our hosted viewer (see {@tutorial url}).

Viewers can be quickly embedded in any HTML document using just two lines of source code (see {@tutorial embeddable}).

####Mouse Controls####

Movement | &nbsp;&nbsp;&nbsp;&nbsp; | Mouse Input | &nbsp;&nbsp;&nbsp;&nbsp; | Touch Input
---------|-|:------------|------------
Rotation | |   Primary Mouse Button | | Single touch
Translation | |    Middle Mouse Button or Ctrl+Primary | | Triple touch
Zoom | | Scroll Wheel or Second Mouse Button or Shift+Primary | | Pinch (double touch)
Slab | | Ctrl+Second | | Not Available


##Developing with 3Dmol.js##

3Dmol.js provides a full-featured API for manipulating and styling molecular data.

###Getting the source code###
The library is available as a single minified JavaScript file (includes jQuery):

``` 
{@lang xml}<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script> 
```

<br>
An un-minified file is also provided for debugging purposes.
``` 
{@lang xml}<script src="http://3Dmol.csb.pitt.edu/build/3Dmol.js"></script> 
```

There is also an unminified version provided that is missing jQuery for use when compiling your own minified libraries.
```
{@lang xml}<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-nojquery.js"></script>
```

The files hosted by 3Dmol.csb.pitt.edu closely track the development version and so
will change frequently.  If you desire more stability you may copy the files into your
own project or, alternatively, we host release snapshots on the <a href="https://cdnjs.com/libraries/3Dmol">CloudFare</a>
content delivery network (note that you must specify the release version).
```
{@lang xml}<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.3.0/3Dmol-min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.3.0/3Dmol-nojquery.js"></script>
```

Using the CDN will likely provide the best network performance to your users, but features will
lag behind the development branch as we only plan to tag new releases a few times a year.


<br>
The full source distribution is available [from github](https://github.com/3dmol/3Dmol.js).

```
git clone https://github.com/3dmol/3Dmol.js.git
``` 
The checked out code can be built using [grunt](http://gruntjs.com/).
```
grunt build
```


Since 3Dmol.js is licensed under the permissive BSD open-source license, you are free
to copy this code and use it any any project, as long the code is properly acknowledged.

###Using the source code###

Every 3Dmol.js viewer canvas corresponds to a {@link $3Dmol.GLViewer} object. The viewer object
includes methods for setting viewer properties and for creating and manipulating molecular models, surfaces
and custom geometric shapes within the view.

A {@link $3Dmol.GLModel} object represents molecular data.  Each model object stores its own
rendering data and provides a convient way to reference a defined part of the scene.

Models are manipulated and styled using {@link AtomSpec} JavaScript objects. 

An example of a viewer that manipulates the styles of the embedded objects is shown below.  View the source code for the implementation details.

<iframe width=800, height=800 src="../tests/example.html"></iframe> 

### Citing 3Dmol.js

If this software is useful in your work, please use the following citation:

> Nicholas Rego and David Koes   
> 3Dmol.js: molecular visualization with WebGL   
> Bioinformatics (2015) 31 (8): 1322-1324 [doi:10.1093/bioinformatics/btu829](http://doi.org/10.1093/bioinformatics/btu829)

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
to integrate with <a href="http://sourceforge.net/projects/jsmol/">other cheminformatics libraries</a> to
automatically provide additional formats and analyses.

**Does 3Dmol.js work with touch devices?**

Yes, as long as they support WebGL.  For example, it runs great in Safari on an iPad running iOS 8.

##Contact

Please address any questions or concerns to [dkoes@pitt.edu](mailto:dkoes@pitt.edu).  
You may also [submit an issue](https://github.com/3dmol/3Dmol.js/issues) on github.

##Funding

3DMol.js is funded  is funded through R01GM108340 from the National Institute of General Medical Sciences. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institute of General Medical Sciences or the National Institutes of Health. 

