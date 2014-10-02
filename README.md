<style>
.page-title {visibility: hidden; height: 0px; width: 0px;} //hack to get rid of Index
</style>
<script src="http://3Dmol.csb.pitt.edu/release/3Dmol-min.js"></script> 
#$3Dmol.js
<div  style="float: right; height: 250px; width: 250px; position: relative;" class='3Dmoljs_viewer' data-pdb='1UBQ' data-backgroundcolor='0xffffff' data-style='{"cartoon":{"color": "spectrum"}}'></div>  
<script>
setInterval(function() {
 if($3Dmol.viewers) if($3Dmol.viewers[0]) {
    var view = $3Dmol.viewers[0];
    view.rotate(1);
 }
}, 50);
</script>
##Overview    

$3Dmol.js is an object-oriented, webGL based JavaScript library for online molecular visualization - No Java required!
With $3Dmol.js, you can add beautifully rendered molecular visualizations to your web applications.  Features include
 * support for pdb, sdf, mol2, xyz, and cube formats
 * parallelized molecular surface computation
 * sphere, stick, line, cross, cartoon, and surface styles
 * atom property based selection and styling
 * labels
 * clickable interactivity with molecular data
 * geometric shapes included spheres and arrows

## Getting Started ##

Molecular data can be shared and visualized without writing any HTML
 using only a declarative URL specification and our hosted viewer (see {@tutorial url}).

Viewers can be quickly embedded in any HTML document using just two lines of source code (see {@tutorial embeddable}).

##Developing with $3Dmol.js##

$3Dmol.js provides a full-featured API for manipulating and styling molecular data.

###Getting the source code###
The library is available as a single minified JavaScript file:

``` 
{@lang xml}<script src="http://3Dmol.csb.pitt.edu/release/3Dmol-min.js"></script> 
```

<br>
An un-minified file is also provided for debugging purposes.
``` 
{@lang xml}<script src="http://3Dmol.csb.pitt.edu/release/3Dmol.js"></script> 
```

<br>
The full source distribution is available [from github](https://github.com/dkoes/$3Dmol).

```
git clone https://github.com/dkoes/$3Dmol.git
``` 

###Using the source code###

Every $3Dmol.js viewer canvas corresponds to a {@link $3Dmol.GLViewer} object. The viewer object
includes methods for setting viewer properties and for creating and manipulating molecular models, surfaces
and custom geometric shapes within the view.

A {@link $3Dmol.GLModel} object represents molecular data.  Each model object stores its own
rendering data and provides a convient way to reference a defined part of the scene.

Models are manipulated and styled using {@link AtomSpec} JavaScript objects. 

An example of a viewer that manipulates the styles of the embedded objects is shown below.  View the source code for the implementation details.

<iframe width=800, height=800 src="http://3Dmol.csb.pitt.edu/doc/example.html"></iframe> 
