<script src="../build/3Dmol-min.js"></script>
<script
  src="https://code.jquery.com/jquery-3.6.4.min.js"
  integrity="sha256-oP6HI9z1XaZNBrJURtCoUT5SUnxFr8s3BzRl+cbzUq8="
  crossorigin="anonymous"></script>
  
<style>
.mol-container {
  width: 100%;
  height: 400px;
  position: relative;
  border: 1px solid #999;
}

.align-center {
  width: 600px;
  margin: 20px auto 10px;
  text-align: center;
}
</style>

The most simple way to display and style a 3Dmol instance is by embedding 3Dmol parameters [within a URL](./tutorial-url.html) or [within HTML tags](./tutorial-embeddable.html). 

However, if you are interacting with this viewer from your own web pages then you may want to take advantage of the 3Dmol API.

Here we will use the API to create, load and style the 3Dmol instance. Note, this tutorial assumes a basic knowledge of HTML, Javascript and jQuery.

### Creating the 3DMol instance

First, make sure your HTML document links to the 3Dmol code.

```
<script src="https://3Dmol.org/build/3Dmol-min.js""></script>
```

Now create a container tag in the HTML document that will hold the vizualisation.

```
<div id="container-01" class="mol-container"></div>
```

Note: 3Dmol will adopt the size of the container so we need to make sure this size has been explicitly set (i.e. use CSS or an inline ```style``` attribute).

```
<style>
.mol-container {
  width: 60%;
  height: 400px;
  position: relative;
}
</style>
```
 
Once the web page has loaded, we can create a new instance of 3Dmol.

To make sure everything works, we are going to add a sphere, set the camera, render the scene, then add a zoom.

```
<script>
$(function() {
  let element = $('#container-01');
  let config = { backgroundColor: 'orange' };
  let viewer = $3Dmol.createViewer( element, config );
  viewer.addSphere({ center: {x:0, y:0, z:0}, radius: 10.0, color: 'green' });
  viewer.zoomTo();
  viewer.render();
  viewer.zoom(0.8, 2000);
});
</script>
```

Documentation: [$3Dmol.createViewer()]($3Dmol.html#createViewer)

<div class="align-center">
  <button id="btn-01" class="btn btn-primary">Try it</button>
  <div id="container-01" class="mol-container"></div>
</div>

<script>
jQuery(function() {
  let viewer = '';
  $('#btn-01').on('click', function() {
    let element = $('#container-01');
    let config = { backgroundColor: 'orange' };
    
    viewer = $3Dmol.createViewer( element, config );
    viewer.addSphere({ center: {x:0, y:0, z:0}, radius: 10.0, color: 'green' });
    viewer.zoomTo();
    viewer.render();
    viewer.zoom(0.8, 2000);
  });
  
  $('#btn-01-alt').on('click', function() {
    viewer.setBackgroundColor('white');
  });
});
</script>

If this has worked, you should see a rather fetching green ball in front of an orange background. If not, then now's a good time to get familiar with the developer console on your favourite browser and check for typos ([Firefox](https://developer.mozilla.org/en/docs/Tools/Web_Console), [Chrome](https://developers.google.com/web/tools/chrome-devtools/), ...)

Note, the ```viewer``` variable now contains an instance of GLViewer and we use the [$3Dmol.GLViewer API]($3Dmol.GLViewer) to change the orange background white.

```
  viewer.setBackgroundColor('white');
```

<div class="align-center">
  <button id="btn-01-alt" class="btn btn-primary">Try it</button>
</div>

Documentation: [$3Dmol.GLViewer API]($3Dmol.GLViewer)

### Loading data dynamically

Since jQuery is already loaded, we'll use ```jQuery.ajax()``` to load PDB data from an external source. If the operation is successful then we can feed the raw PDB data into our viewer, otherwise we can report the problem.

```
  let viewer = $3Dmol.createViewer( element, config );
  let pdbUri = '/path/to/your/pdb/files/1ycr.pdb';
  jQuery.ajax( pdbUri, { 
    success: function(data) {
      let v = viewer;
      v.addModel( data, "pdb" );                       /* load data */
      v.setStyle({}, {cartoon: {color: 'spectrum'}});  /* style all atoms */
      v.zoomTo();                                      /* set camera */
      v.render();                                      /* render scene */
      v.zoom(1.2, 1000);                               /* slight zoom */
    },
    error: function(hdr, status, err) {
      console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
  });
```

All these methods (and more) are covered in detail in the [API documentation for GLViewer]($3Dmol.GLViewer).

<div class="align-center">
  <button class="btn btn-primary" id="btn-02">Try it</button>
  <div id="container-02" class="mol-container"></div>
</div>

<script>
jQuery(function() {
  let element = $('#container-02');
  let config = { defaultcolors: $3Dmol.rasmolElementColors, backgroundColor: 'white' };
  let viewer = $3Dmol.createViewer( element, config );

  $('#btn-02').on('click', function() {
    let pdbUri = '../tutorials/doc-data/1ycr.pdb';
    jQuery.ajax( pdbUri, { 
      success: function(data) {
        let v = viewer;
        v.addModel( data, "pdb" );                       /* load data */
        v.setStyle({}, {cartoon: {color: 'spectrum'}});  /* style all atoms */
        v.zoomTo();                                      /* set camera */
        v.render();                                      /* render scene */
        v.zoom(1.2, 1000);                               /* slight zoom */
      },
      error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
      },
    });
  })
})
</script>


***Note on using external data files (CORS)***

By default, Javascript will only be allowed to load data from the same domain as the web page from which it has been invoked (i.e. if your web page is being served from "my.domain.com" then javascript on that web page will only be able to load data from "my.domain.com"). This is a standard security restriction called "Cross-origin resource sharing" [CORS](https://en.wikipedia.org/wiki/Cross-origin_resource_sharing) - there are ways around this restriction, however for the sake of this tutorial we assume that the external PDB data file resides on your server.


## Dynamic styles

One of the advantages of using the API is that you have greater control of deciding how 3Dmol should interact to other events going on in the browser (user interaction, javascript events, etc). 

Another advantage is that this also gives you more flexibility with styling.

Here we create our own function that will decide on the color of a particular atom based on its own properties:

```
  let colorAsSnake = function(atom) {
    return atom.resi % 2 == 0 ? 'white' : 'green';
  };
```

Then we can apply that colouring function to a particular atom selection (in this case, all atoms in chain 'A')
    
```
  viewer.setStyle({chain: 'A'}, {cartoon: {colorfunc: colorAsSnake}});
```

<div class="align-center">
  <button class="btn btn-primary" id="btn-03">Try it</button>
  <div id="container-03" class="mol-container"></div>
</div>

<script>
jQuery(function() {
  let element = $('#container-03');
  let config = { defaultcolors: $3Dmol.rasmolElementColors, backgroundColor: 'white' };
  let viewer = $3Dmol.createViewer( element, config );
  let colorAsSnake = function(atom) {
    return atom.resi % 2 == 0 ? 'white' : 'green';
  };
  $('#btn-03').on('click', function() {
    let pdbUri = '../tutorials/doc-data/1ycr.pdb';
    jQuery.ajax( pdbUri, { 
      success: function(data) {
        let v = viewer;
        v.addModel( data, "pdb" );
        v.setStyle( {}, {cartoon: { colorfunc: colorAsSnake }} );
        v.zoomTo();
        v.render();
        v.zoom(1.2, 1000);
      },
      error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
      },
    });
  })
})
</script>

