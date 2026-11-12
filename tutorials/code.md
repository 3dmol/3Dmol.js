<script src="/doc/_assets/tutorial.js" defer></script>
  
<style>
.mol-container {
  width: 100%;
  height: 400px;
  position: relative;
  border: 1px solid #999;
}

.tutorial-example {
  width: 100%;
  max-width: 600px;
  margin: 20px auto 10px;
  text-align: center;
}

.tutorial-try-button {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  padding: 0.6rem 1.2rem;
  margin-bottom: 0.75rem;
  border: 1px solid #0056b3;
  border-radius: 0.375rem;
  background: #0069d9;
  color: #fff;
  font: inherit;
  font-weight: 600;
  line-height: 1.25;
  cursor: pointer;
  box-shadow: 0 1px 2px rgb(0 0 0 / 15%);
}

.tutorial-try-button:hover {
  background: #0056b3;
}

.tutorial-try-button:focus-visible {
  outline: 3px solid var(--clean-accent, #007bff);
  outline-offset: 3px;
}

.tutorial-try-button:active {
  background: #004085;
}

.tutorial-status:empty {
  display: none;
}
</style>

The most simple way to display and style a 3Dmol instance is by embedding 3Dmol parameters [within a URL](./tutorial-url.html) or [within HTML tags](./tutorial-embeddable.html). 

However, if you are interacting with this viewer from your own web pages then you may want to take advantage of the 3Dmol API.

Here we will use the API to create, load and style the 3Dmol instance. This tutorial assumes a basic knowledge of HTML and JavaScript.

### Creating the 3DMol instance

First, make sure your HTML document links to the 3Dmol code.

```
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
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
  let element = document.querySelector('#container-01');
  let config = { backgroundColor: 'orange' };
  let viewer = $3Dmol.createViewer( element, config );
  viewer.addSphere({ center: {x:0, y:0, z:0}, radius: 10.0, color: 'green' });
  viewer.zoomTo();
  viewer.render();
  viewer.zoom(0.8, 2000);
</script>
```

Documentation: [$3Dmol.createViewer()](global.html#createViewer)

<div class="tutorial-example">
  <button type="button" id="btn-01" class="tutorial-try-button" aria-controls="container-01">Try it</button>
  <p class="tutorial-status" role="status"></p>
  <div id="container-01" class="mol-container"></div>
</div>

If this has worked, you should see a rather fetching green ball in front of an orange background. If not, then now's a good time to get familiar with the developer console on your favourite browser and check for typos ([Firefox](https://developer.mozilla.org/en/docs/Tools/Web_Console), [Chrome](https://developers.google.com/web/tools/chrome-devtools/), ...)

Documentation: [$3Dmol.GLViewer API](GLViewer.html)

### Loading data dynamically

Use `fetch()` to load PDB data from a URL, then pass the text to the viewer. Check the response before using it so that a missing file produces a useful error.

```
  async function loadProtein(viewer, pdbUri) {
    const response = await fetch(pdbUri);
    if (!response.ok) {
      throw new Error(`Failed to load PDB: ${response.status}`);
    }
    const data = await response.text();
    viewer.removeAllModels();
    viewer.addModel(data, 'pdb');
    viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
    viewer.zoomTo();
    viewer.render();
    viewer.zoom(1.2, 1000);
  }

  const viewer = $3Dmol.createViewer(element, config);
  loadProtein(viewer, '/path/to/your/pdb/files/1ycr.pdb')
    .catch(error => console.error(error));
```

All these methods (and more) are covered in detail in the [API documentation for GLViewer](GLViewer.html).

The interactive examples below include the 1YCR protein, so they also work when this documentation is opened as a local file.

<div class="tutorial-example">
  <button type="button" class="tutorial-try-button" id="btn-02" aria-controls="container-02">Try it</button>
  <p class="tutorial-status" role="status"></p>
  <div id="container-02" class="mol-container"></div>
</div>


***Note on using external data files (CORS)***

When fetching a file from another origin, that server must allow the request through [CORS](https://developer.mozilla.org/en-US/docs/Web/HTTP/Guides/CORS). Browsers also restrict fetching local files through `file://`. Serve your own page and PDB file over HTTP, or use a data server that permits cross-origin requests.


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

<div class="tutorial-example">
  <button type="button" class="tutorial-try-button" id="btn-03" aria-controls="container-03">Try it</button>
  <p class="tutorial-status" role="status"></p>
  <div id="container-03" class="mol-container"></div>
</div>

