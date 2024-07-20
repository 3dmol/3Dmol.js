/** testcontextmenu.js
 *
 * Test updates for context menu.
 * Original behavior supported invoking context menu by right-clicking on an
 * atom in the same offsetParent.
 *
 * Updated behavior:
 * - Support longtouch to open context menu
 * - Support context menu for Shapes
 * - Suppress "click" events (mouseup, touchup) when context menu is invoked
 * - Pass the mouse event to the userContextMenuHandler in case of a different
 *   offset parent.
 *
 * This has two GLViewers, with different offsetParent configurations.
 * 1. View1 tests right-clicks on atoms when the menu has the same offsetParent
 *    as the canvas
 *    - Right-clicking an atom should open the menu with top,left at the click
 *    - The "Click" label should not be seen if the click event is suppressed.
 * 2. View2 tests right-clicks on the background and longtouch on Shapes
 *    when the menu has a different offsetParent from the canvas.
 *    - Background right-clicks should open the menu with top,left at the click
 *    - longtouch on the Sphere should open the menu with top,left at the touch
 *    - The "Click" label should not be seen if click event is suppressed.
 *
 * Note: this test has two 3Dmol viewers and positioned html decorations.
 * In order for glcheck to check it all, exportAllToCanvas gathers everything
 * into resultCanvas, which is the one found by glcheck.
 */

/* skip */ //does not fit into evaluation framework, but apparently works for glcheck?
  
/* @script
$3Dmol.tstctx_testView1 = function(viewer) {
    if (typeof(PointerEvent) === 'undefined') {
        // These tests can be executed by the browser, glcheck, or jest (which uses JSDOM).
        // In the jest/jsdom case, PointerEvents are not defined, so mock them here.
        $3Dmol.tstctx_isJest = true;
        window.PointerEvent = class PointerEvent extends CustomEvent {
            constructor(type, eventInitDict) {
                super(type, eventInitDict);
                for (const field of [
                    'pageX', 'pageY',
                    'screenX', 'screenY',
                    'clientX', 'clientY',
                    'targetTouches'
                ]) {
                    if (eventInitDict[field]) this[field] = eventInitDict[field];
                }
            }
        };
    }

    viewer.userContextMenuHandler = (selected, x, y, intersects, ev) => {
        // Since the menu has the same offsetParent as the viewer,
        // we can use the same x,y coordinates.
        $3Dmol.tstctx_openMenu(viewer, selected, x, y);
    }
    $3Dmol.tstctx_setupViewer(viewer);
    console.log('Testing right-click on atom in View1');
    // Target 8-oclock H in View1
    const { x, y } = $3Dmol.tstctx_getCoords(68, 120, 'View1');
    $3Dmol.tstctx_rightClickAt(viewer, { x, y });
}

$3Dmol.tstctx_testView2 = function(viewer) {
    viewer.userContextMenuHandler = (selected, x, y, intersects, ev) => {
        // Since the menu has a different offsetParent from the canvas,
        // we need to look at the mouseEvent to get preferred menu coordinates.
        $3Dmol.tstctx_openMenu(viewer, selected, ev?.pageX || x, ev?.pageY || y);
        if (!ev) {
            console.warn(
                'event not passed into handler; View2 has offset position'
            );
        }
    }
    viewer.addSphere({
        // Sphere added around 10-oclock H in View2
        center: { x: -1.8, y: 1.67, z: 0 }, radius: 1.0,
        color: 'purple', alpha: 0.4,
        contextMenuEnabled: true,
        clickable: true,
        callback: (shape) => viewer.addLabel("Click", { position: shape }),
    });

    $3Dmol.tstctx_setupViewer(viewer);
    console.log('Testing right-click on background in View2');
    let { x, y } = $3Dmol.tstctx_getCoords(229, 369, 'View2');;
    $3Dmol.tstctx_rightClickAt(viewer, { x, y });
    console.log('Testing longtouch on sphere in View2');
    ({ x, y } = $3Dmol.tstctx_getCoords(270, 270, 'View2'));
    // Add a slight change of position to long touch, because physically, it is
    // hard to have completely zero movement when doing a long touch.
    const { x: endX, y: endY } = $3Dmol.tstctx_getCoords(275, 265, 'View2');
    $3Dmol.tstctx_longTouchAt(viewer, { x, y, endX, endY });

    // Test that touch-dragging too far away (> 5) cancels the contextmenu event
    // Note: this introduces a slight rotation to the benzene. The reference image
    // doesn't have this, but it is subtle enough to pass the pixelmatch threshold.
    setTimeout(() => {
        const [newX, newY] = [x+10, y+10];
        $3Dmol.tstctx_longTouchAt(viewer, { x: newX, y: newY, endX: newX - 0.5, endY: newY + 5.1 });
    }, 1500);

    // Tell glcheck to complete the test after other touch timers have completed
    setTimeout($3Dmol.tstctx_finish, 3000);
}

$3Dmol.tstctx_rightClickAt = function(viewer, { x, y, endX=0, endY=0 }) {
    // Emulate the sequence of events when a user clicks to open context menu:
    // 1. mousedown
    // 2. mousemove (if there is movement)
    // 2. mouseup (GLViewer will maybe interpret this as a "click")
    // 3. contextmenu
    // endX,endY params allow to specify movement between mousedown and mouseup;
    // for mouse right-clicks, contextmenu is only fired if there's no movement.
    const preventDefault = () => {};
    const startEvt = { pageX: x, pageY: y, preventDefault, which: 3 };
    const endEvt = { pageX: endX || x, pageY: endY || y, preventDefault};
    viewer._handleMouseDown(startEvt);
    if (endX || endY) {
        viewer._handleMouseMove(endEvt);
    }
    viewer._handleMouseUp(endEvt);
    viewer._handleContextMenu(endEvt);
}

$3Dmol.tstctx_longTouchAt = function(viewer, { x, y, endX=0, endY=0 }) {
    // Emulate the sequence of events when user does longtouch to open menu:
    // 1. touchstart (starts the longtouch timer)
    // 2. touchmove (allow a small amount of movement for longtouch contextmenu)
    // 3. touchend (fired after longtouch timer completes)
    const preventDefault = () => {};
    const startEvt = {
        type: 'touchstart',
        targetTouches: [{ pageX: x, pageY: y, clientX: x, clientY: y }],
        preventDefault,
    };
    const endEvt = {
        type: 'touchmove',
        targetTouches: [{
            pageX: endX || x, pageY: endY || y,
            clientX: endX || x, clientY: endY || y,
        }],
        preventDefault,
    };
    viewer._handleMouseDown(startEvt);
    viewer._handleMouseMove(endEvt);
    setTimeout(
        () => viewer._handleMouseUp({ ...endEvt, type: 'touchend' }),
        1200 // after longtouch timeout (1000)
    );
}

// HELPER FUNCTIONS
$3Dmol.tstctx_getCoords = function(x, y, viewerId) {
    if (!$3Dmol.tstctx_isJest) return { x, y };
    // In jest/jsdom tests, the offsets are not the same as glcheck / browser
    // Subtract the left,top numbers for View1/2 divs in the html, plus border
    const offsetX = viewerId === 'View1' ? (20 + 1) : (220 + 1);
    const offsetY = viewerId === 'View1' ? (20 + 1) : (220 + 1);
    return {
        x: x - offsetX,
        y: y - offsetY,
    };
}

$3Dmol.tstctx_setupViewer = function(viewer) {
    viewer.setHeight(180);
    viewer.setWidth(180);
    viewer.zoomTo();

    // Note: atoms have to be clickable to support context menu.
    // This is because they only get an insersectObject if they are clickable.
    viewer.setClickable({}, true, (atom) => {
        viewer.addLabel("Click", { position: atom });
        console.log(`Clicked on atom ${atom.atom} at ${atom.x},${atom.y},${atom.z}`);
    });
    viewer.enableContextMenu({},true);
    viewer.render( ); // This space is necessary to prevent "callback" insertion
}

$3Dmol.tstctx_finish = async function() {
    await $3Dmol.tstctx_exportAllToCanvas();
    callback(); // provided and needed by test machinery
}

// glcheck checks the first canvas it sees, so need to gather the two views and
// the HTML elements for the "context menus" into a single canvas for glcheck.
$3Dmol.tstctx_exportAllToCanvas = async function() {
    const canvas = resultCanvas;
    const ctx = canvas.getContext('2d');
    const menusImg = await nodeToImg(workingDiv, canvas.width, canvas.height);
    // All these numbers are doubles of the styles from the html below
    ctx.drawImage(View1Canvas, 40, 40, 360, 360);
    ctx.drawImage(View2Canvas, 440, 440, 360, 360);
    ctx.drawImage(menusImg, 0, 0, 1600, 1600);

    async function nodeToImg(node, width, height) {
        // Helper function to convert the html decorations to an image that can be
        // loaded into resultsCanvas. Not perfect but good enough for pixelmatch!
        // Represent html as svg, load svg into <img>
        const html = node.outerHTML;
        const svg = `
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}">
<foreignObject width="100%" height="100%">
    <div xmlns="http://www.w3.org/1999/xhtml">${html}</div>
</foreignObject>
</svg>`;
        return await svgToImg(svg);
    }

    function svgToImg(svg) {
        return new Promise((resolve, reject) => {
            const loadTimer = setTimeout(() => {
                // jest tests use jsdom which doesn't load images by default.
                // In that case, the onload will never fire to resolve the Promise.
                // This setTimeout ensures that the promise will complete.
                // Since jest only executes render tests for code coverage, not for
                // correctness, and this image is only used by glcheck,
                // we can resolve the promise as success without actually loading the image.
                resolve(img);
            }, 5000);

            const img = new Image();
            img.onload = () => {
                clearTimeout(loadTimer);
                resolve(img);
            };
            img.onerror = reject;
            img.src = 'data:image/svg+xml,' + encodeURIComponent(svg);
        });
    }
}

$3Dmol.tstctx_openMenu = function(viewer, selected, x, y) {
    const type = selected ? (selected.atom ? 'atom' : 'shape') : 'background';
    const menuElt = document.getElementById(`${type}-menu`);
    console.log(`Context menu invoked for ${type} at ${x}, ${y}, target at ${selected?.x}, ${selected?.y}`);
    menuElt.style.left = x;
    menuElt.style.top = y;
    menuElt.style.display = 'block';
}

*/

/* @div
<div>
  <!-- This div holds the resultsCanvas, to combine the two 3Dmol canvases and the HTML decorations -->
  <div id="resultDiv" style="position: absolute; left: 402px; width: 800px; height: 800px;">
    <canvas id="resultCanvas" height=800 width=800>
  </div>

  <!-- This div hold the 3Dmol test canvases -->
  <div id="workingDiv" style="width: 400px; height: 400px;">
    <!-- View1 has canvas and menu inside the same offsetParent element -->
    <div id="View1" style="position: absolute; left: 20px; top: 20px; border: black solid 1px;">
      <div class='viewer_3Dmoljs'  style="width: 180px; height: 180px;" data-backgroundColor="white" data-href="../test_structs/benzene.sdf" data-style="stick" data-callback='$3Dmol.tstctx_testView1' data-config="id:View1Canvas"></div>
      <div id="atom-menu" style="display: none; position: absolute; background-color: rgba(255, 255, 255, .50); color: blue; border: solid blue 2px; z-index: 100;">
        <div>View1 Mouse Atom Menu</span></div>
      </div>
    </div>

    <!--
    View2 has the menu outside the canvas' offsetParent.
    This means that the default x,y coming from GLViewer can't be used
    to position the menu, since they are relative to the canvas's offsetParent.
    -->
    <div id="View2" style="position: absolute; left: 220px; top: 220px; border: black solid 1px;">
      <div class='viewer_3Dmoljs'  style="width: 180px; height: 180px;" data-backgroundColor="white" data-href="../test_structs/benzene.sdf" data-style="stick" data-callback='$3Dmol.tstctx_testView2' data-config="id:View2Canvas"></div>
    </div>
    <div id="shape-menu" style="display: none; position: absolute; background-color: rgba(255, 255, 255, .50); color: blue; border: solid limegreen 2px; z-index: 100;">
      <div>View2 Shape Menu</div>
    </div>
    <div id="background-menu" style="display: none; position: absolute; background-color: rgba(255, 255, 255, .50); color: blue; border: solid magenta 2px; z-index: 100;">
      <div>View2 Background Menu</div>
    </div>
  </div>
</div>
*/
