/** contextmenu.js
 *
 * Test updates for context menu.
 * Original behavior supported invoking context menu by right-clicking on an
 * atom in the same absolute parent.
 *
 * Updated behavior:
 * - Support longtouch to open context menu
 * - Support context menu for Shapes
 * - Suppress "click" events (mouseup, touchup) when context menu is invoked
 * - Pass the mouse event to the userContextMenuHandler in case of a different
 *   offset parent.
 *
 * How this tests the features:
 * For View1, the canvas and the menu have the same absolute parent element, so
 * the click handler can use the default x and y to position the menu.
 * For View2, the canvas and the menu have different absolute parent elements,
 * so the click handler needs other mouseevent coords to position the menu.
 *
 * 1. View1 simulates mouse right-clicks:
 *    - Right-clicking an atom should open the menu at the atom
 *    - The "Click" label should not be seen if the click event is suppressed.
 * 2. View2 simulates right-click on background and longtouch on Shapes
 *    - Right-clicking on the background should open the menu at click location.
 *    - Longpress on the Sphere should open the context menu at the sphere.
 *    - The "Click" label should not be seen if click event is suppressed.
 */

function testView1(viewer) {
    viewer.userContextMenuHandler = (selected, x, y, intersects, ev) => {
        // Since the menu is in the same absolute parent as the viewer,
        // we can use the same x,y coordinates.
        openMenu(viewer, selected, x, y);
    }
    setupViewer(viewer);
    console.log('Testing right-click on atom in View1');
    rightClickAt(viewer, { x: 68, y: 120 });
}

function testView2(viewer) {
    viewer.userContextMenuHandler = (selected, x, y, intersects, ev) => {
        // Since the menu has a different absolute parent as the canvas,
        // we need to look at the mouseEvent to get preferred menu coordinates.
        openMenu(viewer, selected, ev?.pageX || x, ev?.pageY || y);
        if (!ev) {
            console.warn(
                'event not passed into handler; View2 will has offset position'
            );
        }
    }
    viewer.addSphere({
        center: { x: -1.8, y: 1.67, z: 0 }, radius: 1.0,
        color: 'purple', alpha: 0.4,
        contextMenuEnabled: true,
        clickable: true,
        callback: (shape) => viewer.addLabel("Click", { position: shape }),
    });

    setupViewer(viewer);
    console.log('Testing right-click on background in View2');
    rightClickAt(viewer, { x: 229, y: 369 });
    console.log('Testing longtouch on sphere in View2');
    // Add a slight change of position to long touch, because physically, it is
    // hard to have completely zero movement when doing a long touch.
    longTouchAt(viewer, { x: 270, y: 270, endX: 270, endY: 270 });
}

function rightClickAt(viewer, { x, y, endX=0, endY=0 }) {
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

function longTouchAt(viewer, { x, y, endX=0, endY=0 }) {
    // Emulate the sequence of events when a user clicks to open context menu:
    // 1. touchstart (starts the longtouch timer)
    // 2. touchmove (allow a small amount of movement for longtouch contextmenu)
    // 3. touchend (fires after longtouch timer completes)
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
        2000
    );
}

function setupViewer(viewer) {
    viewer.setHeight(180);
    viewer.setWidth(180);
    viewer.zoomTo();

    // Note: atoms have to be clickable to support context menu.
    // This is because they only get an insersectObject if they are clickable.
    viewer.setClickable({}, true, (atom) => {
        viewer.addLabel("Click", { position: atom });
    });
    viewer.enableContextMenu({},true);
    viewer.render( );
}

function openMenu(viewer, selected, x, y) {
    const type = targetType(selected);
    const menuElt = document.getElementById(`${type}-menu`);
    console.log(`${type} Menu ${x}, ${y}, ${selected?.x}, ${selected?.y}`);
    menuElt.style.left = x;
    menuElt.style.top = y;
    menuElt.style.display = 'block';
}

function targetType(selected) {
    if (!selected) return 'background';
    return selected.atom ? 'atom' : 'shape';
}

  /* @div
  <div style="width: 400px; height: 400px;">
    <!-- View1 has canvas and menu inside the same absolute parent element -->
    <div id="View1" style="position: absolute; left: 20px; top: 20px; border: black solid 1px;">
      <div  class='viewer_3Dmoljs'  style="width: 180px; height: 180px;" data-backgroundColor="white" data-href="../test_structs/benzene.sdf" data-style="stick" data-callback='testView1'></div>
      <div id="atom-menu" style="display: none; position: absolute; background-color: rgba(255, 255, 255, .50); color: blue; border: solid blue 2px; z-index: 100;">
        <div>View1 Mouse Atom Menu</span></div>
      </div>
    </div>

    <!--
    View2 has the menu outside the canvas' absolute parent.
    This means that the default x,y coming from GLViewer can't be used
    to position the menu, since they are relative to the absolute parent.
    -->
    <div id="View2" style="position: absolute; left: 220px; top: 220px; border: black solid 1px;">
      <div  class='viewer_3Dmoljs'  style="width: 180px; height: 180px;" data-backgroundColor="white" data-href="../test_structs/benzene.sdf" data-style="stick" data-callback='testView2'></div>
    </div>
    <div id="shape-menu" style="display: none; position: absolute; background-color: rgba(255, 255, 255, .50); color: blue; border: solid limegreen 2px; z-index: 100;">
      <div>View2 Shape Menu</div>
    </div>
    <div id="background-menu" style="display: none; position: absolute; background-color: rgba(255, 255, 255, .50); color: blue; border: solid magenta 2px; z-index: 100;">
      <div>View2 Background Menu</div>
    </div>
  </div>
*/


