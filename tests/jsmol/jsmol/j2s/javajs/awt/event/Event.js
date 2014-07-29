Clazz.declarePackage ("javajs.awt.event");
c$ = Clazz.decorateAsClass (function () {
this.source = null;
Clazz.instantialize (this, arguments);
}, javajs.awt.event, "Event");
Clazz.defineMethod (c$, "getSource", 
function () {
return this.source;
});
Clazz.defineStatics (c$,
"MOUSE_LEFT", 16,
"MOUSE_MIDDLE", 8,
"MOUSE_RIGHT", 4,
"MOUSE_WHEEL", 32,
"MAC_COMMAND", 20,
"BUTTON_MASK", 28,
"MOUSE_DOWN", 501,
"MOUSE_UP", 502,
"MOUSE_MOVE", 503,
"MOUSE_ENTER", 504,
"MOUSE_EXIT", 505,
"MOUSE_DRAG", 506,
"SHIFT_MASK", 1,
"ALT_MASK", 8,
"CTRL_MASK", 2,
"CTRL_ALT", 10,
"CTRL_SHIFT", 3,
"META_MASK", 4,
"VK_SHIFT", 16,
"VK_ALT", 18,
"VK_CONTROL", 17,
"VK_META", 157,
"VK_LEFT", 37,
"VK_RIGHT", 39,
"VK_PERIOD", 46,
"VK_SPACE", 32,
"VK_DOWN", 40,
"VK_UP", 38,
"VK_ESCAPE", 27,
"VK_DELETE", 127,
"VK_BACK_SPACE", 8,
"VK_PAGE_DOWN", 34,
"VK_PAGE_UP", 33,
"MOVED", 0,
"DRAGGED", 1,
"CLICKED", 2,
"WHEELED", 3,
"PRESSED", 4,
"RELEASED", 5);
