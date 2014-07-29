Clazz.declarePackage ("J.c");
Clazz.load (["java.lang.Enum"], "J.c.CBK", ["JU.SB"], function () {
c$ = Clazz.declareType (J.c, "CBK", Enum);
c$.getCallback = Clazz.defineMethod (c$, "getCallback", 
function (name) {
name = name.toUpperCase ();
name = name.substring (0, Math.max (name.indexOf ("CALLBACK"), 0));
for (var item, $item = 0, $$item = J.c.CBK.values (); $item < $$item.length && ((item = $$item[$item]) || true); $item++) if (item.name ().equalsIgnoreCase (name)) return item;

return null;
}, "~S");
c$.getNameList = Clazz.defineMethod (c$, "getNameList", 
function () {
if (J.c.CBK.nameList == null) {
var names =  new JU.SB ();
for (var item, $item = 0, $$item = J.c.CBK.values (); $item < $$item.length && ((item = $$item[$item]) || true); $item++) names.append (item.name ().toLowerCase ()).append ("Callback;");

J.c.CBK.nameList = names.toString ();
}return J.c.CBK.nameList;
});
c$.nameList = null;
Clazz.defineEnumConstant (c$, "ANIMFRAME", 0, []);
Clazz.defineEnumConstant (c$, "APPLETREADY", 1, []);
Clazz.defineEnumConstant (c$, "ATOMMOVED", 2, []);
Clazz.defineEnumConstant (c$, "CLICK", 3, []);
Clazz.defineEnumConstant (c$, "DRAGDROP", 4, []);
Clazz.defineEnumConstant (c$, "ECHO", 5, []);
Clazz.defineEnumConstant (c$, "ERROR", 6, []);
Clazz.defineEnumConstant (c$, "EVAL", 7, []);
Clazz.defineEnumConstant (c$, "HOVER", 8, []);
Clazz.defineEnumConstant (c$, "LOADSTRUCT", 9, []);
Clazz.defineEnumConstant (c$, "MEASURE", 10, []);
Clazz.defineEnumConstant (c$, "MESSAGE", 11, []);
Clazz.defineEnumConstant (c$, "MINIMIZATION", 12, []);
Clazz.defineEnumConstant (c$, "PICK", 13, []);
Clazz.defineEnumConstant (c$, "RESIZE", 14, []);
Clazz.defineEnumConstant (c$, "SCRIPT", 15, []);
Clazz.defineEnumConstant (c$, "SYNC", 16, []);
Clazz.defineEnumConstant (c$, "STRUCTUREMODIFIED", 17, []);
});
