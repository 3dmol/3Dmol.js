/*! jQuery v1.7 jquery.com | jquery.org/license */
(function(a,b){function cA(a){return f.isWindow(a)?a:a.nodeType===9?a.defaultView||a.parentWindow:!1}function cx(a){if(!cm[a]){var b=c.body,d=f("<"+a+">").appendTo(b),e=d.css("display");d.remove();if(e==="none"||e===""){cn||(cn=c.createElement("iframe"),cn.frameBorder=cn.width=cn.height=0),b.appendChild(cn);if(!co||!cn.createElement)co=(cn.contentWindow||cn.contentDocument).document,co.write((c.compatMode==="CSS1Compat"?"<!doctype html>":"")+"<html><body>"),co.close();d=co.createElement(a),co.body.appendChild(d),e=f.css(d,"display"),b.removeChild(cn)}cm[a]=e}return cm[a]}function cw(a,b){var c={};f.each(cs.concat.apply([],cs.slice(0,b)),function(){c[this]=a});return c}function cv(){ct=b}function cu(){setTimeout(cv,0);return ct=f.now()}function cl(){try{return new a.ActiveXObject("Microsoft.XMLHTTP")}catch(b){}}function ck(){try{return new a.XMLHttpRequest}catch(b){}}function ce(a,c){a.dataFilter&&(c=a.dataFilter(c,a.dataType));var d=a.dataTypes,e={},g,h,i=d.length,j,k=d[0],l,m,n,o,p;for(g=1;g<i;g++){if(g===1)for(h in a.converters)typeof h=="string"&&(e[h.toLowerCase()]=a.converters[h]);l=k,k=d[g];if(k==="*")k=l;else if(l!=="*"&&l!==k){m=l+" "+k,n=e[m]||e["* "+k];if(!n){p=b;for(o in e){j=o.split(" ");if(j[0]===l||j[0]==="*"){p=e[j[1]+" "+k];if(p){o=e[o],o===!0?n=p:p===!0&&(n=o);break}}}}!n&&!p&&f.error("No conversion from "+m.replace(" "," to ")),n!==!0&&(c=n?n(c):p(o(c)))}}return c}function cd(a,c,d){var e=a.contents,f=a.dataTypes,g=a.responseFields,h,i,j,k;for(i in g)i in d&&(c[g[i]]=d[i]);while(f[0]==="*")f.shift(),h===b&&(h=a.mimeType||c.getResponseHeader("content-type"));if(h)for(i in e)if(e[i]&&e[i].test(h)){f.unshift(i);break}if(f[0]in d)j=f[0];else{for(i in d){if(!f[0]||a.converters[i+" "+f[0]]){j=i;break}k||(k=i)}j=j||k}if(j){j!==f[0]&&f.unshift(j);return d[j]}}function cc(a,b,c,d){if(f.isArray(b))f.each(b,function(b,e){c||bG.test(a)?d(a,e):cc(a+"["+(typeof e=="object"||f.isArray(e)?b:"")+"]",e,c,d)});else if(!c&&b!=null&&typeof b=="object")for(var e in b)cc(a+"["+e+"]",b[e],c,d);else d(a,b)}function cb(a,c){var d,e,g=f.ajaxSettings.flatOptions||{};for(d in c)c[d]!==b&&((g[d]?a:e||(e={}))[d]=c[d]);e&&f.extend(!0,a,e)}function ca(a,c,d,e,f,g){f=f||c.dataTypes[0],g=g||{},g[f]=!0;var h=a[f],i=0,j=h?h.length:0,k=a===bV,l;for(;i<j&&(k||!l);i++)l=h[i](c,d,e),typeof l=="string"&&(!k||g[l]?l=b:(c.dataTypes.unshift(l),l=ca(a,c,d,e,l,g)));(k||!l)&&!g["*"]&&(l=ca(a,c,d,e,"*",g));return l}function b_(a){return function(b,c){typeof b!="string"&&(c=b,b="*");if(f.isFunction(c)){var d=b.toLowerCase().split(bR),e=0,g=d.length,h,i,j;for(;e<g;e++)h=d[e],j=/^\+/.test(h),j&&(h=h.substr(1)||"*"),i=a[h]=a[h]||[],i[j?"unshift":"push"](c)}}}function bE(a,b,c){var d=b==="width"?a.offsetWidth:a.offsetHeight,e=b==="width"?bz:bA;if(d>0){c!=="border"&&f.each(e,function(){c||(d-=parseFloat(f.css(a,"padding"+this))||0),c==="margin"?d+=parseFloat(f.css(a,c+this))||0:d-=parseFloat(f.css(a,"border"+this+"Width"))||0});return d+"px"}d=bB(a,b,b);if(d<0||d==null)d=a.style[b]||0;d=parseFloat(d)||0,c&&f.each(e,function(){d+=parseFloat(f.css(a,"padding"+this))||0,c!=="padding"&&(d+=parseFloat(f.css(a,"border"+this+"Width"))||0),c==="margin"&&(d+=parseFloat(f.css(a,c+this))||0)});return d+"px"}function br(a,b){b.src?f.ajax({url:b.src,async:!1,dataType:"script"}):f.globalEval((b.text||b.textContent||b.innerHTML||"").replace(bi,"/*$0*/")),b.parentNode&&b.parentNode.removeChild(b)}function bq(a){var b=(a.nodeName||"").toLowerCase();b==="input"?bp(a):b!=="script"&&typeof a.getElementsByTagName!="undefined"&&f.grep(a.getElementsByTagName("input"),bp)}function bp(a){if(a.type==="checkbox"||a.type==="radio")a.defaultChecked=a.checked}function bo(a){return typeof a.getElementsByTagName!="undefined"?a.getElementsByTagName("*"):typeof a.querySelectorAll!="undefined"?a.querySelectorAll("*"):[]}function bn(a,b){var c;if(b.nodeType===1){b.clearAttributes&&b.clearAttributes(),b.mergeAttributes&&b.mergeAttributes(a),c=b.nodeName.toLowerCase();if(c==="object")b.outerHTML=a.outerHTML;else if(c!=="input"||a.type!=="checkbox"&&a.type!=="radio"){if(c==="option")b.selected=a.defaultSelected;else if(c==="input"||c==="textarea")b.defaultValue=a.defaultValue}else a.checked&&(b.defaultChecked=b.checked=a.checked),b.value!==a.value&&(b.value=a.value);b.removeAttribute(f.expando)}}function bm(a,b){if(b.nodeType===1&&!!f.hasData(a)){var c,d,e,g=f._data(a),h=f._data(b,g),i=g.events;if(i){delete h.handle,h.events={};for(c in i)for(d=0,e=i[c].length;d<e;d++)f.event.add(b,c+(i[c][d].namespace?".":"")+i[c][d].namespace,i[c][d],i[c][d].data)}h.data&&(h.data=f.extend({},h.data))}}function bl(a,b){return f.nodeName(a,"table")?a.getElementsByTagName("tbody")[0]||a.appendChild(a.ownerDocument.createElement("tbody")):a}function X(a){var b=Y.split(" "),c=a.createDocumentFragment();if(c.createElement)while(b.length)c.createElement(b.pop());return c}function W(a,b,c){b=b||0;if(f.isFunction(b))return f.grep(a,function(a,d){var e=!!b.call(a,d,a);return e===c});if(b.nodeType)return f.grep(a,function(a,d){return a===b===c});if(typeof b=="string"){var d=f.grep(a,function(a){return a.nodeType===1});if(R.test(b))return f.filter(b,d,!c);b=f.filter(b,d)}return f.grep(a,function(a,d){return f.inArray(a,b)>=0===c})}function V(a){return!a||!a.parentNode||a.parentNode.nodeType===11}function N(){return!0}function M(){return!1}function n(a,b,c){var d=b+"defer",e=b+"queue",g=b+"mark",h=f._data(a,d);h&&(c==="queue"||!f._data(a,e))&&(c==="mark"||!f._data(a,g))&&setTimeout(function(){!f._data(a,e)&&!f._data(a,g)&&(f.removeData(a,d,!0),h.fire())},0)}function m(a){for(var b in a){if(b==="data"&&f.isEmptyObject(a[b]))continue;if(b!=="toJSON")return!1}return!0}function l(a,c,d){if(d===b&&a.nodeType===1){var e="data-"+c.replace(k,"-$1").toLowerCase();d=a.getAttribute(e);if(typeof d=="string"){try{d=d==="true"?!0:d==="false"?!1:d==="null"?null:f.isNumeric(d)?parseFloat(d):j.test(d)?f.parseJSON(d):d}catch(g){}f.data(a,c,d)}else d=b}return d}function h(a){var b=g[a]={},c,d;a=a.split(/\s+/);for(c=0,d=a.length;c<d;c++)b[a[c]]=!0;return b}var c=a.document,d=a.navigator,e=a.location,f=function(){function K(){if(!e.isReady){try{c.documentElement.doScroll("left")}catch(a){setTimeout(K,1);return}e.ready()}}var e=function(a,b){return new e.fn.init(a,b,h)},f=a.jQuery,g=a.$,h,i=/^(?:[^#<]*(<[\w\W]+>)[^>]*$|#([\w\-]*)$)/,j=/\S/,k=/^\s+/,l=/\s+$/,m=/\d/,n=/^<(\w+)\s*\/?>(?:<\/\1>)?$/,o=/^[\],:{}\s]*$/,p=/\\(?:["\\\/bfnrt]|u[0-9a-fA-F]{4})/g,q=/"[^"\\\n\r]*"|true|false|null|-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?/g,r=/(?:^|:|,)(?:\s*\[)+/g,s=/(webkit)[ \/]([\w.]+)/,t=/(opera)(?:.*version)?[ \/]([\w.]+)/,u=/(msie) ([\w.]+)/,v=/(mozilla)(?:.*? rv:([\w.]+))?/,w=/-([a-z]|[0-9])/ig,x=/^-ms-/,y=function(a,b){return(b+"").toUpperCase()},z=d.userAgent,A,B,C,D=Object.prototype.toString,E=Object.prototype.hasOwnProperty,F=Array.prototype.push,G=Array.prototype.slice,H=String.prototype.trim,I=Array.prototype.indexOf,J={};e.fn=e.prototype={constructor:e,init:function(a,d,f){var g,h,j,k;if(!a)return this;if(a.nodeType){this.context=this[0]=a,this.length=1;return this}if(a==="body"&&!d&&c.body){this.context=c,this[0]=c.body,this.selector=a,this.length=1;return this}if(typeof a=="string"){a.charAt(0)!=="<"||a.charAt(a.length-1)!==">"||a.length<3?g=i.exec(a):g=[null,a,null];if(g&&(g[1]||!d)){if(g[1]){d=d instanceof e?d[0]:d,k=d?d.ownerDocument||d:c,j=n.exec(a),j?e.isPlainObject(d)?(a=[c.createElement(j[1])],e.fn.attr.call(a,d,!0)):a=[k.createElement(j[1])]:(j=e.buildFragment([g[1]],[k]),a=(j.cacheable?e.clone(j.fragment):j.fragment).childNodes);return e.merge(this,a)}h=c.getElementById(g[2]);if(h&&h.parentNode){if(h.id!==g[2])return f.find(a);this.length=1,this[0]=h}this.context=c,this.selector=a;return this}return!d||d.jquery?(d||f).find(a):this.constructor(d).find(a)}if(e.isFunction(a))return f.ready(a);a.selector!==b&&(this.selector=a.selector,this.context=a.context);return e.makeArray(a,this)},selector:"",jquery:"1.7",length:0,size:function(){return this.length},toArray:function(){return G.call(this,0)},get:function(a){return a==null?this.toArray():a<0?this[this.length+a]:this[a]},pushStack:function(a,b,c){var d=this.constructor();e.isArray(a)?F.apply(d,a):e.merge(d,a),d.prevObject=this,d.context=this.context,b==="find"?d.selector=this.selector+(this.selector?" ":"")+c:b&&(d.selector=this.selector+"."+b+"("+c+")");return d},each:function(a,b){return e.each(this,a,b)},ready:function(a){e.bindReady(),B.add(a);return this},eq:function(a){return a===-1?this.slice(a):this.slice(a,+a+1)},first:function(){return this.eq(0)},last:function(){return this.eq(-1)},slice:function(){return this.pushStack(G.apply(this,arguments),"slice",G.call(arguments).join(","))},map:function(a){return this.pushStack(e.map(this,function(b,c){return a.call(b,c,b)}))},end:function(){return this.prevObject||this.constructor(null)},push:F,sort:[].sort,splice:[].splice},e.fn.init.prototype=e.fn,e.extend=e.fn.extend=function(){var a,c,d,f,g,h,i=arguments[0]||{},j=1,k=arguments.length,l=!1;typeof i=="boolean"&&(l=i,i=arguments[1]||{},j=2),typeof i!="object"&&!e.isFunction(i)&&(i={}),k===j&&(i=this,--j);for(;j<k;j++)if((a=arguments[j])!=null)for(c in a){d=i[c],f=a[c];if(i===f)continue;l&&f&&(e.isPlainObject(f)||(g=e.isArray(f)))?(g?(g=!1,h=d&&e.isArray(d)?d:[]):h=d&&e.isPlainObject(d)?d:{},i[c]=e.extend(l,h,f)):f!==b&&(i[c]=f)}return i},e.extend({noConflict:function(b){a.$===e&&(a.$=g),b&&a.jQuery===e&&(a.jQuery=f);return e},isReady:!1,readyWait:1,holdReady:function(a){a?e.readyWait++:e.ready(!0)},ready:function(a){if(a===!0&&!--e.readyWait||a!==!0&&!e.isReady){if(!c.body)return setTimeout(e.ready,1);e.isReady=!0;if(a!==!0&&--e.readyWait>0)return;B.fireWith(c,[e]),e.fn.trigger&&e(c).trigger("ready").unbind("ready")}},bindReady:function(){if(!B){B=e.Callbacks("once memory");if(c.readyState==="complete")return setTimeout(e.ready,1);if(c.addEventListener)c.addEventListener("DOMContentLoaded",C,!1),a.addEventListener("load",e.ready,!1);else if(c.attachEvent){c.attachEvent("onreadystatechange",C),a.attachEvent("onload",e.ready);var b=!1;try{b=a.frameElement==null}catch(d){}c.documentElement.doScroll&&b&&K()}}},isFunction:function(a){return e.type(a)==="function"},isArray:Array.isArray||function(a){return e.type(a)==="array"},isWindow:function(a){return a&&typeof a=="object"&&"setInterval"in a},isNumeric:function(a){return a!=null&&m.test(a)&&!isNaN(a)},type:function(a){return a==null?String(a):J[D.call(a)]||"object"},isPlainObject:function(a){if(!a||e.type(a)!=="object"||a.nodeType||e.isWindow(a))return!1;try{if(a.constructor&&!E.call(a,"constructor")&&!E.call(a.constructor.prototype,"isPrototypeOf"))return!1}catch(c){return!1}var d;for(d in a);return d===b||E.call(a,d)},isEmptyObject:function(a){for(var b in a)return!1;return!0},error:function(a){throw a},parseJSON:function(b){if(typeof b!="string"||!b)return null;b=e.trim(b);if(a.JSON&&a.JSON.parse)return a.JSON.parse(b);if(o.test(b.replace(p,"@").replace(q,"]").replace(r,"")))return(new Function("return "+b))();e.error("Invalid JSON: "+b)},parseXML:function(c){var d,f;try{a.DOMParser?(f=new DOMParser,d=f.parseFromString(c,"text/xml")):(d=new ActiveXObject("Microsoft.XMLDOM"),d.async="false",d.loadXML(c))}catch(g){d=b}(!d||!d.documentElement||d.getElementsByTagName("parsererror").length)&&e.error("Invalid XML: "+c);return d},noop:function(){},globalEval:function(b){b&&j.test(b)&&(a.execScript||function(b){a.eval.call(a,b)})(b)},camelCase:function(a){return a.replace(x,"ms-").replace(w,y)},nodeName:function(a,b){return a.nodeName&&a.nodeName.toUpperCase()===b.toUpperCase()},each:function(a,c,d){var f,g=0,h=a.length,i=h===b||e.isFunction(a);if(d){if(i){for(f in a)if(c.apply(a[f],d)===!1)break}else for(;g<h;)if(c.apply(a[g++],d)===!1)break}else if(i){for(f in a)if(c.call(a[f],f,a[f])===!1)break}else for(;g<h;)if(c.call(a[g],g,a[g++])===!1)break;return a},trim:H?function(a){return a==null?"":H.call(a)}:function(a){return a==null?"":(a+"").replace(k,"").replace(l,"")},makeArray:function(a,b){var c=b||[];if(a!=null){var d=e.type(a);a.length==null||d==="string"||d==="function"||d==="regexp"||e.isWindow(a)?F.call(c,a):e.merge(c,a)}return c},inArray:function(a,b,c){var d;if(b){if(I)return I.call(b,a,c);d=b.length,c=c?c<0?Math.max(0,d+c):c:0;for(;c<d;c++)if(c in b&&b[c]===a)return c}return-1},merge:function(a,c){var d=a.length,e=0;if(typeof c.length=="number")for(var f=c.length;e<f;e++)a[d++]=c[e];else while(c[e]!==b)a[d++]=c[e++];a.length=d;return a},grep:function(a,b,c){var d=[],e;c=!!c;for(var f=0,g=a.length;f<g;f++)e=!!b(a[f],f),c!==e&&d.push(a[f]);return d},map:function(a,c,d){var f,g,h=[],i=0,j=a.length,k=a instanceof e||j!==b&&typeof j=="number"&&(j>0&&a[0]&&a[j-1]||j===0||e.isArray(a));if(k)for(;i<j;i++)f=c(a[i],i,d),f!=null&&(h[h.length]=f);else for(g in a)f=c(a[g],g,d),f!=null&&(h[h.length]=f);return h.concat.apply([],h)},guid:1,proxy:function(a,c){if(typeof c=="string"){var d=a[c];c=a,a=d}if(!e.isFunction(a))return b;var f=G.call(arguments,2),g=function(){return a.apply(c,f.concat(G.call(arguments)))};g.guid=a.guid=a.guid||g.guid||e.guid++;return g},access:function(a,c,d,f,g,h){var i=a.length;if(typeof c=="object"){for(var j in c)e.access(a,j,c[j],f,g,d);return a}if(d!==b){f=!h&&f&&e.isFunction(d);for(var k=0;k<i;k++)g(a[k],c,f?d.call(a[k],k,g(a[k],c)):d,h);return a}return i?g(a[0],c):b},now:function(){return(new Date).getTime()},uaMatch:function(a){a=a.toLowerCase();var b=s.exec(a)||t.exec(a)||u.exec(a)||a.indexOf("compatible")<0&&v.exec(a)||[];return{browser:b[1]||"",version:b[2]||"0"}},sub:function(){function a(b,c){return new a.fn.init(b,c)}e.extend(!0,a,this),a.superclass=this,a.fn=a.prototype=this(),a.fn.constructor=a,a.sub=this.sub,a.fn.init=function(d,f){f&&f instanceof e&&!(f instanceof a)&&(f=a(f));return e.fn.init.call(this,d,f,b)},a.fn.init.prototype=a.fn;var b=a(c);return a},browser:{}}),e.each("Boolean Number String Function Array Date RegExp Object".split(" "),function(a,b){J["[object "+b+"]"]=b.toLowerCase()}),A=e.uaMatch(z),A.browser&&(e.browser[A.browser]=!0,e.browser.version=A.version),e.browser.webkit&&(e.browser.safari=!0),j.test(" ")&&(k=/^[\s\xA0]+/,l=/[\s\xA0]+$/),h=e(c),c.addEventListener?C=function(){c.removeEventListener("DOMContentLoaded",C,!1),e.ready()}:c.attachEvent&&(C=function(){c.readyState==="complete"&&(c.detachEvent("onreadystatechange",C),e.ready())}),typeof define=="function"&&define.amd&&define.amd.jQuery&&define("jquery",[],function(){return e});return e}(),g={};f.Callbacks=function(a){a=a?g[a]||h(a):{};var c=[],d=[],e,i,j,k,l,m=function(b){var d,e,g,h,i;for(d=0,e=b.length;d<e;d++)g=b[d],h=f.type(g),h==="array"?m(g):h==="function"&&(!a.unique||!o.has(g))&&c.push(g)},n=function(b,f){f=f||[],e=!a.memory||[b,f],i=!0,l=j||0,j=0,k=c.length;for(;c&&l<k;l++)if(c[l].apply(b,f)===!1&&a.stopOnFalse){e=!0;break}i=!1,c&&(a.once?e===!0?o.disable():c=[]:d&&d.length&&(e=d.shift(),o.fireWith(e[0],e[1])))},o={add:function(){if(c){var a=c.length;m(arguments),i?k=c.length:e&&e!==!0&&(j=a,n(e[0],e[1]))}return this},remove:function(){if(c){var b=arguments,d=0,e=b.length;for(;d<e;d++)for(var f=0;f<c.length;f++)if(b[d]===c[f]){i&&f<=k&&(k--,f<=l&&l--),c.splice(f--,1);if(a.unique)break}}return this},has:function(a){if(c){var b=0,d=c.length;for(;b<d;b++)if(a===c[b])return!0}return!1},empty:function(){c=[];return this},disable:function(){c=d=e=b;return this},disabled:function(){return!c},lock:function(){d=b,(!e||e===!0)&&o.disable();return this},locked:function(){return!d},fireWith:function(b,c){d&&(i?a.once||d.push([b,c]):(!a.once||!e)&&n(b,c));return this},fire:function(){o.fireWith(this,arguments);return this},fired:function(){return!!e}};return o};var i=[].slice;f.extend({Deferred:function(a){var b=f.Callbacks("once memory"),c=f.Callbacks("once memory"),d=f.Callbacks("memory"),e="pending",g={resolve:b,reject:c,notify:d},h={done:b.add,fail:c.add,progress:d.add,state:function(){return e},isResolved:b.fired,isRejected:c.fired,then:function(a,b,c){i.done(a).fail(b).progress(c);return this},always:function(){return i.done.apply(i,arguments).fail.apply(i,arguments)},pipe:function(a,b,c){return f.Deferred(function(d){f.each({done:[a,"resolve"],fail:[b,"reject"],progress:[c,"notify"]},function(a,b){var c=b[0],e=b[1],g;f.isFunction(c)?i[a](function(){g=c.apply(this,arguments),g&&f.isFunction(g.promise)?g.promise().then(d.resolve,d.reject,d.notify):d[e+"With"](this===i?d:this,[g])}):i[a](d[e])})}).promise()},promise:function(a){if(a==null)a=h;else for(var b in h)a[b]=h[b];return a}},i=h.promise({}),j;for(j in g)i[j]=g[j].fire,i[j+"With"]=g[j].fireWith;i.done(function(){e="resolved"},c.disable,d.lock).fail(function(){e="rejected"},b.disable,d.lock),a&&a.call(i,i);return i},when:function(a){function m(a){return function(b){e[a]=arguments.length>1?i.call(arguments,0):b,j.notifyWith(k,e)}}function l(a){return function(c){b[a]=arguments.length>1?i.call(arguments,0):c,--g||j.resolveWith(j,b)}}var b=i.call(arguments,0),c=0,d=b.length,e=Array(d),g=d,h=d,j=d<=1&&a&&f.isFunction(a.promise)?a:f.Deferred(),k=j.promise();if(d>1){for(;c<d;c++)b[c]&&b[c].promise&&f.isFunction(b[c].promise)?b[c].promise().then(l(c),j.reject,m(c)):--g;g||j.resolveWith(j,b)}else j!==a&&j.resolveWith(j,d?[a]:[]);return k}}),f.support=function(){var a=c.createElement("div"),b=c.documentElement,d,e,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u;a.setAttribute("className","t"),a.innerHTML="   <link/><table></table><a href='/a' style='top:1px;float:left;opacity:.55;'>a</a><input type='checkbox'/><nav></nav>",d=a.getElementsByTagName("*"),e=a.getElementsByTagName("a")[0];if(!d||!d.length||!e)return{};g=c.createElement("select"),h=g.appendChild(c.createElement("option")),i=a.getElementsByTagName("input")[0],k={leadingWhitespace:a.firstChild.nodeType===3,tbody:!a.getElementsByTagName("tbody").length,htmlSerialize:!!a.getElementsByTagName("link").length,style:/top/.test(e.getAttribute("style")),hrefNormalized:e.getAttribute("href")==="/a",opacity:/^0.55/.test(e.style.opacity),cssFloat:!!e.style.cssFloat,unknownElems:!!a.getElementsByTagName("nav").length,checkOn:i.value==="on",optSelected:h.selected,getSetAttribute:a.className!=="t",enctype:!!c.createElement("form").enctype,submitBubbles:!0,changeBubbles:!0,focusinBubbles:!1,deleteExpando:!0,noCloneEvent:!0,inlineBlockNeedsLayout:!1,shrinkWrapBlocks:!1,reliableMarginRight:!0},i.checked=!0,k.noCloneChecked=i.cloneNode(!0).checked,g.disabled=!0,k.optDisabled=!h.disabled;try{delete a.test}catch(v){k.deleteExpando=!1}!a.addEventListener&&a.attachEvent&&a.fireEvent&&(a.attachEvent("onclick",function(){k.noCloneEvent=!1}),a.cloneNode(!0).fireEvent("onclick")),i=c.createElement("input"),i.value="t",i.setAttribute("type","radio"),k.radioValue=i.value==="t",i.setAttribute("checked","checked"),a.appendChild(i),l=c.createDocumentFragment(),l.appendChild(a.lastChild),k.checkClone=l.cloneNode(!0).cloneNode(!0).lastChild.checked,a.innerHTML="",a.style.width=a.style.paddingLeft="1px",m=c.getElementsByTagName("body")[0],o=c.createElement(m?"div":"body"),p={visibility:"hidden",width:0,height:0,border:0,margin:0,background:"none"},m&&f.extend(p,{position:"absolute",left:"-999px",top:"-999px"});for(t in p)o.style[t]=p[t];o.appendChild(a),n=m||b,n.insertBefore(o,n.firstChild),k.appendChecked=i.checked,k.boxModel=a.offsetWidth===2,"zoom"in a.style&&(a.style.display="inline",a.style.zoom=1,k.inlineBlockNeedsLayout=a.offsetWidth===2,a.style.display="",a.innerHTML="<div style='width:4px;'></div>",k.shrinkWrapBlocks=a.offsetWidth!==2),a.innerHTML="<table><tr><td style='padding:0;border:0;display:none'></td><td>t</td></tr></table>",q=a.getElementsByTagName("td"),u=q[0].offsetHeight===0,q[0].style.display="",q[1].style.display="none",k.reliableHiddenOffsets=u&&q[0].offsetHeight===0,a.innerHTML="",c.defaultView&&c.defaultView.getComputedStyle&&(j=c.createElement("div"),j.style.width="0",j.style.marginRight="0",a.appendChild(j),k.reliableMarginRight=(parseInt((c.defaultView.getComputedStyle(j,null)||{marginRight:0}).marginRight,10)||0)===0);if(a.attachEvent)for(t in{submit:1,change:1,focusin:1})s="on"+t,u=s in a,u||(a.setAttribute(s,"return;"),u=typeof a[s]=="function"),k[t+"Bubbles"]=u;f(function(){var a,b,d,e,g,h,i=1,j="position:absolute;top:0;left:0;width:1px;height:1px;margin:0;",l="visibility:hidden;border:0;",n="style='"+j+"border:5px solid #000;padding:0;'",p="<div "+n+"><div></div></div>"+"<table "+n+" cellpadding='0' cellspacing='0'>"+"<tr><td></td></tr></table>";m=c.getElementsByTagName("body")[0];!m||(a=c.createElement("div"),a.style.cssText=l+"width:0;height:0;position:static;top:0;margin-top:"+i+"px",m.insertBefore(a,m.firstChild),o=c.createElement("div"),o.style.cssText=j+l,o.innerHTML=p,a.appendChild(o),b=o.firstChild,d=b.firstChild,g=b.nextSibling.firstChild.firstChild,h={doesNotAddBorder:d.offsetTop!==5,doesAddBorderForTableAndCells:g.offsetTop===5},d.style.position="fixed",d.style.top="20px",h.fixedPosition=d.offsetTop===20||d.offsetTop===15,d.style.position=d.style.top="",b.style.overflow="hidden",b.style.position="relative",h.subtractsBorderForOverflowNotVisible=d.offsetTop===-5,h.doesNotIncludeMarginInBodyOffset=m.offsetTop!==i,m.removeChild(a),o=a=null,f.extend(k,h))}),o.innerHTML="",n.removeChild(o),o=l=g=h=m=j=a=i=null;return k}(),f.boxModel=f.support.boxModel;var j=/^(?:\{.*\}|\[.*\])$/,k=/([A-Z])/g;f.extend({cache:{},uuid:0,expando:"jQuery"+(f.fn.jquery+Math.random()).replace(/\D/g,""),noData:{embed:!0,object:"clsid:D27CDB6E-AE6D-11cf-96B8-444553540000",applet:!0},hasData:function(a){a=a.nodeType?f.cache[a[f.expando]]:a[f.expando];return!!a&&!m(a)},data:function(a,c,d,e){if(!!f.acceptData(a)){var g,h,i,j=f.expando,k=typeof c=="string",l=a.nodeType,m=l?f.cache:a,n=l?a[f.expando]:a[f.expando]&&f.expando,o=c==="events";if((!n||!m[n]||!o&&!e&&!m[n].data)&&k&&d===b)return;n||(l?a[f.expando]=n=++f.uuid:n=f.expando),m[n]||(m[n]={},l||(m[n].toJSON=f.noop));if(typeof c=="object"||typeof c=="function")e?m[n]=f.extend(m[n],c):m[n].data=f.extend(m[n].data,c);g=h=m[n],e||(h.data||(h.data={}),h=h.data),d!==b&&(h[f.camelCase(c)]=d);if(o&&!h[c])return g.events;k?(i=h[c],i==null&&(i=h[f.camelCase(c)])):i=h;return i}},removeData:function(a,b,c){if(!!f.acceptData(a)){var d,e,g,h=f.expando,i=a.nodeType,j=i?f.cache:a,k=i?a[f.expando]:f.expando;if(!j[k])return;if(b){d=c?j[k]:j[k].data;if(d){f.isArray(b)?b=b:b in d?b=[b]:(b=f.camelCase(b),b in d?b=[b]:b=b.split(" "));for(e=0,g=b.length;e<g;e++)delete d[b[e]];if(!(c?m:f.isEmptyObject)(d))return}}if(!c){delete j[k].data;if(!m(j[k]))return}f.support.deleteExpando||!j.setInterval?delete j[k]:j[k]=null,i&&(f.support.deleteExpando?delete a[f.expando]:a.removeAttribute?a.removeAttribute(f.expando):a[f.expando]=null)}},_data:function(a,b,c){return f.data(a,b,c,!0)},acceptData:function(a){if(a.nodeName){var b=f.noData[a.nodeName.toLowerCase()];if(b)return b!==!0&&a.getAttribute("classid")===b}return!0}}),f.fn.extend({data:function(a,c){var d,e,g,h=null;if(typeof a=="undefined"){if(this.length){h=f.data(this[0]);if(this[0].nodeType===1&&!f._data(this[0],"parsedAttrs")){e=this[0].attributes;for(var i=0,j=e.length;i<j;i++)g=e[i].name,g.indexOf("data-")===0&&(g=f.camelCase(g.substring(5)),l(this[0],g,h[g]));f._data(this[0],"parsedAttrs",!0)}}return h}if(typeof a=="object")return this.each(function(){f.data(this,a)});d=a.split("."),d[1]=d[1]?"."+d[1]:"";if(c===b){h=this.triggerHandler("getData"+d[1]+"!",[d[0]]),h===b&&this.length&&(h=f.data(this[0],a),h=l(this[0],a,h));return h===b&&d[1]?this.data(d[0]):h}return this.each(function(){var b=f(this),e=[d[0],c];b.triggerHandler("setData"+d[1]+"!",e),f.data(this,a,c),b.triggerHandler("changeData"+d[1]+"!",e)})},removeData:function(a){return this.each(function(){f.removeData(this,a)})}}),f.extend({_mark:function(a,b){a&&(b=(b||"fx")+"mark",f._data(a,b,(f._data(a,b)||0)+1))},_unmark:function(a,b,c){a!==!0&&(c=b,b=a,a=!1);if(b){c=c||"fx";var d=c+"mark",e=a?0:(f._data(b,d)||1)-1;e?f._data(b,d,e):(f.removeData(b,d,!0),n(b,c,"mark"))}},queue:function(a,b,c){var d;if(a){b=(b||"fx")+"queue",d=f._data(a,b),c&&(!d||f.isArray(c)?d=f._data(a,b,f.makeArray(c)):d.push(c));return d||[]}},dequeue:function(a,b){b=b||"fx";var c=f.queue(a,b),d=c.shift(),e={};d==="inprogress"&&(d=c.shift()),d&&(b==="fx"&&c.unshift("inprogress"),f._data(a,b+".run",e),d.call(a,function(){f.dequeue(a,b)},e)),c.length||(f.removeData(a,b+"queue "+b+".run",!0),n(a,b,"queue"))}}),f.fn.extend({queue:function(a,c){typeof a!="string"&&(c=a,a="fx");if(c===b)return f.queue(this[0],a);return this.each(function(){var b=f.queue(this,a,c);a==="fx"&&b[0]!=="inprogress"&&f.dequeue(this,a)})},dequeue:function(a){return this.each(function(){f.dequeue(this,a)})},delay:function(a,b){a=f.fx?f.fx.speeds[a]||a:a,b=b||"fx";return this.queue(b,function(b,c){var d=setTimeout(b,a);c.stop=function(){clearTimeout(d)}})},clearQueue:function(a){return this.queue(a||"fx",[])},promise:function(a,c){function m(){--h||d.resolveWith(e,[e])}typeof a!="string"&&(c=a,a=b),a=a||"fx";var d=f.Deferred(),e=this,g=e.length,h=1,i=a+"defer",j=a+"queue",k=a+"mark",l;while(g--)if(l=f.data(e[g],i,b,!0)||(f.data(e[g],j,b,!0)||f.data(e[g],k,b,!0))&&f.data(e[g],i,f.Callbacks("once memory"),!0))h++,l.add(m);m();return d.promise()}});var o=/[\n\t\r]/g,p=/\s+/,q=/\r/g,r=/^(?:button|input)$/i,s=/^(?:button|input|object|select|textarea)$/i,t=/^a(?:rea)?$/i,u=/^(?:autofocus|autoplay|async|checked|controls|defer|disabled|hidden|loop|multiple|open|readonly|required|scoped|selected)$/i,v=f.support.getSetAttribute,w,x,y;f.fn.extend({attr:function(a,b){return f.access(this,a,b,!0,f.attr)},removeAttr:function(a){return this.each(function(){f.removeAttr(this,a)})},prop:function(a,b){return f.access(this,a,b,!0,f.prop)},removeProp:function(a){a=f.propFix[a]||a;return this.each(function(){try{this[a]=b,delete this[a]}catch(c){}})},addClass:function(a){var b,c,d,e,g,h,i;if(f.isFunction(a))return this.each(function(b){f(this).addClass(a.call(this,b,this.className))});if(a&&typeof a=="string"){b=a.split(p);for(c=0,d=this.length;c<d;c++){e=this[c];if(e.nodeType===1)if(!e.className&&b.length===1)e.className=a;else{g=" "+e.className+" ";for(h=0,i=b.length;h<i;h++)~g.indexOf(" "+b[h]+" ")||(g+=b[h]+" ");e.className=f.trim(g)}}}return this},removeClass:function(a){var c,d,e,g,h,i,j;if(f.isFunction(a))return this.each(function(b){f(this).removeClass(a.call(this,b,this.className))});if(a&&typeof a=="string"||a===b){c=(a||"").split(p);for(d=0,e=this.length;d<e;d++){g=this[d];if(g.nodeType===1&&g.className)if(a){h=(" "+g.className+" ").replace(o," ");for(i=0,j=c.length;i<j;i++)h=h.replace(" "+c[i]+" "," ");g.className=f.trim(h)}else g.className=""}}return this},toggleClass:function(a,b){var c=typeof a,d=typeof b=="boolean";if(f.isFunction(a))return this.each(function(c){f(this).toggleClass(a.call(this,c,this.className,b),b)});return this.each(function(){if(c==="string"){var e,g=0,h=f(this),i=b,j=a.split(p);while(e=j[g++])i=d?i:!h.hasClass(e),h[i?"addClass":"removeClass"](e)}else if(c==="undefined"||c==="boolean")this.className&&f._data(this,"__className__",this.className),this.className=this.className||a===!1?"":f._data(this,"__className__")||""})},hasClass:function(a){var b=" "+a+" ",c=0,d=this.length;for(;c<d;c++)if(this[c].nodeType===1&&(" "+this[c].className+" ").replace(o," ").indexOf(b)>-1)return!0;return!1},val:function(a){var c,d,e,g=this[0];if(!arguments.length){if(g){c=f.valHooks[g.nodeName.toLowerCase()]||f.valHooks[g.type];if(c&&"get"in c&&(d=c.get(g,"value"))!==b)return d;d=g.value;return typeof d=="string"?d.replace(q,""):d==null?"":d}return b}e=f.isFunction(a);return this.each(function(d){var g=f(this),h;if(this.nodeType===1){e?h=a.call(this,d,g.val()):h=a,h==null?h="":typeof h=="number"?h+="":f.isArray(h)&&(h=f.map(h,function(a){return a==null?"":a+""})),c=f.valHooks[this.nodeName.toLowerCase()]||f.valHooks[this.type];if(!c||!("set"in c)||c.set(this,h,"value")===b)this.value=h}})}}),f.extend({valHooks:{option:{get:function(a){var b=a.attributes.value;return!b||b.specified?a.value:a.text}},select:{get:function(a){var b,c,d,e,g=a.selectedIndex,h=[],i=a.options,j=a.type==="select-one";if(g<0)return null;c=j?g:0,d=j?g+1:i.length;for(;c<d;c++){e=i[c];if(e.selected&&(f.support.optDisabled?!e.disabled:e.getAttribute("disabled")===null)&&(!e.parentNode.disabled||!f.nodeName(e.parentNode,"optgroup"))){b=f(e).val();if(j)return b;h.push(b)}}if(j&&!h.length&&i.length)return f(i[g]).val();return h},set:function(a,b){var c=f.makeArray(b);f(a).find("option").each(function(){this.selected=f.inArray(f(this).val(),c)>=0}),c.length||(a.selectedIndex=-1);return c}}},attrFn:{val:!0,css:!0,html:!0,text:!0,data:!0,width:!0,height:!0,offset:!0},attr:function(a,c,d,e){var g,h,i,j=a.nodeType;if(!a||j===3||j===8||j===2)return b;if(e&&c in f.attrFn)return f(a)[c](d);if(!("getAttribute"in a))return f.prop(a,c,d);i=j!==1||!f.isXMLDoc(a),i&&(c=c.toLowerCase(),h=f.attrHooks[c]||(u.test(c)?x:w));if(d!==b){if(d===null){f.removeAttr(a,c);return b}if(h&&"set"in h&&i&&(g=h.set(a,d,c))!==b)return g;a.setAttribute(c,""+d);return d}if(h&&"get"in h&&i&&(g=h.get(a,c))!==null)return g;g=a.getAttribute(c);return g===null?b:g},removeAttr:function(a,b){var c,d,e,g,h=0;if(a.nodeType===1){d=(b||"").split(p),g=d.length;for(;h<g;h++)e=d[h].toLowerCase(),c=f.propFix[e]||e,f.attr(a,e,""),a.removeAttribute(v?e:c),u.test(e)&&c in a&&(a[c]=!1)}},attrHooks:{type:{set:function(a,b){if(r.test(a.nodeName)&&a.parentNode)f.error("type property can't be changed");else if(!f.support.radioValue&&b==="radio"&&f.nodeName(a,"input")){var c=a.value;a.setAttribute("type",b),c&&(a.value=c);return b}}},value:{get:function(a,b){if(w&&f.nodeName(a,"button"))return w.get(a,b);return b in a?a.value:null},set:function(a,b,c){if(w&&f.nodeName(a,"button"))return w.set(a,b,c);a.value=b}}},propFix:{tabindex:"tabIndex",readonly:"readOnly","for":"htmlFor","class":"className",maxlength:"maxLength",cellspacing:"cellSpacing",cellpadding:"cellPadding",rowspan:"rowSpan",colspan:"colSpan",usemap:"useMap",frameborder:"frameBorder",contenteditable:"contentEditable"},prop:function(a,c,d){var e,g,h,i=a.nodeType;if(!a||i===3||i===8||i===2)return b;h=i!==1||!f.isXMLDoc(a),h&&(c=f.propFix[c]||c,g=f.propHooks[c]);return d!==b?g&&"set"in g&&(e=g.set(a,d,c))!==b?e:a[c]=d:g&&"get"in g&&(e=g.get(a,c))!==null?e:a[c]},propHooks:{tabIndex:{get:function(a){var c=a.getAttributeNode("tabindex");return c&&c.specified?parseInt(c.value,10):s.test(a.nodeName)||t.test(a.nodeName)&&a.href?0:b}}}}),f.attrHooks.tabindex=f.propHooks.tabIndex,x={get:function(a,c){var d,e=f.prop(a,c);return e===!0||typeof e!="boolean"&&(d=a.getAttributeNode(c))&&d.nodeValue!==!1?c.toLowerCase():b},set:function(a,b,c){var d;b===!1?f.removeAttr(a,c):(d=f.propFix[c]||c,d in a&&(a[d]=!0),a.setAttribute(c,c.toLowerCase()));return c}},v||(y={name:!0,id:!0},w=f.valHooks.button={get:function(a,c){var d;d=a.getAttributeNode(c);return d&&(y[c]?d.nodeValue!=="":d.specified)?d.nodeValue:b},set:function(a,b,d){var e=a.getAttributeNode(d);e||(e=c.createAttribute(d),a.setAttributeNode(e));return e.nodeValue=b+""}},f.attrHooks.tabindex.set=w.set,f.each(["width","height"],function(a,b){f.attrHooks[b]=f.extend(f.attrHooks[b],{set:function(a,c){if(c===""){a.setAttribute(b,"auto");return c}}})}),f.attrHooks.contenteditable={get:w.get,set:function(a,b,c){b===""&&(b="false"),w.set(a,b,c)}}),f.support.hrefNormalized||f.each(["href","src","width","height"],function(a,c){f.attrHooks[c]=f.extend(f.attrHooks[c],{get:function(a){var d=a.getAttribute(c,2);return d===null?b:d}})}),f.support.style||(f.attrHooks.style={get:function(a){return a.style.cssText.toLowerCase()||b},set:function(a,b){return a.style.cssText=""+b}}),f.support.optSelected||(f.propHooks.selected=f.extend(f.propHooks.selected,{get:function(a){var b=a.parentNode;b&&(b.selectedIndex,b.parentNode&&b.parentNode.selectedIndex);return null}})),f.support.enctype||(f.propFix.enctype="encoding"),f.support.checkOn||f.each(["radio","checkbox"],function(){f.valHooks[this]={get:function(a){return a.getAttribute("value")===null?"on":a.value}}}),f.each(["radio","checkbox"],function(){f.valHooks[this]=f.extend(f.valHooks[this],{set:function(a,b){if(f.isArray(b))return a.checked=f.inArray(f(a).val(),b)>=0}})});var z=/\.(.*)$/,A=/^(?:textarea|input|select)$/i,B=/\./g,C=/ /g,D=/[^\w\s.|`]/g,E=/^([^\.]*)?(?:\.(.+))?$/,F=/\bhover(\.\S+)?/,G=/^key/,H=/^(?:mouse|contextmenu)|click/,I=/^(\w*)(?:#([\w\-]+))?(?:\.([\w\-]+))?$/,J=function(a){var b=I.exec(a);b&&
(b[1]=(b[1]||"").toLowerCase(),b[3]=b[3]&&new RegExp("(?:^|\\s)"+b[3]+"(?:\\s|$)"));return b},K=function(a,b){return(!b[1]||a.nodeName.toLowerCase()===b[1])&&(!b[2]||a.id===b[2])&&(!b[3]||b[3].test(a.className))},L=function(a){return f.event.special.hover?a:a.replace(F,"mouseenter$1 mouseleave$1")};f.event={add:function(a,c,d,e,g){var h,i,j,k,l,m,n,o,p,q,r,s;if(!(a.nodeType===3||a.nodeType===8||!c||!d||!(h=f._data(a)))){d.handler&&(p=d,d=p.handler),d.guid||(d.guid=f.guid++),j=h.events,j||(h.events=j={}),i=h.handle,i||(h.handle=i=function(a){return typeof f!="undefined"&&(!a||f.event.triggered!==a.type)?f.event.dispatch.apply(i.elem,arguments):b},i.elem=a),c=L(c).split(" ");for(k=0;k<c.length;k++){l=E.exec(c[k])||[],m=l[1],n=(l[2]||"").split(".").sort(),s=f.event.special[m]||{},m=(g?s.delegateType:s.bindType)||m,s=f.event.special[m]||{},o=f.extend({type:m,origType:l[1],data:e,handler:d,guid:d.guid,selector:g,namespace:n.join(".")},p),g&&(o.quick=J(g),!o.quick&&f.expr.match.POS.test(g)&&(o.isPositional=!0)),r=j[m];if(!r){r=j[m]=[],r.delegateCount=0;if(!s.setup||s.setup.call(a,e,n,i)===!1)a.addEventListener?a.addEventListener(m,i,!1):a.attachEvent&&a.attachEvent("on"+m,i)}s.add&&(s.add.call(a,o),o.handler.guid||(o.handler.guid=d.guid)),g?r.splice(r.delegateCount++,0,o):r.push(o),f.event.global[m]=!0}a=null}},global:{},remove:function(a,b,c,d){var e=f.hasData(a)&&f._data(a),g,h,i,j,k,l,m,n,o,p,q;if(!!e&&!!(m=e.events)){b=L(b||"").split(" ");for(g=0;g<b.length;g++){h=E.exec(b[g])||[],i=h[1],j=h[2];if(!i){j=j?"."+j:"";for(l in m)f.event.remove(a,l+j,c,d);return}n=f.event.special[i]||{},i=(d?n.delegateType:n.bindType)||i,p=m[i]||[],k=p.length,j=j?new RegExp("(^|\\.)"+j.split(".").sort().join("\\.(?:.*\\.)?")+"(\\.|$)"):null;if(c||j||d||n.remove)for(l=0;l<p.length;l++){q=p[l];if(!c||c.guid===q.guid)if(!j||j.test(q.namespace))if(!d||d===q.selector||d==="**"&&q.selector)p.splice(l--,1),q.selector&&p.delegateCount--,n.remove&&n.remove.call(a,q)}else p.length=0;p.length===0&&k!==p.length&&((!n.teardown||n.teardown.call(a,j)===!1)&&f.removeEvent(a,i,e.handle),delete m[i])}f.isEmptyObject(m)&&(o=e.handle,o&&(o.elem=null),f.removeData(a,["events","handle"],!0))}},customEvent:{getData:!0,setData:!0,changeData:!0},trigger:function(c,d,e,g){if(!e||e.nodeType!==3&&e.nodeType!==8){var h=c.type||c,i=[],j,k,l,m,n,o,p,q,r,s;h.indexOf("!")>=0&&(h=h.slice(0,-1),k=!0),h.indexOf(".")>=0&&(i=h.split("."),h=i.shift(),i.sort());if((!e||f.event.customEvent[h])&&!f.event.global[h])return;c=typeof c=="object"?c[f.expando]?c:new f.Event(h,c):new f.Event(h),c.type=h,c.isTrigger=!0,c.exclusive=k,c.namespace=i.join("."),c.namespace_re=c.namespace?new RegExp("(^|\\.)"+i.join("\\.(?:.*\\.)?")+"(\\.|$)"):null,o=h.indexOf(":")<0?"on"+h:"",(g||!e)&&c.preventDefault();if(!e){j=f.cache;for(l in j)j[l].events&&j[l].events[h]&&f.event.trigger(c,d,j[l].handle.elem,!0);return}c.result=b,c.target||(c.target=e),d=d!=null?f.makeArray(d):[],d.unshift(c),p=f.event.special[h]||{};if(p.trigger&&p.trigger.apply(e,d)===!1)return;r=[[e,p.bindType||h]];if(!g&&!p.noBubble&&!f.isWindow(e)){s=p.delegateType||h,n=null;for(m=e.parentNode;m;m=m.parentNode)r.push([m,s]),n=m;n&&n===e.ownerDocument&&r.push([n.defaultView||n.parentWindow||a,s])}for(l=0;l<r.length;l++){m=r[l][0],c.type=r[l][1],q=(f._data(m,"events")||{})[c.type]&&f._data(m,"handle"),q&&q.apply(m,d),q=o&&m[o],q&&f.acceptData(m)&&q.apply(m,d);if(c.isPropagationStopped())break}c.type=h,c.isDefaultPrevented()||(!p._default||p._default.apply(e.ownerDocument,d)===!1)&&(h!=="click"||!f.nodeName(e,"a"))&&f.acceptData(e)&&o&&e[h]&&(h!=="focus"&&h!=="blur"||c.target.offsetWidth!==0)&&!f.isWindow(e)&&(n=e[o],n&&(e[o]=null),f.event.triggered=h,e[h](),f.event.triggered=b,n&&(e[o]=n));return c.result}},dispatch:function(c){c=f.event.fix(c||a.event);var d=(f._data(this,"events")||{})[c.type]||[],e=d.delegateCount,g=[].slice.call(arguments,0),h=!c.exclusive&&!c.namespace,i=(f.event.special[c.type]||{}).handle,j=[],k,l,m,n,o,p,q,r,s,t,u;g[0]=c,c.delegateTarget=this;if(e&&!c.target.disabled&&(!c.button||c.type!=="click"))for(m=c.target;m!=this;m=m.parentNode||this){o={},q=[];for(k=0;k<e;k++)r=d[k],s=r.selector,t=o[s],r.isPositional?t=(t||(o[s]=f(s))).index(m)>=0:t===b&&(t=o[s]=r.quick?K(m,r.quick):f(m).is(s)),t&&q.push(r);q.length&&j.push({elem:m,matches:q})}d.length>e&&j.push({elem:this,matches:d.slice(e)});for(k=0;k<j.length&&!c.isPropagationStopped();k++){p=j[k],c.currentTarget=p.elem;for(l=0;l<p.matches.length&&!c.isImmediatePropagationStopped();l++){r=p.matches[l];if(h||!c.namespace&&!r.namespace||c.namespace_re&&c.namespace_re.test(r.namespace))c.data=r.data,c.handleObj=r,n=(i||r.handler).apply(p.elem,g),n!==b&&(c.result=n,n===!1&&(c.preventDefault(),c.stopPropagation()))}}return c.result},props:"attrChange attrName relatedNode srcElement altKey bubbles cancelable ctrlKey currentTarget eventPhase metaKey relatedTarget shiftKey target timeStamp view which".split(" "),fixHooks:{},keyHooks:{props:"char charCode key keyCode".split(" "),filter:function(a,b){a.which==null&&(a.which=b.charCode!=null?b.charCode:b.keyCode);return a}},mouseHooks:{props:"button buttons clientX clientY fromElement offsetX offsetY pageX pageY screenX screenY toElement wheelDelta".split(" "),filter:function(a,d){var e,f,g,h=d.button,i=d.fromElement;a.pageX==null&&d.clientX!=null&&(e=a.target.ownerDocument||c,f=e.documentElement,g=e.body,a.pageX=d.clientX+(f&&f.scrollLeft||g&&g.scrollLeft||0)-(f&&f.clientLeft||g&&g.clientLeft||0),a.pageY=d.clientY+(f&&f.scrollTop||g&&g.scrollTop||0)-(f&&f.clientTop||g&&g.clientTop||0)),!a.relatedTarget&&i&&(a.relatedTarget=i===a.target?d.toElement:i),!a.which&&h!==b&&(a.which=h&1?1:h&2?3:h&4?2:0);return a}},fix:function(a){if(a[f.expando])return a;var d,e,g=a,h=f.event.fixHooks[a.type]||{},i=h.props?this.props.concat(h.props):this.props;a=f.Event(g);for(d=i.length;d;)e=i[--d],a[e]=g[e];a.target||(a.target=g.srcElement||c),a.target.nodeType===3&&(a.target=a.target.parentNode),a.metaKey===b&&(a.metaKey=a.ctrlKey);return h.filter?h.filter(a,g):a},special:{ready:{setup:f.bindReady},focus:{delegateType:"focusin",noBubble:!0},blur:{delegateType:"focusout",noBubble:!0},beforeunload:{setup:function(a,b,c){f.isWindow(this)&&(this.onbeforeunload=c)},teardown:function(a,b){this.onbeforeunload===b&&(this.onbeforeunload=null)}}},simulate:function(a,b,c,d){var e=f.extend(new f.Event,c,{type:a,isSimulated:!0,originalEvent:{}});d?f.event.trigger(e,null,b):f.event.dispatch.call(b,e),e.isDefaultPrevented()&&c.preventDefault()}},f.event.handle=f.event.dispatch,f.removeEvent=c.removeEventListener?function(a,b,c){a.removeEventListener&&a.removeEventListener(b,c,!1)}:function(a,b,c){a.detachEvent&&a.detachEvent("on"+b,c)},f.Event=function(a,b){if(!(this instanceof f.Event))return new f.Event(a,b);a&&a.type?(this.originalEvent=a,this.type=a.type,this.isDefaultPrevented=a.defaultPrevented||a.returnValue===!1||a.getPreventDefault&&a.getPreventDefault()?N:M):this.type=a,b&&f.extend(this,b),this.timeStamp=a&&a.timeStamp||f.now(),this[f.expando]=!0},f.Event.prototype={preventDefault:function(){this.isDefaultPrevented=N;var a=this.originalEvent;!a||(a.preventDefault?a.preventDefault():a.returnValue=!1)},stopPropagation:function(){this.isPropagationStopped=N;var a=this.originalEvent;!a||(a.stopPropagation&&a.stopPropagation(),a.cancelBubble=!0)},stopImmediatePropagation:function(){this.isImmediatePropagationStopped=N,this.stopPropagation()},isDefaultPrevented:M,isPropagationStopped:M,isImmediatePropagationStopped:M},f.each({mouseenter:"mouseover",mouseleave:"mouseout"},function(a,b){f.event.special[a]=f.event.special[b]={delegateType:b,bindType:b,handle:function(a){var b=this,c=a.relatedTarget,d=a.handleObj,e=d.selector,g,h;if(!c||d.origType===a.type||c!==b&&!f.contains(b,c))g=a.type,a.type=d.origType,h=d.handler.apply(this,arguments),a.type=g;return h}}}),f.support.submitBubbles||(f.event.special.submit={setup:function(){if(f.nodeName(this,"form"))return!1;f.event.add(this,"click._submit keypress._submit",function(a){var c=a.target,d=f.nodeName(c,"input")||f.nodeName(c,"button")?c.form:b;d&&!d._submit_attached&&(f.event.add(d,"submit._submit",function(a){this.parentNode&&f.event.simulate("submit",this.parentNode,a,!0)}),d._submit_attached=!0)})},teardown:function(){if(f.nodeName(this,"form"))return!1;f.event.remove(this,"._submit")}}),f.support.changeBubbles||(f.event.special.change={setup:function(){if(A.test(this.nodeName)){if(this.type==="checkbox"||this.type==="radio")f.event.add(this,"propertychange._change",function(a){a.originalEvent.propertyName==="checked"&&(this._just_changed=!0)}),f.event.add(this,"click._change",function(a){this._just_changed&&(this._just_changed=!1,f.event.simulate("change",this,a,!0))});return!1}f.event.add(this,"beforeactivate._change",function(a){var b=a.target;A.test(b.nodeName)&&!b._change_attached&&(f.event.add(b,"change._change",function(a){this.parentNode&&!a.isSimulated&&f.event.simulate("change",this.parentNode,a,!0)}),b._change_attached=!0)})},handle:function(a){var b=a.target;if(this!==b||a.isSimulated||a.isTrigger||b.type!=="radio"&&b.type!=="checkbox")return a.handleObj.handler.apply(this,arguments)},teardown:function(){f.event.remove(this,"._change");return A.test(this.nodeName)}}),f.support.focusinBubbles||f.each({focus:"focusin",blur:"focusout"},function(a,b){var d=0,e=function(a){f.event.simulate(b,a.target,f.event.fix(a),!0)};f.event.special[b]={setup:function(){d++===0&&c.addEventListener(a,e,!0)},teardown:function(){--d===0&&c.removeEventListener(a,e,!0)}}}),f.fn.extend({on:function(a,c,d,e,g){var h,i;if(typeof a=="object"){typeof c!="string"&&(d=c,c=b);for(i in a)this.on(i,c,d,a[i],g);return this}d==null&&e==null?(e=c,d=c=b):e==null&&(typeof c=="string"?(e=d,d=b):(e=d,d=c,c=b));if(e===!1)e=M;else if(!e)return this;g===1&&(h=e,e=function(a){f().off(a);return h.apply(this,arguments)},e.guid=h.guid||(h.guid=f.guid++));return this.each(function(){f.event.add(this,a,e,d,c)})},one:function(a,b,c,d){return this.on.call(this,a,b,c,d,1)},off:function(a,c,d){if(a&&a.preventDefault&&a.handleObj){var e=a.handleObj;f(a.delegateTarget).off(e.namespace?e.type+"."+e.namespace:e.type,e.selector,e.handler);return this}if(typeof a=="object"){for(var g in a)this.off(g,c,a[g]);return this}if(c===!1||typeof c=="function")d=c,c=b;d===!1&&(d=M);return this.each(function(){f.event.remove(this,a,d,c)})},bind:function(a,b,c){return this.on(a,null,b,c)},unbind:function(a,b){return this.off(a,null,b)},live:function(a,b,c){f(this.context).on(a,this.selector,b,c);return this},die:function(a,b){f(this.context).off(a,this.selector||"**",b);return this},delegate:function(a,b,c,d){return this.on(b,a,c,d)},undelegate:function(a,b,c){return arguments.length==1?this.off(a,"**"):this.off(b,a,c)},trigger:function(a,b){return this.each(function(){f.event.trigger(a,b,this)})},triggerHandler:function(a,b){if(this[0])return f.event.trigger(a,b,this[0],!0)},toggle:function(a){var b=arguments,c=a.guid||f.guid++,d=0,e=function(c){var e=(f._data(this,"lastToggle"+a.guid)||0)%d;f._data(this,"lastToggle"+a.guid,e+1),c.preventDefault();return b[e].apply(this,arguments)||!1};e.guid=c;while(d<b.length)b[d++].guid=c;return this.click(e)},hover:function(a,b){return this.mouseenter(a).mouseleave(b||a)}}),f.each("blur focus focusin focusout load resize scroll unload click dblclick mousedown mouseup mousemove mouseover mouseout mouseenter mouseleave change select submit keydown keypress keyup error contextmenu".split(" "),function(a,b){f.fn[b]=function(a,c){c==null&&(c=a,a=null);return arguments.length>0?this.bind(b,a,c):this.trigger(b)},f.attrFn&&(f.attrFn[b]=!0),G.test(b)&&(f.event.fixHooks[b]=f.event.keyHooks),H.test(b)&&(f.event.fixHooks[b]=f.event.mouseHooks)}),function(){function x(a,b,c,e,f,g){for(var h=0,i=e.length;h<i;h++){var j=e[h];if(j){var k=!1;j=j[a];while(j){if(j[d]===c){k=e[j.sizset];break}if(j.nodeType===1){g||(j[d]=c,j.sizset=h);if(typeof b!="string"){if(j===b){k=!0;break}}else if(m.filter(b,[j]).length>0){k=j;break}}j=j[a]}e[h]=k}}}function w(a,b,c,e,f,g){for(var h=0,i=e.length;h<i;h++){var j=e[h];if(j){var k=!1;j=j[a];while(j){if(j[d]===c){k=e[j.sizset];break}j.nodeType===1&&!g&&(j[d]=c,j.sizset=h);if(j.nodeName.toLowerCase()===b){k=j;break}j=j[a]}e[h]=k}}}var a=/((?:\((?:\([^()]+\)|[^()]+)+\)|\[(?:\[[^\[\]]*\]|['"][^'"]*['"]|[^\[\]'"]+)+\]|\\.|[^ >+~,(\[\\]+)+|[>+~])(\s*,\s*)?((?:.|\r|\n)*)/g,d="sizcache"+(Math.random()+"").replace(".",""),e=0,g=Object.prototype.toString,h=!1,i=!0,j=/\\/g,k=/\r\n/g,l=/\W/;[0,0].sort(function(){i=!1;return 0});var m=function(b,d,e,f){e=e||[],d=d||c;var h=d;if(d.nodeType!==1&&d.nodeType!==9)return[];if(!b||typeof b!="string")return e;var i,j,k,l,n,q,r,t,u=!0,v=m.isXML(d),w=[],x=b;do{a.exec(""),i=a.exec(x);if(i){x=i[3],w.push(i[1]);if(i[2]){l=i[3];break}}}while(i);if(w.length>1&&p.exec(b))if(w.length===2&&o.relative[w[0]])j=y(w[0]+w[1],d,f);else{j=o.relative[w[0]]?[d]:m(w.shift(),d);while(w.length)b=w.shift(),o.relative[b]&&(b+=w.shift()),j=y(b,j,f)}else{!f&&w.length>1&&d.nodeType===9&&!v&&o.match.ID.test(w[0])&&!o.match.ID.test(w[w.length-1])&&(n=m.find(w.shift(),d,v),d=n.expr?m.filter(n.expr,n.set)[0]:n.set[0]);if(d){n=f?{expr:w.pop(),set:s(f)}:m.find(w.pop(),w.length===1&&(w[0]==="~"||w[0]==="+")&&d.parentNode?d.parentNode:d,v),j=n.expr?m.filter(n.expr,n.set):n.set,w.length>0?k=s(j):u=!1;while(w.length)q=w.pop(),r=q,o.relative[q]?r=w.pop():q="",r==null&&(r=d),o.relative[q](k,r,v)}else k=w=[]}k||(k=j),k||m.error(q||b);if(g.call(k)==="[object Array]")if(!u)e.push.apply(e,k);else if(d&&d.nodeType===1)for(t=0;k[t]!=null;t++)k[t]&&(k[t]===!0||k[t].nodeType===1&&m.contains(d,k[t]))&&e.push(j[t]);else for(t=0;k[t]!=null;t++)k[t]&&k[t].nodeType===1&&e.push(j[t]);else s(k,e);l&&(m(l,h,e,f),m.uniqueSort(e));return e};m.uniqueSort=function(a){if(u){h=i,a.sort(u);if(h)for(var b=1;b<a.length;b++)a[b]===a[b-1]&&a.splice(b--,1)}return a},m.matches=function(a,b){return m(a,null,null,b)},m.matchesSelector=function(a,b){return m(b,null,null,[a]).length>0},m.find=function(a,b,c){var d,e,f,g,h,i;if(!a)return[];for(e=0,f=o.order.length;e<f;e++){h=o.order[e];if(g=o.leftMatch[h].exec(a)){i=g[1],g.splice(1,1);if(i.substr(i.length-1)!=="\\"){g[1]=(g[1]||"").replace(j,""),d=o.find[h](g,b,c);if(d!=null){a=a.replace(o.match[h],"");break}}}}d||(d=typeof b.getElementsByTagName!="undefined"?b.getElementsByTagName("*"):[]);return{set:d,expr:a}},m.filter=function(a,c,d,e){var f,g,h,i,j,k,l,n,p,q=a,r=[],s=c,t=c&&c[0]&&m.isXML(c[0]);while(a&&c.length){for(h in o.filter)if((f=o.leftMatch[h].exec(a))!=null&&f[2]){k=o.filter[h],l=f[1],g=!1,f.splice(1,1);if(l.substr(l.length-1)==="\\")continue;s===r&&(r=[]);if(o.preFilter[h]){f=o.preFilter[h](f,s,d,r,e,t);if(!f)g=i=!0;else if(f===!0)continue}if(f)for(n=0;(j=s[n])!=null;n++)j&&(i=k(j,f,n,s),p=e^i,d&&i!=null?p?g=!0:s[n]=!1:p&&(r.push(j),g=!0));if(i!==b){d||(s=r),a=a.replace(o.match[h],"");if(!g)return[];break}}if(a===q)if(g==null)m.error(a);else break;q=a}return s},m.error=function(a){throw"Syntax error, unrecognized expression: "+a};var n=m.getText=function(a){var b,c,d=a.nodeType,e="";if(d){if(d===1){if(typeof a.textContent=="string")return a.textContent;if(typeof a.innerText=="string")return a.innerText.replace(k,"");for(a=a.firstChild;a;a=a.nextSibling)e+=n(a)}else if(d===3||d===4)return a.nodeValue}else for(b=0;c=a[b];b++)c.nodeType!==8&&(e+=n(c));return e},o=m.selectors={order:["ID","NAME","TAG"],match:{ID:/#((?:[\w\u00c0-\uFFFF\-]|\\.)+)/,CLASS:/\.((?:[\w\u00c0-\uFFFF\-]|\\.)+)/,NAME:/\[name=['"]*((?:[\w\u00c0-\uFFFF\-]|\\.)+)['"]*\]/,ATTR:/\[\s*((?:[\w\u00c0-\uFFFF\-]|\\.)+)\s*(?:(\S?=)\s*(?:(['"])(.*?)\3|(#?(?:[\w\u00c0-\uFFFF\-]|\\.)*)|)|)\s*\]/,TAG:/^((?:[\w\u00c0-\uFFFF\*\-]|\\.)+)/,CHILD:/:(only|nth|last|first)-child(?:\(\s*(even|odd|(?:[+\-]?\d+|(?:[+\-]?\d*)?n\s*(?:[+\-]\s*\d+)?))\s*\))?/,POS:/:(nth|eq|gt|lt|first|last|even|odd)(?:\((\d*)\))?(?=[^\-]|$)/,PSEUDO:/:((?:[\w\u00c0-\uFFFF\-]|\\.)+)(?:\((['"]?)((?:\([^\)]+\)|[^\(\)]*)+)\2\))?/},leftMatch:{},attrMap:{"class":"className","for":"htmlFor"},attrHandle:{href:function(a){return a.getAttribute("href")},type:function(a){return a.getAttribute("type")}},relative:{"+":function(a,b){var c=typeof b=="string",d=c&&!l.test(b),e=c&&!d;d&&(b=b.toLowerCase());for(var f=0,g=a.length,h;f<g;f++)if(h=a[f]){while((h=h.previousSibling)&&h.nodeType!==1);a[f]=e||h&&h.nodeName.toLowerCase()===b?h||!1:h===b}e&&m.filter(b,a,!0)},">":function(a,b){var c,d=typeof b=="string",e=0,f=a.length;if(d&&!l.test(b)){b=b.toLowerCase();for(;e<f;e++){c=a[e];if(c){var g=c.parentNode;a[e]=g.nodeName.toLowerCase()===b?g:!1}}}else{for(;e<f;e++)c=a[e],c&&(a[e]=d?c.parentNode:c.parentNode===b);d&&m.filter(b,a,!0)}},"":function(a,b,c){var d,f=e++,g=x;typeof b=="string"&&!l.test(b)&&(b=b.toLowerCase(),d=b,g=w),g("parentNode",b,f,a,d,c)},"~":function(a,b,c){var d,f=e++,g=x;typeof b=="string"&&!l.test(b)&&(b=b.toLowerCase(),d=b,g=w),g("previousSibling",b,f,a,d,c)}},find:{ID:function(a,b,c){if(typeof b.getElementById!="undefined"&&!c){var d=b.getElementById(a[1]);return d&&d.parentNode?[d]:[]}},NAME:function(a,b){if(typeof b.getElementsByName!="undefined"){var c=[],d=b.getElementsByName(a[1]);for(var e=0,f=d.length;e<f;e++)d[e].getAttribute("name")===a[1]&&c.push(d[e]);return c.length===0?null:c}},TAG:function(a,b){if(typeof b.getElementsByTagName!="undefined")return b.getElementsByTagName(a[1])}},preFilter:{CLASS:function(a,b,c,d,e,f){a=" "+a[1].replace(j,"")+" ";if(f)return a;for(var g=0,h;(h=b[g])!=null;g++)h&&(e^(h.className&&(" "+h.className+" ").replace(/[\t\n\r]/g," ").indexOf(a)>=0)?c||d.push(h):c&&(b[g]=!1));return!1},ID:function(a){return a[1].replace(j,"")},TAG:function(a,b){return a[1].replace(j,"").toLowerCase()},CHILD:function(a){if(a[1]==="nth"){a[2]||m.error(a[0]),a[2]=a[2].replace(/^\+|\s*/g,"");var b=/(-?)(\d*)(?:n([+\-]?\d*))?/.exec(a[2]==="even"&&"2n"||a[2]==="odd"&&"2n+1"||!/\D/.test(a[2])&&"0n+"+a[2]||a[2]);a[2]=b[1]+(b[2]||1)-0,a[3]=b[3]-0}else a[2]&&m.error(a[0]);a[0]=e++;return a},ATTR:function(a,b,c,d,e,f){var g=a[1]=a[1].replace(j,"");!f&&o.attrMap[g]&&(a[1]=o.attrMap[g]),a[4]=(a[4]||a[5]||"").replace(j,""),a[2]==="~="&&(a[4]=" "+a[4]+" ");return a},PSEUDO:function(b,c,d,e,f){if(b[1]==="not")if((a.exec(b[3])||"").length>1||/^\w/.test(b[3]))b[3]=m(b[3],null,null,c);else{var g=m.filter(b[3],c,d,!0^f);d||e.push.apply(e,g);return!1}else if(o.match.POS.test(b[0])||o.match.CHILD.test(b[0]))return!0;return b},POS:function(a){a.unshift(!0);return a}},filters:{enabled:function(a){return a.disabled===!1&&a.type!=="hidden"},disabled:function(a){return a.disabled===!0},checked:function(a){return a.checked===!0},selected:function(a){a.parentNode&&a.parentNode.selectedIndex;return a.selected===!0},parent:function(a){return!!a.firstChild},empty:function(a){return!a.firstChild},has:function(a,b,c){return!!m(c[3],a).length},header:function(a){return/h\d/i.test(a.nodeName)},text:function(a){var b=a.getAttribute("type"),c=a.type;return a.nodeName.toLowerCase()==="input"&&"text"===c&&(b===c||b===null)},radio:function(a){return a.nodeName.toLowerCase()==="input"&&"radio"===a.type},checkbox:function(a){return a.nodeName.toLowerCase()==="input"&&"checkbox"===a.type},file:function(a){return a.nodeName.toLowerCase()==="input"&&"file"===a.type},password:function(a){return a.nodeName.toLowerCase()==="input"&&"password"===a.type},submit:function(a){var b=a.nodeName.toLowerCase();return(b==="input"||b==="button")&&"submit"===a.type},image:function(a){return a.nodeName.toLowerCase()==="input"&&"image"===a.type},reset:function(a){var b=a.nodeName.toLowerCase();return(b==="input"||b==="button")&&"reset"===a.type},button:function(a){var b=a.nodeName.toLowerCase();return b==="input"&&"button"===a.type||b==="button"},input:function(a){return/input|select|textarea|button/i.test(a.nodeName)},focus:function(a){return a===a.ownerDocument.activeElement}},setFilters:{first:function(a,b){return b===0},last:function(a,b,c,d){return b===d.length-1},even:function(a,b){return b%2===0},odd:function(a,b){return b%2===1},lt:function(a,b,c){return b<c[3]-0},gt:function(a,b,c){return b>c[3]-0},nth:function(a,b,c){return c[3]-0===b},eq:function(a,b,c){return c[3]-0===b}},filter:{PSEUDO:function(a,b,c,d){var e=b[1],f=o.filters[e];if(f)return f(a,c,b,d);if(e==="contains")return(a.textContent||a.innerText||n([a])||"").indexOf(b[3])>=0;if(e==="not"){var g=b[3];for(var h=0,i=g.length;h<i;h++)if(g[h]===a)return!1;return!0}m.error(e)},CHILD:function(a,b){var c,e,f,g,h,i,j,k=b[1],l=a;switch(k){case"only":case"first":while(l=l.previousSibling)if(l.nodeType===1)return!1;if(k==="first")return!0;l=a;case"last":while(l=l.nextSibling)if(l.nodeType===1)return!1;return!0;case"nth":c=b[2],e=b[3];if(c===1&&e===0)return!0;f=b[0],g=a.parentNode;if(g&&(g[d]!==f||!a.nodeIndex)){i=0;for(l=g.firstChild;l;l=l.nextSibling)l.nodeType===1&&(l.nodeIndex=++i);g[d]=f}j=a.nodeIndex-e;return c===0?j===0:j%c===0&&j/c>=0}},ID:function(a,b){return a.nodeType===1&&a.getAttribute("id")===b},TAG:function(a,b){return b==="*"&&a.nodeType===1||!!a.nodeName&&a.nodeName.toLowerCase()===b},CLASS:function(a,b){return(" "+(a.className||a.getAttribute("class"))+" ").indexOf(b)>-1},ATTR:function(a,b){var c=b[1],d=m.attr?m.attr(a,c):o.attrHandle[c]?o.attrHandle[c](a):a[c]!=null?a[c]:a.getAttribute(c),e=d+"",f=b[2],g=b[4];return d==null?f==="!=":!f&&m.attr?d!=null:f==="="?e===g:f==="*="?e.indexOf(g)>=0:f==="~="?(" "+e+" ").indexOf(g)>=0:g?f==="!="?e!==g:f==="^="?e.indexOf(g)===0:f==="$="?e.substr(e.length-g.length)===g:f==="|="?e===g||e.substr(0,g.length+1)===g+"-":!1:e&&d!==!1},POS:function(a,b,c,d){var e=b[2],f=o.setFilters[e];if(f)return f(a,c,b,d)}}},p=o.match.POS,q=function(a,b){return"\\"+(b-0+1)};for(var r in o.match)o.match[r]=new RegExp(o.match[r].source+/(?![^\[]*\])(?![^\(]*\))/.source),o.leftMatch[r]=new RegExp(/(^(?:.|\r|\n)*?)/.source+o.match[r].source.replace(/\\(\d+)/g,q));var s=function(a,b){a=Array.prototype.slice.call(a,0);if(b){b.push.apply(b,a);return b}return a};try{Array.prototype.slice.call(c.documentElement.childNodes,0)[0].nodeType}catch(t){s=function(a,b){var c=0,d=b||[];if(g.call(a)==="[object Array]")Array.prototype.push.apply(d,a);else if(typeof a.length=="number")for(var e=a.length;c<e;c++)d.push(a[c]);else for(;a[c];c++)d.push(a[c]);return d}}var u,v;c.documentElement.compareDocumentPosition?u=function(a,b){if(a===b){h=!0;return 0}if(!a.compareDocumentPosition||!b.compareDocumentPosition)return a.compareDocumentPosition?-1:1;return a.compareDocumentPosition(b)&4?-1:1}:(u=function(a,b){if(a===b){h=!0;return 0}if(a.sourceIndex&&b.sourceIndex)return a.sourceIndex-b.sourceIndex;var c,d,e=[],f=[],g=a.parentNode,i=b.parentNode,j=g;if(g===i)return v(a,b);if(!g)return-1;if(!i)return 1;while(j)e.unshift(j),j=j.parentNode;j=i;while(j)f.unshift(j),j=j.parentNode;c=e.length,d=f.length;for(var k=0;k<c&&k<d;k++)if(e[k]!==f[k])return v(e[k],f[k]);return k===c?v(a,f[k],-1):v(e[k],b,1)},v=function(a,b,c){if(a===b)return c;var d=a.nextSibling;while(d){if(d===b)return-1;d=d.nextSibling}return 1}),function(){var a=c.createElement("div"),d="script"+(new Date).getTime(),e=c.documentElement;a.innerHTML="<a name='"+d+"'/>",e.insertBefore(a,e.firstChild),c.getElementById(d)&&(o.find.ID=function(a,c,d){if(typeof c.getElementById!="undefined"&&!d){var e=c.getElementById(a[1]);return e?e.id===a[1]||typeof e.getAttributeNode!="undefined"&&e.getAttributeNode("id").nodeValue===a[1]?[e]:b:[]}},o.filter.ID=function(a,b){var c=typeof a.getAttributeNode!="undefined"&&a.getAttributeNode("id");return a.nodeType===1&&c&&c.nodeValue===b}),e.removeChild(a),e=a=null}(),function(){var a=c.createElement("div");a.appendChild(c.createComment("")),a.getElementsByTagName("*").length>0&&(o.find.TAG=function(a,b){var c=b.getElementsByTagName(a[1]);if(a[1]==="*"){var d=[];for(var e=0;c[e];e++)c[e].nodeType===1&&d.push(c[e]);c=d}return c}),a.innerHTML="<a href='#'></a>",a.firstChild&&typeof a.firstChild.getAttribute!="undefined"&&a.firstChild.getAttribute("href")!=="#"&&(o.attrHandle.href=function(a){return a.getAttribute("href",2)}),a=null}(),c.querySelectorAll&&function(){var a=m,b=c.createElement("div"),d="__sizzle__";b.innerHTML="<p class='TEST'></p>";if(!b.querySelectorAll||b.querySelectorAll(".TEST").length!==0){m=function(b,e,f,g){e=e||c;if(!g&&!m.isXML(e)){var h=/^(\w+$)|^\.([\w\-]+$)|^#([\w\-]+$)/.exec(b);if(h&&(e.nodeType===1||e.nodeType===9)){if(h[1])return s(e.getElementsByTagName(b),f);if(h[2]&&o.find.CLASS&&e.getElementsByClassName)return s(e.getElementsByClassName(h[2]),f)}if(e.nodeType===9){if(b==="body"&&e.body)return s([e.body],f);if(h&&h[3]){var i=e.getElementById(h[3]);if(!i||!i.parentNode)return s([],f);if(i.id===h[3])return s([i],f)}try{return s(e.querySelectorAll(b),f)}catch(j){}}else if(e.nodeType===1&&e.nodeName.toLowerCase()!=="object"){var k=e,l=e.getAttribute("id"),n=l||d,p=e.parentNode,q=/^\s*[+~]/.test(b);l?n=n.replace(/'/g,"\\$&"):e.setAttribute("id",n),q&&p&&(e=e.parentNode);try{if(!q||p)return s(e.querySelectorAll("[id='"+n+"'] "+b),f)}catch(r){}finally{l||k.removeAttribute("id")}}}return a(b,e,f,g)};for(var e in a)m[e]=a[e];b=null}}(),function(){var a=c.documentElement,b=a.matchesSelector||a.mozMatchesSelector||a.webkitMatchesSelector||a.msMatchesSelector;if(b){var d=!b.call(c.createElement("div"),"div"),e=!1;try{b.call(c.documentElement,"[test!='']:sizzle")}catch(f){e=!0}m.matchesSelector=function(a,c){c=c.replace(/\=\s*([^'"\]]*)\s*\]/g,"='$1']");if(!m.isXML(a))try{if(e||!o.match.PSEUDO.test(c)&&!/!=/.test(c)){var f=b.call(a,c);if(f||!d||a.document&&a.document.nodeType!==11)return f}}catch(g){}return m(c,null,null,[a]).length>0}}}(),function(){var a=c.createElement("div");a.innerHTML="<div class='test e'></div><div class='test'></div>";if(!!a.getElementsByClassName&&a.getElementsByClassName("e").length!==0){a.lastChild.className="e";if(a.getElementsByClassName("e").length===1)return;o.order.splice(1,0,"CLASS"),o.find.CLASS=function(a,b,c){if(typeof b.getElementsByClassName!="undefined"&&!c)return b.getElementsByClassName(a[1])},a=null}}(),c.documentElement.contains?m.contains=function(a,b){return a!==b&&(a.contains?a.contains(b):!0)}:c.documentElement.compareDocumentPosition?m.contains=function(a,b){return!!(a.compareDocumentPosition(b)&16)}:m.contains=function(){return!1},m.isXML=function(a){var b=(a?a.ownerDocument||a:0).documentElement;return b?b.nodeName!=="HTML":!1};var y=function(a,b,c){var d,e=[],f="",g=b.nodeType?[b]:b;while(d=o.match.PSEUDO.exec(a))f+=d[0],a=a.replace(o.match.PSEUDO,"");a=o.relative[a]?a+"*":a;for(var h=0,i=g.length;h<i;h++)m(a,g[h],e,c);return m.filter(f,e)};m.attr=f.attr,m.selectors.attrMap={},f.find=m,f.expr=m.selectors,f.expr[":"]=f.expr.filters,f.unique=m.uniqueSort,f.text=m.getText,f.isXMLDoc=m.isXML,f.contains=m.contains}();var O=/Until$/,P=/^(?:parents|prevUntil|prevAll)/,Q=/,/,R=/^.[^:#\[\.,]*$/,S=Array.prototype.slice,T=f.expr.match.POS,U={children:!0,contents:!0,next:!0,prev:!0};f.fn.extend({find:function(a){var b=this,c,d;if(typeof a!="string")return f(a).filter(function(){for(c=0,d=b.length;c<d;c++)if(f.contains(b[c],this))return!0});var e=this.pushStack("","find",a),g,h,i;for(c=0,d=this.length;c<d;c++){g=e.length,f.find(a,this[c],e);if(c>0)for(h=g;h<e.length;h++)for(i=0;i<g;i++)if(e[i]===e[h]){e.splice(h--,1);break}}return e},has:function(a){var b=f(a);return this.filter(function(){for(var a=0,c=b.length;a<c;a++)if(f.contains(this,b[a]))return!0})},not:function(a){return this.pushStack(W(this,a,!1),"not",a)},filter:function(a){return this.pushStack(W(this,a,!0),"filter",a)},is:function(a){return!!a&&(typeof a=="string"?T.test(a)?f(a,this.context).index(this[0])>=0:f.filter(a,this).length>0:this.filter(a).length>0)},closest:function(a,b){var c=[],d,e,g=this[0];if(f.isArray(a)){var h=1;while(g&&g.ownerDocument&&g!==b){for(d=0;d<a.length;d++)f(g).is(a[d])&&c.push({selector:a[d],elem:g,level:h});g=g.parentNode,h++}return c}var i=T.test(a)||typeof a!="string"?f(a,b||this.context):0;for(d=0,e=this.length;d<e;d++){g=this[d];while(g){if(i?i.index(g)>-1:f.find.matchesSelector(g,a)){c.push(g);break}g=g.parentNode;if(!g||!g.ownerDocument||g===b||g.nodeType===11)break}}c=c.length>1?f.unique(c):c;return this.pushStack(c,"closest",a)},index:function(a){if(!a)return this[0]&&this[0].parentNode?this.prevAll().length:-1;if(typeof a=="string")return f.inArray(this[0],f(a));return f.inArray(a.jquery?a[0]:a,this)},add:function(a,b){var c=typeof a=="string"?f(a,b):f.makeArray(a&&a.nodeType?[a]:a),d=f.merge(this.get(),c);return this.pushStack(V(c[0])||V(d[0])?d:f.unique(d))},andSelf:function(){return this.add(this.prevObject)}}),f.each({parent:function(a){var b=a.parentNode;return b&&b.nodeType!==11?b:null},parents:function(a){return f.dir(a,"parentNode")},parentsUntil:function(a,b,c){return f.dir(a,"parentNode",c)},next:function(a){return f.nth(a,2,"nextSibling")},prev:function(a){return f.nth(a,2,"previousSibling")},nextAll:function(a){return f.dir(a,"nextSibling")},prevAll:function(a){return f.dir(a,"previousSibling")},nextUntil:function(a,b,c){return f.dir(a,"nextSibling",c)},prevUntil:function(a,b,c){return f.dir(a,"previousSibling",c)},siblings:function(a){return f.sibling(a.parentNode.firstChild,a)},children:function(a){return f.sibling(a.firstChild)},contents:function(a){return f.nodeName(a,"iframe")?a.contentDocument||a.contentWindow.document:f.makeArray(a.childNodes)}},function(a,b){f.fn[a]=function(c,d){var e=f.map(this,b,c),g=S.call(arguments);O.test(a)||(d=c),d&&typeof d=="string"&&(e=f.filter(d,e)),e=this.length>1&&!U[a]?f.unique(e):e,(this.length>1||Q.test(d))&&P.test(a)&&(e=e.reverse());return this.pushStack(e,a,g.join(","))}}),f.extend({filter:function(a,b,c){c&&(a=":not("+a+")");return b.length===1?f.find.matchesSelector(b[0],a)?[b[0]]:[]:f.find.matches(a,b)},dir:function(a,c,d){var e=[],g=a[c];while(g&&g.nodeType!==9&&(d===b||g.nodeType!==1||!f(g).is(d)))g.nodeType===1&&e.push(g),g=g[c];return e},nth:function(a,b,c,d){b=b||1;var e=0;for(;a;a=a[c])if(a.nodeType===1&&++e===b)break;return a},sibling:function(a,b){var c=[];for(;a;a=a.nextSibling)a.nodeType===1&&a!==b&&c.push(a);return c}});var Y="abbr article aside audio canvas datalist details figcaption figure footer header hgroup mark meter nav output progress section summary time video",Z=/ jQuery\d+="(?:\d+|null)"/g,$=/^\s+/,_=/<(?!area|br|col|embed|hr|img|input|link|meta|param)(([\w:]+)[^>]*)\/>/ig,ba=/<([\w:]+)/,bb=/<tbody/i,bc=/<|&#?\w+;/,bd=/<(?:script|style)/i,be=/<(?:script|object|embed|option|style)/i,bf=new RegExp("<(?:"+Y.replace(" ","|")+")","i"),bg=/checked\s*(?:[^=]|=\s*.checked.)/i,bh=/\/(java|ecma)script/i,bi=/^\s*<!(?:\[CDATA\[|\-\-)/,bj={option:[1,"<select multiple='multiple'>","</select>"],legend:[1,"<fieldset>","</fieldset>"],thead:[1,"<table>","</table>"],tr:[2,"<table><tbody>","</tbody></table>"],td:[3,"<table><tbody><tr>","</tr></tbody></table>"],col:[2,"<table><tbody></tbody><colgroup>","</colgroup></table>"],area:[1,"<map>","</map>"],_default:[0,"",""]},bk=X(c);bj.optgroup=bj.option,bj.tbody=bj.tfoot=bj.colgroup=bj.caption=bj.thead,bj.th=bj.td,f.support.htmlSerialize||(bj._default=[1,"div<div>","</div>"]),f.fn.extend({text:function(a){if(f.isFunction(a))return this.each(function(b){var c=f(this);c.text(a.call(this,b,c.text()))});if(typeof a!="object"&&a!==b)return this.empty().append((this[0]&&this[0].ownerDocument||c).createTextNode(a));return f.text(this)},wrapAll:function(a){if(f.isFunction(a))return this.each(function(b){f(this).wrapAll(a.call(this,b))});if(this[0]){var b=f(a,this[0].ownerDocument).eq(0).clone(!0);this[0].parentNode&&b.insertBefore(this[0]),b.map(function(){var a=this;while(a.firstChild&&a.firstChild.nodeType===1)a=a.firstChild;return a}).append(this)}return this},wrapInner:function(a){if(f.isFunction(a))return this.each(function(b){f(this).wrapInner(a.call(this,b))});return this.each(function(){var b=f(this),c=b.contents();c.length?c.wrapAll(a):b.append(a)})},wrap:function(a){return this.each(function(){f(this).wrapAll(a)})},unwrap:function(){return this.parent().each(function(){f.nodeName(this,"body")||f(this).replaceWith(this.childNodes)}).end()},append:function(){return this.domManip(arguments,!0,function(a){this.nodeType===1&&this.appendChild(a)})},prepend:function(){return this.domManip(arguments,!0,function(a){this.nodeType===1&&this.insertBefore(a,this.firstChild)})},before:function(){if(this[0]&&this[0].parentNode)return this.domManip(arguments,!1,function(a){this.parentNode.insertBefore(a,this)});if(arguments.length){var a=f(arguments[0]);a.push.apply(a,this.toArray());return this.pushStack(a,"before",arguments)}},after:function(){if(this[0]&&this[0].parentNode)return this.domManip(arguments,!1,function(a){this.parentNode.insertBefore(a,this.nextSibling)});if(arguments.length){var a=this.pushStack(this,"after"
,arguments);a.push.apply(a,f(arguments[0]).toArray());return a}},remove:function(a,b){for(var c=0,d;(d=this[c])!=null;c++)if(!a||f.filter(a,[d]).length)!b&&d.nodeType===1&&(f.cleanData(d.getElementsByTagName("*")),f.cleanData([d])),d.parentNode&&d.parentNode.removeChild(d);return this},empty:function(){for(var a=0,b;(b=this[a])!=null;a++){b.nodeType===1&&f.cleanData(b.getElementsByTagName("*"));while(b.firstChild)b.removeChild(b.firstChild)}return this},clone:function(a,b){a=a==null?!1:a,b=b==null?a:b;return this.map(function(){return f.clone(this,a,b)})},html:function(a){if(a===b)return this[0]&&this[0].nodeType===1?this[0].innerHTML.replace(Z,""):null;if(typeof a=="string"&&!bd.test(a)&&(f.support.leadingWhitespace||!$.test(a))&&!bj[(ba.exec(a)||["",""])[1].toLowerCase()]){a=a.replace(_,"<$1></$2>");try{for(var c=0,d=this.length;c<d;c++)this[c].nodeType===1&&(f.cleanData(this[c].getElementsByTagName("*")),this[c].innerHTML=a)}catch(e){this.empty().append(a)}}else f.isFunction(a)?this.each(function(b){var c=f(this);c.html(a.call(this,b,c.html()))}):this.empty().append(a);return this},replaceWith:function(a){if(this[0]&&this[0].parentNode){if(f.isFunction(a))return this.each(function(b){var c=f(this),d=c.html();c.replaceWith(a.call(this,b,d))});typeof a!="string"&&(a=f(a).detach());return this.each(function(){var b=this.nextSibling,c=this.parentNode;f(this).remove(),b?f(b).before(a):f(c).append(a)})}return this.length?this.pushStack(f(f.isFunction(a)?a():a),"replaceWith",a):this},detach:function(a){return this.remove(a,!0)},domManip:function(a,c,d){var e,g,h,i,j=a[0],k=[];if(!f.support.checkClone&&arguments.length===3&&typeof j=="string"&&bg.test(j))return this.each(function(){f(this).domManip(a,c,d,!0)});if(f.isFunction(j))return this.each(function(e){var g=f(this);a[0]=j.call(this,e,c?g.html():b),g.domManip(a,c,d)});if(this[0]){i=j&&j.parentNode,f.support.parentNode&&i&&i.nodeType===11&&i.childNodes.length===this.length?e={fragment:i}:e=f.buildFragment(a,this,k),h=e.fragment,h.childNodes.length===1?g=h=h.firstChild:g=h.firstChild;if(g){c=c&&f.nodeName(g,"tr");for(var l=0,m=this.length,n=m-1;l<m;l++)d.call(c?bl(this[l],g):this[l],e.cacheable||m>1&&l<n?f.clone(h,!0,!0):h)}k.length&&f.each(k,br)}return this}}),f.buildFragment=function(a,b,d){var e,g,h,i,j=a[0];b&&b[0]&&(i=b[0].ownerDocument||b[0]),i.createDocumentFragment||(i=c),a.length===1&&typeof j=="string"&&j.length<512&&i===c&&j.charAt(0)==="<"&&!be.test(j)&&(f.support.checkClone||!bg.test(j))&&!f.support.unknownElems&&bf.test(j)&&(g=!0,h=f.fragments[j],h&&h!==1&&(e=h)),e||(e=i.createDocumentFragment(),f.clean(a,i,e,d)),g&&(f.fragments[j]=h?e:1);return{fragment:e,cacheable:g}},f.fragments={},f.each({appendTo:"append",prependTo:"prepend",insertBefore:"before",insertAfter:"after",replaceAll:"replaceWith"},function(a,b){f.fn[a]=function(c){var d=[],e=f(c),g=this.length===1&&this[0].parentNode;if(g&&g.nodeType===11&&g.childNodes.length===1&&e.length===1){e[b](this[0]);return this}for(var h=0,i=e.length;h<i;h++){var j=(h>0?this.clone(!0):this).get();f(e[h])[b](j),d=d.concat(j)}return this.pushStack(d,a,e.selector)}}),f.extend({clone:function(a,b,c){var d=a.cloneNode(!0),e,g,h;if((!f.support.noCloneEvent||!f.support.noCloneChecked)&&(a.nodeType===1||a.nodeType===11)&&!f.isXMLDoc(a)){bn(a,d),e=bo(a),g=bo(d);for(h=0;e[h];++h)g[h]&&bn(e[h],g[h])}if(b){bm(a,d);if(c){e=bo(a),g=bo(d);for(h=0;e[h];++h)bm(e[h],g[h])}}e=g=null;return d},clean:function(a,b,d,e){var g;b=b||c,typeof b.createElement=="undefined"&&(b=b.ownerDocument||b[0]&&b[0].ownerDocument||c);var h=[],i;for(var j=0,k;(k=a[j])!=null;j++){typeof k=="number"&&(k+="");if(!k)continue;if(typeof k=="string")if(!bc.test(k))k=b.createTextNode(k);else{k=k.replace(_,"<$1></$2>");var l=(ba.exec(k)||["",""])[1].toLowerCase(),m=bj[l]||bj._default,n=m[0],o=b.createElement("div");b===c?bk.appendChild(o):X(b).appendChild(o),o.innerHTML=m[1]+k+m[2];while(n--)o=o.lastChild;if(!f.support.tbody){var p=bb.test(k),q=l==="table"&&!p?o.firstChild&&o.firstChild.childNodes:m[1]==="<table>"&&!p?o.childNodes:[];for(i=q.length-1;i>=0;--i)f.nodeName(q[i],"tbody")&&!q[i].childNodes.length&&q[i].parentNode.removeChild(q[i])}!f.support.leadingWhitespace&&$.test(k)&&o.insertBefore(b.createTextNode($.exec(k)[0]),o.firstChild),k=o.childNodes}var r;if(!f.support.appendChecked)if(k[0]&&typeof (r=k.length)=="number")for(i=0;i<r;i++)bq(k[i]);else bq(k);k.nodeType?h.push(k):h=f.merge(h,k)}if(d){g=function(a){return!a.type||bh.test(a.type)};for(j=0;h[j];j++)if(e&&f.nodeName(h[j],"script")&&(!h[j].type||h[j].type.toLowerCase()==="text/javascript"))e.push(h[j].parentNode?h[j].parentNode.removeChild(h[j]):h[j]);else{if(h[j].nodeType===1){var s=f.grep(h[j].getElementsByTagName("script"),g);h.splice.apply(h,[j+1,0].concat(s))}d.appendChild(h[j])}}return h},cleanData:function(a){var b,c,d=f.cache,e=f.event.special,g=f.support.deleteExpando;for(var h=0,i;(i=a[h])!=null;h++){if(i.nodeName&&f.noData[i.nodeName.toLowerCase()])continue;c=i[f.expando];if(c){b=d[c];if(b&&b.events){for(var j in b.events)e[j]?f.event.remove(i,j):f.removeEvent(i,j,b.handle);b.handle&&(b.handle.elem=null)}g?delete i[f.expando]:i.removeAttribute&&i.removeAttribute(f.expando),delete d[c]}}}});var bs=/alpha\([^)]*\)/i,bt=/opacity=([^)]*)/,bu=/([A-Z]|^ms)/g,bv=/^-?\d+(?:px)?$/i,bw=/^-?\d/,bx=/^([\-+])=([\-+.\de]+)/,by={position:"absolute",visibility:"hidden",display:"block"},bz=["Left","Right"],bA=["Top","Bottom"],bB,bC,bD;f.fn.css=function(a,c){if(arguments.length===2&&c===b)return this;return f.access(this,a,c,!0,function(a,c,d){return d!==b?f.style(a,c,d):f.css(a,c)})},f.extend({cssHooks:{opacity:{get:function(a,b){if(b){var c=bB(a,"opacity","opacity");return c===""?"1":c}return a.style.opacity}}},cssNumber:{fillOpacity:!0,fontWeight:!0,lineHeight:!0,opacity:!0,orphans:!0,widows:!0,zIndex:!0,zoom:!0},cssProps:{"float":f.support.cssFloat?"cssFloat":"styleFloat"},style:function(a,c,d,e){if(!!a&&a.nodeType!==3&&a.nodeType!==8&&!!a.style){var g,h,i=f.camelCase(c),j=a.style,k=f.cssHooks[i];c=f.cssProps[i]||i;if(d===b){if(k&&"get"in k&&(g=k.get(a,!1,e))!==b)return g;return j[c]}h=typeof d,h==="string"&&(g=bx.exec(d))&&(d=+(g[1]+1)*+g[2]+parseFloat(f.css(a,c)),h="number");if(d==null||h==="number"&&isNaN(d))return;h==="number"&&!f.cssNumber[i]&&(d+="px");if(!k||!("set"in k)||(d=k.set(a,d))!==b)try{j[c]=d}catch(l){}}},css:function(a,c,d){var e,g;c=f.camelCase(c),g=f.cssHooks[c],c=f.cssProps[c]||c,c==="cssFloat"&&(c="float");if(g&&"get"in g&&(e=g.get(a,!0,d))!==b)return e;if(bB)return bB(a,c)},swap:function(a,b,c){var d={};for(var e in b)d[e]=a.style[e],a.style[e]=b[e];c.call(a);for(e in b)a.style[e]=d[e]}}),f.curCSS=f.css,f.each(["height","width"],function(a,b){f.cssHooks[b]={get:function(a,c,d){var e;if(c){if(a.offsetWidth!==0)return bE(a,b,d);f.swap(a,by,function(){e=bE(a,b,d)});return e}},set:function(a,b){if(!bv.test(b))return b;b=parseFloat(b);if(b>=0)return b+"px"}}}),f.support.opacity||(f.cssHooks.opacity={get:function(a,b){return bt.test((b&&a.currentStyle?a.currentStyle.filter:a.style.filter)||"")?parseFloat(RegExp.$1)/100+"":b?"1":""},set:function(a,b){var c=a.style,d=a.currentStyle,e=f.isNumeric(b)?"alpha(opacity="+b*100+")":"",g=d&&d.filter||c.filter||"";c.zoom=1;if(b>=1&&f.trim(g.replace(bs,""))===""){c.removeAttribute("filter");if(d&&!d.filter)return}c.filter=bs.test(g)?g.replace(bs,e):g+" "+e}}),f(function(){f.support.reliableMarginRight||(f.cssHooks.marginRight={get:function(a,b){var c;f.swap(a,{display:"inline-block"},function(){b?c=bB(a,"margin-right","marginRight"):c=a.style.marginRight});return c}})}),c.defaultView&&c.defaultView.getComputedStyle&&(bC=function(a,c){var d,e,g;c=c.replace(bu,"-$1").toLowerCase();if(!(e=a.ownerDocument.defaultView))return b;if(g=e.getComputedStyle(a,null))d=g.getPropertyValue(c),d===""&&!f.contains(a.ownerDocument.documentElement,a)&&(d=f.style(a,c));return d}),c.documentElement.currentStyle&&(bD=function(a,b){var c,d,e,f=a.currentStyle&&a.currentStyle[b],g=a.style;f===null&&g&&(e=g[b])&&(f=e),!bv.test(f)&&bw.test(f)&&(c=g.left,d=a.runtimeStyle&&a.runtimeStyle.left,d&&(a.runtimeStyle.left=a.currentStyle.left),g.left=b==="fontSize"?"1em":f||0,f=g.pixelLeft+"px",g.left=c,d&&(a.runtimeStyle.left=d));return f===""?"auto":f}),bB=bC||bD,f.expr&&f.expr.filters&&(f.expr.filters.hidden=function(a){var b=a.offsetWidth,c=a.offsetHeight;return b===0&&c===0||!f.support.reliableHiddenOffsets&&(a.style&&a.style.display||f.css(a,"display"))==="none"},f.expr.filters.visible=function(a){return!f.expr.filters.hidden(a)});var bF=/%20/g,bG=/\[\]$/,bH=/\r?\n/g,bI=/#.*$/,bJ=/^(.*?):[ \t]*([^\r\n]*)\r?$/mg,bK=/^(?:color|date|datetime|datetime-local|email|hidden|month|number|password|range|search|tel|text|time|url|week)$/i,bL=/^(?:about|app|app\-storage|.+\-extension|file|res|widget):$/,bM=/^(?:GET|HEAD)$/,bN=/^\/\//,bO=/\?/,bP=/<script\b[^<]*(?:(?!<\/script>)<[^<]*)*<\/script>/gi,bQ=/^(?:select|textarea)/i,bR=/\s+/,bS=/([?&])_=[^&]*/,bT=/^([\w\+\.\-]+:)(?:\/\/([^\/?#:]*)(?::(\d+))?)?/,bU=f.fn.load,bV={},bW={},bX,bY,bZ=["*/"]+["*"];try{bX=e.href}catch(b$){bX=c.createElement("a"),bX.href="",bX=bX.href}bY=bT.exec(bX.toLowerCase())||[],f.fn.extend({load:function(a,c,d){if(typeof a!="string"&&bU)return bU.apply(this,arguments);if(!this.length)return this;var e=a.indexOf(" ");if(e>=0){var g=a.slice(e,a.length);a=a.slice(0,e)}var h="GET";c&&(f.isFunction(c)?(d=c,c=b):typeof c=="object"&&(c=f.param(c,f.ajaxSettings.traditional),h="POST"));var i=this;f.ajax({url:a,type:h,dataType:"html",data:c,complete:function(a,b,c){c=a.responseText,a.isResolved()&&(a.done(function(a){c=a}),i.html(g?f("<div>").append(c.replace(bP,"")).find(g):c)),d&&i.each(d,[c,b,a])}});return this},serialize:function(){return f.param(this.serializeArray())},serializeArray:function(){return this.map(function(){return this.elements?f.makeArray(this.elements):this}).filter(function(){return this.name&&!this.disabled&&(this.checked||bQ.test(this.nodeName)||bK.test(this.type))}).map(function(a,b){var c=f(this).val();return c==null?null:f.isArray(c)?f.map(c,function(a,c){return{name:b.name,value:a.replace(bH,"\r\n")}}):{name:b.name,value:c.replace(bH,"\r\n")}}).get()}}),f.each("ajaxStart ajaxStop ajaxComplete ajaxError ajaxSuccess ajaxSend".split(" "),function(a,b){f.fn[b]=function(a){return this.bind(b,a)}}),f.each(["get","post"],function(a,c){f[c]=function(a,d,e,g){f.isFunction(d)&&(g=g||e,e=d,d=b);return f.ajax({type:c,url:a,data:d,success:e,dataType:g})}}),f.extend({getScript:function(a,c){return f.get(a,b,c,"script")},getJSON:function(a,b,c){return f.get(a,b,c,"json")},ajaxSetup:function(a,b){b?cb(a,f.ajaxSettings):(b=a,a=f.ajaxSettings),cb(a,b);return a},ajaxSettings:{url:bX,isLocal:bL.test(bY[1]),global:!0,type:"GET",contentType:"application/x-www-form-urlencoded",processData:!0,async:!0,accepts:{xml:"application/xml, text/xml",html:"text/html",text:"text/plain",json:"application/json, text/javascript","*":bZ},contents:{xml:/xml/,html:/html/,json:/json/},responseFields:{xml:"responseXML",text:"responseText"},converters:{"* text":a.String,"text html":!0,"text json":f.parseJSON,"text xml":f.parseXML},flatOptions:{context:!0,url:!0}},ajaxPrefilter:b_(bV),ajaxTransport:b_(bW),ajax:function(a,c){function w(a,c,l,m){if(s!==2){s=2,q&&clearTimeout(q),p=b,n=m||"",v.readyState=a>0?4:0;var o,r,u,w=c,x=l?cd(d,v,l):b,y,z;if(a>=200&&a<300||a===304){if(d.ifModified){if(y=v.getResponseHeader("Last-Modified"))f.lastModified[k]=y;if(z=v.getResponseHeader("Etag"))f.etag[k]=z}if(a===304)w="notmodified",o=!0;else try{r=ce(d,x),w="success",o=!0}catch(A){w="parsererror",u=A}}else{u=w;if(!w||a)w="error",a<0&&(a=0)}v.status=a,v.statusText=""+(c||w),o?h.resolveWith(e,[r,w,v]):h.rejectWith(e,[v,w,u]),v.statusCode(j),j=b,t&&g.trigger("ajax"+(o?"Success":"Error"),[v,d,o?r:u]),i.fireWith(e,[v,w]),t&&(g.trigger("ajaxComplete",[v,d]),--f.active||f.event.trigger("ajaxStop"))}}typeof a=="object"&&(c=a,a=b),c=c||{};var d=f.ajaxSetup({},c),e=d.context||d,g=e!==d&&(e.nodeType||e instanceof f)?f(e):f.event,h=f.Deferred(),i=f.Callbacks("once memory"),j=d.statusCode||{},k,l={},m={},n,o,p,q,r,s=0,t,u,v={readyState:0,setRequestHeader:function(a,b){if(!s){var c=a.toLowerCase();a=m[c]=m[c]||a,l[a]=b}return this},getAllResponseHeaders:function(){return s===2?n:null},getResponseHeader:function(a){var c;if(s===2){if(!o){o={};while(c=bJ.exec(n))o[c[1].toLowerCase()]=c[2]}c=o[a.toLowerCase()]}return c===b?null:c},overrideMimeType:function(a){s||(d.mimeType=a);return this},abort:function(a){a=a||"abort",p&&p.abort(a),w(0,a);return this}};h.promise(v),v.success=v.done,v.error=v.fail,v.complete=i.add,v.statusCode=function(a){if(a){var b;if(s<2)for(b in a)j[b]=[j[b],a[b]];else b=a[v.status],v.then(b,b)}return this},d.url=((a||d.url)+"").replace(bI,"").replace(bN,bY[1]+"//"),d.dataTypes=f.trim(d.dataType||"*").toLowerCase().split(bR),d.crossDomain==null&&(r=bT.exec(d.url.toLowerCase()),d.crossDomain=!(!r||r[1]==bY[1]&&r[2]==bY[2]&&(r[3]||(r[1]==="http:"?80:443))==(bY[3]||(bY[1]==="http:"?80:443)))),d.data&&d.processData&&typeof d.data!="string"&&(d.data=f.param(d.data,d.traditional)),ca(bV,d,c,v);if(s===2)return!1;t=d.global,d.type=d.type.toUpperCase(),d.hasContent=!bM.test(d.type),t&&f.active++===0&&f.event.trigger("ajaxStart");if(!d.hasContent){d.data&&(d.url+=(bO.test(d.url)?"&":"?")+d.data,delete d.data),k=d.url;if(d.cache===!1){var x=f.now(),y=d.url.replace(bS,"$1_="+x);d.url=y+(y===d.url?(bO.test(d.url)?"&":"?")+"_="+x:"")}}(d.data&&d.hasContent&&d.contentType!==!1||c.contentType)&&v.setRequestHeader("Content-Type",d.contentType),d.ifModified&&(k=k||d.url,f.lastModified[k]&&v.setRequestHeader("If-Modified-Since",f.lastModified[k]),f.etag[k]&&v.setRequestHeader("If-None-Match",f.etag[k])),v.setRequestHeader("Accept",d.dataTypes[0]&&d.accepts[d.dataTypes[0]]?d.accepts[d.dataTypes[0]]+(d.dataTypes[0]!=="*"?", "+bZ+"; q=0.01":""):d.accepts["*"]);for(u in d.headers)v.setRequestHeader(u,d.headers[u]);if(d.beforeSend&&(d.beforeSend.call(e,v,d)===!1||s===2)){v.abort();return!1}for(u in{success:1,error:1,complete:1})v[u](d[u]);p=ca(bW,d,c,v);if(!p)w(-1,"No Transport");else{v.readyState=1,t&&g.trigger("ajaxSend",[v,d]),d.async&&d.timeout>0&&(q=setTimeout(function(){v.abort("timeout")},d.timeout));try{s=1,p.send(l,w)}catch(z){s<2?w(-1,z):f.error(z)}}return v},param:function(a,c){var d=[],e=function(a,b){b=f.isFunction(b)?b():b,d[d.length]=encodeURIComponent(a)+"="+encodeURIComponent(b)};c===b&&(c=f.ajaxSettings.traditional);if(f.isArray(a)||a.jquery&&!f.isPlainObject(a))f.each(a,function(){e(this.name,this.value)});else for(var g in a)cc(g,a[g],c,e);return d.join("&").replace(bF,"+")}}),f.extend({active:0,lastModified:{},etag:{}});var cf=f.now(),cg=/(\=)\?(&|$)|\?\?/i;f.ajaxSetup({jsonp:"callback",jsonpCallback:function(){return f.expando+"_"+cf++}}),f.ajaxPrefilter("json jsonp",function(b,c,d){var e=b.contentType==="application/x-www-form-urlencoded"&&typeof b.data=="string";if(b.dataTypes[0]==="jsonp"||b.jsonp!==!1&&(cg.test(b.url)||e&&cg.test(b.data))){var g,h=b.jsonpCallback=f.isFunction(b.jsonpCallback)?b.jsonpCallback():b.jsonpCallback,i=a[h],j=b.url,k=b.data,l="$1"+h+"$2";b.jsonp!==!1&&(j=j.replace(cg,l),b.url===j&&(e&&(k=k.replace(cg,l)),b.data===k&&(j+=(/\?/.test(j)?"&":"?")+b.jsonp+"="+h))),b.url=j,b.data=k,a[h]=function(a){g=[a]},d.always(function(){a[h]=i,g&&f.isFunction(i)&&a[h](g[0])}),b.converters["script json"]=function(){g||f.error(h+" was not called");return g[0]},b.dataTypes[0]="json";return"script"}}),f.ajaxSetup({accepts:{script:"text/javascript, application/javascript, application/ecmascript, application/x-ecmascript"},contents:{script:/javascript|ecmascript/},converters:{"text script":function(a){f.globalEval(a);return a}}}),f.ajaxPrefilter("script",function(a){a.cache===b&&(a.cache=!1),a.crossDomain&&(a.type="GET",a.global=!1)}),f.ajaxTransport("script",function(a){if(a.crossDomain){var d,e=c.head||c.getElementsByTagName("head")[0]||c.documentElement;return{send:function(f,g){d=c.createElement("script"),d.async="async",a.scriptCharset&&(d.charset=a.scriptCharset),d.src=a.url,d.onload=d.onreadystatechange=function(a,c){if(c||!d.readyState||/loaded|complete/.test(d.readyState))d.onload=d.onreadystatechange=null,e&&d.parentNode&&e.removeChild(d),d=b,c||g(200,"success")},e.insertBefore(d,e.firstChild)},abort:function(){d&&d.onload(0,1)}}}});var ch=a.ActiveXObject?function(){for(var a in cj)cj[a](0,1)}:!1,ci=0,cj;f.ajaxSettings.xhr=a.ActiveXObject?function(){return!this.isLocal&&ck()||cl()}:ck,function(a){f.extend(f.support,{ajax:!!a,cors:!!a&&"withCredentials"in a})}(f.ajaxSettings.xhr()),f.support.ajax&&f.ajaxTransport(function(c){if(!c.crossDomain||f.support.cors){var d;return{send:function(e,g){var h=c.xhr(),i,j;c.username?h.open(c.type,c.url,c.async,c.username,c.password):h.open(c.type,c.url,c.async);if(c.xhrFields)for(j in c.xhrFields)h[j]=c.xhrFields[j];c.mimeType&&h.overrideMimeType&&h.overrideMimeType(c.mimeType),!c.crossDomain&&!e["X-Requested-With"]&&(e["X-Requested-With"]="XMLHttpRequest");try{for(j in e)h.setRequestHeader(j,e[j])}catch(k){}h.send(c.hasContent&&c.data||null),d=function(a,e){var j,k,l,m,n;try{if(d&&(e||h.readyState===4)){d=b,i&&(h.onreadystatechange=f.noop,ch&&delete cj[i]);if(e)h.readyState!==4&&h.abort();else{j=h.status,l=h.getAllResponseHeaders(),m={},n=h.responseXML,n&&n.documentElement&&(m.xml=n),m.text=h.responseText;try{k=h.statusText}catch(o){k=""}!j&&c.isLocal&&!c.crossDomain?j=m.text?200:404:j===1223&&(j=204)}}}catch(p){e||g(-1,p)}m&&g(j,k,m,l)},!c.async||h.readyState===4?d():(i=++ci,ch&&(cj||(cj={},f(a).unload(ch)),cj[i]=d),h.onreadystatechange=d)},abort:function(){d&&d(0,1)}}}});var cm={},cn,co,cp=/^(?:toggle|show|hide)$/,cq=/^([+\-]=)?([\d+.\-]+)([a-z%]*)$/i,cr,cs=[["height","marginTop","marginBottom","paddingTop","paddingBottom"],["width","marginLeft","marginRight","paddingLeft","paddingRight"],["opacity"]],ct;f.fn.extend({show:function(a,b,c){var d,e;if(a||a===0)return this.animate(cw("show",3),a,b,c);for(var g=0,h=this.length;g<h;g++)d=this[g],d.style&&(e=d.style.display,!f._data(d,"olddisplay")&&e==="none"&&(e=d.style.display=""),e===""&&f.css(d,"display")==="none"&&f._data(d,"olddisplay",cx(d.nodeName)));for(g=0;g<h;g++){d=this[g];if(d.style){e=d.style.display;if(e===""||e==="none")d.style.display=f._data(d,"olddisplay")||""}}return this},hide:function(a,b,c){if(a||a===0)return this.animate(cw("hide",3),a,b,c);var d,e,g=0,h=this.length;for(;g<h;g++)d=this[g],d.style&&(e=f.css(d,"display"),e!=="none"&&!f._data(d,"olddisplay")&&f._data(d,"olddisplay",e));for(g=0;g<h;g++)this[g].style&&(this[g].style.display="none");return this},_toggle:f.fn.toggle,toggle:function(a,b,c){var d=typeof a=="boolean";f.isFunction(a)&&f.isFunction(b)?this._toggle.apply(this,arguments):a==null||d?this.each(function(){var b=d?a:f(this).is(":hidden");f(this)[b?"show":"hide"]()}):this.animate(cw("toggle",3),a,b,c);return this},fadeTo:function(a,b,c,d){return this.filter(":hidden").css("opacity",0).show().end().animate({opacity:b},a,c,d)},animate:function(a,b,c,d){function g(){e.queue===!1&&f._mark(this);var b=f.extend({},e),c=this.nodeType===1,d=c&&f(this).is(":hidden"),g,h,i,j,k,l,m,n,o;b.animatedProperties={};for(i in a){g=f.camelCase(i),i!==g&&(a[g]=a[i],delete a[i]),h=a[g],f.isArray(h)?(b.animatedProperties[g]=h[1],h=a[g]=h[0]):b.animatedProperties[g]=b.specialEasing&&b.specialEasing[g]||b.easing||"swing";if(h==="hide"&&d||h==="show"&&!d)return b.complete.call(this);c&&(g==="height"||g==="width")&&(b.overflow=[this.style.overflow,this.style.overflowX,this.style.overflowY],f.css(this,"display")==="inline"&&f.css(this,"float")==="none"&&(!f.support.inlineBlockNeedsLayout||cx(this.nodeName)==="inline"?this.style.display="inline-block":this.style.zoom=1))}b.overflow!=null&&(this.style.overflow="hidden");for(i in a)j=new f.fx(this,b,i),h=a[i],cp.test(h)?(o=f._data(this,"toggle"+i)||(h==="toggle"?d?"show":"hide":0),o?(f._data(this,"toggle"+i,o==="show"?"hide":"show"),j[o]()):j[h]()):(k=cq.exec(h),l=j.cur(),k?(m=parseFloat(k[2]),n=k[3]||(f.cssNumber[i]?"":"px"),n!=="px"&&(f.style(this,i,(m||1)+n),l=(m||1)/j.cur()*l,f.style(this,i,l+n)),k[1]&&(m=(k[1]==="-="?-1:1)*m+l),j.custom(l,m,n)):j.custom(l,h,""));return!0}var e=f.speed(b,c,d);if(f.isEmptyObject(a))return this.each(e.complete,[!1]);a=f.extend({},a);return e.queue===!1?this.each(g):this.queue(e.queue,g)},stop:function(a,c,d){typeof a!="string"&&(d=c,c=a,a=b),c&&a!==!1&&this.queue(a||"fx",[]);return this.each(function(){function h(a,b,c){var e=b[c];f.removeData(a,c,!0),e.stop(d)}var b,c=!1,e=f.timers,g=f._data(this);d||f._unmark(!0,this);if(a==null)for(b in g)g[b].stop&&b.indexOf(".run")===b.length-4&&h(this,g,b);else g[b=a+".run"]&&g[b].stop&&h(this,g,b);for(b=e.length;b--;)e[b].elem===this&&(a==null||e[b].queue===a)&&(d?e[b](!0):e[b].saveState(),c=!0,e.splice(b,1));(!d||!c)&&f.dequeue(this,a)})}}),f.each({slideDown:cw("show",1),slideUp:cw("hide",1),slideToggle:cw("toggle",1),fadeIn:{opacity:"show"},fadeOut:{opacity:"hide"},fadeToggle:{opacity:"toggle"}},function(a,b){f.fn[a]=function(a,c,d){return this.animate(b,a,c,d)}}),f.extend({speed:function(a,b,c){var d=a&&typeof a=="object"?f.extend({},a):{complete:c||!c&&b||f.isFunction(a)&&a,duration:a,easing:c&&b||b&&!f.isFunction(b)&&b};d.duration=f.fx.off?0:typeof d.duration=="number"?d.duration:d.duration in f.fx.speeds?f.fx.speeds[d.duration]:f.fx.speeds._default;if(d.queue==null||d.queue===!0)d.queue="fx";d.old=d.complete,d.complete=function(a){f.isFunction(d.old)&&d.old.call(this),d.queue?f.dequeue(this,d.queue):a!==!1&&f._unmark(this)};return d},easing:{linear:function(a,b,c,d){return c+d*a},swing:function(a,b,c,d){return(-Math.cos(a*Math.PI)/2+.5)*d+c}},timers:[],fx:function(a,b,c){this.options=b,this.elem=a,this.prop=c,b.orig=b.orig||{}}}),f.fx.prototype={update:function(){this.options.step&&this.options.step.call(this.elem,this.now,this),(f.fx.step[this.prop]||f.fx.step._default)(this)},cur:function(){if(this.elem[this.prop]!=null&&(!this.elem.style||this.elem.style[this.prop]==null))return this.elem[this.prop];var a,b=f.css(this.elem,this.prop);return isNaN(a=parseFloat(b))?!b||b==="auto"?0:b:a},custom:function(a,c,d){function h(a){return e.step(a)}var e=this,g=f.fx;this.startTime=ct||cu(),this.end=c,this.now=this.start=a,this.pos=this.state=0,this.unit=d||this.unit||(f.cssNumber[this.prop]?"":"px"),h.queue=this.options.queue,h.elem=this.elem,h.saveState=function(){e.options.hide&&f._data(e.elem,"fxshow"+e.prop)===b&&f._data(e.elem,"fxshow"+e.prop,e.start)},h()&&f.timers.push(h)&&!cr&&(cr=setInterval(g.tick,g.interval))},show:function(){var a=f._data(this.elem,"fxshow"+this.prop);this.options.orig[this.prop]=a||f.style(this.elem,this.prop),this.options.show=!0,a!==b?this.custom(this.cur(),a):this.custom(this.prop==="width"||this.prop==="height"?1:0,this.cur()),f(this.elem).show()},hide:function(){this.options.orig[this.prop]=f._data(this.elem,"fxshow"+this.prop)||f.style(this.elem,this.prop),this.options.hide=!0,this.custom(this.cur(),0)},step:function(a){var b,c,d,e=ct||cu(),g=!0,h=this.elem,i=this.options;if(a||e>=i.duration+this.startTime){this.now=this.end,this.pos=this.state=1,this.update(),i.animatedProperties[this.prop]=!0;for(b in i.animatedProperties)i.animatedProperties[b]!==!0&&(g=!1);if(g){i.overflow!=null&&!f.support.shrinkWrapBlocks&&f.each(["","X","Y"],function(a,b){h.style["overflow"+b]=i.overflow[a]}),i.hide&&f(h).hide();if(i.hide||i.show)for(b in i.animatedProperties)f.style(h,b,i.orig[b]),f.removeData(h,"fxshow"+b,!0),f.removeData(h,"toggle"+b,!0);d=i.complete,d&&(i.complete=!1,d.call(h))}return!1}i.duration==Infinity?this.now=e:(c=e-this.startTime,this.state=c/i.duration,this.pos=f.easing[i.animatedProperties[this.prop]](this.state,c,0,1,i.duration),this.now=this.start+(this.end-this.start)*this.pos),this.update();return!0}},f.extend(f.fx,{tick:function(){var a,b=f.timers,c=0;for(;c<b.length;c++)a=b[c],!a()&&b[c]===a&&b.splice(c--,1);b.length||f.fx.stop()},interval:13,stop:function(){clearInterval(cr),cr=null},speeds:{slow:600,fast:200,_default:400},step:{opacity:function(a){f.style(a.elem,"opacity",a.now)},_default:function(a){a.elem.style&&a.elem.style[a.prop]!=null?a.elem.style[a.prop]=a.now+a.unit:a.elem[a.prop]=a.now}}}),f.each(["width","height"],function(a,b){f.fx.step[b]=function(a){f.style(a.elem,b,Math.max(0,a.now))}}),f.expr&&f.expr.filters&&(f.expr.filters.animated=function(a){return f.grep(f.timers,function(b){return a===b.elem}).length});var cy=/^t(?:able|d|h)$/i,cz=/^(?:body|html)$/i;"getBoundingClientRect"in c.documentElement?f.fn.offset=function(a){var b=this[0],c;if(a)return this.each(function(b){f.offset.setOffset(this,a,b)});if(!b||!b.ownerDocument)return null;if(b===b.ownerDocument.body)return f.offset.bodyOffset(b);try{c=b.getBoundingClientRect()}catch(d){}var e=b.ownerDocument,g=e.documentElement;if(!c||!f.contains(g,b))return c?{top:c.top,left:c.left}:{top:0,left:0};var h=e.body,i=cA(e),j=g.clientTop||h.clientTop||0,k=g.clientLeft||h.clientLeft||0,l=i.pageYOffset||f.support.boxModel&&g.scrollTop||h.scrollTop,m=i.pageXOffset||f.support.boxModel&&g.scrollLeft||h.scrollLeft,n=c.top+l-j,o=c.left+m-k;return{top:n,left:o}}:f.fn.offset=function(a){var b=this[0];if(a)return this.each(function(b){f.offset.setOffset(this,a,b)});if(!b||!b.ownerDocument)return null;if(b===b.ownerDocument.body)return f.offset.bodyOffset(b);var c,d=b.offsetParent,e=b,g=b.ownerDocument,h=g.documentElement,i=g.body,j=g.defaultView,k=j?j.getComputedStyle(b,null):b.currentStyle,l=b.offsetTop,m=b.offsetLeft;while((b=b.parentNode)&&b!==i&&b!==h){if(f.support.fixedPosition&&k.position==="fixed")break;c=j?j.getComputedStyle(b,null):b.currentStyle,l-=b.scrollTop,m-=b.scrollLeft,b===d&&(l+=b.offsetTop,m+=b.offsetLeft,f.support.doesNotAddBorder&&(!f.support.doesAddBorderForTableAndCells||!cy.test(b.nodeName))&&(l+=parseFloat(c.borderTopWidth)||0,m+=parseFloat(c.borderLeftWidth)||0),e=d,d=b.offsetParent),f.support.subtractsBorderForOverflowNotVisible&&c.overflow!=="visible"&&(l+=parseFloat(c.borderTopWidth)||0,m+=parseFloat(c.borderLeftWidth)||0),k=c}if(k.position==="relative"||k.position==="static")l+=i.offsetTop,m+=i.offsetLeft;f.support.fixedPosition&&k.position==="fixed"&&(l+=Math.max(h.scrollTop,i.scrollTop),m+=Math.max(h.scrollLeft,i.scrollLeft));return{top:l,left:m}},f.offset={bodyOffset:function(a){var b=a.offsetTop,c=a.offsetLeft;f.support.doesNotIncludeMarginInBodyOffset&&(b+=parseFloat(f.css(a,"marginTop"))||0,c+=parseFloat(f.css(a,"marginLeft"))||0);return{top:b,left:c}},setOffset:function(a,b,c){var d=f.css(a,"position");d==="static"&&(a.style.position="relative");var e=f(a),g=e.offset(),h=f.css(a,"top"),i=f.css(a,"left"),j=(d==="absolute"||d==="fixed")&&f.inArray("auto",[h,i])>-1,k={},l={},m,n;j?(l=e.position(),m=l.top,n=l.left):(m=parseFloat(h)||0,n=parseFloat(i)||0),f.isFunction(b)&&(b=b.call(a,c,g)),b.top!=null&&(k.top=b.top-g.top+m),b.left!=null&&(k.left=b.left-g.left+n),"using"in b?b.using.call(a,k):e.css(k)}},f.fn.extend({position:function(){if(!this[0])return null;var a=this[0],b=this.offsetParent(),c=this.offset(),d=cz.test(b[0].nodeName)?{top:0,left:0}:b.offset();c.top-=parseFloat(f.css(a,"marginTop"))||0,c.left-=parseFloat(f.css(a,"marginLeft"))||0,d.top+=parseFloat(f.css(b[0],"borderTopWidth"))||0,d.left+=parseFloat(f.css(b[0],"borderLeftWidth"))||0;return{top:c.top-d.top,left:c.left-d.left}},offsetParent:function(){return this.map(function(){var a=this.offsetParent||c.body;while(a&&!cz.test(a.nodeName)&&f.css(a,"position")==="static")a=a.offsetParent;return a})}}),f.each(["Left","Top"],function(a,c){var d="scroll"+c;f.fn[d]=function(c){var e,g;if(c===b){e=this[0];if(!e)return null;g=cA(e);return g?"pageXOffset"in g?g[a?"pageYOffset":"pageXOffset"]:f.support.boxModel&&g.document.documentElement[d]||g.document.body[d]:e[d]}return this.each(function(){g=cA(this),g?g.scrollTo(a?f(g).scrollLeft():c,a?c:f(g).scrollTop()):this[d]=c})}}),f.each(["Height","Width"],function(a,c){var d=c.toLowerCase();f.fn["inner"+c]=function(){var a=this[0];return a?a.style?parseFloat(f.css(a,d,"padding")):this[d]():null},f.fn["outer"+c]=function(a){var b=this[0];return b?b.style?parseFloat(f.css(b,d,a?"margin":"border")):this[d]():null},f.fn[d]=function(a){var e=this[0];if(!e)return a==null?null:this;if(f.isFunction(a))return this.each(function(b){var c=f(this);c[d](a.call(this,b,c[d]()))});if(f.isWindow(e)){var g=e.document.documentElement["client"+c],h=e.document.body;return e.document.compatMode==="CSS1Compat"&&g||h&&h["client"+c]||g}if(e.nodeType===9)return Math.max(e.documentElement["client"+c],e.body["scroll"+c],e.documentElement["scroll"+c],e.body["offset"+c],e.documentElement["offset"+c]);if(a===b){var i=f.css(e,d),j=parseFloat(i);return f.isNumeric(j)?j:i}return this.css(d,typeof a=="string"?a:a+"px")}}),a.jQuery=a.$=f})(window);

// Three.js - http://github.com/mrdoob/three.js
'use strict';var THREE=THREE||{REVISION:"49"};self.Int32Array||(self.Int32Array=Array,self.Float32Array=Array);
(function(){for(var a=0,b=["ms","moz","webkit","o"],c=0;c<b.length&&!window.requestAnimationFrame;++c){window.requestAnimationFrame=window[b[c]+"RequestAnimationFrame"];window.cancelAnimationFrame=window[b[c]+"CancelAnimationFrame"]||window[b[c]+"CancelRequestAnimationFrame"]}if(!window.requestAnimationFrame)window.requestAnimationFrame=function(b){var c=Date.now(),f=Math.max(0,16-(c-a)),h=window.setTimeout(function(){b(c+f)},f);a=c+f;return h};if(!window.cancelAnimationFrame)window.cancelAnimationFrame=
function(a){clearTimeout(a)}})();THREE.Clock=function(a){this.autoStart=a!==void 0?a:true;this.elapsedTime=this.oldTime=this.startTime=0;this.running=false};THREE.Clock.prototype.start=function(){this.oldTime=this.startTime=Date.now();this.running=true};THREE.Clock.prototype.stop=function(){this.getElapsedTime();this.running=false};THREE.Clock.prototype.getElapsedTime=function(){return this.elapsedTime=this.elapsedTime+this.getDelta()};
THREE.Clock.prototype.getDelta=function(){var a=0;this.autoStart&&!this.running&&this.start();if(this.running){var b=Date.now(),a=0.001*(b-this.oldTime);this.oldTime=b;this.elapsedTime=this.elapsedTime+a}return a};THREE.Color=function(a){a!==void 0&&this.setHex(a);return this};
THREE.Color.prototype={constructor:THREE.Color,r:1,g:1,b:1,copy:function(a){this.r=a.r;this.g=a.g;this.b=a.b;return this},copyGammaToLinear:function(a){this.r=a.r*a.r;this.g=a.g*a.g;this.b=a.b*a.b;return this},copyLinearToGamma:function(a){this.r=Math.sqrt(a.r);this.g=Math.sqrt(a.g);this.b=Math.sqrt(a.b);return this},convertGammaToLinear:function(){var a=this.r,b=this.g,c=this.b;this.r=a*a;this.g=b*b;this.b=c*c;return this},convertLinearToGamma:function(){this.r=Math.sqrt(this.r);this.g=Math.sqrt(this.g);
this.b=Math.sqrt(this.b);return this},setRGB:function(a,b,c){this.r=a;this.g=b;this.b=c;return this},setHSV:function(a,b,c){var d,e,f;if(c===0)this.r=this.g=this.b=0;else{d=Math.floor(a*6);e=a*6-d;a=c*(1-b);f=c*(1-b*e);b=c*(1-b*(1-e));switch(d){case 1:this.r=f;this.g=c;this.b=a;break;case 2:this.r=a;this.g=c;this.b=b;break;case 3:this.r=a;this.g=f;this.b=c;break;case 4:this.r=b;this.g=a;this.b=c;break;case 5:this.r=c;this.g=a;this.b=f;break;case 6:case 0:this.r=c;this.g=b;this.b=a}}return this},setHex:function(a){a=
Math.floor(a);this.r=(a>>16&255)/255;this.g=(a>>8&255)/255;this.b=(a&255)/255;return this},lerpSelf:function(a,b){this.r=this.r+(a.r-this.r)*b;this.g=this.g+(a.g-this.g)*b;this.b=this.b+(a.b-this.b)*b;return this},getHex:function(){return Math.floor(this.r*255)<<16^Math.floor(this.g*255)<<8^Math.floor(this.b*255)},getContextStyle:function(){return"rgb("+Math.floor(this.r*255)+","+Math.floor(this.g*255)+","+Math.floor(this.b*255)+")"},clone:function(){return(new THREE.Color).setRGB(this.r,this.g,this.b)}};
THREE.Vector2=function(a,b){this.x=a||0;this.y=b||0};
THREE.Vector2.prototype={constructor:THREE.Vector2,set:function(a,b){this.x=a;this.y=b;return this},copy:function(a){this.x=a.x;this.y=a.y;return this},add:function(a,b){this.x=a.x+b.x;this.y=a.y+b.y;return this},addSelf:function(a){this.x=this.x+a.x;this.y=this.y+a.y;return this},sub:function(a,b){this.x=a.x-b.x;this.y=a.y-b.y;return this},subSelf:function(a){this.x=this.x-a.x;this.y=this.y-a.y;return this},multiplyScalar:function(a){this.x=this.x*a;this.y=this.y*a;return this},divideScalar:function(a){if(a){this.x=
this.x/a;this.y=this.y/a}else this.set(0,0);return this},negate:function(){return this.multiplyScalar(-1)},dot:function(a){return this.x*a.x+this.y*a.y},lengthSq:function(){return this.x*this.x+this.y*this.y},length:function(){return Math.sqrt(this.lengthSq())},normalize:function(){return this.divideScalar(this.length())},distanceTo:function(a){return Math.sqrt(this.distanceToSquared(a))},distanceToSquared:function(a){var b=this.x-a.x,a=this.y-a.y;return b*b+a*a},setLength:function(a){return this.normalize().multiplyScalar(a)},
lerpSelf:function(a,b){this.x=this.x+(a.x-this.x)*b;this.y=this.y+(a.y-this.y)*b;return this},equals:function(a){return a.x===this.x&&a.y===this.y},isZero:function(){return this.lengthSq()<1.0E-4},clone:function(){return new THREE.Vector2(this.x,this.y)}};THREE.Vector3=function(a,b,c){this.x=a||0;this.y=b||0;this.z=c||0};
THREE.Vector3.prototype={constructor:THREE.Vector3,set:function(a,b,c){this.x=a;this.y=b;this.z=c;return this},setX:function(a){this.x=a;return this},setY:function(a){this.y=a;return this},setZ:function(a){this.z=a;return this},copy:function(a){this.x=a.x;this.y=a.y;this.z=a.z;return this},add:function(a,b){this.x=a.x+b.x;this.y=a.y+b.y;this.z=a.z+b.z;return this},addSelf:function(a){this.x=this.x+a.x;this.y=this.y+a.y;this.z=this.z+a.z;return this},addScalar:function(a){this.x=this.x+a;this.y=this.y+
a;this.z=this.z+a;return this},sub:function(a,b){this.x=a.x-b.x;this.y=a.y-b.y;this.z=a.z-b.z;return this},subSelf:function(a){this.x=this.x-a.x;this.y=this.y-a.y;this.z=this.z-a.z;return this},multiply:function(a,b){this.x=a.x*b.x;this.y=a.y*b.y;this.z=a.z*b.z;return this},multiplySelf:function(a){this.x=this.x*a.x;this.y=this.y*a.y;this.z=this.z*a.z;return this},multiplyScalar:function(a){this.x=this.x*a;this.y=this.y*a;this.z=this.z*a;return this},divideSelf:function(a){this.x=this.x/a.x;this.y=
this.y/a.y;this.z=this.z/a.z;return this},divideScalar:function(a){if(a){this.x=this.x/a;this.y=this.y/a;this.z=this.z/a}else this.z=this.y=this.x=0;return this},negate:function(){return this.multiplyScalar(-1)},dot:function(a){return this.x*a.x+this.y*a.y+this.z*a.z},lengthSq:function(){return this.x*this.x+this.y*this.y+this.z*this.z},length:function(){return Math.sqrt(this.lengthSq())},lengthManhattan:function(){return Math.abs(this.x)+Math.abs(this.y)+Math.abs(this.z)},normalize:function(){return this.divideScalar(this.length())},
setLength:function(a){return this.normalize().multiplyScalar(a)},lerpSelf:function(a,b){this.x=this.x+(a.x-this.x)*b;this.y=this.y+(a.y-this.y)*b;this.z=this.z+(a.z-this.z)*b;return this},cross:function(a,b){this.x=a.y*b.z-a.z*b.y;this.y=a.z*b.x-a.x*b.z;this.z=a.x*b.y-a.y*b.x;return this},crossSelf:function(a){var b=this.x,c=this.y,d=this.z;this.x=c*a.z-d*a.y;this.y=d*a.x-b*a.z;this.z=b*a.y-c*a.x;return this},distanceTo:function(a){return Math.sqrt(this.distanceToSquared(a))},distanceToSquared:function(a){return(new THREE.Vector3).sub(this,
a).lengthSq()},getPositionFromMatrix:function(a){this.x=a.elements[12];this.y=a.elements[13];this.z=a.elements[14];return this},getRotationFromMatrix:function(a,b){var c=b?b.x:1,d=b?b.y:1,e=b?b.z:1,f=a.elements[0]/c,h=a.elements[4]/d,c=a.elements[1]/c,d=a.elements[5]/d,i=a.elements[9]/e,l=a.elements[10]/e;this.y=Math.asin(a.elements[8]/e);e=Math.cos(this.y);if(Math.abs(e)>1.0E-5){this.x=Math.atan2(-i/e,l/e);this.z=Math.atan2(-h/e,f/e)}else{this.x=0;this.z=Math.atan2(c,d)}return this},getScaleFromMatrix:function(a){var b=
this.set(a.elements[0],a.elements[1],a.elements[2]).length(),c=this.set(a.elements[4],a.elements[5],a.elements[6]).length(),a=this.set(a.elements[8],a.elements[9],a.elements[10]).length();this.x=b;this.y=c;this.z=a},equals:function(a){return a.x===this.x&&a.y===this.y&&a.z===this.z},isZero:function(){return this.lengthSq()<1.0E-4},clone:function(){return new THREE.Vector3(this.x,this.y,this.z)}};THREE.Vector4=function(a,b,c,d){this.x=a||0;this.y=b||0;this.z=c||0;this.w=d!==void 0?d:1};
THREE.Vector4.prototype={constructor:THREE.Vector4,set:function(a,b,c,d){this.x=a;this.y=b;this.z=c;this.w=d;return this},copy:function(a){this.x=a.x;this.y=a.y;this.z=a.z;this.w=a.w!==void 0?a.w:1;return this},add:function(a,b){this.x=a.x+b.x;this.y=a.y+b.y;this.z=a.z+b.z;this.w=a.w+b.w;return this},addSelf:function(a){this.x=this.x+a.x;this.y=this.y+a.y;this.z=this.z+a.z;this.w=this.w+a.w;return this},sub:function(a,b){this.x=a.x-b.x;this.y=a.y-b.y;this.z=a.z-b.z;this.w=a.w-b.w;return this},subSelf:function(a){this.x=
this.x-a.x;this.y=this.y-a.y;this.z=this.z-a.z;this.w=this.w-a.w;return this},multiplyScalar:function(a){this.x=this.x*a;this.y=this.y*a;this.z=this.z*a;this.w=this.w*a;return this},divideScalar:function(a){if(a){this.x=this.x/a;this.y=this.y/a;this.z=this.z/a;this.w=this.w/a}else{this.z=this.y=this.x=0;this.w=1}return this},negate:function(){return this.multiplyScalar(-1)},dot:function(a){return this.x*a.x+this.y*a.y+this.z*a.z+this.w*a.w},lengthSq:function(){return this.dot(this)},length:function(){return Math.sqrt(this.lengthSq())},
normalize:function(){return this.divideScalar(this.length())},setLength:function(a){return this.normalize().multiplyScalar(a)},lerpSelf:function(a,b){this.x=this.x+(a.x-this.x)*b;this.y=this.y+(a.y-this.y)*b;this.z=this.z+(a.z-this.z)*b;this.w=this.w+(a.w-this.w)*b;return this},clone:function(){return new THREE.Vector4(this.x,this.y,this.z,this.w)}};THREE.Frustum=function(){this.planes=[new THREE.Vector4,new THREE.Vector4,new THREE.Vector4,new THREE.Vector4,new THREE.Vector4,new THREE.Vector4]};
THREE.Frustum.prototype.setFromMatrix=function(a){var b,c=this.planes,d=a.elements,a=d[0];b=d[1];var e=d[2],f=d[3],h=d[4],i=d[5],l=d[6],k=d[7],j=d[8],m=d[9],n=d[10],o=d[11],s=d[12],p=d[13],q=d[14],d=d[15];c[0].set(f-a,k-h,o-j,d-s);c[1].set(f+a,k+h,o+j,d+s);c[2].set(f+b,k+i,o+m,d+p);c[3].set(f-b,k-i,o-m,d-p);c[4].set(f-e,k-l,o-n,d-q);c[5].set(f+e,k+l,o+n,d+q);for(a=0;a<6;a++){b=c[a];b.divideScalar(Math.sqrt(b.x*b.x+b.y*b.y+b.z*b.z))}};
THREE.Frustum.prototype.contains=function(a){for(var b=this.planes,c=a.matrixWorld,d=c.elements,c=-a.geometry.boundingSphere.radius*c.getMaxScaleOnAxis(),e=0;e<6;e++){a=b[e].x*d[12]+b[e].y*d[13]+b[e].z*d[14]+b[e].w;if(a<=c)return false}return true};THREE.Frustum.__v1=new THREE.Vector3;
THREE.Ray=function(a,b){function c(a,b,c){s.sub(c,a);w=s.dot(b);A=p.add(a,q.copy(b).multiplyScalar(w));return y=c.distanceTo(A)}function d(a,b,c,d){s.sub(d,b);p.sub(c,b);q.sub(a,b);u=s.dot(s);H=s.dot(p);B=s.dot(q);K=p.dot(p);N=p.dot(q);Y=1/(u*K-H*H);ca=(K*B-H*N)*Y;I=(u*N-H*B)*Y;return ca>=0&&I>=0&&ca+I<1}this.origin=a||new THREE.Vector3;this.direction=b||new THREE.Vector3;var e=1.0E-4;this.setPrecision=function(a){e=a};var f=new THREE.Vector3,h=new THREE.Vector3,i=new THREE.Vector3,l=new THREE.Vector3,
k=new THREE.Vector3,j=new THREE.Vector3,m=new THREE.Vector3,n=new THREE.Vector3,o=new THREE.Vector3;this.intersectObject=function(a){var b,s=[];if(a instanceof THREE.Particle){var p=c(this.origin,this.direction,a.matrixWorld.getPosition());if(p>a.scale.x)return[];b={distance:p,point:a.position,face:null,object:a};s.push(b)}else if(a instanceof THREE.Mesh){var p=c(this.origin,this.direction,a.matrixWorld.getPosition()),g=THREE.Frustum.__v1.set(a.matrixWorld.getColumnX().length(),a.matrixWorld.getColumnY().length(),
a.matrixWorld.getColumnZ().length());if(p>a.geometry.boundingSphere.radius*Math.max(g.x,Math.max(g.y,g.z)))return s;var q,u,y=a.geometry,w=y.vertices,A;a.matrixRotationWorld.extractRotation(a.matrixWorld);p=0;for(g=y.faces.length;p<g;p++){b=y.faces[p];k.copy(this.origin);j.copy(this.direction);A=a.matrixWorld;m=A.multiplyVector3(m.copy(b.centroid)).subSelf(k);n=a.matrixRotationWorld.multiplyVector3(n.copy(b.normal));q=j.dot(n);if(!(Math.abs(q)<e)){u=n.dot(m)/q;if(!(u<0)&&(a.doubleSided||(a.flipSided?
q>0:q<0))){o.add(k,j.multiplyScalar(u));if(b instanceof THREE.Face3){f=A.multiplyVector3(f.copy(w[b.a]));h=A.multiplyVector3(h.copy(w[b.b]));i=A.multiplyVector3(i.copy(w[b.c]));if(d(o,f,h,i)){b={distance:k.distanceTo(o),point:o.clone(),face:b,object:a};s.push(b)}}else if(b instanceof THREE.Face4){f=A.multiplyVector3(f.copy(w[b.a]));h=A.multiplyVector3(h.copy(w[b.b]));i=A.multiplyVector3(i.copy(w[b.c]));l=A.multiplyVector3(l.copy(w[b.d]));if(d(o,f,h,l)||d(o,h,i,l)){b={distance:k.distanceTo(o),point:o.clone(),
face:b,object:a};s.push(b)}}}}}}return s};this.intersectObjects=function(a){for(var b=[],c=0,d=a.length;c<d;c++)Array.prototype.push.apply(b,this.intersectObject(a[c]));b.sort(function(a,b){return a.distance-b.distance});return b};var s=new THREE.Vector3,p=new THREE.Vector3,q=new THREE.Vector3,w,A,y,u,H,B,K,N,Y,ca,I};
THREE.Rectangle=function(){function a(){f=d-b;h=e-c}var b,c,d,e,f,h,i=true;this.getX=function(){return b};this.getY=function(){return c};this.getWidth=function(){return f};this.getHeight=function(){return h};this.getLeft=function(){return b};this.getTop=function(){return c};this.getRight=function(){return d};this.getBottom=function(){return e};this.set=function(f,h,j,m){i=false;b=f;c=h;d=j;e=m;a()};this.addPoint=function(f,h){if(i){i=false;b=f;c=h;d=f;e=h}else{b=b<f?b:f;c=c<h?c:h;d=d>f?d:f;e=e>h?
e:h}a()};this.add3Points=function(f,h,j,m,n,o){if(i){i=false;b=f<j?f<n?f:n:j<n?j:n;c=h<m?h<o?h:o:m<o?m:o;d=f>j?f>n?f:n:j>n?j:n;e=h>m?h>o?h:o:m>o?m:o}else{b=f<j?f<n?f<b?f:b:n<b?n:b:j<n?j<b?j:b:n<b?n:b;c=h<m?h<o?h<c?h:c:o<c?o:c:m<o?m<c?m:c:o<c?o:c;d=f>j?f>n?f>d?f:d:n>d?n:d:j>n?j>d?j:d:n>d?n:d;e=h>m?h>o?h>e?h:e:o>e?o:e:m>o?m>e?m:e:o>e?o:e}a()};this.addRectangle=function(f){if(i){i=false;b=f.getLeft();c=f.getTop();d=f.getRight();e=f.getBottom()}else{b=b<f.getLeft()?b:f.getLeft();c=c<f.getTop()?c:f.getTop();
d=d>f.getRight()?d:f.getRight();e=e>f.getBottom()?e:f.getBottom()}a()};this.inflate=function(f){b=b-f;c=c-f;d=d+f;e=e+f;a()};this.minSelf=function(f){b=b>f.getLeft()?b:f.getLeft();c=c>f.getTop()?c:f.getTop();d=d<f.getRight()?d:f.getRight();e=e<f.getBottom()?e:f.getBottom();a()};this.intersects=function(a){return d<a.getLeft()||b>a.getRight()||e<a.getTop()||c>a.getBottom()?false:true};this.empty=function(){i=true;e=d=c=b=0;a()};this.isEmpty=function(){return i}};
THREE.Math={clamp:function(a,b,c){return a<b?b:a>c?c:a},clampBottom:function(a,b){return a<b?b:a},mapLinear:function(a,b,c,d,e){return d+(a-b)*(e-d)/(c-b)},random16:function(){return(65280*Math.random()+255*Math.random())/65535},randInt:function(a,b){return a+Math.floor(Math.random()*(b-a+1))},randFloat:function(a,b){return a+Math.random()*(b-a)},randFloatSpread:function(a){return a*(0.5-Math.random())},sign:function(a){return a<0?-1:a>0?1:0}};THREE.Matrix3=function(){this.elements=new Float32Array(9)};
THREE.Matrix3.prototype={constructor:THREE.Matrix3,getInverse:function(a){var b=a.elements,a=b[10]*b[5]-b[6]*b[9],c=-b[10]*b[1]+b[2]*b[9],d=b[6]*b[1]-b[2]*b[5],e=-b[10]*b[4]+b[6]*b[8],f=b[10]*b[0]-b[2]*b[8],h=-b[6]*b[0]+b[2]*b[4],i=b[9]*b[4]-b[5]*b[8],l=-b[9]*b[0]+b[1]*b[8],k=b[5]*b[0]-b[1]*b[4],b=b[0]*a+b[1]*e+b[2]*i;b===0&&console.warn("Matrix3.getInverse(): determinant == 0");var b=1/b,j=this.elements;j[0]=b*a;j[1]=b*c;j[2]=b*d;j[3]=b*e;j[4]=b*f;j[5]=b*h;j[6]=b*i;j[7]=b*l;j[8]=b*k;return this},
transpose:function(){var a,b=this.elements;a=b[1];b[1]=b[3];b[3]=a;a=b[2];b[2]=b[6];b[6]=a;a=b[5];b[5]=b[7];b[7]=a;return this},transposeIntoArray:function(a){var b=this.m;a[0]=b[0];a[1]=b[3];a[2]=b[6];a[3]=b[1];a[4]=b[4];a[5]=b[7];a[6]=b[2];a[7]=b[5];a[8]=b[8];return this}};THREE.Matrix4=function(a,b,c,d,e,f,h,i,l,k,j,m,n,o,s,p){this.elements=new Float32Array(16);this.set(a!==void 0?a:1,b||0,c||0,d||0,e||0,f!==void 0?f:1,h||0,i||0,l||0,k||0,j!==void 0?j:1,m||0,n||0,o||0,s||0,p!==void 0?p:1)};
THREE.Matrix4.prototype={constructor:THREE.Matrix4,set:function(a,b,c,d,e,f,h,i,l,k,j,m,n,o,s,p){var q=this.elements;q[0]=a;q[4]=b;q[8]=c;q[12]=d;q[1]=e;q[5]=f;q[9]=h;q[13]=i;q[2]=l;q[6]=k;q[10]=j;q[14]=m;q[3]=n;q[7]=o;q[11]=s;q[15]=p;return this},identity:function(){this.set(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);return this},copy:function(a){a=a.elements;this.set(a[0],a[4],a[8],a[12],a[1],a[5],a[9],a[13],a[2],a[6],a[10],a[14],a[3],a[7],a[11],a[15]);return this},lookAt:function(a,b,c){var d=this.elements,
e=THREE.Matrix4.__v1,f=THREE.Matrix4.__v2,h=THREE.Matrix4.__v3;h.sub(a,b).normalize();if(h.length()===0)h.z=1;e.cross(c,h).normalize();if(e.length()===0){h.x=h.x+1.0E-4;e.cross(c,h).normalize()}f.cross(h,e);d[0]=e.x;d[4]=f.x;d[8]=h.x;d[1]=e.y;d[5]=f.y;d[9]=h.y;d[2]=e.z;d[6]=f.z;d[10]=h.z;return this},multiply:function(a,b){var c=a.elements,d=b.elements,e=this.elements,f=c[0],h=c[4],i=c[8],l=c[12],k=c[1],j=c[5],m=c[9],n=c[13],o=c[2],s=c[6],p=c[10],q=c[14],w=c[3],A=c[7],y=c[11],c=c[15],u=d[0],H=d[4],
B=d[8],K=d[12],N=d[1],Y=d[5],ca=d[9],I=d[13],ba=d[2],ja=d[6],ya=d[10],D=d[14],g=d[3],Na=d[7],za=d[11],d=d[15];e[0]=f*u+h*N+i*ba+l*g;e[4]=f*H+h*Y+i*ja+l*Na;e[8]=f*B+h*ca+i*ya+l*za;e[12]=f*K+h*I+i*D+l*d;e[1]=k*u+j*N+m*ba+n*g;e[5]=k*H+j*Y+m*ja+n*Na;e[9]=k*B+j*ca+m*ya+n*za;e[13]=k*K+j*I+m*D+n*d;e[2]=o*u+s*N+p*ba+q*g;e[6]=o*H+s*Y+p*ja+q*Na;e[10]=o*B+s*ca+p*ya+q*za;e[14]=o*K+s*I+p*D+q*d;e[3]=w*u+A*N+y*ba+c*g;e[7]=w*H+A*Y+y*ja+c*Na;e[11]=w*B+A*ca+y*ya+c*za;e[15]=w*K+A*I+y*D+c*d;return this},multiplySelf:function(a){return this.multiply(this,
a)},multiplyToArray:function(a,b,c){var d=this.elements;this.multiply(a,b);c[0]=d[0];c[1]=d[1];c[2]=d[2];c[3]=d[3];c[4]=d[4];c[5]=d[5];c[6]=d[6];c[7]=d[7];c[8]=d[8];c[9]=d[9];c[10]=d[10];c[11]=d[11];c[12]=d[12];c[13]=d[13];c[14]=d[14];c[15]=d[15];return this},multiplyScalar:function(a){var b=this.elements;b[0]=b[0]*a;b[4]=b[4]*a;b[8]=b[8]*a;b[12]=b[12]*a;b[1]=b[1]*a;b[5]=b[5]*a;b[9]=b[9]*a;b[13]=b[13]*a;b[2]=b[2]*a;b[6]=b[6]*a;b[10]=b[10]*a;b[14]=b[14]*a;b[3]=b[3]*a;b[7]=b[7]*a;b[11]=b[11]*a;b[15]=
b[15]*a;return this},multiplyVector3:function(a){var b=this.elements,c=a.x,d=a.y,e=a.z,f=1/(b[3]*c+b[7]*d+b[11]*e+b[15]);a.x=(b[0]*c+b[4]*d+b[8]*e+b[12])*f;a.y=(b[1]*c+b[5]*d+b[9]*e+b[13])*f;a.z=(b[2]*c+b[6]*d+b[10]*e+b[14])*f;return a},multiplyVector4:function(a){var b=this.elements,c=a.x,d=a.y,e=a.z,f=a.w;a.x=b[0]*c+b[4]*d+b[8]*e+b[12]*f;a.y=b[1]*c+b[5]*d+b[9]*e+b[13]*f;a.z=b[2]*c+b[6]*d+b[10]*e+b[14]*f;a.w=b[3]*c+b[7]*d+b[11]*e+b[15]*f;return a},rotateAxis:function(a){var b=this.elements,c=a.x,
d=a.y,e=a.z;a.x=c*b[0]+d*b[4]+e*b[8];a.y=c*b[1]+d*b[5]+e*b[9];a.z=c*b[2]+d*b[6]+e*b[10];a.normalize();return a},crossVector:function(a){var b=this.elements,c=new THREE.Vector4;c.x=b[0]*a.x+b[4]*a.y+b[8]*a.z+b[12]*a.w;c.y=b[1]*a.x+b[5]*a.y+b[9]*a.z+b[13]*a.w;c.z=b[2]*a.x+b[6]*a.y+b[10]*a.z+b[14]*a.w;c.w=a.w?b[3]*a.x+b[7]*a.y+b[11]*a.z+b[15]*a.w:1;return c},determinant:function(){var a=this.elements,b=a[0],c=a[4],d=a[8],e=a[12],f=a[1],h=a[5],i=a[9],l=a[13],k=a[2],j=a[6],m=a[10],n=a[14],o=a[3],s=a[7],
p=a[11],a=a[15];return e*i*j*o-d*l*j*o-e*h*m*o+c*l*m*o+d*h*n*o-c*i*n*o-e*i*k*s+d*l*k*s+e*f*m*s-b*l*m*s-d*f*n*s+b*i*n*s+e*h*k*p-c*l*k*p-e*f*j*p+b*l*j*p+c*f*n*p-b*h*n*p-d*h*k*a+c*i*k*a+d*f*j*a-b*i*j*a-c*f*m*a+b*h*m*a},transpose:function(){var a=this.elements,b;b=a[1];a[1]=a[4];a[4]=b;b=a[2];a[2]=a[8];a[8]=b;b=a[6];a[6]=a[9];a[9]=b;b=a[3];a[3]=a[12];a[12]=b;b=a[7];a[7]=a[13];a[13]=b;b=a[11];a[11]=a[14];a[14]=b;return this},flattenToArray:function(a){var b=this.elements;a[0]=b[0];a[1]=b[1];a[2]=b[2];
a[3]=b[3];a[4]=b[4];a[5]=b[5];a[6]=b[6];a[7]=b[7];a[8]=b[8];a[9]=b[9];a[10]=b[10];a[11]=b[11];a[12]=b[12];a[13]=b[13];a[14]=b[14];a[15]=b[15];return a},flattenToArrayOffset:function(a,b){var c=this.elements;a[b]=c[0];a[b+1]=c[1];a[b+2]=c[2];a[b+3]=c[3];a[b+4]=c[4];a[b+5]=c[5];a[b+6]=c[6];a[b+7]=c[7];a[b+8]=c[8];a[b+9]=c[9];a[b+10]=c[10];a[b+11]=c[11];a[b+12]=c[12];a[b+13]=c[13];a[b+14]=c[14];a[b+15]=c[15];return a},getPosition:function(){var a=this.elements;return THREE.Matrix4.__v1.set(a[12],a[13],
a[14])},setPosition:function(a){var b=this.elements;b[12]=a.x;b[13]=a.y;b[14]=a.z;return this},getColumnX:function(){var a=this.elements;return THREE.Matrix4.__v1.set(a[0],a[1],a[2])},getColumnY:function(){var a=this.elements;return THREE.Matrix4.__v1.set(a[4],a[5],a[6])},getColumnZ:function(){var a=this.elements;return THREE.Matrix4.__v1.set(a[8],a[9],a[10])},getInverse:function(a){var b=this.elements,c=a.elements,d=c[0],e=c[4],f=c[8],h=c[12],i=c[1],l=c[5],k=c[9],j=c[13],m=c[2],n=c[6],o=c[10],s=
c[14],p=c[3],q=c[7],w=c[11],c=c[15];b[0]=k*s*q-j*o*q+j*n*w-l*s*w-k*n*c+l*o*c;b[4]=h*o*q-f*s*q-h*n*w+e*s*w+f*n*c-e*o*c;b[8]=f*j*q-h*k*q+h*l*w-e*j*w-f*l*c+e*k*c;b[12]=h*k*n-f*j*n-h*l*o+e*j*o+f*l*s-e*k*s;b[1]=j*o*p-k*s*p-j*m*w+i*s*w+k*m*c-i*o*c;b[5]=f*s*p-h*o*p+h*m*w-d*s*w-f*m*c+d*o*c;b[9]=h*k*p-f*j*p-h*i*w+d*j*w+f*i*c-d*k*c;b[13]=f*j*m-h*k*m+h*i*o-d*j*o-f*i*s+d*k*s;b[2]=l*s*p-j*n*p+j*m*q-i*s*q-l*m*c+i*n*c;b[6]=h*n*p-e*s*p-h*m*q+d*s*q+e*m*c-d*n*c;b[10]=e*j*p-h*l*p+h*i*q-d*j*q-e*i*c+d*l*c;b[14]=h*l*m-
e*j*m-h*i*n+d*j*n+e*i*s-d*l*s;b[3]=k*n*p-l*o*p-k*m*q+i*o*q+l*m*w-i*n*w;b[7]=e*o*p-f*n*p+f*m*q-d*o*q-e*m*w+d*n*w;b[11]=f*l*p-e*k*p-f*i*q+d*k*q+e*i*w-d*l*w;b[15]=e*k*m-f*l*m+f*i*n-d*k*n-e*i*o+d*l*o;this.multiplyScalar(1/a.determinant());return this},setRotationFromEuler:function(a,b){var c=this.elements,d=a.x,e=a.y,f=a.z,h=Math.cos(d),d=Math.sin(d),i=Math.cos(e),e=Math.sin(e),l=Math.cos(f),f=Math.sin(f);switch(b){case "YXZ":var k=i*l,j=i*f,m=e*l,n=e*f;c[0]=k+n*d;c[4]=m*d-j;c[8]=h*e;c[1]=h*f;c[5]=h*
l;c[9]=-d;c[2]=j*d-m;c[6]=n+k*d;c[10]=h*i;break;case "ZXY":k=i*l;j=i*f;m=e*l;n=e*f;c[0]=k-n*d;c[4]=-h*f;c[8]=m+j*d;c[1]=j+m*d;c[5]=h*l;c[9]=n-k*d;c[2]=-h*e;c[6]=d;c[10]=h*i;break;case "ZYX":k=h*l;j=h*f;m=d*l;n=d*f;c[0]=i*l;c[4]=m*e-j;c[8]=k*e+n;c[1]=i*f;c[5]=n*e+k;c[9]=j*e-m;c[2]=-e;c[6]=d*i;c[10]=h*i;break;case "YZX":k=h*i;j=h*e;m=d*i;n=d*e;c[0]=i*l;c[4]=n-k*f;c[8]=m*f+j;c[1]=f;c[5]=h*l;c[9]=-d*l;c[2]=-e*l;c[6]=j*f+m;c[10]=k-n*f;break;case "XZY":k=h*i;j=h*e;m=d*i;n=d*e;c[0]=i*l;c[4]=-f;c[8]=e*l;
c[1]=k*f+n;c[5]=h*l;c[9]=j*f-m;c[2]=m*f-j;c[6]=d*l;c[10]=n*f+k;break;default:k=h*l;j=h*f;m=d*l;n=d*f;c[0]=i*l;c[4]=-i*f;c[8]=e;c[1]=j+m*e;c[5]=k-n*e;c[9]=-d*i;c[2]=n-k*e;c[6]=m+j*e;c[10]=h*i}return this},setRotationFromQuaternion:function(a){var b=this.elements,c=a.x,d=a.y,e=a.z,f=a.w,h=c+c,i=d+d,l=e+e,a=c*h,k=c*i,c=c*l,j=d*i,d=d*l,e=e*l,h=f*h,i=f*i,f=f*l;b[0]=1-(j+e);b[4]=k-f;b[8]=c+i;b[1]=k+f;b[5]=1-(a+e);b[9]=d-h;b[2]=c-i;b[6]=d+h;b[10]=1-(a+j);return this},compose:function(a,b,c){var d=this.elements,
e=THREE.Matrix4.__m1,f=THREE.Matrix4.__m2;e.identity();e.setRotationFromQuaternion(b);f.makeScale(c.x,c.y,c.z);this.multiply(e,f);d[12]=a.x;d[13]=a.y;d[14]=a.z;return this},decompose:function(a,b,c){var d=this.elements,e=THREE.Matrix4.__v1,f=THREE.Matrix4.__v2,h=THREE.Matrix4.__v3;e.set(d[0],d[1],d[2]);f.set(d[4],d[5],d[6]);h.set(d[8],d[9],d[10]);a=a instanceof THREE.Vector3?a:new THREE.Vector3;b=b instanceof THREE.Quaternion?b:new THREE.Quaternion;c=c instanceof THREE.Vector3?c:new THREE.Vector3;
c.x=e.length();c.y=f.length();c.z=h.length();a.x=d[12];a.y=d[13];a.z=d[14];d=THREE.Matrix4.__m1;d.copy(this);d.elements[0]=d.elements[0]/c.x;d.elements[1]=d.elements[1]/c.x;d.elements[2]=d.elements[2]/c.x;d.elements[4]=d.elements[4]/c.y;d.elements[5]=d.elements[5]/c.y;d.elements[6]=d.elements[6]/c.y;d.elements[8]=d.elements[8]/c.z;d.elements[9]=d.elements[9]/c.z;d.elements[10]=d.elements[10]/c.z;b.setFromRotationMatrix(d);return[a,b,c]},extractPosition:function(a){var b=this.elements,a=a.elements;
b[12]=a[12];b[13]=a[13];b[14]=a[14];return this},extractRotation:function(a){var b=this.elements,a=a.elements,c=THREE.Matrix4.__v1,d=1/c.set(a[0],a[1],a[2]).length(),e=1/c.set(a[4],a[5],a[6]).length(),c=1/c.set(a[8],a[9],a[10]).length();b[0]=a[0]*d;b[1]=a[1]*d;b[2]=a[2]*d;b[4]=a[4]*e;b[5]=a[5]*e;b[6]=a[6]*e;b[8]=a[8]*c;b[9]=a[9]*c;b[10]=a[10]*c;return this},translate:function(a){var b=this.elements,c=a.x,d=a.y,a=a.z;b[12]=b[0]*c+b[4]*d+b[8]*a+b[12];b[13]=b[1]*c+b[5]*d+b[9]*a+b[13];b[14]=b[2]*c+b[6]*
d+b[10]*a+b[14];b[15]=b[3]*c+b[7]*d+b[11]*a+b[15];return this},rotateX:function(a){var b=this.elements,c=b[4],d=b[5],e=b[6],f=b[7],h=b[8],i=b[9],l=b[10],k=b[11],j=Math.cos(a),a=Math.sin(a);b[4]=j*c+a*h;b[5]=j*d+a*i;b[6]=j*e+a*l;b[7]=j*f+a*k;b[8]=j*h-a*c;b[9]=j*i-a*d;b[10]=j*l-a*e;b[11]=j*k-a*f;return this},rotateY:function(a){var b=this.elements,c=b[0],d=b[1],e=b[2],f=b[3],h=b[8],i=b[9],l=b[10],k=b[11],j=Math.cos(a),a=Math.sin(a);b[0]=j*c-a*h;b[1]=j*d-a*i;b[2]=j*e-a*l;b[3]=j*f-a*k;b[8]=j*h+a*c;b[9]=
j*i+a*d;b[10]=j*l+a*e;b[11]=j*k+a*f;return this},rotateZ:function(a){var b=this.elements,c=b[0],d=b[1],e=b[2],f=b[3],h=b[4],i=b[5],l=b[6],k=b[7],j=Math.cos(a),a=Math.sin(a);b[0]=j*c+a*h;b[1]=j*d+a*i;b[2]=j*e+a*l;b[3]=j*f+a*k;b[4]=j*h-a*c;b[5]=j*i-a*d;b[6]=j*l-a*e;b[7]=j*k-a*f;return this},rotateByAxis:function(a,b){var c=this.elements;if(a.x===1&&a.y===0&&a.z===0)return this.rotateX(b);if(a.x===0&&a.y===1&&a.z===0)return this.rotateY(b);if(a.x===0&&a.y===0&&a.z===1)return this.rotateZ(b);var d=a.x,
e=a.y,f=a.z,h=Math.sqrt(d*d+e*e+f*f),d=d/h,e=e/h,f=f/h,h=d*d,i=e*e,l=f*f,k=Math.cos(b),j=Math.sin(b),m=1-k,n=d*e*m,o=d*f*m,m=e*f*m,d=d*j,s=e*j,j=f*j,f=h+(1-h)*k,h=n+j,e=o-s,n=n-j,i=i+(1-i)*k,j=m+d,o=o+s,m=m-d,l=l+(1-l)*k,k=c[0],d=c[1],s=c[2],p=c[3],q=c[4],w=c[5],A=c[6],y=c[7],u=c[8],H=c[9],B=c[10],K=c[11];c[0]=f*k+h*q+e*u;c[1]=f*d+h*w+e*H;c[2]=f*s+h*A+e*B;c[3]=f*p+h*y+e*K;c[4]=n*k+i*q+j*u;c[5]=n*d+i*w+j*H;c[6]=n*s+i*A+j*B;c[7]=n*p+i*y+j*K;c[8]=o*k+m*q+l*u;c[9]=o*d+m*w+l*H;c[10]=o*s+m*A+l*B;c[11]=
o*p+m*y+l*K;return this},scale:function(a){var b=this.elements,c=a.x,d=a.y,a=a.z;b[0]=b[0]*c;b[4]=b[4]*d;b[8]=b[8]*a;b[1]=b[1]*c;b[5]=b[5]*d;b[9]=b[9]*a;b[2]=b[2]*c;b[6]=b[6]*d;b[10]=b[10]*a;b[3]=b[3]*c;b[7]=b[7]*d;b[11]=b[11]*a;return this},getMaxScaleOnAxis:function(){var a=this.elements;return Math.sqrt(Math.max(a[0]*a[0]+a[1]*a[1]+a[2]*a[2],Math.max(a[4]*a[4]+a[5]*a[5]+a[6]*a[6],a[8]*a[8]+a[9]*a[9]+a[10]*a[10])))},makeTranslation:function(a,b,c){this.set(1,0,0,a,0,1,0,b,0,0,1,c,0,0,0,1);return this},
makeRotationX:function(a){var b=Math.cos(a),a=Math.sin(a);this.set(1,0,0,0,0,b,-a,0,0,a,b,0,0,0,0,1);return this},makeRotationY:function(a){var b=Math.cos(a),a=Math.sin(a);this.set(b,0,a,0,0,1,0,0,-a,0,b,0,0,0,0,1);return this},makeRotationZ:function(a){var b=Math.cos(a),a=Math.sin(a);this.set(b,-a,0,0,a,b,0,0,0,0,1,0,0,0,0,1);return this},makeRotationAxis:function(a,b){var c=Math.cos(b),d=Math.sin(b),e=1-c,f=a.x,h=a.y,i=a.z,l=e*f,k=e*h;this.set(l*f+c,l*h-d*i,l*i+d*h,0,l*h+d*i,k*h+c,k*i-d*f,0,l*i-
d*h,k*i+d*f,e*i*i+c,0,0,0,0,1);return this},makeScale:function(a,b,c){this.set(a,0,0,0,0,b,0,0,0,0,c,0,0,0,0,1);return this},makeFrustum:function(a,b,c,d,e,f){var h=this.elements;h[0]=2*e/(b-a);h[4]=0;h[8]=(b+a)/(b-a);h[12]=0;h[1]=0;h[5]=2*e/(d-c);h[9]=(d+c)/(d-c);h[13]=0;h[2]=0;h[6]=0;h[10]=-(f+e)/(f-e);h[14]=-2*f*e/(f-e);h[3]=0;h[7]=0;h[11]=-1;h[15]=0;return this},makePerspective:function(a,b,c,d){var a=c*Math.tan(a*Math.PI/360),e=-a;return this.makeFrustum(e*b,a*b,e,a,c,d)},makeOrthographic:function(a,
b,c,d,e,f){var h=this.elements,i=b-a,l=c-d,k=f-e;h[0]=2/i;h[4]=0;h[8]=0;h[12]=-((b+a)/i);h[1]=0;h[5]=2/l;h[9]=0;h[13]=-((c+d)/l);h[2]=0;h[6]=0;h[10]=-2/k;h[14]=-((f+e)/k);h[3]=0;h[7]=0;h[11]=0;h[15]=1;return this},clone:function(){var a=this.elements;return new THREE.Matrix4(a[0],a[4],a[8],a[12],a[1],a[5],a[9],a[13],a[2],a[6],a[10],a[14],a[3],a[7],a[11],a[15])}};THREE.Matrix4.__v1=new THREE.Vector3;THREE.Matrix4.__v2=new THREE.Vector3;THREE.Matrix4.__v3=new THREE.Vector3;THREE.Matrix4.__m1=new THREE.Matrix4;
THREE.Matrix4.__m2=new THREE.Matrix4;
THREE.Object3D=function(){this.id=THREE.Object3DCount++;this.name="";this.parent=void 0;this.children=[];this.up=new THREE.Vector3(0,1,0);this.position=new THREE.Vector3;this.rotation=new THREE.Vector3;this.eulerOrder="XYZ";this.scale=new THREE.Vector3(1,1,1);this.flipSided=this.doubleSided=false;this.renderDepth=null;this.rotationAutoUpdate=true;this.matrix=new THREE.Matrix4;this.matrixWorld=new THREE.Matrix4;this.matrixRotationWorld=new THREE.Matrix4;this.matrixWorldNeedsUpdate=this.matrixAutoUpdate=
true;this.quaternion=new THREE.Quaternion;this.useQuaternion=false;this.boundRadius=0;this.boundRadiusScale=1;this.visible=true;this.receiveShadow=this.castShadow=false;this.frustumCulled=true;this._vector=new THREE.Vector3};
THREE.Object3D.prototype={constructor:THREE.Object3D,applyMatrix:function(a){this.matrix.multiply(a,this.matrix);this.scale.getScaleFromMatrix(this.matrix);this.rotation.getRotationFromMatrix(this.matrix,this.scale);this.position.getPositionFromMatrix(this.matrix)},translate:function(a,b){this.matrix.rotateAxis(b);this.position.addSelf(b.multiplyScalar(a))},translateX:function(a){this.translate(a,this._vector.set(1,0,0))},translateY:function(a){this.translate(a,this._vector.set(0,1,0))},translateZ:function(a){this.translate(a,
this._vector.set(0,0,1))},lookAt:function(a){this.matrix.lookAt(a,this.position,this.up);this.rotationAutoUpdate&&this.rotation.getRotationFromMatrix(this.matrix)},add:function(a){if(a===this)console.warn("THREE.Object3D.add: An object can't be added as a child of itself.");else if(a instanceof THREE.Object3D){a.parent!==void 0&&a.parent.remove(a);a.parent=this;this.children.push(a);for(var b=this;b.parent!==void 0;)b=b.parent;b!==void 0&&b instanceof THREE.Scene&&b.__addObject(a)}},remove:function(a){var b=
this.children.indexOf(a);if(b!==-1){a.parent=void 0;this.children.splice(b,1);for(b=this;b.parent!==void 0;)b=b.parent;b!==void 0&&b instanceof THREE.Scene&&b.__removeObject(a)}},getChildByName:function(a,b){var c,d,e;c=0;for(d=this.children.length;c<d;c++){e=this.children[c];if(e.name===a)return e;if(b){e=e.getChildByName(a,b);if(e!==void 0)return e}}},updateMatrix:function(){this.matrix.setPosition(this.position);this.useQuaternion?this.matrix.setRotationFromQuaternion(this.quaternion):this.matrix.setRotationFromEuler(this.rotation,
this.eulerOrder);if(this.scale.x!==1||this.scale.y!==1||this.scale.z!==1){this.matrix.scale(this.scale);this.boundRadiusScale=Math.max(this.scale.x,Math.max(this.scale.y,this.scale.z))}this.matrixWorldNeedsUpdate=true},updateMatrixWorld:function(a){this.matrixAutoUpdate&&this.updateMatrix();if(this.matrixWorldNeedsUpdate||a){this.parent?this.matrixWorld.multiply(this.parent.matrixWorld,this.matrix):this.matrixWorld.copy(this.matrix);this.matrixWorldNeedsUpdate=false;a=true}for(var b=0,c=this.children.length;b<
c;b++)this.children[b].updateMatrixWorld(a)}};THREE.Object3DCount=0;
THREE.Projector=function(){function a(){var a=l[i]=l[i]||new THREE.RenderableVertex;i++;return a}function b(a,b){return b.z-a.z}function c(a,b){var c=0,d=1,g=a.z+a.w,e=b.z+b.w,f=-a.z+a.w,h=-b.z+b.w;if(g>=0&&e>=0&&f>=0&&h>=0)return true;if(g<0&&e<0||f<0&&h<0)return false;g<0?c=Math.max(c,g/(g-e)):e<0&&(d=Math.min(d,g/(g-e)));f<0?c=Math.max(c,f/(f-h)):h<0&&(d=Math.min(d,f/(f-h)));if(d<c)return false;a.lerpSelf(b,c);b.lerpSelf(a,1-d);return true}var d,e,f=[],h,i,l=[],k,j,m=[],n,o=[],s,p,q=[],w,A,y=[],
u={objects:[],sprites:[],lights:[],elements:[]},H=new THREE.Vector3,B=new THREE.Vector4,K=new THREE.Matrix4,N=new THREE.Matrix4,Y=new THREE.Frustum,ca=new THREE.Vector4,I=new THREE.Vector4;this.projectVector=function(a,b){b.matrixWorldInverse.getInverse(b.matrixWorld);K.multiply(b.projectionMatrix,b.matrixWorldInverse);K.multiplyVector3(a);return a};this.unprojectVector=function(a,b){b.projectionMatrixInverse.getInverse(b.projectionMatrix);K.multiply(b.matrixWorld,b.projectionMatrixInverse);K.multiplyVector3(a);
return a};this.pickingRay=function(a,b){var c;a.z=-1;c=new THREE.Vector3(a.x,a.y,1);this.unprojectVector(a,b);this.unprojectVector(c,b);c.subSelf(a).normalize();return new THREE.Ray(a,c)};this.projectGraph=function(a,c){e=0;u.objects.length=0;u.sprites.length=0;u.lights.length=0;var h=function(a){if(a.visible!==false){if((a instanceof THREE.Mesh||a instanceof THREE.Line)&&(a.frustumCulled===false||Y.contains(a))){H.copy(a.matrixWorld.getPosition());K.multiplyVector3(H);var b=f[e]=f[e]||new THREE.RenderableObject;
e++;d=b;d.object=a;d.z=H.z;u.objects.push(d)}else a instanceof THREE.Light&&u.lights.push(a);for(var b=0,c=a.children.length;b<c;b++)h(a.children[b])}};h(a);c&&u.objects.sort(b);return u};this.projectScene=function(d,e,f){var D=e.near,g=e.far,H=false,za,Da,$,R,J,Z,Q,pa,O,sa,Ga,Ha,La,Oa,Ta;A=p=n=j=0;u.elements.length=0;if(e.parent===void 0){console.warn("DEPRECATED: Camera hasn't been added to a Scene. Adding it...");d.add(e)}d.updateMatrixWorld();e.matrixWorldInverse.getInverse(e.matrixWorld);K.multiply(e.projectionMatrix,
e.matrixWorldInverse);Y.setFromMatrix(K);u=this.projectGraph(d,false);d=0;for(za=u.objects.length;d<za;d++){O=u.objects[d].object;sa=O.matrixWorld;i=0;if(O instanceof THREE.Mesh){Ga=O.geometry;Ha=O.geometry.materials;R=Ga.vertices;La=Ga.faces;Oa=Ga.faceVertexUvs;Ga=O.matrixRotationWorld.extractRotation(sa);Da=0;for($=R.length;Da<$;Da++){h=a();h.positionWorld.copy(R[Da]);sa.multiplyVector3(h.positionWorld);h.positionScreen.copy(h.positionWorld);K.multiplyVector4(h.positionScreen);h.positionScreen.x=
h.positionScreen.x/h.positionScreen.w;h.positionScreen.y=h.positionScreen.y/h.positionScreen.w;h.visible=h.positionScreen.z>D&&h.positionScreen.z<g}R=0;for(Da=La.length;R<Da;R++){$=La[R];if($ instanceof THREE.Face3){J=l[$.a];Z=l[$.b];Q=l[$.c];if(J.visible&&Z.visible&&Q.visible){H=(Q.positionScreen.x-J.positionScreen.x)*(Z.positionScreen.y-J.positionScreen.y)-(Q.positionScreen.y-J.positionScreen.y)*(Z.positionScreen.x-J.positionScreen.x)<0;if(O.doubleSided||H!=O.flipSided){pa=m[j]=m[j]||new THREE.RenderableFace3;
j++;k=pa;k.v1.copy(J);k.v2.copy(Z);k.v3.copy(Q)}else continue}else continue}else if($ instanceof THREE.Face4){J=l[$.a];Z=l[$.b];Q=l[$.c];pa=l[$.d];if(J.visible&&Z.visible&&Q.visible&&pa.visible){H=(pa.positionScreen.x-J.positionScreen.x)*(Z.positionScreen.y-J.positionScreen.y)-(pa.positionScreen.y-J.positionScreen.y)*(Z.positionScreen.x-J.positionScreen.x)<0||(Z.positionScreen.x-Q.positionScreen.x)*(pa.positionScreen.y-Q.positionScreen.y)-(Z.positionScreen.y-Q.positionScreen.y)*(pa.positionScreen.x-
Q.positionScreen.x)<0;if(O.doubleSided||H!=O.flipSided){Ta=o[n]=o[n]||new THREE.RenderableFace4;n++;k=Ta;k.v1.copy(J);k.v2.copy(Z);k.v3.copy(Q);k.v4.copy(pa)}else continue}else continue}k.normalWorld.copy($.normal);!H&&(O.flipSided||O.doubleSided)&&k.normalWorld.negate();Ga.multiplyVector3(k.normalWorld);k.centroidWorld.copy($.centroid);sa.multiplyVector3(k.centroidWorld);k.centroidScreen.copy(k.centroidWorld);K.multiplyVector3(k.centroidScreen);Q=$.vertexNormals;J=0;for(Z=Q.length;J<Z;J++){pa=k.vertexNormalsWorld[J];
pa.copy(Q[J]);!H&&(O.flipSided||O.doubleSided)&&pa.negate();Ga.multiplyVector3(pa)}J=0;for(Z=Oa.length;J<Z;J++)if(Ta=Oa[J][R]){Q=0;for(pa=Ta.length;Q<pa;Q++)k.uvs[J][Q]=Ta[Q]}k.material=O.material;k.faceMaterial=$.materialIndex!==null?Ha[$.materialIndex]:null;k.z=k.centroidScreen.z;u.elements.push(k)}}else if(O instanceof THREE.Line){N.multiply(K,sa);R=O.geometry.vertices;J=a();J.positionScreen.copy(R[0]);N.multiplyVector4(J.positionScreen);sa=O.type===THREE.LinePieces?2:1;Da=1;for($=R.length;Da<
$;Da++){J=a();J.positionScreen.copy(R[Da]);N.multiplyVector4(J.positionScreen);if(!((Da+1)%sa>0)){Z=l[i-2];ca.copy(J.positionScreen);I.copy(Z.positionScreen);if(c(ca,I)){ca.multiplyScalar(1/ca.w);I.multiplyScalar(1/I.w);Ha=q[p]=q[p]||new THREE.RenderableLine;p++;s=Ha;s.v1.positionScreen.copy(ca);s.v2.positionScreen.copy(I);s.z=Math.max(ca.z,I.z);s.material=O.material;u.elements.push(s)}}}}}d=0;for(za=u.sprites.length;d<za;d++){O=u.sprites[d].object;sa=O.matrixWorld;if(O instanceof THREE.Particle){B.set(sa.elements[12],
sa.elements[13],sa.elements[14],1);K.multiplyVector4(B);B.z=B.z/B.w;if(B.z>0&&B.z<1){D=y[A]=y[A]||new THREE.RenderableParticle;A++;w=D;w.x=B.x/B.w;w.y=B.y/B.w;w.z=B.z;w.rotation=O.rotation.z;w.scale.x=O.scale.x*Math.abs(w.x-(B.x+e.projectionMatrix.elements[0])/(B.w+e.projectionMatrix.elements[12]));w.scale.y=O.scale.y*Math.abs(w.y-(B.y+e.projectionMatrix.elements[5])/(B.w+e.projectionMatrix.elements[13]));w.material=O.material;u.elements.push(w)}}}f&&u.elements.sort(b);return u}};
THREE.Quaternion=function(a,b,c,d){this.x=a||0;this.y=b||0;this.z=c||0;this.w=d!==void 0?d:1};
THREE.Quaternion.prototype={constructor:THREE.Quaternion,set:function(a,b,c,d){this.x=a;this.y=b;this.z=c;this.w=d;return this},copy:function(a){this.x=a.x;this.y=a.y;this.z=a.z;this.w=a.w;return this},setFromEuler:function(a){var b=Math.PI/360,c=a.x*b,d=a.y*b,e=a.z*b,a=Math.cos(d),d=Math.sin(d),b=Math.cos(-e),e=Math.sin(-e),f=Math.cos(c),c=Math.sin(c),h=a*b,i=d*e;this.w=h*f-i*c;this.x=h*c+i*f;this.y=d*b*f+a*e*c;this.z=a*e*f-d*b*c;return this},setFromAxisAngle:function(a,b){var c=b/2,d=Math.sin(c);
this.x=a.x*d;this.y=a.y*d;this.z=a.z*d;this.w=Math.cos(c);return this},setFromRotationMatrix:function(a){var b=Math.pow(a.determinant(),1/3);this.w=Math.sqrt(Math.max(0,b+a.elements[0]+a.elements[5]+a.elements[10]))/2;this.x=Math.sqrt(Math.max(0,b+a.elements[0]-a.elements[5]-a.elements[10]))/2;this.y=Math.sqrt(Math.max(0,b-a.elements[0]+a.elements[5]-a.elements[10]))/2;this.z=Math.sqrt(Math.max(0,b-a.elements[0]-a.elements[5]+a.elements[10]))/2;this.x=a.elements[6]-a.elements[9]<0?-Math.abs(this.x):
Math.abs(this.x);this.y=a.elements[8]-a.elements[2]<0?-Math.abs(this.y):Math.abs(this.y);this.z=a.elements[1]-a.elements[4]<0?-Math.abs(this.z):Math.abs(this.z);this.normalize();return this},calculateW:function(){this.w=-Math.sqrt(Math.abs(1-this.x*this.x-this.y*this.y-this.z*this.z));return this},inverse:function(){this.x=this.x*-1;this.y=this.y*-1;this.z=this.z*-1;return this},length:function(){return Math.sqrt(this.x*this.x+this.y*this.y+this.z*this.z+this.w*this.w)},normalize:function(){var a=
Math.sqrt(this.x*this.x+this.y*this.y+this.z*this.z+this.w*this.w);if(a===0)this.w=this.z=this.y=this.x=0;else{a=1/a;this.x=this.x*a;this.y=this.y*a;this.z=this.z*a;this.w=this.w*a}return this},multiply:function(a,b){this.x=a.x*b.w+a.y*b.z-a.z*b.y+a.w*b.x;this.y=-a.x*b.z+a.y*b.w+a.z*b.x+a.w*b.y;this.z=a.x*b.y-a.y*b.x+a.z*b.w+a.w*b.z;this.w=-a.x*b.x-a.y*b.y-a.z*b.z+a.w*b.w;return this},multiplySelf:function(a){var b=this.x,c=this.y,d=this.z,e=this.w,f=a.x,h=a.y,i=a.z,a=a.w;this.x=b*a+e*f+c*i-d*h;this.y=
c*a+e*h+d*f-b*i;this.z=d*a+e*i+b*h-c*f;this.w=e*a-b*f-c*h-d*i;return this},multiplyVector3:function(a,b){b||(b=a);var c=a.x,d=a.y,e=a.z,f=this.x,h=this.y,i=this.z,l=this.w,k=l*c+h*e-i*d,j=l*d+i*c-f*e,m=l*e+f*d-h*c,c=-f*c-h*d-i*e;b.x=k*l+c*-f+j*-i-m*-h;b.y=j*l+c*-h+m*-f-k*-i;b.z=m*l+c*-i+k*-h-j*-f;return b},clone:function(){return new THREE.Quaternion(this.x,this.y,this.z,this.w)}};
THREE.Quaternion.slerp=function(a,b,c,d){var e=a.w*b.w+a.x*b.x+a.y*b.y+a.z*b.z;if(e<0){c.w=-b.w;c.x=-b.x;c.y=-b.y;c.z=-b.z;e=-e}else c.copy(b);if(Math.abs(e)>=1){c.w=a.w;c.x=a.x;c.y=a.y;c.z=a.z;return c}var f=Math.acos(e),e=Math.sqrt(1-e*e);if(Math.abs(e)<0.001){c.w=0.5*(a.w+b.w);c.x=0.5*(a.x+b.x);c.y=0.5*(a.y+b.y);c.z=0.5*(a.z+b.z);return c}b=Math.sin((1-d)*f)/e;d=Math.sin(d*f)/e;c.w=a.w*b+c.w*d;c.x=a.x*b+c.x*d;c.y=a.y*b+c.y*d;c.z=a.z*b+c.z*d;return c};THREE.Vertex=function(){console.warn("THREE.Vertex has been DEPRECATED. Use THREE.Vector3 instead.")};
THREE.Face3=function(a,b,c,d,e,f){this.a=a;this.b=b;this.c=c;this.normal=d instanceof THREE.Vector3?d:new THREE.Vector3;this.vertexNormals=d instanceof Array?d:[];this.color=e instanceof THREE.Color?e:new THREE.Color;this.vertexColors=e instanceof Array?e:[];this.vertexTangents=[];this.materialIndex=f;this.centroid=new THREE.Vector3};
THREE.Face3.prototype={constructor:THREE.Face3,clone:function(){var a=new THREE.Face3(this.a,this.b,this.c);a.normal.copy(this.normal);a.color.copy(this.color);a.centroid.copy(this.centroid);a.materialIndex=this.materialIndex;var b,c;b=0;for(c=this.vertexNormals.length;b<c;b++)a.vertexNormals[b]=this.vertexNormals[b].clone();b=0;for(c=this.vertexColors.length;b<c;b++)a.vertexColors[b]=this.vertexColors[b].clone();b=0;for(c=this.vertexTangents.length;b<c;b++)a.vertexTangents[b]=this.vertexTangents[b].clone();
return a}};THREE.Face4=function(a,b,c,d,e,f,h){this.a=a;this.b=b;this.c=c;this.d=d;this.normal=e instanceof THREE.Vector3?e:new THREE.Vector3;this.vertexNormals=e instanceof Array?e:[];this.color=f instanceof THREE.Color?f:new THREE.Color;this.vertexColors=f instanceof Array?f:[];this.vertexTangents=[];this.materialIndex=h;this.centroid=new THREE.Vector3};
THREE.Face4.prototype={constructor:THREE.Face4,clone:function(){var a=new THREE.Face4(this.a,this.b,this.c,this.d);a.normal.copy(this.normal);a.color.copy(this.color);a.centroid.copy(this.centroid);a.materialIndex=this.materialIndex;var b,c;b=0;for(c=this.vertexNormals.length;b<c;b++)a.vertexNormals[b]=this.vertexNormals[b].clone();b=0;for(c=this.vertexColors.length;b<c;b++)a.vertexColors[b]=this.vertexColors[b].clone();b=0;for(c=this.vertexTangents.length;b<c;b++)a.vertexTangents[b]=this.vertexTangents[b].clone();
return a}};THREE.UV=function(a,b){this.u=a||0;this.v=b||0};THREE.UV.prototype={constructor:THREE.UV,set:function(a,b){this.u=a;this.v=b;return this},copy:function(a){this.u=a.u;this.v=a.v;return this},lerpSelf:function(a,b){this.u=this.u+(a.u-this.u)*b;this.v=this.v+(a.v-this.v)*b;return this},clone:function(){return new THREE.UV(this.u,this.v)}};
THREE.Geometry=function(){this.id=THREE.GeometryCount++;this.vertices=[];this.colors=[];this.materials=[];this.faces=[];this.faceUvs=[[]];this.faceVertexUvs=[[]];this.morphTargets=[];this.morphColors=[];this.morphNormals=[];this.skinWeights=[];this.skinIndices=[];this.boundingSphere=this.boundingBox=null;this.dynamic=this.hasTangents=false};
THREE.Geometry.prototype={constructor:THREE.Geometry,applyMatrix:function(a){var b=new THREE.Matrix4;b.extractRotation(a);for(var c=0,d=this.vertices.length;c<d;c++)a.multiplyVector3(this.vertices[c]);c=0;for(d=this.faces.length;c<d;c++){var e=this.faces[c];b.multiplyVector3(e.normal);for(var f=0,h=e.vertexNormals.length;f<h;f++)b.multiplyVector3(e.vertexNormals[f]);a.multiplyVector3(e.centroid)}},computeCentroids:function(){var a,b,c;a=0;for(b=this.faces.length;a<b;a++){c=this.faces[a];c.centroid.set(0,
0,0);if(c instanceof THREE.Face3){c.centroid.addSelf(this.vertices[c.a]);c.centroid.addSelf(this.vertices[c.b]);c.centroid.addSelf(this.vertices[c.c]);c.centroid.divideScalar(3)}else if(c instanceof THREE.Face4){c.centroid.addSelf(this.vertices[c.a]);c.centroid.addSelf(this.vertices[c.b]);c.centroid.addSelf(this.vertices[c.c]);c.centroid.addSelf(this.vertices[c.d]);c.centroid.divideScalar(4)}}},computeFaceNormals:function(){var a,b,c,d,e,f,h=new THREE.Vector3,i=new THREE.Vector3;a=0;for(b=this.faces.length;a<
b;a++){c=this.faces[a];d=this.vertices[c.a];e=this.vertices[c.b];f=this.vertices[c.c];h.sub(f,e);i.sub(d,e);h.crossSelf(i);h.isZero()||h.normalize();c.normal.copy(h)}},computeVertexNormals:function(){var a,b,c,d;if(this.__tmpVertices===void 0){d=this.__tmpVertices=Array(this.vertices.length);a=0;for(b=this.vertices.length;a<b;a++)d[a]=new THREE.Vector3;a=0;for(b=this.faces.length;a<b;a++){c=this.faces[a];if(c instanceof THREE.Face3)c.vertexNormals=[new THREE.Vector3,new THREE.Vector3,new THREE.Vector3];
else if(c instanceof THREE.Face4)c.vertexNormals=[new THREE.Vector3,new THREE.Vector3,new THREE.Vector3,new THREE.Vector3]}}else{d=this.__tmpVertices;a=0;for(b=this.vertices.length;a<b;a++)d[a].set(0,0,0)}a=0;for(b=this.faces.length;a<b;a++){c=this.faces[a];if(c instanceof THREE.Face3){d[c.a].addSelf(c.normal);d[c.b].addSelf(c.normal);d[c.c].addSelf(c.normal)}else if(c instanceof THREE.Face4){d[c.a].addSelf(c.normal);d[c.b].addSelf(c.normal);d[c.c].addSelf(c.normal);d[c.d].addSelf(c.normal)}}a=0;
for(b=this.vertices.length;a<b;a++)d[a].normalize();a=0;for(b=this.faces.length;a<b;a++){c=this.faces[a];if(c instanceof THREE.Face3){c.vertexNormals[0].copy(d[c.a]);c.vertexNormals[1].copy(d[c.b]);c.vertexNormals[2].copy(d[c.c])}else if(c instanceof THREE.Face4){c.vertexNormals[0].copy(d[c.a]);c.vertexNormals[1].copy(d[c.b]);c.vertexNormals[2].copy(d[c.c]);c.vertexNormals[3].copy(d[c.d])}}},computeMorphNormals:function(){var a,b,c,d,e;c=0;for(d=this.faces.length;c<d;c++){e=this.faces[c];e.__originalFaceNormal?
e.__originalFaceNormal.copy(e.normal):e.__originalFaceNormal=e.normal.clone();if(!e.__originalVertexNormals)e.__originalVertexNormals=[];a=0;for(b=e.vertexNormals.length;a<b;a++)e.__originalVertexNormals[a]?e.__originalVertexNormals[a].copy(e.vertexNormals[a]):e.__originalVertexNormals[a]=e.vertexNormals[a].clone()}var f=new THREE.Geometry;f.faces=this.faces;a=0;for(b=this.morphTargets.length;a<b;a++){if(!this.morphNormals[a]){this.morphNormals[a]={};this.morphNormals[a].faceNormals=[];this.morphNormals[a].vertexNormals=
[];var h=this.morphNormals[a].faceNormals,i=this.morphNormals[a].vertexNormals,l,k;c=0;for(d=this.faces.length;c<d;c++){e=this.faces[c];l=new THREE.Vector3;k=e instanceof THREE.Face3?{a:new THREE.Vector3,b:new THREE.Vector3,c:new THREE.Vector3}:{a:new THREE.Vector3,b:new THREE.Vector3,c:new THREE.Vector3,d:new THREE.Vector3};h.push(l);i.push(k)}}h=this.morphNormals[a];f.vertices=this.morphTargets[a].vertices;f.computeFaceNormals();f.computeVertexNormals();c=0;for(d=this.faces.length;c<d;c++){e=this.faces[c];
l=h.faceNormals[c];k=h.vertexNormals[c];l.copy(e.normal);if(e instanceof THREE.Face3){k.a.copy(e.vertexNormals[0]);k.b.copy(e.vertexNormals[1]);k.c.copy(e.vertexNormals[2])}else{k.a.copy(e.vertexNormals[0]);k.b.copy(e.vertexNormals[1]);k.c.copy(e.vertexNormals[2]);k.d.copy(e.vertexNormals[3])}}}c=0;for(d=this.faces.length;c<d;c++){e=this.faces[c];e.normal=e.__originalFaceNormal;e.vertexNormals=e.__originalVertexNormals}},computeTangents:function(){function a(a,b,c,d,g,e,f){i=a.vertices[b];l=a.vertices[c];
k=a.vertices[d];j=h[g];m=h[e];n=h[f];o=l.x-i.x;s=k.x-i.x;p=l.y-i.y;q=k.y-i.y;w=l.z-i.z;A=k.z-i.z;y=m.u-j.u;u=n.u-j.u;H=m.v-j.v;B=n.v-j.v;K=1/(y*B-u*H);I.set((B*o-H*s)*K,(B*p-H*q)*K,(B*w-H*A)*K);ba.set((y*s-u*o)*K,(y*q-u*p)*K,(y*A-u*w)*K);Y[b].addSelf(I);Y[c].addSelf(I);Y[d].addSelf(I);ca[b].addSelf(ba);ca[c].addSelf(ba);ca[d].addSelf(ba)}var b,c,d,e,f,h,i,l,k,j,m,n,o,s,p,q,w,A,y,u,H,B,K,N,Y=[],ca=[],I=new THREE.Vector3,ba=new THREE.Vector3,ja=new THREE.Vector3,ya=new THREE.Vector3,D=new THREE.Vector3;
b=0;for(c=this.vertices.length;b<c;b++){Y[b]=new THREE.Vector3;ca[b]=new THREE.Vector3}b=0;for(c=this.faces.length;b<c;b++){f=this.faces[b];h=this.faceVertexUvs[0][b];if(f instanceof THREE.Face3)a(this,f.a,f.b,f.c,0,1,2);else if(f instanceof THREE.Face4){a(this,f.a,f.b,f.d,0,1,3);a(this,f.b,f.c,f.d,1,2,3)}}var g=["a","b","c","d"];b=0;for(c=this.faces.length;b<c;b++){f=this.faces[b];for(d=0;d<f.vertexNormals.length;d++){D.copy(f.vertexNormals[d]);e=f[g[d]];N=Y[e];ja.copy(N);ja.subSelf(D.multiplyScalar(D.dot(N))).normalize();
ya.cross(f.vertexNormals[d],N);e=ya.dot(ca[e]);e=e<0?-1:1;f.vertexTangents[d]=new THREE.Vector4(ja.x,ja.y,ja.z,e)}}this.hasTangents=true},computeBoundingBox:function(){if(!this.boundingBox)this.boundingBox={min:new THREE.Vector3,max:new THREE.Vector3};if(this.vertices.length>0){var a;a=this.vertices[0];this.boundingBox.min.copy(a);this.boundingBox.max.copy(a);for(var b=this.boundingBox.min,c=this.boundingBox.max,d=1,e=this.vertices.length;d<e;d++){a=this.vertices[d];if(a.x<b.x)b.x=a.x;else if(a.x>
c.x)c.x=a.x;if(a.y<b.y)b.y=a.y;else if(a.y>c.y)c.y=a.y;if(a.z<b.z)b.z=a.z;else if(a.z>c.z)c.z=a.z}}else{this.boundingBox.min.set(0,0,0);this.boundingBox.max.set(0,0,0)}},computeBoundingSphere:function(){if(!this.boundingSphere)this.boundingSphere={radius:0};for(var a,b=0,c=0,d=this.vertices.length;c<d;c++){a=this.vertices[c].length();a>b&&(b=a)}this.boundingSphere.radius=b},mergeVertices:function(){var a={},b=[],c=[],d,e=Math.pow(10,4),f,h,i;f=0;for(h=this.vertices.length;f<h;f++){d=this.vertices[f];
d=[Math.round(d.x*e),Math.round(d.y*e),Math.round(d.z*e)].join("_");if(a[d]===void 0){a[d]=f;b.push(this.vertices[f]);c[f]=b.length-1}else c[f]=c[a[d]]}f=0;for(h=this.faces.length;f<h;f++){e=this.faces[f];if(e instanceof THREE.Face3){e.a=c[e.a];e.b=c[e.b];e.c=c[e.c]}else if(e instanceof THREE.Face4){e.a=c[e.a];e.b=c[e.b];e.c=c[e.c];e.d=c[e.d];d=[e.a,e.b,e.c,e.d];for(a=3;a>0;a--)if(d.indexOf(e["abcd"[a]])!=a){d.splice(a,1);this.faces[f]=new THREE.Face3(d[0],d[1],d[2]);e=0;for(d=this.faceVertexUvs.length;e<
d;e++)(i=this.faceVertexUvs[e][f])&&i.splice(a,1);break}}}c=this.vertices.length-b.length;this.vertices=b;return c}};THREE.GeometryCount=0;
THREE.Spline=function(a){function b(a,b,c,d,e,f,h){a=(c-a)*0.5;d=(d-b)*0.5;return(2*(b-c)+a+d)*h+(-3*(b-c)-2*a-d)*f+a*e+b}this.points=a;var c=[],d={x:0,y:0,z:0},e,f,h,i,l,k,j,m,n;this.initFromArray=function(a){this.points=[];for(var b=0;b<a.length;b++)this.points[b]={x:a[b][0],y:a[b][1],z:a[b][2]}};this.getPoint=function(a){e=(this.points.length-1)*a;f=Math.floor(e);h=e-f;c[0]=f===0?f:f-1;c[1]=f;c[2]=f>this.points.length-2?this.points.length-1:f+1;c[3]=f>this.points.length-3?this.points.length-1:
f+2;k=this.points[c[0]];j=this.points[c[1]];m=this.points[c[2]];n=this.points[c[3]];i=h*h;l=h*i;d.x=b(k.x,j.x,m.x,n.x,h,i,l);d.y=b(k.y,j.y,m.y,n.y,h,i,l);d.z=b(k.z,j.z,m.z,n.z,h,i,l);return d};this.getControlPointsArray=function(){var a,b,c=this.points.length,d=[];for(a=0;a<c;a++){b=this.points[a];d[a]=[b.x,b.y,b.z]}return d};this.getLength=function(a){var b,c,d,e=b=b=0,f=new THREE.Vector3,h=new THREE.Vector3,i=[],k=0;i[0]=0;a||(a=100);c=this.points.length*a;f.copy(this.points[0]);for(a=1;a<c;a++){b=
a/c;d=this.getPoint(b);h.copy(d);k=k+h.distanceTo(f);f.copy(d);b=(this.points.length-1)*b;b=Math.floor(b);if(b!=e){i[b]=k;e=b}}i[i.length]=k;return{chunks:i,total:k}};this.reparametrizeByArcLength=function(a){var b,c,d,e,f,h,i=[],k=new THREE.Vector3,j=this.getLength();i.push(k.copy(this.points[0]).clone());for(b=1;b<this.points.length;b++){c=j.chunks[b]-j.chunks[b-1];h=Math.ceil(a*c/j.total);e=(b-1)/(this.points.length-1);f=b/(this.points.length-1);for(c=1;c<h-1;c++){d=e+c*(1/h)*(f-e);d=this.getPoint(d);
i.push(k.copy(d).clone())}i.push(k.copy(this.points[b]).clone())}this.points=i}};THREE.Camera=function(){THREE.Object3D.call(this);this.matrixWorldInverse=new THREE.Matrix4;this.projectionMatrix=new THREE.Matrix4;this.projectionMatrixInverse=new THREE.Matrix4};THREE.Camera.prototype=new THREE.Object3D;THREE.Camera.prototype.constructor=THREE.Camera;THREE.Camera.prototype.lookAt=function(a){this.matrix.lookAt(this.position,a,this.up);this.rotationAutoUpdate&&this.rotation.getRotationFromMatrix(this.matrix)};
THREE.OrthographicCamera=function(a,b,c,d,e,f){THREE.Camera.call(this);this.left=a;this.right=b;this.top=c;this.bottom=d;this.near=e!==void 0?e:0.1;this.far=f!==void 0?f:2E3;this.updateProjectionMatrix()};THREE.OrthographicCamera.prototype=new THREE.Camera;THREE.OrthographicCamera.prototype.constructor=THREE.OrthographicCamera;THREE.OrthographicCamera.prototype.updateProjectionMatrix=function(){this.projectionMatrix.makeOrthographic(this.left,this.right,this.top,this.bottom,this.near,this.far)};
THREE.PerspectiveCamera=function(a,b,c,d){THREE.Camera.call(this);this.fov=a!==void 0?a:50;this.aspect=b!==void 0?b:1;this.near=c!==void 0?c:0.1;this.far=d!==void 0?d:2E3;this.updateProjectionMatrix()};THREE.PerspectiveCamera.prototype=new THREE.Camera;THREE.PerspectiveCamera.prototype.constructor=THREE.PerspectiveCamera;THREE.PerspectiveCamera.prototype.setLens=function(a,b){this.fov=2*Math.atan((b!==void 0?b:24)/(a*2))*(180/Math.PI);this.updateProjectionMatrix()};
THREE.PerspectiveCamera.prototype.setViewOffset=function(a,b,c,d,e,f){this.fullWidth=a;this.fullHeight=b;this.x=c;this.y=d;this.width=e;this.height=f;this.updateProjectionMatrix()};
THREE.PerspectiveCamera.prototype.updateProjectionMatrix=function(){if(this.fullWidth){var a=this.fullWidth/this.fullHeight,b=Math.tan(this.fov*Math.PI/360)*this.near,c=-b,d=a*c,a=Math.abs(a*b-d),c=Math.abs(b-c);this.projectionMatrix.makeFrustum(d+this.x*a/this.fullWidth,d+(this.x+this.width)*a/this.fullWidth,b-(this.y+this.height)*c/this.fullHeight,b-this.y*c/this.fullHeight,this.near,this.far)}else this.projectionMatrix.makePerspective(this.fov,this.aspect,this.near,this.far)};
THREE.Light=function(a){THREE.Object3D.call(this);this.color=new THREE.Color(a)};THREE.Light.prototype=new THREE.Object3D;THREE.Light.prototype.constructor=THREE.Light;THREE.Light.prototype.supr=THREE.Object3D.prototype;THREE.AmbientLight=function(a){THREE.Light.call(this,a)};THREE.AmbientLight.prototype=new THREE.Light;THREE.AmbientLight.prototype.constructor=THREE.AmbientLight;
THREE.DirectionalLight=function(a,b,c){THREE.Light.call(this,a);this.position=new THREE.Vector3(0,1,0);this.target=new THREE.Object3D;this.intensity=b!==void 0?b:1;this.distance=c!==void 0?c:0;this.onlyShadow=this.castShadow=false;this.shadowCameraNear=50;this.shadowCameraFar=5E3;this.shadowCameraLeft=-500;this.shadowCameraTop=this.shadowCameraRight=500;this.shadowCameraBottom=-500;this.shadowCameraVisible=false;this.shadowBias=0;this.shadowDarkness=0.5;this.shadowMapHeight=this.shadowMapWidth=512;
this.shadowCascade=false;this.shadowCascadeOffset=new THREE.Vector3(0,0,-1E3);this.shadowCascadeCount=2;this.shadowCascadeBias=[0,0,0];this.shadowCascadeWidth=[512,512,512];this.shadowCascadeHeight=[512,512,512];this.shadowCascadeNearZ=[-1,0.99,0.998];this.shadowCascadeFarZ=[0.99,0.998,1];this.shadowCascadeArray=[];this.shadowMatrix=this.shadowCamera=this.shadowMapSize=this.shadowMap=null};THREE.DirectionalLight.prototype=new THREE.Light;THREE.DirectionalLight.prototype.constructor=THREE.DirectionalLight;
THREE.PointLight=function(a,b,c){THREE.Light.call(this,a);this.position=new THREE.Vector3(0,0,0);this.intensity=b!==void 0?b:1;this.distance=c!==void 0?c:0};THREE.PointLight.prototype=new THREE.Light;THREE.PointLight.prototype.constructor=THREE.PointLight;
THREE.SpotLight=function(a,b,c,d,e){THREE.Light.call(this,a);this.position=new THREE.Vector3(0,1,0);this.target=new THREE.Object3D;this.intensity=b!==void 0?b:1;this.distance=c!==void 0?c:0;this.angle=d!==void 0?d:Math.PI/2;this.exponent=e!==void 0?e:10;this.onlyShadow=this.castShadow=false;this.shadowCameraNear=50;this.shadowCameraFar=5E3;this.shadowCameraFov=50;this.shadowCameraVisible=false;this.shadowBias=0;this.shadowDarkness=0.5;this.shadowMapHeight=this.shadowMapWidth=512;this.shadowMatrix=
this.shadowCamera=this.shadowMapSize=this.shadowMap=null};THREE.SpotLight.prototype=new THREE.Light;THREE.SpotLight.prototype.constructor=THREE.SpotLight;
THREE.Material=function(a){a=a||{};this.id=THREE.MaterialCount++;this.name="";this.opacity=a.opacity!==void 0?a.opacity:1;this.transparent=a.transparent!==void 0?a.transparent:false;this.blending=a.blending!==void 0?a.blending:THREE.NormalBlending;this.blendSrc=a.blendSrc!==void 0?a.blendSrc:THREE.SrcAlphaFactor;this.blendDst=a.blendDst!==void 0?a.blendDst:THREE.OneMinusSrcAlphaFactor;this.blendEquation=a.blendEquation!==void 0?a.blendEquation:THREE.AddEquation;this.depthTest=a.depthTest!==void 0?
a.depthTest:true;this.depthWrite=a.depthWrite!==void 0?a.depthWrite:true;this.polygonOffset=a.polygonOffset!==void 0?a.polygonOffset:false;this.polygonOffsetFactor=a.polygonOffsetFactor!==void 0?a.polygonOffsetFactor:0;this.polygonOffsetUnits=a.polygonOffsetUnits!==void 0?a.polygonOffsetUnits:0;this.alphaTest=a.alphaTest!==void 0?a.alphaTest:0;this.overdraw=a.overdraw!==void 0?a.overdraw:false;this.needsUpdate=this.visible=true};THREE.MaterialCount=0;THREE.NoShading=0;THREE.FlatShading=1;
THREE.SmoothShading=2;THREE.NoColors=0;THREE.FaceColors=1;THREE.VertexColors=2;THREE.NoBlending=0;THREE.NormalBlending=1;THREE.AdditiveBlending=2;THREE.SubtractiveBlending=3;THREE.MultiplyBlending=4;THREE.AdditiveAlphaBlending=5;THREE.CustomBlending=6;THREE.AddEquation=100;THREE.SubtractEquation=101;THREE.ReverseSubtractEquation=102;THREE.ZeroFactor=200;THREE.OneFactor=201;THREE.SrcColorFactor=202;THREE.OneMinusSrcColorFactor=203;THREE.SrcAlphaFactor=204;THREE.OneMinusSrcAlphaFactor=205;
THREE.DstAlphaFactor=206;THREE.OneMinusDstAlphaFactor=207;THREE.DstColorFactor=208;THREE.OneMinusDstColorFactor=209;THREE.SrcAlphaSaturateFactor=210;
THREE.LineBasicMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.color=a.color!==void 0?new THREE.Color(a.color):new THREE.Color(16777215);this.linewidth=a.linewidth!==void 0?a.linewidth:1;this.linecap=a.linecap!==void 0?a.linecap:"round";this.linejoin=a.linejoin!==void 0?a.linejoin:"round";this.vertexColors=a.vertexColors?a.vertexColors:false;this.fog=a.fog!==void 0?a.fog:true};THREE.LineBasicMaterial.prototype=new THREE.Material;THREE.LineBasicMaterial.prototype.constructor=THREE.LineBasicMaterial;
THREE.MeshBasicMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.color=a.color!==void 0?new THREE.Color(a.color):new THREE.Color(16777215);this.map=a.map!==void 0?a.map:null;this.lightMap=a.lightMap!==void 0?a.lightMap:null;this.envMap=a.envMap!==void 0?a.envMap:null;this.combine=a.combine!==void 0?a.combine:THREE.MultiplyOperation;this.reflectivity=a.reflectivity!==void 0?a.reflectivity:1;this.refractionRatio=a.refractionRatio!==void 0?a.refractionRatio:0.98;this.fog=a.fog!==void 0?a.fog:
true;this.shading=a.shading!==void 0?a.shading:THREE.SmoothShading;this.wireframe=a.wireframe!==void 0?a.wireframe:false;this.wireframeLinewidth=a.wireframeLinewidth!==void 0?a.wireframeLinewidth:1;this.wireframeLinecap=a.wireframeLinecap!==void 0?a.wireframeLinecap:"round";this.wireframeLinejoin=a.wireframeLinejoin!==void 0?a.wireframeLinejoin:"round";this.vertexColors=a.vertexColors!==void 0?a.vertexColors:THREE.NoColors;this.skinning=a.skinning!==void 0?a.skinning:false;this.morphTargets=a.morphTargets!==
void 0?a.morphTargets:false};THREE.MeshBasicMaterial.prototype=new THREE.Material;THREE.MeshBasicMaterial.prototype.constructor=THREE.MeshBasicMaterial;
THREE.MeshLambertMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.color=a.color!==void 0?new THREE.Color(a.color):new THREE.Color(16777215);this.ambient=a.ambient!==void 0?new THREE.Color(a.ambient):new THREE.Color(16777215);this.emissive=a.emissive!==void 0?new THREE.Color(a.emissive):new THREE.Color(0);this.wrapAround=a.wrapAround!==void 0?a.wrapAround:false;this.wrapRGB=new THREE.Vector3(1,1,1);this.map=a.map!==void 0?a.map:null;this.lightMap=a.lightMap!==void 0?a.lightMap:null;this.envMap=
a.envMap!==void 0?a.envMap:null;this.combine=a.combine!==void 0?a.combine:THREE.MultiplyOperation;this.reflectivity=a.reflectivity!==void 0?a.reflectivity:1;this.refractionRatio=a.refractionRatio!==void 0?a.refractionRatio:0.98;this.fog=a.fog!==void 0?a.fog:true;this.shading=a.shading!==void 0?a.shading:THREE.SmoothShading;this.wireframe=a.wireframe!==void 0?a.wireframe:false;this.wireframeLinewidth=a.wireframeLinewidth!==void 0?a.wireframeLinewidth:1;this.wireframeLinecap=a.wireframeLinecap!==void 0?
a.wireframeLinecap:"round";this.wireframeLinejoin=a.wireframeLinejoin!==void 0?a.wireframeLinejoin:"round";this.vertexColors=a.vertexColors!==void 0?a.vertexColors:THREE.NoColors;this.skinning=a.skinning!==void 0?a.skinning:false;this.morphTargets=a.morphTargets!==void 0?a.morphTargets:false;this.morphNormals=a.morphNormals!==void 0?a.morphNormals:false};THREE.MeshLambertMaterial.prototype=new THREE.Material;THREE.MeshLambertMaterial.prototype.constructor=THREE.MeshLambertMaterial;
THREE.MeshPhongMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.color=a.color!==void 0?new THREE.Color(a.color):new THREE.Color(16777215);this.ambient=a.ambient!==void 0?new THREE.Color(a.ambient):new THREE.Color(16777215);this.emissive=a.emissive!==void 0?new THREE.Color(a.emissive):new THREE.Color(0);this.specular=a.specular!==void 0?new THREE.Color(a.specular):new THREE.Color(1118481);this.shininess=a.shininess!==void 0?a.shininess:30;this.metal=a.metal!==void 0?a.metal:false;this.perPixel=
a.perPixel!==void 0?a.perPixel:false;this.wrapAround=a.wrapAround!==void 0?a.wrapAround:false;this.wrapRGB=new THREE.Vector3(1,1,1);this.map=a.map!==void 0?a.map:null;this.lightMap=a.lightMap!==void 0?a.lightMap:null;this.envMap=a.envMap!==void 0?a.envMap:null;this.combine=a.combine!==void 0?a.combine:THREE.MultiplyOperation;this.reflectivity=a.reflectivity!==void 0?a.reflectivity:1;this.refractionRatio=a.refractionRatio!==void 0?a.refractionRatio:0.98;this.fog=a.fog!==void 0?a.fog:true;this.shading=
a.shading!==void 0?a.shading:THREE.SmoothShading;this.wireframe=a.wireframe!==void 0?a.wireframe:false;this.wireframeLinewidth=a.wireframeLinewidth!==void 0?a.wireframeLinewidth:1;this.wireframeLinecap=a.wireframeLinecap!==void 0?a.wireframeLinecap:"round";this.wireframeLinejoin=a.wireframeLinejoin!==void 0?a.wireframeLinejoin:"round";this.vertexColors=a.vertexColors!==void 0?a.vertexColors:THREE.NoColors;this.skinning=a.skinning!==void 0?a.skinning:false;this.morphTargets=a.morphTargets!==void 0?
a.morphTargets:false;this.morphNormals=a.morphNormals!==void 0?a.morphNormals:false};THREE.MeshPhongMaterial.prototype=new THREE.Material;THREE.MeshPhongMaterial.prototype.constructor=THREE.MeshPhongMaterial;THREE.MeshDepthMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.shading=a.shading!==void 0?a.shading:THREE.SmoothShading;this.wireframe=a.wireframe!==void 0?a.wireframe:false;this.wireframeLinewidth=a.wireframeLinewidth!==void 0?a.wireframeLinewidth:1};
THREE.MeshDepthMaterial.prototype=new THREE.Material;THREE.MeshDepthMaterial.prototype.constructor=THREE.MeshDepthMaterial;THREE.MeshNormalMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.shading=a.shading?a.shading:THREE.FlatShading;this.wireframe=a.wireframe?a.wireframe:false;this.wireframeLinewidth=a.wireframeLinewidth?a.wireframeLinewidth:1};THREE.MeshNormalMaterial.prototype=new THREE.Material;THREE.MeshNormalMaterial.prototype.constructor=THREE.MeshNormalMaterial;
THREE.MeshFaceMaterial=function(){};THREE.ParticleBasicMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.color=a.color!==void 0?new THREE.Color(a.color):new THREE.Color(16777215);this.map=a.map!==void 0?a.map:null;this.size=a.size!==void 0?a.size:1;this.sizeAttenuation=a.sizeAttenuation!==void 0?a.sizeAttenuation:true;this.vertexColors=a.vertexColors!==void 0?a.vertexColors:false;this.fog=a.fog!==void 0?a.fog:true};THREE.ParticleBasicMaterial.prototype=new THREE.Material;
THREE.ParticleBasicMaterial.prototype.constructor=THREE.ParticleBasicMaterial;
THREE.ShaderMaterial=function(a){THREE.Material.call(this,a);a=a||{};this.fragmentShader=a.fragmentShader!==void 0?a.fragmentShader:"void main() {}";this.vertexShader=a.vertexShader!==void 0?a.vertexShader:"void main() {}";this.uniforms=a.uniforms!==void 0?a.uniforms:{};this.attributes=a.attributes;this.shading=a.shading!==void 0?a.shading:THREE.SmoothShading;this.wireframe=a.wireframe!==void 0?a.wireframe:false;this.wireframeLinewidth=a.wireframeLinewidth!==void 0?a.wireframeLinewidth:1;this.fog=
a.fog!==void 0?a.fog:false;this.lights=a.lights!==void 0?a.lights:false;this.vertexColors=a.vertexColors!==void 0?a.vertexColors:THREE.NoColors;this.skinning=a.skinning!==void 0?a.skinning:false;this.morphTargets=a.morphTargets!==void 0?a.morphTargets:false;this.morphNormals=a.morphNormals!==void 0?a.morphNormals:false};THREE.ShaderMaterial.prototype=new THREE.Material;THREE.ShaderMaterial.prototype.constructor=THREE.ShaderMaterial;
THREE.Texture=function(a,b,c,d,e,f,h,i){this.id=THREE.TextureCount++;this.image=a;this.mapping=b!==void 0?b:new THREE.UVMapping;this.wrapS=c!==void 0?c:THREE.ClampToEdgeWrapping;this.wrapT=d!==void 0?d:THREE.ClampToEdgeWrapping;this.magFilter=e!==void 0?e:THREE.LinearFilter;this.minFilter=f!==void 0?f:THREE.LinearMipMapLinearFilter;this.format=h!==void 0?h:THREE.RGBAFormat;this.type=i!==void 0?i:THREE.UnsignedByteType;this.offset=new THREE.Vector2(0,0);this.repeat=new THREE.Vector2(1,1);this.generateMipmaps=
true;this.needsUpdate=this.premultiplyAlpha=false;this.onUpdate=null};THREE.Texture.prototype={constructor:THREE.Texture,clone:function(){var a=new THREE.Texture(this.image,this.mapping,this.wrapS,this.wrapT,this.magFilter,this.minFilter,this.format,this.type);a.offset.copy(this.offset);a.repeat.copy(this.repeat);return a}};THREE.TextureCount=0;THREE.MultiplyOperation=0;THREE.MixOperation=1;THREE.UVMapping=function(){};THREE.CubeReflectionMapping=function(){};THREE.CubeRefractionMapping=function(){};
THREE.SphericalReflectionMapping=function(){};THREE.SphericalRefractionMapping=function(){};THREE.RepeatWrapping=0;THREE.ClampToEdgeWrapping=1;THREE.MirroredRepeatWrapping=2;THREE.NearestFilter=3;THREE.NearestMipMapNearestFilter=4;THREE.NearestMipMapLinearFilter=5;THREE.LinearFilter=6;THREE.LinearMipMapNearestFilter=7;THREE.LinearMipMapLinearFilter=8;THREE.ByteType=9;THREE.UnsignedByteType=10;THREE.ShortType=11;THREE.UnsignedShortType=12;THREE.IntType=13;THREE.UnsignedIntType=14;THREE.FloatType=15;
THREE.AlphaFormat=16;THREE.RGBFormat=17;THREE.RGBAFormat=18;THREE.LuminanceFormat=19;THREE.LuminanceAlphaFormat=20;THREE.DataTexture=function(a,b,c,d,e,f,h,i,l,k){THREE.Texture.call(this,null,f,h,i,l,k,d,e);this.image={data:a,width:b,height:c}};THREE.DataTexture.prototype=new THREE.Texture;THREE.DataTexture.prototype.constructor=THREE.DataTexture;
THREE.DataTexture.prototype.clone=function(){var a=new THREE.DataTexture(this.image.data,this.image.width,this.image.height,this.format,this.type,this.mapping,this.wrapS,this.wrapT,this.magFilter,this.minFilter);a.offset.copy(this.offset);a.repeat.copy(this.repeat);return a};THREE.Particle=function(a){THREE.Object3D.call(this);this.material=a};THREE.Particle.prototype=new THREE.Object3D;THREE.Particle.prototype.constructor=THREE.Particle;
THREE.ParticleSystem=function(a,b){THREE.Object3D.call(this);this.geometry=a;this.material=b!==void 0?b:new THREE.ParticleBasicMaterial({color:Math.random()*16777215});this.sortParticles=false;if(this.geometry){this.geometry.boundingSphere||this.geometry.computeBoundingSphere();this.boundRadius=a.boundingSphere.radius}this.frustumCulled=false};THREE.ParticleSystem.prototype=new THREE.Object3D;THREE.ParticleSystem.prototype.constructor=THREE.ParticleSystem;
THREE.Line=function(a,b,c){THREE.Object3D.call(this);this.geometry=a;this.material=b!==void 0?b:new THREE.LineBasicMaterial({color:Math.random()*16777215});this.type=c!==void 0?c:THREE.LineStrip;this.geometry&&(this.geometry.boundingSphere||this.geometry.computeBoundingSphere())};THREE.LineStrip=0;THREE.LinePieces=1;THREE.Line.prototype=new THREE.Object3D;THREE.Line.prototype.constructor=THREE.Line;
THREE.Mesh=function(a,b){THREE.Object3D.call(this);this.geometry=a;this.material=b!==void 0?b:new THREE.MeshBasicMaterial({color:Math.random()*16777215,wireframe:true});if(this.geometry){this.geometry.boundingSphere||this.geometry.computeBoundingSphere();this.boundRadius=a.boundingSphere.radius;if(this.geometry.morphTargets.length){this.morphTargetBase=-1;this.morphTargetForcedOrder=[];this.morphTargetInfluences=[];this.morphTargetDictionary={};for(var c=0;c<this.geometry.morphTargets.length;c++){this.morphTargetInfluences.push(0);
this.morphTargetDictionary[this.geometry.morphTargets[c].name]=c}}}};THREE.Mesh.prototype=new THREE.Object3D;THREE.Mesh.prototype.constructor=THREE.Mesh;THREE.Mesh.prototype.supr=THREE.Object3D.prototype;THREE.Mesh.prototype.getMorphTargetIndexByName=function(a){if(this.morphTargetDictionary[a]!==void 0)return this.morphTargetDictionary[a];console.log("THREE.Mesh.getMorphTargetIndexByName: morph target "+a+" does not exist. Returning 0.");return 0};
THREE.Ribbon=function(a,b){THREE.Object3D.call(this);this.geometry=a;this.material=b};THREE.Ribbon.prototype=new THREE.Object3D;THREE.Ribbon.prototype.constructor=THREE.Ribbon;THREE.LOD=function(){THREE.Object3D.call(this);this.LODs=[]};THREE.LOD.prototype=new THREE.Object3D;THREE.LOD.prototype.constructor=THREE.LOD;THREE.LOD.prototype.supr=THREE.Object3D.prototype;
THREE.LOD.prototype.addLevel=function(a,b){b===void 0&&(b=0);for(var b=Math.abs(b),c=0;c<this.LODs.length;c++)if(b<this.LODs[c].visibleAtDistance)break;this.LODs.splice(c,0,{visibleAtDistance:b,object3D:a});this.add(a)};
THREE.LOD.prototype.update=function(a){if(this.LODs.length>1){a.matrixWorldInverse.getInverse(a.matrixWorld);a=a.matrixWorldInverse;a=-(a.elements[2]*this.matrixWorld.elements[12]+a.elements[6]*this.matrixWorld.elements[13]+a.elements[10]*this.matrixWorld.elements[14]+a.elements[14]);this.LODs[0].object3D.visible=true;for(var b=1;b<this.LODs.length;b++)if(a>=this.LODs[b].visibleAtDistance){this.LODs[b-1].object3D.visible=false;this.LODs[b].object3D.visible=true}else break;for(;b<this.LODs.length;b++)this.LODs[b].object3D.visible=
false}};
THREE.Sprite=function(a){THREE.Object3D.call(this);this.color=a.color!==void 0?new THREE.Color(a.color):new THREE.Color(16777215);this.map=a.map!==void 0?a.map:new THREE.Texture;this.blending=a.blending!==void 0?a.blending:THREE.NormalBlending;this.blendSrc=a.blendSrc!==void 0?a.blendSrc:THREE.SrcAlphaFactor;this.blendDst=a.blendDst!==void 0?a.blendDst:THREE.OneMinusSrcAlphaFactor;this.blendEquation=a.blendEquation!==void 0?a.blendEquation:THREE.AddEquation;this.useScreenCoordinates=a.useScreenCoordinates!==void 0?
a.useScreenCoordinates:true;this.mergeWith3D=a.mergeWith3D!==void 0?a.mergeWith3D:!this.useScreenCoordinates;this.affectedByDistance=a.affectedByDistance!==void 0?a.affectedByDistance:!this.useScreenCoordinates;this.scaleByViewport=a.scaleByViewport!==void 0?a.scaleByViewport:!this.affectedByDistance;this.alignment=a.alignment instanceof THREE.Vector2?a.alignment:THREE.SpriteAlignment.center;this.rotation3d=this.rotation;this.rotation=0;this.opacity=1;this.uvOffset=new THREE.Vector2(0,0);this.uvScale=
new THREE.Vector2(1,1)};THREE.Sprite.prototype=new THREE.Object3D;THREE.Sprite.prototype.constructor=THREE.Sprite;THREE.Sprite.prototype.updateMatrix=function(){this.matrix.setPosition(this.position);this.rotation3d.set(0,0,this.rotation);this.matrix.setRotationFromEuler(this.rotation3d);if(this.scale.x!==1||this.scale.y!==1){this.matrix.scale(this.scale);this.boundRadiusScale=Math.max(this.scale.x,this.scale.y)}this.matrixWorldNeedsUpdate=true};THREE.SpriteAlignment={};
THREE.SpriteAlignment.topLeft=new THREE.Vector2(1,-1);THREE.SpriteAlignment.topCenter=new THREE.Vector2(0,-1);THREE.SpriteAlignment.topRight=new THREE.Vector2(-1,-1);THREE.SpriteAlignment.centerLeft=new THREE.Vector2(1,0);THREE.SpriteAlignment.center=new THREE.Vector2(0,0);THREE.SpriteAlignment.centerRight=new THREE.Vector2(-1,0);THREE.SpriteAlignment.bottomLeft=new THREE.Vector2(1,1);THREE.SpriteAlignment.bottomCenter=new THREE.Vector2(0,1);
THREE.SpriteAlignment.bottomRight=new THREE.Vector2(-1,1);THREE.Scene=function(){THREE.Object3D.call(this);this.overrideMaterial=this.fog=null;this.matrixAutoUpdate=false;this.__objects=[];this.__lights=[];this.__objectsAdded=[];this.__objectsRemoved=[]};THREE.Scene.prototype=new THREE.Object3D;THREE.Scene.prototype.constructor=THREE.Scene;
THREE.Scene.prototype.__addObject=function(a){if(a instanceof THREE.Light)this.__lights.indexOf(a)===-1&&this.__lights.push(a);else if(!(a instanceof THREE.Camera)&&this.__objects.indexOf(a)===-1){this.__objects.push(a);this.__objectsAdded.push(a);var b=this.__objectsRemoved.indexOf(a);b!==-1&&this.__objectsRemoved.splice(b,1)}for(b=0;b<a.children.length;b++)this.__addObject(a.children[b])};
THREE.Scene.prototype.__removeObject=function(a){if(a instanceof THREE.Light){var b=this.__lights.indexOf(a);b!==-1&&this.__lights.splice(b,1)}else if(!(a instanceof THREE.Camera)){b=this.__objects.indexOf(a);if(b!==-1){this.__objects.splice(b,1);this.__objectsRemoved.push(a);b=this.__objectsAdded.indexOf(a);b!==-1&&this.__objectsAdded.splice(b,1)}}for(b=0;b<a.children.length;b++)this.__removeObject(a.children[b])};
THREE.Fog=function(a,b,c){this.color=new THREE.Color(a);this.near=b!==void 0?b:1;this.far=c!==void 0?c:1E3};THREE.FogExp2=function(a,b){this.color=new THREE.Color(a);this.density=b!==void 0?b:2.5E-4};
THREE.ShaderChunk={fog_pars_fragment:"#ifdef USE_FOG\nuniform vec3 fogColor;\n#ifdef FOG_EXP2\nuniform float fogDensity;\n#else\nuniform float fogNear;\nuniform float fogFar;\n#endif\n#endif",fog_fragment:"#ifdef USE_FOG\nfloat depth = gl_FragCoord.z / gl_FragCoord.w;\n#ifdef FOG_EXP2\nconst float LOG2 = 1.442695;\nfloat fogFactor = exp2( - fogDensity * fogDensity * depth * depth * LOG2 );\nfogFactor = 1.0 - clamp( fogFactor, 0.0, 1.0 );\n#else\nfloat fogFactor = smoothstep( fogNear, fogFar, depth );\n#endif\ngl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );\n#endif",envmap_pars_fragment:"#ifdef USE_ENVMAP\nvarying vec3 vReflect;\nuniform float reflectivity;\nuniform samplerCube envMap;\nuniform float flipEnvMap;\nuniform int combine;\n#endif",
envmap_fragment:"#ifdef USE_ENVMAP\n#ifdef DOUBLE_SIDED\nfloat flipNormal = ( -1.0 + 2.0 * float( gl_FrontFacing ) );\nvec4 cubeColor = textureCube( envMap, flipNormal * vec3( flipEnvMap * vReflect.x, vReflect.yz ) );\n#else\nvec4 cubeColor = textureCube( envMap, vec3( flipEnvMap * vReflect.x, vReflect.yz ) );\n#endif\n#ifdef GAMMA_INPUT\ncubeColor.xyz *= cubeColor.xyz;\n#endif\nif ( combine == 1 ) {\ngl_FragColor.xyz = mix( gl_FragColor.xyz, cubeColor.xyz, reflectivity );\n} else {\ngl_FragColor.xyz = gl_FragColor.xyz * cubeColor.xyz;\n}\n#endif",
envmap_pars_vertex:"#ifdef USE_ENVMAP\nvarying vec3 vReflect;\nuniform float refractionRatio;\nuniform bool useRefract;\n#endif",envmap_vertex:"#ifdef USE_ENVMAP\nvec4 mPosition = objectMatrix * vec4( position, 1.0 );\nvec3 nWorld = mat3( objectMatrix[ 0 ].xyz, objectMatrix[ 1 ].xyz, objectMatrix[ 2 ].xyz ) * normal;\nif ( useRefract ) {\nvReflect = refract( normalize( mPosition.xyz - cameraPosition ), normalize( nWorld.xyz ), refractionRatio );\n} else {\nvReflect = reflect( normalize( mPosition.xyz - cameraPosition ), normalize( nWorld.xyz ) );\n}\n#endif",
map_particle_pars_fragment:"#ifdef USE_MAP\nuniform sampler2D map;\n#endif",map_particle_fragment:"#ifdef USE_MAP\ngl_FragColor = gl_FragColor * texture2D( map, gl_PointCoord );\n#endif",map_pars_vertex:"#ifdef USE_MAP\nvarying vec2 vUv;\nuniform vec4 offsetRepeat;\n#endif",map_pars_fragment:"#ifdef USE_MAP\nvarying vec2 vUv;\nuniform sampler2D map;\n#endif",map_vertex:"#ifdef USE_MAP\nvUv = uv * offsetRepeat.zw + offsetRepeat.xy;\n#endif",map_fragment:"#ifdef USE_MAP\n#ifdef GAMMA_INPUT\nvec4 texelColor = texture2D( map, vUv );\ntexelColor.xyz *= texelColor.xyz;\ngl_FragColor = gl_FragColor * texelColor;\n#else\ngl_FragColor = gl_FragColor * texture2D( map, vUv );\n#endif\n#endif",
lightmap_pars_fragment:"#ifdef USE_LIGHTMAP\nvarying vec2 vUv2;\nuniform sampler2D lightMap;\n#endif",lightmap_pars_vertex:"#ifdef USE_LIGHTMAP\nvarying vec2 vUv2;\n#endif",lightmap_fragment:"#ifdef USE_LIGHTMAP\ngl_FragColor = gl_FragColor * texture2D( lightMap, vUv2 );\n#endif",lightmap_vertex:"#ifdef USE_LIGHTMAP\nvUv2 = uv2;\n#endif",lights_lambert_pars_vertex:"uniform vec3 ambient;\nuniform vec3 diffuse;\nuniform vec3 emissive;\nuniform vec3 ambientLightColor;\n#if MAX_DIR_LIGHTS > 0\nuniform vec3 directionalLightColor[ MAX_DIR_LIGHTS ];\nuniform vec3 directionalLightDirection[ MAX_DIR_LIGHTS ];\n#endif\n#if MAX_POINT_LIGHTS > 0\nuniform vec3 pointLightColor[ MAX_POINT_LIGHTS ];\nuniform vec3 pointLightPosition[ MAX_POINT_LIGHTS ];\nuniform float pointLightDistance[ MAX_POINT_LIGHTS ];\n#endif\n#if MAX_SPOT_LIGHTS > 0\nuniform vec3 spotLightColor[ MAX_SPOT_LIGHTS ];\nuniform vec3 spotLightPosition[ MAX_SPOT_LIGHTS ];\nuniform vec3 spotLightDirection[ MAX_SPOT_LIGHTS ];\nuniform float spotLightDistance[ MAX_SPOT_LIGHTS ];\nuniform float spotLightAngle[ MAX_SPOT_LIGHTS ];\nuniform float spotLightExponent[ MAX_SPOT_LIGHTS ];\n#endif\n#ifdef WRAP_AROUND\nuniform vec3 wrapRGB;\n#endif",
lights_lambert_vertex:"vLightFront = vec3( 0.0 );\n#ifdef DOUBLE_SIDED\nvLightBack = vec3( 0.0 );\n#endif\ntransformedNormal = normalize( transformedNormal );\n#if MAX_DIR_LIGHTS > 0\nfor( int i = 0; i < MAX_DIR_LIGHTS; i ++ ) {\nvec4 lDirection = viewMatrix * vec4( directionalLightDirection[ i ], 0.0 );\nvec3 dirVector = normalize( lDirection.xyz );\nfloat dotProduct = dot( transformedNormal, dirVector );\nvec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );\n#ifdef DOUBLE_SIDED\nvec3 directionalLightWeightingBack = vec3( max( -dotProduct, 0.0 ) );\n#ifdef WRAP_AROUND\nvec3 directionalLightWeightingHalfBack = vec3( max( -0.5 * dotProduct + 0.5, 0.0 ) );\n#endif\n#endif\n#ifdef WRAP_AROUND\nvec3 directionalLightWeightingHalf = vec3( max( 0.5 * dotProduct + 0.5, 0.0 ) );\ndirectionalLightWeighting = mix( directionalLightWeighting, directionalLightWeightingHalf, wrapRGB );\n#ifdef DOUBLE_SIDED\ndirectionalLightWeightingBack = mix( directionalLightWeightingBack, directionalLightWeightingHalfBack, wrapRGB );\n#endif\n#endif\nvLightFront += directionalLightColor[ i ] * directionalLightWeighting;\n#ifdef DOUBLE_SIDED\nvLightBack += directionalLightColor[ i ] * directionalLightWeightingBack;\n#endif\n}\n#endif\n#if MAX_POINT_LIGHTS > 0\nfor( int i = 0; i < MAX_POINT_LIGHTS; i ++ ) {\nvec4 lPosition = viewMatrix * vec4( pointLightPosition[ i ], 1.0 );\nvec3 lVector = lPosition.xyz - mvPosition.xyz;\nfloat lDistance = 1.0;\nif ( pointLightDistance[ i ] > 0.0 )\nlDistance = 1.0 - min( ( length( lVector ) / pointLightDistance[ i ] ), 1.0 );\nlVector = normalize( lVector );\nfloat dotProduct = dot( transformedNormal, lVector );\nvec3 pointLightWeighting = vec3( max( dotProduct, 0.0 ) );\n#ifdef DOUBLE_SIDED\nvec3 pointLightWeightingBack = vec3( max( -dotProduct, 0.0 ) );\n#ifdef WRAP_AROUND\nvec3 pointLightWeightingHalfBack = vec3( max( -0.5 * dotProduct + 0.5, 0.0 ) );\n#endif\n#endif\n#ifdef WRAP_AROUND\nvec3 pointLightWeightingHalf = vec3( max( 0.5 * dotProduct + 0.5, 0.0 ) );\npointLightWeighting = mix( pointLightWeighting, pointLightWeightingHalf, wrapRGB );\n#ifdef DOUBLE_SIDED\npointLightWeightingBack = mix( pointLightWeightingBack, pointLightWeightingHalfBack, wrapRGB );\n#endif\n#endif\nvLightFront += pointLightColor[ i ] * pointLightWeighting * lDistance;\n#ifdef DOUBLE_SIDED\nvLightBack += pointLightColor[ i ] * pointLightWeightingBack * lDistance;\n#endif\n}\n#endif\n#if MAX_SPOT_LIGHTS > 0\nfor( int i = 0; i < MAX_SPOT_LIGHTS; i ++ ) {\nvec4 lPosition = viewMatrix * vec4( spotLightPosition[ i ], 1.0 );\nvec3 lVector = lPosition.xyz - mvPosition.xyz;\nlVector = normalize( lVector );\nfloat spotEffect = dot( spotLightDirection[ i ], normalize( spotLightPosition[ i ] - mPosition.xyz ) );\nif ( spotEffect > spotLightAngle[ i ] ) {\nspotEffect = pow( spotEffect, spotLightExponent[ i ] );\nfloat lDistance = 1.0;\nif ( spotLightDistance[ i ] > 0.0 )\nlDistance = 1.0 - min( ( length( lVector ) / spotLightDistance[ i ] ), 1.0 );\nfloat dotProduct = dot( transformedNormal, lVector );\nvec3 spotLightWeighting = vec3( max( dotProduct, 0.0 ) );\n#ifdef DOUBLE_SIDED\nvec3 spotLightWeightingBack = vec3( max( -dotProduct, 0.0 ) );\n#ifdef WRAP_AROUND\nvec3 spotLightWeightingHalfBack = vec3( max( -0.5 * dotProduct + 0.5, 0.0 ) );\n#endif\n#endif\n#ifdef WRAP_AROUND\nvec3 spotLightWeightingHalf = vec3( max( 0.5 * dotProduct + 0.5, 0.0 ) );\nspotLightWeighting = mix( spotLightWeighting, spotLightWeightingHalf, wrapRGB );\n#ifdef DOUBLE_SIDED\nspotLightWeightingBack = mix( spotLightWeightingBack, spotLightWeightingHalfBack, wrapRGB );\n#endif\n#endif\nvLightFront += spotLightColor[ i ] * spotLightWeighting * lDistance * spotEffect;\n#ifdef DOUBLE_SIDED\nvLightBack += spotLightColor[ i ] * spotLightWeightingBack * lDistance * spotEffect;\n#endif\n}\n}\n#endif\nvLightFront = vLightFront * diffuse + ambient * ambientLightColor + emissive;\n#ifdef DOUBLE_SIDED\nvLightBack = vLightBack * diffuse + ambient * ambientLightColor + emissive;\n#endif",
lights_phong_pars_vertex:"#ifndef PHONG_PER_PIXEL\n#if MAX_POINT_LIGHTS > 0\nuniform vec3 pointLightPosition[ MAX_POINT_LIGHTS ];\nuniform float pointLightDistance[ MAX_POINT_LIGHTS ];\nvarying vec4 vPointLight[ MAX_POINT_LIGHTS ];\n#endif\n#if MAX_SPOT_LIGHTS > 0\nuniform vec3 spotLightPosition[ MAX_SPOT_LIGHTS ];\nuniform float spotLightDistance[ MAX_SPOT_LIGHTS ];\nvarying vec4 vSpotLight[ MAX_SPOT_LIGHTS ];\n#endif\n#endif\n#if MAX_SPOT_LIGHTS > 0\nvarying vec3 vWorldPosition;\n#endif",lights_phong_vertex:"#ifndef PHONG_PER_PIXEL\n#if MAX_POINT_LIGHTS > 0\nfor( int i = 0; i < MAX_POINT_LIGHTS; i ++ ) {\nvec4 lPosition = viewMatrix * vec4( pointLightPosition[ i ], 1.0 );\nvec3 lVector = lPosition.xyz - mvPosition.xyz;\nfloat lDistance = 1.0;\nif ( pointLightDistance[ i ] > 0.0 )\nlDistance = 1.0 - min( ( length( lVector ) / pointLightDistance[ i ] ), 1.0 );\nvPointLight[ i ] = vec4( lVector, lDistance );\n}\n#endif\n#if MAX_SPOT_LIGHTS > 0\nfor( int i = 0; i < MAX_SPOT_LIGHTS; i ++ ) {\nvec4 lPosition = viewMatrix * vec4( spotLightPosition[ i ], 1.0 );\nvec3 lVector = lPosition.xyz - mvPosition.xyz;\nfloat lDistance = 1.0;\nif ( spotLightDistance[ i ] > 0.0 )\nlDistance = 1.0 - min( ( length( lVector ) / spotLightDistance[ i ] ), 1.0 );\nvSpotLight[ i ] = vec4( lVector, lDistance );\n}\n#endif\n#endif\n#if MAX_SPOT_LIGHTS > 0\nvWorldPosition = mPosition.xyz;\n#endif",
lights_phong_pars_fragment:"uniform vec3 ambientLightColor;\n#if MAX_DIR_LIGHTS > 0\nuniform vec3 directionalLightColor[ MAX_DIR_LIGHTS ];\nuniform vec3 directionalLightDirection[ MAX_DIR_LIGHTS ];\n#endif\n#if MAX_POINT_LIGHTS > 0\nuniform vec3 pointLightColor[ MAX_POINT_LIGHTS ];\n#ifdef PHONG_PER_PIXEL\nuniform vec3 pointLightPosition[ MAX_POINT_LIGHTS ];\nuniform float pointLightDistance[ MAX_POINT_LIGHTS ];\n#else\nvarying vec4 vPointLight[ MAX_POINT_LIGHTS ];\n#endif\n#endif\n#if MAX_SPOT_LIGHTS > 0\nuniform vec3 spotLightColor[ MAX_SPOT_LIGHTS ];\nuniform vec3 spotLightPosition[ MAX_SPOT_LIGHTS ];\nuniform vec3 spotLightDirection[ MAX_SPOT_LIGHTS ];\nuniform float spotLightAngle[ MAX_SPOT_LIGHTS ];\nuniform float spotLightExponent[ MAX_SPOT_LIGHTS ];\n#ifdef PHONG_PER_PIXEL\nuniform float spotLightDistance[ MAX_SPOT_LIGHTS ];\n#else\nvarying vec4 vSpotLight[ MAX_SPOT_LIGHTS ];\n#endif\nvarying vec3 vWorldPosition;\n#endif\n#ifdef WRAP_AROUND\nuniform vec3 wrapRGB;\n#endif\nvarying vec3 vViewPosition;\nvarying vec3 vNormal;",
lights_phong_fragment:"vec3 normal = normalize( vNormal );\nvec3 viewPosition = normalize( vViewPosition );\n#ifdef DOUBLE_SIDED\nnormal = normal * ( -1.0 + 2.0 * float( gl_FrontFacing ) );\n#endif\n#if MAX_POINT_LIGHTS > 0\nvec3 pointDiffuse  = vec3( 0.0 );\nvec3 pointSpecular = vec3( 0.0 );\nfor ( int i = 0; i < MAX_POINT_LIGHTS; i ++ ) {\n#ifdef PHONG_PER_PIXEL\nvec4 lPosition = viewMatrix * vec4( pointLightPosition[ i ], 1.0 );\nvec3 lVector = lPosition.xyz + vViewPosition.xyz;\nfloat lDistance = 1.0;\nif ( pointLightDistance[ i ] > 0.0 )\nlDistance = 1.0 - min( ( length( lVector ) / pointLightDistance[ i ] ), 1.0 );\nlVector = normalize( lVector );\n#else\nvec3 lVector = normalize( vPointLight[ i ].xyz );\nfloat lDistance = vPointLight[ i ].w;\n#endif\nfloat dotProduct = dot( normal, lVector );\n#ifdef WRAP_AROUND\nfloat pointDiffuseWeightFull = max( dotProduct, 0.0 );\nfloat pointDiffuseWeightHalf = max( 0.5 * dotProduct + 0.5, 0.0 );\nvec3 pointDiffuseWeight = mix( vec3 ( pointDiffuseWeightFull ), vec3( pointDiffuseWeightHalf ), wrapRGB );\n#else\nfloat pointDiffuseWeight = max( dotProduct, 0.0 );\n#endif\npointDiffuse  += diffuse * pointLightColor[ i ] * pointDiffuseWeight * lDistance;\nvec3 pointHalfVector = normalize( lVector + viewPosition );\nfloat pointDotNormalHalf = max( dot( normal, pointHalfVector ), 0.0 );\nfloat pointSpecularWeight = max( pow( pointDotNormalHalf, shininess ), 0.0 );\n#ifdef PHYSICALLY_BASED_SHADING\nfloat specularNormalization = ( shininess + 2.0001 ) / 8.0;\nvec3 schlick = specular + vec3( 1.0 - specular ) * pow( 1.0 - dot( lVector, pointHalfVector ), 5.0 );\npointSpecular += schlick * pointLightColor[ i ] * pointSpecularWeight * pointDiffuseWeight * lDistance * specularNormalization;\n#else\npointSpecular += specular * pointLightColor[ i ] * pointSpecularWeight * pointDiffuseWeight * lDistance;\n#endif\n}\n#endif\n#if MAX_SPOT_LIGHTS > 0\nvec3 spotDiffuse  = vec3( 0.0 );\nvec3 spotSpecular = vec3( 0.0 );\nfor ( int i = 0; i < MAX_SPOT_LIGHTS; i ++ ) {\n#ifdef PHONG_PER_PIXEL\nvec4 lPosition = viewMatrix * vec4( spotLightPosition[ i ], 1.0 );\nvec3 lVector = lPosition.xyz + vViewPosition.xyz;\nfloat lDistance = 1.0;\nif ( spotLightDistance[ i ] > 0.0 )\nlDistance = 1.0 - min( ( length( lVector ) / spotLightDistance[ i ] ), 1.0 );\nlVector = normalize( lVector );\n#else\nvec3 lVector = normalize( vSpotLight[ i ].xyz );\nfloat lDistance = vSpotLight[ i ].w;\n#endif\nfloat spotEffect = dot( spotLightDirection[ i ], normalize( spotLightPosition[ i ] - vWorldPosition ) );\nif ( spotEffect > spotLightAngle[ i ] ) {\nspotEffect = pow( spotEffect, spotLightExponent[ i ] );\nfloat dotProduct = dot( normal, lVector );\n#ifdef WRAP_AROUND\nfloat spotDiffuseWeightFull = max( dotProduct, 0.0 );\nfloat spotDiffuseWeightHalf = max( 0.5 * dotProduct + 0.5, 0.0 );\nvec3 spotDiffuseWeight = mix( vec3 ( spotDiffuseWeightFull ), vec3( spotDiffuseWeightHalf ), wrapRGB );\n#else\nfloat spotDiffuseWeight = max( dotProduct, 0.0 );\n#endif\nspotDiffuse += diffuse * spotLightColor[ i ] * spotDiffuseWeight * lDistance * spotEffect;\nvec3 spotHalfVector = normalize( lVector + viewPosition );\nfloat spotDotNormalHalf = max( dot( normal, spotHalfVector ), 0.0 );\nfloat spotSpecularWeight = max( pow( spotDotNormalHalf, shininess ), 0.0 );\n#ifdef PHYSICALLY_BASED_SHADING\nfloat specularNormalization = ( shininess + 2.0001 ) / 8.0;\nvec3 schlick = specular + vec3( 1.0 - specular ) * pow( 1.0 - dot( lVector, spotHalfVector ), 5.0 );\nspotSpecular += schlick * spotLightColor[ i ] * spotSpecularWeight * spotDiffuseWeight * lDistance * specularNormalization * spotEffect;\n#else\nspotSpecular += specular * spotLightColor[ i ] * spotSpecularWeight * spotDiffuseWeight * lDistance * spotEffect;\n#endif\n}\n}\n#endif\n#if MAX_DIR_LIGHTS > 0\nvec3 dirDiffuse  = vec3( 0.0 );\nvec3 dirSpecular = vec3( 0.0 );\nfor( int i = 0; i < MAX_DIR_LIGHTS; i ++ ) {\nvec4 lDirection = viewMatrix * vec4( directionalLightDirection[ i ], 0.0 );\nvec3 dirVector = normalize( lDirection.xyz );\nfloat dotProduct = dot( normal, dirVector );\n#ifdef WRAP_AROUND\nfloat dirDiffuseWeightFull = max( dotProduct, 0.0 );\nfloat dirDiffuseWeightHalf = max( 0.5 * dotProduct + 0.5, 0.0 );\nvec3 dirDiffuseWeight = mix( vec3( dirDiffuseWeightFull ), vec3( dirDiffuseWeightHalf ), wrapRGB );\n#else\nfloat dirDiffuseWeight = max( dotProduct, 0.0 );\n#endif\ndirDiffuse  += diffuse * directionalLightColor[ i ] * dirDiffuseWeight;\nvec3 dirHalfVector = normalize( dirVector + viewPosition );\nfloat dirDotNormalHalf = max( dot( normal, dirHalfVector ), 0.0 );\nfloat dirSpecularWeight = max( pow( dirDotNormalHalf, shininess ), 0.0 );\n#ifdef PHYSICALLY_BASED_SHADING\nfloat specularNormalization = ( shininess + 2.0001 ) / 8.0;\nvec3 schlick = specular + vec3( 1.0 - specular ) * pow( 1.0 - dot( dirVector, dirHalfVector ), 5.0 );\ndirSpecular += schlick * directionalLightColor[ i ] * dirSpecularWeight * dirDiffuseWeight * specularNormalization;\n#else\ndirSpecular += specular * directionalLightColor[ i ] * dirSpecularWeight * dirDiffuseWeight;\n#endif\n}\n#endif\nvec3 totalDiffuse = vec3( 0.0 );\nvec3 totalSpecular = vec3( 0.0 );\n#if MAX_DIR_LIGHTS > 0\ntotalDiffuse += dirDiffuse;\ntotalSpecular += dirSpecular;\n#endif\n#if MAX_POINT_LIGHTS > 0\ntotalDiffuse += pointDiffuse;\ntotalSpecular += pointSpecular;\n#endif\n#if MAX_SPOT_LIGHTS > 0\ntotalDiffuse += spotDiffuse;\ntotalSpecular += spotSpecular;\n#endif\n#ifdef METAL\ngl_FragColor.xyz = gl_FragColor.xyz * ( emissive + totalDiffuse + ambientLightColor * ambient + totalSpecular );\n#else\ngl_FragColor.xyz = gl_FragColor.xyz * ( emissive + totalDiffuse + ambientLightColor * ambient ) + totalSpecular;\n#endif",
color_pars_fragment:"#ifdef USE_COLOR\nvarying vec3 vColor;\n#endif",color_fragment:"#ifdef USE_COLOR\ngl_FragColor = gl_FragColor * vec4( vColor, opacity );\n#endif",color_pars_vertex:"#ifdef USE_COLOR\nvarying vec3 vColor;\n#endif",color_vertex:"#ifdef USE_COLOR\n#ifdef GAMMA_INPUT\nvColor = color * color;\n#else\nvColor = color;\n#endif\n#endif",skinning_pars_vertex:"#ifdef USE_SKINNING\nuniform mat4 boneGlobalMatrices[ MAX_BONES ];\n#endif",skinning_vertex:"#ifdef USE_SKINNING\ngl_Position  = ( boneGlobalMatrices[ int( skinIndex.x ) ] * skinVertexA ) * skinWeight.x;\ngl_Position += ( boneGlobalMatrices[ int( skinIndex.y ) ] * skinVertexB ) * skinWeight.y;\ngl_Position  = projectionMatrix * modelViewMatrix * gl_Position;\n#endif",
morphtarget_pars_vertex:"#ifdef USE_MORPHTARGETS\n#ifndef USE_MORPHNORMALS\nuniform float morphTargetInfluences[ 8 ];\n#else\nuniform float morphTargetInfluences[ 4 ];\n#endif\n#endif",morphtarget_vertex:"#ifdef USE_MORPHTARGETS\nvec3 morphed = vec3( 0.0 );\nmorphed += ( morphTarget0 - position ) * morphTargetInfluences[ 0 ];\nmorphed += ( morphTarget1 - position ) * morphTargetInfluences[ 1 ];\nmorphed += ( morphTarget2 - position ) * morphTargetInfluences[ 2 ];\nmorphed += ( morphTarget3 - position ) * morphTargetInfluences[ 3 ];\n#ifndef USE_MORPHNORMALS\nmorphed += ( morphTarget4 - position ) * morphTargetInfluences[ 4 ];\nmorphed += ( morphTarget5 - position ) * morphTargetInfluences[ 5 ];\nmorphed += ( morphTarget6 - position ) * morphTargetInfluences[ 6 ];\nmorphed += ( morphTarget7 - position ) * morphTargetInfluences[ 7 ];\n#endif\nmorphed += position;\ngl_Position = projectionMatrix * modelViewMatrix * vec4( morphed, 1.0 );\n#endif",
default_vertex:"#ifndef USE_MORPHTARGETS\n#ifndef USE_SKINNING\ngl_Position = projectionMatrix * mvPosition;\n#endif\n#endif",morphnormal_vertex:"#ifdef USE_MORPHNORMALS\nvec3 morphedNormal = vec3( 0.0 );\nmorphedNormal +=  ( morphNormal0 - normal ) * morphTargetInfluences[ 0 ];\nmorphedNormal +=  ( morphNormal1 - normal ) * morphTargetInfluences[ 1 ];\nmorphedNormal +=  ( morphNormal2 - normal ) * morphTargetInfluences[ 2 ];\nmorphedNormal +=  ( morphNormal3 - normal ) * morphTargetInfluences[ 3 ];\nmorphedNormal += normal;\nvec3 transformedNormal = normalMatrix * morphedNormal;\n#else\nvec3 transformedNormal = normalMatrix * normal;\n#endif",
shadowmap_pars_fragment:"#ifdef USE_SHADOWMAP\nuniform sampler2D shadowMap[ MAX_SHADOWS ];\nuniform vec2 shadowMapSize[ MAX_SHADOWS ];\nuniform float shadowDarkness[ MAX_SHADOWS ];\nuniform float shadowBias[ MAX_SHADOWS ];\nvarying vec4 vShadowCoord[ MAX_SHADOWS ];\nfloat unpackDepth( const in vec4 rgba_depth ) {\nconst vec4 bit_shift = vec4( 1.0 / ( 256.0 * 256.0 * 256.0 ), 1.0 / ( 256.0 * 256.0 ), 1.0 / 256.0, 1.0 );\nfloat depth = dot( rgba_depth, bit_shift );\nreturn depth;\n}\n#endif",shadowmap_fragment:"#ifdef USE_SHADOWMAP\n#ifdef SHADOWMAP_DEBUG\nvec3 frustumColors[3];\nfrustumColors[0] = vec3( 1.0, 0.5, 0.0 );\nfrustumColors[1] = vec3( 0.0, 1.0, 0.8 );\nfrustumColors[2] = vec3( 0.0, 0.5, 1.0 );\n#endif\n#ifdef SHADOWMAP_CASCADE\nint inFrustumCount = 0;\n#endif\nfloat fDepth;\nvec3 shadowColor = vec3( 1.0 );\nfor( int i = 0; i < MAX_SHADOWS; i ++ ) {\nvec3 shadowCoord = vShadowCoord[ i ].xyz / vShadowCoord[ i ].w;\nbvec4 inFrustumVec = bvec4 ( shadowCoord.x >= 0.0, shadowCoord.x <= 1.0, shadowCoord.y >= 0.0, shadowCoord.y <= 1.0 );\nbool inFrustum = all( inFrustumVec );\n#ifdef SHADOWMAP_CASCADE\ninFrustumCount += int( inFrustum );\nbvec3 frustumTestVec = bvec3( inFrustum, inFrustumCount == 1, shadowCoord.z <= 1.0 );\n#else\nbvec2 frustumTestVec = bvec2( inFrustum, shadowCoord.z <= 1.0 );\n#endif\nbool frustumTest = all( frustumTestVec );\nif ( frustumTest ) {\nshadowCoord.z += shadowBias[ i ];\n#ifdef SHADOWMAP_SOFT\nfloat shadow = 0.0;\nconst float shadowDelta = 1.0 / 9.0;\nfloat xPixelOffset = 1.0 / shadowMapSize[ i ].x;\nfloat yPixelOffset = 1.0 / shadowMapSize[ i ].y;\nfloat dx0 = -1.25 * xPixelOffset;\nfloat dy0 = -1.25 * yPixelOffset;\nfloat dx1 = 1.25 * xPixelOffset;\nfloat dy1 = 1.25 * yPixelOffset;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( dx0, dy0 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( 0.0, dy0 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( dx1, dy0 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( dx0, 0.0 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( dx1, 0.0 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( dx0, dy1 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( 0.0, dy1 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nfDepth = unpackDepth( texture2D( shadowMap[ i ], shadowCoord.xy + vec2( dx1, dy1 ) ) );\nif ( fDepth < shadowCoord.z ) shadow += shadowDelta;\nshadowColor = shadowColor * vec3( ( 1.0 - shadowDarkness[ i ] * shadow ) );\n#else\nvec4 rgbaDepth = texture2D( shadowMap[ i ], shadowCoord.xy );\nfloat fDepth = unpackDepth( rgbaDepth );\nif ( fDepth < shadowCoord.z )\nshadowColor = shadowColor * vec3( 1.0 - shadowDarkness[ i ] );\n#endif\n}\n#ifdef SHADOWMAP_DEBUG\n#ifdef SHADOWMAP_CASCADE\nif ( inFrustum && inFrustumCount == 1 ) gl_FragColor.xyz *= frustumColors[ i ];\n#else\nif ( inFrustum ) gl_FragColor.xyz *= frustumColors[ i ];\n#endif\n#endif\n}\n#ifdef GAMMA_OUTPUT\nshadowColor *= shadowColor;\n#endif\ngl_FragColor.xyz = gl_FragColor.xyz * shadowColor;\n#endif",
shadowmap_pars_vertex:"#ifdef USE_SHADOWMAP\nvarying vec4 vShadowCoord[ MAX_SHADOWS ];\nuniform mat4 shadowMatrix[ MAX_SHADOWS ];\n#endif",shadowmap_vertex:"#ifdef USE_SHADOWMAP\nfor( int i = 0; i < MAX_SHADOWS; i ++ ) {\n#ifdef USE_MORPHTARGETS\nvShadowCoord[ i ] = shadowMatrix[ i ] * objectMatrix * vec4( morphed, 1.0 );\n#else\nvShadowCoord[ i ] = shadowMatrix[ i ] * objectMatrix * vec4( position, 1.0 );\n#endif\n}\n#endif",alphatest_fragment:"#ifdef ALPHATEST\nif ( gl_FragColor.a < ALPHATEST ) discard;\n#endif",
linear_to_gamma_fragment:"#ifdef GAMMA_OUTPUT\ngl_FragColor.xyz = sqrt( gl_FragColor.xyz );\n#endif"};
THREE.UniformsUtils={merge:function(a){var b,c,d,e={};for(b=0;b<a.length;b++){d=this.clone(a[b]);for(c in d)e[c]=d[c]}return e},clone:function(a){var b,c,d,e={};for(b in a){e[b]={};for(c in a[b]){d=a[b][c];e[b][c]=d instanceof THREE.Color||d instanceof THREE.Vector2||d instanceof THREE.Vector3||d instanceof THREE.Vector4||d instanceof THREE.Matrix4||d instanceof THREE.Texture?d.clone():d instanceof Array?d.slice():d}}return e}};
THREE.UniformsLib={common:{diffuse:{type:"c",value:new THREE.Color(15658734)},opacity:{type:"f",value:1},map:{type:"t",value:0,texture:null},offsetRepeat:{type:"v4",value:new THREE.Vector4(0,0,1,1)},lightMap:{type:"t",value:2,texture:null},envMap:{type:"t",value:1,texture:null},flipEnvMap:{type:"f",value:-1},useRefract:{type:"i",value:0},reflectivity:{type:"f",value:1},refractionRatio:{type:"f",value:0.98},combine:{type:"i",value:0},morphTargetInfluences:{type:"f",value:0}},fog:{fogDensity:{type:"f",
value:2.5E-4},fogNear:{type:"f",value:1},fogFar:{type:"f",value:2E3},fogColor:{type:"c",value:new THREE.Color(16777215)}},lights:{ambientLightColor:{type:"fv",value:[]},directionalLightDirection:{type:"fv",value:[]},directionalLightColor:{type:"fv",value:[]},pointLightColor:{type:"fv",value:[]},pointLightPosition:{type:"fv",value:[]},pointLightDistance:{type:"fv1",value:[]},spotLightColor:{type:"fv",value:[]},spotLightPosition:{type:"fv",value:[]},spotLightDirection:{type:"fv",value:[]},spotLightDistance:{type:"fv1",
value:[]},spotLightAngle:{type:"fv1",value:[]},spotLightExponent:{type:"fv1",value:[]}},particle:{psColor:{type:"c",value:new THREE.Color(15658734)},opacity:{type:"f",value:1},size:{type:"f",value:1},scale:{type:"f",value:1},map:{type:"t",value:0,texture:null},fogDensity:{type:"f",value:2.5E-4},fogNear:{type:"f",value:1},fogFar:{type:"f",value:2E3},fogColor:{type:"c",value:new THREE.Color(16777215)}},shadowmap:{shadowMap:{type:"tv",value:6,texture:[]},shadowMapSize:{type:"v2v",value:[]},shadowBias:{type:"fv1",
value:[]},shadowDarkness:{type:"fv1",value:[]},shadowMatrix:{type:"m4v",value:[]}}};
THREE.ShaderLib={depth:{uniforms:{mNear:{type:"f",value:1},mFar:{type:"f",value:2E3},opacity:{type:"f",value:1}},vertexShader:"void main() {\ngl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );\n}",fragmentShader:"uniform float mNear;\nuniform float mFar;\nuniform float opacity;\nvoid main() {\nfloat depth = gl_FragCoord.z / gl_FragCoord.w;\nfloat color = 1.0 - smoothstep( mNear, mFar, depth );\ngl_FragColor = vec4( vec3( color ), opacity );\n}"},normal:{uniforms:{opacity:{type:"f",
value:1}},vertexShader:"varying vec3 vNormal;\nvoid main() {\nvec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );\nvNormal = normalMatrix * normal;\ngl_Position = projectionMatrix * mvPosition;\n}",fragmentShader:"uniform float opacity;\nvarying vec3 vNormal;\nvoid main() {\ngl_FragColor = vec4( 0.5 * normalize( vNormal ) + 0.5, opacity );\n}"},basic:{uniforms:THREE.UniformsUtils.merge([THREE.UniformsLib.common,THREE.UniformsLib.fog,THREE.UniformsLib.shadowmap]),vertexShader:[THREE.ShaderChunk.map_pars_vertex,
THREE.ShaderChunk.lightmap_pars_vertex,THREE.ShaderChunk.envmap_pars_vertex,THREE.ShaderChunk.color_pars_vertex,THREE.ShaderChunk.skinning_pars_vertex,THREE.ShaderChunk.morphtarget_pars_vertex,THREE.ShaderChunk.shadowmap_pars_vertex,"void main() {\nvec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",THREE.ShaderChunk.map_vertex,THREE.ShaderChunk.lightmap_vertex,THREE.ShaderChunk.envmap_vertex,THREE.ShaderChunk.color_vertex,THREE.ShaderChunk.skinning_vertex,THREE.ShaderChunk.morphtarget_vertex,
THREE.ShaderChunk.default_vertex,THREE.ShaderChunk.shadowmap_vertex,"}"].join("\n"),fragmentShader:["uniform vec3 diffuse;\nuniform float opacity;",THREE.ShaderChunk.color_pars_fragment,THREE.ShaderChunk.map_pars_fragment,THREE.ShaderChunk.lightmap_pars_fragment,THREE.ShaderChunk.envmap_pars_fragment,THREE.ShaderChunk.fog_pars_fragment,THREE.ShaderChunk.shadowmap_pars_fragment,"void main() {\ngl_FragColor = vec4( diffuse, opacity );",THREE.ShaderChunk.map_fragment,THREE.ShaderChunk.alphatest_fragment,
THREE.ShaderChunk.lightmap_fragment,THREE.ShaderChunk.color_fragment,THREE.ShaderChunk.envmap_fragment,THREE.ShaderChunk.shadowmap_fragment,THREE.ShaderChunk.linear_to_gamma_fragment,THREE.ShaderChunk.fog_fragment,"}"].join("\n")},lambert:{uniforms:THREE.UniformsUtils.merge([THREE.UniformsLib.common,THREE.UniformsLib.fog,THREE.UniformsLib.lights,THREE.UniformsLib.shadowmap,{ambient:{type:"c",value:new THREE.Color(16777215)},emissive:{type:"c",value:new THREE.Color(0)},wrapRGB:{type:"v3",value:new THREE.Vector3(1,
1,1)}}]),vertexShader:["varying vec3 vLightFront;\n#ifdef DOUBLE_SIDED\nvarying vec3 vLightBack;\n#endif",THREE.ShaderChunk.map_pars_vertex,THREE.ShaderChunk.lightmap_pars_vertex,THREE.ShaderChunk.envmap_pars_vertex,THREE.ShaderChunk.lights_lambert_pars_vertex,THREE.ShaderChunk.color_pars_vertex,THREE.ShaderChunk.skinning_pars_vertex,THREE.ShaderChunk.morphtarget_pars_vertex,THREE.ShaderChunk.shadowmap_pars_vertex,"void main() {\nvec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",THREE.ShaderChunk.map_vertex,
THREE.ShaderChunk.lightmap_vertex,THREE.ShaderChunk.envmap_vertex,THREE.ShaderChunk.color_vertex,THREE.ShaderChunk.morphnormal_vertex,"#ifndef USE_ENVMAP\nvec4 mPosition = objectMatrix * vec4( position, 1.0 );\n#endif",THREE.ShaderChunk.lights_lambert_vertex,THREE.ShaderChunk.skinning_vertex,THREE.ShaderChunk.morphtarget_vertex,THREE.ShaderChunk.default_vertex,THREE.ShaderChunk.shadowmap_vertex,"}"].join("\n"),fragmentShader:["uniform float opacity;\nvarying vec3 vLightFront;\n#ifdef DOUBLE_SIDED\nvarying vec3 vLightBack;\n#endif",
THREE.ShaderChunk.color_pars_fragment,THREE.ShaderChunk.map_pars_fragment,THREE.ShaderChunk.lightmap_pars_fragment,THREE.ShaderChunk.envmap_pars_fragment,THREE.ShaderChunk.fog_pars_fragment,THREE.ShaderChunk.shadowmap_pars_fragment,"void main() {\ngl_FragColor = vec4( vec3 ( 1.0 ), opacity );",THREE.ShaderChunk.map_fragment,THREE.ShaderChunk.alphatest_fragment,"#ifdef DOUBLE_SIDED\nif ( gl_FrontFacing )\ngl_FragColor.xyz *= vLightFront;\nelse\ngl_FragColor.xyz *= vLightBack;\n#else\ngl_FragColor.xyz *= vLightFront;\n#endif",
THREE.ShaderChunk.lightmap_fragment,THREE.ShaderChunk.color_fragment,THREE.ShaderChunk.envmap_fragment,THREE.ShaderChunk.shadowmap_fragment,THREE.ShaderChunk.linear_to_gamma_fragment,THREE.ShaderChunk.fog_fragment,"}"].join("\n")},phong:{uniforms:THREE.UniformsUtils.merge([THREE.UniformsLib.common,THREE.UniformsLib.fog,THREE.UniformsLib.lights,THREE.UniformsLib.shadowmap,{ambient:{type:"c",value:new THREE.Color(16777215)},emissive:{type:"c",value:new THREE.Color(0)},specular:{type:"c",value:new THREE.Color(1118481)},
shininess:{type:"f",value:30},wrapRGB:{type:"v3",value:new THREE.Vector3(1,1,1)}}]),vertexShader:["varying vec3 vViewPosition;\nvarying vec3 vNormal;",THREE.ShaderChunk.map_pars_vertex,THREE.ShaderChunk.lightmap_pars_vertex,THREE.ShaderChunk.envmap_pars_vertex,THREE.ShaderChunk.lights_phong_pars_vertex,THREE.ShaderChunk.color_pars_vertex,THREE.ShaderChunk.skinning_pars_vertex,THREE.ShaderChunk.morphtarget_pars_vertex,THREE.ShaderChunk.shadowmap_pars_vertex,"void main() {\nvec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
THREE.ShaderChunk.map_vertex,THREE.ShaderChunk.lightmap_vertex,THREE.ShaderChunk.envmap_vertex,THREE.ShaderChunk.color_vertex,"#ifndef USE_ENVMAP\nvec4 mPosition = objectMatrix * vec4( position, 1.0 );\n#endif\nvViewPosition = -mvPosition.xyz;",THREE.ShaderChunk.morphnormal_vertex,"vNormal = transformedNormal;",THREE.ShaderChunk.lights_phong_vertex,THREE.ShaderChunk.skinning_vertex,THREE.ShaderChunk.morphtarget_vertex,THREE.ShaderChunk.default_vertex,THREE.ShaderChunk.shadowmap_vertex,"}"].join("\n"),
fragmentShader:["uniform vec3 diffuse;\nuniform float opacity;\nuniform vec3 ambient;\nuniform vec3 emissive;\nuniform vec3 specular;\nuniform float shininess;",THREE.ShaderChunk.color_pars_fragment,THREE.ShaderChunk.map_pars_fragment,THREE.ShaderChunk.lightmap_pars_fragment,THREE.ShaderChunk.envmap_pars_fragment,THREE.ShaderChunk.fog_pars_fragment,THREE.ShaderChunk.lights_phong_pars_fragment,THREE.ShaderChunk.shadowmap_pars_fragment,"void main() {\ngl_FragColor = vec4( vec3 ( 1.0 ), opacity );",
THREE.ShaderChunk.map_fragment,THREE.ShaderChunk.alphatest_fragment,THREE.ShaderChunk.lights_phong_fragment,THREE.ShaderChunk.lightmap_fragment,THREE.ShaderChunk.color_fragment,THREE.ShaderChunk.envmap_fragment,THREE.ShaderChunk.shadowmap_fragment,THREE.ShaderChunk.linear_to_gamma_fragment,THREE.ShaderChunk.fog_fragment,"}"].join("\n")},particle_basic:{uniforms:THREE.UniformsUtils.merge([THREE.UniformsLib.particle,THREE.UniformsLib.shadowmap]),vertexShader:["uniform float size;\nuniform float scale;",
THREE.ShaderChunk.color_pars_vertex,THREE.ShaderChunk.shadowmap_pars_vertex,"void main() {",THREE.ShaderChunk.color_vertex,"vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );\n#ifdef USE_SIZEATTENUATION\ngl_PointSize = size * ( scale / length( mvPosition.xyz ) );\n#else\ngl_PointSize = size;\n#endif\ngl_Position = projectionMatrix * mvPosition;",THREE.ShaderChunk.shadowmap_vertex,"}"].join("\n"),fragmentShader:["uniform vec3 psColor;\nuniform float opacity;",THREE.ShaderChunk.color_pars_fragment,
THREE.ShaderChunk.map_particle_pars_fragment,THREE.ShaderChunk.fog_pars_fragment,THREE.ShaderChunk.shadowmap_pars_fragment,"void main() {\ngl_FragColor = vec4( psColor, opacity );",THREE.ShaderChunk.map_particle_fragment,THREE.ShaderChunk.alphatest_fragment,THREE.ShaderChunk.color_fragment,THREE.ShaderChunk.shadowmap_fragment,THREE.ShaderChunk.fog_fragment,"}"].join("\n")},depthRGBA:{uniforms:{},vertexShader:[THREE.ShaderChunk.morphtarget_pars_vertex,"void main() {\nvec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
THREE.ShaderChunk.morphtarget_vertex,THREE.ShaderChunk.default_vertex,"}"].join("\n"),fragmentShader:"vec4 pack_depth( const in float depth ) {\nconst vec4 bit_shift = vec4( 256.0 * 256.0 * 256.0, 256.0 * 256.0, 256.0, 1.0 );\nconst vec4 bit_mask  = vec4( 0.0, 1.0 / 256.0, 1.0 / 256.0, 1.0 / 256.0 );\nvec4 res = fract( depth * bit_shift );\nres -= res.xxyz * bit_mask;\nreturn res;\n}\nvoid main() {\ngl_FragData[ 0 ] = pack_depth( gl_FragCoord.z );\n}"}};
THREE.WebGLRenderer=function(a){function b(a,b){var c=a.vertices.length,d=b.material;if(d.attributes){if(a.__webglCustomAttributesList===void 0)a.__webglCustomAttributesList=[];for(var e in d.attributes){var f=d.attributes[e];if(!f.__webglInitialized||f.createUniqueBuffers){f.__webglInitialized=true;var h=1;f.type==="v2"?h=2:f.type==="v3"?h=3:f.type==="v4"?h=4:f.type==="c"&&(h=3);f.size=h;f.array=new Float32Array(c*h);f.buffer=g.createBuffer();f.buffer.belongsToAttribute=e;f.needsUpdate=true}a.__webglCustomAttributesList.push(f)}}}
function c(a,b){if(a.material&&!(a.material instanceof THREE.MeshFaceMaterial))return a.material;if(b.materialIndex>=0)return a.geometry.materials[b.materialIndex]}function d(a){return a instanceof THREE.MeshBasicMaterial&&!a.envMap||a instanceof THREE.MeshDepthMaterial?false:a&&a.shading!==void 0&&a.shading===THREE.SmoothShading?THREE.SmoothShading:THREE.FlatShading}function e(a){return a.map||a.lightMap||a instanceof THREE.ShaderMaterial?true:false}function f(a,b,c){var d,e,f,h,i=a.vertices;h=i.length;
var k=a.colors,j=k.length,l=a.__vertexArray,n=a.__colorArray,o=a.__sortArray,m=a.verticesNeedUpdate,p=a.colorsNeedUpdate,s=a.__webglCustomAttributesList;if(c.sortParticles){Ub.copy(qb);Ub.multiplySelf(c.matrixWorld);for(d=0;d<h;d++){e=i[d];Ra.copy(e);Ub.multiplyVector3(Ra);o[d]=[Ra.z,d]}o.sort(function(a,b){return b[0]-a[0]});for(d=0;d<h;d++){e=i[o[d][1]];f=d*3;l[f]=e.x;l[f+1]=e.y;l[f+2]=e.z}for(d=0;d<j;d++){f=d*3;e=k[o[d][1]];n[f]=e.r;n[f+1]=e.g;n[f+2]=e.b}if(s){k=0;for(j=s.length;k<j;k++){i=s[k];
if(i.boundTo===void 0||i.boundTo==="vertices"){f=0;e=i.value.length;if(i.size===1)for(d=0;d<e;d++){h=o[d][1];i.array[d]=i.value[h]}else if(i.size===2)for(d=0;d<e;d++){h=o[d][1];h=i.value[h];i.array[f]=h.x;i.array[f+1]=h.y;f=f+2}else if(i.size===3)if(i.type==="c")for(d=0;d<e;d++){h=o[d][1];h=i.value[h];i.array[f]=h.r;i.array[f+1]=h.g;i.array[f+2]=h.b;f=f+3}else for(d=0;d<e;d++){h=o[d][1];h=i.value[h];i.array[f]=h.x;i.array[f+1]=h.y;i.array[f+2]=h.z;f=f+3}else if(i.size===4)for(d=0;d<e;d++){h=o[d][1];
h=i.value[h];i.array[f]=h.x;i.array[f+1]=h.y;i.array[f+2]=h.z;i.array[f+3]=h.w;f=f+4}}}}}else{if(m)for(d=0;d<h;d++){e=i[d];f=d*3;l[f]=e.x;l[f+1]=e.y;l[f+2]=e.z}if(p)for(d=0;d<j;d++){e=k[d];f=d*3;n[f]=e.r;n[f+1]=e.g;n[f+2]=e.b}if(s){k=0;for(j=s.length;k<j;k++){i=s[k];if(i.needsUpdate&&(i.boundTo===void 0||i.boundTo==="vertices")){e=i.value.length;f=0;if(i.size===1)for(d=0;d<e;d++)i.array[d]=i.value[d];else if(i.size===2)for(d=0;d<e;d++){h=i.value[d];i.array[f]=h.x;i.array[f+1]=h.y;f=f+2}else if(i.size===
3)if(i.type==="c")for(d=0;d<e;d++){h=i.value[d];i.array[f]=h.r;i.array[f+1]=h.g;i.array[f+2]=h.b;f=f+3}else for(d=0;d<e;d++){h=i.value[d];i.array[f]=h.x;i.array[f+1]=h.y;i.array[f+2]=h.z;f=f+3}else if(i.size===4)for(d=0;d<e;d++){h=i.value[d];i.array[f]=h.x;i.array[f+1]=h.y;i.array[f+2]=h.z;i.array[f+3]=h.w;f=f+4}}}}}if(m||c.sortParticles){g.bindBuffer(g.ARRAY_BUFFER,a.__webglVertexBuffer);g.bufferData(g.ARRAY_BUFFER,l,b)}if(p||c.sortParticles){g.bindBuffer(g.ARRAY_BUFFER,a.__webglColorBuffer);g.bufferData(g.ARRAY_BUFFER,
n,b)}if(s){k=0;for(j=s.length;k<j;k++){i=s[k];if(i.needsUpdate||c.sortParticles){g.bindBuffer(g.ARRAY_BUFFER,i.buffer);g.bufferData(g.ARRAY_BUFFER,i.array,b)}}}}function h(a,b){return b.z-a.z}function i(a,b,c){if(a.length)for(var d=0,e=a.length;d<e;d++){J=za=null;$=R=pa=Q=Oa=La=O=-1;rb=true;a[d].render(b,c,jc,kc);J=za=null;$=R=pa=Q=Oa=La=O=-1;rb=true}}function l(a,b,c,d,e,g,f,h){var i,k,j,l;if(b){k=a.length-1;l=b=-1}else{k=0;b=a.length;l=1}for(var n=k;n!==b;n=n+l){i=a[n];if(i.render){k=i.object;j=
i.buffer;if(h)i=h;else{i=i[c];if(!i)continue;f&&D.setBlending(i.blending,i.blendEquation,i.blendSrc,i.blendDst);D.setDepthTest(i.depthTest);D.setDepthWrite(i.depthWrite);p(i.polygonOffset,i.polygonOffsetFactor,i.polygonOffsetUnits)}D.setObjectFaces(k);j instanceof THREE.BufferGeometry?D.renderBufferDirect(d,e,g,i,j,k):D.renderBuffer(d,e,g,i,j,k)}}}function k(a,b,c,d,e,g,f){for(var h,i,k=0,j=a.length;k<j;k++){h=a[k];i=h.object;if(i.visible){if(f)h=f;else{h=h[b];if(!h)continue;g&&D.setBlending(h.blending,
h.blendEquation,h.blendSrc,h.blendDst);D.setDepthTest(h.depthTest);D.setDepthWrite(h.depthWrite);p(h.polygonOffset,h.polygonOffsetFactor,h.polygonOffsetUnits)}D.renderImmediateObject(c,d,e,h,i)}}}function j(a,b,c){a.push({buffer:b,object:c,opaque:null,transparent:null})}function m(a){for(var b in a.attributes)if(a.attributes[b].needsUpdate)return true;return false}function n(a){for(var b in a.attributes)a.attributes[b].needsUpdate=false}function o(a,b,c,d,e){if(!d.program||d.needsUpdate){D.initMaterial(d,
b,c,e);d.needsUpdate=false}if(d.morphTargets&&!e.__webglMorphTargetInfluences){e.__webglMorphTargetInfluences=new Float32Array(D.maxMorphTargets);for(var f=0,h=D.maxMorphTargets;f<h;f++)e.__webglMorphTargetInfluences[f]=0}var i=false,f=d.program,h=f.uniforms,k=d.uniforms;if(f!==za){g.useProgram(f);za=f;i=true}if(d.id!==$){$=d.id;i=true}if(i||a!==J){g.uniformMatrix4fv(h.projectionMatrix,false,a._projectionMatrixArray);a!==J&&(J=a)}if(i){if(c&&d.fog){k.fogColor.value=c.color;if(c instanceof THREE.Fog){k.fogNear.value=
c.near;k.fogFar.value=c.far}else if(c instanceof THREE.FogExp2)k.fogDensity.value=c.density}if(d instanceof THREE.MeshPhongMaterial||d instanceof THREE.MeshLambertMaterial||d.lights){if(rb){for(var j,l=0,n=0,o=0,m,p,s,q=lc,y=q.directional.colors,A=q.directional.positions,B=q.point.colors,K=q.point.positions,N=q.point.distances,Y=q.spot.colors,O=q.spot.positions,ba=q.spot.distances,ja=q.spot.directions,ca=q.spot.angles,Q=q.spot.exponents,Z=0,I=0,R=0,la=s=0,c=la=0,i=b.length;c<i;c++){j=b[c];if(!j.onlyShadow){m=
j.color;p=j.intensity;s=j.distance;if(j instanceof THREE.AmbientLight)if(D.gammaInput){l=l+m.r*m.r;n=n+m.g*m.g;o=o+m.b*m.b}else{l=l+m.r;n=n+m.g;o=o+m.b}else if(j instanceof THREE.DirectionalLight){s=Z*3;if(D.gammaInput){y[s]=m.r*m.r*p*p;y[s+1]=m.g*m.g*p*p;y[s+2]=m.b*m.b*p*p}else{y[s]=m.r*p;y[s+1]=m.g*p;y[s+2]=m.b*p}Ca.copy(j.matrixWorld.getPosition());Ca.subSelf(j.target.matrixWorld.getPosition());Ca.normalize();A[s]=Ca.x;A[s+1]=Ca.y;A[s+2]=Ca.z;Z=Z+1}else if(j instanceof THREE.PointLight){la=I*3;
if(D.gammaInput){B[la]=m.r*m.r*p*p;B[la+1]=m.g*m.g*p*p;B[la+2]=m.b*m.b*p*p}else{B[la]=m.r*p;B[la+1]=m.g*p;B[la+2]=m.b*p}m=j.matrixWorld.getPosition();K[la]=m.x;K[la+1]=m.y;K[la+2]=m.z;N[I]=s;I=I+1}else if(j instanceof THREE.SpotLight){la=R*3;if(D.gammaInput){Y[la]=m.r*m.r*p*p;Y[la+1]=m.g*m.g*p*p;Y[la+2]=m.b*m.b*p*p}else{Y[la]=m.r*p;Y[la+1]=m.g*p;Y[la+2]=m.b*p}m=j.matrixWorld.getPosition();O[la]=m.x;O[la+1]=m.y;O[la+2]=m.z;ba[R]=s;Ca.copy(m);Ca.subSelf(j.target.matrixWorld.getPosition());Ca.normalize();
ja[la]=Ca.x;ja[la+1]=Ca.y;ja[la+2]=Ca.z;ca[R]=Math.cos(j.angle);Q[R]=j.exponent;R=R+1}}}c=Z*3;for(i=y.length;c<i;c++)y[c]=0;c=I*3;for(i=B.length;c<i;c++)B[c]=0;c=R*3;for(i=Y.length;c<i;c++)Y[c]=0;q.directional.length=Z;q.point.length=I;q.spot.length=R;q.ambient[0]=l;q.ambient[1]=n;q.ambient[2]=o;rb=false}c=lc;k.ambientLightColor.value=c.ambient;k.directionalLightColor.value=c.directional.colors;k.directionalLightDirection.value=c.directional.positions;k.pointLightColor.value=c.point.colors;k.pointLightPosition.value=
c.point.positions;k.pointLightDistance.value=c.point.distances;k.spotLightColor.value=c.spot.colors;k.spotLightPosition.value=c.spot.positions;k.spotLightDistance.value=c.spot.distances;k.spotLightDirection.value=c.spot.directions;k.spotLightAngle.value=c.spot.angles;k.spotLightExponent.value=c.spot.exponents}if(d instanceof THREE.MeshBasicMaterial||d instanceof THREE.MeshLambertMaterial||d instanceof THREE.MeshPhongMaterial){k.opacity.value=d.opacity;D.gammaInput?k.diffuse.value.copyGammaToLinear(d.color):
k.diffuse.value=d.color;(k.map.texture=d.map)&&k.offsetRepeat.value.set(d.map.offset.x,d.map.offset.y,d.map.repeat.x,d.map.repeat.y);k.lightMap.texture=d.lightMap;k.envMap.texture=d.envMap;k.flipEnvMap.value=d.envMap instanceof THREE.WebGLRenderTargetCube?1:-1;k.reflectivity.value=d.reflectivity;k.refractionRatio.value=d.refractionRatio;k.combine.value=d.combine;k.useRefract.value=d.envMap&&d.envMap.mapping instanceof THREE.CubeRefractionMapping}if(d instanceof THREE.LineBasicMaterial){k.diffuse.value=
d.color;k.opacity.value=d.opacity}else if(d instanceof THREE.ParticleBasicMaterial){k.psColor.value=d.color;k.opacity.value=d.opacity;k.size.value=d.size;k.scale.value=H.height/2;k.map.texture=d.map}else if(d instanceof THREE.MeshPhongMaterial){k.shininess.value=d.shininess;if(D.gammaInput){k.ambient.value.copyGammaToLinear(d.ambient);k.emissive.value.copyGammaToLinear(d.emissive);k.specular.value.copyGammaToLinear(d.specular)}else{k.ambient.value=d.ambient;k.emissive.value=d.emissive;k.specular.value=
d.specular}d.wrapAround&&k.wrapRGB.value.copy(d.wrapRGB)}else if(d instanceof THREE.MeshLambertMaterial){if(D.gammaInput){k.ambient.value.copyGammaToLinear(d.ambient);k.emissive.value.copyGammaToLinear(d.emissive)}else{k.ambient.value=d.ambient;k.emissive.value=d.emissive}d.wrapAround&&k.wrapRGB.value.copy(d.wrapRGB)}else if(d instanceof THREE.MeshDepthMaterial){k.mNear.value=a.near;k.mFar.value=a.far;k.opacity.value=d.opacity}else if(d instanceof THREE.MeshNormalMaterial)k.opacity.value=d.opacity;
if(e.receiveShadow&&!d._shadowPass&&k.shadowMatrix){i=c=0;for(j=b.length;i<j;i++){l=b[i];if(l.castShadow&&(l instanceof THREE.SpotLight||l instanceof THREE.DirectionalLight&&!l.shadowCascade)){k.shadowMap.texture[c]=l.shadowMap;k.shadowMapSize.value[c]=l.shadowMapSize;k.shadowMatrix.value[c]=l.shadowMatrix;k.shadowDarkness.value[c]=l.shadowDarkness;k.shadowBias.value[c]=l.shadowBias;c++}}}b=d.uniformsList;k=0;for(c=b.length;k<c;k++)if(l=f.uniforms[b[k][1]]){i=b[k][0];n=i.type;j=i.value;switch(n){case "i":g.uniform1i(l,
j);break;case "f":g.uniform1f(l,j);break;case "v2":g.uniform2f(l,j.x,j.y);break;case "v3":g.uniform3f(l,j.x,j.y,j.z);break;case "v4":g.uniform4f(l,j.x,j.y,j.z,j.w);break;case "c":g.uniform3f(l,j.r,j.g,j.b);break;case "fv1":g.uniform1fv(l,j);break;case "fv":g.uniform3fv(l,j);break;case "v2v":if(!i._array)i._array=new Float32Array(2*j.length);n=0;for(o=j.length;n<o;n++){q=n*2;i._array[q]=j[n].x;i._array[q+1]=j[n].y}g.uniform2fv(l,i._array);break;case "v3v":if(!i._array)i._array=new Float32Array(3*j.length);
n=0;for(o=j.length;n<o;n++){q=n*3;i._array[q]=j[n].x;i._array[q+1]=j[n].y;i._array[q+2]=j[n].z}g.uniform3fv(l,i._array);break;case "v4v":if(!i._array)i._array=new Float32Array(4*j.length);n=0;for(o=j.length;n<o;n++){q=n*4;i._array[q]=j[n].x;i._array[q+1]=j[n].y;i._array[q+2]=j[n].z;i._array[q+3]=j[n].w}g.uniform4fv(l,i._array);break;case "m4":if(!i._array)i._array=new Float32Array(16);j.flattenToArray(i._array);g.uniformMatrix4fv(l,false,i._array);break;case "m4v":if(!i._array)i._array=new Float32Array(16*
j.length);n=0;for(o=j.length;n<o;n++)j[n].flattenToArrayOffset(i._array,n*16);g.uniformMatrix4fv(l,false,i._array);break;case "t":g.uniform1i(l,j);l=i.texture;if(!l)continue;if(l.image instanceof Array&&l.image.length===6){i=l;if(i.image.length===6)if(i.needsUpdate){if(!i.image.__webglTextureCube)i.image.__webglTextureCube=g.createTexture();g.activeTexture(g.TEXTURE0+j);g.bindTexture(g.TEXTURE_CUBE_MAP,i.image.__webglTextureCube);j=[];for(l=0;l<6;l++){n=j;o=l;if(D.autoScaleCubemaps){q=i.image[l];
A=yc;if(!(q.width<=A&&q.height<=A)){B=Math.max(q.width,q.height);y=Math.floor(q.width*A/B);A=Math.floor(q.height*A/B);B=document.createElement("canvas");B.width=y;B.height=A;B.getContext("2d").drawImage(q,0,0,q.width,q.height,0,0,y,A);q=B}}else q=i.image[l];n[o]=q}l=j[0];n=(l.width&l.width-1)===0&&(l.height&l.height-1)===0;o=u(i.format);q=u(i.type);w(g.TEXTURE_CUBE_MAP,i,n);for(l=0;l<6;l++)g.texImage2D(g.TEXTURE_CUBE_MAP_POSITIVE_X+l,0,o,o,q,j[l]);i.generateMipmaps&&n&&g.generateMipmap(g.TEXTURE_CUBE_MAP);
i.needsUpdate=false;if(i.onUpdate)i.onUpdate()}else{g.activeTexture(g.TEXTURE0+j);g.bindTexture(g.TEXTURE_CUBE_MAP,i.image.__webglTextureCube)}}else if(l instanceof THREE.WebGLRenderTargetCube){i=l;g.activeTexture(g.TEXTURE0+j);g.bindTexture(g.TEXTURE_CUBE_MAP,i.__webglTexture)}else D.setTexture(l,j);break;case "tv":if(!i._array){i._array=[];n=0;for(o=i.texture.length;n<o;n++)i._array[n]=j+n}g.uniform1iv(l,i._array);n=0;for(o=i.texture.length;n<o;n++)(l=i.texture[n])&&D.setTexture(l,i._array[n])}}if((d instanceof
THREE.ShaderMaterial||d instanceof THREE.MeshPhongMaterial||d.envMap)&&h.cameraPosition!==null){b=a.matrixWorld.getPosition();g.uniform3f(h.cameraPosition,b.x,b.y,b.z)}(d instanceof THREE.MeshPhongMaterial||d instanceof THREE.MeshLambertMaterial||d instanceof THREE.ShaderMaterial||d.skinning)&&h.viewMatrix!==null&&g.uniformMatrix4fv(h.viewMatrix,false,a._viewMatrixArray);d.skinning&&g.uniformMatrix4fv(h.boneGlobalMatrices,false,e.boneMatrices)}g.uniformMatrix4fv(h.modelViewMatrix,false,e._modelViewMatrix.elements);
h.normalMatrix&&g.uniformMatrix3fv(h.normalMatrix,false,e._normalMatrix.elements);h.objectMatrix!==null&&g.uniformMatrix4fv(h.objectMatrix,false,e.matrixWorld.elements);return f}function s(a,b){a._modelViewMatrix.multiply(b.matrixWorldInverse,a.matrixWorld);a._normalMatrix.getInverse(a._modelViewMatrix);a._normalMatrix.transpose()}function p(a,b,c){if(Ta!==a){a?g.enable(g.POLYGON_OFFSET_FILL):g.disable(g.POLYGON_OFFSET_FILL);Ta=a}if(a&&(Vb!==b||Wb!==c)){g.polygonOffset(b,c);Vb=b;Wb=c}}function q(a,
b){var c;a==="fragment"?c=g.createShader(g.FRAGMENT_SHADER):a==="vertex"&&(c=g.createShader(g.VERTEX_SHADER));g.shaderSource(c,b);g.compileShader(c);if(!g.getShaderParameter(c,g.COMPILE_STATUS)){console.error(g.getShaderInfoLog(c));console.error(b);return null}return c}function w(a,b,c){if(c){g.texParameteri(a,g.TEXTURE_WRAP_S,u(b.wrapS));g.texParameteri(a,g.TEXTURE_WRAP_T,u(b.wrapT));g.texParameteri(a,g.TEXTURE_MAG_FILTER,u(b.magFilter));g.texParameteri(a,g.TEXTURE_MIN_FILTER,u(b.minFilter))}else{g.texParameteri(a,
g.TEXTURE_WRAP_S,g.CLAMP_TO_EDGE);g.texParameteri(a,g.TEXTURE_WRAP_T,g.CLAMP_TO_EDGE);g.texParameteri(a,g.TEXTURE_MAG_FILTER,y(b.magFilter));g.texParameteri(a,g.TEXTURE_MIN_FILTER,y(b.minFilter))}}function A(a,b){g.bindRenderbuffer(g.RENDERBUFFER,a);if(b.depthBuffer&&!b.stencilBuffer){g.renderbufferStorage(g.RENDERBUFFER,g.DEPTH_COMPONENT16,b.width,b.height);g.framebufferRenderbuffer(g.FRAMEBUFFER,g.DEPTH_ATTACHMENT,g.RENDERBUFFER,a)}else if(b.depthBuffer&&b.stencilBuffer){g.renderbufferStorage(g.RENDERBUFFER,
g.DEPTH_STENCIL,b.width,b.height);g.framebufferRenderbuffer(g.FRAMEBUFFER,g.DEPTH_STENCIL_ATTACHMENT,g.RENDERBUFFER,a)}else g.renderbufferStorage(g.RENDERBUFFER,g.RGBA4,b.width,b.height)}function y(a){switch(a){case THREE.NearestFilter:case THREE.NearestMipMapNearestFilter:case THREE.NearestMipMapLinearFilter:return g.NEAREST;default:return g.LINEAR}}function u(a){switch(a){case THREE.RepeatWrapping:return g.REPEAT;case THREE.ClampToEdgeWrapping:return g.CLAMP_TO_EDGE;case THREE.MirroredRepeatWrapping:return g.MIRRORED_REPEAT;
case THREE.NearestFilter:return g.NEAREST;case THREE.NearestMipMapNearestFilter:return g.NEAREST_MIPMAP_NEAREST;case THREE.NearestMipMapLinearFilter:return g.NEAREST_MIPMAP_LINEAR;case THREE.LinearFilter:return g.LINEAR;case THREE.LinearMipMapNearestFilter:return g.LINEAR_MIPMAP_NEAREST;case THREE.LinearMipMapLinearFilter:return g.LINEAR_MIPMAP_LINEAR;case THREE.ByteType:return g.BYTE;case THREE.UnsignedByteType:return g.UNSIGNED_BYTE;case THREE.ShortType:return g.SHORT;case THREE.UnsignedShortType:return g.UNSIGNED_SHORT;
case THREE.IntType:return g.INT;case THREE.UnsignedIntType:return g.UNSIGNED_INT;case THREE.FloatType:return g.FLOAT;case THREE.AlphaFormat:return g.ALPHA;case THREE.RGBFormat:return g.RGB;case THREE.RGBAFormat:return g.RGBA;case THREE.LuminanceFormat:return g.LUMINANCE;case THREE.LuminanceAlphaFormat:return g.LUMINANCE_ALPHA;case THREE.AddEquation:return g.FUNC_ADD;case THREE.SubtractEquation:return g.FUNC_SUBTRACT;case THREE.ReverseSubtractEquation:return g.FUNC_REVERSE_SUBTRACT;case THREE.ZeroFactor:return g.ZERO;
case THREE.OneFactor:return g.ONE;case THREE.SrcColorFactor:return g.SRC_COLOR;case THREE.OneMinusSrcColorFactor:return g.ONE_MINUS_SRC_COLOR;case THREE.SrcAlphaFactor:return g.SRC_ALPHA;case THREE.OneMinusSrcAlphaFactor:return g.ONE_MINUS_SRC_ALPHA;case THREE.DstAlphaFactor:return g.DST_ALPHA;case THREE.OneMinusDstAlphaFactor:return g.ONE_MINUS_DST_ALPHA;case THREE.DstColorFactor:return g.DST_COLOR;case THREE.OneMinusDstColorFactor:return g.ONE_MINUS_DST_COLOR;case THREE.SrcAlphaSaturateFactor:return g.SRC_ALPHA_SATURATE}return 0}
console.log("THREE.WebGLRenderer",THREE.REVISION);var a=a||{},H=a.canvas!==void 0?a.canvas:document.createElement("canvas"),B=a.precision!==void 0?a.precision:"highp",K=a.alpha!==void 0?a.alpha:true,N=a.premultipliedAlpha!==void 0?a.premultipliedAlpha:true,Y=a.antialias!==void 0?a.antialias:false,ca=a.stencil!==void 0?a.stencil:true,I=a.preserveDrawingBuffer!==void 0?a.preserveDrawingBuffer:false,ba=a.clearColor!==void 0?new THREE.Color(a.clearColor):new THREE.Color(0),ja=a.clearAlpha!==void 0?a.clearAlpha:
0,ya=a.maxLights!==void 0?a.maxLights:4;this.domElement=H;this.context=null;this.autoUpdateScene=this.autoUpdateObjects=this.sortObjects=this.autoClearStencil=this.autoClearDepth=this.autoClearColor=this.autoClear=true;this.shadowMapEnabled=this.physicallyBasedShading=this.gammaOutput=this.gammaInput=false;this.shadowMapCullFrontFaces=this.shadowMapSoft=this.shadowMapAutoUpdate=true;this.shadowMapCascade=this.shadowMapDebug=false;this.maxMorphTargets=8;this.maxMorphNormals=4;this.autoScaleCubemaps=
true;this.renderPluginsPre=[];this.renderPluginsPost=[];this.info={memory:{programs:0,geometries:0,textures:0},render:{calls:0,vertices:0,faces:0,points:0}};var D=this,g,Na=[],za=null,Da=null,$=-1,R=null,J=null,Z=0,Q=-1,pa=-1,O=-1,sa=-1,Ga=-1,Ha=-1,La=-1,Oa=-1,Ta=null,Vb=null,Wb=null,Db=null,Eb=0,Fb=0,Xb=0,Gb=0,jc=0,kc=0,Yb=new THREE.Frustum,qb=new THREE.Matrix4,Ub=new THREE.Matrix4,Ra=new THREE.Vector4,Ca=new THREE.Vector3,rb=true,lc={ambient:[0,0,0],directional:{length:0,colors:[],positions:[]},
point:{length:0,colors:[],positions:[],distances:[]},spot:{length:0,colors:[],positions:[],distances:[],directions:[],angles:[],exponents:[]}};g=function(){var a;try{if(!(a=H.getContext("experimental-webgl",{alpha:K,premultipliedAlpha:N,antialias:Y,stencil:ca,preserveDrawingBuffer:I})))throw"Error creating WebGL context.";}catch(b){console.error(b)}a.getExtension("OES_texture_float")||console.log("THREE.WebGLRenderer: Float textures not supported.");return a}();g.clearColor(0,0,0,1);g.clearDepth(1);
g.clearStencil(0);g.enable(g.DEPTH_TEST);g.depthFunc(g.LEQUAL);g.frontFace(g.CCW);g.cullFace(g.BACK);g.enable(g.CULL_FACE);g.enable(g.BLEND);g.blendEquation(g.FUNC_ADD);g.blendFunc(g.SRC_ALPHA,g.ONE_MINUS_SRC_ALPHA);g.clearColor(ba.r,ba.g,ba.b,ja);this.context=g;var Zb=g.getParameter(g.MAX_VERTEX_TEXTURE_IMAGE_UNITS);g.getParameter(g.MAX_TEXTURE_SIZE);var yc=g.getParameter(g.MAX_CUBE_MAP_TEXTURE_SIZE);this.getContext=function(){return g};this.supportsVertexTextures=function(){return Zb>0};this.setSize=
function(a,b){H.width=a;H.height=b;this.setViewport(0,0,H.width,H.height)};this.setViewport=function(a,b,c,d){Eb=a;Fb=b;Xb=c;Gb=d;g.viewport(Eb,Fb,Xb,Gb)};this.setScissor=function(a,b,c,d){g.scissor(a,b,c,d)};this.enableScissorTest=function(a){a?g.enable(g.SCISSOR_TEST):g.disable(g.SCISSOR_TEST)};this.setClearColorHex=function(a,b){ba.setHex(a);ja=b;g.clearColor(ba.r,ba.g,ba.b,ja)};this.setClearColor=function(a,b){ba.copy(a);ja=b;g.clearColor(ba.r,ba.g,ba.b,ja)};this.getClearColor=function(){return ba};
this.getClearAlpha=function(){return ja};this.clear=function(a,b,c){var d=0;if(a===void 0||a)d=d|g.COLOR_BUFFER_BIT;if(b===void 0||b)d=d|g.DEPTH_BUFFER_BIT;if(c===void 0||c)d=d|g.STENCIL_BUFFER_BIT;g.clear(d)};this.clearTarget=function(a,b,c,d){this.setRenderTarget(a);this.clear(b,c,d)};this.addPostPlugin=function(a){a.init(this);this.renderPluginsPost.push(a)};this.addPrePlugin=function(a){a.init(this);this.renderPluginsPre.push(a)};this.deallocateObject=function(a){if(a.__webglInit){a.__webglInit=
false;delete a._modelViewMatrix;delete a._normalMatrix;delete a._normalMatrixArray;delete a._modelViewMatrixArray;delete a._objectMatrixArray;if(a instanceof THREE.Mesh)for(var b in a.geometry.geometryGroups){var c=a.geometry.geometryGroups[b];g.deleteBuffer(c.__webglVertexBuffer);g.deleteBuffer(c.__webglNormalBuffer);g.deleteBuffer(c.__webglTangentBuffer);g.deleteBuffer(c.__webglColorBuffer);g.deleteBuffer(c.__webglUVBuffer);g.deleteBuffer(c.__webglUV2Buffer);g.deleteBuffer(c.__webglSkinVertexABuffer);
g.deleteBuffer(c.__webglSkinVertexBBuffer);g.deleteBuffer(c.__webglSkinIndicesBuffer);g.deleteBuffer(c.__webglSkinWeightsBuffer);g.deleteBuffer(c.__webglFaceBuffer);g.deleteBuffer(c.__webglLineBuffer);var d=void 0,e=void 0;if(c.numMorphTargets){d=0;for(e=c.numMorphTargets;d<e;d++)g.deleteBuffer(c.__webglMorphTargetsBuffers[d])}if(c.numMorphNormals){d=0;for(e=c.numMorphNormals;d<e;d++)g.deleteBuffer(c.__webglMorphNormalsBuffers[d])}if(c.__webglCustomAttributesList){d=void 0;for(d in c.__webglCustomAttributesList)g.deleteBuffer(c.__webglCustomAttributesList[d].buffer)}D.info.memory.geometries--}else if(a instanceof
THREE.Line){a=a.geometry;g.deleteBuffer(a.__webglVertexBuffer);g.deleteBuffer(a.__webglColorBuffer);D.info.memory.geometries--}}};this.deallocateTexture=function(a){if(a.__webglInit){a.__webglInit=false;g.deleteTexture(a.__webglTexture);D.info.memory.textures--}};this.deallocateRenderTarget=function(a){if(a&&a.__webglTexture){g.deleteTexture(a.__webglTexture);if(a instanceof THREE.WebGLRenderTargetCube)for(var b=0;b<6;b++){g.deleteFramebuffer(a.__webglFramebuffer[b]);g.deleteRenderbuffer(a.__webglRenderbuffer[b])}else{g.deleteFramebuffer(a.__webglFramebuffer);
g.deleteRenderbuffer(a.__webglRenderbuffer)}}};this.updateShadowMap=function(a,b){za=null;$=R=Oa=La=O=-1;rb=true;pa=Q=-1;this.shadowMapPlugin.update(a,b)};this.renderBufferImmediate=function(a,b,c){if(!a.__webglVertexBuffer)a.__webglVertexBuffer=g.createBuffer();if(!a.__webglNormalBuffer)a.__webglNormalBuffer=g.createBuffer();if(a.hasPos){g.bindBuffer(g.ARRAY_BUFFER,a.__webglVertexBuffer);g.bufferData(g.ARRAY_BUFFER,a.positionArray,g.DYNAMIC_DRAW);g.enableVertexAttribArray(b.attributes.position);
g.vertexAttribPointer(b.attributes.position,3,g.FLOAT,false,0,0)}if(a.hasNormal){g.bindBuffer(g.ARRAY_BUFFER,a.__webglNormalBuffer);if(c===THREE.FlatShading){var d,e,f,h,i,k,j,l,n,m,o=a.count*3;for(m=0;m<o;m=m+9){c=a.normalArray;d=c[m];e=c[m+1];f=c[m+2];h=c[m+3];k=c[m+4];l=c[m+5];i=c[m+6];j=c[m+7];n=c[m+8];d=(d+h+i)/3;e=(e+k+j)/3;f=(f+l+n)/3;c[m]=d;c[m+1]=e;c[m+2]=f;c[m+3]=d;c[m+4]=e;c[m+5]=f;c[m+6]=d;c[m+7]=e;c[m+8]=f}}g.bufferData(g.ARRAY_BUFFER,a.normalArray,g.DYNAMIC_DRAW);g.enableVertexAttribArray(b.attributes.normal);
g.vertexAttribPointer(b.attributes.normal,3,g.FLOAT,false,0,0)}g.drawArrays(g.TRIANGLES,0,a.count);a.count=0};this.renderBufferDirect=function(a,b,c,d,e,f){if(d.visible!==false){c=o(a,b,c,d,f);a=c.attributes;b=false;d=e.id*16777215+c.id*2+(d.wireframe?1:0);if(d!==R){R=d;b=true}if(f instanceof THREE.Mesh){f=e.offsets;d=0;for(c=f.length;d<c;++d){if(b){g.bindBuffer(g.ARRAY_BUFFER,e.vertexPositionBuffer);g.vertexAttribPointer(a.position,e.vertexPositionBuffer.itemSize,g.FLOAT,false,0,f[d].index*12);if(a.normal>=
0&&e.vertexNormalBuffer){g.bindBuffer(g.ARRAY_BUFFER,e.vertexNormalBuffer);g.vertexAttribPointer(a.normal,e.vertexNormalBuffer.itemSize,g.FLOAT,false,0,f[d].index*12)}if(a.uv>=0&&e.vertexUvBuffer)if(e.vertexUvBuffer){g.bindBuffer(g.ARRAY_BUFFER,e.vertexUvBuffer);g.vertexAttribPointer(a.uv,e.vertexUvBuffer.itemSize,g.FLOAT,false,0,f[d].index*8);g.enableVertexAttribArray(a.uv)}else g.disableVertexAttribArray(a.uv);if(a.color>=0&&e.vertexColorBuffer){g.bindBuffer(g.ARRAY_BUFFER,e.vertexColorBuffer);
g.vertexAttribPointer(a.color,e.vertexColorBuffer.itemSize,g.FLOAT,false,0,f[d].index*16)}g.bindBuffer(g.ELEMENT_ARRAY_BUFFER,e.vertexIndexBuffer)}g.drawElements(g.TRIANGLES,f[d].count,g.UNSIGNED_SHORT,f[d].start*2);D.info.render.calls++;D.info.render.vertices=D.info.render.vertices+f[d].count;D.info.render.faces=D.info.render.faces+f[d].count/3}}}};this.renderBuffer=function(a,b,c,d,e,f){if(d.visible!==false){var h,i,c=o(a,b,c,d,f),b=c.attributes,a=false,c=e.id*16777215+c.id*2+(d.wireframe?1:0);
if(c!==R){R=c;a=true}if(!d.morphTargets&&b.position>=0){if(a){g.bindBuffer(g.ARRAY_BUFFER,e.__webglVertexBuffer);g.vertexAttribPointer(b.position,3,g.FLOAT,false,0,0)}}else if(f.morphTargetBase){c=d.program.attributes;if(f.morphTargetBase!==-1){g.bindBuffer(g.ARRAY_BUFFER,e.__webglMorphTargetsBuffers[f.morphTargetBase]);g.vertexAttribPointer(c.position,3,g.FLOAT,false,0,0)}else if(c.position>=0){g.bindBuffer(g.ARRAY_BUFFER,e.__webglVertexBuffer);g.vertexAttribPointer(c.position,3,g.FLOAT,false,0,
0)}if(f.morphTargetForcedOrder.length){h=0;var k=f.morphTargetForcedOrder;for(i=f.morphTargetInfluences;h<d.numSupportedMorphTargets&&h<k.length;){g.bindBuffer(g.ARRAY_BUFFER,e.__webglMorphTargetsBuffers[k[h]]);g.vertexAttribPointer(c["morphTarget"+h],3,g.FLOAT,false,0,0);if(d.morphNormals){g.bindBuffer(g.ARRAY_BUFFER,e.__webglMorphNormalsBuffers[k[h]]);g.vertexAttribPointer(c["morphNormal"+h],3,g.FLOAT,false,0,0)}f.__webglMorphTargetInfluences[h]=i[k[h]];h++}}else{var k=[],j=-1,l=0;i=f.morphTargetInfluences;
var n,m=i.length;h=0;for(f.morphTargetBase!==-1&&(k[f.morphTargetBase]=true);h<d.numSupportedMorphTargets;){for(n=0;n<m;n++)if(!k[n]&&i[n]>j){l=n;j=i[l]}g.bindBuffer(g.ARRAY_BUFFER,e.__webglMorphTargetsBuffers[l]);g.vertexAttribPointer(c["morphTarget"+h],3,g.FLOAT,false,0,0);if(d.morphNormals){g.bindBuffer(g.ARRAY_BUFFER,e.__webglMorphNormalsBuffers[l]);g.vertexAttribPointer(c["morphNormal"+h],3,g.FLOAT,false,0,0)}f.__webglMorphTargetInfluences[h]=j;k[l]=1;j=-1;h++}}d.program.uniforms.morphTargetInfluences!==
null&&g.uniform1fv(d.program.uniforms.morphTargetInfluences,f.__webglMorphTargetInfluences)}if(a){if(e.__webglCustomAttributesList){h=0;for(i=e.__webglCustomAttributesList.length;h<i;h++){c=e.__webglCustomAttributesList[h];if(b[c.buffer.belongsToAttribute]>=0){g.bindBuffer(g.ARRAY_BUFFER,c.buffer);g.vertexAttribPointer(b[c.buffer.belongsToAttribute],c.size,g.FLOAT,false,0,0)}}}if(b.color>=0){g.bindBuffer(g.ARRAY_BUFFER,e.__webglColorBuffer);g.vertexAttribPointer(b.color,3,g.FLOAT,false,0,0)}if(b.normal>=
0){g.bindBuffer(g.ARRAY_BUFFER,e.__webglNormalBuffer);g.vertexAttribPointer(b.normal,3,g.FLOAT,false,0,0)}if(b.tangent>=0){g.bindBuffer(g.ARRAY_BUFFER,e.__webglTangentBuffer);g.vertexAttribPointer(b.tangent,4,g.FLOAT,false,0,0)}if(b.uv>=0)if(e.__webglUVBuffer){g.bindBuffer(g.ARRAY_BUFFER,e.__webglUVBuffer);g.vertexAttribPointer(b.uv,2,g.FLOAT,false,0,0);g.enableVertexAttribArray(b.uv)}else g.disableVertexAttribArray(b.uv);if(b.uv2>=0)if(e.__webglUV2Buffer){g.bindBuffer(g.ARRAY_BUFFER,e.__webglUV2Buffer);
g.vertexAttribPointer(b.uv2,2,g.FLOAT,false,0,0);g.enableVertexAttribArray(b.uv2)}else g.disableVertexAttribArray(b.uv2);if(d.skinning&&b.skinVertexA>=0&&b.skinVertexB>=0&&b.skinIndex>=0&&b.skinWeight>=0){g.bindBuffer(g.ARRAY_BUFFER,e.__webglSkinVertexABuffer);g.vertexAttribPointer(b.skinVertexA,4,g.FLOAT,false,0,0);g.bindBuffer(g.ARRAY_BUFFER,e.__webglSkinVertexBBuffer);g.vertexAttribPointer(b.skinVertexB,4,g.FLOAT,false,0,0);g.bindBuffer(g.ARRAY_BUFFER,e.__webglSkinIndicesBuffer);g.vertexAttribPointer(b.skinIndex,
4,g.FLOAT,false,0,0);g.bindBuffer(g.ARRAY_BUFFER,e.__webglSkinWeightsBuffer);g.vertexAttribPointer(b.skinWeight,4,g.FLOAT,false,0,0)}}if(f instanceof THREE.Mesh){if(d.wireframe){d=d.wireframeLinewidth;if(d!==Db){g.lineWidth(d);Db=d}a&&g.bindBuffer(g.ELEMENT_ARRAY_BUFFER,e.__webglLineBuffer);g.drawElements(g.LINES,e.__webglLineCount,g.UNSIGNED_SHORT,0)}else{a&&g.bindBuffer(g.ELEMENT_ARRAY_BUFFER,e.__webglFaceBuffer);g.drawElements(g.TRIANGLES,e.__webglFaceCount,g.UNSIGNED_SHORT,0)}D.info.render.calls++;
D.info.render.vertices=D.info.render.vertices+e.__webglFaceCount;D.info.render.faces=D.info.render.faces+e.__webglFaceCount/3}else if(f instanceof THREE.Line){f=f.type===THREE.LineStrip?g.LINE_STRIP:g.LINES;d=d.linewidth;if(d!==Db){g.lineWidth(d);Db=d}g.drawArrays(f,0,e.__webglLineCount);D.info.render.calls++}}};this.render=function(a,b,c,d){var e,f,j,n,m=a.__lights,o=a.fog;$=-1;rb=true;if(b.parent===void 0){console.warn("DEPRECATED: Camera hasn't been added to a Scene. Adding it...");a.add(b)}this.autoUpdateScene&&
a.updateMatrixWorld();if(!b._viewMatrixArray)b._viewMatrixArray=new Float32Array(16);if(!b._projectionMatrixArray)b._projectionMatrixArray=new Float32Array(16);b.matrixWorldInverse.getInverse(b.matrixWorld);b.matrixWorldInverse.flattenToArray(b._viewMatrixArray);b.projectionMatrix.flattenToArray(b._projectionMatrixArray);qb.multiply(b.projectionMatrix,b.matrixWorldInverse);Yb.setFromMatrix(qb);this.autoUpdateObjects&&this.initWebGLObjects(a);i(this.renderPluginsPre,a,b);D.info.render.calls=0;D.info.render.vertices=
0;D.info.render.faces=0;D.info.render.points=0;this.setRenderTarget(c);(this.autoClear||d)&&this.clear(this.autoClearColor,this.autoClearDepth,this.autoClearStencil);n=a.__webglObjects;d=0;for(e=n.length;d<e;d++){f=n[d];j=f.object;f.render=false;if(j.visible&&(!(j instanceof THREE.Mesh||j instanceof THREE.ParticleSystem)||!j.frustumCulled||Yb.contains(j))){s(j,b);var q=f,u=q.object,y=q.buffer,w=void 0,w=w=void 0,w=u.material;if(w instanceof THREE.MeshFaceMaterial){w=y.materialIndex;if(w>=0){w=u.geometry.materials[w];
if(w.transparent){q.transparent=w;q.opaque=null}else{q.opaque=w;q.transparent=null}}}else if(w)if(w.transparent){q.transparent=w;q.opaque=null}else{q.opaque=w;q.transparent=null}f.render=true;if(this.sortObjects)if(j.renderDepth)f.z=j.renderDepth;else{Ra.copy(j.matrixWorld.getPosition());qb.multiplyVector3(Ra);f.z=Ra.z}}}this.sortObjects&&n.sort(h);n=a.__webglObjectsImmediate;d=0;for(e=n.length;d<e;d++){f=n[d];j=f.object;if(j.visible){s(j,b);j=f.object.material;if(j.transparent){f.transparent=j;f.opaque=
null}else{f.opaque=j;f.transparent=null}}}if(a.overrideMaterial){d=a.overrideMaterial;this.setBlending(d.blending,d.blendEquation,d.blendSrc,d.blendDst);this.setDepthTest(d.depthTest);this.setDepthWrite(d.depthWrite);p(d.polygonOffset,d.polygonOffsetFactor,d.polygonOffsetUnits);l(a.__webglObjects,false,"",b,m,o,true,d);k(a.__webglObjectsImmediate,"",b,m,o,false,d)}else{this.setBlending(THREE.NormalBlending);l(a.__webglObjects,true,"opaque",b,m,o,false);k(a.__webglObjectsImmediate,"opaque",b,m,o,false);
l(a.__webglObjects,false,"transparent",b,m,o,true);k(a.__webglObjectsImmediate,"transparent",b,m,o,true)}i(this.renderPluginsPost,a,b);if(c&&c.generateMipmaps&&c.minFilter!==THREE.NearestFilter&&c.minFilter!==THREE.LinearFilter)if(c instanceof THREE.WebGLRenderTargetCube){g.bindTexture(g.TEXTURE_CUBE_MAP,c.__webglTexture);g.generateMipmap(g.TEXTURE_CUBE_MAP);g.bindTexture(g.TEXTURE_CUBE_MAP,null)}else{g.bindTexture(g.TEXTURE_2D,c.__webglTexture);g.generateMipmap(g.TEXTURE_2D);g.bindTexture(g.TEXTURE_2D,
null)}this.setDepthTest(true);this.setDepthWrite(true)};this.renderImmediateObject=function(a,b,c,d,e){var f=o(a,b,c,d,e);R=-1;D.setObjectFaces(e);e.immediateRenderCallback?e.immediateRenderCallback(f,g,Yb):e.render(function(a){D.renderBufferImmediate(a,f,d.shading)})};this.initWebGLObjects=function(a){if(!a.__webglObjects){a.__webglObjects=[];a.__webglObjectsImmediate=[];a.__webglSprites=[];a.__webglFlares=[]}for(;a.__objectsAdded.length;){var h=a.__objectsAdded[0],i=a,k=void 0,l=void 0,o=void 0;
if(!h.__webglInit){h.__webglInit=true;h._modelViewMatrix=new THREE.Matrix4;h._normalMatrix=new THREE.Matrix3;if(h instanceof THREE.Mesh){l=h.geometry;if(l instanceof THREE.Geometry){if(l.geometryGroups===void 0){var q=l,p=void 0,s=void 0,u=void 0,w=void 0,y=void 0,A=void 0,B=void 0,H={},K=q.morphTargets.length,Y=q.morphNormals.length;q.geometryGroups={};p=0;for(s=q.faces.length;p<s;p++){u=q.faces[p];w=u.materialIndex;A=w!==void 0?w:-1;H[A]===void 0&&(H[A]={hash:A,counter:0});B=H[A].hash+"_"+H[A].counter;
q.geometryGroups[B]===void 0&&(q.geometryGroups[B]={faces3:[],faces4:[],materialIndex:w,vertices:0,numMorphTargets:K,numMorphNormals:Y});y=u instanceof THREE.Face3?3:4;if(q.geometryGroups[B].vertices+y>65535){H[A].counter=H[A].counter+1;B=H[A].hash+"_"+H[A].counter;q.geometryGroups[B]===void 0&&(q.geometryGroups[B]={faces3:[],faces4:[],materialIndex:w,vertices:0,numMorphTargets:K,numMorphNormals:Y})}u instanceof THREE.Face3?q.geometryGroups[B].faces3.push(p):q.geometryGroups[B].faces4.push(p);q.geometryGroups[B].vertices=
q.geometryGroups[B].vertices+y}q.geometryGroupsList=[];var J=void 0;for(J in q.geometryGroups){q.geometryGroups[J].id=Z++;q.geometryGroupsList.push(q.geometryGroups[J])}}for(k in l.geometryGroups){o=l.geometryGroups[k];if(!o.__webglVertexBuffer){var N=o;N.__webglVertexBuffer=g.createBuffer();N.__webglNormalBuffer=g.createBuffer();N.__webglTangentBuffer=g.createBuffer();N.__webglColorBuffer=g.createBuffer();N.__webglUVBuffer=g.createBuffer();N.__webglUV2Buffer=g.createBuffer();N.__webglSkinVertexABuffer=
g.createBuffer();N.__webglSkinVertexBBuffer=g.createBuffer();N.__webglSkinIndicesBuffer=g.createBuffer();N.__webglSkinWeightsBuffer=g.createBuffer();N.__webglFaceBuffer=g.createBuffer();N.__webglLineBuffer=g.createBuffer();var O=void 0,R=void 0;if(N.numMorphTargets){N.__webglMorphTargetsBuffers=[];O=0;for(R=N.numMorphTargets;O<R;O++)N.__webglMorphTargetsBuffers.push(g.createBuffer())}if(N.numMorphNormals){N.__webglMorphNormalsBuffers=[];O=0;for(R=N.numMorphNormals;O<R;O++)N.__webglMorphNormalsBuffers.push(g.createBuffer())}D.info.memory.geometries++;
var I=o,ba=h,ja=ba.geometry,ca=I.faces3,$=I.faces4,Q=ca.length*3+$.length*4,za=ca.length*1+$.length*2,pa=ca.length*3+$.length*4,ya=c(ba,I),Na=e(ya),la=d(ya),Da=ya.vertexColors?ya.vertexColors:false;I.__vertexArray=new Float32Array(Q*3);if(la)I.__normalArray=new Float32Array(Q*3);if(ja.hasTangents)I.__tangentArray=new Float32Array(Q*4);if(Da)I.__colorArray=new Float32Array(Q*3);if(Na){if(ja.faceUvs.length>0||ja.faceVertexUvs.length>0)I.__uvArray=new Float32Array(Q*2);if(ja.faceUvs.length>1||ja.faceVertexUvs.length>
1)I.__uv2Array=new Float32Array(Q*2)}if(ba.geometry.skinWeights.length&&ba.geometry.skinIndices.length){I.__skinVertexAArray=new Float32Array(Q*4);I.__skinVertexBArray=new Float32Array(Q*4);I.__skinIndexArray=new Float32Array(Q*4);I.__skinWeightArray=new Float32Array(Q*4)}I.__faceArray=new Uint16Array(za*3);I.__lineArray=new Uint16Array(pa*2);var sa=void 0,Ca=void 0;if(I.numMorphTargets){I.__morphTargetsArrays=[];sa=0;for(Ca=I.numMorphTargets;sa<Ca;sa++)I.__morphTargetsArrays.push(new Float32Array(Q*
3))}if(I.numMorphNormals){I.__morphNormalsArrays=[];sa=0;for(Ca=I.numMorphNormals;sa<Ca;sa++)I.__morphNormalsArrays.push(new Float32Array(Q*3))}I.__webglFaceCount=za*3;I.__webglLineCount=pa*2;if(ya.attributes){if(I.__webglCustomAttributesList===void 0)I.__webglCustomAttributesList=[];var Ga=void 0;for(Ga in ya.attributes){var La=ya.attributes[Ga],Ma={},Oa;for(Oa in La)Ma[Oa]=La[Oa];if(!Ma.__webglInitialized||Ma.createUniqueBuffers){Ma.__webglInitialized=true;var Ha=1;Ma.type==="v2"?Ha=2:Ma.type===
"v3"?Ha=3:Ma.type==="v4"?Ha=4:Ma.type==="c"&&(Ha=3);Ma.size=Ha;Ma.array=new Float32Array(Q*Ha);Ma.buffer=g.createBuffer();Ma.buffer.belongsToAttribute=Ga;La.needsUpdate=true;Ma.__original=La}I.__webglCustomAttributesList.push(Ma)}}I.__inittedArrays=true;l.verticesNeedUpdate=true;l.morphTargetsNeedUpdate=true;l.elementsNeedUpdate=true;l.uvsNeedUpdate=true;l.normalsNeedUpdate=true;l.tangetsNeedUpdate=true;l.colorsNeedUpdate=true}}}}else if(h instanceof THREE.Line){l=h.geometry;if(!l.__webglVertexBuffer){var Ta=
l;Ta.__webglVertexBuffer=g.createBuffer();Ta.__webglColorBuffer=g.createBuffer();D.info.memory.geometries++;var Ra=l,rb=h,qb=Ra.vertices.length;Ra.__vertexArray=new Float32Array(qb*3);Ra.__colorArray=new Float32Array(qb*3);Ra.__webglLineCount=qb;b(Ra,rb);l.verticesNeedUpdate=true;l.colorsNeedUpdate=true}}else if(h instanceof THREE.ParticleSystem){l=h.geometry;if(!l.__webglVertexBuffer){var Db=l;Db.__webglVertexBuffer=g.createBuffer();Db.__webglColorBuffer=g.createBuffer();D.info.geometries++;var Mb=
l,Ub=h,Eb=Mb.vertices.length;Mb.__vertexArray=new Float32Array(Eb*3);Mb.__colorArray=new Float32Array(Eb*3);Mb.__sortArray=[];Mb.__webglParticleCount=Eb;b(Mb,Ub);l.verticesNeedUpdate=true;l.colorsNeedUpdate=true}}}if(!h.__webglActive){if(h instanceof THREE.Mesh){l=h.geometry;if(l instanceof THREE.BufferGeometry)j(i.__webglObjects,l,h);else for(k in l.geometryGroups){o=l.geometryGroups[k];j(i.__webglObjects,o,h)}}else if(h instanceof THREE.Line){l=h.geometry;j(i.__webglObjects,l,h)}h.__webglActive=
true}a.__objectsAdded.splice(0,1)}for(;a.__objectsRemoved.length;){var mc=a.__objectsRemoved[0];if(mc instanceof THREE.Mesh||mc instanceof THREE.Line)for(var Fb=a.__webglObjects,Xb=mc,nc=Fb.length-1;nc>=0;nc--)Fb[nc].object===Xb&&Fb.splice(nc,1);mc.__webglActive=false;a.__objectsRemoved.splice(0,1)}for(var Gb=0,Yb=a.__webglObjects.length;Gb<Yb;Gb++){var Va=a.__webglObjects[Gb].object,V=Va.geometry,$b=void 0,Nb=void 0,Ea=void 0;if(Va instanceof THREE.Mesh)if(V instanceof THREE.BufferGeometry){V.verticesNeedUpdate=
false;V.elementsNeedUpdate=false;V.uvsNeedUpdate=false;V.normalsNeedUpdate=false;V.colorsNeedUpdate=false}else{for(var zc=0,jc=V.geometryGroupsList.length;zc<jc;zc++){$b=V.geometryGroupsList[zc];Ea=c(Va,$b);Nb=Ea.attributes&&m(Ea);if(V.verticesNeedUpdate||V.morphTargetsNeedUpdate||V.elementsNeedUpdate||V.uvsNeedUpdate||V.normalsNeedUpdate||V.colorsNeedUpdate||V.tangetsNeedUpdate||Nb){var M=$b,kc=Va,Ia=g.DYNAMIC_DRAW,lc=!V.dynamic,Hb=Ea;if(M.__inittedArrays){var Vb=d(Hb),Ac=Hb.vertexColors?Hb.vertexColors:
false,Wb=e(Hb),oc=Vb===THREE.SmoothShading,v=void 0,C=void 0,Sa=void 0,z=void 0,Ob=void 0,sb=void 0,Ua=void 0,pc=void 0,lb=void 0,Pb=void 0,Qb=void 0,E=void 0,F=void 0,G=void 0,W=void 0,Wa=void 0,Xa=void 0,Ya=void 0,ac=void 0,Za=void 0,$a=void 0,ab=void 0,bc=void 0,bb=void 0,cb=void 0,db=void 0,cc=void 0,eb=void 0,fb=void 0,gb=void 0,dc=void 0,hb=void 0,ib=void 0,jb=void 0,ec=void 0,tb=void 0,ub=void 0,vb=void 0,qc=void 0,wb=void 0,xb=void 0,yb=void 0,rc=void 0,S=void 0,Zb=void 0,zb=void 0,Rb=void 0,
Sb=void 0,ta=void 0,Ic=void 0,qa=void 0,ra=void 0,Ab=void 0,mb=void 0,ka=0,oa=0,nb=0,ob=0,Pa=0,xa=0,X=0,Aa=0,ma=0,x=0,L=0,t=0,Ja=void 0,ua=M.__vertexArray,fc=M.__uvArray,gc=M.__uv2Array,Qa=M.__normalArray,da=M.__tangentArray,va=M.__colorArray,ea=M.__skinVertexAArray,fa=M.__skinVertexBArray,ga=M.__skinIndexArray,ha=M.__skinWeightArray,Bc=M.__morphTargetsArrays,Cc=M.__morphNormalsArrays,Dc=M.__webglCustomAttributesList,r=void 0,kb=M.__faceArray,Ka=M.__lineArray,Ba=kc.geometry,yc=Ba.elementsNeedUpdate,
Jc=Ba.uvsNeedUpdate,Nc=Ba.normalsNeedUpdate,Oc=Ba.tangetsNeedUpdate,Pc=Ba.colorsNeedUpdate,Qc=Ba.morphTargetsNeedUpdate,Ib=Ba.vertices,T=M.faces3,U=M.faces4,na=Ba.faces,Ec=Ba.faceVertexUvs[0],Fc=Ba.faceVertexUvs[1],Jb=Ba.skinVerticesA,Kb=Ba.skinVerticesB,Lb=Ba.skinIndices,Bb=Ba.skinWeights,Cb=Ba.morphTargets,sc=Ba.morphNormals;if(Ba.verticesNeedUpdate){v=0;for(C=T.length;v<C;v++){z=na[T[v]];E=Ib[z.a];F=Ib[z.b];G=Ib[z.c];ua[oa]=E.x;ua[oa+1]=E.y;ua[oa+2]=E.z;ua[oa+3]=F.x;ua[oa+4]=F.y;ua[oa+5]=F.z;ua[oa+
6]=G.x;ua[oa+7]=G.y;ua[oa+8]=G.z;oa=oa+9}v=0;for(C=U.length;v<C;v++){z=na[U[v]];E=Ib[z.a];F=Ib[z.b];G=Ib[z.c];W=Ib[z.d];ua[oa]=E.x;ua[oa+1]=E.y;ua[oa+2]=E.z;ua[oa+3]=F.x;ua[oa+4]=F.y;ua[oa+5]=F.z;ua[oa+6]=G.x;ua[oa+7]=G.y;ua[oa+8]=G.z;ua[oa+9]=W.x;ua[oa+10]=W.y;ua[oa+11]=W.z;oa=oa+12}g.bindBuffer(g.ARRAY_BUFFER,M.__webglVertexBuffer);g.bufferData(g.ARRAY_BUFFER,ua,Ia)}if(Qc){ta=0;for(Ic=Cb.length;ta<Ic;ta++){v=L=0;for(C=T.length;v<C;v++){Ab=T[v];z=na[Ab];E=Cb[ta].vertices[z.a];F=Cb[ta].vertices[z.b];
G=Cb[ta].vertices[z.c];qa=Bc[ta];qa[L]=E.x;qa[L+1]=E.y;qa[L+2]=E.z;qa[L+3]=F.x;qa[L+4]=F.y;qa[L+5]=F.z;qa[L+6]=G.x;qa[L+7]=G.y;qa[L+8]=G.z;if(Hb.morphNormals){if(oc){mb=sc[ta].vertexNormals[Ab];Za=mb.a;$a=mb.b;ab=mb.c}else ab=$a=Za=sc[ta].faceNormals[Ab];ra=Cc[ta];ra[L]=Za.x;ra[L+1]=Za.y;ra[L+2]=Za.z;ra[L+3]=$a.x;ra[L+4]=$a.y;ra[L+5]=$a.z;ra[L+6]=ab.x;ra[L+7]=ab.y;ra[L+8]=ab.z}L=L+9}v=0;for(C=U.length;v<C;v++){Ab=U[v];z=na[Ab];E=Cb[ta].vertices[z.a];F=Cb[ta].vertices[z.b];G=Cb[ta].vertices[z.c];W=
Cb[ta].vertices[z.d];qa=Bc[ta];qa[L]=E.x;qa[L+1]=E.y;qa[L+2]=E.z;qa[L+3]=F.x;qa[L+4]=F.y;qa[L+5]=F.z;qa[L+6]=G.x;qa[L+7]=G.y;qa[L+8]=G.z;qa[L+9]=W.x;qa[L+10]=W.y;qa[L+11]=W.z;if(Hb.morphNormals){if(oc){mb=sc[ta].vertexNormals[Ab];Za=mb.a;$a=mb.b;ab=mb.c;bc=mb.d}else bc=ab=$a=Za=sc[ta].faceNormals[Ab];ra=Cc[ta];ra[L]=Za.x;ra[L+1]=Za.y;ra[L+2]=Za.z;ra[L+3]=$a.x;ra[L+4]=$a.y;ra[L+5]=$a.z;ra[L+6]=ab.x;ra[L+7]=ab.y;ra[L+8]=ab.z;ra[L+9]=bc.x;ra[L+10]=bc.y;ra[L+11]=bc.z}L=L+12}g.bindBuffer(g.ARRAY_BUFFER,
M.__webglMorphTargetsBuffers[ta]);g.bufferData(g.ARRAY_BUFFER,Bc[ta],Ia);if(Hb.morphNormals){g.bindBuffer(g.ARRAY_BUFFER,M.__webglMorphNormalsBuffers[ta]);g.bufferData(g.ARRAY_BUFFER,Cc[ta],Ia)}}}if(Bb.length){v=0;for(C=T.length;v<C;v++){z=na[T[v]];eb=Bb[z.a];fb=Bb[z.b];gb=Bb[z.c];ha[x]=eb.x;ha[x+1]=eb.y;ha[x+2]=eb.z;ha[x+3]=eb.w;ha[x+4]=fb.x;ha[x+5]=fb.y;ha[x+6]=fb.z;ha[x+7]=fb.w;ha[x+8]=gb.x;ha[x+9]=gb.y;ha[x+10]=gb.z;ha[x+11]=gb.w;hb=Lb[z.a];ib=Lb[z.b];jb=Lb[z.c];ga[x]=hb.x;ga[x+1]=hb.y;ga[x+2]=
hb.z;ga[x+3]=hb.w;ga[x+4]=ib.x;ga[x+5]=ib.y;ga[x+6]=ib.z;ga[x+7]=ib.w;ga[x+8]=jb.x;ga[x+9]=jb.y;ga[x+10]=jb.z;ga[x+11]=jb.w;tb=Jb[z.a];ub=Jb[z.b];vb=Jb[z.c];ea[x]=tb.x;ea[x+1]=tb.y;ea[x+2]=tb.z;ea[x+3]=1;ea[x+4]=ub.x;ea[x+5]=ub.y;ea[x+6]=ub.z;ea[x+7]=1;ea[x+8]=vb.x;ea[x+9]=vb.y;ea[x+10]=vb.z;ea[x+11]=1;wb=Kb[z.a];xb=Kb[z.b];yb=Kb[z.c];fa[x]=wb.x;fa[x+1]=wb.y;fa[x+2]=wb.z;fa[x+3]=1;fa[x+4]=xb.x;fa[x+5]=xb.y;fa[x+6]=xb.z;fa[x+7]=1;fa[x+8]=yb.x;fa[x+9]=yb.y;fa[x+10]=yb.z;fa[x+11]=1;x=x+12}v=0;for(C=
U.length;v<C;v++){z=na[U[v]];eb=Bb[z.a];fb=Bb[z.b];gb=Bb[z.c];dc=Bb[z.d];ha[x]=eb.x;ha[x+1]=eb.y;ha[x+2]=eb.z;ha[x+3]=eb.w;ha[x+4]=fb.x;ha[x+5]=fb.y;ha[x+6]=fb.z;ha[x+7]=fb.w;ha[x+8]=gb.x;ha[x+9]=gb.y;ha[x+10]=gb.z;ha[x+11]=gb.w;ha[x+12]=dc.x;ha[x+13]=dc.y;ha[x+14]=dc.z;ha[x+15]=dc.w;hb=Lb[z.a];ib=Lb[z.b];jb=Lb[z.c];ec=Lb[z.d];ga[x]=hb.x;ga[x+1]=hb.y;ga[x+2]=hb.z;ga[x+3]=hb.w;ga[x+4]=ib.x;ga[x+5]=ib.y;ga[x+6]=ib.z;ga[x+7]=ib.w;ga[x+8]=jb.x;ga[x+9]=jb.y;ga[x+10]=jb.z;ga[x+11]=jb.w;ga[x+12]=ec.x;ga[x+
13]=ec.y;ga[x+14]=ec.z;ga[x+15]=ec.w;tb=Jb[z.a];ub=Jb[z.b];vb=Jb[z.c];qc=Jb[z.d];ea[x]=tb.x;ea[x+1]=tb.y;ea[x+2]=tb.z;ea[x+3]=1;ea[x+4]=ub.x;ea[x+5]=ub.y;ea[x+6]=ub.z;ea[x+7]=1;ea[x+8]=vb.x;ea[x+9]=vb.y;ea[x+10]=vb.z;ea[x+11]=1;ea[x+12]=qc.x;ea[x+13]=qc.y;ea[x+14]=qc.z;ea[x+15]=1;wb=Kb[z.a];xb=Kb[z.b];yb=Kb[z.c];rc=Kb[z.d];fa[x]=wb.x;fa[x+1]=wb.y;fa[x+2]=wb.z;fa[x+3]=1;fa[x+4]=xb.x;fa[x+5]=xb.y;fa[x+6]=xb.z;fa[x+7]=1;fa[x+8]=yb.x;fa[x+9]=yb.y;fa[x+10]=yb.z;fa[x+11]=1;fa[x+12]=rc.x;fa[x+13]=rc.y;fa[x+
14]=rc.z;fa[x+15]=1;x=x+16}if(x>0){g.bindBuffer(g.ARRAY_BUFFER,M.__webglSkinVertexABuffer);g.bufferData(g.ARRAY_BUFFER,ea,Ia);g.bindBuffer(g.ARRAY_BUFFER,M.__webglSkinVertexBBuffer);g.bufferData(g.ARRAY_BUFFER,fa,Ia);g.bindBuffer(g.ARRAY_BUFFER,M.__webglSkinIndicesBuffer);g.bufferData(g.ARRAY_BUFFER,ga,Ia);g.bindBuffer(g.ARRAY_BUFFER,M.__webglSkinWeightsBuffer);g.bufferData(g.ARRAY_BUFFER,ha,Ia)}}if(Pc&&Ac){v=0;for(C=T.length;v<C;v++){z=na[T[v]];Ua=z.vertexColors;pc=z.color;if(Ua.length===3&&Ac===
THREE.VertexColors){bb=Ua[0];cb=Ua[1];db=Ua[2]}else db=cb=bb=pc;va[ma]=bb.r;va[ma+1]=bb.g;va[ma+2]=bb.b;va[ma+3]=cb.r;va[ma+4]=cb.g;va[ma+5]=cb.b;va[ma+6]=db.r;va[ma+7]=db.g;va[ma+8]=db.b;ma=ma+9}v=0;for(C=U.length;v<C;v++){z=na[U[v]];Ua=z.vertexColors;pc=z.color;if(Ua.length===4&&Ac===THREE.VertexColors){bb=Ua[0];cb=Ua[1];db=Ua[2];cc=Ua[3]}else cc=db=cb=bb=pc;va[ma]=bb.r;va[ma+1]=bb.g;va[ma+2]=bb.b;va[ma+3]=cb.r;va[ma+4]=cb.g;va[ma+5]=cb.b;va[ma+6]=db.r;va[ma+7]=db.g;va[ma+8]=db.b;va[ma+9]=cc.r;
va[ma+10]=cc.g;va[ma+11]=cc.b;ma=ma+12}if(ma>0){g.bindBuffer(g.ARRAY_BUFFER,M.__webglColorBuffer);g.bufferData(g.ARRAY_BUFFER,va,Ia)}}if(Oc&&Ba.hasTangents){v=0;for(C=T.length;v<C;v++){z=na[T[v]];lb=z.vertexTangents;Wa=lb[0];Xa=lb[1];Ya=lb[2];da[X]=Wa.x;da[X+1]=Wa.y;da[X+2]=Wa.z;da[X+3]=Wa.w;da[X+4]=Xa.x;da[X+5]=Xa.y;da[X+6]=Xa.z;da[X+7]=Xa.w;da[X+8]=Ya.x;da[X+9]=Ya.y;da[X+10]=Ya.z;da[X+11]=Ya.w;X=X+12}v=0;for(C=U.length;v<C;v++){z=na[U[v]];lb=z.vertexTangents;Wa=lb[0];Xa=lb[1];Ya=lb[2];ac=lb[3];
da[X]=Wa.x;da[X+1]=Wa.y;da[X+2]=Wa.z;da[X+3]=Wa.w;da[X+4]=Xa.x;da[X+5]=Xa.y;da[X+6]=Xa.z;da[X+7]=Xa.w;da[X+8]=Ya.x;da[X+9]=Ya.y;da[X+10]=Ya.z;da[X+11]=Ya.w;da[X+12]=ac.x;da[X+13]=ac.y;da[X+14]=ac.z;da[X+15]=ac.w;X=X+16}g.bindBuffer(g.ARRAY_BUFFER,M.__webglTangentBuffer);g.bufferData(g.ARRAY_BUFFER,da,Ia)}if(Nc&&Vb){v=0;for(C=T.length;v<C;v++){z=na[T[v]];Ob=z.vertexNormals;sb=z.normal;if(Ob.length===3&&oc)for(S=0;S<3;S++){zb=Ob[S];Qa[xa]=zb.x;Qa[xa+1]=zb.y;Qa[xa+2]=zb.z;xa=xa+3}else for(S=0;S<3;S++){Qa[xa]=
sb.x;Qa[xa+1]=sb.y;Qa[xa+2]=sb.z;xa=xa+3}}v=0;for(C=U.length;v<C;v++){z=na[U[v]];Ob=z.vertexNormals;sb=z.normal;if(Ob.length===4&&oc)for(S=0;S<4;S++){zb=Ob[S];Qa[xa]=zb.x;Qa[xa+1]=zb.y;Qa[xa+2]=zb.z;xa=xa+3}else for(S=0;S<4;S++){Qa[xa]=sb.x;Qa[xa+1]=sb.y;Qa[xa+2]=sb.z;xa=xa+3}}g.bindBuffer(g.ARRAY_BUFFER,M.__webglNormalBuffer);g.bufferData(g.ARRAY_BUFFER,Qa,Ia)}if(Jc&&Ec&&Wb){v=0;for(C=T.length;v<C;v++){Sa=T[v];z=na[Sa];Pb=Ec[Sa];if(Pb!==void 0)for(S=0;S<3;S++){Rb=Pb[S];fc[nb]=Rb.u;fc[nb+1]=Rb.v;
nb=nb+2}}v=0;for(C=U.length;v<C;v++){Sa=U[v];z=na[Sa];Pb=Ec[Sa];if(Pb!==void 0)for(S=0;S<4;S++){Rb=Pb[S];fc[nb]=Rb.u;fc[nb+1]=Rb.v;nb=nb+2}}if(nb>0){g.bindBuffer(g.ARRAY_BUFFER,M.__webglUVBuffer);g.bufferData(g.ARRAY_BUFFER,fc,Ia)}}if(Jc&&Fc&&Wb){v=0;for(C=T.length;v<C;v++){Sa=T[v];z=na[Sa];Qb=Fc[Sa];if(Qb!==void 0)for(S=0;S<3;S++){Sb=Qb[S];gc[ob]=Sb.u;gc[ob+1]=Sb.v;ob=ob+2}}v=0;for(C=U.length;v<C;v++){Sa=U[v];z=na[Sa];Qb=Fc[Sa];if(Qb!==void 0)for(S=0;S<4;S++){Sb=Qb[S];gc[ob]=Sb.u;gc[ob+1]=Sb.v;ob=
ob+2}}if(ob>0){g.bindBuffer(g.ARRAY_BUFFER,M.__webglUV2Buffer);g.bufferData(g.ARRAY_BUFFER,gc,Ia)}}if(yc){v=0;for(C=T.length;v<C;v++){z=na[T[v]];kb[Pa]=ka;kb[Pa+1]=ka+1;kb[Pa+2]=ka+2;Pa=Pa+3;Ka[Aa]=ka;Ka[Aa+1]=ka+1;Ka[Aa+2]=ka;Ka[Aa+3]=ka+2;Ka[Aa+4]=ka+1;Ka[Aa+5]=ka+2;Aa=Aa+6;ka=ka+3}v=0;for(C=U.length;v<C;v++){z=na[U[v]];kb[Pa]=ka;kb[Pa+1]=ka+1;kb[Pa+2]=ka+3;kb[Pa+3]=ka+1;kb[Pa+4]=ka+2;kb[Pa+5]=ka+3;Pa=Pa+6;Ka[Aa]=ka;Ka[Aa+1]=ka+1;Ka[Aa+2]=ka;Ka[Aa+3]=ka+3;Ka[Aa+4]=ka+1;Ka[Aa+5]=ka+2;Ka[Aa+6]=ka+
2;Ka[Aa+7]=ka+3;Aa=Aa+8;ka=ka+4}g.bindBuffer(g.ELEMENT_ARRAY_BUFFER,M.__webglFaceBuffer);g.bufferData(g.ELEMENT_ARRAY_BUFFER,kb,Ia);g.bindBuffer(g.ELEMENT_ARRAY_BUFFER,M.__webglLineBuffer);g.bufferData(g.ELEMENT_ARRAY_BUFFER,Ka,Ia)}if(Dc){S=0;for(Zb=Dc.length;S<Zb;S++){r=Dc[S];if(r.__original.needsUpdate){t=0;if(r.size===1)if(r.boundTo===void 0||r.boundTo==="vertices"){v=0;for(C=T.length;v<C;v++){z=na[T[v]];r.array[t]=r.value[z.a];r.array[t+1]=r.value[z.b];r.array[t+2]=r.value[z.c];t=t+3}v=0;for(C=
U.length;v<C;v++){z=na[U[v]];r.array[t]=r.value[z.a];r.array[t+1]=r.value[z.b];r.array[t+2]=r.value[z.c];r.array[t+3]=r.value[z.d];t=t+4}}else{if(r.boundTo==="faces"){v=0;for(C=T.length;v<C;v++){Ja=r.value[T[v]];r.array[t]=Ja;r.array[t+1]=Ja;r.array[t+2]=Ja;t=t+3}v=0;for(C=U.length;v<C;v++){Ja=r.value[U[v]];r.array[t]=Ja;r.array[t+1]=Ja;r.array[t+2]=Ja;r.array[t+3]=Ja;t=t+4}}}else if(r.size===2)if(r.boundTo===void 0||r.boundTo==="vertices"){v=0;for(C=T.length;v<C;v++){z=na[T[v]];E=r.value[z.a];F=
r.value[z.b];G=r.value[z.c];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=F.x;r.array[t+3]=F.y;r.array[t+4]=G.x;r.array[t+5]=G.y;t=t+6}v=0;for(C=U.length;v<C;v++){z=na[U[v]];E=r.value[z.a];F=r.value[z.b];G=r.value[z.c];W=r.value[z.d];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=F.x;r.array[t+3]=F.y;r.array[t+4]=G.x;r.array[t+5]=G.y;r.array[t+6]=W.x;r.array[t+7]=W.y;t=t+8}}else{if(r.boundTo==="faces"){v=0;for(C=T.length;v<C;v++){G=F=E=Ja=r.value[T[v]];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=F.x;
r.array[t+3]=F.y;r.array[t+4]=G.x;r.array[t+5]=G.y;t=t+6}v=0;for(C=U.length;v<C;v++){W=G=F=E=Ja=r.value[U[v]];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=F.x;r.array[t+3]=F.y;r.array[t+4]=G.x;r.array[t+5]=G.y;r.array[t+6]=W.x;r.array[t+7]=W.y;t=t+8}}}else if(r.size===3){var P;P=r.type==="c"?["r","g","b"]:["x","y","z"];if(r.boundTo===void 0||r.boundTo==="vertices"){v=0;for(C=T.length;v<C;v++){z=na[T[v]];E=r.value[z.a];F=r.value[z.b];G=r.value[z.c];r.array[t]=E[P[0]];r.array[t+1]=E[P[1]];r.array[t+
2]=E[P[2]];r.array[t+3]=F[P[0]];r.array[t+4]=F[P[1]];r.array[t+5]=F[P[2]];r.array[t+6]=G[P[0]];r.array[t+7]=G[P[1]];r.array[t+8]=G[P[2]];t=t+9}v=0;for(C=U.length;v<C;v++){z=na[U[v]];E=r.value[z.a];F=r.value[z.b];G=r.value[z.c];W=r.value[z.d];r.array[t]=E[P[0]];r.array[t+1]=E[P[1]];r.array[t+2]=E[P[2]];r.array[t+3]=F[P[0]];r.array[t+4]=F[P[1]];r.array[t+5]=F[P[2]];r.array[t+6]=G[P[0]];r.array[t+7]=G[P[1]];r.array[t+8]=G[P[2]];r.array[t+9]=W[P[0]];r.array[t+10]=W[P[1]];r.array[t+11]=W[P[2]];t=t+12}}else if(r.boundTo===
"faces"){v=0;for(C=T.length;v<C;v++){G=F=E=Ja=r.value[T[v]];r.array[t]=E[P[0]];r.array[t+1]=E[P[1]];r.array[t+2]=E[P[2]];r.array[t+3]=F[P[0]];r.array[t+4]=F[P[1]];r.array[t+5]=F[P[2]];r.array[t+6]=G[P[0]];r.array[t+7]=G[P[1]];r.array[t+8]=G[P[2]];t=t+9}v=0;for(C=U.length;v<C;v++){W=G=F=E=Ja=r.value[U[v]];r.array[t]=E[P[0]];r.array[t+1]=E[P[1]];r.array[t+2]=E[P[2]];r.array[t+3]=F[P[0]];r.array[t+4]=F[P[1]];r.array[t+5]=F[P[2]];r.array[t+6]=G[P[0]];r.array[t+7]=G[P[1]];r.array[t+8]=G[P[2]];r.array[t+
9]=W[P[0]];r.array[t+10]=W[P[1]];r.array[t+11]=W[P[2]];t=t+12}}}else if(r.size===4)if(r.boundTo===void 0||r.boundTo==="vertices"){v=0;for(C=T.length;v<C;v++){z=na[T[v]];E=r.value[z.a];F=r.value[z.b];G=r.value[z.c];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=E.z;r.array[t+3]=E.w;r.array[t+4]=F.x;r.array[t+5]=F.y;r.array[t+6]=F.z;r.array[t+7]=F.w;r.array[t+8]=G.x;r.array[t+9]=G.y;r.array[t+10]=G.z;r.array[t+11]=G.w;t=t+12}v=0;for(C=U.length;v<C;v++){z=na[U[v]];E=r.value[z.a];F=r.value[z.b];G=r.value[z.c];
W=r.value[z.d];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=E.z;r.array[t+3]=E.w;r.array[t+4]=F.x;r.array[t+5]=F.y;r.array[t+6]=F.z;r.array[t+7]=F.w;r.array[t+8]=G.x;r.array[t+9]=G.y;r.array[t+10]=G.z;r.array[t+11]=G.w;r.array[t+12]=W.x;r.array[t+13]=W.y;r.array[t+14]=W.z;r.array[t+15]=W.w;t=t+16}}else if(r.boundTo==="faces"){v=0;for(C=T.length;v<C;v++){G=F=E=Ja=r.value[T[v]];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=E.z;r.array[t+3]=E.w;r.array[t+4]=F.x;r.array[t+5]=F.y;r.array[t+6]=F.z;r.array[t+
7]=F.w;r.array[t+8]=G.x;r.array[t+9]=G.y;r.array[t+10]=G.z;r.array[t+11]=G.w;t=t+12}v=0;for(C=U.length;v<C;v++){W=G=F=E=Ja=r.value[U[v]];r.array[t]=E.x;r.array[t+1]=E.y;r.array[t+2]=E.z;r.array[t+3]=E.w;r.array[t+4]=F.x;r.array[t+5]=F.y;r.array[t+6]=F.z;r.array[t+7]=F.w;r.array[t+8]=G.x;r.array[t+9]=G.y;r.array[t+10]=G.z;r.array[t+11]=G.w;r.array[t+12]=W.x;r.array[t+13]=W.y;r.array[t+14]=W.z;r.array[t+15]=W.w;t=t+16}}g.bindBuffer(g.ARRAY_BUFFER,r.buffer);g.bufferData(g.ARRAY_BUFFER,r.array,Ia)}}}if(lc){delete M.__inittedArrays;
delete M.__colorArray;delete M.__normalArray;delete M.__tangentArray;delete M.__uvArray;delete M.__uv2Array;delete M.__faceArray;delete M.__vertexArray;delete M.__lineArray;delete M.__skinVertexAArray;delete M.__skinVertexBArray;delete M.__skinIndexArray;delete M.__skinWeightArray}}}}V.verticesNeedUpdate=false;V.morphTargetsNeedUpdate=false;V.elementsNeedUpdate=false;V.uvsNeedUpdate=false;V.normalsNeedUpdate=false;V.colorsNeedUpdate=false;V.tangetsNeedUpdate=false;Ea.attributes&&n(Ea)}else if(Va instanceof
THREE.Line){Ea=c(Va,$b);Nb=Ea.attributes&&m(Ea);if(V.verticesNeedUpdate||V.colorsNeedUpdate||Nb){var pb=V,Gc=g.DYNAMIC_DRAW,hc=void 0,ic=void 0,tc=void 0,ia=void 0,uc=void 0,Kc=pb.vertices,Lc=pb.colors,Rc=Kc.length,Sc=Lc.length,vc=pb.__vertexArray,wc=pb.__colorArray,Tc=pb.colorsNeedUpdate,Hc=pb.__webglCustomAttributesList,xc=void 0,Mc=void 0,wa=void 0,Tb=void 0,Fa=void 0,aa=void 0;if(pb.verticesNeedUpdate){for(hc=0;hc<Rc;hc++){tc=Kc[hc];ia=hc*3;vc[ia]=tc.x;vc[ia+1]=tc.y;vc[ia+2]=tc.z}g.bindBuffer(g.ARRAY_BUFFER,
pb.__webglVertexBuffer);g.bufferData(g.ARRAY_BUFFER,vc,Gc)}if(Tc){for(ic=0;ic<Sc;ic++){uc=Lc[ic];ia=ic*3;wc[ia]=uc.r;wc[ia+1]=uc.g;wc[ia+2]=uc.b}g.bindBuffer(g.ARRAY_BUFFER,pb.__webglColorBuffer);g.bufferData(g.ARRAY_BUFFER,wc,Gc)}if(Hc){xc=0;for(Mc=Hc.length;xc<Mc;xc++){aa=Hc[xc];if(aa.needsUpdate&&(aa.boundTo===void 0||aa.boundTo==="vertices")){ia=0;Tb=aa.value.length;if(aa.size===1)for(wa=0;wa<Tb;wa++)aa.array[wa]=aa.value[wa];else if(aa.size===2)for(wa=0;wa<Tb;wa++){Fa=aa.value[wa];aa.array[ia]=
Fa.x;aa.array[ia+1]=Fa.y;ia=ia+2}else if(aa.size===3)if(aa.type==="c")for(wa=0;wa<Tb;wa++){Fa=aa.value[wa];aa.array[ia]=Fa.r;aa.array[ia+1]=Fa.g;aa.array[ia+2]=Fa.b;ia=ia+3}else for(wa=0;wa<Tb;wa++){Fa=aa.value[wa];aa.array[ia]=Fa.x;aa.array[ia+1]=Fa.y;aa.array[ia+2]=Fa.z;ia=ia+3}else if(aa.size===4)for(wa=0;wa<Tb;wa++){Fa=aa.value[wa];aa.array[ia]=Fa.x;aa.array[ia+1]=Fa.y;aa.array[ia+2]=Fa.z;aa.array[ia+3]=Fa.w;ia=ia+4}g.bindBuffer(g.ARRAY_BUFFER,aa.buffer);g.bufferData(g.ARRAY_BUFFER,aa.array,Gc)}}}}V.verticesNeedUpdate=
false;V.colorsNeedUpdate=false;Ea.attributes&&n(Ea)}else if(Va instanceof THREE.ParticleSystem){Ea=c(Va,$b);Nb=Ea.attributes&&m(Ea);(V.verticesNeedUpdate||V.colorsNeedUpdate||Va.sortParticles||Nb)&&f(V,g.DYNAMIC_DRAW,Va);V.verticesNeedUpdate=false;V.colorsNeedUpdate=false;Ea.attributes&&n(Ea)}}};this.initMaterial=function(a,b,c,d){var e,f,h;a instanceof THREE.MeshDepthMaterial?h="depth":a instanceof THREE.MeshNormalMaterial?h="normal":a instanceof THREE.MeshBasicMaterial?h="basic":a instanceof THREE.MeshLambertMaterial?
h="lambert":a instanceof THREE.MeshPhongMaterial?h="phong":a instanceof THREE.LineBasicMaterial?h="basic":a instanceof THREE.ParticleBasicMaterial&&(h="particle_basic");if(h){var i=THREE.ShaderLib[h];a.uniforms=THREE.UniformsUtils.clone(i.uniforms);a.vertexShader=i.vertexShader;a.fragmentShader=i.fragmentShader}var k,j,l,n,m;k=n=m=i=0;for(j=b.length;k<j;k++){l=b[k];if(!l.onlyShadow){l instanceof THREE.DirectionalLight&&n++;l instanceof THREE.PointLight&&m++;l instanceof THREE.SpotLight&&i++}}if(m+
i+n<=ya){j=n;l=m;n=i}else{j=Math.ceil(ya*n/(m+n));n=l=ya-j}var o=0,i=0;for(m=b.length;i<m;i++){k=b[i];if(k.castShadow){k instanceof THREE.SpotLight&&o++;k instanceof THREE.DirectionalLight&&!k.shadowCascade&&o++}}var p;a:{m=a.fragmentShader;k=a.vertexShader;var i=a.uniforms,b=a.attributes,c={map:!!a.map,envMap:!!a.envMap,lightMap:!!a.lightMap,vertexColors:a.vertexColors,fog:c,useFog:a.fog,sizeAttenuation:a.sizeAttenuation,skinning:a.skinning,maxBones:50,morphTargets:a.morphTargets,morphNormals:a.morphNormals,
maxMorphTargets:this.maxMorphTargets,maxMorphNormals:this.maxMorphNormals,maxDirLights:j,maxPointLights:l,maxSpotLights:n,maxShadows:o,shadowMapEnabled:this.shadowMapEnabled&&d.receiveShadow,shadowMapSoft:this.shadowMapSoft,shadowMapDebug:this.shadowMapDebug,shadowMapCascade:this.shadowMapCascade,alphaTest:a.alphaTest,metal:a.metal,perPixel:a.perPixel,wrapAround:a.wrapAround,doubleSided:d&&d.doubleSided},s,d=[];if(h)d.push(h);else{d.push(m);d.push(k)}for(s in c){d.push(s);d.push(c[s])}h=d.join();
s=0;for(d=Na.length;s<d;s++)if(Na[s].code===h){p=Na[s].program;break a}s=g.createProgram();d=["precision "+B+" float;",Zb>0?"#define VERTEX_TEXTURES":"",D.gammaInput?"#define GAMMA_INPUT":"",D.gammaOutput?"#define GAMMA_OUTPUT":"",D.physicallyBasedShading?"#define PHYSICALLY_BASED_SHADING":"","#define MAX_DIR_LIGHTS "+c.maxDirLights,"#define MAX_POINT_LIGHTS "+c.maxPointLights,"#define MAX_SPOT_LIGHTS "+c.maxSpotLights,"#define MAX_SHADOWS "+c.maxShadows,"#define MAX_BONES "+c.maxBones,c.map?"#define USE_MAP":
"",c.envMap?"#define USE_ENVMAP":"",c.lightMap?"#define USE_LIGHTMAP":"",c.vertexColors?"#define USE_COLOR":"",c.skinning?"#define USE_SKINNING":"",c.morphTargets?"#define USE_MORPHTARGETS":"",c.morphNormals?"#define USE_MORPHNORMALS":"",c.perPixel?"#define PHONG_PER_PIXEL":"",c.wrapAround?"#define WRAP_AROUND":"",c.doubleSided?"#define DOUBLE_SIDED":"",c.shadowMapEnabled?"#define USE_SHADOWMAP":"",c.shadowMapSoft?"#define SHADOWMAP_SOFT":"",c.shadowMapDebug?"#define SHADOWMAP_DEBUG":"",c.shadowMapCascade?
"#define SHADOWMAP_CASCADE":"",c.sizeAttenuation?"#define USE_SIZEATTENUATION":"","uniform mat4 objectMatrix;\nuniform mat4 modelViewMatrix;\nuniform mat4 projectionMatrix;\nuniform mat4 viewMatrix;\nuniform mat3 normalMatrix;\nuniform vec3 cameraPosition;\nattribute vec3 position;\nattribute vec3 normal;\nattribute vec2 uv;\nattribute vec2 uv2;\n#ifdef USE_COLOR\nattribute vec3 color;\n#endif\n#ifdef USE_MORPHTARGETS\nattribute vec3 morphTarget0;\nattribute vec3 morphTarget1;\nattribute vec3 morphTarget2;\nattribute vec3 morphTarget3;\n#ifdef USE_MORPHNORMALS\nattribute vec3 morphNormal0;\nattribute vec3 morphNormal1;\nattribute vec3 morphNormal2;\nattribute vec3 morphNormal3;\n#else\nattribute vec3 morphTarget4;\nattribute vec3 morphTarget5;\nattribute vec3 morphTarget6;\nattribute vec3 morphTarget7;\n#endif\n#endif\n#ifdef USE_SKINNING\nattribute vec4 skinVertexA;\nattribute vec4 skinVertexB;\nattribute vec4 skinIndex;\nattribute vec4 skinWeight;\n#endif\n"].join("\n");
j=["precision "+B+" float;","#define MAX_DIR_LIGHTS "+c.maxDirLights,"#define MAX_POINT_LIGHTS "+c.maxPointLights,"#define MAX_SPOT_LIGHTS "+c.maxSpotLights,"#define MAX_SHADOWS "+c.maxShadows,c.alphaTest?"#define ALPHATEST "+c.alphaTest:"",D.gammaInput?"#define GAMMA_INPUT":"",D.gammaOutput?"#define GAMMA_OUTPUT":"",D.physicallyBasedShading?"#define PHYSICALLY_BASED_SHADING":"",c.useFog&&c.fog?"#define USE_FOG":"",c.useFog&&c.fog instanceof THREE.FogExp2?"#define FOG_EXP2":"",c.map?"#define USE_MAP":
"",c.envMap?"#define USE_ENVMAP":"",c.lightMap?"#define USE_LIGHTMAP":"",c.vertexColors?"#define USE_COLOR":"",c.metal?"#define METAL":"",c.perPixel?"#define PHONG_PER_PIXEL":"",c.wrapAround?"#define WRAP_AROUND":"",c.doubleSided?"#define DOUBLE_SIDED":"",c.shadowMapEnabled?"#define USE_SHADOWMAP":"",c.shadowMapSoft?"#define SHADOWMAP_SOFT":"",c.shadowMapDebug?"#define SHADOWMAP_DEBUG":"",c.shadowMapCascade?"#define SHADOWMAP_CASCADE":"","uniform mat4 viewMatrix;\nuniform vec3 cameraPosition;\n"].join("\n");
g.attachShader(s,q("fragment",j+m));g.attachShader(s,q("vertex",d+k));g.linkProgram(s);g.getProgramParameter(s,g.LINK_STATUS)||console.error("Could not initialise shader\nVALIDATE_STATUS: "+g.getProgramParameter(s,g.VALIDATE_STATUS)+", gl error ["+g.getError()+"]");s.uniforms={};s.attributes={};var u,d=["viewMatrix","modelViewMatrix","projectionMatrix","normalMatrix","objectMatrix","cameraPosition","boneGlobalMatrices","morphTargetInfluences"];for(u in i)d.push(u);u=d;d=0;for(i=u.length;d<i;d++){m=
u[d];s.uniforms[m]=g.getUniformLocation(s,m)}d=["position","normal","uv","uv2","tangent","color","skinVertexA","skinVertexB","skinIndex","skinWeight"];for(u=0;u<c.maxMorphTargets;u++)d.push("morphTarget"+u);for(u=0;u<c.maxMorphNormals;u++)d.push("morphNormal"+u);for(p in b)d.push(p);p=d;u=0;for(b=p.length;u<b;u++){c=p[u];s.attributes[c]=g.getAttribLocation(s,c)}s.id=Na.length;Na.push({program:s,code:h});D.info.memory.programs=Na.length;p=s}a.program=p;p=a.program.attributes;p.position>=0&&g.enableVertexAttribArray(p.position);
p.color>=0&&g.enableVertexAttribArray(p.color);p.normal>=0&&g.enableVertexAttribArray(p.normal);p.tangent>=0&&g.enableVertexAttribArray(p.tangent);if(a.skinning&&p.skinVertexA>=0&&p.skinVertexB>=0&&p.skinIndex>=0&&p.skinWeight>=0){g.enableVertexAttribArray(p.skinVertexA);g.enableVertexAttribArray(p.skinVertexB);g.enableVertexAttribArray(p.skinIndex);g.enableVertexAttribArray(p.skinWeight)}if(a.attributes)for(f in a.attributes)p[f]!==void 0&&p[f]>=0&&g.enableVertexAttribArray(p[f]);if(a.morphTargets){a.numSupportedMorphTargets=
0;s="morphTarget";for(f=0;f<this.maxMorphTargets;f++){u=s+f;if(p[u]>=0){g.enableVertexAttribArray(p[u]);a.numSupportedMorphTargets++}}}if(a.morphNormals){a.numSupportedMorphNormals=0;s="morphNormal";for(f=0;f<this.maxMorphNormals;f++){u=s+f;if(p[u]>=0){g.enableVertexAttribArray(p[u]);a.numSupportedMorphNormals++}}}a.uniformsList=[];for(e in a.uniforms)a.uniformsList.push([a.uniforms[e],e])};this.setFaceCulling=function(a,b){if(a){!b||b==="ccw"?g.frontFace(g.CCW):g.frontFace(g.CW);a==="back"?g.cullFace(g.BACK):
a==="front"?g.cullFace(g.FRONT):g.cullFace(g.FRONT_AND_BACK);g.enable(g.CULL_FACE)}else g.disable(g.CULL_FACE)};this.setObjectFaces=function(a){if(Q!==a.doubleSided){a.doubleSided?g.disable(g.CULL_FACE):g.enable(g.CULL_FACE);Q=a.doubleSided}if(pa!==a.flipSided){a.flipSided?g.frontFace(g.CW):g.frontFace(g.CCW);pa=a.flipSided}};this.setDepthTest=function(a){if(La!==a){a?g.enable(g.DEPTH_TEST):g.disable(g.DEPTH_TEST);La=a}};this.setDepthWrite=function(a){if(Oa!==a){g.depthMask(a);Oa=a}};this.setBlending=
function(a,b,c,d){if(a!==O){switch(a){case THREE.NoBlending:g.disable(g.BLEND);break;case THREE.AdditiveBlending:g.enable(g.BLEND);g.blendEquation(g.FUNC_ADD);g.blendFunc(g.SRC_ALPHA,g.ONE);break;case THREE.SubtractiveBlending:g.enable(g.BLEND);g.blendEquation(g.FUNC_ADD);g.blendFunc(g.ZERO,g.ONE_MINUS_SRC_COLOR);break;case THREE.MultiplyBlending:g.enable(g.BLEND);g.blendEquation(g.FUNC_ADD);g.blendFunc(g.ZERO,g.SRC_COLOR);break;case THREE.CustomBlending:g.enable(g.BLEND);break;default:g.enable(g.BLEND);
g.blendEquationSeparate(g.FUNC_ADD,g.FUNC_ADD);g.blendFuncSeparate(g.SRC_ALPHA,g.ONE_MINUS_SRC_ALPHA,g.ONE,g.ONE_MINUS_SRC_ALPHA)}O=a}if(a===THREE.CustomBlending){if(b!==sa){g.blendEquation(u(b));sa=b}if(c!==Ga||d!==Ha){g.blendFunc(u(c),u(d));Ga=c;Ha=d}}else Ha=Ga=sa=null};this.setTexture=function(a,b){if(a.needsUpdate){if(!a.__webglInit){a.__webglInit=true;a.__webglTexture=g.createTexture();D.info.memory.textures++}g.activeTexture(g.TEXTURE0+b);g.bindTexture(g.TEXTURE_2D,a.__webglTexture);g.pixelStorei(g.UNPACK_PREMULTIPLY_ALPHA_WEBGL,
a.premultiplyAlpha);var c=a.image,d=(c.width&c.width-1)===0&&(c.height&c.height-1)===0,e=u(a.format),f=u(a.type);w(g.TEXTURE_2D,a,d);a instanceof THREE.DataTexture?g.texImage2D(g.TEXTURE_2D,0,e,c.width,c.height,0,e,f,c.data):g.texImage2D(g.TEXTURE_2D,0,e,e,f,a.image);a.generateMipmaps&&d&&g.generateMipmap(g.TEXTURE_2D);a.needsUpdate=false;if(a.onUpdate)a.onUpdate()}else{g.activeTexture(g.TEXTURE0+b);g.bindTexture(g.TEXTURE_2D,a.__webglTexture)}};this.setRenderTarget=function(a){var b=a instanceof
THREE.WebGLRenderTargetCube;if(a&&!a.__webglFramebuffer){if(a.depthBuffer===void 0)a.depthBuffer=true;if(a.stencilBuffer===void 0)a.stencilBuffer=true;a.__webglTexture=g.createTexture();var c=(a.width&a.width-1)===0&&(a.height&a.height-1)===0,d=u(a.format),e=u(a.type);if(b){a.__webglFramebuffer=[];a.__webglRenderbuffer=[];g.bindTexture(g.TEXTURE_CUBE_MAP,a.__webglTexture);w(g.TEXTURE_CUBE_MAP,a,c);for(var f=0;f<6;f++){a.__webglFramebuffer[f]=g.createFramebuffer();a.__webglRenderbuffer[f]=g.createRenderbuffer();
g.texImage2D(g.TEXTURE_CUBE_MAP_POSITIVE_X+f,0,d,a.width,a.height,0,d,e,null);var h=a,i=g.TEXTURE_CUBE_MAP_POSITIVE_X+f;g.bindFramebuffer(g.FRAMEBUFFER,a.__webglFramebuffer[f]);g.framebufferTexture2D(g.FRAMEBUFFER,g.COLOR_ATTACHMENT0,i,h.__webglTexture,0);A(a.__webglRenderbuffer[f],a)}c&&g.generateMipmap(g.TEXTURE_CUBE_MAP)}else{a.__webglFramebuffer=g.createFramebuffer();a.__webglRenderbuffer=g.createRenderbuffer();g.bindTexture(g.TEXTURE_2D,a.__webglTexture);w(g.TEXTURE_2D,a,c);g.texImage2D(g.TEXTURE_2D,
0,d,a.width,a.height,0,d,e,null);d=g.TEXTURE_2D;g.bindFramebuffer(g.FRAMEBUFFER,a.__webglFramebuffer);g.framebufferTexture2D(g.FRAMEBUFFER,g.COLOR_ATTACHMENT0,d,a.__webglTexture,0);A(a.__webglRenderbuffer,a);c&&g.generateMipmap(g.TEXTURE_2D)}b?g.bindTexture(g.TEXTURE_CUBE_MAP,null):g.bindTexture(g.TEXTURE_2D,null);g.bindRenderbuffer(g.RENDERBUFFER,null);g.bindFramebuffer(g.FRAMEBUFFER,null)}if(a){b=b?a.__webglFramebuffer[a.activeCubeFace]:a.__webglFramebuffer;c=a.width;a=a.height;e=d=0}else{b=null;
c=Xb;a=Gb;d=Eb;e=Fb}if(b!==Da){g.bindFramebuffer(g.FRAMEBUFFER,b);g.viewport(d,e,c,a);Da=b}jc=c;kc=a}};
THREE.WebGLRenderTarget=function(a,b,c){this.width=a;this.height=b;c=c||{};this.wrapS=c.wrapS!==void 0?c.wrapS:THREE.ClampToEdgeWrapping;this.wrapT=c.wrapT!==void 0?c.wrapT:THREE.ClampToEdgeWrapping;this.magFilter=c.magFilter!==void 0?c.magFilter:THREE.LinearFilter;this.minFilter=c.minFilter!==void 0?c.minFilter:THREE.LinearMipMapLinearFilter;this.offset=new THREE.Vector2(0,0);this.repeat=new THREE.Vector2(1,1);this.format=c.format!==void 0?c.format:THREE.RGBAFormat;this.type=c.type!==void 0?c.type:
THREE.UnsignedByteType;this.depthBuffer=c.depthBuffer!==void 0?c.depthBuffer:true;this.stencilBuffer=c.stencilBuffer!==void 0?c.stencilBuffer:true;this.generateMipmaps=true};
THREE.WebGLRenderTarget.prototype.clone=function(){var a=new THREE.WebGLRenderTarget(this.width,this.height);a.wrapS=this.wrapS;a.wrapT=this.wrapT;a.magFilter=this.magFilter;a.minFilter=this.minFilter;a.offset.copy(this.offset);a.repeat.copy(this.repeat);a.format=this.format;a.type=this.type;a.depthBuffer=this.depthBuffer;a.stencilBuffer=this.stencilBuffer;return a};THREE.WebGLRenderTargetCube=function(a,b,c){THREE.WebGLRenderTarget.call(this,a,b,c);this.activeCubeFace=0};
THREE.WebGLRenderTargetCube.prototype=new THREE.WebGLRenderTarget;THREE.WebGLRenderTargetCube.prototype.constructor=THREE.WebGLRenderTargetCube;THREE.RenderableVertex=function(){this.positionWorld=new THREE.Vector3;this.positionScreen=new THREE.Vector4;this.visible=true};THREE.RenderableVertex.prototype.copy=function(a){this.positionWorld.copy(a.positionWorld);this.positionScreen.copy(a.positionScreen)};
THREE.RenderableFace3=function(){this.v1=new THREE.RenderableVertex;this.v2=new THREE.RenderableVertex;this.v3=new THREE.RenderableVertex;this.centroidWorld=new THREE.Vector3;this.centroidScreen=new THREE.Vector3;this.normalWorld=new THREE.Vector3;this.vertexNormalsWorld=[new THREE.Vector3,new THREE.Vector3,new THREE.Vector3];this.faceMaterial=this.material=null;this.uvs=[[]];this.z=null};
THREE.RenderableFace4=function(){this.v1=new THREE.RenderableVertex;this.v2=new THREE.RenderableVertex;this.v3=new THREE.RenderableVertex;this.v4=new THREE.RenderableVertex;this.centroidWorld=new THREE.Vector3;this.centroidScreen=new THREE.Vector3;this.normalWorld=new THREE.Vector3;this.vertexNormalsWorld=[new THREE.Vector3,new THREE.Vector3,new THREE.Vector3,new THREE.Vector3];this.faceMaterial=this.material=null;this.uvs=[[]];this.z=null};THREE.RenderableObject=function(){this.z=this.object=null};
THREE.RenderableLine=function(){this.z=null;this.v1=new THREE.RenderableVertex;this.v2=new THREE.RenderableVertex;this.material=null};
THREE.GeometryUtils={merge:function(a,b){for(var c,d,e=a.vertices.length,f=b instanceof THREE.Mesh?b.geometry:b,h=a.vertices,i=f.vertices,l=a.faces,k=f.faces,j=a.faceVertexUvs[0],m=f.faceVertexUvs[0],n={},o=0;o<a.materials.length;o++)n[a.materials[o].id]=o;if(b instanceof THREE.Mesh){b.matrixAutoUpdate&&b.updateMatrix();c=b.matrix;d=new THREE.Matrix4;d.extractRotation(c,b.scale)}for(var o=0,s=i.length;o<s;o++){var p=i[o].clone();c&&c.multiplyVector3(p);h.push(p)}o=0;for(s=k.length;o<s;o++){var h=
k[o],q,w,A=h.vertexNormals,y=h.vertexColors;h instanceof THREE.Face3?q=new THREE.Face3(h.a+e,h.b+e,h.c+e):h instanceof THREE.Face4&&(q=new THREE.Face4(h.a+e,h.b+e,h.c+e,h.d+e));q.normal.copy(h.normal);d&&d.multiplyVector3(q.normal);i=0;for(p=A.length;i<p;i++){w=A[i].clone();d&&d.multiplyVector3(w);q.vertexNormals.push(w)}q.color.copy(h.color);i=0;for(p=y.length;i<p;i++){w=y[i];q.vertexColors.push(w.clone())}if(h.materialIndex!==void 0){i=f.materials[h.materialIndex];p=i.id;y=n[p];if(y===void 0){y=
a.materials.length;n[p]=y;a.materials.push(i)}q.materialIndex=y}q.centroid.copy(h.centroid);c&&c.multiplyVector3(q.centroid);l.push(q)}o=0;for(s=m.length;o<s;o++){c=m[o];d=[];i=0;for(p=c.length;i<p;i++)d.push(new THREE.UV(c[i].u,c[i].v));j.push(d)}},clone:function(a){var b=new THREE.Geometry,c,d=a.vertices,e=a.faces,f=a.faceVertexUvs[0];if(a.materials)b.materials=a.materials.slice();a=0;for(c=d.length;a<c;a++)b.vertices.push(d[a].clone());a=0;for(c=e.length;a<c;a++)b.faces.push(e[a].clone());a=0;
for(c=f.length;a<c;a++){for(var d=f[a],e=[],h=0,i=d.length;h<i;h++)e.push(new THREE.UV(d[h].u,d[h].v));b.faceVertexUvs[0].push(e)}return b},randomPointInTriangle:function(a,b,c){var d,e,f,h=new THREE.Vector3,i=THREE.GeometryUtils.__v1;d=THREE.GeometryUtils.random();e=THREE.GeometryUtils.random();if(d+e>1){d=1-d;e=1-e}f=1-d-e;h.copy(a);h.multiplyScalar(d);i.copy(b);i.multiplyScalar(e);h.addSelf(i);i.copy(c);i.multiplyScalar(f);h.addSelf(i);return h},randomPointInFace:function(a,b,c){var d,e,f;if(a instanceof
THREE.Face3){d=b.vertices[a.a];e=b.vertices[a.b];f=b.vertices[a.c];return THREE.GeometryUtils.randomPointInTriangle(d,e,f)}if(a instanceof THREE.Face4){d=b.vertices[a.a];e=b.vertices[a.b];f=b.vertices[a.c];var b=b.vertices[a.d],h;if(c)if(a._area1&&a._area2){c=a._area1;h=a._area2}else{c=THREE.GeometryUtils.triangleArea(d,e,b);h=THREE.GeometryUtils.triangleArea(e,f,b);a._area1=c;a._area2=h}else{c=THREE.GeometryUtils.triangleArea(d,e,b);h=THREE.GeometryUtils.triangleArea(e,f,b)}return THREE.GeometryUtils.random()*
(c+h)<c?THREE.GeometryUtils.randomPointInTriangle(d,e,b):THREE.GeometryUtils.randomPointInTriangle(e,f,b)}},randomPointsInGeometry:function(a,b){function c(a){function b(c,d){if(d<c)return c;var e=c+Math.floor((d-c)/2);return k[e]>a?b(c,e-1):k[e]<a?b(e+1,d):e}return b(0,k.length-1)}var d,e,f=a.faces,h=a.vertices,i=f.length,l=0,k=[],j,m,n,o;for(e=0;e<i;e++){d=f[e];if(d instanceof THREE.Face3){j=h[d.a];m=h[d.b];n=h[d.c];d._area=THREE.GeometryUtils.triangleArea(j,m,n)}else if(d instanceof THREE.Face4){j=
h[d.a];m=h[d.b];n=h[d.c];o=h[d.d];d._area1=THREE.GeometryUtils.triangleArea(j,m,o);d._area2=THREE.GeometryUtils.triangleArea(m,n,o);d._area=d._area1+d._area2}l=l+d._area;k[e]=l}d=[];for(e=0;e<b;e++){h=THREE.GeometryUtils.random()*l;h=c(h);d[e]=THREE.GeometryUtils.randomPointInFace(f[h],a,true)}return d},triangleArea:function(a,b,c){var d,e=THREE.GeometryUtils.__v1;e.sub(a,b);d=e.length();e.sub(a,c);a=e.length();e.sub(b,c);c=e.length();b=0.5*(d+a+c);return Math.sqrt(b*(b-d)*(b-a)*(b-c))},center:function(a){a.computeBoundingBox();
var b=a.boundingBox,c=new THREE.Vector3;c.add(b.min,b.max);c.multiplyScalar(-0.5);a.applyMatrix((new THREE.Matrix4).makeTranslation(c.x,c.y,c.z));a.computeBoundingBox();return c},normalizeUVs:function(a){for(var a=a.faceVertexUvs[0],b=0,c=a.length;b<c;b++)for(var d=a[b],e=0,f=d.length;e<f;e++){if(d[e].u!==1)d[e].u=d[e].u-Math.floor(d[e].u);if(d[e].v!==1)d[e].v=d[e].v-Math.floor(d[e].v)}},triangulateQuads:function(a){var b,c,d,e,f=[],h=[],i=[];b=0;for(c=a.faceUvs.length;b<c;b++)h[b]=[];b=0;for(c=a.faceVertexUvs.length;b<
c;b++)i[b]=[];b=0;for(c=a.faces.length;b<c;b++){d=a.faces[b];if(d instanceof THREE.Face4){e=d.a;var l=d.b,k=d.c,j=d.d,m=new THREE.Face3,n=new THREE.Face3;m.color.copy(d.color);n.color.copy(d.color);m.materialIndex=d.materialIndex;n.materialIndex=d.materialIndex;m.a=e;m.b=l;m.c=j;n.a=l;n.b=k;n.c=j;if(d.vertexColors.length===4){m.vertexColors[0]=d.vertexColors[0].clone();m.vertexColors[1]=d.vertexColors[1].clone();m.vertexColors[2]=d.vertexColors[3].clone();n.vertexColors[0]=d.vertexColors[1].clone();
n.vertexColors[1]=d.vertexColors[2].clone();n.vertexColors[2]=d.vertexColors[3].clone()}f.push(m,n);d=0;for(e=a.faceVertexUvs.length;d<e;d++)if(a.faceVertexUvs[d].length){m=a.faceVertexUvs[d][b];l=m[1];k=m[2];j=m[3];m=[m[0].clone(),l.clone(),j.clone()];l=[l.clone(),k.clone(),j.clone()];i[d].push(m,l)}d=0;for(e=a.faceUvs.length;d<e;d++)if(a.faceUvs[d].length){l=a.faceUvs[d][b];h[d].push(l,l)}}else{f.push(d);d=0;for(e=a.faceUvs.length;d<e;d++)h[d].push(a.faceUvs[d]);d=0;for(e=a.faceVertexUvs.length;d<
e;d++)i[d].push(a.faceVertexUvs[d])}}a.faces=f;a.faceUvs=h;a.faceVertexUvs=i;a.computeCentroids();a.computeFaceNormals();a.computeVertexNormals();a.hasTangents&&a.computeTangents()},explode:function(a){for(var b=[],c=0,d=a.faces.length;c<d;c++){var e=b.length,f=a.faces[c];if(f instanceof THREE.Face4){var h=f.a,i=f.b,l=f.c,h=a.vertices[h],i=a.vertices[i],l=a.vertices[l],k=a.vertices[f.d];b.push(h.clone());b.push(i.clone());b.push(l.clone());b.push(k.clone());f.a=e;f.b=e+1;f.c=e+2;f.d=e+3}else{h=f.a;
i=f.b;l=f.c;h=a.vertices[h];i=a.vertices[i];l=a.vertices[l];b.push(h.clone());b.push(i.clone());b.push(l.clone());f.a=e;f.b=e+1;f.c=e+2}}a.vertices=b;delete a.__tmpVertices},tessellate:function(a,b){var c,d,e,f,h,i,l,k,j,m,n,o,s,p,q,w,A,y,u,H=[],B=[];c=0;for(d=a.faceVertexUvs.length;c<d;c++)B[c]=[];c=0;for(d=a.faces.length;c<d;c++){e=a.faces[c];if(e instanceof THREE.Face3){f=e.a;h=e.b;i=e.c;k=a.vertices[f];j=a.vertices[h];m=a.vertices[i];o=k.distanceTo(j);s=j.distanceTo(m);n=k.distanceTo(m);if(o>
b||s>b||n>b){l=a.vertices.length;y=e.clone();u=e.clone();if(o>=s&&o>=n){k=k.clone();k.lerpSelf(j,0.5);y.a=f;y.b=l;y.c=i;u.a=l;u.b=h;u.c=i;if(e.vertexNormals.length===3){f=e.vertexNormals[0].clone();f.lerpSelf(e.vertexNormals[1],0.5);y.vertexNormals[1].copy(f);u.vertexNormals[0].copy(f)}if(e.vertexColors.length===3){f=e.vertexColors[0].clone();f.lerpSelf(e.vertexColors[1],0.5);y.vertexColors[1].copy(f);u.vertexColors[0].copy(f)}e=0}else if(s>=o&&s>=n){k=j.clone();k.lerpSelf(m,0.5);y.a=f;y.b=h;y.c=
l;u.a=l;u.b=i;u.c=f;if(e.vertexNormals.length===3){f=e.vertexNormals[1].clone();f.lerpSelf(e.vertexNormals[2],0.5);y.vertexNormals[2].copy(f);u.vertexNormals[0].copy(f);u.vertexNormals[1].copy(e.vertexNormals[2]);u.vertexNormals[2].copy(e.vertexNormals[0])}if(e.vertexColors.length===3){f=e.vertexColors[1].clone();f.lerpSelf(e.vertexColors[2],0.5);y.vertexColors[2].copy(f);u.vertexColors[0].copy(f);u.vertexColors[1].copy(e.vertexColors[2]);u.vertexColors[2].copy(e.vertexColors[0])}e=1}else{k=k.clone();
k.lerpSelf(m,0.5);y.a=f;y.b=h;y.c=l;u.a=l;u.b=h;u.c=i;if(e.vertexNormals.length===3){f=e.vertexNormals[0].clone();f.lerpSelf(e.vertexNormals[2],0.5);y.vertexNormals[2].copy(f);u.vertexNormals[0].copy(f)}if(e.vertexColors.length===3){f=e.vertexColors[0].clone();f.lerpSelf(e.vertexColors[2],0.5);y.vertexColors[2].copy(f);u.vertexColors[0].copy(f)}e=2}H.push(y,u);a.vertices.push(k);f=0;for(h=a.faceVertexUvs.length;f<h;f++)if(a.faceVertexUvs[f].length){k=a.faceVertexUvs[f][c];u=k[0];i=k[1];y=k[2];if(e===
0){j=u.clone();j.lerpSelf(i,0.5);k=[u.clone(),j.clone(),y.clone()];i=[j.clone(),i.clone(),y.clone()]}else if(e===1){j=i.clone();j.lerpSelf(y,0.5);k=[u.clone(),i.clone(),j.clone()];i=[j.clone(),y.clone(),u.clone()]}else{j=u.clone();j.lerpSelf(y,0.5);k=[u.clone(),i.clone(),j.clone()];i=[j.clone(),i.clone(),y.clone()]}B[f].push(k,i)}}else{H.push(e);f=0;for(h=a.faceVertexUvs.length;f<h;f++)B[f].push(a.faceVertexUvs[f][c])}}else{f=e.a;h=e.b;i=e.c;l=e.d;k=a.vertices[f];j=a.vertices[h];m=a.vertices[i];n=
a.vertices[l];o=k.distanceTo(j);s=j.distanceTo(m);p=m.distanceTo(n);q=k.distanceTo(n);if(o>b||s>b||p>b||q>b){w=a.vertices.length;A=a.vertices.length+1;y=e.clone();u=e.clone();if(o>=s&&o>=p&&o>=q||p>=s&&p>=o&&p>=q){o=k.clone();o.lerpSelf(j,0.5);j=m.clone();j.lerpSelf(n,0.5);y.a=f;y.b=w;y.c=A;y.d=l;u.a=w;u.b=h;u.c=i;u.d=A;if(e.vertexNormals.length===4){f=e.vertexNormals[0].clone();f.lerpSelf(e.vertexNormals[1],0.5);h=e.vertexNormals[2].clone();h.lerpSelf(e.vertexNormals[3],0.5);y.vertexNormals[1].copy(f);
y.vertexNormals[2].copy(h);u.vertexNormals[0].copy(f);u.vertexNormals[3].copy(h)}if(e.vertexColors.length===4){f=e.vertexColors[0].clone();f.lerpSelf(e.vertexColors[1],0.5);h=e.vertexColors[2].clone();h.lerpSelf(e.vertexColors[3],0.5);y.vertexColors[1].copy(f);y.vertexColors[2].copy(h);u.vertexColors[0].copy(f);u.vertexColors[3].copy(h)}e=0}else{o=j.clone();o.lerpSelf(m,0.5);j=n.clone();j.lerpSelf(k,0.5);y.a=f;y.b=h;y.c=w;y.d=A;u.a=A;u.b=w;u.c=i;u.d=l;if(e.vertexNormals.length===4){f=e.vertexNormals[1].clone();
f.lerpSelf(e.vertexNormals[2],0.5);h=e.vertexNormals[3].clone();h.lerpSelf(e.vertexNormals[0],0.5);y.vertexNormals[2].copy(f);y.vertexNormals[3].copy(h);u.vertexNormals[0].copy(h);u.vertexNormals[1].copy(f)}if(e.vertexColors.length===4){f=e.vertexColors[1].clone();f.lerpSelf(e.vertexColors[2],0.5);h=e.vertexColors[3].clone();h.lerpSelf(e.vertexColors[0],0.5);y.vertexColors[2].copy(f);y.vertexColors[3].copy(h);u.vertexColors[0].copy(h);u.vertexColors[1].copy(f)}e=1}H.push(y,u);a.vertices.push(o,j);
f=0;for(h=a.faceVertexUvs.length;f<h;f++)if(a.faceVertexUvs[f].length){k=a.faceVertexUvs[f][c];u=k[0];i=k[1];y=k[2];k=k[3];if(e===0){j=u.clone();j.lerpSelf(i,0.5);m=y.clone();m.lerpSelf(k,0.5);u=[u.clone(),j.clone(),m.clone(),k.clone()];i=[j.clone(),i.clone(),y.clone(),m.clone()]}else{j=i.clone();j.lerpSelf(y,0.5);m=k.clone();m.lerpSelf(u,0.5);u=[u.clone(),i.clone(),j.clone(),m.clone()];i=[m.clone(),j.clone(),y.clone(),k.clone()]}B[f].push(u,i)}}else{H.push(e);f=0;for(h=a.faceVertexUvs.length;f<h;f++)B[f].push(a.faceVertexUvs[f][c])}}}a.faces=
H;a.faceVertexUvs=B}};THREE.GeometryUtils.random=THREE.Math.random16;THREE.GeometryUtils.__v1=new THREE.Vector3;
THREE.ImageUtils={crossOrigin:"anonymous",loadTexture:function(a,b,c){var d=new Image,e=new THREE.Texture(d,b);d.onload=function(){e.needsUpdate=true;c&&c(this)};d.crossOrigin=this.crossOrigin;d.src=a;return e},loadTextureCube:function(a,b,c){var d,e=[],f=new THREE.Texture(e,b),b=e.loadCount=0;for(d=a.length;b<d;++b){e[b]=new Image;e[b].onload=function(){e.loadCount=e.loadCount+1;if(e.loadCount===6)f.needsUpdate=true;c&&c(this)};e[b].crossOrigin=this.crossOrigin;e[b].src=a[b]}return f},getNormalMap:function(a,
b){var c=function(a){var b=Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);return[a[0]/b,a[1]/b,a[2]/b]},b=b|1,d=a.width,e=a.height,f=document.createElement("canvas");f.width=d;f.height=e;var h=f.getContext("2d");h.drawImage(a,0,0);for(var i=h.getImageData(0,0,d,e).data,l=h.createImageData(d,e),k=l.data,j=0;j<d;j++)for(var m=0;m<e;m++){var n=m-1<0?0:m-1,o=m+1>e-1?e-1:m+1,s=j-1<0?0:j-1,p=j+1>d-1?d-1:j+1,q=[],w=[0,0,i[(m*d+j)*4]/255*b];q.push([-1,0,i[(m*d+s)*4]/255*b]);q.push([-1,-1,i[(n*d+s)*4]/255*b]);q.push([0,
-1,i[(n*d+j)*4]/255*b]);q.push([1,-1,i[(n*d+p)*4]/255*b]);q.push([1,0,i[(m*d+p)*4]/255*b]);q.push([1,1,i[(o*d+p)*4]/255*b]);q.push([0,1,i[(o*d+j)*4]/255*b]);q.push([-1,1,i[(o*d+s)*4]/255*b]);n=[];s=q.length;for(o=0;o<s;o++){var p=q[o],A=q[(o+1)%s],p=[p[0]-w[0],p[1]-w[1],p[2]-w[2]],A=[A[0]-w[0],A[1]-w[1],A[2]-w[2]];n.push(c([p[1]*A[2]-p[2]*A[1],p[2]*A[0]-p[0]*A[2],p[0]*A[1]-p[1]*A[0]]))}q=[0,0,0];for(o=0;o<n.length;o++){q[0]=q[0]+n[o][0];q[1]=q[1]+n[o][1];q[2]=q[2]+n[o][2]}q[0]=q[0]/n.length;q[1]=
q[1]/n.length;q[2]=q[2]/n.length;w=(m*d+j)*4;k[w]=(q[0]+1)/2*255|0;k[w+1]=(q[1]+0.5)*255|0;k[w+2]=q[2]*255|0;k[w+3]=255}h.putImageData(l,0,0);return f},generateDataTexture:function(a,b,c){for(var d=a*b,e=new Uint8Array(3*d),f=Math.floor(c.r*255),h=Math.floor(c.g*255),c=Math.floor(c.b*255),i=0;i<d;i++){e[i*3]=f;e[i*3+1]=h;e[i*3+2]=c}a=new THREE.DataTexture(e,a,b,THREE.RGBFormat);a.needsUpdate=true;return a}};
THREE.SceneUtils={showHierarchy:function(a,b){THREE.SceneUtils.traverseHierarchy(a,function(a){a.visible=b})},traverseHierarchy:function(a,b){var c,d,e=a.children.length;for(d=0;d<e;d++){c=a.children[d];b(c);THREE.SceneUtils.traverseHierarchy(c,b)}},createMultiMaterialObject:function(a,b){var c,d=b.length,e=new THREE.Object3D;for(c=0;c<d;c++){var f=new THREE.Mesh(a,b[c]);e.add(f)}return e},cloneObject:function(a){var b;if(a instanceof THREE.Mesh)b=new THREE.Mesh(a.geometry,a.material);else if(a instanceof
THREE.Line)b=new THREE.Line(a.geometry,a.material,a.type);else if(a instanceof THREE.Ribbon)b=new THREE.Ribbon(a.geometry,a.material);else if(a instanceof THREE.ParticleSystem){b=new THREE.ParticleSystem(a.geometry,a.material);b.sortParticles=a.sortParticles}else if(a instanceof THREE.Particle)b=new THREE.Particle(a.material);else if(a instanceof THREE.Sprite){b=new THREE.Sprite({});b.color.copy(a.color);b.map=a.map;b.blending=a.blending;b.useScreenCoordinates=a.useScreenCoordinates;b.mergeWith3D=
a.mergeWith3D;b.affectedByDistance=a.affectedByDistance;b.scaleByViewport=a.scaleByViewport;b.alignment=a.alignment;b.rotation3d.copy(a.rotation3d);b.rotation=a.rotation;b.opacity=a.opacity;b.uvOffset.copy(a.uvOffset);b.uvScale.copy(a.uvScale)}else a instanceof THREE.LOD?b=new THREE.LOD:a instanceof THREE.Object3D&&(b=new THREE.Object3D);b.name=a.name;b.parent=a.parent;b.up.copy(a.up);b.position.copy(a.position);b.rotation instanceof THREE.Vector3&&b.rotation.copy(a.rotation);b.eulerOrder=a.eulerOrder;
b.scale.copy(a.scale);b.dynamic=a.dynamic;b.doubleSided=a.doubleSided;b.flipSided=a.flipSided;b.renderDepth=a.renderDepth;b.rotationAutoUpdate=a.rotationAutoUpdate;b.matrix.copy(a.matrix);b.matrixWorld.copy(a.matrixWorld);b.matrixRotationWorld.copy(a.matrixRotationWorld);b.matrixAutoUpdate=a.matrixAutoUpdate;b.matrixWorldNeedsUpdate=a.matrixWorldNeedsUpdate;b.quaternion.copy(a.quaternion);b.useQuaternion=a.useQuaternion;b.boundRadius=a.boundRadius;b.boundRadiusScale=a.boundRadiusScale;b.visible=a.visible;
b.castShadow=a.castShadow;b.receiveShadow=a.receiveShadow;b.frustumCulled=a.frustumCulled;for(var c=0;c<a.children.length;c++){var d=THREE.SceneUtils.cloneObject(a.children[c]);b.children[c]=d;d.parent=b}if(a instanceof THREE.LOD)for(c=0;c<a.LODs.length;c++)b.LODs[c]={visibleAtDistance:a.LODs[c].visibleAtDistance,object3D:b.children[c]};return b},detach:function(a,b,c){a.applyMatrix(b.matrixWorld);b.remove(a);c.add(a)},attach:function(a,b,c){var d=new THREE.Matrix4;d.getInverse(c.matrixWorld);a.applyMatrix(d);
b.remove(a);c.add(a)}};
THREE.WebGLRenderer&&(THREE.ShaderUtils={lib:{fresnel:{uniforms:{mRefractionRatio:{type:"f",value:1.02},mFresnelBias:{type:"f",value:0.1},mFresnelPower:{type:"f",value:2},mFresnelScale:{type:"f",value:1},tCube:{type:"t",value:1,texture:null}},fragmentShader:"uniform samplerCube tCube;\nvarying vec3 vReflect;\nvarying vec3 vRefract[3];\nvarying float vReflectionFactor;\nvoid main() {\nvec4 reflectedColor = textureCube( tCube, vec3( -vReflect.x, vReflect.yz ) );\nvec4 refractedColor = vec4( 1.0, 1.0, 1.0, 1.0 );\nrefractedColor.r = textureCube( tCube, vec3( -vRefract[0].x, vRefract[0].yz ) ).r;\nrefractedColor.g = textureCube( tCube, vec3( -vRefract[1].x, vRefract[1].yz ) ).g;\nrefractedColor.b = textureCube( tCube, vec3( -vRefract[2].x, vRefract[2].yz ) ).b;\nrefractedColor.a = 1.0;\ngl_FragColor = mix( refractedColor, reflectedColor, clamp( vReflectionFactor, 0.0, 1.0 ) );\n}",vertexShader:"uniform float mRefractionRatio;\nuniform float mFresnelBias;\nuniform float mFresnelScale;\nuniform float mFresnelPower;\nvarying vec3 vReflect;\nvarying vec3 vRefract[3];\nvarying float vReflectionFactor;\nvoid main() {\nvec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );\nvec4 mPosition = objectMatrix * vec4( position, 1.0 );\nvec3 nWorld = normalize ( mat3( objectMatrix[0].xyz, objectMatrix[1].xyz, objectMatrix[2].xyz ) * normal );\nvec3 I = mPosition.xyz - cameraPosition;\nvReflect = reflect( I, nWorld );\nvRefract[0] = refract( normalize( I ), nWorld, mRefractionRatio );\nvRefract[1] = refract( normalize( I ), nWorld, mRefractionRatio * 0.99 );\nvRefract[2] = refract( normalize( I ), nWorld, mRefractionRatio * 0.98 );\nvReflectionFactor = mFresnelBias + mFresnelScale * pow( 1.0 + dot( normalize( I ), nWorld ), mFresnelPower );\ngl_Position = projectionMatrix * mvPosition;\n}"},
normal:{uniforms:THREE.UniformsUtils.merge([THREE.UniformsLib.fog,THREE.UniformsLib.lights,THREE.UniformsLib.shadowmap,{enableAO:{type:"i",value:0},enableDiffuse:{type:"i",value:0},enableSpecular:{type:"i",value:0},enableReflection:{type:"i",value:0},tDiffuse:{type:"t",value:0,texture:null},tCube:{type:"t",value:1,texture:null},tNormal:{type:"t",value:2,texture:null},tSpecular:{type:"t",value:3,texture:null},tAO:{type:"t",value:4,texture:null},tDisplacement:{type:"t",value:5,texture:null},uNormalScale:{type:"f",
value:1},uDisplacementBias:{type:"f",value:0},uDisplacementScale:{type:"f",value:1},uDiffuseColor:{type:"c",value:new THREE.Color(16777215)},uSpecularColor:{type:"c",value:new THREE.Color(1118481)},uAmbientColor:{type:"c",value:new THREE.Color(16777215)},uShininess:{type:"f",value:30},uOpacity:{type:"f",value:1},uReflectivity:{type:"f",value:0.5},uOffset:{type:"v2",value:new THREE.Vector2(0,0)},uRepeat:{type:"v2",value:new THREE.Vector2(1,1)},wrapRGB:{type:"v3",value:new THREE.Vector3(1,1,1)}}]),
fragmentShader:["uniform vec3 uAmbientColor;\nuniform vec3 uDiffuseColor;\nuniform vec3 uSpecularColor;\nuniform float uShininess;\nuniform float uOpacity;\nuniform bool enableDiffuse;\nuniform bool enableSpecular;\nuniform bool enableAO;\nuniform bool enableReflection;\nuniform sampler2D tDiffuse;\nuniform sampler2D tNormal;\nuniform sampler2D tSpecular;\nuniform sampler2D tAO;\nuniform samplerCube tCube;\nuniform float uNormalScale;\nuniform float uReflectivity;\nvarying vec3 vTangent;\nvarying vec3 vBinormal;\nvarying vec3 vNormal;\nvarying vec2 vUv;\nuniform vec3 ambientLightColor;\n#if MAX_DIR_LIGHTS > 0\nuniform vec3 directionalLightColor[ MAX_DIR_LIGHTS ];\nuniform vec3 directionalLightDirection[ MAX_DIR_LIGHTS ];\n#endif\n#if MAX_POINT_LIGHTS > 0\nuniform vec3 pointLightColor[ MAX_POINT_LIGHTS ];\nvarying vec4 vPointLight[ MAX_POINT_LIGHTS ];\n#endif\n#ifdef WRAP_AROUND\nuniform vec3 wrapRGB;\n#endif\nvarying vec3 vViewPosition;",
THREE.ShaderChunk.shadowmap_pars_fragment,THREE.ShaderChunk.fog_pars_fragment,"void main() {\ngl_FragColor = vec4( vec3( 1.0 ), uOpacity );\nvec3 specularTex = vec3( 1.0 );\nvec3 normalTex = texture2D( tNormal, vUv ).xyz * 2.0 - 1.0;\nnormalTex.xy *= uNormalScale;\nnormalTex = normalize( normalTex );\nif( enableDiffuse ) {\n#ifdef GAMMA_INPUT\nvec4 texelColor = texture2D( tDiffuse, vUv );\ntexelColor.xyz *= texelColor.xyz;\ngl_FragColor = gl_FragColor * texelColor;\n#else\ngl_FragColor = gl_FragColor * texture2D( tDiffuse, vUv );\n#endif\n}\nif( enableAO ) {\n#ifdef GAMMA_INPUT\nvec4 aoColor = texture2D( tAO, vUv );\naoColor.xyz *= aoColor.xyz;\ngl_FragColor.xyz = gl_FragColor.xyz * aoColor.xyz;\n#else\ngl_FragColor.xyz = gl_FragColor.xyz * texture2D( tAO, vUv ).xyz;\n#endif\n}\nif( enableSpecular )\nspecularTex = texture2D( tSpecular, vUv ).xyz;\nmat3 tsb = mat3( normalize( vTangent ), normalize( vBinormal ), normalize( vNormal ) );\nvec3 finalNormal = tsb * normalTex;\nvec3 normal = normalize( finalNormal );\nvec3 viewPosition = normalize( vViewPosition );\n#if MAX_POINT_LIGHTS > 0\nvec3 pointDiffuse = vec3( 0.0 );\nvec3 pointSpecular = vec3( 0.0 );\nfor ( int i = 0; i < MAX_POINT_LIGHTS; i ++ ) {\nvec3 pointVector = normalize( vPointLight[ i ].xyz );\nfloat pointDistance = vPointLight[ i ].w;\n#ifdef WRAP_AROUND\nfloat pointDiffuseWeightFull = max( dot( normal, pointVector ), 0.0 );\nfloat pointDiffuseWeightHalf = max( 0.5 * dot( normal, pointVector ) + 0.5, 0.0 );\nvec3 pointDiffuseWeight = mix( vec3 ( pointDiffuseWeightFull ), vec3( pointDiffuseWeightHalf ), wrapRGB );\n#else\nfloat pointDiffuseWeight = max( dot( normal, pointVector ), 0.0 );\n#endif\npointDiffuse += pointDistance * pointLightColor[ i ] * uDiffuseColor * pointDiffuseWeight;\nvec3 pointHalfVector = normalize( pointVector + viewPosition );\nfloat pointDotNormalHalf = max( dot( normal, pointHalfVector ), 0.0 );\nfloat pointSpecularWeight = specularTex.r * max( pow( pointDotNormalHalf, uShininess ), 0.0 );\n#ifdef PHYSICALLY_BASED_SHADING\nfloat specularNormalization = ( uShininess + 2.0001 ) / 8.0;\nvec3 schlick = uSpecularColor + vec3( 1.0 - uSpecularColor ) * pow( 1.0 - dot( pointVector, pointHalfVector ), 5.0 );\npointSpecular += schlick * pointLightColor[ i ] * pointSpecularWeight * pointDiffuseWeight * pointDistance * specularNormalization;\n#else\npointSpecular += pointDistance * pointLightColor[ i ] * uSpecularColor * pointSpecularWeight * pointDiffuseWeight;\n#endif\n}\n#endif\n#if MAX_DIR_LIGHTS > 0\nvec3 dirDiffuse = vec3( 0.0 );\nvec3 dirSpecular = vec3( 0.0 );\nfor( int i = 0; i < MAX_DIR_LIGHTS; i++ ) {\nvec4 lDirection = viewMatrix * vec4( directionalLightDirection[ i ], 0.0 );\nvec3 dirVector = normalize( lDirection.xyz );\n#ifdef WRAP_AROUND\nfloat directionalLightWeightingFull = max( dot( normal, dirVector ), 0.0 );\nfloat directionalLightWeightingHalf = max( 0.5 * dot( normal, dirVector ) + 0.5, 0.0 );\nvec3 dirDiffuseWeight = mix( vec3( directionalLightWeightingFull ), vec3( directionalLightWeightingHalf ), wrapRGB );\n#else\nfloat dirDiffuseWeight = max( dot( normal, dirVector ), 0.0 );\n#endif\ndirDiffuse += directionalLightColor[ i ] * uDiffuseColor * dirDiffuseWeight;\nvec3 dirHalfVector = normalize( dirVector + viewPosition );\nfloat dirDotNormalHalf = max( dot( normal, dirHalfVector ), 0.0 );\nfloat dirSpecularWeight = specularTex.r * max( pow( dirDotNormalHalf, uShininess ), 0.0 );\n#ifdef PHYSICALLY_BASED_SHADING\nfloat specularNormalization = ( uShininess + 2.0001 ) / 8.0;\nvec3 schlick = uSpecularColor + vec3( 1.0 - uSpecularColor ) * pow( 1.0 - dot( dirVector, dirHalfVector ), 5.0 );\ndirSpecular += schlick * directionalLightColor[ i ] * dirSpecularWeight * dirDiffuseWeight * specularNormalization;\n#else\ndirSpecular += directionalLightColor[ i ] * uSpecularColor * dirSpecularWeight * dirDiffuseWeight;\n#endif\n}\n#endif\nvec3 totalDiffuse = vec3( 0.0 );\nvec3 totalSpecular = vec3( 0.0 );\n#if MAX_DIR_LIGHTS > 0\ntotalDiffuse += dirDiffuse;\ntotalSpecular += dirSpecular;\n#endif\n#if MAX_POINT_LIGHTS > 0\ntotalDiffuse += pointDiffuse;\ntotalSpecular += pointSpecular;\n#endif\ngl_FragColor.xyz = gl_FragColor.xyz * ( totalDiffuse + ambientLightColor * uAmbientColor) + totalSpecular;\nif ( enableReflection ) {\nvec3 wPos = cameraPosition - vViewPosition;\nvec3 vReflect = reflect( normalize( wPos ), normal );\nvec4 cubeColor = textureCube( tCube, vec3( -vReflect.x, vReflect.yz ) );\n#ifdef GAMMA_INPUT\ncubeColor.xyz *= cubeColor.xyz;\n#endif\ngl_FragColor.xyz = mix( gl_FragColor.xyz, cubeColor.xyz, specularTex.r * uReflectivity );\n}",
THREE.ShaderChunk.shadowmap_fragment,THREE.ShaderChunk.linear_to_gamma_fragment,THREE.ShaderChunk.fog_fragment,"}"].join("\n"),vertexShader:["attribute vec4 tangent;\nuniform vec2 uOffset;\nuniform vec2 uRepeat;\n#ifdef VERTEX_TEXTURES\nuniform sampler2D tDisplacement;\nuniform float uDisplacementScale;\nuniform float uDisplacementBias;\n#endif\nvarying vec3 vTangent;\nvarying vec3 vBinormal;\nvarying vec3 vNormal;\nvarying vec2 vUv;\n#if MAX_POINT_LIGHTS > 0\nuniform vec3 pointLightPosition[ MAX_POINT_LIGHTS ];\nuniform float pointLightDistance[ MAX_POINT_LIGHTS ];\nvarying vec4 vPointLight[ MAX_POINT_LIGHTS ];\n#endif\nvarying vec3 vViewPosition;",
THREE.ShaderChunk.shadowmap_pars_vertex,"void main() {\nvec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );\nvViewPosition = -mvPosition.xyz;\nvNormal = normalMatrix * normal;\nvTangent = normalMatrix * tangent.xyz;\nvBinormal = cross( vNormal, vTangent ) * tangent.w;\nvUv = uv * uRepeat + uOffset;\n#if MAX_POINT_LIGHTS > 0\nfor( int i = 0; i < MAX_POINT_LIGHTS; i++ ) {\nvec4 lPosition = viewMatrix * vec4( pointLightPosition[ i ], 1.0 );\nvec3 lVector = lPosition.xyz - mvPosition.xyz;\nfloat lDistance = 1.0;\nif ( pointLightDistance[ i ] > 0.0 )\nlDistance = 1.0 - min( ( length( lVector ) / pointLightDistance[ i ] ), 1.0 );\nlVector = normalize( lVector );\nvPointLight[ i ] = vec4( lVector, lDistance );\n}\n#endif\n#ifdef VERTEX_TEXTURES\nvec3 dv = texture2D( tDisplacement, uv ).xyz;\nfloat df = uDisplacementScale * dv.x + uDisplacementBias;\nvec4 displacedPosition = vec4( normalize( vNormal.xyz ) * df, 0.0 ) + mvPosition;\ngl_Position = projectionMatrix * displacedPosition;\n#else\ngl_Position = projectionMatrix * mvPosition;\n#endif",
THREE.ShaderChunk.shadowmap_vertex,"}"].join("\n")},cube:{uniforms:{tCube:{type:"t",value:1,texture:null},tFlip:{type:"f",value:-1}},vertexShader:"varying vec3 vViewPosition;\nvoid main() {\nvec4 mPosition = objectMatrix * vec4( position, 1.0 );\nvViewPosition = cameraPosition - mPosition.xyz;\ngl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );\n}",fragmentShader:"uniform samplerCube tCube;\nuniform float tFlip;\nvarying vec3 vViewPosition;\nvoid main() {\nvec3 wPos = cameraPosition - vViewPosition;\ngl_FragColor = textureCube( tCube, vec3( tFlip * wPos.x, wPos.yz ) );\n}"}}});
THREE.BufferGeometry=function(){this.id=THREE.GeometryCount++;this.vertexColorArray=this.vertexUvArray=this.vertexNormalArray=this.vertexPositionArray=this.vertexIndexArray=this.vertexColorBuffer=this.vertexUvBuffer=this.vertexNormalBuffer=this.vertexPositionBuffer=this.vertexIndexBuffer=null;this.dynamic=false;this.boundingSphere=this.boundingBox=null;this.morphTargets=[]};THREE.BufferGeometry.prototype={constructor:THREE.BufferGeometry,computeBoundingBox:function(){},computeBoundingSphere:function(){}};
THREE.CubeGeometry=function(a,b,c,d,e,f,h,i){function l(a,b,c,h,i,j,l,n){var m,o=d||1,p=e||1,q=i/2,g=j/2,s=k.vertices.length;if(a==="x"&&b==="y"||a==="y"&&b==="x")m="z";else if(a==="x"&&b==="z"||a==="z"&&b==="x"){m="y";p=f||1}else if(a==="z"&&b==="y"||a==="y"&&b==="z"){m="x";o=f||1}var w=o+1,y=p+1,A=i/o,R=j/p,J=new THREE.Vector3;J[m]=l>0?1:-1;for(i=0;i<y;i++)for(j=0;j<w;j++){var Z=new THREE.Vector3;Z[a]=(j*A-q)*c;Z[b]=(i*R-g)*h;Z[m]=l;k.vertices.push(Z)}for(i=0;i<p;i++)for(j=0;j<o;j++){a=new THREE.Face4(j+
w*i+s,j+w*(i+1)+s,j+1+w*(i+1)+s,j+1+w*i+s);a.normal.copy(J);a.vertexNormals.push(J.clone(),J.clone(),J.clone(),J.clone());a.materialIndex=n;k.faces.push(a);k.faceVertexUvs[0].push([new THREE.UV(j/o,i/p),new THREE.UV(j/o,(i+1)/p),new THREE.UV((j+1)/o,(i+1)/p),new THREE.UV((j+1)/o,i/p)])}}THREE.Geometry.call(this);var k=this,j=a/2,m=b/2,n=c/2,o,s,p,q,w,A;if(h!==void 0){if(h instanceof Array)this.materials=h;else{this.materials=[];for(o=0;o<6;o++)this.materials.push(h)}o=0;q=1;s=2;w=3;p=4;A=5}else this.materials=
[];this.sides={px:true,nx:true,py:true,ny:true,pz:true,nz:true};if(i!=void 0)for(var y in i)this.sides[y]!==void 0&&(this.sides[y]=i[y]);this.sides.px&&l("z","y",-1,-1,c,b,j,o);this.sides.nx&&l("z","y",1,-1,c,b,-j,q);this.sides.py&&l("x","z",1,1,a,c,m,s);this.sides.ny&&l("x","z",1,-1,a,c,-m,w);this.sides.pz&&l("x","y",1,-1,a,b,n,p);this.sides.nz&&l("x","y",-1,-1,a,b,-n,A);this.computeCentroids();this.mergeVertices()};THREE.CubeGeometry.prototype=new THREE.Geometry;
THREE.CubeGeometry.prototype.constructor=THREE.CubeGeometry;
THREE.CylinderGeometry=function(a,b,c,d,e,f){THREE.Geometry.call(this);var a=a!==void 0?a:20,b=b!==void 0?b:20,c=c!==void 0?c:100,h=c/2,d=d||8,e=e||1,i,l,k=[],j=[];for(l=0;l<=e;l++){var m=[],n=[],o=l/e,s=o*(b-a)+a;for(i=0;i<=d;i++){var p=i/d,q=new THREE.Vector3;q.x=s*Math.sin(p*Math.PI*2);q.y=-o*c+h;q.z=s*Math.cos(p*Math.PI*2);this.vertices.push(q);m.push(this.vertices.length-1);n.push(new THREE.UV(p,o))}k.push(m);j.push(n)}c=(b-a)/c;for(i=0;i<d;i++){if(a!==0){m=this.vertices[k[0][i]].clone();n=this.vertices[k[0][i+
1]].clone()}else{m=this.vertices[k[1][i]].clone();n=this.vertices[k[1][i+1]].clone()}m.setY(Math.sqrt(m.x*m.x+m.z*m.z)*c).normalize();n.setY(Math.sqrt(n.x*n.x+n.z*n.z)*c).normalize();for(l=0;l<e;l++){var o=k[l][i],s=k[l+1][i],p=k[l+1][i+1],q=k[l][i+1],w=m.clone(),A=m.clone(),y=n.clone(),u=n.clone(),H=j[l][i].clone(),B=j[l+1][i].clone(),K=j[l+1][i+1].clone(),N=j[l][i+1].clone();this.faces.push(new THREE.Face4(o,s,p,q,[w,A,y,u]));this.faceVertexUvs[0].push([H,B,K,N])}}if(!f&&a>0){this.vertices.push(new THREE.Vector3(0,
h,0));for(i=0;i<d;i++){o=k[0][i];s=k[0][i+1];p=this.vertices.length-1;w=new THREE.Vector3(0,1,0);A=new THREE.Vector3(0,1,0);y=new THREE.Vector3(0,1,0);H=j[0][i].clone();B=j[0][i+1].clone();K=new THREE.UV(B.u,0);this.faces.push(new THREE.Face3(o,s,p,[w,A,y]));this.faceVertexUvs[0].push([H,B,K])}}if(!f&&b>0){this.vertices.push(new THREE.Vector3(0,-h,0));for(i=0;i<d;i++){o=k[l][i+1];s=k[l][i];p=this.vertices.length-1;w=new THREE.Vector3(0,-1,0);A=new THREE.Vector3(0,-1,0);y=new THREE.Vector3(0,-1,0);
H=j[l][i+1].clone();B=j[l][i].clone();K=new THREE.UV(B.u,1);this.faces.push(new THREE.Face3(o,s,p,[w,A,y]));this.faceVertexUvs[0].push([H,B,K])}}this.computeCentroids();this.computeFaceNormals()};THREE.CylinderGeometry.prototype=new THREE.Geometry;THREE.CylinderGeometry.prototype.constructor=THREE.CylinderGeometry;
THREE.PlaneGeometry=function(a,b,c,d){THREE.Geometry.call(this);for(var e=a/2,f=b/2,c=c||1,d=d||1,h=c+1,i=d+1,l=a/c,k=b/d,j=new THREE.Vector3(0,1,0),a=0;a<i;a++)for(b=0;b<h;b++)this.vertices.push(new THREE.Vector3(b*l-e,0,a*k-f));for(a=0;a<d;a++)for(b=0;b<c;b++){e=new THREE.Face4(b+h*a,b+h*(a+1),b+1+h*(a+1),b+1+h*a);e.normal.copy(j);e.vertexNormals.push(j.clone(),j.clone(),j.clone(),j.clone());this.faces.push(e);this.faceVertexUvs[0].push([new THREE.UV(b/c,a/d),new THREE.UV(b/c,(a+1)/d),new THREE.UV((b+
1)/c,(a+1)/d),new THREE.UV((b+1)/c,a/d)])}this.computeCentroids()};THREE.PlaneGeometry.prototype=new THREE.Geometry;THREE.PlaneGeometry.prototype.constructor=THREE.PlaneGeometry;
THREE.SphereGeometry=function(a,b,c,d,e,f,h){THREE.Geometry.call(this);var a=a||50,d=d!==void 0?d:0,e=e!==void 0?e:Math.PI*2,f=f!==void 0?f:0,h=h!==void 0?h:Math.PI,b=Math.max(3,Math.floor(b)||8),c=Math.max(2,Math.floor(c)||6),i,l,k=[],j=[];for(l=0;l<=c;l++){var m=[],n=[];for(i=0;i<=b;i++){var o=i/b,s=l/c,p=new THREE.Vector3;p.x=-a*Math.cos(d+o*e)*Math.sin(f+s*h);p.y=a*Math.cos(f+s*h);p.z=a*Math.sin(d+o*e)*Math.sin(f+s*h);this.vertices.push(p);m.push(this.vertices.length-1);n.push(new THREE.UV(o,
s))}k.push(m);j.push(n)}for(l=0;l<c;l++)for(i=0;i<b;i++){var d=k[l][i+1],e=k[l][i],f=k[l+1][i],h=k[l+1][i+1],m=this.vertices[d].clone().normalize(),n=this.vertices[e].clone().normalize(),o=this.vertices[f].clone().normalize(),s=this.vertices[h].clone().normalize(),p=j[l][i+1].clone(),q=j[l][i].clone(),w=j[l+1][i].clone(),A=j[l+1][i+1].clone();if(Math.abs(this.vertices[d].y)==a){this.faces.push(new THREE.Face3(d,f,h,[m,o,s]));this.faceVertexUvs[0].push([p,w,A])}else if(Math.abs(this.vertices[f].y)==
a){this.faces.push(new THREE.Face3(d,e,f,[m,n,o]));this.faceVertexUvs[0].push([p,q,w])}else{this.faces.push(new THREE.Face4(d,e,f,h,[m,n,o,s]));this.faceVertexUvs[0].push([p,q,w,A])}}this.computeCentroids();this.computeFaceNormals();this.boundingSphere={radius:a}};THREE.SphereGeometry.prototype=new THREE.Geometry;THREE.SphereGeometry.prototype.constructor=THREE.SphereGeometry;
THREE.PolyhedronGeometry=function(a,b,c,d){function e(a){var b=a.normalize().clone();b.index=l.vertices.push(b)-1;var c=Math.atan2(a.z,-a.x)/2/Math.PI+0.5,a=Math.atan2(-a.y,Math.sqrt(a.x*a.x+a.z*a.z))/Math.PI+0.5;b.uv=new THREE.UV(c,a);return b}function f(a,b,c,d){if(d<1){d=new THREE.Face3(a.index,b.index,c.index,[a.clone(),b.clone(),c.clone()]);d.centroid.addSelf(a).addSelf(b).addSelf(c).divideScalar(3);d.normal=d.centroid.clone().normalize();l.faces.push(d);d=Math.atan2(d.centroid.z,-d.centroid.x);
l.faceVertexUvs[0].push([i(a.uv,a,d),i(b.uv,b,d),i(c.uv,c,d)])}else{d=d-1;f(a,h(a,b),h(a,c),d);f(h(a,b),b,h(b,c),d);f(h(a,c),h(b,c),c,d);f(h(a,b),h(b,c),h(a,c),d)}}function h(a,b){m[a.index]||(m[a.index]=[]);m[b.index]||(m[b.index]=[]);var c=m[a.index][b.index];c===void 0&&(m[a.index][b.index]=m[b.index][a.index]=c=e((new THREE.Vector3).add(a,b).divideScalar(2)));return c}function i(a,b,c){c<0&&a.u===1&&(a=new THREE.UV(a.u-1,a.v));b.x===0&&b.z===0&&(a=new THREE.UV(c/2/Math.PI+0.5,a.v));return a}THREE.Geometry.call(this);
for(var c=c||1,d=d||0,l=this,k=0,j=a.length;k<j;k++)e(new THREE.Vector3(a[k][0],a[k][1],a[k][2]));for(var m=[],a=this.vertices,k=0,j=b.length;k<j;k++)f(a[b[k][0]],a[b[k][1]],a[b[k][2]],d);this.mergeVertices();k=0;for(j=this.vertices.length;k<j;k++)this.vertices[k].multiplyScalar(c);this.computeCentroids();this.boundingSphere={radius:c}};THREE.PolyhedronGeometry.prototype=new THREE.Geometry;THREE.PolyhedronGeometry.prototype.constructor=THREE.PolyhedronGeometry;
THREE.IcosahedronGeometry=function(a,b){var c=(1+Math.sqrt(5))/2;THREE.PolyhedronGeometry.call(this,[[-1,c,0],[1,c,0],[-1,-c,0],[1,-c,0],[0,-1,c],[0,1,c],[0,-1,-c],[0,1,-c],[c,0,-1],[c,0,1],[-c,0,-1],[-c,0,1]],[[0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],[1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],[3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],[4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1]],a,b)};THREE.IcosahedronGeometry.prototype=new THREE.Geometry;THREE.IcosahedronGeometry.prototype.constructor=THREE.IcosahedronGeometry;

/*
 GLmol - Molecular Viewer on WebGL/Javascript (0.47)
  (C) Copyright 2011-2012, biochem_fan
      License: dual license of MIT or LGPL3

  Contributors:
    Robert Hanson for parseXYZ, deferred instantiation

  This program uses
      Three.js 
         https://github.com/mrdoob/three.js
         Copyright (c) 2010-2012 three.js Authors. All rights reserved.
      jQuery
         http://jquery.org/
         Copyright (c) 2011 John Resig
 */

// Workaround for Intel GMA series (gl_FrontFacing causes compilation error)
THREE.ShaderLib.lambert.fragmentShader = THREE.ShaderLib.lambert.fragmentShader.replace("gl_FrontFacing", "true");
THREE.ShaderLib.lambert.vertexShader = THREE.ShaderLib.lambert.vertexShader.replace(/\}$/, "#ifdef DOUBLE_SIDED\n if (transformedNormal.z < 0.0) vLightFront = vLightBack;\n #endif\n }");

var TV3 = THREE.Vector3, TF3 = THREE.Face3, TCo = THREE.Color;

THREE.Geometry.prototype.colorAll = function (color) {
   for (var i = 0; i < this.faces.length; i++) {
      this.faces[i].color = color;
   }
};

THREE.Matrix4.prototype.isIdentity = function() {
   for (var i = 0; i < 4; i++)
      for (var j = 0; j < 4; j++) 
         if (this.elements[i * 4 + j] != (i == j) ? 1 : 0) return false;
   return true;
};

var GLmol = (function() {
function GLmol(id, suppressAutoload) {
   if (id) this.create(id, suppressAutoload);
   return true;
}

GLmol.prototype.create = function(id, suppressAutoload) {
   this.Nucleotides = ['  G', '  A', '  T', '  C', '  U', ' DG', ' DA', ' DT', ' DC', ' DU'];
   this.ElementColors = {"H": 0xCCCCCC, "C": 0xAAAAAA, "O": 0xCC0000, "N": 0x0000CC, "S": 0xCCCC00, "P": 0x6622CC,
                         "F": 0x00CC00, "CL": 0x00CC00, "BR": 0x882200, "I": 0x6600AA,
                         "FE": 0xCC6600, "CA": 0x8888AA};
// Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
   this.vdwRadii = {"H": 1.2, "Li": 1.82, "Na": 2.27, "K": 2.75, "C": 1.7, "N": 1.55, "O": 1.52,
                   "F": 1.47, "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "SE": 1.90,
                   "ZN": 1.39, "CU": 1.4, "NI": 1.63};

   this.id = id;
   this.aaScale = 1; // or 2

   this.container = $('#' + this.id);
   this.WIDTH = this.container.width() * this.aaScale, this.HEIGHT = this.container.height() * this.aaScale;
   this.ASPECT = this.WIDTH / this.HEIGHT;
   this.NEAR = 1, FAR = 800;
   this.CAMERA_Z = -150;
   this.renderer = new THREE.WebGLRenderer({antialias: true});
   this.renderer.sortObjects = false; // hopefully improve performance
   // 'antialias: true' now works in Firefox too!
   // setting this.aaScale = 2 will enable antialias in older Firefox but GPU load increases.
   this.renderer.domElement.style.width = "100%";
   this.renderer.domElement.style.height = "100%";
   this.container.append(this.renderer.domElement);
   this.renderer.setSize(this.WIDTH, this.HEIGHT);

   this.camera = new THREE.PerspectiveCamera(20, this.ASPECT, 1, 800); // will be updated anyway
   this.camera.position = new TV3(0, 0, this.CAMERA_Z);
   this.camera.lookAt(new TV3(0, 0, 0));
   this.perspectiveCamera = this.camera;
   this.orthoscopicCamera = new THREE.OrthographicCamera();
   this.orthoscopicCamera.position.z = this.CAMERA_Z;
   this.orthoscopicCamera.lookAt(new TV3(0, 0, 0));

   var self = this;
   $(window).resize(function() { // only window can capture resize event
      self.WIDTH = self.container.width() * self.aaScale;
      self.HEIGHT = self.container.height() * self.aaScale;
      self.ASPECT = self.WIDTH / self.HEIGHT;
      self.renderer.setSize(self.WIDTH, self.HEIGHT);
      self.camera.aspect = self.ASPECT;
      self.camera.updateProjectionMatrix();
      self.show();
   });

   this.scene = null;
   this.rotationGroup = null; // which contains modelGroup
   this.modelGroup = null;

   this.bgColor = 0x000000;
   this.fov = 20;
   this.fogStart = 0.4;
   this.slabNear = -50; // relative to the center of rotationGroup
   this.slabFar = +50;

   // Default values
   this.sphereRadius = 1.5; 
   this.cylinderRadius = 0.4;
   this.lineWidth = 1.5 * this.aaScale;
   this.curveWidth = 3 * this.aaScale;
   this.defaultColor = 0xCCCCCC;
   this.sphereQuality = 16; //16;
   this.cylinderQuality = 16; //8;
   this.axisDIV = 5; // 3 still gives acceptable quality
   this.strandDIV = 6;
   this.nucleicAcidStrandDIV = 4;
   this.tubeDIV = 8;
   this.coilWidth = 0.3;
   this.helixSheetWidth = 1.3;
   this.nucleicAcidWidth = 0.8;
   this.thickness = 0.4;
 
   // UI variables
   this.cq = new THREE.Quaternion(1, 0, 0, 0);
   this.dq = new THREE.Quaternion(1, 0, 0, 0);
   this.isDragging = false;
   this.mouseStartX = 0;
   this.mouseStartY = 0;
   this.currentModelPos = 0;
   this.cz = 0;
   this.enableMouse();

   if (suppressAutoload) return;
   this.loadMolecule();
}

GLmol.prototype.setupLights = function(scene) {
   var directionalLight =  new THREE.DirectionalLight(0xFFFFFF);
   directionalLight.position = new TV3(0.2, 0.2, -1).normalize();
   directionalLight.intensity = 1.2;
   scene.add(directionalLight);
   var ambientLight = new THREE.AmbientLight(0x202020);
   scene.add(ambientLight);
};

GLmol.prototype.parseSDF = function(str) {
   var atoms = this.atoms;
   var protein = this.protein;

   var lines = str.split("\n");
   if (lines.length < 4) return;
   var atomCount = parseInt(lines[3].substr(0, 3));
   if (isNaN(atomCount) || atomCount <= 0) return;
   var bondCount = parseInt(lines[3].substr(3, 3));
   var offset = 4;
   if (lines.length < 4 + atomCount + bondCount) return;
   for (var i = 1; i <= atomCount; i++) {
      var line = lines[offset];
      offset++;
      var atom = {};
      atom.serial = i;
      atom.x = parseFloat(line.substr(0, 10));
      atom.y = parseFloat(line.substr(10, 10));
      atom.z = parseFloat(line.substr(20, 10));
      atom.hetflag = true;
      atom.atom = atom.elem = line.substr(31, 3).replace(/ /g, "");
      atom.bonds = [];
      atom.bondOrder = [];
      atoms[i] = atom;
   }
   for (i = 1; i <= bondCount; i++) {
      var line = lines[offset];
      offset++;
      var from = parseInt(line.substr(0, 3));
      var to = parseInt(line.substr(3, 3));
      var order = parseInt(line.substr(6, 3));
      atoms[from].bonds.push(to);
      atoms[from].bondOrder.push(order);
      atoms[to].bonds.push(from);
      atoms[to].bondOrder.push(order);
   }

   protein.smallMolecule = true;
   return true;
};

GLmol.prototype.parseXYZ = function(str) {
   var atoms = this.atoms;
   var protein = this.protein;

   var lines = str.split("\n");
   if (lines.length < 3) return;
   var atomCount = parseInt(lines[0].substr(0, 3));
   if (isNaN(atomCount) || atomCount <= 0) return;
   if (lines.length < atomCount + 2) return;
   var offset = 2;
   for (var i = 1; i <= atomCount; i++) {
      var line = lines[offset++];
      var tokens = line.replace(/^\s+/, "").replace(/\s+/g," ").split(" ");
      console.log(tokens);
      var atom = {};
      atom.serial = i;
      atom.atom = atom.elem = tokens[0];
      atom.x = parseFloat(tokens[1]);
      atom.y = parseFloat(tokens[2]);
      atom.z = parseFloat(tokens[3]);
      atom.hetflag = true;
      atom.bonds = [];
      atom.bondOrder = [];
      atoms[i] = atom;
   }
   for (var i = 1; i < atomCount; i++) // hopefully XYZ is small enough
      for (var j = i + 1; j <= atomCount; j++)
         if (this.isConnected(atoms[i], atoms[j])) {
	    atoms[i].bonds.push(j);
	    atoms[i].bondOrder.push(1);
    	    atoms[j].bonds.push(i);
     	    atoms[j].bondOrder.push(1);
         }
   protein.smallMolecule = true;
   return true;
};

GLmol.prototype.parsePDB2 = function(str) {
   var atoms = this.atoms;
   var protein = this.protein;
   var molID;

   var atoms_cnt = 0;
   lines = str.split("\n");
   for (var i = 0; i < lines.length; i++) {
      line = lines[i].replace(/^\s*/, ''); // remove indent
      var recordName = line.substr(0, 6);
      if (recordName == 'ATOM  ' || recordName == 'HETATM') {
         var atom, resn, chain, resi, x, y, z, hetflag, elem, serial, altLoc, b;
         altLoc = line.substr(16, 1);
         if (altLoc != ' ' && altLoc != 'A') continue; // FIXME: ad hoc
         serial = parseInt(line.substr(6, 5));
         atom = line.substr(12, 4).replace(/ /g, "");
         resn = line.substr(17, 3);
         chain = line.substr(21, 1);
         resi = parseInt(line.substr(22, 5)); 
         x = parseFloat(line.substr(30, 8));
         y = parseFloat(line.substr(38, 8));
         z = parseFloat(line.substr(46, 8));
         b = parseFloat(line.substr(60, 8));
         elem = line.substr(76, 2).replace(/ /g, "");
         if (elem == '') { // for some incorrect PDB files
            elem = line.substr(12, 4).replace(/ /g,"");
         }
         if (line[0] == 'H') hetflag = true;
         else hetflag = false;
         atoms[serial] = {'resn': resn, 'x': x, 'y': y, 'z': z, 'elem': elem,
  'hetflag': hetflag, 'chain': chain, 'resi': resi, 'serial': serial, 'atom': atom,
  'bonds': [], 'ss': 'c', 'color': 0xFFFFFF, 'bonds': [], 'bondOrder': [], 'b': b /*', altLoc': altLoc*/};
      } else if (recordName == 'SHEET ') {
         var startChain = line.substr(21, 1);
         var startResi = parseInt(line.substr(22, 4));
         var endChain = line.substr(32, 1);
         var endResi = parseInt(line.substr(33, 4));
         protein.sheet.push([startChain, startResi, endChain, endResi]);
     } else if (recordName == 'CONECT') {
// MEMO: We don't have to parse SSBOND, LINK because both are also 
// described in CONECT. But what about 2JYT???
         var from = parseInt(line.substr(6, 5));
         for (var j = 0; j < 4; j++) {
            var to = parseInt(line.substr([11, 16, 21, 26][j], 5));
            if (isNaN(to)) continue;
            if (atoms[from] != undefined) {
               atoms[from].bonds.push(to);
               atoms[from].bondOrder.push(1);
            }
         }
     } else if (recordName == 'HELIX ') {
         var startChain = line.substr(19, 1);
         var startResi = parseInt(line.substr(21, 4));
         var endChain = line.substr(31, 1);
         var endResi = parseInt(line.substr(33, 4));
         protein.helix.push([startChain, startResi, endChain, endResi]);
     } else if (recordName == 'CRYST1') {
         protein.a = parseFloat(line.substr(6, 9));
         protein.b = parseFloat(line.substr(15, 9));
         protein.c = parseFloat(line.substr(24, 9));
         protein.alpha = parseFloat(line.substr(33, 7));
         protein.beta = parseFloat(line.substr(40, 7));
         protein.gamma = parseFloat(line.substr(47, 7));
         protein.spacegroup = line.substr(55, 11);
         this.defineCell();
      } else if (recordName == 'REMARK') {
         var type = parseInt(line.substr(7, 3));
         if (type == 290 && line.substr(13, 5) == 'SMTRY') {
            var n = parseInt(line[18]) - 1;
            var m = parseInt(line.substr(21, 2));
            if (protein.symMat[m] == undefined) protein.symMat[m] = new THREE.Matrix4().identity();
            protein.symMat[m].elements[n] = parseFloat(line.substr(24, 9));
            protein.symMat[m].elements[n + 4] = parseFloat(line.substr(34, 9));
            protein.symMat[m].elements[n + 8] = parseFloat(line.substr(44, 9));
            protein.symMat[m].elements[n + 12] = parseFloat(line.substr(54, 10));
         } else if (type == 350 && line.substr(13, 5) == 'BIOMT') {
            var n = parseInt(line[18]) - 1;
            var m = parseInt(line.substr(21, 2));
            if (protein.biomtMatrices[m] == undefined) protein.biomtMatrices[m] = new THREE.Matrix4().identity();
            protein.biomtMatrices[m].elements[n] = parseFloat(line.substr(24, 9));
            protein.biomtMatrices[m].elements[n + 4] = parseFloat(line.substr(34, 9));
            protein.biomtMatrices[m].elements[n + 8] = parseFloat(line.substr(44, 9));
            protein.biomtMatrices[m].elements[n + 12] = parseFloat(line.substr(54, 10));
         } else if (type == 350 && line.substr(11, 11) == 'BIOMOLECULE') {
             protein.biomtMatrices = []; protein.biomtChains = '';
         } else if (type == 350 && line.substr(34, 6) == 'CHAINS') {
             protein.biomtChains += line.substr(41, 40);
         }
      } else if (recordName == 'HEADER') {
         protein.pdbID = line.substr(62, 4);
      } else if (recordName == 'TITLE ') {
         if (protein.title == undefined) protein.title = "";
            protein.title += line.substr(10, 70) + "\n"; // CHECK: why 60 is not enough???
      } else if (recordName == 'COMPND') {
              // TODO: Implement me!
      }
   }

   // Assign secondary structures 
   for (i = 0; i < atoms.length; i++) {
      atom = atoms[i]; if (atom == undefined) continue;

      var found = false;
      // MEMO: Can start chain and end chain differ?
      for (j = 0; j < protein.sheet.length; j++) {
         if (atom.chain != protein.sheet[j][0]) continue;
         if (atom.resi < protein.sheet[j][1]) continue;
         if (atom.resi > protein.sheet[j][3]) continue;
         atom.ss = 's';
         if (atom.resi == protein.sheet[j][1]) atom.ssbegin = true;
         if (atom.resi == protein.sheet[j][3]) atom.ssend = true;
      }
      for (j = 0; j < protein.helix.length; j++) {
         if (atom.chain != protein.helix[j][0]) continue;
         if (atom.resi < protein.helix[j][1]) continue;
         if (atom.resi > protein.helix[j][3]) continue;
         atom.ss = 'h';
         if (atom.resi == protein.helix[j][1]) atom.ssbegin = true;
         else if (atom.resi == protein.helix[j][3]) atom.ssend = true;
      }
   }
   protein.smallMolecule = false;
   return true;
};

// Catmull-Rom subdivision
GLmol.prototype.subdivide = function(_points, DIV) { // points as Vector3
   var ret = [];
   var points = _points;
   points = new Array(); // Smoothing test
   points.push(_points[0]);
   for (var i = 1, lim = _points.length - 1; i < lim; i++) {
      var p1 = _points[i], p2 = _points[i + 1];
      if (p1.smoothen) points.push(new TV3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2));
      else points.push(p1);
   }
   points.push(_points[_points.length - 1]);

   for (var i = -1, size = points.length; i <= size - 3; i++) {
      var p0 = points[(i == -1) ? 0 : i];
      var p1 = points[i + 1], p2 = points[i + 2];
      var p3 = points[(i == size - 3) ? size - 1 : i + 3];
      var v0 = new TV3().sub(p2, p0).multiplyScalar(0.5);
      var v1 = new TV3().sub(p3, p1).multiplyScalar(0.5);
      for (var j = 0; j < DIV; j++) {
         var t = 1.0 / DIV * j;
         var x = p1.x + t * v0.x 
                  + t * t * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x)
                  + t * t * t * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
         var y = p1.y + t * v0.y 
                  + t * t * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y)
                  + t * t * t * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
         var z = p1.z + t * v0.z 
                  + t * t * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z)
                  + t * t * t * (2 * p1.z - 2 * p2.z + v0.z + v1.z);
         ret.push(new TV3(x, y, z));
      }
   }
   ret.push(points[points.length - 1]);
   return ret;
};

GLmol.prototype.drawAtomsAsSphere = function(group, atomlist, defaultRadius, forceDefault, scale) {
   var sphereGeometry = new THREE.SphereGeometry(1, this.sphereQuality, this.sphereQuality); // r, seg, ring
   for (var i = 0; i < atomlist.length; i++) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      var sphereMaterial = new THREE.MeshLambertMaterial({color: atom.color});
      var sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
      group.add(sphere);
      var r = (!forceDefault && this.vdwRadii[atom.elem] != undefined) ? this.vdwRadii[atom.elem] : defaultRadius;
      if (!forceDefault && scale) r *= scale;
      sphere.scale.x = sphere.scale.y = sphere.scale.z = r;
      sphere.position.x = atom.x;
      sphere.position.y = atom.y;
      sphere.position.z = atom.z;
   }
};

// about two times faster than sphere when div = 2
GLmol.prototype.drawAtomsAsIcosahedron = function(group, atomlist, defaultRadius, forceDefault) {
   var geo = this.IcosahedronGeometry();
   for (var i = 0; i < atomlist.length; i++) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      var mat = new THREE.MeshLambertMaterial({color: atom.color});
      var sphere = new THREE.Mesh(geo, mat);
      sphere.scale.x = sphere.scale.y = sphere.scale.z = (!forceDefault && this.vdwRadii[atom.elem] != undefined) ? this.vdwRadii[atom.elem] : defaultRadius;
      group.add(sphere);
      sphere.position.x = atom.x;
      sphere.position.y = atom.y;
      sphere.position.z = atom.z;
   }
};

GLmol.prototype.isConnected = function(atom1, atom2) {
   var s = atom1.bonds.indexOf(atom2.serial);
   if (s != -1) return atom1.bondOrder[s];

   if (this.protein.smallMolecule && (atom1.hetflag || atom2.hetflag)) return 0; // CHECK: or should I ?

   var distSquared = (atom1.x - atom2.x) * (atom1.x - atom2.x) + 
                     (atom1.y - atom2.y) * (atom1.y - atom2.y) + 
                     (atom1.z - atom2.z) * (atom1.z - atom2.z);

//   if (atom1.altLoc != atom2.altLoc) return false;
   if (isNaN(distSquared)) return 0;
   if (distSquared < 0.5) return 0; // maybe duplicate position.

   if (distSquared > 1.3 && (atom1.elem == 'H' || atom2.elem == 'H' || atom1.elem == 'D' || atom2.elem == 'D')) return 0;
   if (distSquared < 3.42 && (atom1.elem == 'S' || atom2.elem == 'S')) return 1;
   if (distSquared > 2.78) return 0;
   return 1;
};

GLmol.prototype.drawBondAsStickSub = function(group, atom1, atom2, bondR, order) {
   var delta, tmp;
   if (order > 1) delta = this.calcBondDelta(atom1, atom2, bondR * 2.3);
   var p1 = new TV3(atom1.x, atom1.y, atom1.z);
   var p2 = new TV3(atom2.x, atom2.y, atom2.z);
   var mp = p1.clone().addSelf(p2).multiplyScalar(0.5);

   var c1 = new TCo(atom1.color), c2 = new TCo(atom2.color);
   if (order == 1 || order == 3) {
      this.drawCylinder(group, p1, mp, bondR, atom1.color);
      this.drawCylinder(group, p2, mp, bondR, atom2.color);
   }
   if (order > 1) {
      tmp = mp.clone().addSelf(delta);
      this.drawCylinder(group, p1.clone().addSelf(delta), tmp, bondR, atom1.color);
      this.drawCylinder(group, p2.clone().addSelf(delta), tmp, bondR, atom2.color);
      tmp = mp.clone().subSelf(delta);
      this.drawCylinder(group, p1.clone().subSelf(delta), tmp, bondR, atom1.color);
      this.drawCylinder(group, p2.clone().subSelf(delta), tmp, bondR, atom2.color);
   }
};

GLmol.prototype.drawBondsAsStick = function(group, atomlist, bondR, atomR, ignoreNonbonded, multipleBonds, scale) {
   var sphereGeometry = new THREE.SphereGeometry(1, this.sphereQuality, this.sphereQuality);
   var nAtoms = atomlist.length, mp;
   var forSpheres = [];
   if (!!multipleBonds) bondR /= 2.5;
   for (var _i = 0; _i < nAtoms; _i++) {
      var i = atomlist[_i];
      var atom1 = this.atoms[i];
      if (atom1 == undefined) continue;
      for (var _j = _i + 1; _j < _i + 30 && _j < nAtoms; _j++) {
         var j = atomlist[_j];
         var atom2 = this.atoms[j];
         if (atom2 == undefined) continue;
         var order = this.isConnected(atom1, atom2);
         if (order == 0) continue;
         atom1.connected = atom2.connected = true;
         this.drawBondAsStickSub(group, atom1, atom2, bondR, (!!multipleBonds) ? order : 1);
      }
      for (var _j = 0; _j < atom1.bonds.length; _j++) {
         var j = atom1.bonds[_j];
         if (j < i + 30) continue; // be conservative!
         if (atomlist.indexOf(j) == -1) continue;
         var atom2 = this.atoms[j];
         if (atom2 == undefined) continue;
         atom1.connected = atom2.connected = true;
         this.drawBondAsStickSub(group, atom1, atom2, bondR, (!!multipleBonds) ? atom1.bondOrder[_j] : 1);
      }
      if (atom1.connected) forSpheres.push(i);
   }
   this.drawAtomsAsSphere(group, forSpheres, atomR, !scale, scale);
};

GLmol.prototype.defineCell = function() {
    var p = this.protein;
    if (p.a == undefined) return;

    p.ax = p.a;
    p.ay = 0;
    p.az = 0;
    p.bx = p.b * Math.cos(Math.PI / 180.0 * p.gamma);
    p.by = p.b * Math.sin(Math.PI / 180.0 * p.gamma);
    p.bz = 0;
    p.cx = p.c * Math.cos(Math.PI / 180.0 * p.beta);
    p.cy = p.c * (Math.cos(Math.PI / 180.0 * p.alpha) - 
               Math.cos(Math.PI / 180.0 * p.gamma) 
             * Math.cos(Math.PI / 180.0 * p.beta)
             / Math.sin(Math.PI / 180.0 * p.gamma));
    p.cz = Math.sqrt(p.c * p.c * Math.sin(Math.PI / 180.0 * p.beta)
               * Math.sin(Math.PI / 180.0 * p.beta) - p.cy * p.cy);
};

GLmol.prototype.drawUnitcell = function(group) {
    var p = this.protein;
    if (p.a == undefined) return;

    var vertices = [[0, 0, 0], [p.ax, p.ay, p.az], [p.bx, p.by, p.bz], [p.ax + p.bx, p.ay + p.by, p.az + p.bz],
          [p.cx, p.cy, p.cz], [p.cx + p.ax, p.cy + p.ay,  p.cz + p.az], [p.cx + p.bx, p.cy + p.by, p.cz + p.bz], [p.cx + p.ax + p.bx, p.cy + p.ay + p.by, p.cz + p.az + p.bz]];
    var edges = [0, 1, 0, 2, 1, 3, 2, 3, 4, 5, 4, 6, 5, 7, 6, 7, 0, 4, 1, 5, 2, 6, 3, 7];    

    var geo = new THREE.Geometry();
    for (var i = 0; i < edges.length; i++) {
       geo.vertices.push(new TV3(vertices[edges[i]][0], vertices[edges[i]][1], vertices[edges[i]][2]));
    }
   var lineMaterial = new THREE.LineBasicMaterial({linewidth: 1, color: 0xcccccc});
   var line = new THREE.Line(geo, lineMaterial);
   line.type = THREE.LinePieces;
   group.add(line);
};

// TODO: Find inner side of a ring
GLmol.prototype.calcBondDelta = function(atom1, atom2, sep) {
   var dot;
   var axis = new TV3(atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z).normalize();
   var found = null;
   for (var i = 0; i < atom1.bonds.length && !found; i++) {
      var atom = this.atoms[atom1.bonds[i]]; if (!atom) continue;
      if (atom.serial != atom2.serial && atom.elem != 'H') found = atom;
   }
   for (var i = 0; i < atom2.bonds.length && !found; i++) {
      var atom = this.atoms[atom2.bonds[i]]; if (!atom) continue;
      if (atom.serial != atom1.serial && atom.elem != 'H') found = atom;
   }
   if (found) {
      var tmp = new TV3(atom1.x - found.x, atom1.y - found.y, atom1.z - found.z).normalize();
      dot = tmp.dot(axis);
      delta = new TV3(tmp.x - axis.x * dot, tmp.y - axis.y * dot, tmp.z - axis.z * dot);
   }
   if (!found || Math.abs(dot - 1) < 0.001 || Math.abs(dot + 1) < 0.001) {
      if (axis.x < 0.01 && axis.y < 0.01) {
         delta = new TV3(0, -axis.z, axis.y);
      } else {
         delta = new TV3(-axis.y, axis.x, 0);
      }
   }
   delta.normalize().multiplyScalar(sep);
   return delta;
};

GLmol.prototype.drawBondsAsLineSub = function(geo, atom1, atom2, order) {
   var delta, tmp, vs = geo.vertices, cs = geo.colors;
   if (order > 1) delta = this.calcBondDelta(atom1, atom2, 0.15);
   var p1 = new TV3(atom1.x, atom1.y, atom1.z);
   var p2 = new TV3(atom2.x, atom2.y, atom2.z);
   var mp = p1.clone().addSelf(p2).multiplyScalar(0.5);

   var c1 = new TCo(atom1.color), c2 = new TCo(atom2.color);
   if (order == 1 || order == 3) {
      vs.push(p1); cs.push(c1); vs.push(mp); cs.push(c1);
      vs.push(p2); cs.push(c2); vs.push(mp); cs.push(c2);
   }
   if (order > 1) {
      vs.push(p1.clone().addSelf(delta)); cs.push(c1);
      vs.push(tmp = mp.clone().addSelf(delta)); cs.push(c1);
      vs.push(p2.clone().addSelf(delta)); cs.push(c2);
      vs.push(tmp); cs.push(c2);
      vs.push(p1.clone().subSelf(delta)); cs.push(c1);
      vs.push(tmp = mp.clone().subSelf(delta)); cs.push(c1);
      vs.push(p2.clone().subSelf(delta)); cs.push(c2);
      vs.push(tmp); cs.push(c2);
   }
};

GLmol.prototype.drawBondsAsLine = function(group, atomlist, lineWidth) {
   var geo = new THREE.Geometry();   
   var nAtoms = atomlist.length;

   for (var _i = 0; _i < nAtoms; _i++) {
      var i = atomlist[_i];
      var  atom1 = this.atoms[i];
      if (atom1 == undefined) continue;
      for (var _j = _i + 1; _j < _i + 30 && _j < nAtoms; _j++) {
         var j = atomlist[_j];
         var atom2 = this.atoms[j];
         if (atom2 == undefined) continue;
         var order = this.isConnected(atom1, atom2);
         if (order == 0) continue;

         this.drawBondsAsLineSub(geo, atom1, atom2, order);
      }
      for (var _j = 0; _j < atom1.bonds.length; _j++) {
          var j = atom1.bonds[_j];
          if (j < i + 30) continue; // be conservative!
          if (atomlist.indexOf(j) == -1) continue;
          var atom2 = this.atoms[j];
          if (atom2 == undefined) continue;
          this.drawBondsAsLineSub(geo, atom1, atom2, atom1.bondOrder[_j]);
      }
    }
   var lineMaterial = new THREE.LineBasicMaterial({linewidth: lineWidth});
   lineMaterial.vertexColors = true;

   var line = new THREE.Line(geo, lineMaterial);
   line.type = THREE.LinePieces;
   group.add(line);
};

GLmol.prototype.drawSmoothCurve = function(group, _points, width, colors, div) {
   if (_points.length == 0) return;

   div = (div == undefined) ? 5 : div;

   var geo = new THREE.Geometry();
   var points = this.subdivide(_points, div);

   for (var i = 0; i < points.length; i++) {
      geo.vertices.push(points[i]);
      geo.colors.push(new TCo(colors[(i == 0) ? 0 : Math.round((i - 1) / div)]));
  }
  var lineMaterial = new THREE.LineBasicMaterial({linewidth: width});
  lineMaterial.vertexColors = true;
  var line = new THREE.Line(geo, lineMaterial);
  line.type = THREE.LineStrip;
  group.add(line);
};

GLmol.prototype.drawAsCross = function(group, atomlist, delta) {
   var geo = new THREE.Geometry();
   var points = [[delta, 0, 0], [-delta, 0, 0], [0, delta, 0], [0, -delta, 0], [0, 0, delta], [0, 0, -delta]];
 
   for (var i = 0, lim = atomlist.length; i < lim; i++) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      var c = new TCo(atom.color);
      for (var j = 0; j < 6; j++) {
         geo.vertices.push(new TV3(atom.x + points[j][0], atom.y + points[j][1], atom.z + points[j][2]));
         geo.colors.push(c);
      }
  }
  var lineMaterial = new THREE.LineBasicMaterial({linewidth: this.lineWidth});
  lineMaterial.vertexColors = true;
  var line = new THREE.Line(geo, lineMaterial, THREE.LinePieces);
  group.add(line);
};

// FIXME: Winkled...
GLmol.prototype.drawSmoothTube = function(group, _points, colors, radii) {
   if (_points.length < 2) return;

   var circleDiv = this.tubeDIV, axisDiv = this.axisDIV;
   var geo = new THREE.Geometry();
   var points = this.subdivide(_points, axisDiv);
   var prevAxis1 = new TV3(), prevAxis2;

   for (var i = 0, lim = points.length; i < lim; i++) {
      var r, idx = (i - 1) / axisDiv;
      if (i == 0) r = radii[0];
      else { 
         if (idx % 1 == 0) r = radii[idx];
         else {
            var floored = Math.floor(idx);
            var tmp = idx - floored;
            r = radii[floored] * tmp + radii[floored + 1] * (1 - tmp);
         }
      }
      var delta, axis1, axis2;

      if (i < lim - 1) {
         delta = new TV3().sub(points[i], points[i + 1]);
         axis1 = new TV3(0, - delta.z, delta.y).normalize().multiplyScalar(r);
         axis2 = new TV3().cross(delta, axis1).normalize().multiplyScalar(r);
//      var dir = 1, offset = 0;
         if (prevAxis1.dot(axis1) < 0) {
                 axis1.negate(); axis2.negate();  //dir = -1;//offset = 2 * Math.PI / axisDiv;
         }
         prevAxis1 = axis1; prevAxis2 = axis2;
      } else {
         axis1 = prevAxis1; axis2 = prevAxis2;
      }

      for (var j = 0; j < circleDiv; j++) {
         var angle = 2 * Math.PI / circleDiv * j; //* dir  + offset;
         var c = Math.cos(angle), s = Math.sin(angle);
         geo.vertices.push(new TV3(
         points[i].x + c * axis1.x + s * axis2.x,
         points[i].y + c * axis1.y + s * axis2.y, 
         points[i].z + c * axis1.z + s * axis2.z));
      }
   }

   var offset = 0;
   for (var i = 0, lim = points.length - 1; i < lim; i++) {
      var c =  new TCo(colors[Math.round((i - 1)/ axisDiv)]);

      var reg = 0;
      var r1 = new TV3().sub(geo.vertices[offset], geo.vertices[offset + circleDiv]).lengthSq();
      var r2 = new TV3().sub(geo.vertices[offset], geo.vertices[offset + circleDiv + 1]).lengthSq();
      if (r1 > r2) {r1 = r2; reg = 1;};
      for (var j = 0; j < circleDiv; j++) {
          geo.faces.push(new TF3(offset + j, offset + (j + reg) % circleDiv + circleDiv, offset + (j + 1) % circleDiv));
          geo.faces.push(new TF3(offset + (j + 1) % circleDiv, offset + (j + reg) % circleDiv + circleDiv, offset + (j + reg + 1) % circleDiv + circleDiv));
          geo.faces[geo.faces.length -2].color = c;
          geo.faces[geo.faces.length -1].color = c;
      }
      offset += circleDiv;
   }
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var mat = new THREE.MeshLambertMaterial();
   mat.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, mat);
   mesh.doubleSided = true;
   group.add(mesh);
};


GLmol.prototype.drawMainchainCurve = function(group, atomlist, curveWidth, atomName, div) {
   var points = [], colors = [];
   var currentChain, currentResi;
   if (div == undefined) div = 5;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == atomName) && !atom.hetflag) {
         if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
            this.drawSmoothCurve(group, points, curveWidth, colors, div);
            points = [];
            colors = [];
         }
         points.push(new TV3(atom.x, atom.y, atom.z));
         colors.push(atom.color);
         currentChain = atom.chain;
         currentResi = atom.resi;
      }
   }
    this.drawSmoothCurve(group, points, curveWidth, colors, div);
};

GLmol.prototype.drawMainchainTube = function(group, atomlist, atomName, radius) {
   var points = [], colors = [], radii = [];
   var currentChain, currentResi;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == atomName) && !atom.hetflag) {
         if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
            this.drawSmoothTube(group, points, colors, radii);
            points = []; colors = []; radii = [];
         }
         points.push(new TV3(atom.x, atom.y, atom.z));
         if (radius == undefined) {
            radii.push((atom.b > 0) ? atom.b / 100 : 0.3);
         } else {
            radii.push(radius);
         }
         colors.push(atom.color);
         currentChain = atom.chain;
         currentResi = atom.resi;
      }
   }
   this.drawSmoothTube(group, points, colors, radii);
};

GLmol.prototype.drawStrip = function(group, p1, p2, colors, div, thickness) {
   if ((p1.length) < 2) return;
   div = div || this.axisDIV;
   p1 = this.subdivide(p1, div);
   p2 = this.subdivide(p2, div);
   if (!thickness) return this.drawThinStrip(group, p1, p2, colors, div);

   var geo = new THREE.Geometry();
   var vs = geo.vertices, fs = geo.faces;
   var axis, p1v, p2v, a1v, a2v;
   for (var i = 0, lim = p1.length; i < lim; i++) {
      vs.push(p1v = p1[i]); // 0
      vs.push(p1v); // 1
      vs.push(p2v = p2[i]); // 2
      vs.push(p2v); // 3
      if (i < lim - 1) {
         var toNext = p1[i + 1].clone().subSelf(p1[i]);
         var toSide = p2[i].clone().subSelf(p1[i]);
         axis = toSide.crossSelf(toNext).normalize().multiplyScalar(thickness);
      }
      vs.push(a1v = p1[i].clone().addSelf(axis)); // 4
      vs.push(a1v); // 5
      vs.push(a2v = p2[i].clone().addSelf(axis)); // 6
      vs.push(a2v); // 7
   }
   var faces = [[0, 2, -6, -8], [-4, -2, 6, 4], [7, 3, -5, -1], [-3, -7, 1, 5]];
   for (var i = 1, lim = p1.length; i < lim; i++) {
      var offset = 8 * i, color = new TCo(colors[Math.round((i - 1)/ div)]);
      for (var j = 0; j < 4; j++) {
         var f = new THREE.Face4(offset + faces[j][0], offset + faces[j][1], offset + faces[j][2], offset + faces[j][3], undefined, color);
         fs.push(f);
      }
   }
   var vsize = vs.length - 8; // Cap
   for (var i = 0; i < 4; i++) {vs.push(vs[i * 2]); vs.push(vs[vsize + i * 2])};
   vsize += 8;
   fs.push(new THREE.Face4(vsize, vsize + 2, vsize + 6, vsize + 4, undefined, fs[0].color));
   fs.push(new THREE.Face4(vsize + 1, vsize + 5, vsize + 7, vsize + 3, undefined, fs[fs.length - 3].color));
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var material =  new THREE.MeshLambertMaterial();
   material.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, material);
   mesh.doubleSided = true;
   group.add(mesh);
};


GLmol.prototype.drawThinStrip = function(group, p1, p2, colors, div) {
   var geo = new THREE.Geometry();
   for (var i = 0, lim = p1.length; i < lim; i++) {
      geo.vertices.push(p1[i]); // 2i
      geo.vertices.push(p2[i]); // 2i + 1
   }
   for (var i = 1, lim = p1.length; i < lim; i++) {
      var f = new THREE.Face4(2 * i, 2 * i + 1, 2 * i - 1, 2 * i - 2);
      f.color = new TCo(colors[Math.round((i - 1)/ div)]);
      geo.faces.push(f);
   }
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var material =  new THREE.MeshLambertMaterial();
   material.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, material);
   mesh.doubleSided = true;
   group.add(mesh);
};


GLmol.prototype.IcosahedronGeometry = function() {
   if (!this.icosahedron) this.icosahedron = new THREE.IcosahedronGeometry(1);
   return this.icosahedron;
};

GLmol.prototype.drawCylinder = function(group, from, to, radius, color, cap) {
   if (!from || !to) return;

   var midpoint = new TV3().add(from, to).multiplyScalar(0.5);
   var color = new TCo(color);

   if (!this.cylinderGeometry) {
      this.cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, this.cylinderQuality, 1, !cap);
      this.cylinderGeometry.faceUvs = [];
      this.faceVertexUvs = [];
   }
   var cylinderMaterial = new THREE.MeshLambertMaterial({color: color.getHex()});
   var cylinder = new THREE.Mesh(this.cylinderGeometry, cylinderMaterial);
   cylinder.position = midpoint;
   cylinder.lookAt(from);
   cylinder.updateMatrix();
   cylinder.matrixAutoUpdate = false;
   var m = new THREE.Matrix4().makeScale(radius, radius, from.distanceTo(to));
   m.rotateX(Math.PI / 2);
   cylinder.matrix.multiplySelf(m);
   group.add(cylinder);
};

// FIXME: transition!
GLmol.prototype.drawHelixAsCylinder = function(group, atomlist, radius) {
   var start = null;
   var currentChain, currentResi;

   var others = [], beta = [];

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;
      if ((atom.ss != 'h' && atom.ss != 's') || atom.ssend || atom.ssbegin) others.push(atom.serial);
      if (atom.ss == 's') beta.push(atom.serial);
      if (atom.atom != 'CA') continue;

      if (atom.ss == 'h' && atom.ssend) {
         if (start != null) this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(atom.x, atom.y, atom.z), radius, atom.color, true);
         start = null;
      }
      currentChain = atom.chain;
      currentResi = atom.resi;
      if (start == null && atom.ss == 'h' && atom.ssbegin) start = atom;
   }
   if (start != null) this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(atom.x, atom.y, atom.z), radius, atom.color);
   this.drawMainchainTube(group, others, "CA", 0.3);
   this.drawStrand(group, beta, undefined, undefined, true,  0, this.helixSheetWidth, false, this.thickness * 2);
};

GLmol.prototype.drawCartoon = function(group, atomlist, doNotSmoothen, thickness) {
   this.drawStrand(group, atomlist, 2, undefined, true, undefined, undefined, doNotSmoothen, thickness);
};

GLmol.prototype.drawStrand = function(group, atomlist, num, div, fill, coilWidth, helixSheetWidth, doNotSmoothen, thickness) {
   num = num || this.strandDIV;
   div = div || this.axisDIV;
   coilWidth = coilWidth || this.coilWidth;
   doNotSmoothen == (doNotSmoothen == undefined) ? false : doNotSmoothen;
   helixSheetWidth = helixSheetWidth || this.helixSheetWidth;
   var points = []; for (var k = 0; k < num; k++) points[k] = [];
   var colors = [];
   var currentChain, currentResi, currentCA;
   var prevCO = null, ss=null, ssborder = false;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == 'O' || atom.atom == 'CA') && !atom.hetflag) {
         if (atom.atom == 'CA') {
            if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
               for (var j = 0; !thickness && j < num; j++)
                  this.drawSmoothCurve(group, points[j], 1 ,colors, div);
               if (fill) this.drawStrip(group, points[0], points[num - 1], colors, div, thickness);
               var points = []; for (var k = 0; k < num; k++) points[k] = [];
               colors = [];
               prevCO = null; ss = null; ssborder = false;
            }
            currentCA = new TV3(atom.x, atom.y, atom.z);
            currentChain = atom.chain;
            currentResi = atom.resi;
            ss = atom.ss; ssborder = atom.ssstart || atom.ssend;
            colors.push(atom.color);
         } else { // O
            var O = new TV3(atom.x, atom.y, atom.z);
            O.subSelf(currentCA);
            O.normalize(); // can be omitted for performance
            O.multiplyScalar((ss == 'c') ? coilWidth : helixSheetWidth); 
            if (prevCO != undefined && O.dot(prevCO) < 0) O.negate();
            prevCO = O;
            for (var j = 0; j < num; j++) {
               var delta = -1 + 2 / (num - 1) * j;
               var v = new TV3(currentCA.x + prevCO.x * delta, 
                               currentCA.y + prevCO.y * delta, currentCA.z + prevCO.z * delta);
               if (!doNotSmoothen && ss == 's') v.smoothen = true;
               points[j].push(v);
            }                         
         }
      }
   }
   for (var j = 0; !thickness && j < num; j++)
      this.drawSmoothCurve(group, points[j], 1 ,colors, div);
   if (fill) this.drawStrip(group, points[0], points[num - 1], colors, div, thickness);
};

GLmol.prototype.drawNucleicAcidLadderSub = function(geo, lineGeo, atoms, color) {
//        color.r *= 0.9; color.g *= 0.9; color.b *= 0.9;
   if (atoms[0] != undefined && atoms[1] != undefined && atoms[2] != undefined &&
       atoms[3] != undefined && atoms[4] != undefined && atoms[5] != undefined) {
      var baseFaceId = geo.vertices.length;
      for (var i = 0; i <= 5; i++) geo.vertices.push(atoms[i]);
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 4, baseFaceId + 5));
          for (var j = geo.faces.length - 4, lim = geo.faces.length; j < lim; j++) geo.faces[j].color = color;
    }
    if (atoms[4] != undefined && atoms[3] != undefined && atoms[6] != undefined &&
       atoms[7] != undefined && atoms[8] != undefined) {
       var baseFaceId = geo.vertices.length;
       geo.vertices.push(atoms[4]);
       geo.vertices.push(atoms[3]);
       geo.vertices.push(atoms[6]);
       geo.vertices.push(atoms[7]);
       geo.vertices.push(atoms[8]);
       for (var i = 0; i <= 4; i++) geo.colors.push(color);
       geo.faces.push(new TF3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
       geo.faces.push(new TF3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
       geo.faces.push(new TF3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
       for (var j = geo.faces.length - 3, lim = geo.faces.length; j < lim; j++) geo.faces[j].color = color;
    }
};

GLmol.prototype.drawNucleicAcidLadder = function(group, atomlist) {
   var geo = new THREE.Geometry();
   var lineGeo = new THREE.Geometry();
   var baseAtoms = ["N1", "C2", "N3", "C4", "C5", "C6", "N9", "C8", "N7"];
   var currentChain, currentResi, currentComponent = new Array(baseAtoms.length);
   var color = new TCo(0xcc0000);
   
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);
         currentComponent = new Array(baseAtoms.length);
      }
      var pos = baseAtoms.indexOf(atom.atom);
      if (pos != -1) currentComponent[pos] = new TV3(atom.x, atom.y, atom.z);
      if (atom.atom == 'O3\'') color = new TCo(atom.color);
      currentResi = atom.resi; currentChain = atom.chain;
   }
   this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);
   geo.computeFaceNormals();
   var mat = new THREE.MeshLambertMaterial();
   mat.vertexColors = THREE.VertexColors;
   var mesh = new THREE.Mesh(geo, mat);
   mesh.doubleSided = true;
   group.add(mesh);
};

GLmol.prototype.drawNucleicAcidStick = function(group, atomlist) {
   var currentChain, currentResi, start = null, end = null;
   
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         if (start != null && end != null)
            this.drawCylinder(group, new TV3(start.x, start.y, start.z), 
                              new TV3(end.x, end.y, end.z), 0.3, start.color, true);
         start = null; end = null;
      }
      if (atom.atom == 'O3\'') start = atom;
      if (atom.resn == '  A' || atom.resn == '  G' || atom.resn == ' DA' || atom.resn == ' DG') {
         if (atom.atom == 'N1')  end = atom; //  N1(AG), N3(CTU)
      } else if (atom.atom == 'N3') {
         end = atom;
      }
      currentResi = atom.resi; currentChain = atom.chain;
   }
   if (start != null && end != null)
      this.drawCylinder(group, new TV3(start.x, start.y, start.z), 
                        new TV3(end.x, end.y, end.z), 0.3, start.color, true);
};

GLmol.prototype.drawNucleicAcidLine = function(group, atomlist) {
   var currentChain, currentResi, start = null, end = null;
   var geo = new THREE.Geometry();

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         if (start != null && end != null) {
            geo.vertices.push(new TV3(start.x, start.y, start.z));
            geo.colors.push(new TCo(start.color));
            geo.vertices.push(new TV3(end.x, end.y, end.z));
            geo.colors.push(new TCo(start.color));
         }
         start = null; end = null;
      }
      if (atom.atom == 'O3\'') start = atom;
      if (atom.resn == '  A' || atom.resn == '  G' || atom.resn == ' DA' || atom.resn == ' DG') {
         if (atom.atom == 'N1')  end = atom; //  N1(AG), N3(CTU)
      } else if (atom.atom == 'N3') {
         end = atom;
      }
      currentResi = atom.resi; currentChain = atom.chain;
   }
   if (start != null && end != null) {
      geo.vertices.push(new TV3(start.x, start.y, start.z));
      geo.colors.push(new TCo(start.color));
      geo.vertices.push(new TV3(end.x, end.y, end.z));
      geo.colors.push(new TCo(start.color));
    }
   var mat =  new THREE.LineBasicMaterial({linewidth: 1, linejoin: false});
   mat.linewidth = 1.5; mat.vertexColors = true;
   var line = new THREE.Line(geo, mat, THREE.LinePieces);
   group.add(line);
};

GLmol.prototype.drawCartoonNucleicAcid = function(group, atomlist, div, thickness) {
        this.drawStrandNucleicAcid(group, atomlist, 2, div, true, undefined, thickness);
};

GLmol.prototype.drawStrandNucleicAcid = function(group, atomlist, num, div, fill, nucleicAcidWidth, thickness) {
   nucleicAcidWidth = nucleicAcidWidth || this.nucleicAcidWidth;
   div = div || this.axisDIV;
   num = num || this.nucleicAcidStrandDIV;
   var points = []; for (var k = 0; k < num; k++) points[k] = [];
   var colors = [];
   var currentChain, currentResi, currentO3;
   var prevOO = null;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == 'O3\'' || atom.atom == 'OP2') && !atom.hetflag) {
         if (atom.atom == 'O3\'') { // to connect 3' end. FIXME: better way to do?
            if (currentChain != atom.chain || currentResi + 1 != atom.resi) {               
               if (currentO3) {
                  for (var j = 0; j < num; j++) {
                     var delta = -1 + 2 / (num - 1) * j;
                     points[j].push(new TV3(currentO3.x + prevOO.x * delta, 
                      currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
                  }
               }
               if (fill) this.drawStrip(group, points[0], points[1], colors, div, thickness);
               for (var j = 0; !thickness && j < num; j++)
                  this.drawSmoothCurve(group, points[j], 1 ,colors, div);
               var points = []; for (var k = 0; k < num; k++) points[k] = [];
               colors = [];
               prevOO = null;
            }
            currentO3 = new TV3(atom.x, atom.y, atom.z);
            currentChain = atom.chain;
            currentResi = atom.resi;
            colors.push(atom.color);
         } else { // OP2
            if (!currentO3) {prevOO = null; continue;} // for 5' phosphate (e.g. 3QX3)
            var O = new TV3(atom.x, atom.y, atom.z);
            O.subSelf(currentO3);
            O.normalize().multiplyScalar(nucleicAcidWidth);  // TODO: refactor
            if (prevOO != undefined && O.dot(prevOO) < 0) {
               O.negate();
            }
            prevOO = O;
            for (var j = 0; j < num; j++) {
               var delta = -1 + 2 / (num - 1) * j;
               points[j].push(new TV3(currentO3.x + prevOO.x * delta, 
                 currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
            }
            currentO3 = null;
         }
      }
   }
   if (currentO3) {
      for (var j = 0; j < num; j++) {
         var delta = -1 + 2 / (num - 1) * j;
         points[j].push(new TV3(currentO3.x + prevOO.x * delta, 
           currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
      }
   }
   if (fill) this.drawStrip(group, points[0], points[1], colors, div, thickness); 
   for (var j = 0; !thickness && j < num; j++)
      this.drawSmoothCurve(group, points[j], 1 ,colors, div);
};

GLmol.prototype.drawDottedLines = function(group, points, color) {
    var geo = new THREE.Geometry();
    var step = 0.3;

    for (var i = 0, lim = Math.floor(points.length / 2); i < lim; i++) {
        var p1 = points[2 * i], p2 = points[2 * i + 1];
        var delta = p2.clone().subSelf(p1);
        var dist = delta.length();
        delta.normalize().multiplyScalar(step);
        var jlim =  Math.floor(dist / step);
        for (var j = 0; j < jlim; j++) {
           var p = new TV3(p1.x + delta.x * j, p1.y + delta.y * j, p1.z + delta.z * j);
           geo.vertices.push(p);
        }
        if (jlim % 2 == 1) geo.vertices.push(p2);
    }

    var mat = new THREE.LineBasicMaterial({'color': color.getHex()});
    mat.linewidth = 2;
    var line = new THREE.Line(geo, mat, THREE.LinePieces);
    group.add(line);
};

GLmol.prototype.getAllAtoms = function() {
   var ret = [];
   for (var i in this.atoms) {
      ret.push(this.atoms[i].serial);
   }
   return ret;
};

// Probably I can refactor using higher-order functions.
GLmol.prototype.getHetatms = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.removeSolvents = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.resn != 'HOH') ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getProteins = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (!atom.hetflag) ret.push(atom.serial);
   }
   return ret;
};

// TODO: Test
GLmol.prototype.excludeAtoms = function(atomlist, deleteList) {
   var ret = [];
   var blackList = new Object();
   for (var _i in deleteList) blackList[deleteList[_i]] = true;

   for (var _i in atomlist) {
      var i = atomlist[_i];

      if (!blackList[i]) ret.push(i);
   }
   return ret;
};

GLmol.prototype.getSidechains = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (atom.atom == 'C' || atom.atom == 'O' || (atom.atom == 'N' && atom.resn != "PRO")) continue;
      ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getAtomsWithin = function(atomlist, extent) {
   var ret = [];

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.x < extent[0][0] || atom.x > extent[1][0]) continue;
      if (atom.y < extent[0][1] || atom.y > extent[1][1]) continue;
      if (atom.z < extent[0][2] || atom.z > extent[1][2]) continue;
      ret.push(atom.serial);      
   }
   return ret;
};

GLmol.prototype.getExtent = function(atomlist) {
   var xmin = ymin = zmin = 9999;
   var xmax = ymax = zmax = -9999;
   var xsum = ysum = zsum = cnt = 0;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;
      cnt++;
      xsum += atom.x; ysum += atom.y; zsum += atom.z;

      xmin = (xmin < atom.x) ? xmin : atom.x;
      ymin = (ymin < atom.y) ? ymin : atom.y;
      zmin = (zmin < atom.z) ? zmin : atom.z;
      xmax = (xmax > atom.x) ? xmax : atom.x;
      ymax = (ymax > atom.y) ? ymax : atom.y;
      zmax = (zmax > atom.z) ? zmax : atom.z;
   }
   return [[xmin, ymin, zmin], [xmax, ymax, zmax], [xsum / cnt, ysum / cnt, zsum / cnt]];
};

GLmol.prototype.getResiduesById = function(atomlist, resi) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (resi.indexOf(atom.resi) != -1) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getResidueBySS = function(atomlist, ss) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (ss.indexOf(atom.ss) != -1) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getChain = function(atomlist, chain) {
   var ret = [], chains = {};
   chain = chain.toString(); // concat if Array
   for (var i = 0, lim = chain.length; i < lim; i++) chains[chain.substr(i, 1)] = true;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (chains[atom.chain]) ret.push(atom.serial);
   }
   return ret;
};

// for HETATM only
GLmol.prototype.getNonbonded = function(atomlist, chain) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag && atom.bonds.length == 0) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.colorByAtom = function(atomlist, colors) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      var c = colors[atom.elem];
      if (c == undefined) c = this.ElementColors[atom.elem];
      if (c == undefined) c = this.defaultColor;
      atom.color = c;
   }
};


// MEMO: Color only CA. maybe I should add atom.cartoonColor.
GLmol.prototype.colorByStructure = function(atomlist, helixColor, sheetColor, colorSidechains) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (!colorSidechains && (atom.atom != 'CA' || atom.hetflag)) continue;
      if (atom.ss[0] == 's') atom.color = sheetColor;
      else if (atom.ss[0] == 'h') atom.color = helixColor;
   }
};

GLmol.prototype.colorByBFactor = function(atomlist, colorSidechains) {
   var minB = 1000, maxB = -1000;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         if (minB > atom.b) minB = atom.b;
         if (maxB < atom.b) maxB = atom.b;
      }
   }

   var mid = (maxB + minB) / 2;

   var range = (maxB - minB) / 2;
   if (range < 0.01 && range > -0.01) return;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         var color = new TCo(0);
         if (atom.b < mid)
            color.setHSV(0.667, (mid - atom.b) / range, 1);
         else
            color.setHSV(0, (atom.b - mid) / range, 1);
         atom.color = color.getHex();
      }
   }
};

GLmol.prototype.colorByChain = function(atomlist, colorSidechains) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         var color = new TCo(0);
         color.setHSV((atom.chain.charCodeAt(0) * 5) % 17 / 17.0, 1, 0.9);
         atom.color = color.getHex();
      }
   }
};

GLmol.prototype.colorByResidue = function(atomlist, residueColors) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      c = residueColors[atom.resn]
      if (c != undefined) atom.color = c;
   }
};

GLmol.prototype.colorAtoms = function(atomlist, c) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      atom.color = c;
   }
};

GLmol.prototype.colorByPolarity = function(atomlist, polar, nonpolar) {
   var polarResidues = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS'];
   var nonPolarResidues = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP'];
   var colorMap = {};
   for (var i in polarResidues) colorMap[polarResidues[i]] = polar;
   for (i in nonPolarResidues) colorMap[nonPolarResidues[i]] = nonpolar;
   this.colorByResidue(atomlist, colorMap);   
};

// TODO: Add near(atomlist, neighbor, distanceCutoff)
// TODO: Add expandToResidue(atomlist)

GLmol.prototype.colorChainbow = function(atomlist, colorSidechains) {
   var cnt = 0;
   var atom, i;
   for (i in atomlist) {
      atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if ((colorSidechains || atom.atom != 'CA' || atom.atom != 'O3\'') && !atom.hetflag)
         cnt++;
   }

   var total = cnt;
   cnt = 0;
   for (i in atomlist) {
      atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if ((colorSidechains || atom.atom != 'CA' || atom.atom != 'O3\'') && !atom.hetflag) {
         var color = new TCo(0);
         color.setHSV(240.0 / 360 * (1 - cnt / total), 1, 0.9);
         atom.color = color.getHex();
         cnt++;
      }
   }
};

GLmol.prototype.drawSymmetryMates2 = function(group, asu, matrices) {
   if (matrices == undefined) return;
   asu.matrixAutoUpdate = false;

   var cnt = 1;
   this.protein.appliedMatrix = new THREE.Matrix4();
   for (var i = 0; i < matrices.length; i++) {
      var mat = matrices[i];
      if (mat == undefined || mat.isIdentity()) continue;
      console.log(mat);
      var symmetryMate = THREE.SceneUtils.cloneObject(asu);
      symmetryMate.matrix = mat;
      group.add(symmetryMate);
      for (var j = 0; j < 16; j++) this.protein.appliedMatrix.elements[j] += mat.elements[j];
      cnt++;
   }
   this.protein.appliedMatrix.multiplyScalar(cnt);
};


GLmol.prototype.drawSymmetryMatesWithTranslation2 = function(group, asu, matrices) {
   if (matrices == undefined) return;
   var p = this.protein;
   asu.matrixAutoUpdate = false;

   for (var i = 0; i < matrices.length; i++) {
      var mat = matrices[i];
      if (mat == undefined) continue;

      for (var a = -1; a <=0; a++) {
         for (var b = -1; b <= 0; b++) {
             for (var c = -1; c <= 0; c++) {
                var translationMat = new THREE.Matrix4().makeTranslation(
                   p.ax * a + p.bx * b + p.cx * c,
                   p.ay * a + p.by * b + p.cy * c,
                   p.az * a + p.bz * b + p.cz * c);
                var symop = mat.clone().multiplySelf(translationMat);
                if (symop.isIdentity()) continue;
                var symmetryMate = THREE.SceneUtils.cloneObject(asu);
                symmetryMate.matrix = symop;
                group.add(symmetryMate);
             }
         }
      }
   }
};

GLmol.prototype.defineRepresentation = function() {
   var all = this.getAllAtoms();
   var hetatm = this.removeSolvents(this.getHetatms(all));
   this.colorByAtom(all, {});
   this.colorByChain(all);

   this.drawAtomsAsSphere(this.modelGroup, hetatm, this.sphereRadius); 
   this.drawMainchainCurve(this.modelGroup, all, this.curveWidth, 'P');
   this.drawCartoon(this.modelGroup, all, this.curveWidth);
};

GLmol.prototype.getView = function() {
   if (!this.modelGroup) return [0, 0, 0, 0, 0, 0, 0, 1];
   var pos = this.modelGroup.position;
   var q = this.rotationGroup.quaternion;
   return [pos.x, pos.y, pos.z, this.rotationGroup.position.z, q.x, q.y, q.z, q.w];
};

GLmol.prototype.setView = function(arg) {
   if (!this.modelGroup || !this.rotationGroup) return;
   this.modelGroup.position.x = arg[0];
   this.modelGroup.position.y = arg[1];
   this.modelGroup.position.z = arg[2];
   this.rotationGroup.position.z = arg[3];
   this.rotationGroup.quaternion.x = arg[4];
   this.rotationGroup.quaternion.y = arg[5];
   this.rotationGroup.quaternion.z = arg[6];
   this.rotationGroup.quaternion.w = arg[7];
   this.show();
};

GLmol.prototype.setBackground = function(hex, a) {
   a = a | 1.0;
   this.bgColor = hex;
   this.renderer.setClearColorHex(hex, a);
   this.scene.fog.color = new TCo(hex);
};

GLmol.prototype.initializeScene = function() {
   // CHECK: Should I explicitly call scene.deallocateObject?
   this.scene = new THREE.Scene();
   this.scene.fog = new THREE.Fog(this.bgColor, 100, 200);

   this.modelGroup = new THREE.Object3D();
   this.rotationGroup = new THREE.Object3D();
   this.rotationGroup.useQuaternion = true;
   this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
   this.rotationGroup.add(this.modelGroup);

   this.scene.add(this.rotationGroup);
   this.setupLights(this.scene);
};

GLmol.prototype.zoomInto = function(atomlist, keepSlab) {
   var tmp = this.getExtent(atomlist);
   var center = new TV3(tmp[2][0], tmp[2][1], tmp[2][2]);//(tmp[0][0] + tmp[1][0]) / 2, (tmp[0][1] + tmp[1][1]) / 2, (tmp[0][2] + tmp[1][2]) / 2);
   if (this.protein.appliedMatrix) {center = this.protein.appliedMatrix.multiplyVector3(center);}
   this.modelGroup.position = center.multiplyScalar(-1);
   var x = tmp[1][0] - tmp[0][0], y = tmp[1][1] - tmp[0][1], z = tmp[1][2] - tmp[0][2];

   var maxD = Math.sqrt(x * x + y * y + z * z);
   if (maxD < 25) maxD = 25;

   if (!keepSlab) {
      this.slabNear = -maxD / 1.9;
      this.slabFar = maxD / 3;
   }

   this.rotationGroup.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * this.camera.fov / 2) - 150;
   this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
};

GLmol.prototype.rebuildScene = function() {
   time = new Date();

   var view = this.getView();
   this.initializeScene();
   this.defineRepresentation();
   this.setView(view);

   console.log("builded scene in " + (+new Date() - time) + "ms");
};

GLmol.prototype.loadMolecule = function(repressZoom) {
   this.loadMoleculeStr(repressZoom, $('#' + this.id + '_src').val());
};

GLmol.prototype.loadMoleculeStr = function(repressZoom, source) {
   var time = new Date();

   this.protein = {sheet: [], helix: [], biomtChains: '', biomtMatrices: [], symMat: [], pdbID: '', title: ''};
   this.atoms = [];

   this.parsePDB2(source);
   if (!this.parseSDF(source)) this.parseXYZ(source);
   console.log("parsed in " + (+new Date() - time) + "ms");
   
   var title = $('#' + this.id + '_pdbTitle');
   var titleStr = '';
   if (this.protein.pdbID != '') titleStr += '<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' + this.protein.pdbID + '">' + this.protein.pdbID + '</a>';
   if (this.protein.title != '') titleStr += '<br>' + this.protein.title;
   title.html(titleStr);

   this.rebuildScene(true);
   if (repressZoom == undefined || !repressZoom) this.zoomInto(this.getAllAtoms());

   this.show();
 };

GLmol.prototype.setSlabAndFog = function() {
   var center = this.rotationGroup.position.z - this.camera.position.z;
   if (center < 1) center = 1;
   this.camera.near = center + this.slabNear;
   if (this.camera.near < 1) this.camera.near = 1;
   this.camera.far = center + this.slabFar;
   if (this.camera.near + 1 > this.camera.far) this.camera.far = this.camera.near + 1;
   if (this.camera instanceof THREE.PerspectiveCamera) {
      this.camera.fov = this.fov;
   } else {
      this.camera.right = center * Math.tan(Math.PI / 180 * this.fov);
      this.camera.left = - this.camera.right;
      this.camera.top = this.camera.right / this.ASPECT;
      this.camera.bottom = - this.camera.top;
   }
   this.camera.updateProjectionMatrix();
   this.scene.fog.near = this.camera.near + this.fogStart * (this.camera.far - this.camera.near);
//   if (this.scene.fog.near > center) this.scene.fog.near = center;
   this.scene.fog.far = this.camera.far;
};

GLmol.prototype.enableMouse = function() {
   var me = this, glDOM = $(this.renderer.domElement); 

   // TODO: Better touch panel support. 
   // Contribution is needed as I don't own any iOS or Android device with WebGL support.
   glDOM.bind('mousedown touchstart', function(ev) {
      ev.preventDefault();
      if (!me.scene) return;
      var x = ev.pageX, y = ev.pageY;
      if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
         x = ev.originalEvent.targetTouches[0].pageX;
         y = ev.originalEvent.targetTouches[0].pageY;
      }
      if (x == undefined) return;
      me.isDragging = true;
      me.mouseButton = ev.which;
      me.mouseStartX = x;
      me.mouseStartY = y;
      me.cq = me.rotationGroup.quaternion;
      me.cz = me.rotationGroup.position.z;
      me.currentModelPos = me.modelGroup.position.clone();
      me.cslabNear = me.slabNear;
      me.cslabFar = me.slabFar;
    });

   glDOM.bind('DOMMouseScroll mousewheel', function(ev) { // Zoom
      ev.preventDefault();
      if (!me.scene) return;
      var scaleFactor = (me.rotationGroup.position.z - me.CAMERA_Z) * 0.85;
      if (ev.originalEvent.detail) { // Webkit
         me.rotationGroup.position.z += scaleFactor * ev.originalEvent.detail / 10;
      } else if (ev.originalEvent.wheelDelta) { // Firefox
         me.rotationGroup.position.z -= scaleFactor * ev.originalEvent.wheelDelta / 400;
      }
      console.log(ev.originalEvent.wheelDelta, ev.originalEvent.detail, me.rotationGroup.position.z);
      me.show();
   });
   glDOM.bind("contextmenu", function(ev) {ev.preventDefault();});
   $('body').bind('mouseup touchend', function(ev) {
      me.isDragging = false;
   });

   glDOM.bind('mousemove touchmove', function(ev) { // touchmove
      ev.preventDefault();
      if (!me.scene) return;
      if (!me.isDragging) return;
      var mode = 0;
      var modeRadio = $('input[name=' + me.id + '_mouseMode]:checked');
      if (modeRadio.length > 0) mode = parseInt(modeRadio.val());

      var x = ev.pageX, y = ev.pageY;
      if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
         x = ev.originalEvent.targetTouches[0].pageX;
         y = ev.originalEvent.targetTouches[0].pageY;
      }
      if (x == undefined) return;
      var dx = (x - me.mouseStartX) / me.WIDTH;
      var dy = (y - me.mouseStartY) / me.HEIGHT;
      var r = Math.sqrt(dx * dx + dy * dy);
      if (mode == 3 || (me.mouseButton == 3 && ev.ctrlKey)) { // Slab
          me.slabNear = me.cslabNear + dx * 100;
          me.slabFar = me.cslabFar + dy * 100;
      } else if (mode == 2 || me.mouseButton == 3 || ev.shiftKey) { // Zoom
         var scaleFactor = (me.rotationGroup.position.z - me.CAMERA_Z) * 0.85; 
         if (scaleFactor < 80) scaleFactor = 80;
         me.rotationGroup.position.z = me.cz - dy * scaleFactor;
      } else if (mode == 1 || me.mouseButton == 2 || ev.ctrlKey) { // Translate
         var scaleFactor = (me.rotationGroup.position.z - me.CAMERA_Z) * 0.85;
         if (scaleFactor < 20) scaleFactor = 20;
         var translationByScreen = new TV3(- dx * scaleFactor, - dy * scaleFactor, 0);
         var q = me.rotationGroup.quaternion;
         var qinv = new THREE.Quaternion(q.x, q.y, q.z, q.w).inverse().normalize(); 
         var translation = qinv.multiplyVector3(translationByScreen);
         me.modelGroup.position.x = me.currentModelPos.x + translation.x;
         me.modelGroup.position.y = me.currentModelPos.y + translation.y;
         me.modelGroup.position.z = me.currentModelPos.z + translation.z;
      } else if ((mode == 0 || me.mouseButton == 1) && r != 0) { // Rotate
         var rs = Math.sin(r * Math.PI) / r;
         me.dq.x = Math.cos(r * Math.PI); 
         me.dq.y = 0;
         me.dq.z =  rs * dx; 
         me.dq.w =  rs * dy;
         me.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0); 
         me.rotationGroup.quaternion.multiplySelf(me.dq);
         me.rotationGroup.quaternion.multiplySelf(me.cq);
      }
      me.show();
   });
};


GLmol.prototype.show = function() {
   if (!this.scene) return;

   var time = new Date();
   this.setSlabAndFog();
   this.renderer.render(this.scene, this.camera);
   console.log("rendered in " + (+new Date() - time) + "ms");
};

// For scripting
GLmol.prototype.doFunc = function(func) {
    func(this);
};

return GLmol;
}());
