(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}
	
	var functionBoundary = "[\\s\\(\\)]";
	
	function createFunctionRule(words, tokenName) {
		var map = sunlight.util.createHashMap(words, functionBoundary, false);
		return function(context) {
			//must be right after a "("
			var token = sunlight.util.getPreviousNonWsToken(context.getAllTokens(), context.count());
			if (!token) {
				return null;
			}
			
			if (token.name !== "punctuation" || token.value !== "(" && (token.name !== "operator" || token.value !== "#'")) {
				return null;
			}
			
			return sunlight.util.matchWord(context, map, tokenName);
		};
	}
	
	sunlight.registerLanguage("lisp", {
		scopes: {
			string: [ ["\"", "\""] ],
			comment: [ [";", "\n", null, true] ],
			keywordArgument: [ [":", { regex: /[\s\)\(]/, length: 1 }, null, true] ]
		},
		
		customTokens: {
			//lisp allows "-" as part of the ident, which is matched by \b so we have to set the boundary manually
			keyword: {
				values: [
					"always","appending","append","as","collecting","collect","counting","count","doing","do","finally","for","if",
					"initially","maximize","maximizing","minimize","minimizing","named","nconcing","nconc","never","repeat",
					"return","summing","sum","thereis","unless","until","when","while","with",
			
					//other loop keywords
					"into", "in"
				],
				boundary: "[^\\w-]"
			},
			
			globalVariable: {
				values: [
					"*applyhook*","*break-on-signals*","*break-on-warnings*","*compile-file-pathname*","*compile-file-truename*","*compile-print*",
					"*compile-verbose*","*debug-io*","*debugger-hook*","*gensym-counter*","*terminal-io*","*default-pathname-defaults*","*error-output*",
					"*evalhook*","*features*","*load-pathname*","*load-print*","*load-truename*","*load-verbose*","*macroexpand-hook*","*modules*",
					"*package*","*print-array*","*print-base*","*print-case*","*print-circle*","*print-escape*","*print-gensym*","*print-length*",
					"*print-level*","*print-lines*","*print-miser-width*","*print-pprint-dispatch*","*print-pretty*","*print-radix*","*print-readably*",
					"*print-right-margin*","*query-io*","*random-state*","*read-base*","*read-default-float-format*","*read-eval*","*read-suppress*",
					"*readtable*","*standard-input*","*standard-output*","*suppress-series-warnings*","*trace-output*",
					
					"***", "**", "*",
					"+++", "++", "+",
					"-",
					"///", "//", "/"
				],
				boundary: functionBoundary
			},
			
			constant: {
				values: [
					"array-dimension-limit","array-rank-limit","array-total-size-limit","call-arguments-limit","char-bits-limit","char-code-limit",
					"char-control-bit","char-font-limit","char-hyper-bit","char-meta-bit","char-super-bit","double-float-epsilon",
					"double-float-negative-epsilon","internal-time-units-per-second","lambda-list-keywords","lambda-parameters-limit",
					"least-negative-double-float","least-negative-long-float","least-negative-normalized-double-float",
					"least-negative-normalized-long-float","least-negative-normalized-short-float","least-negative-normalized-single-float",
					"least-negative-short-float","least-negative-single-float","least-positive-double-float","least-positive-long-float",
					"least-positive-normalized-double-float","least-positive-normalized-long-float","least-positive-normalized-short-float",
					"least-positive-normalized-single-float","least-positive-short-float","least-positive-single-float","long-float-epsilon",
					"long-float-negative-epsilon","most-negative-double-float","most-negative-fixnum","most-negative-long-float",
					"most-negative-short-float","most-negative-single-float","most-positive-double-float","most-positive-fixnum",
					"most-positive-long-float","most-positive-short-float","most-positive-single-float","multiple-values-limit","nil","pi",
					"short-float-epsilon","short-float-negative-epsilon","single-float-epsilon","single-float-negative-epsilon","t"
				],
				boundary: functionBoundary
			},
			
			declarationSpecifier: {
				values: ["off-line-port","optimizable-series-function","propagate-alterability"],
				boundary: functionBoundary
			}
		},
		
		customParseRules: [
			//#<n>r operator where n is an integer
			function(context) {
				var peek, 
					count = 0,
					token;
				
				if (context.reader.current() !== "#") {
					return null;
				}
				
				while ((peek = context.reader.peek(++count)) !== context.reader.EOF && peek.length === count) {
					if (!/\d$/.test(peek)) {
						if (/[AR]$/i.test(peek)) {
							break;
						}
						
						return null;
					}
				}
				
				if (peek.length !== count) {
					return null;
				}
				
				token = context.createToken("operator", "#" + peek, context.reader.getLine(), context.reader.getColumn());
				context.reader.read(peek.length);
				return token;
			},
			
			//characters prepended by the #\ operator are read as idents
			function(context) {
				var prevToken,
					peek,
					value,
					line = context.reader.getLine(), 
					column = context.reader.getColumn();
					
				if (context.defaultData.text !== "" || /\s/.test(context.reader.current())) {
					//whitespace is not allowed
					return null;
				}
				
				prevToken = context.getAllTokens()[context.count() - 1];
				if (!prevToken || prevToken.name !== "operator" || prevToken.value !== "#\\") {
					return null;
				}
				
				//the next characters up until whitespace or ( or ) are part of the ident
				value = context.reader.current();
				while ((peek = context.reader.peek()) !== context.reader.EOF) {
					if (/[\s\(\)]/.test(peek)) {
						break;
					}
					
					value += context.reader.read();
				}
				
				return context.createToken("ident", value, line, column);
			},
			
			//variables
			function(context) {
				var token,
					value = "*",
					peek, 
					line = context.reader.getLine(), 
					column = context.reader.getColumn();
				
				if (context.reader.current() !== "*") {
					return null;
				}
				
				token = sunlight.util.getPreviousNonWsToken(context.getAllTokens(), context.count());
				if (token && token.name === "punctuation" && token.value === "(") {
					//function that starts with "*"
					return null;
				}
				
				if (/[\s\*\)\(]/.test(context.reader.peek())) {
					return null;
				}
				
				//read until *
				while ((peek = context.reader.peek()) !== context.reader.EOF) {
					value += context.reader.read();
					
					if (peek === "*") {
						break;
					}
				}
				
				return context.createToken("variable", value, line, column);
			},
			
			//function after #' operator
			function() {
				var boundary = new RegExp(functionBoundary);
			
				return function(context) {
					var token,
						peek,
						value,
						line = context.reader.getLine(), 
						column = context.reader.getColumn();
					
					if (context.defaultData.text !== "" || boundary.test(context.reader.current())) {
						//whitespace is not allowed or we're already at the boundary
						return null;
					}
					
					token = context.getAllTokens()[context.count() - 1];
					if (!token || token.name !== "operator" || token.value !== "#'") {
						return null;
					}
					
					//read until function boundary
					value = context.reader.current();
					while ((peek = context.reader.peek()) !== context.reader.EOF) {
						if (boundary.test(peek)) {
							break;
						}
						
						value += context.reader.read();
					}
					
					return context.createToken("function", value, line, column);
				};
			}(),
			
			//types/specifiers
			createFunctionRule([
				"arithmetic-error","cell-error","condition","control-error","division-by-zero","end-of-file","error","file-error",
				"floating-point-overflow","floating-point-underflow","package-error","program-error","restart","series","series-element-type",
				"serious-condition","simple-condition","simple-error","simple-type-error","simple-warning","storage-condition","stream-error",
				"type-error","unbound-variable","undefined-function","warning"
			], "type"),
			
			//special forms
			createFunctionRule([
				"block","catch","compiler-let","declare","eval-when","flet","function","generic-flet","generic-labels","go","if","let*",
				"let","load-time-value","locally","multiple-value-call","multiple-value-prog1","progn","progv","quote","return-from","setq",
				"symbol-macrolet","tagbody","the","throw","unwind-protect","with-added-methods"
			], "specialForm"),
			
			//macros
			createFunctionRule([
				"and","assert","call-method","case","ccase","check-type","compiler-let","cond","ctypecase","decf","declaim","defclass","defgeneric",
				"define-compiler-macro","define-condition","define-declaration","define-method-combination","define-modify-macro","define-setf-method",
				"defmacro","defmethod","defpackage","defstruct","deftype","defun","defvar","destructuring-bind","do*","do-all-symbols","do-external-symbols",
				"do-symbols","dolist","dotimes","do","ecase","encapsulated","etypecase","formatter","gathering","generic-function","handler-bind",
				"handler-case","ignore-errors","in-package","incf","iterate","locally","loop-finish","loop","mapping","multiple-value-bind","multiple-value-list",
				"multiple-value-setq","next-in","nth-value","or","pop","pprint-exit-if-list-exhausted","pprint-logical-block","pprint-pop",
				"print-unreadable-object","producing","prog*","prog1","prog2","prog","psetf","psetq","pushnew","push","remf","restart-bind","restart-case",
				"return","rotatef","setf","shiftf","step","terminate-producing","time","trace","typecase","unless","untrace","when","with-accessors",
				"with-compilation-unit","with-condition-restarts","with-hash-table-iterator","with-input-from-string","with-open-file","with-open-stream",
				"with-output-to-string","with-package-iterator","with-simple-restart","with-slots","with-standard-io-syntax",
				
				//others
				"defparameter"
			], "macro"),
			
			//functions
			createFunctionRule([
				"*", "+", "-", "/", "1+", "1-", "<=", "<",  ">=", ">", "=",
				
				//i guess this is a function
				"lambda",
				
				"abort","abs","acons","acosh","acos","add-method","adjoin","adjust-array","adjustable-array-p","alpha-char-p","alphanumericp",
				"alter","append","applyhook","apply","apropos-list","apropos","aref","arithmetic-error-operands","arithmetic-error-operation",
				"array-dimensions","array-dimension","array-element-type","array-has-fill-pointer-p","array-in-bounds-p","array-rank",
				"array-row-major-index","array-total-size","arrayp","ash","asinh","asin","assoc-if","assoc-if-not","assoc","atanh","atan","atom",
				"augment-environment","bit-andc1","bit-andc2","bit-and","bit-eqv","bit-ior","bit-nand","bit-nor","bit-not","bit-orc1","bit-orc2",
				"bit-vector-p","bit-xor","bit","boole","both-case-p","boundp","break","broadcast-stream-streams","butlast","byte-position",
				"byte-size","byte","caaaar","caaadr","caaar","caadar","caaddr","caadr","caar","cadaar","cadadr","cadar","caddar","cadddr","caddr",
				"cadr","call-next-method","car","catenate","cdaaar","cdaadr","cdaar","cdadar","cdaddr","cdadr","cdar","cddaar","cddadr","cddar",
				"cdddar","cddddr","cdddr","cddr","cdr","ceiling","cell-error-name","cerror","change-class","char-bits","char-bit","char-code",
				"char-downcase","char-equal","char-font","char-greaterp","char-int","char-lessp","char-name","char-not-equal","char-not-greaterp",
				"char-not-lessp","char-upcase","char/=","char<=","char<","char=","char>=","char>","char","characterp","character","choose-if",
				"choose","chunk","cis","class-name","class-of","clear-input","close","clrhash","code-char","coerce","collect-alist","collect-and",
				"collect-append","collect-file","collect-first","collect-fn","collect-hash","collect-last","collect-length","collect-max","collect-min",
				"collect-nconc","collect-nth","collect-or","collect-plist","collect-sum","collecting-fn","collect","commonp","compile-file",
				"compile-file-pathname","compiled-function-p","compiler-macro-function","compiler-macroexpand","compiler-macroexpand-1","compile",
				"complement","complexp","complex","compute-applicable-methods","compute-restarts","concatenated-stream-streams","concatenate",
				"conjugate","consp","constantp","cons","continue","copy-alist","copy-list","copy-pprint-dispatch","copy-readtable","copy-seq",
				"copy-symbol","copy-tree","cosh","cos","cotruncate","count-if","count-if-not","count","declaration-information","decode-float",
				"decode-universal-time","delete-duplicates","delete-file","delete-if","delete-if-not","delete-package","delete","denominator",
				"deposit-field","describe-object","describe","digit-char-p","digit-char","directory-namestring","directory","disassemble","documentation",
				"dpb","dribble","echo-stream-input-stream","echo-stream-output-stream","ed","eighth","elt","enclose","encode-universal-time","endp",
				"enough-namestring","ensure-generic-function","eql","eq","equalp","equal","error","evalhook","eval","evenp","every","expand","export",
				"expt","exp","fboundp","fdefinition","ffloor","fifth","file-author","file-error-pathname","file-length","file-namestring","file-position",
				"file-string-length","file-write-date","fill-pointer","fill","find-all-symbols","find-class","find-if-not","find-if","find-method",
				"find-package","find-restart","find-symbol","find","finish-output","first","float-digits","float-precision","float-radix","float-sign",
				"floatp","float","floor","format","fourth","funcall","function-information","function-keywords","function-lambda-expression","functionp",
				"f","gatherer","gcd","generator","gensym","gentemp","get-decoded-time","get-internal-real-time","get-internal-run-time",
				"get-output-stream-string","get-properties","get-setf-method-multiple-value","get-setf-method","get-universal-time","getf","gethash",
				"get","graphic-char-p","hash-table-count","hash-table-p","hash-table-rehash-size","hash-table-rehash-threshold","hash-table-size",
				"hash-table-test","host-namestring","identity","imagpart","import","in-package","initialize-instance","input-stream-p","inspect",
				"int-char","integer-decode-float","integer-length","integerp","interactive-stream-p","intern","intersection","invalid-method-error",
				"invoke-debugger","invoke-restart","isqrt","keywordp","last","latch","lcm","ldb-test","ldb","ldiff","length","lisp-implementation-type",
				"lisp-implementation-version","list*","list-all-packages","list-length","listen","listp","list","load-logical-pathname-translations",
				"load","logandc1","logandc2","logand","logbitp","logcount","logeqv","logical-pathname-translations","logical-pathname","logior",
				"lognand","lognor","lognot","logorc1","logorc2","logtest","logxor","log","long-site-name","lower-case-p","machine-instance","machine-type",
				"machine-version","macro-function","macroexpand-1","macroexpand","make-array","make-broadcast-stream","make-char","make-concatenated-stream",
				"make-condition","make-dispatch-macro-character","make-echo-stream","make-hash-table","make-instances-obsolete","make-instance","make-list",
				"make-load-form-saving-slots","make-load-form","make-package","make-pathname","make-random-state","make-sequence","make-string-input-stream",
				"make-string-output-stream","make-string","make-symbol","make-synonym-stream","make-two-way-stream","makunbound","map-fn","map-into","mapcan",
				"mapcar","mapcon","mapc","maphash","maplist","mapl","map","mask-field","mask","max","member-if","member-if-not","member","merge-pathnames",
				"merge","method-combination-error","method-qualifiers","mingle","minusp","min","mismatch","mod","muffle-warning","name-char","namestring",
				"nbutlast","nconc","next-method-p","next-out","nintersection","ninth","no-applicable-method","no-next-method","notany","notevery","not",
				"nreconc","nreverse","nset-difference","nset-exclusive-or","nstring-capitalize","nstring-downcase","nstring-upcase","nsublis","nsubst",
				"nsubst-if-not","nsubst-if","nsubstitute-if-not","nsubstitute-if","nsubstitute","nthcdr","nth","null","numberp","numerator","nunion","oddp",
				"open-stream-p","open","output-stream-p","package-error-package","package-name","package-nicknames","package-shadowing-symbols",
				"package-use-list","package-used-by-list","packagep","pairlis","parse-integer","parse-macro","parse-namestring","pathname-device",
				"pathname-directory","pathname-host","pathname-match-p","pathname-name","pathname-type","pathname-version","pathnamep","pathname","peek-char",
				"phase","plusp","position-if-not","position-if","positions","position","pprint-dispatch","pprint-fill","pprint-indent","pprint-linear",
				"pprint-newline","pprint-tabular","pprint-tab","previous","prin1","print-object","probe-file","proclaim","provide","random-state-p",
				"random","rassoc-if-not","rassoc-if","rassoc","rationalize","rationalp","rational","read-byte","read-char-no-hang","read-char",
				"read-delimited-list","read-from-string","read-line","read","read-preserving-whitespace","readtable-case","readtablep","realpart","realp",
				"reduce","reinitialize-instance","remhash","remove-duplicates","remove-method","remove","remprop","rem","rename-file","rename-package",
				"replace","require","restart-name","rest","result-of","revappend","reverse","room","round","row-major-aref","rplaca","rplacd","sbit",
				"scale-float","scan-alist","scan-file","scan-fn-inclusive","scan-fn","scan-hash","scan-lists-of-lists-fringe","scan-lists-of-lists",
				"scan-multiple","scan-plist","scan-range","scan-sublists","scan-symbols","scan","schar","search","second","series","set-char-bit",
				"set-difference","set-dispatch-macro-character","set-exclusive-or","set-macro-character","set-pprint-dispatch","set-syntax-from-char",
				"set","seventh","shadow","shadowing-import","shared-initialize","short-site-name","signal","signum","simple-bit-vector-p",
				"simple-condition-format-arguments","simple-condition-format-string","simple-string-p","simple-vector-p","sinh","sin","sixth","sleep",
				"slot-boundp","slot-exists-p","slot-makunbound","slot-missing","slot-unbound","slot-value","software-type","software-version","some",
				"sort","special-form-p","split-if","split","sqrt","stable-sort","standard-char-p","store-value","stream-element-type","stream-error-stream",
				"stream-external-format","streamp","string-capitalize","string-char-p","string-downcase","string-equal","string-greaterp","string-left-trim",
				"string-lessp","string-not-equal","string-not-greaterp","string-not-lessp","string-right-trim","string-trim","string-upcase","string/=",
				"string<","string<=","string=","string>","string>=","stringp","string","sublis","subseq","subseries","subsetp","subst-if-not","subst-if",
				"substitute-if-not","substitute-if","substitute","subst","subtypep","svref","sxhash","symbol-function","symbol-name","symbol-package",
				"symbol-plist","symbol-value","symbolp","synonym-stream-symbol","tailp","tanh","tan","tenth","terpri","third","to-alter","translate-logical-pathname",
				"translate-pathname","tree-equal","truename","truncate","two-way-stream-input-stream","two-way-stream-output-stream","type-error-datum",
				"type-error-expected-type","type-of","typep","unexport","unintern","union","unread-char","until-if","until","unuse-package",
				"update-instance-for-different-class","update-instance-for-redefined-class","upgraded-array-element-type","upgraded-complex-part-type","upper-case-p",
				"use-package","use-value","user-homedir-pathname","values-list","values","variable-information","vector-pop","vector-push-extend","vector-push",
				"vectorp","vector","warn","wild-pathname-p","write-byte","write-char","write-string","write-to-string","write","y-or-n-p","yes-or-no-p","zerop"
			], "function")
		],
		
		identFirstLetter: /[A-Za-z]/,
		identAfterFirstLetter: /[\w-]/,
		
		namedIdentRules: {
			custom: [
				function() {
					var defMacros = [
						"defmacro","defmethod","defpackage","defstruct","deftype","defun","defvar","define-compiler-macro",
						"define-condition","define-declaration","define-method-combination","define-modify-macro",
						"define-setf-method"
					];
					
					return function(context) {
						//if the previous token is "defun" function, then this is a user-defined function
						var prevToken = sunlight.util.getPreviousNonWsToken(context.tokens, context.index),
							identValue = context.tokens[context.index].value;
						
						if (prevToken && prevToken.name === "macro" && sunlight.util.contains(defMacros, prevToken.value)) {
							context.items.userDefinedFunctions.push(identValue);
						}
						
						return sunlight.util.contains(context.items.userDefinedFunctions, identValue);
					};
				}()
			]
		},
		
		contextItems: {
			userDefinedFunctions: []
		},

		operators: [
			"=>",
			"#B", "#b", "#O", "#o", "#X", "#x",
			"#C", "#c", "#S", "#s", "#P", "#p",
			
			"#.", "#:", "#\\", "#'", "#",
			"'",
			"...", "..", "."
		]
	});
}(this["Sunlight"]));