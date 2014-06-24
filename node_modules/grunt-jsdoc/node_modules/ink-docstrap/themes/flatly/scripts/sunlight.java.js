(function(sunlight, undefined){

	//this is remarkably similar to the C# language definition
	//ergo, there is a lot of copypasting going on here...

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}
	
	function createNamedIdentFunction(func) {
		var typeDefinitionRegex = /^T([A-Z0-9]\w*)?$/;
		return function(context) {
			return !typeDefinitionRegex.test(context.tokens[context.index].value) && func(context);
		};
	};
	
	var primitives = ["boolean", "byte", "char", "double", "float", "int", "long", "short"],
		acceptableKeywords = primitives.concat(["extends"]);

	sunlight.registerLanguage("java", {
		keywords: [
			//http://download.oracle.com/javase/tutorial/java/nutsandbolts/_keywords.html
			"abstract", "assert", "boolean", "break", "byte", "case", "catch", "char", "class", "const",
			"continue", "default", "do", "double", "else", "enum", "extends", "final", "finally", "float",
			"for", "goto", "if", "implements", "import", "instanceof", "int", "interface", "long", "native",
			"new", "package", "private", "protected", "public", "return", "short", "static", "strictfp" /* wtf? */, "super",
			"switch", "synchronized", "this", "throw", "throws", "transient", "try", "void", "volatile", "while",
			
			//literals
			"null", "false", "true"
		],
		
		scopes: {
			string: [ ["\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])], ["'", "'", ["\'", "\\\\"]] ],
			comment: [ ["//", "\n", null, true], ["/*", "*/"] ],
			annotation: [ ["@", { length: 1, regex: /[\s\(]/ }, null, true] ]
		},

		identFirstLetter: /[A-Za-z_]/,
		identAfterFirstLetter: /\w/,

		//these are mostly stolen from the C# lang file
		namedIdentRules: {
			custom: [
				//generic definitions/params between "<" and ">"
				createNamedIdentFunction(function(context) {
					//between < and > and preceded by an ident and not preceded by "class"
					var index = context.index, 
						token,
						foundIdent = false,
						bracketCountLeft = [0, 0];
						
					//look for "<" preceded by an ident but not "class"
					//if we run into ">" before "," or "<" then it's a big fail
					while ((token = context.tokens[--index]) !== undefined) {
						if (token.name === "keyword" && token.value === "class") {
							//this must be a generic class type definition, e.g. Foo<T>, and we don't want to color the "T"
							return false;
						}
						
						if (token.name === "operator") {
							switch (token.value) {
								case "<":
								case "<<":
									bracketCountLeft[0] += token.value.length;
									break;
								case ">":
								case ">>":
									if (bracketCountLeft[0] === 0) {
										return false;
									}
									
									bracketCountLeft[1] += token.value.length;
									break;
								case ".":
								case "?":
								case "&":
									//allows generic method invocations, like "Foo" in "foo.Resolve<Foo>()"
									//and <? extends Bar & Baz>
									break;
								default:
									return false;
							}
							
							continue;
						}
						
						if (
							(token.name === "keyword" && sunlight.util.contains(acceptableKeywords, token.value))
							|| token.name === "default"
							|| (token.name === "punctuation" && token.value === ",")
						) {
							continue;
						}
						
						if (token.name === "ident") {
							foundIdent = true;
							continue;
						}
						
						//anything else means we're no longer in a generic definition
						break;
					}
					
					if (!foundIdent || bracketCountLeft[0] === 0) {
						//not inside a generic definition
						return false;
					}
					
					//now look forward to make sure the generic definition is closed
					//this avoids false positives like "foo < bar"
					index = context.index;
					while ((token = context.tokens[++index]) !== undefined) {
						if (token.name === "operator" && (token.value === ">" || token.value === ">>")) {
							return true;
						}
						
						if (
							(token.name === "keyword" && sunlight.util.contains(acceptableKeywords, token.value))
							|| (token.name === "operator" && sunlight.util.contains(["<", "<<", ">", ">>"], token.value))
							|| (token.name === "punctuation" && token.value === ",")
							|| token.name === "ident"
							|| token.name === "default"
						) {
							continue;
						}
						
						return false;
					}
					
					return false;
				}),
				
				//generic declarations and return values (ident preceding a generic definition)
				//this finds "Foo" in "Foo<Bar> foo"
				createNamedIdentFunction(function(context) {
					//if it's preceded by an ident or a primitive/alias keyword then it's no good (i.e. a generic method definition like "public void Foo<T>")
					//also a big fail if it is preceded by a ., i.e. a generic method invocation like container.Resolve()
					var token = sunlight.util.getPreviousNonWsToken(context.tokens, context.index),
						index,
						bracketCount;
					
					if (token !== undefined) {
						if (
							token.name === "ident" 
							|| (token.name === "keyword" && sunlight.util.contains(primitives.concat(["void"]), token.value))
							|| (token.name === "operator" && token.value === ".")
						) {
							return false;
						}
					}
					
					//needs to be immediately followed by <, then by idents or ?, acceptable keywords and ",", and then closed by >, then immediately followed by an ident
					token = sunlight.util.getNextNonWsToken(context.tokens, context.index);
					if (!token || token.name !== "operator" || token.value !== "<") {
						return false;
					}
					
					index = context.index;
					bracketCount = [0, 0]; //open (<), close (>)
					while ((token = context.tokens[++index]) !== undefined) {
						if (token.name === "operator") {
							switch (token.value) {
								case "<":
									bracketCount[0]++;
									break;
								case "<<":
									bracketCount[0] += 2;
									break;
								case ">":
									bracketCount[1]++;
									break;
								case ">>":
									bracketCount[1] += 2;
									break;
								case "?":
								case "&":
									//e.g. Foo<? extends Bar & Baz>
									break;
								default:
									return false;
							}
							
							//if bracket counts match, get the f out
							if (bracketCount[0] > 0 && bracketCount[0] === bracketCount[1]) {
								break;
							}
							
							continue;
						}
						
						if (
							token.name === "default"
							|| token.name === "ident"
							|| (token.name === "keyword" && sunlight.util.contains(acceptableKeywords, token.value))
							|| (token.name === "punctuation" && token.value === ",")
						) {
							continue;
						}
						
						return false;
					}
					
					//verify bracket count
					if (bracketCount[0] !== bracketCount[1]) {
						//mismatched generics, could be something scary
						return false;
					}
					
					//next token should be optional whitespace followed by an ident
					token = context.tokens[++index];
					if (!token || (token.name !== "default" && token.name !== "ident")) {
						return false;
					}
					
					if (token.name === "default") {
						token = context.tokens[++index];
						if (!token || token.name !== "ident") {
							return false;
						}
					}
					
					return true;
				}),
				
				//fully qualified type names
				createNamedIdentFunction(function(context) {
					//next token is not "."
					var nextToken = sunlight.util.getNextNonWsToken(context.tokens, context.index),
						token,
						index,
						previous;
					
					if (nextToken && nextToken.name === "operator" && nextToken.value === ".") {
						return false;
					}
					
					//go backward and make sure that there are only idents and dots before the new keyword
					index = context.index;
					previous = context.tokens[index];
					while ((token = context.tokens[--index]) !== undefined) {
						if (token.name === "keyword" && (token.value === "new" || token.value === "import" || token.value === "instanceof")) {
							return true;
						}
						
						if (token.name === "default") {
							continue;
						}
						
						if (token.name === "ident") {
							if (previous && previous.name === "ident") {
								return false;
							}
							
							previous = token;
							continue;
						}
						
						if (token.name === "operator" && token.value === ".") {
							if (previous && previous.name !== "ident") {
								return false;
							}
							
							previous = token;
							continue;
						}
						
						break;
					}
					
					return false;
				}),
			
				//casting
				function() {
					var precedes = [
						[sunlight.util.whitespace, { token: "punctuation", values: [")"] }, sunlight.util.whitespace, { token: "ident" }],
						[sunlight.util.whitespace, { token: "punctuation", values: [")"] }, sunlight.util.whitespace, { token: "keyword", values: ["this"] }]
					];
				
					return createNamedIdentFunction(function(context) {
						var token,
							index,
							prevToken,
							precedesIsSatisfied = function(tokens) {
								for (var i = 0; i < precedes.length; i++) {
									if (sunlight.util.createProceduralRule(context.index + 1, 1, precedes[i], false)(tokens)) {
										return true;
									}
								}
								
								return false;
							}(context.tokens);
						
						if (!precedesIsSatisfied) {
							return false;
						}
						
						//make sure the previous tokens are "(" and then not a keyword
						//this'll make sure that things like "if (foo) doSomething();" won't color "foo"
						
						index = context.index;
						while (token = context.tokens[--index]) {
							if (token.name === "punctuation" && token.value === "(") {
								prevToken = sunlight.util.getPreviousNonWsToken(context.tokens, index);
								if (prevToken && prevToken.name === "keyword") {
									return false;
								}
								
								return true;
							}
						}
						
						return false;
					});
				}(),
			
				//can't use the follows/precedes/between utilities since we need to verify that it doesn't match the type definition naming convention
				createNamedIdentFunction(function(context) {
					var i,
						follows = [
							[{ token: "ident" }, sunlight.util.whitespace, { token: "keyword", values: ["extends", "implements"] }, sunlight.util.whitespace],
					
							//method/property return values
							//class/interface names
							[{ token: "keyword", values: ["class", "interface", "enum", "public", "private", "protected", "static", "final"] }, sunlight.util.whitespace],
							
							//bounded generic interface constraints
							//this matches "MyInterface" in "public <T extends Object & MyInterface>..."
							[
								{ token: "keyword", values: ["extends"] },
								{ token: "default" }, 
								{ token: "ident" }, 
								{ token: "default" }, 
								{ token: "operator", values: ["&"] },
								{ token: "default" }
							]
						],
						
						precedes = [
							//arrays
							[sunlight.util.whitespace, { token: "punctuation", values: ["["] }, sunlight.util.whitespace, { token: "punctuation", values: ["]"] }], //in method parameters

							//assignment: Object o = new object();
							//method parameters: public int Foo(Foo foo, Bar b, Object o) { }
							[{ token: "default" }, { token: "ident" }]
						],
						
						between = [
							{ opener: { token: "keyword", values: ["implements", "throws"] }, closer: { token: "punctuation", values: ["{"] } }
						];
					
					for (i = 0; i < follows.length; i++) {
						if (sunlight.util.createProceduralRule(context.index - 1, -1, follows[i], false)(context.tokens)) {
							return true;
						}
					}
					
					for (i = 0; i < precedes.length; i++) {
						if (sunlight.util.createProceduralRule(context.index + 1, 1, precedes[i], false)(context.tokens)) {
							return true;
						}
					}
					
					for (i = 0; i < between.length; i++) {
						if (sunlight.util.createBetweenRule(context.index, between[i].opener, between[i].closer, false)(context.tokens)) {
							return true;
						}
					}
					
					return false;
				})
			]
		},

		operators: [
			//arithmetic
			"++", "+=", "+",
			"--", "-=", "-",
			      "*=", "*",
			      "/=", "/",
			      "%=", "%",

			//boolean
			"&&", "||",

			//bitwise
			"|=",   "|",
			"&=",   "&",
			"^=",   "^",
			">>>=", ">>>", ">>=", ">>",
			"<<=", "<<",

			//inequality
			"<=", "<",
			">=", ">",
			"==", "!=",

			//unary
			"!", "~",

			//other
			"?", ":", ".", "="
		]
	});
}(this["Sunlight"]));