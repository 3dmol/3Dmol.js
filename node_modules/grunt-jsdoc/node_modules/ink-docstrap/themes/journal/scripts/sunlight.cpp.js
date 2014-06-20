(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}
	
	var primitives = ["int", "char", "void", "long", "signed", "unsigned", "double", "bool", "typename", "class", "short", "wchar_t", "struct"],
		acceptableKeywords = ["int", "char", "void", "long", "signed", "unsigned", "double", "bool", "typename", "class", "short", "wchar_t"];

	sunlight.registerLanguage("cpp", {
		//http://www.cppreference.com/wiki/keywords/start
		keywords: [
			"and","default","noexcept","template","and_eq","delete","not","this","alignof","double",
			"not_eq","thread_local","asm","dynamic_cast","nullptr","throw","auto","else","operator",
			"true","bitand","enum","or","try","bitor","explicittodo","or_eq","typedef","bool","export",
			"private","typeid","break","externtodo","protected","typename","case","false","public","union",
			"catch","float","register","using","char","for","reinterpret_cast","unsigned","char16_t",
			"friend","return","void","char32_t","goto","short","wchar_t","class","if","signed","virtual",
			"compl","inline","sizeof","volatile","const","int","static","while","constexpr","long",
			"static_assert","xor","const_cast","mutable","static_cast","xor_eq","continue","namespace",
			"struct","decltype","new","switch"
		],
		
		customTokens: {
			constant: {
				values: [
					"EXIT_SUCCESS", "EXIT_FAILURE",
					"SIG_DFL", "SIG_IGN", "SIG_ERR", "SIGABRT", "SIGFPE", "SIGILL", "SIGINT", "SIGSEGV", "SIGTERM"
				],
				boundary: "\\b"
			},
			
			//http://www.cppreference.com/wiki/utility/types/start
			basicType: {
				values: ["ptrdiff_t", "size_t", "nullptr_t", "max_align_t"],
				boundary: "\\b"
			},
			
			ellipsis: {
				values: ["..."],
				boundary: ""
			}
		},

		scopes: {
			string: [ ["\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])] ],
			"char": [ ["'", "'", ["\\\'", "\\\\"]] ],
			comment: [ ["//", "\n", null, true], ["/*", "*/"] ],
			preprocessorDirective: [ ["#", "\n", null, true] ]
		},
		
		customParseRules: [	
		],

		identFirstLetter: /[A-Za-z_]/,
		identAfterFirstLetter: /\w/,

		namedIdentRules: {
			custom: [
				//pointer default declarations, e.g. pointer* myPointer;
				function() {
					var precedes = [[
							sunlight.util.whitespace, 
							{ token: "operator", values: ["*", "**"] }, 
							{ token: "default" },
							{ token: "ident" }, 
							sunlight.util.whitespace, 
							{ token: "punctuation", values: [";"] }
						], [
							//function parameters
							{ token: "default" },
							{ token: "operator", values: ["&"] },
							sunlight.util.whitespace,
							{ token: "ident" }
							
						]
					];
					
					return function(context) {
						//basically, can't be on the right hand side of an equals sign
						//so we traverse the tokens backward, and if we run into a "=" before a ";" or a "{", it's no good
						
						var token,
							index,
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
						
						//make sure we're not on the left side of the equals sign
						index = context.index;
						while (token = context.tokens[--index]) {
							if (token.name === "punctuation" && (token.value === ";" || token.value === "{")) {
								return true;
							}
							
							if (token.name === "operator" && token.value === "=") {
								return false;
							}
						}
						
						return false;
					};
				}(),
				
				//casting
				function() {
					var precedes = [
						[sunlight.util.whitespace, { token: "punctuation", values: [")"] }, sunlight.util.whitespace, { token: "ident" }],
						[{ token: "operator", values: ["*", "**"] }, sunlight.util.whitespace, { token: "punctuation", values: [")"] }, sunlight.util.whitespace, { token: "operator", values: ["&"], optional: true }, { token: "ident" }]
					];
				
					return function(context) {
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
					};
				}(),
			
				//generic definitions/params between "<" and ">"
				function(context) {
					//between < and > and preceded by an ident and not preceded by "class"
					var index = context.index,
						token,
						prevToken = sunlight.util.getPreviousNonWsToken(context.tokens, context.index),
						foundIdent = false,
						bracketCountLeft;
					
					//if the previous token is a keyword, then we don't care about it
					if (!prevToken || prevToken.name === "keyword") {
						return false;
					}
					
					//look for "<" preceded by an ident but not "class"
					//if we run into ">" before "," or "<" then it's a big fail
					bracketCountLeft = [0, 0];
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
									//allows generic method invocations, like "Foo" in "foo.Resolve<Foo>()"
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
				},
				
				//ident before generic definitions, e.g. "foo" in "foo<bar>"
				function(context) {
					//if it's preceded by an ident or a primitive/alias keyword then it's no good (i.e. a generic method definition like "public void Foo<T>")
					//also a big fail if it is preceded by a ., i.e. a generic method invocation like container.Resolve()
					var token = sunlight.util.getPreviousNonWsToken(context.tokens, context.index),
						index,
						bracketCount;
					
					if (token !== undefined) {
						if (
							token.name === "ident" 
							|| (token.name === "keyword" && sunlight.util.contains(primitives.concat(["string", "object", "void"]), token.value))
							|| (token.name === "operator" && token.value === ".")
						) {
							return false;
						}
					}
					
					//needs to be immediately followed by <, then by idents, acceptable keywords and ",", and then closed by >, then immediately followed by an ident
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
								default:
									return false;
							}
							
							//if bracket counts match, get the f out
							if (bracketCount[0] === bracketCount[1]) {
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
				},
				
				//after class keyword but inside <>
				function(context) {
					var prevToken = sunlight.util.getPreviousNonWsToken(context.tokens, context.index),
						token,
						index;
					
					if (!prevToken || prevToken.name !== "keyword" || prevToken.value !== "class") {
						return false;
					}
					
					//make sure we're not inside <>
					//easiest way is to go forward and verify that we hit a "{" before a ">"
					index = context.index;
					while (token = context.tokens[++index]) {
						if (token.name === "punctuation" && token.value === "{") {
							return true;
						}
						
						if (token.name === "operator" && sunlight.util.contains([">", ">>"], token.value)) {
							return false;
						}
					}
					
					return false;
				}
			],
			
			follows: [
				[{ token: "keyword", values: ["enum", "struct", "union"] }, sunlight.util.whitespace]
			],
			
			precedes: [
				//normal parameters/declarations
				[{ token: "default" }, { token: "ident" }],
				
				[sunlight.util.whitespace, { token: "operator", values: ["*", "**"] }, { token: "default" }, { token: "ident" }, sunlight.util.whitespace, { token: "operator", values: ["=", ","] }],
				[sunlight.util.whitespace, { token: "operator", values: ["*", "**"] }, { token: "default" }, { token: "operator", values: ["&"] }, sunlight.util.whitespace, { token: "ident" }, sunlight.util.whitespace, { token: "operator", values: ["=", ","] }],
				
				//e.g. "std" in "std::char_traits<CharT>"
				[sunlight.util.whitespace, { token: "operator", values: ["::"] }]
			]
		},
		
		//http://www.cppreference.com/wiki/language/operator_precedence
		operators: [
			"==", "=",
			"+=", "++", "+",
			"->*", "->", "-=", "--", "-",
			"**", "*=", "*", //added ** for double pointer convenience
			"/=", "/",
			"%=", "%",
			"!=", "!",
			">>=", ">>", ">=", ">",
			"<<=", "<<", "<=", "<",
			"&=", "&&", "&",
			"|=", "||", "|",
			"~",
			"^=", "^",
			".*", ".",
			"?",
			"::", ":",
			","
		]
		
	});
}(this["Sunlight"]));