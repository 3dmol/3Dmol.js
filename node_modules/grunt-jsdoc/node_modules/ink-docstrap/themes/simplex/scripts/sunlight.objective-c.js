(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}
	
	sunlight.registerLanguage("objective-c", {
		keywords: [
			//c++ keywords
			"and","default","noexcept","template","and_eq","delete","not","this","alignof","double",
			"not_eq","thread_local","asm","dynamic_cast","nullptr","throw","auto","else","operator",
			"true","bitand","enum","or","try","bitor","explicittodo","or_eq","typedef","bool","export",
			"private","typeid","break","externtodo","protected","typename","case","false","public","union",
			"catch","float","register","using","char","for","reinterpret_cast","unsigned","char16_t",
			"friend","return","void","char32_t","goto","short","wchar_t","if","signed","virtual",
			"compl","inline","sizeof","volatile","const","int","static","while","constexpr","long",
			"static_assert","xor","const_cast","mutable","static_cast","xor_eq","continue","namespace",
			"struct","decltype","new","switch",
			
			//objective c keywords
			"id", "self", "nil", "super", "in", "out", "inout", "bycopy", "byval", "oneway",
			
			"SEL", "BOOL", "YES", "NO",
			
			"@interface", "@implementation", "@end", "@class",
			"@private", "@public", "@package", "@protected",
			"@protocol", "@optional", "@required",
			"@property", "@synthesize", "@dynamic",
			"@selector",
			"@try", "@catch", "@finally", "@throw",
			"@synchronized",
			"@encode",
			
			"__attribute__", 
			
			//these seem to be conditional, somehow...
			"__weak", "__strong"
		],
		
		customTokens: {
			constant: {
				values: [
					"EXIT_SUCCESS", "EXIT_FAILURE",
					"SIG_DFL", "SIG_IGN", "SIG_ERR", "SIGABRT", "SIGFPE", "SIGILL", "SIGINT", "SIGSEGV", "SIGTERM"
				],
				boundary: "\\b"
			}
		},

		scopes: {
			string: [ ["\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])], ["@\"", "\"", ["\\\\", "\\\""]] ],
			"char": [ ["'", "'", ["\\\'", "\\\\"]] ],
			comment: [ ["//", "\n", null, true], ["/*", "*/"] ],
			preprocessorDirective: [ ["#", "\n", null, true] ]
		},
		
		customParseRules: [
			//message destination (e.g. method calls)
			function(context) {
				var peek,
					count,
					ident,
					possibleMessageArgument,
					match,
					parenCount,
					bracketCount,
					token,
					index,
					exprCount;
				
				//read the ident first
				if (!context.language.identFirstLetter.test(context.reader.current())) {
					return null;
				}
				
				count = 0;
				while ((peek = context.reader.peek(++count)) && peek.length === count) {
					if (!context.language.identAfterFirstLetter.test(sunlight.util.last(peek))) {
						break;
					}
				}
				
				ident = context.reader.current() + peek.substring(0, peek.length - 1);
				count = count - 1;
				possibleMessageArgument = false;
				while ((peek = context.reader.peek(++count)) && peek.length === count) {
					if (!/\s$/.test(peek)) {
						match = /([\]:])$/.exec(peek);
						if (match === null) {
							//not a message destination
							return null;
						}
						
						possibleMessageArgument = match[1] === ":" && !/::$/.test(context.reader.peek(count + 1));
						break;
					}
				}
				
				//must be the second expression after "["
				parenCount = 0;
				bracketCount = 0;
				index = context.count();
				exprCount = 1;
				while (token = context.token(--index)) {
					if (exprCount > 1 && !possibleMessageArgument) {
						return null;
					}
					
					if (token.name === "punctuation") {
						switch (token.value) {
							case ";":
							case "{":
							case "}":
								//short circuit rules
								return null;
							case "(":
								parenCount--;
								break;
							case ")":
								parenCount++;
								break;
							case "[":
								if (bracketCount === 0 && parenCount === 0) {
									if (exprCount >= 1) {
										token = context.createToken(possibleMessageArgument && exprCount > 1 ? "messageArgumentName" : "messageDestination", ident, context.reader.getLine(), context.reader.getColumn());
										context.reader.read(ident.length - 1);
										return token;
									}
									
									return null;
								}
								
								bracketCount--;
								break;
							case "]":
								bracketCount++;
								break;
						}
					}
					
					if (bracketCount === 0 && parenCount === 0 && token.name === "default") {
						exprCount++;
					}
				}
				
				return null;
			},
			
			//@property attributes
			function() {
				var attributes = sunlight.util.createHashMap([
					"getter", "setter", 
					"readonly", "readwrite",
					"assign", "retain", "copy",
					"nonatomic"
				], "\\b");
				
				return function(context) {
					var token = sunlight.util.matchWord(context, attributes, "keyword", true),
						prevToken,
						index;
					if (!token) {
						return null;
					}
					
					//must be inside () after @property
					
					//look backward for "("
					//if we find a ";" before a "(" then that's no good
					index = context.count();
					while (prevToken = context.token(--index)) {
						if (prevToken.name === "punctuation") {
							if (prevToken.value === "(") {
								//previous token must be @property
								prevToken = sunlight.util.getPreviousNonWsToken(context.getAllTokens(), index);
								if (!prevToken || prevToken.name !== "keyword" || prevToken.value !== "@property") {
									return null;
								}
								
								token.line = context.reader.getLine();
								token.column = context.reader.getColumn();
								context.reader.read(token.value.length - 1);
								return token;
							} else if (prevToken.value === ";") {
								return null;
							}
						}
					}
					
					return null;
				};
			}()
		],

		identFirstLetter: /[A-Za-z_]/,
		identAfterFirstLetter: /\w/,

		//after classname in () (categories)
		namedIdentRules: {
			custom: [
				//naming convention: NS.+, CG.+ are assumed to be built in objects
				function(context) {
					var regex = /^(NS|CG).+$/,
						nextToken = sunlight.util.getNextNonWsToken(context.tokens, context.index);
					return regex.test(context.tokens[context.index].value) && (!nextToken || nextToken.name !== "punctuation" || nextToken.value !== "(");
				},
				
				//call to class or alloc
				function(context) {
					var nextToken = sunlight.util.getNextNonWsToken(context.tokens, context.index);
					return nextToken && nextToken.name === "messageDestination" && (nextToken.value === "class" || nextToken.value === "alloc");
				},
				
				//ident followed by an ident, but not inside [] and not followed by ":"
				function(context) {
					var token,
						index,
						parenCount;

					if (!sunlight.util.createProceduralRule(context.index + 1, 1, [
						{ token: "default" },
						{ token: "ident" }
					])(context.tokens)) {
						return false;
					}

					if (sunlight.util.createProceduralRule(context.index + 1, 1, [
						{ token: "default" },
						{ token: "ident" },
						sunlight.util.whitespace,
						{ token: "operator", value: ":" }
					])(context.tokens)) {
						//should not be followed by a colon, as that indicates this is an argument definition
						return false;
					}
					
					//must be between []
					index = context.index;
					parenCount = 0;
					while (token = context.tokens[--index]) {
						if (token.name === "punctuation") {
							switch (token.value) {
								case "[":
									return false;
								case "{":
								case ",":
									return true;
								case "(":
									if (parenCount === 0) {
										return true;
									}

									parenCount++;
									break;
								case ")":
									parenCount--;
									break;
							}
						}
					}
					
					return true;
				},
				
				//pointer default declarations, e.g. pointer* myPointer;
				function() {
					var precedes = [[
							sunlight.util.whitespace,
							{ token: "operator", values: ["*", "**"] }, 
							sunlight.util.whitespace,
							{ token: "ident" }, 
							sunlight.util.whitespace, 
							{ token: "punctuation", values: [";"] }
						], [
							//function parameters
							sunlight.util.whitespace,
							{ token: "operator", values: ["&", "*", "**"] },
							sunlight.util.whitespace,
							{ token: "ident" }
							
						]
					];
					
					return function(context) {
						//basically, can't be on the right hand side of an equals sign
						//so we traverse the tokens backward, and if we run into a "=" before a ";" or a "{", it's no good
						
						var precedesIsSatisfied,
							isPartOfProperty,
							foundEquals,
							token,
							index;
						
						precedesIsSatisfied = function(tokens) {
							var i;
							for (i = 0; i < precedes.length; i++) {
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
						//objc addition: okay if part of a @property statement
						isPartOfProperty = false;
						foundEquals = false;
						index = context.index;
						while (token = context.tokens[--index]) {
							if (token.name === "punctuation" && (token.value === ";" || token.value === "{")) {
								return isPartOfProperty || !foundEquals;
							}
							
							if (token.name === "operator" && token.value === "=") {
								foundEquals = true;
							} else if (token.name === "keyword" && token.value === "@property") {
								isPartOfProperty = true;
							}
						}
						
						return false;
					};
				}(),
				
				//casting
				function() {
					var precedes = [
						[sunlight.util.whitespace, { token: "punctuation", values: [")"] }, sunlight.util.whitespace, { token: "ident" }],
						[sunlight.util.whitespace, { token: "punctuation", values: [")"] }, sunlight.util.whitespace, { token: "punctuation", values: ["["] }],
						[
							sunlight.util.whitespace, 
							{ token: "operator", values: ["*", "**"] }, 
							sunlight.util.whitespace, 
							{ token: "punctuation", values: [")"] }, 
							sunlight.util.whitespace, 
							{ token: "operator", values: ["&"], optional: true }, 
							{ token: "ident" }
						],
						
						[
							sunlight.util.whitespace, 
							{ token: "operator", values: ["*", "**"] }, 
							sunlight.util.whitespace, 
							{ token: "punctuation", values: [")"] }, 
							sunlight.util.whitespace, 
							{ token: "operator", values: ["&"], optional: true }, 
							{ token: "punctuation", values: ["["] }
						]
					];
				
					return function(context) {
						var token, 
							index,
							prevToken,
							precedesIsSatisfied;
							
						precedesIsSatisfied = function(tokens) {
							var i;
							for (i = 0; i < precedes.length; i++) {
								if (sunlight.util.createProceduralRule(context.index + 1, 1, precedes[i], false)(tokens)) {
									return true;
								}
							}
							
							return false;
						}(context.tokens);
						
						if (!precedesIsSatisfied) {
							return false;
						}
						
						//make sure the previous tokens are "(" and then not a keyword or an ident
						//this'll make sure that things like "if (foo) doSomething();" and "bar(foo)" won't color "foo"
						
						index = context.index;
						while (token = context.tokens[--index]) {
							if (token.name === "punctuation" && token.value === "(") {
								prevToken = sunlight.util.getPreviousNonWsToken(context.tokens, index);
								if (prevToken) {
									if (prevToken.name === "ident") {
										return false;
									}
									
									if (prevToken.name === "keyword" && sunlight.util.contains(["if", "while"], prevToken.value)) {
										return false;
									}
								}
								
								return true;
							}
						}
						
						return false;
					};
				}(),
				
				//generic definitions/params between "<" and ">"
				//stolen and slightly modified from cpp, this is actually for protocols, since objective-c doesn't have generics
				function(context) {
					//between < and > and preceded by an ident and not preceded by "class"
					var index = context.index, 
						token, 
						foundIdent, 
						bracketCountLeft, 
						prevToken;
					
					//if the previous token is a keyword, then we don't care about it
					prevToken = sunlight.util.getPreviousNonWsToken(context.tokens, context.index);
					if (!prevToken || prevToken.name === "keyword") {
						return false;
					}
					
					//look for "<" preceded by an ident but not "class"
					//if we run into ">" before "," or "<" then it's a big fail
					foundIdent = false;
					bracketCountLeft = [0, 0];
					while ((token = context.tokens[--index]) !== undefined) {
						if (token.name === "operator") {
							switch (token.value) {
								case "<":
								case "<<":
									bracketCountLeft[0] += token.value.length;
									continue;
								case ">":
								case ">>":
									if (bracketCountLeft[0] === 0) {
										return false;
									}
									
									bracketCountLeft[1] += token.value.length;
									continue;
								case ".":
								case "::":
									//allows generic method invocations, like "Foo" in "foo.Resolve<Foo>()"
								case "*":
									//allows pointers
									continue;
							}
						}
						
						if (
							//(token.name === "keyword" && sunlight.util.contains(acceptableKeywords, token.value))
							token.name === "default"
							|| (token.name === "punctuation" && token.value === ",")
						) {
							continue;
						}
						
						if (token.name === "ident" || (token.name === "keyword" && sunlight.util.contains(["id", "static_cast"], token.value))) {
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
							//(token.name === "keyword" && sunlight.util.contains(acceptableKeywords, token.value))
							(token.name === "operator" && sunlight.util.contains(["<", "<<", ">", ">>", "::", "*"], token.value))
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
				
				//ident before <>
				//stolen from c++/java/c#
				function(context) {
					//if it's preceded by an ident or a primitive/alias keyword then it's no good (i.e. a generic method definition like "public void Foo<T>")
					//also a big fail if it is preceded by a ., i.e. a generic method invocation like container.Resolve()
					var token = sunlight.util.getPreviousNonWsToken(context.tokens, context.index),
						index,
						bracketCount;
						
					if (token !== undefined) {
						if (
							token.name === "ident" 
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
				}
			],
			
			follows: [
				[{ token: "keyword", values: ["@interface", "@protocol", "@implementation"] }, { token: "default" }]
			],
			
			precedes: [
				[{ token: "operator", values: ["::"] }],
				
				[
					sunlight.util.whitespace, 
					{ token: "operator", values: ["*", "**"] },
					{ token: "default" }, 
					{ token: "ident" }, 
					sunlight.util.whitespace, 
					{ token: "operator", values: ["=", ","] }
				],
				
				[
					sunlight.util.whitespace, 
					{ token: "operator", values: ["*", "**"] }, 
					sunlight.util.whitespace, 
					{ token: "operator", values: ["&"] }, 
					sunlight.util.whitespace, 
					{ token: "ident" }, 
					sunlight.util.whitespace, 
					{ token: "operator", values: ["=", ","] }
				]
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