(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}
	
	if (!sunlight.isRegistered("xml")) {
		throw "Scala requires the XML language to be registered";
	}

	sunlight.registerLanguage("scala", {
		keywords: [
			"abstract","case","catch","class","def","do","else","extends","false","final","finally","forSome","for",
			"if","implicit","import","lazy","match","new","null","object","override","package","private","protected",
			"return","sealed","super","this","throw","trait","try","true","type","val","var","while","with","yield"
		],
		
		embeddedLanguages: {
			xml: {
				switchTo: function(context) {
					var prevToken;
					if (context.reader.current() !== "<" || !/[\w!?]/.test(context.reader.peek())) {
						return false;
					}
					
					if (context.defaultData.text !== "") {
						//preceded by whitespace
						return true;
					}
					
					prevToken = context.token(context.count() - 1);
					return prevToken && prevToken.name === "punctuation" && sunlight.util.contains(["(", "{"], prevToken.value);
				},
				
				switchBack: function(context) {
					var prevToken = context.token(context.count() - 1);
					if (!prevToken) {
						return false;
					}
					
					if (prevToken.name === "tagName") {
						if (!context.items.literalXmlOpenTag) {
							context.items.literalXmlOpenTag = prevToken.value;
						}
					} else if (prevToken.name === "operator") {
						switch (prevToken.value) {
							case "<":
								context.items.literalXmlNestingLevel++;
								break;
							case "</":
							case "/>":
								context.items.literalXmlNestingLevel--;
								break;
						}
					}
					
					if (context.items.literalXmlOpenTag && context.items.literalXmlNestingLevel === 0 && (prevToken.value === ">" || prevToken.value === "/>")) {
						return true;
					}
					
					return false;
				}
			}
		},
		
		scopes: {
			string: [ ["\"\"\"", "\"\"\""], ["\"", "\"", ["\\\\", "\\\""]] ],
			"char": [ ["'", "'", ["\\\\", "\\'"]] ],
			quotedIdent: [ ["`", "`", ["\\`", "\\\\"]] ],
			comment: [ ["//", "\n", null, true], ["/*", "*/"] ],
			annotation: [ ["@", { length: 1, regex: /\W/ }, null, true] ]
		},
		
		identFirstLetter: /[A-Za-z]/,
		identAfterFirstLetter: /\w/,
		
		customParseRules: [
			//symbol literals
			function(context) {
				var line = context.reader.getLine(),
					column = context.reader.getColumn();

				if (context.reader.current() !== "'") {
					return false;
				}

				var match = /^(\w+)(?!')/i.exec(context.reader.peekSubstring());
				if (!match) {
					return false;
				}

				context.reader.read(match[1].length);
				return context.createToken("symbolLiteral", "'" + match[1], line, column);
			},

			//case classes: can't distinguish between a case class and a function call so we need to keep track of them
			function(context) {
				var prevToken,
					peek,
					ident = context.reader.current(), 
					line = context.reader.getLine(), 
					column = context.reader.getColumn();
				
				if (context.defaultData.text === "") {
					return false;
				}
				
				if (!/[A-Za-z]/.test(context.reader.current())) {
					return false;
				}
				
				prevToken = context.token(context.count() - 1);
				if (context.defaultData.text === "" || !prevToken || prevToken.name !== "keyword" || !sunlight.util.contains(["class", "type", "trait", "object"], prevToken.value)) {
					return false;
				}
				
				//read the ident
				while (peek = context.reader.peek()) {
					if (!/\w/.test(peek)) {
						break;
					}
					
					ident += context.reader.read();
				}
				
				context.items.userDefinedTypes.push(ident);
				return context.createToken("ident", ident, line, column);
			}
		],

		namedIdentRules: {
			custom: [	
				//some built in types
				function() {
					var builtInTypes = [
						"Nil", "Nothing", "Unit", "Pair", "Map", "String", "List", "Int",
						
						"Seq", "Option", "Double", "AnyRef", "AnyVal", "Any", "ScalaObject",
						"Float", "Long", "Short", "Byte", "Char", "Boolean"
					];
					
					return function(context) {
						//next token is not "."
						var nextToken = sunlight.util.getNextNonWsToken(context.tokens, context.index);
						if (nextToken && nextToken.name === "operator" && nextToken.value === ".") {
							return false;
						}
					
						return sunlight.util.contains(builtInTypes, context.tokens[context.index].value);
					};
				}(),
				
				//user-defined types
				function(context) {
					return sunlight.util.contains(context.items.userDefinedTypes, context.tokens[context.index].value);
				},
				
				//fully qualified type names after "new"
				function(context) {
					//next token is not "."
					var nextToken = sunlight.util.getNextNonWsToken(context.tokens, context.index),
						token,
						index,
						previous;
					
					if (nextToken && nextToken.name === "operator" && nextToken.value === ".") {
						return false;
					}
					
					//go backward and make sure that there are only idents and dots after the new keyword
					index = context.index;
					previous = context.tokens[index];
					while ((token = context.tokens[--index]) !== undefined) {
						if (token.name === "keyword" && token.value === "new") {
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
				},
				
				function() {
					var follows = [
							[{ token: "keyword", values: ["class", "object", "extends", "new", "type", "trait"] }, { token: "default" }],
							[{ token: "operator", values: [":"] }, sunlight.util.whitespace],
							[{ token: "operator", values: ["#"] }],
							[{ token: "keyword", values: ["type"] }, { token: "default" }, { token: "ident" }, sunlight.util.whitespace, { token: "operator", values: ["="] }, sunlight.util.whitespace]
						],
						between = [
							//generics
							{ opener: { token: "punctuation", values: ["["] }, closer: { token: "punctuation", values: ["]"] } }
						];
					
					return function(context) {
						var i;
						if (/^[A-Z]([A-Z0-9]\w*)?$/.test(context.tokens[context.index].value)) {
							//generic type names are assumed to start with a capital letter optionally followed by a number or another capital letter
							//e.g. A, T1, TFrom, etc.
							return false;
						}
						
						for (i = 0; i < follows.length; i++) {
							if (sunlight.util.createProceduralRule(context.index - 1, -1, follows[i], false)(context.tokens)) {
								return true;
							}
						}
						
						for (i = 0; i < between.length; i++) {
							if (sunlight.util.createBetweenRule(context.index, between[i].opener, between[i].closer, false)(context.tokens)) {
								return true;
							}
						}
						
						return false;
					};
				}()
				
			]
		},
		
		contextItems: {
			literalXmlOpenTag: null,
			literalXmlNestingLevel: 0,
			userDefinedTypes: []
		},

		operators: [
			"++", "+=", "+",
			"--", "-=", "->", "-",
			"*=", "*",
			"^=", "^^", "^",
			"~>", "~",
			"!=", "!",
			"&&", "&=", "&",
			"||", "|=", "|",
			
			">>>", ">>=", ">>", ">",
			"<<=", "<<", "<~", "<=", "<%", "<",
			
			"%>", "%=", "%",
			
			"::", ":<", ":>", ":",
			
			"==", "=",
		
			"@", "#", "_", "."
		]

	});
	
	sunlight.globalOptions.enableScalaXmlInterpolation = false;
	
	sunlight.bind("beforeHighlight", function(context) {
		if (context.language.name === "scala") {
			this.options.enableScalaXmlInterpolation = true;
		}
	});
	sunlight.bind("afterHighlight", function(context) {
		this.options.enableScalaXmlInterpolation = false;
	});
	
}(this["Sunlight"]));