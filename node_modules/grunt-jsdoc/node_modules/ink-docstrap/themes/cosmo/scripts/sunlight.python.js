(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}

	sunlight.registerLanguage("python", {
		keywords: [
			//http://docs.python.org/py3k/reference/lexical_analysis.html#keywords
			"False","class","finally","is","return","None","continue","for","lambda","try","True",
			"def","from","nonlocal","while","and","del","global","not","with","as","elif","if","or",
			"yield","assert","else","import","pass","break","except","in","raise"
		],
		
		customTokens: {
			ellipsis: {
				values: ["..."],
				boundary: ""
			},
			
			//http://docs.python.org/py3k/reference/lexical_analysis.html#delimiters
			delimiter: {
				values: ["(", ")", "[", "]", "{", "}", ",", ":", ".", ";", "@", "=", "+=", "-=", "*=", "/=", "//=", "%=", "&=", "|=", "^=", ">>=", "<<=", "**="],
				boundary: ""
			},
			
			//http://docs.python.org/py3k/library/constants.html
			constant: {
				values: ["NotImplemented", "Ellipsis", "False", "True", "None", "__debug__"],
				boundary: "\\b"
			},
			
			attribute: {
				values: [
					"__doc__", "__name__", "__module__", "__defaults__", "__code__", "__globals__", "__dict__", "__closure__", "__annotations__", "__kwedefaults__",
					"__self__", "__func__", "__bases__"
				],
				boundary: "\\b"
			},
			
			//http://docs.python.org/py3k/reference/datamodel.html#specialnames
			specialMethod: {
				values: [
					"__next__", "__new__", "__init__", "__del__", "__repr__", "__str__", "__format__",
					"__lt__", "__le__", "__eq__", "__ne__", "__gt__", "__ge__",
					"__hash__", "__bool__",
					"__call__", "__prepare__",
					"__getattr__", "__getattribute__", "__setattr__", "__setattribute__", "__delattr__", "__dir__",
					"__get__", "__set__", "__delete__",
					"__slots__",
					"__instancecheck__", "__subclasscheck__",
					"__getitem__", "__setitem__", "__delitem__", "__iter__", "__reversed__", "__contains__",
					"__add__", "__sub__", "__mul__", "__truediv__", "__floordiv__", "__mod__", "__divmod__", "__pow__", "__lshift__", "__rshift__", "__and__", "__xor__", "__or__",
					"__radd__", "__rsub__", "__rmul__", "__rtruediv__", "__rfloordiv__", "__rmod__", "__rdivmod__", "__rpow__", "__rlshift__", "__rrshift__", "__rand__", "__xror__","__ror__",
					"__iadd__", "__isub__", "__imul__", "__itruediv__", "__ifloordiv__", "__imod__", "__idivmod__", "__ipow__", "__ilshift__", "__irshift__", "__iand__", "__xror__", "__ior__",
					"__neg__", "__pos__", "__abs__", "__invert__", "__complex__", "__int__", "__float__", "__round__", "__index__",
					"__enter__", "__exit__"
				],
				boundary: "\\b"
			},
			
			//http://docs.python.org/py3k/library/functions.html
			"function": {
				values: [
					"abs","dict","help","min","setattr","all","dir","hex","next","slice","any","divmod","id","object","sorted","ascii","enumerate","input",
					"oct","staticmethod","bin","eval","int","open","str","bool","exec","isinstance","ord","sum","bytearray","filter","issubclass","pow",
					"super","bytes","float","iter","print","tuple","callable","format","len","property","type","chr","frozenset","list","range","vars",
					"classmethod","getattr","locals","repr","zip","compile","globals","map","reversed","__import__","complex","hasattr","max","round",
					"delattr","hash","memoryview","set"
				],
				boundary: "\\b"
			}
		},

		scopes: {
			longString: [ 
				["\"\"\"", "\"\"\"", sunlight.util.escapeSequences.concat(["\\\""])], 
				["'''", "'''", sunlight.util.escapeSequences.concat(["\\'"])] 
			],
			rawLongString: [ 
				["r\"\"\"", "\"\"\""], 
				["R\"\"\"", "\"\"\""], 
				["r'''", "'''"], 
				["R'''", "'''"]
			],
			binaryLongString: [
				//raw binary
				["br\"\"\"", "\"\"\""], 
				["bR\"\"\"", "\"\"\""], 
				["Br\"\"\"", "\"\"\""], 
				["BR\"\"\"", "\"\"\""],
				["br'''", "'''"], 
				["bR'''", "'''"], 
				["Br'''", "'''"], 
				["BR'''", "'''"],
				
				//just binary
				["b\"\"\"", "\"\"\"", sunlight.util.escapeSequences.concat(["\\\""])], 
				["B\"\"\"", "\"\"\"", sunlight.util.escapeSequences.concat(["\\\""])],
				["b'''", "'''", sunlight.util.escapeSequences.concat(["\\'"])], 
				["B'''", "'''", sunlight.util.escapeSequences.concat(["\\'"])]
			],
			
			string: [ ["\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])], ["'", "'", sunlight.util.escapeSequences.concat(["\\\'", "\\\\"])] ],
			rawString: [ ["r\"", "\""], ["R\"", "\""], ["r'", "'"], ["R'", "'"] ],
			binaryString: [ 
				//just binary
				["b\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])],
				["b'", "'", sunlight.util.escapeSequences.concat(["\\'"])],
				["B\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])],
				["B'", "'", sunlight.util.escapeSequences.concat(["\\'"])],
				
				//raw binary
				["br\"", "\""], ["bR\"", "\""], ["Br\"", "\""], ["BR\"", "\""],
				["br'", "'"], ["bR'", "'"], ["Br'", "'"], ["BR'", "'"]
			],
			
			comment: [ ["#", "\n", null, true] ]
		},
		
		identFirstLetter: /[A-Za-z_]/,
		identAfterFirstLetter: /\w/,
		
		namedIdentRules: {
			follows: [
				//class names
				//function names
				//exception names
				[{ token: "keyword", values: ["class", "def", "raise", "except"] }, sunlight.util.whitespace]
			]
		},

		operators: [
			//http://docs.python.org/py3k/reference/lexical_analysis.html#operators
			"+", 
			"-", 
			"**", "*",
			"//",  "/",  
			"%", 
			"&", 
			"|",
			"^", 
			"~", 
			"<<", "<=", "<", 
			">>", ">=", ">",   
			"==", "!="
		]
	});
}(this["Sunlight"]));