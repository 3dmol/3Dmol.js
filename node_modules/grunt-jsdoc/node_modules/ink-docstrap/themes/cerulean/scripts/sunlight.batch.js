(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}

	sunlight.registerLanguage("batch", {
		caseInsensitive: true,
		
		scopes: {
			string: [ ["\"", "\"", ["\\\"", "\\\\"]] ],
			comment: [ ["REM", "\n", null, true], ["::", "\n", null, true] ],
			variable: [["%", { regex: /[^\w%]/, length: 1}, null, true]]
		},
		
		customParseRules: [
			//labels
			function(context) {
				var colon,
					peek,
					value = "",
					line,
					column;
				
				if (!context.reader.isSolWs() || context.reader.current() !== ":" || context.reader.peek() === ":") {
					return null;
				}
				
				colon = context.createToken("operator", ":", context.reader.getLine(), context.reader.getColumn());
				
				//label, read until whitespace
				while (peek = context.reader.peek()) {
					if (/\s/.test(peek)) {
						break;
					}
					
					value += context.reader.read();
					if (!line) {
						line = context.reader.getLine();
						column = context.reader.getColumn();
					}
				}
				
				if (value === "") {
					return null;
				}
				
				
				return [colon, context.createToken("label", value, line, column)];
			},
			
			//label after goto statements
			function(context) {
				var matches = sunlight.util.createProceduralRule(context.count() - 1, -1, [{ token: "keyword", values: ["goto"] }, { token: "operator", values: [":"], optional: true }], true),
					peek,
					value,
					line = context.reader.getLine(),
					column = context.reader.getColumn();
				
				if (!matches(context.getAllTokens())) {
					return null;
				}
				
				value = context.reader.current();
				while (peek = context.reader.peek()) {
					if (/[\W]/.test(peek)) {
						break;
					}
					
					value += context.reader.read();
				}
				
				return context.createToken("label", value, line, column);
			},
			
			//keywords have to be handled manually because strings don't have to be quoted
			//e.g. we don't want to highlight "do" in "echo do you have the time?"
			function() {
				var keywords = sunlight.util.createHashMap([
					//commands
					"assoc","attrib","break","bcdedit","cacls","call","cd","chcp","chdir","chkdsk","chkntfs","cls","cmd",
					"color","comp","compact","convertfcopy","date","del","dir","diskcomp","diskcopy","diskpart","doskey",
					"driverquery","echo","endlocal","erase","exit","fc","findstr","find","format","for","fsutil","ftype",
					"goto","gpresult","graftabl","help","icacls","if","label","md","mkdir","mklink","mode","more","move",
					"openfiles","path","pause","popd","print","prompt","pushd","rd","recover",/*"rem",*/"rename","ren",
					"replace","rmdir","robocopy","setlocal","set","schtasks","sc","shift","shutdown","sort","start",
					"subst","systeminfo","tasklist","taskkill","time","title","tree","type","verify","ver","vol","xcopy",
					"wmic",
					
					"lfnfor",
					
					//keywords
					"do", "else", "errorlevel", "exist", "in", "not",
					"choice",
					"com1", "con", "prn", "aux", "nul", "lpt1",
					"exit", "eof", "off", "on",
					
					"equ","neq","lss","leq","gtr","geq"
				], "\\b", true);
				
				return function(context) {
					var token = sunlight.util.matchWord(context, keywords, "keyword", true),
						prevToken,
						index;
					
					if (!token) {
						return null;
					}
					
					//look backward for "echo" or "title" or "set" or "|" or beginning of line
					//if we find "echo" or "set" or "title" or "=" before "|" or sol then it's a fail
					
					if (!context.reader.isSolWs()) {
						index = context.count();
						while (prevToken = context.token(--index)) {
							if (prevToken.name === "keyword" && sunlight.util.contains(["echo", "title", "set"], prevToken.value)) {
								return null;
							}
							if (prevToken.name === "operator" && prevToken.value === "=") {
								return null;
							}
							
							//pipe
							if (prevToken.name === "operator" && prevToken.value === "|") {
								break;
							}
							
							//sol
							if (prevToken.name === "default" && prevToken.value.indexOf("\n") >= 0) {
								break;
							}
						}
					}
					
					context.reader.read(token.value.length - 1);
					return token;
				};
			}()
		],
		
		identFirstLetter: /[A-Za-z_\.]/,
		identAfterFirstLetter: /[\w-]/,
		
		operators: [
			"&&", "||", "&", ":", "/", "==", "|", "@", "*", ">>", ">", "<", "==!", "!", "=", "+"
		]

	});
}(this["Sunlight"]));