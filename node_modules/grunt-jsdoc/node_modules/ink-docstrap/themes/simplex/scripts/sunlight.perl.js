(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}
	
	function isValidForRegexLiteral(context) {
		var previousNonWsToken = context.token(context.count() - 1),
			previousToken = null;
		
		if (context.defaultData.text !== "") {
			previousToken = context.createToken("default", context.defaultData.text); 
		}
		
		if (!previousToken) {
			previousToken = previousNonWsToken;
		}
		
		//first token of the string
		if (previousToken === undefined) {
			return true;
		}
		
		if (sunlight.util.contains(["keyword", "ident", "number", "variable", "specialVariable"], previousNonWsToken.name)) {
			//this is not valid context for a regex literal
			return false;
		}
		if (previousNonWsToken.name === "punctuation" && !sunlight.util.contains(["(", "{", "[", ",", ";"], previousNonWsToken.value)) {
			return false;
		}
		
		return true;
	}
	
	function readWhileWhitespace(context) {
		var value = "", 
			peek;
		
		while ((peek = context.reader.peek()) !== context.reader.EOF && /\s/.test(peek)) {
			value += context.reader.read();
		}
		
		return value;
	}
	
	function readBetweenDelimiters(context, delimiter, appendCloser) {
		var opener = delimiter,
			closer = delimiter,
			trackCloser = sunlight.util.contains(["[", "(", "{"], opener),
			peek2,
			value = opener,
			next,
			closerCount = 1
		
		switch (delimiter) {
			case "[":
				closer = "]";
				break;
			case "(":
				closer = ")";
				break;
			case "{":
				closer = "}";
				break;
		}
		
		while (context.reader.peek() !== context.reader.EOF) {
			peek2 = context.reader.peek(2);
			if (peek2 === "\\" + closer || peek2 === "\\\\") {
				//escaped backslash or escaped closer
				value += context.reader.read(2);
				continue;
			}
			
			next = context.reader.read();
			
			if (next === opener && trackCloser) {
				closerCount++;
			} else if (next === closer && --closerCount <= 0) {
				if (trackCloser || appendCloser) {
					value += next;
				}
				break;
			}

			value += next;
		}
		
		return value;
	}
	
	//perl allows whitespace before delimiters (wtf?)
	function getDelimiterIfValid(context, peekCount) {
		var peek = context.reader.peek(peekCount);
		while (peek.length === peekCount && /\s$/.test(peek)) {
			peek = context.reader.peek(++peekCount);
		}
		
		if (/[\w]$/.test(peek)) {
			return false;
		}
		
		context.reader.read(peekCount);
		return {
			delimiter: context.reader.current(),
			value: peek.substring(0, peek.length - 1)
		};
	}

	sunlight.registerLanguage("perl", {
		keywords: [
			"caller","die","dump","eval","exit","goto","last","next","redo","return","sub","wantarray",
			"break","continue","given","when","default",
			"import","local","my","our","state",
			"do","no","package","require","use",
			"bless","dbmclose","dbmopen","ref","tied","untie","tie",
			
			"if", "elsif", "else", "unless", "while", "foreach", "for", "until",
			
			"not", "or", "and"
		],
		
		customTokens: {
			"function": {
				values: [
					"chomp","chop","chr","crypt","hex","index","length","oct","ord","rindex","sprintf","substr",
					"pos","quotemeta","split","study",
					"abs","atan2","cos","exp","hex","int","log","oct","rand","sin","sqrt","srand",
					"pop","push","shift","splice","unshift",
					"grep","join","map","reverse","sort",
					"delete","each","exists","keys","values",
					"binmode","closedir","close","eof","fileno","flock","format","getc","print","printf","readdir","rewinddir",
					"say","seekdir","seek","select","syscall","sysread","sysseek","tell","telldir","truncate","warn","write",
					"pack","syswrite","unpack","vec",

					"chdir","chmod","chown","chroot","fcntl","glob","ioctl","link","lstat","mkdir","open","opendir","readlink","rename","rmdir","stat","symlink","sysopen","umask","unlink","utime",
					"defined","dump","eval","formline","reset","scalar","undef",
					"alarm","exec","fork","getpgrp","getppid","getpriority","kill","pipe","setpgrp","setpriority","sleep","system","wait","waitpid",
					"accept","bind","connect","getpeername","getsockname","getsockopt","listen","recv","send","setsockopt","shutdown","socket","socketpair",
					"msgctl","msgget","msgrcv","msgsnd","semctl","semget","semop","shmctl","shmget","shmread","shmwrite",
					"endgrent","endhostent","endnetent","endpwent","getgrent","getgrgid","getgrnam","getlogin","getpwent","getpwnam","getpwuid","setgrent","setpwent",
					"endprotoent","endservent","gethostbyaddr","gethostbyname","gethostent","getnetbyaddr","getnetbyname","getnetent","getprotobyname","getprotobynumber",
					"getprotoent","getservbyname","getservbyport","getservent","sethostent","setnetent","setprotoent","setservent",
					"gmtime","localtime","times","time",

					"lcfirst","lc","lock","prototype","readline",
					"readpipe","read","ucfirst","uc"
				],
				boundary: "\\b"
			},
			
			//http://perldoc.perl.org/perlvar.html
			//jesus, perl...
			specialVariable: {
				values: [
					"$.", "$<", "$_", "$/", "$!", "$ARG", "$&", "$a", "$b", "$MATCH", "$PREMATCH", "${^MATCH}", "${^PREMATCH}", "$POSTMATCH", "$'", 
					"$LAST_PAREN_MATCH", "$+", "$LAST_SUBMATCH_RESULT", "$^N", "$INPUT_LINE_NUMBER", "$NR", "$.", "$INPUT_RECORD_SEPARATOR", "$RS",
					"$OUTPUT_AUTOFLUSH", "$OFS", "$,", 
					
					"@LAST_MATCH_END", "@+", 
					
					"%LAST_PAREN_MATCH", "%+",
					
					"$OUTPUT_RECORD_SEPARATOR","$ORS","$LIST_SEPARATOR","$\"","$SUBSCRIPT_SEPARATOR","$SUBSEP","$;","$FORMAT_PAGE_NUMBER","$%",
					"$FORMAT_LINES_PER_PAGE","$=","$FORMAT_LINES_LEFT","$-","@LAST_MATCH_START","@-","%-","$FORMAT_NAME","$~","$FORMAT_TOP_NAME",
					"$FORMAT_LINE_BREAK_CHARACTERS","$:","$FORMAT_FORMFEED","$^L","$ACCUMULATOR","$^A","$CHILD_ERROR","$?","${^CHILD_ERROR_NATIVE}",
					"${^ENCODING}","$OS_ERROR","$ERRNO","$!","%OS_ERROR","%ERRNO","%!","$EXTENDED_OS_ERROR","$^E","$EVAL_ERROR","$@","$PROCESS_ID","$PID",
					"$$","$REAL_USER_ID","$UID","$<","$EFFECTIVE_USER_ID","$EUID","$>","$REAL_GROUP_ID","$GID","$(","$EFFECTIVE_GROUP_ID","$EGID","$)",
					"$PROGRAM_NAME","$0","$[","$]","$COMPILING","$^C","$DEBUGGING","$^D","${^RE_DEBUG_FLAGS}","${^RE_TRIE_MAXBUF}","$SYSTEM_FD_MAX","$^F",
					"$^H","%^H","$INPLACE_EDIT","$^I","$^M","$OSNAME","$^O","${^OPEN}","$PERLDB","$^P","$LAST_REGEXP_CODE_RESULT","$^R","$EXCEPTIONS_BEING_CAUGHT",
					"$^S","$BASETIME","$^T","${^TAINT}","${^UNICODE}","${^UTF8CACHE}","${^UTF8LOCALE}","$PERL_VERSION","$^V","$WARNING","$^W","${^WARNING_BITS}",
					"${^WIN32_SLOPPY_STAT}","$EXECUTABLE_NAME","$^X","ARGV","$ARGV","@ARGV","ARGVOUT","@F","@INC","@ARG","@_","%INC","%ENV","$ENV","%SIG","$SIG",
					
					"$^",
					
					"$#array"
				],
				boundary: "\\W"
			}
		},

		scopes: {
			string: [ ["\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])], ["'", "'", ["\\\'", "\\\\"]] ],
			comment: [ ["#", "\n", null, true] ],
			variable: [ 
				["$#", { length: 1, regex: /[\W]/ }, null, true], //array count
				["$", { length: 1, regex: /[\W]/ }, null, true], 
				["@", { length: 1, regex: /[\W]/ }, null, true], 
				["%", { length: 1, regex: /[\W]/ }, null, true] 
			]
		},
		
		customParseRules: [
			//qr/STRING/msixpo, m/PATTERN/msixpogc, /PATTERN/msixpogc, // (empty pattern), ?pattern?
			//y///, tr///, s/PATTERN/REPLACEMENT/msixpogce 
			function(context) {
				var current,
					peek,
					delimiter,
					hasReplace = false,
					value = "",
					delimiterInfo,
					line = context.reader.getLine(),
					column = context.reader.getColumn(),
					fetchSecondDelimiter;
				
				if (!isValidForRegexLiteral(context)) {
					return null;
				}
				
				current = context.reader.current();
				peek = context.reader.peek();
				
				if (current === "/" || current === "?") {
					delimiter = current;
				} else if (current === "m" || current === "y" || current === "s") {
					if (!(delimiterInfo = getDelimiterIfValid(context, 1))) {
						return null;
					}
					
					value = current + delimiterInfo.value;
					delimiter = delimiterInfo.delimiter;
					hasReplace = current === "y" || current === "s";
				} else if ((current === "t" || current === "q") && peek === "r") {
					if (!(delimiterInfo = getDelimiterIfValid(context, 2))) {
						return null;
					}
					
					hasReplace = current === "t";
					value = current + delimiterInfo.value;
					delimiter = delimiterInfo.delimiter;
				} else {
					return null;
				}
				
				//read the regex literal
				value += readBetweenDelimiters(context, delimiter, !hasReplace);
				if (hasReplace) {
					//apparently whitespace between search and replace is allowed, so read the whitespace, if it exists
					value += readWhileWhitespace(context);
					//new delimiter
					
					fetchSecondDelimiter = sunlight.util.contains(["[", "(", "{"], delimiter);
					if (fetchSecondDelimiter) {
						delimiterInfo = getDelimiterIfValid(context, 1);
						if (delimiterInfo) {
							value += delimiterInfo.value;
							delimiter = delimiterInfo.delimiter;
						}
					}
					
					value += readBetweenDelimiters(context, delimiter, true);
				}
				
				//read the regex modifiers (we just assume any character is valid)
				while (context.reader.peek() !== context.reader.EOF) {
					if (!/[A-Za-z]/.test(context.reader.peek())) {
						break;
					}
					
					value += context.reader.read();
				}
				
				return context.createToken("regexLiteral", value, line, column);
			},
			
			//raw strings
			function(context) {
				var value = "q",
					readCount = 1,
					peek,
					line = context.reader.getLine(), 
					column = context.reader.getColumn();
					
				//begin with q, qw, qx, or qq  with a non-alphanumeric delimiter (opening bracket/paren are closed by corresponding closing bracket/paren)
				if (context.reader.current() !== "q") {
					return null;
				}
				
				peek = context.reader.peek();
				if (peek === "q" || peek === "w" || peek == "x") {
					readCount++;
				}
				
				if (/[A-Za-z0-9]$/.test(context.reader.peek(readCount))) {
					//potential % operator
					return null;
				}
				
				value += context.reader.read(readCount - 1) + readBetweenDelimiters(context, context.reader.read(), true);
				return context.createToken("rawString", value, line, column);
			},
			
			//heredoc declaration (stolen from ruby)
			function(context) {
				var prevToken,
					line = context.reader.getLine(), 
					column = context.reader.getColumn(),
					value = "<<",
					ident = "",
					current,
					delimiter = "",
					peek,
					peek2;
					
				if (context.reader.current() !== "<" || context.reader.peek() !== "<") {
					return null;
				}
				
				//cannot be preceded by an ident or a number or a string
				prevToken = sunlight.util.getPreviousNonWsToken(context.getAllTokens(), context.count() - 1);
				if (prevToken && (prevToken.name === "ident" || prevToken.name === "number" || prevToken.name === "string")) {
					return null;
				}
				
				//can be between quotes (double, single or back) or not, or preceded by a hyphen
				
				context.reader.read(2);
				
				current = context.reader.current();
				if (current === "-") {
					context.reader.read();
					value += current;
					current = context.reader.current();
				}
				
				if (sunlight.util.contains(["\"", "'", "`"], current)) {
					delimiter = current;
				} else {
					ident = current;
				}
				
				value += current;
				
				while ((peek = context.reader.peek()) !== context.reader.EOF) {
					if (peek === "\n" || (delimiter === "" && /\W/.test(peek))) {
						break;
					}
					
					if (peek === "\\") {
						peek2 = context.reader.peek(2);
						if (delimiter !== "" && sunlight.util.contains(["\\" + delimiter, "\\\\"], peek2)) {
							value += peek2;
							ident += context.reader.read(2);
							continue;
						}
					}
					
					value += context.reader.read();
					
					if (delimiter !== "" && peek === delimiter) {
						break;
					}
					
					ident += peek;
				}
				
				context.items.heredocQueue.push(ident);
				
				return context.createToken("heredocDeclaration", value, line, column);
			},
			
			//heredoc
			function(context) {
				var tokens = [],
					declaration,
					line,
					column,
					value,
					peekIdent;
				
				if (context.items.heredocQueue.length === 0) {
					return null;
				}
				
				//there must have been at least one line break since the heredoc declaration(s)
				if (context.defaultData.text.replace(/[^\n]/g, "").length === 0) {
					return null;
				}
				
				//we're confirmed to be in the heredoc body, so read until all of the heredoc declarations have been satisfied
				
				value = context.reader.current();
				while (context.items.heredocQueue.length > 0 && context.reader.peek() !== context.reader.EOF) {
					declaration = context.items.heredocQueue.shift();
					line = context.reader.getLine(), column = context.reader.getColumn();
					
					//read until "\n{declaration}\n"
					while (context.reader.peek() !== context.reader.EOF) {
						peekIdent = context.reader.peek(declaration.length + 2);
						if (peekIdent === "\n" + declaration || peekIdent === "\n" + declaration + "\n") {
							value += context.reader.read(declaration.length + 2);
							break;
						}
						
						value += context.reader.read();
					}
					
					tokens.push(context.createToken("heredoc", value, line, column));
					value = "";
				}
				
				return tokens.length > 0 ? tokens : null;
			},
			
			//pod: http://perldoc.perl.org/perlpod.html
			//stolen from ruby
			function(context) {
				var value = "=",
					line = context.reader.getLine(),
					column = context.reader.getColumn(),
					foundEnd = false,
					peek;
				
				//these begin on with a line that starts with "=begin" and end with a line that starts with "=end"
				//apparently stuff on the same line as "=end" is also part of the comment
				
				if (context.reader.current() !== "=" || !context.reader.isSol()) {
					return null;
				}
				
				//read until "\n=cut" and then everything until the end of that line
				while ((peek = context.reader.peek()) !== context.reader.EOF) {
					if (!foundEnd && context.reader.peek(5) === "\n=cut") {
						foundEnd = true;
						value += context.reader.read(5);
						continue;
					}
					
					if (foundEnd && peek === "\n") {
						break;
					}
					
					value += context.reader.read();
				}
				
				return context.createToken("docComment", value, line, column);
			}
		],

		identFirstLetter: /[A-Za-z_]/,
		identAfterFirstLetter: /\w/,

		namedIdentRules: {
			follows: [
				[{ token: "keyword", values: ["sub"] }, { token: "default" }],
				[{ token: "operator", values: ["\\&"] }, sunlight.util.whitespace]
			]
		},
		
		operators: [
			"++", "+=", "+",
			"--", "-=", "-",
			"**=", "**", "*=", "*",
			"//=", "/=", "//", "/",
			"%=", "%",
			"=>", "=~", "==", "=",
			"!", "!~", "!=",
			"~", "~~",
			"\\&", "\\",
			
			"&&=", "&=", "&&", "&", 
			"||=", "||", "|=", "|",
			
			"<<=", "<=>", "<<", "<=", "<",
			">>=", ">>", ">=", ">",
			"^=", "^",
			
			"?", "::", ":",
			
			"...", ".=", "..", ".",
			
			",",
			
			"x=", "x", //seriously, perl?
			
			"lt", "gt", "le", "ge", "eq", "ne", "cmp"
		],
		
		contextItems: {
			heredocQueue: []
		}
		
	});
}(this["Sunlight"]));