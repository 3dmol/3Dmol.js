(function(sunlight, undefined){

	if (sunlight === undefined || sunlight["registerLanguage"] === undefined) {
		throw "Include sunlight.js before including language files";
	}
	
	sunlight.registerLanguage("nginx", {
		scopes: {
			string: [ ["\"", "\"", sunlight.util.escapeSequences.concat(["\\\""])], ["'", "'", ["\\'", "\\\\"]] ],
			comment: [ ["#", "\n", null, true] ],
			variable: [ ["$", { length: 1, regex: /[\W]/ }, null, true] ],
			label: [ ["@", { length: 1, regex: /[\W]/ }, null, true] ],
			ssiCommand: [ ["<!--#", "-->"] ] //http://wiki.nginx.org/HttpSsiModule#SSI_Commands
		},
		
		customParseRules: [
			//context detection (server, location, upstream, http, mail, etc.)
			function() {
				var words = sunlight.util.createHashMap(["server", "location", "upstream", "http", "mail", "types", "map", "split-clients", "geo", "limit_except"], "\\b", false);
				
				return function(context) {
					var token = sunlight.util.matchWord(context, words, "context", true),
						peek,
						count;
					
					if (!token) {
						return null;
					}
					
					//if we encounter a "{" before a ";", we're good to go
					count = token.value.length;
					while ((peek = context.reader.peek(count++)) !== context.reader.EOF) {
						if (/\{$/.test(peek)) {
							context.reader.read(token.value.length - 1); //already read the first letter
							return token;
						}
						
						if (/;$/.test(peek)) {
							break;
						}
					}
					
					return null;
				};
			}(),
			
			//regex literals
			function(context) {
				//can come after server_name if it starts with ~ and ends with [\s;]
				//or after "location ^~"
				//or after "location ~"
				//or after "location ~*"
				
				var current = context.reader.current(),
					token,
					index,
					isRegexLiteral = false,
					regexLiteral = current,
					peek,
					prevToken,
					line = context.reader.getLine(), 
					column = context.reader.getColumn();
					
				if (/[\s\n]/.test(current)) {
					return null;
				}
				
				index = context.count() - 1;
				while ((token = context.token(index--)) !== context.reader.EOF) {
					if (token.name === "regexLiteral" || (token.name === "punctuation" && (token.value === "{" || token.value === "}" || token.value === ";"))) {
						return null;
					}
					
					if (token.name === "operator" && sunlight.util.contains(["^~", "~", "~*"], token.value)) {
						//check previous token for "location" keyword
						prevToken = sunlight.util.getPreviousWhile(context.getAllTokens(), index, function(token) {
							return token.name === "default" || token.name === "comment";
						});
						
						if (prevToken.name === "context" && prevToken.value === "location") {
							isRegexLiteral = true;
							break;
						}
					} else if (token.name === "keyword" && token.value === "server_name") {
						//server_name check
						if (current !== "~") {
							return null;
						}
						
						isRegexLiteral = true;
						break;
					}
				}
				
				if (!isRegexLiteral) {
					return null;
				}
				
				//read to the end of the regex literal
				while ((peek = context.reader.peek()) !== context.reader.EOF) {
					if (/[\s;\n]/.test(peek)) {
						break;
					}
					
					regexLiteral += context.reader.read();
				}
				
				return context.createToken("regexLiteral", regexLiteral, line, column);
			},
			
			//directives, can't just make them keywords because then values (not directives) will possibly be highlighted
			//e.g. "index index.html" will have both "index"es highlighted when we just want the first
			function() {
				var directives = sunlight.util.createHashMap([
					//core: http://wiki.nginx.org/CoreModule
					"daemon", "env", "debug_points", "error_log", "include", "lock_file", "master_process", "pid", "ssl_engine", "timer_resolution",
					"user", "worker_cpu_affinity", "worker_priority", "worker_processes", "worker_rlimit_core", "worker_rlimit_nofile",
					"worker_rlimit_sigpending", "working_directory",
					
					//events: http://wiki.nginx.org/EventsModule
					"accept_mutext_delay", "accept_mutex", "debug_connection", "devpoll_changes", "devpoll_events", "kqueue_changes",
					"kqueue_events", "epoll_events", "multi_accept", "rtsig_signo", "rtsig_overflow_events", "rtsig_overflow_text",
					"rtsig_overflow_threshold","use", "worker_connections",
					
					//httpCore: http://wiki.nginx.org/HttpCoreModule
					"aio","alias","chunked_transfer_encoding","client_body_in_file_only","client_body_in_single_buffer","client_body_buffer_size",
					"client_body_temp_path","client_body_timeout","client_header_buffer_size","client_header_timeout","client_max_body_size",
					"default_type","directio","error_page","if_modified_since","internal","keepalive_timeout","keepalive_requests",
					"large_client_header_buffers","limit_except","limit_rate_after","limit_rate","listen","location","log_not_found",
					"log_subrequest","msie_padding","msie_refresh","open_file_cache_errors","open_file_cache_min_uses",
					"open_file_cache_valid","open_file_cache","optimize_server_names","port_in_redirect","post_action","recursive_error_pages",
					"resolver_timeout","resolver","root","satisfy_any","satisfy","send_timeout","sendfile","server_name","server_name_in_redirect",
					"server_names_hash_max_size","server_names_hash_bucket_size","server_tokens","server","tcp_nodelay","tcp_nopush","try_files","types",
					"underscores_in_headers",
					
					//httpAccess: http://wiki.nginx.org/HttpAccessModule
					"allow", "deny",
					
					//httpAuthBasic: http://wiki.nginx.org/HttpAuthBasicModule
					"auth_basic_user_file", "auth_basic",
					
					//httpAutoindex: http://wiki.nginx.org/HttpAutoindexModule
					"autoindex_exact_size", "autoindex_localtime", "autoindex",
					
					//httpBrowser: http://wiki.nginx.org/HttpBrowserModule
					"ancient_browser_value", "ancient_browser", "modern_browser_value","modern_browser",
					
					//httpCharset: http://wiki.nginx.org/HttpCharsetModule
					"charset_map", "override_charset", "source_charset", "charset",
					
					//emptyGif: http://wiki.nginx.org/HttpEmptyGifModule
					"empty_gif",
					
					//fCgi: http://wiki.nginx.org/HttpFcgiModule
					"fastcgi_bind","fastcgi_buffer_size","fastcgi_buffers","fastcgi_cache_key","fastcgi_cache_path",
					"fastcgi_cache_methods","fastcgi_cache_min_uses","fastcgi_cache_use_stale","fastcgi_cache_valid","fastcgi_cache","fastcgi_connect_timeout",
					"fastcgi_index","fastcgi_hide_header","fastcgi_ignore_client_abort","fastcgi_ignore_headers","fastcgi_intercept_errors",
					"fastcgi_max_temp_file_size","fastcgi_no_cache","fastcgi_next_upstream","fastcgi_param","fastcgi_pass_header","fastcgi_pass",
					"fastcgi_read_timeout","fastcgi_redirect_errors","fastcgi_send_timeout","fastcgi_split_path_info",
					"fastcgi_store_access","fastcgi_store","fastcgi_temp_path",
					
					//geo: http://wiki.nginx.org/HttpGeoModule
					"geo",
					
					//gzip: http://wiki.nginx.org/HttpGzipModule
					"gzip_buffers", "gzip_comp_level", "gzip_disable", "gzip_http_version", "gzip_min_length", "gzip_proxied", "gzip_types", "gzip_vary", "gzip",
					
					//headers: http://wiki.nginx.org/HttpHeadersModule
					"add_header", "expires",
					
					//index: http://wiki.nginx.org/HttpIndexModule
					"index",
					
					//limitRequests: http://wiki.nginx.org/HttpLimitReqModule
					"limit_req_log_level", "limit_req_zone", "limit_req",
					
					//limitZone: http://wiki.nginx.org/HttpLimitZoneModule
					"limit_zone", "limit_conn_log_level", "limit_conn",
					
					//log: http://wiki.nginx.org/HttpLogModule
					"access_log", "log_format", "open_log_file_cache",
					
					//map: http://wiki.nginx.org/HttpMapModule
					"map_hash_max_size", "map_hash_bucket_size", "map",
					
					//memcached: http://wiki.nginx.org/HttpMemcachedModule
					"memcached_pass", "memcached_connect_timeout", "memcached_read_timeout", "memcached_send_timeout", "memcached_buffer_size", "memcached_next_upstream",
					
					//proxy: http://wiki.nginx.org/HttpProxyModule
					"proxy_bind","proxy_buffer_size","proxy_buffering","proxy_buffers","proxy_busy_buffers_size","proxy_cache_bypass","proxy_cache_key",
					"proxy_cache_methods","proxy_cache_min_uses","proxy_cache_path","proxy_cache_use_stale","proxy_cache_valid","proxy_cache","proxy_connect_timeout",
					"proxy_headers_hash_bucket_size","proxy_headers_hash_max_size","proxy_hide_header","proxy_ignore_client_abort","proxy_ignore_headers",
					"proxy_intercept_errors","proxy_max_temp_file_size","proxy_method","proxy_next_upstream","proxy_no_cache","proxy_pass_header",
					"proxy_pass_request_body","proxy_pass_request_headers","proxy_pass","proxy_read_timeout","proxy_redirect_errors","proxy_redirect","proxy_send_lowat",
					"proxy_send_timeout","proxy_set_body","proxy_set_header","proxy_ssl_session_reuse","proxy_store_access","proxy_store","proxy_temp_file_write_size",
					"proxy_temp_path","proxy_upstream_fail_timeout","proxy_upstream_max_fails",
					
					//referer: http://wiki.nginx.org/HttpRefererModule
					"valid_referers",
					
					//rewrite: http://wiki.nginx.org/HttpRewriteModule
					"break", "if", "return", "rewrite", "set", "uninitialized_variable_warn",
					
					//scgi: http://wiki.nginx.org/HttpScgiModule
					"scgi_bind","scgi_buffer_size","scgi_buffers","scgi_busy_buffers_size","scgi_cache_bypass","scgi_cache_key","scgi_cache_methods",
					"scgi_cache_min_uses","scgi_cache_path","scgi_cache_use_stale","scgi_cache_valid","scgi_cache","scgi_connect_timeout","scgi_hide_header","scgi_ignore_client_abort",
					"scgi_ignore_headers","scgi_intercept_errors","scgi_max_temp_file_size","scgi_next_upstream","scgi_no_cache","scgi_param","scgi_pass_header",
					"scgi_pass_request_body","scgi_pass_request_headers","scgi_pass","scgi_read_timeout","scgi_send_timeout","scgi_store_access","scgi_store","scgi_temp_file_write_size",
					"scgi_temp_path",
					
					//clients: http://wiki.nginx.org/HttpSplitClientsModule
					"split-clients", //is this a typo in the docs? should it be an underscore?
					
					//ssi: http://wiki.nginx.org/HttpSsiModule
					"ssi", "ssi_silent_errors", "ssi_types", "ssi_value_length",
					
					//upstream: http://wiki.nginx.org/HttpUpstreamModule
					"ip_hash", "server", "upstream",
					
					//userid: http://wiki.nginx.org/HttpUserIdModule
					"userid_domain", "userid_expires", "userid_name", "userid_p3p", "userid_path", "userid_service", "userid",
					
					//uwsgi: http://wiki.nginx.org/HttpUwsgiModule
					"uwsgi_bind","uwsgi_buffer_size","uwsgi_buffers","uwsgi_busy_buffers_size","uwsgi_cache_bypass","uwsgi_cache_key",
					"uwsgi_cache_methods","uwsgi_cache_min_uses","uwsgi_cache_path","uwsgi_cache_use_stale","uwsgi_cache_valid","uwsgi_cache","uwsgi_connect_timeout",
					"uwsgi_hide_header","uwsgi_ignore_client_abort","uwsgi_ignore_headers","uwsgi_intercept_errors","uwsgi_max_temp_file_size","uwsgi_modifier1",
					"uwsgi_modifier2","uwsgi_next_upstream","uwsgi_no_cache","uwsgi_param","uwsgi_pass_header","uwsgi_pass_request_body",
					"uwsgi_pass_request_headers","uwsgi_pass","uwsgi_read_timeout","uwsgi_send_timeout","uwsgi_store_access","uwsgi_store","uwsgi_string",
					"uwsgi_temp_file_write_size","uwsgi_temp_path",
					
					//optional http modules
					
					//addition: http://wiki.nginx.org/HttpAdditionModule
					"add_before_body", "add_after_body", "addition_types",
					
					//perl: http://wiki.nginx.org/EmbeddedPerlModule
					"perl_modules", "perl_require", "perl_set", "perl",
					
					//flv: http://wiki.nginx.org/HttpFlvStreamModule
					"flv",
					
					//geoIp: http://wiki.nginx.org/HttpGeoIPModule
					"geoip_country", "geoip_city",
					
					//google performance tools: http://wiki.nginx.org/GooglePerftoolsModule
					"google_perftools_profiles",
					
					//gzip static: http://wiki.nginx.org/HttpGzipStaticModule
					"gzip_static", "gzip_http_version", "gzip_proxied", /*"gzip_disable", "gzip_vary",*/ //<-- already accounted for in gzip module
					
					//image filter: http://wiki.nginx.org/HttpImageFilterModule
					"image_filter_buffer", "image_filter_jpeg_quality", "image_filter_transparency", "image_filter",
					
					//random index: http://wiki.nginx.org/HttpRandomIndexModule
					"random_index",
					
					//realip: http://wiki.nginx.org/HttpRealIpModule
					"set_real_ip_from", "real_ip_header",
					
					//secure link: http://wiki.nginx.org/HttpSecureLinkModule
					"secure_link_secret", "secure_link_md5", "secure_link",
					
					//ssl: http://wiki.nginx.org/HttpSslModule
					"ssl_certificate_key","ssl_client_certificate","ssl_certificate","ssl_dhparam","ssl_ciphers","ssl_crl",
					"ssl_prefer_server_ciphers","ssl_protocols","ssl_verify_client","ssl_verify_depth","ssl_session_cache",
					"ssl_session_timeout","ssl_engine", "ssl",
					
					//stub stats: http://wiki.nginx.org/HttpStubStatusModule
					"stub_status",
					
					//sub: http://wiki.nginx.org/HttpSubModule
					"sub_filter_once", "sub_filter_types", "sub_filter",
					
					//dav: http://wiki.nginx.org/HttpDavModule
					"dav_access", "dav_methods", "create_full_put_path",
					
					//xslt: http://wiki.nginx.org/HttpXsltModule
					"xml_entities", "xslt_stylesheet", "xslt_types",
					
					//mail modules
					
					//core: http://wiki.nginx.org/MailCoreModule
					"auth","imap_capabilities","imap_client_buffer","listen","pop3_auth","pop3_capabilities",
					"protocol","server","server_name","smtp_auth","smtp_capabilities","so_keepalive","timeout",
					
					//auth: http://wiki.nginx.org/MailAuthModule
					"auth_http", "auth_http_header","auth_http_timeout",
					
					//proxy: http://wiki.nginx.org/MailProxyModule
					"proxy_buffer", "proxy_pass_error_message", "proxy_timeout", "proxy", "xclient",
					
					//ssl: http://wiki.nginx.org/MailSslModule
					/*"ssl","ssl_certificate","ssl_certificate_key","ssl_ciphers","ssl_prefer_server_ciphers","ssl_protocols","ssl_session_cache","ssl_session_timeout",*/ //<- covered in http ssl module
					"starttls",

					//echo: http://wiki.nginx.org/HttpEchoModule
					"echo_duplicate", "echo_flush", "echo_sleep", "echo_blocking_sleep", "echo_reset_timer",
					"echo_read_request_body", "echo_location_async", "echo_location", "echo_subrequest_async",
					"echo_subrequest", "echo_foreach_split", "echo_end", "echo_request_body", "echo_exec", 
					"echo_before_body", "echo_after_body", "echo",
					
					//others
					"default", "output_buffers"
				], "[\\s;]", false);
				
				return function(context) {
					var token = sunlight.util.matchWord(context, directives, "keyword", true),
						prevToken;
					
					if (!token) {
						return null;
					}
					
					//must be the first word in a statement
					//which means first token in the string, or first non-ws token after "{", "}" or ";"
					prevToken = sunlight.util.getPreviousWhile(context.getAllTokens(), context.count(), function(token) {
						return token.name === "default" || token.name === "comment";
					});
					
					if (!prevToken || (prevToken.name === "punctuation" && (sunlight.util.contains(["{", "}", ";"], prevToken.value)))) {
						context.reader.read(token.value.length - 1); //already read the first character
						return token;
					}
					
					return null;
				};
			}()
		],

		identFirstLetter: /[A-Za-z_-]/,
		identAfterFirstLetter: /[\w-]/,

		operators: [
			"~*", "~", 
			"^~",
			"=",
			"::", ":"
		]
		
	});
}(this["Sunlight"]));
