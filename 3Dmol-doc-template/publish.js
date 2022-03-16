

/**
 * @module template/publish
 * @type {*}
 */
/* global env: true */

const template = require( 'jsdoc/template' );
	const fs = require( 'jsdoc/fs' );
	const _ = require( 'underscore' );
	const path = require( 'jsdoc/path' );

	const {taffy} = require( 'taffydb' );
	const {handle} = require( 'jsdoc/util/error' );
	const helper = require( 'jsdoc/util/templateHelper' );
// jsdoc node support is still a bit odd
	// moment = require( "./moment" ),
	const {htmlsafe} = helper;
	const {linkto} = helper;
	const {resolveAuthorLinks} = helper;
	const {scopeToPunc} = helper;
	const hasOwnProp = Object.prototype.hasOwnProperty;
	const conf = env.conf.templates || {};
	let data;
	let view;
	let outdir = env.opts.destination;

const globalUrl = helper.getUniqueFilename( 'types' );
const indexUrl = helper.getUniqueFilename( 'index' );

const navOptions = {
	systemName            : conf.systemName || "Documentation",
	navType               : conf.navType || "vertical",
	footer                : conf.footer || "",
	copyright             : conf.copyright || "",
	theme                 : conf.theme || "simplex",
	syntaxTheme           : conf.syntaxTheme || "default",
	linenums              : conf.linenums,
	collapseSymbols       : conf.collapseSymbols || false,
	inverseNav            : conf.inverseNav,
	outputSourceFiles     : conf.outputSourceFiles === true,
	sourceRootPath        : conf.sourceRootPath,
	outputSourcePath      : conf.outputSourcePath,
	dateFormat            : conf.dateFormat,
	analytics             : conf.analytics || null,
	highlightTutorialCode : conf.highlightTutorialCode
};

const navigationMaster = {
	index     : {
		title   : navOptions.systemName,
		link    : indexUrl,
		members : []
	},
	namespace : {
		title   : "Namespaces",
		link    : helper.getUniqueFilename( "namespaces.list" ),
		members : []
	},
	module    : {
		title   : "Modules",
		link    : helper.getUniqueFilename( "modules.list" ),
		members : []
	},
	class     : {
		title   : "Classes",
		link    : helper.getUniqueFilename( 'classes.list' ),
		members : []
	},

	mixin    : {
		title   : "Mixins",
		link    : helper.getUniqueFilename( "mixins.list" ),
		members : []
	},
	event    : {
		title   : "Events",
		link    : helper.getUniqueFilename( "events.list" ),
		members : []
	},
	tutorial : {
		title   : "Tutorials",
		link    : helper.getUniqueFilename( "tutorials.list" ),
		members : []
	},
	global   : {
		title   : "Type Specifications",
		link    : globalUrl,
		members : []

	},
	external : {
		title   : "Externals",
		link    : helper.getUniqueFilename( "externals.list" ),
		members : []
	}
};

function find( spec ) {
	return helper.find( data, spec );
}

function tutoriallink( tutorial ) {
	return helper.toTutorial( tutorial, null, { tag : 'em', classname : 'disabled', prefix : 'Tutorial: ' } );
}

function getAncestorLinks( doclet ) {
	return helper.getAncestorLinks( data, doclet );
}

function hashToLink( doclet, hash ) {
	if ( !/^(#.+)/.test( hash ) ) { return hash; }

	let url = helper.createLink( doclet );

	url = url.replace( /(#.+|$)/, hash );
	return `<a href="${  url  }">${  hash  }</a>`;
}

function needsSignature( doclet ) {
	let needsSig = false;

	// function and class definitions always get a signature
	if ( doclet.kind === 'function' || doclet.kind === 'class' ) {
		needsSig = true;
	}
	// typedefs that contain functions get a signature, too
	else if ( doclet.kind === 'typedef' && doclet.type && doclet.type.names &&
		doclet.type.names.length ) {
		for ( let i = 0, l = doclet.type.names.length; i < l; i++ ) {
			if ( doclet.type.names[i].toLowerCase() === 'function' ) {
				needsSig = true;
				break;
			}
		}
	}

	return needsSig;
}

function addSignatureParams( f ) {
	const params = helper.getSignatureParams( f, 'optional' );

	f.signature = `${f.signature || ''  }(${  params.join( ', ' )  })`;
}

function addSignatureReturns( f ) {
	const returnTypes = helper.getSignatureReturns( f );

	f.signature = `<span class="signature">${  f.signature || ''  }</span>` + 
	`<span class="type-signature">${  returnTypes.length ? ` &rarr; {${  returnTypes.join( '|' )  }}` : ''  }</span>`;
}

function addSignatureTypes( f ) {
	const types = helper.getSignatureTypes( f );

	f.signature = `${f.signature || ''  }<span class="type-signature">${  types.length ? ` :${  types.join( '|' )}` : ''  }</span>`;
}

function addAttribs( f ) {
	const attribs = helper.getAttribs( f );

	f.attribs = `<span class="type-signature">${  htmlsafe( attribs.length ? `<${  attribs.join( ', ' )  }> ` : '' )  }</span>`;
}

function shortenPaths( files, commonPrefix ) {
//	// always use forward slashes
//	var regexp = new RegExp( '\\\\', 'g' );
//
//	var prefix = commonPrefix.toLowerCase().replace( regexp, "/" );
//
//	Object.keys( files ).forEach( function ( file ) {
//		files[file].shortened = files[file]
//			.resolved
//			.toLowerCase()
//			.replace( regexp, '/' )
//			.replace( prefix, '' );
//	} );

	Object.keys(files).forEach((file) => {
		files[file].shortened = files[file].resolved.replace(commonPrefix, '')
			// always use forward slashes
			.replace(/\\/g, '/');
	});


	return files;
}

function getPathFromDoclet( doclet ) {
	if ( !doclet.meta ) {
		return;
	}

	return path.normalize(doclet.meta.path && doclet.meta.path !== 'null' ?
		`${doclet.meta.path  }/${  doclet.meta.filename}` :
		doclet.meta.filename);
}

function generate( docType, title, docs, filename, resolveLinks ) {
	resolveLinks = resolveLinks !== false;

	const docData = {
		title,
		docs,
		docType
	};

	const outpath = path.join( outdir, filename );
		let html = view.render( 'container.tmpl', docData );

	if ( resolveLinks ) {
		html = helper.resolveLinks( html ); // turn {@link foo} into <a href="foodoc.html">foo</a>
	}

	fs.writeFileSync( outpath, html, 'utf8' );
}

function generateSourceFiles( sourceFiles ) {
	Object.keys( sourceFiles ).forEach( ( file ) => {
		let source;
		// links are keyed to the shortened path in each doclet's `meta.shortpath` property
		const sourceOutfile = helper.getUniqueFilename( sourceFiles[file].shortened );
		helper.registerLink( sourceFiles[file].shortened, sourceOutfile );

		try {
			source = {
				kind : 'source',
				code : helper.htmlsafe( fs.readFileSync( sourceFiles[file].resolved, 'utf8' ) )
			};
		}
		catch ( e ) {
			handle( e );
		}

		generate( 'source', `Source: ${  sourceFiles[file].shortened}`, [source], sourceOutfile,
			false );
	} );
}

/**
 * Look for classes or functions with the same name as modules (which indicates that the module
 * exports only that class or function), then attach the classes or functions to the `module`
 * property of the appropriate module doclets. The name of each class or function is also updated
 * for display purposes. This function mutates the original arrays.
 *
 * @private
 * @param {Array.<module:jsdoc/doclet.Doclet>} doclets - The array of classes and functions to
 * check.
 * @param {Array.<module:jsdoc/doclet.Doclet>} modules - The array of module doclets to search.
 */
function attachModuleSymbols( doclets, modules ) {
	const symbols = {};

	// build a lookup table
	doclets.forEach( ( symbol ) => {
		symbols[symbol.longname] = symbol;
	} );

	return modules.map( ( module ) => {
		if ( symbols[module.longname] ) {
			module.module = symbols[module.longname];
			module.module.name = `${module.module.name.replace( 'module:', 'require("' )  }")`;
		}
	} );
}

/**
 * Create the navigation sidebar.
 * @param {object} members The members that will be used to create the sidebar.
 * @param {array<object>} members.classes
 * @param {array<object>} members.externals
 * @param {array<object>} members.globals
 * @param {array<object>} members.mixins
 * @param {array<object>} members.modules
 * @param {array<object>} members.namespaces
 * @param {array<object>} members.tutorials
 * @param {array<object>} members.events
 * @return {string} The HTML for the navigation sidebar.
 */
function buildNav( members ) {

	const seen = {};
	const nav = navigationMaster;
	if ( members.modules.length ) {

		members.modules.forEach( ( m ) => {
			if ( !hasOwnProp.call( seen, m.longname ) ) {

				nav.module.members.push( linkto( m.longname, m.longname.replace("module:", "") ) );
			}
			seen[m.longname] = true;
		} );
	}

	if ( members.externals.length ) {

		members.externals.forEach( ( e ) => {
			if ( !hasOwnProp.call( seen, e.longname ) ) {

				nav.external.members.push( linkto( e.longname, e.name.replace( /(^"|"$)/g, '' ) ) );
			}
			seen[e.longname] = true;
		} );
	}

	if ( members.classes.length ) {

		members.classes.forEach( ( c ) => {
			if ( !hasOwnProp.call( seen, c.longname ) ) {

				nav.class.members.push( linkto( c.longname, c.longname.replace("module:", "") ) );
			}
			seen[c.longname] = true;
		} );

	}

	if ( members.events.length ) {

		members.events.forEach( ( e ) => {
			if ( !hasOwnProp.call( seen, e.longname ) ) {

				nav.event.members.push( linkto( e.longname, e.longname.replace("module:", "") ) );
			}
			seen[e.longname] = true;
		} );

	}

	if ( members.namespaces.length ) {

		members.namespaces.forEach( ( n ) => {
			if ( !hasOwnProp.call( seen, n.longname ) ) {

				nav.namespace.members.push( linkto( n.longname, n.longname.replace("module:", "") ) );
			}
			seen[n.longname] = true;
		} );

	}

	if ( members.mixins.length ) {

		members.mixins.forEach( ( m ) => {
			if ( !hasOwnProp.call( seen, m.longname ) ) {

				nav.mixin.members.push( linkto( m.longname, m.longname.replace("module:", "") ) );
			}
			seen[m.longname] = true;
		} );

	}

	if ( members.tutorials.length ) {

		members.tutorials.forEach( ( t ) => {

			nav.tutorial.members.push( tutoriallink( t.name ) );
		} );

	}

	if ( members.globals.length ) {
		members.globals.forEach( ( g ) => {
			if ( !hasOwnProp.call( seen, g.longname ) ) {

				nav.global.members.push( linkto( g.longname, g.longname.replace("module:", "") ) );
			}
			seen[g.longname] = true;
		} );
	}

	const topLevelNav = [];
	_.each( nav, ( entry, name ) => {
		if ( entry.members.length > 0 && name !== "index" ) {
			topLevelNav.push( {
				title   : entry.title,
				link    : entry.link,
				members : entry.members
			} );
		}
	} );
	nav.topLevelNav = topLevelNav;
}

/**
 @param {TAFFY} taffyData See <http://taffydb.com/>.
 @param {object} opts
 @param {Tutorial} tutorials
 */
exports.publish = function ( taffyData, opts, tutorials ) {
	data = taffyData;

	conf.default = conf.default || {};

	const templatePath = opts.template;
	view = new template.Template( `${templatePath  }/tmpl` );

	// claim some special filenames in advance, so the All-Powerful Overseer of Filename Uniqueness
	// doesn't try to hand them out later
//	var indexUrl = helper.getUniqueFilename( 'index' );
	// don't call registerLink() on this one! 'index' is also a valid longname

//	var globalUrl = helper.getUniqueFilename( 'global' );
	helper.registerLink( 'global', globalUrl );

	// set up templating
	view.layout = 'layout.tmpl';

	// set up tutorials for helper
	helper.setTutorials( tutorials );

	data = helper.prune( data );
	data.sort( 'longname, version, since' );
	helper.addEventListeners( data );

	let sourceFiles = {};
	const sourceFilePaths = [];
	data().each( ( doclet ) => {
		doclet.attribs = '';

		if ( doclet.examples ) {
			doclet.examples = doclet.examples.map( ( example ) => {
				let caption; let code; let lang;

				if ( example.match( /^\s*<caption>([\s\S]+?)<\/caption>(\s*[\n\r])([\s\S]+)$/i ) ) {
					caption = RegExp.$1;
					code = RegExp.$3;
				}

				lang = /{@lang (.*?)}/.exec( example );

				if ( lang && lang[1] ) {
					example = example.replace( lang[0], "" );
					lang = lang[1];

				} else {
					lang = null;
				}

				return {
					caption : caption || '',
					code    : code || example,
					lang    : (lang && navOptions.highlightTutorialCode) ? lang : "javascript"
				};
			} );
		}
		if ( doclet.see ) {
			doclet.see.forEach( ( seeItem, i ) => {
				doclet.see[i] = hashToLink( doclet, seeItem );
			} );
		}

		// build a list of source files
		let sourcePath;
		if ( doclet.meta ) {
			sourcePath = getPathFromDoclet( doclet );
			sourceFiles[sourcePath] = {
				resolved  : sourcePath,
				shortened : null
			};

			sourceFilePaths.push( sourcePath );

		}
	} );

	// update outdir if necessary, then create outdir
	const packageInfo = ( find( {kind : 'package'} ) || [] ) [0];
	if ( packageInfo && packageInfo.name ) {
		outdir = path.join( outdir, packageInfo.name, packageInfo.version );
	}
	fs.mkPath( outdir );

	// copy static files to outdir
	const fromDir = path.join( templatePath, 'static' );
		const staticFiles = fs.ls( fromDir, 3 );

	staticFiles.forEach( ( fileName ) => {
		const toDir = fs.toDir( fileName.replace( fromDir, outdir ) );
		fs.mkPath( toDir );
		fs.copyFileSync( fileName, toDir );
	} );

	if ( sourceFilePaths.length ) {
		let payload = navOptions.sourceRootPath;
		if ( !payload ) {
			payload = path.commonPrefix( sourceFilePaths );
		}
		sourceFiles = shortenPaths( sourceFiles, payload );
	}
	data().each( ( doclet ) => {
		const url = helper.createLink( doclet );
		helper.registerLink( doclet.longname, url );

		// add a shortened version of the full path
		let docletPath;
		if ( doclet.meta ) {
			docletPath = getPathFromDoclet( doclet );
			if ( !_.isEmpty( sourceFiles[docletPath] ) ) {
				docletPath = sourceFiles[docletPath].shortened;
				if ( docletPath ) {
					doclet.meta.shortpath = docletPath;
				}
			}
		}
	} );

	data().each( ( doclet ) => {
		const url = helper.longnameToUrl[doclet.longname];

		if ( url.indexOf( '#' ) > -1 ) {
			doclet.id = helper.longnameToUrl[doclet.longname].split( /#/ ).pop();
		}
		else {
			doclet.id = doclet.name;
		}

		if ( needsSignature( doclet ) ) {
			addSignatureParams( doclet );
			addSignatureReturns( doclet );
			addAttribs( doclet );
		}
	} );

	// do this after the urls have all been generated
	data().each( ( doclet ) => {
		doclet.ancestors = getAncestorLinks( doclet );

		if ( doclet.kind === 'member' ) {
			addSignatureTypes( doclet );
			addAttribs( doclet );
		}

		if ( doclet.kind === 'constant' ) {
			addSignatureTypes( doclet );
			addAttribs( doclet );
			doclet.kind = 'member';
		}
	} );

	const members = helper.getMembers( data );
	members.tutorials = tutorials.children;

	// add template helpers
	view.find = find;
	view.linkto = linkto;
	view.resolveAuthorLinks = resolveAuthorLinks;
	view.tutoriallink = tutoriallink;
	view.htmlsafe = htmlsafe;
	// view.moment = moment;

	// once for all
	buildNav( members );
	view.nav = navigationMaster;
	view.navOptions = navOptions;
	attachModuleSymbols( find( { kind : ['class', 'function'], longname : {left : 'module:'} } ),
		members.modules );

	// only output pretty-printed source files if requested; do this before generating any other
	// pages, so the other pages can link to the source files
	if ( navOptions.outputSourceFiles ) {
		generateSourceFiles( sourceFiles );
	}

	if ( members.globals.length ) {
		generate( 'global', 'Global', [
			{kind : 'globalobj'}
		], globalUrl );
	}

	// some browsers can't make the dropdown work
	if ( view.nav.module && view.nav.module.members.length ) {
		generate( 'module', view.nav.module.title, [
			{kind : 'sectionIndex', contents : view.nav.module}
		], navigationMaster.module.link );
	}

	if ( view.nav.class && view.nav.class.members.length ) {
		generate( 'class', view.nav.class.title, [
			{kind : 'sectionIndex', contents : view.nav.class}
		], navigationMaster.class.link );
	}

	if ( view.nav.namespace && view.nav.namespace.members.length ) {
		generate( 'namespace', view.nav.namespace.title, [
			{kind : 'sectionIndex', contents : view.nav.namespace}
		], navigationMaster.namespace.link );
	}

	if ( view.nav.mixin && view.nav.mixin.members.length ) {
		generate( 'mixin', view.nav.mixin.title, [
			{kind : 'sectionIndex', contents : view.nav.mixin}
		], navigationMaster.mixin.link );
	}

	if ( view.nav.external && view.nav.external.members.length ) {
		generate( 'external', view.nav.external.title, [
			{kind : 'sectionIndex', contents : view.nav.external}
		], navigationMaster.external.link );
	}

	if ( view.nav.tutorial && view.nav.tutorial.members.length ) {
		generate( 'tutorial', view.nav.tutorial.title, [
			{kind : 'sectionIndex', contents : view.nav.tutorial}
		], navigationMaster.tutorial.link );
	}

	// index page displays information from package.json and lists files
	const files = find( {kind : 'file'} );
		const packages = find( {kind : 'package'} );

	generate( 'index', 'Index',
		packages.concat(
			[
				{kind : 'mainpage', readme : opts.readme, longname : (opts.mainpagetitle) ? opts.mainpagetitle : 'Main Page'}
			]
		).concat( files ),
		indexUrl );

	// set up the lists that we'll use to generate pages
	const classes = taffy( members.classes );
	const modules = taffy( members.modules );
	const namespaces = taffy( members.namespaces );
	const mixins = taffy( members.mixins );
	const externals = taffy( members.externals );

	for ( const longname in helper.longnameToUrl ) {
		if ( hasOwnProp.call( helper.longnameToUrl, longname ) ) {
			const myClasses = helper.find( classes, {longname} );
			if ( myClasses.length ) {
				generate( 'class', `Class: ${  myClasses[0].name}`, myClasses, helper.longnameToUrl[longname] );
			}

			const myModules = helper.find( modules, {longname} );
			if ( myModules.length ) {
				generate( 'module', `Module: ${  myModules[0].name}`, myModules, helper.longnameToUrl[longname] );
			}

			const myNamespaces = helper.find( namespaces, {longname} );
			if ( myNamespaces.length ) {
				generate( 'namespace', `Namespace: ${  myNamespaces[0].name}`, myNamespaces, helper.longnameToUrl[longname] );
			}

			const myMixins = helper.find( mixins, {longname} );
			if ( myMixins.length ) {
				generate( 'mixin', `Mixin: ${  myMixins[0].name}`, myMixins, helper.longnameToUrl[longname] );
			}

			const myExternals = helper.find( externals, {longname} );
			if ( myExternals.length ) {
				generate( 'external', `External: ${  myExternals[0].name}`, myExternals, helper.longnameToUrl[longname] );
			}
		}
	}

	// TODO: move the tutorial functions to templateHelper.js
	function generateTutorial( title, tutorial, filename ) {
		const tutorialData = {
			title,
			header   : tutorial.title,
			content  : tutorial.parse(),
			children : tutorial.children,
			docs     : null
		};

		const tutorialPath = path.join( outdir, filename );
			let html = view.render( 'tutorial.tmpl', tutorialData );

		// yes, you can use {@link} in tutorials too!
		html = helper.resolveLinks( html ); // turn {@link foo} into <a href="foodoc.html">foo</a>

		fs.writeFileSync( tutorialPath, html, 'utf8' );
	}

	// tutorials can have only one parent so there is no risk for loops
	function saveChildren( node ) {
		node.children.forEach( ( child ) => {
			generateTutorial( `tutorial${  child.title}`, child, helper.tutorialToUrl( child.name ) );
			saveChildren( child );
		} );
	}

	saveChildren( tutorials );
};
