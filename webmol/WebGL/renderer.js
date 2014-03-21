/**
Simplified webGL renderer 
 */

WebMol.Renderer = function ( parameters ) {
    
    parameters = parameters || {};
    
    var _canvas = parameters.canvas !== undefined ? parameters.canvas : document.createElement( 'canvas' ),

    _precision = parameters.precision !== undefined ? parameters.precision : 'highp',

    _alpha = parameters.alpha !== undefined ? parameters.alpha : true,
    _premultipliedAlpha = parameters.premultipliedAlpha !== undefined ? parameters.premultipliedAlpha : true,
    _antialias = parameters.antialias !== undefined ? parameters.antialias : false,
    _stencil = parameters.stencil !== undefined ? parameters.stencil : true,
    _preserveDrawingBuffer = parameters.preserveDrawingBuffer !== undefined ? parameters.preserveDrawingBuffer : false,

    _clearColor = parameters.clearColor !== undefined ? new WebMol.Color( parameters.clearColor ) : new WebMol.Color( 0x000000 ),
    _clearAlpha = parameters.clearAlpha !== undefined ? parameters.clearAlpha : 0;
    
    this.domElement = _canvas;
    this.context = null;
    this.devicePixelRatio = parameters.devicePixelRatio !== undefined
    ? parameters.devicePixelRatio
    : self.devicePixelRatio !== undefined
    ? self.devicePixelRatio
    : 1;

    // clearing

    this.autoClear = true;
    this.autoClearColor = true;
    this.autoClearDepth = true;
    this.autoClearStencil = true;

    // scene graph

    this.sortObjects = true;

    this.autoUpdateObjects = true;
    this.autoUpdateScene = true;
    
    this.renderPluginsPost = [];
    
    // info

    this.info = {

        memory: {
    
        programs: 0,
        geometries: 0,
        textures: 0
    
        },
    
        render: {
    
        calls: 0,
        vertices: 0,
        faces: 0,
        points: 0
    
        }

    };

    // internal properties

    var _this = this,

    _programs = [],
    _programs_counter = 0,

    // internal state cache

    _currentProgram = null,
    _currentFramebuffer = null,
    _currentMaterialId = -1,
    _currentGeometryGroupHash = null,
    _currentCamera = null,
    _geometryGroupCounter = 0,

    _usedTextureUnits = 0,

    // GL state cache

    _oldDoubleSided = -1,
    _oldFlipSided = -1,

    _oldBlending = -1,

    _oldBlendEquation = -1,
    _oldBlendSrc = -1,
    _oldBlendDst = -1,

    _oldDepthTest = -1,
    _oldDepthWrite = -1,

    _oldPolygonOffset = null,
    _oldPolygonOffsetFactor = null,
    _oldPolygonOffsetUnits = null,

    _oldLineWidth = null,

    _viewportX = 0,
    _viewportY = 0,
    _viewportWidth = 0,
    _viewportHeight = 0,
    _currentWidth = 0,
    _currentHeight = 0,

    _enabledAttributes = {},

     // camera matrices cache

    _projScreenMatrix = new WebMol.Matrix4(),

    _vector3 = new WebMol.Vector3(),

    // light arrays cache

    _direction = new WebMol.Vector3(),

    _lightsNeedUpdate = true,

    _lights = {

            ambient: [0,0,0],
            directional: { length: 0, colors: new Array(), positions: new Array() },
            point: { length: 0, colors: new Array(), positions: new Array(), distances: new Array() },
            spot: { length: 0, colors: new Array(), positions: new Array(), distances: new Array(), directions: new Array(), anglesCos: new Array(), exponents: new Array() },
            hemi: { length: 0, skyColors: new Array(), groundColors: new Array(), positions: new Array() }

    };

    // initialize

    var _gl;

    initGL();

    setDefaultGLState();

    this.context = _gl;    

    // API

    this.getContext = function () {

            return _gl;

    };

    this.getPrecision = function () {

            return _precision;

    };
    
    this.setClearColorHex = function ( hex, alpha ) {

            _clearColor.setHex( hex );
            _clearAlpha = alpha;

            _gl.clearColor( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha );

    };

    this.setSize = function ( width, height ) {

            _canvas.width = width * this.devicePixelRatio;
            _canvas.height = height * this.devicePixelRatio;

            _canvas.style.width = width + 'px';
            _canvas.style.height = height + 'px';

            this.setViewport( 0, 0, _canvas.width, _canvas.height );

    };

    this.setViewport = function ( x, y, width, height ) {

            _viewportX = x !== undefined ? x : 0;
            _viewportY = y !== undefined ? y : 0;

            _viewportWidth = width !== undefined ? width : _canvas.width;
            _viewportHeight = height !== undefined ? height : _canvas.height;

            _gl.viewport( _viewportX, _viewportY, _viewportWidth, _viewportHeight );

    };

    this.clear = function ( color, depth, stencil ) {

            var bits = 0;

            if ( color === undefined || color ) bits |= _gl.COLOR_BUFFER_BIT;
            if ( depth === undefined || depth ) bits |= _gl.DEPTH_BUFFER_BIT;
            if ( stencil === undefined || stencil ) bits |= _gl.STENCIL_BUFFER_BIT;

            _gl.clear( bits );

    };

    this.clearTarget = function ( renderTarget, color, depth, stencil ) {

            this.setRenderTarget( renderTarget );
            this.clear( color, depth, stencil );

    };

    this.setMaterialFaces = function ( material ) {

            var doubleSided = material.side === WebMol.DoubleSide;
            var flipSided = material.side === WebMol.BackSide;

            if ( _oldDoubleSided !== doubleSided ) {

                    if ( doubleSided ) {

                            _gl.disable( _gl.CULL_FACE );

                    } else {

                            _gl.enable( _gl.CULL_FACE );

                    }

                    _oldDoubleSided = doubleSided;

            }

            if ( _oldFlipSided !== flipSided ) {

                    if ( flipSided ) {

                            _gl.frontFace( _gl.CW );

                    } else {

                            _gl.frontFace( _gl.CCW );

                    }

                    _oldFlipSided = flipSided;

            }    

    };
    
    this.setDepthTest = function ( depthTest ) {

            if ( _oldDepthTest !== depthTest ) {

                    if ( depthTest ) {

                            _gl.enable( _gl.DEPTH_TEST );

                    } else {

                            _gl.disable( _gl.DEPTH_TEST );

                    }

                    _oldDepthTest = depthTest;

            }

    };

    this.setDepthWrite = function ( depthWrite ) {

            if ( _oldDepthWrite !== depthWrite ) {

                    _gl.depthMask( depthWrite );
                    _oldDepthWrite = depthWrite;

            }

    };

    this.setBlending = function( blending ) {

            if (blending === WebMol.NoBlending) 
                    _gl.disable( _gl.BLEND );

            else {
                    _gl.enable( _gl.BLEND );
                    _gl.blendEquationSeparate( _gl.FUNC_ADD, _gl.FUNC_ADD );
                    _gl.blendFuncSeparate( _gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA, _gl.ONE, _gl.ONE_MINUS_SRC_ALPHA );
            }

            _oldBlending = blending;
    };
    
    // Plugins
    
    this.addPostPlugin = function(plugin) {

        plugin.init(this);
        this.renderPluginsPost.push(plugin);

    };

    // Sorting

    function numericalSort ( a, b ) {

            return b[ 0 ] - a[ 0 ];

    };

    function enableAttribute( attribute ) {

            if ( ! _enabledAttributes[ attribute ] ) {

                    _gl.enableVertexAttribArray( attribute );
                    _enabledAttributes[ attribute ] = true;

            }

    };

    function disableAttributes() {

            for ( var attribute in _enabledAttributes ) {

                    if ( _enabledAttributes[ attribute ] ) {

                            _gl.disableVertexAttribArray( attribute );
                            _enabledAttributes[ attribute ] = false;

                    }

            }

    };    

    function setPolygonOffset ( polygonOffset, factor, units) {

            if ( _oldPolygonOffset !== polygonOffset ) {

                    if (polygonOffset)
                            _gl.enable( _gl.POLYGON_OFFSET_FILL );
                    else
                            _gl.disable( _gl.POLYGON_OFFSET_FILL );
            }
    };

    function setLineWidth ( width ) {

        if ( width !== _oldLineWidth ) {
                _gl.lineWidth(width);
                _oldLineWidth = width;
        }

    };
    
    var onGeometryDispose = function(event) {
        
        var geometry = event.target;
        geometry.removeEventListener('dispose', onGeometryDispose);
        
        deallocateGeometry(geometry);
        
        _this.info.memory.geometries--;
        
    };
    
    var onTextureDispose = function(event) {

        var texture = event.target;

        texture.removeEventListener('dispose', onTextureDispose);

        deallocateTexture(texture);

        _this.info.memory.textures--;


    };
    
    var onMaterialDispose = function(event) {
        
        var material = event.target;
        material.removeEventListener('dispose', onMaterialDispose);
        
        deallocateMaterial(material);
        
    };
    
    var deallocateGeometry = function(geometry) {
        
        geometry.__webglInit = undefined;
        
        if (geometry.__webglVertexBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglVertexBuffer);
        
        if (geometry.__webglColorBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglColorBuffer);
        
        if (geometry.geometryChunks !== undefined) {
            
            for (var g in geometry.geometryChunks) {  
                
                var geometryChunk = geometry.geometryChunks[g];

                if (geometryChunk.__webglVertexBuffer !== undefined)
                    _gl.deleteBuffer(geometryChunk.__webglVertexBuffer);

                if (geometryChunk.__webglColorBuffer !== undefined)
                    _gl.deleteBuffer(geometryChunk.__webglColorBuffer);
                
                if (geometryChunk.__webglNormalBuffer !== undefined)
                    _gl.deleteBuffer(geometryChunk.__webglNormalBuffer);  
                
                if (geometryChunk.__webglFaceBuffer !== undefined)
                    _gl.deleteBuffer(geometryChunk.__webglFaceBuffer);
            }
        }
    };
    
    var deallocateMaterial = function (material) {

        var program = material.program;

        if ( program === undefined ) return;

        material.program = undefined;

        // only deallocate GL program if this was the last use of shared program
        // assumed there is only single copy of any program in the _programs list
        // (that's how it's constructed)

        var i, il, programInfo;
        var deleteProgram = false;

        for ( i = 0, il = _programs.length; i < il; i ++ ) {

                programInfo = _programs[ i ];

                if ( programInfo.program === program ) {

                        programInfo.usedTimes --;

                        if ( programInfo.usedTimes === 0 ) {

                                deleteProgram = true;

                        }

                        break;

                }

        }

        if ( deleteProgram === true ) {

                // avoid using array.splice, this is costlier than creating new array from scratch

                var newPrograms = [];

                for ( i = 0, il = _programs.length; i < il; i ++ ) {

                        programInfo = _programs[ i ];

                        if ( programInfo.program !== program ) {

                                newPrograms.push( programInfo );

                        }

                }

                _programs = newPrograms;

                _gl.deleteProgram( program );

                _this.info.memory.programs --;

        }

    };
    
    var deallocateTexture = function(texture) {

        if (texture.image && texture.image.__webglTextureCube) {

            // cube texture

            _gl.deleteTexture(texture.image.__webglTextureCube);

        } 
        
        else {

            // 2D texture

            if ( ! texture.__webglInit ) return;

            texture.__webglInit = false;
            _gl.deleteTexture( texture.__webglTexture );

        }

    };

    //Compile and return shader
    function getShader (type, str) {

        var shader;

        if (type === "fragment")
                shader = _gl.createShader( _gl.FRAGMENT_SHADER );
        else if (type === "vertex")
                shader = _gl.createShader( _gl.VERTEX_SHADER );

        _gl.shaderSource(shader, str);
        _gl.compileShader(shader);

        if ( ! _gl.getShaderParameter(shader, _gl.COMPILE_STATUS) ) {

                console.error(_gl.getShaderInfoLog(shader));
                console.error("could not initialize shader");
                return null;

        }

        return shader;

    };    


    //Compile appropriate shaders (if necessary) from source code and attach to gl program.
    function buildProgram(fragmentShader, vertexShader, uniforms) {

        var p, pl, d, program, code;
        var chunks = [];

        chunks.push(fragmentShader);
        chunks.push(vertexShader);

        code = chunks.join();

        //check if program has already been compiled

        for (p = 0, pl = _programs.length; p < pl; p++) {

                var programInfo = _programs[p];

                if (programInfo.code === code) {

                        programInfo.usedTimes++;

                        return programInfo.program;
                }
        }

        //Set up new program and compile shaders

        program = _gl.createProgram();
        
        //set up precision
        var precision = _precision;
        var prefix = "precision " + precision + " float;\n";
        
        var glFragmentShader = getShader("fragment", prefix + fragmentShader);
        var glVertexShader = getShader("vertex", prefix + vertexShader);

        _gl.attachShader(program, glVertexShader);
        _gl.attachShader(program, glFragmentShader);

        _gl.linkProgram(program);

        if (! _gl.getProgramParameter(program, _gl.LINK_STATUS) )
                console.error("Could not initialize shader");

        //gather and cache uniform variables and attributes

        program.uniforms = {};
        program.attributes = {};

        var identifiers, u, a, i;

        //uniform vars
        identifiers = 
                [ 'viewMatrix', 'modelViewMatrix', 'projectionMatrix', 'normalMatrix', 'modelMatrix', 'cameraPosition' ];

        //custom uniform vars
        for (u in uniforms) 
                identifiers.push(u);

        for (i = 0; i < identifiers.length; i++) {

                var uniformVar = identifiers[i];
                program.uniforms[uniformVar] = _gl.getUniformLocation(program, uniformVar);

        }

        //attributes
        identifiers = 
                [ 'position', 'normal', 'color', 'lineDistance' ];

        /*
        for (a in attributes)
                identifiers.push(a);
        */

        for (i = 0; i < identifiers.length; i++) {

                var attributeVar = identifiers[i];
                program.attributes[attributeVar] = _gl.getAttribLocation(program, attributeVar);
        }

        program.id = _programs_counter++;
        _programs.push( {program: program, code: code, usedTimes: 1} );
        _this.info.memory.programs = _programs.length;

        return program;
    };

    //TODO: need to set up shader attributes and uniforms as attributes on material object after attaching prgm
    //We need to attach appropriate uniform variables to material after shaders have been chosen
    this.initMaterial = function ( material, lights, fog, object ) {

        material.addEventListener('dispose', onMaterialDispose);

        var u, a, identifiers, i, parameters, maxLightCount, maxBones, maxShadows, shaderID;

        if (material instanceof WebMol.LineBasicMaterial)
                shaderID = "basic";
        else if (material instanceof WebMol.MeshLambertMaterial)
                shaderID = "lambert";

        if (shaderID) {

                var shader = WebMol.ShaderLib[shaderID];
                material.shaderType = shaderID;
                material.vertexShader = shader.vertexShader;
                material.fragmentShader = shader.fragmentShader;
                material.uniforms = WebMol.ShaderUtils.clone(shader.uniforms);
                //TODO: set material uniforms to shader uniform variables

        }

        material.program = buildProgram(material.fragmentShader, material.vertexShader, material.uniforms);

    };

    function setProgram( camera, lights, fog, material, object ) {

        if ( material.needsUpdate ) {

                if (material.program)
                        deallocateMaterial(material);

                        _this.initMaterial( material, lights, fog, object );
                        material.needsUpdate = false;
        }

        var refreshMaterial = false;

        //p_uniforms: uniformVarName => uniformLocation
        //m_uniforms: uniformVarName => uniformJsVal
        var program = material.program,
                p_uniforms = program.uniforms,
                m_uniforms = material.uniforms;

        if (program != _currentProgram) {        
                _gl.useProgram(program);
                _currentProgram = program;

                refreshMaterial = true;
        }

        if (material.id != _currentMaterialId) {
                _currentMaterialId = material.id;
                refreshMaterial = true;
        }

        if (camera != _currentCamera) {    
                _currentCamera = camera;
                refreshMaterial = true;
        }

        //Send projection matrix to uniform variable in shader
        if (refreshMaterial) {

                //Load projection, model-view matrices for perspective
                _gl.uniformMatrix4fv(p_uniforms.projectionMatrix, false, camera.projectionMatrix.elements);
                _gl.uniformMatrix4fv(p_uniforms.modelViewMatrix, false, object._modelViewMatrix.elements);

                //Set up correct fog uniform vals
                m_uniforms.fogColor.value = fog.color;
                m_uniforms.fogNear.value = fog.near;
                m_uniforms.fogFar.value = fog.far;

                //Set up lights for lambert shader
                if (material.shaderType === "lambert") {

                        //load view and normal matrices for directional and object lighting
                        _gl.uniformMatrix4fv(p_uniforms.viewMatrix, false, camera.matrixWorldInverse.elements);
                        _gl.uniformMatrix3fv(p_uniforms.normalMatrix, false, object._normalMatrix.elements);
                        //_gl.uniformMatrix4fv(p_uniforms.modelMatrix, false, object.matrixWorld.elements);

                        if (_lightsNeedUpdate) {
                                setupLights(program, lights);
                                _lightsNeedUpdate = false;
                        }

                        //Set up correct light uniform var vals
                        m_uniforms.ambientLightColor.value = _lights.ambient;
                        m_uniforms.directionalLightColor.value = _lights.directional.colors;
                        m_uniforms.directionalLightDirection.value = _lights.directional.positions;
                        m_uniforms.ambient.value = material.ambient;
                        m_uniforms.emissive.value = material.emissive;

                }

                //opacity, diffuse, emissive, etc
                m_uniforms.opacity.value = material.opacity;
                m_uniforms.diffuse.value = material.color;

                //Load any other material specific uniform variables to gl shaders
                loadMaterialUniforms(p_uniforms, m_uniforms);

        }

        return program;

    };

    function loadMaterialUniforms(p_uniforms, m_uniforms) {
            var uniformVar, type, uniformVal, uniformLoc;

            for (uniformVar in m_uniforms) {
                    if (! p_uniforms[uniformVar])
                            continue;

                    type = m_uniforms[uniformVar].type;
                    uniformVal = m_uniforms[uniformVar].value;
                    uniformLoc = p_uniforms[uniformVar];

                    //single float
                    if (type === 'f')
                            _gl.uniform1f(uniformLoc, uniformVal);
                    //array of floats
                    else if (type === 'fv')
                            _gl.uniform3fv(uniformLoc, uniformVal);
                    //color - r,g,b floats
                    else if (type === 'c')
                            _gl.uniform3f(uniformLoc, uniformVal.r, uniformVal.g, uniformVal.b);

            }

    };

    this.renderBuffer = function ( camera, lights, fog, material, geometryGroup, object ) {

            if ( ! material.visible )
                    return;

            var program, attributes, linewidth, primitives, a, attribute, i, il;

            //Sets up proper vertex and fragment shaders and attaches them to webGL program
            //Also sets appropriate uniform variables 
            program = setProgram(camera, lights, fog, material, object);

            attributes = program.attributes;

            var updateBuffers = false,
                    geometryGroupHash = (geometryGroup.id * 0xffffff) + (program.id * 2);

            if (geometryGroupHash !== _currentGeometryGroupHash) {
                    _currentGeometryGroupHash = geometryGroupHash;
                    updateBuffers = true;
            }

            //rebind shader attributes to appropriate (and already initialized) gl buffers
            if (updateBuffers) {

                    disableAttributes();

                    // Vertices
                    if (attributes.position >= 0) {            
                            _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglVertexBuffer );
                            enableAttribute( attributes.position );
                            _gl.vertexAttribPointer( attributes.position, 3, _gl.FLOAT, false, 0, 0 );    
                    }

                    // Colors
                    if (attributes.color >= 0) {
                            _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglColorBuffer);
                            enableAttribute( attributes.color );
                            _gl.vertexAttribPointer( attributes.color, 3, _gl.FLOAT, false, 0, 0 );
                    }

                    // Normals (lambert shader only)
                    if (attributes.normal >=0) {
                            _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryGroup.__webglNormalBuffer );
                            enableAttribute( attributes.normal );
                            _gl.vertexAttribPointer( attributes.normal, 3, _gl.FLOAT, false, 0, 0 );
                    }

            }

            //Render

            //lambert shaders - draw triangles
            //TODO: make sure geometryGroup's face count is setup correctly
            if (object instanceof WebMol.Mesh) {

                    var faceCount = geometryGroup.__faceArray.length;

                    if (updateBuffers)
                            _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryGroup.__webglFaceBuffer );
                    _gl.drawElements( _gl.TRIANGLES, faceCount, _gl.UNSIGNED_SHORT, 0 );

                    _this.info.render.calls++;
                    _this.info.render.vertices += faceCount;
                    _this.info.render.faces += faceCount / 3;
            }

            //basic shaders - draw lines
            else if (object instanceof WebMol.Line) {

                    var lineCount = geometryGroup.vertices;

                    setLineWidth(material.lineWidth);
                    _gl.drawArrays( _gl.LINES, 0, lineCount );

                    _this.info.render.calls++;
            }

    };

    //rendering
    function renderObjects ( renderList, reverse, materialType, camera, lights, fog, useBlending, overrideMaterial)  {

            var webglObject, object, buffer, material, start, end, delta;

            //Forward or backward render

            if (reverse) {
                    start = renderList.length - 1;
                    end = -1;
                    delta = -1;
            }

            else {
                    start = 0;
                    end = renderList.length;
                    delta = 1;
            }

            for (var i = start; i !== end; i += delta) {

                    webglObject = renderList[i];

                    if (webglObject.render) {

                            object = webglObject.object;
                            buffer = webglObject.buffer;
                            material = webglObject[materialType];

                            if ( ! material )
                                    continue;

                            if (useBlending)
                                    _this.setBlending(material.blending);

                            _this.setDepthTest(material.depthTest);
                            _this.setDepthWrite(material.depthWrite);
                            setPolygonOffset(material.polygonOffset, material.polygonOffsetFactor, material.polygonOffsetUnits);

                            _this.setMaterialFaces(material);

                            _this.renderBuffer(camera, lights, fog, material, buffer, object);
                    }
            }

    };
    
    this.render = function ( scene, camera, renderTarget, forceClear ) {

            if ( camera instanceof WebMol.Camera === false )  {

                    console.error( 'WebMol.Renderer.render: camera is not an instance of WebMol.Camera.' );
                    return;

            }

            var i, il,

            webglObject, object,
            renderList,

            lights = scene.__lights,
            fog = scene.fog;

            // reset caching for this frame

            _currentMaterialId = -1;
            _lightsNeedUpdate = true;

            // update scene graph

            if ( this.autoUpdateScene ) scene.updateMatrixWorld();

            // update camera matrices
            //Pretty sure camera's parent is always going to be undefined for our purposes...
            if ( camera.parent === undefined ) camera.updateMatrixWorld();

            camera.matrixWorldInverse.getInverse( camera.matrixWorld );

            _projScreenMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

            // update WebGL objects

            if ( this.autoUpdateObjects ) this.initWebGLObjects( scene );


            _this.info.render.calls = 0;
            _this.info.render.vertices = 0;
            _this.info.render.faces = 0;
            _this.info.render.points = 0;

            _currentWidth = _viewportWidth;
            _currentHeight = _viewportHeight;

            if ( this.autoClear || forceClear ) {

                    this.clear( this.autoClearColor, this.autoClearDepth, this.autoClearStencil );

            }

            // set matrices for regular objects (frustum culled)

            renderList = scene.__webglObjects;

            for ( i = 0, il = renderList.length; i < il; i ++ ) {

                    webglObject = renderList[ i ];
                    object = webglObject.object;

                    webglObject.render = false;

                    if ( object.visible ) {        
                            setupMatrices( object, camera );
                            unrollBufferMaterial( webglObject );
                            webglObject.render = true;
                    }
            }

            // set matrices for immediate objects

            var material = null;

            // opaque pass (front-to-back order)

            this.setBlending( WebMol.NoBlending );

            renderObjects( scene.__webglObjects, true, "opaque", camera, lights, fog, false, material );

            // transparent pass (back-to-front order)

            renderObjects( scene.__webglObjects, false, "transparent", camera, lights, fog, true, material );

            // Render plugins (e.g. sprites), and reset state
            
            renderPlugins(this.renderPluginsPost, scene, camera);

            // Ensure depth buffer writing is enabled so it can be cleared on next render

            this.setDepthTest( true );
            this.setDepthWrite( true );

            //_gl.finish();

    };
    
    function renderPlugins(plugins, scene, camera) {
        
        //Reset state once regardless
        //This should also fix cartoon render bug (after transparent surface render)
        
        _currentGeometryGroupHash = -1;
        _currentProgram = null;
        _currentCamera = null;
        _oldBlending = -1;
        _oldDepthWrite = -1;
        _oldDepthTest = -1;
        _oldDoubleSided = -1;
        _currentMaterialId = -1;
        _oldFlipSided = -1;
        
        
        if (!plugins.length)
            return;
        
        for (var i = 0, il = plugins.length; i < il; i++) {
            
            _lightsNeedUpdate = true;
            
            plugins[i].render(scene, camera, _currentWidth, _currentHeight);
            
            //Reset state after plugin render
            _currentGeometryGroupHash = -1;
            _currentProgram = null;
            _currentCamera = null;
            _oldBlending = -1;
            _oldDepthWrite = -1;
            _oldDepthTest = -1;
            _oldDoubleSided = -1;
            _currentMaterialId = -1;
            _oldFlipSided = -1;       
                
        }  
        
    };

    this.initWebGLObjects = function ( scene ) {

            if ( !scene.__webglObjects ) {

                    scene.__webglObjects = [];
                    scene.__webglObjectsImmediate = [];
                    scene.__webglSprites = [];
                    scene.__webglFlares = [];

            }

            //Add objects; this sets up buffers for each geometryChunk
            if (scene.__objectsAdded.length) {
                
                while(scene.__objectsAdded.length){
                    addObject(scene.__objectsAdded[0], scene);
                    scene.__objectsAdded.splice(0, 1);
                }
                
                //Force buffer update during render
                //Hackish fix for initial cartoon-render-then-transparent-surface bug
                _currentGeometryGroupHash = -1;
                
            }

            while (scene.__objectsRemoved.length) {

                    removeObject(scene.__objectsRemoved[ 0 ], scene);
                    scene.__objectsRemoved.splice(0, 1);

            }

            // update must be called after objects adding / removal
            //This sends typed arrays to GL buffers for each geometryChunk
            for ( var o = 0, ol = scene.__webglObjects.length; o < ol; o ++ ) {

                    updateObject(scene.__webglObjects[ o ].object);

            }

    };
    
    // Objects adding

    function addObject (object, scene) {

        var g, geometry, material, geometryGroup;

        if ( !object.__webglInit ) {

            object.__webglInit = true;

            object._modelViewMatrix = new WebMol.Matrix4();
            object._normalMatrix = new WebMol.Matrix3();

            if (object.geometry !== undefined && object.geometry.__webglInit === undefined) {

                    object.geometry.__webglInit = true;
                    object.geometry.addEventListener('dispose', onGeometryDispose);

            }

            if (object instanceof WebMol.Mesh){
                geometry = object.geometry;
                material = object.material;

                for (g in geometry.geometryChunks) {

                        var geometryChunk = geometry.geometryChunks[ g ];
                        
                        geometryChunk.id = _geometryGroupCounter++;

                        // initialise VBO on the first access

                        if ( !geometryChunk.__webglVertexBuffer ) {

                                createMeshBuffers(geometryChunk);

                                geometry.verticesNeedUpdate = true;
                                geometry.elementsNeedUpdate = true;
                                geometry.normalsNeedUpdate = true;
                                geometry.colorsNeedUpdate = true;

                        }
                }
            } 
            
            else if (object instanceof WebMol.Line) {
                
                geometry = object.geometry;
                if ( !geometry.__webglVertexBuffer ) {

                    createLineBuffers( geometry );
                    geometry.verticesNeedUpdate = true;
                    geometry.colorsNeedUpdate = true;

                }
            }
        
        }
        
        if ( ! object.__webglActive ) {
            
            if (object instanceof WebMol.Mesh) {
                
                geometry = object.geometry;

                for ( g in geometry.geometryChunks ) {
                        geometryChunk = geometry.geometryChunks[g];

                        addBuffer(scene.__webglObjects, geometryChunk, object);
                }
                
            }
            
            else if (object instanceof WebMol.Line) {
                geometry = object.geometry;
                addBuffer(scene.__webglObjects, geometry, object);
            }
            
            
            //Sprite
            else if (object instanceof WebMol.Sprite) 
                scene.__webglSprites.push(object);
         
            
            object.__webglActive = true;
            
        }

    };

    function updateObject ( object ) {

        var geometry = object.geometry, material = object.material,
                geometryGroup, geometryChunk, customAttributesDirty;

        if ( object instanceof WebMol.Mesh ) {
            
            for (g in geometry.geometryChunks ) {
                
                geometryChunk = geometry.geometryChunks[ g ];

                if ( geometry.verticesNeedUpdate || geometry.elementsNeedUpdate ||
                         geometry.colorsNeedUpdate || geometry.normalsNeedUpdate) {
                        setBuffers( geometryChunk, _gl.DYNAMIC_DRAW );
                }
            }
            
            geometry.verticesNeedUpdate = false;
            geometry.elementsNeedUpdate = false;
            geometry.normalsNeedUpdate = false;
            geometry.colorsNeedUpdate = false;

            geometry.buffersNeedUpdate = false;

        }
        
        else if ( object instanceof WebMol.Line ) {
            
            if ( geometry.verticesNeedUpdate || geometry.colorsNeedUpdate) 
                setBuffers( geometry, _gl.DYNAMIC_DRAW, true);
                
                geometry.veticesNeedUpdate = false;
                geometry.colorsNeedUpdate = false;
                
        }

    };
    
    function removeObject( object, scene ) {

        if (object instanceof WebMol.Mesh || object instanceof WebMol.Line )
            removeInstances(scene.__webglObjects, object);

        else if (object instanceof WebMol.Sprite)
            removeInstancesDirect(scene.__webglSprites, object);
            
        object.__webglActive = false;

    };

    function removeInstances( objList, object ) {

            for (var o = objList.length - 1; o >= 0; --o) {

                    if (objList[o].object === object) 
                            objList.splice(o, 1);

            }
    };

    function removeInstancesDirect( objList, object ) {

            for (var o = objList.length - 1; o >= 0; --o) {

                    if (objList[o] === object) 
                            objList.splice(o, 1);

            }
    };

    function unrollBufferMaterial( globject ) {

            var object = globject.object;
            var material = object.material;

            if ( material.transparent) {                    
                    globject.opaque = null;
                    globject.transparent = material;
            }

            else {
                    globject.opaque = material;
                    globject.transparent = null;
            }

    };

    function setBuffers( geometryChunk, hint, line ) {

        var vertexArray = geometryChunk.__vertexArray;
        var colorArray = geometryChunk.__colorArray;
         
        //vertex buffers
        _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryChunk.__webglVertexBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, vertexArray, hint );        

        //color buffers
        _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryChunk.__webglColorBuffer );
        _gl.bufferData( _gl.ARRAY_BUFFER, colorArray, hint );    
              
        //set buffers for line render
        
            //set buffers for mesh render
            if (line === undefined) {

                var normalArray = geometryChunk.__normalArray;
                var faceArray = geometryChunk.__faceArray;

                //normal buffers
                _gl.bindBuffer( _gl.ARRAY_BUFFER, geometryChunk.__webglNormalBuffer );
                _gl.bufferData( _gl.ARRAY_BUFFER, normalArray, hint );        

                //face (index) buffers
                _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, geometryChunk.__webglFaceBuffer );
                _gl.bufferData( _gl.ELEMENT_ARRAY_BUFFER, faceArray, hint );    

            }

    };
    
    //Creates appropriate gl buffers for geometry chunk
    //TODO: do we need line buffer for mesh objects?
    //Also, can we integrate this with createLineBuffers?
    function createMeshBuffers ( geometryChunk ) {

        geometryChunk.__webglVertexBuffer = _gl.createBuffer();
        geometryChunk.__webglNormalBuffer = _gl.createBuffer();
        geometryChunk.__webglColorBuffer = _gl.createBuffer();

        geometryChunk.__webglFaceBuffer = _gl.createBuffer();

        _this.info.memory.geometries++;
    };
    
    function createLineBuffers ( geometry ) {
        
        geometry.__webglVertexBuffer = _gl.createBuffer();
        geometry.__webglColorBuffer = _gl.createBuffer();
        
        _this.info.memory.geometries++;
    };

    function addBuffer (objlist, buffer, object) {

        objlist.push(
                {
                        buffer: buffer,
                        object: object,
                        opaque: null,
                        transparent: null
                }
        );

    };

    function setupMatrices (object, camera) {

        object._modelViewMatrix.multiplyMatrices( camera.matrixWorldInverse, object.matrixWorld );

        object._normalMatrix.getInverse( object._modelViewMatrix );
        object._normalMatrix.transpose();

    };

    function isPowerOfTwo ( value ) {

        return ( value & ( value - 1 ) ) === 0;

    };
    
    // Fallback filters for non-power-of-2 textures

    function filterFallback ( f ) {

        return _gl.LINEAR;

    };

    function setTextureParameters ( textureType, texture, isImagePowerOfTwo ) {

        if ( isImagePowerOfTwo ) {

            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, paramToGL( texture.wrapS ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, paramToGL( texture.wrapT ) );

            _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, paramToGL( texture.magFilter ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, paramToGL( texture.minFilter ) );

        } else {

            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE );
            _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE );

            _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, filterFallback( texture.magFilter ) );
            _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, filterFallback( texture.minFilter ) );

        }

    };
    
    this.setTexture = function (texture, slot) {

        if (texture.needsUpdate) {

            if ( !texture.__webglInit ) {

                texture.__webglInit = true;

                texture.addEventListener('dispose', onTextureDispose);

                texture.__webglTexture = _gl.createTexture();

                _this.info.memory.textures++;

            }

            _gl.activeTexture(_gl.TEXTURE0 + slot);
            _gl.bindTexture(_gl.TEXTURE_2D, texture.__webglTexture);

            _gl.pixelStorei(_gl.UNPACK_FLIP_Y_WEBGL, texture.flipY);
            _gl.pixelStorei(_gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, texture.premultiplyAlpha);
            _gl.pixelStorei(_gl.UNPACK_ALIGNMENT, texture.unpackAlignment);

            var image = texture.image,
            isImagePowerOfTwo = isPowerOfTwo(image.width) && isPowerOfTwo(image.height),
            glFormat = paramToGL(texture.format),
            glType = paramToGL(texture.type);

            setTextureParameters(_gl.TEXTURE_2D, texture, isImagePowerOfTwo);

            var mipmap, mipmaps = texture.mipmaps;

            // regular Texture (image, video, canvas)

            // use manually created mipmaps if available
            // if there are no manual mipmaps
            // set 0 level mipmap and then use GL to generate other mipmap levels

            if ( mipmaps.length > 0 && isImagePowerOfTwo ) {

                for ( var i = 0, il = mipmaps.length; i < il; i ++ ) {
                    mipmap = mipmaps[ i ];
                    _gl.texImage2D( _gl.TEXTURE_2D, i, glFormat, glFormat, glType, mipmap );
                }
                
                texture.generateMipmaps = false;
            } 
            
            else 
                _gl.texImage2D( _gl.TEXTURE_2D, 0, glFormat, glFormat, glType, texture.image );

            
            if ( texture.generateMipmaps && isImagePowerOfTwo ) _gl.generateMipmap( _gl.TEXTURE_2D );

            texture.needsUpdate = false;

            if ( texture.onUpdate ) texture.onUpdate();

        } else {

            _gl.activeTexture( _gl.TEXTURE0 + slot );
            _gl.bindTexture( _gl.TEXTURE_2D, texture.__webglTexture );

        }

    };
    
    // Map constants to WebGL constants

    function paramToGL ( p ) {

        if ( p === WebMol.UnsignedByteType ) return _gl.UNSIGNED_BYTE;
        if ( p === WebMol.RGBAFormat ) return _gl.RGBA;

        return 0;

    };
    
    function setupLights ( program, lights ) {
        var l, ll, light, n,
        r = 0, g = 0, b = 0,
        color,
        position,
        intensity,
        distance,
        
        zlights = _lights,
        
        dirColors = zlights.directional.colors,
        dirPositions = zlights.directional.positions,
        
        dirCount = 0,
        dirLength = 0,
        dirOffset = 0;
        
        for ( l = 0, ll = lights.length; l < ll; l++) {
            
            light = lights[l];
            
            color = light.color;
            intensity = light.intensity;
            distance = light.distance;
            
            if (light instanceof WebMol.Light) {
                
                dirCount++;
                
                _direction.getPositionFromMatrix(light.matrixWorld);
                _vector3.getPositionFromMatrix(light.target.matrixWorld);
                _direction.sub(_vector3);
                _direction.normalize();
                
                if (_direction.x === 0 && _direction.y === 0 && _direction.z === 0)
                    continue;
                
                dirPositions[dirOffset] = _direction.x;
                dirPositions[dirOffset + 1] = _direction.y;
                dirPositions[dirOffset + 2] = _direction.z;

                dirColors[dirOffset] = color.r * intensity;
                dirColors[dirOffset + 1] = color.g * intensity;
                dirColors[dirOffset + 2] = color.b * intensity;
                
                dirOffset += 3;
                
                dirLength++;    
            }
        
        }

        zlights.ambient[0] = r;
        zlights.ambient[1] = g;
        zlights.ambient[2] = b;
        zlights.directional.length = dirLength;
    };

    function initGL () {

            try {

                    if ( ! ( _gl = _canvas.getContext( 'experimental-webgl', { alpha: _alpha, premultipliedAlpha: _premultipliedAlpha, antialias: _antialias, stencil: _stencil, preserveDrawingBuffer: _preserveDrawingBuffer } ) ) ) {

                            throw 'Error creating WebGL context.';

                    }

            } catch ( error ) {

                    console.error( error );

            }

    };

    function setDefaultGLState () {

            _gl.clearColor( 0, 0, 0, 1 );
            _gl.clearDepth( 1 );
            _gl.clearStencil( 0 );

            _gl.enable( _gl.DEPTH_TEST );
            _gl.depthFunc( _gl.LEQUAL );

            _gl.frontFace( _gl.CCW );
            _gl.cullFace( _gl.BACK );
            _gl.enable( _gl.CULL_FACE );

            _gl.enable( _gl.BLEND );
            _gl.blendEquation( _gl.FUNC_ADD );
            _gl.blendFunc( _gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA );

            _gl.clearColor( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha );

    };
    
    this.addPostPlugin(new WebMol.SpritePlugin());
        
};

