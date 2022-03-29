/* eslint-disable no-cond-assign */
/**
 * Simplified webGL renderer
 */

$3Dmol.Renderer = function(parameters) {

    parameters = parameters || {};
    this.row = parameters.row;
    this.col = parameters.col;
    this.rows = parameters.rows;
    this.cols = parameters.cols;
    const _canvas = parameters.canvas !== undefined ? parameters.canvas
            : document.createElement('canvas');


    const _precision = parameters.precision !== undefined ? parameters.precision
            : 'highp'; const _alpha = parameters.alpha !== undefined ? parameters.alpha
            : true; const _premultipliedAlpha = parameters.premultipliedAlpha !== undefined ? parameters.premultipliedAlpha
            : true; const _antialias = parameters.antialias !== undefined ? parameters.antialias
            : false; const _stencil = parameters.stencil !== undefined ? parameters.stencil
            : true; const _preserveDrawingBuffer = parameters.preserveDrawingBuffer !== undefined ? parameters.preserveDrawingBuffer
            : false; const _clearColor = parameters.clearColor !== undefined ? new $3Dmol.Color(
            parameters.clearColor) : new $3Dmol.Color(0x000000);
             let _clearAlpha = parameters.clearAlpha !== undefined ? parameters.clearAlpha : 0;
            let _outlineMaterial = new $3Dmol.MeshOutlineMaterial(parameters.outline);
            let _outlineSphereImposterMaterial = new $3Dmol.SphereImposterOutlineMaterial(parameters.outline);
            let _outlineStickImposterMaterial = new $3Dmol.StickImposterOutlineMaterial(parameters.outline);
            let _outlineEnabled = !!parameters.outline
            ;
    this.domElement = _canvas;
    this.context = null;
    this.devicePixelRatio = 1.0; // set in setSize

    _canvas.id=parameters.id;
    this.autoClear = true;
    this.autoClearColor = true;
    this.autoClearDepth = true;
    this.autoClearStencil = true;

    // scene graph

    this.sortObjects = true;

    this.autoUpdateObjects = true;
    this.autoUpdateScene = true;

    // info

    this.info = {
        memory : {

            programs : 0,
            geometries : 0,
            textures : 0

        },
        render : {

            calls : 0,
            vertices : 0,
            faces : 0,
            points : 0

        }
    };

    // internal properties
    const _this = this;
    let _programs = []; let _programsCounter = 0;
    let _webglversion = 1;
    // internal state cache
    let _currentProgram = null;
    let _currentMaterialId = -1; let _currentGeometryGroupHash = null; let _currentCamera = null; let _geometryGroupCounter = 0;

    // GL state cache
    let _oldDoubleSided = -1; let _oldFlipSided = -1;
    let _oldBlending = -1;
    let _oldDepthTest = -1; let _oldDepthWrite = -1;
    const _oldPolygonOffset = null;
    let _oldLineWidth = null;

    let _viewportWidth = 0; let _viewportHeight = 0; let _currentWidth = 0; let _currentHeight = 0;
    const _enabledAttributes = {};

    // camera matrices cache
    const _projScreenMatrix = new $3Dmol.Matrix4();
    const _vector3 = new $3Dmol.Vector3();

    const _worldInverse = new $3Dmol.Matrix4();
    const _projInverse = new $3Dmol.Matrix4();
    const _textureMatrix = new $3Dmol.Matrix4();
    // light arrays cache
    const _direction = new $3Dmol.Vector3();
    let _lightsNeedUpdate = true;

    const _lights = {
        ambient : [ 0, 0, 0 ],
        directional : {
            length : 0,
            colors : [],
            positions : []
        },
        point : {
            length : 0,
            colors : [],
            positions : [],
            distances : []
        },
        spot : {
            length : 0,
            colors : [],
            positions : [],
            distances : [],
            directions : [],
            anglesCos : [],
            exponents : []
        },
        hemi : {
            length : 0,
            skyColors : [],
            groundColors : [],
            positions : []
        }

    };

    const sprites = new $3Dmol.SpritePlugin();

    // screensshader related variables
    let _screenshader = null;
    let _vertexattribpos = null;
    let _screenQuadVBO = null;

    // framebuffer variables
    let _fb = null;
    let _targetTexture = null;
    let _depthTexture = null;

    // initialize
    let _gl;

    initGL();
    setDefaultGLState();

    this.context = _gl;
    const _extInstanced = _gl.getExtension("ANGLE_instanced_arrays");
    const _extFragDepth = _gl.getExtension("EXT_frag_depth");
    const _extFloatLinear = _gl.getExtension('OES_texture_float_linear');
    const _extColorBufferFloat = _gl.getExtension('EXT_color_buffer_float');

    // API

    this.supportedExtensions = function() {
        return {supportsAIA: Boolean(_extInstanced),
            supportsImposters:  Boolean(_extFragDepth) || !isWebGL1()
            };
    };

    this.getContext = function() {
        return _gl;
    };

    this.isLost = function() {
        return _gl.isContextLost();
    };
    
    this.getPrecision = function() {
        return _precision;
    };

    this.setClearColorHex = function(hex, alpha) {
        _clearColor.setHex(hex);
        _clearAlpha = alpha;

        _gl.clearColor(_clearColor.r, _clearColor.g, _clearColor.b,
                        _clearAlpha);
    };

    this.enableOutline = function(parameters) {
        _outlineMaterial = new $3Dmol.MeshOutlineMaterial(parameters);
        _outlineSphereImposterMaterial = new $3Dmol.SphereImposterOutlineMaterial(parameters);
        _outlineStickImposterMaterial = new $3Dmol.StickImposterOutlineMaterial(parameters);
        _outlineEnabled = true;
    };

    this.disableOutline = function() {
        _outlineEnabled = false;
    };
    this.setViewport = function(){
        if(this.rows !== undefined && this.cols !== undefined && this.row !== undefined && this.col !== undefined){

            const wid = _canvas.width/this.cols;
            const hei = _canvas.height/this.rows;

            _viewportWidth =  wid;
            _viewportHeight = hei;

            _gl.enable(_gl.SCISSOR_TEST);
            _gl.scissor(wid*this.col,hei * this.row, wid, hei);
            _gl.viewport(wid * this.col , hei * this.row, wid, hei);

        }
    };

    this.setSize = function(width, height) {
        // zooming (in the browser) changes the pixel ratio and width/height
        this.devicePixelRatio = (window.devicePixelRatio !== undefined) ? window.devicePixelRatio : 1;
        // with antialiasing on (which doesn't seem to do much), render at double rsolution to eliminate jaggies
	      // my iphone crashes if we do though, so as a hacky workaround, don't do it with retina displays
        if(_antialias && this.devicePixelRatio < 2.0) this.devicePixelRatio *= 2.0;


        if(this.rows !== undefined && this.cols !== undefined && this.row !== undefined && this.col !== undefined){
            const wid = width/this.cols;
            const hei = height/this.rows;
            _canvas.width = width* this.devicePixelRatio;
            _canvas.height = height*this.devicePixelRatio;

            _viewportWidth =  wid * this.devicePixelRatio;
            _viewportHeight = hei * this.devicePixelRatio;

            _canvas.style.width = `${width  }px`;
            _canvas.style.height = `${height  }px`;

            this.setViewport();
        }else{
            _viewportWidth = _canvas.width = width * this.devicePixelRatio;
            _viewportHeight =  _canvas.height = height * this.devicePixelRatio;

            _canvas.style.width = `${width  }px`;
            _canvas.style.height = `${height  }px`;

            _gl.viewport(0, 0, _gl.drawingBufferWidth, _gl.drawingBufferHeight);
        }
        this.initFrameBuffer();
    };

    this.clear = function(color, depth, stencil) {

        let bits = 0;
        if (color === undefined || color)
            bits |= _gl.COLOR_BUFFER_BIT;
        if (depth === undefined || depth)
            bits |= _gl.DEPTH_BUFFER_BIT;
        if (stencil === undefined || stencil)
            bits |= _gl.STENCIL_BUFFER_BIT;
        _gl.clear(bits);

    };

    this.clearTarget = function(color, depth, stencil) {

        this.clear(color, depth, stencil);

    };

    this.setMaterialFaces = function(material, reflected) {

        const doubleSided = material.side === $3Dmol.DoubleSide;
        let flipSided = material.side === $3Dmol.BackSide;

        if(!material.imposter) // ignore reflection with imposters
            flipSided = reflected ? !flipSided : flipSided;

        if (_oldDoubleSided !== doubleSided) {

            if (doubleSided) {

                _gl.disable(_gl.CULL_FACE);

            } else {

                _gl.enable(_gl.CULL_FACE);

            }

            _oldDoubleSided = doubleSided;

        }

        if (_oldFlipSided !== flipSided) {

            if (flipSided) {

                _gl.frontFace(_gl.CW);

            } else {

                _gl.frontFace(_gl.CCW);

            }

            _oldFlipSided = flipSided;

        }

        _gl.cullFace(_gl.BACK);

    };

    this.setDepthTest = function(depthTest) {

        if (_oldDepthTest !== depthTest) {

            if (depthTest) {

                _gl.enable(_gl.DEPTH_TEST);

            } else {

                _gl.disable(_gl.DEPTH_TEST);

            }

            _oldDepthTest = depthTest;
        }
    };

    this.setDepthWrite = function(depthWrite) {

        if (_oldDepthWrite !== depthWrite) {

            _gl.depthMask(depthWrite);
            _oldDepthWrite = depthWrite;
        }
    };

    this.setBlending = function(blending) {

        if (!blending) {
            _gl.disable(_gl.BLEND);

        } else {
            _gl.enable(_gl.BLEND);
            _gl.blendEquationSeparate(_gl.FUNC_ADD, _gl.FUNC_ADD);
            _gl.blendFuncSeparate(_gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA,
                    _gl.ONE, _gl.ONE_MINUS_SRC_ALPHA);
        }
        _oldBlending = blending;
    };

    function enableAttribute(attribute) {

        if (!_enabledAttributes[attribute]) {
            _gl.enableVertexAttribArray(attribute);
            _enabledAttributes[attribute] = true;
        }
    }

    function disableAttributes() {

        for ( const attribute in _enabledAttributes) {

            if (_enabledAttributes[attribute]) {

                _gl.disableVertexAttribArray(attribute);
                _enabledAttributes[attribute] = false;
            }
        }
    }

    function setPolygonOffset(polygonOffset) {

        if (_oldPolygonOffset !== polygonOffset) {

            if (polygonOffset)
                _gl.enable(_gl.POLYGON_OFFSET_FILL);
            else
                _gl.disable(_gl.POLYGON_OFFSET_FILL);
        }
    }

    function setLineWidth(width) {

        if (width !== _oldLineWidth) {
            _gl.lineWidth(width);
            _oldLineWidth = width;
        }
    }

    const deallocateGeometry = function(geometry) {

        geometry.__webglInit = undefined;

        if (geometry.__webglVertexBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglVertexBuffer);

        if (geometry.__webglColorBuffer !== undefined)
            _gl.deleteBuffer(geometry.__webglColorBuffer);

        if (geometry.geometryGroups !== undefined) {

            for (let g = 0, gl = geometry.groups; g < gl; g++) {

                const geometryGroup = geometry.geometryGroups[g];

                if (geometryGroup.__webglVertexBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglVertexBuffer);

                if (geometryGroup.__webglColorBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglColorBuffer);

                if (geometryGroup.__webglNormalBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglNormalBuffer);

                if (geometryGroup.__webglFaceBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglFaceBuffer);

                if (geometryGroup.__webglLineBuffer !== undefined)
                    _gl.deleteBuffer(geometryGroup.__webglLineBuffer);

            }
        }
    };

    const deallocateMaterial = function(material) {

        const {program} = material;

        if (program === undefined)
            return;

        material.program = undefined;

        // only deallocate GL program if this was the last use of shared program
        // assumed there is only single copy of any program in the _programs
        // list
        // (that's how it's constructed)

        let i; let il; let programInfo;
        let deleteProgram = false;

        for (i = 0, il = _programs.length; i < il; i++) {

            programInfo = _programs[i];

            if (programInfo.program === program) {

                programInfo.usedTimes-=1;

                if (programInfo.usedTimes === 0) {

                    deleteProgram = true;

                }

                break;

            }

        }

        if (deleteProgram === true) {

            // avoid using array.splice, this is costlier than creating new
            // array from scratch

            const newPrograms = [];

            for (i = 0, il = _programs.length; i < il; i++) {

                programInfo = _programs[i];

                if (programInfo.program !== program) {

                    newPrograms.push(programInfo);

                }

            }

            _programs = newPrograms;

            _gl.deleteProgram(program);

            _this.info.memory.programs-=1;

        }

    };

    const deallocateTexture = function(texture) {

        if (texture.image && texture.image.__webglTextureCube) {

            // cube texture

            _gl.deleteTexture(texture.image.__webglTextureCube);

        }

        else {

            // 2D texture

            if (!texture.__webglInit)
                return;

            texture.__webglInit = false;
            _gl.deleteTexture(texture.__webglTexture);

        }

    };


    const onGeometryDispose = function(event) {

        const geometry = event.target;
        geometry.removeEventListener('dispose', onGeometryDispose);

        deallocateGeometry(geometry);

        _this.info.memory.geometries-=1;

    };

    const onTextureDispose = function(event) {

        const texture = event.target;

        texture.removeEventListener('dispose', onTextureDispose);

        deallocateTexture(texture);

        _this.info.memory.textures-=1;

    };

    const onMaterialDispose = function(event) {

        const material = event.target;
        material.removeEventListener('dispose', onMaterialDispose);

        deallocateMaterial(material);

    };

    // Compile and return shader
    function getShader(type, str) {

        let shader;

        if(!isWebGL1() && !str.startsWith("#version")) {
            // convert webgl1 to webgl2, unless a version is already explicit
            str = str.replace(/gl_FragDepthEXT/g,"gl_FragDepth");
            if(type === "fragment") {
                str = str.replace(/varying/g,"in");
            } else {
                str = str.replace(/varying/g,"out");
            }
            str = str.replace(/attribute/g,"in");
            str = str.replace(/texture2D/g,"texture");
            str = str.replace(/\/\/DEFINEFRAGCOLOR/g,'out vec4 glFragColor;');
            str = str.replace(/gl_FragColor/g,"glFragColor");
            str = `#version 300 es\n${str}`;
        }
        if (type === "fragment")
            shader = _gl.createShader(_gl.FRAGMENT_SHADER);
        else if (type === "vertex")
            shader = _gl.createShader(_gl.VERTEX_SHADER);

        _gl.shaderSource(shader, str);
        _gl.compileShader(shader);

        if (!_gl.getShaderParameter(shader, _gl.COMPILE_STATUS)) {

            console.error(_gl.getShaderInfoLog(shader));
            console.error("could not initialize shader");
            return null;

        }

        return shader;

    }

    // Compile appropriate shaders (if necessary) from source code and attach to
    // gl program.
    function buildProgram(fragmentShader, vertexShader, uniforms, parameters) {

        let p; let pl; let program; let code;
        const chunks = [];

        chunks.push(fragmentShader);
        chunks.push(vertexShader);

        for (p in parameters) {
            chunks.push(p);
            chunks.push(parameters[p]);
        }

        code = chunks.join();

        // check if program has already been compiled

        for (p = 0, pl = _programs.length; p < pl; p++) {

            const programInfo = _programs[p];

            if (programInfo.code === code) {

                programInfo.usedTimes+=1;

                return programInfo.program;
            }
        }

        // check if program requires webgl2
        if (isWebGL1()){
            if (parameters.volumetric)
                throw new Error("Volumetric rendering requires webgl2 which is not supported by your hardware.");
        }


        // Set up new program and compile shaders

        program = _gl.createProgram();

        // set up precision
        const precision = _precision;
        const prefix = `precision ${  precision  } float;`;

        const prefixVertex = [
                parameters.volumetric ? "#version 300 es" : "", prefix ]
                .join("\n");

        const prefixFragment = [
                parameters.volumetric ? "#version 300 es" : "",
                parameters.fragdepth && isWebGL1() ? "#extension GL_EXT_frag_depth: enable"
                        : "",
                parameters.wireframe ? "#define WIREFRAME 1" : "", prefix ]
                .join("\n");

        const glFragmentShader = getShader("fragment", prefixFragment
                + fragmentShader);
        const glVertexShader = getShader("vertex", prefixVertex + vertexShader);

        _gl.attachShader(program, glVertexShader);
        _gl.attachShader(program, glFragmentShader);

        _gl.linkProgram(program);

        if (!_gl.getProgramParameter(program, _gl.LINK_STATUS))
            console.error("Could not initialize shader");

        // gather and cache uniform variables and attributes

        program.uniforms = {};
        program.attributes = {};

        let identifiers; let u; let i;

        // uniform vars
        identifiers = [ 'viewMatrix', 'modelViewMatrix', 'projectionMatrix',
                'normalMatrix'];

        // custom uniform vars
        for (u in uniforms)
            identifiers.push(u);

        for (i = 0; i < identifiers.length; i++) {

            const uniformVar = identifiers[i];
            program.uniforms[uniformVar] = _gl.getUniformLocation(program,
                    uniformVar);

        }

        // attributes
        identifiers = [ 'position', 'normal', 'color', 'lineDistance',
                'offset', 'radius' ];

        /*
         * for (a in attributes) identifiers.push(a);
         */

        for (i = 0; i < identifiers.length; i++) {

            const attributeVar = identifiers[i];
            program.attributes[attributeVar] = _gl.getAttribLocation(program,
                    attributeVar);
        }

        program.id = _programsCounter+=1;
        _programs.push({
            program,
            code,
            usedTimes : 1
        });
        _this.info.memory.programs = _programs.length;

        return program;
    }

    // TODO: need to set up shader attributes and uniforms as attributes on
    // material object after attaching prgm
    // We need to attach appropriate uniform variables to material after shaders
    // have been chosen
    this.initMaterial = function(material) {

        material.addEventListener('dispose', onMaterialDispose);

        let parameters; let shaderID;

        shaderID = material.shaderID;

        if (shaderID) {

            const shader = $3Dmol.ShaderLib[shaderID];
            material.vertexShader = shader.vertexShader;
            material.fragmentShader = shader.fragmentShader;
            material.uniforms = $3Dmol.ShaderUtils.clone(shader.uniforms);
            // TODO: set material uniforms to shader uniform variables

        }

        parameters = {
            wireframe : material.wireframe,
            fragdepth : material.imposter,
            volumetric : material.volumetric
        };

        material.program = buildProgram(material.fragmentShader,
                material.vertexShader, material.uniforms, parameters);

    };

    function setProgram(camera, lights, fog, material, object, renderer) {

        if (material.needsUpdate) {

            if (material.program)
                deallocateMaterial(material);

            _this.initMaterial(material, lights, fog, object);
            material.needsUpdate = false;
        }

        let refreshMaterial = false;

        // pUniform: uniformVarName => uniformLocation
        // mUniform: uniformVarName => uniformJsVal
        const {program} = material; const pUniform = program.uniforms; const mUniform = material.uniforms;

        if (program !== _currentProgram) {
            _gl.useProgram(program);
            _currentProgram = program;

            refreshMaterial = true;
        }

        if (material.id !== _currentMaterialId) {
            _currentMaterialId = material.id;
            refreshMaterial = true;
        }

        if (camera !== _currentCamera) {
            _currentCamera = camera;
            refreshMaterial = true;
        }

        _gl.uniformMatrix4fv(pUniform.projectionMatrix, false,
                camera.projectionMatrix.elements);
        _gl.uniformMatrix4fv(pUniform.modelViewMatrix, false,
                object._modelViewMatrix.elements);
        _gl.uniformMatrix3fv(pUniform.normalMatrix, false,
                object._normalMatrix.elements);

        // Send projection matrix to uniform variable in shader
        if (refreshMaterial) {

            // Load projection, model-view matrices for perspective

            // Set up correct fog uniform vals
            mUniform.fogColor.value = fog.color;
            mUniform.fogNear.value = fog.near;
            mUniform.fogFar.value = fog.far;

            // Set up lights for lambert shader
            if (material.shaderID.startsWith("lambert")
                    || material.shaderID === "instanced"
                    || material.shaderID.endsWith("imposter")) {

                // load view and normal matrices for directional and object
                // lighting
                _gl.uniformMatrix4fv(pUniform.viewMatrix, false,
                        camera.matrixWorldInverse.elements);

                if (_lightsNeedUpdate) {
                    setupLights(program, lights);
                    _lightsNeedUpdate = false;
                }

                // Set up correct light uniform var vals
                mUniform.directionalLightColor.value = _lights.directional.colors;
                mUniform.directionalLightDirection.value = _lights.directional.positions;

            } else if (material.shaderID.endsWith("outline")) {
                mUniform.outlineColor.value = material.outlineColor;
                mUniform.outlineWidth.value = material.outlineWidth;
                mUniform.outlinePushback.value = material.outlinePushback;
            }  else if (material.shaderID === "volumetric") {

                // need a matrix that maps back from model coordinates to texture coordinates
                //  textureMat*modelInv*position
                object._modelViewMatrix.getScale(_direction); // scale factor of conversion
                _worldInverse.getInverse(object._modelViewMatrix);
                _projInverse.getInverse(camera.projectionMatrix);
                _textureMatrix.multiplyMatrices(object.material.texmatrix,_worldInverse);
                _gl.uniformMatrix4fv(pUniform.textmat, false, _textureMatrix.elements);
                _gl.uniformMatrix4fv(pUniform.projinv, false, _projInverse.elements);

                //  need the resolution (step size of ray in viewer coordinates)
                const invscale = Math.min(Math.min(_direction.x,_direction.y),_direction.z);
                mUniform.step.value = object.material.unit*invscale;
                mUniform.maxdepth.value = object.material.maxdepth*invscale;
                mUniform.transfermax.value = object.material.transfermax;
                mUniform.transfermin.value = object.material.transfermin;
                mUniform.subsamples.value = object.material.subsamples;



                renderer.setTexture(object.material.transferfn, 4, false);
                renderer.setTexture(object.material.map, 3, true);
                // depth texture from the renderbuffer, for volumetric integration with surfaces
                _gl.activeTexture(_gl.TEXTURE5);
                _gl.bindTexture(_gl.TEXTURE_2D, _depthTexture);

            }

            // opacity, diffuse, emissive, etc
            mUniform.opacity.value = material.opacity;

            // Load any other material specific uniform variables to gl shaders
            loadMaterialUniforms(pUniform, mUniform);

        }

        return program;

    }

    function loadMaterialUniforms(pUniform, mUniform) {
        let uniformVar; let type; let uniformVal; let uniformLoc;

        for (uniformVar in mUniform) {
            if (!pUniform[uniformVar])
                continue;

            type = mUniform[uniformVar].type;
            uniformVal = mUniform[uniformVar].value;
            uniformLoc = pUniform[uniformVar];

            // single float
            if (type === 'f')
                _gl.uniform1f(uniformLoc, uniformVal);
            // single integer
            else if (type === 'i')
                _gl.uniform1i(uniformLoc, uniformVal);
            // array of floats
            else if (type === 'fv')
                _gl.uniform3fv(uniformLoc, uniformVal);
            // color - r,g,b floats
            else if (type === 'c')
                _gl.uniform3f(uniformLoc, uniformVal.r, uniformVal.g,
                        uniformVal.b);
            else if (type === 'f4')
                _gl.uniform4f(uniformLoc, uniformVal[0], uniformVal[1],
                        uniformVal[2],uniformVal[3]);

        }

    }

    this.renderBuffer = function(camera, lights, fog, material, geometryGroup,
            object) {

        if (!material.visible)
            return;

        let program; let attributes;

        // Sets up proper vertex and fragment shaders and attaches them to webGL
        // program
        // Also sets appropriate uniform variables
        program = setProgram(camera, lights, fog, material, object, this);

        attributes = program.attributes;

        let updateBuffers = false; const wireframeBit = material.wireframe ? 1 : 0; const geometryGroupHash = (geometryGroup.id * 0xffffff)
                + (program.id * 2) + wireframeBit;

        if (geometryGroupHash !== _currentGeometryGroupHash) {
            _currentGeometryGroupHash = geometryGroupHash;
            updateBuffers = true;
        }

        // rebind shader attributes to appropriate (and already initialized) gl
        // buffers
        if (updateBuffers) {

            disableAttributes();

            // Vertices
            if (attributes.position >= 0) {
                _gl.bindBuffer(_gl.ARRAY_BUFFER,
                        geometryGroup.__webglVertexBuffer);
                enableAttribute(attributes.position);
                _gl.vertexAttribPointer(attributes.position, 3, _gl.FLOAT,
                        false, 0, 0);
            }

            // Colors
            if (attributes.color >= 0) {
                _gl.bindBuffer(_gl.ARRAY_BUFFER,
                        geometryGroup.__webglColorBuffer);
                enableAttribute(attributes.color);
                _gl.vertexAttribPointer(attributes.color, 3, _gl.FLOAT, false,
                        0, 0);
            }

            // Normals
            if (attributes.normal >= 0) {
                _gl.bindBuffer(_gl.ARRAY_BUFFER,
                        geometryGroup.__webglNormalBuffer);
                enableAttribute(attributes.normal);
                _gl.vertexAttribPointer(attributes.normal, 3, _gl.FLOAT, false,
                        0, 0);
            }

            // Offsets (Instanced only)
            if (attributes.offset >= 0) {
                _gl.bindBuffer(_gl.ARRAY_BUFFER,
                        geometryGroup.__webglOffsetBuffer);
                enableAttribute(attributes.offset);
                _gl.vertexAttribPointer(attributes.offset, 3, _gl.FLOAT, false,
                        0, 0);
            }

            // Radii (Instanced only)
            if (attributes.radius >= 0) {
                _gl.bindBuffer(_gl.ARRAY_BUFFER,
                        geometryGroup.__webglRadiusBuffer);
                enableAttribute(attributes.radius);
                _gl.vertexAttribPointer(attributes.radius, 1, _gl.FLOAT, false,
                        0, 0);
            }

        }

        // Render
        let faceCount; let lineCount;
        // lambert shaders - draw triangles
        // TODO: make sure geometryGroup's face count is setup correctly
        if (object instanceof $3Dmol.Mesh) {

            if (material.shaderID === "instanced") {
                const sphereGeometryGroup = material.sphere.geometryGroups[0];
                if (updateBuffers) {
                    _gl.bindBuffer(_gl.ARRAY_BUFFER,
                            geometryGroup.__webglVertexBuffer);
                    _gl.bufferData(_gl.ARRAY_BUFFER,
                            sphereGeometryGroup.vertexArray, _gl.STATIC_DRAW);
                    _gl.bindBuffer(_gl.ARRAY_BUFFER,
                            geometryGroup.__webglNormalBuffer);
                    _gl.bufferData(_gl.ARRAY_BUFFER,
                            sphereGeometryGroup.normalArray, _gl.STATIC_DRAW);
                    _gl.bindBuffer(_gl.ELEMENT_ARRAY_BUFFER,
                            geometryGroup.__webglFaceBuffer);
                    _gl.bufferData(_gl.ELEMENT_ARRAY_BUFFER,
                            sphereGeometryGroup.faceArray, _gl.STATIC_DRAW);
                }

                faceCount = sphereGeometryGroup.faceidx;

                _extInstanced.vertexAttribDivisorANGLE(attributes.offset, 1);
                _extInstanced.vertexAttribDivisorANGLE(attributes.radius, 1);
                _extInstanced.vertexAttribDivisorANGLE(attributes.color, 1);

                _extInstanced.drawElementsInstancedANGLE(_gl.TRIANGLES,
                        faceCount, _gl.UNSIGNED_SHORT, 0,
                        geometryGroup.radiusArray.length);

                _extInstanced.vertexAttribDivisorANGLE(attributes.offset, 0);
                _extInstanced.vertexAttribDivisorANGLE(attributes.radius, 0);
                _extInstanced.vertexAttribDivisorANGLE(attributes.color, 0);

            }

            else if (material.wireframe) {
                lineCount = geometryGroup.lineidx;
                setLineWidth(material.wireframeLinewidth);

                if (updateBuffers)
                    _gl.bindBuffer(_gl.ELEMENT_ARRAY_BUFFER,
                            geometryGroup.__webglLineBuffer);

                _gl.drawElements(_gl.LINES, lineCount, _gl.UNSIGNED_SHORT, 0);
            }

            else {
                faceCount = geometryGroup.faceidx;

                if (updateBuffers)
                    _gl.bindBuffer(_gl.ELEMENT_ARRAY_BUFFER,
                            geometryGroup.__webglFaceBuffer);
                _gl.drawElements(_gl.TRIANGLES, faceCount, _gl.UNSIGNED_SHORT,
                        0);

            }

            _this.info.render.calls+=1;
            _this.info.render.vertices += faceCount;
            _this.info.render.faces += faceCount / 3;
        }

        // basic shaders - draw lines
        else if (object instanceof $3Dmol.Line) {
            lineCount = geometryGroup.vertices;

            setLineWidth(material.linewidth);
            _gl.drawArrays(_gl.LINES, 0, lineCount);

            _this.info.render.calls+=1;
        }

    };

    // rendering
    function renderObjects(renderList, reverse, materialType, camera, lights,
            fog, useBlending) {

        let webglObject; let object; let buffer; let material; let start; let end; let delta;

        // Forward or backward render

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

        for (let i = start; i !== end; i += delta) {

            webglObject = renderList[i];

            if (webglObject.render) {

                object = webglObject.object;
                buffer = webglObject.buffer;
                material = webglObject[materialType];

                if (!material)
                    continue;

                if (useBlending)
                    _this.setBlending(true);

                _this.setDepthTest(material.depthTest);
                _this.setDepthWrite(material.depthWrite);
                setPolygonOffset(material.polygonOffset,
                        material.polygonOffsetFactor,
                        material.polygonOffsetUnits);

                const reflected = object._modelViewMatrix.isReflected();

                _this.setMaterialFaces(material, reflected);

                _this.renderBuffer(camera, lights, fog, material, buffer,
                        object);
                if (_outlineEnabled || material.outline) {
                    if(material.shaderID === 'sphereimposter') {
                        _this.renderBuffer(camera, lights, fog, _outlineSphereImposterMaterial,
                                buffer, object);
                    }
                    else if(material.shaderID === 'stickimposter') {
                        _this.renderBuffer(camera, lights, fog, _outlineStickImposterMaterial,
                                buffer, object);
                    }
                    else if(!material.wireframe
                        && material.shaderID !== 'basic'
                        && material.opacity !== 0.0) {
                        _this.renderBuffer(camera, lights, fog, _outlineMaterial,
                            buffer, object);
                    }
                }
            }
        }

    }

    this.render = function(scene, camera, forceClear) {

        if (camera instanceof $3Dmol.Camera === false) {

            console
                    .error('$3Dmol.Renderer.render: camera is not an instance of $3Dmol.Camera.');
            return;

        }

        let i; let il;

        let webglObject; let object; let renderList;

        const lights = scene.__lights; const {fog} = scene;

        // reset caching for this frame

        _currentMaterialId = -1;
        _lightsNeedUpdate = true;

        // update scene graph

        if (this.autoUpdateScene)
            scene.updateMatrixWorld();

        // update camera matrices
        // Pretty sure camera's parent is always going to be undefined for our
        // purposes...
        if (camera.parent === undefined)
            camera.updateMatrixWorld();

        camera.matrixWorldInverse.getInverse(camera.matrixWorld);

        _projScreenMatrix.multiplyMatrices(camera.projectionMatrix,
                camera.matrixWorldInverse);

        // update WebGL objects

        if (this.autoUpdateObjects)
            this.initWebGLObjects(scene);

        _this.info.render.calls = 0;
        _this.info.render.vertices = 0;
        _this.info.render.faces = 0;
        _this.info.render.points = 0;

        _currentWidth = _viewportWidth;
        _currentHeight = _viewportHeight;
        this.setViewport();
        this.setFrameBuffer();

        if (this.autoClear || forceClear) {
            _gl.clearColor(_clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha);
            this.clear(this.autoClearColor, this.autoClearDepth, this.autoClearStencil);
        }


        // set matrices for regular objects (frustum culled)

        renderList = scene.__webglObjects;
        let hasvolumetric = false;
        for (i = 0, il = renderList.length; i < il; i++) {

            webglObject = renderList[i];
            object = webglObject.object;

            webglObject.render = false;

            if (object.visible) {
                setupMatrices(object, camera);
                unrollBufferMaterial(webglObject);
                webglObject.render = true;
                if(webglObject.volumetric) hasvolumetric = true;
            }
        }

        // set matrices for immediate objects

        const material = null;

        // opaque pass (front-to-back order)

        this.setBlending(false);

        renderObjects(scene.__webglObjects, true, "opaque", camera, lights,
                fog, false, material);

        // Render embedded labels (sprites)
        renderSprites(scene, camera, false);

        // prime depth buffer
        renderObjects(scene.__webglObjects, true, "blank", camera, lights, fog,
                true, material);

        // transparent pass (back-to-front order)

        renderObjects(scene.__webglObjects, false, "transparent", camera,
                lights, fog, true, material);

        // volumetric is separate
        if(hasvolumetric && _fb) {
            // disconnect framebuffer to get depth texture
            this.reinitFrameBuffer();
            renderObjects(scene.__webglObjects, false, "volumetric", camera,
                    lights, fog, true, material);
        }


        this.renderFrameBuffertoScreen();
        this.setDepthTest(true);
        this.setDepthWrite(true);

        // Render floating labels (sprites)
        renderSprites(scene, camera, true);
    };

    function renderSprites(scene, camera, inFront) {

        // Reset state once regardless
        // This should also fix cartoon render bug (after transparent surface
        // render)

        _currentGeometryGroupHash = -1;
        _currentProgram = null;
        _currentCamera = null;
        _oldBlending = -1;
        _oldDepthWrite = -1;
        _oldDepthTest = -1;
        _oldDoubleSided = -1;
        _currentMaterialId = -1;
        _oldFlipSided = -1;
        _lightsNeedUpdate = true;

        sprites.render(scene, camera, _currentWidth, _currentHeight, inFront);

        // Reset state a
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

    // reinitialize framebuffer without the depth texture attached so we can read to it
    // do not allocate new textures
    this.reinitFrameBuffer = function() {
        // only needed/works with webgl2
        if (isWebGL1()) return;

        // Create and bind the framebuffer
        _fb = _gl.createFramebuffer();
        _gl.bindFramebuffer(_gl.FRAMEBUFFER, _fb);
        _gl.framebufferTexture2D(_gl.FRAMEBUFFER, _gl.COLOR_ATTACHMENT0, _gl.TEXTURE_2D, _targetTexture, 0);
    };

    // setup framebuffer for drawing into, assumes buffers already allocated
    this.setFrameBuffer = function() {
        if (isWebGL1() || !_fb) return;
        const width = _viewportWidth;
        const height = _viewportHeight;

        // when using framebuffer, always draw from origin, will shift the viewport when we render
        _gl.enable(_gl.SCISSOR_TEST);
        _gl.scissor(0,0, width, height);
        _gl.viewport(0,0, width, height);

        // color texture
        _gl.bindTexture(_gl.TEXTURE_2D, _targetTexture);
        _gl.texImage2D(_gl.TEXTURE_2D, 0, _gl.RGBA, width, height, 0, _gl.RGBA, _gl.UNSIGNED_BYTE, null);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_MIN_FILTER, _gl.LINEAR);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_MAG_FILTER, _gl.LINEAR);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE);

        // depth texture
        _gl.bindTexture(_gl.TEXTURE_2D, _depthTexture);
        _gl.texImage2D(_gl.TEXTURE_2D, 0, _gl.DEPTH_COMPONENT32F, width, height, 0, _gl.DEPTH_COMPONENT, _gl.FLOAT, null);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_MIN_FILTER, _gl.NEAREST);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_MAG_FILTER, _gl.NEAREST);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE);
        _gl.texParameteri(_gl.TEXTURE_2D, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE);

        // bind fb
        _gl.bindFramebuffer(_gl.FRAMEBUFFER, _fb);
        _gl.framebufferTexture2D(_gl.FRAMEBUFFER, _gl.COLOR_ATTACHMENT0, _gl.TEXTURE_2D, _targetTexture, 0);
        _gl.framebufferTexture2D(_gl.FRAMEBUFFER, _gl.DEPTH_ATTACHMENT,  _gl.TEXTURE_2D, _depthTexture, 0);

    };

    // allocate buffers for framebuffer, needs to be called with every resize
    this.initFrameBuffer = function() {
        // only needed/works with webgl2
        if (isWebGL1()) return;


        const width = _viewportWidth;
        const height = _viewportHeight;

        // when using framebuffer, always draw from origin, will shift the viewport when we render
        _gl.enable(_gl.SCISSOR_TEST);
        _gl.scissor(0,0, width, height);
        _gl.viewport(0,0, width, height);

        // create textures and frame buffer, will be initialized in setFrameBuffer
        _targetTexture = _gl.createTexture();
        _depthTexture = _gl.createTexture();
        _fb = _gl.createFramebuffer();

        // build screenshader
        const screenshader = $3Dmol.ShaderLib.screen;

        _screenshader = buildProgram(screenshader.fragmentShader,
            screenshader.vertexShader, screenshader.uniforms, {});
        _vertexattribpos = _gl.getAttribLocation(_screenshader, 'vertexPosition');
        // create the vertex array and attrib array for the full screenquad
        const verts = [
            // First triangle:
             1.0,  1.0,
            -1.0,  1.0,
            -1.0, -1.0,
            // Second triangle:
            -1.0, -1.0,
             1.0, -1.0,
             1.0,  1.0
        ];
        _screenQuadVBO = _gl.createBuffer();
        _gl.bindBuffer(_gl.ARRAY_BUFFER, _screenQuadVBO);
        _gl.bufferData(_gl.ARRAY_BUFFER, new Float32Array(verts), _gl.STATIC_DRAW);

    };

    this.renderFrameBuffertoScreen = function(){
        // only needed/works with webgl2
        if (isWebGL1() || _fb === null) return;

        this.setViewport(); // draw texture in correct viewport

        // bind default framebuffer
        _gl.bindFramebuffer(_gl.FRAMEBUFFER, null);
        _gl.clear(_gl.COLOR_BUFFER_BIT | _gl.DEPTH_BUFFER_BIT);
        _gl.frontFace(_gl.CCW);
        _gl.cullFace(_gl.BACK);

        // set screen shader and use it
        _gl.useProgram(_screenshader);
        _currentProgram = _screenshader;
        // disable depth test
        this.setDepthTest(-1);
        this.setDepthWrite(-1);

        // bind vertexarray buffer and texture
        _gl.bindBuffer(_gl.ARRAY_BUFFER, _screenQuadVBO);
        _gl.enableVertexAttribArray(_vertexattribpos);
        _gl.vertexAttribPointer(_vertexattribpos, 2, _gl.FLOAT, false, 0, 0);

        _gl.activeTexture(_gl.TEXTURE0);
        _gl.bindTexture(_gl.TEXTURE_2D, _targetTexture);

        // Draw 6 vertexes => 2 triangles:
        _gl.drawArrays(_gl.TRIANGLES, 0, 6);

    };


    this.initWebGLObjects = function(scene) {

        if (!scene.__webglObjects) {

            scene.__webglObjects = [];
            scene.__webglObjectsImmediate = [];
            scene.__webglSprites = [];
            scene.__webglFlares = [];

        }

        // Add objects; this sets up buffers for each geometryGroup
        if (scene.__objectsAdded.length) {

            while (scene.__objectsAdded.length) {
                addObject(scene.__objectsAdded[0], scene);
                scene.__objectsAdded.splice(0, 1);
            }

            // Force buffer update during render
            // Hackish fix for initial cartoon-render-then-transparent-surface
            // bug
            _currentGeometryGroupHash = -1;

        }

        while (scene.__objectsRemoved.length) {

            removeObject(scene.__objectsRemoved[0], scene);
            scene.__objectsRemoved.splice(0, 1);

        }

        // update must be called after objects adding / removal
        // This sends typed arrays to GL buffers for each geometryGroup
        for (let o = 0, ol = scene.__webglObjects.length; o < ol; o++) {

            updateObject(scene.__webglObjects[o].object);

        }

    };

    // Objects adding

    function addObject(object, scene) {

        let g; let gl; let geometry; let material; let geometryGroup;

        if (!object.__webglInit) {

            object.__webglInit = true;

            object._modelViewMatrix = new $3Dmol.Matrix4();
            object._normalMatrix = new $3Dmol.Matrix3();

            if (object.geometry !== undefined
                    && object.geometry.__webglInit === undefined) {

                object.geometry.__webglInit = true;
                object.geometry.addEventListener('dispose', onGeometryDispose);

            }

            if (object instanceof $3Dmol.Mesh || object instanceof $3Dmol.Line) {
                geometry = object.geometry;
                material = object.material;

                for (g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {

                    geometryGroup = geometry.geometryGroups[g];

                    geometryGroup.id = _geometryGroupCounter+=1;

                    // initialise VBO on the first access

                    if (!geometryGroup.__webglVertexBuffer) {

                        if (object instanceof $3Dmol.Mesh) {
                            createMeshBuffers(geometryGroup);
                            geometry.elementsNeedUpdate = true;
                            geometry.normalsNeedUpdate = true;
                        }

                        else if (object instanceof $3Dmol.Line)
                            createLineBuffers(geometryGroup);

                        geometry.verticesNeedUpdate = true;
                        geometry.colorsNeedUpdate = true;

                    }

                }

            }

        }

        if (!object.__webglActive) {

            if (object instanceof $3Dmol.Mesh || object instanceof $3Dmol.Line) {

                geometry = object.geometry;

                for (g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {
                    geometryGroup = geometry.geometryGroups[g];

                    addBuffer(scene.__webglObjects, geometryGroup, object);
                }

            }

            // Sprite
            else if (object instanceof $3Dmol.Sprite)
                scene.__webglSprites.push(object);

            object.__webglActive = true;

        }

    }

    function updateObject(object) {

        const {geometry} = object; let geometryGroup;

        if (object instanceof $3Dmol.Mesh || object instanceof $3Dmol.Line) {

            for (let g = 0, gl = geometry.geometryGroups.length; g < gl; g++) {

                geometryGroup = geometry.geometryGroups[g];

                if (geometry.verticesNeedUpdate || geometry.elementsNeedUpdate
                        || geometry.colorsNeedUpdate
                        || geometry.normalsNeedUpdate) {
                    setBuffers(geometryGroup, _gl.STATIC_DRAW);
                }
            }

            geometry.verticesNeedUpdate = false;
            geometry.elementsNeedUpdate = false;
            geometry.normalsNeedUpdate = false;
            geometry.colorsNeedUpdate = false;

            geometry.buffersNeedUpdate = false;

        }

    }

    function removeObject(object, scene) {

        if (object instanceof $3Dmol.Mesh || object instanceof $3Dmol.Line)
            removeInstances(scene.__webglObjects, object);

        else if (object instanceof $3Dmol.Sprite)
            removeInstancesDirect(scene.__webglSprites, object);

        object.__webglActive = false;

    }

    function removeInstances(objList, object) {

        for (let o = objList.length - 1; o >= 0; --o) {

            if (objList[o].object === object)
                objList.splice(o, 1);

        }
    }

    function removeInstancesDirect(objList, object) {

        for (let o = objList.length - 1; o >= 0; --o) {

            if (objList[o] === object)
                objList.splice(o, 1);

        }
    }

    function unrollBufferMaterial(globject) {

        const {object} = globject;
        const {material} = object;

        if (material.volumetric) {
            globject.opaque = null;
            globject.transparent = null;
            globject.volumetric = material;
        }
        else if (material.transparent) {
            globject.opaque = null;
            globject.volumetric = null;
            globject.transparent = material;
            if (!material.wireframe) {
                const blankMaterial = material.clone();
                blankMaterial.opacity = 0.0;
                globject.blank = blankMaterial;
            }
        } else {
            globject.opaque = material;
            globject.transparent = null;
            globject.volumetric = null;
        }

    }

    function setBuffers(geometryGroup, hint) {

        const {vertexArray} = geometryGroup;
        const {colorArray} = geometryGroup;

        // offset buffers
        if (geometryGroup.__webglOffsetBuffer !== undefined ) {
            _gl.bindBuffer(_gl.ARRAY_BUFFER, geometryGroup.__webglOffsetBuffer);
            _gl.bufferData(_gl.ARRAY_BUFFER, vertexArray, hint);
        }
        else {
            // normal, non-instanced case
            _gl.bindBuffer(_gl.ARRAY_BUFFER, geometryGroup.__webglVertexBuffer);
            _gl.bufferData(_gl.ARRAY_BUFFER, vertexArray, hint);
        }
        // color buffers
        _gl.bindBuffer(_gl.ARRAY_BUFFER, geometryGroup.__webglColorBuffer);
        _gl.bufferData(_gl.ARRAY_BUFFER, colorArray, hint);

        // normal buffers
        if (geometryGroup.normalArray
                && geometryGroup.__webglNormalBuffer !== undefined) {
            const {normalArray} = geometryGroup;
            _gl.bindBuffer(_gl.ARRAY_BUFFER, geometryGroup.__webglNormalBuffer);
            _gl.bufferData(_gl.ARRAY_BUFFER, normalArray, hint);

        }



        // radius buffers
        if (geometryGroup.radiusArray
                && geometryGroup.__webglRadiusBuffer !== undefined) {
            _gl.bindBuffer(_gl.ARRAY_BUFFER, geometryGroup.__webglRadiusBuffer);
            _gl.bufferData(_gl.ARRAY_BUFFER, geometryGroup.radiusArray, hint);
        }

        // face (index) buffers
        if (geometryGroup.faceArray
                && geometryGroup.__webglFaceBuffer !== undefined) {
            const {faceArray} = geometryGroup;
            _gl.bindBuffer(_gl.ELEMENT_ARRAY_BUFFER,
                    geometryGroup.__webglFaceBuffer);
            _gl.bufferData(_gl.ELEMENT_ARRAY_BUFFER, faceArray, hint);

        }

        // line (index) buffers (for wireframe)
        if (geometryGroup.lineArray
                && geometryGroup.__webglLineBuffer !== undefined) {
            const {lineArray} = geometryGroup;
            _gl.bindBuffer(_gl.ELEMENT_ARRAY_BUFFER,
                    geometryGroup.__webglLineBuffer);
            _gl.bufferData(_gl.ELEMENT_ARRAY_BUFFER, lineArray, hint);
        }

    }

    // Creates appropriate gl buffers for geometry chunk
    // TODO: do we need line buffer for mesh objects?
    // Also, can we integrate this with createLineBuffers?
    function createMeshBuffers(geometryGroup) {

        if (geometryGroup.radiusArray) {
            geometryGroup.__webglRadiusBuffer = _gl.createBuffer();
        }
        if(geometryGroup.useOffset) {
            geometryGroup.__webglOffsetBuffer = _gl.createBuffer();
        }
        geometryGroup.__webglVertexBuffer = _gl.createBuffer();
        geometryGroup.__webglNormalBuffer = _gl.createBuffer();
        geometryGroup.__webglColorBuffer = _gl.createBuffer();

        geometryGroup.__webglFaceBuffer = _gl.createBuffer();
        geometryGroup.__webglLineBuffer = _gl.createBuffer();

        _this.info.memory.geometries+=1;
    }

    function createLineBuffers(geometry) {

        geometry.__webglVertexBuffer = _gl.createBuffer();
        geometry.__webglColorBuffer = _gl.createBuffer();

        _this.info.memory.geometries+=1;
    }

    function addBuffer(objlist, buffer, object) {

        objlist.push({
            buffer,
            object,
            opaque : null,
            transparent : null
        });

    }

    function setupMatrices(object, camera) {

        object._modelViewMatrix.multiplyMatrices(camera.matrixWorldInverse,
                object.matrixWorld);

        object._normalMatrix.getInverse(object._modelViewMatrix);
        object._normalMatrix.transpose();

    }

    // Fallback filters for non-power-of-2 textures
    function filterFallback() {

        return _gl.LINEAR;
    }

    function setTextureParameters(textureType, texture) {

        if (textureType === _gl.TEXTURE_2D){
            _gl.texParameteri(textureType, _gl.TEXTURE_WRAP_S,
                    _gl.CLAMP_TO_EDGE);
            _gl.texParameteri(textureType, _gl.TEXTURE_WRAP_T,
                    _gl.CLAMP_TO_EDGE);
            _gl.texParameteri(textureType, _gl.TEXTURE_MAG_FILTER,
                    filterFallback(texture.magFilter));
            _gl.texParameteri(textureType, _gl.TEXTURE_MIN_FILTER,
                    filterFallback(texture.minFilter));

        } else { // 3Dtexture
            _gl.texParameteri(textureType, _gl.TEXTURE_WRAP_S,
                _gl.CLAMP_TO_EDGE);
            _gl.texParameteri(textureType, _gl.TEXTURE_WRAP_T,
                _gl.CLAMP_TO_EDGE);
            _gl.texParameteri(textureType, _gl.TEXTURE_WRAP_R,
                _gl.CLAMP_TO_EDGE);

            if(_extColorBufferFloat && _extFloatLinear) {
              // linear interpolation isn't supported by default (despite being the default??)
              _gl.texParameteri(textureType, _gl.TEXTURE_MAG_FILTER,_gl.LINEAR);
              _gl.texParameteri(textureType, _gl.TEXTURE_MIN_FILTER,_gl.LINEAR);
            } else {
              _gl.texParameteri(textureType, _gl.TEXTURE_MAG_FILTER,_gl.NEAREST);
              _gl.texParameteri(textureType, _gl.TEXTURE_MIN_FILTER,_gl.NEAREST);
            }

        }

    }

    this.getYRatio = function() {
      if(this.rows !== undefined && this.row !== undefined) return this.rows;
      return 1;
    };

    this.getXRatio = function() {
        if(this.cols !== undefined && this.col !== undefined) return this.cols;
        return 1;
    };

    this.getAspect = function(width,height){
        if(width === undefined || height === undefined){
            width = _canvas.width;
            height = _canvas.height;
        }
        let aspect = width/height;
        if(this.rows !== undefined && this.cols !== undefined && this.row !== undefined && this.col !== undefined){
            const wid = width/this.cols;
            const hei = height/this.rows;
            aspect = wid/hei;
        }
        return aspect;
    };

    this.setTexture = function(texture, slot, is3D) {

        if (texture.needsUpdate) {

            if (!texture.__webglInit) {
                texture.__webglInit = true;
                texture.addEventListener('dispose', onTextureDispose);
                texture.__webglTexture = _gl.createTexture();
                _this.info.memory.textures+=1;

            }

            _gl.activeTexture(_gl.TEXTURE0 + slot);
            const gltextureType = is3D ? _gl.TEXTURE_3D : _gl.TEXTURE_2D;
            _gl.bindTexture(gltextureType, texture.__webglTexture);
            _gl.pixelStorei(_gl.UNPACK_FLIP_Y_WEBGL, texture.flipY);
            _gl.pixelStorei(_gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, texture.premultiplyAlpha);
            _gl.pixelStorei(_gl.UNPACK_ALIGNMENT, texture.unpackAlignment);
            _gl.pixelStorei(_gl.PACK_ALIGNMENT, texture.unpackAlignment);

            const glFormat = paramToGL(texture.format); const glType = paramToGL(texture.type);

            if(!is3D) {
              const {image} = texture;
              let {width} = image; // might not be defined
              let {height} = image;
              if(typeof(width) === 'undefined') { // if no width,
                width = image.length;
                if(glFormat === _gl.RGBA) {
                    width /= 4; // each element takes up 4 bytes
                }
                height = 1;
              }
              setTextureParameters(_gl.TEXTURE_2D, texture);
              if(!isWebGL1()) { // webgl2
                _gl.texImage2D(_gl.TEXTURE_2D, 0, glFormat, width, height, 0, glFormat, glType, texture.image);
              } else {
                _gl.texImage2D(_gl.TEXTURE_2D, 0, glFormat, glFormat, glType, texture.image);
              }
            } else { // 3D
              setTextureParameters(_gl.TEXTURE_3D, texture);
              _gl.texImage3D(_gl.TEXTURE_3D, 0, _gl.R32F, texture.image.size.z, texture.image.size.y,
                  texture.image.size.x, 0, _gl.RED, _gl.FLOAT, texture.image.data);
            }

            texture.needsUpdate = false;

            if (texture.onUpdate)
                texture.onUpdate();

        } else {

            _gl.activeTexture(_gl.TEXTURE0 + slot);
            if (is3D)
                _gl.bindTexture(_gl.TEXTURE_3D, texture.__webglTexture);
            else
                _gl.bindTexture(_gl.TEXTURE_2D, texture.__webglTexture);

        }

    };

    // Map constants to WebGL constants

    function paramToGL(p) {

        if (p === $3Dmol.UnsignedByteType)
            return _gl.UNSIGNED_BYTE;
        if (p === $3Dmol.RGBAFormat)
            return _gl.RGBA;
        if (p === $3Dmol.NearestFilter)
            return _gl.NEAREST;

        return 0;

    }

    function setupLights(program, lights) {
        let l; let ll; let light; const r = 0; const g = 0; const b = 0; let color; let intensity; let distance;

        const zlights = _lights;

        const dirColors = zlights.directional.colors; const dirPositions = zlights.directional.positions;

        let dirCount = 0; let dirLength = 0; let dirOffset = 0;

        for (l = 0, ll = lights.length; l < ll; l++) {

            light = lights[l];

            color = light.color;
            intensity = light.intensity;
            distance = light.distance;

            if (light instanceof $3Dmol.Light) {

                dirCount+=1;

                _direction.getPositionFromMatrix(light.matrixWorld);
                _vector3.getPositionFromMatrix(light.target.matrixWorld);
                _direction.sub(_vector3);
                _direction.normalize();

                if (_direction.x === 0 && _direction.y === 0
                        && _direction.z === 0)
                    continue;

                dirPositions[dirOffset] = _direction.x;
                dirPositions[dirOffset + 1] = _direction.y;
                dirPositions[dirOffset + 2] = _direction.z;

                dirColors[dirOffset] = color.r * intensity;
                dirColors[dirOffset + 1] = color.g * intensity;
                dirColors[dirOffset + 2] = color.b * intensity;

                dirOffset += 3;

                dirLength+=1;
            }

        }

        zlights.ambient[0] = r;
        zlights.ambient[1] = g;
        zlights.ambient[2] = b;
        zlights.directional.length = dirLength;
    }

    function initGL() {
        // note setting antialis to true doesn't seem to do much and
        // causes problems on iOS Safari

        try {
            if (!(_gl = _canvas.getContext('webgl2', {
                alpha : _alpha,
                premultipliedAlpha : _premultipliedAlpha,
                antialias : false,
                stencil : _stencil,
                preserveDrawingBuffer : _preserveDrawingBuffer
            }))) {
                if (!(_gl = _canvas.getContext('experimental-webgl', {
                    alpha : _alpha,
                    premultipliedAlpha : _premultipliedAlpha,
                    antialias : false,
                    stencil : _stencil,
                    preserveDrawingBuffer : _preserveDrawingBuffer
                }))) {
                    if (!(_gl = _canvas.getContext('webgl', {
                        alpha : _alpha,
                        premultipliedAlpha : _premultipliedAlpha,
                        antialias : false,
                        stencil : _stencil,
                        preserveDrawingBuffer : _preserveDrawingBuffer
                    }))) {
                        throw new Error("Error creating WebGL context.");
                    }
                }
            }
            const vers = _gl.getParameter(_gl.VERSION);
            _webglversion = parseInt(vers[6]);
        } catch (error) {

            console.error(error);
        }

    }

    function isWebGL1() {
      return _webglversion === 1;
    }

    this.supportsVolumetric = function() { return !isWebGL1(); };

    function setDefaultGLState() {

        _gl.clearColor(0, 0, 0, 1);
        _gl.clearDepth(1);
        _gl.clearStencil(0);

        _gl.enable(_gl.DEPTH_TEST);
        _gl.depthFunc(_gl.LEQUAL);

        _gl.frontFace(_gl.CCW);
        _gl.cullFace(_gl.BACK);
        _gl.enable(_gl.CULL_FACE);

        _gl.enable(_gl.BLEND);
        _gl.blendEquation(_gl.FUNC_ADD);
        _gl.blendFunc(_gl.SRC_ALPHA, _gl.ONE_MINUS_SRC_ALPHA);

        _gl.clearColor(_clearColor.r, _clearColor.g, _clearColor.b,
                        _clearAlpha);
    }

    sprites.init(this);

};
