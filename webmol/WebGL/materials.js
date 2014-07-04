/**
 * Line and Mesh material types
 * @constructor
 */
WebMol.Material = function () {

    WebMol.EventDispatcher.call( this );

    this.id = WebMol.MaterialIdCount ++;

    this.name = '';
    
    //TODO: Which of these instance variables can I remove??
    this.side = WebMol.FrontSide;

    this.opacity = 1;
    this.transparent = false;

    this.blending = WebMol.NormalBlending;

    this.depthTest = true;
    this.depthWrite = true;

    this.polygonOffset = false;
    this.polygonOffsetFactor = 0;
    this.polygonOffsetUnits = 0;

    this.alphaTest = 0;

    this.visible = true;

    this.needsUpdate = true;

};


WebMol.Material.prototype.setValues = function ( values ) {

    if ( values === undefined ) return;

    for ( var key in values ) {

        var newValue = values[ key ];

        if ( newValue === undefined ) {

            console.warn( 'WebMol.Material: \'' + key + '\' parameter is undefined.' );
            continue;

        }

        if ( key in this ) {

            var currentValue = this[ key ];

            if ( currentValue instanceof WebMol.Color && newValue instanceof WebMol.Color ) {

                currentValue.copy( newValue );

            } else if ( currentValue instanceof WebMol.Color ) {

                currentValue.set( newValue );

            } else if ( currentValue instanceof WebMol.Vector3 && newValue instanceof WebMol.Vector3 ) {

                currentValue.copy( newValue );

            } else {

                this[ key ] = newValue;

            }

        }

    }

};
//TODO: might want to look into blending equations
WebMol.Material.prototype.clone = function ( material ) {

    if ( material === undefined ) material = new WebMol.Material();

    material.name = this.name;

    material.side = this.side;

    material.opacity = this.opacity;
    material.transparent = this.transparent;

    material.blending = this.blending;

    material.depthTest = this.depthTest;
    material.depthWrite = this.depthWrite;

    material.polygonOffset = this.polygonOffset;
    material.polygonOffsetFactor = this.polygonOffsetFactor;
    material.polygonOffsetUnits = this.polygonOffsetUnits;

    material.alphaTest = this.alphaTest;

    material.overdraw = this.overdraw;

    material.visible = this.visible;

    return material;

};

WebMol.Material.prototype.dispose = function () {

    this.dispatchEvent( { type: 'dispose' } );

};

WebMol.MaterialIdCount = 0;

//Line basic material
/** @constructor */
WebMol.LineBasicMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    
    this.linewidth = 1;
    this.linecap = 'round';
    this.linejoin = 'round';
    
    this.vertexColors = false;
    
    this.fog = true;
    
    this.setValues(parameters);
    
};

WebMol.LineBasicMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.LineBasicMaterial.prototype.clone = function() {
  
    var material = new WebMol.LineBasicMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy();
    
};

//Mesh Lambert material
/** @constructor */
WebMol.MeshLambertMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    this.ambient = new WebMol.Color(0xfffff);
    this.emissive = new WebMol.Color(0x000000);
    
    //TODO: Which of these instance variables do I really need?
    this.wrapAround = false;
    this.wrapRGB = new WebMol.Vector3(1,1,1);
    
    this.map = null;
    
    this.lightMap = null;
    
    this.specularMap = null;
    
    this.envMap = null;
    this.reflectivity = 1;
    this.refractionRatio = 0.98;
    
    this.fog = true;
    
    this.wireframe = false;
    this.wireframeLinewidth = 1;
    this.wireframeLinecap = 'round';
    this.wireframeLinejoin = 'round';
    
    this.shading = WebMol.SmoothShading;
    
    this.vertexColors = WebMol.NoColors;
    
    this.skinning = false;
    
    this.setValues(parameters);
    
};

WebMol.MeshLambertMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.MeshLambertMaterial.prototype.clone = function() {
  
    var material = new WebMol.MeshLambertMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy(this.color);
    material.ambient.copy(this.ambient);
    material.emissive.copy(this.emissive);
    
    material.wrapAround = this.wrapAround;
    material.wrapRGB.copy(this.wrapRGB);
    
    material.map = this.map;
    
    material.lightMap = this.lightMap;
    
    material.specularMap = this.specularMap;
    
    material.envMap = this.envMap;
    material.combine = this.combine;
    material.reflectivity = this.reflectivity;
    material.refractionRatio = this.refractionRatio;
    
    material.fog = this.fog;
    
    material.shading = this.shading;
    
    material.vertexColors = this.vertexColors;
    
    material.skinning = this.skinning;
    material.morphTargets = this.morphTargets;
    material.morphNormals = this.morphNormals;
    
    return material;
    
};


//Sprite material
/** @constructor */
WebMol.SpriteMaterial = function(parameters) {
    
    WebMol.Material.call(this);
    
    this.color = new WebMol.Color(0xffffff);
    this.map = new WebMol.Texture();
    
    this.useScreenCoordinates = true;
    this.depthTest = !this.useScreenCoordinates;
    this.sizeAttenuation = !this.useScreenCoordinates;
    this.scaleByViewPort = !this.sizeAttenuation;
    this.alignment = WebMol.SpriteAlignment.center.clone();
    
    this.fog = false; // use scene fog
    
    this.uvOffset = new WebMol.Vector2(0, 0);
    this.uvScale = new WebMol.Vector2(1, 1);
    
    this.setValues(parameters);
    
    parameters = parameters || {};
    
    if (parameters.depthTest === undefined)
        this.depthTest = !this.useScreenCoordinates;
    if (parameters.sizeAttenuation === undefined)
        this.sizeAttenuation = !this.useScreenCoordinates;
    if (parameters.scaleByViewPort === undefined)
        this.scaleByViewPort = !this.sizeAttenuation;
    
};

WebMol.SpriteMaterial.prototype = Object.create(WebMol.Material.prototype);

WebMol.SpriteMaterial.prototype.clone = function() {
    
    var material = new WebMol.SpriteMaterial();
    
    WebMol.Material.prototype.clone.call(this, material);
    
    material.color.copy(this.color);
    material.map = this.map;
    
    material.useScreenCoordinates = useScreenCoordinates;
    material.sizeAttenuation = this.sizeAttenuation;
    material.scaleByViewport = this.scaleByViewPort;
    material.alignment.copy(this.alignment);
    
    material.uvOffset.copy(this.uvOffset);
    
    return material;
    
};

//Alignment for Sprites

WebMol.SpriteAlignment = {};
WebMol.SpriteAlignment.topLeft = new WebMol.Vector2(1, -1);
WebMol.SpriteAlignment.topCenter = new WebMol.Vector2(0, -1);
WebMol.SpriteAlignment.topRight = new WebMol.Vector2(-1, -1);
WebMol.SpriteAlignment.centerLeft = new WebMol.Vector2(1, 0);
WebMol.SpriteAlignment.center = new WebMol.Vector2(0, 0);
WebMol.SpriteAlignment.centerRight = new WebMol.Vector2(-1, 0);
WebMol.SpriteAlignment.bottomLeft = new WebMol.Vector2(1, 1);
WebMol.SpriteAlignment.bottomCenter = new WebMol.Vector2(0, 1);
WebMol.SpriteAlignment.bottomRight = new WebMol.Vector2(-1, 1);


//Texture
//We really only create textures from 2d rendering contexts (to display text labels)
/** @constructor */
WebMol.Texture = function(image) {

    WebMol.EventDispatcher.call(this);
    
    this.id = WebMol.TextureIdCount++;
    
    this.name = "";
    
    this.image = image;
    this.mipmaps = [];
    
    this.mapping = new WebMol.UVMapping();
    
    this.wrapS = WebMol.ClampToEdgeWrapping;
    this.wrapT = WebMol.ClampToEdgeWrapping;
    
    this.magFilter = WebMol.LinearFilter;
    this.minFilter = WebMol.LinearMipMapLinearFilter;
    
    this.anisotropy = 1;
    
    this.format = WebMol.RGBAFormat;
    this.type = WebMol.UnsignedByteType;
    
    this.offset = new WebMol.Vector2(0, 0);
    this.repeat = new WebMol.Vector2(1, 1);
    
    this.generateMipmaps = true;
    this.premultiplyAlpha = false;
    this.flipY = true;
    this.unpackAlignment = 4;
    
    this.needsUpdate = false;
    this.onUpdate = null;
    
};

WebMol.Texture.prototype = {

    constructor : WebMol.Texture,
    
    clone : function(texture) {
        
        if (texture === undefined)
            texture = new WebMol.Texture();
        
        texture.image = this.image;
        texture.mipmaps = this.mipmaps.slice(0);
        
        texture.mapping = this.mapping;
        
        texture.wrapS = this.wrapS;
        texture.wrapT = this.wrapT;
        
        texture.magFilter = this.magFilter;
        texture.minFilter = this.minFilter;
        
        texture.anisotropy = this.anisotropy;
        
        texture.format = this.format;
        texture.type = this.type;
        
        texture.offset.copy(this.offset);
        texture.repeat.copy(this.repeat);
        
        texture.generateMipmaps = this.generateMipmaps;
        texture.premultiplyAlpha = this.premultiplyAlpha;
        texture.flipY = this.flipY;
        texture.unpackAlignment = this.unpackAlignment;
        
        return texture;
        
    },
    
    dispose : function() {
        
        this.dispatchEvent( {type: 'dispose'});
        
    }    
    
};

WebMol.TextureIdCount = 0;


