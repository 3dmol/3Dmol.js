/* 
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

$3Dmol.ShaderUtils = {
    
    clone: function ( uniforms_src ) {
        
        var u, p, parameter, parameter_src, uniforms_clone = {};
        
        for (u in uniforms_src) {
            uniforms_clone[u] = {};
            uniforms_clone[u].type = uniforms_src[u].type;
            
            var srcValue = uniforms_src[u].value;
            
            if (srcValue instanceof $3Dmol.Color)
                uniforms_clone[u].value = srcValue.clone();
            else if (typeof srcValue === "number")
                uniforms_clone[u].value = srcValue;
            else if (srcValue instanceof Array) 
                uniforms_clone[u].value = [];
            else
                console.error("Error copying shader uniforms from ShaderLib: unknown type for uniform");
            
        }
        
        return uniforms_clone;
    }
};

$3Dmol.ShaderLib = { 
    'basic' : {
        fragmentShader : [                    
"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vColor;",

"void main() {",
    
"    gl_FragColor = vec4( vColor, opacity );",
    
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",    
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
    
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"
                                                     
].join("\n"),
        
        vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",

"attribute vec3 position;",
"attribute vec3 color;",

"varying vec3 vColor;",

"void main() {",

"    vColor = color;",
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
"    gl_Position = projectionMatrix * mvPosition;",

"}"
        
].join("\n"),
    
        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000}
        }

    },
    
 'sphereimposter' : {
        fragmentShader : [
"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",
"uniform mat4 projectionMatrix;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",
"uniform float uDepth;",
"uniform vec3 directionalLightColor[ 1 ];",

"varying vec3 vColor;",
"varying vec2 mapping;",
"varying float rval;",
"varying vec3 vLight;",
"varying vec3 center;",


"void main() {",
"    float lensqr = dot(mapping,mapping);",
"    float rsqr = rval*rval;",
"    if(lensqr > rsqr)",
"       discard;",
"    float z = sqrt(rsqr-lensqr);",
"    vec3 cameraPos = center+ vec3(mapping.x,mapping.y,z);",
"    vec4 clipPos = projectionMatrix * vec4(cameraPos, 1.0);",
"    float ndcDepth = clipPos.z / clipPos.w;",
"    gl_FragDepthEXT = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;",
"    vec3 norm = normalize(vec3(mapping.x,mapping.y,z));",
"    float dotProduct = dot( norm, vLight );",
"    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );",    
"    vec3 vLight = directionalLightColor[ 0 ] * directionalLightWeighting;",
"    gl_FragColor = vec4(vLight*vColor, opacity*opacity );", 
"    float fogFactor = smoothstep( fogNear, fogFar, gl_FragDepthEXT/gl_FragCoord.w );",
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",


"}"
                                                     
].join("\n"),
        
        vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 directionalLightColor[ 1 ];",
"uniform vec3 directionalLightDirection[ 1 ];",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"varying vec2 mapping;",
"varying vec3 vColor;",
"varying float rval;",
"varying vec3 vLight;",
"varying vec3 center;",

"void main() {",

"    vColor = color;",
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
"    center = mvPosition.xyz;",
"    vec4 projPosition = projectionMatrix * mvPosition;",
"    vec4 adjust = projectionMatrix* vec4(normal,0.0); adjust.z = 0.0; adjust.w = 0.0;",
"    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );",
"    vLight = normalize( lDirection.xyz );",
"    mapping = normal.xy;",
"    rval = abs(normal.x);",
"    gl_Position = projPosition+adjust;",

"}"
        
].join("\n"),
    
        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000},
            directionalLightColor: { type: 'fv', value: [] },
            directionalLightDirection: { type: 'fv', value: [] }
        }

    },
    
    
    'lambert' : { 
        fragmentShader : [

"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vLightFront;",
"varying vec3 vColor;",

"void main() {",
    
"    gl_FragColor = vec4( vec3 ( 1.0 ), opacity );",
    
"    #ifndef WIREFRAME",
"    gl_FragColor.xyz *= vLightFront;",
"    #endif",
    
"    gl_FragColor = gl_FragColor * vec4( vColor, opacity );",
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",
    
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
    
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"


].join("\n"),
       
       vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 directionalLightColor[ 1 ];",
"uniform vec3 directionalLightDirection[ 1 ];",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"varying vec3 vColor;",
"varying vec3 vLightFront;",

"void main() {",
    
"    vColor = color;",
    
"    vec3 objectNormal = normal;",  
"    vec3 transformedNormal = normalMatrix * objectNormal;",    
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
    
"    vLightFront = vec3( 0.0 );",
    
"    transformedNormal = normalize( transformedNormal );",
    
"    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );",
"    vec3 dirVector = normalize( lDirection.xyz );",
"    float dotProduct = dot( transformedNormal, dirVector );",
"    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );",
    
"    vLightFront += directionalLightColor[ 0 ] * directionalLightWeighting;",
    
"    gl_Position = projectionMatrix * mvPosition;",
"}"
           
].join("\n"),

        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000},
            directionalLightColor: { type: 'fv', value: [] },
            directionalLightDirection: { type: 'fv', value: [] }
        }

    },


    'instanced' : {
        fragmentShader : [

"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vLightFront;",
"varying vec3 vColor;",

"void main() {",

"    gl_FragColor = vec4( vec3 ( 1.0 ), opacity );",

"    #ifndef WIREFRAME",
"    gl_FragColor.xyz *= vLightFront;",
"    #endif",

"    gl_FragColor = gl_FragColor * vec4( vColor, opacity );",
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",

"    float fogFactor = smoothstep( fogNear, fogFar, depth );",

"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"


].join("\n"),

       vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 directionalLightColor[ 1 ];",
"uniform vec3 directionalLightDirection[ 1 ];",

"attribute vec3 offset;",
"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",
"attribute float radius;",

"varying vec3 vColor;",
"varying vec3 vLightFront;",

"void main() {",

"    vColor = color;",

"    vec3 objectNormal = normal;",
"    vec3 transformedNormal = normalMatrix * objectNormal;",
"    vec4 mvPosition = modelViewMatrix * vec4( position * radius + offset, 1.0 );",

"    vLightFront = vec3( 0.0 );",

"    transformedNormal = normalize( transformedNormal );",

"    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );",
"    vec3 dirVector = normalize( lDirection.xyz );",
"    float dotProduct = dot( transformedNormal, dirVector );",
"    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );",

"    vLightFront += directionalLightColor[ 0 ] * directionalLightWeighting;",

"    gl_Position = projectionMatrix * mvPosition;",
"}"

].join("\n"),

        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000},
            directionalLightColor: { type: 'fv', value: [] },
            directionalLightDirection: { type: 'fv', value: [] }
        }

    },
 
//for outline
     'outline' : { 
        fragmentShader : [

"uniform float opacity;",
"uniform vec3 outlineColor;",
"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"void main() {",
    
"    gl_FragColor = vec4( outlineColor, 1 );",
"}"


].join("\n"),
       
       vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform float outlineWidth;",
"uniform float outlinePushback;",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"void main() {",

"    vec4 norm = modelViewMatrix*vec4(normalize(normal),0.0);",
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
"    mvPosition.xy += norm.xy*outlineWidth;",
"    gl_Position = projectionMatrix * mvPosition;",
"    mvPosition.z -= outlinePushback;", //go backwards in model space
"    vec4 pushpos = projectionMatrix*mvPosition;", //project to get z in projection space, I'm probably missing some simple math to do the same thing..
"    gl_Position.z = gl_Position.w*pushpos.z/pushpos.w;",
"}"
           
].join("\n"),

        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            outlineColor: { type: 'c', value: new $3Dmol.Color(0.0, 0.0, 0.0) },
            fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000},           
            outlineWidth: { type: 'f', value: 0.1 },
            outlinePushback: { type: 'f', value: 1.0 },
        }

    },
//for outlining sphere imposter
    'sphereimposteroutline' : { 
       fragmentShader : [

"uniform float opacity;",
"uniform vec3 outlineColor;",
"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",
"uniform mat4 projectionMatrix;",
"varying vec2 mapping;",
"varying float rval;",
"varying vec3 center;",

"uniform float outlinePushback;",


"void main() {",
"    float lensqr = dot(mapping,mapping);",
"    float rsqr = rval*rval;",
"    if(lensqr > rsqr)",
"       discard;",
"    float z = sqrt(rsqr-lensqr);",
"    vec3 cameraPos = center+ vec3(mapping.x,mapping.y,z-outlinePushback);",
"    vec4 clipPos = projectionMatrix * vec4(cameraPos, 1.0);",
"    float ndcDepth = clipPos.z / clipPos.w;",
"    gl_FragDepthEXT = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;",
"    gl_FragColor = vec4(outlineColor, 1 );",
"}"


].join("\n"),
      
      vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform float outlineWidth;",
"uniform float outlinePushback;",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"varying vec2 mapping;",
"varying float rval;",
"varying vec3 center;",

"void main() {",

"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
"    center = mvPosition.xyz;",
"    vec4 projPosition = projectionMatrix * mvPosition;",
"    vec2 norm = normal.xy + vec2(sign(normal.x)*outlineWidth,sign(normal.y)*outlineWidth);", 
"    vec4 adjust = projectionMatrix* vec4(norm,normal.z,0.0); adjust.z = 0.0; adjust.w = 0.0;",
"    mapping = norm.xy;",
"    rval = abs(norm.x);",
"    gl_Position = projPosition+adjust;",
"}"
          
].join("\n"),

       uniforms : {
           opacity: { type: 'f', value: 1.0 },
           outlineColor: { type: 'c', value: new $3Dmol.Color(0.0, 0.0, 0.0) },
           fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
           fogNear: { type: 'f', value: 1.0 },
           fogFar: { type: 'f', value: 2000},           
           outlineWidth: { type: 'f', value: 0.1 },
           outlinePushback: { type: 'f', value: 1.0 },
       }

   },
   //for double sided lighting
   'stickimposter' : { 
       fragmentShader : [

"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vLightFront;",
"varying vec3 vLightBack;",
"varying vec3 cposition;",
"varying vec3 cnormal;",
"varying vec3 vColor;",
"varying float currR;",
"varying float radius;",

"void main() {",   
   
"    float mult = sqrt(radius*radius-currR*currR)/radius;",
"    float vLight = dot(normalize(-cnormal),vec3(0,0,1));",
"    gl_FragColor = vec4( vColor, opacity );",
"    gl_FragColor.xyz *= vLight;",
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",
   
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
   
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"


].join("\n"),
      
      vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 directionalLightColor[ 1 ];",
"uniform vec3 directionalLightDirection[ 1 ];",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"varying vec3 vColor;",
"varying vec3 vLightFront;",
"varying vec3 vLightBack;",
"varying vec3 cposition;",
"varying vec3 cnormal;",
"varying float currR;",
"varying float radius;",

"void main() {",
   
"    vColor = color;",
"    radius = length(normal);",
"    vec3 projpt = position+normal;",
"    vec3 transformedNormal = normalMatrix * normal;",    
"    vec3 zray = vec3(0.0,0.0,1.0);",
"    vec3 cr = cross(zray,transformedNormal);",
"    vec3 adjust = normalize(cr) * radius;",
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 )+vec4(adjust,0.0);",
"    vec4 position = projectionMatrix * mvPosition;",
"    currR = (cr.x < 0.0) ? -radius : radius;",
"    cnormal = (transformedNormal.z < 0.0) ? transformedNormal : -transformedNormal;",
"    gl_Position = position;",
"}"
          
].join("\n"),

       uniforms : {
           opacity: { type: 'f', value: 1.0 },
           fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
           fogNear: { type: 'f', value: 1.0 },
           fogFar: { type: 'f', value: 2000},           
           directionalLightColor: { type: 'fv', value: [] },
           directionalLightDirection: { type: 'fv', value: [] }
       }

   },
   
    //for double sided lighting
    'lambertdouble' : { 
        fragmentShader : [

"uniform mat4 viewMatrix;",
"uniform vec3 cameraPosition;",
"uniform float opacity;",

"uniform vec3 fogColor;",
"uniform float fogNear;",
"uniform float fogFar;",

"varying vec3 vLightFront;",
"varying vec3 vLightBack;",

"varying vec3 vColor;",

"void main() {",
    
"    gl_FragColor = vec4( vec3 ( 1.0 ), opacity );",
    
"    #ifndef WIREFRAME",
"    if ( gl_FrontFacing )",
"       gl_FragColor.xyz *= vLightFront;",
"    else",
"       gl_FragColor.xyz *= vLightBack;",
"    #endif",
    
"    gl_FragColor = gl_FragColor * vec4( vColor, opacity );",
"    float depth = gl_FragCoord.z / gl_FragCoord.w;",
    
"    float fogFactor = smoothstep( fogNear, fogFar, depth );",
    
"    gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );",

"}"


].join("\n"),
       
       vertexShader : [

"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform mat4 viewMatrix;",
"uniform mat3 normalMatrix;",
"uniform vec3 cameraPosition;",
"uniform vec3 directionalLightColor[ 1 ];",
"uniform vec3 directionalLightDirection[ 1 ];",

"attribute vec3 position;",
"attribute vec3 normal;",
"attribute vec3 color;",

"varying vec3 vColor;",
"varying vec3 vLightFront;",
"varying vec3 vLightBack;",

"void main() {",
    
"    vColor = color;",
    
"    vec3 objectNormal = normal;",  
"    vec3 transformedNormal = normalMatrix * objectNormal;",    
"    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );",
    
"    vLightFront = vec3( 0.0 );",
"    vLightBack = vec3( 0.0 );",
    
"    transformedNormal = normalize( transformedNormal );",
    
"    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );",
"    vec3 dirVector = normalize( lDirection.xyz );",
"    float dotProduct = dot( transformedNormal, dirVector );",
"    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );",
"    vec3 directionalLightWeightingBack = vec3( max( -dotProduct, 0.0 ) );",

"    vLightFront += directionalLightColor[ 0 ] * directionalLightWeighting;",
"    vLightBack += directionalLightColor[ 0 ] * directionalLightWeightingBack;",

"    gl_Position = projectionMatrix * mvPosition;",
"}"
           
].join("\n"),

        uniforms : {
            opacity: { type: 'f', value: 1.0 },
            fogColor: { type: 'c', value: new $3Dmol.Color(1.0, 1.0, 1.0) },
            fogNear: { type: 'f', value: 1.0 },
            fogFar: { type: 'f', value: 2000},           
            directionalLightColor: { type: 'fv', value: [] },
            directionalLightDirection: { type: 'fv', value: [] }
        }

    },
    
    
    'sprite': {
        
        fragmentShader: [
                                                         
"uniform vec3 color;",
"uniform sampler2D map;",
"uniform float opacity;",

"uniform int fogType;",
"uniform vec3 fogColor;",
"uniform float fogDensity;",
"uniform float fogNear;",
"uniform float fogFar;",
"uniform float alphaTest;",

"varying vec2 vUV;",

"void main() {",
    
"    vec4 texture = texture2D(map, vUV);",
    
"    if (texture.a < alphaTest) discard;",
    
"    gl_FragColor = vec4(color * texture.xyz, texture.a * opacity);",
    
"    if (fogType > 0) {",
        
"        float depth = gl_FragCoord.z / gl_FragCoord.w;",
"        float fogFactor = 0.0;",
        
"        if (fogType == 1) {",
"            fogFactor = smoothstep(fogNear, fogFar, depth);",
"        }",
        
"        else {",
"            const float LOG2 = 1.442695;",
"            float fogFactor = exp2(- fogDensity * fogDensity * depth * depth * LOG2);",
"            fogFactor = 1.0 - clamp(fogFactor, 0.0, 1.0);",
"        }",
        
"        gl_FragColor = mix(gl_FragColor, vec4(fogColor, gl_FragColor.w), fogFactor);",
        
"    }",
"}"                                              
            
].join("\n"),
        
        vertexShader: [

"uniform int useScreenCoordinates;",
"uniform vec3 screenPosition;",
"uniform mat4 modelViewMatrix;",
"uniform mat4 projectionMatrix;",
"uniform float rotation;",
"uniform vec2 scale;",
"uniform vec2 alignment;",
"uniform vec2 uvOffset;",
"uniform vec2 uvScale;",

"attribute vec2 position;",
"attribute vec2 uv;",

"varying vec2 vUV;",

"void main() {",
    
"    vUV = uvOffset + uv * uvScale;",
    
"    vec2 alignedPosition = position + alignment;",
    
"    vec2 rotatedPosition;",
"    rotatedPosition.x = ( cos(rotation) * alignedPosition.x - sin(rotation) * alignedPosition.y ) * scale.x;",
"    rotatedPosition.y = ( sin(rotation) * alignedPosition.x + cos(rotation) * alignedPosition.y ) * scale.y;",
    
"    vec4 finalPosition;",
    
"    if(useScreenCoordinates != 0) {",
"        finalPosition = vec4(screenPosition.xy + rotatedPosition, screenPosition.z, 1.0);",
"    }",
    
"    else {",
"        finalPosition = projectionMatrix * modelViewMatrix * vec4(0.0, 0.0, 0.0, 1.0); finalPosition /= finalPosition.w;",
"        finalPosition.xy += rotatedPosition; ",
"    }",
    
"    gl_Position = finalPosition;",
    
"}"
       
].join("\n"),

        uniforms : {
            
        }
        
    }
    
};
