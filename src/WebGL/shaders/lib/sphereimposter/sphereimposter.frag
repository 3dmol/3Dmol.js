
uniform mat4 viewMatrix;
uniform float opacity;
uniform mat4 projectionMatrix;

uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
uniform float uDepth;
uniform vec3 directionalLightColor[ 1 ];

varying vec3 vColor;
varying float vAlpha;
varying vec2 mapping;
varying float rval;
varying vec3 vLight;
varying vec3 center;

#ifdef SHADED
uniform highp sampler2D shading;
#endif

//DEFINEFRAGCOLOR

void main() {
    float lensqr = dot(mapping,mapping);
    float rsqr = rval*rval;
    if(lensqr > rsqr)
       discard;
    float z = sqrt(rsqr-lensqr);
    vec3 cameraPos = center+ vec3(mapping.x,mapping.y,z);
    vec4 clipPos = projectionMatrix * vec4(cameraPos, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepthEXT = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;
    vec3 norm = normalize(vec3(mapping.x,mapping.y,z));
    float dotProduct = dot( norm, vLight );
    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );
    vec3 vLight = directionalLightColor[ 0 ] * directionalLightWeighting;
    vec3 color = vLight*vColor;
#ifdef SHADED
    ivec2 dim = textureSize(shading,0);
    float shadowFactor = texture2D(shading,vec2(gl_FragCoord.x/float(dim.x),gl_FragCoord.y/float(dim.y))).r;
    color *= shadowFactor;
#endif    
    // Per-vertex alpha (issue #166) applies LINEARLY: a user styling opacity:0.25 on an atom
    // expects a quarter-strength ghost. Only the material opacity keeps this shader's legacy
    // squared response curve, so whole-model opacity renders bit-identical to before
    // (vAlpha defaults to 1.0 via the attribute's generic value when no alphaArray exists).
    gl_FragColor = vec4(color, opacity*opacity*vAlpha );

    if(fogNear != fogFar) {
        float depth = -cameraPos.z;
        float fogFactor = smoothstep( fogNear, fogFar, depth );
        gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );
    }
     
}

