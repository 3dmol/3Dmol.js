
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
    // Per-vertex alpha (issue #166) applies LINEARLY: an atom styled opacity:0.5 renders at
    // 0.5. The material opacity keeps this shader's historical squared response, so existing
    // whole-model translucency is bit-identical to before -- vAlpha is 1.0 there (the
    // attribute's generic value, since geometries without an alphaArray never bind one), and
    // this collapses to the original opacity*opacity.
    //
    // The two therefore disagree: whole-model 0.5 renders like 0.25, per-atom 0.5 renders at
    // 0.5. That is deliberate. Squaring vAlpha to match would make the new API inherit a
    // legacy quirk, and unsquaring the material would change every committed reference image.
    gl_FragColor = vec4(color, opacity*opacity*vAlpha );

    if(fogNear != fogFar) {
        float depth = -cameraPos.z;
        float fogFactor = smoothstep( fogNear, fogFar, depth );
        gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );
    }
     
}

