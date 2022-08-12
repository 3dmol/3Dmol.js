

uniform float opacity;
uniform vec3 outlineColor;
uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
uniform mat4 projectionMatrix;
varying vec2 mapping;
varying float rval;
varying vec3 center;

uniform float outlinePushback;

//DEFINEFRAGCOLOR

void main() {
    float lensqr = dot(mapping,mapping);
    float rsqr = rval*rval;
    if(lensqr > rsqr)
       discard;
    float z = sqrt(rsqr-lensqr);
    vec3 cameraPos = center+ vec3(mapping.x,mapping.y,z-outlinePushback);
    vec4 clipPos = projectionMatrix * vec4(cameraPos, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepthEXT = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;
    gl_FragColor = vec4(outlineColor, 1 );
}


