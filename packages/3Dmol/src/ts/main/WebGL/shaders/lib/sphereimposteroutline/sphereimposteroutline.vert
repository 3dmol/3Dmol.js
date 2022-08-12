

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform float outlineWidth;
uniform float outlinePushback;

attribute vec3 position;
attribute vec3 normal;
attribute vec3 color;

varying vec2 mapping;
varying float rval;
varying vec3 center;

void main() {

    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );
    center = mvPosition.xyz;
    vec4 projPosition = projectionMatrix * mvPosition;
    vec2 norm = normal.xy + vec2(sign(normal.x)*outlineWidth,sign(normal.y)*outlineWidth);
    vec4 adjust = projectionMatrix* vec4(norm,normal.z,0.0); adjust.z = 0.0; adjust.w = 0.0;
    mapping = norm.xy;
    rval = abs(norm.x);
    gl_Position = projPosition+adjust;
}

