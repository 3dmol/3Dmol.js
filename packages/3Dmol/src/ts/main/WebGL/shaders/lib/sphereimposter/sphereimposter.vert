uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;
uniform mat3 normalMatrix;
uniform vec3 directionalLightColor[ 1 ];
uniform vec3 directionalLightDirection[ 1 ];

attribute vec3 position;
attribute vec3 normal;
attribute vec3 color;

varying vec2 mapping;
varying vec3 vColor;
varying float rval;
varying vec3 vLight;
varying vec3 center;

void main() {

    vColor = color;
    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );
    center = mvPosition.xyz;
    vec4 projPosition = projectionMatrix * mvPosition;
    vec4 adjust = projectionMatrix* vec4(normal,0.0); adjust.z = 0.0; adjust.w = 0.0;
    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );
    vLight = normalize( lDirection.xyz );
    mapping = normal.xy;
    rval = abs(normal.x);
    gl_Position = projPosition+adjust;

}
