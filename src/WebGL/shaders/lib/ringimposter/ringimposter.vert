uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;
uniform mat3 normalMatrix;
uniform vec3 directionalLightColor[ 1 ];
uniform vec3 directionalLightDirection[ 1 ];

attribute vec3 position;
attribute vec3 normal;
attribute vec3 color;
attribute float radius;

varying vec2 mapping;
varying vec3 vColor;
varying float vOuterR;
varying float vMinorR;
varying vec3 vLight;
varying vec3 center;
varying vec3 vRingNormal;

void main() {

    // Decode: normal direction = ring plane normal, length = minorR
    vMinorR = length(normal);
    vec3 ringNormalObj = normal / max(vMinorR, 0.0001);

    // Decode: |radius| = outerR, sign(radius) = X billboard corner
    vOuterR = abs(radius);
    float xSign = sign(radius);

    // Decode: sign(color.b) = Y billboard corner, actual color = (r, g, |b|)
    float ySign = sign(color.b);
    vColor = vec3(color.r, color.g, abs(color.b));

    // Transform ring center to view space
    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );
    center = mvPosition.xyz;

    // Transform ring normal to view space
    vRingNormal = normalize( normalMatrix * ringNormalObj );

    // Billboard expansion
    mapping = vec2(xSign * vOuterR, ySign * vOuterR);
    vec4 projPosition = projectionMatrix * mvPosition;
    vec4 adjust = projectionMatrix * vec4(mapping, 0.0, 0.0);
    adjust.z = 0.0; adjust.w = 0.0;

    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );
    vLight = normalize( lDirection.xyz );
    gl_Position = projPosition + adjust;

}
