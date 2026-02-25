
uniform float opacity;
uniform mat4 projectionMatrix;

uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;

varying vec3 vColor;
varying vec2 mapping;
varying float vOuterR;
varying float vMinorR;
varying vec3 vLight;
varying vec3 center;
varying vec3 vRingNormal;

//DEFINEFRAGCOLOR

void main() {
    float majorR = vOuterR - vMinorR;

    // The fragment represents a view ray at screen position (center.xy + mapping).
    // Intersect this ray with the ring plane to find the in-plane distance.
    // Ring plane: dot(P - center, vRingNormal) = 0
    // Ray: P = (center.x + mapping.x, center.y + mapping.y, center.z + t)
    // Solution: mapping.x*nx + mapping.y*ny + t*nz = 0  =>  t = -a/nz
    float a = mapping.x * vRingNormal.x + mapping.y * vRingNormal.y;
    float nz = vRingNormal.z;

    // Clamp |nz| to avoid singularity for edge-on rings
    float nzSafe = nz;
    if (abs(nzSafe) < 0.01) nzSafe = (nzSafe >= 0.0) ? 0.01 : -0.01;
    float zHit = -a / nzSafe;

    // In-ring-plane distance from center
    // hitPoint = (mapping.x, mapping.y, zHit) is in the ring plane by construction
    float dFromCenterSq = mapping.x * mapping.x + mapping.y * mapping.y + zHit * zHit;
    float dFromCenter = sqrt(dFromCenterSq);

    // Distance from ring centerline (in ring plane)
    float dRadial = dFromCenter - majorR;

    // Tube test
    if (dRadial * dRadial > vMinorR * vMinorR)
        discard;

    // Tube cross-section height
    float tubeZ = sqrt(max(0.0, vMinorR * vMinorR - dRadial * dRadial));

    // Always render the front surface (closest to camera).
    // The tube extends ±tubeZ along vRingNormal from the ring plane.
    // When vRingNormal.z > 0, +tubeZ is toward camera; when < 0, -tubeZ is.
    float tubeSign = (vRingNormal.z >= 0.0) ? 1.0 : -1.0;

    // Direction from center to the ring-plane hit point (normalized)
    vec3 hitPoint = vec3(mapping.x, mapping.y, zHit);
    vec3 radDir = (dFromCenter > 0.0001) ? hitPoint / dFromCenter : vec3(1.0, 0.0, 0.0);

    // Surface normal: radial component + tube bulge toward viewer
    vec3 norm = normalize(dRadial * radDir + tubeSign * tubeZ * vRingNormal);

    // Surface position: ring-plane intersection + tube bulge toward viewer
    vec3 surfacePos = center + hitPoint + tubeSign * tubeZ * vRingNormal;

    vec4 clipPos = projectionMatrix * vec4(surfacePos, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepthEXT = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

    float dotProduct = dot( norm, vLight );
    vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );
    vec3 vLightResult = directionalLightWeighting;
    vec3 color = vLightResult * vColor;
    gl_FragColor = vec4(color, opacity * opacity);

    if(fogNear != fogFar) {
        float depth = -surfacePos.z;
        float fogFactor = smoothstep( fogNear, fogFar, depth );
        gl_FragColor = mix( gl_FragColor, vec4( fogColor, gl_FragColor.w ), fogFactor );
    }

}
