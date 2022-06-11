

uniform int useScreenCoordinates;
uniform vec3 screenPosition;
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform float rotation;
uniform vec2 scale;
uniform vec2 alignment;
uniform vec2 uvOffset;
uniform vec2 uvScale;

attribute vec2 position;
attribute vec2 uv;

varying vec2 vUV;

void main() {

    vUV = uvOffset + uv * uvScale;

    vec2 alignedPosition = position + alignment;

    vec2 rotatedPosition;
    rotatedPosition.x = ( cos(rotation) * alignedPosition.x - sin(rotation) * alignedPosition.y ) * scale.x;
    rotatedPosition.y = ( sin(rotation) * alignedPosition.x + cos(rotation) * alignedPosition.y ) * scale.y;

    vec4 finalPosition;

    if(useScreenCoordinates != 0) {
        finalPosition = vec4(screenPosition.xy + rotatedPosition, screenPosition.z, 1.0);
    }

    else {
        finalPosition = projectionMatrix * modelViewMatrix * vec4(0.0, 0.0, 0.0, 1.0); finalPosition /= finalPosition.w;
        finalPosition.xy += rotatedPosition; 
    }

    gl_Position = finalPosition;

}

