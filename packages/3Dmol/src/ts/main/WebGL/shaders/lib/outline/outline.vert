

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform float outlineWidth;
uniform float outlinePushback;

attribute vec3 position;
attribute vec3 normal;
attribute vec3 color;

void main() {

    vec4 norm = modelViewMatrix*vec4(normalize(normal),0.0);
    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );
    mvPosition.xy += norm.xy*outlineWidth;
    gl_Position = projectionMatrix * mvPosition;
    mvPosition.z -= outlinePushback; //go backwards in model space
    vec4 pushpos = projectionMatrix*mvPosition; //project to get z in projection space, I'm probably missing some simple math to do the same thing..
    gl_Position.z = gl_Position.w*pushpos.z/pushpos.w;
}

