uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;

in vec3 position;
out vec4 mvPosition;
void main() {

    mvPosition = modelViewMatrix * vec4( position, 1.0 );
    gl_Position = projectionMatrix*mvPosition;
}
