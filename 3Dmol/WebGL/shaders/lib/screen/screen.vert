attribute vec2 vertexPosition;
varying highp vec2 vTexCoords;
const vec2 scale = vec2(0.5, 0.5);

void main() {
   vTexCoords  = vertexPosition * scale + scale; // scale vertex attribute to [0,1] range
   gl_Position = vec4(vertexPosition, 0.0, 1.0);
}
        