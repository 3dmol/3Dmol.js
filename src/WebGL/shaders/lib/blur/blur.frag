const float INV_TOTAL_SAMPLES_FACTOR = 1.0 / 9.0;
uniform highp sampler2D inTex;
varying highp vec2 vTexCoords;
    
void main() {
 
 vec2 texelSize = 1.0 / vec2(textureSize(inTex,0));
 float blurred_visibility_factor = 0.0f;
 for (int t = -1; t <= 1; ++t) {
  for (int s = -1; s <= 1; ++s) {
   vec2 offset = vec2(float(s), float(t)) * texelSize;
   blurred_visibility_factor += texture2D(inTex, vTexCoords + offset).r;
  }
 }
    
 gl_FragDepthEXT = blurred_visibility_factor * INV_TOTAL_SAMPLES_FACTOR;
 
}