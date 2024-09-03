

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform float outlineWidth;
uniform float outlinePushback;
uniform float outlineMaxPixels;
uniform float vWidth;
uniform float vHeight;

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

    vec4 adjust = projectionMatrix* vec4(norm,normal.z,1.0); 
    mapping = norm.xy;
    rval = abs(norm.x);
    gl_Position = projPosition+vec4(adjust.xy,0.0,0.0);

    if(outlineMaxPixels > 0.0) {
        vec4 unadjusted = projectionMatrix*vec4(center.x+normal.x, center.y,center.z,1.0); 
        vec4 ccoord = projectionMatrix*vec4(center.xyz,1.0);
        adjust = projectionMatrix* vec4(center.x+norm.x,center.y,center.z,1.0); 
        //subtract center 
        unadjusted.xyz -= ccoord.xyz;
        adjust.xyz -= ccoord.xyz;
        unadjusted /= unadjusted.w;
        adjust /= adjust.w;
        float diff = abs(adjust.x-unadjusted.x);
        diff *= vWidth;
        if(diff > outlineMaxPixels) {
            
            float fixlen = abs(unadjusted.x) + outlineMaxPixels/vWidth;
            //adjsut reval by ratio of lengths
            rval *= fixlen/abs(adjust.x);
        }

    }
}

