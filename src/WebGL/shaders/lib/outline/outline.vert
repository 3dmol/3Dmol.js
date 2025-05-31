

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform float outlineWidth;
uniform float outlinePushback;
uniform float vWidth;
uniform float vHeight;
uniform float outlineMaxPixels;

attribute vec3 position;
attribute vec3 normal;
attribute vec3 color;

varying vec4 mvPosition;

void main() {

    vec4 norm = modelViewMatrix*vec4(normalize(normal),0.0);
    mvPosition = modelViewMatrix * vec4( position, 1.0 );
    mvPosition.xy += norm.xy*outlineWidth;
    vec4 outpos = projectionMatrix * mvPosition;

    if(outlineMaxPixels > 0.0) {
        vec4 unadjusted = projectionMatrix*modelViewMatrix * vec4( position, 1.0 );
        float w = outpos.w;
        //normalize homogeneous coords
        unadjusted /= unadjusted.w;
        outpos /= outpos.w;
        vec2 diff = outpos.xy-unadjusted.xy;
        //put into pixels
        diff.x *= vWidth;
        diff.y *= vHeight;
        if ( length(diff) > outlineMaxPixels) {
            vec2 ndiff = normalize(diff)*outlineMaxPixels;
            ndiff.x /= vWidth;
            ndiff.y /= vHeight;
            outpos.xy = unadjusted.xy;
            outpos.xy += ndiff;
        }
        outpos *= w; //if I don't do this things blow up
    }
    gl_Position = outpos;
    mvPosition.z -= outlinePushback; //go backwards in model space
    vec4 pushpos = projectionMatrix*mvPosition; //project to get z in projection space, I'm probably missing some simple math to do the same thing..
    gl_Position.z = gl_Position.w*pushpos.z/pushpos.w;
}

