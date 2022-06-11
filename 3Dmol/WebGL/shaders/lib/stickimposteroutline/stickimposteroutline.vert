

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;
uniform mat3 normalMatrix;
uniform vec3 directionalLightColor[ 1 ];
uniform vec3 directionalLightDirection[ 1 ];
uniform vec3 outlineColor;
uniform float outlineWidth;
uniform float outlinePushback;


attribute vec3 position;
attribute vec3 normal;
attribute vec3 color;
attribute float radius;

varying vec3 vColor;
varying vec3 vLight;
varying vec3 cposition;
varying vec3 p1;
varying vec3 p2;
varying float r;

void main() {

    vColor = outlineColor;
    float rad = radius+sign(radius)*outlineWidth;
    r = abs(rad);
    vec4 to = modelViewMatrix*vec4(normal, 1.0); //normal is other point of cylinder
    vec4 pt = modelViewMatrix*vec4(position, 1.0);
//pushback
    to.xyz += normalize(to.xyz)*outlinePushback;
    pt.xyz += normalize(pt.xyz)*outlinePushback;

    vec4 mvPosition = pt;
    p1 = pt.xyz; p2 = to.xyz;
    vec3 norm = to.xyz-pt.xyz;
    float mult = 1.1; //slop to account for perspective of sphere
    if(length(p1) > length(p2)) { //billboard at level of closest point
       mvPosition = to;
    }
    vec3 n = normalize(mvPosition.xyz);
//intersect with the plane defined by the camera looking at the billboard point
    if(color.z >= 0.0) { //p1
       vec3 pnorm = normalize(p1);
       float t = dot(mvPosition.xyz-p1,n)/dot(pnorm,n);
       mvPosition.xyz = p1+t*pnorm;
    } else {
       vec3 pnorm = normalize(p2);
       float t = dot(mvPosition.xyz-p2,n)/dot(pnorm,n);
       mvPosition.xyz = p2+t*pnorm;
       mult *= -1.0;
    }
    vec3 cr = normalize(cross(mvPosition.xyz,norm))*rad;
    vec3 doublecr = normalize(cross(mvPosition.xyz,cr))*rad;
    mvPosition.xy +=  mult*(cr + doublecr).xy;
    cposition = mvPosition.xyz;
    gl_Position = projectionMatrix * mvPosition;
    vLight = vec3(1.0,1.0,1.0);
}

