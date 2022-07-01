

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

varying vec3 vColor;
varying vec3 vLight;
varying vec3 cposition;
varying vec3 p1;
varying vec3 p2;
varying float r;

void main() {

    vColor = color; vColor.z = abs(vColor.z); //z indicates which vertex and so would vary
    r = abs(radius);
    vec4 to = modelViewMatrix*vec4(normal, 1.0); //normal is other point of cylinder
    vec4 pt = modelViewMatrix*vec4(position, 1.0);
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
       if(projectionMatrix[3][3] == 0.0) { //perspective
         vec3 pnorm = normalize(p1);
         float t = dot(mvPosition.xyz-p1,n)/dot(pnorm,n);
         mvPosition.xyz = p1+t*pnorm; 
       } else { //orthographic
         mvPosition.xyz = p1;
       }
    } else {
      if(projectionMatrix[3][3] == 0.0) { //perspective
         vec3 pnorm = normalize(p2);
         float t = dot(mvPosition.xyz-p2,n)/dot(pnorm,n);
         mvPosition.xyz = p2+t*pnorm;
       } else { //orthographic
         mvPosition.xyz = p2;
       } 
       mult *= -1.0;
    }
    vec3 cr = normalize(cross(mvPosition.xyz,norm))*radius;
    vec3 doublecr = normalize(cross(mvPosition.xyz,cr))*radius;
    mvPosition.xyz +=  mult*(cr + doublecr).xyz;
    cposition = mvPosition.xyz;
    gl_Position = projectionMatrix * mvPosition;
    vec4 lDirection = viewMatrix * vec4( directionalLightDirection[ 0 ], 0.0 );
    vLight = normalize( lDirection.xyz )*directionalLightColor[0]; //not really sure this is right, but color is always white so..
}

