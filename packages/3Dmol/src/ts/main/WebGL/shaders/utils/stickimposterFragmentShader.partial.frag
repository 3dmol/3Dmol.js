uniform float opacity;
uniform mat4 projectionMatrix;

uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;

varying vec3 vLight;
varying vec3 vColor;
varying vec3 cposition;
varying vec3 p1;
varying vec3 p2;
varying float r;

//DEFINEFRAGCOLOR

//cylinder-ray intersection testing taken from http://mrl.nyu.edu/~dzorin/cg05/lecture12.pdf
//also useful: http://stackoverflow.com/questions/9595300/cylinder-impostor-in-glsl
//with a bit more care (caps) this could be a general cylinder imposter (see also outline)
void main() {
    vec3 color = abs(vColor);
    vec3 pos = cposition;
    vec3 p = pos; //ray point
    vec3 v = vec3(0.0,0.0,-1.0); //ray normal - orthographic
    if(projectionMatrix[3][3] == 0.0) v = normalize(pos); //ray normal - perspective
    vec3 pa = p1; //cyl start
    vec3 va = normalize(p2-p1); //cyl norm
    vec3 tmp1 = v-(dot(v,va)*va);
    vec3 deltap = p-pa;
    float A = dot(tmp1,tmp1);
    if(A == 0.0) discard;
    vec3 tmp2 = deltap-(dot(deltap,va)*va);
    float B = 2.0*dot(tmp1, tmp2);
    float C = dot(tmp2,tmp2)-r*r;
//quadratic equation!
    float det = (B*B) - (4.0*A*C);
    if(det < 0.0) discard;
    float sqrtDet = sqrt(det);
    float posT = (-B+sqrtDet)/(2.0*A);
    float negT = (-B-sqrtDet)/(2.0*A);
    float intersectionT = min(posT,negT);
    vec3 qi = p+v*intersectionT;
    float dotp1 = dot(va,qi-p1);
    float dotp2 = dot(va,qi-p2);
    vec3 norm;
    if( dotp1 < 0.0 || dotp2 > 0.0) { //(p-c)^2 + 2(p-c)vt +v^2+t^2 - r^2 = 0
       vec3 cp;
       if( dotp1 < 0.0) {  
//        if(vColor.x < 0.0 ) discard; //color sign bit indicates if we should cap or not
        cp = p1;
       } else {
//          if(vColor.y < 0.0 ) discard;
          cp = p2;
       }
       vec3 diff = p-cp;
       A = dot(v,v);
       B = dot(diff,v)*2.0;
       C = dot(diff,diff)-r*r;
       det = (B*B) - (4.0*C);
       if(det < 0.0) discard;
       sqrtDet = sqrt(det);
       posT = (-B+sqrtDet)/(2.0);
       negT = (-B-sqrtDet)/(2.0);
       float t = min(posT,negT);
       qi = p+v*t; 
       norm = normalize(qi-cp); 
    } else {
       norm = normalize(qi-(dotp1*va + p1));
    }
    vec4 clipPos = projectionMatrix * vec4(qi, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    float depth = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;
    gl_FragDepthEXT = depth;