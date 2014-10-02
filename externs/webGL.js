

/** @typeDef {(string|number)} */
var numlike;


var $3Dmol = {};

/**
 * @constructor
 * @param {...number} color
 */
$3Dmol.Color = function(color) {};
$3Dmol.Color.r;
$3Dmol.Color.g;
$3Dmol.Color.b;

$3Dmol.Math = {};
$3Dmol.Math.clamp = function(x, min, max) {};
$3Dmol.Math.degToRad = function(deg) {};

// Quaternion
/**
 * @constructor
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} w
 */
$3Dmol.Quaternion = function(x, y, z, w) {};
/**
 * 
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} w
 * @return {undefined}
 */
$3Dmol.Quaternion.x;
$3Dmol.Quaternion.y;
$3Dmol.Quaternion.z;
$3Dmol.Quaternion.w;
/** @return {$3Dmol.Quaternion} */
$3Dmol.Quaternion.prototype.set = function(x, y, z, w) {};
/**
 * 
 * @param {$3Dmol.Quaternion} q
 * @return {$3Dmol.Quaternion}
 */
$3Dmol.Quaternion.prototype.copy = function(q) {};
/** @return {$3Dmol.Quaternion} */
$3Dmol.Quaternion.prototype.conjugate = function() {};
/** @return {$3Dmol.Quaternion} */
$3Dmol.Quaternion.prototype.inverse = function() {};
$3Dmol.Quaternion.prototype.length = function() {};
/** @return {$3Dmol.Quaternion} */
$3Dmol.Quaternion.prototype.normalize = function () {};
/**
 * 
 * @param {$3Dmol.Quaternion} g
 * @return {$3Dmol.Quaternion}
 */
$3Dmol.Quaternion.prototype.multiply = function(g) {};
/** 
 * @param {$3Dmol.Quaternion} a
 * @param {$3Dmol.Quaternion} b
 * @return {$3Dmol.Quaternion} 
 */
$3Dmol.Quaternion.prototype.multiplyQuaternions = function(a, b) {};

/** @constructor 
 *  @param {...number} args
 */
$3Dmol.Vector2 = function(args) {};
$3Dmol.Vector2.x;
$3Dmol.Vector2.y;
/** @return {$3Dmol.Vector2} */
$3Dmol.Vector2.prototype.set = function(x, y) {};
/**
 * 
 * @param {$3Dmol.Vector2} a
 * @param {$3Dmol.Vector2} b
 * @return {$3Dmol.Vector2}
 */
$3Dmol.Vector2.protoype.subVectors = function(a, b) {};
/**
 * 
 * @param {$3Dmol.Vector2} v
 * @return {$3Dmol.Vector2}
 */
$3Dmol.Vector2.prototype.copy = function(v) {};
/** @return {$3Dmol.Vector2} */
$3Dmol.Vecto2.prototype.clone = function() {};

//A 3 Vector
/** @constructor 
 *  @param {...number} args
 */
$3Dmol.Vector3 = function(args) {};
$3Dmol.Vector3.x;
$3Dmol.Vector3.y;
$3Dmol.Vector3.z;
/** @return {$3Dmol.Vector3} */
$3Dmol.Vector3.prototype.set = function(x, y, z) {};
/** 
 * @param {$3Dmol.Vector3} v
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.copy = function(v) {};
/**
 * @param {$3Dmol.Vector3} v
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.add = function(v) {};
/**
 * 
 * @param {$3Dmol.Vector3} a
 * @param {$3Dmol.Vector3} b
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.addVectors = function(a, b) {};
/**
 * @param {$3Dmol.Vector3} v
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.sub = function(v) {};
/**
 * 
 * @param {$3Dmol.Vector3} a
 * @param {$3Dmol.Vector3} b
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.subVectors = function(a, b) {};
/**
 * @param {number} s
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.multiplyScalar = function(s) {};
/**
 * @param {number} s
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.divideScalar = function(s) {};
/**
 * @param {$3Dmol.Vector3} v
 * @return {number}
 */
$3Dmol.Vector3.prototype.distanceTo = function(v) {};
/**
 * @param {$3Dmol.Vector3} v
 * @return {number}
 */
$3Dmol.Vector3.prototype.distanceToSquared = function(v) {};
/**
 * @param {$3Dmol.Matrix4} m
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.applyMatrix4 = function(m) {};
/**
 * @param {$3Dmol.Matrix4} m
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.applyProjection = function(m) {};
/**
 * @param {$3Dmol.Quaternion} q
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.applyQuaternion = function(q) {};
/** @return {$3Dmol.Vector3} */
$3Dmol.Vector3.prototype.negate = function() {};
/**
 * @param {$3Dmol.Vector3} v
 * @return {number}
 */
$3Dmol.Vector3.prototype.dot = function(v) {};
$3Dmol.Vector3.prototype.length = function() {};
$3Dmol.Vector3.prototype.lengthSq = function() {};
/** @return {$3Dmol.Vector3} */
$3Dmol.Vector3.prototype.normalize = function() {};
/**
 * @param {$3Dmol.Vector3} v
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.cross = function(v) {};
/**
 * 
 * @param {$3Dmol.Vector3} a
 * @param {$3Dmol.Vector3} b
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.crossVectors = function(a, b) {};
/**
 * @param {$3Dmol.Matrix4} m
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.getPositionFromMatrix = function(m) {};
/**
 * 
 * @param {$3Dmol.Matrix4} m
 * @param {string=} order
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Vector3.prototype.setEulerFromRotationMatrix = function(m, order) {};
/** @return {$3Dmol.Vector3} */
$3Dmol.Vector3.prototype.clone = function() {};

//Matrices

//Matrix3
/** 
 * @constructor 
 * @param {..number} args
 */
$3Dmol.Matrix3 = function(args) {};
$3Dmol.Matrix3.elements;
/** @return {$3Dmol.Matrix3} */
$3Dmol.Matrix3.prototype.set = function(n11, n12, n13, n21, n22, n23, n31, n32, n33) {};
/** @return {$3Dmol.Matrix3} */
$3Dmol.Matrix3.prototype.identity = function() {};
/** 
 * @param {$3Dmol.Matrix3} 
 * @return {$3Dmol.Matrix3} 
 */
$3Dmol.Matrix3.prototype.copy = function(m) {};
/** 
 * @param {number} s
 * @return {$3Dmol.Matrix3} */
$3Dmol.Matrix3.prototype.multiplyScalar = function ( s ) {};
/** 
 * @param {$3Dmol.Matrix3} matrix
 * @param {boolean=} throwOnInvertible
 * @return {$3Dmol.Matrix3} */
$3Dmol.Matrix3.prototype.getInverse = function ( matrix, throwOnInvertible ) {};
/** @return {$3Dmol.Matrix3} */
$3Dmol.Matrix3.prototype.transpose = function () {};
/** @return {$3Dmol.Matrix3} */
$3Dmol.Matrix3.prototype.clone = function () {};

//Matrix 4
/** 
 * @constructor 
 * @param {...number} args
 */
$3Dmol.Matrix4 = function(args) {};
$3Dmol.Matrix4.elements;
/** @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.set = function ( n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44 ) {};
/** @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.identity = function () {};
/** 
 * @param {$3Dmol.Matrix4} m
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.copy = function (m) {};
/** 
 * @param {$3Dmol.Vector3}
 * @param {string=} order
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.setRotationFromEuler = function (v, order) {};
/** 
 * @param {$3Dmol.Quaternion} q
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.setRotationFromQuaternion = function (q) {};
/** 
 * @param {$3Dmol.Vector3} eye
 * @param {$3Dmol.Vector3} target
 * @param {$3Dmol.Vector3} up
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.lookAt = function ( eye, target, up ) {};
/** 
 * @param {$3Dmol.Matrix4} a
 * @param {$3Dmol.Matrix4} b
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.multiplyMatrices = function (a, b) {};
/** 
 * @param {number} s
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.multiplyScalar = function (s) {};
/** @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.transpose = function () {};
/** @return {$3Dmol.Vector3} */
$3Dmol.Matrix4.prototype.getPosition = function() {};
/** 
 * @param {$3Dmol.Vector3} v
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.setPosition = function (v) {};
/** 
 * @param {$3Dmol.Matrix4} m
 * @param {boolean=} throwOnInvertible
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.getInverse = function (m, throwOnInvertible) {};
/** 
 * @param {$3Dmol.Vector3} translation
 * @param {$3Dmol.Quaternion} rotation
 * @param {$3Dmol.Vector3} scale
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.compose = function (translation, rotation, scale) {};
/** 
 * @param {$3Dmol.Vector3=} translation
 * @param {$3Dmol.Quaternion=} rotation
 * @param {$3Dmol.Vector3=} scale
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.decompose = function (translation, rotation, scale) {};
/** 
 * @param {number} v
 * @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.scale = function (v) {};
$3Dmol.Matrix4.prototype.getMaxScaleOnAxis = function() {};
/** @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.makeFrustum = function (left, right, bottom, top, near, far) {};
/** @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.makePerspective = function (fov, aspect, near, far) {};
/** @return {$3Dmol.Matrix4} */
$3Dmol.Matrix4.prototype.clone = function () {};

/**
 * @constructor
 * @param {$3Dmol.Vector3=} origin
 * @param {$3Dmol.Vector3=} direction
 */
$3Dmol.Ray = function(origin, direction) {};
/**
 * @param {$3Dmol.Vector3} origin
 * @param {$3Dmol.Vector3} direction
 * @return {$3Dmol.Ray}
 */
$3Dmol.Ray.prototype.set = function(origin, direction){};
/** 
 * @param {$3Dmol.Ray} ray
 * @return {$3Dmol.Ray} */
$3Dmol.Ray.prototype.copy = function(ray) {};
/**
 * @param {number} t
 * @param {$3Dmol.Vector3=} optionalTarget
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Ray.prototype.at = function(t, optionalTarget) {};
/** @return {$3Dmol.Ray} */
$3Dmol.Ray.prototype.recast = function(t) {};
/**
 * 
 * @param {$3Dmol.Vector3} point
 * @param {$3Dmol.Vector3=} optionalTarget
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Ray.prototype.closestPointToPoint = function(point, optionalTarget) {};
/**
 * 
 * @param {$3Dmol.Vector3} point
 * @return {number}
 */
$3Dmol.Ray.prototype.distanceToPoint = function(point) {};
$3Dmol.Ray.prototype.isIntersectionCylinder = function() {};
/**
 * 
 * @param {$3Dmol.Sphere} sphere
 * @return {boolean}
 */
$3Dmol.Ray.prototype.isIntersectionSphere = function(sphere) {};
$3Dmol.Ray.prototype.isIntersectionPlane = function(plane) {};
$3Dmol.Ray.prototype.distanceToPlane = function(plane) {};
$3Dmol.Ray.prototype.intersectPlane = function(plane, optionalTarget) {};
/** 
 * @param {$3Dmol.Matrix4} matrix4
 * @return {$3Dmol.Ray} */
$3Dmol.Ray.prototype.applyMatrix4 = function(matrix4) {};
/**
 * 
 * @param {$3Dmol.Ray} ray
 * @return {boolean}
 */
$3Dmol.Ray.prototype.equals = function(ray) {};
/** @return {$3Dmol.Ray} */
$3Dmol.Ray.prototype.clone = function() {};

//Intersection sphere for sphere, stick render
/**
 * @constructor
 * @param {$3Dmol.Vector3=} center
 * @param {number=} radius
 */
$3Dmol.Sphere = function(center, radius) {};
/** 
 * @param {$3Dmol.Vector3} center
 * @param {number} radius
 * @return {$3Dmol.Sphere} */
$3Dmol.Sphere.prototype.set = function(center, radius) {};
/** 
 * @param {$3Dmol.Sphere} sphere
 * @return {$3Dmol.Sphere} */
$3Dmol.Sphere.prototype.copy = function(sphere) {};
/** 
 * @param {$3Dmol.Matrix4} matrix
 * @return {$3Dmol.Sphere} */
$3Dmol.Sphere.prototype.applyMatrix4 = function(matrix) {};
/** @return {$3Dmol.Sphere} */
$3Dmol.Sphere.prototype.translate = function(offset) {};
/**
 * @param {$3Dmol.Sphere} sphere
 * @return {boolean}
 */
$3Dmol.Sphere.prototype.equals = function(sphere) {};
/** @return {$3Dmol.Sphere} */
$3Dmol.Sphere.prototype.clone = function() {};

/**
 * @constructor
 * @param {$3Dmol.Vector3=} c1
 * @param {$3Dmol.Vector3=} c2
 * @param {number=} radius
 */
$3Dmol.Cylinder = function(c1, c2, radius) {};
/** 
 * @param {$3Dmol.Cylinder} cylinder
 * @return {$3Dmol.Cylinder} */
$3Dmol.Cylinder.prototype.copy = function(cylinder) {};
$3Dmol.Cylinder.prototype.lengthSq = function() {};
/** 
 * @param {$3Dmol.Matrix4} matrix
 * @return {$3Dmol.Cylinder} */
$3Dmol.Cylinder.prototype.applyMatrix4 = function(matrix) {};

/**
 * @constructor 
 * @param {$3Dmol.Vector3=} a
 * @param {$3Dmol.Vector3=} b
 * @param {$3Dmol.Vector3=} c
 */
$3Dmol.Triangle = function(a, b, c){};
/**
 * @param {$3Dmol.Triangle} triangle
 * @return {$3Dmol.Triangle}
 */
$3Dmol.Triangle.prototype.copy = function(triangle) {};
/**
 * @param {$3Dmol.Matrix4} matrix
 * @return {$3Dmol.Triangle} */
$3Dmol.Triangle.prototype.applyMatrix4 = function(matrix) {};
/** @return {$3Dmol.Vector3} */
$3Dmol.Triangle.prototype.getNormal = function() {};

//Event Handling
$3Dmol.EventDispatcher = function() {};

//Object3D base constructor function
/** @constructor */
$3Dmol.Object3D = function() {};
/** @param {$3Dmol.Vector3} vector */
$3Dmol.Object3D.prototype.lookAt = function(vector) {};
/** @param {$3Dmol.Object3D} object */
$3Dmol.Object3D.prototype.add = function(object) {};
/** @param {$3Dmol.Object3D} object */
$3Dmol.Object3D.prototype.remove = function(object) {};
$3Dmol.Object3D.prototype.updateMatrix = function() {};
/** @param {boolean=} force */
$3Dmol.Object3D.prototype.updateMatrixWorld = function(force) {};
/** 
 * @param {$3Dmol.Object3D} object
 * @return {$3Dmol.Object3D}
 */
$3Dmol.Object3D.prototype.clone = function(object) {};
$3Dmol.Object3D.id;   
$3Dmol.Object3D.name;

$3Dmol.Object3D.parent;
$3Dmol.Object3D.children;
/** @type {$3Dmol.Vector3} */
$3Dmol.Object3D.position;
/** @type {$3Dmol.Vector3} */
$3Dmol.Object3D.rotation;
/** @type {$3Dmol.Matrix4} */
$3Dmol.Object3D.matrix;
/** @type {$3Dmol.Matrix4} */
$3Dmol.Object3D.matrixWorld;
/** @type {$3Dmol.Quaternion} */
$3Dmol.Object3D.quaternion;
/** @type {string} */
$3Dmol.Object3D.eulerOrder;
/** @type {$3Dmol.Vector3} */
$3Dmol.Object3D.up;
/** @type {$3Dmol.Vector3} */
$3Dmol.Object3D.scale;

$3Dmol.Object3D.matrixAutoUpdate;
$3Dmol.Object3D.matrixWorldNeedsUpdate;
$3Dmol.Object3D.rotationAutoUpdate;
$3Dmol.Object3D.useQuaternion;

$3Dmol.Object3D.visible;
$3Dmol.Object3DIDCount;

/**
 * @constructor
 * @param {boolean=} mesh
 */
$3Dmol.Geometry = function(mesh) {};

/**
 * @constructor
 * @param {number=} id
 */
var geometryGroup = function(id) {};
/** @return {$3Dmol.Vector3} */
geometryGroup.prototype.getCentroid = function() {}; 
geometryGroup.prototype.setNormals = function() {};
geometryGroup.prototype.setLineIndices = function() {};
geometryGroup.prototype.truncateArrayBuffers = function() {};
geometryGroup.vertexArray;
geometryGroup.colorArray;
geometryGroup.normalArray;
geometryGroup.faceArray;
geometryGroup.lineArray;
geometryGroup.vertices;
geometryGroup.faceidx;
geometryGroup.lineidx;
geometryGroup.id;

/** @return {geometryGroup} */
$3Dmol.Geometry.prototype.addGeoGroup = function() {};
/** 
 * @param {number=} addVertices
 * @return {geometryGroup} */
$3Dmol.Geometry.prototype.updateGeoGroup = function(addVertices) {};
/**
 * @param {boolean} three
 */
$3Dmol.Geometry.prototype.setUpNormals = function(three) {};

/** @return {undefined} */
$3Dmol.Geometry.prototype.setUpWireframe = function() {};
/** @return {undefined} */
$3Dmol.Geometry.prototype.initTypedArrays = function() {};
/** @return {undefined} */
$3Dmol.Geometry.prototype.dispose = function() {};
/** @type {number} */
$3Dmol.Geometry.prototype.vertices;

$3Dmol.Geometry.id;

$3Dmol.Geometry.name;

$3Dmol.Geometry.hasTangents;

$3Dmol.Geometry.dynamic;
$3Dmol.Geometry.mesh;
$3Dmol.Geometry.verticesNeedUpdate;
$3Dmol.Geometry.elementsNeedUpdate;
$3Dmol.Geometry.normalsNeedUpdate;
$3Dmol.Geometry.colorsNeedUpdate;
$3Dmol.Geometry.buffersNeedUpdate;
/** @type {Array.<geometryGroup>} */
$3Dmol.Geometry.geometryGroups;
$3Dmol.Geometry.groups;

$3Dmol.GeometryIDCount;

/**
 * @constructor
 * @param {$3Dmol.Vector3} origin
 * @param {$3Dmol.Vector3} direction
 * @param {number=} far
 * @param {number=} near
 */
$3Dmol.Raycaster = function(origin, direction, far, near) {};
//Raycaster
$3Dmol.Raycaster.prototype.precision;
$3Dmol.Raycaster.prototype.linePrecision;
/**
 * @param {$3Dmol.Vector3} origin
 * @param {$3Dmol.Vector3} direction
 */
$3Dmol.Raycaster.prototype.set = function(origin, direction) {};
/**
 * @param {$3Dmol.Object3D} group
 * @param {Array.<Object>} objects
 * @return {Array.<Object>} 
 */
$3Dmol.Raycaster.prototype.intersectObjects = function(group, objects) {};

/** @struct */
var IntersectionShapes = {};
/** @type {Array.<$3Dmol.Sphere>} */
IntersectionShapes.sphere;
/** @type {Array.<$3Dmol.Vector3>} */
IntersectionShapes.line;
/** @type {Array.<$3Dmol.Triangle> */
IntersectionShapes.triangle;
/** @type {Array.<$3Dmol.Cylinder> */
IntersectionShapes.cylinder;

/** @constructor */
$3Dmol.Projector = function () {};
/**
 * @param {$3Dmol.Vector3} vector
 * @param {$3Dmol.Camera} camera
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Projector.projectVector = function(vector, camera) {};
/**
 * 
 * @param {$3Dmol.Vector3} vector
 * @param {$3Dmol.Camera} camera
 * @return {$3Dmol.Vector3}
 */
$3Dmol.Projector.unprojectVector = function(vector, camera) {};

/**
 * @constructor
 * @extends $3Dmol.Object3D
 * @param {number=} fov
 * @param {number=} aspect
 * @param {number=} near
 * @param {number=} far
 * 
 */
$3Dmol.Camera = function(fov, aspect, near, far) {};
    
/** @override */
$3Dmol.Camera.prototype.lookAt = function(vector){};

$3Dmol.Camera.prototype.updateProjectionMatrix = function () {};

/** @type {number} */
$3Dmol.Camera.fov;
/** @type {number} */
$3Dmol.Camera.aspect;
/** @type {number} */
$3Dmol.Camera.near;
/** @type {number} */
$3Dmol.Camera.far;

/** @type {$3Dmol.Matrix4} */
$3Dmol.Camera.projectionMatrix;
/** @type {$3Dmol.Matrix4} */
$3Dmol.Camera.projectionMatrixInverse;
/** @type {$3Dmol.Matrix4} */
$3Dmol.Camera.matrixWorldInverse;

$3Dmol.SpritePlugin = function () {};

/** 
 * @constructor
 * @extends {$3Dmol.Object3D}
 * @param {number} hex
 * @param {number=} intensity
 */
$3Dmol.Light = function(hex, intensity) {};
/** @type {$3Dmol.Color} */
$3Dmol.Light.color;
/** @type {$3Dmol.Vector3} */
$3Dmol.Light.position;
/** @type {$3Dmol.Object3D} */
$3Dmol.Light.target;
/** @type {number} */
$3Dmol.Light.intensity;
/** @type {boolean} */
$3Dmol.Light.castShadow;
/** @type {boolean} */
$3Dmol.Light.onlyShadow;

/** 
 * @constructor 
 *  @extends $3Dmol.EventDispatcher
 */
$3Dmol.Material = function () {};

$3Dmol.Material.id;
$3Dmol.Material.name;
$3Dmol.Material.side;
$3Dmol.Material.opacity;
$3Dmol.Material.transparent;
$3Dmol.Material.blending;
$3Dmol.Material.depthTest;
$3Dmol.Material.depthWrite;
$3Dmol.Material.polygonOffset;
$3Dmol.Material.polygonOffsetFactor;
$3Dmol.Material.polygonOffsetUnits;
$3Dmol.Material.alphaTest;
$3Dmol.Material.visible;
$3Dmol.Material.needsUpdate;

/**
 * @param {matSpec} values
 */
$3Dmol.Material.prototype.setValues = function(values) {};
/**
 * @param {$3Dmol.Material=} material
 * @return {$3Dmol.Material}
 */
$3Dmol.Material.prototype.clone = function(material) {};

$3Dmol.Material.prototype.dispose = function () {};

/** 
 * Since we can instantiate materials with a matSpec
 * @constructor
 * @struct
 * @extends {$3Dmol.Material}
 */
var matSpec = {};
matSpec.color;
/** @type {numlike} */
matSpec.linewidth;
matSpec.linecap;
matSpec.linejoin;
matSpec.vertexColors;
matSpec.fog;
matSpec.ambient;
matSpec.emissive;
matSpec.wrapAround;
matSpec.wrapRGB;
matSpec.map;
matSpec.lightMap;
matSpec.specularMap;
matSpec.envMap;
matSpec.reflectivity;
matSpec.refractionRatio;
matSpec.wireframe;
matSpec.wireframeLinewidth;
matSpec.wireframeLinecap;
matSpec.wireframeLinejoin;
matSpec.shading;
matSpec.skinning;
matSpec.useScreenCoordinates;
matSpec.depthTest;
matSpec.sizeAttenuation;
matSpec.scaleByViewPort;
matSpec.alignment;
matSpec.uvOffset;
matSpec.uvScale;

$3Dmol.MaterialIdCount;

//material constants
// sides

$3Dmol.FrontSide;
$3Dmol.BackSide;
$3Dmol.DoubleSide;

// blending modes
$3Dmol.NoBlending;
$3Dmol.NormalBlending;
$3Dmol.AdditiveBlending;
$3Dmol.SubtractiveBlending;
$3Dmol.MultiplyBlending;
$3Dmol.CustomBlending;

// shading
$3Dmol.NoShading;
$3Dmol.FlatShading;
$3Dmol.SmoothShading;

// colors
$3Dmol.NoColors;
$3Dmol.FaceColors;
$3Dmol.VertexColors;

//Texture constants
//TODO: Which of these do I need (since I only use textures to display label sprites) ?
$3Dmol.MultiplyOperation;
$3Dmol.MixOperation;
$3Dmol.AddOperation;

// mapping modes

$3Dmol.UVMapping = function() {};

// wrapping modes
$3Dmol.ClampToEdgeWrapping;

//Filters
$3Dmol.LinearFilter;
$3Dmol.LinearMipMapLinearFilter;

//Data types
$3Dmol.UnsignedByteType;

//Pixel formats
$3Dmol.RGBAFormat;

/**
 * @constructor
 * @extends {$3Dmol.Material}
 * @param {matSpec} parameters
 */
$3Dmol.LineBasicMaterial = function(parameters) {};
/**
 * @override
 */
$3Dmol.LineBasicMaterial.prototype.clone = function() {};

/**
 * @construtor
 * @extends {$3Dmol.Material}
 * @param {matSpec} parameters
 */
$3Dmol.MeshLambertMaterial = function(parameters) {};
/** @override */
$3Dmol.MeshLambertMaterial.prototype.clone = function() {};

/**
 * @constructor
 * @struct
 * @extends {$3Dmol.Material}
 * @param {matSpec} parameters
 */
$3Dmol.SpriteMaterial = function(parameters) {};
/** @override */
$3Dmol.SpriteMaterial.prototype.clone = function() {};

//Alignment for Sprites
/** @enum {$3Dmol.Vector2} */
$3Dmol.SpriteAlignment = {};
$3Dmol.SpriteAlignment.topLeft;
$3Dmol.SpriteAlignment.topCenter;
$3Dmol.SpriteAlignment.topRight;
$3Dmol.SpriteAlignment.centerLeft;
$3Dmol.SpriteAlignment.center;
$3Dmol.SpriteAlignment.centerRight;
$3Dmol.SpriteAlignment.bottomLeft;
$3Dmol.SpriteAlignment.bottomCenter;
$3Dmol.SpriteAlignment.bottomRight;

/**
 * @constructor
 * @extends {$3Dmol.EventDispatcher}
 * @param {Object} image - An html canvas element
 */
$3Dmol.Texture = function(image) {};
$3Dmol.Texture.needsUpdate;
/**
 * @param {$3Dmol.Texture=} texture
 * @return {$3Dmol.Texture}
 */
$3Dmol.Texture.prototype.clone = function(texture) {};
$3Dmol.Texture.prototype.dispose = function() {};

$3Dmol.TextureIdCount;

/**
 * @constructor
 * @extends {$3Dmol.Object3D}
 * @param {$3Dmol.Geometry} geometry
 * @param {$3Dmol.Material=} material
 * @param {numlike} type
 */
$3Dmol.Line = function (geometry, material, type) {};

$3Dmol.LineStrip;
$3Dmol.LinePieces;

/** 
 * @param {$3Dmol.Line)
 * @return {$3Dmol.Line}
 */
$3Dmol.Line.prototype.clone = function(object) {};

/**
 * @constructor
 * @extends {$3Dmol.Object3D}
 * @param {$3Dmol.Geometry} geometry
 * @param {$3Dmol.Material} material
 */
$3Dmol.Mesh = function(geometry, material) {};

/**
 * @param {$3Dmol.Mesh} object
 * @return {$3Dmol.Mesh}
 */
$3Dmol.Mesh.prototype.clone = function (object) {};

/** @type {$3Dmol.Geometry} */
$3Dmol.Mesh.geometry;
/** @type {$3Dmol.Material} */
$3Dmol.Mesh.material;

/**
 * 
 * @constructor
 * @extends {$3Dmol.Object3D}
 * @param {$3Dmol.Material}
 */
$3Dmol.Sprite = function(material) {};
$3Dmol.Sprite.rotation3d;
$3Dmol.Sprite.rotation;
/** @override */
$3Dmol.Sprite.prototype.updateMatrix = function() {};
/** 
 * 
 * @param {type} object
 * @return {$3Dmol.Sprite|$3Dmol.Sprite.prototype.clone.object|$3Dmol.Object3D|_L16.my.Sprite.prototype.clone.object|$3Dmol.Object3D.prototype.clone.object}@param {$3Dmol.Sprite}
 * @return {$3Dmol.Sprite}
 */
$3Dmol.Sprite.prototype.clone = function(object) {};

/** 
 * @constructor
 * @struct
 */
var rendererSpec = {};
rendererSpec.canvas;
rendererSpec.precision;
rendererSpec.alpha;
rendererSpec.premultipliedAlpha;
rendererSpec.antialias;
rendererSpec.stencil;
rendererSpec.preserveDrawingBuffer;
rendererSpec.clearColor;
rendererSpec.clearAlpha;
rendererSpec.devicePixelRatio;

/**
 * @constructor
 * @param {rendererSpec} parameters
 */
$3Dmol.Renderer = function(parameters) {};
$3Dmol.Renderer.domElement;
/**
 * @param {number} width
 * @param {number} height
 */
$3Dmol.Renderer.setSize = function(width, height) {};
/**
 * @param {number} hex
 * @param {number} alpha
 */
$3Dmol.Renderer.setClearColorHex = function(hex, alpha) {};
/**
 * 
 * @param {$3Dmol.Scene} scene
 * @param {$3Dmol.Camera} camera
 * @param {boolean=} forceClear
 */
$3Dmol.Renderer.render = function(scene, camera, forceClear) {};

/**
 * @constructor
 * @extends {$3Dmol.Object3D}
 */
$3Dmol.Scene = function() {};
/** @type {$3Dmol.Fog} */
$3Dmol.Scene.fog;    
$3Dmol.Scene.overrideMaterial;
$3Dmol.Scene.matrixAutoUpdate;
/** @type {Array.<$3Dmol.Object3D>} */
$3Dmol.Scene.__objects;
/** @type {Array.<$3Dmol.Light>} */
$3Dmol.Scene.__lights;
/** @type {Array.<$3Dmol.Object3D>} */
$3Dmol.Scene.__objectsAdded;
/** @type {Array.<$3Dmol.Object3D>} */
$3Dmol.Scene.__objectsRemoved;

/**
 * @param {$3Dmol.Object3D} object
 */
$3Dmol.Scene.prototype.__addObject = function(object) {};

/**
 * @param {$3Dmol.Object3D} object
 */
$3Dmol.Scene.prototype.__removeObject = function(object) {};

/**
 * @constructor
 * @param {number} hex
 * @param {number=} near
 * @param {number=} far
 */
$3Dmol.Fog = function(hex, near, far) {};
$3Dmol.Fog.name;
/** @type {$3Dmol.Color} */
$3Dmol.Fog.color;
$3Dmol.Fog.near;
$3Dmol.Fog.far;
/** @return {$3Dmol.Fog} */
$3Dmol.Fog.prototype.clone = function() {};

/**
 * @struct
 */
$3Dmol.ShaderUtils = {};

/**
 * @param {uniformsList} uniforms_src
 * @return {uniformsList}
 */    
$3Dmol.ShaderUtils.clone = function(uniforms_src) {};

/**
 * @dict
 */
$3Dmol.ShaderLib = {};
$3Dmol.ShaderLib['basic'];
$3Dmol.ShaderLib['lambert'];
$3Dmol.ShaderLib['sprite'];

/** @struct */
var shader = {};
/** @type {string} */
shader.fragmentShader;
/** @type {string} **/
shader.vertexShader;
/** @type {uniformsList} */
shader.uniforms;

/** @typedef {{type: string, value}} */
var uniform;

/** @type {Object.<?,{uniform}>} */
var uniformsList = {};
uniformsList.opacity;
uniformsList.diffuse;
uniformsList.fogColor;
uniformsList.fogNear;
uniformsList.fogFar;
uniformsList.ambient;           
uniformsList.emissive;
uniformsList.ambientLightColor;
uniformsList.directionalLightColor;
uniformsList.directionalLightDirection;
