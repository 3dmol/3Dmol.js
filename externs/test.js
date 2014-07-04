var WebMol = {};
WebMol.Math = {};
WebMol.Math.clamp = function(x, min, max) {};
WebMol.Math.degToRad = function(deg) {};

// Quaternion
/**
 * @constructor
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} w
 */
WebMol.Quaternion = function(x, y, z, w) {};
/**
 * 
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} w
 * @returns {undefined}
 */
WebMol.Quaternion.x;
WebMol.Quaternion.y;
WebMol.Quaternion.z;
WebMol.Quaternion.w;
/** @return {WebMol.Quaternion} */
WebMol.Quaternion.prototype.set = function(x, y, z, w) {};
/**
 * 
 * @param {WebMol.Quaternion} q
 * @returns {WebMol.Quaternion}
 */
WebMol.Quaternion.prototype.copy = function(q) {};
/** @return {WebMol.Quaternion} */
WebMol.Quaternion.prototype.conjugate = function() {};
/** @return {WebMol.Quaternion} */
WebMol.Quaternion.prototype.inverse = function() {};
WebMol.Quaternion.prototype.length = function() {};
/** @return {WebMol.Quaternion} */
WebMol.Quaternion.prototype.normalize = function () {};
/**
 * 
 * @param {WebMol.Quaternion} g
 * @returns {WebMol.Quaternion}
 */
WebMol.Quaternion.prototype.multiply = function(g) {};
/** 
 * @param {WebMol.Quaternion} a
 * @param {WebMol.Quaternion} b
 * @return {WebMol.Quaternion} 
 */
WebMol.Quaternion.prototype.multiplyQuaternions = function(a, b) {};

//Matrix 4
/** 
 * @constructor 
 * @param {...number} args
 */
WebMol.Matrix4 = function(args) {};
WebMol.Matrix4.elements;
/** @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.set = function ( n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44 ) {};
/** @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.identity = function () {};
/** 
 * @param {WebMol.Matrix4} m
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.copy = function (m) {};
/** 
 * @param {WebMol.Vector3}
 * @param {string=} order
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.setRotationFromEuler = function (v, order) {};
/** 
 * @param {WebMol.Quaternion} q
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.setRotationFromQuaternion = function (q) {};
/** 
 * @param {WebMol.Vector3} eye
 * @param {WebMol.Vector3} target
 * @param {WebMol.Vector3} up
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.lookAt = function ( eye, target, up ) {};
/** 
 * @param {WebMol.Matrix4} a
 * @param {WebMol.Matrix4} b
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.multiplyMatrices = function (a, b) {};
/** 
 * @param {number} s
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.multiplyScalar = function (s) {};
/** @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.transpose = function () {};
/** @return {WebMol.Vector3} */
WebMol.Matrix4.prototype.getPosition = function() {};
/** 
 * @param {WebMol.Vector3} v
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.setPosition = function (v) {};
/** 
 * @param {WebMol.Matrix4} m
 * @param {boolean=} throwOnInvertible
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.getInverse = function (m, throwOnInvertible) {};
/** 
 * @param {WebMol.Vector3} translation
 * @param {WebMol.Quaternion} rotation
 * @param {WebMol.Vector3} scale
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.compose = function (translation, rotation, scale) {};
/** 
 * @param {WebMol.Vector3=} translation
 * @param {WebMol.Quaternion=} rotation
 * @param {WebMol.Vector3=} scale
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.decompose = function (translation, rotation, scale) {};
/** 
 * @param {number} v
 * @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.scale = function (v) {};
WebMol.Matrix4.prototype.getMaxScaleOnAxis = function() {};
/** @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.makeFrustum = function (left, right, bottom, top, near, far) {};
/** @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.makePerspective = function (fov, aspect, near, far) {};
/** @return {WebMol.Matrix4} */
WebMol.Matrix4.prototype.clone = function () {};

/** @constructor */
WebMol.Vector3 = function(x, y, z) {};

WebMol.Vector3.x;
WebMol.Vector3.y;
WebMol.Vector3.z;
/** @return {WebMol.Vector3} */
WebMol.Vector3.prototype.set = function(x, y, z) {};
/** 
 * @param {WebMol.Vector3} v
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.copy = function(v) {};
/**
 * @param {WebMol.Vector3} v
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.add = function(v) {};
/**
 * 
 * @param {WebMol.Vector3} a
 * @param {WebMol.Vector3} b
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.addVectors = function(a, b) {};
/**
 * @param {WebMol.Vector3} v
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.sub = function(v) {};
/**
 * 
 * @param {WebMol.Vector3} a
 * @param {WebMol.Vector3} b
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.subVectors = function(a, b) {};
/**
 * @param {number} s
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.multiplyScalar = function(s) {};
/**
 * @param {number} s
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.divideScalar = function(s) {};
/**
 * @param {WebMol.Vector3} v
 * @return {number}
 */
WebMol.Vector3.prototype.distanceTo = function(v) {};
/**
 * @param {WebMol.Vector3} v
 * @return {number}
 */
WebMol.Vector3.prototype.distanceToSquared = function(v) {};
/**
 * @param {WebMol.Matrix4} m
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.applyMatrix4 = function(m) {};
/**
 * @param {WebMol.Matrix4} m
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.applyProjection = function(m) {};
/**
 * @param {WebMol.Quaternion} q
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.applyQuaternion = function(q) {};
/** @return {WebMol.Vector3} */
WebMol.Vector3.prototype.negate = function() {};
/**
 * @param {WebMol.Vector3} v
 * @return {number}
 */
WebMol.Vector3.prototype.dot = function(v) {};
WebMol.Vector3.prototype.length = function() {};
WebMol.Vector3.prototype.lengthSq = function() {};
/** @return {WebMol.Vector3} */
WebMol.Vector3.prototype.normalize = function() {};
/**
 * @param {WebMol.Vector3} v
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.cross = function(v) {};
/**
 * 
 * @param {WebMol.Vector3} a
 * @param {WebMol.Vector3} b
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.crossVectors = function(a, b) {};
/**
 * @param {WebMol.Matrix4} m
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.getPositionFromMatrix = function(m) {};
/**
 * 
 * @param {WebMol.Matrix4} m
 * @param {string=} order
 * @return {WebMol.Vector3}
 */
WebMol.Vector3.prototype.setEulerFromRotationMatrix = function(m, order) {};
/** @return {WebMol.Vector3} */
WebMol.Vector3.prototype.clone = function() {};