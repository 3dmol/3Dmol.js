# Math Routines for 3Dmol.js

This is the documentation for the essential math involved in 3Dmol.js. 
Link to the JavaScript file : [Link](https://github.com/3dmol/3Dmol.js/blob/master/3Dmol/WebGL/math.js)
## Classes:
    
 - [Quarternion](#Quarternion-Class)
 - [Vector2](#Vector2-Class)
 - [Vector3](#Vector3-Class)
 - [Matrix3](#Matrix3-Class)
 - [Matrix4](#Matrix4-Class)
 - [Ray](#Ray-Class)
 

### Quarternion Class

#### Short Descriptiton
The quarternion is a set of 4 numbers used to represent rotations. The set is [z,y,z,w] and they are:-
```
x = RotationAxis.x * sin(RotationAngle / 2)
y = RotationAxis.y * sin(RotationAngle / 2)
z = RotationAxis.z * sin(RotationAngle / 2)
w = cos(RotationAngle / 2)
```
where 
- RotationAxis is the axis around which we want to make our rotation. 
- RotationAngle is the angle of rotation to be performed on this axis.

Also see why Quarternions are preferred over Euler Angles.

#### Utility Functions



