# Math Routines for 3Dmol.js

This is the documentation for the essential math involved in 3Dmol.js.
Link to the JavaScript file : [Link](https://github.com/3dmol/3Dmol.js/blob/master/3Dmol/WebGL/math.js)

## Classes

- [Quarternion](#Quarternion-Class)
- [Vector2](#Vector2-Class)
- [Vector3](#Vector3-Class)
- [Matrix3](#Matrix3-Class)
- [Matrix4](#Matrix4-Class)
- [Ray](#Ray-Class)

### Quarternion Class

#### Short Descriptiton

The quarternion is a set of 4 numbers used to represent rotations. The set is [z,y,z,w] and they are:-

```javascript
x = RotationAxis.x * sin(RotationAngle / 2)
y = RotationAxis.y * sin(RotationAngle / 2)
z = RotationAxis.z * sin(RotationAngle / 2)
w = cos(RotationAngle / 2)
```

where

- RotationAxis is the axis around which we want to make our rotation. 
- RotationAngle is the angle of rotation to be performed on this axis.

Also see why Quarternions are preferred over Euler Angles.

#### Constructor for Quarternion Class

**Quarternion(x,y,z,w)**

Parameters:

| Name 	|  Type  	|                 Description                 	|
|:----:	|:------:	|:-------------------------------------------:	|
|   x  	| Number 	| Multiplier of the imaginary basis vector i. 	|
|   y  	| Number 	| Multiplier of the imaginary basis vector j. 	|
|   z  	| Number 	| Multiplier of the imaginary basis vector k. 	|
|   w  	| Number 	|         Multiplier of the real part.        	|

#### Methods

##### set(x,y,z,w)

Sets the parameters x, y, z, and w.

Parameters

| Name 	|  Type  	|                 Description                 	|
|:----:	|:------:	|:-------------------------------------------:	|
|   x  	| Number 	| Multiplier of the imaginary basis vector i. 	|
|   y  	| Number 	| Multiplier of the imaginary basis vector j. 	|
|   z  	| Number 	| Multiplier of the imaginary basis vector k. 	|
|   w  	| Number 	|         Multiplier of the real part.        	|

Returns:

Type

*quarternion*

##### copy(q)

Copies one quarternion to another.

Parameters

| Name  	| Type        	| Description                                        	|
|-------	|-------------	|----------------------------------------------------	|
| q     	| Quarternion 	| The quarternion whose parameters have to be copied 	|

Returns:
 
Type

Quarternion

##### conjugate()

Returns the conjugate of the current quarternion.

Parameters: none

Returns:
 
Type

Quarternion

##### inverse()

Returns the inverse quarternion rotation.

Parameters: none

Returns:
 
Type

Quarternion

##### length()

Returns length from Origin (includes the rotation angle, w)

Parameters: none

Returns:
 
Type

Number

##### lengthxyz()

Returns length from Origin (does not includes the rotation angle, w)

Parameters: none

Returns:
 
Type

Number

##### normalize()

Normalize the quaternion. 
Note : Values of the quarternion are changed.
Returns:
 
Type

Quarternion

##### multiply(q)

Multiplies current quarternion with q.

Parameters

| Name  	| Type        	| Description                                        	|
|-------	|-------------	|----------------------------------------------------	|
| q     	| Quarternion 	| The quarternion whose parameters have to be multiplied	|

Returns:
 
Type

Quarternion

##### multiplyScalar(s)

Multiplies every parameter of the quarternion with the scalar s.

Parameters

| Name  	| Type        	| Description                                        	|
|-------	|-------------	|----------------------------------------------------	|
| s     	| Number 	| The number with which the parameters should be multiplied	|

Returns:
 
Type

Quarternion

##### sub(q)

Subtracts the parameters of q from the current quarternion.

Parameters

| Name  	| Type        	| Description                                        	|
|-------	|-------------	|----------------------------------------------------	|
| q     	| Quarternion 	| The quarternion whose parameters have to be subtracted 	|

Returns:
 
Type

Quarternion

##### clone()

Returns the clone of the current quarternion.

Parameters

None

Returns:
 
Type

Quarternion

##### setfromEuler(e)

Converts and sets the values of the parameters of the current quarternion from Euler Angles.


Parameters

| Name  	| Type        	| Description                                        	|
|-------	|-------------	|----------------------------------------------------	|
| e     	| Euler Angle Vector	| The Euler Angle Vector from whom the values are converted. 	|

Returns:
 
Type

Quarternion









