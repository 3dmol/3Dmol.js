// @ts-check

import { Matrix4 } from "./Matrix4";

// module level declaration of helper matrix to avoid reallocation
var matrix = new Matrix4();

export class Matrix3 {
    constructor(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
        this.elements = new Float32Array(9);

        this.set(
            n11 !== undefined ? n11 : 1,
            n12 || 0,
            n13 || 0,
            n21 || 0,
            n22 !== undefined ? n22 : 1,
            n23 || 0,
            n31 || 0,
            n32 || 0,
            n33 !== undefined ? n33 : 1
        );
    }

    set(n11, n12, n13, n21, n22, n23, n31, n32, n33) {
        var te = this.elements;

        te[0] = n11;
        te[3] = n12;
        te[6] = n13;
        te[1] = n21;
        te[4] = n22;
        te[7] = n23;
        te[2] = n31;
        te[5] = n32;
        te[8] = n33;

        return this;
    }

    identity() {
        this.set(1, 0, 0, 0, 1, 0, 0, 0, 1);

        return this;
    }

    copy(m) {
        var me = m.elements;

        this.set(me[0], me[3], me[6], me[1], me[4], me[7], me[2], me[5], me[8]);
    }

    multiplyScalar(s) {
        var te = this.elements;

        te[0] *= s;
        te[3] *= s;
        te[6] *= s;
        te[1] *= s;
        te[4] *= s;
        te[7] *= s;
        te[2] *= s;
        te[5] *= s;
        te[8] *= s;

        return this;
    }

    getInverse3(matrix) {
        // input: Matrix3

        let me = matrix.elements;
        let te = this.elements;

        te[0] = me[4] * me[8] - me[5] * me[7];
        te[3] = me[6] * me[5] - me[3] * me[8];
        te[6] = me[3] * me[7] - me[6] * me[4];
        te[1] = me[7] * me[2] - me[1] * me[8];
        te[4] = me[0] * me[8] - me[6] * me[2];
        te[7] = me[1] * me[6] - me[0] * me[7];
        te[2] = me[1] * me[5] - me[2] * me[4];
        te[5] = me[2] * me[3] - me[0] * me[5];
        te[8] = me[0] * me[4] - me[1] * me[3];

        let det = me[0] * te[0] + me[3] * te[1] + me[6] * te[2];
        this.multiplyScalar(1.0 / det);

        return this;
    }

    getInverse(matrix, throwOnInvertible) {
        // input: Matrix4

        var me = matrix.elements;
        var te = this.elements;

        te[0] = me[10] * me[5] - me[6] * me[9];
        te[1] = -me[10] * me[1] + me[2] * me[9];
        te[2] = me[6] * me[1] - me[2] * me[5];
        te[3] = -me[10] * me[4] + me[6] * me[8];
        te[4] = me[10] * me[0] - me[2] * me[8];
        te[5] = -me[6] * me[0] + me[2] * me[4];
        te[6] = me[9] * me[4] - me[5] * me[8];
        te[7] = -me[9] * me[0] + me[1] * me[8];
        te[8] = me[5] * me[0] - me[1] * me[4];

        var det = me[0] * te[0] + me[1] * te[3] + me[2] * te[6];

        // no inverse

        if (det === 0) {
            var msg = "Matrix3.getInverse(): can't invert matrix, determinant is 0";

            if (throwOnInvertible || false) {
                throw new Error(msg);
            } else {
                console.warn(msg);
            }

            this.identity();

            return this;
        }

        this.multiplyScalar(1.0 / det);

        return this;
    }

    getDeterminant() {
        var m = this.elements;

        /*
         * |a b c| |d e f| |g h i|
         */

        var determinant =
            m[0] * m[4] * m[8] + // +aei
            m[1] * m[5] * m[6] + // +bfg
            m[2] * m[3] * m[7] - // +cdh
            m[2] * m[4] * m[6] - // -ceg
            m[1] * m[3] * m[8] - // -bdi
            m[0] * m[5] * m[7]; // -afh
        return determinant;
    }

    getMatrix4() {
        var m = this.elements;
        return new Matrix4(
            m[0],
            m[3],
            m[6],
            0,
            m[1],
            m[4],
            m[7],
            0,
            m[2],
            m[5],
            m[8],
            0
        );
    }

    transpose() {
        var tmp,
            m = this.elements;

        tmp = m[1];
        m[1] = m[3];
        m[3] = tmp;
        tmp = m[2];
        m[2] = m[6];
        m[6] = tmp;
        tmp = m[5];
        m[5] = m[7];
        m[7] = tmp;

        return this;
    }

    clone() {
        var te = this.elements;

        return new Matrix3(
            te[0],
            te[3],
            te[6],
            te[1],
            te[4],
            te[7],
            te[2],
            te[5],
            te[8]
        );
    }
    unproject(camera) {
        matrix.multiplyMatrices(
            camera.matrixWorld,
            matrix.getInverse(camera.projectionMatrix)
        );
        return this.applyMatrix4(matrix);
    }
}
