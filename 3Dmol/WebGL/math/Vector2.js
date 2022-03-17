// A 2 Vector
/** @constructor */
export class Vector2 {
  constructor(x, y) {
    this.x = x || 0.0;
    this.y = y || 0.0;
  }

  set(x, y) {
    this.x = x;
    this.y = y;

    return this;
  }

  subVectors(a, b) {
    this.x = a.x - b.x;
    this.y = a.y - b.y;

    return this;
  }

  copy(v) {
    this.x = v.x;
    this.y = v.y;

    return this;
  }

  clone() {
    return new Vector2(this.x, this.y);
  }
}
