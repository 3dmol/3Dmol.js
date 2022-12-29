/** @class 
 *  @subcategory  Math
 * */ 
export class Vector2 {
    x: number;
    y: number;

    constructor(x: number, y: number) {
      this.x = x || 0.0;
      this.y = y || 0.0;
    }
  
    set(x: any, y: any) {
      this.x = x;
      this.y = y;
  
      return this;
    }
  
    subVectors(a: { x: number; y: number; }, b: { x: number; y: number; }) {
      this.x = a.x - b.x;
      this.y = a.y - b.y;
  
      return this;
    }
  
    copy(v: { x: any; y: any; }) {
      this.x = v.x;
      this.y = v.y;
  
      return this;
    }
  
    clone() {
      return new Vector2(this.x, this.y);
    }
  }