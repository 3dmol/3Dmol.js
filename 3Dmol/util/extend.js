// simplified version of $.extend
const extend = (obj1, src1) => {
  for (const key in src1) {
    if (src1[key] && src1[key] !== undefined) {
      // if(Object.prototype.hasOwnProperty.call(src1,key) && src1[key] !== undefined){ // use Object.prototype
      obj1[key] = src1[key];
    }
  }
  return obj1;
};

export default extend;
