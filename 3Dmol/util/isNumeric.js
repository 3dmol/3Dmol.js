const isNumeric = obj => {
  const type = typeof obj;
  return (type === 'number' || type === 'string') && !Number.isNaN(obj - Number.parseFloat(obj));
};

export default isNumeric;
