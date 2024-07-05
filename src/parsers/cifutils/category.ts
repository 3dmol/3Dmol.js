/* Note - all code in this directory is adapted from molstar -dkoes */
import { EncodedCategory } from "./encoding";
import { decode } from './decoder';

export function Category(data: EncodedCategory) {
  const map = Object.create(null);
  const cache = Object.create(null);
  for (const col of data.columns) map[col.name] = col;
  return {
      rowCount: data.rowCount,
      name: data.name.substring(1),
      fieldNames: data.columns.map(c => c.name),
      getField(name) {
          const col = map[name];
          if (!col) return void 0;
          if (!!cache[name]) return cache[name];
          cache[name] = decode(col.data);
          return cache[name];
      }
  };
}