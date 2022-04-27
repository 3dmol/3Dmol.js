export default function multiLineString(f) {
  return f
    .toString()
    .replace(/^[^/]+\/\*!?/, '')
    .replace(/\*\/[^/]+$/, '');
}
