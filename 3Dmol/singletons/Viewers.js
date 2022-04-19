/**
 * @type {Record<string|symbol, import("../GLViewer").default>}
 */
const viewers = {};

export default new Proxy(viewers, {
    set: (target, key, value) => {
        target[key] = value;
        return true;
    },

    get: (target, key) =>  target[key]
})