const $3Dmol = require('../../3Dmol');
const properties = require("./properties.json")

test(() => {
    expect(JSON.stringify(properties)).toEqual(JSON.stringify($3Dmol))
})