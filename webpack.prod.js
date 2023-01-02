/* eslint-disable no-undef*/

const { merge } = require('webpack-merge');
const common = require('./webpack.common.js');
const TerserPlugin = require('terser-webpack-plugin');

module.exports = merge(common, {
  mode: 'production',
  optimization: {
    usedExports: false, //essential for pulling in $3Dmol namespace
    minimizer: [
      new TerserPlugin({
        parallel: true,
        terserOptions: {
          keep_classnames: true,
          compress: true,
          mangle:
          {
            reserved: ["MarchingCube"],
          }
        },
      }),
    ],
  },
  output: {
    filename: `[name]-min.js`,
  }
});
