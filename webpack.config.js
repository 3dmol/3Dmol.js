// tslint:disable:object-literal-sort-keys
const path = require("path");
const webpack = require("webpack");
const pkg = require("./package.json");

const banner = `${pkg.name} v${pkg.version}
${pkg.description}
Author: ${pkg.author}`;

module.exports = {
  target: "web",
  mode: "production",
  entry: "./3Dmol/index.ts",
  output: {
    filename: `index.js`,
    path: path.resolve(__dirname, "tmp"),
    globalObject: "this",
    library: pkg.name,
    libraryTarget: "umd",
  },
  resolve: {
    extensions: [".ts", ".tsx", ".js", ".json"],
  },

  module: {
    rules: [
      { test: /\.tsx?$/, loader: "ts-loader"},
      { test: /\.frag/, loader: "raw-loader" },
      { test: /\.vert/, loader: "raw-loader" },
    ],
  },

  plugins: [new webpack.BannerPlugin({ banner })],
};
