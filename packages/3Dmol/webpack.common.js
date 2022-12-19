// tslint:disable:object-literal-sort-keys
const path = require("path");
const webpack = require("webpack");
const pkg = require("./package.json");
const ESLintPlugin = require('eslint-webpack-plugin');


const banner = `${pkg.name} v${pkg.version}
${pkg.description}
Author: ${pkg.author}`;

module.exports =  {
  target: "web",
  mode: "production",
  entry: ["./src/index.ts", 
        "./src/3dmol.js",
        "./src/autoload.js", 
        "./src/glcartoon.js",
        "./src/gldraw.js",
        "./src/glmodel.js",
        "./src/glshape.js", 
        "./src/glviewer.js",
        "./src/SurfaceWorker.js",
        "./src/ui/ui.js",
        "./src/ui/state.js",
        "./src/ui/icon.js",
        "./src/ui/form.js",        
        "./src/ui/defaultValues.js",
        "./src/exporter.js"
        ],
  output: {
    path: path.resolve(__dirname, "build"),
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
      { test: /\.vert/, loader: "raw-loader" }
    ],
  },

  plugins: [
    new webpack.ProvidePlugin({
        $: "jquery",
        jQuery: "jquery",
        pako: "pako",
        UPNG: "upng-js",
        netcdfjs: "netcdfjs",
        MMTF: path.resolve(__dirname, "./src/vendor/mmtf.js")
    }),
    new webpack.BannerPlugin({ banner }), 
    new ESLintPlugin()],
};
