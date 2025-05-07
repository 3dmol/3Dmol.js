/* eslint-disable no-undef*/
const path = require("path");
const webpack = require("webpack");
const { merge } = require("webpack-merge");
const pkg = require("./package.json");
const ESLintPlugin = require("eslint-webpack-plugin");
const TerserPlugin = require('terser-webpack-plugin');

const banner = `${pkg.name} v${pkg.version}
${pkg.description}
Author: ${pkg.author}`;

base_config = {
  target: "web",
  output: {
    path: path.resolve(__dirname, "build"),
  },
  experiments: {
    outputModule: true, // Required for ES modules
  },
  resolve: {
    extensions: [".ts", ".tsx", ".js", ".json"],
  },
  module: {
    rules: [
      { test: /\.tsx?$/, loader: "ts-loader" },
      { test: /\.frag$/, loader: "raw-loader" },
      { test: /\.vert$/, loader: "raw-loader" },
    ],
  },
  plugins: [
    new webpack.ProvidePlugin({
      MMTF: path.resolve(__dirname, "./src/vendor/mmtf.js"),
      $: "jquery",
    }),
    new webpack.BannerPlugin({ banner }),
    new ESLintPlugin({ exclude: ["node_modules"] }),
  ],
  stats: {
    errorDetails: true,
  },
};

umd_config = {
  entry: {
    "3Dmol": ["./src/index.ts", "./src/SurfaceWorker.js", "./src/exporter.js"],
    "3Dmol.ui": [
      "./src/ui/ui.js",
      "./src/ui/state.js",
      "./src/ui/icon.js",
      "./src/ui/form.js",
      "./src/ui/defaultValues.js",
    ],
  },
  output: {
    library: "[name]",
    libraryTarget: "umd",
    globalObject: "this",
  },
};

es6_config = {
  entry: {
    "3Dmol.es6": ["./src/3Dmol.ts"],
  },
  output: {
    libraryTarget: "module",
  },
}

env_config = {
  development: {
    mode: "development",
    devtool: "inline-source-map",
    output: {
      filename: "[name].js",
    },
  },
  production: {
    mode: "production",
    optimization: {
      usedExports: false,
      minimizer: [
        new TerserPlugin({
          parallel: true,
          terserOptions: {
            keep_classnames: true,
            compress: true,
            mangle: {
              reserved: ["MarchingCube"],
            },
          },
        }),
      ],
    },
    output: {
      filename: `[name]-min.js`,
    },
  },
};

module.exports = function (env, argv) {
  return [
    merge(base_config, umd_config, env_config[argv["mode"]]),
    merge(base_config, es6_config, env_config[argv["mode"]]),
  ];
};
