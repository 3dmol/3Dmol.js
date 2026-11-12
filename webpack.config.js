/* eslint-disable no-undef*/
const path = require("path");
const webpack = require("webpack");
const { merge } = require("webpack-merge");
const pkg = require("./package.json");
const ESLintPlugin = require("eslint-webpack-plugin");
const TerserPlugin = require('terser-webpack-plugin');
const { getBuildInfo } = require('./scripts/build-info.cjs');

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
    new ESLintPlugin({ exclude: ["node_modules"] }),
  ],
  stats: {
    errorDetails: true,
  },
};

umd_config = {
  entry: {
    "3Dmol": ["./src/index.ts", "./src/SurfaceWorker.js", "./src/exporter.js"],
    "3Dmol.ui": ["./src/ui/index.ts"],
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
  const buildInfo = getBuildInfo(__dirname, pkg.version);
  const banner = `${pkg.name} v${buildInfo.version}
Git: ${buildInfo.gitDescribe || 'unavailable'}
Commit: ${buildInfo.commit || 'unavailable'}
${pkg.description}
Author: ${pkg.author}`;
  const build_config = {
    plugins: [
      new webpack.DefinePlugin({
        __3DMOL_BUILD_INFO__: JSON.stringify(buildInfo),
      }),
      new webpack.BannerPlugin({
        banner,
        // Retain the banner after minification, before source maps are generated.
        stage: webpack.Compilation.PROCESS_ASSETS_STAGE_OPTIMIZE_SIZE + 1,
      }),
    ],
  };
  return [
    merge(base_config, umd_config, env_config[argv["mode"]], build_config),
    merge(base_config, es6_config, env_config[argv["mode"]], build_config),
  ];
};
