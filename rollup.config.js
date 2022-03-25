import resolve from "@rollup/plugin-node-resolve";
import typescript from "@rollup/plugin-typescript";
import multi from "@rollup/plugin-multi-entry";
import { terser } from "rollup-plugin-terser";

const mainsrc = ["3Dmol/3dmol.js", "3Dmol/colors.js", "3Dmol/*.js"];
const webglsrc = [
  "3Dmol/WebGL/math.js",
  "3Dmol/WebGL/core.js",
  "3Dmol/WebGL/shapes.js",
  "3Dmol/WebGL/**/*.js",
];
const extsrc = ["js/mmtf.js"];
const uisrc = [
  "3Dmol/ui/ui.js",
  "3Dmol/ui/state.js",
  "3Dmol/ui/icon.js",
  "3Dmol/ui/form.js",
  "3Dmol/ui/defaultValues.js",
  "3Dmol/ui/**/*.js",
];
const ignoresrc = ["3Dmol/SurfaceWorker.js"];

const commonOutputConfig = {
  name: "$3Dmol",
  format: "umd",
  sourcemap: true,
  strict: true,
  exports: "none",
  externalLiveBindings: true,
  extend: true,
  intro: "var $3Dmol = window.$3Dmol || {};",
  outro: "window.$3Dmol = $3Dmol;",
};

const include = [...webglsrc, ...mainsrc, ...extsrc, ...uisrc]
console.log(include)

export default {
  input: {
    include,
    exclude: ignoresrc,
  },
  output: [
    {
      file: "build/3Dmol.js",
      ...commonOutputConfig,
    },
    {
      file: "build/3Dmol.min.js",
      plugins: [terser()],
      ...commonOutputConfig,
    },
    {
      file: "build/3Dmol-nojquery.js",
      globals: {
        jquery: "jQuery",
        jquery: "$",
      },
      ...commonOutputConfig,
    },
    {
      file: "build/3Dmol-nojquery-min.js",
      globals: {
        jquery: "jQuery",
        jquery: "$",
      },
      plugins: [terser()],
      ...commonOutputConfig,
    },
  ],
  plugins: [
    resolve(),
    //typescript(),
    multi(),
  ],
};
