import resolve from "@rollup/plugin-node-resolve";
import commonjs from "@rollup/plugin-commonjs";
import typescript from "@rollup/plugin-typescript";
//import { terser } from "rollup-plugin-terser";

export default {
    input: ["3Dmol/WebGL/index.ts"],
    output: [
        {
            dir: "tmp/",
            name: "$3Dmol",
            format: "umd",
            sourcemap: true,
        }
    ],
    
    plugins:[
        typescript(),
        resolve(),
        commonjs(),
    ]
}