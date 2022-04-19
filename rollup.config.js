import resolve from "@rollup/plugin-node-resolve";
import commonjs from "@rollup/plugin-commonjs";
import { terser} from "rollup-plugin-terser";

export default {
    input: "3Dmol/index.js",
    output: [
        {
            file: "build/3Dmol.js",
            name: "$3Dmol",
            format: "umd",
            sourcemap: true,
        },
        {
            file: "build/3Dmol-min.js",
            name: "$3Dmol",
            format: "umd",
            sourcemap: true,
            plugins: [terser()],
        },
        {
            dir: "build",
            format: "es",
            sourcemap: true,
            preserveModules: true,
        }
    ],
    
    plugins:[
        resolve(),
        commonjs(),
    ]
}