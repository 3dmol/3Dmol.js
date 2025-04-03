const {
    defineConfig,
    globalIgnores,
} = require("eslint/config");

const babelParser = require("@babel/eslint-parser");
const js = require("@eslint/js");

const {
    FlatCompat,
} = require("@eslint/eslintrc");

const compat = new FlatCompat({
    baseDirectory: __dirname,
    recommendedConfig: js.configs.recommended,
    allConfig: js.configs.all
});

module.exports = defineConfig([globalIgnores(["**/vendor/*.js"]), {
    extends: compat.extends("eslint:recommended"),

    languageOptions: {
        globals: {
            $3Dmol: "writable",
            console: "writable",
            document: "readonly",
            $: "readonly",
            window: "readonly",
            self: "writable",
            module: "readonly",
            Blob: "readonly",
            pako: "readonly",
            UPNG: "readonly",
            netcdfjs: "readonly",
            XMLHttpRequest: "readonly",
            alert: "readonly",
            setTimeout: "readonly",
            clearTimeout: "readonly",
            setInterval: "readonly",
            clearInterval: "readonly",
            Worker: "readonly",
            MMTF: "readonly",
            TextDecoder: "readonly",
            FileReader: "readonly",
            ProteinSurface: "readonly",
        },

        parser: babelParser,
        ecmaVersion: 5,
        sourceType: "script",

        parserOptions: {
            requireConfigFile: false,
        },
    },

    rules: {
        "no-prototype-builtins": "off",
    },
}]);