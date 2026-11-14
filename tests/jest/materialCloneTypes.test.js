/** @jest-environment node */
const path = require("path");
const ts = require("typescript");

test("published material clone types match constructed and supplied destinations", () => {
    const names = [
        "Material", "ImposterMaterial", "InstancedMaterial", "LineBasicMaterial",
        "MeshDoubleLambertMaterial", "MeshLambertMaterial", "MeshOutlineMaterial",
        "RingImposterMaterial", "SphereImposterMaterial", "SphereImposterOutlineMaterial",
        "SpriteMaterial", "StickImposterMaterial", "StickImposterOutlineMaterial", "VolumetricMaterial"
    ];
    const filename = path.join(__dirname, "material-clone-types-consumer.ts");
    const source = `
        ${names.map(name => `import { ${name} } from '../../build/types/WebGL/materials/${name}';`).join("\n")}
        const materials: Material[] = [${names.map(name => `new ${name}()`).join(", ")}];
        const baseClone: Material = materials[0].clone();

        ${names.map(name => `
            const ${name}Clone: ${name} = new ${name}().clone();
            const ${name}Undefined: ${name} = new ${name}().clone(undefined);
            class ${name}Child extends ${name} { customValue = 123; }
            const ${name}Destination = new ${name}Child();
            const ${name}Typed: ${name}Child = new ${name}().clone(${name}Destination);
            const ${name}Explicit: ${name}Child = new ${name}().clone<${name}Child>(${name}Destination);
            const ${name}Optional: ${name} = new ${name}().clone(Math.random() ? ${name}Destination : undefined);

            // @ts-expect-error inherited clone constructs the builtin parent, not the custom subclass
            const ${name}Wrong: ${name}Child = new ${name}Child().clone();
            // @ts-expect-error a generic subclass result requires a supplied destination
            new ${name}().clone<${name}Child>();
        `).join("\n")}

        // @ts-expect-error a Lambert destination must have the fields that its clone copies
        new MeshLambertMaterial().clone(new Material());
        // @ts-expect-error a line destination must have the fields that its clone copies
        new LineBasicMaterial().clone(new Material());

        // Existing custom implementations can keep their generic-this override.
        class LegacyCustomMaterial extends MeshLambertMaterial {
            custom = true;
            clone<T extends this>(material: T = new LegacyCustomMaterial() as T): T {
                super.clone.call(this, material);
                material.custom = this.custom;
                return material;
            }
        }
        const legacy: Material = new LegacyCustomMaterial();
        const legacyClone: LegacyCustomMaterial = new LegacyCustomMaterial().clone();
    `;
    const options = {
        noEmit: true,
        strict: true,
        skipLibCheck: true,
        target: ts.ScriptTarget.ES2015,
        module: ts.ModuleKind.ESNext,
        moduleResolution: ts.ModuleResolutionKind.Node10,
        baseUrl: path.resolve(__dirname, "../../build/types"),
        types: []
    };
    const host = ts.createCompilerHost(options);
    const readFile = host.readFile;
    host.readFile = file => file === filename ? source : readFile(file);
    const program = ts.createProgram([filename], options, host);
    const diagnostics = ts.getPreEmitDiagnostics(program);
    expect(ts.formatDiagnostics(diagnostics, host)).toBe("");
});
