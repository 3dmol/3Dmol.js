import type { Object3D } from "./core/Object3D";
import type { Geometry } from "./core/Geometry";
import type { Texture } from "./core/Texture";
import type { Material } from "./materials/Material";
import { VolumetricMaterial } from "./materials/VolumetricMaterial";

const pendingMaterials = new Set<Material>();
let materialFlushScheduled = false;

// Internal scheduler shared by scene owners and renderer-created pass materials.
export function disposeMaterials(materials: Iterable<Material>): void {
    for (const material of materials) pendingMaterials.add(material);
    if (pendingMaterials.size && !materialFlushScheduled) {
        materialFlushScheduled = true;
        // Retire materials at the next Promise microtask, including removals during
        // viewer.render(), so new frame materials can reuse their shader programs.
        Promise.resolve().then(() => {
            const pending = Array.from(pendingMaterials);
            pendingMaterials.clear();
            materialFlushScheduled = false;
            for (const material of pending) material.dispose();
        });
    }
}

/** Release resources owned by a render tree; clones may share those resources. */
export function disposeObjectResources(root: Object3D, options: {
    geometries?: boolean;
    volumeTextures?: boolean;
} = {}): void {
    const geometries = new Set<Geometry>();
    const materials = new Set<Material>();
    const textures = new Set<Texture>();
    const objects = [root];
    while (objects.length) {
        const object = objects.pop();
        if (options.geometries !== false && object.geometry) geometries.add(object.geometry);
        const material = object.material;
        if (material) {
            materials.add(material);
            if (options.volumeTextures && material instanceof VolumetricMaterial) {
                if (material.map) textures.add(material.map);
                if (material.transferfn) textures.add(material.transferfn);
            }
        }
        for (const child of object.children) objects.push(child);
    }

    for (const geometry of geometries) geometry.dispose();
    for (const texture of textures) texture.dispose();
    disposeMaterials(materials);
}
