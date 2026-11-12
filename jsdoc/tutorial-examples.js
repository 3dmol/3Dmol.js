'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');

// Serialized into the browser bundle with the example molecule. Keeping the
// data beside the renderer avoids AJAX restrictions when viewing local docs.
function initializeExamples(pdb) {
    function example(buttonId, draw, backgroundColor, zoom) {
        const button = document.getElementById(buttonId);
        if (!button) return;
        const container = document.getElementById(button.getAttribute('aria-controls'));
        const status = button.parentElement.querySelector('.tutorial-status');
        let viewer;
        button.addEventListener('click', () => {
            try {
                status.textContent = '';
                viewer ||= $3Dmol.createViewer(container, { backgroundColor });
                viewer.clear();
                draw(viewer);
                viewer.zoomTo();
                viewer.zoom(zoom);
                viewer.render();
            } catch (error) {
                status.textContent = 'Unable to display the example. Please reload the page and try again.';
                console.error(error);
            }
        });
    }

    example('btn-01', viewer => {
        viewer.addSphere({ center: { x: 0, y: 0, z: 0 }, radius: 10, color: 'green' });
    }, 'orange', 0.8);

    example('btn-02', viewer => {
        viewer.addModel(pdb, 'pdb');
        viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
    }, 'white', 1.2);

    example('btn-03', viewer => {
        viewer.addModel(pdb, 'pdb');
        viewer.setStyle({}, { cartoon: {
            colorfunc: atom => atom.resi % 2 === 0 ? 'white' : 'green'
        } });
    }, 'white', 1.2);
}

async function prepareTutorialExamples(destination) {
    const root = path.resolve(__dirname, '..');
    const assets = path.resolve(destination, '_assets');
    const [library, pdb] = await Promise.all([
        fs.readFile(path.join(root, 'build/3Dmol-min.js'), 'utf8'),
        fs.readFile(path.join(root, 'tutorials/doc-data/1ycr.pdb'), 'utf8')
    ]);
    await fs.mkdir(assets, { recursive: true });
    await fs.writeFile(path.join(assets, 'tutorial.js'),
        library + `\n;(${initializeExamples.toString()})(${JSON.stringify(pdb)});\n`);
    await fs.copyFile(path.join(root, 'build/3Dmol-min.js.LICENSE.txt'),
        path.join(assets, '3Dmol-min.js.LICENSE.txt'));
}

module.exports = { prepareTutorialExamples };
