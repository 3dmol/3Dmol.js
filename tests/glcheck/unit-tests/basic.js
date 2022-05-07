import $3Dmol from '../../../build/3Dmol/index.js';

glcheck('createViewer', (t, canvas) => {
  const viewer = $3Dmol.createViewer(
    $('<div>', {
      canvas,
    })
  );
  t.ok(viewer, 'viewer created');
  t.done();
});

glcheck('GLViewer.addLabel', (t, canvas) => {
  const viewer = $3Dmol.createViewer(
    $('<div>', {
      canvas,
    })
  );
  const label = viewer.addLabel('Hello world', {
    position: {x: 0, y: 0, z: 0},
    backgroundColor: 'orange',
  });
  t.ok(label, 'label created');
  viewer.render(t.done);
});

glcheck('GLViewer.addLabel animated', async (t, canvas) => {
  const viewer = $3Dmol.createViewer(
    $('<div>', {
      canvas,
    })
  );
  for (var i = 0; i < 10; i++) {
    viewer.addLabel('Hello ' + i, {
      position: {x: i, y: 0, z: i},
      backgroundColor: 'orange',
      frame: i,
    });
  }
  let done = false;
  viewer.animate({loop: 'forward', reps: 1});
  viewer.render();
  await t.loopUntil(() => viewer.surfacesFinished() && !viewer.isAnimated() && !$3Dmol.processing_autoinit);
  t.done();
});
