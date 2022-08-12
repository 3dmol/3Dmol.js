var freq_data = `14
                        Demo
                        C 1.4503 -0.3344 -0.0999 -0.02 -0.04 0.05
                        C 1.6983 -0.7959 1.1842 -0.33 0.26 0.02
                        C 2.4164 0.4823 -0.7083 0.42 -0.27 -0.05
                        O 3.8281 -0.6667 -0.9197 -0.31 0.22 0.01
                        C 4.1148 -1.1750 0.2217 -0.03 -0.04 0.07
                        C 3.3387 -2.1631 0.8005 0.35 -0.20 -0.10
                        H 0.7647 -0.8960 -0.7408 -0.02 0.00 0.02
                        H 2.2562 -0.1690 1.8842 0.08 -0.03 -0.06
                        H 1.0370 -1.5313 1.6497 -0.15 0.08 -0.00
                        H 2.9505 1.2096 -0.0914 -0.06 0.08 -0.02
                        H 2.2916 0.7816 -1.7506 0.22 -0.13 0.01
                        H 4.8526 -0.6497 0.8593 -0.01 0.02 -0.05
                        H 2.7188 -2.7970 0.1653 -0.09 0.08 0.06
                        H 3.6193 -2.5617 1.7790 0.27 -0.14 -0.06`


viewer.addModel(freq_data, "xyz");
viewer.vibrate(10, 0.8, true, null);
viewer.animate({'loop': 'backAndForth'});
viewer.setStyle({}, {stick:{}, sphere: {radius: 0.5}});
viewer.zoomTo();
viewer.render(/* no callback */ );

function demonstration(ind) {
  viewer.stopAnimate(); 
  viewer.vibrate(10, 0.8, true, null);
  //Without the above line, there does not seem to be a problem
  //However, the amplitude cannot be updated otherwise

  viewer.animate({'loop': 'backAndForth',reps:1});

  if(ind == 0) {
    viewer.render(); 
    return;
   }

  setTimeout(function() { demonstration(ind-1)}, 100);
}
demonstration(11); //Call the above function 10 times in 1 second
