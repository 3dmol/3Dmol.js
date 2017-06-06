/*
@data
<textarea style="display: none;" id="test">
4
* (null), Energy   -1000.0000000
N     0.000005    0.019779   -0.000003   -0.157114    0.000052   -0.012746
H     0.931955   -0.364989    0.000003    1.507100   -0.601158   -0.004108
H    -0.465975   -0.364992    0.807088    0.283368    0.257996   -0.583024
H    -0.465979   -0.364991   -0.807088    0.392764    0.342436    0.764260
</textarea>
*/
            var data = $("#test").val();
            viewer.setBackgroundColor(0xffffff);    
            var m = viewer.addModel(data, "xyz");
            m.setStyle({},{stick:{}});
            m.vibrate(10, 1);
            viewer.animate({loop: "forward",reps:1});
            viewer.zoomTo();
            viewer.render();
