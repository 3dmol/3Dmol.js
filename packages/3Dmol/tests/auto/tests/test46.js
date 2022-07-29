var $scope = {}
 $scope.volumetricFilesPath="chgcar.list";

  $scope.CHGCARS=[];
  $scope.MODELS=[];

  $scope.done = 0; //for testing.. don't call test callback until finished

  $scope.init = function() {
    $scope.MAIN_VIEWER=viewer;
    $scope.addModelObject("../test_structs/CONTCAR", true); // Maybe erase in the future
    $scope.addChgcarObject("../test_structs/CHGCAR"); // Maybe erase in the future
    $scope.MAIN_VIEWER.setBackgroundColor(0xffffff);
  }

  $scope.addModelObject = function (name, value) {
    var model = {};
      //this would suit me
      model.format="vasp";
    
    model.name = name;
    model.value = value? true: false;
    $scope.MODELS.push(model);
  }

  $scope.addChgcarObject = function(name) {
    var chgcarObject={};
    var format = name.match(/\..*/);
      chgcarObject.format="vasp"; //It suits my needs
    
    chgcarObject.name = name;
    chgcarObject.data = false;
    chgcarObject.value = true;
    chgcarObject.isovalue = 0.01;
    chgcarObject.opacity = 0.95;
    chgcarObject.alpha = 0.5;
    chgcarObject.smoothness = 1;
    chgcarObject.voxel = false;
    chgcarObject.color = "blue";
    $scope.CHGCARS.push(chgcarObject);
  }

  $scope.getChgcarNames = function() {
    //console.log("Reading chgcars from  "+$scope.volumetricFilesPath);
    $.get($scope.volumetricFilesPath, function(data){
      //console.log(data);
      data = data.split(/[\n\r]/);
      data.forEach(function(name){
        if (name) {
          $scope.addChgcarObject(name);
          //console.log($scope.CHGCARS);
        }
      });
    });
  }


  $scope.clear = function() {
    //console.log("Clearing..");
    $scope.MAIN_VIEWER.clear();
  }



  $scope.render = function() {
    $scope.clear();
    $scope.renderModels();
    $scope.renderChgcar();
    $scope.MAIN_VIEWER.setSlab(-50,50);
    $scope.MAIN_VIEWER.render();
  }

  $scope.renderIso = function() {
    $scope.MAIN_VIEWER.removeAllShapes();
    $scope.MAIN_VIEWER.removeAllSurfaces();
    $scope.renderChgcar();
  }

  $scope.renderVolumetricData = function (chgcarObject) {
    var isovalue   = chgcarObject.isovalue;
    var opacity    = chgcarObject.opacity;
    var alpha      = chgcarObject.alpha;
    var smoothness = chgcarObject.smoothness;
    var voxel      = chgcarObject.voxel;
    var format     = chgcarObject.format;
    var color     = chgcarObject.color;
    var volumetric_path = chgcarObject.name;
    if (chgcarObject.data) {
      $scope.MAIN_VIEWER.addIsosurface(chgcarObject.data , {voxel:voxel , isoval: isovalue  , color: color, opacity:opacity , smoothness:smoothness , alpha: alpha});
      $scope.MAIN_VIEWER.render();
    } else {
      //console.log("Loading volumetric_data from "+volumetric_path);
      $.get(volumetric_path, function (data) {
        //console.log("Volumetric data received");
        var voldata    = new $3Dmol.VolumeData(data, format);
        chgcarObject.data = voldata;
        $scope.MAIN_VIEWER.addIsosurface(voldata , {voxel:voxel , isoval: isovalue  , color: color, opacity:opacity , smoothness:smoothness , alpha: alpha});
        $scope.MAIN_VIEWER.render();
        $scope.done++;
      }).fail(function(err) {console.log(err);});
    }
  }

  $scope.renderChgcar = function() {
    //console.log("Rendering Chgcar");
    $scope.CHGCARS.forEach(function(chgcarObject, index){
      if (chgcarObject.value) {
        $scope.renderVolumetricData(chgcarObject);
      }
    });

  }

  $scope.renderModel = function (model) {
    var modelPath=model.name;
    var format=model.format;
    $.get(modelPath,function(data){
      //console.log("Structural data received");
      var model = $scope.MAIN_VIEWER.addModel(data, format);
      model.setStyle({}, {sphere:{scale: 0.2}, stick:{radius:0.1}});
      $scope.MAIN_VIEWER.zoomTo();
      $scope.MAIN_VIEWER.render();
      $scope.done++;
    });
  }
  $scope.renderModels = function() {
    $scope.MODELS.forEach(function(model){
      $scope.renderModel(model);
    });
  }

  $scope.init();
  $scope.render();

  function checkdone() {
        // Check if condition met. If not, re-check later (msec).
        if ($scope.done !== 2) {
            setTimeout(checkdone, 500);            
        } else {
            viewer.render(callback);
        }
  }
  checkdone();
