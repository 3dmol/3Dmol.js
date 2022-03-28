var $scope = {}
	
$scope.MODELS=[];

$scope.init = function() {
  $scope.MAIN_VIEWER=viewer;
  $scope.addModelObject("../test_structs/CONTCAR", true);
  $scope.MAIN_VIEWER.setBackgroundColor(0xffffff);
}

$scope.addModelObject = function (name, value) {
  var model = {};
model.format="vasp";    
  model.name = name;
  model.value = value? true: false;
  $scope.MODELS.push(model);
}

$scope.clear = function() {
  $scope.MAIN_VIEWER.clear();
}

$scope.render = function() {
  $scope.clear();
  $scope.renderModels();
  $scope.MAIN_VIEWER.setSlab(-50,50);
  $scope.MAIN_VIEWER.render();
}

$scope.renderModel = function (model) {
  var modelPath=model.name;
  var format=model.format;
  $.get(modelPath,function(data){
    var model = $scope.MAIN_VIEWER.addModel(data, format);
  $scope.MAIN_VIEWER.addUnitCell(model);
    model.setStyle({}, {sphere:{scale: 0.2}, stick:{radius:0.1}});	  
    $scope.MAIN_VIEWER.zoomTo();
    $scope.MAIN_VIEWER.render();
    viewer.render();
  });
}

$scope.renderModels = function() {
  $scope.MODELS.forEach(function(model){
    $scope.renderModel(model);
  });
}

$scope.init();
$scope.renderModels();
