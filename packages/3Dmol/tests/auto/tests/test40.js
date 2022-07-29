function triangle(viewer) {
    var vertices = [];
    var normals = [];
    var colors = [];
    var r = 20;
    //triangle
    vertices.push(new $3Dmol.Vector3(0,0,0));
    vertices.push(new $3Dmol.Vector3(r,0,0));
    vertices.push(new $3Dmol.Vector3(0,r,0));
    
    normals.push(new $3Dmol.Vector3(0,0,1));
    normals.push(new $3Dmol.Vector3(0,0,1));
    normals.push(new $3Dmol.Vector3(0,0,1));
    
    colors.push({r:1,g:0,b:0});
    colors.push({r:0,g:1,b:0});
    colors.push({r:0,g:0,b:1});

    var faces = [ 0,1,2 ];
    
    var spec = {vertexArr:vertices, normalArr: normals, faceArr:faces,color:colors};
    viewer.addCustom(spec);
}

function cylinder(viewer) {
    var vertices = [];
    var normals = [];
    var colors = [];
    var r = 10;

    // "Cylinder" with 4 points (aka a open box with weird normals)
    vertices.push(new $3Dmol.Vector3(r,-2*r,0));
    vertices.push(new $3Dmol.Vector3(0,-2*r,r));
    vertices.push(new $3Dmol.Vector3(-r,-2*r,0));
    vertices.push(new $3Dmol.Vector3(0,-2*r,-r));
    vertices.push(new $3Dmol.Vector3(r,-r,0));
    vertices.push(new $3Dmol.Vector3(0,-r,r));
    vertices.push(new $3Dmol.Vector3(-r,-r,0));
    vertices.push(new $3Dmol.Vector3(0,-r,-r));
    
    normals.push(new $3Dmol.Vector3(1,0,0));
    normals.push(new $3Dmol.Vector3(0,0,1));
    normals.push(new $3Dmol.Vector3(-1,0,0));
    normals.push(new $3Dmol.Vector3(0,0,-1));
    normals.push(new $3Dmol.Vector3(1,0,0));
    normals.push(new $3Dmol.Vector3(0,0,1));
    normals.push(new $3Dmol.Vector3(-1,0,0));
    normals.push(new $3Dmol.Vector3(0,0,-1));

    colors.push({r:1,g:0,b:0});
    colors.push({r:1,g:1,b:0});
    colors.push({r:1,g:1,b:1});
    colors.push({r:1,g:0,b:1});
    colors.push({r:.5,g:0,b:0});
    colors.push({r:.5,g:.5,b:0});
    colors.push({r:.5,g:.5,b:.5});
    colors.push({r:.5,g:0,b:.5});

    var faces = [ 0,4,1, 4,5,1,
            1,5,2, 5,6,2,
            2,6,3, 6,7,3,
            3,7,0, 7,4,0 ];

    var spec = {vertexArr:vertices, normalArr: normals, faceArr:faces,color:colors};
    viewer.addCustom(spec);
}


            

            triangle(viewer);
            cylinder(viewer);
            viewer.render();
