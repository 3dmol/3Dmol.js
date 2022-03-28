/* eslint-disable eqeqeq */
/**
 * A GLShape is a collection of user specified shapes.
 *
 * @constructor $3Dmol.GLShape
 * @extends {ShapeSpec}
 * @param {number} sid - Unique identifier
 * @param {ShapeSpec} stylespec - shape style specification
 */

$3Dmol.GLShape = (function glshapeConstructor() {
  // Marching cube, to match with protein surface generation
  const ISDONE = 2;

  const finalizeGeo = function finalizeGeo(geo) {
    // to avoid creating a bunch of geometries, we leave geoGroup untruncated
    // until render is called, at which point we truncate;
    // successive called up updateGeo will return a new geometry
    const geoGroup = geo.updateGeoGroup(0);
    if (geoGroup.vertices > 0) {
      geoGroup.truncateArrayBuffers(true, true);
    }
  };

  /**
   *
   * @param {$3Dmol.Geometry}
   *            geo
   * @param {$3Dmol.Color |
   *            colorlike} color
   */
  const updateColor = function updateColor(geo, color) {
    color = color || $3Dmol.CC.color(color);
    geo.colorsNeedUpdate = true;

    let r;
    let g;
    let b;
    if (color.constructor !== Array) {
      r = color.r;
      g = color.g;
      b = color.b;
    }

    for (const gg in geo.geometryGroups) {
      const geoGroup = geo.geometryGroups[gg];
      const colorArr = geoGroup.colorArray;

      for (let i = 0, il = geoGroup.vertices; i < il; ++i) {
        if (color.constructor === Array) {
          const c = color[i];
          r = c.r;
          g = c.g;
          b = c.b;
        }

        colorArr[i * 3] = r;
        colorArr[i * 3 + 1] = g;
        colorArr[i * 3 + 2] = b;
      }
    }
  };

  /**
   * @param {$3Dmol.GLShape}
   *            shape
   * @param {geometryGroup}
   *            geoGroup
   * @param {ArrowSpec}
   *            spec
   */
  const drawArrow = function drawArrow(shape, geo, spec) {
    const from = spec.start;
    const {end} = spec;
    const {radius} = spec;
    const {radiusRatio} = spec;
    let {mid} = spec;
    const midoffset = spec.midpos;

    if (!(from && end)) return;

    const geoGroup = geo.updateGeoGroup(51);

    // vertices

    const dir = end.clone().sub(from);
    if (midoffset) {
      // absolute offset, convert to relative
      const length = dir.length();
      if (midoffset > 0) mid = midoffset / length;
      else mid = (length + midoffset) / length;
    }

    dir.multiplyScalar(mid);

    const to = from.clone().add(dir);
    const negDir = dir.clone().negate();

    shape.intersectionShape.cylinder.push(new $3Dmol.Cylinder(from.clone(), to.clone(), radius));
    shape.intersectionShape.sphere.push(new $3Dmol.Sphere(from.clone(), radius));

    // get orthonormal vector
    const nvecs = [];
    nvecs[0] = dir.clone();
    if (Math.abs(nvecs[0].x) > 0.0001) nvecs[0].y += 1;
    else nvecs[0].x += 1;
    nvecs[0].cross(dir);
    nvecs[0].normalize();

    nvecs[0] = nvecs[0];
    // another orth vector
    nvecs[4] = nvecs[0].clone();
    nvecs[4].crossVectors(nvecs[0], dir);
    nvecs[4].normalize();
    nvecs[8] = nvecs[0].clone().negate();
    nvecs[12] = nvecs[4].clone().negate();

    // now quarter positions
    nvecs[2] = nvecs[0].clone().add(nvecs[4]).normalize();
    nvecs[6] = nvecs[4].clone().add(nvecs[8]).normalize();
    nvecs[10] = nvecs[8].clone().add(nvecs[12]).normalize();
    nvecs[14] = nvecs[12].clone().add(nvecs[0]).normalize();

    // eights
    nvecs[1] = nvecs[0].clone().add(nvecs[2]).normalize();
    nvecs[3] = nvecs[2].clone().add(nvecs[4]).normalize();
    nvecs[5] = nvecs[4].clone().add(nvecs[6]).normalize();
    nvecs[7] = nvecs[6].clone().add(nvecs[8]).normalize();
    nvecs[9] = nvecs[8].clone().add(nvecs[10]).normalize();
    nvecs[11] = nvecs[10].clone().add(nvecs[12]).normalize();
    nvecs[13] = nvecs[12].clone().add(nvecs[14]).normalize();
    nvecs[15] = nvecs[14].clone().add(nvecs[0]).normalize();

    const start = geoGroup.vertices;
    const {vertexArray} = geoGroup;
    const {faceArray} = geoGroup;
    const {normalArray} = geoGroup;
    const {lineArray} = geoGroup;

    let offset;
    let i;
    let n;
    // add vertices, opposing vertices paired together
    for (i = 0, n = nvecs.length; i < n; ++i) {
      offset = 3 * (start + 3 * i);
      const bottom = nvecs[i].clone().multiplyScalar(radius).add(from);
      const top = nvecs[i].clone().multiplyScalar(radius).add(to);
      const conebase = nvecs[i]
        .clone()
        .multiplyScalar(radius * radiusRatio)
        .add(to);

      vertexArray[offset] = bottom.x;
      vertexArray[offset + 1] = bottom.y;
      vertexArray[offset + 2] = bottom.z;

      vertexArray[offset + 3] = top.x;
      vertexArray[offset + 4] = top.y;
      vertexArray[offset + 5] = top.z;

      vertexArray[offset + 6] = conebase.x;
      vertexArray[offset + 7] = conebase.y;
      vertexArray[offset + 8] = conebase.z;

      if (i > 0) {
        const prev_x = vertexArray[offset - 3];
        const prev_y = vertexArray[offset - 2];
        const prev_z = vertexArray[offset - 1];

        const c = new $3Dmol.Vector3(prev_x, prev_y, prev_z);
        const b = end.clone();
        const b2 = to.clone();
        const a = new $3Dmol.Vector3(conebase.x, conebase.y, conebase.z);

        shape.intersectionShape.triangle.push(new $3Dmol.Triangle(a, b, c));
        shape.intersectionShape.triangle.push(new $3Dmol.Triangle(c.clone(), b2, a.clone()));
      }
    }

    geoGroup.vertices += 48;
    offset = geoGroup.vertices * 3;

    // caps
    vertexArray[offset] = from.x;
    vertexArray[offset + 1] = from.y;
    vertexArray[offset + 2] = from.z;

    vertexArray[offset + 3] = to.x;
    vertexArray[offset + 4] = to.y;
    vertexArray[offset + 5] = to.z;

    vertexArray[offset + 6] = end.x;
    vertexArray[offset + 7] = end.y;
    vertexArray[offset + 8] = end.z;

    geoGroup.vertices += 3;

    // now faces
    let face;
    let norm;
    let faceoffset;
    let lineoffset;
    let t1;
    let t2;
    let t2b;
    let t3;
    let t3b;
    let t4;
    let t1offset;
    let t2offset;
    let t2boffset;
    let t3offset;
    let t3boffset;
    let t4offset;
    let n1;
    let n2;
    let n3;
    let n4;
    const fromi = geoGroup.vertices - 3;
    const toi = geoGroup.vertices - 2;
    const endi = geoGroup.vertices - 1;
    const fromoffset = fromi * 3;
    const tooffset = toi * 3;
    const endoffset = endi * 3;
    for (i = 0, n = nvecs.length - 1; i < n; ++i) {
      const ti = start + 3 * i;
      offset = ti * 3;
      faceoffset = geoGroup.faceidx;
      lineoffset = geoGroup.lineidx;

      t1 = ti;
      t1offset = t1 * 3;
      t2 = ti + 1;
      t2offset = t2 * 3;
      t2b = ti + 2;
      t2boffset = t2b * 3;
      t3 = ti + 4;
      t3offset = t3 * 3;
      t3b = ti + 5;
      t3boffset = t3b * 3;
      t4 = ti + 3;
      t4offset = t4 * 3;

      // face = [t1, t2, t4], [t2, t3, t4];
      // face = [t1, t2, t3, t4];

      norm = [nvecs[i], nvecs[i], nvecs[i + 1], nvecs[i + 1]];

      n1 = n2 = nvecs[i];
      n3 = n4 = nvecs[i + 1];

      normalArray[t1offset] = n1.x;
      normalArray[t2offset] = n2.x;
      normalArray[t4offset] = n4.x;
      normalArray[t1offset + 1] = n1.y;
      normalArray[t2offset + 1] = n2.y;
      normalArray[t4offset + 1] = n4.y;
      normalArray[t1offset + 2] = n1.z;
      normalArray[t2offset + 2] = n2.z;
      normalArray[t4offset + 2] = n4.z;

      normalArray[t2offset] = n2.x;
      normalArray[t3offset] = n3.x;
      normalArray[t4offset] = n4.x;
      normalArray[t2offset + 1] = n2.y;
      normalArray[t3offset + 1] = n3.y;
      normalArray[t4offset + 1] = n4.y;
      normalArray[t2offset + 2] = n2.z;
      normalArray[t3offset + 2] = n3.z;
      normalArray[t4offset + 2] = n4.z;

      normalArray[t2boffset] = n2.x;
      normalArray[t3boffset] = n3.x;
      normalArray[t2boffset + 1] = n2.y;
      normalArray[t3boffset + 1] = n3.y;
      normalArray[t2boffset + 2] = n2.z;
      normalArray[t3boffset + 2] = n3.z;

      // sides
      faceArray[faceoffset] = t1;
      faceArray[faceoffset + 1] = t2;
      faceArray[faceoffset + 2] = t4;
      faceArray[faceoffset + 3] = t2;
      faceArray[faceoffset + 4] = t3;
      faceArray[faceoffset + 5] = t4;
      // caps
      faceArray[faceoffset + 6] = t1;
      faceArray[faceoffset + 7] = t4;
      faceArray[faceoffset + 8] = fromi;
      faceArray[faceoffset + 9] = t2b;
      faceArray[faceoffset + 10] = toi;
      faceArray[faceoffset + 11] = t3b;
      // arrowhead
      faceArray[faceoffset + 12] = t2b;
      faceArray[faceoffset + 13] = endi;
      faceArray[faceoffset + 14] = t3b;

      // sides
      lineArray[lineoffset] = t1;
      lineArray[lineoffset + 1] = t2;
      lineArray[lineoffset + 2] = t1;
      lineArray[lineoffset + 3] = t4;
      // lineArray[lineoffset+4] = t2, lineArray[lineoffset+5] = t3;
      lineArray[lineoffset + 4] = t3;
      lineArray[lineoffset + 5] = t4;
      // caps
      lineArray[lineoffset + 6] = t1;
      lineArray[lineoffset + 7] = t4;
      // lineArray[lineoffset+10] = t1, lineArray[lineoffset+11] = fromi;
      // lineArray[lineoffset+12] = t4, lineArray[lineoffset+13] = fromi;

      lineArray[lineoffset + 8] = t2b;
      lineArray[lineoffset + 9] = t2; // toi
      lineArray[lineoffset + 10] = t2b;
      lineArray[lineoffset + 11] = t3b;
      lineArray[lineoffset + 12] = t3;
      lineArray[lineoffset + 13] = t3b; // toi
      // arrowhead
      lineArray[lineoffset + 14] = t2b;
      lineArray[lineoffset + 15] = endi;
      lineArray[lineoffset + 16] = t2b;
      lineArray[lineoffset + 17] = t3b;
      lineArray[lineoffset + 18] = endi;
      lineArray[lineoffset + 19] = t3b;

      geoGroup.faceidx += 15;
      geoGroup.lineidx += 20;
    }
    // final face

    face = [start + 45, start + 46, start + 1, start, start + 47, start + 2];
    norm = [nvecs[15], nvecs[15], nvecs[0], nvecs[0]];

    faceoffset = geoGroup.faceidx;
    lineoffset = geoGroup.lineidx;

    t1 = face[0];
    t1offset = t1 * 3;
    t2 = face[1];
    t2offset = t2 * 3;
    t2b = face[4];
    t2boffset = t2b * 3;
    t3 = face[2];
    t3offset = t3 * 3;
    t3b = face[5];
    t3boffset = t3b * 3;
    t4 = face[3];
    t4offset = t4 * 3;

    n1 = n2 = nvecs[15];
    n3 = n4 = nvecs[0];

    normalArray[t1offset] = n1.x;
    normalArray[t2offset] = n2.x;
    normalArray[t4offset] = n4.x;
    normalArray[t1offset + 1] = n1.y;
    normalArray[t2offset + 1] = n2.y;
    normalArray[t4offset + 1] = n4.y;
    normalArray[t1offset + 2] = n1.z;
    normalArray[t2offset + 2] = n2.z;
    normalArray[t4offset + 2] = n4.z;

    normalArray[t2offset] = n2.x;
    normalArray[t3offset] = n3.x;
    normalArray[t4offset] = n4.x;
    normalArray[t2offset + 1] = n2.y;
    normalArray[t3offset + 1] = n3.y;
    normalArray[t4offset + 1] = n4.y;
    normalArray[t2offset + 2] = n2.z;
    normalArray[t3offset + 2] = n3.z;
    normalArray[t4offset + 2] = n4.z;

    normalArray[t2boffset] = n2.x;
    normalArray[t3boffset] = n3.x;
    normalArray[t2boffset + 1] = n2.y;
    normalArray[t3boffset + 1] = n3.y;
    normalArray[t2boffset + 2] = n2.z;
    normalArray[t3boffset + 2] = n3.z;

    // Cap normals
    dir.normalize();
    negDir.normalize();
    normalArray[fromoffset] = negDir.x;
    normalArray[tooffset] = normalArray[endoffset] = dir.x;
    normalArray[fromoffset + 1] = negDir.y;
    normalArray[tooffset + 1] = normalArray[endoffset + 1] = dir.y;
    normalArray[fromoffset + 2] = negDir.z;
    normalArray[tooffset + 2] = normalArray[endoffset + 2] = dir.z;

    // Final side
    faceArray[faceoffset] = t1;
    faceArray[faceoffset + 1] = t2;
    faceArray[faceoffset + 2] = t4;
    faceArray[faceoffset + 3] = t2;
    faceArray[faceoffset + 4] = t3;
    faceArray[faceoffset + 5] = t4;
    // final caps
    faceArray[faceoffset + 6] = t1;
    faceArray[faceoffset + 7] = t4;
    faceArray[faceoffset + 8] = fromi;
    faceArray[faceoffset + 9] = t2b;
    faceArray[faceoffset + 10] = toi;
    faceArray[faceoffset + 11] = t3b;
    // final arrowhead
    faceArray[faceoffset + 12] = t2b;
    faceArray[faceoffset + 13] = endi;
    faceArray[faceoffset + 14] = t3b;

    // sides
    lineArray[lineoffset] = t1;
    lineArray[lineoffset + 1] = t2;
    lineArray[lineoffset + 2] = t1;
    lineArray[lineoffset + 3] = t4;
    // lineArray[lineoffset+4] = t2, lineArray[lineoffset+5] = t3;
    lineArray[lineoffset + 4] = t3;
    lineArray[lineoffset + 5] = t4;
    // caps
    lineArray[lineoffset + 6] = t1;
    lineArray[lineoffset + 7] = t4;
    // lineArray[lineoffset+10] = t1, lineArray[lineoffset+11] = fromi;
    // lineArray[lineoffset+12] = t4, lineArray[lineoffset+13] = fromi;

    lineArray[lineoffset + 8] = t2b;
    lineArray[lineoffset + 9] = t2; // toi
    lineArray[lineoffset + 10] = t2b;
    lineArray[lineoffset + 11] = t3b;
    lineArray[lineoffset + 12] = t3;
    lineArray[lineoffset + 13] = t3b; // toi
    // arrowhead
    lineArray[lineoffset + 14] = t2b;
    lineArray[lineoffset + 15] = endi;
    lineArray[lineoffset + 16] = t2b;
    lineArray[lineoffset + 17] = t3b;
    lineArray[lineoffset + 18] = endi;
    lineArray[lineoffset + 19] = t3b;

    geoGroup.faceidx += 15;
    geoGroup.lineidx += 20;
  };

  // Update a bounding sphere's position and radius
  // from list of centroids and new points
  /**
   * @param {$3Dmol.Sphere}
   *            sphere
   * @param {Object}
   *            components, centroid of all objects in shape
   * @param {Array}
   *            points, flat array of all points in shape
   * @param {int} numPoints, number of valid poitns in points
   */
  const updateBoundingFromPoints = function updateBoundingFromPoints(sphere, components, points, numPoints) {
    sphere.center.set(0, 0, 0);

    // previously I weighted each component's center equally, but I think
    // it is better to use all points
    let xmin = Infinity;
    let ymin = Infinity;
    let zmin = Infinity;
    let xmax = -Infinity;
    let ymax = -Infinity;
    let zmax = -Infinity;
    if (sphere.box) {
      xmin = sphere.box.min.x;
      xmax = sphere.box.max.x;
      ymin = sphere.box.min.y;
      ymax = sphere.box.max.y;
      zmin = sphere.box.min.z;
      zmax = sphere.box.max.z;
    }

    for (let i = 0, il = numPoints; i < il; i++) {
      const x = points[i * 3];
      const y = points[i * 3 + 1];
      const z = points[i * 3 + 2];
      if (x < xmin) xmin = x;
      if (y < ymin) ymin = y;
      if (z < zmin) zmin = z;
      if (x > xmax) xmax = x;
      if (y > ymax) ymax = y;
      if (z > zmax) zmax = z;
    }

    sphere.center.set((xmax + xmin) / 2, (ymax + ymin) / 2, (zmax + zmin) / 2);
    sphere.radius = sphere.center.distanceTo({x: xmax, y: ymax, z: zmax});
    sphere.box = {min: {x: xmin, y: ymin, z: zmin}, max: {x: xmax, y: ymax, z: zmax}};
  };

  // helper function for adding an appropriately sized mesh
  const addCustomGeo = function addCustomGeo(shape, geo, mesh, color, clickable) {
    const geoGroup = geo.addGeoGroup();
    const {vertexArr} = mesh;
    const {normalArr} = mesh;
    const {faceArr} = mesh;

    geoGroup.vertices = vertexArr.length;
    geoGroup.faceidx = faceArr.length;

    let offset;
    let v;
    let a;
    let b;
    let c;
    let i;
    let il;
    let r;
    let g;
    const {vertexArray} = geoGroup;
    const {colorArray} = geoGroup;

    if (color.constructor !== Array) {
      r = color.r;
      g = color.g;
      b = color.b;
    }
    for (i = 0, il = geoGroup.vertices; i < il; ++i) {
      offset = i * 3;
      v = vertexArr[i];
      vertexArray[offset] = v.x;
      vertexArray[offset + 1] = v.y;
      vertexArray[offset + 2] = v.z;

      if (color.constructor === Array) {
        c = color[i];
        r = c.r;
        g = c.g;
        b = c.b;
      }

      colorArray[offset] = r;
      colorArray[offset + 1] = g;
      colorArray[offset + 2] = b;
    }

    if (clickable) {
      for (i = 0, il = geoGroup.faceidx / 3; i < il; ++i) {
        offset = i * 3;
        a = faceArr[offset];
        b = faceArr[offset + 1];
        c = faceArr[offset + 2];
        const vA = new $3Dmol.Vector3();
        const vB = new $3Dmol.Vector3();
        const vC = new $3Dmol.Vector3();
        shape.intersectionShape.triangle.push(
          new $3Dmol.Triangle(vA.copy(vertexArr[a]), vB.copy(vertexArr[b]), vC.copy(vertexArr[c]))
        );
      }
    }

    if (clickable) {
      const center = new $3Dmol.Vector3(0, 0, 0);
      let cnt = 0;
      for (let g = 0; g < geo.geometryGroups.length; g++) {
        center.add(geo.geometryGroups[g].getCentroid());
        cnt++;
      }
      center.divideScalar(cnt);

      updateBoundingFromPoints(
        shape.boundingSphere,
        {centroid: center},
        vertexArray,
        geoGroup.vertices
      );
    }

    geoGroup.faceArray = new Uint16Array(faceArr);

    geoGroup.truncateArrayBuffers(true, true);

    if (normalArr.length < geoGroup.vertices) geoGroup.setNormals();
    else {
      const normalArray = (geoGroup.normalArray = new Float32Array(geoGroup.vertices * 3));
      let n;
      for (i = 0, il = geoGroup.vertices; i < il; ++i) {
        offset = i * 3;
        n = normalArr[i];
        normalArray[offset] = n.x;
        normalArray[offset + 1] = n.y;
        normalArray[offset + 2] = n.z;
      }
    }

    geoGroup.setLineIndices();
    geoGroup.lineidx = geoGroup.lineArray.length;
  };

  // handles custom shape generation from user supplied arrays
  // May need to generate normal and/or line indices
  /**
   * @param {$3Dmol.GLShape}
   *            shape
   * @param {geometry}
   *            geo
   * @param {CustomSpec}
   *            customSpec
   */
  const drawCustom = function drawCustom(shape, geo, customSpec) {
    const mesh = customSpec;
    const {vertexArr} = mesh;
    const {faceArr} = mesh;
    if (vertexArr.length === 0 || faceArr.length === 0) {
      console.warn(
        'Error adding custom shape component: No vertices and/or face indices supplied!'
      );
    }

    let {color} = customSpec;
    if (typeof color == 'undefined') {
      color = shape.color;
    }
    color = $3Dmol.CC.color(color);

    // var firstgeo = geo.geometryGroups.length;
    const splits = $3Dmol.splitMesh(mesh);
    for (let i = 0, n = splits.length; i < n; i++) {
      addCustomGeo(
        shape,
        geo,
        splits[i],
        splits[i].colorArr ? splits[i].colorArr : color,
        customSpec.clickable
      );
    }
  };

  /**
   *
   * @param {$3Dmol.GLShape}
   *            shape
   * @param {ShapeSpec}
   *            stylespec
   * @returns {undefined}
   */
  const updateFromStyle = function updateFromStyle(shape, stylespec) {
    if (typeof stylespec.color != 'undefined') {
      shape.color = stylespec.color || new $3Dmol.Color();
      if (!(stylespec.color instanceof $3Dmol.Color))
        shape.color = $3Dmol.CC.color(stylespec.color);
    } else {
      shape.color = $3Dmol.CC.color(0);
    }
    shape.wireframe = !!stylespec.wireframe;
    // opacity is the preferred nomenclature, support alpha for backwards compat
    shape.opacity = stylespec.alpha ? $3Dmol.Math.clamp(stylespec.alpha, 0.0, 1.0) : 1.0;
    if (typeof stylespec.opacity != 'undefined') {
      shape.opacity = $3Dmol.Math.clamp(stylespec.opacity, 0.0, 1.0);
    }
    shape.side = stylespec.side !== undefined ? stylespec.side : $3Dmol.DoubleSide;

    shape.linewidth = typeof stylespec.linewidth == 'undefined' ? 1 : stylespec.linewidth;
    // Click handling
    shape.clickable = !!stylespec.clickable;
    shape.callback = $3Dmol.makeFunction(stylespec.callback);
    shape.hoverable = !!stylespec.hoverable;
    shape.hover_callback = $3Dmol.makeFunction(stylespec.hover_callback);
    shape.unhover_callback = $3Dmol.makeFunction(stylespec.unhover_callback);

    shape.hidden = stylespec.hidden;
    shape.frame = stylespec.frame;
  };

  /**
   * Custom renderable shape
   *
   * @constructor $3Dmol.GLShape
   *
   * @param {Object}
   *            stylespec
   * @returns {$3Dmol.GLShape}
   */
  function GLShape(stylespec) {
    stylespec = stylespec || {};
    $3Dmol.ShapeIDCount++;

    this.boundingSphere = new $3Dmol.Sphere();
    /** @type {IntersectionShapes} */
    this.intersectionShape = {
      sphere: [],
      cylinder: [],
      line: [],
      triangle: [],
    };

    updateFromStyle(this, stylespec);

    // Keep track of shape components and their centroids
    const components = [];
    let shapeObj = null;
    let renderedShapeObj = null;

    const geo = new $3Dmol.Geometry(true);
    const linegeo = new $3Dmol.Geometry(true);

    /** Update shape with new style specification
         * @function $3Dmol.GLShape#updateStyle
         * @param {ShapeSpec} newspec
         * @return {$3Dmol.GLShape}
           @example 
            let sphere = viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
            sphere.updateStyle({color:'yellow',opacity:0.5});
            viewer.render();
         */
    this.updateStyle = function (newspec) {
      for (const prop in newspec) {
        stylespec[prop] = newspec[prop];
      }

      updateFromStyle(this, stylespec);

      if (newspec.voldata && newspec.volscheme) {
        $3Dmol.adjustVolumeStyle(newspec);

        // convert volumetric data into colors
        const scheme = newspec.volscheme;
        const {voldata} = newspec;
        const cc = $3Dmol.CC;
        const range = scheme.range() || [-1, 1];
        geo.setColors((x, y, z) => {
          const val = voldata.getVal(x, y, z);
          const col = cc.color(scheme.valueToHex(val, range));
          return col;
        });
        delete this.color;
      }
    };

    /**
         * Creates a custom shape from supplied vertex and face arrays
         * @function $3Dmol.GLShape#addCustom
         * @param {CustomShapeSpec} customSpec
         * @return {$3Dmol.GLShape}
         
         */
    this.addCustom = function (customSpec) {
      customSpec.vertexArr = customSpec.vertexArr || [];
      customSpec.faceArr = customSpec.faceArr || [];
      customSpec.normalArr = customSpec.normalArr || [];

      // will split mesh as needed
      drawCustom(this, geo, customSpec);
    };

    /**
         * Creates a sphere shape
         * @function $3Dmol.GLShape#addSphere
         * @param {SphereSpec} sphereSpec
         * @return {$3Dmol.GLShape}
         @example 
         viewer.addSphere({center:{x:0,y:0,z:0},radius:10.0,color:'red'});
         
         viewer.render();
         */
    this.addSphere = function (sphereSpec) {
      sphereSpec.center = sphereSpec.center || {
        x: 0,
        y: 0,
        z: 0,
      };
      sphereSpec.radius = sphereSpec.radius
        ? $3Dmol.Math.clamp(sphereSpec.radius, 0, Infinity)
        : 1.5;
      sphereSpec.color = $3Dmol.CC.color(sphereSpec.color);

      this.intersectionShape.sphere.push(new $3Dmol.Sphere(sphereSpec.center, sphereSpec.radius));

      $3Dmol.GLDraw.drawSphere(geo, sphereSpec.center, sphereSpec.radius, sphereSpec.color);

      components.push({
        centroid: new $3Dmol.Vector3(sphereSpec.center.x, sphereSpec.center.y, sphereSpec.center.z),
      });
      const geoGroup = geo.updateGeoGroup(0);

      updateBoundingFromPoints(
        this.boundingSphere,
        components,
        geoGroup.vertexArray,
        geoGroup.vertices
      );
    };

    /**
         * Creates a box
         * @function $3Dmol.GLShape#addBox
         * @param {BoxSpec} boxSpec
         * @return {$3Dmol.GLShape}
         @example 
         var shape = viewer.addShape({color:'red'});
         shape.addBox({corner: {x:1,y:2,z:0}, dimensions: {w: 4, h: 2, d: 6}});
         shape.addBox({corner: {x:-5,y:-3,z:0},
                       dimensions: { w: {x:1,y:1,z:0},
                                     h: {x:-1,y:1,z:0},
                                     d: {x:0,y:0,z:1} }});
         viewer.zoomTo();
         viewer.rotate(30);
         viewer.render();
         */
    this.addBox = function (boxSpec) {
      const dim = boxSpec.dimensions || {w: 1, h: 1, d: 1};

      // dimensions may be scalar or vector quantities
      let {w} = dim;
      if (typeof dim.w == 'number') {
        w = {x: dim.w, y: 0, z: 0};
      }
      let {h} = dim;
      if (typeof dim.h == 'number') {
        h = {x: 0, y: dim.h, z: 0};
      }
      let {d} = dim;
      if (typeof dim.d == 'number') {
        d = {x: 0, y: 0, z: dim.d};
      }

      // can position using corner OR center
      let c = boxSpec.corner;
      if (c == undefined) {
        if (boxSpec.center !== undefined) {
          c = {
            x: boxSpec.center.x - 0.5 * (w.x + h.x + d.x),
            y: boxSpec.center.y - 0.5 * (w.y + h.y + d.y),
            z: boxSpec.center.z - 0.5 * (w.z + h.z + d.z),
          };
        } else {
          // default to origin
          c = {x: 0, y: 0, z: 0};
        }
      }

      // 8 vertices
      const uv = [
        {x: c.x, y: c.y, z: c.z},
        {x: c.x + w.x, y: c.y + w.y, z: c.z + w.z},
        {x: c.x + h.x, y: c.y + h.y, z: c.z + h.z},
        {x: c.x + w.x + h.x, y: c.y + w.y + h.y, z: c.z + w.z + h.z},
        {x: c.x + d.x, y: c.y + d.y, z: c.z + d.z},
        {x: c.x + w.x + d.x, y: c.y + w.y + d.y, z: c.z + w.z + d.z},
        {x: c.x + h.x + d.x, y: c.y + h.y + d.y, z: c.z + h.z + d.z},
        {x: c.x + w.x + h.x + d.x, y: c.y + w.y + h.y + d.y, z: c.z + w.z + h.z + d.z},
      ];

      // but.. so that we can have sharp issues, we want a unique normal
      // for each face - since normals are associated with vertices, need to duplicate

      // bottom
      // 0 1
      // 2 3
      // top
      // 4 5
      // 6 7
      const verts = [];
      const faces = [];
      // bottom
      verts.splice(verts.length, 0, uv[0], uv[1], uv[2], uv[3]);
      faces.splice(faces.length, 0, 0, 2, 1, 1, 2, 3);
      let foff = 4;
      // front
      verts.splice(verts.length, 0, uv[2], uv[3], uv[6], uv[7]);
      faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
      foff += 4;
      // back
      verts.splice(verts.length, 0, uv[4], uv[5], uv[0], uv[1]);
      faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
      foff += 4;
      // top
      verts.splice(verts.length, 0, uv[6], uv[7], uv[4], uv[5]);
      faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
      foff += 4;
      // right
      verts.splice(verts.length, 0, uv[3], uv[1], uv[7], uv[5]);
      faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
      foff += 4;
      // left
      verts.splice(verts.length, 0, uv[2], uv[6], uv[0], uv[4]); // fix: was 2 0 6 4 , was flipped! will this ruin anything?
      // and is this the reason for having double sided lambert shading? the box had a flipped face
      faces.splice(faces.length, 0, foff + 0, foff + 2, foff + 1, foff + 1, foff + 2, foff + 3);
      foff += 4;

      const spec = $3Dmol.extend({}, boxSpec);
      spec.vertexArr = verts;
      spec.faceArr = faces;
      spec.normalArr = [];
      drawCustom(this, geo, spec);

      const centroid = new $3Dmol.Vector3();
      components.push({
        centroid: centroid.addVectors(uv[0], uv[7]).multiplyScalar(0.5),
      });
      const geoGroup = geo.updateGeoGroup(0);
      updateBoundingFromPoints(
        this.boundingSphere,
        components,
        geoGroup.vertexArray,
        geoGroup.vertices
      );
    };

    /**
         * Creates a cylinder shape
         * @function $3Dmol.GLShape#addCylinder
         * @param {CylinderSpec} cylinderSpec
         * @return {$3Dmol.GLShape}
         @example
              viewer.addCylinder({start:{x:0.0,y:0.0,z:0.0},
                                  end:{x:10.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  fromCap:1,
                                  toCap:2,
                                  color:'red',
                                  hoverable:true,
                                  clickable:true,
                                  callback:function(){ this.color.setHex(0x00FFFF00);viewer.render( );},
                                  hover_callback: function(){ viewer.render( );},
                                  unhover_callback: function(){ this.color.setHex(0xFF000000);viewer.render( );}
                                 });
              viewer.addCylinder({start:{x:0.0,y:2.0,z:0.0},
                                  end:{x:0.0,y:10.0,z:0.0},
                                  radius:0.5,
                                  fromCap:false,
                                  toCap:true,
                                  color:'teal'});
              viewer.addCylinder({start:{x:15.0,y:0.0,z:0.0},
                                  end:{x:20.0,y:0.0,z:0.0},
                                  radius:1.0,
                                  color:'black',
                                  fromCap:false,
                                  toCap:false});
              viewer.render();
         */
    this.addCylinder = function (cylinderSpec) {
      cylinderSpec.start = cylinderSpec.start || {};
      cylinderSpec.end = cylinderSpec.end || {};

      const start = new $3Dmol.Vector3(
        cylinderSpec.start.x || 0,
        cylinderSpec.start.y || 0,
        cylinderSpec.start.z || 0
      );
      const end = new $3Dmol.Vector3(
        cylinderSpec.end.x,
        cylinderSpec.end.y || 0,
        cylinderSpec.end.z || 0
      );
      if (typeof end.x == 'undefined') end.x = 3; // show something even if undefined

      const radius = cylinderSpec.radius || 0.1;
      const color = $3Dmol.CC.color(cylinderSpec.color);

      this.intersectionShape.cylinder.push(new $3Dmol.Cylinder(start, end, radius));

      $3Dmol.GLDraw.drawCylinder(
        geo,
        start,
        end,
        radius,
        color,
        cylinderSpec.fromCap,
        cylinderSpec.toCap
      );

      const centroid = new $3Dmol.Vector3();
      components.push({
        centroid: centroid.addVectors(start, end).multiplyScalar(0.5),
      });
      const geoGroup = geo.updateGeoGroup(0);
      updateBoundingFromPoints(
        this.boundingSphere,
        components,
        geoGroup.vertexArray,
        geoGroup.vertices
      );
    };

    /**
     * Creates a dashed cylinder shape
     * @function $3Dmol.GLShape#addDashedCylinder
     * @param {CylinderSpec} cylinderSpec
     * @return {$3Dmol.GLShape}
     */
    this.addDashedCylinder = function (cylinderSpec) {
      cylinderSpec.start = cylinderSpec.start || {};
      cylinderSpec.end = cylinderSpec.end || {};
      cylinderSpec.dashLength = cylinderSpec.dashLength || 0.25;
      cylinderSpec.gapLength = cylinderSpec.gapLength || 0.25;

      const start = new $3Dmol.Vector3(
        cylinderSpec.start.x || 0,
        cylinderSpec.start.y || 0,
        cylinderSpec.start.z || 0
      );
      const end = new $3Dmol.Vector3(
        cylinderSpec.end.x,
        cylinderSpec.end.y || 0,
        cylinderSpec.end.z || 0
      );
      if (typeof end.x == 'undefined') end.x = 3; // show something even if undefined

      const radius = cylinderSpec.radius || 0.1;
      const color = $3Dmol.CC.color(cylinderSpec.color);

      const cylinderLength = Math.sqrt(
        (start.x - end.x) ** 2 + (start.y - end.y) ** 2 + (start.z - end.z) ** 2
      );

      const count = cylinderLength / (cylinderSpec.gapLength + cylinderSpec.dashLength);

      let new_start = new $3Dmol.Vector3(
        cylinderSpec.start.x || 0,
        cylinderSpec.start.y || 0,
        cylinderSpec.start.z || 0
      );
      let new_end = new $3Dmol.Vector3(
        cylinderSpec.end.x,
        cylinderSpec.end.y || 0,
        cylinderSpec.end.z || 0
      );

      const gapVector = new $3Dmol.Vector3(
        (end.x - start.x) / (cylinderLength / cylinderSpec.gapLength),
        (end.y - start.y) / (cylinderLength / cylinderSpec.gapLength),
        (end.z - start.z) / (cylinderLength / cylinderSpec.gapLength)
      );
      const dashVector = new $3Dmol.Vector3(
        (end.x - start.x) / (cylinderLength / cylinderSpec.dashLength),
        (end.y - start.y) / (cylinderLength / cylinderSpec.dashLength),
        (end.z - start.z) / (cylinderLength / cylinderSpec.dashLength)
      );

      for (let place = 0; place < count; place++) {
        new_end = new $3Dmol.Vector3(
          new_start.x + dashVector.x,
          new_start.y + dashVector.y,
          new_start.z + dashVector.z
        );

        this.intersectionShape.cylinder.push(new $3Dmol.Cylinder(new_start, new_end, radius));

        $3Dmol.GLDraw.drawCylinder(
          geo,
          new_start,
          new_end,
          radius,
          color,
          cylinderSpec.fromCap,
          cylinderSpec.toCap
        );

        new_start = new $3Dmol.Vector3(
          new_end.x + gapVector.x,
          new_end.y + gapVector.y,
          new_end.z + gapVector.z
        );
      }
      const centroid = new $3Dmol.Vector3();
      components.push({
        centroid: centroid.addVectors(start, end).multiplyScalar(0.5),
      });
      const geoGroup = geo.updateGeoGroup(0);
      updateBoundingFromPoints(
        this.boundingSphere,
        components,
        geoGroup.vertexArray,
        geoGroup.vertices
      );
    };

    /**
     * Creates a curved shape
     * @function $3Dmol.GLShape#addCurve
     * @param {CurveSpec} curveSpec
     * @return {$3Dmol.GLShape}
     */
    this.addCurve = function (curveSpec) {
      curveSpec.points = curveSpec.points || [];
      curveSpec.smooth = curveSpec.smooth || 10;
      if (typeof curveSpec.fromCap == 'undefined') curveSpec.fromCap = 2;
      if (typeof curveSpec.toCap == 'undefined') curveSpec.toCap = 2;

      // subdivide into smoothed spline points
      const points = $3Dmol.subdivide_spline(curveSpec.points, curveSpec.smooth);

      if (points.length < 3) {
        console.log('Too few points in addCurve');
        return;
      }

      const radius = curveSpec.radius || 0.1;
      const color = $3Dmol.CC.color(curveSpec.color);
      // TODO TODO - this is very inefficient, should create our
      // own water tight model with proper normals...

      // if arrows are requested, peel off enough points to fit
      // at least 2*r of arrowness
      let start = 0;
      let end = points.length - 1;
      const segmentlen = points[0].distanceTo(points[1]);
      const npts = Math.ceil((2 * radius) / segmentlen);
      if (curveSpec.toArrow) {
        end -= npts;
        const arrowspec = {
          start: points[end],
          end: points[points.length - 1],
          radius,
          color,
          mid: 0.0001,
        };
        this.addArrow(arrowspec);
      }
      if (curveSpec.fromArrow) {
        start += npts;
        const arrowspec = {
          start: points[start],
          end: points[0],
          radius,
          color,
          mid: 0.0001,
        };
        this.addArrow(arrowspec);
      }

      const midway = Math.ceil(points.length / 2);
      const middleSpec = {radius, color, fromCap: 2, toCap: 2};
      for (let i = start; i < end; i++) {
        middleSpec.start = points[i];
        middleSpec.end = points[i + 1];
        middleSpec.fromCap = 2;
        middleSpec.toCap = 2;
        if (i < midway) {
          middleSpec.fromCap = 2;
          middleSpec.toCap = 0;
        } else if (i > midway) {
          middleSpec.fromCap = 0;
          middleSpec.toCap = 2;
        } else {
          middleSpec.fromCap = 2;
          middleSpec.toCap = 2;
        }

        this.addCylinder(middleSpec);
      }
    };

    /**
         * Creates a line shape
         * @function $3Dmol.GLShape#addLine         
         * @param {LineSpec} lineSpec
         * @return {$3Dmol.GLShape}
         @example
         $3Dmol.download("pdb:2ABJ",viewer,{},function(){
                  viewer.addLine({dashed:true,start:{x:0,y:0,z:0},end:{x:100,y:100,z:100}});
                  viewer.render(callback);
              });

         */
    this.addLine = function (lineSpec) {
      lineSpec.start = lineSpec.start || {};
      lineSpec.end = lineSpec.end || {};

      const start = new $3Dmol.Vector3(
        lineSpec.start.x || 0,
        lineSpec.start.y || 0,
        lineSpec.start.z || 0
      );
      const end = new $3Dmol.Vector3(lineSpec.end.x, lineSpec.end.y || 0, lineSpec.end.z || 0);
      if (typeof end.x == 'undefined') end.x = 3; // show something even if undefined

      let geoGroup = geo.updateGeoGroup(2);

      // make line from start to end
      // for consistency with rest of shapes, uses vertices and lines rather
      // than a separate line geometry
      const vstart = geoGroup.vertices;
      const i = vstart * 3;
      const {vertexArray} = geoGroup;
      vertexArray[i] = start.x;
      vertexArray[i + 1] = start.y;
      vertexArray[i + 2] = start.z;
      vertexArray[i + 3] = end.x;
      vertexArray[i + 4] = end.y;
      vertexArray[i + 5] = end.z;
      geoGroup.vertices += 2;

      const {lineArray} = geoGroup;
      const li = geoGroup.lineidx;
      lineArray[li] = vstart;
      lineArray[li + 1] = vstart + 1;
      geoGroup.lineidx += 2;

      const centroid = new $3Dmol.Vector3();
      components.push({
        centroid: centroid.addVectors(start, end).multiplyScalar(0.5),
      });
      geoGroup = geo.updateGeoGroup(0);
      updateBoundingFromPoints(
        this.boundingSphere,
        components,
        geoGroup.vertexArray,
        geoGroup.vertices
      );
    };

    /**
         * Creates an arrow shape
         * @function $3Dmol.GLShape#addArrow        
         * @param {ArrowSpec} arrowSpec
         * @return {$3Dmol.GLShape}
         @example
          $3Dmol.download("pdb:4DM7",viewer,{},function(){
                  viewer.setBackgroundColor(0xffffffff);
                  viewer.addArrow({
                      start: {x:-10.0, y:0.0, z:0.0},
                      end: {x:0.0, y:-10.0, z:0.0},
                      radius: 1.0,
                      radiusRadio:1.0,
                      mid:1.0,
                      clickable:true,
                      callback:function(){
                          this.color.setHex(0xFF0000FF);
                          viewer.render( );
                      }
                  });
                  viewer.render();
                });
         */
    this.addArrow = function (arrowSpec) {
      arrowSpec.start = arrowSpec.start || {};
      arrowSpec.end = arrowSpec.end || {};

      arrowSpec.start = new $3Dmol.Vector3(
        arrowSpec.start.x || 0,
        arrowSpec.start.y || 0,
        arrowSpec.start.z || 0
      );

      if (arrowSpec.dir instanceof $3Dmol.Vector3 && typeof arrowSpec.length === 'number') {
        const end = arrowSpec.dir.clone().multiplyScalar(arrowSpec.length).add(arrowSpec.start);
        arrowSpec.end = end;
      } else {
        arrowSpec.end = new $3Dmol.Vector3(
          arrowSpec.end.x,
          arrowSpec.end.y || 0,
          arrowSpec.end.z || 0
        );
        if (typeof arrowSpec.end.x == 'undefined') arrowSpec.end.x = 3; // show something even if undefined
      }

      arrowSpec.radius = arrowSpec.radius || 0.1;

      arrowSpec.radiusRatio = arrowSpec.radiusRatio || 1.618034;

      arrowSpec.mid = arrowSpec.mid > 0 && arrowSpec.mid < 1 ? arrowSpec.mid : 0.618034;

      drawArrow(this, geo, arrowSpec);

      const centroid = new $3Dmol.Vector3();
      components.push({
        centroid: centroid.addVectors(arrowSpec.start, arrowSpec.end).multiplyScalar(0.5),
      });
      const geoGroup = geo.updateGeoGroup(0);
      updateBoundingFromPoints(
        this.boundingSphere,
        components,
        geoGroup.vertexArray,
        geoGroup.vertices
      );
    };

    const distance_from = function (c1, c2) {
      return Math.sqrt((c1.x - c2.x) ** 2 + (c1.y - c2.y) ** 2 + (c1.z - c2.z) ** 2);
    };

    const inSelectedRegion = function (coordinate, selectedRegion, radius) {
      for (let i = 0; i < selectedRegion.length; i++) {
        if (distance_from(selectedRegion[i], coordinate) <= radius) return true;
      }
      return false;
    };

    /**
         * Create isosurface from voluemetric data.
         * @function $3Dmol.GLShape#addIsosurface         
         * @param {$3Dmol.VolumeData} data - volumetric input data
         * @param {IsoSurfaceSpec} isoSpec - volumetric data shape specification
         * @example //the user can specify a selected region for the isosurface 
         $.get('../test_structs/benzene-homo.cube', function(data){
                  var voldata = new $3Dmol.VolumeData(data, "cube");
                  viewer.addIsosurface(voldata, {isoval: 0.01,
                                                 color: "blue",
                                                 alpha: 0.5,
                                                 smoothness: 10});
                  viewer.addIsosurface(voldata, {isoval: -0.01,
                                                 color: "red",
                                                 smoothness: 5,
                                                 opacity:0.5,
                                                 wireframe:true,
                                                 clickable:true,
                                                 callback:
                                                 function() {
                                                     this.opacity = 0.0;
                                                     viewer.render( );
                                                 }});
                  viewer.setStyle({}, {stick:{}});
                  viewer.zoomTo();
                  viewer.render();
                });
         */
    this.addIsosurface = function (data, volSpec, callback) {
      // may want to cache the arrays geneerated when selectedRegion ==true

      const isoval =
        volSpec.isoval !== undefined && typeof volSpec.isoval === 'number' ? volSpec.isoval : 0.0;
      const voxel = !!volSpec.voxel;
      const smoothness = volSpec.smoothness === undefined ? 1 : volSpec.smoothness;

      const nX = data.size.x;
      const nY = data.size.y;
      const nZ = data.size.z;
      const vertnums = new Int16Array(nX * nY * nZ);
      const vals = data.data;

      let i;
      let il;

      for (i = 0, il = vertnums.length; i < il; ++i) vertnums[i] = -1;

      const bitdata = new Uint8Array(nX * nY * nZ);

      // mark locations partitioned by isoval
      for (i = 0, il = vals.length; i < il; ++i) {
        const val = isoval >= 0 ? vals[i] - isoval : isoval - vals[i];
        if (val > 0) bitdata[i] |= ISDONE;
      }

      let verts = [];
      let faces = [];

      $3Dmol.MarchingCube.march(bitdata, verts, faces, {
        fulltable: true,
        voxel,
        unitCube: data.unit,
        origin: data.origin,
        matrix: data.matrix,
        nX,
        nY,
        nZ,
      });

      if (!voxel && smoothness > 0) $3Dmol.MarchingCube.laplacianSmooth(smoothness, verts, faces);
      const vertexmapping = [];
      const newvertices = [];
      const newfaces = [];

      if (volSpec.selectedRegion && volSpec.coords === undefined) {
        volSpec.coords = volSpec.selectedRegion; // backwards compat for incorrectly documented feature
      }
      if (volSpec.coords !== undefined) {
        let xmax = volSpec.coords[0].x;
        let ymax = volSpec.coords[0].y;
        let zmax = volSpec.coords[0].z;
        let xmin = volSpec.coords[0].x;
        let ymin = volSpec.coords[0].y;
        let zmin = volSpec.coords[0].z;

        for (let i = 0; i < volSpec.coords.length; i++) {
          if (volSpec.coords[i].x > xmax) xmax = volSpec.coords[i].x;
          else if (volSpec.coords[i].x < xmin) xmin = volSpec.coords[i].x;
          if (volSpec.coords[i].y > ymax) ymax = volSpec.coords[i].y;
          else if (volSpec.coords[i].y < ymin) ymin = volSpec.coords[i].y;
          if (volSpec.coords[i].z > zmax) zmax = volSpec.coords[i].z;
          else if (volSpec.coords[i].z < zmin) zmin = volSpec.coords[i].z;
        }

        let rad = 2;
        if (volSpec.radius !== undefined) {
          rad = volSpec.radius; // backwards compat
        }
        if (volSpec.selectedOffset !== undefined) {
          // backwards compat
          rad = volSpec.selectedOffset;
        }
        if (volSpec.seldist !== undefined) {
          rad = volSpec.seldist;
        }

        xmin -= rad;
        xmax += rad;
        ymin -= rad;
        ymax += rad;
        zmin -= rad;
        zmax += rad;

        // accounts for radius
        for (let i = 0; i < verts.length; i++) {
          if (
            verts[i].x > xmin &&
            verts[i].x < xmax &&
            verts[i].y > ymin &&
            verts[i].y < ymax &&
            verts[i].z > zmin &&
            verts[i].z < zmax &&
            inSelectedRegion(verts[i], volSpec.coords, rad)
          ) {
            vertexmapping.push(newvertices.length);
            newvertices.push(verts[i]);
          } else {
            vertexmapping.push(-1);
          }
        }
        for (let i = 0; i + 2 < faces.length; i += 3) {
          if (
            vertexmapping[faces[i]] !== -1 &&
            vertexmapping[faces[i + 1]] !== -1 &&
            vertexmapping[faces[i + 2]] !== -1
          ) {
            newfaces.push(faces[i] - (faces[i] - vertexmapping[faces[i]]));
            newfaces.push(faces[i + 1] - (faces[i + 1] - vertexmapping[faces[i + 1]]));
            newfaces.push(faces[i + 2] - (faces[i + 2] - vertexmapping[faces[i + 2]]));
          }
        }
        verts = newvertices;
        faces = newfaces;
      }

      drawCustom(this, geo, {
        vertexArr: verts,
        faceArr: faces,
        normalArr: [],
        clickable: volSpec.clickable,
        hoverable: volSpec.hoverable,
      });

      this.updateStyle(volSpec);

      // computing bounding sphere from vertices
      const origin = new $3Dmol.Vector3(data.origin.x, data.origin.y, data.origin.z);
      const size = new $3Dmol.Vector3(
        data.size.x * data.unit.x,
        data.size.y * data.unit.y,
        data.size.z * data.unit.z
      );

      const total = new $3Dmol.Vector3(0, 0, 0);
      const maxv = origin.clone();
      const minv = origin.clone().add(size);
      for (let i = 0; i < verts.length; i++) {
        total.add(verts[i]);
        maxv.max(verts[i]);
        minv.min(verts[i]);
      }
      total.divideScalar(verts.length);
      const len1 = total.distanceTo(minv);
      const len2 = total.distanceTo(maxv);
      this.boundingSphere.center = total;
      this.boundingSphere.radius = Math.max(len1, len2);
      if (typeof callback == 'function') callback();
    };

    /**
     * @deprecated Use addIsosurface instead
     * Creates custom shape from volumetric data
     * @param {string} data - Volumetric input data
     * @param {string} fmt - Input data format (e.g. 'cube' for cube file format)
     * @param {IsoSurfaceSpec} isoSpec - Volumetric data shape specification
     * @return {$3Dmol.GLShape}
     */
    this.addVolumetricData = function (data, fmt, volSpec) {
      data = new $3Dmol.VolumeData(data, fmt);
      this.addIsosurface(data, volSpec);
    };

    // for internal use, truncate buffers to save memory
    this.finalize = function () {
      finalizeGeo(geo);
      geo.initTypedArrays();
      return geo;
    };

    /**
     * Initialize webgl objects for rendering
     * @param {$3Dmol.Object3D} group
     *
     */
    this.globj = function (group) {
      if (renderedShapeObj) {
        group.remove(renderedShapeObj);
        renderedShapeObj = null;
      }

      if (this.hidden) return;
      finalizeGeo(geo);
      geo.initTypedArrays();

      if (this.wireframe) {
        geo.setUpWireframe();
      }

      if (typeof this.color != 'undefined') updateColor(geo, this.color);

      shapeObj = new $3Dmol.Object3D();
      let material = null;

      if (this.side == $3Dmol.DoubleSide) {
        material = new $3Dmol.MeshDoubleLambertMaterial({
          wireframe: this.wireframe,
          side: this.side,
          transparent: this.opacity < 1,
          opacity: this.opacity,
          wireframeLinewidth: this.linewidth,
          vertexColors: $3Dmol.VertexColors,
        });
      } else {
        material = new $3Dmol.MeshLambertMaterial({
          wireframe: this.wireframe,
          side: this.side,
          transparent: this.opacity < 1,
          opacity: this.opacity,
          wireframeLinewidth: this.linewidth,
          vertexColors: $3Dmol.VertexColors,
        });
      }

      const mesh = new $3Dmol.Mesh(geo, material);

      shapeObj.add(mesh);

      const lineMaterial = new $3Dmol.LineBasicMaterial({
        linewidth: this.linewidth,
        color: this.color,
      });
      const line = new $3Dmol.Line(linegeo, lineMaterial, $3Dmol.LinePieces);
      shapeObj.add(line);

      renderedShapeObj = shapeObj.clone();
      group.add(renderedShapeObj);
    };

    this.removegl = function (group) {
      if (renderedShapeObj) {
        // dispose of geos and materials
        if (renderedShapeObj.geometry !== undefined) renderedShapeObj.geometry.dispose();
        if (renderedShapeObj.material !== undefined) renderedShapeObj.material.dispose();
        group.remove(renderedShapeObj);
        renderedShapeObj = null;
      }
      shapeObj = null;
    };
  }

  Object.defineProperty(GLShape.prototype, 'position', {
    get() {
      return this.boundingSphere.center;
    },
  });

  Object.defineProperty(GLShape.prototype, 'x', {
    get() {
      return this.boundingSphere.center.x;
    },
  });

  Object.defineProperty(GLShape.prototype, 'y', {
    get() {
      return this.boundingSphere.center.y;
    },
  });

  Object.defineProperty(GLShape.prototype, 'z', {
    get() {
      return this.boundingSphere.center.z;
    },
  });

  return GLShape;
})();

$3Dmol.ShapeIDCount = 0;

$3Dmol.splitMesh = function splitMesh(mesh) {
  const MAXVERT = 64000; // webgl only supports 2^16 elements, leave a little breathing room (require at least 2)
  // peel off 64k vertices rsvh into their own mesh
  // duplicating vertices and normals as necessary to preserve faces and lines

  if (mesh.vertexArr.length < MAXVERT) return [mesh]; // typical case

  const slices = [{vertexArr: [], normalArr: [], faceArr: []}];
  if (mesh.colorArr) slices.colorArr = [];
  const vertSlice = []; // indexed by original vertex to get current slice
  const vertIndex = []; // indexed by original vertex to get index within slice
  let currentSlice = 0;

  // for each face, make sure all three vertices (or copies) are in the same slice
  const faces = mesh.faceArr;
  for (let i = 0, nf = faces.length; i < nf; i += 3) {
    const slice = slices[currentSlice];
    for (let j = 0; j < 3; j++) {
      // process each vertex to make sure it is assigned a slice
      // all vertices of a face must belong to the same slice
      const v = faces[i + j];
      if (vertSlice[v] !== currentSlice) {
        // true if undefined
        vertSlice[v] = currentSlice;
        vertIndex[v] = slice.vertexArr.length;
        slice.vertexArr.push(mesh.vertexArr[v]);
        if (mesh.normalArr && mesh.normalArr[v]) slice.normalArr.push(mesh.normalArr[v]);
        if (mesh.colorArr && mesh.colorArr[v]) slice.colorArr.push(mesh.colorArr[v]);
      }
      slice.faceArr.push(vertIndex[v]);
    }

    if (slice.vertexArr.length >= MAXVERT) {
      // new slice
      slices.push({vertexArr: [], normalArr: [], faceArr: []});
      if (mesh.colorArr) slices.colorArr = [];
      currentSlice++;
    }
  }
  return slices;
};
