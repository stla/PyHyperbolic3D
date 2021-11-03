# PyHyperbolic3D

Python stuff for drawing 3D hyperbolic polyhedra.

![](https://github.com/stla/PyHyperbolic3D/raw/main/examples/icosahedron.png)

![](https://github.com/stla/PyHyperbolic3D/raw/main/examples/icosahedron_slider.gif)

![](https://github.com/stla/PyHyperbolic3D/raw/main/examples/icosahedron_colored.gif)

![](https://github.com/stla/PyHyperbolic3D/raw/main/examples/BarthHyperbolicpolyhedron.gif)


#### `gyrotube(A, B, s, r, npoints=300):`

Tubular hyperbolic segment.

##### Parameters
- **`A,B`** points (lists or arrays)

  The two endpoints of the segment.

- **`s`** positive float

   Curvature parameter.
   
- **`r`** positive float

   Radius of the tube.
   
- **`npoints`** integer

   Number of points along the segment. The default is 300.

##### Returns
A PyVista mesh ready for inclusion in a plotting region.

___

#### `gyrotriangle(A, B, C, s, depth=5, tol=1e-6):`

Hyperbolic triangle.

##### Parameters
- **`A,B,C`** points (lists or arrays)

  The vertices of the triangle.

- **`s`** positive float

   Curvature parameter.
   
- **`depth`** integer

   The number of recursive subdivions. The default is 5.

- **`tol`** small positive float

   The tolerance used to merge duplicated points in the mesh.
The default is 1e-6.

##### Returns
A PyVista mesh ready for inclusion in a plotting region.


