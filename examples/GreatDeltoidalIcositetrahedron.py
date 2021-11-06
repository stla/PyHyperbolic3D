# http://dmccooey.com/polyhedra/GreatDeltoidalIcositetrahedron.html
from math import sqrt
import pyvista as pv
import numpy as np
from hyperbolic3d.gyro import gyrotube, gyrotriangle

##~~ Great deltoidal icositetrahedron ~~##

C0 = (4 - sqrt(2)) / 7
C1 = sqrt(2)
vertices = [
  [ 0.0,  0.0,   C1],
  [ 0.0,  0.0,  -C1],
  [  C1,  0.0,  0.0],
  [ -C1,  0.0,  0.0],
  [ 0.0,   C1,  0.0],
  [ 0.0,  -C1,  0.0],
  [-1.0,  0.0, -1.0],
  [-1.0,  0.0,  1.0],
  [ 1.0,  0.0, -1.0],
  [ 1.0,  0.0,  1.0],
  [-1.0, -1.0,  0.0],
  [-1.0,  1.0,  0.0],
  [ 1.0, -1.0,  0.0],
  [ 1.0,  1.0,  0.0],
  [ 0.0, -1.0, -1.0],
  [ 0.0, -1.0,  1.0],
  [ 0.0,  1.0, -1.0],
  [ 0.0,  1.0,  1.0],
  [ -C0,  -C0,  -C0],
  [ -C0,  -C0,   C0],
  [ -C0,   C0,  -C0],
  [ -C0,   C0,   C0],
  [  C0,  -C0,  -C0],
  [  C0,  -C0,   C0],
  [  C0,   C0,  -C0],
  [  C0,   C0,   C0]
]
faces = [
  [  0,  6, 18, 14 ],
  [  0, 14, 22,  8 ],
  [  0,  8, 24, 16 ],
  [  0, 16, 20,  6 ],
  [  1,  9, 23, 15 ],
  [  1, 15, 19,  7 ],
  [  1,  7, 21, 17 ],
  [  1, 17, 25,  9 ],
  [  2,  7, 19, 10 ],
  [  2, 10, 18,  6 ],
  [  2,  6, 20, 11 ],
  [  2, 11, 21,  7 ],
  [  3,  8, 22, 12 ],
  [  3, 12, 23,  9 ],
  [  3,  9, 25, 13 ],
  [  3, 13, 24,  8 ],
  [  4, 10, 19, 15 ],
  [  4, 15, 23, 12 ],
  [  4, 12, 22, 14 ],
  [  4, 14, 18, 10 ],
  [  5, 11, 20, 16 ],
  [  5, 16, 24, 13 ],
  [  5, 13, 25, 17 ],
  [  5, 17, 21, 11 ]
]


# edges
edges = []
for f in faces:
    edges = edges + [[f[0], f[1]], [f[0], f[2]], [f[1], f[2]]]
edges = np.array(edges, dtype=int)
edges = np.unique(edges, axis=0)


# 
def make_meshes(s, col1, col2):
    meshes = []
    for face in faces:
        pts = [vertices[i] for i in face]
        m1 = gyrotriangle(pts[0], pts[1], pts[2], s, depth=4)
        m2 = gyrotriangle(pts[0], pts[2], pts[3], s, depth=4)
        mesh = m1.merge(m2)
        meshes.append((mesh, col1))
    for idx in edges:
        A = vertices[idx[0]]
        B = vertices[idx[1]]
        hyperbola = gyrotube(A, B, s=s, r=0.01)
        meshes.append((hyperbola, col2))
    for vertex in vertices:
        sphere = pv.Sphere(0.035, center=vertex)
        meshes.append((sphere, col2))
    return meshes

meshes = make_meshes(1, "yellow", "darkmagenta")
pltr = pv.Plotter()
pltr.set_background("#363940")
for mesh in meshes:
    pltr.add_mesh(mesh[0], color=mesh[1], smooth_shading=True, specular=15)
pltr.show()
