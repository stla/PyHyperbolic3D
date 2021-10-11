# -*- coding: utf-8 -*-
import pyvista as pv
import numpy as np
from hyperbolic3d.gyro import gyrotube, gyrotriangle

vertices = [
    [0, 0.618033988749895, 1],
    [0, 0.618033988749895, -1],
    [0, -0.618033988749895, 1],
    [0, -0.618033988749895, -1],
    [0.618033988749895, 1, 0],
    [0.618033988749895, -1, 0],
    [-0.618033988749895, 1, 0],
    [-0.618033988749895, -1, 0],
    [1, 0, 0.618033988749895],
    [-1, 0, 0.618033988749895],
    [1, 0, -0.618033988749895],
    [-1, 0, -0.618033988749895]
]

edges = [
    [0, 2],
    [0, 4],
    [0, 6],
    [0, 8],
    [0, 9],
    [1, 3],
    [1, 4],
    [1, 6],
    [1, 10],
    [1, 11],
    [2, 5],
    [2, 7],
    [2, 8],
    [2, 9],
    [3, 5],
    [3, 7],
    [3, 10],
    [3, 11],
    [4, 6],
    [4, 8],
    [4, 10],
    [5, 7],
    [5, 8],
    [5, 10],
    [6, 9],
    [6, 11],
    [7, 9],
    [7, 11],
    [8, 10],
    [9, 11]
]

faces = [
    [0, 2, 8],
    [0, 8, 4],
    [0, 4, 6],
    [0, 6, 9],
    [0, 9, 2],
    [3, 11, 1],
    [3, 1, 10],
    [3, 10, 5],
    [3, 5, 7],
    [3, 7, 11],
    [8, 2, 5],
    [4, 8, 10],
    [6, 4, 1],
    [9, 6, 11],
    [2, 9, 7],
    [1, 11, 6],
    [10, 1, 4],
    [5, 10, 8],
    [7, 5, 2],
    [11, 7, 9]
]

s_    = np.linspace(0.3, 2, 20)
r     = 0.03
depth = 5

for i,s in enumerate(s_):
    pl = pv.Plotter(window_size=[500,500], off_screen=True)
    pl.camera.zoom(0.3)
    pl.camera.view_angle = 100.0
    for vertex in vertices:
        V = pv.Sphere(2 * r, center=vertex)
        pl.add_mesh(V, smooth_shading=True, color="deeppink", specular=5)
    for edge in edges:
        A = vertices[edge[0]]
        B = vertices[edge[1]]
        AB = gyrotube(A, B, s, r)
        pl.add_mesh(AB, smooth_shading=True, color="darkorange", specular=5)
    for face in faces:
        A = vertices[face[0]]
        B = vertices[face[1]]
        C = vertices[face[2]]
        tmesh = gyrotriangle(A, B, C, s, depth)
        tmesh.point_data["distance"] = np.linalg.norm(tmesh.points, axis=1)
        pl.add_mesh(
            tmesh, smooth_shading=True, specular=50, cmap="turbo",
            show_scalar_bar=False, pbr=False, log_scale=False, flip_scalars=True
        )
    num = str(i) if i>=10 else "0"+str(i)
    pl.show(screenshot="pngs/png"+num+".png")

# magick convert -dispose previous -loop 0 -delay 10 -duplicate 1,-2-1 png*.png icosahedron_colored.gif