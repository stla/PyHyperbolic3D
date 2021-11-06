# -*- coding: utf-8 -*-
from math import pi, cos, sin
import pyvista as pv
import numpy as np
from scipy.optimize import NonlinearConstraint, BFGS, SR1, minimize, Bounds
from hyperbolic3d.gyro import gyrotube, gyrotriangle, gyroABt


rectangles = [
 [0, 1, 6, 5],
 [1, 2, 7, 6],
 [2, 3, 8, 7],
 [3, 4, 9, 8],
 [4, 0, 5, 9]
]

edges = [
 [0, 1],
 [1, 2],
 [2, 3],
 [3, 4],
 [4, 0],
 [5, 6],
 [6, 7],
 [7, 8],
 [8, 9],
 [9, 5],
 [1, 6],
 [2, 7],
 [3, 8],
 [4, 9],
 [0, 5]
]

vs1 = np.array([[cos(2*i*pi/5), sin(2*i*pi/5), 0.3] for i in [0,2,4,1,3]])
vs2 = np.hstack((vs1[:, 0:2], np.full((5,1), -0.3)))
vertices = np.vstack((vs1, vs2))

pentagramms = np.arange(0, 10).reshape(2, 5, order="C")

# intersections of edges
def tintersect(s):
  pts = vertices[pentagramms[0, :], :]
  A = pts[0, :]
  B = pts[1, :]
  C = pts[3, :]
  D = pts[2, :] 
  f = lambda t: np.linalg.norm(gyroABt(A, B, t, s) - gyroABt(C, D, t, s))
  nonlinear_constraint = NonlinearConstraint(
      f, 0, 0.5, jac="2-point", hess=BFGS()
  )
  bounds = Bounds(0, 0.5)
  res = minimize(f, 0.25, method="trust-constr", jac="2-point", hess=SR1(),
               constraints=[nonlinear_constraint],
               options={'verbose': 1}, bounds=bounds)
  return res.x[0]

# draw 
def drawP(pltr, s, col1, col2, col3):
    for rect in rectangles:
        pts = vertices[rect, :]
        m1 = gyrotriangle(pts[0, :], pts[1, :], pts[2, :], s)
        m2 = gyrotriangle(pts[0, :], pts[2, :], pts[3, :], s)
        mesh = m1.merge(m2)
        pltr.add_mesh(mesh, color = col1)
    t = tintersect(s)
    for i in range(2):
        pts = vertices[pentagramms[i,],]
        P0 = gyroABt(pts[0,], pts[1,], t, s)
        P1 = gyroABt(pts[1,], pts[2,], t, s)
        P2 = gyroABt(pts[2,], pts[3,], t, s)
        P3 = gyroABt(pts[3,], pts[4,], t, s)
        P4 = gyroABt(pts[4,], pts[0,], t, s)
        m0 = gyrotriangle(pts[0,], P0, P2, s)
        m1 = gyrotriangle(pts[1,], P1, P3, s)
        m2 = gyrotriangle(pts[2,], P2, P4, s)
        m3 = gyrotriangle(pts[3,], P3, P0, s)
        m4 = gyrotriangle(pts[4,], P4, P1, s)
        mm0 = gyrotriangle(P0, P2, P4, s)
        mm1 = gyrotriangle(P0, P4, P1, s)
        mm2 = gyrotriangle(P0, P1, P3, s)
        mesh = m0.merge(m1).merge(m2).merge(m3).merge(m4).merge(mm0).merge(mm1).merge(mm2)
        pltr.add_mesh(mesh, color = col2, specular=10)
    for idx in edges:
        A = vertices[idx[0], :]
        B = vertices[idx[1], :]
        edge = gyrotube(A, B, s=s, r=0.025)
        pltr.add_mesh(edge, color = col3, specular=10)
    for vertex in vertices:
        sphere = pv.Sphere(0.035, center=vertex)
        pltr.add_mesh(sphere, color = col3, specular=10)

pltr = pv.Plotter()
pltr.set_background("#363940")
drawP(pltr, 0.5, "yellow", "darkmagenta", "orangered")
pltr.show()