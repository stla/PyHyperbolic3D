# -*- coding: utf-8 -*-
import math
import numpy as np
import pyvista as pv


def opposite(A):
    return -A if isinstance(A, np.ndarray) else [-x for x in A]


def dotprod(X, Y=None):
    return np.dot(X, (X if Y is None else Y))


def betaF(A, s=1.0):
    return 1.0 / math.sqrt(1.0 + dotprod(A) / s ** 2)


def gyroadd(A, B, s=1.0):
    A_ = np.array(A)
    B_ = np.array(B)
    betaA = betaF(A_, s)
    betaB = betaF(B_, s)
    return (
        1.0 + (betaA / (1.0 + betaA)) * dotprod(A_, B_) / s ** 2 + 
        (1.0 - betaB) / betaB 
    ) * A_ + B_


def gyroscalar(r, A, s=1.0):
    A_ = np.array(A)
    Anorm = math.sqrt(dotprod(A_))
    return (s / Anorm) * np.sinh(r * np.arcsinh(Anorm / s)) * A_


def gyroABt(A, B, t, s=1.0):
    return gyroadd(A, gyroscalar(t, gyroadd(opposite(A), B, s), s), s)


def gyromidpoint(A, B, s=1.0):
    return gyroABt(A, B, 0.5, s).tolist()


def gyrosegment(A, B, s=1.0, n=50):
    t = [i / n for i in range(n)]
    t.append(1.0)
    return [np.transpose(gyroABt(A, B, t_, s)).tolist() for t_ in t]


def gyrotube(A, B, s, r, npoints=300):
    """
    Tubular hyperbolic segment.

    Parameters
    ----------
    A,B : points (lists or arrays)
        The two endpoints of the segment.
    s : positive float
        Curvature parameter.
    r : positive float
        Radius of the tube.
    npoints : integer
        Number of points along the segment. The default is 300.

    Returns
    -------
    PyVista mesh
        A PyVista mesh ready for inclusion in a plotting region.

    """
    AB = gyrosegment(A, B, s)
    splineAB = pv.Spline(AB, npoints)
    splineAB["scalars"] = np.arange(splineAB.n_points)
    return splineAB.tube(radius=r).extract_geometry()


def gyrosubdiv(A1, A2, A3, s=1.0):
    M12 = gyromidpoint(A1, A2, s)
    M13 = gyromidpoint(A1, A3, s)
    M23 = gyromidpoint(A2, A3, s)
    return [[A1, M12, M13], [A2, M23, M12], [A3, M13, M23], [M12, M13, M23]]


def gyrotriangle(A, B, C, s, depth=5, tol=1e-6):
    """
    Hyperbolic triangle.

    Parameters
    ----------
    A,B,C : points (lists or arrays)
        The vertices of the triangle.
    s : positive float
        Curvature parameter.
    depth : integer
        The number of recursive subdivions. The default is 5.
    tol : small positive float
        The tolerance used to merge duplicated points in the mesh.
        The default is 1e-6.

    Returns
    -------
    PyVista mesh
        A PyVista mesh ready for inclusion in a plotting region.

    """
    subd = gyrosubdiv(A, B, C, s)
    for _ in range(depth - 1):
        lst = list(map(lambda t: gyrosubdiv(*t, s), subd))
        subd = np.concatenate(lst)
    tsubd = list(map(np.transpose, subd))
    vertices = np.transpose(np.concatenate(tsubd, axis=1))
    nvertices = vertices.shape[0]
    ntriangles = nvertices // 3
    repeats3 = np.full((ntriangles, 1), 3)
    indices = np.array(range(nvertices)).reshape(ntriangles, 3)
    faces = np.hstack((repeats3, indices)).flatten()
    return pv.PolyData(vertices, faces).clean(tolerance=tol)
