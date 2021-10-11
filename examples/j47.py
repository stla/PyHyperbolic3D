# -*- coding: utf-8 -*-
import pyvista as pv
from hyperbolic3d.gyro import gyrotriangle, gyrotube

s = 0.1
depth = 5
r = 0.007

vertices = [
 [-0.29686046437836, -0.171145846782417, -0.0472798406950809],
 [-0.292390323045577, -0.120125848492049, 0.140463870911089],
 [-0.278799092206407, 0.0157275970704114, -0.0984933785009629],
 [-0.274328950873624, 0.06674759536078, 0.0892503331052067],
 [-0.231585747542868, -0.124733213852467, -0.224637678084076],
 [-0.219882765599125, 0.00883887577767885, 0.266881740074936],
 [-0.169988681168102, -0.255822195021454, 0.0735679574259245],
 [-0.168712240427772, 0.141844812948812, -0.197720675961363],
 [-0.161479399816813, 0.22439690428859, 0.106055030591479],
 [-0.13802182406713, -0.258244430336326, -0.118377285963722],
 [-0.129519109980523, -0.161198626665482, 0.23873247466925],
 [-0.121498895764232, 0.00138400202593346, -0.323864975544475],
 [-0.107033214542315, 0.166488184705489, 0.283686437561207],
 [-0.0962046829813201, 0.270809537218541, -0.0713028067975163],
 [-0.0458287915766691, -0.167540121048566, -0.263786696503713],
 [-0.0320711111879154, -0.0105167122435899, 0.314029033914774],
 [-0.00864934470742099, 0.159033310953743, -0.307060278058203],
 [-0.0035773491929405, -0.328453827319776, 0.00354911555690258],
 [0.0030536372363208, 0.292605400583889, 0.184459140100808],
 [0.0619037924986957, -0.175349677583013, 0.270790918192074],
 [0.0638582127390304, 0.287998035223471, -0.180642408894355],
 [0.0683283540718134, 0.339018033513839, 0.00710130271181298],
 [0.0713758114178494, -0.0183552303800157, -0.307118821651861],
 [0.0851334918066029, 0.138668178424961, 0.270696908766624],
 [0.139740081268431, -0.27871907638255, 0.12543960867227],
 [0.145594110902549, -0.181691171966097, -0.231728252980892],
 [0.168823810210457, 0.132326684041877, -0.231822262406339],
 [0.171706938369401, -0.281141311697422, -0.0665056347173764],
 [0.177326524297064, 0.229372487712722, 0.125287498226633],
 [0.209293381398035, 0.226950252397849, -0.0666577451630137],
 [0.237188080061038, -0.128037161960658, 0.200736167917793],
 [0.251544823781764, 0.0660365461266406, 0.200678066897604],
 [0.28891154136392, -0.131956421028871, -0.109837759865526],
 [0.303268285084646, 0.0621172870584274, -0.109895860885714],
 [0.329381112551499, -0.0373328526728975, 0.0553267573778008]
]

edges = [
 [2, 3],
 [3, 8],
 [8, 13],
 [7, 13],
 [2, 7],
 [8, 12],
 [3, 5],
 [5, 12],
 [8, 18],
 [12, 18],
 [13, 21],
 [18, 21],
 [13, 20],
 [20, 21],
 [7, 16],
 [16, 20],
 [7, 11],
 [11, 16],
 [2, 4],
 [4, 11],
 [0, 2],
 [0, 4],
 [1, 3],
 [0, 1],
 [1, 5],
 [30, 34],
 [24, 30],
 [24, 27],
 [27, 32],
 [32, 34],
 [19, 30],
 [30, 31],
 [23, 31],
 [15, 23],
 [15, 19],
 [19, 24],
 [10, 19],
 [10, 15],
 [17, 24],
 [6, 10],
 [6, 17],
 [17, 27],
 [9, 17],
 [6, 9],
 [25, 27],
 [9, 14],
 [14, 25],
 [25, 32],
 [22, 25],
 [14, 22],
 [32, 33],
 [22, 26],
 [26, 33],
 [33, 34],
 [29, 33],
 [26, 29],
 [31, 34],
 [28, 29],
 [28, 31],
 [23, 28],
 [16, 22],
 [11, 22],
 [11, 14],
 [4, 14],
 [4, 9],
 [0, 9],
 [0, 6],
 [1, 6],
 [1, 10],
 [5, 10],
 [5, 15],
 [12, 15],
 [12, 23],
 [18, 23],
 [18, 28],
 [21, 28],
 [21, 29],
 [20, 29],
 [20, 26],
 [16, 26]
]

faces = [
    [3, 8, 13, 7, 2],
    [8, 3, 5, 12],
    [8, 12, 18],
    [13, 8, 18, 21],
    [13, 21, 20],
    [7, 13, 20, 16],
    [7, 16, 11],
    [2, 7, 11, 4],
    [2, 4, 0],
    [3, 2, 0, 1],
    [3, 1, 5],
    [30, 24, 27, 32, 34],
    [30, 31, 23, 15, 19],
    [30, 19, 24],
    [19, 15, 10],
    [24, 19, 10, 6, 17],
    [24, 17, 27],
    [17, 6, 9],
    [27, 17, 9, 14, 25],
    [27, 25, 32],
    [25, 14, 22],
    [32, 25, 22, 26, 33],
    [32, 33, 34],
    [33, 26, 29],
    [34, 33, 29, 28, 31],
    [34, 31, 30],
    [31, 28, 23],
    [16, 22, 11],
    [11, 22, 14],
    [11, 14, 4],
    [4, 14, 9],
    [4, 9, 0],
    [0, 9, 6],
    [0, 6, 1],
    [1, 6, 10],
    [1, 10, 5],
    [5, 10, 15],
    [5, 15, 12],
    [12, 15, 23],
    [12, 23, 18],
    [18, 23, 28],
    [18, 28, 21],
    [21, 28, 29],
    [21, 29, 20],
    [20, 29, 26],
    [20, 26, 16],
    [16, 26, 22]
]



pl = pv.Plotter()
for vertex in vertices:
    V = pv.Sphere(2 * r, center=vertex)
    pl.add_mesh(V, smooth_shading=True, color="yellow")
for edge in edges:
    A = vertices[edge[0]]
    B = vertices[edge[1]]
    AB = gyrotube(A, B, s, r)
    pl.add_mesh(AB, smooth_shading=True, color="yellow")
for face in faces:
    A = vertices[face[0]]
    B = vertices[face[1]]
    C = vertices[face[2]]
    tmesh1 = gyrotriangle(A, B, C, s, depth)
    if len(face) == 3:
        pl.add_mesh(tmesh1, smooth_shading=True, specular=5, color="navy")
    else:
        D = vertices[face[3]]
        tmesh2 = gyrotriangle(A, C, D, s, depth)
        tmesh12 = tmesh1.merge(tmesh2)
        if len(face) == 4:
            pl.add_mesh(tmesh12, smooth_shading=True, specular=5, color="navy")
        else:
            E = vertices[face[4]]
            tmesh3 = gyrotriangle(A, D, E, s, depth)
            tmesh123 = tmesh12.merge(tmesh3)
            pl.add_mesh(tmesh123, smooth_shading=True, specular=5, color="navy")
pl.show()

