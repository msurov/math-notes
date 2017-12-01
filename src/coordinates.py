import numpy  as np
from numpy import sin,cos,array
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def Rz(a):
    return array([
        [cos(a), -sin(a), 0, 0],
        [sin(a),  cos(a), 0, 0],
        [     0,       0, 1, 0],
        [     0,       0, 0, 1]
    ])


def Ry(a):
    return array([
        [ cos(a),  0, sin(a), 0],
        [      0,  1,      0, 0],
        [-sin(a),  0, cos(a), 0],
        [      0,  0,      0, 1]
    ])


def Dz(d):
    return array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, d],
        [0, 0, 0, 1]
    ])


def draw_line(ax, p1, p2):
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 'o-', label='parametric curve')


def origin(P):
    return P[0:3,3]


def draw_robot(q, l):
    p = np.array([0,0,0,1])

    P0 = np.eye(4)
    P1 = P0.dot(Dz(l[0]).dot(Rz(q[0])))
    P2 = P1.dot(Dz(l[1]).dot(Ry(q[1])))
    P3 = P2.dot(Dz(l[2]).dot(Ry(q[2])))
    P4 = P3.dot(Dz(l[3]).dot(Rz(q[3])))
    P5 = P4.dot(Dz(l[4]).dot(Ry(q[4])))
    P6 = P5.dot(Dz(l[5]).dot(Rz(q[5])))
    P7 = P6.dot(Dz(l[6]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    draw_line(ax, origin(P0), origin(P1))
    draw_line(ax, origin(P1), origin(P2))
    draw_line(ax, origin(P2), origin(P3))
    draw_line(ax, origin(P3), origin(P4))
    draw_line(ax, origin(P4), origin(P5))
    draw_line(ax, origin(P5), origin(P6))
    draw_line(ax, origin(P6), origin(P7))
    plt.show()


q = [0,0.1,0.2,0.3,0.4,0.5]
l = [1,1,1,1,1,1,1]
draw_robot(q, l)
