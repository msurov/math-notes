import numpy as np
from misc.math import integrate
import matplotlib.pyplot as plt


A = np.array([
    [0., 1.],
    [0., 0.]
])
J = 0.008
B = np.array([
    [0.],
    [1./J]
])
Bd = np.array([
    [0.],
    [1./J*1.005]
])
C = np.array([1., 0.])
L = np.array([10., 200.])

K = A - np.outer(L, C)
l,_ = np.linalg.eig(K)
print l
exit()

def feedback(t, x):
    # return 0.1*t + 2.1 * np.sin(1.212345124*t)
    return 0.005*np.sin(t)


def dynamics(t, state):
    x = state[0:2]
    x_ = state[2:4]
    u = feedback(t, x_)
    dx = A.dot(x) + B.dot([u])
    y = C.dot(x)
    y_ = C.dot(x_)
    dx_ = A.dot(x_) + Bd.dot([u]) + L.dot(y - y_)
    return np.concatenate((dx, dx_))


x0 = [5., 0., 1., 0.1]
tspan = [0, 20]
t,x = integrate(dynamics, x0, tspan, step=1e-3)

e = np.linalg.norm(x[:,0:2] - x[:,2:4], axis=1)
# print e
# plt.plot(t, e)
plt.plot(t, x[:,2])
plt.plot(t, x[:,0])
plt.show()
