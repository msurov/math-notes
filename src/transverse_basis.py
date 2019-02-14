import numpy as np
import matplotlib.pyplot as plt
from misc.math import integrate



class Curve:
    def __init__(self, dim):
        self.nfreq = 2
        self.ndim = dim

        self.a = np.random.random((self.nfreq, self.ndim)) - 0.5
        self.b0 = np.random.random(self.ndim) - 0.5
        self.b = np.random.random((self.nfreq, self.ndim)) - 0.5

        self.n = np.arange(1, self.nfreq+1)
        self.an = np.reshape(self.n, (-1,1)) * self.a
        self.bn = np.reshape(self.n, (-1,1)) * self.b
        self.an2 = np.reshape(self.n**2, (-1,1)) * self.a
        self.bn2 = np.reshape(self.n**2, (-1,1)) * self.b


    def val(self, t):
        sin_nt = np.sin(np.outer(t, self.n))
        cos_nt = np.cos(np.outer(t, self.n))
        return sin_nt.dot(self.a) + cos_nt.dot(self.b) + self.b0


    def der1(self, t):
        '''
            d gamma / dt
        '''
        sin_nt = np.sin(np.outer(t, self.n))
        cos_nt = np.cos(np.outer(t, self.n))
        return cos_nt.dot(self.an) - sin_nt.dot(self.bn)


    def der2(self, t):
        '''
            d2 gamma / dt2
        '''
        sin_nt = np.sin(np.outer(t, self.n))
        cos_nt = np.cos(np.outer(t, self.n))
        return -sin_nt.dot(self.an2) - cos_nt.dot(self.bn2)


    def tangent(self, t):
        '''
            d gamma / ds
        '''
        d = self.der1(t)
        k = 1. / np.linalg.norm(d, axis=1)
        k = np.reshape(k, (len(d), 1))
        return k * d


    def curvature(self, t):
        '''
            d2gamma/ds2
        '''
        d1 = self.der1(t)
        d2 = self.der2(t)
        k = 1. / np.linalg.norm(d1, axis=1)
        k = np.reshape(k, (len(d1), 1))
        d1d2 = np.reshape(np.sum(d1 * d2, axis=1), (-1,1))
        return k**2 * d2 - k**4 * d1 * d1d2


    def speed(self, t):
        '''
            ds/dt
        '''
        d = self.der1(t)
        v = np.linalg.norm(d, axis=1)
        v = np.reshape(v, (len(d), 1))
        return v


    def __call__(self, t):
        return self.val(t)



def get_perp(v):
    v = np.reshape(v, (-1,1))
    U,l,Vt = np.linalg.svd(v.T)
    perp = Vt[1:,:] / l[0]
    return perp.T


def get_basis(v):
    v = np.reshape(v, (-1,1))
    U,l,Vt = np.linalg.svd(v.T)
    d = len(v)
    basis = Vt / np.linalg.det(Vt)**d
    return basis


def test():
    np.random.seed(3)
    t = np.linspace(0, 4*np.pi, 10000)
    ndim = 5
    gamma = Curve(ndim)


    def A(t):
        cur, = gamma.curvature(t)
        d, = gamma.der1(t)
        k = 0.
        A = (k * np.eye(ndim) - np.outer(d, cur))
        return A

    def rhs(t,E):
        return A(t).dot(E)


    tau0, = gamma.tangent(0.)
    E0 = get_perp(tau0)
    t, E = integrate(rhs, E0, [0, 1], step=1e-3)
    assert np.all(np.abs([gamma.tangent(ti).dot(Ei) for ti,Ei in zip(t,E)]) < 1e-10), 'E is not perpendicular to tangent'
    assert np.all([np.allclose(np.dot(Ei.T, Ei), np.eye(ndim-1)) for Ei in E]), 'E is not orthonormal' 


if __name__ == '__main__':
    test()
