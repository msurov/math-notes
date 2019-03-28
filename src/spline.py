import numpy as np
import matplotlib.pyplot as plt


def arr_from_str(s):
    arr = s.split()
    return np.array([float(e) for e in arr])


def _a(t, x, i, p, r):
    return (x - t[i]) / (t[i+1+p-r] - t[i])


def _d(t, c, x, i, p, r):
    if (r == 0):
        return c[i]

    a = _a(t, x, i, p, r)
    return (1. - a) * _d(t, c, x, i - 1, p, r - 1) + a * _d(t, c, x, i, p, r-1)


def get_idx(t, x, p):
    if x < t[0]:
        l = p
    elif x >= t[-1]:
        l = len(t) - p - 2
    else:
        idx, = np.nonzero(~(t <= x))
        l = idx[0] - 1
    return l


def deBoor(t, c, x, p=3):
    l = get_idx(t, x, p)
    return _d(t, c, x, l, p, p)


def deBoor_der(t, c, x, p=3):
    l = get_idx(t, x, p)
    return _d(t, c, x, l, p, p)

'''
def _B(t, x, i, l, p):
    if p == 0:
        return 1 if i == l else 0
    if i < l-p:
        return 0
    if i > l:
        return 0
    if i == l-p:
        return (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * _B(t, x, i+1, l, p-1)
    if i == l:
        return (x - t[i]) / (t[i+p] - t[i]) * _B(t, x, i, l, p-1)

    return (x - t[i]) / (t[i+p] - t[i]) * _B(t, x, i, l, p-1) + \
        (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * _B(t, x, i+1, l, p-1)
'''

def spline_B(t, x, i, p=3):
    l = get_idx(t, x, p)

    def B(i,p):
        if p == 0:
            return 1 if i == l else 0
        if i < l-p:
            return 0
        if i > l:
            return 0
        if i == l-p:
            return (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)
        if i == l:
            return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1)

        return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1) + \
            (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)

    return B(i,p)


def spline_B_der(t, x, i, p=3):
    l = get_idx(t, x, p)
    N = len(t) - p - 1

    def B(i,p):
        if p == 0:
            return 1 if i == l else 0
        if i < l-p:
            return 0
        if i > l:
            return 0
        if i == l-p:
            return (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)
        if i == l:
            return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1)

        return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1) + \
            (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)

    def B_der(i, p=3):
        if i == N-1:
            return p / (t[i+p] - t[i]) * B(i, p-1)
        elif i == 0:
            return - p / (t[i+p+1] - t[i+1]) * B(i+1, p-1)
        else:
            return p / (t[i+p] - t[i]) * B(i, p-1) - p / (t[i+p+1] - t[i+1]) * B(i+1, p-1)

    return B_der(i,p)


def spline_basis(t, x, p=3):
    l = get_idx(t, x, p)

    def B(i,p):
        if p == 0:
            return 1 if i == l else 0
        if i < l-p:
            return 0
        if i > l:
            return 0
        if i == l-p:
            return (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)
        if i == l:
            return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1)

        return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1) + \
            (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)

    N = len(t) - p - 1
    B_arr = np.zeros(N, dtype='float')

    for i in range(0, N):
        B_arr[i] = B(i,p)

    return B_arr


def spline_basis_der(t, x, p=3):
    l = get_idx(t, x, p)
    n = len(t) - p - 1
    D_B_arr = np.zeros(n, dtype='float')

    def B(i,p):
        if p == 0:
            return 1 if i == l else 0
        if i < l-p:
            return 0
        if i > l:
            return 0
        if i == l-p:
            return (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)
        if i == l:
            return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1)

        return (x - t[i]) / (t[i+p] - t[i]) * B(i, p-1) + \
            (t[i+p+1] - x) / (t[i+p+1] - t[i+1]) * B(i+1, p-1)

    def B_der(i, p=3):
        if i == n-1:
            return p / (t[i+p] - t[i]) * B(i, p-1)
        elif i == 0:
            return - p / (t[i+p+1] - t[i+1]) * B(i+1, p-1)
        else:
            return p / (t[i+p] - t[i]) * B(i, p-1) - p / (t[i+p+1] - t[i+1]) * B(i+1, p-1)

    for i in range(0, n):
        D_B_arr[i] = B_der(i,p)

    return D_B_arr


def interpolate2(x, y, p):
    K = len(x)
    N = K + p - 1

    t = np.zeros(K + 2 * p, dtype=float)
    t[0:p] = x[0]
    t[-p:] = x[-1]
    t[p:-p] = x

    c = np.zeros(N, dtype=float)
    A = np.zeros((N, N), dtype='float')
    b = np.zeros(N, dtype='float')

    for k in range(0,K):
        for i in range(0,N):
            A[k,i] = spline_B(t, x[k], i, p)

    for i in range(0,N):
        A[K,i] = spline_B_der(t, t[0], i, p)
        A[K+1,i] = spline_B_der(t, t[-1], i, p)


    b[0:K] = y
    b[K:K+2] = 0

    c = np.linalg.inv(A).dot(b)
    return t,c


# t = arr_from_str("0        0        0        0 0.105263 0.157895 0.210526 0.263158 0.315789 0.368421 0.421053 0.473684 0.526316 0.578947 0.631579 0.684211 0.736842 0.789474 0.842105 0.894737        1    1        1        1")
# c = arr_from_str("-1.291e-18  0.0350878  0.0877192   0.157312   0.209071   0.260251    0.31071   0.360309    0.40891   0.456379   0.502583   0.547396   0.590693   0.632353   0.672263   0.710311   0.746391   0.791743   0.822514   0.841471")
# x = arr_from_str("0 0.0526316  0.105263  0.157895  0.210526  0.263158  0.315789  0.368421  0.421053  0.473684  0.526316  0.578947  0.631579  0.684211  0.736842  0.789474  0.842105  0.894737  0.947368         1")
# y = arr_from_str("0 0.0526073  0.105069  0.157239  0.208975  0.260131  0.310567  0.360143  0.408721  0.456168  0.502351  0.547143   0.59042  0.632061  0.671953  0.709983  0.746047  0.780044  0.811882  0.841471")

x = np.linspace(0,1,100)
y = np.exp(x)
plt.plot(x, y, 'o')

t,c = interpolate2(x, y, 3)

print t
print c

x = np.linspace(-0.1, 1.1, 1000)
y = np.array([deBoor(t, c, xi) for xi in x])
plt.plot(x, y)

y = np.exp(x)
plt.plot(x, y)

plt.grid()
plt.show()
