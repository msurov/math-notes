from sympy.diffgeom import Manifold, Patch, CoordSystem
from sympy import cos, sin

import sympy as sy
import sympy.diffgeom as dif
from sympy.diffgeom.rn import R2
from sympy.diffgeom import CovarDerivativeOp
from sympy.diffgeom import metric_to_Christoffel_2nd, metric_to_Christoffel_1st, metric_to_Riemann_components, metric_to_Ricci_components
from sympy.diffgeom import TensorProduct as TP


a,g,l,E = sy.symbols('a g l E', positive=True)

Mfd = Manifold('M', 2)
Ptch = Patch('P', Mfd)
CS = CoordSystem('cyl', Ptch, ['x', 'phi'])
x,phi = CS.coord_functions()
dx,dphi = CS.base_oneforms()
ex,ephi = CS.base_vectors()

_g = a * TP(dx,dx) + -sin(phi) * TP(dx,dphi) - sin(phi) * TP(dphi,dx) + l * TP(dphi,dphi)
mu = (E - g * sin(phi)) * _g
G = metric_to_Christoffel_2nd(mu)

print G[0,:,:]
print G[1,:,:]
