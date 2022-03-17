import numpy as np
import sympy as sy
from ortho_poly import find_orthonormal_basis


ncoefs = 10
x = sy.symbols('x', real=True)
basis = find_orthonormal_basis(x, ncoefs)
basis_coefs = [np.array(e.all_coeffs(), float) for e in basis]

# 1. Collocation points
roots = np.roots(basis_coefs[-1])
roots = np.sort(roots)

# 2. Derivative matrix
