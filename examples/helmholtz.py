#! /usr/bin/env python3
#
# In this script we solve the 1-D helmholtz equation :math:`u_{,ii} + k^2 u = -f` on a unit
# interval :math:`I`, subject to boundary
# conditions:
#
# .. math:: u(0) &= 0
#
#           u'(1) - i k u(1) &= 0
#

from nutils import mesh, function, solver, export, cli, testing
import numpy, treelog

def main(nelems:int, btype:str, degree:int):
  '''
  Helmholtz problem on a unit interval.

  .. arguments::

     nelems [10]
       Number of elements along edge.
     btype [std]
       Type of basis function (std/spline), availability depending on the
       selected element type.
     degree [1]
       Polynomial degree.
  '''

  domain, geom = mesh.rectilinear([numpy.linspace(0, 1, nelems + 1)])

  ns = function.Namespace()
  ns.x = geom
  ns.basis = domain.basis('std', degree=degree)
  ns.lhs = function.Argument('lhs', ns.basis.shape, dtype=complex)
  ns.u = 'basis_n lhs_n'
  ns.k = 0.2
  ns.I = 1j

  # .. math:: ∀ v: ∫_Ω (u' cong(v)' - k^2 u cong(v) d:x - i k u(1) cong(v(1)) = ∫_Ω f cong(u) d:x

  res = domain.integral('(d(u, x_i) conj(d(basis_n, x_i)) - k^2 u conj(basis_n)) J(x)' @ ns, degree=degree * 2)
  res -= domain.boundary["right"].integral("I k u conj(basis_n) J(x)" @ ns, degree=degree * 2)

  A_ij, b_i = domain.boundary["left"].integrate([ns.eval_ij('basis_i basis_j J(x)'), ns.eval_i('0 basis_i J(x)')], degree=3)
  bcons = A_ij.rowsupp(1e-15)
  cons = numpy.full_like(b_i, numpy.nan, dtype=float)
  cons[bcons] = A_ij.submatrix(bcons, bcons).solve(b_i[bcons])

  lhs = solver.solve_linear('lhs', res, constrain=cons)

  return cons, lhs

if __name__ == '__main__':
  cli.run(main)

# class test(testing.TestCase):

#   @testing.requires('matplotlib')
#   def test_default(self):
#     cons, lhs = main(nelems=4, btype='std', degree=1)
