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

def main(nelems: int, btype: str, degree: int, k: int = 3, f: str = 'sin(x_0) - cos(x_0)'):
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
     k [3]
       Wave number.
     f ['sin(x_0) - cos(x_0)']
       Inhomogeneous function in u'' + k^2 u = -f, a a function of 'x_0'.
  '''

  domain, geom = mesh.rectilinear([numpy.linspace(0, 1, nelems + 1)])

  ns = function.Namespace()
  ns.x = geom
  ns.basis = domain.basis('std', degree=degree)
  ns.lhs = function.Argument('lhs', ns.basis.shape, dtype=complex)
  ns.u = 'basis_n lhs_n'
  ns.I = 1j
  ns.k = k
  ns.f = f

  # .. math:: ∀ v: ∫_Ω (u' cong(v)' - k^2 u cong(v) d:x - i k u(1) cong(v(1)) = ∫_Ω f cong(v) d:x
  # by substitution: :math:`\bar{v} -> v`, this is equivalent to
  # .. math:: ∀ v: ∫_Ω (u' v' - k^2 u v d:x - i k u(1) v(1) = ∫_Ω f v d:x
  res = domain.integral('(d(u, x_i) conj(d(basis_n, x_i)) - k^2 u conj(basis_n)) J(x)' @ ns, degree=degree * 2)
  res -= domain.boundary["right"].integral("I k u conj(basis_n) J(x)" @ ns, degree=degree * 2)
  res -= domain.integral('f conj(basis_n) J(x)' @ ns, degree=degree * 2)

  A_ij, b_i = domain.boundary["left"].integrate(
      [ns.eval_ij('basis_i basis_j J(x)'), ns.eval_i('0 basis_i J(x)')], degree=3)
  bcons = A_ij.rowsupp(1e-15)
  cons = numpy.full_like(b_i, numpy.nan, dtype=float)
  cons[bcons] = A_ij.submatrix(bcons, bcons).solve(b_i[bcons])

  lhs = solver.solve_linear('lhs', res, constrain=cons)

  ns.sol = '''
    (1 / (-k + k^3)) (k cos(x_0) - k cos(k x_0) - k sin(x_0) - I k sin(k x_0) +
    cos(1) cos(k) sin(k x_0) + I k cos(1) cos(k) sin(k x_0) +
    cos(k) sin(1) sin(k x_0) - I k cos(k) sin(1) sin(k x_0) +
    I cos(1) sin(k) sin(k x_0) - k cos(1) sin(k) sin(k x_0) +
    I sin(1) sin(k) sin(k x_0) + k sin(1) sin(k) sin(k x_0))
  '''

  err = domain.integral('(abs(u - sol)^2 J(x)' @ ns, degree=degree * 2).eval(lhs=lhs)**.5
  treelog.user('L2 error: {:.2e}'.format(err))

  return cons, lhs, err

if __name__ == '__main__':
  cli.run(main)

class test(testing.TestCase):
  @testing.requires('matplotlib')
  def test_default(self):
    cons, lhs, err = main(nelems=4, btype='std', degree=1)
