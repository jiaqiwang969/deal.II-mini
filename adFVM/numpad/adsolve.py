# solve nonlinear systems, and differentiate through implicit relations
# that are established through nonlinear solvers
# Copyright (C) 2014
# Qiqi Wang  qiqi.wang@gmail.com
# engineer-chaos.blogspot.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division, print_function, absolute_import

import os
import sys
import unittest
import numbers
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg

sys.path.append(os.path.realpath('..')) # for running unittest

from numpad.adstate import _add_ops
from numpad.adarray import *
from numpad.adarray import __DEBUG_MODE__, _DEBUG_perturb_new

class ResidualState(IntermediateState):
    def __init__(self, prev_state):
        host = prev_state.host()
        IntermediateState.__init__(self, host, prev_state, 1, None)

    def tos(self):
        if hasattr(self, 'solution') and self.solution():
            yield self.solution()
        for state in IntermediateState.tos(self):
            yield state

    # --------- adjoint differentiation -------- #

    def diff_adjoint(self, f_diff_dependers):
        f_diff_self = 0
        iter_f_diff_dependers = iter(f_diff_dependers)

        if hasattr(self, 'solution') and self.solution():
            f_diff_soln = next(iter_f_diff_dependers)

            if f_diff_soln is not 0:
                # inverse of Jacobian matrix
                self_diff_soln = self.solution().jacobian.T
                # check if diagonal matrix
                n = self_diff_soln.shape[0]
                try:
                    is_diag = self_diff_soln.indices.size == n and \
                              self_diff_soln.indptr.size == n+1 and \
                              all(self_diff_soln.indices == np.arange(n)) and \
                              all(self_diff_soln.indptr == np.arange(n+1))
                except TypeError:
                    is_diag = False
                if is_diag:
                    assert f_diff_soln.shape[-1] == n
                    # inverse of diagonal matrix
                    soln_diff_self = self_diff_soln.copy()
                    soln_diff_self.data = 1. / np.array(soln_diff_self.data)
                    f_diff_self = -f_diff_soln * soln_diff_self
                else:
                    if hasattr(f_diff_soln, 'todense'):
                        f_diff_soln = f_diff_soln.todense()
                    f_diff_soln = np.array(f_diff_soln)
                    soln_diff_self = splinalg.factorized(self_diff_soln.tocsc())
                    f_diff_self = np.array([-soln_diff_self(b) \
                                            for b in f_diff_soln])
                    f_diff_self = f_diff_self.reshape(f_diff_soln.shape)
                    f_diff_self = sp.csr_matrix(f_diff_self)

        f_diff_self_1 = IntermediateState.diff_adjoint(self,
                                                       iter_f_diff_dependers)
        return _add_ops(f_diff_self, f_diff_self_1)


class SolutionState(IntermediateState):
    def __init__(self, host, residual_state, jacobian):
        IntermediateState.__init__(self, host, None, None, None)
        assert isinstance(residual_state, ResidualState)
        assert residual_state.size == self.size
        residual_state.solution = weakref.ref(self)
        self.residual = residual_state
        if not isinstance(jacobian, numbers.Number):
            assert jacobian.shape == (self.size, self.size)
        self.jacobian = jacobian

    def obliviate(self):
        IntermediateState.obliviate(self)
        self.residual = None
        self.jacobian = None

    def froms(self):
        if hasattr(self, 'residual') and self.residual:
            yield self.residual
        # it should have no other dependees
        assert self.prev is None and self.other is None

    # --------- tangent differentiation -------- #

    def diff_tangent(self, dependees_diff_u):
        if hasattr(self, 'residual') and self.residual:
            resid_diff_u, = dependees_diff_u
            if resid_diff_u is 0:
                return 0
            else:
                # inverse of Jacobian matrix
                resid_diff_self = self.jacobian
                # check if diagonal matrix
                n = resid_diff_self.shape[0]
                try:
                    is_diag = resid_diff_self.indices.size == n and \
                              resid_diff_self.indptr.size == n+1 and \
                              all(resid_diff_self.indices == np.arange(n)) and \
                              all(resid_diff_self.indptr == np.arange(n+1))
                except TypeError:
                    is_diag = False
                if is_diag:
                    # inverse of diagonal matrix
                    self_diff_resid = resid_diff_self.copy()
                    self_diff_resid.data = 1 / np.array(self_diff_resid.data)

                    self_diff_u = -self_diff_resid * resid_diff_u
                    return self_diff_u
                else:
                    if hasattr(resid_diff_u, 'todense'):
                        resid_diff_u = resid_diff_u.todense()
                    resid_diff_u = np.array(resid_diff_u)
                    self_diff_resid = splinalg.factorized(resid_diff_self.tocsc())
                    self_diff_u = np.transpose([-self_diff_resid(b) \
                                            for b in resid_diff_u.T])
                    self_diff_u = self_diff_u.reshape(resid_diff_u.shape)
                    return sp.csr_matrix(self_diff_u)
        else:
            return 0


class adsolution(adarray):
    def __init__(self, solution, residual, n_Newton):
        assert isinstance(solution, adarray)
        assert isinstance(residual, adarray)

        residual._current_state = ResidualState(residual._current_state)

        adarray.__init__(self, solution._value)
        self._current_state = SolutionState(self, residual._current_state,
                                            residual.diff(solution))
        self._n_Newton = n_Newton
        self._res_norm = np.linalg.norm(residual._value.reshape(residual.size))

        _DEBUG_perturb_new(self)

    def obliviate(self):
        adarray.obliviate(self)
        del self._n_Newton
        del self._res_norm


def solve_newton_with_dt(func, u0, args, kargs, dt,
          max_iter, abs_tol, rel_tol, verbose):
    u = adarray(value(u0).copy())
    _DEBUG_perturb_new(u)
    for i_Newton in range(max_iter):
        if dt == np.inf:
            res = func(u, *args, **kargs)
        else:
            res = (u - u0) / dt + func(u, *args, **kargs)
        res_norm = np.linalg.norm(res._value.reshape(res.size), np.inf)
        if verbose:
            print('    ', i_Newton, res_norm)

        if i_Newton == 0:
            res_norm0 = res_norm
        if res_norm < max(abs_tol, rel_tol * res_norm0):
            return adsolution(u, res, i_Newton + 1)
        if not np.isfinite(res_norm) or res_norm > res_norm0 * 1E6:
            break

        # Newton update
        J = res.diff(u).tocsr()
        if J.shape[0] > 1:
            minus_du = splinalg.spsolve(J, np.ravel(res._value),
                                        use_umfpack=False)
        else:
            minus_du = res._value / J.toarray()[0,0]
        u._value -= minus_du.reshape(u.shape)
        u = adarray(u._value)  # unlink operation history if any
        _DEBUG_perturb_new(u)

    # not converged
    u = adarray(value(u0).copy())
    res = func(u, *args, **kargs)
    return adsolution(u, res, np.inf)

def psuedo_time_continuation(func, u0, args, kargs, dt_min, dt_max, dt_ratio,
          max_iter, abs_tol, rel_tol, verbose):
    dt = dt_min
    while np.isfinite(dt):
        if np.abs(dt) > np.abs(dt_max):
            dt = np.inf
        if verbose:
            print('continuation step with dt = ', dt)
        u = solve_newton_with_dt(func, u0, args, kargs, dt,
                max_iter, abs_tol, rel_tol, verbose)
        if not np.isfinite(u._n_Newton):
            return dt, None
        else:
            dt *= dt_ratio
            u0 = u
    if np.isfinite(u._n_Newton):
        return np.inf, u

def solve(func, u0, args=(), kargs={},
          max_iter=10, abs_tol=1E-6, rel_tol=1E-6, verbose=True):
    func = replace__globals__(func)

    # try straight Newton first
    u = solve_newton_with_dt(func, u0, args, kargs, np.inf,
                max_iter, abs_tol, rel_tol, verbose)
    if np.isfinite(u._n_Newton):
        return u

    if verbose:
        print('Newton failed to converge, starting psuedo time continuation')

    dt_min, dt_max, dt_ratio = 1, 4, 2
    u = None
    while u is None:
        dt_diverged_p, u = psuedo_time_continuation(
                func, u0, args, kargs, dt_min, dt_max, dt_ratio,
                max_iter, abs_tol, rel_tol, verbose)
        if u is not None: return u

        dt_diverged_m, u = psuedo_time_continuation(
                func, u0, args, kargs, -dt_min, -dt_max, dt_ratio,
                max_iter, abs_tol, rel_tol, verbose)
        if u is not None: return u

        if np.abs(dt_diverged_p) == dt_min and np.abs(dt_diverged_m) == dt_min:
            if verbose:
                print('Continuation failed at start, halfing min dt')
            dt_min /= 2
        elif np.isinf(dt_diverged_p) or np.isinf(dt_diverged_m):
            if verbose:
                print('Continuation failed at end, doubling max dt')
            dt_max *= 2
        else:
            if verbose:
                print('Continuation failed in middle')
                print('Decreasing step to ', sqrt(dt_ratio))
            dt_min /= dt_ratio
            dt_ratio = sqrt(dt_ratio)
        sys.stdout.flush()



# =========================================================== #
#                                                             #
#                         unittests                           #
#                                                             #
# =========================================================== #

class _Poisson1dTest(unittest.TestCase):
    def residual(self, u, f, dx):
        res = -2 * u
        res[1:] += u[:-1]
        res[:-1] += u[1:]
        return res / dx**2 + f

    def testPoisson1d(self):
        N = 256
        dx = adarray(1. / N)

        f = ones(N-1)
        u = zeros(N-1)

        u = solve(self.residual, u, (f, dx))

        x = np.linspace(0, 1, N+1)[1:-1]
        self.assertAlmostEqual(0, np.abs(u._value - 0.5 * x * (1 - x)).max())

        # solve tangent equation
        dudx = np.array(u.diff(dx).todense()).reshape(u.shape)
        self.assertAlmostEqual(0, np.abs(dudx - 2 * u._value / dx._value).max())

        # solve adjoint equation
        J = u.sum()
        dJdf = J.diff(f)
        self.assertAlmostEqual(0, np.abs(dJdf - u._value).max())


class _Poisson2dTest(unittest.TestCase):
    def residual(self, u, f, dx, dy):
        res = -(2 / dx**2 + 2 / dy**2) * u
        res[1:,:] += u[:-1,:] / dx**2
        res[:-1,:] += u[1:,:] / dx**2
        res[:,1:] += u[:,:-1] / dy**2
        res[:,:-1] += u[:,1:] / dy**2
        res += f
        return res

    def testPoisson2d(self):
        #N, M = 256, 512
        N, M = 256, 64
        dx, dy = adarray([1. / N, 1. / M])

        f = ones((N-1, M-1))
        u = ones((N-1, M-1))

        u = solve(self.residual, u, (f, dx, dy))

        x = np.linspace(0, 1, N+1)[1:-1]
        y = np.linspace(0, 1, M+1)[1:-1]

        # solve tangent equation
        dudx = np.array(u.diff(dx).todense()).reshape(u.shape)
        dudy = np.array(u.diff(dy).todense()).reshape(u.shape)

        self.assertAlmostEqual(0,
            abs(2 * u._value - (dudx * dx._value + dudy * dy._value)).max())

        # solve adjoint equation
        J = u.sum()
        dJdf = J.diff(f)

        self.assertAlmostEqual(0, abs(np.ravel(u._value) - dJdf).max())


class _HardSolveTest(unittest.TestCase):
    def testSin2xPlusX(self):
        def sin2xPlusX(x):
            return sin(2 * x) + x
        x = solve(sin2xPlusX, array(100), rel_tol=1E-12, abs_tol=1E-12,
                verbose=False)
        self.assertAlmostEqual(0, value(x))

    def testSinAxPlusX(self):
        N = 10
        A = linspace(0, 2, N)
        def sinAxPlusX(x):
            return sin(A * x) + x
        x = solve(sinAxPlusX, 100 * ones(N), rel_tol=1E-12, abs_tol=1E-12,
                verbose=False)
        self.assertAlmostEqual(0, np.linalg.norm(value(x)))

if __name__ == '__main__':
    # _Poisson2dTest().testPoisson2d()
    unittest.main()
