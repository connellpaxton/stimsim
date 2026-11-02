"""Numerical checks for the passive fallback, independent of NEURON."""

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from analyze import load_swc
from passive_demo import factor_tree, simulate, solve_tree


class PassiveTests(unittest.TestCase):
    def test_tree_solver_matches_dense_linear_solve(self):
        parent = np.array([-1, 0, 0, 1, 1, 2])
        order = np.array([0, 1, 2, 3, 4, 5])
        coupling = np.array([0., .3, 1.1, .7, 2., .4])
        diagonal = np.full(6, .5) + coupling
        np.add.at(diagonal, parent[1:], coupling[1:])
        matrix = np.diag(diagonal)
        for i in order[1:]:
            matrix[i, parent[i]] = matrix[parent[i], i] = -coupling[i]
        rhs = np.array([1., 2., -3., .5, 0., 4.])
        actual = solve_tree(rhs, factor_tree(diagonal, coupling, parent, order), coupling, parent, order)
        np.testing.assert_allclose(actual, np.linalg.solve(matrix, rhs), rtol=1e-12, atol=1e-12)

    def test_single_soma_matches_exact_rc_and_converges(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'soma.swc'
            path.write_text('0 1 0 0 0 10 -1\n')
            errors = []
            for dt in [.1, .05]:
                times, traces, info = simulate(path, duration=3, dt=dt, amplitude=.01)
                leak = info['area_um2'] * 6.3e-7
                tau = info['passive_time_constant_ms']
                expected = -65 + .01 / leak * (1 - np.exp(-1 / tau)) * np.exp(-1 / tau)
                errors.append(abs(traces[-1, 0] - expected))
                self.assertAlmostEqual(traces[0, 0], -65)
            self.assertLess(errors[1], errors[0])
            self.assertLess(errors[1], .005)

    def test_unstimulated_branched_tree_stays_at_rest(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'tree.swc'
            path.write_text('3 3 20 10 0 1 1\n0 1 0 0 0 5 -1\n1 3 10 0 0 1 0\n2 2 20 -10 0 1 1\n')
            _, traces, _ = simulate(path, duration=1, amplitude=0)
            np.testing.assert_array_equal(traces, np.full_like(traces, -65))

    def test_invalid_morphology_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.swc'
            for content in ['0 1 0 0 0 5 -1\n1 3 1 0 0 1 9\n', '0 1 0 0 0 5 -1\n1 3 1 0 0 1 2\n2 3 2 0 0 1 1\n']:
                path.write_text(content)
                with self.assertRaises(ValueError):
                    load_swc(path)


if __name__ == '__main__':
    unittest.main()
