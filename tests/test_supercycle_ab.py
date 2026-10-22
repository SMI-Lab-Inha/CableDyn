# SPDX-License-Identifier: Apache-2.0
"""The rainflow count and damage-equivalent load of validation/scripts/openfast_supercycle_ab.py."""

from __future__ import annotations

import importlib.util
import unittest
from collections import Counter
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "validation" / "scripts" / "openfast_supercycle_ab.py"
_SPEC = importlib.util.spec_from_file_location("openfast_supercycle_ab", _SCRIPT)
ab = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ab)

# ASTM E1049-85 (2017), Fig. 6: the rainflow counts of this history
_ASTM_HISTORY = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0]
_ASTM_COUNTS = {3.0: 0.5, 4.0: 1.5, 6.0: 0.5, 8.0: 1.0, 9.0: 0.5}


def _counts(history):
    c = Counter()
    for r, n in ab.rainflow(history):
        c[r] += n
    return dict(c)


class RainflowTest(unittest.TestCase):
    def test_astm_example(self):
        self.assertEqual(_counts(_ASTM_HISTORY), _ASTM_COUNTS)

    def test_plateaus_count_once(self):
        # equal neighbours (a tension written to seven digits) neither add nor drop a turning point
        plateaued = [-2.0, 1.0, 1.0, -3.0, -3.0, 5.0, 5.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0]
        self.assertEqual(_counts(plateaued), _ASTM_COUNTS)
        self.assertEqual(ab.turning_points([1.0, 2.0, 2.0, 1.0]), [1.0, 2.0, 1.0])

    def test_damage_equivalent_load(self):
        # a sine of amplitude A over N whole periods in T seconds: N full cycles of range 2A
        a, n, t = 3.0, 40, 200.0
        history = [a * ((-1.0) ** k) for k in range(2 * n + 1)]
        expected = (n * (2.0 * a) ** 3 / t) ** (1.0 / 3.0)
        self.assertAlmostEqual(ab.damage_equivalent_load(history, t, 3.0), expected, places=12)


if __name__ == "__main__":
    unittest.main()
