# This code is part of a Qiskit project.
#
# (C) Copyright IBM 2020, 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

""" Test Driver OpenMolcas Extra """

import unittest

from test import QiskitNatureTestCase
import numpy as np
from qiskit_nature.second_q.drivers import MolcasDriver
from qiskit_nature import QiskitNatureError
import qiskit_nature.optionals as _optionals


class TestDriverPsi4Extra(QiskitNatureTestCase):
    """OpenMolcas Driver extra tests for driver specifics, errors etc"""

    @unittest.skipIf(not _optionals.HAS_MOLCAS, "OpenMolcas not available.")
    def setUp(self):
        super().setUp()

    def test_input_format_list(self):
        """input as a list"""
        driver = MolcasDriver(
            [
                "&GATEWAY",
                "coord",
                "2",
                "Angstrom",
                "H  0.0 0.0 0.0",
                "H  0.0 0.0 0.7",
                "basis=sto-3g",
                "Group=Nosym",
                "&SEWARD",
                "&SCF",
                "&RASSCF",
                "Spin= 1; Nactel= 2 0 0; Inactive= 0; Ras2= 2",
                "DMPO",
            ]
        )
        driver_result = driver.run()
        with self.subTest("energy"):
            self.assertAlmostEqual(driver_result.reference_energy, -1.117, places=3)

    def test_input_format_string(self):
        """input as a multi line string"""
        cfg = """
&GATEWAY
coord
2
Angstrom
H 0.0 0.0 0.0
H 0.0 0.0 0.7
basis=sto-3g
Group=Nosym
&SEWARD
&SCF
&RASSCF
Spin= 1; Nactel= 2 0 0; Inactive= 0; Ras2= 2
DMPO
"""
        driver = MolcasDriver(cfg)
        driver_result = driver.run()
        with self.subTest("energy"):
            self.assertAlmostEqual(driver_result.reference_energy, -1.117, places=3)

    def test_input_format_fail(self):
        """input type failure"""
        with self.assertRaises(QiskitNatureError):
            _ = MolcasDriver(1.000)


if __name__ == "__main__":
    unittest.main()
