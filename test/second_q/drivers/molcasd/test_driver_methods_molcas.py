# This code is part of a Qiskit project.
#
# (C) Copyright IBM 2019, 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

""" Test Driver Methods OpenMolcas """

import unittest

from test.second_q.drivers.test_driver_methods_gsc import TestDriverMethods
from qiskit_nature.second_q.drivers import MolcasDriver
import qiskit_nature.optionals as _optionals


class TestDriverMethodsMolcas(TestDriverMethods):
    """Driver Methods OpenMolcas tests"""

    molcas_h2_config = """
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

    @unittest.skipIf(not _optionals.HAS_MOLCAS, "OpenMolcas not available.")
    def setUp(self):
        super().setUp()
        self.h2 = "H 0 0 0; H 0 0 0.7"
        self.ref_energies = {"h2": -1.137675}  # FCI energy for H2 at 0.7 A with sto-3g basis
        MolcasDriver(config=self.molcas_h2_config)

    def test_h2_rhf(self):
        """H2 rhf test"""
        driver = MolcasDriver(config=self.molcas_h2_config)
        result = self._run_driver(driver)
        self._assert_energy(result, "h2")


if __name__ == "__main__":
    unittest.main()
