# This code is part of a Qiskit project.
#
# (C) Copyright IBM 2018, 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

""" Test Driver OpenMolcas """

import unittest

from test import QiskitNatureTestCase
from test.second_q.drivers.test_driver import TestDriver
from qiskit_nature.second_q.drivers import MolcasDriver
import qiskit_nature.optionals as _optionals


class TestDriverMolcas(QiskitNatureTestCase, TestDriver):
    """Molcas Driver tests."""

    @unittest.skipIf(not _optionals.HAS_MOLCAS, "OpenMolcas not available.")
    def setUp(self):
        super().setUp()
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
        self.driver_result = driver.run()


class TestDriverMolcasMolecule(QiskitNatureTestCase, TestDriver):
    """Molcas Driver molecule tests."""

    @unittest.skipIf(not _optionals.HAS_MOLCAS, "OpenMolcas not available.")
    def setUp(self):
        super().setUp()
        driver = MolcasDriver.from_molecule(TestDriver.MOLECULE)
        self.driver_result = driver.run()


if __name__ == "__main__":
    unittest.main()
