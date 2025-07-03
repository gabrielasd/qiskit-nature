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

""" OpenMolcas Driver."""

from __future__ import annotations

import io
import logging
import os
import sys
import subprocess
import tempfile
from typing import Any, TYPE_CHECKING

import numpy as np

from qiskit_nature import QiskitNatureError
from qiskit_nature.constants import PERIODIC_TABLE
from qiskit_nature.units import DistanceUnit
from qiskit_nature.exceptions import UnsupportMethodError
import qiskit_nature.optionals as _optionals
from qiskit_nature.settings import settings
from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo
from qiskit_nature.second_q.formats.qcschema import QCSchema
from qiskit_nature.second_q.formats.qcschema_translator import qcschema_to_problem
from qiskit_nature.second_q.formats.fcidump import FCIDump
from qiskit_nature.second_q.problems import ElectronicBasis, ElectronicStructureProblem
from qiskit_nature.utils import get_einsum

from ..electronic_structure_driver import ElectronicStructureDriver, MethodType, _QCSchemaData

logger = logging.getLogger(__name__)


@_optionals.HAS_MOLCAS.require_in_instance
class MolcasDriver(ElectronicStructureDriver):
    """
    Qiskit Nature driver using the OpenMolcas program.
    """

    def __init__(
        self,
        config: str | list[str] = "&GATEWAY\ncoord\n2\nangstrom\nH  0.0 0.0 0.0\nH  0.0 0.0 0.735\n"
        "basis=sto-3g\nGroup=Nosym\n&SEWARD\n&SCF\n"
        "&RASSCF\nSpin= 1; Nactel= 2 0 0; Inactive= 0; Ras2= 2\nDMPO\n\n",
    ) -> None:
        """
        Args:
            config: A molecular configuration conforming to Molcas format.

        Raises:
            QiskitNatureError: Invalid Input
        """
        super().__init__()
        if not isinstance(config, str) and not isinstance(config, list):
            raise QiskitNatureError(f"Invalid config for Molcas Driver '{config}'")

        if isinstance(config, list):
            config = "\n".join(config)

        # TODO: Check that the config contains the required sections specially the RASSCF section
        self._config = config
        self._qcschemadata = _QCSchemaData()

    @staticmethod
    @_optionals.HAS_MOLCAS.require_in_call
    def from_molecule(
        molecule: MoleculeInfo,
        *,
        basis: str = "sto-3g",
        method: MethodType = MethodType.RHF,
        driver_kwargs: dict[str, Any] | None = None,
    ) -> "MolcasDriver":
        """Creates a driver from a molecule.

        Args:
            molecule: the molecular information.
            basis: the basis set.
            method: the SCF method type.
            driver_kwargs: keyword arguments to be passed to driver.

        Returns:
            The constructed driver instance.

        Raises:
            QiskitNatureError: when an unknown unit is encountered.
        """
        # TODO: Figure out how to specify the number of orbitals in the basis set.
        # This is needed for the RASSCF calculation that generates the FCIDUMP.

        # # Ignore kwargs parameter for this driver
        # del driver_kwargs
        # MolcasDriver.check_method_supported(method)
        # basis = MolcasDriver.to_driver_basis(basis)

        # if molecule.units == DistanceUnit.ANGSTROM:
        #     units = "Angstrom"
        # elif molecule.units == DistanceUnit.BOHR:
        #     units = "Bohr"
        # else:
        #     raise QiskitNatureError(f"Unknown unit '{molecule.units.value}'")
        
        # name = "".join(molecule.symbols)
        # geom = "\n".join(
        #     [
        #         name + " " + " ".join(map(str, coord))
        #         for (name, coord) in zip(molecule.symbols, molecule.coords)
        #     ]
        # )
        # cfg1 = "&GATEWAY\n"
        # cfg2 = f"coord\n{len(molecule.symbols)}\n{units}\n"
        # cfg2 += f"{geom}\nbasis={basis}\nGroup=Nosym\n"
        # cfg3 = f"&SEWARD\n&SCF\n\n"

        # nbasis = None
        # nelec = sum([PERIODIC_TABLE.index(symbol) for symbol in molecule.symbols])
        # nelec -= molecule.charge
        # cfg4 = f"&RASSCF\nSpin= 1; Nactel= {nelec} 0 0; Inactive=0; Ras2={nbasis}\nDMPO\n\n"

        # return MolcasDriver(cfg1 + cfg2 + cfg3)
        raise NotImplementedError()

    @staticmethod
    def to_driver_basis(basis: str) -> str:
        """Converts basis to a driver acceptable basis.

        Args:
            basis: The basis set to be used.

        Returns:
            A driver acceptable basis.
        """
        if basis == "sto3g":
            return "sto-3g"
        return basis

    @staticmethod
    def check_method_supported(method: MethodType) -> None:
        """Checks that Gaussian supports this method.

        Args:
            method: the SCF method type.

        Raises:
            UnsupportMethodError: If the method is not supported.
        """
        if method not in [MethodType.RHF]:
            raise UnsupportMethodError(f"Invalid Molcas method {method.value}.")

    def run(self) -> ElectronicStructureProblem:
        cfg = self._config

        logger.debug(
            "User supplied configuration raw: '%s'",
            cfg.replace("\r", "\\r").replace("\n", "\\n"),
        )

        logger.debug("User supplied configuration\n%s", cfg)

        run_directory = os.getcwd()
        file_fd, input_file = tempfile.mkstemp(suffix=".inp", prefix="molcasf", dir=run_directory)
        os.close(file_fd)
        with open(input_file, "w", encoding="utf8") as stream:
            stream.write(cfg)
        
        output_file = input_file.replace(".inp", ".out") 
        MolcasDriver._run_pymolcas(input_file, output_file)
        if logger.isEnabledFor(logging.DEBUG):
            with open(output_file, "r", encoding="utf8") as file:
                logger.debug("OpenMolcas output file:\n%s", file.read())

        file_name = input_file.rpartition('.')[0]        
        self._qcschemadata = MolcasDriver._parse_molcas_files(file_name)
        try:
            base_name = os.path.basename(output_file).rpartition('.')[0]
            for local_file in os.listdir(run_directory):
                if local_file.startswith(f"{base_name}"):
                    os.remove(run_directory + "/" + local_file)
                elif local_file.startswith("xmldump"):
                    os.remove(local_file)
        except Exception:  # pylint: disable=broad-except
            logger.warning("Failed to remove OpenMolcas files starting with %s", base_name)

        return self.to_problem()

    def to_qcschema(self, *, include_dipole: bool = False) -> QCSchema:
        return MolcasDriver._to_qcschema(self._qcschemadata, include_dipole=include_dipole)

    def to_problem(
        self,
        *,
        basis: ElectronicBasis = ElectronicBasis.MO,
        include_dipole: bool = False,
    ) -> ElectronicStructureProblem:
        return qcschema_to_problem(
            self.to_qcschema(include_dipole=include_dipole),
            basis=basis,
            include_dipole=include_dipole,
        )

    @staticmethod
    def _run_pymolcas(input_file: str, output_file: str) -> str:
        process = None
        try:
            with subprocess.Popen(
                [_optionals.MOLCAS, input_file, "-o", output_file],
                stdout=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                stdout, _ = process.communicate()
                process.wait()
        except Exception as ex:
            if process is not None:
                process.kill()

            raise QiskitNatureError(f"{_optionals.MOLCAS_DESC} run has failed") from ex

        if process.returncode != 0:
            errmsg = ""
            if stdout is not None:
                lines = stdout.splitlines()
                for line in lines:
                    logger.error(line)
                    errmsg += line + "\n"
            raise QiskitNatureError(
                f"{_optionals.MOLCAS_DESC} process return code {process.returncode}: {errmsg}"
            )
    
    @staticmethod
    def _parse_molcas_files(fname: str, *, include_dipole: bool = False) -> QCSchema:
        from .utils import parse_molden, parse_output
        data = _QCSchemaData()

        if not os.path.exists(f"{fname}.FciDmp"):
            raise QiskitNatureError(f"Missing file {fname}.FciDmp")
        if not os.path.exists(f"{fname}.scf.molden"):
            raise QiskitNatureError(f"Missing file {fname}.scf.molden")
        if not os.path.exists(f"{fname}.out"):
            raise QiskitNatureError(f"Missing file {fname}.out")
        
        # TODO: Support adding results from (UHF, ROHF) calculation into data
        fcidump = FCIDump.from_file(f"{fname}.FciDmp")
        data.hij = None     # TODO: transform MO integrals to AO integrals
        data.hij_b = None
        data.hij_mo = fcidump.hij
        data.hij_mo_b = fcidump.hij_b
        # if _has_B:
        #     _q_h1b.transform(_q_hf_wavefn.Cb())
        #     data.hij_mo_b = np.asarray(_q_h1b)

        # # TODO: add support for symmetry-reduced integrals
        # data.eri = np.asarray(_q_mints.ao_eri())
        data.eri_mo = fcidump.hijkl
        data.eri_mo_ba = None
        data.eri_mo_bb = None
        # if _has_B:
        #     data.eri_mo_bb = np.asarray(_q_mints.mo_eri(_q_hf_wavefn.Cb(), _q_hf_wavefn.Cb(),
        #                                                 _q_hf_wavefn.Cb(), _q_hf_wavefn.Cb()))
        #     data.eri_mo_ba = np.asarray(_q_mints.mo_eri(_q_hf_wavefn.Cb(), _q_hf_wavefn.Cb(),
        #                                                 _q_hf_wavefn.Ca(), _q_hf_wavefn.Ca()))

        molden_data = parse_molden(f"{fname}.scf.molden")
        data.overlap = None
        data.mo_coeff = np.asarray([mo["coeffs"] for mo in molden_data["orbitals"]])
        data.mo_coeff_b = None
        data.mo_energy = np.asarray([mo["energy"] for mo in molden_data["orbitals"]])
        data.mo_energy_b = None
        data.mo_occ = np.asarray([mo["occup"] for mo in molden_data["orbitals"]])
        data.mo_occ_b = None
        data.symbols = []
        data.coords  = np.empty([molden_data["natoms"], 3])
        for _n, atom in enumerate(molden_data["atoms"]):
            data.symbols.append(PERIODIC_TABLE[atom["atnum"]])
            data.coords[_n][0] = atom['x']
            data.coords[_n][1] = atom['y']
            data.coords[_n][2] = atom['z']        
        data.coords = data.coords.flatten()

        logdata = parse_output(f"{fname}.out")
        data.e_nuc = fcidump.constant_energy
        data.e_ref = logdata["etotal"]
        data.multiplicity = fcidump.multiplicity
        data.charge = None      # TODO: add mass and charge information
        data.masses = None     
        data.method = None      # TODO: determine method from output
        data.basis = logdata["basis"]
        data.creator = "OpenMolcas"
        data.version = logdata["version"]
        data.routine = None
        data.nbasis = None
        data.nmo = data.mo_coeff.shape[0]
        data.nalpha = (fcidump.num_electrons + data.multiplicity - 1) // 2
        data.nbeta = (fcidump.num_electrons - data.multiplicity + 1) // 2
        data.keywords = None

        # TODO: add dipole moment information

        return data
