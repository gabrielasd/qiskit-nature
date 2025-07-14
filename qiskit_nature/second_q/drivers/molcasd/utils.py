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

""" OpenMolcas Utility Methods """

import re

from qiskit_nature.constants import PERIODIC_TABLE
from ..electronic_structure_driver import MethodType


def make_molcas_cofig(
        geom: str | list[str],
        charge: int,
        basis: str,
        nbasis: int,
        method: str, 
        units: str,
        wfn_tol: float = 5e-2
    ) -> str:
    """Make a Molcas configuration string from the given parameters.

    Args:
        geom: The elements and coordinates of all atoms in the system.
        charge: The charge of the molecule.
        basis: A basis set name as recognized by OpenMolcas.
        nbasis: The number of spatial orbitals. This parameter is required to generate and FCIDUMP file.
        method: The SCF method type to be used for the OpenMolcas calculation. For now the MolcasDriver
        only supports RHF and CASSCF methods.
        units: Denotes the unit of coordinates. Valid values are "Angstrom" and "Bohr".
        wfn_tol: This keyword is used to specify the threshold for printing the coefficients and Slater 
        determinant expansion of a CASSCF wavefunction to a ``VecDet`` file.
    """

    if isinstance(geom, list):
        geom = "\n".join(geom)

    atoms = geom.split('\n')
    natom = len(atoms)
    nelec = sum([PERIODIC_TABLE.index(atom[0][0]) for atom in atoms])
    nelec -= charge
    
    # OpenMolcas input instructions
    cfg1 = "&GATEWAY\n"
    cfg2 = f"coord\n{natom}\n{units}\n"
    cfg2 += f"{geom}\nbasis={basis}\nGroup=Nosym\n"
    cfg3 = f"&SEWARD\n&SCF\n"
    cfg4 = f"&RASSCF\nSpin= 1; Nactel= {nelec} 0 0; Inactive= 0; Ras2= {nbasis}\n"
    cfg4 += f"DMPO\n\n"

    if method == "casscf":
        cfg4 += f"&RASSCF\nprwf = {wfn_tol}\nPRSD\n\n"
    
    return cfg1 + cfg2 + cfg3 + cfg4


def parse_molden(molden: str):
    """Parse a Molden file to extract atomic and molecular orbital data.

    Parameters
    ----------
    fname : str
        The path to the Molden file.

    Returns
    -------
    data : dict
        A dictionary containing the following keys:
        - 'natoms': Number of atoms in the system.
        - 'atoms': A list of dictionaries, each containing the following keys:
            - 'label': Atomic symbol.
            - 'atnum': Atomic number.
            - 'x', 'y', 'z': Cartesian coordinates of the atom.
        - 'orbitals': A list of dictionaries, each containing the following keys:
            - 'energy': Energy of the molecular orbital.
            - 'occup': Occupation number of the molecular orbital.
            - 'coeffs': Coefficients of the molecular orbital in the basis set.
    """
    # TODO: add charge and basis set attributes
    data = {
        "natoms": 0,
        "atoms": [],
        "orbitals": []
    }

    with open(molden, 'r') as f:
        lines = f.readlines()

    # section = None
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Get number of atoms
        if line.startswith('[N_Atoms]'):
            # section = 'natoms'
            i += 1
            data['natoms'] = int(lines[i].strip())

        # Get atomic data: Symbol, atomic number and coordinates
        elif line.startswith('[Atoms]'):
            # section = 'atoms'
            i += 1
            while i < len(lines) and not lines[i].startswith('['):
                parts = lines[i].split()
                if len(parts) >= 6:
                    atom = {
                        'label': parts[0],
                        'atnum': int(parts[2]),
                        'x': float(parts[3]),
                        'y': float(parts[4]),
                        'z': float(parts[5])
                    }
                    data['atoms'].append(atom)
                i += 1
            continue  # skip normal i += 1 to avoid missing next section header
        
        # Get molecular orbitals: energies, occupations and coefficients
        # TODO: check that L81-L101 are correct
        elif line.startswith('[MO]'):
            # section = 'mo'
            i += 1
            orbital = {}
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith('Sym='):
                    if orbital:
                        data['orbitals'].append(orbital)
                    orbital = {'coeffs': []}
                elif line.startswith('Ene='):
                    orbital['energy'] = float(line.split('=')[1])
                elif line.startswith('Occup='):
                    orbital['occup'] = float(line.split('=')[1])
                elif re.match(r'^\d+', line):  # coefficient line
                    parts = line.split()
                    if len(parts) == 2:
                        orbital['coeffs'].append(float(parts[1]))
                i += 1
            if orbital:
                data['orbitals'].append(orbital)

        else:
            i += 1

    return data


def parse_output(out: str):
    """Parse the output of a Molcas calculation to extract relevant information."""
    ver = None
    basis = None
    nbasis = None
    etotal = None
    
    with open(out, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if "version:" in line:
            ver = line.split()[-1].strip()
        elif "basis=" in line:
            basis = line.split("=")[-1].strip()
        elif "Basis functions" in line:
            nbasis = int(line.split()[-1].strip())
        elif line.startswith("::") and "Total SCF energy" in line:
            etotal = float(line.split()[-1].strip())
        
    return {
        "version": ver,
        "basis": basis,
        "nbasis": nbasis,
        "etotal": etotal,
    }
