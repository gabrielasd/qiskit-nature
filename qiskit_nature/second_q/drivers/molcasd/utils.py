import re


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
        elif line.startswith("::"):
            etotal = float(line.split()[-1].strip())
        
    return {
        "version": ver,
        "basis": basis,
        "nbasis": nbasis,
        "etotal": etotal,
    }
