#!/usr/bin/python3
"""A program to calculate the system energies of FCC CU with differnt
cutoff energies and kpoints. Changing these parameters affects the
accuracy of the calculations
"""
import sys
sys.path.append("/home/costrouc/work/projects/python-qe/")
from pyqe import QE

import numpy as np

def create_input():
    """
    As a rough starting template we will use dictionaries
    """
    qe = QE()
    qe.control.add_keypairs({
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'pseudo_dir': './',
        'prefix': 'cu'
    })
    qe.system.add_keypairs({
        'ibrav': 2,
        'celldm(1)': 6.73, 
        'nat': 1,
        'ntyp': 1,
        'ecutwfc': 25.0,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.02
    })
    qe.electrons.add_keypairs({
        'diagonalization': 'david',
    })
    qe.atomic_species.add_atom_type('Cu', 63.55, 'Cu.blyp-d-hgh.UPF')
    qe.atomic_positions.add_atom_position('Cu', [0.0, 0.0, 0.0])
    qe.k_points.automatic([4, 4, 4], [0, 0, 0])

    return qe


if __name__ == "__main__":


    kpoints = [[1, 1, 1], [2, 2, 2], [4, 4, 4], [8, 8, 8]]
    lattice_params = np.linspace(6.25, 7.25, 10)
    ecutwfc_values = [15.0, 20.0, 25.0, 30.0]

    i = 0
    data = []
    qe = create_input()
    for kpoint in kpoints:
        for lattice_param in lattice_params:
            for ecutwfc_value in ecutwfc_values:

                prefix = "cu.{0}".format(i)
                qe.control.add_keypairs({"prefix": prefix})
                qe.system.add_keypairs({"ecutwfc": ecutwfc_value,
                                        "celldm(1)": lattice_param})
                qe.k_points.automatic(kpoint, [0, 0, 0])
                qe.validate()
                results = qe.run(infile="in", outfile="out", errfile="err")

                data.append(["cu.{0}".format(i),
                             ecutwfc_value,
                             lattice_param,
                             kpoint,
                             results['total energy']])
                i += 1

    print(data)
