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

    qe_namelists = {
        "control": {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'outdir': './saves/',
            'pseudo_dir': './',
            'prefix': 'cu'
        },
        "system": {
            'ibrav': 2,
            'celldm(1)': 6.73,
            'nat': 1,
            'ntyp': 1,
            'ecutwfc': 25.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02
        },
        "electrons": {
            'diagonalization': 'cg'
        }
    }

    qe = QE(qe_namelists)

    qe.atomic_species.add_atom_type('Cu', 63.55, 'Cu.pz-d-rrkjus.UPF')
    qe.atomic_positions.add_atom_position('Cu', [0.0, 0.0, 0.0])
    qe.k_points.automatic([4, 4, 4], [0, 0, 0])

    return qe

def bmeos(V0, E0, V, B0, Bprime):
    eta = np.power(V/V0, 1.0/3.0)
    term = 9*B0*V0/16.0 * (eta**2 - 1)**2
    term *= (6 + Bprime * (eta**2 - 1) - 4*eta**2)
    return E0 + term

# Fitting to Birch–Murnaghan EOS(equation of state)
def solve_bmeos(volume, energy):
    index, E0 = min(enumerate(energy), key=lambda row: row[1])
    V0 = volume[index]

    volume = np.array(volume)
    energy = np.array(energy)
    
    error_func = lambda param, V, E: bmeos(V0, E0, V, param[0], param[1]) - E
    init_guess = [80.0, 4.0]

    from scipy.optimize import leastsq
    return leastsq(error_func, init_guess, args=(volume, energy))


if __name__ == "__main__":


    kpoints = [[8, 8, 8]]
    lattice_params = np.linspace(5.5, 8.0, 30)
    ecutwfc_values = [30.0]

    data = []
    from itertools import product
    qe = create_input()
    for i, [kpoint, celldm, ecutwfc] in enumerate(product(
            kpoints, lattice_params, ecutwfc_values)):

        prefix = "cu.{0}".format(i)
        qe.control.add_keypairs({"prefix": prefix})
        qe.system.add_keypairs({"ecutwfc": ecutwfc,
                                "celldm(1)": celldm})
        qe.k_points.automatic(kpoint, [0, 0, 0])
        qe.validate()

        results = qe.run()

        data.append(["cu.{0}".format(i),
                     kpoint,
                     ecutwfc,
                     results['volume'],
                     results['total energy']])


    # Write results to files
    import json
    outfilename = "results.json"
    json.dump({"header": ["prefix", "kpoint", "ecutwfc", "volume", "total energy"],
               "data": data},
              open(outfilename, "w"), indent=4)

    # Plot data
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    ax.set_title("Energy vs.Volume Cu fcc")
    ax.set_xlabel("Volume a.u.^3")
    ax.set_ylabel("Total Energy")

    for kpoint, ecutwfc in product(kpoints, ecutwfc_values):
        filter_data = filter(lambda row: row[1] == kpoint and row[2] == ecutwfc,
                             data)
        volume, energy = zip(*[[row[3], row[4]] for row in filter_data])

        print(solve_bmeos(volume, energy))

        label_str = "kpoint: {0}".format(kpoint)
        label_str += " ecutwfc: {0}".format(ecutwfc)
        ax.plot(volume, energy, '.', label=label_str)

    ax.legend()
    fig.savefig("hw2.png")

