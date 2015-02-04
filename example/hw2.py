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
        'outdir': './saves/',
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
    lattice_params = np.linspace(5.0, 13.0, 10)
    #ecutwfc_values = [15.0, 20.0, 25.0, 30.0]
    ecutwfc_values = [70.0]

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
                     celldm,
                     results['total energy'][0]])


    # Write results to files
    outfilename = "results.txt"
    with open(outfilename, "w") as outfile:
        line_str = "{0} | {1} | {2} | {3} | {4}\n"
        outfile.write(line_str.format(
            "prefix", "kpoint", "ecutwfc", "celldm", "total energy"))
        for row in data:
            outfile.write(line_str.format(
                row[0], row[1], row[2], row[3], row[4]))

    # Plot results
    plot_data = []
    ecutwfc = ecutwfc_values[0]
    for kpoint in kpoints:
        kpoint_data = [d for d in data if d[1] == kpoint and d[2] == ecutwfc]
        kpoint_data.sort(key=lambda d: d[3])
        kpoint_x = [d[3] for d in kpoint_data]
        kpoint_y = [d[4] for d in kpoint_data]
        plot_data.append([kpoint_x, kpoint_y])

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    title_str = "Energy vs. Lattice Parameter ecutwfc = {0}"
    ax.set_title(title_str.format(ecutwfc))
    ax.set_xlabel("Lattice Parameter")
    ax.set_ylabel("Total Energy")
    for kpoint, [celldm, tot_energy] in zip(kpoints, plot_data):
        ax.plot(celldm, tot_energy, label="{0}".format(kpoint))
    ax.legend()
    fig.savefig("hw2.png")
    fig.show()
