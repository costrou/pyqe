#!/usr/bin/python3
"""A program to calculate the system energies of FCC CU with differnt
cutoff energies and kpoints. Changing these parameters affects the
accuracy of the calculations
"""
import sys
sys.path.append("/home/costrouc/work/projects/python-qe/")
from pyqe import QE

import numpy as np

# Burch Murnahan EOS
def bmeos(V0, E0, B0, BP, V):
    eta = np.power(V0 / V, 1.0/3.0)
    term1 = (eta**2 - 1)**3 * BP
    term2 = (eta**2 - 1)**2 * (6.0 - 4*eta**2)
    return E0 + 9*B0*V0/16.0 * (term1 + term2)


# Fitting to Birch–Murnaghan EOS(equation of state)
def solve_bmeos(volume, energy):
    volume = np.array(volume)
    energy = np.array(energy)

    a, b, c = np.polyfit(volume, energy, 2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    BP = 3.5
    init_guess = (V0, E0, B0, BP)

    def error_func(param, V, E):
        V0, E0, B0, BP = param
        return E - bmeos(V0, E0, B0, BP, V)

    from scipy.optimize import leastsq
    return leastsq(error_func, init_guess, args=(volume, energy))


if __name__ == "__main__":

    qe_namelist = {
        "control": {
            'calculation': 'scf',
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

    kpoints = [[8, 8, 8]]
    lattice_params = np.linspace(6.00, 8.00, 10)
    ecutwfc_values = [20.0, 40.0, 60.0, 80.0]

    # Convert a.u.^3 to angstroms

    bhor = 0.529177249 # Angstroms/bhor
    Ry = 13.605698066 # eV/Ry
    AeVtoGPa = 160.2176487 #GPa/(eV/A^3)

    data = []
    from itertools import product
    for i, [kpoint, celldm, ecutwfc] in enumerate(product(
            kpoints, lattice_params, ecutwfc_values)):
        print("iter {0} celldm{1}".format(i, celldm))

        prefix = "cu.{0}".format(i)
        qe_namelist['control']['prefix'] = prefix
        qe_namelist['system']['ecutwfc'] = ecutwfc
        qe_namelist['system']['celldm(1)'] = celldm

        qe = QE(qe_namelist)

        qe.atomic_species.add_atom_type('Cu', 63.55, 'Cu.pbe-mt_fhi.UPF')
        qe.atomic_positions.add_atom_position('Cu', [0.0, 0.0, 0.0])
        qe.k_points.automatic(kpoint, [0, 0, 0])

        qe.validate()

        results = qe.run(infile="in", outfile="out", errfile="err")

        volume = results['volume'] * bhor**3
        energy = results['total energy'] * Ry

        data.append(["cu.{0}".format(i),
                     kpoint,
                     ecutwfc,
                     volume,
                     energy])


    # Write results to files
    import json
    outfilename = "results.json"
    json.dump({"header": ["prefix", "kpoint", "ecutwfc", "volume", "total energy"],
               "data": data},
              open(outfilename, "w"), indent=4)

    # Plot data
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for kpoint, ecutwfc in product(kpoints, ecutwfc_values):
        volume, energy = zip(*[[row[3], row[4]] for row in data if \
                               row[1] == kpoint and row[2] == ecutwfc])

        # Data Points
        label_str = "kpoint: {0}".format(kpoint)
        label_str += " ecutwfc: {0}".format(ecutwfc)
        ax.plot(volume, energy, '.', label=label_str)

        # Fitting
        result, iteration = solve_bmeos(volume, energy)

        x_pts = np.linspace(min(volume), max(volume), 1000)
        y_pts = bmeos(result[0], result[1], result[2], result[3], x_pts)
        
        bulk_modulus = result[2] * AeVtoGPa
        label_str = "Bulk Modulus: {0:.2f} [GPa]"
        ax.plot(x_pts, y_pts, '-', label=label_str.format(bulk_modulus))

    ax.legend()
    ax.set_title("Energy vs.Volume Cu fcc")
    ax.set_xlabel("Volume [A^3]")
    ax.set_ylabel("Total Energy [eV]")
    fig.set_size_inches(20.0, 12.0)
    fig.savefig("hw2.png", dpi=200)

