"""A program to calculate the system energies of FCC CU with differnt
cutoff energies and kpoints. Changing these parameters affects the
accuracy of the calculations
"""

def create_input():
    """
    As a rough starting template we will use dictionaries
    """
    import sys
    sys.path.append("../")
    from espresso import QE
    qe = QE()
    qe.control.addKeypairs({
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'pseudo_dir': './',
        'prefix': 'cu'
    })
    qe.system.addKeypairs({
        'ibrav': 2,
        'celldm(i)': 6.73, #Yeah I had to lie for now
        'nat': 1,
        'ntyp': 1,
        'ecutwfc': 25.0,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.02
    })
    qe.electrons.update({
        'diagonalization': 'david',
    })

    qe.atomic_species.add_atom_type('Cu', 63.55, 'Cu.blyp-d-hgh.UPF')
    qe.atomic_positions.add_atom_position('Cu', [0.0, 0.0, 0.0])
    qe.k_points.automatic([4, 4, 4], [0, 0, 0])

    return qe


if __name__ == "__main__":
    qe = create_input()
    qe.validate()
    qe.to_file("hello.in")
    qe.run()
