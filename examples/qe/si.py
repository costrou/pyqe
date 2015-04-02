from pyqe.espresso import PWBase

#Data Setup
qe_namelist = {
    'control': {
        'calculation': 'scf',
        'prefix': 'silicon',
        'tstress': True,
        'tprnfor': True,
        'pseudo_dir': './',
        'outdir': './',
    },
    'system': {
        'ibrav': 2,
        'celldm(1)': 10.20,
        'nat': 2,
        'ntyp': 1,
        'ecutwfc': 18.0,
    },
    "electrons": {
        'diagonalization': 'davidson',
        'mixing_mode': 'plain',
        'mixing_beta': 0.7,
        'conv_thr': 1E-8,
    }
}

atom_types = [
    ['Si', 28.086, 'Si.pz-vbc.UPF']
]

atom_positions = [
    ['Si', [0.0, 0.0, 0.0]],
    ['Si', [0.25, 0.25, 0.25]]
]

kpoints = (4, 4, 4)

# Create QE Run
qe = PWBase()
qe.add_keypairs_to_namelist(qe_namelist)

for symbol, mass, pseudo in atom_types:
    qe.atomic_species.add_atom_type(symbol, mass, pseudo)

for symbol, position in atom_positions:
    qe.atomic_positions.add_atom_position(symbol, position)

qe.k_points.from_monkhorst_pack([8, 8, 8], [0, 0, 0])

# Validate and Run quantum espresso
qe.validate()
results = qe.run()
