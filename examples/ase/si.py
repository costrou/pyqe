from pyqe.ase.qe import QE
from ase.lattice.spacegroup import crystal

a = 5.43
material = crystal(('Si'), basis=[0.0, 0.0, 0.0], spacegroup=227,
                   cellpar=[a, a, a, 90, 90, 90])

material.calc = QE(calculation='scf',
                   prefix='silicon',
                   ecutwfc=16.0,
                   kpts={'size':[4, 4, 4]}, 
                   pseudo={'Si': 'Si.pz-vbc.UPF'}, 
                   occupations='fixed',
                   pseudo_dir='./',
                   outdir='./')

print("Potential Energy: {0}".format(material.get_potential_energy()))

print("\nForces on Atoms")
for i, f in material.get_forces():
    print("{0} {1}: <{2} {3} {4}>".format(material[i].symbol, i, f[0], f[1], f[2]))

s = material.get_stress()
print("\nStress")
print("xx {0}, yy {1}, zz {2}, yz {3}, xz {4}, xy {5}".format(s[0], s[1], s[2], s[3], s[4], x[5]))
