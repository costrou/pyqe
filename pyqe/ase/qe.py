"""
ASE Calculator for Quantum Espresso

 - potential_energy
 - force
 - stress
"""

from ase.calculators.calculator import Calculator, all_changes, all_properties
from pyqe.espresso import PWBase
import numpy as np

class QE(Calculator):
    """Quantum Espresso ASE Calculator
    
    """
    implemented_properties = ['energy', 'forces', 'stress']

    default_parameters = {
        'calculation': 'scf',
        'convergence': {'energy': 1E-6}
    }

    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                 atoms=None, **kwargs):
        """Accepted Keypairs
        calculation:
            str - only scf[default] supported at the moment
        ecutwfc: 
            plane wave cuttoff
        nbands:
            number of bands for calculation (see pw.x -> nbnd)
        usesymm:
            whether or not to use symmetry in calculation (True/False)
        maxiter:
            maximum number of iterations in an scf step
        convergence:
            {'energy', <value>} - only one implemented
        kpts:
            (1, 1, 1) - gamma point
            (n1, n2, n3) - Morstead Packing
            [[k1, k2, k3] ... ] - List of kpoints
        prefix: (self.label is meaningless at moment)
            str - string to append to output files
        pseudo:
            {'atom symbol': 'name of pseudo potential', ...} - dictionary of pseudo per atom
        outdir:
            str - relative or absolute path to output directory.
                  Will create it if it does not exist.
        pseudo_dir:
            str - relative or absolute path to pseudo directory.
        occupations:
            str - 'smearing', 'tetrahedra', 'fixed', 'from_input'
        input_dft:
            str - functional to use
        fft_mesh:
            [nr1, nr2, nr3] - fft mesh to use for calculation
        keypairs:
            way to initialize via the base class PWBase. 
            (DO NOT SET CELL OR ATOM POSITIONS via PWBase this is done 
            automagically via ASE Atoms)
        debug:
            if True will print input and output QE files

        AUTOMATICALLY SET KEYPAIRS:
        Set so we can always extract stress and forces:
            control.tstress=True
            control.tprnfor=True
        From atoms object:
            CARD 'CELL PARAMETERS (angstroms)' from atoms.cell
            CARD 'ATOMIC POSITIONS (angstroms)' from atoms[i]
            system.nat from number of unique atoms in atoms[i]
            system.ntyp from len(atoms)
        Thus DO NOT SET:
            A, B, C, cosAB, cosBC, cosAC, celldm(1-6)
            or any of these previously mentioned. Be smart
        """
        Calculator.__init__(self, restart, ignore_bad_restart_file, label, atoms, **kwargs)
        self._pw = None

    def set(self, **kwargs):
        """Sets parameters list and determines if calculation needs to be reset

        *for now it is reset if any parameters change*
        """
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()

    def initialize(self):
        """Updates PWBase to reflect self.parameters and atoms to be ready for a run

        """
        self._pw = PWBase()

        # Setup Unit Cell
        cell = self.atoms.cell
        self._pw.system.add_keypair(('ibrav', 0))
        self._pw.cell_parameters.add_lattice_vec(cell[0], cell[1], cell[2], option='angstrom')

        # Setup Atom Positions
        self._pw.system.add_keypair(('nat', len(self.atoms)))
        for atom in self.atoms:
            self._pw.atomic_positions.add_atom_position(atom.symbol, atom.position, option='angstrom')
        
        # Setup Atom Psuedo Potentials
        from itertools import groupby
        for symbol, atoms in groupby(self.atoms, lambda atom: atom.symbol):
            pseudo_file = self.parameters['pseudo'].get(symbol)
            if pseudo_file:
                self._pw.atomic_species.add_atom_type(symbol, list(atoms)[0].mass, pseudo_file)
            else:
                error_str = "Must provide pseudo potential for all atom types {0}"
                raise Exception(error_str.format(symbol))
        self._pw.system.add_keypair(('ntyp', len(self._pw.atomic_species.atoms)))

        
        # Set User Input Parameters [easy ones first, hard last]
        if self.parameters.get('calculation'):
            self._pw.control.add_keypair(('calculation', self.parameters.get('calculation')))

        if self.parameters.get('ecutwfc'):
            self._pw.system.add_keypair(('ecutwfc', self.parameters['ecutwfc']))

        if self.parameters.get('nbands'):
            self._pw.system.add_keypair(('nbnd', self.parameters['nbands']))

        if self.parameters.get('usesymm'):
            self._pw.system.add_keypair(('nosym', not self.parameters['usesymm']))

        if self.parameters.get('maxiter'):
            self._pw.electrons.add_keypair(('electron_maxstep', self.parameters['maxiter']))

        if self.parameters.get('prefix'):
            self._pw.control.add_keypair(('prefix', self.parameters['prefix']))

        if self.parameters.get('outdir'):
            self._pw.control.add_keypair(('outdir', self.parameters['outdir']))

        if self.parameters.get('pseudo_dir'):
            self._pw.control.add_keypair(('pseudo_dir', self.parameters['pseudo_dir']))

        if self.parameters.get('occupations'):
            self._pw.system.add_keypair(('occupations', self.parameters['occupations']))

        if self.parameters.get('input_dft'):
            self._pw.system.add_keypair(('input_dft', self.parameters['input_dft']))

        if self.parameters.get('fft_mesh'):
            fft_mesh = self.parameters['fft_mesh']
            if len(fft_mesh) == 3 and isinstance(fft_mesh[0], int):
                self._pw.system.add_keypairs({'nr1': fft_mesh[0],
                                               'nr2': fft_mesh[1],
                                               'nr3': fft_mesh[2]})
            else:
                error_str = 'FFT MESH = [nr1, nr2, nr3] each is integer'
                raise Exception(error_str)

        if self.parameters.get('convergence'):
            conv = self.parameters['convergence']
            if isinstance(conv, dict):
                if conv.get('energy'):
                    self._pw.electrons.add_keypair(('conv_thr', conv['energy']))
                else:
                    raise Exception('Only energy convergence implemented currently')
            else:
                raise Exception("Unknown convergence format declared (dict required)")

        if self.parameters.get('kpts'):
            kpts = self.parameters['kpts']
            if isinstance(kpts, dict):
                if kpts.get('density'):
                    from ase.calculator.calculator import kptdensity2monkhorstpack
                    self._pw.k_points.from_list(kptdensity2monkhorstpack(
                        kpts['density'], kpts.get('even', False)))
                elif kpts.get('size'):
                    self._pw.k_points.from_monkhorst_pack(
                        kpts['size'], kpts.get('offset', [0, 0, 0]))
            elif isinstance(kpts, list):
                self._pw.k_points.from_list(kpts)
            else:
                raise Exception("Unknown kpts format declared")

        # Add all PWBase initializer keypairs [READ VALUES TO NOT SET]
        if self.parameters.get('keypairs'):
            self._pw.add_keypairs_to_namelist(self.parameters['keypairs'])

        # Always calculate stress and strain
        self._pw.control.add_keypairs({
            'tstress': True,
            'tprnfor': True})

        self._pw.validate()


    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """Determines if a calculation needs to be preformed

        """
        if system_changes or not self.results:
            if atoms:
                self.set_atoms(atoms)

            self.initialize()

            # Run pw.x
            if self.parameters.get('debug'):
                results = self._pw.run(infile="in", outfile="out", errfile="err")
            else:
                results = self._pw.run()

            # Extract Results and convert to appropriate units
            from ase.units import Bohr, Ry

            energy = results['calculation']['total energy'] / Ry

            forces = np.array([_[2] for _ in results['calculation']['forces']]) * Bohr/Ry

            
            stress = np.array(results['calculation']['stress']) * Bohr**3/Ry
            # xx, yy, zz, yz, xz, xy
            stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                               stress[1, 2], stress[0, 2], stress[0, 1]])

            self.results = {'energy': energy,
                            'forces': forces,
                            'stress': stress}

        return self.results[properties[0]]
