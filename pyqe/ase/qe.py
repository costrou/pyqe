"""
ASE Calculator for Quantum Espresso

 - potential_energy
 - force
 - stress
"""

from ase.calculators.calculator import Calculator, equal
from ase.units import Bohr, Ry, Hartree

from pyqe.espresso import PWBase
import numpy as np

def calculation(property_name):
    def calculation_decorator(property_function):
        def calculation_wrapper(self, atoms=None, *args, **kwargs):
            if property_name not in self.implemented_properties:
                raise NotImplementedError(property_name)

            if self.atoms is None and atoms is None:
                raise Exception("Atoms object required for calculation")

            if self._calculation_required(atoms, property_name):
                if atoms is not None:
                    self.set_atoms(atoms)

                self._initialize()

                results = self._calculate(property_name)
                self._set_results(results)

            return property_function(self, atoms, *args, **kwargs)
        return calculation_wrapper
    return calculation_decorator


class QE(Calculator):
    """Quantum Espresso ASE Calculator
    
    """
    implemented_properties = ['energy', 'forces', 'stress', 
                              'ibz_kpoint_eigenvalues', 'ibz_kpoints_position', 'ibz_kpoints_weight', 
                              'nspins', 'nbands', 'xc_functional', 'fermi_energy', 'dos']

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
            plane wave cuttoff [eV] notice not in Ry!
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
        # Read parameters from file
        if 'parameters' in kwargs:
            filename = kwargs.pop('parameters')
            parameters = Parameters.read(filename)
            parameters.update(kwargs)
            kwargs = parameters

        changed_parameters = {}

        for key, value in kwargs.items():
            oldvalue = self.parameters.get(key)
            if key not in self.parameters or not equal(value, oldvalue):
                changed_parameters[key] = value
                self.parameters[key] = value

        if changed_parameters:
            self.reset()

        return changed_parameters

    def reset(self):
        """Clear all information from old calculation."""
        self.results = {}

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()

    def _initialize(self):
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
            self._pw.system.add_keypair(('ecutwfc', self.parameters['ecutwfc'] / Ry))

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

        if self.parameters.get('kpts') is not None:
            kpts = self.parameters['kpts']

            if isinstance(kpts, dict):
                if kpts.get('density'):
                    from ase.calculator.calculator import kptdensity2monkhorstpack
                    self._pw.k_points.from_list(kptdensity2monkhorstpack(
                        kpts['density'], kpts.get('even', False)))
                elif kpts.get('size'):
                    self._pw.k_points.from_monkhorst_pack(
                        kpts['size'], kpts.get('offset', [0, 0, 0]))
                else:
                    raise Exception("Unknown kpts format declared")
            elif hasattr(kpts, "__iter__"):
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

    def _calculate(self, property_name):
        """Preforms calculation.

        """
        if self.parameters.get('debug'):
            return self._pw.run(infile="in", outfile="out", errfile="err")
        return self._pw.run()

    def _calculation_required(self, atoms, property_name):
        # TODO need to do better check on results (using cached idea)
        if self.results == {}:
            return True

        if self.atoms and atoms:
            system_changes = self.check_state(atoms)
            if system_changes:
                return True
        
        return False

    def _set_results(self, results):
        # Extract Results from outfile and convert to appropriate units


        # Updates from output file (some values may not be output!!)
        if self.parameters.get('calculation') != "bands":
            energy = results['calculation']['total energy'] * Ry

            forces = np.array([_[2] for _ in results['calculation']['forces']]) * Ry / Bohr

            stress = np.array(results['calculation']['stress']) * Ry / (Bohr**3)
            # xx, yy, zz, yz, xz, xy
            stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                               stress[1, 2], stress[0, 2], stress[0, 1]])

            self.results.update(
                {'energy': energy,
                 'forces': forces,
                 'stress': stress})

        # Update from data file (guarenteed to be in output)
        fermi_energy = results['data-file']['band-structure-info']['fermi-energy'] * Hartree

        for kpt in results['data-file']['kpoints']:
            kpt['eigenvalues'] = np.array(kpt['eigenvalues']) * Hartree

        self.results.update(
            {'fermi-energy': fermi_energy,
             'xc-functional': results['data-file']['exchange-correlation'],
             'nspins': results['data-file']['band-structure-info']['number spin-components'],
             'nbands': results['data-file']['band-structure-info']['number bands'],
             'charge-density': results['data-file']['charge-density'],
             'ibz-kpoints': results['data-file']['kpoints']})


    @calculation("energy")
    def get_potential_energy(self, atoms=None):
        return self.results['energy']

    @calculation("forces")
    def get_forces(self, atoms=None):
        return self.results['forces']

    @calculation("stress")
    def get_stress(self, atoms=None):
        return self.results['stress']

    @calculation("ibz_kpoint_eigenvalues")
    def get_eigenvalues(self, atoms=None, kpt=0, spin=0):
        if spin == 1:
            raise Exception("Spin systems not implemented yet!")
        return self.results['ibz-kpoints'][kpt]['eigenvalues']

    @calculation("ibz_kpoints_position")
    def get_ibz_k_points(self, atoms=None):
        return np.array([kpoint['coordinate'] for kpoint in self.results['ibz-kpoints']])

    @calculation("ibz_kpoints_weight")
    def get_k_point_weights(self, atoms=None):
        return np.array([kpoint['weight'] for kpoint in self.results['ibz-kpoints']])

    @calculation("nspins")
    def get_number_of_spins(self, atoms=None):
        return self.results['nspins']
    
    @calculation("nbands")
    def get_number_of_bands(self, atoms=None):
        return self.results['nbands']

    @calculation("xc_functional")
    def get_xc_functional(self, atoms=None):
        return self.results['xc_functional']

    @calculation("fermi_energy")
    def get_fermi_level(self, atoms=None):
        return self.results['fermi-energy']
    
    @calculation("pseudo_density")
    def get_pseudo_density(self, atoms=None, spin=None, pad=True):
        if spin != None:
            raise Exception("Spin systems not implemented yet!")
        return self.results['charge-density']

    @calculation("dos")
    def get_dos(self, atoms=None, spin=None, width=0.1, npts=201):
        from ase.dft import DOS
        dos = DOS(self, width=width, npts=npts)
        return dos.get_energies(), dos.get_dos(spin)
        
