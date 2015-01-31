"""
All the Cards for Quantum Espresso

FULLY IMPLEMENTED:
  ATOMIC_SPECIES
  CELL_PARAMETERS

PARTIAL IMPLEMENTATION:
  K_POINTS
  ATOMIC_POSITIONS

NOT IMPLEMENTED:
  OCCUPATIONS
  CONSTRAINTS
  ATOMIC_FORCES
"""
import numpy as np

# Global Variables
QE_TAB = "   "

class AtomicSpecies:
    """
    Card: ATOMIC_SPECIES

    Atomic Species defines the atoms to be used in the dft calculation.
    For each line the symbol, mass, and pseudo potential to be used for the
    atom is supplied.
    """
    def __init__(self):
        self.name = "ATOMIC_SPECIES"
        self.atoms = []

    def num_atoms(self):
        """
        returns the numbers of atom types defined
        """
        return len(self.atoms)

    def add_atom_type(self, symbol, mass, pseudopot):
        """
        Add atom:
        symbol    (1-2) characters
        mass      double
        pseudopot path to pseudo potential file
        """
        if len(symbol) > 2:
            error_str = "Chemical Symbol 1-2 Characters [{0}]".format(symbol)
            raise Exception(error_str)

        self.atoms.append([symbol, mass, pseudopot])

    def __str__(self):
        atoms_str = "{0}\n".format(self.name)
        for atom in self.atoms:
            atoms_str += QE_TAB + "{0} {1} {2}\n".format(atom[0], atom[1], atom[2])
        return atoms_str

class AtomicPositions:
    """
    ATOMIC_POSITIONS Card

    options:
    alat       : atomic positions are in cartesian coordinates, in
                 units of the lattice parameter (either celldm(1)
                 or A). If no option is specified, 'alat' is assumed;
                 not specifying units is DEPRECATED and will no
                 longer be allowed in the future

    bohr       : atomic positions are in cartesian coordinate,
                 in atomic units (i.e. Bohr radii)

    angstrom   : atomic positions are in cartesian coordinates,
                 in Angstrom

    crystal    : atomic positions are in crystal coordinates, i.e.
                 in relative coordinates of the primitive lattice
                 vectors as defined either in card CELL_PARAMETERS
                 or via the ibrav + celldm / a,b,c... variables

    crystal_sg : atomic positions are in crystal coordinates, i.e.
                 in relative coordinates of the primitive lattice.
                 This option differs from the previous one because
                 in this case only the symmetry inequivalent atoms
                 are given. The variable space_group must indicate
                 the space group number used to find the symmetry
                 equivalent atoms. The other variables that control
                 this option are uniqueb, origin_choice, and
                 rhombohedral.
    """
    options = ["alat", "bhor", "crystal", "crystal_sg", "angstrom"]

    def __init__(self):
        self.name = "ATOMIC_POSITIONS"
        self.option = None
        self.atom_positions = []

    def num_atom_positions(self):
        """
        returns the numbers of atom positions defined
        """
        return len(self.atom_positions)

    def add_atom_position(self, symbol, position, option="alat"):
        """
        Add atom position:
        symbol    (1-2) characters
        position  [x, y, z]
        """
        if len(symbol) > 2:
            error_str = "Chemical Symbol 1-2 Characters [{0}]".format(symbol)
            raise Exception(error_str)

        if not len(position) == 3:
            raise Exception("Must be vector len 3")

        if option not in AtomicPositions.options:
            error_str = "CELL_PARAMETERS {0} not valid option".format(self.option)
            raise Exception(error_str)

        self.option = option
        self.atom_positions.append([symbol, np.array(position)])

    def validate(self):
        """
        Validate the class ATOMIC_POSITIONS is properly setup
        (to the best of my knowledge)
        """
        if self.option not in AtomicPositions.options:
            error_str = "ATOMIC_POSITIONS {0} not valid (should never happen)".format(self.option)
            raise Exception(error_str)

    def __str__(self):
        atomicpositions_str = "{0} ({1})\n".format(self.name, self.option)
        for atom_position in self.atom_positions:
            symbol = atom_position[0]
            position_str = " ".join(map(str, atom_position[1]))
            atomicpositions_str += QE_TAB + symbol + " " + position_str + "\n"
        return atomicpositions_str


class KPoints:
    """
    KPoint options (see function for description):
     - tpiba
     - tpiba_b
     - tpiba_c
     - crystal
     - crystal_b
     - crystal_c
     - gamma

    DEFAULT: tpiba
    """

    options = ["tpiba", "tpiba_b", "tpiba_c",
               "crystal", "crystal_b", "crystal_c",
               "gamma", "automatic"]

    def __init__(self):
        self.name = "K_POINTS"
        self.option = None
        self.config = None

    def automatic(self, grid, offset):
        """
        Automatically generated uniform grid of k-points, i.e,
        generates ( nk1, nk2, nk3 ) grid with ( sk1, sk2, sk3 )
        offset.  nk1, nk2, nk3 as in Monkhorst-Pack grids k1, k2, k3
        must be 0 ( no offset ) or 1 ( grid displaced by half a grid
        step in the corresponding direction ) BEWARE: only grids
        having the full symmetry of the crystal work with
        tetrahedra. Some grids with offset may not work.

        FORMAT:
        K_POINTS (automatic)
          nk1 nk2 nk3 sk1 sk2 sk3
        """
        # Validate Input
        if not len(grid) == len(offset) == 3:
            raise Exception("Must be vector len 3")

        self.option = "automatic"
        self.config = [np.array(grid), np.array(offset)]

    def validate(self):
        """
        Validate that class K_POINTS is properly setup
        (to the best of my knowledge)
        """
        if self.option not in KPoints.options:
            error_str = "K_POINT {0} not valid (should never happen)".format(self.option)
            raise Exception(error_str)

    def __str__(self):
        kpoint_str = "{0} ({1})\n".format(self.name, self.option)

        if self.option == "automatic":
            # A little trick to convert list of int to delimited string
            grid_str = " ".join(map(str, self.config[0]))
            offset_str = " ".join(map(str, self.config[1]))
            kpoint_str += QE_TAB + grid_str + " " + offset_str + "\n"
        else:
            error_str = "K_POINT ({0}) not implemeted yet!".format(self.option)
            raise Exception(error_str)

        return kpoint_str


class CellParameters:
    """
    Card: CELL_PARAMETERS
    Defines the vectors for lattice if ibrav = 0. This Card is optional.

    bohr     : latticle vectors in bhor radii
               alat = sqrt(vec1*vec1)

    angstrom : lattice vectors in angstroms
               alat = sqrt(vec1*vec1)

    alat     : lattice vectors in units of the lattice parameter
               either celldm(1) or a.

    DEFAULT: alat
    """
    options = ["alat", "bhor", "angstrom"]

    def __init__(self):
        self.name = "CELL_PARAMETERS"
        self.option = None
        self.lattice_vec = None

    def add_lattice_vec(self, vec1, vec2, vec3, option="alat"):
        """
        Sets the lattice vectors [vec1, vec2, vec3].

        option - specifies the units to use:
         - bohr
         - angstrom
         - alat
        """
        if not len(vec1) == len(vec2) == len(vec3):
            raise Exception("Must be vector len 3")

        if option not in CellParameters.options:
            error_str = "CELL_PARAMETERS {0} not valid option".format(self.option)
            raise Exception(error_str)

        self.option = option
        self.lattice_vec = np.array([vec1, vec2, vec3])

    def validate(self):
        """
        Validate that class K_POINTS is properly setup
        (to the best of my knowledge)
        """
        if self.option not in CellParameters.options:
            error_str = "CELL_PARAMETER {0} not valid (should never happen)".format(self.option)
            raise Exception(error_str)

    def __str__(self):
        cellparameter_str = "{0} ({1})\n".format(self.name, self.option)
        v1_str = " ".join(map(str, self.lattice_vec[0]))
        cellparameter_str += QE_TAB + v1_str + "\n"
        v2_str = " ".join(map(str, self.lattice_vec[1]))
        cellparameter_str += QE_TAB + v2_str + "\n"
        v3_str = " ".join(map(str, self.lattice_vec[2]))
        cellparameter_str += QE_TAB + v3_str + "\n"

        return cellparameter_str

